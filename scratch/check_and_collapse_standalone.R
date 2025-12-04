# Unified function to check and collapse both random fields and temporal processes
check_and_collapse <- function(
    tmb_obj,
    report_vals,
    spatial,
    spatiotemporal,
    time_varying_type,
    n_m,
    n_t,
    omit_spatial_intercept,
    collapse = TRUE,                         # Control what to check
    collapse_variance_threshold = 0.01,      # For sigma_O and sigma_E
    collapse_correlation_threshold = 0.01,   # Distance from 0 or 1 for rho
    delta,
    silent,
    sd_report = NULL) {                      # Add sd_report parameter
  
  # Manual testing block - set values here for interactive debugging
  if(FALSE){  
    tmb_obj <- fit8_adj$tmb_obj
    report_vals <- tmb_obj$report()
    spatial = "off"
    spatiotemporal = "off"
    time_varying_type = "ar1"
    n_m = 1
    n_t = 10
    omit_spatial_intercept = FALSE
    collapse = TRUE
    collapse_variance_threshold = 0.01
    collapse_correlation_threshold = 0.01
    delta = FALSE
    silent = FALSE
    sd_report = TMB::sdreport(tmb_obj, getJointPrecision = FALSE)
  }
  
  # Determine what to check based on collapse parameter
  check_variance <- FALSE
  check_correlation <- FALSE
  
  if (is.logical(collapse)) {
    if (collapse) {
      check_variance <- TRUE
      check_correlation <- TRUE
    }
    # If FALSE, both stay FALSE (no checking)
  } else if (collapse == "variance") {
    check_variance <- TRUE
  } else if (collapse == "correlation") {
    check_correlation <- TRUE
  }
  
  # Create sd_report if not provided and we need it for correlation checks
  if (check_correlation && is.null(sd_report)) {
    sd_report <- tryCatch({
      # Suppress TMB output
      capture.output({
        sd_rep <- TMB::sdreport(tmb_obj, getJointPrecision = FALSE)
      }, type = "message")
      sd_rep
    }, error = function(e) {
      if (!silent) {
        cli_inform(c(
          "i" = "Could not compute standard errors for correlation checks",
          ">" = "Skipping correlation-based collapse checks"
        ))
      }
      NULL
    })
  }
  
  do_refit <- FALSE
  spatial_updated <- spatial
  spatiotemporal_updated <- spatiotemporal
  time_varying_type_updated <- time_varying_type
  warning_info <- NULL  # Store warnings to display after fit
  
  # ============================================================================
  # PART 1: VARIANCE COLLAPSE (sigma_O and sigma_E)
  # ============================================================================
  # Only run if check_variance is TRUE
  
  if (check_variance) {
    # Check spatial field variance
    if (any(spatial == "on") && !omit_spatial_intercept) {
      est_sigma_O <- report_vals$sigma_O
      if (length(est_sigma_O) > 0 && any(est_sigma_O < collapse_variance_threshold)) {
        which_sigma <- which(est_sigma_O < collapse_variance_threshold)
        if (!silent) {
          cli_inform(c(
            "!" = "Spatial variance below threshold ({collapse_variance_threshold}) detected",
            "i" = "Affected model(s): {paste(which_sigma, collapse = ', ')}",
            ">" = "Refitting with spatial field(s) disabled"
          ))
        }
        for (i in which_sigma) {
          spatial_updated[i] <- "off"
        }
        do_refit <- TRUE
      }
    }
    
    # Check spatiotemporal field variance
    if (!all(spatiotemporal == "off")) {
      est_sigma_E <- report_vals$sigma_E
      if (length(est_sigma_E) > 0) {
        est_sigma_E_by_model <- numeric(n_m)
        for (m in seq_len(n_m)) {
          idx_start <- (m - 1) * n_t + 1
          idx_end <- m * n_t
          if (idx_end <= length(est_sigma_E)) {
            est_sigma_E_by_model[m] <- min(est_sigma_E[idx_start:idx_end])
          }
        }
        if (any(est_sigma_E_by_model < collapse_variance_threshold, na.rm = TRUE)) {
          which_sigma <- which(est_sigma_E_by_model < collapse_variance_threshold)
          if (!silent) {
            cli_inform(c(
              "!" = "Spatiotemporal variance below threshold ({collapse_variance_threshold}) detected",
              "i" = "Affected model(s): {paste(which_sigma, collapse = ', ')}",
              ">" = "Refitting with spatiotemporal field(s) disabled"
            ))
          }
          for (m in which_sigma) {
            spatiotemporal_updated[m] <- "off"
          }
          do_refit <- TRUE
        }
      }
    }
  }
  
  # ============================================================================
  # PART 2: CORRELATION COLLAPSE (rho parameters for AR1 processes)
  # ============================================================================
  # Only run if check_correlation is TRUE
  
  if (check_correlation) {
    # Helper function to get SE for a parameter
    get_param_se <- function(tmb_obj, param_name, index = 1) {
      if (is.null(tmb_obj$sd_report)) return(NULL)
      idx <- which(names(tmb_obj$sd_report$value) == param_name)
      if (length(idx) == 0) return(NULL)
      if (index > length(idx)) return(NULL)
      
      # Try to get SE from appropriate covariance matrix
      se_val <- tryCatch({
        # First try cov.fixed (for fixed effect parameters like rho for spatiotemporal)
        if (!is.null(tmb_obj$sd_report$cov.fixed) && idx[index] <= nrow(tmb_obj$sd_report$cov.fixed)) {
          sqrt(tmb_obj$sd_report$cov.fixed[idx[index], idx[index]])
        } else if (!is.null(tmb_obj$sd_report$diag.cov.random) && idx[index] <= length(tmb_obj$sd_report$diag.cov.random)) {
          # Try diag.cov.random (for random effect parameters like rho_time)
          sqrt(tmb_obj$sd_report$diag.cov.random[idx[index]])
        } else {
          NULL
        }
      }, error = function(e) NULL)
      
      return(se_val)
    }
    
    # Check spatiotemporal AR1 correlation
    if (any(spatiotemporal_updated == "ar1")) {
      est_rho <- report_vals$rho
      if (length(est_rho) > 0) {
        
        # HIGH CORRELATION: rho near 1 → switch to RW
        # Collapse if EITHER:
        # 1. Point estimate very close to 1
        # 2. OR CI overlaps 1 (uncertain but possibly near 1)
        which_high_rho <- integer(0)
        for (i in which(spatiotemporal_updated == "ar1")) {
          rho_val <- est_rho[i]
          rho_se <- get_param_se(tmb_obj, "rho", i)
          
          collapse_to_rw <- FALSE
          
          # Check 1: Point estimate near 1
          if (rho_val > (1 - collapse_correlation_threshold)) {
            collapse_to_rw <- TRUE
          }
          
          # Check 2: OR CI overlaps 1
          if (!is.null(rho_se)) {
            rho_lower <- rho_val - 1.96 * rho_se
            rho_upper <- rho_val + 1.96 * rho_se
            if (rho_lower < 1 && rho_upper > 1) {
              collapse_to_rw <- TRUE
            }
          }
          
          if (collapse_to_rw) {
            which_high_rho <- c(which_high_rho, i)
          }
        }
        
        if (length(which_high_rho) > 0) {
          if (!silent) {
            cli_inform(c(
              "!" = "Spatiotemporal AR1 correlation near 1 detected",
              "i" = "Affected model(s): {paste(which_high_rho, collapse = ', ')}",
              ">" = "Refitting with spatiotemporal = 'rw' (random walk)"
            ))
          }
          for (i in which_high_rho) {
            spatiotemporal_updated[i] <- "rw"
          }
          do_refit <- TRUE
        }
        
        # LOW CORRELATION: rho near 0 → switch to IID
        # Collapse if EITHER:
        # 1. Point estimate very close to 0
        # 2. OR CI is wide and overlaps 0 (uncertain correlation)
        which_low_rho <- integer(0)
        for (i in which(spatiotemporal_updated == "ar1")) {
          rho_val <- est_rho[i]
          rho_se <- get_param_se(tmb_obj, "rho", i)
          
          collapse_to_iid <- FALSE
          
          # Check 1: Point estimate near 0
          if (abs(rho_val) < collapse_correlation_threshold) {
            collapse_to_iid <- TRUE
          }
          
          # Check 2: OR CI overlaps 0 (uncertain)
          if (!is.null(rho_se)) {
            rho_lower <- rho_val - 1.96 * rho_se
            rho_upper <- rho_val + 1.96 * rho_se
            if (rho_lower < 0 && rho_upper > 0) {
              collapse_to_iid <- TRUE
            }
          }
          
          if (collapse_to_iid) {
            which_low_rho <- c(which_low_rho, i)
          }
        }
        
        if (length(which_low_rho) > 0) {
          if (!silent) {
            cli_inform(c(
              "!" = "Spatiotemporal AR1 correlation near 0 detected",
              "i" = "Affected model(s): {paste(which_low_rho, collapse = ', ')}",
              ">" = "Refitting with spatiotemporal = 'iid' (independent)"
            ))
          }
          for (i in which_low_rho) {
            spatiotemporal_updated[i] <- "iid"
          }
          do_refit <- TRUE
        }
      }
    }
    
    # Check time-varying AR1 correlation
    if (!is.null(time_varying_type) && any(time_varying_type == "ar1")) {
      # Get unscaled rho_time from sd_report
      est_rho_time_unscaled <- NULL
      rho_time_se_unscaled <- NULL
      
      if (!is.null(sd_report)) {
        # Get estimate
        pl <- as.list(sd_report, "Estimate")
        if (!is.null(pl$rho_time_unscaled)) {
          est_rho_time_unscaled <- pl$rho_time_unscaled
          # Handle potential matrix/array structure
          if (is.matrix(est_rho_time_unscaled) || is.array(est_rho_time_unscaled)) {
            est_rho_time_unscaled <- est_rho_time_unscaled[1, 1]
          }
        }
        
        # Get SE
        pls <- as.list(sd_report, "Std. Error")
        if (!is.null(pls$rho_time_unscaled)) {
          rho_time_se_unscaled <- pls$rho_time_unscaled
          # Handle potential matrix/array structure
          if (is.matrix(rho_time_se_unscaled) || is.array(rho_time_se_unscaled)) {
            rho_time_se_unscaled <- rho_time_se_unscaled[1, 1]
          }
        }
      }
      
      if (!is.null(est_rho_time_unscaled)) {
        # Transform to bounded scale [-1, 1]
        bound <- function(x) 2 / (1 + exp(-x)) - 1
        est_rho_time <- bound(est_rho_time_unscaled)
        
        # HIGH CORRELATION: rho_time near 1 → switch to RW
        # Collapse if EITHER:
        # 1. Point estimate very close to 1
        # 2. OR CI overlaps 1 (uncertain but possibly near 1)
        collapse_to_rw <- FALSE
        
        # Check 1: Point estimate near 1
        if (est_rho_time > (1 - collapse_correlation_threshold)) {
          collapse_to_rw <- TRUE
        }
        
        # Check 2: OR CI overlaps 1
        if (!is.null(rho_time_se_unscaled)) {
          rho_time_lower <- bound(est_rho_time_unscaled - 1.96 * rho_time_se_unscaled)
          rho_time_upper <- bound(est_rho_time_unscaled + 1.96 * rho_time_se_unscaled)
          if (rho_time_lower < 1 && rho_time_upper > 1) {
            collapse_to_rw <- TRUE
          }
        }
        
        if (collapse_to_rw) {
          if (!silent) {
            cli_inform(c(
              "!" = "Time-varying AR1 correlation near 1 (rho_time = {round(est_rho_time, 3)}) detected",
              ">" = "Refitting with time_varying_type = 'rw' (random walk)"
            ))
          }
          time_varying_type_updated <- "rw"
          do_refit <- TRUE
        }
        
        # LOW CORRELATION: rho_time near 0 → WARNING ONLY
        # Warn if EITHER:
        # 1. Point estimate very close to 0
        # 2. OR CI is wide and overlaps 0 (uncertain correlation)
        warn_low_rho <- FALSE
        
        # Check 1: Point estimate near 0
        if (abs(est_rho_time) < collapse_correlation_threshold) {
          warn_low_rho <- TRUE
        }
        
        # Check 2: OR CI overlaps 0 (uncertain)
        if (!is.null(rho_time_se_unscaled)) {
          rho_time_lower <- bound(est_rho_time_unscaled - 1.96 * rho_time_se_unscaled)
          rho_time_upper <- bound(est_rho_time_unscaled + 1.96 * rho_time_se_unscaled)
          if (rho_time_lower < 0 && rho_time_upper > 0) {
            warn_low_rho <- TRUE
          }
        }
        
        if (warn_low_rho) {
          # Store warning info to display after fit completes
          warning_info <- list(
            type = "time_varying_low_rho",
            rho_time = est_rho_time,
            message = paste0(
              "Time-varying AR1 correlation near 0 (rho_time = ", round(est_rho_time, 3), ") detected. ",
              "95% CI overlaps zero - time-varying effects show no temporal correlation. ",
              "Consider removing time_varying or using (1 | time_column) in formula instead."
            )
          )
          # NOTE: Do NOT set do_refit = TRUE
          # time_varying_type = "iid" is not allowed - requires manual formula change
        } else {
          warning_info <- NULL
        }
      }
    }
  }
  
  # ============================================================================
  # Prepare arguments for update()
  # ============================================================================
  
  if (delta) {
    spatial_arg <- as.list(spatial_updated)
    spatiotemporal_arg <- as.list(spatiotemporal_updated)
  } else {
    spatial_arg <- spatial_updated[1]
    spatiotemporal_arg <- spatiotemporal_updated[1]
  }
  
  list(
    do_refit = do_refit,
    spatial_arg = spatial_arg,
    spatiotemporal_arg = spatiotemporal_arg,
    time_varying_type_arg = time_varying_type_updated,
    warning_info = warning_info  # Return warning to display later
  )
}