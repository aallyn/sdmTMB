#' Turn sdmTMB model output into a tidy data frame
#'
#' @param x Output from [sdmTMB()].
#' @param effects A character value. One of `"fixed"` ('fixed' or main-effect
#'   parameters), `"ran_pars"` (standard deviations, spatial range, and other
#'   random effect and dispersion-related terms), or `"ran_vals"` (individual
#'   random intercepts, if included; behaves like `ranef()`).
#' @param conf.int Include a confidence interval?
#' @param conf.level Confidence level for CI.
#' @param exponentiate Whether to exponentiate the fixed-effect coefficient
#'   estimates and confidence intervals.
#' @param model Which model to tidy if a delta model (1 or 2).
#' @param silent Omit any messages?
#' @param ... Extra arguments (not used).
#'
#' @return A data frame
#' @details
#' Follows the conventions of the \pkg{broom} and \pkg{broom.mixed} packages.
#'
#' Currently, `effects = "ran_pars"` also includes dispersion-related terms
#' (e.g., `phi`), which are not actually associated with random effects.
#'
#' Standard errors for spatial variance terms fit in log space (e.g., variance
#' terms, range, or parameters associated with the observation error) are
#' omitted to avoid confusion. Confidence intervals are still available.
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom stats plogis
#' @examples
#' fit <- sdmTMB(density ~ poly(depth_scaled, 2, raw = TRUE),
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie()
#' )
#' tidy(fit)
#' tidy(fit, conf.int = TRUE)
#' tidy(fit, "ran_pars", conf.int = TRUE)
#'
#' pcod_2011$fyear <- as.factor(pcod_2011$year)
#' fit <- sdmTMB(density ~ poly(depth_scaled, 2, raw = TRUE) + (1 | fyear),
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie()
#' )
#' tidy(fit, "ran_vals")

tidy.sdmTMB <- function(x, effects = c("fixed", "ran_pars", "ran_vals"), model = 1,
                 conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE,
                 silent = FALSE, ...) {
  effects <- match.arg(effects)
  assert_that(is.logical(exponentiate))
  assert_that(is.logical(conf.int))
  if (conf.int) {
    assert_that(is.numeric(conf.level),
      conf.level > 0, conf.level < 1,
      length(conf.level) == 1,
      msg = "`conf.level` must be length 1 and between 0 and 1")
  }

  crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  if (exponentiate) trans <- exp else trans <- I

  reinitialize(x)

  delta <- isTRUE(x$family$delta)
  assert_that(is.numeric(model))
  assert_that(length(model) == 1L)
  if (delta) assert_that(model %in% c(1, 2), msg = "`model` must be 1 or 2.")
  if (!delta) assert_that(model == 1, msg = "Only one model: `model` must be 1.")

  se_rep <- as.list(x$sd_report, "Std. Error", report = TRUE)
  est_rep <- as.list(x$sd_report, "Estimate", report = TRUE)
  se <- as.list(x$sd_report, "Std. Error", report = FALSE)
  est <- as.list(x$sd_report, "Estimate", report = FALSE)

  se <- c(se, se_rep)
  est <- c(est, est_rep)
  # cleanup:
  est$epsilon_st <- NULL
  est$zeta_s <- NULL
  est$omega_s <- NULL
  est$ln_H_input <- NULL

  se$epsilon_st <- NULL
  se$zeta_s <- NULL
  se$omega_s <- NULL
  se$ln_H_input <- NULL

  subset_pars <- function(p, model) {
    p$b_j <- if (model == 1) p$b_j else p$b_j2
    p$ln_tau_O <- p$ln_tau_O[model]
    p$ln_tau_Z <- p$ln_tau_Z[model]
    p$ln_tau_E <- p$ln_tau_E[model]
    p$ln_kappa <- as.numeric(p$ln_kappa[,model])
    p$ln_phi <- p$ln_phi[model]
    p$ln_tau_V <- as.numeric(p$ln_tau_V[,model])
    p$ar1_phi <- as.numeric(p$ar1_phi[model])
    p$ln_tau_G <- as.numeric(p$ln_tau_G[,model])
    p$log_sigma_O <- as.numeric(p$log_sigma_O[1,model])
    p$log_sigma_E <- as.numeric(p$log_sigma_E[1,model])
    p$log_sigma_Z <- as.numeric(p$log_sigma_Z[,model])
    p$log_range <- as.numeric(p$log_range[,model])

    p$phi <- p$phi[model]
    p$range <- as.numeric(p$range[,model])
    p$sigma_E <- as.numeric(p$sigma_E[1,model])
    p$sigma_O <- as.numeric(p$sigma_O[1,model])
    p$sigma_Z <- as.numeric(p$sigma_Z[,model])
    p$sigma_G <- as.numeric(p$sigma_G[,model])
    p
  }
  est <- subset_pars(est, model)
  se <- subset_pars(se, model)

  if (x$family$family[[model]] %in% c("binomial", "poisson")) {
    se$ln_phi <- NULL
    est$ln_phi <- NULL
    se$phi <- NULL
    est$phi <- NULL
  }

  ii <- 1

  # grab fixed effects:
  .formula <- x$split_formula[[model]]$form_no_bars
  .formula <- remove_s_and_t2(.formula)
  if (!"mgcv" %in% names(x)) x[["mgcv"]] <- FALSE
  fe_names <- colnames(model.matrix(.formula, x$data))

  b_j <- est$b_j[!fe_names == "offset", drop = TRUE]
  b_j_se <- se$b_j[!fe_names == "offset", drop = TRUE]
  fe_names <- fe_names[!fe_names == "offset"]
  out <- data.frame(term = fe_names, estimate = b_j, std.error = b_j_se, stringsAsFactors = FALSE)

  if (x$tmb_data$threshold_func > 0) {
    if (x$threshold_function == 1L) {
      par_name <- paste0(x$threshold_parameter, c("-slope", "-breakpt"))
    } else {
      par_name <- paste0(x$threshold_parameter, c("-s50", "-s95", "-smax"))
    }
    out <- rbind(
      out,
      data.frame(
        term = par_name, estimate = est$b_threshold[,model,drop=TRUE],
        std.error = se$b_threshold[,model,drop=TRUE], stringsAsFactors = FALSE
      )
    )
  }

  if (conf.int) {
    out$conf.low <- as.numeric(trans(out$estimate - crit * out$std.error))
    out$conf.high <- as.numeric(trans(out$estimate + crit * out$std.error))
  }
  # must wrap in as.numeric() otherwise I() leaves 'AsIs' class that affects emmeans package
  out$estimate <- as.numeric(trans(out$estimate))
  if (exponentiate) out$std.error <- NULL

  out_re <- list()
  log_name <- c("log_range")
  name <- c("range")
  if (!isTRUE(is.na(x$tmb_map$ln_phi))) {
    log_name <- c(log_name, "ln_phi")
    name <- c(name, "phi")
  }
  if (x$tmb_data$include_spatial[model]) {
    log_name <- c(log_name, "log_sigma_O")
    name <- c(name, "sigma_O")
  }
  if (!x$tmb_data$spatial_only[model]) {
    log_name <- c(log_name, "log_sigma_E")
    name <- c(name, "sigma_E")
  }
  if (x$tmb_data$spatial_covariate) {
    log_name <- c(log_name, "log_sigma_Z")
    name <- c(name, "sigma_Z")
  }
  if (x$tmb_data$random_walk) {
    log_name <- c(log_name, "ln_tau_V")
    name <- c(name, "tau_V")
  }
  if (length(est$ln_tau_G) > 0L) {
    log_name <- c(log_name, "ln_tau_G")
    name <- c(name, "sigma_G")
  }

  j <- 0
  if (!"log_range" %in% names(est)) {
    cli_warn("This model was fit with an old version of sdmTMB. Some parameters may not be available to the tidy() method. Re-fit the model with the current version of sdmTMB if you need access to any missing parameters.")
  }

  for (i in name) {
    j <- j + 1
    if (i %in% names(est)) {
      .e <- est[[log_name[j]]]
      .se <- se[[log_name[j]]]
      .e <- if (is.null(.e)) NA else .e
      .se <- if (is.null(.se)) NA else .se

      non_log_name <- gsub("ln_", "", gsub("log_", "", log_name))
      this <- non_log_name[j]
      if (this == "tau_G") this <- "sigma_G"
      if (this == "tau_V") this <- "sigma_V"
      this_se <- as.numeric(se[[this]])
      this_est <- as.numeric(est[[this]])
      if (length(this_est) && !(all(this_se == 0) && all(this_est == 0))) {
        out_re[[i]] <- data.frame(
          term = i, estimate = this_est, std.error = this_se,
          conf.low = exp(.e - crit * .se),
          conf.high = exp(.e + crit * .se),
          stringsAsFactors = FALSE
        )
      }
      ii <- ii + 1
    }
  }
  discard <- unlist(lapply(out_re, function(x) length(x) == 1L)) # e.g. old models and phi
  out_re[discard] <- NULL

  if ("tweedie" %in% x$family$family) {
    out_re$tweedie_p <- data.frame(
      term = "tweedie_p", estimate = plogis(est$thetaf) + 1,
      std.error = se$tweedie_p, stringsAsFactors = FALSE)
    out_re$tweedie_p$conf.low <- plogis(est$thetaf - crit * se$thetaf) + 1
    out_re$tweedie_p$conf.high <- plogis(est$thetaf + crit * se$thetaf) + 1
    ii <- ii + 1
  }

  if ("ar1_phi" %in% names(est)) {
    ar_phi <- est$ar1_phi
    ar_phi_se <- se$ar1_phi
    rho_est <- 2 * stats::plogis(ar_phi) - 1
    rho_lwr <- 2 * stats::plogis(ar_phi - crit * ar_phi_se) - 1
    rho_upr <- 2 * stats::plogis(ar_phi + crit * ar_phi_se) - 1
    out_re[[ii]] <- data.frame(
      term = "rho", estimate = rho_est, std.error = NA,
      conf.low = rho_lwr, conf.high = rho_upr, stringsAsFactors = FALSE
    )
    ii <- ii + 1
  }

  if (all(!x$tmb_data$include_spatial) && all(x$tmb_data$spatial_only)) out_re$range <- NULL

  out_re <- do.call("rbind", out_re)
  row.names(out_re) <- NULL

  if (identical(est$ln_tau_E, 0)) out_re <- out_re[out_re$term != "sigma_E", ]
  if (identical(est$ln_tau_V, 0)) out_re <- out_re[out_re$term != "sigma_V", ]
  if (identical(est$ln_tau_G, 0)) out_re <- out_re[out_re$term != "sigma_G", ]
  if (identical(est$ln_tau_O, 0)) out_re <- out_re[out_re$term != "sigma_O", ]
  if (identical(est$ln_tau_Z, 0)) out_re <- out_re[out_re$term != "sigma_Z", ]
  if (is.na(x$tmb_map$ar1_phi[model])) out_re <- out_re[out_re$term != "rho", ]

  if (!conf.int) {
    out_re[["conf.low"]] <- NULL
    out_re[["conf.high"]] <- NULL
  }

  # random intercepts
  n_re_int <- x$split_formula[[model]]$n_bars
  if (n_re_int == 0 && effects == "ran_vals") {
    cli::cli_abort("effects = 'ran_vals' currently only works with random intercepts (e.g., `+ (1 | g)`).")
  }
  if (n_re_int > 0) {
    out_ranef <- list()
    re_est <- as.list(x$sd_report, "Estimate")$RE
    re_ses <- as.list(x$sd_report, "Std. Error")$RE
    for(jj in 1:n_re_int) {
      level_names <- levels(x$data[[x$split_formula[[model]]$barnames[jj]]])
      n_levels <- length(level_names)
      re_name <- x$split_formula[[model]]$barnames[jj]

      if(jj==1) {
        start_pos <- 1
        end_pos <- n_levels
      } else {
        start_pos <- end_pos + 1
        end_pos <- start_pos + n_levels - 1
      }
      out_ranef[[jj]] <- data.frame(
        term = paste0(re_name,"_",level_names),
        estimate = re_est[start_pos:end_pos,model],
        std.error = re_ses[start_pos:end_pos,model],
        conf.low = re_est[start_pos:end_pos,model] - crit * re_ses[start_pos:end_pos,model],
        conf.high = re_est[start_pos:end_pos,model] + crit * re_ses[start_pos:end_pos,model],
        stringsAsFactors = FALSE
      )
      if (!conf.int) {
        out_ranef[[jj]][["conf.low"]] <- NULL
        out_ranef[[jj]][["conf.high"]] <- NULL
      }
    }
    out_ranef <- do.call("rbind", out_ranef)
    row.names(out_ranef) <- NULL
  }

  out <- unique(out) # range can be duplicated
  out_re <- unique(out_re)

  if (requireNamespace("tibble", quietly = TRUE)) {
    frm <- tibble::as_tibble
  } else {
    frm <- as.data.frame
  }

  if (effects == "fixed") {
    return(frm(out))
  } else if (effects == "ran_vals") {
    return(frm(out_ranef))
  } else if (effects == "ran_pars") {
    return(frm(out_re))
  } else {
    cli_abort("The specified 'effects' type is not available.")
  }
}

#' @importFrom generics tidy
#' @export
generics::tidy
