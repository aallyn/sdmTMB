#!/usr/bin/env Rscript
# Standalone Test Script for check_and_collapse()
# Run this with: Rscript standalone_test.R
#
# NOTE: In the actual package, check_and_collapse() is an internal function
# within fit.R. This standalone test sources it separately for testing purposes.

# ============================================================================
# SETUP
# ============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("  Testing check_and_collapse() Function\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("\n")

# Source the function
if (!file.exists(here::here("scratch/check_and_collapse_standalone.R"))) {
  stop("ERROR: check_and_collapse_UNIFIED.R not found in current directory!\n",
       "Please make sure you're running this from the directory containing the function file.")
}

cat("NOTE: In the package, this function is internal within fit.R\n")
cat("      This test sources it separately for standalone testing.\n\n")

source(here::here("scratch/check_and_collapse_standalone.R"))

# Mock cli_inform if cli package not available
if (!requireNamespace("cli", quietly = TRUE)) {
  cli_inform <- function(x) {
    if (is.list(x)) {
      for (msg in x) {
        cat("INFO:", msg, "\n")
      }
    } else {
      cat("INFO:", x, "\n")
    }
  }
  cat("Note: cli package not found, using simple message output\n\n")
} else {
  library(cli)
}

# Test counter
test_num <- 0
passed <- 0
failed <- 0

# Helper function to run a test
run_test <- function(test_name, test_func) {
  test_num <<- test_num + 1
  cat("\n")
  cat("────────────────────────────────────────────────────────────────\n")
  cat(sprintf("Test %d: %s\n", test_num, test_name))
  cat("────────────────────────────────────────────────────────────────\n")
  
  tryCatch({
    test_func()
    passed <<- passed + 1
    cat("✓ PASSED\n")
  }, error = function(e) {
    failed <<- failed + 1
    cat("✗ FAILED:", conditionMessage(e), "\n")
  })
}

# ============================================================================
# TEST 1: Low variance (sigma_O) - should turn spatial off
# ============================================================================

run_test("Low spatial variance (sigma_O < threshold)", function() {
  cat("\nScenario: sigma_O = 0.005 (below 0.01 threshold)\n")
  cat("Expected: Turn spatial field OFF\n\n")
  
  report_vals <- list(
    sigma_O = 0.005,  # Below threshold
    rho = 0.5         # Middle range, shouldn't trigger
  )
  
  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals,
    spatial = "on",
    spatiotemporal = "ar1",
    time_varying_type = NULL,
    n_m = 1,
    n_t = 10,
    omit_spatial_intercept = FALSE,
    collapse_threshold = 0.01,
    rho_threshold_upper = 0.99,
    rho_threshold_lower = 0.01,
    delta = FALSE,
    silent = FALSE
  )
  
  cat("\nResults:\n")
  cat("  do_refit:", result$do_refit, "\n")
  cat("  spatial_arg:", result$spatial_arg, "\n")
  cat("  spatiotemporal_arg:", result$spatiotemporal_arg, "\n")
  
  stopifnot(result$do_refit == TRUE)
  stopifnot(result$spatial_arg == "off")
  stopifnot(result$spatiotemporal_arg == "ar1")  # Should not change
})

# ============================================================================
# TEST 2: High correlation (rho) - should switch AR1 to RW
# ============================================================================

run_test("High spatiotemporal correlation (rho > 0.99)", function() {
  cat("\nScenario: rho = 0.995 (above 0.99 threshold)\n")
  cat("Expected: Switch spatiotemporal from AR1 to RW\n\n")
  
  report_vals <- list(
    sigma_O = 0.5,   # Fine
    rho = 0.995      # Above threshold
  )
  
  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals,
    spatial = "on",
    spatiotemporal = "ar1",
    time_varying_type = NULL,
    n_m = 1,
    n_t = 10,
    omit_spatial_intercept = FALSE,
    collapse_threshold = 0.01,
    rho_threshold_upper = 0.99,
    rho_threshold_lower = 0.01,
    delta = FALSE,
    silent = FALSE
  )
  
  cat("\nResults:\n")
  cat("  do_refit:", result$do_refit, "\n")
  cat("  spatial_arg:", result$spatial_arg, "\n")
  cat("  spatiotemporal_arg:", result$spatiotemporal_arg, "\n")
  
  stopifnot(result$do_refit == TRUE)
  stopifnot(result$spatial_arg == "on")  # Should not change
  stopifnot(result$spatiotemporal_arg == "rw")
})

# ============================================================================
# TEST 3: Low correlation (rho) - should switch AR1 to IID
# ============================================================================

run_test("Low spatiotemporal correlation (|rho| < 0.01)", function() {
  cat("\nScenario: rho = 0.005 (below 0.01 threshold)\n")
  cat("Expected: Switch spatiotemporal from AR1 to IID\n\n")
  
  report_vals <- list(
    sigma_O = 0.5,   # Fine
    rho = 0.005      # Below threshold
  )
  
  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals,
    spatial = "on",
    spatiotemporal = "ar1",
    time_varying_type = NULL,
    n_m = 1,
    n_t = 10,
    omit_spatial_intercept = FALSE,
    collapse_threshold = 0.01,
    rho_threshold_upper = 0.99,
    rho_threshold_lower = 0.01,
    delta = FALSE,
    silent = FALSE
  )
  
  cat("\nResults:\n")
  cat("  do_refit:", result$do_refit, "\n")
  cat("  spatial_arg:", result$spatial_arg, "\n")
  cat("  spatiotemporal_arg:", result$spatiotemporal_arg, "\n")
  
  stopifnot(result$do_refit == TRUE)
  stopifnot(result$spatial_arg == "on")
  stopifnot(result$spatiotemporal_arg == "iid")
})

# ============================================================================
# TEST 4: Both variance AND correlation issues
# ============================================================================

run_test("Both variance AND correlation issues", function() {
  cat("\nScenario: sigma_O = 0.005 AND rho = 0.995\n")
  cat("Expected: Turn spatial OFF and change spatiotemporal AR1 -> ... (depends on sigma_E)\n\n")
  
  report_vals <- list(
    sigma_O = 0.005,  # Below threshold
    sigma_E = rep(0.5, 10),  # Fine
    rho = 0.995       # Above threshold
  )
  
  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals,
    spatial = "on",
    spatiotemporal = "ar1",
    time_varying_type = NULL,
    n_m = 1,
    n_t = 10,
    omit_spatial_intercept = FALSE,
    collapse_threshold = 0.01,
    rho_threshold_upper = 0.99,
    rho_threshold_lower = 0.01,
    delta = FALSE,
    silent = FALSE
  )
  
  cat("\nResults:\n")
  cat("  do_refit:", result$do_refit, "\n")
  cat("  spatial_arg:", result$spatial_arg, "\n")
  cat("  spatiotemporal_arg:", result$spatiotemporal_arg, "\n")
  
  stopifnot(result$do_refit == TRUE)
  stopifnot(result$spatial_arg == "off")      # Variance issue
  stopifnot(result$spatiotemporal_arg == "rw") # Correlation issue
})

# ============================================================================
# TEST 5: Time-varying with high rho_time
# ============================================================================

run_test("High time-varying correlation (rho_time > 0.99)", function() {
  cat("\nScenario: rho_time = 0.996 (above 0.99 threshold)\n")
  cat("Expected: Switch time_varying_type from AR1 to RW\n\n")
  
  report_vals <- list(
    rho_time = matrix(0.996, 1, 1)  # High correlation
  )
  
  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals,
    spatial = "off",
    spatiotemporal = "off",
    time_varying_type = "ar1",
    n_m = 1,
    n_t = 10,
    omit_spatial_intercept = FALSE,
    collapse_threshold = 0.01,
    rho_threshold_upper = 0.99,
    rho_threshold_lower = 0.01,
    delta = FALSE,
    silent = FALSE
  )
  
  cat("\nResults:\n")
  cat("  do_refit:", result$do_refit, "\n")
  cat("  time_varying_type_arg:", result$time_varying_type_arg, "\n")
  
  stopifnot(result$do_refit == TRUE)
  stopifnot(result$time_varying_type_arg == "rw")
})

# ============================================================================
# TEST 6: No issues (should not refit)
# ============================================================================

run_test("No issues detected", function() {
  cat("\nScenario: All parameters in acceptable ranges\n")
  cat("Expected: No changes, do_refit = FALSE\n\n")
  
  report_vals <- list(
    sigma_O = 0.5,
    sigma_E = rep(0.5, 10),
    rho = 0.5,
    rho_time = matrix(0.6, 1, 1)
  )
  
  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals,
    spatial = "on",
    spatiotemporal = "ar1",
    time_varying_type = "ar1",
    n_m = 1,
    n_t = 10,
    omit_spatial_intercept = FALSE,
    collapse_threshold = 0.01,
    rho_threshold_upper = 0.99,
    rho_threshold_lower = 0.01,
    delta = FALSE,
    silent = TRUE
  )
  
  cat("\nResults:\n")
  cat("  do_refit:", result$do_refit, "\n")
  
  stopifnot(result$do_refit == FALSE)
  stopifnot(result$spatial_arg == "on")
  stopifnot(result$spatiotemporal_arg == "ar1")
  stopifnot(result$time_varying_type_arg == "ar1")
})

# ============================================================================
# TEST 7: Delta model with different issues per component
# ============================================================================

run_test("Delta model (different issues per component)", function() {
  cat("\nScenario: Model 1 has low sigma_O, Model 2 has high rho\n")
  cat("Expected: Different changes for each component\n\n")
  
  report_vals <- list(
    sigma_O = c(0.005, 0.5),   # Model 1 low, Model 2 fine
    sigma_E = rep(c(0.5, 0.5), each = 10),
    rho = c(0.5, 0.995)        # Model 1 fine, Model 2 high
  )
  
  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals,
    spatial = c("on", "on"),
    spatiotemporal = c("ar1", "ar1"),
    time_varying_type = NULL,
    n_m = 2,
    n_t = 10,
    omit_spatial_intercept = FALSE,
    collapse_threshold = 0.01,
    rho_threshold_upper = 0.99,
    rho_threshold_lower = 0.01,
    delta = TRUE,
    silent = FALSE
  )
  
  cat("\nResults:\n")
  cat("  do_refit:", result$do_refit, "\n")
  cat("  spatial_arg:", unlist(result$spatial_arg), "\n")
  cat("  spatiotemporal_arg:", unlist(result$spatiotemporal_arg), "\n")
  
  stopifnot(result$do_refit == TRUE)
  stopifnot(is.list(result$spatial_arg))
  stopifnot(result$spatial_arg[[1]] == "off")    # Model 1: low sigma_O
  stopifnot(result$spatial_arg[[2]] == "on")     # Model 2: fine
  stopifnot(result$spatiotemporal_arg[[1]] == "ar1")  # Model 1: fine
  stopifnot(result$spatiotemporal_arg[[2]] == "rw")   # Model 2: high rho
})

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("  Test Results\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("\n")
cat(sprintf("  Total tests: %d\n", test_num))
cat(sprintf("  Passed: %d\n", passed))
cat(sprintf("  Failed: %d\n", failed))
cat("\n")

if (failed == 0) {
  cat("✓ ALL TESTS PASSED!\n")
  cat("\n")
  cat("The check_and_collapse() function is working correctly.\n")
  cat("\n")
  quit(status = 0)
} else {
  cat("✗ SOME TESTS FAILED\n")
  cat("\n")
  cat("Please review the failed tests above.\n")
  cat("\n")
  quit(status = 1)
}