devtools::load_all()

# Manual testing
# ============================================================================
# Setup: Create test data and mesh
# ============================================================================
set.seed(123)

predictor_dat <- data.frame(
  X = runif(1000), 
  Y = runif(1000),
  a1 = rnorm(1000), 
  year = sample(1:5, size = 1000, replace = TRUE)
)

mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

# Simulate data with NO spatial or spatiotemporal variation
sim_dat <- sdmTMB_simulate(
  formula = ~ 1 + a1,
  data = predictor_dat,
  time = "year",
  mesh = mesh,
  family = tweedie(),
  range = 2,
  sigma_E = 0,      # No spatiotemporal variation
  phi = 0.2,
  sigma_O = 0,      # No spatial variation
  seed = 42,
  B = c(0.2, -0.4)
)

# ============================================================================
# TEST 1: Default behavior (collapse = TRUE)
# ============================================================================
cat("Fitting model with spatial and spatiotemporal fields...\n")
cat("Expected: Both should collapse to 'off' (low variance)\n\n")

fit1 <- sdmTMB(
  observed ~ a1,
  data = sim_dat, 
  mesh = mesh, 
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  control = sdmTMBcontrol(collapse = TRUE)  # Default
)

fit1$spatial
fit1$spatiotemporal

# ============================================================================
# TEST 2: Only variance checking on spatial/spatio-temporal pieces (collapse = "variance")
# ============================================================================
fit2 <- sdmTMB(
  observed ~ a1,
  data = sim_dat, 
  mesh = mesh, 
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  control = sdmTMBcontrol(collapse = "variance")
)
fit2$spatial
fit2$spatiotemporal

# ============================================================================
# TEST 3: Disabled collapse (collapse = FALSE)
# ============================================================================
fit3 <- sdmTMB(
  observed ~ a1,
  data = sim_dat, 
  mesh = mesh, 
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  control = sdmTMBcontrol(collapse = FALSE)
)
fit3$spatial
fit3$spatiotemporal

# ============================================================================
# TEST 4: Custom variance threshold
# ============================================================================
set.seed(456)
predictor_dat2 <- data.frame(
  X = runif(1000), 
  Y = runif(1000),
  a1 = rnorm(1000), 
  year = sample(1:5, size = 1000, replace = TRUE)
)

sim_dat2 <- sdmTMB_simulate(
  formula = ~ 1 + a1,
  data = predictor_dat2,
  time = "year",
  mesh = mesh,
  family = tweedie(),
  range = 2,
  sigma_E = 0.2,    # Small spatiotemporal variance
  phi = 0.2,
  sigma_O = 0,
  seed = 43,
  B = c(0.2, -0.4)
)

fit4a <- sdmTMB(
  observed ~ a1,
  data = sim_dat2, 
  mesh = mesh, 
  time = "year",
  spatial = "off",
  spatiotemporal = "iid",
  control = sdmTMBcontrol(collapse = TRUE, variance_threshold = 0.001)
)
fit4a$spatial
fit4a$spatiotemporal

# Adjusting threshold, should collapse
fit4b <- sdmTMB(
  observed ~ a1,
  data = sim_dat2, 
  mesh = mesh, 
  time = "year",
  spatial = "off",
  spatiotemporal = "iid",
  control = sdmTMBcontrol(collapse = TRUE, variance_threshold = 0.3)
)
fit4b$spatial
fit4b$spatiotemporal

# ============================================================================
# TEST 5: Backward compatibility with deprecated arguments
# ============================================================================

# This should give a deprecation warning but still work
fit5 <- suppressWarnings(
  sdmTMB(
    observed ~ a1,
    data = sim_dat, 
    mesh = mesh, 
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(collapse_spatial_variance = TRUE)
  )
)
fit5$spatial
fit5$spatio_temporal


# ============================================================================
# TEST 6: Delta model
# ============================================================================

# Add some zeros
sim_dat$observed[sample(1:500, size = 50, replace = FALSE)] <- 0

fit6 <- sdmTMB(
  observed ~ a1,
  data = sim_dat, 
  mesh = mesh, 
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  family = delta_gamma(),
  control = sdmTMBcontrol(collapse = TRUE)
)
fit6$spatial
fit6$spatiotemporal


# ============================================================================
# Setup for time-varying tests
# ============================================================================

set.seed(789)
predictor_dat_tv <- data.frame(
  X = runif(1000), 
  Y = runif(1000),
  a1 = rnorm(1000), 
  year = rep(1:10, each = 100)  # 10 time steps for AR1
)

mesh_tv <- make_mesh(predictor_dat_tv, xy_cols = c("X", "Y"), cutoff = 0.1)

sim_dat_tv <- sdmTMB_simulate(
  formula = ~ 1,           # No fixed effects
  data = predictor_dat_tv,
  time = "year",
  mesh = mesh_tv,
  family = gaussian(),
  range = 0.5,
  time_varying = ~ 1,      # Time-varying intercept
  time_varying_type = "ar1",  # ← ADD THIS (default is "iid")
  rho_time = 0.99,             
  sigma_V = 0.5,           # ← Change to sigma_V (not sigma_E for time-varying)
  sigma_E = 0,             # No spatiotemporal
  sigma_O = 0,             # No spatial
  phi = 0.1,
  seed = 99,
  B = c(1)
)

# ============================================================================
# TEST 7: AR1 with HIGH correlation (rho near 1) → Should collapse to RW
# ============================================================================
fit7_base <- sdmTMB(
  observed ~ 0,
  data = sim_dat_tv,
  mesh = mesh_tv,
  time = "year",
  spatial = "off",
  spatiotemporal = "off",
  time_varying = ~ 1,
  time_varying_type = "ar1",
  control = sdmTMBcontrol(
    collapse = FALSE
  )
)
fit7_base$call['time_varying_type']
pl <- as.list(fit7_base$sd_report, "Estimate")
pls <- as.list(fit7_base$sd_report, "Std. Error")
.rho <- pl$rho_time_unscaled[1,1]
.se <- pls$rho_time_unscaled[1,1]

bound <- function(x) 2 * plogis(x) - 1

bound(.rho)
bound(.rho + 1.96 * .se)
bound(.rho - 1.96 * .se)

# Expected: Should detect high rho and switch to 'rw' (random walk)
fit7_adj <- sdmTMB(
  observed ~ 0,
  data = sim_dat_tv,
  mesh = mesh_tv,
  time = "year",
  spatial = "off",
  spatiotemporal = "off",
  time_varying = ~ 1,
  time_varying_type = "ar1",
  control = sdmTMBcontrol(
    collapse = TRUE,
    rho_threshold_upper = 0.99
  )
)
fit7_adj$call['time_varying_type']  # Check result
fit7_adj

# ============================================================================
# TEST 8: AR1 with LOW correlation (rho near 0) → Should collapse to IID
# ============================================================================
sim_dat_iid<- sdmTMB_simulate(
  formula = ~ 1,           # No fixed effects
  data = predictor_dat_tv,
  time = "year",
  mesh = mesh_tv,
  family = gaussian(),
  range = 0.5,
  time_varying = ~ 1,      # Time-varying intercept
  time_varying_type = "ar1", 
  rho_time = 0,             
  sigma_V = 0.5,           # ← Change to sigma_V (not sigma_E for time-varying)
  sigma_E = 0,             # No spatiotemporal
  sigma_O = 0,             # No spatial
  phi = 0.1,
  seed = 99,
  B = c(1)
)

fit8_base <- sdmTMB(
  observed ~ 1,
  data = sim_dat_iid,
  mesh = mesh_tv,
  time = "year",
  spatial = "off",
  spatiotemporal = "off",
  time_varying = ~ 1,
  time_varying_type = "ar1",  # Fitting AR1 to IID data
  control = sdmTMBcontrol(
    collapse = FALSE
  ),
  silent = FALSE
)

fit8_base$call['time_varying_type']
pl <- as.list(fit8_base$sd_report, "Estimate")
pls <- as.list(fit8_base$sd_report, "Std. Error")
.rho <- pl$rho_time_unscaled[1,1]
.se <- pls$rho_time_unscaled[1,1]

bound <- function(x) 2 * plogis(x) - 1

bound(.rho)
bound(.rho + 1.96 * .se)
bound(.rho - 1.96 * .se)

fit8_adj <- sdmTMB(
  observed ~ 1,
  data = sim_dat_iid,
  mesh = mesh_tv,
  time = "year",
  spatial = "off",
  spatiotemporal = "off",
  time_varying = ~ 1,
  time_varying_type = "ar1",  # Fitting AR1 to IID data
  control = sdmTMBcontrol(
    collapse = TRUE,
    collapse_variance_threshold = 0.01,
    collapse_correlation_threshold = 0.01
  ),
  silent = FALSE
)
fit8_adj$call
pl <- as.list(fit8_adj$sd_report, "Estimate")
pls <- as.list(fit8_adj$sd_report, "Std. Error")
.rho <- pl$rho_time_unscaled[1,1]
.se <- pls$rho_time_unscaled[1,1]

bound <- function(x) 2 * plogis(x) - 1

bound(.rho)
bound(.rho + 1.96 * .se)
bound(.rho - 1.96 * .se)


source(here::here("scratch/check_and_collapse_standalone.R"))
result <- check_and_collapse(
  tmb_obj = fit8_adj$tmb_obj,
  report_vals = fit8_adj$tmb_obj$report(),
  spatial = "off",
  spatiotemporal = "off",
  time_varying_type = "ar1",
  n_m = 1,
  n_t = 10,
  omit_spatial_intercept = FALSE,
  collapse = TRUE,
  collapse_variance_threshold = 0.01,
  collapse_correlation_threshold = 0.01,
  delta = FALSE,
  silent = FALSE
)

tmb_obj <- fit8_adj$tmb_obj
report_vals <- tmb_obj$report()

# What's rho_time in the report?
print(paste("rho_time from report:", report_vals$rho_time))

# What would the checks be?
est_rho_time <- report_vals$rho_time
if (is.matrix(est_rho_time) || is.array(est_rho_time)) {
  est_rho_time <- est_rho_time[1, 1]
}

print(paste("abs(rho_time):", abs(est_rho_time)))
print(paste("< 0.01?", abs(est_rho_time) < 0.01))


# ============================================================================
# TEST 9: Only check correlation (collapse = "correlation")
# ============================================================================
# Using collapse = 'correlation' to ONLY check AR1 correlation
# Expected: Checks time_varying rho but NOT spatial/spatiotemporal variance

# Add spatial field with low variance
sim_dat_tv2 <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat_tv,
  time = "year",
  mesh = mesh_tv,
  family = gaussian(),
  range = 0.5,
  sigma_E = 0,    # Low spatiotemporal - but we won't check this
  sigma_O = 0,    # Low spatial - but we won't check this
  phi = 0.1,
  seed = 101,
  B = c(1)
)

fit9 <- sdmTMB(
  observed ~ 1,
  data = sim_dat_tv2,
  mesh = mesh_tv,
  time = "year",
  spatial = "on",              # Low variance but won't collapse
  spatiotemporal = "iid",      # Low variance but won't collapse
  time_varying = ~ 1,
  time_varying_type = "ar1",   # This WILL be checked
  control = sdmTMBcontrol(collapse = "correlation")  # Only check correlation
)

fit9$spatial            # Should stay 'on' - variance not checked
fit9$spatiotemporal     # Should stay 'iid' - variance not checked
fit9$time_varying_type  # May change based on rho

# ============================================================================
# TEST 10: Both variance AND correlation checking (collapse = TRUE)
# ============================================================================
# Using collapse = TRUE (default) to check EVERYTHING
# Expected: Checks spatial/spatiotemporal variance AND time_varying rho

fit10 <- sdmTMB(
  observed ~ 1,
  data = sim_dat_tv2,
  mesh = mesh_tv,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  time_varying = ~ 1,
  time_varying_type = "ar1",
  control = sdmTMBcontrol(collapse = TRUE)  # Check BOTH
)

fit10$spatial            # Should collapse to 'off' if variance low
fit10$spatiotemporal     # Should collapse to 'off' if variance low
fit10$time_varying_type  # May change based on rho

# ============================================================================
# SUMMARY
# ============================================================================
# VARIANCE CHECKS (spatial/spatiotemporal):
#   - Detects when sigma_O or sigma_E is very small
#   - Collapses 'on' → 'off' or 'iid' → 'off'
#
# CORRELATION CHECKS (time_varying AR1):
#   - Detects when rho is near 1 (high correlation)
#   - Detects when rho is near 0 (no correlation)
#   - Collapses 'ar1' → 'rw' (if rho near 1)
#   - Collapses 'ar1' → 'iid' (if |rho| near 0)
#
# CONTROL OPTIONS:
#   collapse = TRUE:          Check both variance & correlation
#   collapse = 'variance':    Check variance only
#   collapse = 'correlation': Check correlation only
#   collapse = FALSE:         Check nothing

rmarkdown::render(
  input = here::here("vignettes/articles/collapse_vignette.rmd"),
  output_file = "collapse_vignette.html",
  output_dir = here::here("vignettes/articles/")
)

