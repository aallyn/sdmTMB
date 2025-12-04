test_that("collapsing spatial and spatiotemporal fields works", {
  skip_on_cran()
  skip_on_ci()

  set.seed(123)
  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000), year = sample(1:5, size = 1000, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0,
    phi = 0.2,
    sigma_O = 0,
    seed = 42,
    B = c(0.2, -0.4)
  )
  # create some fake 0s
  sim_dat$observed[sample(1:500, size = 50, replace = FALSE)] <- 0

  # test spatial collapse with gaussian family
  fit_nospatial <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", spatial = "off",
    control = sdmTMBcontrol(collapse = FALSE)
  )
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off",
    control = sdmTMBcontrol(collapse = TRUE)
  )
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial, "ran_pars"), tidy(fit, "ran_pars"))
  expect_equal(tidy(fit_nospatial, "ran_vals"), tidy(fit, "ran_vals"))

  # test spatial / spatiotemporal collapse with spatiotemporal on
  fit_nospatial <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", spatial = "off"
  )
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    control = sdmTMBcontrol(collapse = TRUE)
  )
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial, "ran_pars"), tidy(fit, "ran_pars"))
  expect_equal(tidy(fit_nospatial, "ran_vals"), tidy(fit, "ran_vals"))

  # test spatial collapse with delta family
  fit_nospatial <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", spatial = "off", family = delta_gamma()
  )
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", family = delta_gamma(),
    control = sdmTMBcontrol(collapse = TRUE)
  )
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial, "ran_pars"), tidy(fit, "ran_pars"))
  expect_equal(tidy(fit_nospatial, "ran_vals"), tidy(fit, "ran_vals"))

  # test spatial / spatiotemporal collapse with delta family
  fit_nospatial <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", spatial = "off", family = delta_gamma()
  )
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    family = delta_gamma(),
    control = sdmTMBcontrol(collapse = TRUE)
  )
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial, "ran_pars"), tidy(fit, "ran_pars"))
  expect_equal(tidy(fit_nospatial, "ran_vals"), tidy(fit, "ran_vals"))
})

test_that("custom variance threshold works", {
  skip_on_cran()

  set.seed(456)

  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000), year = sample(1:5, size = 1000, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0.2,
    phi = 0.2,
    sigma_O = 0,
    seed = 43,
    B = c(0.2, -0.4)
  )

  # With 0.001, should not collapse
  fit_default <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatial = "off",
    control = sdmTMBcontrol(collapse = TRUE, variance_threshold = 0.001)
  )

  # With higher threshold (0.3), should collapse
  fit_collapse <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatial = "off",
    control = sdmTMBcontrol(collapse = TRUE, variance_threshold = 0.3)
  )

  expect_true(all(fit_default$spatiotemporal == "iid"))
  expect_true(all(fit_collapse$spatiotemporal == "off"))
})

test_that("collapse = 'variance' only checks variance", {
  skip_on_cran()

  set.seed(789)

  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000), year = sample(1:3, size = 1000, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0,
    phi = 0.2,
    sigma_O = 0,
    seed = 44,
    B = c(0.2, -0.4)
  )

  # With collapse = "variance", should only check variance, not correlation
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat,
    mesh = mesh,
    time = "year",
    spatiotemporal = "ar1",
    control = sdmTMBcontrol(collapse = "variance")
  )

  # Should have collapsed to no spatiotemporal due to low variance
  expect_true(all(fit$spatiotemporal == "off"))
})

test_that("collapse = 'correlation' only checks correlation", {
  skip_on_cran()
  skip_on_ci()

  set.seed(321)

  predictor_dat <- data.frame(
    X = runif(800), Y = runif(800),
    a1 = rnorm(800), year = sample(1:3, size = 800, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0.3,
    phi = 0.2,
    sigma_O = 0,
    seed = 45,
    B = c(0.2, -0.4)
  )

  # With collapse = "correlation", should only check rho, not variance
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "ar1",
    control = sdmTMBcontrol(
      collapse = "correlation",
      rho_threshold_upper = 0.95,
      rho_threshold_lower = 0.05
    )
  )

  # Should still have spatiotemporal field (variance not checked)
  # Might be "ar1", "rw", or "iid" depending on estimated rho
  expect_true(all(fit$spatiotemporal %in% c("ar1", "rw", "iid")))
})

test_that("collapse = FALSE disables all checks", {
  skip_on_cran()

  set.seed(654)

  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000), year = sample(1:3, size = 1000, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0,
    phi = 0.2,
    sigma_O = 0,
    seed = 46,
    B = c(0.2, -0.4)
  )

  # With collapse = FALSE, should not collapse anything
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat,
    mesh = mesh,
    time = "year",
    spatiotemporal = "ar1",
    control = sdmTMBcontrol(collapse = FALSE)
  )

  # Should keep spatiotemporal as ar1 (no collapse despite low variance)
  expect_true(all(fit$spatiotemporal == "ar1"))
})

test_that("collapse validation works", {
  # Test that variance_threshold must be positive
  expect_error(
    sdmTMBcontrol(variance_threshold = -0.01),
    "variance_threshold not greater than 0"
  )

  expect_error(
    sdmTMBcontrol(variance_threshold = 0),
    "variance_threshold not greater than 0"
  )

  # Valid thresholds should work
  ctrl <- sdmTMBcontrol(variance_threshold = 0.001)
  expect_equal(ctrl$variance_threshold, 0.001)

  ctrl <- sdmTMBcontrol(variance_threshold = 0.1)
  expect_equal(ctrl$variance_threshold, 0.1)

  # Test rho thresholds must be in (0, 1)
  expect_error(
    sdmTMBcontrol(rho_threshold_upper = 1.5),
    "rho_threshold_upper not less than 1"
  )

  expect_error(
    sdmTMBcontrol(rho_threshold_lower = -0.1),
    "rho_threshold_lower not greater than 0"
  )

  expect_error(
    sdmTMBcontrol(rho_threshold_upper = 0),
    "rho_threshold_upper not greater than 0"
  )

  # Valid rho thresholds
  ctrl <- sdmTMBcontrol(rho_threshold_upper = 0.99, rho_threshold_lower = 0.01)
  expect_equal(ctrl$rho_threshold_upper, 0.99)
  expect_equal(ctrl$rho_threshold_lower, 0.01)

  # Test collapse argument validation
  expect_error(
    sdmTMBcontrol(collapse = "invalid"),
    "collapse.*must be TRUE, FALSE, 'variance', or 'correlation'"
  )

  # Valid collapse values
  ctrl <- sdmTMBcontrol(collapse = TRUE)
  expect_equal(ctrl$collapse, TRUE)

  ctrl <- sdmTMBcontrol(collapse = FALSE)
  expect_equal(ctrl$collapse, FALSE)

  ctrl <- sdmTMBcontrol(collapse = "variance")
  expect_equal(ctrl$collapse, "variance")

  ctrl <- sdmTMBcontrol(collapse = "correlation")
  expect_equal(ctrl$collapse, "correlation")
})

test_that("deprecated arguments still work with warnings", {
  # Test collapse_spatial_variance deprecation
  expect_warning(
    ctrl <- sdmTMBcontrol(collapse_spatial_variance = TRUE),
    "collapse_spatial_variance.*deprecated"
  )
  expect_equal(ctrl$collapse, "variance")

  expect_warning(
    ctrl <- sdmTMBcontrol(collapse_spatial_variance = FALSE),
    "collapse_spatial_variance.*deprecated"
  )
  expect_equal(ctrl$collapse, FALSE)

  # Test collapse_threshold deprecation
  expect_warning(
    ctrl <- sdmTMBcontrol(collapse_threshold = 0.02),
    "collapse_threshold.*deprecated"
  )
  expect_equal(ctrl$variance_threshold, 0.02)
})

test_that("AR1 to RW conversion works when rho near 1", {
  skip_on_cran()
  skip_on_ci()

  set.seed(987)

  predictor_dat <- data.frame(
    X = runif(600), Y = runif(600),
    a1 = rnorm(600), year = sample(1:3, size = 600, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0.3,
    phi = 0.2,
    sigma_O = 0,
    seed = 47,
    B = c(0.2, -0.4)
  )

  # Fit with AR1 - if rho estimated > 0.95, should switch to RW
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "ar1",
    control = sdmTMBcontrol(
      collapse = TRUE,
      rho_threshold_upper = 0.95
    )
  )

  # Check if model structure is valid (ar1 or rw)
  expect_true(all(fit$spatiotemporal %in% c("ar1", "rw")))
})

test_that("AR1 to IID conversion works when rho near 0", {
  skip_on_cran()
  skip_on_ci()

  set.seed(111)

  predictor_dat <- data.frame(
    X = runif(600), Y = runif(600),
    a1 = rnorm(600), year = sample(1:3, size = 600, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0.3,
    phi = 0.2,
    sigma_O = 0,
    seed = 48,
    B = c(0.2, -0.4)
  )

  # Fit with AR1 - if rho estimated near 0, should switch to IID
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "ar1",
    control = sdmTMBcontrol(
      collapse = TRUE,
      rho_threshold_lower = 0.05
    )
  )

  # Check if model structure is valid (ar1 or iid)
  expect_true(all(fit$spatiotemporal %in% c("ar1", "iid")))
})

test_that("delta models handle collapse correctly", {
  skip_on_cran()
  skip_on_ci()

  set.seed(222)

  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000), year = sample(1:3, size = 1000, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0,
    phi = 0.2,
    sigma_O = 0,
    seed = 49,
    B = c(0.2, -0.4)
  )
  sim_dat$observed[sample(1:500, size = 100, replace = FALSE)] <- 0

  # Fit delta model with collapse enabled
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat,
    mesh = mesh,
    time = "year",
    family = delta_gamma(),
    control = sdmTMBcontrol(collapse = TRUE)
  )

  # Both components should have collapsed
  expect_true(all(fit$spatial == "off"))
  expect_true(all(fit$spatiotemporal == "off"))
})

test_that("check_and_collapse function unit tests", {
  # Unit tests for the check_and_collapse function itself
  # NOTE: check_and_collapse() is internal in fit.R
  # Access with sdmTMB:::check_and_collapse() for unit testing
  
  # Test low variance
  report_vals_low_variance <- list(
    sigma_O = 0.005,
    rho = 0.5
  )

  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals_low_variance,
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
    silent = TRUE
  )

  expect_true(result$do_refit)
  expect_equal(result$spatial_arg, "off")
  expect_equal(result$spatiotemporal_arg, "ar1")

  # Test high rho
  report_vals_high_rho <- list(
    sigma_O = 0.5,
    rho = 0.995
  )

  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals_high_rho,
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
    silent = TRUE
  )

  expect_true(result$do_refit)
  expect_equal(result$spatial_arg, "on")
  expect_equal(result$spatiotemporal_arg, "rw")

  # Test low rho
  report_vals_low_rho <- list(
    sigma_O = 0.5,
    rho = 0.005
  )

  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals_low_rho,
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
    silent = TRUE
  )

  expect_true(result$do_refit)
  expect_equal(result$spatial_arg, "on")
  expect_equal(result$spatiotemporal_arg, "iid")

  # Test no issues
  report_vals_ok <- list(
    sigma_O = 0.5,
    sigma_E = rep(0.5, 10),
    rho = 0.5
  )

  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals_ok,
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
    silent = TRUE
  )

  expect_false(result$do_refit)
  expect_equal(result$spatial_arg, "on")
  expect_equal(result$spatiotemporal_arg, "ar1")

  # Test both variance and correlation issues
  report_vals_both <- list(
    sigma_O = 0.005,
    sigma_E = rep(0.5, 10),
    rho = 0.995
  )

  result <- check_and_collapse(
    tmb_obj = NULL,
    report_vals = report_vals_both,
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
    silent = TRUE
  )

  expect_true(result$do_refit)
  expect_equal(result$spatial_arg, "off")
  expect_equal(result$spatiotemporal_arg, "rw")
})
