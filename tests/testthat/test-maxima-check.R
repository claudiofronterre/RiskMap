test_that("the Gauss-Hermite rule integrates polynomials exactly", {
  rule <- gauss_hermite_rule(60L)
  expect_equal(sum(rule$w), sqrt(pi), tolerance = 1e-12)
  expect_equal(sum(rule$w * rule$x^2), sqrt(pi) / 2, tolerance = 1e-12)
  expect_equal(sum(rule$w * rule$x^4), 3 * sqrt(pi) / 4, tolerance = 1e-12)
})

# Log-likelihood with independent location effects by one-dimensional
# numerical integration per location, on the engine's scale.
integrate_independent_sites <- function(D, y, units_m, location, offset,
                                        family, beta, log_variance) {
  variance <- exp(log_variance)
  sum(vapply(unique(location), function(site) {
    rows <- which(location == site)
    log_density <- function(effects) vapply(effects, function(effect) {
      eta <- drop(D[rows, , drop = FALSE] %*% beta) + offset[rows] + effect
      sum(if (family == "binomial") {
        y[rows] * eta - units_m[rows] * log1p(exp(eta))
      } else {
        y[rows] * eta - units_m[rows] * exp(eta)
      }) - effect^2 / (2 * variance)
    }, numeric(1))
    top <- stats::optimize(log_density, c(-10, 10), maximum = TRUE)$objective
    log(stats::integrate(function(effects) exp(log_density(effects) - top),
                         -Inf, Inf, rel.tol = 1e-12)$value) + top -
      log(2 * pi * variance) / 2
  }, numeric(1)))
}

test_that("the independent-site likelihood matches numerical integration", {
  set.seed(3)
  location <- c(1:8, 2L, 5L, 5L)
  D <- cbind(1, stats::rnorm(length(location)))
  units_m <- rep(c(10, 25), length.out = length(location))
  offset <- stats::rnorm(length(location), sd = 0.1)
  eta <- -1 + 0.4 * D[, 2] + stats::rnorm(8)[location] + offset
  responses <- list(
    binomial = stats::rbinom(length(location), units_m, stats::plogis(eta)),
    poisson = stats::rpois(length(location), units_m * exp(eta))
  )
  for (family in names(responses)) {
    y <- responses[[family]]
    sites <- make_independent_site_likelihood(D, y, units_m, location, offset,
                                              family)
    beta <- c(-0.8, 0.3)
    at <- sites$log_likelihood(beta, log(0.7))
    expect_equal(at$value,
                 integrate_independent_sites(D, y, units_m, location, offset,
                                             family, beta, log(0.7)),
                 tolerance = 1e-10)
    expect_equal(at$gradient, numDeriv::grad(function(parameters) {
      sites$log_likelihood(parameters[1:2], parameters[3])$value
    }, c(beta, log(0.7))), tolerance = 1e-6)
    fit <- sites$maximise(beta)
    expect_lt(max(abs(sites$log_likelihood(fit$beta,
                                           fit$log_variance)$gradient)), 1e-3)
    expect_lt(fit$quadrature_error, 1e-6)
  }
})

test_that("the Jensen bound is below the likelihood and tight at zero scale", {
  y <- c(3, 7)
  units_m <- c(10, 12)
  beta <- -0.2
  variance <- 0.8
  sites <- make_independent_site_likelihood(cbind(c(1, 1)), y, units_m, 1:2,
                                            c(0, 0), "binomial")
  at <- sites$log_likelihood(beta, log(variance))
  independent <- list(log_likelihood = at$value,
                      log_variance = log(variance), mean = at$mean,
                      variance = at$variance)
  # The exact log-likelihood of two correlated effects by tensor quadrature.
  exact <- function(correlation) {
    root <- t(chol(variance * matrix(c(1, correlation, correlation, 1), 2)))
    rule <- gauss_hermite_rule(80L)
    grid <- expand.grid(first = seq_along(rule$x), second = seq_along(rule$x))
    effects <- sqrt(2) * cbind(rule$x[grid$first], rule$x[grid$second]) %*%
      t(root)
    eta <- effects + beta
    log_likelihood <- y[1] * eta[, 1] - units_m[1] * log1p(exp(eta[, 1])) +
      y[2] * eta[, 2] - units_m[2] * log1p(exp(eta[, 2]))
    log(sum(rule$w[grid$first] * rule$w[grid$second] *
              exp(log_likelihood)) / pi)
  }
  bound <- function(correlation) {
    compute_jensen_bound(independent,
                         variance * matrix(c(1, correlation, correlation, 1), 2))
  }
  expect_equal(bound(1e-10), at$value, tolerance = 1e-10)
  for (correlation in c(0.1, 0.5, 0.9)) {
    expect_lt(bound(correlation), exact(correlation))
  }
  expect_true(is.na(compute_jensen_bound(independent, diag(c(1, -1)))))
})

test_that("the three-point step finds the best log variance", {
  value_at <- function(point) -(point[1] - 1.2)^2
  step <- step_log_variance(value_at, c(0, 5), value_at(c(0, 5)), 1L, 0.5, 2,
                            c(-10, -10), c(10, 10))
  expect_equal(step$point, c(1.2, 5))
  expect_equal(step$value, 0)
  # A step never returns a worse point than the one it starts from.
  flat <- step_log_variance(function(point) if (point[1] == 0) 1 else 0,
                            c(0, 5), 1, 1L, 0.5, 2, c(-10, -10), c(10, 10))
  expect_identical(flat$point, c(0, 5))
})

# A stand-in engine and screen for the check: the objective has a maximum at
# a long scale (log phi 2) and a short-scale part `short()`, the screen plugs
# in coefficient and log variance 0, and the independent-site log-likelihood
# is `limit`.
make_toy_check <- function(short, limit = -50, best_log_variance = function(t) 0) {
  objective <- function(parameters, z) {
    value <- -parameters[1]^2 -
      (parameters[2] - best_log_variance(parameters[3]))^2 +
      log(exp(-(parameters[3] - 2)^2) + short(parameters[3]))
    if (is.finite(value)) value else NaN
  }
  engine <- list(latent_dimension = 1L,
                 evaluate = function(parameters, z, log_weights_only) {
                   rep(objective(parameters, z), nrow(z))
                 })
  screen <- list(
    beta_index = 1L, sigma2_index = 2L, phi_index = 3L,
    log_phi = c(-6, -5, -4, -3), margin = 0.25, window = 2, step = 0.5,
    step_limit = 2, alarm_multiplier = 2.58, comparison_tolerance = 1e-4,
    independent = function(beta) {
      list(beta = 0, log_variance = 0, log_likelihood = limit,
           quadrature_error = 0)
    },
    plug_in = function(estimate, fit, log_phi) {
      c(fit$beta, fit$log_variance, log_phi)
    },
    bound = function(point, fit) c(value = limit - 1, error = 0),
    at_limit = function(log_phi) log_phi <= -9.99,
    announce = function(...) NULL
  )
  limits <- list(variance_index = integer(), linked_index = list(),
                 boundary_log_variance = -9, tolerance = 0.1,
                 dispersed = list(sigma2_index = 2L, phi_index = 3L,
                                  log_sigma2 = log(c(0.25, 4, 0.25, 4)),
                                  log_phi = log(c(1, 1, 50, 50)),
                                  labels = c("a", "b", "c", "d"),
                                  trigger = log(0.5), screen = screen))
  lower <- rep(-10, 3)
  upper <- rep(10, 3)
  fitted <- optimise_on_bank(objective, c(0, 0, 1.8), NULL, lower, upper,
                             TRUE, limits)
  list(run = function() {
    check_other_maxima(engine, objective, NULL, fitted, NULL, lower, upper,
                       limits, 100L, 1L)
  }, fitted = fitted)
}

test_that("the short screen finds a higher maximum at a short scale", {
  toy <- make_toy_check(function(t) 3 * exp(-(t + 5)^2))
  expect_equal(toy$fitted$fit$par[3], 2, tolerance = 1e-3)
  checked <- toy$run()
  expect_identical(checked$screen$gate, "opened")
  expect_identical(checked$check$start, c("fitted", "short-scale screen"))
  expect_true(checked$separate)
  expect_equal(checked$optimised$fit$par[3], -5, tolerance = 1e-3)
  expect_false(checked$numerically_independent)
  expect_false(checked$screen$alarm$fired)
})

test_that("the short screen stays closed when short scales fit badly", {
  checked <- make_toy_check(function(t) 1e-8 * exp(-(t + 5)^2))$run()
  expect_identical(checked$screen$gate, "closed")
  expect_null(checked$check)
  expect_true(all(is.na(checked$screen$screen$stepped)))
  expect_equal(checked$optimised$fit$par[3], 2, tolerance = 1e-3)
})

test_that("the variance step opens a gate that the plug-in leaves closed", {
  # At short scales the best log variance is 1.2, not the plug-in's 0, which
  # leaves the plug-in about 0.44 below the fit and the step reaches it.
  toy <- make_toy_check(function(t) exp(1) * exp(-(t + 5)^2),
                        best_log_variance = function(t) 1.2 * (t < 0))
  checked <- toy$run()
  table <- checked$screen$screen
  expect_lt(max(table$plug_in, na.rm = TRUE), -toy$fitted$fit$objective - 0.25)
  expect_true(any(is.finite(table$stepped)))
  expect_identical(checked$screen$gate, "opened")
  expect_equal(checked$optimised$fit$par[3], -5, tolerance = 1e-3)
})

test_that("a failed evaluation in the short screen opens the gate", {
  checked <- make_toy_check(function(t) if (t < -5.5) NA else 0)$run()
  expect_true(checked$screen$failed_evaluation)
  expect_identical(checked$screen$gate, "opened")
  expect_true(all(is.finite(checked$optimised$fit$par)))
})

test_that("an exact alarm forces a search and persists if it finds nothing", {
  checked <- make_toy_check(function(t) 1e-8 * exp(-(t + 5)^2),
                            limit = 5)$run()
  expect_identical(checked$screen$gate, "closed")
  expect_true(checked$screen$alarm$fired)
  expect_identical(checked$screen$alarm$log_phi, -Inf)
  expect_true("exact alarm" %in% checked$check$start)
  expect_true(checked$screen$alarm$persists)
})

test_that("a short scale not clearly above the limit is kept with a warning", {
  # The likelihood rises towards log(3) as the scale goes to zero, the exact
  # limit is 0.001 above that, and the optimiser stops on the flat approach
  # before reaching the limit.
  checked <- make_toy_check(function(t) 3 * stats::plogis(-3 * (t + 4)),
                            limit = log(3) + 0.001)$run()
  expect_identical(checked$screen$gate, "opened")
  expect_lt(checked$optimised$fit$par[3], -6)
  expect_false(checked$numerically_independent)
  expect_false(is.null(checked$screen$near_limit))
  # An unresolved comparison with the limit does not suppress the alarm.
  expect_true(checked$screen$alarm$persists)
  # A higher maximum at a short positive scale is neither.
  interior <- make_toy_check(function(t) 3 * exp(-(t + 7)^2),
                             limit = 0.5)$run()
  expect_equal(interior$optimised$fit$par[3], -7, tolerance = 1e-3)
  expect_false(interior$numerically_independent)
  expect_null(interior$screen$near_limit)
})

test_that("a fit with negligible spatial correlation is numerically independent", {
  # The likelihood keeps rising as the scale goes to zero, and the refit runs
  # to the lower bound of the log scale.
  checked <- make_toy_check(function(t) exp(-0.5 * t), limit = 6)$run()
  expect_equal(checked$optimised$fit$par[3], -10)
  expect_true(checked$numerically_independent)
  expect_null(checked$screen$near_limit)
})

test_that("a numerically independent stage has no information or certificate", {
  engine <- list(latent_dimension = 1L,
                 proposal_type = function() "laplace",
                 proposal_failures = function() 0L,
                 evaluate = function(...) rep(0, 1))
  kept <- list(fit = list(par = c(0, 0, -10), objective = 1, convergence = 0L),
               boundary = integer(), restarted = FALSE,
               restart_improved = FALSE)
  local_mocked_bindings(
    optimise_on_bank = function(...) kept,
    check_other_maxima = function(engine, objective, penalty, optimised,
                                  ...) {
      list(optimised = optimised, check = NULL, separate = FALSE,
           ambiguous = FALSE, screen = list(gate = "opened"),
           numerically_independent = TRUE)
    },
    compute_observed_information = function(...) stop("information computed")
  )
  stage <- run_certified_stage(engine, function(...) -1, NULL, c(0, 0, 0),
                               10L, 10L, 1:4, rep(-10, 3), rep(10, 3), TRUE,
                               list(dispersed = list()))
  expect_identical(stage$outcome, "numerically independent")
  expect_false(stage$certified)
  expect_identical(stage$maxima_screen$gate, "opened")
})

test_that("the outcome of the screen is reported to the user", {
  layout <- cov_par_layout(1L, fix_tau2 = FALSE)
  persisting <- list(maxima_screen = list(alarm = list(persists = TRUE)),
                     maxima_check = NULL, better_maximum = FALSE,
                     ambiguous_maxima = FALSE)
  expect_warning(report_maxima_check(persisting, layout, 1),
                 "not the maximum likelihood estimate")
  zero <- list(outcome = "numerically independent", certified = FALSE,
               maxima_screen = list(independent = list(log_variance = log(2))),
               boundary_parameters = integer(), history = list(1))
  expect_warning(report_mcml_outcome(zero, character(), layout, FALSE, 0.1,
                                     FALSE),
                 "numerically one of independent.*variance 2.*unresolved")
  near <- list(near_limit = list(gap = 0.0011, se = 0.0005, draws = 5000L,
                                 log_phi = log(0.023)))
  # The warning quotes the scale of the returned estimate, not that of the
  # fit compared with the limit.
  expect_warning(report_near_limit(near, scale = 0.2271),
                 paste0("phi = 0.227\\).*not clearly distinguishable.*",
                        "difference 0.0011.*profile\\(fit, parameters = \"phi\"\\)"))
  expect_silent(report_near_limit(list(), scale = 0.2271))
})

test_that("the screen is recorded for fits with and without random effects", {
  screened <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
    control_mcml = fixture_control_mcml, return_samples = FALSE,
    messages = FALSE
  ))
  record <- attr(screened, "optimizer")$mcml_maxima_screen
  expect_true(is.finite(record$independent$log_likelihood))
  expect_named(record$independent$beta, c("(Intercept)", "cov"))
  expect_contains(c("opened", "closed", "not run"), record$gate)
  expect_true(is.finite(record$alarm$statistic))
  unavailable <- attr(binomial_model, "optimizer")$mcml_maxima_screen
  expect_identical(unavailable$gate, "unavailable")
  expect_identical(unavailable$reason, "unstructured random effects")
})

test_that("the records of the screen are restored to user units", {
  history <- list(list(
    estimate = c(0, 0, 1, log(2)),
    maxima_screen = list(
      screen = data.frame(log_phi = log(c(0.1, 0.2))),
      alarm = list(log_phi = log(0.1)),
      near_limit = list(log_phi = log(0.05)),
      independent = list(beta = c(1, 2), log_variance = 0)
    )
  ))
  transform <- matrix(c(1, 0, -1, 2), 2, dimnames = list(NULL, c("a", "b")))
  restored <- restore_mcml_history_scale(history, transform,
                                         distance_scale = 10)
  screen <- restored[[1]]$maxima_screen
  expect_equal(screen$screen$log_phi, log(c(1, 2)))
  expect_equal(screen$alarm$log_phi, log(1))
  expect_equal(screen$near_limit$log_phi, log(0.5))
  expect_equal(screen$independent$beta, c(a = -1, b = 4))
})

test_that("an exact value equal to the fit to numerical precision is no alarm", {
  # The fit is at the limit itself, as with a variance at zero, and the
  # comparison has no Monte Carlo error.
  toy <- make_toy_check(function(t) 1e-8 * exp(-(t + 5)^2),
                        limit = log(1 + 1e-8 * exp(-49)) + 2e-9)
  expect_false(toy$run()$screen$alarm$fired)
})

test_that("an ambiguous choice is judged again on a larger bank", {
  # Two separated points whose likelihoods differ by 0.004, the second with
  # its own noise of sd 0.05 per draw, so that 50 draws cannot tell them
  # apart. The second bank is the one recorded.
  engine <- list(latent_dimension = 1L,
                 evaluate = function(parameters, z, log_weights_only) {
                   0.004 * parameters[1] + parameters[1] * z[, 1] * 0.05
                 })
  information <- diag(100, 2)
  once <- choose_maximum(engine, NULL, list(c(0, 0), c(1, 0)), c(-10, -10),
                         information, 50L, 1L, 0.1)
  twice <- choose_maximum(engine, NULL, list(c(0, 0), c(1, 0)), c(-10, -10),
                          information, 50L, 1L, 0.1, enlarged_seed = 2L)
  expect_true(once$ambiguous)
  expect_identical(unique(once$table$independent_draws), 50L)
  expect_identical(unique(twice$table$independent_draws), 200L)
})

test_that("a failed evaluation in the variance step opens the gate", {
  screen <- list(sigma2_index = 2L, margin = 0.25, window = 2, step = 0.5,
                 step_limit = 2, log_phi = c(-2, -1),
                 plug_in = function(estimate, fit, log_phi) c(0, 0, log_phi))
  # The plug-in values are finite and within the window, and every
  # evaluation of the variance step fails.
  objective <- function(parameters, z) if (parameters[2] == 0) -1 else NA
  screened <- screen_short_scales(objective, NULL, c(0, 0, 2), 0, list(),
                                  screen, rep(-10, 3), rep(10, 3))
  expect_true(screened$failed)
  expect_true(screened$opened)
})

test_that("the Jensen bound carries its own quadrature error", {
  y <- c(3, 7, 0, 12)
  units_m <- c(10, 12, 8, 20)
  sites <- make_independent_site_likelihood(cbind(rep(1, 4)), y, units_m, 1:4,
                                            numeric(4), "binomial")
  fit <- sites$maximise(0)
  expect_named(fit$coarse, c("log_likelihood", "mean", "variance"))
  engine <- list(layout = cov_par_layout(1L, fix_tau2 = FALSE),
                 unpack = function(parameters) list(sigma2 = 1, nugget_ratio = 0),
                 spatial_covariance = function(point) {
                   exp(point[2]) * exp(-as.matrix(dist(1:4)) / exp(point[3]))
                 })
  screen <- make_short_scale_screen(engine, sites, as.matrix(dist(1:4)))
  bound <- screen$bound(c(fit$beta, fit$log_variance, log(0.5)), fit)
  expect_named(bound, c("value", "error"))
  expect_gte(bound[["error"]], 0)
  expect_lt(bound[["error"]], 1e-6)
})
