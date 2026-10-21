test_that("adaptive importance likelihood agrees with one-dimensional quadrature", {
  y <- 2
  trials <- 5
  beta <- -0.4
  sigma2 <- 0.8
  parameters <- c(beta, log(sigma2), log(1))
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 1), y = y, units_m = trials,
    distances = matrix(0), family = "binomial"
  )
  exact <- log(integrate(
    function(spatial_effect) {
      linear_predictor <- beta + spatial_effect
      exp(y * linear_predictor - trials * softplus(linear_predictor)) *
        dnorm(spatial_effect, sd = sqrt(sigma2))
    },
    lower = -Inf, upper = Inf, rel.tol = 1e-12
  )$value)

  set.seed(8128)
  bank <- matrix(rnorm(20000), ncol = 1)
  observed <- engine$evaluate(parameters, bank)

  expect_equal(observed, exact, tolerance = 0.005)
})

test_that("adaptive diagnostics use independent banks and absolute ESS", {
  coordinate <- seq(0, 1, length.out = 30)
  covariate <- seq(-1, 1, length.out = 30)
  trials <- rep(20, 30)
  y <- round(trials * plogis(-0.4 + 0.5 * covariate +
                               0.4 * sin(2 * pi * coordinate)))
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y,
    units_m = trials, distances = as.matrix(dist(coordinate)),
    family = "binomial"
  )
  fit <- fit_adaptive_mcml(
    engine = engine$evaluate,
    start = c(-0.4, 0.5, log(0.8), log(0.4)),
    latent_dimension = 30,
    initial_budget = 1000,
    maximum_budget = 1000,
    maximum_stages = 1,
    min_absolute_ess = 100,
    max_pareto_k = 10,
    seed = 2718,
    lower = c(-10, -10, -8, -8),
    upper = c(10, 10, 4, 4)
  )
  stage <- fit$history[[1L]]

  expect_true(stage$fit_effective_samples >= 100)
  expect_true(stage$validation_effective_samples >= 100)
  expect_true(is.finite(stage$max_parameter_mcse_ratio))
  expect_true(is.finite(stage$score_agreement_statistic))
  expect_identical(fit$accepted, stage$accepted)
})

test_that("boundary variance components and their linked range are detected", {
  estimate <- c(0.1, 0.2, -12, 1.5, -0.5)
  expect_identical(
    adaptive_boundary_parameters(estimate, c(3L, 5L), -9),
    3L
  )
  expect_identical(
    adaptive_boundary_parameters(estimate, c(3L, 5L), -9,
                                 linked_index = list(`3` = 4L)),
    c(3L, 4L)
  )
  expect_identical(
    adaptive_boundary_parameters(c(0.1, 0.2, -1, 1.5), 3L, -9,
                                 linked_index = list(`3` = 4L)),
    integer()
  )
})

test_that("variance parameters are described in plain language", {
  layout <- cov_par_layout(2L, fix_tau2 = TRUE, n_re = 1L)
  expect_identical(
    describe_variance_parameters(c(5L, 6L), layout, re_names = "village"),
    c("nugget variance (tau2)", "random-effect variance for 'village'")
  )
  expect_identical(describe_variance_parameters(integer(), layout),
                   character())
})

test_that("Laplace start candidates include data-based and interior values", {
  layout <- cov_par_layout(2L, fix_tau2 = TRUE)
  distances <- as.matrix(dist(c(0, 1, 2, 3)))
  start <- c(0, 0, log(2), log(5), -15)
  candidates <- laplace_start_candidates(start, layout, distances)
  expect_identical(names(candidates),
                   c("supplied", "data_range", "interior_variances"))
  expect_identical(candidates$supplied, start)
  expect_equal(candidates$data_range[4L],
               log(median(distances[upper.tri(distances)]) / 3))
  expect_equal(candidates$interior_variances[c(3L, 5L)],
               c(log(2), log(0.25)))
})

test_that("Laplace multistart returns the Laplace value at its estimate", {
  coordinate <- seq(0, 1, length.out = 25)
  covariate <- seq(-1, 1, length.out = 25)
  trials <- rep(15, 25)
  y <- round(trials * plogis(-0.3 + 0.4 * covariate +
                               0.5 * sin(2 * pi * coordinate)))
  distances <- as.matrix(dist(coordinate))
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y, units_m = trials,
    distances = distances, family = "binomial"
  )
  lower <- c(-10, -10, -8, -8)
  upper <- c(10, 10, 4, 4)
  candidates <- laplace_start_candidates(c(0, 0, log(0.5), log(0.3)),
                                         engine$layout, distances)
  laplace <- laplace_multistart(engine$evaluate, candidates, 25, lower, upper)
  zero <- matrix(0, nrow = 1, ncol = 25)

  expect_equal(laplace$log_likelihood,
               engine$evaluate(laplace$estimate, zero))
  expect_true(all(laplace$log_likelihood >=
                    vapply(candidates, function(candidate) {
                      engine$evaluate(candidate, zero)
                    }, numeric(1))))
  expect_contains(names(candidates), laplace$candidate)
})

test_that("a zero spatial variance is held at the boundary without escalation", {
  set.seed(31)
  coordinate <- seq(0, 1, length.out = 60)
  covariate <- rnorm(60)
  trials <- rep(25, 60)
  # No spatial effect: the maximum likelihood spatial variance is zero.
  y <- rbinom(60, trials, plogis(-0.5 + 0.6 * covariate))
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y, units_m = trials,
    distances = as.matrix(dist(coordinate)), family = "binomial"
  )
  fit <- fit_adaptive_mcml(
    engine = engine$evaluate, start = c(-0.5, 0.6, -14, log(0.2)),
    latent_dimension = 60, initial_budget = 1000, maximum_budget = 8000,
    maximum_stages = 4, min_absolute_ess = 100, seed = 99,
    lower = c(-10, -10, -20, -8), upper = c(10, 10, 4, 4),
    variance_index = 3L, linked_index = list(`3` = 4L)
  )

  expect_identical(fit$boundary_parameters, c(3L, 4L))
  expect_false(fit$identifiability_failure)
  expect_true(fit$accepted)
  expect_identical(fit$budget, 1000L)
  expect_true(all(is.finite(diag(fit$covariance)[1:2])))
  expect_true(all(is.na(diag(fit$covariance)[3:4])))
})

test_that("set_control_mcml validates laplace_start", {
  expect_false(set_control_mcml()$laplace_start)
  expect_true(set_control_mcml(laplace_start = TRUE)$laplace_start)
  expect_error(set_control_mcml(laplace_start = NA), "laplace_start")
  expect_error(set_control_mcml(laplace_start = "yes"), "laplace_start")
})

test_that("set_control_mcml validates hold_short_range", {
  expect_true(set_control_mcml()$hold_short_range)
  expect_false(set_control_mcml(hold_short_range = FALSE)$hold_short_range)
  expect_error(set_control_mcml(hold_short_range = NA), "hold_short_range")
})

test_that("median nearest-neighbour distance ignores duplicate locations", {
  distances <- as.matrix(dist(c(0, 0, 1, 3)))
  expect_identical(median_nearest_neighbour_distance(distances), 1)
})

test_that("a range shorter than the location spacing is held at the spacing", {
  set.seed(7)
  coordinate <- seq(0, 10, length.out = 60)
  covariate <- rnorm(60)
  trials <- rep(20, 60)
  # Independent location effects give a short unconstrained range (about
  # 0.3 here). A minimum of 1 forces the hold, which is what is tested.
  y <- rbinom(60, trials, plogis(-0.3 + 0.5 * covariate + rnorm(60, sd = 0.8)))
  distances <- as.matrix(dist(coordinate))
  minimum <- log(1)
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y, units_m = trials,
    distances = distances, family = "binomial"
  )
  fit <- fit_adaptive_mcml(
    engine = engine$evaluate, start = c(-0.3, 0.5, log(0.5), log(0.01)),
    latent_dimension = 60, initial_budget = 1000, maximum_budget = 4000,
    maximum_stages = 3, min_absolute_ess = 100, seed = 5,
    lower = c(-10, -10, -20, -20), upper = c(10, 10, 4, 4),
    variance_index = 3L, linked_index = list(`3` = 4L),
    range_index = 4L, minimum_log_range = minimum
  )

  expect_true(fit$range_held)
  expect_equal(fit$estimate[4L], minimum, tolerance = 1e-6)
  expect_contains(fit$boundary_parameters, 4L)
  expect_true(is.na(fit$covariance[4L, 4L]))
  expect_true(all(is.finite(diag(fit$covariance)[1:3])))
})

test_that("the PC penalty is a proper density on the working scale", {
  layout <- cov_par_layout(2L, fix_tau2 = FALSE)
  penalty <- make_pc_penalty(layout, range_threshold = 2,
                             sigma_threshold = 1, sigma_probability = 0.05,
                             range_probability = 0.05)
  at <- function(log_sigma2, log_phi) penalty(c(0, 0, log_sigma2, log_phi))
  density_sigma <- function(a) exp(vapply(a, at, 0, log_phi = 0))
  density_range <- function(b) exp(vapply(b, function(x) at(0, x), 0))
  sigma_total <- integrate(density_sigma, -40, 15)$value
  range_total <- integrate(density_range, -30, 30)$value
  # The penalty is separable, log pi = g(log sigma^2) + h(log phi), so it
  # integrates to one exactly when sigma_total * range_total = exp(g(0) + h(0)).
  expect_equal(sigma_total * range_total, exp(at(0, 0)), tolerance = 1e-6)
  # Tail probabilities match the specification.
  expect_equal(integrate(density_sigma, log(1), 15)$value / sigma_total, 0.05,
               tolerance = 1e-6)
  expect_equal(integrate(density_range, -30, log(2))$value / range_total, 0.05,
               tolerance = 1e-6)
})

test_that("penalising the engine leaves importance weights unchanged", {
  coordinate <- seq(0, 1, length.out = 20)
  covariate <- seq(-1, 1, length.out = 20)
  trials <- rep(10, 20)
  y <- round(trials * plogis(-0.2 + 0.3 * covariate))
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y, units_m = trials,
    distances = as.matrix(dist(coordinate)), family = "binomial"
  )
  penalty <- function(parameters) -sum(parameters^2)
  penalised <- penalise_engine(engine$evaluate, penalty)
  theta <- c(-0.2, 0.3, log(0.5), log(0.3))
  set.seed(3)
  z <- matrix(rnorm(200 * 20), 200)
  expect_identical(penalised(theta, z, log_weights_only = TRUE),
                   engine$evaluate(theta, z, log_weights_only = TRUE))
  expect_equal(penalised(theta, z), engine$evaluate(theta, z) + penalty(theta))
  expect_equal(penalised(theta, z, diagnostics = TRUE)$value,
               engine$evaluate(theta, z) + penalty(theta))
})

test_that("set_control_mcml validates the PC prior options", {
  expect_false(set_control_mcml()$pc_prior)
  control <- set_control_mcml(pc_prior = TRUE, pc_range_threshold = 5)
  expect_true(control$pc_prior)
  expect_identical(control$pc_range_threshold, 5)
  expect_error(set_control_mcml(pc_prior = NA), "pc_prior")
  expect_error(set_control_mcml(pc_sigma_probability = 1),
               "pc_sigma_probability")
})

test_that("glgpm fits with the PC prior and reports the unpenalised likelihood", {
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    control_mcml = set_control_mcml(
      initial_samples = 500, max_samples = 500, max_stages = 1,
      min_effective_samples = 1, max_pareto_k = 100,
      max_mcse_ratio = 0.99, score_probability = 0.999999, pc_prior = TRUE
    ),
    messages = FALSE
  ))
  expect_true(attr(fit, "optimizer")$mcml_pc_prior)
  expect_true(all(is.finite(unlist(fit$estimate))))
  expect_true(is.finite(fit$log_lik))
})
