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
  expect_true(set_control_mcml()$laplace_start)
  expect_false(set_control_mcml(laplace_start = FALSE)$laplace_start)
  expect_error(set_control_mcml(laplace_start = NA), "laplace_start")
  expect_error(set_control_mcml(laplace_start = "yes"), "laplace_start")
})
