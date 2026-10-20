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
