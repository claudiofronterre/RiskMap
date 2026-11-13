# Negative binomial with its size estimated (#253, step 6b), through the
# family block of the layout, with the records of the binomial kernel taken
# at the size of each parameter value.

size_survey <- function(n, field_variance, theta, seed, repeats = 1L) {
  set.seed(seed)
  locations <- cbind(runif(n), runif(n)) * 100
  coords <- locations[rep(seq_len(n), each = repeats), , drop = FALSE]
  x <- rnorm(n * repeats)
  m <- runif(n * repeats, 0.5, 3)
  field <- if (field_variance > 0) {
    covariance <- field_variance * exp(-as.matrix(dist(locations)) / 20)
    drop(t(chol(covariance)) %*% rnorm(n))[rep(seq_len(n), each = repeats)]
  } else {
    numeric(n * repeats)
  }
  y <- rnbinom(n * repeats, size = theta,
               mu = m * exp(0.5 + 0.4 * x + field))
  sf::st_as_sf(data.frame(y = y, x = x, m = m, X = coords[, 1] * 1000,
                          Y = coords[, 2] * 1000),
               coords = c("X", "Y"), crs = 32636)
}

test_that("the engine takes the records and their terms at the estimated size", {
  set.seed(2)
  n <- 12
  coords <- cbind(runif(n), runif(n))
  y <- rnbinom(n, size = 2, mu = 3)
  m <- runif(n, 0.5, 2)
  D <- cbind(1, rnorm(n))
  distances <- as.matrix(dist(coords))
  engine_at <- function(...) {
    make_adaptive_mcml_engine(D, y, m, distances, family = "negative_binomial",
                              proposal = "laplace", ...)
  }
  held <- engine_at(family_parameters = list(theta = 2))
  estimated <- engine_at(estimated_family = "theta")
  expect_identical(estimated$layout$index$theta, 5L)
  expect_identical(estimated$layout$family, "theta")
  parameters <- c(0.3, -0.2, log(0.5), log(0.3))
  draws <- matrix(rnorm(40 * n), 40)
  # The same parameter value, with the size on the scale of log(1 / theta),
  # gives the same weights plus the terms the records leave out.
  difference <- estimated$evaluate(c(parameters, -log(2)), draws) -
    held$evaluate(parameters, draws)
  expect_equal(difference,
               log_likelihood_constant("negative_binomial", y, m,
                                       list(theta = 2)),
               tolerance = 1e-10)
  expect_identical(estimated$records_at(estimated$unpack(c(parameters,
                                                            -log(2))))$units_m,
                   y + 2)
  # The size held at its Poisson limit has a slope, and so is tested.
  expect_identical(estimated$testable_variances(5L), 5L)
  expect_identical(held$testable_variances(5L), integer())
  expect_error(engine_at(estimated_family = "theta", integrate_noise = TRUE),
               class = "riskmap_error_internal")
})

test_that("the slope at the Poisson limit is the derivative in 1 / theta", {
  traits <- family_traits("negative_binomial")
  y <- c(0, 2, 7, 1)
  mean <- c(0.5, 2.5, 4, 1.2)
  # A forward difference of second order from the Poisson limit, at sizes
  # near 10^4, where dnbinom() is accurate.
  alpha <- 1e-4
  log_density <- function(alpha) {
    if (alpha == 0) return(sum(dpois(y, mean, log = TRUE)))
    sum(dnbinom(y, size = 1 / alpha, mu = mean, log = TRUE))
  }
  numerical <- (4 * log_density(alpha) - log_density(2 * alpha) -
                  3 * log_density(0)) / (2 * alpha)
  expect_equal(sum(traits$limit_slope(y, mean)), numerical, tolerance = 1e-4)
})

test_that("as the field vanishes the estimated size is the regression's", {
  data <- size_survey(80, 0, 1.5, seed = 21)
  regression <- stats::optim(
    c(0, 0, 0),
    function(p) {
      -sum(dnbinom(data$y, size = exp(p[3]),
                   mu = data$m * exp(p[1] + p[2] * data$x), log = TRUE))
    },
    method = "BFGS", control = list(reltol = 1e-14)
  )
  expect_warning(
    fit <- glgpm(y ~ x + gp(nugget = FALSE), data = data,
                 family = negative_binomial(), denominator = m,
                 control_mcmc = control_mcmc, messages = FALSE, seed = 3),
    "effectively zero"
  )
  estimates <- coef(fit)
  expect_equal(unname(estimates$beta), regression$par[1:2], tolerance = 1e-5)
  expect_equal(estimates$theta, exp(regression$par[3]), tolerance = 1e-5)
  expect_equal(fit$log_lik + sum(data$y * log(data$m)) -
                 sum(lgamma(data$y + 1)), -regression$value,
               tolerance = 1e-8)
  # The start is the maximum of the same regression.
  expect_equal(negative_binomial_start_size(data$y, data$m,
                                            cbind(1, data$x),
                                            numeric(nrow(data))),
               exp(regression$par[3]), tolerance = 1e-4)
})

test_that("a fit that estimates the size reports, predicts and refits with it", {
  data <- size_survey(100, 0.6, 2, seed = 31)
  fit <- glgpm(y ~ x + gp(nugget = FALSE), data = data,
               family = negative_binomial(), denominator = m,
               control_mcmc = control_mcmc, messages = FALSE, seed = 4)
  theta <- coef(fit)$theta
  expect_identical(estimated_family_parameters(fit), "theta")
  expect_identical(held_family_parameters(fit), list())
  expect_null(fit$family_parameters)
  expect_identical(theta, exp(-fit$estimate$theta))
  expect_identical(fit_family_values(fit), list(theta = theta))
  expect_identical(fit_family_argument(fit), negative_binomial())
  expect_identical(fit_family_argument(fit, hold = TRUE),
                   negative_binomial(theta = theta))
  fit_summary <- summary(fit)
  expect_identical(rownames(fit_summary$family_parameters),
                   "Negative binomial size")
  expect_equal(fit_summary$family_parameters[, "Estimate"], theta,
               ignore_attr = TRUE)
  expect_output(print(fit), "Negative binomial size")
  expect_output(print(fit_summary), "Dispersion")
  # The profile rebuilds the likelihood with the terms of the records in its
  # weights, so its value at the estimate is the fit's.
  problem <- make_profile_problem(fit, 2000L, 1L)
  expect_lt(abs(problem$value(problem$bank_estimate) - fit$log_lik), 0.01)
  # The bootstrap refits the size from the estimate back to the estimate.
  refit <- make_bootstrap_refit(fit)
  expect_equal(refit(fit$y, 4L), unlist(coef(fit)), tolerance = 1e-5)
  # A model that simulates from the fit holds the size at its estimate.
  expect_identical(simulation_model(fit)$family_parameters,
                   list(theta = theta))
  sampler <- parameter_conditioning_sampler(fit)
  set.seed(5)
  draws <- sampler(fit$estimate, 1000L)
  expect_lt(max(abs(colMeans(draws$samples) - colMeans(fit$S_samples)) /
                  apply(fit$S_samples, 2L, stats::sd)), 0.6)
})

test_that("an estimated size takes its start from start_pars", {
  data <- size_survey(30, 0, 2, seed = 7)
  expect_input_error(
    glgpm(y ~ x + gp(nugget = FALSE), data = data,
          family = negative_binomial(), denominator = m,
          start_pars = list(theta = -1), messages = FALSE),
    "start_pars"
  )
  expect_input_error(
    glgpm(y ~ x + gp(nugget = FALSE), data = data,
          family = negative_binomial(theta = 2), denominator = m,
          start_pars = list(theta = 3), messages = FALSE),
    "start_pars"
  )
})

test_that("the size and a nugget need more than one record per location", {
  data <- size_survey(30, 0.5, 2, seed = 8)
  expect_input_error(
    glgpm(y ~ x + gp(nugget = TRUE), data = data,
          family = negative_binomial(), denominator = m, messages = FALSE),
    "formula", "one observation per location"
  )
  # With two records per location they can be told apart, and the fit is
  # sampled, since the tables of the integrated form are for one size.
  data <- size_survey(40, 0.5, 2, seed = 9, repeats = 2L)
  fit <- glgpm(y ~ x + gp(nugget = TRUE), data = data,
               family = negative_binomial(), denominator = m,
               control_mcml = set_control_mcml(integrate_noise = TRUE),
               return_samples = FALSE, messages = FALSE, seed = 2)
  expect_false(attr(fit, "optimizer")$mcml_integrate_noise)
  expect_named(fit$estimate, c("beta", "sigma2", "phi", "nu2", "theta"))
})
