# Negative binomial at a known size (#253, step 6a), fitted by the engine as
# the binomial model of its records (engine_records()).

negative_binomial_survey <- function(n, field_variance, theta, seed) {
  set.seed(seed)
  coords <- cbind(runif(n, 0, 100), runif(n, 0, 100))
  x <- rnorm(n)
  m <- runif(n, 0.5, 3)
  field <- if (field_variance > 0) {
    covariance <- field_variance * exp(-as.matrix(dist(coords)) / 20)
    drop(t(chol(covariance)) %*% rnorm(n))
  } else {
    numeric(n)
  }
  y <- rnbinom(n, size = theta, mu = m * exp(0.5 + 0.4 * x + field))
  sf::st_as_sf(data.frame(y = y, x = x, m = m, X = coords[, 1],
                          Y = coords[, 2]),
               coords = c("X", "Y"), crs = 32636)
}

test_that("the records are the binomial kernel, with the constant they leave out", {
  theta <- 2.5
  y <- c(0, 3, 12, 1, 7, 40)
  m <- c(1, 2.5, 4, 0.5, 3, 0.2)
  eta <- c(-1, 0.3, 1.2, -0.4, 0.8, 5)
  parameters <- list(theta = theta)
  records <- engine_records("negative_binomial", y, m, eta, parameters)
  expect_identical(records$family, "binomial")
  expect_identical(records$y, y)
  expect_identical(records$units_m, y + theta)
  expect_identical(records$offset, eta + log(m / theta))
  direct <- dnbinom(y, size = theta, mu = m * exp(eta), log = TRUE)
  constant <- log_likelihood_constant("negative_binomial", y, parameters)
  expect_identical(constant, sum(lgamma(y + theta) - lgamma(theta)))
  # The compiled kernel of the engine's weights, and the terms of its mode
  # search, at the records' linear predictors.
  kernel <- make_family("binomial")
  compiled <- kernel$log_likelihood(matrix(records$offset, 1L), y,
                                    records$units_m)
  expect_equal(compiled + constant - sum(lgamma(y + 1)), sum(direct),
               tolerance = 1e-13)
  terms <- kernel$terms(records$offset, y, records$units_m)
  expect_equal(terms$log_likelihood + constant - sum(lgamma(y + 1)),
               sum(direct), tolerance = 1e-13)
  # The score and curvature of each record in eta, as the derivatives of the
  # negative binomial log-likelihood written directly.
  mean <- m * exp(eta)
  expect_equal(terms$score, theta * (y - mean) / (theta + mean),
               tolerance = 1e-13)
  expect_equal(terms$curvature,
               theta * mean * (y + theta) / (theta + mean)^2,
               tolerance = 1e-13)
  # The engine fits the other families' records as they are.
  expect_identical(engine_records("poisson", y, m, eta),
                   list(family = "poisson", y = y, units_m = m, offset = eta))
  expect_identical(log_likelihood_constant("binomial", y), 0)
  # Without exposures or offsets, each is taken as 1 or 0.
  expect_identical(
    engine_records("negative_binomial", y, NULL, NULL, parameters)$offset,
    rep(log(1 / theta), length(y))
  )
})

test_that("the negative binomial family is given by its function only", {
  expect_identical(
    read_family(negative_binomial(theta = 2)),
    list(family = "negative_binomial", invlink = NULL,
         parameters = list(theta = 2))
  )
  expect_identical(
    read_family(negative_binomial(theta = 2), families = engine_families()),
    read_family(negative_binomial(theta = 2))
  )
  expect_input_error(negative_binomial(theta = 0), "theta")
  expect_input_error(negative_binomial(theta = c(1, 2)), "theta")
  expect_input_error(negative_binomial(), "theta")
  expect_input_error(read_family(negative_binomial(2), invlink = exp),
                     "invlink", "log link only")
  # A string has no size, and errors name only the families users can give.
  error <- expect_input_error(read_family("negative_binomial"), "family")
  expect_false(grepl("negative", conditionMessage(error)))
  expect_input_error(
    read_family("negative_binomial", families = engine_families()), "family"
  )
  # The engine sees its records as binomial ones, so the family has no
  # object, and a caller that loses its size fails.
  expect_error(make_family("negative_binomial"),
               class = "riskmap_error_internal")
  expect_error(engine_records("negative_binomial", 1, 1, 0),
               class = "riskmap_error_internal")
  expect_error(engine_records("negative_binomial", 1, 1, 0,
                              list(theta = 2), invlink = exp),
               class = "riskmap_error_internal")
  expect_error(read_family_link("negative_binomial", "log"),
               class = "riskmap_error_internal")
  expect_identical(check_engine_family("negative_binomial")$kernel_family,
                   "binomial")
  link <- make_inverse_link_functions("negative_binomial", NULL, 3L)
  expect_identical(link[c("name", "link")],
                   list(name = "canonical", link = "log"))
  expect_identical(link$inv(c(-1, 0, 2)), exp(c(-1, 0, 2)))
  expect_identical(
    fit_family_argument(list(family = "negative_binomial",
                             family_parameters = list(theta = 3))),
    negative_binomial(theta = 3)
  )
  expect_identical(fit_family_argument(binomial_model), "binomial")
})

test_that("the traits of the negative binomial family give its closed forms", {
  traits <- family_traits("negative_binomial")
  parameters <- list(theta = 1.5)
  means <- c(0.2, 3, 1.4)
  units_m <- c(2, 1, 0.5)
  set.seed(8)
  expected <- rnbinom(3, size = 1.5, mu = means * units_m)
  set.seed(8)
  expect_identical(traits$simulate(means, units_m, parameters), expected)
  expect_error(traits$simulate(means, units_m, NULL),
               class = "riskmap_error_internal")
  expect_identical(traits$fitted_mean(means, units_m), units_m * means)
  # The support reaches the draws, the response and the 0.999 quantile.
  draws <- c(0, 2, 30)
  expect_identical(traits$score_support(draws, 4, 2, means, parameters), 0:30)
  expect_identical(traits$score_support(draws, 41, 2, means, parameters),
                   0:41)
  rates <- c(20, 30, 40)
  quantile <- qnbinom(0.999, size = 1.5, mu = mean(2 * rates))
  expect_gt(quantile, 41)
  expect_identical(traits$score_support(draws, 4, 2, rates, parameters),
                   0:quantile)
  expect_true(traits$validate(c(0, 5), c(0.5, 2)))
  expect_input_error(traits$validate(c(0, 1.5), c(1, 1)), "data",
                     "whole numbers")
  expect_input_error(traits$validate(c(0, 1), c(1, 0)), "denominator",
                     "positive")
  expect_false(traits$upper_limit)
  expect_true(traits$discrete)
  expect_identical(traits$response_quantity, "rate")
  expect_identical(model_title(traits),
                   "Negative binomial geostatistical model")
  expect_null(traits$start_coefficients)
})

test_that("a negative binomial fit is the binomial fit of its records", {
  data <- negative_binomial_survey(40, 0.5, 2, seed = 34)
  fit <- glgpm(y ~ x + gp(nugget = FALSE), data = data,
               family = negative_binomial(theta = 2), denominator = m,
               distance_units = "m", control_mcmc = control_mcmc,
               messages = FALSE, seed = 4)
  data$trials <- data$y + 2
  data$records_offset <- log(data$m / 2)
  binomial_fit <- glgpm(y ~ x + offset(records_offset) + gp(nugget = FALSE),
                        data = data, family = "binomial",
                        denominator = trials, distance_units = "m",
                        control_mcmc = control_mcmc, messages = FALSE,
                        seed = 4)
  expect_identical(fit$estimate, binomial_fit$estimate)
  expect_identical(fit$covariance, binomial_fit$covariance)
  expect_identical(fit$S_samples, binomial_fit$S_samples)
  expect_identical(attr(fit, "starting_values"),
                   attr(binomial_fit, "starting_values"))
  expect_equal(fit$log_lik - binomial_fit$log_lik,
               sum(lgamma(data$y + 2) - lgamma(2)), tolerance = 1e-14)
  # The fit keeps the model's records, link and size.
  expect_identical(fit$family, "negative_binomial")
  expect_identical(fit$family_parameters, list(theta = 2))
  expect_identical(fit$units_m, data$m)
  expect_identical(fit$link_function[c("name", "link")],
                   list(name = "canonical", link = "log"))
  expect_equal(fit$linear_predictors + log(data$m / 2),
               binomial_fit$linear_predictors, tolerance = 1e-12)
  # The fitted values are expected counts, under the log link of the mean.
  eta <- sweep(fit$S_samples, 2L, drop(fit$D %*% fit$estimate$beta), "+")
  expect_equal(fit$fitted_values, data$m * colMeans(exp(eta)),
               tolerance = 1e-10)
  expect_output(print(fit), "Negative binomial geostatistical model")
  expect_output(print(summary(fit)), "canonical log link")
})

test_that("as the field vanishes the fit is the negative binomial regression", {
  data <- negative_binomial_survey(40, 0, 1.5, seed = 21)
  regression <- stats::optim(
    c(0, 0),
    function(beta) {
      -sum(dnbinom(data$y, size = 1.5,
                   mu = data$m * exp(beta[1] + beta[2] * data$x), log = TRUE))
    },
    method = "BFGS", control = list(reltol = 1e-14)
  )
  expect_warning(
    fit <- glgpm(y ~ x + gp(nugget = FALSE), data = data,
                 family = negative_binomial(theta = 1.5), denominator = m,
                 distance_units = "m", control_mcmc = control_mcmc,
                 messages = FALSE, seed = 3),
    "effectively zero"
  )
  expect_identical(attr(fit, "optimizer")$mcml_boundary,
                   "spatial variance (sigma2)")
  expect_equal(unname(fit$estimate$beta), regression$par, tolerance = 1e-6)
  # The reported log-likelihood leaves out -lgamma(y + 1), as the Poisson
  # family's does.
  expect_equal(fit$log_lik - sum(lgamma(data$y + 1)), -regression$value,
               tolerance = 1e-8)
})

test_that("a negative binomial fit predicts, simulates and refits with its size", {
  data <- negative_binomial_survey(30, 0.5, 2, seed = 41)
  fit <- glgpm(y ~ x + gp(nugget = FALSE), data = data,
               family = negative_binomial(theta = 2), denominator = m,
               distance_units = "m", control_mcmc = control_mcmc,
               messages = FALSE, seed = 5)
  # The predicted rate is the mean of exp() of the linear predictor.
  prediction <- predict(fit, data[1:3, ], control_mcmc = control_mcmc,
                        keep_samples = TRUE, messages = FALSE, seed = 1)
  table <- prediction$predictions
  link <- table$fixed + table$offset + prediction$samples$spatial
  expect_equal(table$link_mean, rowMeans(link), tolerance = 1e-12)
  expect_equal(table$response_mean, rowMeans(exp(link)), tolerance = 1e-12)
  model <- simulation_model(fit)
  expect_identical(model$family, "negative_binomial")
  expect_identical(model$family_parameters, list(theta = 2))
  # The responses have the negative binomial variance mu + mu^2 / theta
  # about their means, well above the Poisson variance.
  simulation <- simulate_glgpm(fit, nsim = 300, seed = 6)
  counts <- data$m * simulation$samples$data[, , "mean"]
  squares <- (simulation$samples$data[, , "response"] - counts)^2
  expect_equal(mean(squares / (counts + counts^2 / 2)), 1, tolerance = 0.1)
  expect_gt(mean(squares / counts), 1.5)
  specified <- specify_glgpm(
    y ~ x + gp(nugget = FALSE), data, negative_binomial(theta = 2),
    list(beta = c(0.5, 0.4), sigma2 = 0.5, phi = 20), denominator = m,
    distance_units = "m"
  )
  expect_identical(specified$family_parameters, list(theta = 2))
  arguments <- assessment_refit_args(fit, fit$formula, data, control_mcmc)
  expect_identical(arguments$family, negative_binomial(theta = 2))
  expect_null(arguments$invlink)
})

test_that("the integrated form tabulates records whose trials are not whole", {
  set.seed(61)
  n <- 80
  coords <- cbind(runif(n, 0, 100), runif(n, 0, 100))
  x <- rnorm(n)
  m <- runif(n, 2, 10)
  field <- drop(t(chol(0.4 * exp(-as.matrix(dist(coords)) / 25))) %*% rnorm(n))
  nugget <- rnorm(n, sd = sqrt(0.6))
  y <- rnbinom(n, size = 3.5, mu = m * exp(-0.5 + 0.4 * x + field + nugget))
  data <- sf::st_as_sf(data.frame(y = y, x = x, m = m, X = coords[, 1],
                                  Y = coords[, 2]),
                       coords = c("X", "Y"), crs = 32636)
  fit_in <- function(integrate_noise) {
    glgpm(y ~ x + gp(nugget = TRUE), data = data,
          family = negative_binomial(theta = 3.5), denominator = m,
          distance_units = "m",
          control_mcml = set_control_mcml(integrate_noise = integrate_noise),
          return_samples = FALSE, messages = FALSE, seed = 4)
  }
  sampled <- fit_in(FALSE)
  integrated <- fit_in(TRUE)
  expect_false(attr(sampled, "optimizer")$mcml_integrate_noise)
  expect_identical(attr(integrated, "optimizer")$mcml_integrated_locations,
                   80L)
  # The two forms of the same likelihood agree to their Monte Carlo error.
  expect_equal(unlist(integrated$estimate), unlist(sampled$estimate),
               tolerance = 0.01)
  expect_equal(integrated$log_lik, sampled$log_lik, tolerance = 1e-4)
})
