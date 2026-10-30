test_that("Gaussian weights agree with a direct covariance solve", {
  # Repeated locations with both nested and crossed random-effect groups.
  location <- rep(1:3, c(3, 2, 2))
  covariance <- diag(c(1, 2, 3, 0.4, 0.7))
  for (group in list(c(1, 1, 1, 2, 2, 1, 1), c(1, 2, 1, 2, 1, 1, 2))) {
    incidence <- cbind(diag(3)[location, ], diag(2)[group, ])
    cross_covariance <- covariance %*% t(incidence)
    observation_covariance <- incidence %*% covariance %*% t(incidence) +
      diag(length(location)) * 0.2
    weights <- gaussian_prediction_weights(Matrix::Matrix(incidence, sparse = TRUE),
                                            covariance, cbind(location, group), 0.2)
    expect_equal(weights(cross_covariance),
                 cross_covariance %*% solve(observation_covariance),
                 tolerance = 1e-12)
  }
})

test_that("Gaussian predictions remain finite with tiny measurement error", {
  model <- gaussian_intercept_model
  for (variance in c(1e-10, 1e-12)) {
    model$estimate$sigma2_me <- log(variance)
    for (type in c("marginal", "joint")) {
      expect_warning(prediction <- predict(
        model,
        newdata = gaussian_data,
        type = type,
        control_mcmc = control_mcmc,
        keep_samples = TRUE,
        messages = FALSE
      ), NA)
      expect_true(all(is.finite(prediction$samples$spatial)))
    }
  }
})

test_that("Gaussian predictions use a fixed measurement error variance", {
  # coef() omits sigma2_me when it is fixed, so the fixed value must be
  # taken from the fit. A copy that treats the same value as estimated
  # must give identical predictions.
  fixed_model <- glgpm(y ~ cov + gp(), data = gaussian_data,
                       family = "gaussian", fix_var_me = 0.2,
                       messages = FALSE)
  estimated_model <- fixed_model
  estimated_model$fix_var_me <- NULL
  estimated_model$estimate$sigma2_me <- log(0.2)
  for (type in c("marginal", "joint")) {
    fixed <- predict(fixed_model, newdata = gaussian_data, type = type,
                     control_mcmc = control_mcmc, keep_samples = TRUE,
                     messages = FALSE)
    estimated <- predict(estimated_model, newdata = gaussian_data,
                         type = type, control_mcmc = control_mcmc,
                         keep_samples = TRUE, messages = FALSE)
    expect_equal(fixed$samples$spatial, estimated$samples$spatial,
                 tolerance = 1e-12)
  }
})
