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
      expect_warning(prediction <- setup_prediction(model,
                                                    grid_pred = st_geometry(gaussian_data),
                                                    type = type,
                                                    control_sim = control_mcmc,
                                                    messages = FALSE), NA)
      expect_true(all(is.finite(prediction$S_samples)))
    }
  }
})
