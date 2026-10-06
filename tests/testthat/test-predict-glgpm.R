test_that("predict_glgpm keeps prediction data and summaries aligned", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2, seed = 41
  )

  result <- predict_glgpm(
    gaussian_model,
    newdata = prediction_data,
    control_mcmc = prediction_control,
    keep_samples = TRUE,
    messages = FALSE
  )

  expect_s3_class(result, "RiskMap_prediction")
  expect_s3_class(result$data, "sf")
  expect_identical(st_geometry(result$data), st_geometry(prediction_data))
  expect_equal(nrow(result$data), nrow(prediction_data))
  expect_true(all(c(
    "fixed", "offset", "spatial_mean", "spatial_lower", "spatial_upper",
    "link_mean", "link_lower", "link_upper", "response_mean",
    "response_lower", "response_upper"
  ) %in% names(result$data)))
  expect_equal(result$data$link_mean, result$data$response_mean)
  expect_equal(nrow(result$samples$spatial), nrow(prediction_data))
})

test_that("predict_glgpm derives predictors and offsets from one sf object", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_data$offset <- seq_len(nrow(prediction_data)) / 100
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2, seed = 42
  )

  result <- predict_glgpm(
    gaussian_offset_model,
    newdata = prediction_data,
    control_mcmc = prediction_control,
    messages = FALSE
  )

  expect_equal(result$data$offset, prediction_data$offset)
  expect_equal(result$data$fixed,
               as.numeric(cbind(1, prediction_data$cov) %*%
                            coef(gaussian_offset_model)$beta))
})

test_that("stored final samples are reused when MCMC controls are omitted", {
  prediction_setup <- setup_prediction(binomial_model, messages = FALSE)
  spatial_columns <- seq_len(nrow(binomial_model$coords))

  expect_equal(
    prediction_setup$S_samples,
    t(binomial_model$S_samples[, spatial_columns, drop = FALSE])
  )
})

test_that("newdata NULL returns aligned predictions at observed rows", {
  result <- predict_glgpm(
    binomial_model,
    newdata = NULL,
    keep_samples = TRUE,
    messages = FALSE
  )

  expect_true(result$observed_locations)
  expect_equal(nrow(result$data), nrow(binomial_model$data))
  expect_identical(st_geometry(result$data), st_geometry(binomial_model$data))
  expect_equal(nrow(result$samples$spatial), nrow(binomial_model$data))
  expect_true(all(result$data$response_mean >= 0 &
                    result$data$response_mean <= 1))
})

test_that("optional random-effect and nugget components are explicit", {
  conditional <- predict_glgpm(
    binomial_model,
    prediction = "joint",
    components = c("fixed", "spatial", "random_effects"),
    keep_samples = TRUE,
    messages = FALSE
  )
  expect_equal(dim(conditional$samples$random_effects),
               dim(conditional$samples$spatial))

  observation_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2, seed = 43
  )
  observation <- predict_glgpm(
    gaussian_offset_model,
    components = c("fixed", "spatial", "nugget"),
    control_mcmc = observation_control,
    keep_samples = TRUE,
    messages = FALSE
  )
  expect_equal(dim(observation$samples$nugget),
               dim(observation$samples$spatial))
})

test_that("predict_glgpm validates its concise public interface", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10

  expect_error(
    predict_glgpm(gaussian_model, prediction_data, components = "spatial"),
    "must include both 'fixed' and 'spatial'"
  )
  expect_error(
    predict_glgpm(gaussian_model, prediction_data, summaries = "mode"),
    "'summaries' must contain unique values"
  )
  expect_error(
    predict_glgpm(gaussian_model, prediction_data, level = 1),
    "strictly between 0 and 1"
  )
})

test_that("coefficient of variation is reported only where interpretable", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2, seed = 44
  )

  result <- predict_glgpm(
    gaussian_model,
    prediction_data,
    summaries = c("mean", "cv"),
    control_mcmc = prediction_control,
    messages = FALSE
  )

  expect_true("response_cv" %in% names(result$data))
  expect_false(any(c("spatial_cv", "link_cv") %in% names(result$data)))
})
