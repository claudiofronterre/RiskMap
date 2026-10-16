test_that("predict.RiskMap keeps prediction data and summaries aligned", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2, seed = 41
  )

  result <- predict(
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

test_that("predict.RiskMap derives predictors and offsets from one sf object", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_data$offset <- seq_len(nrow(prediction_data)) / 100
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2, seed = 42
  )

  result <- predict(
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
  prediction <- predict(
    binomial_model,
    keep_samples = TRUE,
    messages = FALSE
  )
  spatial_columns <- seq_len(nrow(binomial_model$coords))

  expect_equal(
    prediction$samples$spatial,
    t(binomial_model$S_samples[, spatial_columns, drop = FALSE])
  )
})

test_that("newdata NULL returns aligned predictions at observed rows", {
  result <- predict(
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
  conditional <- predict(
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
  observation <- predict(
    gaussian_offset_model,
    components = c("fixed", "spatial", "nugget"),
    control_mcmc = observation_control,
    keep_samples = TRUE,
    messages = FALSE
  )
  expect_equal(dim(observation$samples$nugget),
               dim(observation$samples$spatial))
})

test_that("predict.RiskMap validates its concise public interface", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10

  expect_error(
    predict(gaussian_model, prediction_data, components = "spatial"),
    "must include both 'fixed' and 'spatial'"
  )
  expect_error(
    predict(gaussian_model, prediction_data, summaries = "mode"),
    "'summaries' must contain unique values"
  )
  expect_error(
    predict(gaussian_model, prediction_data, level = 1),
    "strictly between 0 and 1"
  )
})

test_that("coefficient of variation is reported only where interpretable", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2, seed = 44
  )

  result <- predict(
    gaussian_model,
    prediction_data,
    summaries = c("mean", "cv"),
    control_mcmc = prediction_control,
    messages = FALSE
  )

  expect_true("response_cv" %in% names(result$data))
  expect_false(any(c("spatial_cv", "link_cv") %in% names(result$data)))
})

test_that("prediction transforms coordinates internally and preserves output geometry", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  geographic_data <- st_transform(prediction_data, 4326)
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2
  )

  projected <- predict(
    gaussian_model, prediction_data, control_mcmc = prediction_control,
    seed = 45, messages = FALSE
  )
  geographic <- predict(
    gaussian_model, geographic_data, control_mcmc = prediction_control,
    seed = 45, messages = FALSE
  )

  expect_identical(st_geometry(geographic$data), st_geometry(geographic_data))
  expect_equal(geographic$data$response_mean, projected$data$response_mean)
})

test_that("missing predictors preserve rows and return missing predictions", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_data$cov[2] <- NA_real_
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2
  )

  expect_warning(
    result <- predict(
      gaussian_model, prediction_data, control_mcmc = prediction_control,
      seed = 46, keep_samples = TRUE, messages = FALSE
    ),
    "1 prediction row"
  )

  expect_equal(nrow(result$data), nrow(prediction_data))
  expect_true(is.na(result$data$response_mean[2]))
  expect_true(all(is.na(result$samples$spatial[2, ])))
  expect_true(all(is.finite(result$data$response_mean[-2])))
})

test_that("boundaries trigger blockwise joint prediction and retain polygons", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_data$population <- seq_len(nrow(prediction_data))
  boundaries <- areal
  boundaries$district <- paste0("district_", seq_len(nrow(boundaries)))
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2
  )

  result <- predict(
    gaussian_model,
    newdata = prediction_data,
    boundaries = boundaries,
    weights = "population",
    control_mcmc = prediction_control,
    seed = 47,
    keep_samples = TRUE,
    messages = FALSE
  )

  expect_s3_class(result, "RiskMap_prediction")
  expect_true(result$aggregated)
  expect_identical(result$prediction, "joint")
  expect_identical(st_geometry(result$data), st_geometry(boundaries))
  expect_identical(result$data$district, boundaries$district)
  expect_equal(nrow(result$data), nrow(boundaries))
  expect_equal(nrow(result$samples$response), nrow(boundaries))
  expect_true(any(is.finite(result$data$response_mean)))
})

test_that("prediction seed makes generated samples reproducible", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2
  )

  first <- predict(
    gaussian_model, prediction_data, control_mcmc = prediction_control,
    seed = 48, keep_samples = TRUE, messages = FALSE
  )
  second <- predict(
    gaussian_model, prediction_data, control_mcmc = prediction_control,
    seed = 48, keep_samples = TRUE, messages = FALSE
  )

  expect_equal(first$samples$spatial, second$samples$spatial)
})

test_that("polygon blocks share one non-Gaussian conditioning run", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  model_without_samples <- binomial_model
  model_without_samples$S_samples <- NULL
  original_sampler <- getFromNamespace("laplace_sampling_mcmc", "RiskMap")
  sampler_calls <- 0L
  local_mocked_bindings(
    laplace_sampling_mcmc = function(...) {
      sampler_calls <<- sampler_calls + 1L
      original_sampler(...)
    },
    .package = "RiskMap"
  )
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2
  )

  predict(
    model_without_samples,
    newdata = prediction_data,
    boundaries = areal,
    control_mcmc = prediction_control,
    seed = 49,
    messages = FALSE
  )

  expect_identical(sampler_calls, 1L)
})
