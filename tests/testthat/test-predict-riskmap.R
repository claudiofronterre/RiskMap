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
    include = "random_effects",
    keep_samples = TRUE,
    messages = FALSE
  )
  expect_equal(dim(conditional$samples$random_effects),
               dim(conditional$samples$spatial))
  expect_true(all(c(
    "random_effects_mean", "random_effects_lower", "random_effects_upper"
  ) %in% names(conditional$data)))

  observation_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2, seed = 43
  )
  observation <- predict(
    gaussian_offset_model,
    include = "nugget",
    control_mcmc = observation_control,
    keep_samples = TRUE,
    messages = FALSE
  )
  expect_equal(dim(observation$samples$nugget),
               dim(observation$samples$spatial))
  expect_true(all(c("nugget_mean", "nugget_lower", "nugget_upper") %in%
                    names(observation$data)))
})

test_that("prediction distinguishes estimated and fixed nugget variances", {
  expect_equal(
    fitted_nugget_variance(gaussian_offset_model),
    coef(gaussian_offset_model)$tau2
  )

  fixed_nugget_model <- gaussian_model
  fixed_nugget_model$fix_tau2 <- 0.25
  fixed_nugget_model$estimate$nu2 <- NULL
  expect_equal(fitted_nugget_variance(fixed_nugget_model), 0.25)

  setup <- list(par_hat = coef(fixed_nugget_model))
  set.seed(430)
  samples <- prediction_nugget_samples(
    fixed_nugget_model, setup, n_locations = 2000, n_samples = 20,
    include = TRUE
  )
  expect_equal(stats::sd(as.numeric(samples)), 0.5, tolerance = 0.01)
})

test_that("predict.RiskMap validates its concise public interface", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10

  expect_error(
    predict(gaussian_model, prediction_data, include = "spatial"),
    "'include' must contain unique values"
  )
  expect_error(
    predict(gaussian_model, prediction_data, summaries = "mode"),
    "'summaries' must contain unique values"
  )
  expect_error(
    predict(gaussian_model, prediction_data, level = 1),
    "strictly between 0 and 1"
  )

  prediction_data$i <- -1L
  expect_error(
    predict(
      gaussian_model, prediction_data, prediction = "joint",
      include = "random_effects", messages = FALSE
    ),
    "unknown levels for random effect"
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
    aggregate = c("mean", "median", "sum", "min", "max"),
    weights = "population",
    control_mcmc = prediction_control,
    seed = 47,
    keep_samples = TRUE,
    messages = FALSE
  )

  expect_s3_class(result, "RiskMap_prediction")
  expect_true(result$aggregated)
  expect_identical(result$prediction, "joint")
  expect_identical(result$coverage, "exact")
  expect_identical(st_geometry(result$data), st_geometry(boundaries))
  expect_identical(result$data$district, boundaries$district)
  expect_equal(nrow(result$data), nrow(boundaries))
  expect_true(all(c(
    "boundary_id", "aggregate", "quantity", "summary", "threshold", "value"
  ) %in% names(result$summary)))
  expect_equal(nrow(result$samples$mean$response), nrow(boundaries))
  expect_true(any(is.finite(result$summary$value)))
})

test_that("polygon aggregation applies weights on the correct scale", {
  values <- matrix(c(1, 3, 2, 6), nrow = 2)
  weights <- c(1, 3)

  expect_equal(
    aggregate_prediction_values(values, "mean", weights),
    matrix(c(2.5, 5), nrow = 1)
  )
  expect_equal(
    aggregate_prediction_values(values, "sum", weights),
    matrix(c(10, 20), nrow = 1)
  )
  expect_equal(
    aggregate_prediction_values(values, "median", weights),
    matrix(c(3, 6), nrow = 1)
  )
  expect_equal(
    aggregate_prediction_values(values, "median", rep(1, 2)),
    matrix(c(2, 4), nrow = 1)
  )
  expect_equal(
    aggregate_prediction_values(values, "min", weights),
    matrix(c(1, 2), nrow = 1)
  )
  expect_equal(
    aggregate_prediction_values(values, "max", weights),
    matrix(c(3, 6), nrow = 1)
  )
})

test_that("RiskMap grids use exact fractional polygon coverage", {
  outer <- st_sf(
    geometry = st_sfc(st_polygon(list(matrix(
      c(0, 0, 2, 0, 2, 2, 0, 2, 0, 0), ncol = 2, byrow = TRUE
    )))),
    crs = 32631
  )
  prediction_grid <- create_grid(outer, spacing = 1, distance_units = "m")
  boundaries <- st_sf(
    area = c("left", "right"),
    geometry = st_sfc(
      st_polygon(list(matrix(
        c(0, 0, 0.75, 0, 0.75, 2, 0, 2, 0, 0),
        ncol = 2, byrow = TRUE
      ))),
      st_polygon(list(matrix(
        c(0.75, 0, 2, 0, 2, 2, 0.75, 2, 0.75, 0),
        ncol = 2, byrow = TRUE
      )))
    ),
    crs = 32631
  )

  assignment <- assign_prediction_locations(prediction_grid, boundaries)

  expect_identical(assignment$method, "exact")
  expect_length(assignment$area_id, 6)
  expect_equal(
    as.numeric(tapply(assignment$coverage, assignment$area_id, sum)),
    c(1.5, 2.5),
    tolerance = 1e-6
  )
  expect_length(unique(assignment$data$grid_id), 4)

  transformed_assignment <- assign_prediction_locations(
    st_transform(prediction_grid, 4326), st_transform(boundaries, 4326)
  )
  expect_identical(transformed_assignment$method, "exact")
  expect_equal(
    as.numeric(tapply(
      transformed_assignment$coverage,
      transformed_assignment$area_id,
      sum
    )),
    c(1.5, 2.5),
    tolerance = 1e-6
  )
})

test_that("ordinary sf points retain point-in-polygon assignment", {
  points <- st_sf(
    value = 1:2,
    geometry = st_sfc(st_point(c(0.25, 0.5)), st_point(c(1.25, 0.5))),
    crs = 32631
  )
  boundaries <- st_sf(
    geometry = st_sfc(st_polygon(list(matrix(
      c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE
    )))),
    crs = 32631
  )

  assignment <- assign_prediction_locations(points, boundaries)

  expect_identical(assignment$method, "point")
  expect_identical(assignment$data$value, 1L)
  expect_identical(assignment$coverage, 1)
  expect_identical(assignment$outside, 2L)
})

test_that("areal exceedance summaries are tidy and threshold-specific", {
  summary_operations <- normalise_prediction_operations(
    c("mean", "exceedance"),
    c("mean", "median", "sd", "cv", "lower", "upper", "exceedance"),
    "summaries"
  )
  summaries <- tidy_prediction_summaries(
    values = matrix(c(0.1, 0.2, 0.4, 0.8), nrow = 1),
    boundary_id = 2L,
    aggregate = "mean",
    quantity = "response",
    summaries = summary_operations,
    thresholds = c(0.15, 0.5),
    level = 0.95
  )

  expect_equal(summaries$summary, c("mean", "exceedance", "exceedance"))
  expect_equal(summaries$threshold, c(NA, 0.15, 0.5))
  expect_equal(summaries$value, c(0.375, 0.75, 0.25))
})

test_that("prediction interval width uses the requested probability level", {
  summary_operations <- normalise_prediction_operations(
    c("lower", "upper", "interval_width"),
    c("mean", "median", "sd", "cv", "lower", "upper", "interval_width",
      "exceedance"),
    "summaries"
  )
  result <- tidy_prediction_summaries(
    matrix(1:100, nrow = 1), 1L, "mean", "response",
    summary_operations, NULL, 0.8
  )

  expect_equal(result$value[3], result$value[2] - result$value[1])
})

test_that("built-in and custom aggregation and summary functions can be mixed", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_data$population <- seq_len(nrow(prediction_data))
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2
  )

  result <- predict(
    gaussian_model,
    newdata = prediction_data,
    boundaries = areal,
    aggregate = list(
      mean = "mean",
      range = function(values, weights) {
        apply(values, 2L, max) - apply(values, 2L, min)
      }
    ),
    summaries = list(mean = "mean", iqr = stats::IQR),
    weights = "population",
    control_mcmc = prediction_control,
    seed = 471,
    messages = FALSE
  )

  expect_equal(sort(unique(result$summary$aggregate)), c("mean", "range"))
  expect_true(all(c("mean", "iqr") %in%
                    result$summary$summary[
                      result$summary$quantity == "response"
                    ]))
  expect_identical(result$aggregate, c("mean", "range"))
  expect_identical(result$summaries, c("mean", "iqr"))
})

test_that("custom prediction operations have strict output contracts", {
  expect_error(
    aggregate_prediction_values(
      matrix(1:6, nrow = 2), function(values, weights) 1, c(1, 1)
    ),
    "one numeric value per predictive draw"
  )
  expect_error(
    normalise_prediction_operations(
      list(function(values, weights) colMeans(values)),
      c("mean", "sum"), "aggregate"
    ),
    "unique, non-empty name"
  )
})

test_that("areal aggregation must be chosen explicitly", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10

  expect_error(
    predict(
      gaussian_model, newdata = prediction_data, boundaries = areal,
      messages = FALSE
    ),
    "'aggregate' must be supplied"
  )
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
    aggregate = "mean",
    control_mcmc = prediction_control,
    seed = 49,
    messages = FALSE
  )

  expect_identical(sampler_calls, 1L)
})

test_that("non-Gaussian boundary predictions reuse shared random effects", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_data$i <- rep(binomial_model$re[[1]], length.out = nrow(prediction_data))

  result <- predict(
    binomial_model,
    newdata = prediction_data,
    boundaries = areal,
    aggregate = "mean",
    include = "random_effects",
    keep_samples = TRUE,
    messages = FALSE
  )

  expect_true(all(c("mean", "lower", "upper") %in%
                    result$summary$summary[
                      result$summary$quantity == "random_effects"
                    ]))
  expect_equal(nrow(result$samples$mean$random_effects), nrow(areal))
})

test_that("Gaussian boundary predictions retain one shared random-effect draw", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  prediction_data$i <- rep(gaussian_model$re[[1]],
                           length.out = nrow(prediction_data))
  prediction_control <- set_control_mcmc(
    n_sim = 40, burnin = 20, thin = 2
  )

  result <- predict(
    gaussian_model,
    newdata = prediction_data,
    boundaries = areal,
    aggregate = "mean",
    include = "random_effects",
    control_mcmc = prediction_control,
    seed = 50,
    keep_samples = TRUE,
    messages = FALSE
  )

  expect_true(all(c("mean", "lower", "upper") %in%
                    result$summary$summary[
                      result$summary$quantity == "random_effects"
                    ]))
  expect_equal(nrow(result$samples$mean$random_effects), nrow(areal))
})

test_that("location predictions report exceedance probabilities", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  result <- predict(
    gaussian_model,
    newdata = prediction_data,
    summaries = c("mean", "exceedance"),
    thresholds = c(0, 1),
    keep_samples = TRUE,
    seed = 501,
    messages = FALSE
  )
  link <- prediction_link_samples(result)

  expect_equal(result$data$response_exceedance_0, rowMeans(link > 0))
  expect_equal(result$data$response_exceedance_1, rowMeans(link > 1))
  expect_equal(result$thresholds, c(0, 1))
})

test_that("derived targets are transformed before summarising and aggregating", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
  location_result <- predict(
    gaussian_model,
    newdata = prediction_data,
    transform = list(square = function(link) link^2),
    summaries = c("mean", "exceedance"),
    thresholds = 1,
    keep_samples = TRUE,
    seed = 502,
    messages = FALSE
  )
  link <- prediction_link_samples(location_result)
  expect_equal(location_result$data$square_mean, rowMeans(link^2))
  expect_equal(location_result$data$square_exceedance_1,
               rowMeans(link^2 > 1))

  area_result <- predict(
    gaussian_model,
    newdata = prediction_data,
    boundaries = areal,
    aggregate = "mean",
    transform = list(square = function(link) link^2),
    summaries = "mean",
    seed = 503,
    messages = FALSE
  )
  expect_true("square" %in% area_result$summary$quantity)
})

test_that("derived targets validate names and returned dimensions", {
  prediction_data <- grid
  prediction_data$cov <- seq_len(nrow(prediction_data)) / 10

  expect_error(
    predict(gaussian_model, prediction_data,
            transform = list(response = function(link) link),
            messages = FALSE),
    "reserved prediction quantities"
  )
  expect_error(
    predict(gaussian_model, prediction_data,
            transform = list(bad = function(link) rowMeans(link)),
            messages = FALSE),
    "same dimensions"
  )
})

test_that("Gaussian block conditioning matches direct multivariate conditioning", {
  model <- gaussian_model
  model$estimate$sigma2_re[] <- log(0.35)
  model$estimate$sigma2_me <- log(0.2)
  pars <- coef(model)
  distance_scale <- attr(model, "distance_scale") %||% 1
  fitted_coords <- model$coords / distance_scale
  fitted_phi <- pars$phi / distance_scale
  spatial_covariance <- pars$sigma2 * matern_correlation(
    pairwise_distances(fitted_coords), phi = fitted_phi,
    kappa = model$kappa, return_sym_matrix = TRUE
  )
  diag(spatial_covariance) <- diag(spatial_covariance) +
    (pars$tau2 %||% model$fix_tau2 %||% 0)

  spatial_design <- matrix(0, nrow(model$data), nrow(spatial_covariance))
  spatial_design[cbind(seq_len(nrow(model$data)), model$ID_coords)] <- 1
  random_design <- matrix(0, nrow(model$data), length(model$re[[1]]))
  random_design[cbind(seq_len(nrow(model$data)), model$ID_re[, 1])] <- 1
  incidence <- cbind(spatial_design, random_design)
  random_covariance <- diag(pars$sigma2_re[[1]], ncol(random_design))
  latent_covariance <- as.matrix(Matrix::bdiag(
    spatial_covariance, random_covariance
  ))
  residual <- model$y -
    as.numeric(model$D %*% pars$beta + model$cov_offset)
  observation_covariance <- incidence %*% latent_covariance %*% t(incidence) +
    diag(pars$sigma2_me, nrow(incidence))
  direct_gain <- latent_covariance %*% t(incidence) %*%
    solve(observation_covariance)
  expected_mean <- as.numeric(direct_gain %*% residual)
  expected_covariance <- latent_covariance -
    direct_gain %*% incidence %*% latent_covariance

  conditioned <- gaussian_conditioning_distribution(
    model, spatial_covariance, pars, residual
  )
  expect_equal(conditioned$mean, expected_mean, tolerance = 1e-9)
  expect_equal(conditioned$covariance, expected_covariance, tolerance = 1e-9)

  prediction_coords <- coordinates_in_units(
    st_geometry(grid[1:2, ]), model$distance_units
  ) / distance_scale
  cross_covariance <- pars$sigma2 * matern_correlation(
    cross_distances(prediction_coords, fitted_coords), phi = fitted_phi,
    kappa = model$kappa
  )
  prediction_covariance <- pars$sigma2 * matern_correlation(
    pairwise_distances(prediction_coords), phi = fitted_phi,
    kappa = model$kappa, return_sym_matrix = TRUE
  )
  spatial_weights <- cross_covariance %*% solve(spatial_covariance)
  n_spatial <- nrow(spatial_covariance)
  spatial_rows <- seq_len(n_spatial)
  random_rows <- n_spatial + seq_len(ncol(random_design))
  block_mean <- c(
    as.numeric(spatial_weights %*% conditioned$mean[spatial_rows]),
    conditioned$mean[random_rows]
  )
  block_covariance <- rbind(
    cbind(
      prediction_covariance - spatial_weights %*% t(cross_covariance) +
        spatial_weights %*%
        conditioned$covariance[spatial_rows, spatial_rows] %*%
        t(spatial_weights),
      spatial_weights %*%
        conditioned$covariance[spatial_rows, random_rows, drop = FALSE]
    ),
    cbind(
      conditioned$covariance[random_rows, spatial_rows, drop = FALSE] %*%
        t(spatial_weights),
      conditioned$covariance[random_rows, random_rows, drop = FALSE]
    )
  )
  target_to_latent <- as.matrix(Matrix::bdiag(
    cross_covariance, random_covariance
  ))
  target_covariance <- as.matrix(Matrix::bdiag(
    prediction_covariance, random_covariance
  ))
  target_to_observation <- target_to_latent %*% t(incidence)
  direct_target_gain <- target_to_observation %*%
    solve(observation_covariance)
  expect_equal(block_mean,
               as.numeric(direct_target_gain %*% residual), tolerance = 1e-9)
  expect_equal(
    block_covariance,
    target_covariance - direct_target_gain %*% t(target_to_observation),
    tolerance = 1e-9
  )
})
