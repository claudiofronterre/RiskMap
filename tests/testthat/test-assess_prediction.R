test_that("assess_prediction produces errors", {

  expect_error(
    assess_prediction("not list"),
    "'object' must be a list of fitted models of class 'RiskMap'"
  )

  expect_error(
    assess_prediction(list(a = 1)),
    "'object' must be a list of fitted models of class 'RiskMap'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "random"),
    "'method' must be either 'cluster' or 'regularized'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "cluster"),
    "when 'method' is 'cluster' you must supply 'fold'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "cluster",
                      fold = 0.1),
    "'fold' must be a single positive integer"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized"),
    "when 'method' is 'regularized' you must supply 'min_dist'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1),
    "when 'method' is 'regularized' you must supply 'n_size'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1,
                      iter = 2.1),
    "'iter' must be a single positive integer"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1,
                      control_sim = "not sim"),
    "'control_sim' must come from 'set_control_mcmc"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = -1,
                      n_size = 1
    ),
    "'min_dist' must be a single positive"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1.1
    ),
    "'n_size' must be a single positive integer"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1,
                      keep_par_fixed = "not logical"
    ),
    "'keep_par_fixed' must be either TRUE or FALSE"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1,
                      control_sim = "not mcmc"
    ),
    "'control_sim' must come from 'set_control_mcmc"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1,
                      plot_fold = "not true"
    ),
    "'plot_fold' must be either TRUE or FALSE"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1,
                      messages = "not true"
    ),
    "'messages' must be either TRUE or FALSE"
  )

  different_rows <- gaussian_model
  different_rows$data_sf <- different_rows$data_sf[1:9,]

  expect_error(
    assess_prediction(list(gaussian_model,
                           different_rows),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1
    ),
    "All models in 'object' supplied must have the same number of observations"
  )

  different_order <- gaussian_model
  different_order$data_sf <- different_order$data_sf[c(6:10, 1:5),]

  expect_error(
    assess_prediction(list(gaussian_model,
                           different_order),
                      method = "regularized",
                      min_dist = 1,
                      n_size = 1
    ),
    "All models in 'object' must have data in the same row order and geometry"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      user_split = matrix(1:2, ncol = 1)),
    "'user_split' matrix must have the same number of rows as the data in the model"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      user_split = matrix(1:20, ncol = 2)),
    "'user_split' matrix must have a number of columns equal to 'iter'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      user_split = list(1, 2)),
    "'user_split' list must have the same length as 'iter'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      user_split = "not list or matrix"),
    "'user_split' must be a matrix or a list"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      user_split = list(sample(10))),
    "The length of values in 'user_split' to create the test set must be less than the number of rows in the data"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      user_split = list(rep(0,5))),
    "The values in 'user_split' must be row indices of the data"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      user_split = list(c(1.1, 2, 3.3))),
    "The values in 'user_split' must be row indices of the data"
  )

})

make_assess_prediction_fit <- function(data, covariate) {
  coords <- st_coordinates(data)
  fit <- list(
    formula = as.formula(paste("y ~", covariate)),
    data_sf = data,
    family = "gaussian",
    estimate = c(0, 0, 0, 0),
    D = matrix(1, nrow = nrow(data), ncol = 2),
    re = list(),
    cov_offset = NULL,
    fix_tau2 = 0,
    fix_var_me = 0,
    sst = FALSE,
    units_m = rep(1, nrow(data)),
    y = data$y,
    ID_coords = seq_len(nrow(data)),
    coords = coords,
    crs = 4326,
    scale_to_km = FALSE,
    call = list(den = quote(units_m)),
    model_id = covariate
  )
  class(fit) <- "RiskMap"
  fit
}

test_that("assess_prediction uses each model's own data_sf for held-out predictors", {
  geom <- st_sfc(
    st_point(c(0, 0)),
    st_point(c(1, 1)),
    st_point(c(2, 2)),
    crs = 4326
  )
  data_x1 <- st_sf(y = c(1, 2, 3), x1 = c(10, 20, 30), geometry = geom)
  data_x2 <- st_sf(y = c(1, 2, 3), x2 = c(100, 200, 300), geometry = geom)

  seen_predictors <- list()
  local_mocked_bindings(
    setup_prediction = function(object, grid_pred, predictors, ...) {
      seen_predictors[[object$model_id]] <<- names(predictors)
      list(predictors = predictors)
    },
    predict_grid_target = function(object, ...) {
      list(lp_samples = matrix(rep(1, nrow(object$predictors) * 2),
                               nrow = nrow(object$predictors)))
    },
    .package = "RiskMap"
  )

  out <- assess_prediction(
    list(model_x1 = make_assess_prediction_fit(data_x1, "x1"),
         model_x2 = make_assess_prediction_fit(data_x2, "x2")),
    user_split = matrix(c(0, 1, 1), ncol = 1),
    plot_fold = FALSE,
    messages = FALSE,
    which_metric = "CRPS"
  )

  expect_s3_class(out, "RiskMap_cross_validation")

  expect_setequal(names(out), c("test_set", "model"))

  expect_true("x1" %in% seen_predictors$x1)
  expect_false("x2" %in% seen_predictors$x1)
  expect_true("x2" %in% seen_predictors$x2)
  expect_false("x1" %in% seen_predictors$x2)
})

test_that("assess_prediction requires aligned model data", {
  geom_a <- st_sfc(st_point(c(0, 0)), st_point(c(1, 1)), crs = 4326)
  geom_b <- st_sfc(st_point(c(0, 0)), st_point(c(2, 2)), crs = 4326)
  data_a <- st_sf(y = c(1, 2), x1 = c(10, 20), geometry = geom_a)
  data_b <- st_sf(y = c(1, 2), x2 = c(100, 200), geometry = geom_b)

  expect_error(
    assess_prediction(
      list(model_x1 = make_assess_prediction_fit(data_a, "x1"),
           model_x2 = make_assess_prediction_fit(data_b, "x2")),
      user_split = matrix(c(0, 1), ncol = 1),
      plot_fold = FALSE,
      messages = FALSE,
      which_metric = "CRPS"
    ),
    "same row order and geometry"
  )
})

test_that("assess_prediction re-encodes random effects after subsetting", {
  user_split <- matrix(c(1, 1, rep(0, nrow(gaussian_data) - 2)), ncol = 1)

  expect_no_warning(
    out <- assess_prediction(
      list(model = gaussian_model),
      user_split = user_split,
      control_sim = control_mcmc,
      plot_fold = FALSE,
      messages = FALSE,
      which_metric = "CRPS"
    )
  )

  expect_s3_class(out, "RiskMap_cross_validation")
})

test_that("assess_prediction splits test data correctly", {

  n_folds <- 2

  result <- assess_prediction(
    list(intercept_only = gaussian_intercept_model,
         with_covariate = gaussian_model),
    method = "cluster",
    fold = n_folds,
    messages = FALSE)

  expect_length(result$test_set, n_folds)
  expect_equal(sum(unlist(lapply(result$test_set, nrow))), n)

  combined <- do.call(rbind, result$test_set)
  expect_true(all(!duplicated(combined)))

  n_folds <- 3

  result <- assess_prediction(
    list(intercept_only = gaussian_intercept_model,
         with_covariate = gaussian_model),
    method = "cluster",
    fold = n_folds,
    messages = FALSE)

  expect_length(result$test_set, n_folds)
  expect_equal(sum(unlist(lapply(result$test_set, nrow))), n)

  combined <- do.call(rbind, result$test_set)
  expect_true(all(!duplicated(combined)))

  n_size <- 4

  result <- assess_prediction(
    list(intercept_only = gaussian_intercept_model,
         with_covariate = gaussian_model),
    method = "regularized",
    n_size = n_size,
    min_dist = 1,
    messages = FALSE)

  expect_length(result$test_set, 1)
  expect_equal(nrow(result$test_set[[1]]), n_size)


  result <- assess_prediction(
    list(gaussian_model),
    user_split = matrix(
      sample(c(rep(1, n/2), rep(0, n/2))),
      ncol = 1),
    messages = FALSE)

  expect_length(result$test_set, 1)
  expect_equal(nrow(result$test_set[[1]]), n/2)

})

test_that("assess_prediction can refit correctly for all model families", {

  result <- assess_prediction(
    list(gaussian_model),
    method = "regularized",
    min_dist = 1,
    n_size = 1,
    keep_par_fixed = FALSE,
    messages = FALSE)

  expect_setequal(names(result), c("test_set", "model"))

  result <- assess_prediction(
    list(binomial_model),
    method = "regularized",
    min_dist = 1,
    n_size = 1,
    keep_par_fixed = FALSE,
    control_sim = control_mcmc,
    messages = FALSE)

  expect_setequal(names(result), c("test_set", "model"))

  result <- assess_prediction(
    list(poisson_model),
    method = "regularized",
    min_dist = 1,
    n_size = 1,
    keep_par_fixed = FALSE,
    control_sim = control_mcmc,
    messages = FALSE)

  expect_setequal(names(result), c("test_set", "model"))

})


test_that("AnPIT area computes trapezoidal absolute distance", {
  u <- seq(0, 1, length.out = 1001)

  expect_equal(.anpit_area(u, u), 0)
  expect_equal(.anpit_area(rep(0, length(u)), u), 0.5)
  expect_equal(.anpit_area(u^2, u), 1 / 6, tolerance = 1e-6)
})

test_that("assess_prediction reports AnPIT area as a scalar score", {
  geom <- st_sfc(
    st_point(c(0, 0)),
    st_point(c(1, 1)),
    st_point(c(2, 2)),
    crs = 4326
  )
  data <- st_sf(y = c(0, 0, 0), x1 = c(10, 20, 30), geometry = geom)

  testthat::local_mocked_bindings(
    setup_prediction = function(object, grid_pred, predictors, ...) {
      list(predictors = predictors)
    },
    predict_grid_target = function(object, ...) {
      n <- nrow(object$predictors)
      list(lp_samples = matrix(rep(c(-1, 0, 1), each = n), nrow = n))
    },
    .package = "RiskMap"
  )

  out <- assess_prediction(
    list(model_x1 = make_assess_prediction_fit(data, "x1")),
    user_split = matrix(c(0, 1, 1), ncol = 1),
    plot_fold = FALSE,
    messages = FALSE,
    which_metric = "AnPIT"
  )

  expect_named(out$model$model_x1$score, "AnPIT_area")
  expect_length(out$model$model_x1$score$AnPIT_area, 1)
  expect_type(out$model$model_x1$score$AnPIT_area[[1]], "double")
  expect_true(is.finite(out$model$model_x1$score$AnPIT_area[[1]]))
  expect_true(out$model$model_x1$score$AnPIT_area[[1]] >= 0)
  expect_true(out$model$model_x1$score$AnPIT_area[[1]] <= 0.5)
  expect_length(out$model$model_x1$PIT[[1]], 2)

  summary_out <- summary(out)
  expect_true("AnPIT_area" %in% colnames(summary_out))
  expect_equal(summary_out["model_x1", "AnPIT_area"],
               out$model$model_x1$score$AnPIT_area[[1]])
})
