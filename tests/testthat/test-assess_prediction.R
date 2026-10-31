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
    "'method' must be either 'cluster', 'regularized' or 'user'"
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
    "when 'method' is 'regularized' you must supply 'size'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      iter = 2.1),
    "'iter' must be a single positive integer"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      control_mcmc = "not sim"),
    "'control_mcmc' must come from 'set_control_mcmc"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = -1,
                      size = 1
    ),
    "'min_dist' must be a single positive"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1.1
    ),
    "'size' must be a single positive integer"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      keep_par_fixed = "not logical"
    ),
    "'keep_par_fixed' must be either TRUE or FALSE"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      control_mcmc = "not mcmc"
    ),
    "'control_mcmc' must come from 'set_control_mcmc"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      plot_fold = "not true"
    ),
    "'plot_fold' must be either TRUE or FALSE"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      messages = "not true"
    ),
    "'messages' must be either TRUE or FALSE"
  )

  different_rows <- gaussian_model
  different_rows$data <- different_rows$data[1:9,]

  expect_error(
    assess_prediction(list(gaussian_model,
                           different_rows),
                      method = "regularized",
                      min_dist = 1,
                      size = 1
    ),
    "All models in 'object' supplied must have the same number of observations"
  )

  different_order <- gaussian_model
  different_order$data <- different_order$data[c(6:10, 1:5),]

  expect_error(
    assess_prediction(list(gaussian_model,
                           different_order),
                      method = "regularized",
                      min_dist = 1,
                      size = 1
    ),
    "All models in 'object' must have data in the same row order and geometry"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(1:2, ncol = 1)),
    "'user_split' matrix must have the same number of rows as the data in the model"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(1:20, ncol = 2)),
    "'user_split' matrix must have a number of columns equal to 'iter'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(1, 2)),
    "'user_split' list must have the same length as 'iter'"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = "not list or matrix"),
    "'user_split' must be a matrix or a list"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(sample(10))),
    "must leave at least one observation for training"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(integer(0))),
    "test indices must be non-empty"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(0)),
    "test indices must be row indices between 1 and the number of observations"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(c(1.1, 2, 3.3))),
    "test indices must contain whole numbers"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(c(1, NA, 3))),
    "test indices must be a non-missing numeric vector of row indices"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(c(1, 1, 2))),
    "test indices must not contain duplicate indices"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(c(NA, rep(0, n - 1)), ncol = 1)),
    "'user_split' matrix must not contain missing values"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(c(2, rep(0, n - 1)), ncol = 1)),
    "'user_split' matrix must only contain 0s \\(training\\) and 1s \\(test\\)"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(rep(0, n), ncol = 1)),
    "must contain at least one training \\(0\\) and one test \\(1\\) observation"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(rep(1, n), ncol = 1)),
    "must contain at least one training \\(0\\) and one test \\(1\\) observation"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c(1, 2, 3), out_id = c(3, 4)))),
    "'in_id' and 'out_id' must not overlap"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c(1, 1, 2), out_id = 3:n))),
    "'in_id' must not contain duplicate indices"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = 1:2, out_id = n + 1))),
    "'out_id' must be row indices between 1 and the number of observations"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = integer(0), out_id = 1:n))),
    "'in_id' must be non-empty"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c(1, 2.5), out_id = 3:n))),
    "'in_id' must contain whole numbers"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c("1", "2"), out_id = 3:n))),
    "'in_id' must be a non-missing numeric vector of row indices"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c(1, NA), out_id = 3:n))),
    "'in_id' must be a non-missing numeric vector of row indices"
  )

})

make_assess_prediction_fit <- function(data, covariate) {
  coords <- sf::st_coordinates(data)
  fit <- list(
    formula = stats::as.formula(paste("y ~", covariate)),
    data = data,
    input_crs = sf::st_crs(data),
    family = "gaussian",
    estimate = list(beta = c(0, 0), sigma2 = 0, phi = 0),
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
    distance_units = "m",
    call = list(denominator = quote(units_m)),
    model_id = covariate
  )
  class(fit) <- "RiskMap"
  fit
}

test_that("assess_prediction uses each model's own data for held-out predictors", {
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
    predict.RiskMap = function(object, newdata, ...) {
      seen_predictors[[object$model_id]] <<- names(newdata)
      n <- nrow(newdata)
      structure(list(
        predictions = transform(newdata, fixed = 0, offset = 0),
        samples = list(spatial = matrix(rep(1, n * 2), nrow = n))
      ), class = "RiskMap_prediction")
    },
    .package = "RiskMap"
  )

  out <- assess_prediction(
    list(model_x1 = make_assess_prediction_fit(data_x1, "x1"),
         model_x2 = make_assess_prediction_fit(data_x2, "x2")),
    method = "user",
    user_split = matrix(c(0, 1, 1), ncol = 1),
    plot_fold = FALSE,
    messages = FALSE,
    metrics = "CRPS"
  )

  expect_s3_class(out, "RiskMap_cross_validation")

  expect_setequal(names(out), c("test_set", "model", "keep_par_fixed", "failures"))

  expect_true("x1" %in% seen_predictors$x1)
  expect_false("x2" %in% seen_predictors$x1)
  expect_true("x2" %in% seen_predictors$x2)
  expect_false("x1" %in% seen_predictors$x2)
})

test_that("assess_prediction preserves a custom model CRS when refitting", {

  custom_crs <- st_crs(
    "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
  )
  custom_data <- st_transform(gaussian_data, custom_crs)
  custom_model <- glgpm(y ~ cov + gp(),
                        data = custom_data,
                        family = "gaussian",
                        messages = FALSE)
  user_split <- matrix(c(1, 1, rep(0, nrow(custom_data) - 2)), ncol = 1)

  expect_true(is.na(st_crs(custom_model$data)$epsg))
  expect_no_error(
    assess_prediction(
      list(model = custom_model),
      method = "user",
      keep_par_fixed = FALSE,
      user_split = user_split,
      control_mcmc = control_mcmc, seed = fixture_seed,
      plot_fold = FALSE,
      messages = FALSE,
      metrics = "CRPS"
    )
  )
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
      method = "user",
      user_split = matrix(c(0, 1), ncol = 1),
      plot_fold = FALSE,
      messages = FALSE,
      metrics = "CRPS"
    ),
    "same row order and geometry"
  )
})

test_that("assess_prediction re-encodes random effects after subsetting", {
  user_split <- matrix(c(1, 1, rep(0, nrow(gaussian_data) - 2)), ncol = 1)

  expect_warning(
    out <- assess_prediction(
      list(model = gaussian_model),
      method = "user",
      user_split = user_split,
      control_mcmc = control_mcmc, seed = fixture_seed,
      plot_fold = FALSE,
      messages = FALSE,
      metrics = "CRPS"
    ),
    "do not include the random effects"
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
  expect_identical(sum(unlist(lapply(result$test_set, nrow))), n)

  combined <- do.call(rbind, result$test_set)
  expect_true(!all(duplicated(combined)))

  n_folds <- 3

  result <- assess_prediction(
    list(intercept_only = gaussian_intercept_model,
         with_covariate = gaussian_model),
    method = "cluster",
    fold = n_folds,
    messages = FALSE)

  expect_length(result$test_set, n_folds)
  expect_identical(sum(unlist(lapply(result$test_set, nrow))), n)

  combined <- do.call(rbind, result$test_set)
  expect_true(!all(duplicated(combined)))

  size <- 4L

  result <- assess_prediction(
    list(intercept_only = gaussian_intercept_model,
         with_covariate = gaussian_model),
    method = "regularized",
    size = size,
    min_dist = 1,
    messages = FALSE)

  expect_length(result$test_set, 1)
  expect_identical(nrow(result$test_set[[1]]), size)


  result <- assess_prediction(
    list(gaussian_model),
    method = "user",
    user_split = matrix(
      sample(c(rep(1, n/2), rep(0, n/2))),
      ncol = 1),
    messages = FALSE)

  expect_length(result$test_set, 1)
  expect_identical(nrow(result$test_set[[1]]), as.integer(n/2))

  result <- assess_prediction(
    list(gaussian_model),
    method = "user",
    user_split = list(list(in_id = 1:(n - 2), out_id = (n - 1):n)),
    messages = FALSE)

  expect_length(result$test_set, 1)
  expect_identical(nrow(result$test_set[[1]]), 2L)

  result <- assess_prediction(
    list(gaussian_model),
    method = "user",
    user_split = list((n - 1):n),
    messages = FALSE)

  expect_length(result$test_set, 1)
  expect_identical(nrow(result$test_set[[1]]), 2L)

})

test_that("assess_prediction can refit correctly for all model families", {

  # 'regularized' holds out a single random point from the tiny (n = 10)
  # fixtures; refitting on whichever 9 points remain is occasionally
  # numerically singular for an unlucky holdout, so a seed is used here to
  # keep the test deterministic (was intermittently failing on CI). The
  # fixture seed (see helper-data.R) is verified stable for all three
  # families.
  result <- assess_prediction(
    list(gaussian_model),
    method = "regularized",
    min_dist = 1,
    size = 1,
    keep_par_fixed = FALSE,
    control_mcmc = control_mcmc, seed = fixture_seed,
    messages = FALSE)

  expect_setequal(names(result), c("test_set", "model", "keep_par_fixed", "failures"))

  result <- assess_prediction(
    list(binomial_model),
    method = "regularized",
    min_dist = 1,
    size = 1,
    keep_par_fixed = FALSE,
    control_mcmc = control_mcmc, seed = fixture_seed,
    messages = FALSE)

  expect_setequal(names(result), c("test_set", "model", "keep_par_fixed", "failures"))

  expect_warning(
    result <- assess_prediction(
      list(poisson_model),
      method = "regularized",
      min_dist = 1,
      size = 1,
      keep_par_fixed = FALSE,
      control_mcmc = control_mcmc, seed = fixture_seed,
      messages = FALSE),
    "estimated as effectively zero")

  expect_setequal(names(result), c("test_set", "model", "keep_par_fixed", "failures"))

})

test_that("assess_prediction's seed makes regularized/cluster splits reproducible and restores the RNG state", {

  rng_before <- .Random.seed

  result1 <- assess_prediction(
    list(gaussian_model),
    method = "regularized",
    min_dist = 1,
    size = 1,
    messages = FALSE,
    plot_fold = FALSE,
    control_mcmc = control_mcmc, seed = fixture_seed)

  expect_identical(rng_before, .Random.seed)

  result2 <- assess_prediction(
    list(gaussian_model),
    method = "regularized",
    min_dist = 1,
    size = 1,
    messages = FALSE,
    plot_fold = FALSE,
    control_mcmc = control_mcmc, seed = fixture_seed)

  expect_identical(result1$test_set, result2$test_set)
})


test_that("PIT area computes trapezoidal absolute distance", {
  u <- seq(0, 1, length.out = 1001)

  expect_identical(.anpit_area(u, u), 0)
  expect_equal(.anpit_area(rep(0, length(u)), u), 0.5)
  expect_equal(.anpit_area(u^2, u), 1 / 6, tolerance = 1e-6)
})

test_that("assess_prediction reports PIT area as a scalar score", {
  geom <- st_sfc(
    st_point(c(0, 0)),
    st_point(c(1, 1)),
    st_point(c(2, 2)),
    crs = 4326
  )
  data <- st_sf(y = c(0, 0, 0), x1 = c(10, 20, 30), geometry = geom)

  testthat::local_mocked_bindings(
    predict.RiskMap = function(object, newdata, ...) {
      n <- nrow(newdata)
      structure(list(
        predictions = transform(newdata, fixed = 0, offset = 0),
        samples = list(
          spatial = matrix(rep(c(-1, 0, 1), each = n), nrow = n)
        )
      ), class = "RiskMap_prediction")
    },
    .package = "RiskMap"
  )

  out <- assess_prediction(
    list(model_x1 = make_assess_prediction_fit(data, "x1")),
    method = "user",
    user_split = matrix(c(0, 1, 1), ncol = 1),
    plot_fold = FALSE,
    messages = FALSE,
    metrics = "PIT"
  )

  expect_named(out$model$model_x1$metric, "PIT_area")
  expect_length(out$model$model_x1$metric$PIT_area, 1)
  expect_type(out$model$model_x1$metric$PIT_area[[1]], "double")
  expect_true(is.finite(out$model$model_x1$metric$PIT_area[[1]]))
  expect_gte(out$model$model_x1$metric$PIT_area[[1]], 0)
  expect_lte(out$model$model_x1$metric$PIT_area[[1]], 0.5)
  expect_length(out$model$model_x1$PIT[[1]], 2)

  summary_out <- summary(out)
  expect_true("PIT_area" %in% colnames(summary_out))
  expect_identical(summary_out["model_x1", "PIT_area"],
               out$model$model_x1$metric$PIT_area[[1]])

})

test_that("plot_folds() plots every observation, coloured by training/test set (#157)", {
  n <- nrow(gaussian_data)
  data_split <- list(splits = list(
    list(data = gaussian_data[1:7, ], data_test = gaussian_data[8:10, ]),
    list(data = gaussian_data[c(1:5, 9:10), ], data_test = gaussian_data[6:8, ])
  ))

  p <- plot_folds(data_split)
  expect_s3_class(p, "ggplot")

  built <- ggplot2::ggplot_build(p)
  expect_identical(nrow(built$data[[1]]), 2L * n)
  expect_length(unique(built$data[[1]]$PANEL), 2)
  expect_length(unique(built$data[[1]]$colour), 2)
})

test_that("plot_folds() doesn't facet when there is only one iteration", {
  data_split <- list(splits = list(
    list(data = gaussian_data[1:7, ], data_test = gaussian_data[8:10, ])
  ))

  p <- plot_folds(data_split)
  built <- ggplot2::ggplot_build(p)
  expect_length(unique(built$data[[1]]$PANEL), 1)
  expect_identical(nrow(built$data[[1]]), nrow(gaussian_data))
})

test_that("assess_prediction(plot_fold = TRUE) uses plot_folds() for 'user' and 'regularized', autoplot() for 'cluster' (#157)", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  user_split <- matrix(c(1, 1, rep(0, nrow(gaussian_data) - 2)), ncol = 1)
  expect_no_error(
    assess_prediction(list(model = gaussian_model), method = "user",
                      user_split = user_split, control_mcmc = control_mcmc,
                      seed = fixture_seed,
                      plot_fold = TRUE, messages = FALSE, metrics = "CRPS")
  )

  expect_no_error(
    assess_prediction(list(model = gaussian_model), method = "regularized",
                      min_dist = 1, size = 2, control_mcmc = control_mcmc,
                      seed = fixture_seed,
                      plot_fold = TRUE, messages = FALSE, metrics = "CRPS")
  )

  expect_no_error(
    assess_prediction(list(model = gaussian_model), method = "cluster",
                      fold = 2, control_mcmc = control_mcmc,
                      seed = fixture_seed,
                      plot_fold = TRUE, messages = FALSE, metrics = "CRPS")
  )
})

# A minimal binomial fit with a probit link, scored through a mocked
# prediction whose linear predictor is 1 in every draw
make_probit_fit <- function() {
  geom <- st_sfc(st_point(c(0, 0)), st_point(c(1000, 0)),
                 st_point(c(0, 1000)), crs = 32637)
  data <- st_sf(y = c(3, 9, 6), n = c(10, 10, 10), geometry = geom)
  fit <- make_assess_prediction_fit(data, "1")
  fit$family <- "binomial"
  fit$units_m <- data$n
  fit$estimate <- list(beta = 0, sigma2 = 0, phi = 0)
  fit$link_function <- list(inv = function(x) pnorm(x), name = "custom")
  fit
}

mock_constant_prediction <- function(value) {
  function(object, newdata, ...) {
    n <- nrow(newdata)
    structure(list(
      predictions = transform(newdata, fixed = 0, offset = 0),
      samples = list(spatial = matrix(value, nrow = n, ncol = 4000))
    ), class = "RiskMap_prediction")
  }
}

test_that("assess_prediction scores with the model's inverse link", {
  local_mocked_bindings(predict.RiskMap = mock_constant_prediction(1),
                        .package = "RiskMap")
  out <- assess_prediction(
    list(probit = make_probit_fit()),
    method = "user",
    user_split = matrix(c(0, 1, 1), ncol = 1),
    seed = 1,
    plot_fold = FALSE,
    messages = FALSE,
    metrics = "CRPS"
  )

  exact_crps <- function(y, p) {
    F_k <- pbinom(0:10, 10, p)
    sum((F_k - as.numeric(0:10 >= y))^2)
  }
  expected <- vapply(c(9, 6), exact_crps, numeric(1), p = pnorm(1))
  expect_equal(out$model$probit$metric$CRPS[[1]], expected, tolerance = 0.05)
  logit_crps <- vapply(c(9, 6), exact_crps, numeric(1), p = plogis(1))
  expect_gt(abs(out$model$probit$metric$CRPS[[1]][1] - logit_crps[1]), 0.2)
})

test_that("assess_prediction refits with the model's link and records failed refits", {
  captured <- new.env()
  local_mocked_bindings(
    glgpm = function(..., invlink = NULL) {
      captured$invlink <- invlink
      stop("refit failed")
    },
    .package = "RiskMap"
  )
  expect_warning(
    out <- assess_prediction(
      list(probit = make_probit_fit()),
      method = "user",
      user_split = matrix(c(0, 1, 1), ncol = 1),
      keep_par_fixed = FALSE,
      plot_fold = FALSE,
      messages = FALSE,
      metrics = c("CRPS", "PIT")
    ),
    "failed and their scores are missing"
  )

  expect_equal(captured$invlink$inv(1), pnorm(1))
  expect_identical(out$failures$message, "refit failed")
  expect_identical(out$failures$test_set, 1L)
  expect_true(all(is.na(out$model$probit$metric$CRPS[[1]])))
  expect_true(is.na(out$model$probit$metric$PIT_area[[1]]))
  expect_true(is.na(summary(out)["probit", "CRPS"]))
})

test_that("assess_prediction keeps test rows with missing values in unused columns", {
  geom <- st_sfc(st_point(c(0, 0)), st_point(c(1, 1)), st_point(c(2, 2)),
                 crs = 4326)
  data <- st_sf(y = c(1, 2, 3), x1 = c(10, 20, 30), unused = c(1, NA, NA),
                geometry = geom)
  local_mocked_bindings(predict.RiskMap = mock_constant_prediction(0),
                        .package = "RiskMap")

  out <- assess_prediction(
    list(model_x1 = make_assess_prediction_fit(data, "x1")),
    method = "user",
    user_split = matrix(c(0, 1, 1), ncol = 1),
    plot_fold = FALSE,
    messages = FALSE,
    metrics = "CRPS"
  )

  expect_identical(nrow(out$test_set[[1]]), 2L)
  expect_length(out$model$model_x1$metric$CRPS[[1]], 2)
})

test_that("assess_prediction gives the same model the same scores in any position", {
  # Each test set has one seed shared by all models, and the Gaussian
  # measurement error is drawn from a separate stream
  user_split <- matrix(c(1, 1, 1, rep(0, nrow(gaussian_data) - 3)), ncol = 1)
  out <- assess_prediction(
    list(a = gaussian_intercept_model, b = gaussian_intercept_model),
    method = "user",
    user_split = user_split,
    control_mcmc = control_mcmc, seed = fixture_seed,
    plot_fold = FALSE,
    messages = FALSE
  )
  expect_identical(out$model$a$metric, out$model$b$metric)
  expect_identical(out$model$a$PIT, out$model$b$PIT)
})

test_that("assess_prediction reports how the parameters are treated", {
  user_split <- matrix(c(1, 1, rep(0, nrow(gaussian_data) - 2)), ncol = 1)
  expect_message(
    out <- assess_prediction(
      list(m0 = gaussian_intercept_model, m1 = gaussian_offset_model),
      method = "user",
      user_split = user_split,
      control_mcmc = control_mcmc, seed = fixture_seed,
      plot_fold = FALSE,
      messages = TRUE,
      metrics = "CRPS"
    ),
    "different numbers of parameters"
  )
  expect_true(out$keep_par_fixed)
  expect_output(print(summary(out)),
                "Parameters: fixed at the estimates from all the data")
})

test_that("assess_prediction warns when regularized test sets are smaller than 'size'", {
  expect_warning(
    assess_prediction(
      list(gaussian_intercept_model),
      method = "regularized",
      min_dist = 50,
      size = 3,
      control_mcmc = control_mcmc, seed = fixture_seed,
      plot_fold = FALSE,
      messages = FALSE,
      metrics = "CRPS"
    ),
    "fewer than 'size'"
  )
})
