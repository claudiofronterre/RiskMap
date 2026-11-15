test_that("assess_prediction produces errors", {

  # A method is given, so that the check of `models` is the one reached.
  expect_input_error(
    assess_prediction("not list", method = "cluster"),
    "models"
  )

  expect_input_error(
    assess_prediction(list(a = 1), method = "cluster"),
    "models"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "random"),
    "method"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "cluster"),
    "fold", class = "riskmap_error_input"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "cluster",
                      fold = 0.1),
    "fold", class = "riskmap_error_input"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized"),
    "min_dist", class = "riskmap_error_input"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1),
    "size", class = "riskmap_error_input"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      iter = 2.1),
    "iter", class = "riskmap_error_input"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      control_mcmc = "not sim"),
    "control_mcmc"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = -1,
                      size = 1
    ),
    "min_dist", class = "riskmap_error_input"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1.1
    ),
    "size", class = "riskmap_error_input"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      keep_par_fixed = "not logical"
    ),
    "keep_par_fixed", class = "riskmap_error_input"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      control_mcmc = "not mcmc"
    ),
    "control_mcmc"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      plot_fold = "not true"
    ),
    "plot_fold", class = "riskmap_error_input"
  )

  expect_error(
    assess_prediction(list(gaussian_model),
                      method = "regularized",
                      min_dist = 1,
                      size = 1,
                      messages = "not true"
    ),
    "messages", class = "riskmap_error_input"
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
    "must be fitted to the same number of observations"
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
    "same rows in the same order"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(1:2, ncol = 1)),
    "user_split", "must have one row per observation"
  )

  # The number of test sets of a user split is its own (A4).
  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(c(1, rep(0, n - 1)), ncol = 1),
                      iter = 2),
    "iter", "number of test sets in `user_split` (1)"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(1, 2),
                      iter = 1),
    "iter", "number of test sets in `user_split` (2)"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = "not list or matrix"),
    "user_split", "must be a matrix or a list"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(sample(10))),
    "user_split", "must leave at least one observation for training"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(integer(0))),
    "user_split", "must not be empty"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(0)),
    "user_split", "must be row numbers between 1 and"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(c(1.1, 2, 3.3))),
    "user_split", "must be whole numbers"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(c(1, NA, 3))),
    "user_split", "must be numeric row numbers without missing values"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(c(1, 1, 2))),
    "user_split", "must not repeat a row"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(c(NA, rep(0, n - 1)), ncol = 1)),
    "user_split", "must not have missing values"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(c(2, rep(0, n - 1)), ncol = 1)),
    "user_split", "must hold only 0 for training and 1 for test"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(rep(0, n), ncol = 1)),
    "user_split", "one training observation (0) and one test observation (1)"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = matrix(rep(1, n), ncol = 1)),
    "user_split", "one training observation (0) and one test observation (1)"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c(1, 2, 3), out_id = c(3, 4)))),
    "user_split", "must not overlap"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c(1, 1, 2), out_id = 3:n))),
    "user_split", "must not repeat a row"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = 1:2, out_id = n + 1))),
    "user_split", "must be row numbers between 1 and"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = integer(0), out_id = 1:n))),
    "user_split", "must not be empty"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c(1, 2.5), out_id = 3:n))),
    "user_split", "must be whole numbers"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c("1", "2"), out_id = 3:n))),
    "user_split", "must be numeric row numbers without missing values"
  )

  expect_input_error(
    assess_prediction(list(gaussian_model),
                      method = "user",
                      user_split = list(list(in_id = c(1, NA), out_id = 3:n))),
    "user_split", "must be numeric row numbers without missing values"
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
    "same rows in the same order"
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
    "don't include the random effects"
  )

  expect_s3_class(out, "RiskMap_cross_validation")
})

test_that("assess_prediction splits test data correctly", {
  skip_if_not_installed("spatialsample")
  skip_if_not_installed("spatialEco")

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
  skip_if_not_installed("spatialEco")

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
  skip_if_not_installed("spatialEco")

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

test_that("assess_prediction pools PIT_area over all test observations (#230)", {
  geom <- st_sfc(
    st_point(c(0, 0)),
    st_point(c(1, 1)),
    st_point(c(2, 2)),
    st_point(c(3, 3)),
    crs = 4326
  )
  data <- st_sf(y = c(0, 0.5, -1, 2), x1 = c(10, 20, 30, 40), geometry = geom)

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
    user_split = list(1, 2:4),
    seed = 3,
    messages = FALSE
  )

  model <- out$model$model_x1
  # PIT_area is one number for all test observations, not one per test set.
  expect_named(model$metric, c("CRPS", "SCRPS"))
  expect_length(model$PIT, 2)
  u <- seq(0, 1, length.out = 1000)
  expect_identical(model$PIT_area,
                   .anpit_area(ecdf(unlist(model$PIT))(u), u))
  expect_gte(model$PIT_area_expected, 0)
  expect_lte(model$PIT_area_expected, 0.5)

  # Every predictive distribution has mean 0 and standard deviation 1, so the
  # SCRPS is the negative of Bolin and Wallin's, a loss (#230).
  crps <- function(y) 2 * dnorm(y) + y * (2 * pnorm(y) - 1) - 1 / sqrt(pi)
  y <- data$y
  expect_equal(unlist(model$metric$CRPS), crps(y))
  expected_scrps <- 0.5 * (1 + crps(y) * sqrt(pi) + log(2 / sqrt(pi)))
  expect_equal(unlist(model$metric$SCRPS), expected_scrps)

  summary_out <- summary(out)
  expect_identical(colnames(summary_out),
                   c("CRPS", "SCRPS", "PIT_area", "PIT_area_expected"))
  expect_identical(summary_out["model_x1", "PIT_area"], model$PIT_area)
  expect_identical(summary_out["model_x1", "PIT_area_expected"],
                   model$PIT_area_expected)

  # The expected value is simulated from a stream of its own, and the
  # caller's random numbers are unchanged.
  again <- assess_prediction(
    list(model_x1 = make_assess_prediction_fit(data, "x1")),
    method = "user", user_split = list(1, 2:4), seed = 3, messages = FALSE
  )
  expect_identical(again$model$model_x1$PIT_area_expected,
                   model$PIT_area_expected)
})

test_that("a sharper calibrated prediction has the lower SCRPS (#230)", {
  geom <- st_sfc(st_point(c(0, 0)), st_point(c(1, 1)), crs = 4326)
  data <- st_sf(y = c(0, 0), x1 = c(10, 20), geometry = geom)
  scrps_with_spread <- function(spread) {
    local_mocked_bindings(
      predict.RiskMap = function(object, newdata, ...) {
        n <- nrow(newdata)
        structure(list(
          predictions = transform(newdata, fixed = 0, offset = 0),
          samples = list(spatial = matrix(rep(spread * c(-1, 0, 1), each = n),
                                          nrow = n))
        ), class = "RiskMap_prediction")
      },
      .package = "RiskMap"
    )
    out <- assess_prediction(
      list(model = make_assess_prediction_fit(data, "x1")),
      method = "user", user_split = list(2), messages = FALSE,
      metrics = "SCRPS"
    )
    out$model$model$metric$SCRPS[[1]]
  }
  expect_lt(scrps_with_spread(0.5), scrps_with_spread(2))
})

test_that("the expected PIT_area shrinks with the test observations and allows for repeats (#230)", {
  gaussian_sets <- function(rows) {
    list(list(rows = rows, values = list(mean = rep(0, length(rows)),
                                         sd = rep(1, length(rows)))))
  }
  u <- seq(0, 1, length.out = 1000)
  set.seed(1)
  many <- expected_pit_area(gaussian_sets(1:100), FALSE, u)
  set.seed(1)
  few <- expected_pit_area(gaussian_sets(1:25), FALSE, u)
  # About 0.31 / sqrt(n) for n test observations
  expect_equal(many, 0.31 / sqrt(100), tolerance = 0.15)
  expect_equal(few, 0.31 / sqrt(25), tolerance = 0.15)
  # An observation in two test sets is drawn once, so 25 observations each
  # scored four times are no better than 25.
  set.seed(1)
  repeated <- expected_pit_area(gaussian_sets(rep(1:25, 4)), FALSE, u)
  expect_equal(repeated, few, tolerance = 0.05)

  # A discrete family uses the distribution functions of its predictions. Its
  # value agrees with drawing each count from its distribution and averaging
  # the areas of the AnPIT curves, with some counts of zero.
  set.seed(1)
  lambda <- rexp(100, 1 / 5)
  cdf <- lapply(lambda, function(l) cummax(ppois(0:60, l)))
  npit <- function(y, Fk) {
    below <- if (y == 0) 0 else Fk[y]
    at <- Fk[y + 1]
    if (at > below) pmin(pmax((u - below) / (at - below), 0), 1) else
      as.numeric(u > at)
  }
  brute <- replicate(400, {
    y <- vapply(cdf, function(Fk) findInterval(runif(1), Fk, left.open = TRUE),
                numeric(1))
    .anpit_area(rowMeans(vapply(seq_along(cdf), function(e) npit(y[e], cdf[[e]]),
                                numeric(length(u)))), u)
  })
  set.seed(2)
  discrete <- expected_pit_area(list(list(rows = seq_along(cdf), values = cdf)),
                                TRUE, u, n_replicates = 400L)
  expect_equal(discrete, mean(brute), tolerance = 0.08)
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
  skip_if_not_installed("spatialsample")
  skip_if_not_installed("spatialEco")
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
    "failed, and their scores are missing"
  )

  expect_equal(captured$invlink$inv(1), pnorm(1))
  expect_identical(out$failures$message, "refit failed")
  expect_identical(out$failures$test_set, 1L)
  expect_true(all(is.na(out$model$probit$metric$CRPS[[1]])))
  expect_true(is.na(out$model$probit$PIT_area))
  expect_true(is.na(out$model$probit$PIT_area_expected))
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
                "Parameters are fixed at the estimates from all the data.",
                fixed = TRUE)
})

test_that("assess_prediction warns when regularized test sets are smaller than 'size'", {
  skip_if_not_installed("spatialEco")
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
    "fewer than `size`"
  )
})

test_that("assess_prediction takes min_dist in the distance units of the fits (A5)", {
  skip_if_not_installed("spatialEco")
  in_metres <- glgpm(y ~ gp(), data = gaussian_data, family = "gaussian",
                     distance_units = "m", messages = FALSE)
  expect_identical(gaussian_intercept_model$distance_units, "km")
  kilometres <- assess_prediction(
    list(gaussian_intercept_model), method = "regularized", min_dist = 1,
    size = 3, control_mcmc = control_mcmc, seed = fixture_seed,
    messages = FALSE, metrics = "CRPS"
  )
  metres <- assess_prediction(
    list(in_metres), method = "regularized", min_dist = 1000,
    size = 3, control_mcmc = control_mcmc, seed = fixture_seed,
    messages = FALSE, metrics = "CRPS"
  )
  expect_identical(lapply(metres$test_set, st_geometry),
                   lapply(kilometres$test_set, st_geometry))
  expect_input_error(
    assess_prediction(list(gaussian_intercept_model, in_metres),
                      method = "regularized", min_dist = 1, size = 3,
                      messages = FALSE),
    "models"
  )
})

test_that("min_dist reaches subsample.distance() in metres whatever the CRS (A5)", {
  skip_if_not_installed("spatialEco")
  captured <- new.env()
  local_mocked_bindings(
    subsample.distance = function(x, size, d, ...) {
      captured$d <- d
      stop("captured")
    },
    .package = "spatialEco"
  )
  # spatialEco compares distances in metres, also on a CRS in kilometres.
  data <- st_transform(gaussian_data, st_crs("+proj=utm +zone=37 +units=km"))
  fit <- make_assess_prediction_fit(data, "1")
  fit$distance_units <- "km"
  expect_error(assess_prediction(list(fit), method = "regularized",
                                 min_dist = 2, size = 3, messages = FALSE),
               "captured")
  expect_identical(captured$d, 2000)
  fit$distance_units <- "m"
  expect_error(assess_prediction(list(fit), method = "regularized",
                                 min_dist = 2, size = 3, messages = FALSE),
               "captured")
  expect_identical(captured$d, 2)
})

test_that("assess_prediction checks the shape of a user split (A4)", {
  expect_input_error(
    assess_prediction(list(gaussian_intercept_model), method = "user",
                      user_split = list(in_id = 1:8, out_id = 9:10)),
    "user_split", "in a list of its own"
  )
  expect_input_error(
    assess_prediction(list(gaussian_intercept_model), method = "user",
                      user_split = matrix(numeric(0), nrow = 10, ncol = 0)),
    "user_split", "at least one test set"
  )
  expect_input_error(
    assess_prediction(list(gaussian_intercept_model), method = "user",
                      user_split = data.frame(a = 1:10)),
    "user_split", "not a data frame"
  )
})

test_that("assess_prediction draws nothing unless asked and counts the test sets of a list split (A4)", {
  expect_false(formals(assess_prediction)$plot_fold)
  out <- assess_prediction(
    list(gaussian_intercept_model), method = "user",
    user_split = list(1:2, 3:4, 5:6), control_mcmc = control_mcmc,
    seed = fixture_seed, messages = FALSE, metrics = "CRPS"
  )
  expect_length(out$test_set, 3)
  expect_length(out$model[[1]]$metric$CRPS, 3)
})
