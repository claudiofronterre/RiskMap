make_assess_pp_fit <- function(data, covariate) {
  coords <- sf::st_coordinates(data)
  fit <- list(
    formula = stats::as.formula(paste("y ~", covariate)),
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

test_that("assess_pp uses each model's own data_sf for held-out predictors", {
  geom <- sf::st_sfc(
    sf::st_point(c(0, 0)),
    sf::st_point(c(1, 1)),
    sf::st_point(c(2, 2)),
    crs = 4326
  )
  data_x1 <- sf::st_sf(y = c(1, 2, 3), x1 = c(10, 20, 30), geometry = geom)
  data_x2 <- sf::st_sf(y = c(1, 2, 3), x2 = c(100, 200, 300), geometry = geom)

  seen_predictors <- list()
  testthat::local_mocked_bindings(
    pred_over_grid = function(object, grid_pred, predictors, ...) {
      seen_predictors[[object$model_id]] <<- names(predictors)
      list(predictors = predictors)
    },
    pred_target_grid = function(object, ...) {
      list(lp_samples = matrix(rep(1, nrow(object$predictors) * 2),
                               nrow = nrow(object$predictors)))
    },
    .package = "RiskMap"
  )

  out <- assess_pp(
    list(model_x1 = make_assess_pp_fit(data_x1, "x1"),
         model_x2 = make_assess_pp_fit(data_x2, "x2")),
    user_split = matrix(c(0, 1, 1), ncol = 1),
    plot_fold = FALSE,
    messages = FALSE,
    which_metric = "CRPS"
  )

  expect_s3_class(out, "RiskMap.spatial.cv")
  expect_true("x1" %in% seen_predictors$x1)
  expect_false("x2" %in% seen_predictors$x1)
  expect_true("x2" %in% seen_predictors$x2)
  expect_false("x1" %in% seen_predictors$x2)
})

test_that("assess_pp requires aligned model data", {
  geom_a <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 4326)
  geom_b <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(2, 2)), crs = 4326)
  data_a <- sf::st_sf(y = c(1, 2), x1 = c(10, 20), geometry = geom_a)
  data_b <- sf::st_sf(y = c(1, 2), x2 = c(100, 200), geometry = geom_b)

  expect_error(
    assess_pp(
      list(model_x1 = make_assess_pp_fit(data_a, "x1"),
           model_x2 = make_assess_pp_fit(data_b, "x2")),
      user_split = matrix(c(0, 1), ncol = 1),
      plot_fold = FALSE,
      messages = FALSE,
      which_metric = "CRPS"
    ),
    "same row order and geometry"
  )
})
