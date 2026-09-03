test_that("pred_target_shp produces expected output with default arguments", {

  expected_output <- c("lp_samples", "target", "shp", "f_target", "pd_summary", "grid_pred")

  gaussian_grid <- setup_prediction(gaussian_model, type = "joint")
  gaussian_offset_grid <- setup_prediction(gaussian_offset_model, type = "joint")
  binomial_grid <- setup_prediction(binomial_model, control_sim = control_mcmc, type = "joint")
  poisson_grid <- setup_prediction(poisson_model, control_sim = control_mcmc, type = "joint")

  result <- pred_target_shp(gaussian_grid, areal)
  expect_setequal(names(result), expected_output)

  result <- pred_target_shp(gaussian_offset_grid, areal)
  expect_setequal(names(result), expected_output)

  result <- pred_target_shp(binomial_grid, areal)
  expect_setequal(names(result), expected_output)

  result <- pred_target_shp(poisson_grid, areal)
  expect_setequal(names(result), expected_output)

})

test_that("pred_target_shp preserves posterior samples for one-pixel list-mode regions", {
  grid_pred <- list(
    group_one = sf::st_as_sf(data.frame(x = 0, y = 0), coords = c("x", "y"), crs = 4326),
    group_two = sf::st_as_sf(data.frame(x = c(1, 2), y = c(1, 2)), coords = c("x", "y"), crs = 4326)
  )

  polygon_1 <- sf::st_polygon(list(matrix(c(0,0, 0.1,0, 0.1,0.1, 0,0.1, 0,0), ncol = 2, byrow = TRUE)))
  polygon_2 <- sf::st_polygon(list(matrix(c(1,1, 1.5,1, 1.5,1.5, 1,1.5, 1,1), ncol = 2, byrow = TRUE)))
  shp <- sf::st_sf(region = c("group_one", "group_two"),
                   geometry = sf::st_sfc(polygon_1, polygon_2), crs = sf::st_crs(4326))

  object <- list(
    type = "joint",
    grid_pred = grid_pred,
    S_samples = list(
      matrix(c(0.01, 0.02, 0.03), nrow = 1),
      matrix(c(0.10, 0.10, 0.10, 0.40, 0.40, 0.40), nrow = 2, byrow = TRUE)
    ),
    mu_pred = list(c(0), c(0, 0)),
    cov_offset = list(c(0), c(0, 0)),
    re = list(samples = list()),
    par_hat = list()
  )
  class(object) <- "RiskMap_pred"

  out <- pred_target_shp(
    object,
    shp = shp,
    shp_target = sum,
    weights = list(1, c(0.25, 0.75)),
    standardize_weights = FALSE,
    col_names = "region",
    f_target = list(identity_target = identity),
    pd_summary = list(mean = mean),
    return_shp = FALSE,
    return_target_samples = TRUE,
    messages = FALSE
  )

  expect_equal(out$target_samples$group_one$identity_target, c(0.01, 0.02, 0.03))
  expect_equal(out$target$group_one$identity_target$mean, 0.02)
  expect_equal(out$target_samples$group_two$identity_target, rep(0.325, 3))
})

test_that("pred_target_shp errors on wrong list-mode target orientation", {
  grid_pred <- list(
    group_one = sf::st_as_sf(data.frame(x = c(0, 1), y = c(0, 1)), coords = c("x", "y"), crs = 4326)
  )
  shp <- sf::st_sf(
    region = "group_one",
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 0.1,0, 0.1,0.1, 0,0.1, 0,0), ncol = 2, byrow = TRUE))), crs = 4326)
  )
  object <- list(
    type = "joint",
    grid_pred = grid_pred,
    S_samples = list(matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06), nrow = 2)),
    mu_pred = list(c(0, 0)),
    cov_offset = list(c(0, 0)),
    re = list(samples = list()),
    par_hat = list()
  )
  class(object) <- "RiskMap_pred"

  expect_error(
    pred_target_shp(
      object,
      shp = shp,
      weights = list(c(0.5, 0.5)),
      col_names = "region",
      f_target = list(bad_target = function(x) t(x)),
      pd_summary = list(mean = mean),
      return_shp = FALSE,
      messages = FALSE
    ),
    "expected a 2 x 3 matrix"
  )
})
