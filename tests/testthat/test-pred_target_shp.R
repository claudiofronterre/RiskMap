test_that("predict_areal_target produces expected output with default arguments", {

  expected_output <- c("lp_samples", "target", "boundaries", "f_target", "pd_summary", "grid_pred",
                       "grid_summary_range")

  gaussian_grid <- setup_prediction(gaussian_model, type = "joint")
  gaussian_offset_grid <- setup_prediction(gaussian_offset_model, type = "joint")
  binomial_grid <- setup_prediction(binomial_model, control_mcmc = control_mcmc, type = "joint")
  poisson_grid <- setup_prediction(poisson_model, control_mcmc = control_mcmc, type = "joint")

  result <- predict_areal_target(gaussian_grid, areal)
  expect_setequal(names(result), expected_output)

  result <- predict_areal_target(gaussian_offset_grid, areal)
  expect_setequal(names(result), expected_output)

  result <- predict_areal_target(binomial_grid, areal)
  expect_setequal(names(result), expected_output)

  result <- predict_areal_target(poisson_grid, areal)
  expect_setequal(names(result), expected_output)

})

test_that("predict_areal_target preserves posterior samples for one-pixel list-mode regions", {
  grid_pred <- list(
    group_one = sf::st_as_sf(data.frame(x = 0, y = 0), coords = c("x", "y"), crs = 4326),
    group_two = sf::st_as_sf(data.frame(x = c(1, 2), y = c(1, 2)), coords = c("x", "y"), crs = 4326)
  )

  polygon_1 <- sf::st_polygon(list(matrix(c(0,0, 0.1,0, 0.1,0.1, 0,0.1, 0,0), ncol = 2, byrow = TRUE)))
  polygon_2 <- sf::st_polygon(list(matrix(c(1,1, 1.5,1, 1.5,1.5, 1,1.5, 1,1), ncol = 2, byrow = TRUE)))
  boundaries <- sf::st_sf(region = c("group_one", "group_two"),
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

  out <- predict_areal_target(
    object,
    boundaries = boundaries,
    areal_target = sum,
    weights = list(1, c(0.25, 0.75)),
    standardize_weights = FALSE,
    col_names = "region",
    f_target = list(identity_target = identity),
    pd_summary = list(mean = mean),
    return_boundaries = FALSE,
    return_target_samples = TRUE,
    messages = FALSE
  )

  expect_equal(out$target_samples$group_one$identity_target, c(0.01, 0.02, 0.03))
  expect_equal(out$target$group_one$identity_target$mean, 0.02)
  expect_equal(out$target_samples$group_two$identity_target, rep(0.325, 3))
  expect_false("boundaries" %in% names(out))
})

test_that("predict_areal_target errors on wrong list-mode target orientation", {
  grid_pred <- list(
    group_one = sf::st_as_sf(data.frame(x = c(0, 1), y = c(0, 1)), coords = c("x", "y"), crs = 4326)
  )
  boundaries <- sf::st_sf(
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
    predict_areal_target(
      object,
      boundaries = boundaries,
      weights = list(c(0.5, 0.5)),
      col_names = "region",
      f_target = list(bad_target = function(x) t(x)),
      pd_summary = list(mean = mean),
      return_boundaries = FALSE,
      messages = FALSE
    ),
    "expected a 2 x 3 matrix"
  )
})

test_that("plot.RiskMap_predict_areal_target defaults to the first target and validates target/summary #143", {

  gaussian_grid <- setup_prediction(gaussian_model, type = "joint")
  result <- predict_areal_target(gaussian_grid, areal, messages = FALSE)

  expect_no_error(ggplot2::ggplot_build(plot(result)))
  expect_no_error(ggplot2::ggplot_build(plot(result, target = result$f_target[1])))

  expect_error(
    plot(result, target = "not_a_target"),
    "'target' must be one of"
  )
  expect_error(
    plot(result, summary = "not_a_summary"),
    "'summary' must be one of"
  )
})

test_that("plot.RiskMap_predict_areal_target shares palette and map annotations with the grid plot #149", {

  gaussian_grid <- setup_prediction(gaussian_model, type = "joint")
  result <- predict_areal_target(gaussian_grid, areal, messages = FALSE)

  out <- plot(result)
  expect_s3_class(out, "ggplot")
  expect_no_error(ggplot2::ggplot_build(out))
  layer_classes <- vapply(out$layers, function(layer) class(layer$geom)[1], character(1))
  expect_true("GeomSf" %in% layer_classes)
  expect_true("GeomNorthArrow" %in% layer_classes)
  expect_true("GeomScaleBar" %in% layer_classes)
  expect_equal(out$labels$fill, "linear_target_mean")

  out_plain <- plot(result, north_arrow = FALSE, scale_bar = FALSE)
  plain_classes <- vapply(out_plain$layers, function(layer) class(layer$geom)[1], character(1))
  expect_equal(unname(plain_classes), "GeomSf")

  expect_no_error(ggplot2::ggplot_build(plot(result, palette = "Spectral",
                                             reverse_palette = TRUE)))
  expect_error(plot(result, palette = "not_a_palette"), "'palette' must be one of")
})

test_that("plot.RiskMap_predict_areal_target errors informatively without boundaries #149", {

  gaussian_grid <- setup_prediction(gaussian_model, type = "joint")
  result <- predict_areal_target(gaussian_grid, areal, return_boundaries = FALSE,
                                 messages = FALSE)

  expect_error(plot(result), "return_boundaries = TRUE")
})

test_that("plot.RiskMap_predict_areal_target shares its colour range with the grid plot by default #149", {

  gaussian_grid <- setup_prediction(gaussian_model, type = "joint")
  grid_result <- predict_grid_target(gaussian_grid)
  areal_result <- predict_areal_target(gaussian_grid, areal, messages = FALSE)

  fill_limits <- function(map) {
    ggplot2::ggplot_build(map)$plot$scales$get_scales("fill")$get_limits()
  }

  expect_equal(areal_result$grid_summary_range$linear_target$mean,
               range(grid_result$target$linear_target$mean))
  expect_equal(areal_result$grid_summary_range$linear_target$sd,
               range(grid_result$target$linear_target$sd))
  expect_equal(fill_limits(plot(areal_result)), fill_limits(plot(grid_result)))

  # areal sd is smaller than cell-level sd, so the range must still cover it
  sd_limits <- fill_limits(plot(areal_result, summary = "sd"))
  expect_true(all(areal_result$boundaries$linear_target_sd >= sd_limits[1] &
                    areal_result$boundaries$linear_target_sd <= sd_limits[2]))

  expect_equal(fill_limits(plot(areal_result, limits = "shared")), fill_limits(plot(areal_result)))
})

test_that("prediction maps accept independent and custom colour limits #149", {

  gaussian_grid <- setup_prediction(gaussian_model, type = "joint")
  grid_result <- predict_grid_target(gaussian_grid)
  areal_result <- predict_areal_target(gaussian_grid, areal, messages = FALSE)

  fill_limits <- function(map) {
    ggplot2::ggplot_build(map)$plot$scales$get_scales("fill")$get_limits()
  }

  expect_equal(fill_limits(plot(areal_result, limits = "independent")),
               range(areal_result$boundaries$linear_target_mean))
  expect_equal(fill_limits(plot(grid_result, limits = "independent")),
               range(grid_result$target$linear_target$mean))
  expect_equal(fill_limits(plot(grid_result, limits = "shared")),
               range(grid_result$target$linear_target$mean))

  expect_equal(fill_limits(plot(areal_result, limits = c(-10, 10))), c(-10, 10))
  expect_equal(fill_limits(plot(grid_result, limits = c(-10, 10))), c(-10, 10))

  for (result in list(areal_result, grid_result)) {
    expect_error(plot(result, limits = "other"),
                 "'limits' must be \"shared\", \"independent\" or a numeric vector of length two")
    expect_error(plot(result, limits = NULL),
                 "'limits' must be \"shared\", \"independent\" or a numeric vector of length two")
    expect_error(plot(result, limits = 1), "Custom 'limits' must be a numeric vector of length two")
    expect_error(plot(result, limits = c(0, NA)), "Custom 'limits' must be a numeric vector of length two")
  }
})

test_that("prediction maps add space for the north arrow and scale bar only when shown #149", {

  gaussian_grid <- setup_prediction(gaussian_model, type = "joint")
  areal_result <- predict_areal_target(gaussian_grid, areal, messages = FALSE)
  grid_result <- predict_grid_target(gaussian_grid)

  y_range <- function(map) ggplot2::ggplot_build(map)$layout$panel_params[[1]]$y_range

  for (result in list(areal_result, grid_result)) {
    plain <- y_range(plot(result, north_arrow = FALSE, scale_bar = FALSE))
    arrow_only <- y_range(plot(result, north_arrow = TRUE, scale_bar = FALSE))
    scale_only <- y_range(plot(result, north_arrow = FALSE, scale_bar = TRUE))

    expect_gt(arrow_only[2], plain[2])
    expect_lt(scale_only[1], plain[1])
  }

  extent <- c(xmin = 0, xmax = 10, ymin = 0, ymax = 100)
  expect_equal(pad_map_extent(extent, north_arrow = FALSE, scale_bar = FALSE), c(0, 100))
  expect_equal(pad_map_extent(extent, north_arrow = TRUE, scale_bar = FALSE), c(0, 115))
  expect_equal(pad_map_extent(extent, north_arrow = FALSE, scale_bar = TRUE), c(-8, 100))
})
