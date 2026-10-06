test_that("predict_grid_target produces expected output with default arguments", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- setup_prediction(gaussian_model)
  gaussian_offset_grid <- setup_prediction(gaussian_offset_model)
  binomial_grid <- setup_prediction(binomial_model, control_mcmc = control_mcmc)
  poisson_grid <- setup_prediction(poisson_model, control_mcmc = control_mcmc)

  result <- predict_grid_target(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(gaussian_offset_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(poisson_grid)
  expect_setequal(names(result), expected_output)

})

test_that("predict_grid_target produces expected output when grid is provided", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- setup_prediction(gaussian_model,
                                  grid_pred = grid,
                                  predictors = data.frame(cov = rnorm(nrow(grid))),
                                  re_predictors = data.frame(i = 1:5),
                                  type = "joint")

  binomial_grid <- setup_prediction(binomial_model,
                                  grid_pred = grid,
                                  predictors = data.frame(cov = rnorm(nrow(grid))),
                                  re_predictors = data.frame(i = 1:5),
                                  control_mcmc = control_mcmc,
                                  type = "joint")

  poisson_grid <- setup_prediction(poisson_model,
                                 grid_pred = grid,
                                 predictors = data.frame(cov = rnorm(nrow(grid))),
                                 re_predictors = data.frame(i = 1:5),
                                 control_mcmc = control_mcmc,
                                 type = "joint")

  result <- predict_grid_target(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(poisson_grid)
  expect_setequal(names(result), expected_output)

})


test_that("predict_grid_target produces expected output when in list mode", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- setup_prediction(gaussian_model,
                                  grid_pred = list(grid, grid),
                                  predictors = list(data.frame(cov = rnorm(nrow(grid))),
                                                    data.frame(cov = rnorm(nrow(grid)))),
                                  type = "joint")

  binomial_grid <- setup_prediction(binomial_model,
                                  grid_pred = list(grid, grid),
                                  predictors = list(data.frame(cov = rnorm(nrow(grid))),
                                                    data.frame(cov = rnorm(nrow(grid)))),
                                  control_mcmc = control_mcmc,
                                  type = "joint")

  poisson_grid <- setup_prediction(poisson_model,
                                 grid_pred = list(grid, grid),
                                 predictors = list(data.frame(cov = rnorm(nrow(grid))),
                                                   data.frame(cov = rnorm(nrow(grid)))),
                                 control_mcmc = control_mcmc,
                                 type = "joint")

  result <- predict_grid_target(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(poisson_grid)
  expect_setequal(names(result), expected_output)
})

test_that("predict_grid_target handles one-pixel groups in list mode", {
  grid_pred <- list(
    group_one = sf::st_as_sf(data.frame(x = 0, y = 0), coords = c("x", "y"), crs = 4326),
    group_two = sf::st_as_sf(data.frame(x = c(1, 2), y = c(1, 2)), coords = c("x", "y"), crs = 4326)
  )

  object <- list(
    grid_pred = grid_pred,
    S_samples = list(
      matrix(c(1, 2, 3), nrow = 1),
      matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
    ),
    mu_pred = list(c(0), c(0, 0)),
    cov_offset = list(c(0), c(0, 0)),
    re = list(samples = list()),
    par_hat = list(),
    family = "gaussian"
  )
  class(object) <- "RiskMap_pred"

  out <- predict_grid_target(
    object,
    f_target = list(identity_target = function(x) x),
    pd_summary = list(mean = mean, sd = sd)
  )

  expect_equal(out$target$group_one$identity_target$mean, 2)
  expect_equal(out$target$group_one$identity_target$sd, sd(c(1, 2, 3)))
  expect_equal(out$target$group_two$identity_target$mean, c(2, 5))
  expect_equal(out$target$group_two$identity_target$sd, c(sd(c(1, 2, 3)),
                                                          sd(c(4, 5, 6))))
  expect_equal(dim(out$lp_samples[[1]]), c(1, 3))
  expect_equal(dim(out$lp_samples[[2]]), c(2, 3))
})

test_that("plot.RiskMap_predict_grid_target defaults to the first target and validates target/summary #143", {

  gaussian_grid <- setup_prediction(gaussian_model,
                                  grid_pred = grid,
                                  predictors = data.frame(cov = rnorm(nrow(grid))),
                                  re_predictors = data.frame(i = 1:5),
                                  type = "joint")
  result <- predict_grid_target(gaussian_grid)

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

test_that("plot.RiskMap_predict_grid_target draws tiles with lat/lon axes, legend and map annotations #149", {

  gaussian_grid <- setup_prediction(gaussian_model, grid_pred = grid,
                                    predictors = data.frame(cov = rnorm(nrow(grid))),
                                    type = "joint")
  result <- predict_grid_target(gaussian_grid)

  out <- plot(result)
  expect_s3_class(out, "ggplot")
  expect_no_error(ggplot2::ggplot_build(out))
  layer_classes <- vapply(out$layers, function(layer) class(layer$geom)[1], character(1))
  expect_true("GeomTile" %in% layer_classes)
  expect_true("GeomNorthArrow" %in% layer_classes)
  expect_true("GeomScaleBar" %in% layer_classes)
  expect_s3_class(out$coordinates, "CoordSf")
  expect_equal(out$labels$fill, "linear_target_mean")

  out_plain <- plot(result, north_arrow = FALSE, scale_bar = FALSE)
  plain_classes <- vapply(out_plain$layers, function(layer) class(layer$geom)[1], character(1))
  expect_equal(unname(plain_classes), "GeomTile")
})

test_that("plot.RiskMap_predict_grid_target accepts named palettes, colour vectors and scale arguments #149", {

  gaussian_grid <- setup_prediction(gaussian_model, grid_pred = grid,
                                    predictors = data.frame(cov = rnorm(nrow(grid))),
                                    type = "joint")
  result <- predict_grid_target(gaussian_grid)

  expect_no_error(ggplot2::ggplot_build(plot(result, palette = "Blues")))
  expect_no_error(ggplot2::ggplot_build(plot(result, palette = c("white", "darkred"),
                                             name = "Custom title")))

  forward <- ggplot2::ggplot_build(plot(result, palette = c("white", "black")))
  reverse <- ggplot2::ggplot_build(plot(result, palette = c("white", "black"),
                                        reverse_palette = TRUE))
  lowest <- which.min(result$target$linear_target$mean)
  expect_equal(toupper(forward$data[[1]]$fill[lowest]), "#FFFFFF")
  expect_equal(toupper(reverse$data[[1]]$fill[lowest]), "#000000")

  expect_error(plot(result, palette = "not_a_palette"), "'palette' must be one of")
  expect_error(plot(result, palette = c("white", "not_a_colour")), "invalid colours")
  expect_error(plot(result, palette = 1), "'palette' must be a character vector")
  expect_error(plot(result, north_arrow = "yes"), "'north_arrow' must be either TRUE or FALSE")
  expect_error(plot(result, scale_bar = NA), "'scale_bar' must be either TRUE or FALSE")
  expect_error(plot(result, reverse_palette = 1), "'reverse_palette' must be either TRUE or FALSE")
})

test_that("plot.RiskMap_predict_grid_target handles grids with missing columns without warning #149", {

  # regular 1 km grid with a whole column removed
  grid_coordinates <- expand.grid(x = c(0, 1000, 3000, 4000), y = c(0, 1000, 2000))
  gappy_grid <- st_as_sf(grid_coordinates, coords = c("x", "y"), crs = 32637)
  result <- structure(
    list(target = list(linear_target = list(mean = seq_len(nrow(gappy_grid)))),
         grid_pred = gappy_grid,
         f_target = "linear_target",
         pd_summary = "mean"),
    class = "RiskMap_predict_grid_target"
  )

  expect_no_warning(built <- ggplot2::ggplot_build(plot(result)))
  expect_equal(unique(built$data[[1]]$xmax - built$data[[1]]$xmin), 1000)
})
