test_that("extrapolation indicators separate spatial and covariate novelty", {
  locations <- st_coordinates(gaussian_data)
  cov_range <- range(gaussian_data$cov)
  centre <- colMeans(locations)
  newdata <- st_as_sf(
    data.frame(
      x = c(locations[1, 1], centre[[1]], locations[1, 1], 1e7, 5e6),
      z = c(locations[1, 2], centre[[2]], locations[1, 2], 1e7, 5e6),
      cov = c(gaussian_data$cov[1], mean(cov_range),
              cov_range[2] + diff(cov_range), cov_range[1], NA)
    ),
    coords = c("x", "z"), crs = 32637
  )
  expect_warning(
    extrapolation <- assess_extrapolation(gaussian_model, newdata),
    "missing required values"
  )

  expect_s3_class(extrapolation, c("RiskMap_extrapolation", "sf"))
  expect_identical(nrow(extrapolation), 5L)
  expect_equal(extrapolation$covariates_outside[1:4], c(0L, 0L, 1L, 0L))
  expect_equal(extrapolation$covariate_novelty[3], 1)
  expect_identical(extrapolation$most_novel_covariate[3], "cov")
  expect_true(is.na(extrapolation$most_novel_covariate[1]))
  expect_lt(extrapolation$spatial_unexplained[1],
            extrapolation$spatial_unexplained[2])
  expect_equal(extrapolation$spatial_unexplained[4], 1)
  expect_identical(
    as.character(extrapolation$extrapolation[c(1, 3, 4)]),
    c("Supported", "Covariate", "Spatial")
  )
  expect_true(all(is.na(st_drop_geometry(extrapolation)[5, c(
    "covariates_outside", "covariate_novelty", "covariate_leverage",
    "spatial_unexplained", "extrapolation"
  )])))
  expect_lte(max(extrapolation$covariate_leverage[1:2]), 1 + 1e-8)
  expect_gt(extrapolation$covariate_leverage[3], 1)
})

test_that("unexplained spatial variance matches the kriging formula", {
  parameters <- coef(gaussian_intercept_model)
  observed <- gaussian_intercept_model$coords
  ratio <- parameters$sigma2_me / parameters$sigma2
  correlation <- matern_correlation(as.matrix(dist(observed)),
                                    phi = parameters$phi,
                                    kappa = gaussian_intercept_model$kappa)
  diag(correlation) <- 1 + ratio
  point <- observed[1, , drop = FALSE] + c(0.5, -0.3)
  cross <- matern_correlation(
    sqrt(colSums((t(observed) - as.numeric(point))^2)),
    phi = parameters$phi, kappa = gaussian_intercept_model$kappa
  )
  expected <- 1 - drop(cross %*% solve(correlation, cross))

  factor <- crs_to_distance_factor(gaussian_data,
                                   gaussian_intercept_model$distance_units)
  newdata <- st_as_sf(data.frame(x = point[1] / factor, z = point[2] / factor),
                      coords = c("x", "z"), crs = 32637)
  extrapolation <- assess_extrapolation(gaussian_intercept_model, newdata)
  expect_equal(extrapolation$spatial_unexplained, expected, tolerance = 1e-8)
  expect_identical(extrapolation$covariates_outside, 0L)
})

test_that("extrapolation maps and argument checks", {
  boundary <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
    c(0, 0), c(10000, 0), c(10000, 10000), c(0, 10000), c(0, 0)
  ))), crs = 32637))
  grid <- create_grid(boundary, spacing = 2, distance_units = "km")
  grid$cov <- seq(-3, 3, length.out = nrow(grid))
  extrapolation <- assess_extrapolation(gaussian_model, grid)
  expect_false(is.null(attr(extrapolation, "grid_metadata")))
  for (indicator in c("extrapolation", "spatial_unexplained",
                      "covariate_novelty", "covariate_leverage")) {
    map <- plot(extrapolation, indicator = indicator)
    expect_s3_class(map, "ggplot")
    expect_s3_class(map$layers[[1L]]$geom, "GeomSf")
  }
  expect_identical(ggplot2::get_labs(plot(extrapolation))$fill,
                   "Extrapolation")
  with_missing <- grid
  with_missing$cov[1] <- NA
  expect_warning(missing_map <- plot(assess_extrapolation(gaussian_model,
                                                          with_missing)),
                 "missing required values")
  expect_identical(ggplot2::get_guide_data(missing_map, "fill")$.label,
                   c("Supported", "Spatial", "Covariate",
                     "Spatial and covariate"))
  expect_true("grey90" %in% ggplot2::ggplot_build(missing_map)$data[[1]]$fill)
  skip_if_not_installed("patchwork")
  both <- plot(extrapolation,
               indicator = c("spatial_unexplained", "extrapolation"))
  expect_length(c(both$patches$plots, list(both)), 2L)
  expect_identical(both$patches$layout$ncol, 2L)
  expect_error(plot(extrapolation, indicator = "distance"))
  expect_error(assess_extrapolation(gaussian_model, grid,
                                    spatial_threshold = 2))
  expect_error(assess_extrapolation(list(), grid), "fitted RiskMap model")
})
