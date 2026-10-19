test_that("location predictions produce an adaptable prediction map", {
  points <- st_sf(
    response_mean = c(0.1, 0.4, 0.8),
    response_exceedance_0.2 = c(0.1, 0.55, 0.95),
    geometry = st_sfc(
      st_point(c(0, 0)), st_point(c(1, 0)), st_point(c(2, 0)), crs = 32631
    )
  )
  prediction <- structure(list(predictions = points, boundaries = NULL,
                               boundary_summary = NULL),
                          class = "RiskMap_prediction")

  map <- plot(prediction)
  exceedance_map <- plot(prediction, summary = "exceedance", threshold = 0.2)

  expect_s3_class(map, "ggplot")
  expect_no_error(ggplot_build(map))
  expect_null(map$scales$scales[[1L]]$limits)
  expect_equal(exceedance_map$scales$scales[[1L]]$limits, c(0, 1))
  expect_equal(length(prediction_plot_palette(NULL)), 256L)
  expect_equal(
    prediction_plot_palette("-Viridis"),
    rev(grDevices::hcl.colors(256L, "Viridis"))
  )
})

test_that("polygon prediction maps select tidy summaries without reordering", {
  polygons <- st_sf(
    district = c("a", "b"),
    geometry = st_sfc(
      st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))),
      st_polygon(list(rbind(c(1, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 0)))),
      crs = 32631
    )
  )
  summaries <- data.frame(
    boundary_id = rep(1:2, 2),
    aggregate = rep(c("mean", "sum"), each = 2),
    quantity = "response", summary = "mean", threshold = NA_real_,
    value = c(0.2, 0.7, 20, 70)
  )
  point_predictions <- st_sf(
    response_mean = c(0.25, 0.65),
    geometry = st_centroid(st_geometry(polygons))
  )
  prediction <- structure(
    list(
      predictions = point_predictions, boundaries = polygons,
      boundary_summary = summaries
    ),
    class = "RiskMap_prediction"
  )

  point_map <- plot(prediction)
  map <- plot(prediction, aggregate = "mean", limits = c(0, 1))

  expect_s3_class(point_map, "ggplot")
  expect_no_error(ggplot_build(point_map))
  expect_s3_class(map, "ggplot")
  expect_no_error(ggplot_build(map))
  expect_equal(map$layers[[1L]]$data$.prediction_value, c(0.2, 0.7))
  expect_equal(map$scales$scales[[1L]]$limits, c(0, 1))

  expect_error(plot(prediction, aggregate = "median"),
               "must be one of: mean, sum")
})

test_that("plot arguments reject unavailable or ambiguous requests", {
  points <- st_sf(
    response_mean = 0.3,
    response_exceedance_0.2 = 0.7,
    response_exceedance_0.5 = 0.2,
    geometry = st_sfc(st_point(c(0, 0)), crs = 32631)
  )
  prediction <- structure(list(predictions = points, boundaries = NULL,
                               boundary_summary = NULL),
                          class = "RiskMap_prediction")

  expect_error(plot(prediction, quantity = "spatial"), "not available")
  expect_error(plot(prediction, aggregate = "mean"),
               "does not contain boundary summaries")
  expect_error(plot(prediction, summary = "exceedance"),
               "threshold.*must be supplied")
  expect_error(plot(prediction, limits = c(1, 0)), "increasing")
  expect_warning(
    map <- plot(prediction, summary = "exceedance", threshold = 0.2,
                limits = c(0.1, 0.9)),
    "always use limits"
  )
  expect_equal(map$scales$scales[[1L]]$limits, c(0, 1))
})
