legend_title <- function(map) ggplot2::get_labs(map)$fill

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
  # Points with missing aesthetics are only dropped when the map is drawn.
  expect_no_warning(ggplot_gtable(ggplot_build(map)))
  expect_equal(map$scales$scales[[1L]]$limits, c(0, 0.8))
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

test_that("map defaults depend on what the map shows", {
  points <- st_sf(
    response_mean = c(0.05, 0.2, 0.35),
    response_cv = c(0.2, 0.3, 0.5),
    response_exceedance_0.2 = c(0.1, 0.5, 0.9),
    spatial_mean = c(-0.4, 0.1, 0.8),
    geometry = st_sfc(
      st_point(c(0, 0)), st_point(c(1, 0)), st_point(c(2, 0)), crs = 32631
    )
  )
  prediction <- structure(
    list(predictions = points, boundaries = NULL, boundary_summary = NULL,
         family = "binomial", summaries = c("mean", "cv", "exceedance"),
         thresholds = 0.2, level = 0.95),
    class = "RiskMap_prediction"
  )

  mean_scale <- plot(prediction)$scales$scales[[1L]]
  cv_scale <- plot(prediction, summary = "cv")$scales$scales[[1L]]
  exceedance_scale <- plot(prediction,
                           summary = "exceedance")$scales$scales[[1L]]
  spatial_scale <- plot(prediction, quantity = "spatial")$scales$scales[[1L]]

  expect_identical(legend_title(plot(prediction)), "Mean prevalence")
  expect_identical(legend_title(plot(prediction, summary = "cv")),
                   "CV of prevalence")
  expect_identical(legend_title(plot(prediction, summary = "exceedance")),
                   "Probability prevalence > 20%")
  expect_identical(legend_title(plot(prediction, quantity = "spatial")),
                   "Mean spatial effect")
  expect_identical(mean_scale$get_labels(), c("0%", "10%", "20%", "30%",
                                              "40%"))
  expect_identical(exceedance_scale$get_labels(),
                   c("0%", "25%", "50%", "75%", "100%"))
  expect_identical(cv_scale$get_labels(), c("0.2", "0.3", "0.4", "0.5"))
  expect_identical(spatial_scale$limits, c(-1, 1))
  expect_identical(spatial_scale$get_labels()[[1L]], "\u22121")
  expect_identical(prediction_plot_role("response", "interval_width"),
                   "uncertainty")
  expect_false(identical(prediction_role_palette("level"),
                         prediction_role_palette("uncertainty")))
})

test_that("legend titles name the response of each family", {
  base <- list(level = 0.9)
  poisson <- structure(c(base, family = "poisson"),
                       class = "RiskMap_prediction")
  gaussian <- structure(c(base, family = "gaussian", response_name = "haz"),
                        class = "RiskMap_prediction")

  expect_identical(prediction_legend_title(poisson, "response", "mean", NULL),
                   "Mean rate")
  expect_identical(prediction_legend_title(gaussian, "response", "upper",
                                           NULL),
                   "Upper 90% limit of haz")
  expect_identical(prediction_legend_title(gaussian, "response",
                                           "exceedance", 1.5),
                   "Probability haz > 1.5")
  expect_identical(prediction_legend_title(gaussian, "stunting", "mean",
                                           NULL, aggregate = "sum"),
                   "Mean stunting (area sum)")
})

test_that("both spatial scales share one colour scale", {
  polygons <- st_sf(
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
    value = c(0.3, 0.5, 30, 50)
  )
  points <- st_sf(
    response_mean = c(0.1, 0.9),
    geometry = st_centroid(st_geometry(polygons))
  )
  prediction <- structure(
    list(predictions = points, boundaries = polygons,
         boundary_summary = summaries, family = "binomial",
         summaries = "mean", aggregate = c("mean", "sum")),
    class = "RiskMap_prediction"
  )

  default_map <- plot(prediction)
  location_map <- plot(prediction, spatial_scale = "location")
  both_map <- plot(prediction, spatial_scale = "both")

  expect_identical(default_map$layers[[1L]]$data$.prediction_value,
                   c(0.3, 0.5))
  expect_length(location_map$layers, 2L)
  expect_s3_class(both_map$facet, "FacetWrap")
  expect_identical(levels(both_map$layers[[1L]]$data$.spatial_scale),
                   c("Locations", "Area mean"))
  expect_no_error(ggplot_build(both_map))
  expect_identical(both_map$scales$scales[[1L]]$limits, c(0, 1))
  expect_setequal(
    unlist(lapply(both_map$layers, function(layer) {
      unique(as.character(layer$data$.spatial_scale))
    })),
    c("Locations", "Area mean")
  )
  expect_identical(legend_title(plot(prediction, aggregate = "sum")),
                   "Mean prevalence (area sum)")
  expect_error(plot(prediction, spatial_scale = "both", aggregate = "sum"),
               "same scale as the locations")
})

test_that("breaks give a binned map and several summaries give one row each", {
  points <- st_sf(
    response_mean = c(0.05, 0.2, 0.35),
    response_cv = c(0.2, 0.3, 0.5),
    geometry = st_sfc(
      st_point(c(0, 0)), st_point(c(1, 0)), st_point(c(2, 0)), crs = 32631
    )
  )
  prediction <- structure(
    list(predictions = points, boundaries = NULL, boundary_summary = NULL,
         family = "binomial", summaries = c("mean", "cv")),
    class = "RiskMap_prediction"
  )

  binned <- plot(prediction, breaks = c(0, 0.1, 0.2, 0.4))

  expect_s3_class(binned$scales$scales[[1L]], "ScaleBinned")
  expect_no_error(ggplot_build(binned))
  expect_error(plot(prediction, breaks = c(0, 0.2)), "at least three")
  expect_error(plot(prediction, breaks = c(0, 0.3, 0.2)), "increasing")

  skip_if_not_installed("patchwork")
  both <- plot(prediction, summary = c("mean", "cv"))
  expect_s3_class(both, "patchwork")
  expect_length(both$patches$plots, 1L)
})

test_that("predictions from a fitted model print and plot at both scales", {
  areas <- st_sf(
    geometry = st_make_grid(gaussian_data, n = c(2, 1))
  )
  prediction <- predict(gaussian_intercept_model, newdata = gaussian_data,
                        boundaries = areas, aggregate = "mean",
                        seed = 1, messages = FALSE)

  expect_identical(prediction$response_name, "y")
  expect_output(print(prediction), "Map with plot\\(x, summary = \"mean\"\\)")
  expect_identical(legend_title(plot(prediction)), "Mean y")
  expect_no_error(ggplot_build(plot(prediction, spatial_scale = "both")))
})

test_that("grid maps draw clipped cells and accept extra sf layers", {
  boundary <- st_sf(geometry = st_sfc(
    st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 0)))),
    crs = 32631
  ))
  grid <- create_grid(boundary, spacing = 2, distance_units = "m",
                      keep = "overlap")
  grid$response_mean <- seq(0.1, 0.9, length.out = nrow(grid))
  grid$fixed <- seq(-2, -1, length.out = nrow(grid))
  prediction <- structure(
    list(predictions = grid, boundaries = NULL, boundary_summary = NULL,
         family = "binomial", summaries = "mean"),
    class = "RiskMap_prediction"
  )

  map <- plot(prediction)
  cells <- map$layers[[1L]]$data
  expect_s3_class(map$layers[[1L]]$geom, "GeomSf")
  expect_true(all(as.numeric(st_area(cells)) <= 4 + 1e-8))
  expect_true(all(lengths(st_within(cells, st_buffer(boundary, 1e-6))) > 0))
  expect_length(map$layers, 2L)
  outline <- map$layers[[2L]]$data
  expect_equal(as.numeric(st_area(outline)), as.numeric(st_area(boundary)))
  expect_true(is.na(map$layers[[2L]]$aes_params$fill))

  points <- st_transform(st_sf(geometry = st_sfc(st_point(c(5, 2)),
                                                 crs = 32631)), 4326)
  with_points <- map + geom_sf(data = points, inherit.aes = FALSE)
  built <- ggplot_build(with_points)
  cell_box <- st_bbox(built$data[[1L]]$geometry)
  point_box <- st_bbox(built$data[[3L]]$geometry)
  expect_true(point_box[["xmin"]] >= cell_box[["xmin"]] &&
                point_box[["xmax"]] <= cell_box[["xmax"]])
  expect_true(all(st_geometry_type(cells) %in%
                    c("POLYGON", "MULTIPOLYGON")))

  relabelled <- map + labs(fill = "Prevalence in children")
  expect_identical(legend_title(relabelled), "Prevalence in children")

  fixed_map <- plot(prediction, quantity = "fixed")
  expect_identical(legend_title(fixed_map), "Fixed effects")
  expect_identical(fixed_map$scales$scales[[1L]]$limits, c(-2, 2))
  fixed_limited <- plot(prediction, quantity = "fixed", limits = c(-3, 3))
  expect_identical(fixed_limited$scales$scales[[1L]]$limits, c(-3, 3))
})

test_that("several summaries at one scale are arranged in a grid", {
  skip_if_not_installed("patchwork")
  points <- st_sf(
    response_mean = c(0.05, 0.2, 0.35),
    response_cv = c(0.2, 0.3, 0.5),
    response_interval_width = c(0.1, 0.2, 0.3),
    geometry = st_sfc(
      st_point(c(0, 0)), st_point(c(1, 0)), st_point(c(2, 0)), crs = 32631
    )
  )
  prediction <- structure(
    list(predictions = points, boundaries = NULL, boundary_summary = NULL,
         family = "binomial", summaries = c("mean", "cv")),
    class = "RiskMap_prediction"
  )
  maps <- plot(prediction, summary = c("mean", "cv", "interval_width"))
  expect_identical(maps$patches$layout$ncol, 2L)
  side_margins <- vapply(c(maps$patches$plots, list(maps)), function(map) {
    as.numeric(map$theme$plot.margin[2L])
  }, numeric(1))
  expect_true(all(side_margins == 14))
})

test_that("several quantities are drawn in one call", {
  skip_if_not_installed("patchwork")
  points <- st_sf(
    response_mean = c(0.05, 0.2, 0.35),
    spatial_mean = c(-0.5, 0.1, 0.6),
    fixed = c(-0.3, 0, 0.3),
    geometry = st_sfc(
      st_point(c(0, 0)), st_point(c(1, 0)), st_point(c(2, 0)), crs = 32631
    )
  )
  prediction <- structure(
    list(predictions = points, boundaries = NULL, boundary_summary = NULL,
         family = "binomial", summaries = "mean"),
    class = "RiskMap_prediction"
  )
  maps <- plot(prediction, quantity = c("fixed", "spatial"),
               limits = c(-1, 1))
  panels <- c(maps$patches$plots, list(maps))
  expect_length(panels, 2L)
  expect_identical(maps$patches$layout$ncol, 2L)
  expect_identical(vapply(panels, legend_title, character(1)),
                   c("Fixed effects", "Mean spatial effect"))
  for (panel in panels) {
    expect_identical(panel$scales$scales[[1L]]$limits, c(-1, 1))
  }
  expect_error(plot(prediction, quantity = c("fixed", "spatial"),
                    summary = c("mean", "sd")),
               "several quantities or several summaries")
  expect_error(plot(prediction, quantity = c("spatial", "spatial")),
               "unique quantity names")
})
