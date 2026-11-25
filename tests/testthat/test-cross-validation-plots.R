# plot_cross_validation() and plot_test_sets().

layer_data_of <- function(plot, geom) {
  layers <- vapply(plot$layers, function(layer) class(layer$geom)[1],
                   character(1))
  ggplot2::layer_data(plot, match(geom, layers))
}

test_that("plot_cross_validation() draws a calibration curve for each model", {
  p <- plot_cross_validation(cross_validation)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetWrap")
  expect_identical(levels(p$data$model), c("model_a", "model_b"))
  # The plot() method draws the same.
  expect_identical(ggplot2::layer_data(plot(cross_validation), 3),
                   ggplot2::layer_data(p, 3))
  curves <- cross_validation$calibration$curves
  line <- layer_data_of(p, "GeomLine")
  expect_equal(line$y, curves$curve)
  expect_identical(as.integer(line$PANEL), rep(1:2, each = 101))
  band <- layer_data_of(p, "GeomRibbon")
  expect_equal(band$ymin, curves$lower)
  expect_equal(band$ymax, curves$upper)
  # The areas, with the limit of a correct model.
  area <- cross_validation$calibration$area
  text <- layer_data_of(p, "GeomText")$label
  expect_match(text[1], formatC(area$PIT_area[1], format = "f", digits = 3),
               fixed = TRUE)
  expect_match(text[2], paste("Below", formatC(area$PIT_area_limit[2],
                                               format = "f", digits = 3),
                              "if calibrated"), fixed = TRUE)
  without <- cross_validation
  without$calibration$area$PIT_area_limit <- NA_real_
  expect_no_match(layer_data_of(plot_cross_validation(without), "GeomText")$label,
                  "calibrated")
  # One label for the PIT of Gaussian models and of counts.
  expect_identical(p$labels$y, "Cumulative distribution of the PIT")
  # One model has one panel.
  one <- plot_cross_validation(cross_validation, model = "model_b")
  expect_s3_class(one$facet, "FacetNull")
  expect_identical(nrow(one$data), 101L)
})

test_that("plot_cross_validation() draws PIT histograms", {
  for (bins in c(5, 10, 20)) {
    p <- plot_cross_validation(cross_validation, style = "histogram",
                               bins = bins)
    bars <- layer_data_of(p, "GeomRect")
    expect_identical(nrow(bars), as.integer(2 * bins))
    # The bars of each model have an area of 1.
    expect_equal(as.numeric(tapply(bars$ymax, bars$PANEL, sum)) / bins, c(1, 1))
    # They touch, from 0 to 1.
    expect_equal(bars$xmin[bars$PANEL == 1], seq(0, bins - 1) / bins)
    expect_equal(bars$xmax, bars$xmin + 1 / bins)
    height <- as.numeric(pit_histogram(
      t(cross_validation$calibration$curves$curve[102:202]), bins
    ))
    band <- calibration_band(height, pit_histogram(
      cross_validation$calibration$simulated$model_b, bins
    ))
    steps <- layer_data_of(p, "GeomRibbon")
    expect_equal(steps$ymax[steps$PANEL == 2], rep(band$upper, each = 2))
    expect_equal(steps$ymin[steps$PANEL == 2], rep(band$lower, each = 2))
    expect_equal(
      bars$ymax[bars$PANEL == 2],
      as.numeric(pit_histogram(
        t(cross_validation$calibration$curves$curve[102:202]), bins
      ))
    )
  }
  expect_identical(
    nrow(layer_data_of(plot_cross_validation(cross_validation,
                                             style = "histogram"), "GeomRect")),
    20L
  )
})

test_that("plot_cross_validation() maps the scores on one scale", {
  p <- plot_cross_validation(cross_validation, type = "score")
  expect_identical(levels(p$data$model), c("model_a", "model_b"))
  expect_s3_class(p$facet, "FacetWrap")
  scale <- p$scales$get_scales("fill")
  expect_identical(scale$name, "CRPS (lower is better)")
  scores <- cross_validation$scores
  for (model in c("model_a", "model_b")) {
    expect_equal(sort(p$data$value[p$data$model == model]),
                 sort(scores$CRPS[scores$model == model]))
  }
  # Each value at its own location.
  row <- which(p$data$model == "model_a")[3]
  at <- sf::st_equals(p$data[row, ], gaussian_data)[[1]]
  expect_equal(p$data$value[row],
               scores$CRPS[scores$model == "model_a" & scores$row == at])
  one <- plot_cross_validation(cross_validation, type = "score",
                               score = "SCRPS", model = "model_a")
  expect_s3_class(one$facet, "FacetNull")
  expect_identical(one$scales$get_scales("fill")$name, "SCRPS (lower is better)")
})

test_that("plot_cross_validation() maps the differences from a reference", {
  p <- plot_cross_validation(cross_validation, type = "score",
                             reference = "model_a")
  expect_identical(levels(p$data$model), "model_b")
  scores <- split(cross_validation$scores, cross_validation$scores$model)
  expect_equal(sort(p$data$value),
               sort(scores$model_b$CRPS - scores$model_a$CRPS))
  scale <- p$scales$get_scales("fill")
  expect_identical(scale$name, "Difference in CRPS from model_a")
  expect_equal(sum(scale$limits), 0)
  labels <- scale$labels(scale$breaks)
  expect_match(labels[1], "better than model_a")
  expect_match(labels[5], "worse than model_a")
  expect_match(labels[1], "^\u2212")
  expect_identical(labels[3], "0")
  # Without the reference among the models shown, it is the same map.
  same <- plot_cross_validation(cross_validation, type = "score",
                                reference = "model_a", model = "model_b")
  expect_identical(same$data, p$data)
})

test_that("the scales of the maps end at the 98th percentile", {
  values <- c(seq(0, 1, length.out = 99), 40)
  scale <- score_fill_scale(values, "CRPS")
  expect_lt(max(scale$limits), 1.1)
  expect_match(scale$labels(scale$breaks)[length(scale$breaks)], "or more$")
  # What is beyond has the colour of the end.
  expect_equal(scale$oob(40, scale$limits), max(scale$limits))
  signed <- difference_fill_scale(c(-values, values), "CRPS", "a")
  expect_lte(max(signed$limits), 1)
  expect_match(signed$labels(signed$breaks)[1], "or less\nbetter than a")
  # Equal values and zeros have a scale too.
  expect_no_error(score_fill_scale(rep(2, 4), "CRPS"))
  expect_gt(max(difference_fill_scale(rep(0, 4), "CRPS", "a")$limits), 0)
})

test_that("a location with several scores shows their mean", {
  x <- cross_validation
  x$test_sets$location[2] <- x$test_sets$location[1]
  p <- plot_cross_validation(x, type = "score", model = "model_a")
  expect_identical(nrow(p$data), 9L)
  scores <- x$scores[x$scores$model == "model_a", ]
  expect_true(any(abs(p$data$value -
                        mean(scores$CRPS[scores$row %in% 1:2])) < 1e-12))
})

test_that("the colour of the plots can be changed", {
  expect_identical(formals(plot_cross_validation)$colour, assessment_colour)
  expect_identical(formals(plot_test_sets)$colour, assessment_colour)
  line <- function(...) {
    layer_data_of(plot_cross_validation(cross_validation, ...), "GeomLine")$colour
  }
  expect_identical(unique(line()), assessment_colour)
  expect_identical(unique(line(colour = "firebrick")), "firebrick")
  expect_identical(
    unique(layer_data_of(plot(cross_validation, style = "histogram",
                              colour = "#336699"), "GeomRect")$fill),
    validation_palette("#336699")$fill
  )
  # The paler colours keep the hue, and their lightness is always the same.
  hue <- function(colour) {
    luv <- grDevices::convertColor(t(grDevices::col2rgb(colour)) / 255, "sRGB",
                                   "Luv")
    atan2(luv[3], luv[2])
  }
  palette <- validation_palette("#336699")
  expect_equal(hue(palette$fill), hue("#336699"), tolerance = 0.02)
  lightness <- function(colour) {
    grDevices::convertColor(t(grDevices::col2rgb(colour)) / 255, "sRGB",
                            "Luv")[1]
  }
  for (colour in c("white", "gold", "black", "navy", assessment_colour)) {
    palette <- validation_palette(colour)
    expect_gt(lightness(palette$band), lightness(palette$fill) + 12)
    expect_gt(lightness(palette$fill), lightness(palette$dark) + 12)
  }
  sets <- create_test_sets(gaussian_data)
  expect_identical(
    unique(ggplot2::layer_data(plot_test_sets(sets, colour = "navy"))$colour),
    "navy"
  )
  for (colour in list("no colour", 2, NA, c("red", "blue"), NULL)) {
    expect_input_error(plot_cross_validation(cross_validation, colour = colour),
                       "colour")
    expect_input_error(plot_test_sets(sets, colour = colour), "colour")
  }
  expect_input_error(plot_cross_validation(cross_validation, type = "score",
                                           colour = "red"), "colour",
                     "type = \"calibration\"")
})

test_that("plot_cross_validation() draws without a warning", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (arguments in list(list(), list(style = "histogram"),
                         list(type = "score"),
                         list(type = "score", reference = "model_a"))) {
    expect_no_warning(print(do.call(plot_cross_validation,
                                    c(list(cross_validation), arguments))))
  }
})

test_that("plot_cross_validation() checks its arguments", {
  draw <- function(...) plot_cross_validation(cross_validation, ...)
  expect_input_error(plot_cross_validation(1), "x")
  expect_error(draw(type = "PIT"), class = "riskmap_error_input")
  expect_error(draw(style = "bars"), class = "riskmap_error_input")
  expect_error(draw(type = "score", score = "PIT"),
               class = "riskmap_error_input")
  expect_input_error(draw(model = "model_c"), "model", "\"model_c\"")
  expect_input_error(draw(type = "score", reference = "model_c"), "reference")
  expect_input_error(draw(type = "score", reference = "model_a",
                          model = "model_a"), "model", "other than")
  expect_input_error(draw(score = "CRPS"), "score", "type = \"score\"")
  expect_input_error(draw(reference = "model_a"), "reference",
                     "type = \"score\"")
  expect_input_error(draw(type = "score", style = "curve"), "style",
                     "type = \"calibration\"")
  expect_input_error(draw(bins = 10), "bins", "style = \"histogram\"")
  for (bins in list(7, "ten", c(5, 10), NULL, NA)) {
    expect_input_error(draw(style = "histogram", bins = bins), "bins",
                       "5, 10 or 20")
  }
  # The former arguments.
  expect_input_error(draw(metric = "CRPS"), "metric", "type = \"score\"")
  for (name in c("pit_mode", "pit_test_set", "test_set")) {
    call <- c(list(cross_validation), stats::setNames(list("all"), name))
    expect_input_error(do.call(plot, call), name, "all test observations")
  }
  expect_input_error(draw(combine_pit = TRUE), "combine_pit", "own panel")
  expect_input_error(draw(linewidth = 2), "linewidth")
})

plot_test_data <- function() {
  set.seed(60)
  st_as_sf(data.frame(x = runif(40, 0, 1e5), y = runif(40, 0, 1e5)),
           coords = c("x", "y"), crs = 32637)
}

test_that("plot_test_sets() maps the test sets of every method", {
  data <- plot_test_data()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw <- function(sets) {
    expect_no_warning(print(p <- plot_test_sets(sets)))
    expect_s3_class(p, "ggplot")
    p
  }
  # The clusters together, with their numbers.
  clusters <- create_test_sets(data, "cluster", n_sets = 4, seed = 1)
  p <- draw(clusters)
  expect_identical(nrow(p$layers[[1]]$data), 4L)
  expect_identical(nrow(p$layers[[2]]$data), 40L)
  expect_identical(p$layers[[4]]$data$label, 1:4)
  expect_s3_class(p$facet, "FacetNull")
  expect_identical(p$labels, plot(clusters)$labels)
  # A panel for each repeat.
  twice <- create_test_sets(data, "cluster", n_sets = 4, repeats = 2, seed = 1)
  p <- draw(twice)
  expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
  expect_identical(p$layers[[4]]$data$label, rep(1:4, 2))
  # Clusters of one or two locations, or on a line, have no area.
  draw(create_test_sets(data[1:6, ], "cluster", n_sets = 5, seed = 1))
  draw(create_test_sets(data[1:6, ], "cluster", n_sets = 6, seed = 1))
  on_line <- st_as_sf(data.frame(x = 1:8 * 1e4, y = 0), coords = c("x", "y"),
                      crs = 32637)
  draw(create_test_sets(on_line, "cluster", n_sets = 2, seed = 1))
  # Leave-one-out, the locations.
  p <- draw(create_test_sets(data))
  expect_identical(nrow(p$data), 40L)
  expect_match(p$labels$subtitle, "Each of the 40 locations")
  # A panel for each of the others, the first twelve.
  p <- draw(create_test_sets(data, "random", n_sets = 4, seed = 1))
  expect_identical(nrow(p$data), 160L)
  expect_identical(as.integer(table(p$data$role)), c(120L, 40L))
  p <- draw(create_test_sets(data, "spaced", n_sets = 15, size = 3,
                             min_dist = 10, seed = 1))
  expect_identical(nlevels(p$data$set), 12L)
  expect_identical(p$labels$subtitle, "The first 12 of 15 test sets")
  p <- draw(create_test_sets(data, "spaced", n_sets = 1, size = 3,
                             min_dist = 10, seed = 1))
  expect_s3_class(p$facet, "FacetNull")
  draw(read_test_sets(list(1:3, 4:6), gaussian_intercept_model))
  # In longitude and latitude.
  draw(suppressMessages(create_test_sets(st_transform(data, 4326), "cluster",
                                         n_sets = 4, seed = 1)))
})

test_that("neighbouring clusters have different colours", {
  sets <- create_test_sets(plot_test_data(), "cluster", n_sets = 9, seed = 1)
  p <- plot_test_sets(sets)
  colour <- p$layers[[2]]$data$colour
  set <- rep(seq_along(sets$sets), lengths(lapply(sets$sets, `[[`, "test")))
  rows <- unlist(lapply(sets$sets, `[[`, "test"))
  colour <- tapply(colour, set, unique)
  centres <- t(vapply(seq_along(sets$sets), function(k) {
    colMeans(sets$coordinates[rows[set == k], , drop = FALSE])
  }, numeric(2)))
  distance <- as.matrix(dist(centres))
  diag(distance) <- Inf
  nearest <- apply(distance, 1, which.min)
  expect_true(all(colour != colour[nearest]))
})

test_that("plot_test_sets() draws the distances from the training data", {
  data <- plot_test_data()
  loo <- create_test_sets(data)
  clusters <- create_test_sets(data, "cluster", n_sets = 4, seed = 1)
  p <- plot_test_sets(clusters, type = "distance")
  expect_identical(levels(p$data$what), "Test locations")
  expect_equal(p$data$distance[-1], sort(test_set_distances(clusters)$distance))
  expect_equal(p$data$below, (0:40) / 40)
  expect_match(p$labels$x, "(km)", fixed = TRUE)

  grid <- create_grid(create_convex_hull(data), spacing = 10)
  p <- plot_test_sets(list("Leave-one-out" = loo, Clusters = clusters),
                      type = "distance", grid = grid)
  expect_identical(levels(p$data$what),
                   c("Prediction locations", "Leave-one-out", "Clusters"))
  direct <- st_distance(grid, data)
  expected <- sort(as.numeric(apply(direct, 1, min))) / 1000
  expect_equal(p$data$distance[p$data$what == "Prediction locations"][-1],
               expected)
  # The locations in another CRS are at the same distances.
  elsewhere <- plot_test_sets(loo, type = "distance",
                              grid = st_transform(grid, 4326))
  expect_equal(
    elsewhere$data$distance[elsewhere$data$what == "Prediction locations"],
    p$data$distance[p$data$what == "Prediction locations"], tolerance = 1e-5
  )
  # Two with one name would be one curve.
  expect_input_error(plot_test_sets(list(A = loo, A = clusters),
                                    type = "distance"), "x",
                     "\"A\" is repeated")
  expect_input_error(plot_test_sets(list("Prediction locations" = loo),
                                    type = "distance", grid = grid), "x",
                     "The curve of `grid` is called")
  expect_input_error(plot_test_sets(list("Test sets 2" = loo, clusters),
                                    type = "distance"), "x", "is repeated")
  # Unnamed sets have numbers.
  p <- plot_test_sets(list(loo, clusters), type = "distance")
  expect_identical(levels(p$data$what), c("Test sets 1", "Test sets 2"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(print(p))
})

test_that("plot_test_sets() checks its arguments", {
  data <- plot_test_data()
  loo <- create_test_sets(data)
  expect_input_error(plot_test_sets(), "x")
  expect_input_error(plot_test_sets(1), "x")
  expect_input_error(plot_test_sets(list(loo, 1), type = "distance"), "x")
  expect_input_error(plot_test_sets(list(), type = "distance"), "x")
  expect_input_error(plot_test_sets(rep(list(loo), 6), type = "distance"), "x",
                     "at most 5")
  expect_input_error(plot_test_sets(loo, grid = data), "grid",
                     "type = \"distance\"")
  expect_input_error(plot_test_sets(loo, type = "distance", grid = 1), "grid")
  expect_input_error(plot_test_sets(loo, type = "distance",
                                    grid = st_buffer(data, 10)), "grid")
  expect_error(plot_test_sets(loo, type = "histogram"),
               class = "riskmap_error_input")
  expect_input_error(plot_test_sets(loo, size = 2), "size")
  metres <- loo
  metres$distance_units <- "m"
  expect_input_error(plot_test_sets(list(loo, metres), type = "distance"), "x",
                     "same units of distance")
})
