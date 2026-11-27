# Plots of test sets and of cross-validations.

##' The default colour of the plots that assess a model
##' @noRd
assessment_colour <- "#00857F"

##' Colours shared by the plots that assess a model
##'
##' `line` is the colour of what is assessed, the data or the test locations, and
##' `fill` and `band` are paler and greyer versions of it, for bars and for the
##' band of a correct model around a curve. `limits` is the grey between the
##' limits of bars, `contrast` a second colour, and `sets` colours for groups. Readers with any of the usual colour
##' deficiencies can tell the default ones apart.
##' @param colour A colour, as a name or a hexadecimal code.
##' @noRd
validation_palette <- function(colour = assessment_colour) {
  list(line = colour, fill = paler_colour(colour, 61, 0.7),
       band = paler_colour(colour, 93, 0.25), limits = "grey95",
       contrast = "#D55E00", reference = "#7a7f87", dark = "#1d232b",
       text = "#3b4350", frame = "#d5dae0", training = "#c3cad2",
       sets = c(colour, "#D55E00", "#56B4E9", "#3B4A9B", "#E69F00"))
}

##' A paler and greyer version of a colour
##'
##' Its lightness doesn't depend on that of the colour, so that bars and bands
##' can be told apart whether the colour is light or dark.
##' @param lightness The lightness of the result, from 0 to 100.
##' @param chroma The colourfulness of the result, relative to `colour`.
##' @noRd
paler_colour <- function(colour, lightness, chroma) {
  luv <- grDevices::convertColor(t(grDevices::col2rgb(colour)) / 255,
                                 "sRGB", "Luv")
  luv <- cbind(lightness, chroma * luv[, 2:3, drop = FALSE])
  grDevices::rgb(grDevices::convertColor(luv, "Luv", "sRGB", clip = TRUE))
}

##' Check a colour given for a plot
##' @noRd
check_colour <- function(colour, call = rlang::caller_env()) {
  valid <- length(colour) == 1L && is.character(colour) && !is.na(colour) &&
    tryCatch(is.matrix(grDevices::col2rgb(colour)),
             error = function(error) FALSE)
  if (!valid) {
    abort_input("{.arg colour} must be a colour, as a name or a hexadecimal code, not {describe_value(colour)}.", arg = "colour", call = call)
  }
  invisible(colour)
}

##' The theme of the calibration and distance panels
##' @noRd
validation_theme <- function() {
  profile_theme() +
    theme(strip.text = element_text(face = "bold", size = 10, hjust = 0,
                                    margin = margin(2, 0, 4, 0)),
          panel.spacing = unit(14, "pt"),
          plot.margin = margin(4, 10, 4, 4),
          axis.title = element_text(size = 9.5),
          axis.text = element_text(size = 8),
          legend.position = "top")
}

##' The theme of the maps, without axes and with a light frame
##' @noRd
validation_map_theme <- function() {
  theme_void(base_size = 11) +
    theme(strip.text = element_text(face = "bold", size = 10, hjust = 0,
                                    margin = margin(2, 0, 3, 0)),
          panel.spacing = unit(8, "pt"),
          panel.border = element_rect(colour = validation_palette()$frame,
                                      fill = NA, linewidth = 0.4),
          plot.subtitle = element_text(colour = "#586271", size = 9,
                                       margin = margin(0, 0, 4, 0)),
          legend.position = "top",
          plot.margin = margin(4, 6, 4, 6))
}

##' @title Plot a cross-validation
##'
##' @description
##' Draws the calibration of the models of a cross-validation by
##' [assess_prediction()], or maps their scores, as `plot()` of the result
##' does.
##'
##' @param x The result of [assess_prediction()].
##' @param type What to draw, the `"calibration"` of each model (the default)
##'   or a map of a `"score"`.
##' @param style For the calibration, a `"curve"` (the default) or a
##'   `"histogram"`.
##' @param score For a map, the score to show, `"CRPS"` (the default) or
##'   `"SCRPS"`.
##' @param reference For a map, the name of a model. The other models are then
##'   mapped as differences from its scores.
##' @param model Names of the models to include, all of them by default.
##' @param bins The number of bins of a histogram, 5, 10 (the default) or 20.
##' @param colour For the calibration, the colour of the curves. The histograms
##'   and the bands have paler versions of it. The maps of scores have scales of
##'   their own and don't take it.
##' @param ... For `plot()`, the arguments of `plot_cross_validation()`. For
##'   `plot_cross_validation()`, not used, and an error if given, so that a
##'   misspelt argument is not ignored. Add \pkg{ggplot2} components to the
##'   plot with `+`.
##'
##' @details
##' # Calibration
##'
##' The curve is the cumulative distribution of the PIT values. For Gaussian
##' models it is their empirical distribution function, and for counts the
##' average nonrandomised PIT (AnPIT), which is the average of the empirical
##' distribution functions over the randomisation of the PIT. The curve of a
##' calibrated model stays close to the dashed identity line. An S-shaped curve
##' means that the predictive distributions are too narrow or too wide, and a
##' curve on one side of the line that they are too high or too low.
##'
##' The histogram shows the kind of departure. Its bars are the increases of
##' the curve over each bin, divided by the width of the bin, which is the
##' nonrandomised PIT histogram of Czado *et al.* (2009). A U shape means that
##' the predictive distributions are too narrow, a hump that they are too
##' wide, and a slope that they are too high or too low. The line at 1 is a
##' perfectly calibrated model.
##'
##' The shaded band of the curves, and the area between the dashed lines of the
##' histograms, are the 95% joint band of a correct model with the same test
##' sets. The whole curve, or all the bars together, of a correct model is inside
##' it with a probability of 95%.
##'
##' Each panel gives the area between the curve and the identity line, `PIT_area`,
##' and, after "Below", `PIT_area_limit`. It is a 95% limit and not a certainty. The
##' area of a correct model is below it with a probability of 95%, and above it
##' for one dataset in twenty. [assess_prediction()] says how the band and the
##' limit are computed, how to read the area, and where they are not a check of
##' calibration, as for binary outcomes.
##'
##' # Scores
##'
##' A location in several test sets, or with several observations, shows the
##' mean of its scores. The panels of the models share one colour scale, on
##' which darker colours are worse scores. Few locations often have scores far
##' from the others, and then the maps of the models look alike, and their
##' differences from a `reference` tell more. Negative differences are
##' locations where the model of the panel predicts better than the reference,
##' in green, and positive ones, in orange, where it predicts worse.
##' Both scales end at the 98th percentile, of the scores or of the sizes of
##' the differences, and what is beyond it has the colour of the end.
##'
##' @return A \pkg{ggplot2} object, to which components can be added.
##'
##' @references
##' Czado, C., Gneiting, T., & Held, L. (2009). Predictive model assessment
##' for count data. *Biometrics*, 65(4), 1254–1261.
##' \doi{10.1111/j.1541-0420.2009.01191.x}.
##'
##' @seealso [assess_prediction()], and [plot_test_sets()] for the test sets.
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:100, ]
##' models <- list(
##'   intercept = glgpm(y ~ gp(), data = example_data, family = "gaussian",
##'                     messages = FALSE),
##'   covariate = glgpm(y ~ log(pop_dens) + gp(), data = example_data,
##'                     family = "gaussian", messages = FALSE)
##' )
##' test_sets <- create_test_sets(example_data, method = "cluster",
##'                               n_sets = 5, seed = 1)
##' result <- assess_prediction(models, test_sets, seed = 1, messages = FALSE)
##' # The calibration of each model, with the 95% band of a correct model. The
##' # area of a correct model is below the number given with a probability of 95%.
##' plot_cross_validation(result)
##' plot_cross_validation(result, style = "histogram")
##' # Where the model with the covariate predicts better or worse
##' plot_cross_validation(result, type = "score", reference = "intercept")
##' @export
plot_cross_validation <- function(x, type = c("calibration", "score"),
                                  style = c("curve", "histogram"),
                                  score = c("CRPS", "SCRPS"),
                                  reference = NULL, model = NULL, bins = 10,
                                  colour = "#00857F", ...) {
  check_supplied(x)
  check_plot_object(x, "RiskMap_cross_validation", "plot_cross_validation",
                    "assess_prediction")
  check_dots_unused(..., .hints = c(
    metric = "Use {.code type = \"calibration\"} for the calibration curve, or {.code type = \"score\"} with {.arg score}.",
    pit_mode = "The curve is that of all test observations together.",
    pit_test_set = "The curve is that of all test observations together.",
    combine_pit = "Each model has its own panel.",
    test_set = "The curve is that of all test observations together.",
    combine_panels = "Each model has its own panel."
  ))
  given <- c(style = !missing(style), score = !missing(score),
             bins = !missing(bins), colour = !missing(colour))
  type <- match_choice(type)
  style <- match_choice(style)
  score <- match_choice(score)
  if (!identical(bins %in% c(5, 10, 20), TRUE)) {
    abort_input("{.arg bins} must be 5, 10 or 20, not {describe_value(bins)}.", arg = "bins")
  }
  check_colour(colour)
  models <- x$calibration$area$model
  check_assessed_model(model, models, "model", several = TRUE)
  check_assessed_model(reference, models, "reference")
  check_plot_type(type, style, given, reference)
  shown <- if (is.null(model)) models else models[models %in% model]

  if (type == "calibration") {
    plot_calibration(x, shown, style, bins, colour)
  } else {
    plot_scores(x, shown, score, reference)
  }
}

##' @rdname plot_cross_validation
##' @method plot RiskMap_cross_validation
##' @export
plot.RiskMap_cross_validation <- function(x, ...) plot_cross_validation(x, ...)

##' Stop at arguments that the type of plot doesn't have
##' @noRd
check_plot_type <- function(type, style, given, reference,
                            call = rlang::caller_env()) {
  if (type == "calibration" && (given[["score"]] || !is.null(reference))) {
    arg <- if (given[["score"]]) "score" else "reference"
    abort_input("{.arg {arg}} is for maps, with {.code type = \"score\"}.", arg = arg, call = call)
  }
  if (type == "score" && (given[["style"]] || given[["colour"]])) {
    arg <- if (given[["style"]]) "style" else "colour"
    abort_input("{.arg {arg}} is for {.code type = \"calibration\"}.", arg = arg, call = call)
  }
  if (given[["bins"]] && (type == "score" || style == "curve")) {
    abort_input("{.arg bins} is for {.code style = \"histogram\"}.", arg = "bins", call = call)
  }
  invisible(TRUE)
}

##' The calibration panels of some models of a cross-validation
##' @noRd
plot_calibration <- function(x, models, style, bins, colour) {
  area <- x$calibration$area
  area <- area[match(models, area$model), ]
  area <- area[is.finite(area$PIT_area), ]
  if (!nrow(area)) {
    abort_input("{.arg x} has no scores to assess the calibration with.", arg = "x")
  }
  area$label <- paste0("Area ", format_area(area$PIT_area))
  # A limit needs at least 19 simulated curves.
  limited <- is.finite(area$PIT_area_limit)
  area$label[limited] <- paste0(
    area$label[limited], "\nBelow ", format_area(area$PIT_area_limit[limited]),
    " if calibrated"
  )
  area$model <- factor(area$model, levels = models)
  colours <- validation_palette(colour)
  out <- if (style == "curve") {
    plot_calibration_curves(x, models, area, colours)
  } else {
    plot_calibration_histograms(x, models, area, bins, colours)
  }
  out <- out +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    labs(x = "PIT") + validation_theme() +
    theme(panel.border = element_rect(colour = colours$frame, fill = NA,
                                      linewidth = 0.4))
  if (style == "histogram") {
    out <- out + theme(aspect.ratio = 1, panel.grid = element_blank())
  }
  if (length(models) > 1L) out <- out + facet_wrap(vars(.data$model))
  out
}

##' An area of a calibration curve with three decimals
##' @noRd
format_area <- function(area) formatC(area, format = "f", digits = 3L)

##' The areas of the models, in small text in a corner of their panels
##' @noRd
area_labels <- function(area, colours, y, vjust) {
  geom_text(data = area, aes(label = .data$label), x = 0.98, y = y, hjust = 1,
            vjust = vjust, size = 2.8, lineheight = 0.95, colour = colours$text)
}

##' The calibration curves of the models, with the band of a correct model
##' @noRd
plot_calibration_curves <- function(x, models, area, colours) {
  curves <- x$calibration$curves
  curves <- curves[curves$model %in% models, ]
  curves$model <- factor(curves$model, levels = models)
  ggplot(curves, aes(x = .data$u)) +
    geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper),
                fill = colours$band) +
    annotate("segment", x = 0, y = 0, xend = 1, yend = 1, linetype = "22",
             colour = colours$reference, linewidth = 0.4) +
    geom_line(aes(y = .data$curve), colour = colours$line, linewidth = 0.7) +
    area_labels(area, colours, y = 0.03, vjust = 0) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(y = "Cumulative distribution of the PIT")
}

##' The histograms of the PIT of the models, with the band of a correct model
##' @noRd
plot_calibration_histograms <- function(x, models, area, bins, colours) {
  bars <- do.call(rbind, lapply(as.character(area$model), function(name) {
    curve <- x$calibration$curves
    curve <- t(curve$curve[curve$model == name])
    height <- as.numeric(pit_histogram(curve, bins))
    band <- calibration_band(
      height, pit_histogram(x$calibration$simulated[[name]], bins)
    )
    data.frame(model = name, bin = seq_len(bins), height = height,
               lower = band$lower, upper = band$upper)
  }))
  bars$model <- factor(bars$model, levels = models)
  # The band as steps over the bins.
  steps <- rbind(transform(bars, u = (bars$bin - 1) / bins),
                 transform(bars, u = bars$bin / bins))
  steps <- steps[order(steps$model, steps$bin, steps$u), ]
  # Room above the bars for the areas.
  top <- 1.25 * max(bars$height, bars$upper)
  ggplot(bars) +
    geom_ribbon(data = steps, aes(x = .data$u, ymin = .data$lower,
                                  ymax = .data$upper), fill = colours$limits) +
    geom_rect(aes(xmin = (.data$bin - 1) / bins, xmax = .data$bin / bins,
                  ymin = 0, ymax = .data$height),
              fill = colours$fill, colour = "white", linewidth = 0.25) +
    geom_hline(yintercept = 1, colour = colours$reference, linewidth = 0.3) +
    lapply(c("lower", "upper"), function(limit) {
      geom_path(data = steps, aes(x = .data$u, y = .data[[limit]]),
                colour = colours$dark, linewidth = 0.35, linetype = "42")
    }) +
    area_labels(area, colours, y = 0.98 * top, vjust = 1) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, top), expand = FALSE) +
    labs(y = "Density")
}

##' Maps of a score of some models, or of its differences from a reference
##' @noRd
plot_scores <- function(x, models, score, reference) {
  values <- observation_scores(x, score)
  shown <- models
  if (!is.null(reference)) {
    values <- values - values[, reference]
    shown <- setdiff(models, reference)
    if (!length(shown)) {
      abort_input("{.arg model} must name a model other than {.arg reference}.", arg = "model")
    }
  }
  locations <- score_locations(x, values[, shown, drop = FALSE])
  if (!nrow(locations)) {
    abort_input("{.arg x} has no {score} to map.", arg = "x")
  }
  fill <- if (is.null(reference)) {
    score_fill_scale(locations$value, score)
  } else {
    difference_fill_scale(locations$value, score, reference)
  }
  out <- ggplot(locations) +
    geom_sf(aes(fill = .data$value), shape = 21, size = 2.1, stroke = 0.3,
            colour = "grey35") +
    fill + validation_map_theme()
  if (length(shown) > 1L || !is.null(reference)) {
    out <- out + facet_wrap(vars(.data$model))
  }
  out
}

##' The scores at the test locations, one row for each model and location
##'
##' @param values A matrix from `observation_scores()`, or differences of its
##'   columns.
##' @noRd
score_locations <- function(x, values) {
  rows <- as.integer(rownames(values))
  location <- x$test_sets$location[rows]
  first <- rows[!duplicated(location)]
  pieces <- lapply(colnames(values), function(name) {
    at_location <- tapply(values[, name], factor(location, unique(location)),
                          function(v) mean(v, na.rm = TRUE))
    st_sf(model = name, value = as.numeric(at_location),
          geometry = x$test_sets$geometry[first])
  })
  out <- do.call(rbind, pieces)
  out$model <- factor(out$model, levels = colnames(values))
  out <- out[is.finite(out$value), ]
  # The largest values on top.
  out[order(abs(out$value)), ]
}

##' The horizontal colour bar of the maps of scores, above them
##' @noRd
score_colourbar <- function(width) {
  guide_colourbar(
    title.position = "top", direction = "horizontal",
    theme = theme(
      legend.key.width = unit(width, "pt"), legend.key.height = unit(7, "pt"),
      legend.ticks = element_blank(),
      legend.frame = element_rect(colour = "grey45", linewidth = 0.35),
      legend.title = element_text(hjust = 0.5, size = 9.5),
      legend.text = element_text(size = 8)
    )
  )
}

##' The colours of a score, one scale for every model
##'
##' It ends at the 98th percentile of the scores, so that a few large ones
##' don't leave the others alike.
##' @noRd
score_fill_scale <- function(values, score) {
  limits <- c(min(values), as.numeric(stats::quantile(values, 0.98)))
  if (!diff(limits) > 0) limits <- range(values)
  breaks <- scales::breaks_extended(5)(limits)
  breaks <- breaks[breaks >= limits[1L] & breaks <= limits[2L]]
  beyond <- any(values > limits[2L])
  labels <- function(breaks) {
    text <- minus_labels(breaks)
    if (beyond) text[which.max(breaks)] <- paste(text[which.max(breaks)],
                                                 "or more")
    text
  }
  limits <- if (beyond && length(breaks)) c(limits[1L], max(breaks)) else limits
  scale_fill_gradientn(
    colours = prediction_role_palette("uncertainty"), limits = limits,
    breaks = breaks, labels = labels, oob = scales::oob_squish,
    name = paste0(score, " (lower is better)"), guide = score_colourbar(150)
  )
}

##' The colours of the differences of a score from a reference model
##'
##' A diverging scale with zero in the middle, whose ends are the 98th
##' percentile of the sizes of the differences and say which side is better.
##' @noRd
difference_fill_scale <- function(values, score, reference) {
  limit <- as.numeric(stats::quantile(abs(values), 0.98))
  if (!limit > 0) limit <- max(abs(values), 1e-8)
  limit <- max(pretty(c(0, limit), n = 3L))
  beyond <- any(abs(values) > limit)
  labels <- function(breaks) {
    text <- minus_labels(breaks)
    ends <- c(which.min(breaks), which.max(breaks))
    text[ends] <- paste0(text[ends], if (beyond) c(" or less", " or more"),
                         c("\nbetter than ", "\nworse than "), reference)
    text
  }
  scale_fill_gradientn(
    colours = prediction_role_palette("signed"),
    limits = c(-limit, limit), breaks = c(-1, -0.5, 0, 0.5, 1) * limit,
    labels = labels, oob = scales::oob_squish,
    name = paste0("Difference in ", score, " from ", reference),
    guide = score_colourbar(210)
  )
}

##' Numbers of a legend, with a minus sign instead of a hyphen
##' @noRd
minus_labels <- function(values) {
  sub("^-", "\u2212", format(values, trim = TRUE, drop0trailing = TRUE))
}

##' @title Plot test sets
##'
##' @description
##' Maps the test sets of a cross-validation, or draws the distances of the
##' test locations from their training data. `plot()` of the test sets does
##' the same.
##'
##' @param x Test sets from [create_test_sets()]. For `type = "distance"`, also
##'   a named list of them, to compare.
##' @param type What to draw, a `"map"` of the test sets (the default) or the
##'   distribution of the `"distance"` from a test location to its nearest
##'   training location.
##' @param grid For `type = "distance"`, optional locations at which the model
##'   will predict, such as a grid from [create_grid()] or any `sf` object of
##'   points. Their distances to the nearest location of the data are drawn
##'   too.
##' @param colour The colour of the test locations and of the first line of the
##'   distances. The other clusters and lines have colours of their own.
##' @param ... For `plot()`, the arguments of `plot_test_sets()`. For
##'   `plot_test_sets()`, not used, and an error if given, so that a misspelt
##'   argument is not ignored.
##'
##' @details
##' Clusters are mapped together, with the number of each test set at its
##' centre. Other test sets have a panel each, the first twelve of them, with
##' the test locations in `colour` and the training locations in grey. Leave-one-out
##' test sets are the locations themselves.
##'
##' The distances are drawn as cumulative distributions. When the curve of the
##' test sets is close to that of the prediction grid, the cross-validation
##' predicts about as far from the data as the map does, and its scores tell
##' how good the map is. Test locations that are further from their training
##' data make the model look worse, and test locations that are closer make
##' it look better.
##'
##' @return A \pkg{ggplot2} object, to which components can be added.
##'
##' @seealso [create_test_sets()] and [assess_prediction()].
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:100, ]
##' clusters <- create_test_sets(example_data, method = "cluster", n_sets = 5,
##'                              seed = 1)
##' plot_test_sets(clusters)
##'
##' # Compare the distances with those of a prediction grid
##' grid <- create_grid(create_convex_hull(example_data), spacing = 20)
##' plot_test_sets(list("Leave-one-out" = create_test_sets(example_data),
##'                     "Clusters" = clusters),
##'                type = "distance", grid = grid)
##' @export
plot_test_sets <- function(x, type = c("map", "distance"), grid = NULL,
                           colour = "#00857F", ...) {
  check_supplied(x)
  check_dots_unused(...)
  type <- match_choice(type)
  check_colour(colour)
  colours <- validation_palette(colour)
  if (type == "distance") {
    test_sets <- if (inherits(x, "RiskMap_test_sets")) list(x) else x
    if (!is.list(test_sets) || !length(test_sets) ||
        !all(vapply(test_sets, inherits, logical(1), "RiskMap_test_sets"))) {
      abort_input("{.arg x} must be test sets from {.fn create_test_sets}, or a list of them, not {.obj_type_friendly {x}}.", arg = "x")
    }
    return(plot_test_distances(test_sets, grid, colours))
  }
  check_plot_object(x, "RiskMap_test_sets", "plot_test_sets",
                    "create_test_sets")
  if (!is.null(grid)) {
    abort_input("{.arg grid} is for {.code type = \"distance\"}.", arg = "grid")
  }
  if (x$method == "cluster") {
    plot_cluster_map(x, colours)
  } else if (x$method == "leave_one_out") {
    plot_leave_one_out_map(x, colours)
  } else {
    plot_test_set_panels(x, colours)
  }
}

##' @rdname plot_test_sets
##' @method plot RiskMap_test_sets
##' @export
plot.RiskMap_test_sets <- function(x, ...) plot_test_sets(x, ...)

##' The distinct locations of some rows of the data of test sets
##' @noRd
distinct_test_rows <- function(x, rows) rows[!duplicated(x$location[rows])]

##' A map of the clusters, with their numbers
##' @noRd
plot_cluster_map <- function(x, colours) {
  set <- rep(seq_along(x$sets), lengths(lapply(x$sets, `[[`, "test")))
  rows <- unlist(lapply(x$sets, `[[`, "test"))
  first <- !duplicated(paste(set, x$location[rows]))
  set <- set[first]
  rows <- rows[first]
  number <- stats::ave(seq_along(x$sets), x$repetition, FUN = seq_along)
  centres <- cluster_centres(x, set, rows, colours$sets)
  points <- st_sf(colour = centres$colour[set],
                  repetition = paste("Repeat", x$repetition[set]),
                  geometry = x$geometry[rows])
  labels <- st_sf(label = number, colour = centres$colour,
                  repetition = paste("Repeat", x$repetition),
                  geometry = centres$geometry)
  hulls <- cluster_hulls(x, set, rows)
  out <- ggplot()
  if (length(hulls$set)) {
    hulls <- st_sf(colour = centres$colour[hulls$set],
                   repetition = paste("Repeat", x$repetition[hulls$set]),
                   geometry = hulls$geometry)
    out <- out +
      geom_sf(data = hulls, aes(fill = .data$colour, colour = .data$colour),
              alpha = 0.12, linewidth = 0.3)
  }
  out <- out +
    geom_sf(data = points, aes(colour = .data$colour), size = 1.5) +
    geom_sf(data = labels, colour = "white", size = 4.6, alpha = 0.85) +
    geom_sf_text(data = labels, aes(label = .data$label), size = 2.9,
                 colour = colours$dark, fontface = "bold",
                 fun.geometry = identity) +
    scale_colour_identity() + scale_fill_identity() +
    labs(x = NULL, y = NULL) +
    validation_map_theme()
  if (max(x$repetition) > 1L) {
    order <- paste("Repeat", seq_len(max(x$repetition)))
    out <- out + facet_wrap(vars(factor(.data$repetition, levels = order)))
  }
  out
}

##' The centres of clusters, with colours that differ between neighbours
##'
##' Each cluster takes, among the five colours, the one whose nearest cluster
##' is furthest, so that five colours are enough for any number of clusters.
##' @noRd
cluster_centres <- function(x, set, rows, palette) {
  n_sets <- length(x$sets)
  coordinates <- x$coordinates[rows, , drop = FALSE]
  centres <- cbind(tapply(coordinates[, 1L], set, mean),
                   tapply(coordinates[, 2L], set, mean))
  colour <- integer(n_sets)
  for (k in seq_len(n_sets)) {
    earlier <- which(x$repetition == x$repetition[k] & colour > 0L)
    nearest <- vapply(seq_along(palette), function(candidate) {
      same <- earlier[colour[earlier] == candidate]
      if (!length(same)) return(Inf)
      min((centres[same, 1L] - centres[k, 1L])^2 +
            (centres[same, 2L] - centres[k, 2L])^2)
    }, numeric(1))
    colour[k] <- which.max(nearest)
  }
  list(colour = palette[colour], geometry = plane_points(x, centres))
}

##' Points of the projected plane, in the CRS of the data
##' @noRd
plane_points <- function(x, coordinates) {
  factor <- crs_to_distance_factor(st_sfc(crs = x$crs), x$distance_units)
  points <- st_as_sf(data.frame(x = coordinates[, 1L] / factor,
                                y = coordinates[, 2L] / factor),
                     coords = c("x", "y"), crs = x$crs)
  st_transform(st_geometry(points), st_crs(x$geometry))
}

##' The convex hulls of the clusters with at least three locations
##' @return A list with the clusters that have one, `set`, and their `geometry`.
##' @noRd
cluster_hulls <- function(x, set, rows) {
  large <- which(tabulate(set, length(x$sets)) >= 3L)
  if (!length(large)) return(list(set = integer(0)))
  geometry <- do.call(c, lapply(large, function(k) {
    st_convex_hull(st_union(x$geometry[rows[set == k]]))
  }))
  # The hull of locations on one line isn't a polygon.
  polygon <- grepl("POLYGON", st_geometry_type(geometry))
  list(set = large[polygon], geometry = geometry[polygon])
}

##' A map of the locations of leave-one-out test sets
##' @noRd
plot_leave_one_out_map <- function(x, colours) {
  rows <- distinct_test_rows(x, seq_along(x$location))
  ggplot(st_sf(geometry = x$geometry[rows])) +
    geom_sf(colour = colours$line, size = 1.7) +
    labs(subtitle = paste0("Each of the ", format_count(length(rows)),
                           " locations is a test set")) +
    validation_map_theme()
}

##' Maps of the test and training locations, one for each test set
##' @noRd
plot_test_set_panels <- function(x, colours, n_shown = 12L) {
  n_sets <- length(x$sets)
  shown <- seq_len(min(n_sets, n_shown))
  pieces <- lapply(shown, function(k) {
    training <- distinct_test_rows(x, x$sets[[k]]$train)
    tested <- distinct_test_rows(x, x$sets[[k]]$test)
    st_sf(set = k,
          role = rep(c("Training", "Test"), c(length(training), length(tested))),
          geometry = x$geometry[c(training, tested)])
  })
  data <- do.call(rbind, pieces)
  data$set <- factor(data$set, levels = shown,
                     labels = if (n_sets > 1L) paste("Test set", shown) else "")
  data$role <- factor(data$role, levels = c("Training", "Test"))
  out <- ggplot(data) +
    geom_sf(aes(colour = .data$role, size = .data$role)) +
    scale_colour_manual(values = c(Training = colours$training,
                                   Test = colours$line), name = NULL) +
    scale_size_manual(values = c(Training = 0.8, Test = 1.4), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    validation_map_theme()
  if (n_sets > length(shown)) {
    out <- out + labs(subtitle = paste0("The first ", length(shown), " of ",
                                        n_sets, " test sets"))
  }
  if (n_sets > 1L) out <- out + facet_wrap(vars(.data$set))
  out
}

##' The distances of test locations from their training data
##' @noRd
plot_test_distances <- function(test_sets, grid, colours) {
  names <- names(test_sets)
  if (is.null(names)) names <- rep("", length(test_sets))
  unnamed <- names == ""
  names[unnamed] <- if (length(test_sets) == 1L) "Test locations" else
    paste("Test sets", which(unnamed))
  grid_name <- if (!is.null(grid)) "Prediction locations"
  repeated <- unique(c(grid_name, names)[duplicated(c(grid_name, names))])
  if (length(repeated)) {
    hint <- if (any(repeated %in% grid_name)) {
      "The curve of {.arg grid} is called {.val {grid_name}}."
    }
    abort_input(c("The names of {.arg x} must differ, and {.val {repeated}} {?is/are} repeated.", "i" = hint), arg = "x")
  }
  units <- unique(vapply(test_sets, `[[`, character(1), "distance_units"))
  if (length(units) > 1L) {
    abort_input("The test sets in {.arg x} must use the same units of distance, not {.val {units}}.", arg = "x")
  }
  if (length(test_sets) > length(colours$sets)) {
    abort_input("{.arg x} must have at most {length(colours$sets)} sets of test sets to compare.", arg = "x")
  }
  # A cumulative distribution, from zero at the smallest distance.
  cumulative <- function(distance, what) {
    data.frame(distance = sort(distance)[c(1L, seq_along(distance))],
               below = c(0, seq_along(distance)) / length(distance),
               what = what)
  }
  pieces <- Map(function(sets, name) {
    cumulative(test_set_distances(sets)$distance, name)
  }, test_sets, names)
  if (!is.null(grid)) {
    pieces <- c(list(cumulative(grid_distances(test_sets[[1L]], grid),
                                grid_name)), pieces)
  }
  data <- do.call(rbind, unname(pieces))
  levels <- unique(data$what)
  data$what <- factor(data$what, levels = levels)
  line_colours <- c(if (!is.null(grid)) colours$dark,
                    colours$sets[seq_along(test_sets)])
  ggplot(data, aes(x = .data$distance, y = .data$below,
                   colour = .data$what, linetype = .data$what)) +
    geom_step(linewidth = 0.7) +
    scale_colour_manual(values = stats::setNames(line_colours, levels),
                        name = NULL) +
    scale_linetype_manual(
      values = stats::setNames(c(if (!is.null(grid)) "22",
                                 rep("solid", length(test_sets))), levels),
      name = NULL
    ) +
    scale_y_continuous(labels = scales::label_percent(),
                       limits = c(0, 1), expand = expansion(c(0, 0.02))) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(c(0, 0.03))) +
    labs(x = paste0("Distance to the nearest location with data (", units, ")"),
         y = "Locations within the distance") +
    validation_theme()
}

##' The distance of prediction locations to the nearest location of the data
##' @noRd
grid_distances <- function(test_sets, grid, call = rlang::caller_env()) {
  if (!inherits(grid, c("sf", "sfc"))) {
    abort_input("{.arg grid} must be an {.cls sf} object of points, not {.obj_type_friendly {grid}}.", arg = "grid", call = call)
  }
  check_data(grid, type = "sfc", arg = "grid")
  projected <- st_transform(st_geometry(grid), test_sets$crs)
  nearest_distances(coordinates_in_units(projected, test_sets$distance_units),
                    test_sets$coordinates)
}
