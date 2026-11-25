##' Plot a RiskMap prediction
##'
##' Draws a publication-ready map from the location-level or polygon-level
##' output returned by [predict.RiskMap()].
##'
##' @param x A `RiskMap_prediction` object.
##' @param quantity Predictive quantity to display, such as `"response"`,
##'   `"spatial"`, `"link"`, `"fixed"` (the fixed-effect contribution), or the
##'   name of a target supplied through `transform` in [predict.RiskMap()].
##'   Several quantities can be supplied, for example
##'   `c("fixed", "spatial")`, and are arranged like several summaries. Supply
##'   several quantities or several summaries, not both.
##' @param summary Summary to display, such as `"mean"`, `"median"`, `"lower"`,
##'   `"upper"`, `"sd"`, `"cv"`, `"interval_width"`, or `"exceedance"`. When
##'   `NULL` (the default), the first summary requested in
##'   [predict.RiskMap()] is shown. Several summaries can be supplied; they are
##'   then arranged in two columns, or one per row when
##'   `spatial_scale = "both"`.
##'   Ignored for `"fixed"` and `"offset"`, which have a single value.
##' @param spatial_scale Spatial scale to display. `"auto"` (the default) shows
##'   the polygon summaries when [predict.RiskMap()] was given `boundaries`,
##'   and the location-level predictions otherwise. `"location"` and `"area"`
##'   select one scale, and `"both"` draws them side by side with a shared
##'   colour scale.
##' @param aggregate Within-polygon aggregation to display at the area scale,
##'   such as `"mean"` or `"sum"`. When `NULL`, the first aggregation
##'   calculated is used.
##' @param threshold For an exceedance map, the response-scale threshold to
##'   display. It may be omitted when only one threshold was calculated.
##' @param palette Colour palette. `NULL` selects a default that depends on
##'   what the map shows: a calm-to-alarm sequential palette for predicted
##'   levels, a single-hue palette for uncertainty (`"sd"`, `"cv"`,
##'   `"interval_width"`), a blue-to-red palette centred at 0.5 for exceedance
##'   probabilities, and a palette centred at zero for signed effects such as
##'   `"spatial"` and `"link"`. A character vector of colours or the name of an
##'   HCL palette accepted by [grDevices::hcl.colors()] can be supplied. Prefix
##'   a named palette with `"-"` to reverse it, for example `"-Viridis"`.
##' @param limits Optional two-element numeric vector fixing the colour limits.
##'   By default, limits adapt to the values displayed. Supply the same limits
##'   when separate maps must be compared directly. Exceedance probabilities
##'   always use limits from 0 to 1.
##' @param breaks Optional increasing numeric vector of class boundaries. When
##'   supplied, values are shown in classes, for example programme thresholds
##'   such as `c(0, 0.05, 0.1, 0.2, 1)`.
##' @param north_arrow Logical. Add a north arrow. Defaults to `FALSE`.
##' @param scale_bar Logical. Add a scale bar. Defaults to `FALSE`.
##' @param ... Additional arguments passed to the ggplot2 fill scale, for
##'   example `transform = "sqrt"` for skewed values.
##'
##' @details
##' The legend title names the mapped quantity and is the only title. For
##' binomial models the response is labelled as prevalence and for Poisson
##' models as a rate. Probabilities are labelled as percentages. The title can
##' be replaced with `labs(fill = ...)`.
##'
##' Grid cells created by [create_grid()] are drawn as squares clipped to the
##' grid's study boundary, and the boundary itself is outlined. Maps are ordinary `ggplot` objects drawn with
##' [ggplot2::geom_sf()], so further `sf` layers, such as survey locations,
##' can be added in any coordinate reference system.
##'
##' When `spatial_scale = "both"`, the two scales are drawn side by side,
##' locations first, and share one colour scale so that colours mean the same
##' in both panels. This requires an aggregation on
##' the same scale as the locations (`"mean"`, `"median"`, `"min"` or
##' `"max"`).
##'
##' Outlines of the `boundaries` supplied to [predict.RiskMap()] are drawn
##' automatically on location-level maps.
##'
##' @return A `ggplot` object. When several quantities or summaries are
##'   requested, a `patchwork` object combining one `ggplot` per map.
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:50, ]
##' fit <- glgpm(y ~ gp(), data = example_data, family = "gaussian",
##'              messages = FALSE)
##' areas <- sf::st_sf(geometry = sf::st_make_grid(example_data, n = c(2, 1)))
##' prediction <- predict(fit, newdata = example_data, boundaries = areas,
##'                       aggregate = "mean", messages = FALSE)
##' plot(prediction)
##' plot(prediction, spatial_scale = "both")
##' @method plot RiskMap_prediction
##' @export
##' @importFrom ggspatial annotation_north_arrow annotation_scale
##' @importFrom grid unit
plot.RiskMap_prediction <- function(
    x,
    quantity = "response",
    summary = NULL,
    spatial_scale = c("auto", "location", "area", "both"),
    aggregate = NULL,
    threshold = NULL,
    palette = NULL,
    limits = NULL,
    breaks = NULL,
    north_arrow = FALSE,
    scale_bar = FALSE,
    ...) {
  if (!inherits(x, "RiskMap_prediction")) {
    stop("'x' must be a RiskMap_prediction object", call. = FALSE)
  }
  spatial_scale <- match.arg(spatial_scale)
  check_plot_flag(north_arrow, "north_arrow")
  check_plot_flag(scale_bar, "scale_bar")
  if (!is.character(quantity) || !length(quantity) || anyNA(quantity) ||
      anyDuplicated(quantity)) {
    stop("'quantity' must contain unique quantity names", call. = FALSE)
  }
  if (length(quantity) > 1L && length(summary) > 1L) {
    stop("Supply several quantities or several summaries, not both",
         call. = FALSE)
  }
  if (length(quantity) == 1L) {
    if (quantity %in% c("fixed", "offset")) {
      summary <- "estimate"
    } else if (is.null(summary)) {
      summary <- default_prediction_summary(x)
    }
    if (!is.character(summary) || !length(summary) || anyNA(summary) ||
        anyDuplicated(summary)) {
      stop("'summary' must contain unique summary names", call. = FALSE)
    }
  }

  if (length(quantity) > 1L || length(summary) > 1L) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Plotting several maps requires the 'patchwork' package",
           call. = FALSE)
    }
    if (length(quantity) > 1L) {
      maps <- lapply(quantity, function(one_quantity) {
        plot.RiskMap_prediction(
          x, quantity = one_quantity, summary = summary,
          spatial_scale = spatial_scale, aggregate = aggregate,
          threshold = threshold, palette = palette, limits = limits,
          breaks = breaks, north_arrow = north_arrow, scale_bar = scale_bar,
          ...
        )
      })
    } else {
      maps <- lapply(summary, function(one_summary) {
        plot.RiskMap_prediction(
          x, quantity = quantity, summary = one_summary,
          spatial_scale = spatial_scale, aggregate = aggregate,
          threshold = if (identical(one_summary, "exceedance")) threshold,
          palette = palette, limits = limits, breaks = breaks,
          north_arrow = north_arrow, scale_bar = scale_bar, ...
        )
      })
    }
    # Maps carry their legend on top, so two columns keep them readable.
    layout_columns <- if (spatial_scale == "both") 1L else
      min(2L, length(maps))
    # Side margins keep neighbouring legends apart at book and report widths.
    return(patchwork::wrap_plots(maps, ncol = layout_columns) &
             theme(plot.margin = margin(4, 14, 4, 14)))
  }

  has_areas <- !is.null(x$boundaries) && !is.null(x$boundary_summary)
  if (spatial_scale == "auto") {
    spatial_scale <- if (has_areas || !is.null(aggregate)) "area" else
      "location"
  }
  if (spatial_scale %in% c("area", "both") && !has_areas) {
    stop("This prediction object does not contain boundary summaries",
         call. = FALSE)
  }
  if (spatial_scale %in% c("area", "both")) {
    aggregate <- select_prediction_plot_value(
      aggregate %||% unique(x$boundary_summary$aggregate)[[1L]],
      unique(x$boundary_summary$aggregate), "aggregate"
    )
  }
  if (spatial_scale == "both" &&
      !aggregate %in% c("mean", "median", "min", "max")) {
    stop("spatial_scale = \"both\" needs an aggregation on the same scale as ",
         "the locations: \"mean\", \"median\", \"min\" or \"max\"",
         call. = FALSE)
  }

  panels <- list()
  if (spatial_scale %in% c("location", "both")) {
    panels$location <- prediction_plot_data(x, quantity, summary, NULL,
                                            threshold)
  }
  if (spatial_scale %in% c("area", "both")) {
    panels$area <- prediction_plot_data(x, quantity, summary, aggregate,
                                        threshold)
  }
  threshold <- panels[[1L]]$threshold
  role <- prediction_plot_role(quantity, summary)
  probability <- prediction_is_probability(x, quantity, summary) &&
    (spatial_scale == "location" ||
       aggregate %in% c("mean", "median", "min", "max") ||
       identical(summary, "exceedance"))

  values <- unlist(lapply(panels, function(panel) {
    panel$data$.prediction_value
  }), use.names = FALSE)
  if (role == "evidence" && !is.null(limits) &&
      !identical(as.numeric(limits), c(0, 1))) {
    warning("Exceedance probabilities always use limits from 0 to 1",
            call. = FALSE)
  }
  limits <- check_plot_limits(limits)
  breaks <- check_plot_breaks(breaks)
  # Default limits end on round values so the legend labels both ends.
  colour_limits <- switch(
    role,
    evidence = c(0, 1),
    signed = limits %||%
      (c(-1, 1) * max(abs(pretty(range(values, na.rm = TRUE), n = 4)))),
    limits %||% range(pretty(range(values, na.rm = TRUE), n = 4))
  )
  colours <- if (is.null(palette)) prediction_role_palette(role) else
    prediction_plot_palette(palette)

  facet_labels <- c(location = "Locations",
                    area = paste("Area", aggregate %||% ""))
  map <- ggplot()
  facet_levels <- unname(facet_labels[names(panels)])
  for (scale_name in names(panels)) {
    map <- add_prediction_layer(map, x, panels[[scale_name]],
                                facet_labels[[scale_name]], facet_levels)
  }
  if (length(panels) > 1L) {
    map <- map + facet_wrap(~ .spatial_scale, nrow = 1L)
  }

  map <- map + prediction_fill_scale(
    colours = colours, limits = colour_limits, values = values,
    breaks = breaks, role = role, probability = probability,
    ...
  ) +
    # Set through labs() so that users can replace it with labs(fill = ...).
    labs(fill = prediction_legend_title(
      x, quantity, summary, threshold,
      if (spatial_scale == "area") aggregate
    )) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.justification = "center",
      strip.text = element_blank()
    )
  if (north_arrow) map <- map + annotation_north_arrow(location = "tr")
  if (scale_bar) map <- map + annotation_scale(location = "bl")
  map
}

##' Print a RiskMap prediction
##'
##' Lists the predicted quantities, summaries and spatial scales stored in the
##' object, and how to map them.
##'
##' @param x A `RiskMap_prediction` object.
##' @param ... Not used.
##' @return `x`, invisibly.
##' @method print RiskMap_prediction
##' @export
print.RiskMap_prediction <- function(x, ...) {
  locations <- nrow(x$predictions)
  cat("RiskMap prediction: ", locations, " location",
      if (locations == 1L) "" else "s", sep = "")
  if (!is.null(x$boundaries)) {
    cat(" and", nrow(x$boundaries), "areas")
  }
  cat("\n")
  if (!is.null(x$family)) {
    cat("Response: ", prediction_quantity_label(x, "response"), " (",
        x$family, " model)\n", sep = "")
  }
  if (!is.null(x$type)) cat("Draws:", x$type, "\n")
  quantities <- unique(c("response", x$transform, "link", "spatial"))
  cat("Quantities:", paste(quantities, collapse = ", "), "\n")
  if (length(x$summaries)) {
    cat("Summaries:", paste(x$summaries, collapse = ", "), "\n")
  }
  if (length(x$thresholds)) {
    cat("Exceedance thresholds:", paste(x$thresholds, collapse = ", "), "\n")
  }
  if (length(x$aggregate)) {
    cat("Area aggregations:", paste(x$aggregate, collapse = ", "), "\n")
  }
  cat("\nMap with plot(x, summary = \"", default_prediction_summary(x), "\")",
      sep = "")
  if (!is.null(x$boundaries)) cat(", adding spatial_scale = \"both\"")
  cat(".\n")
  invisible(x)
}

#' First summary requested in predict(), used as the default map
#' @noRd
default_prediction_summary <- function(x) {
  if (length(x$summaries)) return(x$summaries[[1L]])
  "mean"
}

#' Add the layer for one spatial scale to a prediction map
#' @noRd
add_prediction_layer <- function(map, x, panel, facet_label, facet_levels) {
  facet_label <- factor(facet_label, levels = facet_levels)
  data <- panel$data
  if (isTRUE(panel$polygon)) {
    data$.spatial_scale <- facet_label
    return(map + geom_sf(
      data = data, aes(fill = .data$.prediction_value),
      colour = "black", linewidth = 0.16, inherit.aes = FALSE
    ))
  }
  metadata <- attr(data, "grid_metadata", exact = TRUE)
  if (!is.null(metadata) && identical(metadata$cell_shape, "square")) {
    map <- map + geom_sf(
      data = prediction_cells(data, metadata, facet_label),
      aes(fill = .data$.prediction_value), colour = NA, inherit.aes = FALSE
    )
    if (!is.null(metadata$boundary)) {
      # The study boundary marks where the clipped cells end.
      study_area <- st_sf(.spatial_scale = facet_label,
                          geometry = st_geometry(metadata$boundary))
      map <- map + geom_sf(data = study_area, fill = NA, colour = "grey20",
                           linewidth = 0.3, inherit.aes = FALSE)
    }
  } else {
    data$.spatial_scale <- facet_label
    map <- map + geom_sf(
      data = data, aes(fill = .data$.prediction_value),
      # A transparent rather than NA colour, which drops every point.
      shape = 21, colour = "transparent", size = 1.8, inherit.aes = FALSE
    )
  }
  if (!is.null(x$boundaries)) {
    outlines <- st_sf(.spatial_scale = facet_label,
                      geometry = st_geometry(x$boundaries))
    map <- map + geom_sf(data = outlines, fill = NA, colour = "black",
                         linewidth = 0.16, inherit.aes = FALSE)
  }
  map
}

#' Square grid cells as polygons, clipped to the grid's study boundary
#' @noRd
prediction_cells <- function(data, metadata, facet_label) {
  grid_crs <- st_crs(metadata$grid_crs)
  centres <- st_transform(st_geometry(data), grid_crs)
  cells <- st_buffer(centres, dist = metadata$cellsize / 2,
                     endCapStyle = "SQUARE")
  cells <- st_sf(.prediction_value = data$.prediction_value,
                 .spatial_scale = facet_label, geometry = cells)
  if (!is.null(metadata$boundary)) {
    boundary <- st_transform(metadata$boundary, grid_crs)
    cells <- suppressWarnings(st_intersection(cells, boundary))
    # Cells that only touch the boundary leave lines or points; keep areas.
    cells <- suppressWarnings(st_collection_extract(cells, "POLYGON"))
    cells <- cells[!st_is_empty(cells), ]
  }
  cells
}

#' Fill scale for a prediction map
#' @noRd
prediction_fill_scale <- function(colours, limits, values, breaks, role,
                                  probability, ...) {
  range_values <- limits %||% range(values, na.rm = TRUE)
  labels <- if (probability) {
    accuracy <- if (diff(range_values) >= 0.05) 1 else 0.1
    scales::label_percent(accuracy = accuracy)
  } else {
    function(value) {
      sub("^-", "\u2212", format(value, trim = TRUE, drop0trailing = TRUE))
    }
  }
  guide_theme <- theme(
    legend.key.width = unit(105, "pt"),
    legend.key.height = unit(7, "pt"),
    legend.ticks = element_blank(),
    legend.frame = element_rect(colour = "grey45", linewidth = 0.35),
    legend.title = element_text(hjust = 0.5, size = 9),
    legend.text = element_text(size = 7.5)
  )

  if (!is.null(breaks)) {
    n_classes <- length(breaks) - 1L
    class_colours <- colours[round(seq(1, length(colours),
                                       length.out = n_classes))]
    midpoints <- (utils::head(breaks, -1L) + utils::tail(breaks, -1L)) / 2
    positions <- (midpoints - breaks[[1L]]) / diff(range(breaks))
    return(scale_fill_stepsn(
      colours = c(class_colours[[1L]], class_colours,
                  class_colours[[n_classes]]),
      values = c(0, positions, 1),
      limits = range(breaks), breaks = breaks, labels = labels,
      oob = scales::oob_squish, na.value = "grey90",
      guide = guide_coloursteps(title.position = "top",
                                direction = "horizontal",
                                theme = guide_theme),
      ...
    ))
  }

  scale_breaks <- if (role == "evidence") seq(0, 1, by = 0.25) else {
    candidates <- pretty(range_values, n = 4)
    candidates[candidates >= range_values[[1L]] &
                 candidates <= range_values[[2L]]]
  }
  scale_fill_gradientn(
    colours = colours, limits = limits, breaks = scale_breaks,
    labels = labels, oob = scales::oob_squish, na.value = "grey90",
    guide = guide_colourbar(title.position = "top", direction = "horizontal",
                            theme = guide_theme),
    ...
  )
}

#' Select and align the values requested from a prediction object
#' @noRd
prediction_plot_data <- function(x, quantity, summary, aggregate,
                                 threshold) {
  if (!is.character(quantity) || length(quantity) != 1L || is.na(quantity) ||
      !nzchar(quantity)) {
    stop("'quantity' must be one non-empty character string", call. = FALSE)
  }
  if (!is.character(summary) || length(summary) != 1L || is.na(summary) ||
      !nzchar(summary)) {
    stop("'summary' must be one non-empty character string", call. = FALSE)
  }

  if (!is.null(aggregate)) {
    candidates <- x$boundary_summary[
      x$boundary_summary$aggregate == aggregate &
        x$boundary_summary$quantity == quantity &
        x$boundary_summary$summary == summary, , drop = FALSE
    ]
    if (!nrow(candidates)) {
      stop("The requested quantity, summary and aggregation are not available",
           call. = FALSE)
    }
    if (identical(summary, "exceedance")) {
      threshold <- select_prediction_plot_value(
        threshold, unique(candidates$threshold), "threshold"
      )
      candidates <- candidates[candidates$threshold == threshold, ,
                               drop = FALSE]
    } else if (!is.null(threshold)) {
      stop("'threshold' is used only when summary = \"exceedance\"",
           call. = FALSE)
    }
    values <- rep(NA_real_, nrow(x$boundaries))
    values[candidates$boundary_id] <- candidates$value
    data <- x$boundaries
    data$.prediction_value <- values
    return(list(data = data, polygon = TRUE, threshold = threshold))
  }

  column <- if (quantity %in% c("fixed", "offset")) {
    quantity
  } else if (identical(summary, "exceedance")) {
    available <- grep(
      paste0("^", quantity, "_exceedance_"), names(x$predictions), value = TRUE
    )
    if (!length(available)) {
      stop("The requested exceedance quantity is not available", call. = FALSE)
    }
    available_thresholds <- as.numeric(sub(
      paste0("^", quantity, "_exceedance_"), "", available
    ))
    threshold <- select_prediction_plot_value(
      threshold, available_thresholds, "threshold"
    )
    available[[which.min(abs(available_thresholds - threshold))]]
  } else {
    if (!is.null(threshold)) {
      stop("'threshold' is used only when summary = \"exceedance\"",
           call. = FALSE)
    }
    paste(quantity, summary, sep = "_")
  }
  if (!column %in% names(x$predictions)) {
    stop("Prediction column ", shQuote(column), " is not available",
         call. = FALSE)
  }
  data <- x$predictions
  data$.prediction_value <- data[[column]]
  list(data = data, polygon = FALSE, threshold = threshold)
}

#' Choose an unambiguous aggregation or threshold
#' @noRd
select_prediction_plot_value <- function(value, available, argument) {
  available <- available[!is.na(available)]
  if (is.null(value)) {
    if (length(available) != 1L) {
      stop("'", argument, "' must be supplied; available values are: ",
           paste(available, collapse = ", "), call. = FALSE)
    }
    return(available[[1L]])
  }
  if (length(value) != 1L || is.na(value) || !value %in% available) {
    stop("'", argument, "' must be one of: ",
         paste(available, collapse = ", "), call. = FALSE)
  }
  value
}

#' What a map shows, which determines its default colours
#' @noRd
prediction_plot_role <- function(quantity, summary) {
  if (identical(summary, "exceedance")) return("evidence")
  if (summary %in% c("sd", "cv", "interval_width")) return("uncertainty")
  if (quantity %in% c("spatial", "link", "nugget", "random_effects",
                      "fixed", "offset")) {
    return("signed")
  }
  "level"
}

#' Default colours for each kind of map
#' @noRd
prediction_role_palette <- function(role) {
  neutral <- "#EFEBE4"
  switch(
    role,
    level = .riskmap_mist_ember,
    evidence = grDevices::colorRampPalette(
      c("#2166AC", "#8FB3D6", neutral, "#DE9585", "#B2182B")
    )(256L),
    uncertainty = rev(grDevices::hcl.colors(256L, "Purples 3"))[24:256],
    signed = grDevices::colorRampPalette(
      c("#00857F", "#7DBFB9", neutral, "#D69BC2", "#B5408F")
    )(256L)
  )
}

#' Resolve user-supplied prediction-map colours
#' @noRd
prediction_plot_palette <- function(palette) {
  if (is.null(palette)) return(.riskmap_mist_ember)
  if (!is.character(palette) || !length(palette) || anyNA(palette)) {
    stop("'palette' must be a named HCL palette or a vector of colours",
         call. = FALSE)
  }
  if (length(palette) > 1L) {
    tryCatch(grDevices::col2rgb(palette), error = function(error) {
      stop("'palette' contains an invalid colour", call. = FALSE)
    })
    return(palette)
  }
  reverse <- startsWith(palette, "-")
  name <- if (reverse) substring(palette, 2L) else palette
  colours <- tryCatch(
    grDevices::hcl.colors(256L, palette = name),
    error = function(error) stop("Unknown HCL palette ", shQuote(name),
                                 call. = FALSE)
  )
  if (reverse) rev(colours) else colours
}

#' Scientific name of a predicted quantity
#' @noRd
prediction_quantity_label <- function(x, quantity) {
  if (identical(quantity, "response")) {
    family <- x$family %||% ""
    if (identical(family, "binomial")) return("prevalence")
    if (identical(family, "poisson")) return("rate")
    response_name <- x$response_name
    if (length(response_name) == 1L && !is.na(response_name)) {
      return(response_name)
    }
    return("response")
  }
  switch(quantity,
         spatial = "spatial effect",
         fixed = "fixed effects",
         offset = "offset",
         link = "linear predictor",
         nugget = "nugget effect",
         random_effects = "random effects",
         gsub("_", " ", quantity, fixed = TRUE))
}

#' Whether a summary is a probability, labelled as a percentage
#' @noRd
prediction_is_probability <- function(x, quantity, summary) {
  identical(summary, "exceedance") ||
    (identical(quantity, "response") && identical(x$family, "binomial") &&
       !identical(summary, "cv"))
}

#' Legend title naming the mapped quantity; the map has no other title
#' @noRd
prediction_legend_title <- function(x, quantity, summary, threshold,
                                    aggregate = NULL) {
  label <- prediction_quantity_label(x, quantity)
  level <- paste0(format(100 * (x$level %||% 0.95)), "%")
  title <- switch(
    summary,
    estimate = label,
    mean = paste("Mean", label),
    median = paste("Median", label),
    lower = paste("Lower", level, "limit of", label),
    upper = paste("Upper", level, "limit of", label),
    sd = paste("SD of", label),
    cv = paste("CV of", label),
    interval_width = paste("Width of", level, "interval for", label),
    exceedance = paste0(
      "Probability ", label, " > ",
      if (prediction_is_probability(x, quantity, "mean")) {
        paste0(format(100 * threshold, drop0trailing = TRUE), "%")
      } else {
        format(threshold)
      }
    ),
    paste(summary, label)
  )
  if (!is.null(aggregate) && !identical(aggregate, "mean")) {
    title <- paste0(title, " (area ", aggregate, ")")
  }
  paste0(toupper(substring(title, 1L, 1L)), substring(title, 2L))
}

#' Validate a plot switch
#' @noRd
check_plot_flag <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", argument, "' must be TRUE or FALSE", call. = FALSE)
  }
}

#' Validate colour limits
#' @noRd
check_plot_limits <- function(limits) {
  if (is.null(limits)) return(NULL)
  if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits) ||
      any(!is.finite(limits)) || limits[[1L]] >= limits[[2L]]) {
    stop("'limits' must contain two increasing finite numbers", call. = FALSE)
  }
  limits
}

#' Validate class boundaries for a binned map
#' @noRd
check_plot_breaks <- function(breaks) {
  if (is.null(breaks)) return(NULL)
  if (!is.numeric(breaks) || length(breaks) < 3L || anyNA(breaks) ||
      any(!is.finite(breaks)) || is.unsorted(breaks, strictly = TRUE)) {
    stop("'breaks' must contain at least three increasing finite numbers",
         call. = FALSE)
  }
  breaks
}
