##' Plot a RiskMap prediction
##'
##' Draws a publication-ready map from the location-level or polygon-level
##' output returned by [predict.RiskMap()].
##'
##' @param x A `RiskMap_prediction` object.
##' @param quantity Predictive quantity to display, such as `"response"`,
##'   `"spatial"`, or `"link"`.
##' @param summary Summary to display, such as `"mean"`, `"lower"`, `"upper"`,
##'   `"interval_width"`, or `"exceedance"`.
##' @param aggregate For polygon predictions, the within-polygon aggregation to
##'   display. It may be omitted when the object contains only one aggregation.
##' @param threshold For an exceedance map, the response-scale threshold to
##'   display. It may be omitted when only one threshold was calculated.
##' @param palette Colour palette. `NULL` selects a suitable default: the
##'   sequential Mist-ember palette for ordinary prediction summaries and a
##'   blue-to-red palette centred at 0.5 for exceedance probabilities. A
##'   character vector of colours or the name of an HCL palette accepted by
##'   [grDevices::hcl.colors()] can be supplied. Prefix a named palette with
##'   `"-"` to reverse it, for example `"-Viridis"`.
##' @param limits Optional two-element numeric vector fixing the colour limits.
##'   By default, limits adapt to the values in the displayed map. Supply the
##'   same limits when maps must be compared directly. Exceedance probabilities
##'   always use limits from 0 to 1.
##' @param north_arrow Logical. Add a north arrow. Defaults to `FALSE`.
##' @param scale_bar Logical. Add a scale bar. Defaults to `FALSE`.
##' @param ... Additional arguments passed to [ggplot2::scale_fill_gradientn()].
##'
##' @details
##' The default limits emphasize spatial variation within the displayed map;
##' colours therefore need not represent the same values in two different maps.
##' Use `limits` to make maps directly comparable. Exceedance maps are an
##' exception: their scale is always fixed from 0 to 1, with 0.5 at the neutral
##' midpoint.
##'
##' @return A `ggplot` object.
##' @method plot RiskMap_prediction
##' @export
##' @importFrom ggspatial annotation_north_arrow annotation_scale
##' @importFrom grid unit
plot.RiskMap_prediction <- function(
    x,
    quantity = "response",
    summary = "mean",
    aggregate = NULL,
    threshold = NULL,
    palette = NULL,
    limits = NULL,
    north_arrow = FALSE,
    scale_bar = FALSE,
    ...) {
  if (!inherits(x, "RiskMap_prediction")) {
    stop("'x' must be a RiskMap_prediction object", call. = FALSE)
  }
  check_plot_flag(north_arrow, "north_arrow")
  check_plot_flag(scale_bar, "scale_bar")

  plot_data <- prediction_plot_data(
    x, quantity, summary, aggregate, threshold
  )
  is_exceedance <- identical(summary, "exceedance")
  if (is_exceedance && !is.null(limits) && !identical(as.numeric(limits), c(0, 1))) {
    warning("Exceedance probabilities always use limits from 0 to 1",
            call. = FALSE)
  }
  colour_limits <- if (is_exceedance) c(0, 1) else check_plot_limits(limits)
  colours <- prediction_plot_palette(palette, is_exceedance)
  legend_title <- prediction_legend_title(quantity, summary, threshold)

  map <- ggplot()
  if (isTRUE(plot_data$polygon)) {
    map <- map + geom_sf(
      data = plot_data$data,
      aes(fill = .data$.prediction_value),
      colour = "black", linewidth = 0.16
    )
  } else {
    metadata <- attr(plot_data$data, "grid_metadata", exact = TRUE)
    if (!is.null(metadata) && identical(metadata$cell_shape, "square")) {
      coords <- st_coordinates(plot_data$data)
      tile_data <- data.frame(
        x = coords[, 1L], y = coords[, 2L],
        .prediction_value = plot_data$data$.prediction_value
      )
      map <- map + geom_tile(
        data = tile_data,
        aes(x = .data$x, y = .data$y, fill = .data$.prediction_value),
        width = metadata$cellsize, height = metadata$cellsize
      ) + coord_sf(crs = st_crs(plot_data$data), datum = NA)
    } else {
      map <- map + geom_sf(
        data = plot_data$data,
        aes(colour = .data$.prediction_value), size = 1.8
      )
    }
  }

  scale_arguments <- c(list(
    colours = colours,
    limits = colour_limits,
    name = legend_title,
    na.value = "grey90",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(72, "pt"),
      barheight = unit(7, "pt"),
      ticks.colour = "black",
      frame.colour = "grey45",
      frame.linewidth = 0.35
    )
  ), list(...))
  scale_function <- if (isTRUE(plot_data$polygon) ||
                        !is.null(attr(plot_data$data, "grid_metadata", exact = TRUE))) {
    scale_fill_gradientn
  } else {
    scale_colour_gradientn
  }
  map <- map + do.call(scale_function, scale_arguments) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.justification = "center",
      legend.title = element_text(hjust = 0.5),
      plot.title = element_blank()
    )
  if (north_arrow) map <- map + annotation_north_arrow(location = "tr")
  if (scale_bar) map <- map + annotation_scale(location = "bl")
  map
}

#' Select and align the values requested from a prediction object
#' @noRd
prediction_plot_data <- function(x, quantity, summary, aggregate, threshold) {
  if (!is.character(quantity) || length(quantity) != 1L || is.na(quantity) ||
      !nzchar(quantity)) {
    stop("'quantity' must be one non-empty character string", call. = FALSE)
  }
  if (!is.character(summary) || length(summary) != 1L || is.na(summary) ||
      !nzchar(summary)) {
    stop("'summary' must be one non-empty character string", call. = FALSE)
  }

  if (isTRUE(x$aggregated)) {
    aggregate <- select_prediction_plot_value(
      aggregate, unique(x$summary$aggregate), "aggregate"
    )
    candidates <- x$summary[
      x$summary$aggregate == aggregate &
        x$summary$quantity == quantity &
        x$summary$summary == summary, , drop = FALSE
    ]
    if (!nrow(candidates)) {
      stop("The requested quantity, summary and aggregation are not available",
           call. = FALSE)
    }
    if (identical(summary, "exceedance")) {
      threshold <- select_prediction_plot_value(
        threshold, unique(candidates$threshold), "threshold"
      )
      candidates <- candidates[candidates$threshold == threshold, , drop = FALSE]
    } else if (!is.null(threshold)) {
      stop("'threshold' is used only when summary = \"exceedance\"",
           call. = FALSE)
    }
    values <- rep(NA_real_, nrow(x$data))
    values[candidates$boundary_id] <- candidates$value
    data <- x$data
    data$.prediction_value <- values
    return(list(data = data, polygon = TRUE))
  }

  if (!is.null(aggregate)) {
    stop("'aggregate' is used only for polygon predictions", call. = FALSE)
  }
  column <- if (identical(summary, "exceedance")) {
    available <- grep(
      paste0("^", quantity, "_exceedance_"), names(x$data), value = TRUE
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
  if (!column %in% names(x$data)) {
    stop("Prediction column ", shQuote(column), " is not available",
         call. = FALSE)
  }
  data <- x$data
  data$.prediction_value <- data[[column]]
  list(data = data, polygon = FALSE)
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

#' Resolve prediction-map colours
#' @noRd
prediction_plot_palette <- function(palette, exceedance = FALSE) {
  if (is.null(palette)) {
    if (exceedance) return(c("#2166AC", "#F7F7F7", "#B2182B"))
    # Designed as a perceptually ordered calm-to-alarm prevalence scale.
    return(.riskmap_mist_ember)
  }
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

#' Construct a concise legend title without duplicating a plot title
#' @noRd
prediction_legend_title <- function(quantity, summary, threshold) {
  quantity <- gsub("_", " ", quantity, fixed = TRUE)
  if (identical(summary, "exceedance")) {
    return(paste0("Probability ", quantity, " > ", format(threshold)))
  }
  paste(tools::toTitleCase(gsub("_", " ", summary, fixed = TRUE)), quantity)
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
