#' Resolve the CRS and compute pairwise distances for distance-based computations
#'
#' Projected data retain their CRS; longitude/latitude data are automatically
#' reprojected to an appropriate UTM CRS, with an informative message. This is
#' the same convention used by `glgpm()`'s `model_crs`, applied here to
#' functions with no fitted model of their own (`summarise_distance()`,
#' `variogram()`).
#'
#' @param data An `sf` object.
#' @param crs `NULL` to retain an existing projected CRS or auto-select UTM
#'   for longitude/latitude data, or a CRS (must be projected) to reproject
#'   `data` to.
#' @param crs_arg The name of the caller's CRS argument, used in messages and
#'   errors.
#' @param purpose A short description of what the CRS is used for, e.g.
#'   `"computing distances"`, used in the auto-reprojection message.
#' @param distance_units The requested coordinate units, either `"m"` or `"km"`.
#' @param dedupe Whether to drop duplicate coordinates before computing
#'   distances. Defaults to `FALSE`.
#' @return A numeric vector of pairwise distances, in `distance_units`, as
#'   returned by `dist()`.
#' @param call The call named in input errors, so that they name the exported
#'   function that called this one.
#' @noRd
extract_distances <- function(data, crs, crs_arg, purpose, distance_units,
                              dedupe = FALSE, call = rlang::caller_env()) {
  if (!is.null(crs)) {
    check_crs(crs, arg = crs_arg, call = call)
    data <- st_transform(data, crs = crs)
    if (st_is_longlat(data)) {
      abort_input("{.arg {crs_arg}} must be a projected coordinate reference system, not longitude and latitude.", arg = crs_arg, call = call)
    }
  } else if (st_is_longlat(data)) {
    auto_crs <- propose_utm(data)
    data <- st_transform(data, crs = auto_crs)
    riskmap_inform(c("{.arg data} are in longitude and latitude, so they are projected to EPSG:{auto_crs} for {purpose}.", "i" = "Set {.arg {crs_arg}} to choose the projection."))
  }
  coords <- coordinates_in_units(data, distance_units,
                                 if (is.null(crs)) "data" else crs_arg)
  if (dedupe) coords <- unique(coords)
  as.numeric(pairwise_distances(coords))
}

##' @title Summaries of the distances
##' @description Computes the distances between the unique locations in the dataset and returns summary statistics.
##'
##' @param data an object of class `sf` containing point geometries.
##' @param distance_crs `NULL` to retain an existing projected CRS, or
##' automatically reproject longitude/latitude data to an appropriate UTM CRS
##' (with a message reporting the choice). Alternatively, a CRS to reproject
##' `data` to, which must itself be projected.
##' @param distance_units Character string, either `"km"` or `"m"`, indicating whether the
##' distances used are expressed in kilometers or meters. Defaults to `"km"`.
##'
##' @return a named numeric vector of class `RiskMap_distance_summary`
##' containing the following components
##' \describe{
##'   \item{`min`}{the minimum distance}
##'   \item{`max`}{the maximum distance}
##'   \item{`mean`}{the mean distance}
##'   \item{`median`}{the median distance}
##' }
##' The values are stored in full. Its `print()` method shows them to
##' `digits` significant digits, with the number of unique locations and the
##' distance units.
##'
##' @examples
##' data(italy_sim)
##'
##' summarise_distance(italy_sim)
##'
##' @export
summarise_distance <- function(data,
                           distance_crs = NULL,
                           distance_units = c("km", "m")) {
  check_supplied(data)

  check_data(data)

  distance_units <- match_choice(distance_units)

  d <- extract_distances(data, distance_crs, "distance_crs", "computing distances",
                         distance_units, dedupe = TRUE)

  out <- c(min(d), max(d), mean(d), median(d))
  names(out) <- c("min", "max", "mean", "median")

  # n unique locations give n (n - 1) / 2 distances
  n_locations <- round((1 + sqrt(1 + 8 * length(d))) / 2)
  structure(out, class = "RiskMap_distance_summary",
            n_locations = n_locations, distance_units = distance_units)
}

##' @rdname summarise_distance
##' @param x an object returned by `summarise_distance()`.
##' @param digits the number of significant digits to print. Defaults to 3.
##' @param ... further arguments, currently ignored.
##' @export
print.RiskMap_distance_summary <- function(x, digits = 3, ...) {
  check_positive_integer(digits, "digits")
  values <- formatC(as.numeric(x), digits = digits, format = "fg",
                    big.mark = ",")
  names(values) <- names(x)
  out_print(c(
    out_rule("Distances between locations"),
    out_text(paste0(format_count(attr(x, "n_locations")),
                    " unique locations, with distances in ",
                    attr(x, "distance_units"), ".")),
    out_capture(print(noquote(values), right = TRUE))
  ))
  invisible(x)
}

#' Construct an extreme-rank-length global envelope
#'
#' @param curves Matrix with distance classes in rows and curves in columns;
#'   the observed curve must be the first column.
#' @param level Simultaneous coverage probability.
#' @return The envelope bounds.
#' @noRd
global_rank_envelope <- function(curves, level) {
  number_curves <- ncol(curves)

  # A curve is globally extreme when its sorted vector of two-sided pointwise
  # ranks is lexicographically small. ERL resolves most ties left by the
  # minimum-rank measure while preserving the exchangeability of the curves.
  # apply() drops to a vector when there is a single distance class, so the
  # results are reshaped to one row per curve explicitly
  pointwise_ranks <- matrix(apply(curves, 1L, function(values) {
    pmin(rank(values, ties.method = "average"),
         rank(-values, ties.method = "average"))
  }), nrow = number_curves)
  sorted_ranks <- matrix(apply(pointwise_ranks, 1L, sort),
                         nrow = number_curves, byrow = TRUE)
  ordering <- do.call(order, c(as.data.frame(sorted_ranks),
                               list(method = "radix")))
  ordered_ranks <- sorted_ranks[ordering, , drop = FALSE]

  tied_with_previous <- c(FALSE, apply(ordered_ranks[-1L, , drop = FALSE] ==
                                        ordered_ranks[-number_curves, , drop = FALSE],
                                      1L, all))
  group <- cumsum(!tied_with_previous)
  ordered_upper_position <- ave(seq_len(number_curves), group, FUN = max)
  upper_position <- integer(number_curves)
  upper_position[ordering] <- ordered_upper_position

  # The tolerance stops rounding error from excluding one curve too few, as
  # (1 - 0.9) * 1000 is slightly below 100
  number_excluded <- floor((1 - level) * number_curves + 1e-8)
  retained <- upper_position > number_excluded

  list(lower = apply(curves[, retained, drop = FALSE], 1L, min),
       upper = apply(curves[, retained, drop = FALSE], 1L, max))
}


##' @title Empirical variogram with a global permutation envelope
##' @description Computes the empirical semivariogram using lag-distance
##' breakpoints and, optionally, an extreme-rank-length global envelope for
##' spatial independence.
##' @param data an object of class \code{sf} containing the variable for which the variogram
##' is to be computed and the coordinates
##' @param variable a character indicating the name of variable for which the variogram is to be computed.
##' @param breaks an optional numeric vector of lag-distance breakpoints.
##' If supplied, these breakpoints are used directly to define the distance classes.
##' A warning is raised when some classes lie entirely beyond the largest
##' distance between locations, since they contain no pairs.
##' @param n_bins the number of lag-distance classes to generate when \code{breaks = NULL}.
##' By default \code{n_bins = 14}.
##' @param max_dist an optional maximum lag distance. When \code{breaks = NULL},
##' the breakpoints are generated as \code{seq(0, max_dist, length.out = n_bins + 1)}.
##' By default \code{max_dist = NULL} and the upper lag distance is set to \code{d_max/3},
##' where \code{d_max} is the maximum observed distance in the data.
##' @param n_permutations a non-negative integer indicating the number of random
##' permutations used with the observed curve to construct the global envelope.
##' By default \code{n_permutations = 999}; set it to zero to omit the envelope.
##' @param level the simultaneous coverage probability of the global envelope.
##' By default \code{level = 0.95}.
##' @param seed an optional non-negative integer seed for the permutations. The
##' caller's random-number state is restored before the function returns.
##' @param distance_crs `NULL` to retain an existing projected CRS, or
##' automatically reproject longitude/latitude data to an appropriate UTM CRS
##' (with a message reporting the choice). Alternatively, a CRS to reproject
##' `data` to, which must itself be projected.
##' @param distance_units Character string, either \code{"km"} or \code{"m"}, indicating whether
##' the distances used in the variogram are expressed in kilometers or meters.
##' By default \code{distance_units = "km"}, as in \code{\link{glgpm}}.
##' @details The observed variogram is included among the permuted curves so
##' that they are exchangeable under the null hypothesis of spatial
##' independence. Curves are ordered by their extreme-rank-length ordering.
##' Unlike separate pointwise intervals, the resulting envelope controls the
##' probability that the variogram leaves the envelope anywhere across the lag
##' distances. The envelope is intended as an exploratory diagnostic rather
##' than a formal hypothesis test.
##'
##' @return an object of class `RiskMap_variogram` which is a list containing the following components:
##'   \describe{
##'   \item{variogram}{a data-frame containing the following columns:
##'   \describe{
##'     \item{distance}{the mean pair distance within each lag class}
##'     \item{semivariance}{the observed empirical semivariance}
##'     \item{n_pairs}{the number of pairs in the lag class}}
##'   If \code{n_permutations > 0}, the data-frame also contains the following columns:
##'   \describe{
##'     \item{lower_envelope}{the lower bound of the simultaneous envelope}
##'     \item{upper_envelope}{the upper bound of the simultaneous envelope}
##'   }}
##'   \item{distance_units}{the value passed to \code{distance_units}}
##'   \item{n_permutations}{the number of permutations}
##'   \item{breaks}{the calculated breaks}
##'   \item{level}{the simultaneous coverage probability}
##'   \item{envelope_method}{\code{"global_extreme_rank_length"}, or
##'   \code{NULL} if no permutations were requested}
##'   }
##'
##' @examples
##' data(italy_sim)
##'
##' italy_variogram <- variogram(
##'                      data = italy_sim[1:200,],
##'                      variable = "y",
##'                      n_bins = 10,
##'                      n_permutations = 199,
##'                      seed = 123)
##'
##' plot(italy_variogram)
##'
##' @export
variogram <- function(data,
                      variable,
                      breaks = NULL,
                      n_bins = 14L,
                      max_dist = NULL,
                      n_permutations = 999L,
                      level = 0.95,
                      seed = NULL,
                      distance_crs = NULL,
                      distance_units = c("km", "m")) {
  check_supplied(data)
  check_supplied(variable)

  check_data(data)

  if (!is.character(variable) || length(variable) != 1L ||
      is.na(variable) || !nzchar(variable)) {
    supplied <- describe_value(variable)
    abort_input("{.arg variable} must be the name of one column of {.arg data}, not {supplied}.", arg = "variable")
  }
  if (!variable %in% names(data)) {
    abort_input("{.arg variable} must name a column of {.arg data}, and {.val {variable}} doesn't.", arg = "variable")
  }
  if (!is.null(breaks) && !missing(n_bins)) {
    abort_input("{.arg breaks} and {.arg n_bins} can't both be given.", arg = "breaks")
  }
  if (!is.null(breaks) && !is.null(max_dist)) {
    abort_input("{.arg breaks} and {.arg max_dist} can't both be given.", arg = "breaks")
  }
  values <- data[[variable]]
  if (!is.numeric(values) || any(!is.finite(values))) {
    abort_input("The column {.val {variable}} must hold only finite numbers.", arg = "variable")
  }
  if (nrow(data) < 2L) {
    abort_input("{.arg data} must have at least two locations.", arg = "data")
  }

  check_positive_integer(n_bins, "n_bins")
  if (!is.null(max_dist)) {
    check_range(max_dist, min = 0, max = Inf, allow_equal = FALSE,
                arg = "max_dist")
  }
  check_positive_integer(n_permutations, "n_permutations", allow_zero = TRUE)
  check_range(level, min = 0, max = 1, allow_equal = FALSE, arg = "level")
  if (n_permutations > 0L &&
      n_permutations + 1L < ceiling(1 / (1 - level))) {
    needed <- ceiling(1 / (1 - level)) - 1L
    abort_input("{.arg n_permutations} must be at least {needed} for a global envelope at {.arg level} = {level}, not {n_permutations}.", arg = "n_permutations")
  }
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)
  distance_units <- match_choice(distance_units)

  d <- extract_distances(data, distance_crs, "distance_crs", "computing distances",
                        distance_units)
  if (!length(d) || !is.finite(max(d)) || max(d) <= 0) {
    abort_input("{.arg data} must have at least two distinct locations.", arg = "data")
  }

  if (is.null(breaks)) {
    upper_dist <- if (is.null(max_dist)) max(d) / 3 else max_dist
    breaks <- seq(0, upper_dist, length.out = n_bins + 1)
  } else {
    if (!is.numeric(breaks) || length(breaks) < 2L || any(!is.finite(breaks))) {
      abort_input("{.arg breaks} must be at least two finite numbers.", arg = "breaks")
    }
    if (any(diff(breaks) <= 0)) {
      abort_input("{.arg breaks} must be strictly increasing.", arg = "breaks")
    }
    if (min(breaks) < 0) {
      abort_input("{.arg breaks} must be non-negative.", arg = "breaks")
    }
  }
  distance_class <- as.integer(cut(d, breaks = breaks,
                                   include.lowest = TRUE, right = TRUE))
  included <- !is.na(distance_class)
  if (!any(included)) {
    range_d <- signif(range(d), 3)
    abort_input(c("No pair of locations falls within the distance classes.", "i" = "Distances between locations range from {range_d[1]} to {range_d[2]} {distance_units}.", "i" = "Check {.arg breaks}, {.arg max_dist} and {.arg distance_units}."), arg = "breaks")
  }
  # Classes are closed on the right, so a class whose lower limit is at or
  # beyond the largest distance cannot contain any pair
  number_beyond <- sum(breaks[-length(breaks)] >= max(d))
  if (number_beyond > 0L) {
    largest <- signif(max(d), 3)
    riskmap_warn(c("{number_beyond} distance class{?es} {?lies/lie} beyond the largest distance between locations, {largest} {distance_units}, and {cli::qty(number_beyond)}{?has/have} no pairs.", "i" = "Check {.arg breaks}, {.arg max_dist} and {.arg distance_units}."))
  }

  d <- d[included]
  distance_class <- distance_class[included]
  number_bins <- length(breaks) - 1L
  number_pairs <- length(d)
  n_pairs <- tabulate(distance_class, nbins = number_bins)
  distance_sum <- numeric(number_bins)
  summed_distances <- rowsum(d, distance_class, reorder = FALSE)
  distance_sum[as.integer(rownames(summed_distances))] <- summed_distances[, 1L]
  mean_distance <- distance_sum / n_pairs
  mean_distance[n_pairs == 0L] <- NA_real_

  n <- nrow(data)
  first_index <- integer(n * (n - 1L) / 2L)
  second_index <- integer(length(first_index))
  position <- 1L
  for (second in seq_len(n - 1L)) {
    length_block <- n - second
    indices <- position:(position + length_block - 1L)
    first_index[indices] <- (second + 1L):n
    second_index[indices] <- second
    position <- position + length_block
  }

  first_index <- first_index[included] - 1L
  second_index <- second_index[included] - 1L
  bin_index <- distance_class - 1L

  if (!is.null(seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(seed)
  }
  permutations <- matrix(seq_len(n), ncol = 1L)
  if (n_permutations > 0L) {
    permutations <- cbind(
      permutations,
      replicate(n_permutations, sample.int(n), simplify = "matrix")
    )
  }
  permuted_values <- matrix(values[permutations], nrow = n)
  curves <- cpp_binned_semivariances(permuted_values, first_index,
                                     second_index, bin_index, number_bins)

  variogram_data <- data.frame(
    distance = mean_distance,
    semivariance = curves[, 1L],
    n_pairs = n_pairs
  )
  envelope_method <- NULL
  if (n_permutations > 0L) {
    nonempty <- n_pairs > 0L
    envelope <- global_rank_envelope(curves[nonempty, , drop = FALSE], level)
    variogram_data$lower_envelope <- NA_real_
    variogram_data$upper_envelope <- NA_real_
    variogram_data$lower_envelope[nonempty] <- envelope$lower
    variogram_data$upper_envelope[nonempty] <- envelope$upper
    envelope_method <- "global_extreme_rank_length"
  }

  result <- list(
    variogram = variogram_data,
    distance_units = distance_units,
    n_permutations = n_permutations,
    breaks = breaks,
    level = level,
    envelope_method = envelope_method
  )

  class(result) <- "RiskMap_variogram"
  return(result)
}

##' @title Print Method for RiskMap_variogram Objects
##' @description Prints the distance classes and the envelope settings of an
##' empirical variogram computed by \code{\link{variogram}}, followed by its
##' table.
##' @param x an object of class \code{RiskMap_variogram}.
##' @param digits the number of significant digits to print. By default
##' \code{digits = 3}.
##' @param ... further arguments passed to \code{\link{print.data.frame}}.
##' @return \code{x}, invisibly.
##' @seealso \code{\link{variogram}}
##' @method print RiskMap_variogram
##' @export
print.RiskMap_variogram <- function(x, digits = 3, ...) {
  out_print(format_variogram(x, digits, ...))
  invisible(x)
}

##' The lines that print an empirical variogram
##' @noRd
format_variogram <- function(x, digits, ...) {
  classes <- nrow(x$variogram)
  envelope <- if (x$n_permutations > 0L) {
    paste0("The ", 100 * x$level, "% global envelope for spatial ",
           "independence is from ", format_count(x$n_permutations),
           " permutations.")
  } else {
    "There is no envelope, since `n_permutations` is 0."
  }
  outside <- if (x$n_permutations > 0L) {
    # A class without pairs has no semivariance and no envelope, so it is
    # neither inside nor outside.
    table <- x$variogram
    tested <- is.finite(table$semivariance) &
      is.finite(table$lower_envelope) & is.finite(table$upper_envelope)
    beyond <- sum(tested & (table$semivariance < table$lower_envelope |
                              table$semivariance > table$upper_envelope))
    untested <- classes - sum(tested)
    c(i = paste0("The variogram is outside the envelope in ", beyond, " of ",
                 sum(tested), " distance classes",
                 if (untested > 0L) {
                   paste0(", and ", untested,
                          if (untested == 1L) " class has" else
                            " classes have", " no envelope")
                 }, "."))
  }
  # Counts are printed in full, whatever the number of digits
  variogram_table <- x$variogram
  variogram_table$n_pairs <- as.integer(variogram_table$n_pairs)
  c(out_rule("Empirical variogram"),
    out_text(paste0(classes, " distance classes from ",
                    signif(min(x$breaks), digits), " to ",
                    signif(max(x$breaks), digits), " ", x$distance_units,
                    ". ", envelope)),
    "",
    out_capture(print(variogram_table, digits = digits, row.names = FALSE,
                      ...)),
    if (!is.null(outside)) "",
    out_bullets(outside))
}

##' @title Plot Method for RiskMap_variogram Objects
##' @description Plots the empirical variogram generated by \code{\link{variogram}}
##' @param x The output generated by the function \code{\link{variogram}}.
##' @param plot_envelope A logical value indicating if the global envelope of spatial independence
##' generated using the permutation test must be displayed (\code{plot_envelope = TRUE}) or not
##' (\code{plot_envelope = FALSE}). By default \code{plot_envelope = TRUE}. Note: if
##' \code{n_permutations} was 0 when running \code{\link{variogram}}, no envelope
##' can be generated; a warning is raised and the envelope is skipped.
##' @param color If \code{plot_envelope = TRUE}, it sets the fill colour of the envelope; run \code{vignette("ggplot2-specs")} for more details on this argument.
##' @param show_pairs A logical value indicating if each point should be
##' labelled with the number of pairs of locations in its distance class. A
##' semivariance estimated from few pairs is unstable, so the labels help to
##' judge how much each point can be trusted, and to choose the distance
##' classes in \code{\link{variogram}}. By default \code{show_pairs = FALSE}.
##' @param ... Currently unused.
##' @return A \code{ggplot} object representing the empirical variogram plot, optionally including the envelope of spatial independence.
##' @details This function plots the empirical variogram, which shows the spatial dependence structure of the data. Each point is the semivariance of one distance class, and classes without pairs are omitted. Points are not joined by lines, because the empirical variogram is only defined for the distance classes. Both axes start at zero. If \code{plot_envelope} is set to \code{TRUE}, the plot also includes a simultaneous extreme-rank-length envelope for spatial independence, and points that fall outside it are drawn in red.
##' @seealso \code{\link{variogram}}
##' @examples
##' data(galicia)
##' galicia_variogram <- variogram(galicia, variable = "lead",
##'                                distance_units = "km", seed = 1)
##' plot(galicia_variogram, show_pairs = TRUE)
##' @method plot RiskMap_variogram
##' @export
plot.RiskMap_variogram <- function(x,
                                   plot_envelope = TRUE,
                                   color = "royalblue1",
                                   show_pairs = FALSE, ...) {
  check_logical(plot_envelope)
  check_logical(show_pairs)

  if (plot_envelope && x$n_permutations == 0L){
    riskmap_warn(c("The variogram is plotted without an envelope for spatial independence, because {.fn variogram} was run with {.code n_permutations = 0}.", "i" = "Increase {.arg n_permutations}, or set {.code plot_envelope = FALSE}."))
    plot_envelope <- FALSE
  }

  # Classes without pairs have no semivariance to show
  variogram_data <- x$variogram[x$variogram$n_pairs > 0L, , drop = FALSE]

  basic_plot <- ggplot(data = variogram_data,
                       aes(x = .data$distance, y = .data$semivariance))

  point_colours <- rep("black", nrow(variogram_data))
  if (plot_envelope) {
    basic_plot <- basic_plot +
      geom_ribbon(aes(ymin = .data$lower_envelope,
                      ymax = .data$upper_envelope),
                  fill = color, alpha = 0.3)
    outside <- variogram_data$semivariance < variogram_data$lower_envelope |
      variogram_data$semivariance > variogram_data$upper_envelope
    point_colours[outside] <- "firebrick"
  }

  basic_plot <- basic_plot + geom_point(colour = point_colours, size = 2)

  if (show_pairs) {
    # Labels alternate above and below the points so that neighbouring
    # labels do not overlap, and the top of the axis gets extra room.
    labels <- variogram_data
    labels$label <- format(labels$n_pairs, big.mark = ",", trim = TRUE)
    labels$vjust <- ifelse(seq_len(nrow(labels)) %% 2L == 1L, -0.9, 1.9)
    basic_plot <- basic_plot +
      geom_text(data = labels,
                aes(label = .data$label, vjust = .data$vjust),
                size = 2.4, colour = "grey30") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
  }

  x_label <- sprintf("Distance (%s)", x$distance_units)

  basic_plot +
    expand_limits(x = 0, y = 0) +
    labs(x = x_label, y = "Semivariance") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
}
