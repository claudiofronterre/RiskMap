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
#' @noRd
extract_distances <- function(data, crs, crs_arg, purpose, distance_units,
                              dedupe = FALSE) {
  if (!is.null(crs)) {
    check_crs(crs, name = crs_arg)
    data <- st_transform(data, crs = crs)
    if (st_is_longlat(data)) {
      stop("'", crs_arg, "' must be a projected CRS, not longitude/latitude", call. = FALSE)
    }
  } else if (st_is_longlat(data)) {
    auto_crs <- propose_utm(data)
    data <- st_transform(data, crs = auto_crs)
    message("'data' are in longitude/latitude and '", crs_arg, "' was not provided; ",
            "automatically reprojecting to EPSG:", auto_crs,
            " for ", purpose, ". Set '", crs_arg, "' to override.")
  }
  coords <- coordinates_in_units(data, distance_units)
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
##' @return a named vector containing the following components
##' \describe{
##'   \item{`min`}{the minimum distance}
##'   \item{`max`}{the maximum distance}
##'   \item{`mean`}{the mean distance}
##'   \item{`median`}{the minimum distance}
##' }
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

  check_data(data)

  stopifnot("'distance_units' must be either 'km' or 'm'" =
              is.character(distance_units) && all(distance_units %in% c("km", "m")))
  distance_units <- match.arg(distance_units)

  d <- extract_distances(data, distance_crs, "distance_crs", "computing distances",
                         distance_units, dedupe = TRUE)

  out <- c(min(d), max(d), mean(d), median(d))
  names(out) <- c("min", "max", "mean", "median")

  return(out)
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
  pointwise_ranks <- apply(curves, 1L, function(values) {
    pmin(rank(values, ties.method = "average"),
         rank(-values, ties.method = "average"))
  })
  if (is.null(dim(pointwise_ranks))) {
    pointwise_ranks <- matrix(pointwise_ranks, ncol = 1L)
  }
  sorted_ranks <- t(apply(pointwise_ranks, 1L, sort))
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

  number_excluded <- floor((1 - level) * number_curves)
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
##' By default \code{distance_units = "m"}
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
##'.  \describe{
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
##' plot_variogram(italy_variogram,
##'                plot_envelope = TRUE)
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
                      distance_units = c("m", "km")) {

  check_data(data)

  if (!is.character(variable) || length(variable) != 1L ||
      is.na(variable) || !nzchar(variable)) {
    stop("'variable' must be a single object of class 'character'")
  }
  if (!variable %in% names(data)){
    stop("'variable' must be one of the columns in 'data'")
  }
  if (!is.null(breaks) && !missing(n_bins)){
    stop("'breaks' and 'n_bins' cannot both be supplied")
  }
  if (!is.null(breaks) && !missing(max_dist)){
    stop("'breaks' and 'max_dist' cannot both be supplied")
  }
  values <- data[[variable]]
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop("'variable' must contain only finite numeric values", call. = FALSE)
  }
  if (nrow(data) < 2L) {
    stop("'data' must contain at least two locations", call. = FALSE)
  }

  check_positive_integer(n_bins, "n_bins")
  if (!is.null(max_dist)) {
    check_range(max_dist, min = 0, max = Inf, allow_equal = FALSE,
                name = "max_dist")
  }
  check_positive_integer(n_permutations, "n_permutations", allow_zero = TRUE)
  check_range(level, min = 0, max = 1, allow_equal = FALSE, name = "level")
  if (n_permutations > 0L &&
      n_permutations + 1L < ceiling(1 / (1 - level))) {
    stop("'n_permutations' is too small to construct a global envelope at ",
         "the requested 'level'", call. = FALSE)
  }
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)
  stopifnot("'distance_units' must be either 'km' or 'm'" =
              is.character(distance_units) && all(distance_units %in% c("km", "m")))
  distance_units <- match.arg(distance_units)

  d <- extract_distances(data, distance_crs, "distance_crs", "computing distances",
                        distance_units)
  if (!length(d) || !is.finite(max(d)) || max(d) <= 0) {
    stop("'data' must contain at least two distinct locations", call. = FALSE)
  }

  if (is.null(breaks)) {
    upper_dist <- ifelse(is.null(max_dist), max(d) / 3, max_dist)
    breaks <- seq(0, upper_dist, length.out = n_bins + 1)
  } else {
    if (!is.numeric(breaks) || length(breaks) < 2L || any(!is.finite(breaks))) {
      stop("'breaks' must be a numeric vector with at least two values")
    }
    if (any(diff(breaks) <= 0)) {
      stop("'breaks' must be strictly increasing")
    }
    if (min(breaks) < 0) {
      stop("'breaks' must be non-negative")
    }
    upper_dist <- max(breaks)
  }
  if (upper_dist > max(d)){
    stop("the provided lag distances go beyond the maximum observed distance")
  }
  distance_class <- as.integer(cut(d, breaks = breaks,
                                   include.lowest = TRUE, right = TRUE))
  included <- !is.na(distance_class)
  if (!any(included)) {
    stop("the provided lag distances do not match the
          scale of the observed distances; consider setting distance_units = 'km'")
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

##' @title Plotting the empirical variogram
##' @description Plots the empirical variogram generated by \code{\link{variogram}}
##' @param variogram_output The output generated by the function \code{\link{variogram}}.
##' @param plot_envelope A logical value indicating if the global envelope of spatial independence
##' generated using the permutation test must be displayed (\code{plot_envelope = TRUE}) or not
##' (\code{plot_envelope = FALSE}). By default \code{plot_envelope = TRUE}. Note: if
##' \code{n_permutations} was 0 when running \code{\link{variogram}}, no envelope
##' can be generated; a warning is raised and the envelope is skipped.
##' @param color If \code{plot_envelope = TRUE}, it sets the colour of the envelope; run \code{vignette("ggplot2-specs")} for more details on this argument.
##' @return A \code{ggplot} object representing the empirical variogram plot, optionally including the envelope of spatial independence.
##' @details This function plots the empirical variogram, which shows the spatial dependence structure of the data. If \code{plot_envelope} is set to \code{TRUE}, the plot also includes a simultaneous extreme-rank-length envelope for spatial independence.
##' @seealso \code{\link{variogram}}
##' @export
plot_variogram <- function(variogram_output,
                           plot_envelope = TRUE,
                           color = "royalblue1") {

  if (!inherits(variogram_output, "RiskMap_variogram")){
    stop("'variogram' must be an object of class 'RiskMap_variogram'")
  }

  if (plot_envelope && variogram_output$n_permutations == 0L){
    warning("No envelope for spatial independence can be plotted because 'n_permutations' ",
            "was ", variogram_output$n_permutations, " when 'variogram()' was run; ",
            "plotting without the envelope. Increase 'n_permutations' or set ",
            "plot_envelope = FALSE to silence this warning.")
    plot_envelope <- FALSE
  }

  basic_plot <- ggplot(data = variogram_output$variogram,
                       aes(x = .data$distance, y = .data$semivariance))

  if (plot_envelope) {
    basic_plot <- basic_plot +
      geom_ribbon(aes(ymin = .data$lower_envelope,
                      ymax = .data$upper_envelope),
                  fill = color, alpha = 0.3)
  }

  basic_plot <- basic_plot + geom_point() + geom_line()

  x_label <- sprintf("Distance (%s)", variogram_output$distance_units)

  basic_plot + labs(x = x_label, y = "Semivariance")
}
