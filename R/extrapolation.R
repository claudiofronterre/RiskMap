##' Assess where predictions extrapolate beyond the data
##'
##' Computes, for each prediction location, indicators of whether a prediction
##' from a fitted geostatistical model relies on extrapolation, either in space
##' (where the data say little about the spatial effect) or in the covariates
##' (outside the conditions observed in the data).
##'
##' @param object A fitted `RiskMap` model.
##' @param newdata An `sf` object of prediction points with the covariates used
##'   by the model, typically a grid from [create_grid()].
##' @param spatial_threshold Threshold of the share of the spatial variance
##'   left unexplained above which a location is classed as spatial
##'   extrapolation. `NULL` (the default) uses the 95th percentile of the
##'   shares left at the sampled locations when their own data are left out
##'   (see Details). A number strictly between 0 and 1 gives a fixed
##'   threshold.
##' @param sampling_noise Whether the share of the spatial variance left
##'   unexplained includes the sampling noise of the outcomes (the default).
##'   With `FALSE`, it depends only on where the data were collected and on
##'   the nugget, as if the spatial effect had been observed without noise at
##'   every sampled location.
##'
##' @details
##' Covariates are assessed on the model's scale, that is, as the columns of
##' the fixed-effect design matrix. A model with `log(elevation)` is therefore
##' assessed on log elevation, and a spline in age through each of its basis
##' columns. For each column, the distance outside the range observed in the
##' data is expressed as a share of that range.
##'
##' The leverage \eqn{h(x) = d(x)^\top (D^\top D)^{-1} d(x)}, where \eqn{D} is
##' the fitted design matrix, is reported relative to its maximum over the
##' data. Values above 1 indicate a combination of covariates that is more
##' unusual than any in the data, even when each covariate is within its range.
##'
##' The share of the spatial variance left unexplained is
##' \deqn{v(x) = 1 - \rho(x)^\top (R + \sigma^{-2} W)^{-1} \rho(x),}
##' where \eqn{\rho(x)} holds the estimated correlations between \eqn{x} and
##' the sampled locations, \eqn{R} their correlation matrix, \eqn{\sigma^2}
##' the spatial variance, and \eqn{W} a diagonal matrix with
##' \eqn{w_j = \tau^2 + c_j}, the nugget variance plus the sampling variance
##' of the linear predictor at the sampled location \eqn{j}. For a Gaussian
##' model, \eqn{c_j} is the measurement-error variance divided by the number
##' of records at the location, and \eqn{v(x)} is
##' \eqn{\mathrm{Var}[S(x) \mid y] / \sigma^2} at the estimated parameters. For
##' binomial, Poisson and negative binomial models, \eqn{c_j} is the inverse of
##' the information \eqn{-\partial^2 \log f(y \mid \eta) / \partial \eta^2}
##' of the records at the location, summed over them and evaluated at the
##' conditional mode of the spatial effect given the data at the estimated
##' parameters, so that \eqn{v(x)} is the Laplace approximation of
##' \eqn{\mathrm{Var}[S(x) \mid y] / \sigma^2}. With `sampling_noise = FALSE`,
##' \eqn{c_j = 0}.
##'
##' \eqn{v(x)} is close to 0 where the data determine the spatial effect and
##' close to 1 where they carry no information about it. Unlike the distance
##' to the nearest sampled location, it accounts for the estimated scale of
##' spatial correlation, for clustering of the sampled locations and for how
##' much the data at each location say about the spatial effect, which is
##' little where few individuals were sampled or where a prevalence is close
##' to 0 or 1. The variances of `re()` terms are not included in \eqn{W}, so
##' for models with `re()` terms \eqn{v(x)} is lower than the share they
##' leave.
##'
##' A location is classed as covariate extrapolation when any covariate lies
##' outside its observed range or its relative leverage exceeds 1, and as
##' spatial extrapolation when \eqn{v(x)} exceeds the threshold. By default
##' the threshold is the 95th percentile over the sampled locations
##' \eqn{s_j} of
##' \deqn{v_{-j} = \mathrm{Var}[S(s_j) \mid y_{-j}] / \sigma^2
##'   = 1 / [(R + \sigma^{-2} W)^{-1}]_{jj} - w_j / \sigma^2,}
##' the share left at \eqn{s_j} by the data of the other locations.
##' Leave-one-out cross-validation predicts at each sampled location from the
##' others, so the \eqn{v_{-j}} are the amounts of spatial information under
##' which predictions can be checked against held-out data, and a location
##' above their 95th percentile has less than nearly all of them. A percentile
##' rather than the largest value keeps a single remote sampled location from
##' setting the threshold. The percentile is 1 only when the data of the
##' other locations say nothing about the spatial effect at more than 5% of
##' the sampled locations, as when the locations share almost no spatial
##' correlation. No location can then be classed as spatial extrapolation, so
##' `assess_extrapolation()` warns, and a fixed `spatial_threshold` is more
##' informative. The two
##' classes are kept separate because they call for different responses: more
##' sampling where data are sparse, and caution about the form of the
##' covariate effects outside the observed conditions.
##'
##' @return An `sf` object of class `RiskMap_extrapolation` with the rows and
##'   geometry of `newdata` and the columns
##'   \describe{
##'     \item{covariates_outside}{Number of covariates outside their observed
##'       range.}
##'     \item{covariate_novelty}{Largest distance outside an observed range,
##'       as a share of that range (0 when all are within range).}
##'     \item{most_novel_covariate}{Name of the design-matrix column with the
##'       largest distance outside its range, or `NA`.}
##'     \item{covariate_leverage}{Leverage relative to its maximum over the
##'       data.}
##'     \item{spatial_unexplained}{Share of the spatial variance left
##'       unexplained by the data, \eqn{v(x)}.}
##'     \item{extrapolation}{Factor with levels `"Supported"`, `"Spatial"`,
##'       `"Covariate"` and `"Spatial and covariate"`.}
##'   }
##'   Rows with missing covariates have missing indicators. The attribute
##'   `spatial_threshold` holds the threshold used, and `covariate_ranges` the
##'   range of each covariate in the data. Map the result with
##'   [plot.RiskMap_extrapolation()].
##' @seealso [plot.RiskMap_extrapolation()], [predict.RiskMap()]
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:150, ]
##' fit <- glgpm(y ~ log(pop_dens) + gp(), data = example_data,
##'              family = "gaussian", messages = FALSE)
##' extrapolation <- assess_extrapolation(fit, italy_sim[151:200, ])
##' extrapolation
##' table(extrapolation$extrapolation)
##' @export
assess_extrapolation <- function(object, newdata, spatial_threshold = NULL,
                                 sampling_noise = TRUE) {
  check_supplied(object)
  check_supplied(newdata)
  if (!inherits(object, "RiskMap")) {
    abort_input("{.arg object} must be a {.cls RiskMap} fit, not {.obj_type_friendly {object}}.", arg = "object")
  }
  if (!is.null(spatial_threshold)) {
    check_range(spatial_threshold, min = 0, max = 1, allow_equal = FALSE,
                arg = "spatial_threshold")
  }
  check_logical(sampling_noise, arg = "sampling_noise")
  prediction <- validate_prediction_data(object, newdata, include = character())
  model_data <- prediction$model_data
  valid_rows <- prediction$valid_rows

  covariates <- extrapolation_covariates(
    object, prediction_design_matrix(object, st_drop_geometry(model_data))
  )
  spatial <- extrapolation_spatial_unexplained(object, model_data,
                                               sampling_noise)
  spatial_unexplained <- spatial$unexplained
  threshold_given <- !is.null(spatial_threshold)
  if (!threshold_given) {
    spatial_threshold <- stats::quantile(spatial$left_out, 0.95,
                                         names = FALSE)
    if (spatial_threshold >= 1 - 1e-8) {
      riskmap_warn(c(
        "The data leave all of the spatial variance unexplained at more than 5% of the sampled locations when their own data are left out, so the threshold is 1 and no location can be classed as spatial extrapolation.",
        "i" = "Give a number in {.arg spatial_threshold} for a fixed threshold."
      ))
    }
  }
  # Both comparisons allow for rounding, so that a location with no spatial
  # information is not flagged against a threshold of 1 on one platform and
  # supported on another.
  covariate_flag <- covariates$outside > 0 | covariates$leverage > 1 + 1e-8
  spatial_flag <- spatial_unexplained > spatial_threshold + 1e-8
  classes <- c("Supported", "Spatial", "Covariate", "Spatial and covariate")
  extrapolation <- factor(
    classes[1L + spatial_flag + 2L * covariate_flag], levels = classes
  )

  restore <- function(values) {
    out <- rep(values[1L][NA], length(valid_rows))
    out[valid_rows] <- values
    out
  }
  out <- prediction$output_data
  out$covariates_outside <- restore(covariates$outside)
  out$covariate_novelty <- restore(covariates$novelty)
  out$most_novel_covariate <- restore(covariates$most_novel)
  out$covariate_leverage <- restore(covariates$leverage)
  out$spatial_unexplained <- restore(spatial_unexplained)
  out$extrapolation <- factor(restore(as.character(extrapolation)),
                              levels = classes)
  attr(out, "grid_metadata") <- attr(newdata, "grid_metadata", exact = TRUE)
  attr(out, "spatial_threshold") <- spatial_threshold
  attr(out, "spatial_rule") <- if (threshold_given) "given" else
    "leave_one_out"
  attr(out, "sampling_noise") <- sampling_noise
  attr(out, "covariate_ranges") <- covariates$ranges
  class(out) <- c("RiskMap_extrapolation", class(out))
  stamp_version(out)
}

#' Covariate range checks and leverage on the model's scale
#' @noRd
extrapolation_covariates <- function(object, design) {
  observed <- as.matrix(object$D)
  varying <- apply(observed, 2L, function(column) diff(range(column)) > 0)
  ranges <- apply(observed[, varying, drop = FALSE], 2L, range)
  n <- nrow(design)
  novelty <- matrix(0, n, sum(varying))
  if (any(varying)) {
    values <- design[, varying, drop = FALSE]
    width <- ranges[2L, ] - ranges[1L, ]
    below <- sweep(-values, 2L, ranges[1L, ], "+")
    above <- sweep(values, 2L, ranges[2L, ], "-")
    novelty <- sweep(pmax(below, above, 0), 2L, width, "/")
  }
  most <- if (any(varying)) max.col(novelty, ties.method = "first") else
    rep(NA_integer_, n)
  largest <- if (any(varying)) novelty[cbind(seq_len(n), most)] else
    rep(0, n)
  most_novel <- ifelse(largest > 0, colnames(observed)[varying][most],
                       NA_character_)

  inverse <- solve(crossprod(observed))
  leverage <- rowSums((design %*% inverse) * design)
  observed_leverage <- rowSums((observed %*% inverse) * observed)
  list(
    outside = as.integer(rowSums(novelty > 0)),
    novelty = largest,
    most_novel = most_novel,
    leverage = leverage / max(observed_leverage),
    ranges = ranges
  )
}

#' Share of the spatial variance left unexplained by the data
#'
#' Returns `unexplained`, v(x) at each prediction location, and `left_out`,
#' v_{-j} at each sampled location from the data of the others, both from the
#' factor of R + W / sigma2 (`assess_extrapolation()`). The ratios W / sigma2
#' are capped at 1e8, so that a location whose data carry almost no
#' information leaves 1 / (K^{-1})_jj - w_j / sigma2 computable.
#' @noRd
extrapolation_spatial_unexplained <- function(object, model_data,
                                              sampling_noise = TRUE) {
  parameters <- coef(object)
  noise <- if (sampling_noise) {
    extrapolation_sampling_variance(object, parameters)
  } else {
    0
  }
  ratio <- pmin((fitted_nugget_variance(object, parameters) + noise) /
                  parameters$sigma2, 1e8)

  observed <- as.matrix(object$coords)
  correlation <- matern_correlation(as.matrix(stats::dist(observed)),
                                    phi = parameters$phi, kappa = object$kappa)
  diag(correlation) <- 1 + ratio
  factor <- factor_covariance(correlation,
                              "correlation matrix of the sampled locations")
  left_out <- 1 / diag(chol2inv(factor)) - ratio - attr(factor, "jitter")

  locations <- coordinates_in_units(st_geometry(model_data),
                                    object$distance_units)
  unexplained <- numeric(nrow(locations))
  blocks <- split(seq_len(nrow(locations)),
                  ceiling(seq_len(nrow(locations)) / 2000))
  for (rows in blocks) {
    distances <- sqrt(pmax(
      outer(rowSums(locations[rows, , drop = FALSE]^2),
            rowSums(observed^2), "+") -
        2 * locations[rows, , drop = FALSE] %*% t(observed),
      0
    ))
    cross <- matern_correlation(distances, phi = parameters$phi,
                                kappa = object$kappa)
    weights <- backsolve(factor, t(cross), transpose = TRUE)
    unexplained[rows] <- 1 - colSums(weights^2)
  }
  list(unexplained = pmin(pmax(unexplained, 0), 1),
       left_out = pmin(pmax(left_out, 0), 1))
}

#' Sampling variance of the linear predictor at each sampled location
#'
#' For a Gaussian fit, the measurement-error variance divided by the number of
#' records at the location. Otherwise the inverse of the summed curvature
#' -d^2 log f(y | eta) / d eta^2 of the records at the location, at the
#' conditional mode of the location and `re()` effects given the data at the
#' estimates. The mode is found as `laplace_sampling_mcmc()` finds the centre
#' of its sampler, on the covariance that prediction conditions on
#' (`conditioning_covariance()`), with the family's records
#' (`engine_records()`).
#' @noRd
extrapolation_sampling_variance <- function(object, parameters) {
  n_locations <- nrow(object$coords)
  id_coords <- object$ID_coords
  if (fit_is_exact(object)) {
    measurement <- parameters$sigma2_me %||% fit_fix_var_me(object) %||% 0
    return(measurement / tabulate(id_coords, n_locations))
  }
  distance_scale <- attr(object, "distance_scale") %||% 1
  covariance <- conditioning_covariance(
    object, parameters, pairwise_distances(object$coords / distance_scale),
    distance_scale
  )
  invlink <- if (identical(object$link_function$name, "canonical")) {
    NULL
  } else {
    object$link_function[c("inv", "d1", "d2")]
  }
  id_re <- if (!is.null(object$ID_re)) as.matrix(object$ID_re)
  sigma2_re <- if (!is.null(id_re)) parameters$sigma2_re
  mu <- as.numeric(as.matrix(object$D) %*% parameters$beta + object$cov_offset)
  records <- engine_records(object$family, object$y, object$units_m, mu,
                            fit_family_values(object), invlink)
  linkf <- make_inverse_link_functions(records$family, invlink,
                                       length(records$y))
  structure <- sampler_latent_structure(covariance, id_coords, id_re,
                                        sigma2_re)
  found <- tryCatch(
    find_sampler_mode(structure, records$y, records$units_m, records$offset,
                      records$family, linkf),
    error = function(error) NULL
  )
  if (is.null(found)) {
    found <- fallback_sampler_mode(structure, records$y, records$units_m,
                                   records$offset, covariance, id_coords,
                                   id_re, sigma2_re, records$family, invlink,
                                   linkf)
  }
  information <- numeric(n_locations)
  sums <- rowsum(found$weights, id_coords)
  information[as.integer(rownames(sums))] <- sums[, 1L]
  1 / pmax(information, 0)
}

##' @rdname assess_extrapolation
##' @param x A `RiskMap_extrapolation` object, whose `print()` method counts
##'   the locations in each class and states the thresholds of the classes,
##'   returning `x` invisibly.
##' @param ... Not used.
##' @method print RiskMap_extrapolation
##' @export
print.RiskMap_extrapolation <- function(x, ...) {
  threshold <- attr(x, "spatial_threshold", exact = TRUE)
  # A subset that lost the classes or their thresholds prints as an sf object.
  if (!is.factor(x[["extrapolation"]]) || !is.numeric(threshold)) {
    return(NextMethod())
  }
  out_print(format_extrapolation(x, threshold))
  invisible(x)
}

#' The lines that print indicators of extrapolation
#' @noRd
format_extrapolation <- function(x, threshold) {
  classes <- x[["extrapolation"]]
  counts <- table(classes)
  missing <- sum(is.na(classes))
  values <- format(c(as.integer(counts), if (missing) missing),
                   big.mark = ",")
  labels <- c(names(counts), if (missing) "Missing covariates")
  share <- format_extrapolation_share(threshold)
  given <- identical(attr(x, "spatial_rule"), "given")
  isolated <- !given && threshold >= 1 - 1e-8
  spatial <- if (given) {
    paste0("Spatial: more than ", share, " of the spatial variance left ",
           "unexplained, the `spatial_threshold` given.")
  } else if (isolated) {
    NULL
  } else {
    paste0("Spatial: more than ", share, " of the spatial variance left ",
           "unexplained, the 95th percentile of the shares at the sampled ",
           "locations when their own data are left out.")
  }
  isolated_note <- if (isolated) {
    c("!" = paste0("Spatial: none, since the threshold is ", share, ", the ",
                   "95th percentile of the shares at the sampled locations ",
                   "when their own data are left out. The data leave all of ",
                   "the spatial variance unexplained at more than 5% of ",
                   "them. A number in `spatial_threshold` gives a fixed ",
                   "threshold."))
  }
  noise <- if (isFALSE(attr(x, "sampling_noise"))) {
    c(i = "The share leaves out the sampling noise of the outcomes.")
  }
  c(out_rule("Extrapolation"),
    out_text(paste0("Indicators at ", out_count(length(classes), "location"),
                    ".")),
    "",
    out_pairs(labels, values),
    "",
    out_bullets(c(
      if (!is.null(spatial)) c(i = spatial),
      isolated_note,
      i = paste0("Covariate: a covariate outside its range in the data, or ",
                 "a leverage above its largest value there."),
      noise,
      i = paste0("`plot(x)` maps the classes, and ",
                 "`indicator = \"spatial_unexplained\"` the share.")
    )))
}

#' A threshold as a percentage, with the digits that keep a value below 1
#' from showing as 100%
#' @noRd
format_extrapolation_share <- function(threshold) {
  for (digits in 3:8) {
    shown <- signif(100 * threshold, digits)
    if (shown < 100 || threshold >= 1) break
  }
  paste0(format(shown, digits = digits), "%")
}

##' Map indicators of extrapolation
##'
##' Maps indicators returned by [assess_extrapolation()], in the style of
##' [plot.RiskMap_prediction()].
##'
##' @param x A `RiskMap_extrapolation` object.
##' @param indicator Indicator to map: `"extrapolation"` (the default, the
##'   four classes, with the thresholds that [assess_extrapolation()] used),
##'   `"spatial_unexplained"` (the share of the spatial variance left
##'   unexplained by the data, on a scale from 0 to 1), `"covariate_novelty"`
##'   or `"covariate_leverage"`. Several indicators are arranged in two
##'   columns.
##' @param ... Not used.
##' @return A `ggplot` object, to which further layers can be added. When
##'   several indicators are requested, a `patchwork` object combining one
##'   `ggplot` per indicator.
##' @method plot RiskMap_extrapolation
##' @export
plot.RiskMap_extrapolation <- function(x, indicator = "extrapolation", ...) {
  choices <- c("extrapolation", "spatial_unexplained", "covariate_novelty",
               "covariate_leverage")
  matched <- if (is.character(indicator) && length(indicator)) {
    choices[pmatch(indicator, choices, duplicates.ok = TRUE)]
  } else {
    NA_character_
  }
  if (anyNA(matched)) {
    abort_input("{.arg indicator} must be among {.val {choices}}.", arg = "indicator")
  }
  indicator <- matched
  if (length(indicator) > 1L) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      riskmap_abort("Plotting several maps needs the {.pkg patchwork} package.")
    }
    maps <- lapply(indicator, function(one) plot(x, indicator = one))
    return(patchwork::wrap_plots(maps, ncol = min(2L, length(maps))) &
             theme(plot.margin = margin(4, 14, 4, 14)))
  }
  data <- x
  class(data) <- setdiff(class(data), "RiskMap_extrapolation")
  data$.prediction_value <- data[[indicator]]
  attr(data, "grid_metadata") <- attr(x, "grid_metadata", exact = TRUE)
  # The class legend lists all four classes, including absent ones.
  map <- add_prediction_layer(ggplot(), list(boundaries = NULL),
                              list(data = data, polygon = FALSE),
                              facet_label = "", facet_levels = "",
                              show_legend = if (indicator == "extrapolation")
                                TRUE else NA)

  if (indicator == "extrapolation") {
    class_colours <- c("Supported" = "#EFEBE4", "Spatial" = "#8C7BBF",
                       "Covariate" = "#E3A34F",
                       "Spatial and covariate" = "#B2182B")
    # Cells with missing covariates are drawn in grey but get no legend key.
    fill_scale <- scale_fill_manual(
      values = class_colours, breaks = names(class_colours),
      drop = FALSE, na.value = "grey90",
      guide = guide_legend(title.position = "top", nrow = 2L)
    )
    title <- "Extrapolation"
  } else {
    values <- data$.prediction_value
    limits <- switch(
      indicator,
      spatial_unexplained = c(0, 1),
      covariate_novelty = c(0, max(pretty(c(0, values), n = 4))),
      covariate_leverage = range(pretty(c(0, values), n = 4))
    )
    fill_scale <- prediction_fill_scale(
      colours = prediction_role_palette("uncertainty"), limits = limits,
      values = values, breaks = NULL, role = "uncertainty",
      probability = indicator != "covariate_leverage"
    )
    title <- switch(
      indicator,
      spatial_unexplained = "Spatial variance unexplained by the data",
      covariate_novelty = "Distance outside covariate range",
      covariate_leverage = "Leverage relative to the data"
    )
  }
  map + fill_scale + labs(fill = title) + theme_void() +
    theme(legend.position = "top", legend.justification = "center",
          legend.title = element_text(hjust = 0.5, size = 9),
          legend.text = element_text(size = 7.5),
          strip.text = element_blank())
}
