##' Assess where predictions extrapolate beyond the data
##'
##' Computes, for each prediction location, indicators of whether a prediction
##' from a fitted geostatistical model relies on extrapolation, either in space
##' (far from the data relative to the estimated spatial correlation) or in the
##' covariates (outside the conditions observed in the data).
##'
##' @param object A fitted `RiskMap` model.
##' @param newdata An `sf` object of prediction points with the covariates used
##'   by the model, typically a grid from [create_grid()].
##' @param spatial_threshold Share of the spatial variance left unexplained by
##'   the data above which a location is classed as spatial extrapolation.
##'   Defaults to 0.8.
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
##' The unexplained spatial variance is
##' \eqn{1 - \rho(x)^\top (R + \nu^2 I)^{-1} \rho(x)}, where \eqn{\rho(x)}
##' holds the estimated correlations between \eqn{x} and the sampled
##' locations, \eqn{R} their correlation matrix, and \eqn{\nu^2} the ratio of
##' the nugget (and, for Gaussian models, the measurement-error variance) to the
##' spatial variance. It is close to 0 next to well-sampled locations and close
##' to 1 where the data carry no information about the spatial effect. Unlike
##' the distance to the nearest sampled location, it accounts for the estimated
##' range of correlation and for clustering of the sampled locations. It depends
##' only on where the data were collected, not on the outcomes or on the number
##' of individuals at each location. For binomial and Poisson outcomes, which
##' carry less information than a noise-free observation of the linear
##' predictor, it is therefore an optimistic (lower) value.
##'
##' A location is classed as covariate extrapolation when any covariate lies
##' outside its observed range or its relative leverage exceeds 1, and as
##' spatial extrapolation when its unexplained spatial variance exceeds
##' `spatial_threshold`. The two are kept separate because they call for
##' different responses: more sampling where data are sparse, and caution about
##' the form of the covariate effects outside the observed conditions.
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
##'     \item{spatial_unexplained}{Share of the spatial variance not explained
##'       by the data.}
##'     \item{extrapolation}{Factor with levels `"Supported"`, `"Spatial"`,
##'       `"Covariate"` and `"Spatial and covariate"`.}
##'   }
##'   Rows with missing covariates have missing indicators. Map the result with
##'   [plot.RiskMap_extrapolation()].
##' @seealso [plot.RiskMap_extrapolation()], [predict.RiskMap()]
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:50, ]
##' fit <- glgpm(y ~ log(pop_dens) + gp(), data = example_data,
##'              family = "gaussian", messages = FALSE)
##' extrapolation <- assess_extrapolation(fit, italy_sim[51:100, ])
##' table(extrapolation$extrapolation)
##' @export
assess_extrapolation <- function(object, newdata, spatial_threshold = 0.8) {
  if (!inherits(object, "RiskMap")) {
    stop("'object' must be a fitted RiskMap model", call. = FALSE)
  }
  check_range(spatial_threshold, min = 0, max = 1,
              name = "spatial_threshold")
  prediction <- validate_prediction_data(object, newdata, include = character())
  model_data <- prediction$model_data
  valid_rows <- prediction$valid_rows

  covariates <- extrapolation_covariates(
    object, prediction_design_matrix(object, st_drop_geometry(model_data))
  )
  spatial_unexplained <- extrapolation_spatial_unexplained(object, model_data)
  covariate_flag <- covariates$outside > 0 | covariates$leverage > 1 + 1e-8
  spatial_flag <- spatial_unexplained > spatial_threshold
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
  attr(out, "covariate_ranges") <- covariates$ranges
  class(out) <- c("RiskMap_extrapolation", class(out))
  out
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

#' Share of the prior spatial variance not explained by the sampled locations
#' @noRd
extrapolation_spatial_unexplained <- function(object, model_data) {
  parameters <- coef(object)
  ratio <- 0
  if (!is.null(parameters$tau2)) ratio <- ratio + parameters$tau2
  if (!is.null(parameters$sigma2_me)) ratio <- ratio + parameters$sigma2_me
  ratio <- ratio / parameters$sigma2

  observed <- unique(as.matrix(object$coords))
  correlation <- matern_correlation(as.matrix(stats::dist(observed)),
                                    phi = parameters$phi, kappa = object$kappa)
  diag(correlation) <- 1 + ratio
  factor <- chol(correlation)

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
  pmin(pmax(unexplained, 0), 1)
}

##' Map indicators of extrapolation
##'
##' Maps one of the indicators returned by [assess_extrapolation()], in the
##' style of [plot.RiskMap_prediction()].
##'
##' @param x A `RiskMap_extrapolation` object.
##' @param indicator Indicator to map: `"extrapolation"` (the default, the
##'   four classes), `"spatial_unexplained"`, `"covariate_novelty"` or
##'   `"covariate_leverage"`.
##' @param ... Not used.
##' @return A `ggplot` object, to which further layers can be added.
##' @method plot RiskMap_extrapolation
##' @export
plot.RiskMap_extrapolation <- function(
    x,
    indicator = c("extrapolation", "spatial_unexplained", "covariate_novelty",
                  "covariate_leverage"),
    ...) {
  indicator <- match.arg(indicator)
  data <- x
  class(data) <- setdiff(class(data), "RiskMap_extrapolation")
  data$.prediction_value <- data[[indicator]]
  attr(data, "grid_metadata") <- attr(x, "grid_metadata", exact = TRUE)
  map <- add_prediction_layer(ggplot(), list(boundaries = NULL),
                              list(data = data, polygon = FALSE),
                              facet_label = "", facet_levels = "")

  if (indicator == "extrapolation") {
    fill_scale <- scale_fill_manual(
      values = c("Supported" = "#EFEBE4", "Spatial" = "#8C7BBF",
                 "Covariate" = "#E3A34F", "Spatial and covariate" = "#B2182B"),
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
