#' Fast pairwise distances
#'
#' Compute Euclidean distances using the same ordering and attributes as
#' [stats::dist()].
#'
#' @param coordinates Numeric matrix with locations in rows.
#' @return A `dist` object.
#' @noRd
pairwise_distances <- function(coordinates) {
  coordinates <- as.matrix(coordinates)
  distances <- cpp_pairwise_distances(coordinates)

  structure(distances,
            Size = nrow(coordinates),
            Diag = FALSE,
            Upper = FALSE,
            method = "euclidean",
            call = match.call(),
            class = "dist")
}

#' Fast distances between two sets of locations
#'
#' @param first Numeric matrix with the first set of locations in rows.
#' @param second Numeric matrix with the second set of locations in rows.
#' @return A matrix with one row per location in `first` and one column per
#'   location in `second`.
#' @noRd
cross_distances <- function(first, second) {
  cpp_cross_distances(as.matrix(first), as.matrix(second))
}

#' Scale spatial coordinates and range
#'
#' The maximum observed pairwise distance is used for both fitting and future
#' prediction locations. Scaling the coordinates and range by the same value
#' leaves every ratio of distance to range unchanged.
#'
#' @param coordinates Numeric matrix containing observed locations.
#' @param phi Optional spatial range in the coordinate units.
#' @return A list containing scaled coordinates, scaled `phi` when supplied,
#'   and the distance scale.
#' @noRd
scale_spatial_coordinates <- function(coordinates, phi = NULL) {
  coordinates <- as.matrix(coordinates)
  if (nrow(coordinates) < 2) {
    stop("Spatial coordinates must contain at least two distinct locations.",
         call. = FALSE)
  }
  distance_scale <- max(pairwise_distances(coordinates))

  if (!is.finite(distance_scale) || distance_scale <= 0) {
    stop("Spatial coordinates must contain at least two distinct locations.",
         call. = FALSE)
  }

  list(coordinates = coordinates / distance_scale,
       phi = if (is.null(phi)) NULL else phi / distance_scale,
       distance_scale = distance_scale)
}

#' Restore a spatial range to the original coordinate units
#'
#' @param phi_scaled Spatial range expressed on the internally scaled geometry.
#' @param distance_scale Maximum observed pairwise distance used for scaling.
#' @return The spatial range in the original coordinate units.
#' @noRd
restore_spatial_range <- function(phi_scaled, distance_scale) {
  phi_scaled * distance_scale
}

#' Stable Cholesky factorisation
#'
#' Symmetrise a covariance matrix and factorise it. Scale-aware diagonal jitter
#' is attempted only after the unmodified matrix fails.
#'
#' @param covariance Numeric square covariance matrix.
#' @param context Short description used in diagnostics.
#' @param allow_jitter Whether to try reported scale-aware diagonal jitter after
#'   an unmodified factorisation fails.
#' @return An upper-triangular Cholesky factor. The applied jitter is stored in
#'   the `jitter` attribute.
#' @noRd
factor_covariance <- function(covariance, context = "covariance matrix",
                              allow_jitter = TRUE) {
  covariance <- as.matrix(covariance)
  if (nrow(covariance) != ncol(covariance)) {
    stop("The ", context, " must be square.", call. = FALSE)
  }

  covariance <- (covariance + t(covariance)) / 2
  root <- tryCatch(chol(covariance), error = function(error) NULL)
  if (!is.null(root)) {
    attr(root, "jitter") <- 0
    return(root)
  }

  if (!allow_jitter) {
    stop("The ", context, " is not positive definite.", call. = FALSE)
  }

  covariance_scale <- max(abs(diag(covariance)), 1)
  relative_jitter <- 10^seq(-12, -6)
  for (relative_value in relative_jitter) {
    jitter <- relative_value * covariance_scale
    adjusted <- covariance
    diag(adjusted) <- diag(adjusted) + jitter
    root <- tryCatch(chol(adjusted), error = function(error) NULL)
    if (!is.null(root)) {
      attr(root, "jitter") <- jitter
      warning("The ", context, " required diagonal jitter of ",
              format(jitter, scientific = TRUE), " for Cholesky factorisation.",
              call. = FALSE)
      return(root)
    }
  }

  stop("The ", context, " is not positive definite, even after adaptive jitter.",
       call. = FALSE)
}

#' Solve a positive-definite system from its Cholesky factor
#'
#' @param root Upper-triangular factor returned by [factor_covariance()].
#' @param right_hand_side Numeric vector or matrix.
#' @return Solution to `crossprod(root) %*% x = right_hand_side`.
#' @noRd
solve_from_cholesky <- function(root, right_hand_side) {
  backsolve(root, forwardsolve(t(root), right_hand_side))
}

#' Compute a positive-definite log determinant from its Cholesky factor
#'
#' @param root Upper-triangular factor returned by [factor_covariance()].
#' @return Log determinant of the original covariance matrix.
#' @noRd
log_determinant_from_cholesky <- function(root) {
  2 * sum(log(diag(root)))
}

#' Compute prediction weights without forming a covariance inverse
#'
#' @param cross_covariance Prediction-by-observation cross-covariance matrix.
#' @param root Upper-triangular Cholesky factor of the observation covariance.
#' @return `cross_covariance %*% solve(covariance)`.
#' @noRd
cholesky_prediction_weights <- function(cross_covariance, root) {
  t(solve_from_cholesky(root, t(cross_covariance)))
}

#' Validate conditional marginal variances
#'
#' @param marginal_variance Unconditional marginal variance.
#' @param weights Prediction weights.
#' @param cross_covariance Prediction-by-observation cross-covariance matrix.
#' @return Non-negative conditional variances.
#' @noRd
conditional_variances <- function(marginal_variance, weights,
                                  cross_covariance) {
  variance <- marginal_variance - rowSums(weights * cross_covariance)
  tolerance <- 100 * .Machine$double.eps *
    max(1, abs(marginal_variance), abs(variance))

  if (any(variance < -tolerance)) {
    stop("The conditional covariance produced materially negative variances.",
         call. = FALSE)
  }

  pmax(variance, 0)
}

#' Draw independent conditional Gaussian samples
#'
#' @param mean Numeric vector or matrix of conditional means.
#' @param standard_deviation Numeric vector of conditional standard deviations.
#' @param n_samples Number of samples.
#' @return Matrix with prediction locations in rows and samples in columns.
#' @noRd
sample_independent_gaussian <- function(mean, standard_deviation, n_samples) {
  n_prediction <- length(standard_deviation)
  noise <- matrix(rnorm(n_prediction * n_samples), nrow = n_prediction)
  mean + standard_deviation * noise
}

#' Draw correlated conditional Gaussian samples
#'
#' @param mean Numeric vector or matrix of conditional means.
#' @param lower_root Lower-triangular covariance factor.
#' @param n_samples Number of samples.
#' @return Matrix with prediction locations in rows and samples in columns.
#' @noRd
sample_correlated_gaussian <- function(mean, lower_root, n_samples) {
  n_prediction <- nrow(lower_root)
  noise <- matrix(rnorm(n_prediction * n_samples), nrow = n_prediction)
  mean + lower_root %*% noise
}

#' Select an internal marginal-prediction batch size
#'
#' Keep each dense prediction-by-observation intermediate near 64 MiB. Small
#' problems remain unbatched.
#'
#' @param n_prediction Number of prediction locations.
#' @param n_conditioning Number of conditioning locations or observations.
#' @param target_bytes Target size of one dense intermediate matrix.
#' @return Integer batch size.
#' @noRd
marginal_prediction_batch_size <- function(
    n_prediction, n_conditioning, target_bytes = 64 * 1024^2) {
  if (n_prediction == 0 || n_conditioning == 0) {
    return(n_prediction)
  }

  max(1L, min(n_prediction, floor(target_bytes / (8 * n_conditioning))))
}

#' Compute marginal spatial predictions in location batches
#'
#' @param prediction_coordinates Prediction coordinates.
#' @param conditioning_coordinates Coordinates defining cross-covariances.
#' @param weight_function Function mapping cross-covariances to weights.
#' @param conditional_signal Vector or matrix multiplied by the weights.
#' @param marginal_variance Unconditional spatial variance.
#' @param phi Matérn range on the internal coordinate scale.
#' @param kappa Matérn smoothness.
#' @param n_samples Number of predictive samples.
#' @param batch_size Number of prediction locations per batch.
#' @return Matrix with prediction locations in rows and samples in columns.
#' @noRd
batched_marginal_prediction <- function(
    prediction_coordinates, conditioning_coordinates, weight_function,
    conditional_signal, marginal_variance, phi, kappa, n_samples,
    batch_size) {
  n_prediction <- nrow(prediction_coordinates)
  samples <- matrix(rnorm(n_prediction * n_samples), nrow = n_prediction)

  for (start in seq.int(1L, n_prediction, by = batch_size)) {
    index <- start:min(start + batch_size - 1L, n_prediction)
    distance <- cross_distances(
      prediction_coordinates[index, , drop = FALSE],
      conditioning_coordinates
    )
    cross_covariance <- marginal_variance * matern_correlation(
      distance, phi = phi, kappa = kappa
    )
    weights <- weight_function(cross_covariance)
    conditional_mean <- weights %*% conditional_signal
    if (is.null(dim(conditional_signal))) {
      conditional_mean <- as.numeric(conditional_mean)
    }
    conditional_sd <- sqrt(conditional_variances(
      marginal_variance, weights, cross_covariance
    ))
    samples[index, ] <- conditional_mean +
      conditional_sd * samples[index, , drop = FALSE]
  }

  samples
}
