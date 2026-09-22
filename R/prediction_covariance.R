#' Stable Gaussian prediction weights
#'
#' For observation incidence matrix H and latent covariance G, the observation
#' covariance is V = H G H' + omega2 I. Every cross-covariance used for latent
#' prediction has its rows in the column space of H. If Q is an orthonormal
#' basis for that space, the required product is therefore
#' C V^-1 = (C Q) (Q' V Q)^-1 Q'.
#'
#' This avoids the subtraction of terms of order 1 / omega2 in the Woodbury
#' precision formula. It also avoids constructing an observation-sized inverse.
#' When group membership is constant within each location, scaled location
#' indicators provide Q directly. Otherwise a thin, rank-revealing QR is used.
#'
#' @param incidence Observation-by-latent-effect incidence matrix H.
#' @param covariance Latent covariance G.
#' @param ids Matrix containing location and random-effect group indices.
#' @param variance Measurement-error variance omega2.
#' @return A function mapping a latent-to-observation cross-covariance C to
#'   prediction weights C V^-1. This is not a general-purpose inverse of V.
#' @importFrom Matrix sparseMatrix
#' @noRd
gaussian_prediction_weights <- function(incidence, covariance, ids, variance) {
  location <- match(ids[, 1], unique(ids[, 1]))
  first <- match(seq_len(max(location)), location)
  counts <- tabulate(location)
  nested <- all(ids == ids[first[location], , drop = FALSE])

  if (nested) {
    # The group averages are sufficient for predicting all latent effects.
    basis <- sparseMatrix(i = seq_along(location), j = location,
                          x = 1 / sqrt(counts[location]),
                          dims = c(length(location), length(counts)))
    reduced_design <- incidence[first, , drop = FALSE] * sqrt(counts)
  } else {
    decomposition <- qr(as.matrix(incidence))
    basis <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
    reduced_design <- crossprod(basis, as.matrix(incidence))
  }

  reduced_design <- as.matrix(reduced_design)
  reduced_covariance <- reduced_design %*% covariance %*% t(reduced_design)
  diag(reduced_covariance) <- diag(reduced_covariance) + variance
  root <- chol(reduced_covariance)
  basis_transpose <- Matrix::t(basis)

  function(cross_covariance) {
    projected <- as.matrix(cross_covariance %*% basis)
    weights <- backsolve(root, forwardsolve(t(root), t(projected)))
    as.matrix(t(weights) %*% basis_transpose)
  }
}
