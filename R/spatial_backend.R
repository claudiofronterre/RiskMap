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
