##' @title Create Grid of Points Within Boundaries
##'
##' @description
##' Generates regularly spaced point centres within a polygon boundary.
##'
##' @param boundaries An object of class `sf` containing POLYGONS or MULTIPOLYGONS within which the grid of points will be created.
##' @param spacing A single positive number specifying the distance between grid-point centres.
##' @param distance_units Character string, either `"km"` or `"m"`, giving the units of `spacing`. Defaults to `"km"`.
##'
##' @details
##' This function creates point centres within `boundaries`; it does not create polygon cells or define an areal prediction target. The CRS is inherited from `boundaries`. If `boundaries` is in longitude/latitude, it is automatically reprojected to an appropriate UTM zone (see [propose_utm()]) and a message reports the conversion used; transform `boundaries` to a projected CRS yourself first to use a different one.
##'
##' @return
##' An `sf` object containing the generated grid points within the boundaries.
##'
##' @export
##'
##' @examples
##' library(sf)
##'
##' # Example boundary data
##' nc <- st_read(system.file("shape/nc.shp", package="sf"))
##' nc <- st_transform(nc, crs = 32617)
##'
##' # Create grid with 10 km spacing
##' grid <- create_grid(nc, spacing = 10)
##'
##' # Plot the grid
##' plot(st_geometry(nc))
##' plot(grid, add = TRUE, col = 'red')
##'
##' @seealso
##' \code{\link[sf]{st_make_grid}}, \code{\link[sf]{st_intersection}}
##'
##'
create_grid <- function(boundaries,
                        spacing,
                        distance_units = c("km", "m")) {

  check_data(boundaries, "polygon")
  check_positive_number(spacing, "")
  stopifnot("'distance_units' must be either 'km' or 'm'" =
              is.character(distance_units) &&
              all(distance_units %in% c("km", "m")))
  distance_units <- match.arg(distance_units)

  if (st_is_longlat(boundaries)) {
    auto_crs <- propose_utm(boundaries)
    boundaries <- st_transform(boundaries, crs = auto_crs)
    message("'boundaries' is in longitude/latitude; automatically reprojecting to EPSG:",
            auto_crs, " to create the grid. Transform 'boundaries' to a projected CRS ",
            "yourself to override.")
  }

  cellsize <- spacing / crs_to_distance_factor(boundaries, distance_units)
  grid_box <- st_sf(
    geometry = st_make_grid(boundaries,
                            cellsize = cellsize,
                            what = "centers")
  )

  study_boundary <- st_union(st_geometry(boundaries))
  grid_out <- st_filter(grid_box, study_boundary)

  if (nrow(grid_out) == 0) {
    stop(
      "No grid-point centres fall within 'boundaries'; try decreasing 'spacing' and check 'distance_units'.",
      call. = FALSE
    )
  }

  grid_out
}
