##' @title Create Grid of Points Within Shapefile
##'
##' @description
##' Generates regularly spaced point centres within a polygon boundary.
##'
##' @param shp An object of class 'sf' containing POLYGONS or MULTIPOLYGONS within which the grid of points will be created.
##' @param spacing A single positive number specifying the distance between grid-point centres.
##' @param distance_units Character string, either `"km"` or `"m"`, giving the units of `spacing`. Defaults to `"km"`.
##'
##' @details
##' This function creates point centres within the boundaries of `shp`; it does not create polygon cells or define an areal prediction target. The CRS is inherited from `shp`, which must use a projected CRS with recognised linear units.
##'
##' @return
##' An 'sf' object containing the generated grid points within the shapefile.
##'
##' @export
##'
##' @examples
##' library(sf)
##'
##' # Example shapefile data
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
##' \code{\link[sf]{st_make_grid}}, \code{\link[sf]{st_intersects}}, \code{\link[sf]{st_transform}}, \code{\link[sf]{st_crs}}
##'
##'
create_grid <- function(shp,
                        spacing,
                        distance_units = c("km", "m")) {

  check_data(shp, "polygon")
  check_positive_number(spacing, "")
  stopifnot("'distance_units' must be either 'km' or 'm'" =
              is.character(distance_units) &&
              all(distance_units %in% c("km", "m")))
  distance_units <- match.arg(distance_units)

  if (st_is_longlat(shp)) {
    stop(
      "The coordinates of 'shp' are in longitude and latitude; ",
      "transform 'shp' to a suitable projected CRS before creating the grid.",
      call. = FALSE
    )
  }

  cellsize <- spacing / crs_to_distance_factor(shp, distance_units)
  grid_box <- st_sf(
    geometry = st_make_grid(shp,
                            cellsize = cellsize,
                            what = "centers")
  )

  study_boundary <- st_union(st_geometry(shp))
  grid_out <- st_filter(grid_box, study_boundary)

  if (nrow(grid_out) == 0) {
    stop(
      "No grid-point centres fall within 'shp'; try decreasing 'spacing'.",
      call. = FALSE
    )
  }

  grid_out
}
