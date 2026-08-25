##' @title Create Grid of Points Within Shapefile
##'
##' @description
##' Generates a grid of points within a given shapefile. The grid points are created based on a specified spatial resolution.
##'
##' @param shp An object of class 'sf' containing POLYGONS or MULTIPOLYGONS within which the grid of points will be created.
##' @param spat_res Numeric value specifying the spatial resolution in kilometers for the grid.
##' @param grid_crs Coordinate reference system for the grid. If NULL, the CRS of 'shp' is used. The shapefile 'shp' will be transformed to this CRS if specified.
##'
##' @details
##' This function creates a grid of points within the boundaries of the provided shapefile ('shp'). The grid points are generated using the specified spatial resolution ('spat_res'). If a coordinate reference system ('grid_crs') is provided, the shapefile is transformed to this CRS before creating the grid.
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
##' # Create grid with 10 km spatial resolution
##' grid <- create_grid(nc, spat_res = 10)
##'
##' # Plot the grid
##' plot(st_geometry(nc))
##' plot(grid, add = TRUE, col = 'red')
##'
##' @seealso
##' \code{\link[sf]{st_make_grid}}, \code{\link[sf]{st_intersects}}, \code{\link[sf]{st_transform}}, \code{\link[sf]{st_crs}}
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterre@@lancaster.ac.uk}
##'
create_grid <- function(shp,
                        spat_res,
                        grid_crs = NULL) {

  check_data(shp, "polygon")
  check_positive_integer(spat_res, "spat_res")

  if(is.null(grid_crs)) {
    grid_crs <- st_crs(shp)
    if (st_is_longlat(shp)) stop("The coordinates of 'shp' are in longitude and latitude - please set 'grid_crs'")
  } else {
    check_positive_integer(grid_crs, "grid_crs")
    crs <- tryCatch(
            st_crs(grid_crs),
            warning = function(w) {
              stop("The 'grid_crs' provided is not a valid CRS")
            }
          )
    shp <- st_transform(shp, crs = crs)
  }

  grid_box <- st_make_grid(shp,
                           cellsize = spat_res*1000,
                           what="centers")

  intersect <- st_intersects(grid_box,
                                 shp,
                                 sparse = FALSE)

  filter <- rowSums(intersect) > 0

  grid_out <- grid_box[filter]

  if (length(grid_out) == 0){
    stop("No points intersect with the 'shp' - try decreasing the 'spat_res'")
  }

  return(grid_out)
}
