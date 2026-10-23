##' @title Create Grid of Points Within Boundaries
##'
##' @description
##' Generates regularly spaced point centres within a polygon boundary and
##' retains the information needed to reconstruct their square grid cells.
##'
##' @param boundaries An object of class `sf` containing POLYGONS or MULTIPOLYGONS within which the grid of points will be created.
##' @param spacing A single positive number specifying the distance between grid-point centres.
##' @param distance_units Character string, either `"km"` or `"m"`, giving the units of `spacing`. Defaults to `"km"`.
##' @param keep Which cells to keep. `"centre"` (the default) keeps cells
##'   whose centre lies within `boundaries`. `"overlap"` keeps every cell that
##'   overlaps `boundaries`, including cells along the edge whose centre lies
##'   outside. Use `"overlap"` when totals over the study area matter, for
##'   example population counts along a coastline.
##'
##' @details
##' This function returns point centres rather than polygon cells. The returned
##' object has class `RiskMap_grid` and carries lightweight metadata describing
##' the square cells. This allows other RiskMap functions to reconstruct their
##' footprints when exact cell coverage is needed without storing a second copy
##' of the grid as polygons. The metadata also records the study boundary, so
##' that maps can clip edge cells to it.
##'
##' With `keep = "overlap"`, cells along the edge extend beyond `boundaries`.
##' Areal predictions from [predict.RiskMap()] weight each cell by the fraction
##' lying within each polygon, so these cells contribute only their covered
##' part.
##'
##' The CRS is inherited from `boundaries`. If `boundaries` is in
##' longitude/latitude, it is automatically reprojected to an appropriate UTM
##' zone (see [propose_utm()]) and a message reports the conversion used;
##' transform `boundaries` to a projected CRS yourself first to use a different
##' one.
##'
##' @return
##' A `RiskMap_grid` and `sf` object containing the generated grid points. The
##' `grid_metadata` attribute records the spacing, units, cell size, cell area,
##' cell shape, grid CRS, alignment origin, the cells kept and the study
##' boundary. The `grid_id` column provides a
##' stable identifier for each retained cell centre.
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
##' # Also keep the edge cells whose centre lies outside the boundary
##' grid_edges <- create_grid(nc, spacing = 10, keep = "overlap")
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
                        distance_units = c("km", "m"),
                        keep = c("centre", "overlap")) {

  check_data(boundaries, "polygon")
  check_positive_number(spacing, "")
  stopifnot("'distance_units' must be either 'km' or 'm'" =
              is.character(distance_units) &&
              all(distance_units %in% c("km", "m")))
  distance_units <- match.arg(distance_units)
  keep <- match.arg(keep)

  input_crs <- st_crs(boundaries)
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
  grid_out <- if (keep == "centre") {
    st_filter(grid_box, study_boundary)
  } else {
    cells <- st_buffer(st_geometry(grid_box), dist = cellsize / 2,
                       endCapStyle = "SQUARE")
    grid_box[lengths(st_intersects(cells, study_boundary)) > 0, ]
  }

  if (nrow(grid_out) == 0) {
    stop(
      "No grid cells fall within 'boundaries'; try decreasing 'spacing' and check 'distance_units'.",
      call. = FALSE
    )
  }

  grid_out$grid_id <- seq_len(nrow(grid_out))
  attr(grid_out, "grid_metadata") <- list(
    spacing = spacing,
    distance_units = distance_units,
    cellsize = cellsize,
    cell_area = cellsize^2,
    cell_shape = "square",
    grid_crs = st_crs(boundaries),
    input_crs = input_crs,
    origin = unname(c(st_bbox(boundaries)["xmin"],
                      st_bbox(boundaries)["ymin"])),
    keep = keep,
    boundary = study_boundary
  )
  class(grid_out) <- c("RiskMap_grid", class(grid_out))
  grid_out
}
