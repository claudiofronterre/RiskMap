##' @title Create Grid of Points Within Boundaries
##'
##' @description
##' Generates regularly spaced point centres within a polygon boundary, or
##' takes them from the cells of a raster, and retains the information needed
##' to reconstruct their square grid cells.
##'
##' @param boundaries An object of class `sf` containing POLYGONS or
##'   MULTIPOLYGONS within which the grid of points will be created. Optional
##'   when `raster` is supplied.
##' @param spacing A single positive number specifying the distance between
##'   grid-point centres. Taken from the resolution of `raster` when that is
##'   supplied.
##' @param distance_units Character string, either `"km"` or `"m"`, giving the units of `spacing`. Defaults to `"km"`.
##' @param keep Which cells to keep. `"centre"` (the default) keeps cells
##'   whose centre lies within `boundaries`. `"overlap"` keeps every cell that
##'   overlaps `boundaries`, including cells along the edge whose centre lies
##'   outside. Use `"overlap"` when totals over the study area matter, for
##'   example population counts along a coastline.
##' @param raster Optional `SpatRaster` from the terra package, in a projected
##'   coordinate reference system with square cells. When supplied, the grid
##'   consists of the raster's cells: their centres are the grid points and
##'   the raster's layers become columns of the grid.
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
##' A grid built from `raster` is aligned with that raster by construction, so
##' its values need no extraction. This suits population counts: aggregate the
##' counts to the resolution of the grid with [terra::aggregate()] and
##' `fun = "sum"`, after projecting them with [terra::project()] and
##' `method = "sum"` if needed, and use the result as `raster`. With
##' `boundaries`, cells are kept according to `keep`; without, the cells with
##' at least one non-missing value are kept. Cells kept inside `boundaries`
##' may have missing values, for example over water.
##'
##' The CRS is inherited from `boundaries`, or from `raster` when it is
##' supplied. If `boundaries` is in longitude/latitude and no raster is given,
##' it is automatically reprojected to the CRS from [propose_crs()] and a
##' message names the CRS and the units of distance;
##' transform `boundaries` to a projected CRS yourself first to use a
##' different one. A warning of class `riskmap_warning_crs` is raised when the
##' CRS used distorts distances by more than 1% over `boundaries`, or over the
##' grid built from `raster`.
##'
##' @return
##' A `RiskMap_grid` and `sf` object containing the generated grid points. The
##' `grid_metadata` attribute records the spacing, units, cell size, cell area,
##' cell shape, grid CRS, alignment origin, the cells kept and the study
##' boundary. The `grid_id` column provides a stable identifier for each
##' retained cell centre. A grid built from `raster` also has one column for
##' each raster layer.
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
##' # Use the cells of a 20 km raster, with its values as a column
##' cells <- terra::align(terra::ext(nc), 20000)
##' counts <- terra::rast(cells, resolution = 20000, crs = "EPSG:32617")
##' terra::values(counts) <- seq_len(terra::ncell(counts))
##' names(counts) <- "count"
##' grid_raster <- create_grid(nc, raster = counts)
##'
##' @seealso
##' \code{\link[sf]{st_make_grid}}, \code{\link[sf]{st_intersection}}
##'
##'
create_grid <- function(boundaries = NULL,
                        spacing = NULL,
                        distance_units = c("km", "m"),
                        keep = c("centre", "overlap"),
                        raster = NULL) {

  distance_units <- match_choice(distance_units)
  keep <- match_choice(keep)
  if (!is.null(raster)) {
    return(create_grid_from_raster(raster, boundaries, spacing,
                                   distance_units, keep))
  }
  if (is.null(boundaries)) {
    abort_input("{.arg boundaries}, {.arg raster} or both must be given.", arg = "boundaries")
  }
  check_data(boundaries, "polygon")
  check_positive_number(spacing, "")

  input_crs <- st_crs(boundaries)
  boundaries <- project_for_distances(boundaries, NULL, NULL,
                                      "to create the grid",
                                      data_arg = "boundaries",
                                      units = distance_units)

  cellsize <- spacing / crs_to_distance_factor(boundaries, distance_units,
                                               "boundaries")
  grid_box <- st_sf(
    geometry = st_make_grid(boundaries,
                            cellsize = cellsize,
                            what = "centers")
  )

  study_boundary <- st_union(st_geometry(boundaries))
  grid_out <- keep_grid_cells(grid_box, study_boundary, cellsize, keep)

  if (nrow(grid_out) == 0) {
    abort_input(c("No grid cells fall within {.arg boundaries}.", "i" = "Try a smaller {.arg spacing}, and check {.arg distance_units}."), arg = "spacing")
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
  as_riskmap_grid(grid_out)
}

#' Add the RiskMap_grid class to an sf grid
#'
#' vctrs, and therefore dplyr, reject a data frame whose subclasses do not all
#' come before "data.frame", so the class goes just before it, after "sf".
#' @noRd
as_riskmap_grid <- function(grid) {
  classes <- setdiff(class(grid), "RiskMap_grid")
  class(grid) <- append(classes, "RiskMap_grid",
                        after = match("data.frame", classes) - 1L)
  grid
}

#' Keep the grid cells that lie in, or overlap, a study boundary
#' @noRd
keep_grid_cells <- function(grid_box, study_boundary, cellsize, keep) {
  if (keep == "centre") {
    return(st_filter(grid_box, study_boundary))
  }
  cells <- st_buffer(st_geometry(grid_box), dist = cellsize / 2,
                     endCapStyle = "SQUARE")
  grid_box[lengths(st_intersects(cells, study_boundary)) > 0, ]
}

#' Grid points at the centres of the cells of a raster
#' @noRd
create_grid_from_raster <- function(raster, boundaries, spacing,
                                    distance_units, keep) {
  if (!inherits(raster, "SpatRaster")) {
    abort_input("{.arg raster} must be a {.cls SpatRaster} from the {.pkg terra} package, not {.obj_type_friendly {raster}}.", arg = "raster")
  }
  if (terra::is.lonlat(raster)) {
    abort_input(c("{.arg raster} must be in a projected coordinate reference system.", "i" = "Project it first, for example with {.fn terra::project}, using {.code method = \"sum\"} for counts."), arg = "raster")
  }
  resolution <- terra::res(raster)
  if (abs(resolution[1] - resolution[2]) > 1e-8 * max(resolution)) {
    abort_input("{.arg raster} must have square cells.", arg = "raster")
  }
  cellsize <- resolution[1]
  grid_crs <- st_crs(terra::crs(raster))
  centres <- terra::xyFromCell(raster, seq_len(terra::ncell(raster)))
  grid_box <- st_as_sf(
    cbind(as.data.frame(centres), terra::values(raster, dataframe = TRUE)),
    coords = c("x", "y"), crs = grid_crs
  )
  factor <- crs_to_distance_factor(grid_box, distance_units, "raster")
  raster_spacing <- cellsize * factor
  if (!is.null(spacing)) {
    check_positive_number(spacing, "")
    if (abs(spacing - raster_spacing) > 1e-6 * raster_spacing) {
      abort_input("{.arg spacing} ({spacing} {distance_units}) must match the resolution of {.arg raster} ({raster_spacing} {distance_units}).", arg = "spacing")
    }
  }

  if (is.null(boundaries)) {
    values <- st_drop_geometry(grid_box)
    grid_out <- grid_box[rowSums(!is.na(values)) > 0, ]
    study_boundary <- NULL
    keep <- "raster"
  } else {
    check_data(boundaries, "polygon")
    boundaries <- st_transform(boundaries, grid_crs)
    study_boundary <- st_union(st_geometry(boundaries))
    grid_out <- keep_grid_cells(grid_box, study_boundary, cellsize, keep)
  }
  if (nrow(grid_out) == 0) {
    abort_input("No cells of {.arg raster} fall within {.arg boundaries} or have values.", arg = "raster")
  }

  check_crs_distortion(grid_out, NULL, "raster")

  grid_out$grid_id <- seq_len(nrow(grid_out))
  attr(grid_out, "grid_metadata") <- list(
    spacing = raster_spacing,
    distance_units = distance_units,
    cellsize = cellsize,
    cell_area = cellsize^2,
    cell_shape = "square",
    grid_crs = grid_crs,
    input_crs = grid_crs,
    origin = unname(as.vector(terra::ext(raster))[c(1, 3)]),
    keep = keep,
    boundary = study_boundary
  )
  as_riskmap_grid(grid_out)
}
