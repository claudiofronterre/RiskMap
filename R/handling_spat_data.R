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
##' @param distance_crs `NULL` to keep the CRS of projected `boundaries`, or
##'   to project `boundaries` in longitude and latitude to the CRS from
##'   [propose_crs()] (with a message naming the CRS and the units of
##'   distance). Alternatively, a projected CRS to project `boundaries` to,
##'   as in [variogram()] and [summarise_distance()]. A grid built from
##'   `raster` takes the CRS of the raster, so `distance_crs` must then be
##'   `NULL` or that CRS.
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
##' The CRS of the grid is `distance_crs` when it is given, and otherwise that
##' of `boundaries`, or of `raster` when it is supplied. If `boundaries` is in
##' longitude and latitude and neither `distance_crs` nor a raster is given,
##' it is projected to the CRS from [propose_crs()] and a message names the
##' CRS and the units of distance. A warning of class
##' `riskmap_warning_crs` is raised when the CRS used distorts distances by
##' more than 1% over `boundaries`, or over the grid built from `raster`.
##'
##' The grid keeps its cells when rows or columns are selected with `[` or
##' with dplyr verbs such as `select()` and `filter()`, unless only the
##' geometry column is kept. Functions that build a new `sf` object, such as
##' `cbind()`, `merge()` and [sf::st_as_sf()], drop the class `RiskMap_grid`,
##' and [predict.RiskMap()] then assigns each point to the polygon that
##' contains it. Add columns with `$<-`, `dplyr::mutate()` or a dplyr join
##' instead. `rbind()` keeps the cells of its first grid, so use it only to
##' join parts of one grid.
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
                        distance_crs = NULL,
                        distance_units = c("km", "m"),
                        keep = c("centre", "overlap"),
                        raster = NULL) {

  distance_units <- match_choice(distance_units)
  keep <- match_choice(keep)
  if (!is.null(raster)) {
    return(create_grid_from_raster(raster, boundaries, spacing,
                                   distance_crs, distance_units, keep))
  }
  if (is.null(boundaries)) {
    abort_input("{.arg boundaries}, {.arg raster} or both must be given.", arg = "boundaries")
  }
  check_data(boundaries, "polygon")
  check_positive_number(spacing, "")

  input_crs <- st_crs(boundaries)
  boundaries <- project_for_distances(boundaries, distance_crs,
                                      "distance_crs", "to create the grid",
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

#' Keep the cells of a grid when its rows or columns are selected
#'
#' `[.sf` comes first, since "sf" comes before "RiskMap_grid". It removes the
#' class "sf", selects with `[`, which reaches this method, and then adds the
#' geometry back with `st_sf()`, which keeps the class and attributes of the
#' data frame it is given. `[.data.frame` drops the metadata when columns are
#' selected, so it is restored here (#214). sf's dplyr verbs select with `[`
#' too, so they need no method of their own. With `drop = TRUE`, `[.sf` adds
#' no geometry back, so nothing is restored.
#' @noRd
#' @export
`[.RiskMap_grid` <- function(x, ..., drop = FALSE) {
  metadata <- attr(x, "grid_metadata", exact = TRUE)
  out <- NextMethod()
  if (drop || !is.data.frame(out) || is.null(metadata)) return(out)
  attr(out, "grid_metadata") <- metadata
  as_riskmap_grid(out)
}

#' Say that points on a regular grid are assigned to polygons one by one
#'
#' Called by `predict()` when it aggregates over polygons. Points that are
#' not a grid from `create_grid()` with its cells but lie on a square lattice
#' are probably such a grid that lost its cells, for example through
#' `cbind()`, or one made with `sf::st_make_grid()`, so the message says how
#' to have each cell weighted by its coverage instead (#214).
#' @param points The prediction locations as the user gave them, in their own
#'   CRS, since a lattice is no longer one once projected to another.
#' @return Whether the message was given, invisibly.
#' @noRd
inform_point_assignment <- function(points) {
  if (inherits(points, "RiskMap_grid") &&
      !is.null(attr(points, "grid_metadata", exact = TRUE))) {
    return(invisible(FALSE))
  }
  if (!on_square_lattice(points)) return(invisible(FALSE))
  riskmap_inform(c(
    "{.arg newdata} lie on a regular grid that doesn't come from {.fn create_grid}, so each point counts fully in the polygon that contains it.",
    "i" = "Make the grid with {.fn create_grid} to weight each cell by the fraction of it within each polygon. Add columns to it with {.code $<-} or {.fn dplyr::mutate}, since {.fn cbind} and {.fn merge} drop its cells."
  ))
  invisible(TRUE)
}

#' Whether points lie on a square lattice
#'
#' The points are distinct, at least four, in at least two rows and two
#' columns, and their coordinates are whole multiples of one spacing from the
#' first row and column, as the centres of the cells of a square grid, of
#' which some may be missing. At least half of the gaps between neighbouring
#' rows, and between neighbouring columns, are one spacing, as in a grid.
#' Scattered locations whose coordinates are whole metres are on a lattice of
#' 1 m, but their rows and columns are mostly further apart.
#' @param points An `sf` object.
#' @noRd
on_square_lattice <- function(points) {
  geometry <- st_geometry(points)
  if (length(geometry) < 4L ||
      !all(st_geometry_type(geometry) == "POINT")) {
    return(FALSE)
  }
  coordinates <- unique(st_coordinates(geometry)[, 1:2, drop = FALSE])
  if (nrow(coordinates) < 4L || !all(is.finite(coordinates))) return(FALSE)
  steps <- apply(coordinates, 2L, function(values) {
    values <- sort(unique(values))
    gaps <- diff(values)
    gaps <- gaps[gaps > 1e-8 * max(abs(values), 1)]
    if (!length(gaps)) return(NA_real_)
    step <- min(gaps)
    positions <- (values - values[1L]) / step
    on_lattice <- all(abs(positions - round(positions)) <= 1e-6)
    adjacent <- mean(abs(gaps / step - 1) <= 1e-6) >= 0.5
    if (on_lattice && adjacent) step else NA_real_
  })
  all(is.finite(steps)) && abs(steps[1L] - steps[2L]) <= 1e-6 * max(steps)
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
                                    distance_crs, distance_units, keep) {
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
  if (!is.null(distance_crs)) {
    check_crs(distance_crs, arg = "distance_crs")
    if (st_crs(distance_crs) != grid_crs) {
      abort_input(c("{.arg distance_crs} must be the coordinate reference system of {.arg raster}, or {.code NULL}.", "i" = "A grid built from {.arg raster} has the cells of the raster, so it takes its CRS. Project the raster first, for example with {.fn terra::project}, to build the grid in another CRS."), arg = "distance_crs")
    }
  }
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

##' @title Convex Hull of an sf Object
##'
##' @description Computes the convex hull of an `sf` object, returning the boundaries of the smallest polygon that can enclose all geometries in the input.
##'
##' @param data An `sf` data frame object containing geometries.
##'
##' @return An `sf` object representing the convex hull of the input geometries.
##'
##' @details The convex hull is the smallest convex polygon that encloses all points in the input `sf` object. This function computes the convex hull by first uniting all geometries in the input using `st_union()`, and then applying `st_convex_hull()` to obtain the polygonal boundary. The result is returned as an `sf` object containing the convex hull geometry.
##'
##' @seealso \code{\link[sf]{st_convex_hull}}, \code{\link[sf]{st_union}}
##'
##' @examples
##' library(sf)
##'
##' # Create example sf object
##' points <- st_sfc(st_point(c(0,0)), st_point(c(1,1)), st_point(c(2,2)), st_point(c(0,2)))
##' sf_points <- st_sf(geometry = points)
##'
##' # Calculate the convex hull
##' convex_hull_result <- create_convex_hull(sf_points)
##'
##' # Plot the result
##' plot(sf_points, col = 'blue', pch = 19)
##' plot(convex_hull_result, add = TRUE, border = 'red')
##' @export
create_convex_hull <- function(data) {
  check_supplied(data)
  # Check if the input is an sf object
  if (!inherits(data, "sf")) {
    abort_input("{.arg data} must be an {.cls sf} object, not {.obj_type_friendly {data}}.", arg = "data")
  }

  # Get the geometry from the sf object
  geometry <- st_geometry(data)

  # Compute the convex hull
  convex_hull <- st_convex_hull(st_union(geometry))

  # Return the convex hull as an sf object
  return(st_sf(geometry = convex_hull))
}
