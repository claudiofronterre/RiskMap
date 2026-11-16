# Create Grid of Points Within Boundaries

Generates regularly spaced point centres within a polygon boundary, or
takes them from the cells of a raster, and retains the information
needed to reconstruct their square grid cells.

## Usage

``` r
create_grid(
  boundaries = NULL,
  spacing = NULL,
  distance_crs = NULL,
  distance_units = c("km", "m"),
  keep = c("centre", "overlap"),
  raster = NULL
)
```

## Arguments

- boundaries:

  An object of class `sf` containing POLYGONS or MULTIPOLYGONS within
  which the grid of points will be created. Optional when `raster` is
  supplied.

- spacing:

  A single positive number specifying the distance between grid-point
  centres. Taken from the resolution of `raster` when that is supplied.

- distance_crs:

  `NULL` to keep the CRS of projected `boundaries`, or to project
  `boundaries` in longitude and latitude to the CRS from
  [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
  (with a message naming the CRS and the units of distance).
  Alternatively, a projected CRS to project `boundaries` to, as in
  [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md)
  and
  [`summarise_distance()`](https://claudiofronterre.github.io/RiskMap/reference/summarise_distance.md).
  A grid built from `raster` takes the CRS of the raster, so
  `distance_crs` must then be `NULL` or that CRS.

- distance_units:

  Character string, either `"km"` or `"m"`, giving the units of
  `spacing`. Defaults to `"km"`.

- keep:

  Which cells to keep. `"centre"` (the default) keeps cells whose centre
  lies within `boundaries`. `"overlap"` keeps every cell that overlaps
  `boundaries`, including cells along the edge whose centre lies
  outside. Use `"overlap"` when totals over the study area matter, for
  example population counts along a coastline.

- raster:

  Optional `SpatRaster` from the terra package, in a projected
  coordinate reference system with square cells. When supplied, the grid
  consists of the raster's cells: their centres are the grid points and
  the raster's layers become columns of the grid.

## Value

A `RiskMap_grid` and `sf` object containing the generated grid points.
The `grid_metadata` attribute records the spacing, units, cell size,
cell area, cell shape, grid CRS, alignment origin, the cells kept and
the study boundary. The `grid_id` column provides a stable identifier
for each retained cell centre. A grid built from `raster` also has one
column for each raster layer.

## Details

This function returns point centres rather than polygon cells. The
returned object has class `RiskMap_grid` and carries lightweight
metadata describing the square cells. This allows other RiskMap
functions to reconstruct their footprints when exact cell coverage is
needed without storing a second copy of the grid as polygons. The
metadata also records the study boundary, so that maps can clip edge
cells to it.

With `keep = "overlap"`, cells along the edge extend beyond
`boundaries`. Areal predictions from
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
weight each cell by the fraction lying within each polygon, so these
cells contribute only their covered part.

A grid built from `raster` is aligned with that raster by construction,
so its values need no extraction. This suits population counts:
aggregate the counts to the resolution of the grid with
[`terra::aggregate()`](https://rspatial.github.io/terra/reference/aggregate.html)
and `fun = "sum"`, after projecting them with
[`terra::project()`](https://rspatial.github.io/terra/reference/project.html)
and `method = "sum"` if needed, and use the result as `raster`. With
`boundaries`, cells are kept according to `keep`; without, the cells
with at least one non-missing value are kept. Cells kept inside
`boundaries` may have missing values, for example over water.

The CRS of the grid is `distance_crs` when it is given, and otherwise
that of `boundaries`, or of `raster` when it is supplied. If
`boundaries` is in longitude and latitude and neither `distance_crs` nor
a raster is given, it is projected to the CRS from
[`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
and a message names the CRS and the units of distance. A warning of
class `riskmap_warning_crs` is raised when the CRS used distorts
distances by more than 1% over `boundaries`, or over the grid built from
`raster`.

The grid keeps its cells when rows or columns are selected with `[` or
with dplyr verbs such as `select()` and
[`filter()`](https://rdrr.io/r/stats/filter.html), unless only the
geometry column is kept. Functions that build a new `sf` object, such as
[`cbind()`](https://rdrr.io/r/base/cbind.html),
[`merge()`](https://rdrr.io/r/base/merge.html) and
[`sf::st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html),
drop the class `RiskMap_grid`, and
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
then assigns each point to the polygon that contains it. Add columns
with `$<-`,
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
or a dplyr join instead. [`rbind()`](https://rdrr.io/r/base/cbind.html)
keeps the cells of its first grid, so use it only to join parts of one
grid.

## See also

[`sf::st_make_grid()`](https://r-spatial.github.io/sf/reference/st_make_grid.html),
[`sf::st_intersection()`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html)

## Examples

``` r
library(sf)

# Example boundary data
nc <- st_read(system.file("shape/nc.shp", package="sf"))
#> Reading layer `nc' from data source 
#>   `/home/runner/work/_temp/Library/sf/shape/nc.shp' using driver `ESRI Shapefile'
#> Simple feature collection with 100 features and 14 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -84.32385 ymin: 33.88199 xmax: -75.45698 ymax: 36.58965
#> Geodetic CRS:  NAD27
nc <- st_transform(nc, crs = 32617)

# Create grid with 10 km spacing
grid <- create_grid(nc, spacing = 10)

# Also keep the edge cells whose centre lies outside the boundary
grid_edges <- create_grid(nc, spacing = 10, keep = "overlap")

# Plot the grid
plot(st_geometry(nc))
plot(grid, add = TRUE, col = 'red')


# Use the cells of a 20 km raster, with its values as a column
cells <- terra::align(terra::ext(nc), 20000)
counts <- terra::rast(cells, resolution = 20000, crs = "EPSG:32617")
terra::values(counts) <- seq_len(terra::ncell(counts))
names(counts) <- "count"
grid_raster <- create_grid(nc, raster = counts)
```
