# Distortion of distances by a projected CRS

Measures how much a projected coordinate reference system (CRS)
stretches or shrinks distances over the region the data cover.

## Usage

``` r
crs_distortion(data, crs = NULL)
```

## Arguments

- data:

  An `sf` or `sfc` object of points or polygons, in any CRS.

- crs:

  The projected CRS to assess, as an EPSG code, a string or an `sf`
  `crs` object. With `NULL`, the default, the CRS of `data`, which must
  then be projected.

## Value

A data frame with one row and the columns

- `max_error_pct`:

  the distortion of distances, in percent

- `max_anisotropy_pct`:

  the anisotropy, in percent

- `min_scale`, `max_scale`:

  the smallest and largest scale factors over the region

- `n_locations`:

  the number of locations at which they were computed

## Details

At each location of the region, the CRS maps a small circle on the WGS
84 ellipsoid to an ellipse whose semi-axes, relative to the circle, are
the largest and smallest scale factors there (Tissot's indicatrix). They
are computed from the derivatives of the projection by central
differences over 10 m. The distortion is the largest of the differences
between these scale factors and 1 over the region. It also bounds the
relative error of the distance between any two locations of the region,
because the straight line between them in the plane and the shortest
path between them on the Earth both stay within the region, or very
close to it for a region thousands of kilometres across.

The region is the convex hull of the data, whether points or polygons,
so that it holds the paths between them, for example across the sea
between islands. It is represented by up to 2,000 vertices of the data,
its boundary and a grid of points inside it.

The anisotropy is the largest ratio of the two scale factors at a
location, minus 1. It is zero, up to rounding, for a conformal CRS, such
as UTM or any CRS from
[`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md),
so that distances from a location are changed by the same factor in
every direction. An equal-area CRS changes them more in some directions
than in others.

## See also

[`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)

## Examples

``` r
data("infect_sma")
# The UTM zone proposed for the data
crs_distortion(infect_sma, crs = propose_crs(infect_sma, messages = FALSE))
#>   max_error_pct max_anisotropy_pct min_scale max_scale n_locations
#> 1    0.01482633       1.490337e-06 0.9998517 0.9999534        1970
# Pseudo-Mercator, which is not conformal
crs_distortion(infect_sma, crs = 3857)
#>   max_error_pct max_anisotropy_pct min_scale max_scale n_locations
#> 1      28.80057          0.4108479  1.279101  1.288006        1970
```
