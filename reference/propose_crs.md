# Propose a projected CRS for the data

Proposes the projected coordinate reference system (CRS) that keeps
distances between locations closest to their values on the Earth, for
computing distances and fitting geostatistical models. A WGS 84 UTM zone
is proposed when it is nearly as good as the best projection, and a
conformal projection centred on the data otherwise.

## Usage

``` r
propose_crs(
  data,
  epsg_tolerance = 0.005,
  messages = getOption("RiskMap.messages", TRUE)
)
```

## Arguments

- data:

  An `sf` or `sfc` object of points or polygons, in any CRS.

- epsg_tolerance:

  The largest distortion of distances, as a fraction, for which a UTM
  zone is proposed. The default of 0.005 accepts a UTM zone that
  distorts distances by at most 0.5% over the data.

- messages:

  Whether to report the CRS proposed and its distortion. Defaults to the
  option `RiskMap.messages`, or `TRUE` if it is unset.

## Value

The EPSG code of the UTM zone as an integer, or an `sf` `crs` object for
a centred projection. Its name begins with "RiskMap" and gives the
family, the centre (the standard parallel for the Lambert conformal
conic projection), the azimuth of the oblique Mercator projection and
the scale factor `k0`, for example "RiskMap oblique stereographic 23E
4.5S k0 0.9952". Either can be passed to
[`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
and to the `model_crs` argument of
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

## Details

No planar map keeps every distance, so a projected CRS stretches or
shrinks them by a factor that varies over the map. The distortion of a
CRS over the data is the largest relative error of short distances over
the region the data cover, measured on the WGS 84 ellipsoid (see
[`crs_distortion()`](https://claudiofronterre.github.io/RiskMap/reference/crs_distortion.md)).
It also bounds the relative error of every distance between two
locations of the data.

The region is the convex hull of the data, whether points or polygons.
Two kinds of CRS are compared.

- The WGS 84 UTM zone at the centre of the region, the centre of the
  smallest spherical cap enclosing it, and its two neighbours. The
  hemisphere is that of the centre, with the Equator counted as north,
  and does not change distances.

- Four conformal projections centred on the data: oblique stereographic,
  transverse Mercator, oblique Mercator and Lambert conformal conic.
  Their parameters minimise the largest scale error over the region
  (Chebyshev's criterion), and their scale factor balances the smallest
  and largest scale about 1. Their centres are rounded to 0.5 degrees
  and their scale factors to four decimals, so that data with similar
  extents usually get the same CRS.

Conformal projections keep angles, so a covariance that is isotropic on
the Earth stays isotropic locally in the plane, and only its scale
changes slightly over the map.

The best UTM zone is proposed when its distortion is at most
`epsg_tolerance`, and the best centred projection otherwise. A region
about 1,000 km across or less usually gets a UTM zone. The best centred
projection keeps distances within 0.5% over every African country. No
planar CRS keeps them within 1% over a region that holds a disc about
3,600 km across, such as sub-Saharan Africa, and a warning of class
`riskmap_warning_crs` then says so.

Functions that compute distances, such as
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md)
and
[`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md),
call `propose_crs()` for data in longitude and latitude. When several
analyses should share one CRS, for example analyses of different surveys
in one country, call `propose_crs()` once on the country's boundary and
pass the result as `model_crs` or `distance_crs`.

## References

Milnor J (1969). A problem in cartography. *American Mathematical
Monthly* 76, 1101-1112.

Snyder JP (1987). *Map Projections: A Working Manual*. USGS Professional
Paper 1395.

## See also

[`crs_distortion()`](https://claudiofronterre.github.io/RiskMap/reference/crs_distortion.md)
for the distortion of any CRS over the data.

## Examples

``` r
data("infect_sma")
propose_crs(infect_sma)
#> ℹ Proposed EPSG:32610 (WGS 84 / UTM zone 10N), which distorts distances by at
#>   most 0.015% over `data`.
#> [1] 32610

# A region too wide for a UTM zone gets a centred projection
library(sf)
wide <- st_as_sf(data.frame(lon = c(12, 31, 12, 31), lat = c(-13, -13, 5, 5)),
                 coords = c("lon", "lat"), crs = 4326)
crs <- propose_crs(wide)
#> ℹ Proposed a custom Lambert conformal conic projection, "RiskMap Lambert
#>   conformal conic 21.5E 4S k0 0.9937", which distorts distances by at most
#>   0.63% over `data`.
#> ℹ The best UTM zone, EPSG:32734 (WGS 84 / UTM zone 34S), distorts them by up to
#>   1.5%, more than `epsg_tolerance` allows.
crs$Name
#> [1] "RiskMap Lambert conformal conic 21.5E 4S k0 0.9937"
crs_distortion(wide, crs)
#>   max_error_pct max_anisotropy_pct min_scale max_scale n_locations
#> 1     0.6318523       2.107198e-06 0.9937025  1.006319        1883
```
