# Get EPSG code for appropriate UTM zone

Calculates the UTM zone and corresponding EPSG code for a spatial object
with WGS84 coordinates. This is useful for reprojecting spatial data
into a projected coordinate system suitable for distance-based
geostatistical analyses.

## Usage

``` r
get_epsg_utm(sf_object)
```

## Arguments

- sf_object:

  An \`sf\` object with a WGS84 CRS (EPSG:4326).

## Value

An integer indicating the EPSG code for the appropriate UTM zone (e.g.,
32629 for UTM Zone 29N).

## Details

The function computes the centroid of the spatial object, uses its
longitude to identify the UTM zone, and then assigns the corresponding
EPSG code. EPSG codes from 32601 to 32660 correspond to UTM zones in the
Northern Hemisphere, while 32701 to 32760 are used for the Southern
Hemisphere.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(RiskMap)
  data("liberia")
  liberia_sf <- st_as_sf(liberia,
                         coords = c("long", "lat"),
                         crs = 4326)
  get_epsg_utm(liberia_sf)
} # }
```
