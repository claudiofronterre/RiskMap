# Summaries of the distances

Computes the distances between the locations in the data-set and returns
summary statistics of these.

## Usage

``` r
dist_summaries(data, convert_to_utm = TRUE, scale_to_km = FALSE)
```

## Arguments

- data:

  an object of class `sf` containing the variable for which the
  variogram is to be computed and the coordinates

- convert_to_utm:

  a logical value, indicating if the conversion to UTM shuold be
  performed (`convert_to_utm = TRUE`) or the coordinate reference system
  of the data must be used without any conversion
  (`convert_to_utm = FALSE`). By default `convert_to_utm = TRUE`. Note:
  if `convert_to_utm = TRUE` the conversion to UTM is performed using
  the epsg provided by
  [`propose_utm`](https://claudiofronterre.github.io/RiskMap/reference/propose_utm.md).

- scale_to_km:

  a logical value, indicating if the distances used in the variogram
  must be scaled to kilometers (`scale_to_km = TRUE`) or left in meters
  (`scale_to_km = FALSE`). By default `scale_to_km = FALSE`

## Value

a list containing the following components

`min` the minimum distance

`max` the maximum distance

`mean` the mean distance

`median` the minimum distance
