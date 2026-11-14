# Summaries of the distances

Computes the distances between the unique locations in the dataset and
returns summary statistics.

## Usage

``` r
summarise_distance(data, distance_crs = NULL, distance_units = c("km", "m"))

# S3 method for class 'RiskMap_distance_summary'
print(x, digits = 3, ...)
```

## Arguments

- data:

  an object of class `sf` containing point geometries.

- distance_crs:

  `NULL` to retain an existing projected CRS, or automatically reproject
  longitude/latitude data to the CRS from
  [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
  (with a message naming the CRS and the units of distance).
  Alternatively, a CRS to reproject `data` to, which must itself be
  projected. A warning of class `riskmap_warning_crs` is raised when the
  CRS used distorts distances by more than 1% over the data.

- distance_units:

  Character string, either `"km"` or `"m"`, indicating whether the
  distances used are expressed in kilometers or meters. Defaults to
  `"km"`.

- x:

  an object returned by `summarise_distance()`.

- digits:

  the number of significant digits to print. Defaults to 3.

- ...:

  further arguments, currently ignored.

## Value

a named numeric vector of class `RiskMap_distance_summary` containing
the following components

- `min`:

  the minimum distance

- `max`:

  the maximum distance

- `mean`:

  the mean distance

- `median`:

  the median distance

The values are stored in full. Its
[`print()`](https://rdrr.io/r/base/print.html) method shows them to
`digits` significant digits, with the number of unique locations and the
distance units.

## Examples

``` r
data(italy_sim)

summarise_distance(italy_sim)
#> ── Distances between locations ─────────────────────────────────────────────────
#> 200 unique locations, with distances in km.
#>    min    max   mean median 
#>   4.95  1,144    440    411 
```
