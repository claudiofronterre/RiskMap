# Map indicators of extrapolation

Maps indicators returned by
[`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md),
in the style of
[`plot.RiskMap_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_prediction.md).

## Usage

``` r
# S3 method for class 'RiskMap_extrapolation'
plot(x, indicator = "extrapolation", ...)
```

## Arguments

- x:

  A `RiskMap_extrapolation` object.

- indicator:

  Indicator to map: `"extrapolation"` (the default, the four classes),
  `"spatial_unexplained"`, `"covariate_novelty"` or
  `"covariate_leverage"`. Several indicators are arranged in two
  columns.

- ...:

  Not used.

## Value

A `ggplot` object, to which further layers can be added. When several
indicators are requested, a `patchwork` object combining one `ggplot`
per indicator.
