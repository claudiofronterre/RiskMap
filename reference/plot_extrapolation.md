# Map indicators of extrapolation

Maps indicators returned by
[`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md),
in the style of
[`plot_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot_prediction.md).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the result
does the same.

## Usage

``` r
plot_extrapolation(x, indicator = "extrapolation", ...)

# S3 method for class 'RiskMap_extrapolation'
plot(x, ...)
```

## Arguments

- x:

  The result of
  [`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md).

- indicator:

  Indicator to map: `"extrapolation"` (the default, the four classes,
  with the thresholds that
  [`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md)
  used), `"spatial_unexplained"` (the share of the spatial variance left
  unexplained by the data, on a scale from 0 to 1),
  `"covariate_novelty"` or `"covariate_leverage"`. Several indicators
  are arranged in two columns.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  arguments of `plot_extrapolation()`. For `plot_extrapolation()`, not
  used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

A `ggplot` object, to which further layers can be added. When several
indicators are requested, a `patchwork` object combining one `ggplot`
per indicator.

## Examples

``` r
data(italy_sim)
example_data <- italy_sim[1:150, ]
fit <- glgpm(y ~ log(pop_dens) + gp(), data = example_data,
             family = "gaussian", messages = FALSE)
extrapolation <- assess_extrapolation(fit, italy_sim[151:200, ])
plot_extrapolation(extrapolation)
#> Warning: Removed 50 rows containing missing values or values outside the scale range
#> (`geom_sf()`).

# The same map
plot(extrapolation)
#> Warning: Removed 50 rows containing missing values or values outside the scale range
#> (`geom_sf()`).
```
