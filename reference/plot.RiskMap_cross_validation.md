# Plot Method for RiskMap_cross_validation Objects

Draws one metric of a `RiskMap_cross_validation` object returned by
[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md):
a map of the CRPS or SCRPS at the test locations, or the calibration
curve (`"PIT"`), the average nonrandomized PIT (AnPIT) curve for
discrete families or the empirical distribution function of the PIT
values for Gaussian models. Several models are panels of one plot that
share one colour scale.

## Usage

``` r
# S3 method for class 'RiskMap_cross_validation'
plot(
  x,
  metric = NULL,
  model = NULL,
  pit_mode = "average",
  pit_test_set = NULL,
  combine_pit = FALSE,
  ...
)
```

## Arguments

- x:

  A `RiskMap_cross_validation` object.

- metric:

  The metric to draw, one of `"PIT"` (the calibration curve), `"CRPS"`
  and `"SCRPS"` (maps of the score of each test observation), among
  those computed. When `NULL` (the default), the calibration curve if it
  was computed and the first score otherwise. `PIT_area` is one number
  for each model, which
  [`summary()`](https://rdrr.io/r/base/summary.html) gives, so it has no
  map.

- model:

  Character vector of model names to include. Defaults to every model in
  `x$model`.

- pit_mode:

  For the calibration curve only: one of `"average"` (the curve of all
  test observations together, the default), `"single"` (one specific
  test set) or `"all"` (every test set separately).

- pit_test_set:

  Integer; required when `pit_mode = "single"`.

- combine_pit:

  Logical; when `pit_mode = "average"`, draw the calibration curves for
  every model in a single panel (`TRUE`) rather than one panel per model
  (`FALSE`, default).

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored. Add ggplot2 components to the plot with `+`.

## Value

A ggplot2 object, drawn when it is printed, to which ggplot2 components
such as
[`ggplot2::scale_colour_viridis_c()`](https://ggplot2.tidyverse.org/reference/scale_viridis.html)
can be added.

## Details

Lower is better for both scores, and a location in several test sets
shows the mean of its scores. The average calibration curve weights each
test set by its size, so it is the curve of all test observations
together, from which `PIT_area` is computed (see
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)).

## See also

[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
