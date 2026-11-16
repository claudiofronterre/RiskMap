# Plot a cross-validation

Draws one metric of a cross-validation by
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md),
as [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the
result does: a map of the CRPS or SCRPS at the test locations, or the
calibration curve (`"PIT"`), the average nonrandomized PIT (AnPIT) curve
for discrete families or the empirical distribution function of the PIT
values for Gaussian models. Several models are panels of one plot that
share one colour scale.

## Usage

``` r
plot_cross_validation(
  x,
  metric = NULL,
  model = NULL,
  pit_mode = "average",
  pit_test_set = NULL,
  combine_pit = FALSE,
  ...
)

# S3 method for class 'RiskMap_cross_validation'
plot(x, ...)
```

## Arguments

- x:

  The result of
  [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).

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

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  arguments of `plot_cross_validation()`. For `plot_cross_validation()`,
  not used, and an error if given, so that a misspelt argument is not
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

[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)

## Examples

``` r
data(italy_sim)
fit <- glgpm(y ~ gp(), data = italy_sim[1:100, ], family = "gaussian",
             messages = FALSE)
cross_validation <- assess_prediction(list(fit), method = "user",
                                      user_split = list(1:50, 51:100),
                                      seed = 1, messages = FALSE)
plot_cross_validation(cross_validation)

plot_cross_validation(cross_validation, metric = "CRPS")

# The same calibration curve
plot(cross_validation)
```
