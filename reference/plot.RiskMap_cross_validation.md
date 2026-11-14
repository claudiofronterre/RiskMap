# Plot Method for RiskMap_cross_validation Objects

Plots whatever predictive-performance output is available in a
`RiskMap_cross_validation` object returned by
[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md):
a spatial map for each per-location metric (`"CRPS"`, `"SCRPS"`,
`"PIT_area"`) that was requested via `metrics` in that call, and/or a
calibration curve (`"PIT"`) - the Aggregated normalised Probability
Integral Transform (AnPIT) curve for discrete families, or the empirical
PIT curve for Gaussian models - when that was requested. Plot titles and
legends use AnPIT or PIT according to each model's family. By default
every available plot is produced, one per requested model, so you don't
need to know in advance which metrics were computed. Each metric gets
its own plot (or, with more than one model, its own small grid of one
panel per model) rather than everything being squeezed into one combined
grid.

## Usage

``` r
# S3 method for class 'RiskMap_cross_validation'
plot(
  x,
  metric = NULL,
  model = NULL,
  ...,
  pit_mode = "average",
  pit_test_set = NULL,
  combine_pit = FALSE
)
```

## Arguments

- x:

  A `RiskMap_cross_validation` object.

- metric:

  Character vector restricting which plot(s) to produce; one or more of
  `"CRPS"`, `"SCRPS"`, `"PIT_area"` (a spatial map of that score) and
  `"PIT"` (the calibration curve). Defaults to every metric available in
  `x`.

- model:

  Character vector of model names to include. Defaults to every model in
  `x$model`.

- ...:

  Additional ggplot2 components (e.g. `scale_color_gradient()`), layered
  onto every spatial-map plot. Must be passed by position only after
  `model`, since the remaining arguments must be named.

- pit_mode:

  For the calibration curve only: one of `"average"` (average curve
  across test sets, the default), `"single"` (one specific test set) or
  `"all"` (every test set separately).

- pit_test_set:

  Integer; required when `pit_mode = "single"`.

- combine_pit:

  Logical; when `pit_mode = "average"`, draw the calibration curves for
  every model in a single panel (`TRUE`) rather than one panel per model
  (`FALSE`, default).

## Value

When exactly one metric is produced: a single ggplot2 object (one model)
or a grid object from gridExtra (more than one model). When more than
one metric is produced, each metric's plot/grid is drawn in turn as a
side effect (so each one renders as its own figure), and a named list of
them - one per metric - is returned invisibly for later reuse.

## See also

[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
