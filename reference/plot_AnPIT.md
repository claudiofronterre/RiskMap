# Plot Calibration Curves (AnPIT / PIT) from Spatial Cross-Validation

Produce calibration plots from a `RiskMap.spatial.cv` object returned by
[`assess_pp`](https://claudiofronterre.github.io/RiskMap/reference/assess_pp.md).
\* For Binomial or Poisson models the function visualises the
*Aggregated normalised Probability Integral Transform* (AnPIT) curves
stored in `$AnPIT`. \* For Gaussian models it detects the list `$PIT`
and instead plots the empirical *Probability Integral Transform* curve
(ECDF of PIT values) on the same \\u\\-grid.

A 45° dashed red line indicates perfect calibration.

## Usage

``` r
plot_AnPIT(
  object,
  mode = "average",
  test_set = NULL,
  model_name = NULL,
  combine_panels = FALSE
)
```

## Arguments

- object:

  A `RiskMap.spatial.cv` object.

- mode:

  One of `"average"` (average curve across test sets), `"single"` (a
  specific test set), or `"all"` (every test set separately).

- test_set:

  Integer; required when `mode = "single"`.

- model_name:

  Optional character string; if supplied, only that model is plotted.

- combine_panels:

  Logical; when `mode = "average"`, draw all models in a single panel
  (`TRUE`) or one panel per model (`FALSE`, default).

## Value

A ggplot2 object (single plot) or a ggpubr grid.
