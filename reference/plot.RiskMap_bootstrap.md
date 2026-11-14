# Plot a parametric bootstrap

Draws a histogram of the bootstrap estimates of each parameter, with the
estimate, the percentile interval and the bias-corrected interval (see
[`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)).
Every parameter is shown on its natural axis, unless `log_scale = TRUE`.

## Usage

``` r
# S3 method for class 'RiskMap_bootstrap'
plot(x, parameters = NULL, bins = 30, log_scale = FALSE, ...)
```

## Arguments

- x:

  A `RiskMap_bootstrap` object from
  [`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md).

- parameters:

  Names of the parameters to draw, as in the columns of `x$estimates`.
  The default is all of them.

- bins:

  Number of bins of each histogram.

- log_scale:

  Logical, whether to show the variances and the scale on a logarithmic
  axis, which spreads out bootstrap estimates near zero. Estimates at
  zero are then left out. The regression coefficients stay on their
  natural axis.

- ...:

  Not used.

## Value

A `ggplot` object for one parameter, or a `patchwork` object combining
one `ggplot` per parameter.

## See also

[`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)
