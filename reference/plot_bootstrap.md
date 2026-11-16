# Plot a parametric bootstrap

Draws a histogram of the bootstrap estimates of each parameter, with the
estimate, the percentile interval and the bias-corrected interval (see
[`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)).
Every parameter is shown on its natural axis, unless `log_scale = TRUE`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the result
does the same.

## Usage

``` r
plot_bootstrap(x, parameters = NULL, bins = 30, log_scale = FALSE, ...)

# S3 method for class 'RiskMap_bootstrap'
plot(x, ...)
```

## Arguments

- x:

  The result of
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

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  arguments of `plot_bootstrap()`. For `plot_bootstrap()`, not used, and
  an error if given, so that a misspelt argument is not ignored.

## Value

A `ggplot` object for one parameter, or a `patchwork` object combining
one `ggplot` per parameter. Nothing is drawn until the object is
printed.

## See also

[`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)

## Examples

``` r
data(italy_sim)
fit <- glgpm(y ~ gp(), data = italy_sim[1:50, ], family = "gaussian",
             messages = FALSE)
bootstrap <- parametric_bootstrap(fit, n_sim = 10, seed = 1,
                                  messages = FALSE)
plot_bootstrap(bootstrap, parameters = c("sigma2", "phi"),
               log_scale = TRUE)

# The same histograms
plot(bootstrap, parameters = c("sigma2", "phi"), log_scale = TRUE)
```
