# Plot an empirical variogram

Plots the empirical variogram computed by
[`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
[`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md)
result does the same.

## Usage

``` r
plot_variogram(
  x,
  envelope = TRUE,
  colour = "royalblue1",
  show_pairs = FALSE,
  ...
)

# S3 method for class 'RiskMap_variogram'
plot(x, ...)
```

## Arguments

- x:

  The result of
  [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md).

- envelope:

  Whether to draw the global envelope for spatial independence from the
  permutations of
  [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md).
  Defaults to `TRUE`. If
  [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md)
  was run with `n_permutations = 0`, there is no envelope, so a warning
  is raised and the variogram is drawn without it.

- colour:

  The fill colour of the envelope. See `vignette("ggplot2-specs")` for
  the ways a colour can be given.

- show_pairs:

  A logical value indicating if each point should be labelled with the
  number of pairs of locations in its distance class. A semivariance
  estimated from few pairs is unstable, so the labels help to judge how
  much each point can be trusted, and to choose the distance classes in
  [`variogram`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md).
  By default `show_pairs = FALSE`.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  arguments of `plot_variogram()`. For `plot_variogram()`, not used, and
  an error if given, so that a misspelt argument is not ignored.

## Value

A `ggplot` object representing the empirical variogram plot, optionally
including the envelope of spatial independence.

## Details

This function plots the empirical variogram, which shows the spatial
dependence structure of the data. Each point is the semivariance of one
distance class, and classes without pairs are omitted. Points are not
joined by lines, because the empirical variogram is only defined for the
distance classes. Both axes start at zero. With `envelope = TRUE`, the
plot also includes a simultaneous extreme-rank-length envelope for
spatial independence, and points that fall outside it are drawn in red.

## See also

[`variogram`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md)

## Examples

``` r
data(galicia)
galicia_variogram <- variogram(galicia, variable = "lead",
                               distance_units = "km", seed = 1)
plot_variogram(galicia_variogram, show_pairs = TRUE)

# The same plot
plot(galicia_variogram, show_pairs = TRUE)
```
