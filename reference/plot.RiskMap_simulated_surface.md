# Plot Method for RiskMap_simulated_surface Objects

Maps the linear predictor of a surface extracted with
[`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md)
over its prediction grid.

## Usage

``` r
# S3 method for class 'RiskMap_simulated_surface'
plot(x, palette = "viridis", reverse_palette = FALSE, limits = NULL, ...)
```

## Arguments

- x:

  Output from
  [`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md).

- palette:

  Either the name of a palette from
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html)
  (e.g. `"viridis"`, `"Blues"`, `"Spectral"`) or a character vector of
  two or more colours to interpolate between. Defaults to `"viridis"`.

- reverse_palette:

  Logical; if `TRUE`, reverse the order of the palette colours. Defaults
  to `FALSE`.

- limits:

  `NULL` (the default) or a numeric vector of length two giving the
  range of the colour scale, for example to map several simulations on a
  common scale.

- ...:

  Additional arguments passed to
  [`ggplot2::scale_fill_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html),
  e.g. `name` or `breaks`.

## Value

A `ggplot` object.

## Details

The returned object can be modified further with standard ggplot2
functions.

## See also

[`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md),
[`plot.RiskMap_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_simulation.md)
