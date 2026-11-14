# Plot Method for RiskMap_simulation Objects

Maps the linear predictor of one simulation from
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
over its prediction grid, with the sampling locations overlaid when data
were also simulated.

## Usage

``` r
# S3 method for class 'RiskMap_simulation'
plot(x, simulation = 1, palette = "viridis", reverse_palette = FALSE, ...)
```

## Arguments

- x:

  Output from
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
  including a surface.

- simulation:

  A single simulation number to plot, defaulting to one.

- palette:

  Either the name of a palette from
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html)
  (e.g. `"viridis"`, `"Blues"`, `"Spectral"`) or a character vector of
  two or more colours to interpolate between. Defaults to `"viridis"`.

- reverse_palette:

  Logical; if `TRUE`, reverse the order of the palette colours. Defaults
  to `FALSE`.

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
[`plot.RiskMap_simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_simulated_surface.md)
