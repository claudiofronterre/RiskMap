# Plot Method for RiskMap_simulation Objects

Maps one component of one simulation from
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
over its prediction grid, by default the mean, with the sampling
locations overlaid when data were also simulated.

## Usage

``` r
# S3 method for class 'RiskMap_simulation'
plot(x, simulation = 1, quantity = "mean", palette = NULL, limits = NULL, ...)
```

## Arguments

- x:

  Output from
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
  including a surface.

- simulation:

  A single simulation number to plot, defaulting to one.

- quantity:

  The simulated component to map, one of `"mean"` (the default),
  `"linear_predictor"`, `"spatial_effect"`, `"nugget"` and
  `"group_effect"`, as in
  [`simulated_values()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_values.md).
  The mean is the prevalence for a binomial model, the mean per unit
  exposure for a Poisson or negative binomial model and the mean of the
  response for a Gaussian model. The nugget and the group effect are
  zero when the model has none.

- palette:

  Colour palette. `NULL` selects the defaults of
  [`plot.RiskMap_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_prediction.md),
  a calm-to-alarm sequential palette for the mean and a palette centred
  at zero for the linear predictor and the effects. A character vector
  of colours or the name of an HCL palette accepted by
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)
  can be supplied. Prefix a named palette with `"-"` to reverse it, for
  example `"-Viridis"`.

- limits:

  `NULL` (the default) or a numeric vector of two increasing numbers
  giving the range of the colour scale, for example to map several
  simulations on a common scale.

- ...:

  Additional arguments passed to the ggplot2 fill scale, for example
  `transform = "sqrt"` for skewed values. The map sets `colours`,
  `limits`, `breaks`, `labels`, `oob`, `na.value` and `guide` itself. To
  change them, add another fill scale to the map with `+`.

## Value

A `ggplot` object.

## Details

The colours and the legend follow
[`plot.RiskMap_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_prediction.md).
The returned object can be modified further with standard ggplot2
functions.

## See also

[`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md),
[`plot.RiskMap_simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_simulated_surface.md)
