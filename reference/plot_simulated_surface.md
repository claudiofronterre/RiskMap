# Map a simulated surface

Maps one component of a surface extracted with
[`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md)
over its prediction grid, by default the mean.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the surface
does the same.

## Usage

``` r
plot_simulated_surface(
  x,
  quantity = "mean",
  palette = NULL,
  limits = NULL,
  ...
)

# S3 method for class 'RiskMap_simulated_surface'
plot(x, ...)
```

## Arguments

- x:

  The result of
  [`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md).

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
  [`plot_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot_prediction.md),
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

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  arguments of `plot_simulated_surface()`. For
  `plot_simulated_surface()`, arguments passed to the ggplot2 fill
  scale, for example `transform = "sqrt"` for skewed values. The map
  sets `colours`, `limits`, `breaks`, `labels`, `oob`, `na.value` and
  `guide` itself. To change them, add another fill scale to the map with
  `+`.

## Value

A `ggplot` object.

## Details

The colours and the legend follow
[`plot_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot_prediction.md),
and the legend title names the mapped component. A prediction grid that
is regular in the CRS of the model is drawn with tiles, and other
prediction locations as points. The returned object can be modified
further with standard ggplot2 functions.

## See also

[`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md),
[`plot_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_simulation.md)

## Examples

``` r
data(italy_sim)
example_data <- italy_sim[1:50, ]
fit <- glgpm(y ~ gp(), data = example_data, family = "gaussian",
             messages = FALSE)
grid <- create_grid(create_convex_hull(example_data), spacing = 20)
simulation <- simulate_glgpm(fit, nsim = 2, what = c("data", "surface"),
                             prediction_grid = grid, seed = 1)
surface <- simulated_surface(simulation, simulation = 2)
plot_simulated_surface(surface, quantity = "spatial_effect")

# The same map
plot(surface, quantity = "spatial_effect")
```
