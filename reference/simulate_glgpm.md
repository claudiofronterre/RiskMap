# Simulate replicated data and spatial surfaces

Generates unconditional realizations from a fitted or specified model,
keeping its parameters fixed. This supports parametric bootstrap,
replicated data checks and simulation studies with a known spatial
surface.

## Usage

``` r
simulate_glgpm(
  object,
  nsim = 1,
  what = "data",
  sample_locations = NULL,
  prediction_grid = NULL,
  seed = NULL
)
```

## Arguments

- object:

  A fitted `RiskMap` object or a model from
  [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md).

- nsim:

  Number of independent simulations. One is useful for illustrating the
  mechanism. Reliable bootstrap intervals and performance estimates need
  many simulations and an assessment of Monte Carlo uncertainty.

- what:

  `"data"`, `"surface"`, or `c("data", "surface")`.

- sample_locations:

  Optional `sf` point data for simulated responses. Defaults to the
  model's original data. New data must contain the model's covariates,
  grouping variables, offsets and any denominator column.

- prediction_grid:

  `sf` point data containing the covariates, grouping variables and
  offsets for the requested surface. Required for `"surface"`.

- seed:

  Optional non-negative integer seed. When supplied, the caller's
  random-number state is restored on exit.

## Value

A `RiskMap_simulation` object containing `locations`, the transformed
`sf` data used for each requested location set; `samples`, numeric
arrays indexed by location, simulation and component; `nsim`; the
simulation `model`; the requested `what`; and `seed`. Use
[`simulated_data()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md)
to obtain a dataset ready for refitting,
[`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md)
for a surface, and
[`simulated_values()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_values.md)
for a tidy table of selected components.

## Details

The spatial process is drawn jointly at the union of the requested
locations using the model's Matern covariance. There is no nearest-grid
approximation and no conditioning on the observed responses. This
differs from conditional spatial prediction with
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md).

The returned components distinguish the spatial effect `S`, the nugget
and grouped effects, the full linear predictor (including offsets), its
inverse-link `mean`, and the `response`. For binomial models `mean` is a
probability, and for Poisson models it is a rate per unit exposure.
Denominators or exposures enter response generation. Gaussian
measurement error enters only the response. A surface does not include
sampled responses.

Repeated coordinates share the spatial effect and location-level nugget.
Matching group labels share a newly simulated group effect, also across
sample and grid locations. Measurement error and response sampling are
independent across observations conditional on the latent effects.

Locations are transformed to the model CRS when necessary, with a
message reporting the conversion. Exact joint simulation uses a dense
covariance matrix, so memory grows quadratically with the number of
unique locations.

## Examples

``` r
library(sf)
locations <- st_as_sf(data.frame(x = c(0, 1000, 2000), y = 0),
                      coords = c("x", "y"), crs = 32629)

grid <- st_as_sf(data.frame(x = c(500, 1500), y = 500),
                 coords = c("x", "y"), crs = 32629)

model <- specify_glgpm(response ~ gp(),
                       data = locations,
                       family = "gaussian",
                       parameters = list(beta = 1,
                                         sigma2 = 1,
                                         phi = 2,
                                         sigma2_me = 0.1))

sim <- simulate_glgpm(model,
                      nsim = 2,
                      what = c("data", "surface"),
                      prediction_grid = grid,
                      seed = 1)

simulated_data(sim, simulation = 1)
#> Simple feature collection with 3 features and 1 field
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 0 xmax: 2000 ymax: 0
#> Projected CRS: WGS 84 / UTM zone 29N
#>         geometry     response
#> 1    POINT (0 0)  0.851613371
#> 2 POINT (1000 0)  0.889323292
#> 3 POINT (2000 0) -0.002729952

simulated_values(sim, component = "spatial_effect")
#> # A tibble: 10 × 5
#>    simulation location_set location_id component        value
#>         <int> <chr>              <int> <chr>            <dbl>
#>  1          1 data                   1 spatial_effect -0.626 
#>  2          1 data                   2 spatial_effect -0.234 
#>  3          1 data                   3 spatial_effect -0.806 
#>  4          2 data                   1 spatial_effect -0.820 
#>  5          2 data                   2 spatial_effect -0.110 
#>  6          2 data                   3 spatial_effect  0.520 
#>  7          1 surface                1 spatial_effect  0.585 
#>  8          1 surface                2 spatial_effect -0.0475
#>  9          2 surface                1 spatial_effect -0.0238
#> 10          2 surface                2 spatial_effect  0.0448
```
