# Extract simulated values in tidy form

Extract simulated values in tidy form

## Usage

``` r
simulated_values(
  object,
  component = NULL,
  simulation = NULL,
  location_set = NULL
)
```

## Arguments

- object:

  Output from
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md).

- component:

  Components to extract. Choose from `"spatial_effect"`, `"nugget"`,
  `"group_effect"`, `"linear_predictor"`, `"mean"` and, for data
  simulations, `"response"`. Defaults to all available components.

- simulation:

  Simulation numbers to extract. Defaults to all simulations.

- location_set:

  `"data"`, `"surface"`, or both. Defaults to available sets.

## Value

A tibble with `simulation`, `location_set`, `location_id`, `component`
and `value`. Location identifiers refer to rows in `object$locations`.

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
