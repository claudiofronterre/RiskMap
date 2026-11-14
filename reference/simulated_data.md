# Extract a simulated dataset or surface

Extract a simulated dataset or surface

## Usage

``` r
simulated_data(object, simulation = 1)

simulated_surface(object, simulation = 1)
```

## Arguments

- object:

  Output from
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md).

- simulation:

  A single simulation number, defaulting to one. These accessors return
  one `sf` object at a time; use
  [`simulated_values()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_values.md)
  to extract multiple simulations in one tidy table.

## Value

`simulated_data()` returns the sample `sf` data with its response column
replaced by the selected simulation, ready for refitting.
`simulated_surface()` returns the grid `sf` data with columns for the
spatial effect, nugget, group effect, linear predictor and mean, with
class `RiskMap_simulated_surface` so that it can be mapped with
[plot()](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_simulated_surface.md).

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
```
