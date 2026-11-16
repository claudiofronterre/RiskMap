# Assess Simulations

Evaluates models using joint data and surface simulations.

## Usage

``` r
assess_simulation(
  simulation,
  models,
  spatial_scale,
  target = identity,
  aggregate = NULL,
  boundaries = NULL,
  label = NULL,
  objective = c("mse", "classify"),
  breaks = NULL,
  control_mcmc = set_control_mcmc(),
  seed = NULL,
  messages = getOption("RiskMap.messages", TRUE)
)
```

## Arguments

- simulation:

  Output from
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
  with both data and surface, from a generating model without grouped
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  effects. The generating model can have any link.

- models:

  A named list of fitted models returned by
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  or model specifications returned by
  [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md).
  Each object defines a candidate model specification and is refitted to
  every simulated dataset; its estimated parameter values are not
  reused.

- spatial_scale:

  The scale(s) at which predictions are assessed: `"location"`, the
  locations of the prediction grid, `"area"`, or `c("location", "area")`
  to compute both from a single fit-and-predict pass over the
  simulations.

- target:

  A function that converts a numeric matrix of means to a numeric matrix
  of the scientific target, with the same dimensions. The mean is the
  probability for a binomial model, the mean per unit exposure for a
  Poisson model and the value itself for a Gaussian model. The default,
  [`identity()`](https://rdrr.io/r/base/identity.html), scores the mean
  itself, such as prevalence. Other targets are functions of the mean,
  such as [`qlogis()`](https://rdrr.io/r/stats/Logistic.html) to score
  prevalence on the logit scale, or `function(mu) pnorm((-2 - mu) / sd)`
  for the proportion of children stunted from a Gaussian model of
  height-for-age z-scores with residual standard deviation `sd`. See
  Details.

- aggregate:

  How the targets at the locations within an area are combined into the
  target of the area, as in
  [`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md):
  one of `"mean"`, `"median"`, `"sum"`, `"min"` and `"max"`, or a
  function that receives a location-by-draw matrix and the location
  weights and returns one number for each draw. Locations have equal
  weights. Required only when `spatial_scale` includes `"area"`.

- boundaries:

  An `sf` object containing only POLYGON or MULTIPOLYGON geometries,
  required if `spatial_scale` includes `"area"`.

- label:

  The name of the column of `boundaries` that holds distinct names of
  the areas. If `NULL`, the areas are named `"reg1"`, `"reg2"` and so
  on, in a column `"region"`.

- objective:

  A character vector specifying objectives, either `"mse"`,
  `"classify"`, or both.

- breaks:

  A numeric vector of at least three increasing class boundaries for
  classification, as in
  [`plot_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot_prediction.md).
  Required if `objective` includes `"classify"`.

- control_mcmc:

  A control object for MCMC sampling, created with
  [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md).
  Default is
  [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md).

- seed:

  Optional non-negative integer that makes the assessment reproducible;
  see Details. The caller's random number state is restored on exit.
  When `NULL` (the default), the seeds of the simulations are drawn from
  the session's random number generator.

- messages:

  Logical, whether to report on the assessment, with a progress bar as
  described in
  [riskmap_progress](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md).
  Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.

## Value

A list of class `RiskMap_assess_simulation`:

- objective:

  One element per requested `spatial_scale` (`"location"` and/or
  `"area"`), each holding `mse`, a matrix with one row per model and one
  column per simulation, and/or `classify`, the classification metrics
  of each model and simulation.

- fits:

  A data frame with one row per model and simulation, giving the
  `status` of the fit (`"ok"` or `"failed"`), the error `message` of a
  failed fit and any `warnings`.

- estimates:

  A list with, for each model, a matrix of the parameter estimates, one
  row per simulation.

- n_sim, spatial_scale:

  The number of simulations and the spatial scales assessed.

## Details

Every model is scored on the scale of the mean, each through its own
inverse link, so the links of the generating model and of the candidates
need not agree. The true target at each grid location is `target()`
applied to the simulated mean, which is the generating model's inverse
link of the simulated linear predictor. That linear predictor includes
the covariates, any offset and, when the generating model has one, the
realised nugget. Each candidate is refitted to every simulated dataset,
its predictive draws of the linear predictor are mapped to the mean
through the inverse link of the refitted model, and `target()` is then
applied to those draws. A probit candidate on simulations from a logit
model is thus scored on the same prevalence as the truth. Candidates
predict with their own covariates and offset, and those with a nugget
include it in their predictions when the generating model has one.

The target receives the mean rather than the linear predictor, which is
what the `transform` argument of
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
receives, because models with different links share no linear-predictor
scale. Every link is monotone, so a target defined on the linear
predictor of a given link can be written as a function of the mean, as
[`qlogis()`](https://rdrr.io/r/stats/Logistic.html) gives the logit
scale.

The target of an area is `aggregate` applied to the targets at the grid
locations that fall within it, with equal weight, unlike the exact
fractional coverage used by
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
with `boundaries`. Areas without grid locations are left out of the
area-level assessment, with a warning. Predictive draws are joint within
each area, which is all that area-level targets need.

Each simulation has its own seed, derived from `seed` or, without one,
from the current random number generator state, and shared by all models
so that their fits and predictions use common random numbers. A fit or
prediction that fails does not stop the assessment: its results are
missing, and the failure is recorded in `fits` and reported in a
warning, as are warnings from the fits.

## Examples

``` r
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
data(italy_sim)
italy_subset <- italy_sim[1:30, ]
italy_grid <- italy_subset[!duplicated(st_coordinates(italy_subset)), ]
model <- specify_glgpm(
  y ~ pop_dens + gp(), italy_subset, family = "gaussian",
  parameters = list(beta = c(1, 0.001), sigma2 = 1, phi = 20,
                    sigma2_me = 0.1)
)
simulations <- simulate_glgpm(
  model, nsim = 1, what = c("data", "surface"),
  prediction_grid = italy_grid, seed = 1
)
boundary <- st_sf(
  region = "study_area",
  geometry = st_convex_hull(st_union(italy_grid))
)
assessment <- assess_simulation(
  simulations,
  models = list(candidate = model),
  spatial_scale = c("location", "area"),
  aggregate = "mean",
  boundaries = boundary,
  objective = "mse",
  messages = FALSE
)
```
