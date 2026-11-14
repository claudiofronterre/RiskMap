# Spatial prediction from a fitted RiskMap model

Predicts the latent spatial process, linear predictor and response-scale
distribution at observed or new locations. New locations and every
variable required by the fitted model are supplied together in one `sf`
object, keeping their rows aligned.

## Usage

``` r
# S3 method for class 'RiskMap'
predict(
  object,
  newdata = NULL,
  boundaries = NULL,
  aggregate = NULL,
  type = c("auto", "marginal", "joint"),
  include = NULL,
  transform = NULL,
  standardize = NULL,
  summaries = c("mean", "lower", "upper"),
  thresholds = NULL,
  level = 0.95,
  weights = NULL,
  control_mcmc = NULL,
  parameter_uncertainty = FALSE,
  parameter_draws = 50L,
  seed = NULL,
  keep_samples = FALSE,
  messages = getOption("RiskMap.messages", TRUE),
  ...
)
```

## Arguments

- object:

  A fitted `RiskMap` object returned by
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

- newdata:

  An `sf` object with point geometries and the variables needed by the
  fitted model. When `NULL`, predictions are made at the observed
  locations.

- boundaries:

  Optional `sf` polygon object. When supplied, prediction locations are
  assigned to polygons and their joint samples are aggregated within
  each polygon. Grids returned by
  [`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md)
  use exact fractional cell coverage; other point objects use
  point-in-polygon assignment.

- aggregate:

  Character vector naming one or more within-polygon operations:
  `"mean"`, `"median"`, `"sum"`, `"min"`, or `"max"`. It must be
  supplied when `boundaries` is used. Several operations reuse the same
  joint predictive samples. Custom operations can be supplied in a named
  list, optionally mixed with built-ins, for example
  `list(mean = "mean", q90 = function(values, weights) ...)`. A custom
  function receives a location-by-draw matrix and the location weights
  and must return one numeric value per draw.

- type:

  Character string specifying `"auto"`, `"marginal"`, or `"joint"`
  prediction. `"auto"` uses marginal prediction for locations and joint
  prediction when `boundaries` are supplied.

- include:

  `NULL`, or a character vector specifying optional components to add to
  the complete linear predictor and response. Supported values are
  `"random_effects"`, for conditional predictions at known groups, and
  `"nugget"`, for a new realization of observation-level predictive
  variation at locations supplied through `newdata`. Fixed effects, the
  spatial process and any offset are always returned separately.

- transform:

  Optional named list of functions defining derived predictive targets
  from the complete linear-predictor draws. Each function receives a
  location-by-draw matrix and must return a numeric matrix with the same
  dimensions. The ordinary response-scale prediction is always
  calculated, so this argument is needed only for additional scientific
  targets, for example
  `list(stunting = function(link) pnorm((-2 - link) / sd_child))`.

- standardize:

  Optional data frame describing a common reference population over
  which response-scale predictions are averaged. Its columns must be
  fixed-effect predictors from the fitted model and override the
  corresponding columns in `newdata`. Rows receive equal weight unless
  an optional numeric `.weight` column is supplied. Predictor values
  kept together on each row are standardized jointly, preserving their
  empirical association. This answers what would be predicted if every
  location had the same reference-population distribution. `.weight`
  weights reference people or profiles and is distinct from the polygon
  `weights` argument.

- summaries:

  Character vector naming the summaries calculated from each predictive
  distribution. Supported values are `"mean"`, `"median"`, `"sd"`,
  `"cv"`, `"lower"`, `"upper"`, `"interval_width"`, and `"exceedance"`.
  The coefficient of variation is calculated only on the response scale
  and is returned as `NA` when the predictive mean is non-positive or
  numerically too close to zero for a stable ratio. Custom summaries can
  be supplied in a named list, optionally mixed with built-ins, for
  example `list(mean = "mean", iqr = IQR)`. Each custom function
  receives the predictive draws and must return one numeric value.

- thresholds:

  Numeric thresholds used when `summaries` includes `"exceedance"`. One
  exceedance probability is returned for each threshold on the response
  scale and for every target supplied through `transform`.

- level:

  Number strictly between zero and one giving the probability covered by
  the `"lower"` and `"upper"` interval limits and used to calculate
  `"interval_width"`.

- weights:

  `NULL` for equally weighted polygon summaries, or the name of a
  numeric column in `newdata` containing non-negative location weights.
  For `aggregate = "mean"` or `"median"`, weights are interpreted within
  each polygon. A weighted median is returned when weights are supplied.
  For `aggregate = "sum"`, their original scale is retained. Weights do
  not affect `"min"` or `"max"`.

- control_mcmc:

  Optional output from
  [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md).
  When omitted, conditional samples stored in a non-Gaussian fit are
  reused. Supplying a control object requests a fresh conditional run.

- parameter_uncertainty:

  Logical. If `FALSE` (the default), predictions are made at the
  estimated parameters (plug-in prediction). If `TRUE`, they also
  reflect the uncertainty in the estimated parameters; see Details.

- parameter_draws:

  Positive integer giving the number of parameter values drawn when
  `parameter_uncertainty = TRUE`. Default 50. Conditional samples stored
  in the fit are not used in that case, because they condition on the
  estimate only; `control_mcmc`, if supplied, sets only the total number
  of predictive draws.

- seed:

  Optional integer seed for random draws made during prediction,
  including a fresh conditional MCMC run. The caller's random-number
  state is restored on exit.

- keep_samples:

  Logical. If `TRUE`, retain the minimal latent samples needed for
  custom post-processing. Predictive samples are generated to calculate
  summaries regardless of this setting, but are discarded by default to
  limit memory use.

- messages:

  Logical. Display messages, and a progress bar while predicting at many
  locations, as described in
  [riskmap_progress](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md).
  Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.

- ...:

  Additional arguments reserved for future prediction methods.

## Value

A `RiskMap_prediction` object containing:

- predictions:

  The aligned prediction points or grid cells with their summaries.

- boundaries:

  The supplied polygons, or `NULL` when no polygon aggregation was
  requested.

- boundary_summary:

  For polygon prediction, a tidy table with one row per boundary,
  aggregation, component, summary and threshold. `boundary_id` is the
  row number of the corresponding polygon in `boundaries`.

- samples:

  `NULL` when `keep_samples = FALSE`; otherwise a list containing
  spatial samples and any requested random-effect or nugget samples for
  location-level prediction, and, with parameter uncertainty, the
  fixed-effect contribution of each draw (`fixed`). For polygon
  aggregation it is nested first by aggregation and then by component,
  with one matrix row per polygon.

- include:

  Optional components included in the link and response.

- transform:

  Names of additional derived predictive targets.

- standardized:

  Whether response predictions were averaged over a reference population
  supplied through `standardize`.

- aggregate:

  Names of the within-polygon aggregations.

- summaries:

  Names of the calculated summaries.

- thresholds:

  Thresholds used for exceedance probabilities.

- level:

  Probability covered by the interval limits.

- weights:

  Name of the weights column used for the polygon summaries, or `NULL`.

- type:

  Whether prediction was marginal or joint.

- parameter_draws:

  The number of parameter draws when parameter uncertainty was
  propagated, otherwise `NULL`.

- family:

  The fitted response family.

- link:

  The link of the fit when it is not the canonical one, its name or
  `"custom"`, and `NULL` for the canonical link.

- response_name:

  Name of the response variable, used to label maps of Gaussian
  predictions.

- observed_locations:

  Whether `newdata` was omitted.

- coverage:

  `NULL` for location-only prediction. For polygon aggregation,
  `"exact"` when fractional grid-cell coverage was used and `"point"`
  for point-in-polygon assignment. This is output metadata, not a
  user-supplied argument.

## Details

At new locations, the default prediction represents the underlying
spatial surface, with linear predictor \$\$d(x)^T \hat\beta + S(x) +
\mathrm{offset}(x).\$\$ The offset is included automatically when it was
specified in the fitted model. Nugget and unstructured random effects
are excluded unless requested explicitly because they generally
represent observation-level or group-specific variation rather than the
underlying spatial surface. The observed-location distinction for a
fitted nugget is described below.

For non-Gaussian fits, stored final-parameter conditional samples are
reused when `control_mcmc = NULL`. If the model was fitted with
`return_samples = FALSE`, the function draws a new conditional sample
using
[`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md)
defaults. Supplying `control_mcmc` always requests a fresh conditional
run. For Gaussian models, prediction is based on their exact conditional
Gaussian distribution.

Plug-in prediction treats the estimated parameters as known, so its
intervals can be too narrow when the parameters, in particular the
spatial variance and scale, are poorly estimated, as in small or
low-prevalence surveys. With `parameter_uncertainty = TRUE`,
`parameter_draws` parameter values are drawn from the approximate
sampling distribution of the estimate: a multivariate normal
distribution with the estimated covariance matrix, on the scale on which
the model is estimated (log variances and log scale). For each value,
conditional samples of the latent effects are generated at that value
and the prediction is repeated; the draws are then pooled, so every
summary, aggregation and derived target reflects both the latent and the
parameter uncertainty. For binomial, Poisson and negative binomial
models the conditional samples at each parameter value come from a short
run of the same Langevin (MALA) sampler used for plug-in prediction, in
coordinates whitened by the Laplace approximation at that value, started
from a draw of that approximation and using the step size tuned at the
estimate (100 burn-in iterations, then thinning by 10). The total number
of predictive draws equals that of plug-in prediction, 2,000 by default.
Parameters without a finite positive variance, such as those held at a
boundary, keep their estimated value. When the covariance matrix of the
estimates is not available, an error is returned; refitting with
`set_control_mcml(pc_prior = TRUE)` usually resolves this. The location
column `fixed` then reports the average fixed-effect contribution.

At the observed locations, the conditioned location-level Gaussian
effect already includes the fitted nugget when the model has one.
Consequently, `include = "nugget"` is available only for predictions
supplied through `newdata`, where it represents a new independent nugget
realization.

`newdata` and `boundaries` are transformed internally to the fitted
model CRS without changing row order; returned geometries remain in
their input CRS. Rows with missing required predictors are retained with
`NA` predictions and an informative warning.

When `boundaries` are supplied, one prediction run produces both the
location-level output in `predictions` and the areal output in
`boundary_summary`. Joint samples are generated in the smallest
independent blocks that preserve dependence within every polygon.
Polygons sharing a fractionally covered grid cell are placed in the same
block so that the cell's draw is reused rather than simulated twice.
This avoids constructing a full joint covariance over unrelated
prediction locations. Locations outside all polygons remain in
`predictions` but do not contribute to areal aggregation; polygons with
no usable locations remain in `boundaries` with `NA` summaries.

With location weights \\w_i\\, `"mean"` computes \$\$\sum_i w_i T_i /
\sum_i w_i,\$\$ so population values are converted internally to
within-polygon fractions. `"sum"` instead computes \\\sum_i w_i T_i\\
and never normalizes the weights. Thus, a prevalence prediction weighted
by population counts gives population-weighted prevalence under `"mean"`
and expected cases under `"sum"`. Population densities must first be
multiplied by the represented cell areas. When `newdata` was created by
[`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md),
each supplied weight is multiplied by the exact fraction of its square
cell covered by each polygon. Cells crossing a boundary can therefore
contribute to more than one polygon. For an arbitrary `sf` point object,
cell footprints cannot be inferred safely and each point is assigned to
the polygon containing it.

For location-level prediction, output columns use `<quantity>_<summary>`
names. For example, `response_mean`, `response_lower`, and
`response_upper` describe the response-scale prediction, while the
corresponding `spatial_*` and `link_*` columns describe the spatial
process and complete linear predictor. `fixed` is the estimated
fixed-effects contribution and `offset` is reported separately. When
requested through `include`, `random_effects_*` and `nugget_*` report
those contributions separately as well as incorporating them into
`link_*` and `response_*`. Binomial response predictions are
probabilities and Poisson response predictions are rates per unit
exposure.

Derived targets supplied through `transform` are evaluated draw by draw
before any location-level summary or polygon aggregation. Their output
uses the target name as the quantity prefix. With multiple thresholds,
location-level exceedance columns are named, for example,
`response_exceedance_0.2` and `response_exceedance_0.5`; polygon output
stores thresholds in the tidy `threshold` column.

When `standardize` is supplied, every conditional spatial draw is reused
across all rows of the reference population. The response and any
derived targets are evaluated for each reference profile and then
averaged using `.weight`. This is direct standardization (empirical or
counterfactual averaging), not simulation of uncertain predictors. The
reported `fixed` and `link_*` quantities use the weighted mean
fixed-effect contribution; with a nonlinear link, transforming that mean
is not the same as the standardized `response_*` quantity.

## See also

[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md)

## Examples

``` r
data(italy_sim)
example_data <- italy_sim[1:50, ]
fit <- glgpm(y ~ pop_dens + gp(), data = example_data, family = "gaussian",
             messages = FALSE)

# Calculate an additional predictive target from every linear-predictor draw
prediction <- predict(
  fit,
  newdata = example_data[1:5, ],
  transform = list(squared_link = function(link) link^2),
  messages = FALSE
)
prediction$predictions
#> Simple feature collection with 5 features and 19 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 165557.8 ymin: 4845668 xmax: 165557.8 ymax: 4845668
#> Projected CRS: WGS 84 / UTM zone 33N
#> # A tibble: 5 × 20
#>       y region  province pop_dens ID_loc           geometry fixed offset
#> * <dbl> <fct>   <fct>       <dbl>  <int>        <POINT [m]> <dbl>  <dbl>
#> 1  6.74 Toscana Pisa         536.      1 (165557.8 4845668)  8.63      0
#> 2  8.56 Toscana Pisa         536.      1 (165557.8 4845668)  8.63      0
#> 3  8.05 Toscana Pisa         536.      1 (165557.8 4845668)  8.63      0
#> 4  7.53 Toscana Pisa         536.      1 (165557.8 4845668)  8.63      0
#> 5  7.16 Toscana Pisa         536.      1 (165557.8 4845668)  8.63      0
#> # ℹ 12 more variables: spatial_mean <dbl>, spatial_lower <dbl>,
#> #   spatial_upper <dbl>, link_mean <dbl>, link_lower <dbl>, link_upper <dbl>,
#> #   response_mean <dbl>, response_lower <dbl>, response_upper <dbl>,
#> #   squared_link_mean <dbl>, squared_link_lower <dbl>, squared_link_upper <dbl>

# Also reflect the uncertainty in the estimated parameters
uncertain_prediction <- predict(
  fit,
  newdata = example_data[1:5, ],
  parameter_uncertainty = TRUE,
  parameter_draws = 20,
  messages = FALSE
)
uncertain_prediction$predictions
#> Simple feature collection with 5 features and 16 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 165557.8 ymin: 4845668 xmax: 165557.8 ymax: 4845668
#> Projected CRS: WGS 84 / UTM zone 33N
#> # A tibble: 5 × 17
#>       y region  province pop_dens ID_loc           geometry fixed offset
#> * <dbl> <fct>   <fct>       <dbl>  <int>        <POINT [m]> <dbl>  <dbl>
#> 1  6.74 Toscana Pisa         536.      1 (165557.8 4845668)  9.00      0
#> 2  8.56 Toscana Pisa         536.      1 (165557.8 4845668)  9.00      0
#> 3  8.05 Toscana Pisa         536.      1 (165557.8 4845668)  9.00      0
#> 4  7.53 Toscana Pisa         536.      1 (165557.8 4845668)  9.00      0
#> 5  7.16 Toscana Pisa         536.      1 (165557.8 4845668)  9.00      0
#> # ℹ 9 more variables: spatial_mean <dbl>, spatial_lower <dbl>,
#> #   spatial_upper <dbl>, link_mean <dbl>, link_lower <dbl>, link_upper <dbl>,
#> #   response_mean <dbl>, response_lower <dbl>, response_upper <dbl>

# Use the same reference distribution of population density at each location
reference_population <- data.frame(pop_dens = example_data$pop_dens)
standardized_prediction <- predict(
  fit,
  newdata = example_data[1:5, ],
  standardize = reference_population,
  messages = FALSE
)
standardized_prediction$predictions
#> Simple feature collection with 5 features and 16 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 165557.8 ymin: 4845668 xmax: 165557.8 ymax: 4845668
#> Projected CRS: WGS 84 / UTM zone 33N
#> # A tibble: 5 × 17
#>       y region  province pop_dens ID_loc           geometry fixed offset
#> * <dbl> <fct>   <fct>       <dbl>  <int>        <POINT [m]> <dbl>  <dbl>
#> 1  6.74 Toscana Pisa         536.      1 (165557.8 4845668)  9.93      0
#> 2  8.56 Toscana Pisa         536.      1 (165557.8 4845668)  9.93      0
#> 3  8.05 Toscana Pisa         536.      1 (165557.8 4845668)  9.93      0
#> 4  7.53 Toscana Pisa         536.      1 (165557.8 4845668)  9.93      0
#> 5  7.16 Toscana Pisa         536.      1 (165557.8 4845668)  9.93      0
#> # ℹ 9 more variables: spatial_mean <dbl>, spatial_lower <dbl>,
#> #   spatial_upper <dbl>, link_mean <dbl>, link_lower <dbl>, link_upper <dbl>,
#> #   response_mean <dbl>, response_lower <dbl>, response_upper <dbl>

# Aggregate joint predictions within two example areas
areas <- sf::st_sf(
  area = 1:2,
  geometry = sf::st_make_grid(example_data, n = c(2, 1))
)
areal_prediction <- predict(
  fit,
  newdata = example_data,
  boundaries = areas,
  aggregate = "mean",
  messages = FALSE
)
#> Warning: The conditional prediction covariance of 2 of 2 blocks needed a diagonal
#> jitter, of up to 1e-12, to be factorised.
areal_prediction$boundary_summary
#>    boundary_id aggregate quantity  summary threshold       value
#> 1            1      mean    fixed estimate        NA 10.30974195
#> 2            1      mean   offset estimate        NA  0.00000000
#> 3            1      mean  spatial     mean        NA -0.35294890
#> 4            1      mean  spatial    lower        NA -0.51861707
#> 5            1      mean  spatial    upper        NA -0.17903651
#> 6            1      mean     link     mean        NA  9.95679304
#> 7            1      mean     link    lower        NA  9.79112488
#> 8            1      mean     link    upper        NA 10.13070543
#> 9            1      mean response     mean        NA  9.95679304
#> 10           1      mean response    lower        NA  9.79112488
#> 11           1      mean response    upper        NA 10.13070543
#> 12           2      mean    fixed estimate        NA  9.36720129
#> 13           2      mean   offset estimate        NA  0.00000000
#> 14           2      mean  spatial     mean        NA  0.24243285
#> 15           2      mean  spatial    lower        NA  0.03567343
#> 16           2      mean  spatial    upper        NA  0.44124090
#> 17           2      mean     link     mean        NA  9.60963414
#> 18           2      mean     link    lower        NA  9.40287471
#> 19           2      mean     link    upper        NA  9.80844219
#> 20           2      mean response     mean        NA  9.60963414
#> 21           2      mean response    lower        NA  9.40287471
#> 22           2      mean response    upper        NA  9.80844219
```
