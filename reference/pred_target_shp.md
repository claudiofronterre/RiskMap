# Predictive Targets over a Shapefile (grid-aggregated)

Computes predictive targets over polygon features using joint prediction
samples from
[`pred_over_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_over_grid.md).
Targets can incorporate covariates, offsets, optional unstructured
random effects, and (if fitted) mass drug administration (MDA) effects
from a DAST model.

## Usage

``` r
pred_target_shp(
  object,
  shp,
  shp_target = mean,
  weights = NULL,
  standardize_weights = FALSE,
  col_names = NULL,
  include_covariates = TRUE,
  include_nugget = FALSE,
  include_cov_offset = FALSE,
  include_mda_effect = TRUE,
  return_shp = TRUE,
  time_pred = NULL,
  mda_grid = NULL,
  include_re = FALSE,
  f_target = NULL,
  pd_summary = NULL,
  messages = TRUE,
  return_target_samples = FALSE
)
```

## Arguments

- object:

  Output from
  [`pred_over_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_over_grid.md)
  (class `RiskMap.pred.re`), typically fitted with `type = "joint"` so
  that linear predictor samples are available.

- shp:

  An sf polygon object (preferred) or a `data.frame` with an attached
  geometry column, representing regions over which predictions are
  aggregated.

- shp_target:

  A function that aggregates grid-cell values within each polygon to a
  single regional value (default `mean`). Examples: `mean`, `sum`, a
  custom weighted mean, etc.

- weights:

  Optional numeric vector of weights used inside `shp_target`. If
  supplied with `standardize_weights = TRUE`, weights are normalized
  within each region.

- standardize_weights:

  Logical; standardize `weights` within each region (`FALSE` by
  default).

- col_names:

  Name or column index in `shp` containing region identifiers to use in
  outputs.

- include_covariates:

  Logical; include fitted covariate effects in the linear predictor
  (default `TRUE`).

- include_nugget:

  Logical; include the nugget (unstructured measurement error) in the
  linear predictor (default `FALSE`).

- include_cov_offset:

  Logical; include any covariate offset term (default `FALSE`).

- include_mda_effect:

  Logical; include the MDA effect as defined by the fitted DAST model
  (default `TRUE`). Requires `time_pred` and, when applicable,
  `mda_grid`.

- return_shp:

  Logical; if `TRUE`, return the shapefile with appended summary columns
  defined by `pd_summary` (default `TRUE`).

- time_pred:

  Optional numeric scalar (or time index) at which to evaluate the
  predictive target

- mda_grid:

  Optional structure describing MDA schedules aligned with prediction
  grid cells (e.g., a `data.frame`/matrix/list). Used only when
  `include_mda_effect = TRUE`.

- include_re:

  Logical; include unstructured random effects (RE) in the linear
  predictor (default `FALSE`).

- f_target:

  List of target functions applied to linear predictor samples (e.g.,
  `list(prev = plogis)` for prevalence on the probability scale). If
  `NULL`, the identity is used.

- pd_summary:

  Named list of summary functions applied to each region's target
  samples (e.g.,
  `list(mean = mean, sd = sd, q025 = function(x) quantile(x, 0.025), q975 = function(x) quantile(x, 0.975))`).
  Names are used as column suffixes in the outputs.

- messages:

  Logical; if `TRUE`, print progress messages while computing regional
  targets.

- return_target_samples:

  Logical; if `TRUE`, also return the raw target samples per region
  (default `FALSE`).

## Value

An object of class `RiskMap_pred_target_shp` with components:

- `target`: `data.frame` of region-level summaries (one row per region).

- `target_samples`: (optional) `list` with one element per region; each
  contains a `data.frame`/matrix of raw samples for each named target in
  `f_target`, if `return_target_samples = TRUE`.

- `shp`: (optional) the input `sf` object with appended summary columns,
  included if `return_shp = TRUE`.

- `f_target`, `pd_summary`, `grid_pred`: inputs echoed for
  reproducibility.

## Details

For each polygon in `shp`, grid-cell samples of the linear predictor are
transformed with `f_target`, optionally adjusted for covariates, offset,
nugget, MDA effects and/or REs, and then aggregated via `shp_target`
(optionally weighted). The list `pd_summary` is applied to each region's
target samples to produce summary statistics.

## See also

[`pred_over_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_over_grid.md),
[`pred_target_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_target_grid.md)
