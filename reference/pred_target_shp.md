# Predictive Target over a Shapefile

Computes predictions over a shapefile using outputs from the
[`pred_over_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_over_grid.md)
function. This function allows for incorporating covariates, offsets,
and optional unstructured random effects into the predictive target.

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

  Output from \`pred_over_grid\`, a RiskMap.pred.re object.

- shp:

  Spatial dataset (sf or data.frame) representing the shapefile over
  which predictions are computed.

- shp_target:

  Function defining the aggregation method for shapefile targets
  (default is mean).

- weights:

  Optional numeric vector of weights for spatial predictions.

- standardize_weights:

  Logical indicating whether to standardize weights (default is FALSE).

- col_names:

  Column name or index in 'shp' containing region names.

- include_covariates:

  Logical indicating whether to include covariates in predictions
  (default is TRUE).

- include_nugget:

  Logical indicating whether to include the nugget effect (default is
  FALSE).

- include_cov_offset:

  Logical indicating whether to include covariate offset in predictions
  (default is FALSE).

- include_mda_effect:

  Logical indicating whether to include the mass drug administration
  (MDA) effect as defined by the fitted DAST model (default is TRUE).

- return_shp:

  Logical indicating whether to return the shape file with the added
  predictive distribution summaries as defined through `pd_summary`.

- include_re:

  Logical indicating whether to include random effects in predictions
  (default is FALSE).

- f_target:

  List of target functions to apply to the linear predictor samples.

- pd_summary:

  List of summary functions (e.g., mean, sd) to summarize target
  samples.

- return_target_samples:

  Logical indicating whether to return raw samples of the predictive
  targets for each region (default is FALSE).

## Value

An object of class 'RiskMap_pred_target_shp' containing:

- `target` – summaries of predictive targets by region.

- `target_samples` – raw samples of predictive targets (if
  `return_target_samples=TRUE`).

- `shp` – spatial object with appended summary statistics.

- `f_target`, `pd_summary`, `grid_pred`.

## Details

This function computes predictive targets or summaries over a spatial
shapefile using outputs from 'pred_S'. It requires the 'terra' package
for spatial data manipulation and should be used with 'sf' or
'data.frame' objects representing the shapefile.

## See also

[`pred_target_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_target_grid.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
