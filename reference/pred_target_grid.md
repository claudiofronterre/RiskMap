# Predictive Target Over a Regular Spatial Grid

Computes predictions over a regular spatial grid using outputs from the
[`pred_over_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_over_grid.md)
function. This function allows for incorporating covariates, offsets,
MDA effects, and optional unstructured random effects into the
predictive target.

## Usage

``` r
pred_target_grid(
  object,
  include_covariates = TRUE,
  include_nugget = FALSE,
  include_cov_offset = FALSE,
  include_mda_effect = TRUE,
  mda_grid = NULL,
  time_pred = NULL,
  include_re = FALSE,
  f_target = NULL,
  pd_summary = NULL
)
```

## Arguments

- object:

  Output from \`pred_over_grid\`, a RiskMap.pred.re object.

- include_covariates:

  Logical. Include covariates in the predictive target.

- include_nugget:

  Logical. Include the nugget effect in the predictive target.

- include_cov_offset:

  Logical. Include the covariate offset in the predictive target.

- include_mda_effect:

  Logical. Include the MDA effect in the predictive target using a DAST
  model; see
  [`dast`](https://claudiofronterre.github.io/RiskMap/reference/dast.md).

- mda_grid:

  Optional. Grid of MDA coverage values required for predictions using a
  DAST model; see
  [`dast`](https://claudiofronterre.github.io/RiskMap/reference/dast.md).

- time_pred:

  Optional. Time point for prediction required for predictions using a
  DAST model; see
  [`dast`](https://claudiofronterre.github.io/RiskMap/reference/dast.md).

- include_re:

  Logical. Include unstructured random effects in the predictive target.

- f_target:

  Optional. List of functions to apply on the linear predictor samples.

- pd_summary:

  Optional. List of summary functions to apply on the predicted values.

## Value

An object of class 'RiskMap_pred_target_grid' containing predicted
values and summaries over the regular spatial grid.

## See also

[`pred_over_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_over_grid.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
