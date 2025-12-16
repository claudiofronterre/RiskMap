# Prediction of the random effects components and covariates effects over a spatial grid using a fitted generalized linear Gaussian process model

This function computes predictions over a spatial grid using a fitted
model obtained from the
[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
function. It provides point predictions and uncertainty estimates for
the specified locations for each component of the model separately: the
spatial random effects; the unstructured random effects (if included);
and the covariates effects.

## Usage

``` r
pred_over_grid(
  object,
  grid_pred = NULL,
  predictors = NULL,
  re_predictors = NULL,
  pred_cov_offset = NULL,
  control_sim = set_control_sim(),
  type = "marginal",
  messages = TRUE
)
```

## Arguments

- object:

  A RiskMap object obtained from the \`glgpm\` function.

- grid_pred:

  An object of class 'sfc', representing the spatial grid over which
  predictions are to be made. Must be in the same coordinate reference
  system (CRS) as the object passed to 'object'.

- predictors:

  Optional. A data frame containing predictor variables used for
  prediction.

- re_predictors:

  Optional. A data frame containing predictors for unstructured random
  effects, if applicable.

- pred_cov_offset:

  Optional. A numeric vector specifying covariate offsets at prediction
  locations.

- control_sim:

  Control parameters for MCMC sampling. Must be an object of class
  "mcmc.RiskMap" as returned by
  [`set_control_sim`](https://claudiofronterre.github.io/RiskMap/reference/set_control_sim.md).

- type:

  Type of prediction. "marginal" for marginal predictions, "joint" for
  joint predictions.

- messages:

  Logical. If TRUE, display progress messages. Default is TRUE.

## Value

An object of class 'RiskMap.pred.re' containing predicted values,
uncertainty estimates, and additional information.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
