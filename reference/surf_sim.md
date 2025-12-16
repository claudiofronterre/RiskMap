# Simulate surface data based on a spatial model

This function simulates surface data based on a user-defined formula and
other parameters. It allows for simulation of spatial data with various
model families (Gaussian, Binomial, or Poisson). The simulation involves
creating spatially correlated random fields and generating outcomes for
data points in a given prediction grid.

## Usage

``` r
surf_sim(
  n_sim,
  pred_grid,
  formula,
  sampling_f,
  family,
  scale_to_km = TRUE,
  control_mcmc = set_control_sim(),
  par0,
  nugget_over_grid = FALSE,
  include_covariates = TRUE,
  fix_var_me = NULL,
  messages = TRUE
)
```

## Arguments

- n_sim:

  The number of simulations to run.

- pred_grid:

  A spatial object (either \`sf\` or \`data.frame\`) representing the
  prediction grid where the simulation will take place.

- formula:

  A formula object specifying the model to be fitted. It should include
  both fixed effects and random effects if applicable.

- sampling_f:

  A function that returns a sampled dataset (of class \`sf\` or
  \`data.frame\`) to simulate data from.

- family:

  A character string specifying the family of the model. Must be one of
  "gaussian", "binomial", or "poisson".

- scale_to_km:

  A logical indicating whether the coordinates should be scaled to
  kilometers. Defaults to \`TRUE\`.

- control_mcmc:

  A list of control parameters for MCMC (not used in this implementation
  but can be expanded later).

- par0:

  A list containing initial parameter values for the simulation,
  including \`beta\`, \`sigma2\`, \`phi\`, \`tau2\`, and \`sigma2_me\`.

- nugget_over_grid:

  A logical indicating whether to include a nugget effect over the
  entire prediction grid.

- include_covariates:

  A logical indicateing if the covariates (or the intercept if no
  covariates are used) should be included in the linear predictor. By
  default `include_covariates = TRUE`

- fix_var_me:

  A parameter to fix the variance of the random effects for the
  measurement error. Defaults to \`NULL\`.

- messages:

  A logical value indicating whether to print messages during the
  simulation. Defaults to \`TRUE\`.

## Value

A list containing the simulated data (`data_sim`), the linear predictors
(`lp_grid_sim`), a logical value indicating if covariates have been
included in the linear predictor (`include_covariates`), a logical value
indicating if the nugget has been included into the simulations of the
linear predictor over the grid (`nugget_over_grid`), a logical
indicating if a covariate offset has been included in the linear
predictor (`include_cov_offset`), the model parameters set for the
simulation (`par0`) and the family used in the model (`family`).

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>
