# Assess Simulations

This function evaluates the performance of models based on simulation
results from the \`surf_sim\` function.

## Usage

``` r
assess_sim(
  obj_sim,
  models,
  control_mcmc = set_control_sim(),
  spatial_scale,
  messages = TRUE,
  f_grid_target = NULL,
  f_area_target = NULL,
  shp = NULL,
  col_names = NULL,
  pred_objective = c("mse", "classify"),
  categories = NULL
)
```

## Arguments

- obj_sim:

  An object of class \`RiskMap.sim\`, obtained as an output from the
  \`surf_sim\` function.

- models:

  A named list of models to be evaluated.

- control_mcmc:

  A control object for MCMC sampling, created with
  \`set_control_sim()\`. Default is \`set_control_sim()\`.

- spatial_scale:

  The scale at which predictions are assessed, either \`"grid"\` or
  \`"area"\`.

- messages:

  Logical, if \`TRUE\` messages will be displayed during processing.
  Default is \`TRUE\`.

- f_grid_target:

  A function for processing grid-level predictions.

- f_area_target:

  A function for processing area-level predictions.

- shp:

  A shapefile of class \`sf\` or \`data.frame\` for area-level analysis,
  required if \`spatial_scale = "area"\`.

- col_names:

  Column name in \`shp\` containing unique region names. If \`NULL\`,
  defaults to \`"region"\`.

- pred_objective:

  A character vector specifying objectives, either \`"mse"\`,
  \`"classify"\`, or both.

- categories:

  A numeric vector of thresholds defining categories for classification.
  Required if \`pred_objective = "classify"\`.

## Value

A list of class \`RiskMap.sim.res\` containing model evaluation results.
