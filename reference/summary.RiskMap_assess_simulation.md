# Summarize Simulation Results

Summarizes the results of model evaluations from a
`RiskMap_assess_simulation` object. Provides average metrics for
classification by category and overall correct classification (CC)
summary.

## Usage

``` r
# S3 method for class 'RiskMap_assess_simulation'
summary(object, ...)
```

## Arguments

- object:

  An object of class `RiskMap_assess_simulation`, as returned by
  `assess_simulation`.

- ...:

  Additional arguments (not used).

## Value

A list of class `summary.RiskMap_assess_simulation` with one element per
assessed spatial scale (`grid` and/or `area`), each holding:

- `mse`: a data frame with, for each model, the number of simulations
  (`n_sim`), the number with a result (`n_valid`), and the mean and
  standard deviation of the MSE across simulations (`MSE_mean`,
  `MSE_sd`);

- `classify`: for each model, `classify_res`, the average sensitivity,
  specificity, PPV, NPV and CC by category across simulations,
  `n_valid`, the number of simulations with a defined value of each
  metric, `n_sim`, and `cc_summary`, the mean, standard deviation and
  2.5% and 97.5% quantiles of the overall proportion of correct
  classification.

It also holds `fits`, the number of model fits and of failed fits.
