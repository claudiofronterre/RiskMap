# Assess Predictive Performance via Spatial Cross-Validation

This function evaluates the predictive performance of spatial models
fitted to \`RiskMap\` objects using cross-validation. It supports two
classes of diagnostic tools:

\- \*\*Scoring rules\*\*, including the Continuous Ranked Probability
Score (CRPS) and its scaled version (SCRPS), which quantify the
sharpness and calibration of probabilistic forecasts; - \*\*Calibration
diagnostics\*\*, based on the Probability Integral Transform (PIT) for
Gaussian outcomes and Aggregated nonparametric PIT (AnPIT) curves for
discrete outcomes (e.g., Poisson or Binomial).

Cross-validation can be performed using either spatial clustering or
regularized subsampling with a minimum inter-point distance. For each
fold or subset, models can be refitted or evaluated with fixed
parameters, offering flexibility in model validation. The function also
provides visualizations of the spatial distribution of test folds.

## Usage

``` r
assess_pp(
  object,
  keep_par_fixed = TRUE,
  iter = 1,
  fold = NULL,
  n_size = NULL,
  control_sim = set_control_sim(),
  method,
  min_dist = NULL,
  plot_fold = TRUE,
  messages = TRUE,
  which_metric = c("AnPIT", "CRPS", "SCRPS"),
  user_split = NULL,
  ...
)
```

## Arguments

- object:

  A list of \`RiskMap\` objects, each representing a model fitted with
  \`glgpm\`.

- keep_par_fixed:

  Logical; if \`TRUE\`, parameters are kept fixed across folds,
  otherwise the model is re-estimated for each fold.

- iter:

  Integer; number of times to repeat the cross-validation.

- fold:

  Integer; number of folds for cross-validation (required if \`method =
  "cluster"\`).

- n_size:

  Optional; the size of the test set, required if \`method =
  "regularized"\`.

- control_sim:

  Control settings for simulation, an output from \`set_control_sim\`.

- method:

  Character; either \`"cluster"\` or \`"regularized"\` for the
  cross-validation method. The \`"cluster"\` method uses spatial
  clustering as implemented by the `spatial_clustering_cv` function from
  the \`spatialEco\` package, while the \`"regularized"\` method selects
  a subsample of the dataset by imposing a minimum distance, set by the
  \`min_dist\` argument, for a randomly selected subset of locations.

- min_dist:

  Optional; minimum distance for regularized subsampling (required if
  \`method = "regularized"\`).

- plot_fold:

  Logical; if \`TRUE\`, plots each fold's test set.

- messages:

  Logical; if \`TRUE\`, displays progress messages.

- which_metric:

  Character vector; one or more of \`"CRPS"\`, \`"SCRPS"\`, or
  \`"AnPIT"\`, to specify the predictive performance metrics to compute.

- user_split:

  A user-defined cross-validation split. Either: \* a matrix with
  `nrow = n` (number of observations) and `ncol = iter` (number of
  iterations), where entries of `1` indicate membership in the test set
  for that iteration and `0` indicate training set; or \* a list of
  length `iter`, where each element is either a vector of test indices,
  or a list with components `in_id` (training indices) and `out_id`
  (test indices). When supplied, `user_split` overrides the automatic
  clustering or regularized distance splitting defined by `method`.

- ...:

  Additional arguments passed to clustering or subsampling functions.

## Value

A list of class \`RiskMap.spatial.cv\`, containing:

- test_set:

  A list of test sets used for validation, each of class \`'sf'\`.

- model:

  A named list, one per model, each containing:

  score

  :   A list with CRPS and/or SCRPS scores for each fold if requested.

  PIT

  :   (if \`family = "gaussian"\` and \`which_metric\` includes
      \`"AnPIT"\`) A list of PIT values for test data.

  AnPIT

  :   (if \`family\` is discrete and \`which_metric\` includes
      \`"AnPIT"\`) A list of AnPIT curves for test data.

## References

Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of
proper scoring rules. \*Statistical Science\*, 38(1), 140–159.
[doi:10.1214/22-STS864](https://doi.org/10.1214/22-STS864) .

## See also

`spatial_clustering_cv`, `subsample.distance`,
[`plot_AnPIT`](https://claudiofronterre.github.io/RiskMap/reference/plot_AnPIT.md)

## Author

Emanuele Giorgi
