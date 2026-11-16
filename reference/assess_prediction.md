# Assess Predictive Performance via Spatial Cross-Validation

This function evaluates the predictive performance of spatial models
fitted to `RiskMap` objects using cross-validation. It supports two
classes of diagnostic tools:

- **Scoring rules**, the Continuous Ranked Probability Score (CRPS) and
  its scaled version (SCRPS) of Bolin and Wallin (2023), which measure
  both the sharpness and the calibration of the predictive
  distributions. Both are losses, so lower values are better.

- **Calibration diagnostics**, based on the Probability Integral
  Transform (PIT) for Gaussian outcomes, the average nonrandomized PIT
  (AnPIT) curve for discrete outcomes (e.g., Poisson or Binomial), and
  the area between the PIT/AnPIT curve and the diagonal (Giorgi *et al.*
  2026).

Cross-validation can be performed using either spatial clustering,
regularized subsampling with a minimum inter-point distance or a
user-defined test set. For each fold or subset, models can be refitted
or evaluated with fixed parameters, offering flexibility in model
validation.
[`plot_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_cross_validation.md),
or [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the
result, maps the scores and draws the calibration curves.

## Usage

``` r
assess_prediction(
  models,
  method,
  fold = NULL,
  min_dist = NULL,
  size = NULL,
  user_split = NULL,
  iter = 1,
  metrics = c("PIT", "CRPS", "SCRPS"),
  keep_par_fixed = TRUE,
  control_mcmc = set_control_mcmc(),
  plot_fold = FALSE,
  seed = NULL,
  messages = getOption("RiskMap.messages", TRUE),
  ...
)
```

## Arguments

- models:

  A list of `RiskMap` objects, each a model fitted with
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  to the same data. Names of the list name the models in the results.

- method:

  Character; either `"cluster"`, `"regularized"` or `"user"` for the
  cross-validation method:

  The `"cluster"` method uses spatial clustering as implemented by the
  `spatial_clustering_cv` function from the `spatialsample` package,
  which must be installed. The `"regularized"` method selects a
  subsample of the dataset by imposing a minimum distance, set by the
  `min_dist` argument, for a randomly selected subset of locations using
  the `subsample.distance` function from the `spatialEco` package, which
  must be installed. The `"user"` method takes a user-defined test set
  defined by `user_split`

- fold:

  Integer; required when `method = "cluster"` - number of folds for
  cross-validation.

- min_dist:

  Numeric; required when `method = "regularized"`. The minimum distance
  between the locations of a test set, in the distance units of the fits
  (`distance_units` in
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)),
  for example kilometres for fits in kilometres.

- size:

  Integer; the size of the test set, required when
  `method = "regularized"`.

- user_split:

  Required when `method = "user"`. A user-defined cross-validation
  split. Either:

  - a matrix with `nrow = n` (number of observations) and one column for
    each test set, where entries of `1` indicate membership in the test
    set and `0` in the training set; or

  - a list with one element for each test set, each either a vector of
    indices of the dataset to use as the test set, or a list with
    components `in_id` (training indices) and `out_id` (test indices).

- iter:

  Integer; number of times to repeat the cross-validation with
  `method = "cluster"` or `"regularized"`. Defaults to `1`. With
  `method = "user"` the number of test sets is that of `user_split`, so
  `iter` need not be given.

- metrics:

  Character vector; one or more of `"CRPS"`, `"SCRPS"` and `"PIT"`, to
  specify the predictive performance metrics to compute. `"PIT"`
  computes the PIT for Gaussian models or the AnPIT curve for discrete
  models, and also `PIT_area`, the area between the calibration curve of
  all test observations and the diagonal, with its expected value under
  perfect calibration. Defaults to all metrics.

- keep_par_fixed:

  Logical. If `TRUE` (the default), the parameters are fixed at the
  estimates of the fitted model, which used all the data, and only the
  predictions are conditioned on the training set. If `FALSE`, each
  model is re-estimated on each training set, which costs one fit per
  test set and model. Fixed parameters are much faster, and with many
  locations and few parameters the scores change little. The test
  observations still contribute to the fixed estimates, however, so the
  scores are optimistic, and more so for models with more parameters.
  Use `FALSE` to compare models with different numbers of parameters.
  The choice is stored in the result and printed by
  [`summary()`](https://rdrr.io/r/base/summary.html).

- control_mcmc:

  Control settings for simulation, an output from
  [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md).

- plot_fold:

  Logical; whether to draw the training and test sets of each split
  while the assessment runs. Defaults to `FALSE`.

- seed:

  Optional non-negative integer that makes the assessment reproducible.
  It sets the random splits for `method = "cluster"` and
  `method = "regularized"`, and a seed for each test set, which is
  passed to
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  when the models are refitted and to
  [`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md).
  All models share the seed of a test set, so that their predictions use
  common random numbers. The caller's random number state is restored on
  exit. When `NULL` (the default), the seeds are drawn from the
  session's random number generator.

- messages:

  Logical, whether to report on the assessment, with a progress bar as
  described in
  [riskmap_progress](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md).
  Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.

- ...:

  Additional arguments passed to clustering or subsampling functions.

## Value

A list of class `RiskMap_cross_validation`, containing:

- test_set:

  A list of test sets used for validation, each of class `'sf'`.

- model:

  A named list, one per model, each containing:

  metric

  :   A list with the CRPS and SCRPS of each test observation, one
      element per test set, if requested.

  PIT

  :   (if `family = "gaussian"` and `metrics` includes `"PIT"`) A list
      of PIT values for test data.

  AnPIT

  :   (if `family` is discrete and `metrics` includes `"PIT"`) A list of
      AnPIT curves for test data.

  PIT_area, PIT_area_expected

  :   (if `metrics` includes `"PIT"`) The area between the calibration
      curve of all test observations and the diagonal, and its expected
      value under perfect calibration.

- keep_par_fixed:

  The value of `keep_par_fixed`.

- failures:

  A data frame with one row for each model and test set whose refit or
  prediction failed, with columns `model`, `test_set` and `message`. The
  scores of these test sets are `NA`.

## Details

Every score is a loss, so lower values are better. The CRPS is in the
units of the response, which for binomial, Poisson and negative binomial
models are counts (for example numbers of positive individuals), not
prevalences or rates. The SCRPS is reported as \$\$\frac{1}{2}\left\[1 +
\frac{\mathrm{CRPS}(F, y)}{c} + \log(2c)\right\], \quad c = \frac{1}{2}
E\|X - X'\|,\$\$ with \\X\\ and \\X'\\ independent draws from the
predictive distribution \\F\\. This is the negative of the scaled CRPS
of Bolin and Wallin (2023), whose scores are higher for better
predictions, so that the CRPS and the SCRPS share one orientation. The
SCRPS can be negative.

`PIT_area` is the area between the calibration curve and the diagonal,
from 0 for a perfectly calibrated curve to at most 0.5. It is computed
once, from the curve of all test observations together, which is the
average of the curves of the test sets weighted by their sizes. Even a
perfectly calibrated model has a positive area, which shrinks as the
number of test observations grows (about \\0.31/\sqrt{n}\\ for \\n\\
observations). `PIT_area_expected` gives that area for the test
observations at hand, as the mean over 200 sets of observations
simulated from their own predictive distributions. An observation in
several test sets is simulated once and scored under each, as the
observed value is. A `PIT_area` well above `PIT_area_expected` points to
miscalibration.

A failed refit or prediction does not stop the assessment: the scores of
that test set are missing and a warning reports the failures. Test
observations with missing values in a variable used by any of the models
are left out. The predictive distributions of the test observations do
not yet include random effects, so the scores of models with
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
terms are not reliable, and a warning says so.

## References

Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of
proper scoring rules. *Statistical Science*, 38(1), 140–159.
[doi:10.1214/22-STS864](https://doi.org/10.1214/22-STS864) .

Giorgi, E., Fronterre, C. & Diggle, P. J. (2026). A decay-adjusted
spatio-temporal model to account for the impact of mass drug
administration on neglected tropical disease prevalence. *Journal of the
Royal Statistical Society Series A: Statistics in
Society*.[doi:10.1093/jrsssa/qnag100](https://doi.org/10.1093/jrsssa/qnag100)
.

## See also

[`plot_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_cross_validation.md),
which [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the
result also calls, and
[`summary.RiskMap_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md).

## Examples

``` r

data(italy_sim)

fit <- glgpm(
  formula = y ~ gp(),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

# user method with a matrix
cross_validation <-
  assess_prediction(
  list(fit),
  method = "user",
  user_split = matrix(
    sample(c(rep(1, 50), rep(0, 50))),
    ncol = 1)
  )
#> ℹ Parameters are fixed at the estimates from all the data.
#> ℹ Set `keep_par_fixed = FALSE` to re-estimate them in each training set.

summary(cross_validation)
#> ── Cross-validation scores ─────────────────────────────────────────────────────
#> Parameters are fixed at the estimates from all the data.
#> 
#> Scores over all test observations
#>    Model   CRPS  SCRPS PIT_area PIT_area_expected
#>  Model 1 0.3271 0.8130   0.0441            0.0427
#> ℹ Lower is better for every score. The CRPS is in the units of the response,
#>   counts for binomial, Poisson and negative binomial models. PIT_area_expected
#>   is the PIT_area of a perfectly calibrated model on these test observations.
#> 
#> By test set
#>    Model Test set   CRPS  SCRPS
#>  Model 1        1 0.3271 0.8130

# user method with a list
cross_validation <-
  assess_prediction(
  list(fit),
  method = "user",
  user_split = list(
    sample(100, 50))
  )
#> ℹ Parameters are fixed at the estimates from all the data.
#> ℹ Set `keep_par_fixed = FALSE` to re-estimate them in each training set.

summary(cross_validation)
#> ── Cross-validation scores ─────────────────────────────────────────────────────
#> Parameters are fixed at the estimates from all the data.
#> 
#> Scores over all test observations
#>    Model   CRPS  SCRPS PIT_area PIT_area_expected
#>  Model 1 0.2481 0.6520   0.0424            0.0455
#> ℹ Lower is better for every score. The CRPS is in the units of the response,
#>   counts for binomial, Poisson and negative binomial models. PIT_area_expected
#>   is the PIT_area of a perfectly calibrated model on these test observations.
#> 
#> By test set
#>    Model Test set   CRPS  SCRPS
#>  Model 1        1 0.2481 0.6520

# cluster method, which needs spatialsample
if (requireNamespace("spatialsample", quietly = TRUE)) {
  cross_validation <-
    assess_prediction(
      list(fit),
      method = "cluster",
      fold = 2
    )
  summary(cross_validation)
}
#> ℹ Parameters are fixed at the estimates from all the data.
#> ℹ Set `keep_par_fixed = FALSE` to re-estimate them in each training set.
#> ── Cross-validation scores ─────────────────────────────────────────────────────
#> Parameters are fixed at the estimates from all the data.
#> 
#> Scores over all test observations
#>    Model   CRPS  SCRPS PIT_area PIT_area_expected
#>  Model 1 1.2393 1.4537   0.0695            0.0318
#> ℹ Lower is better for every score. The CRPS is in the units of the response,
#>   counts for binomial, Poisson and negative binomial models. PIT_area_expected
#>   is the PIT_area of a perfectly calibrated model on these test observations.
#> 
#> By test set
#>    Model Test set   CRPS  SCRPS
#>  Model 1        1 2.1992 1.8316
#>  Model 1        2 0.9994 1.3592

# regularized method, which needs spatialEco
if (requireNamespace("spatialEco", quietly = TRUE)) {
  cross_validation <-
    assess_prediction(
      list(fit),
      method = "regularized",
      size = 5,
      min_dist = 1
    )
  summary(cross_validation)
}
#> ℹ Parameters are fixed at the estimates from all the data.
#> ℹ Set `keep_par_fixed = FALSE` to re-estimate them in each training set.
#> ── Cross-validation scores ─────────────────────────────────────────────────────
#> Parameters are fixed at the estimates from all the data.
#> 
#> Scores over all test observations
#>    Model   CRPS  SCRPS PIT_area PIT_area_expected
#>  Model 1 0.9860 1.3465   0.1312            0.0477
#> ℹ Lower is better for every score. The CRPS is in the units of the response,
#>   counts for binomial, Poisson and negative binomial models. PIT_area_expected
#>   is the PIT_area of a perfectly calibrated model on these test observations.
#> 
#> By test set
#>    Model Test set   CRPS  SCRPS
#>  Model 1        1 0.9860 1.3465
```
