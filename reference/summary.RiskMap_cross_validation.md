# Summarize Cross-Validation Scores for Spatial RiskMap Models

This function summarizes cross-validation scores for different spatial
models obtained from
[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).

## Usage

``` r
# S3 method for class 'RiskMap_cross_validation'
summary(object, ...)
```

## Arguments

- object:

  A `RiskMap_cross_validation` object containing cross-validation scores
  for each model, as obtained from
  [`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

A matrix of summary scores with models as rows and metrics as columns,
with class `"summary.RiskMap_cross_validation"`.

## Details

The function returns a matrix whose rows are the models and whose
columns are the scores, `CRPS` and `SCRPS` averaged over all test
observations that have a score, and `PIT_area` with `PIT_area_expected`,
its value under perfect calibration, when calibration was assessed (see
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)).
Lower is better for every score. Attributes of the returned object
include:

- `test_set_means`: A list of average scores for each test set and
  model.

- `keep_par_fixed`: Whether the parameters were fixed at the estimates
  from all the data.

- `n_failed`: The number of failed test-set assessments.

Its print gives the averages and the scores of each test set.

## See also

[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
