# Summarize Cross-Validation Scores for Spatial RiskMap Models

This function summarizes cross-validation scores for different spatial
models obtained from
[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).

## Usage

``` r
# S3 method for class 'RiskMap_cross_validation'
summary(object, view_all = TRUE, ...)
```

## Arguments

- object:

  A `RiskMap_cross_validation` object containing cross-validation scores
  for each model, as obtained from
  [`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).

- view_all:

  Logical. If `TRUE`, stores the average scores across test sets for
  each model alongside the overall average across all models. Defaults
  to `TRUE`.

- ...:

  Additional arguments passed to or from other methods.

## Value

A matrix of summary scores with models as rows and metrics as columns,
with class `"summary.RiskMap_cross_validation"`.

## Details

The function computes and returns a matrix where rows correspond to
models and columns correspond to performance metrics (e.g., CRPS,
SCRPS). Scores are weighted by subset sizes to compute averages.
Attributes of the returned object include:

- `test_set_means`: A list of average scores for each test set and
  model.

- `overall_averages`: Overall averages for each metric across all
  models.

- `view_all`: Indicates whether averages across test sets are available
  for visualization.

## See also

[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
