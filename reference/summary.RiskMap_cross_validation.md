# Summarise a cross-validation

Gives the average scores of each model over all test observations, their
differences from a reference model, and the areas of the calibration
curves.

## Usage

``` r
# S3 method for class 'RiskMap_cross_validation'
summary(object, reference = NULL, ...)
```

## Arguments

- object:

  The result of
  [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).

- reference:

  Optional name of the model that the others are compared with.

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

A data frame of class `summary.RiskMap_cross_validation` with one row
for each model, named after it, and the columns `CRPS` and `SCRPS`,
`PIT_area`, `PIT_area_expected` and `PIT_area_limit`, and, when
`reference` is given, `CRPS_difference`, `CRPS_se`, `SCRPS_difference`
and `SCRPS_se`. Its attributes include `reference`, `refit`, and for
each model `n_observations`, the number of observations with scores,
`n_failed`, the number of failed test-set assessments, and `binary`,
whether it is a model of binary outcomes.

## Details

Lower is better for both scores, and a negative difference means a
better score than that of the reference model. A `PIT_area` below
`PIT_area_limit` is one that a correct model would give (see
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)).
For a model of binary outcomes the area is not a check of calibration,
so the print of the summary and
[`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)
give it without its limit, and say why.

An observation in several test sets is first given the average of its
scores, so that it counts once. The scores of a model are then averaged
over the observations. The differences from the reference are those of
the observations that both models have scores for, and their standard
error is the standard deviation of these differences divided by the
square root of their number. It treats the observations as independent,
so it is only approximate, and it is too small when the differences at
neighbouring locations are alike.

The scores of each test set are in the `scores` element of the result of
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).

## See also

[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md),
and
[`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)
for a table to render.
