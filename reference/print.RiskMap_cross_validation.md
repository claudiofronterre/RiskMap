# Print the Results of a Spatial Cross-Validation

Print method for objects of class `"RiskMap_cross_validation"` from
[`assess_prediction`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).
It gives the models, the test sets, how the parameters were set and the
scores of each model averaged over the test sets.
[`summary.RiskMap_cross_validation`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md)
gives the scores of each test set, and
[`plot.RiskMap_cross_validation`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_cross_validation.md)
maps the scores at the test locations and draws the calibration curve.

## Usage

``` r
# S3 method for class 'RiskMap_cross_validation'
print(x, ...)
```

## Arguments

- x:

  An object of class `"RiskMap_cross_validation"`.

- ...:

  Not used.

## Value

Invisibly returns `x`.
