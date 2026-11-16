# Print the Results of a Spatial Cross-Validation

Print method for objects of class `"RiskMap_cross_validation"` from
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).
It gives the models, the test sets, how the parameters were set and the
scores of each model over all test observations.
[`summary.RiskMap_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md)
gives the scores of each test set, and
[`plot_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_cross_validation.md)
maps the scores at the test locations or draws the calibration curve.

## Usage

``` r
# S3 method for class 'RiskMap_cross_validation'
print(x, ...)
```

## Arguments

- x:

  An object of class `"RiskMap_cross_validation"`.

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

Invisibly returns `x`.
