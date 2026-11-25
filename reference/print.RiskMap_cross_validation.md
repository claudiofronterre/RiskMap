# Print a cross-validation

Gives the models, the test sets, how the parameters were set and the
scores of each model over all test observations.

## Usage

``` r
# S3 method for class 'RiskMap_cross_validation'
print(x, ...)
```

## Arguments

- x:

  The result of
  [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

Invisibly returns `x`.

## See also

[summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md)
compares the models with one of them, and
[`plot_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_cross_validation.md)
draws the calibration and the scores.
