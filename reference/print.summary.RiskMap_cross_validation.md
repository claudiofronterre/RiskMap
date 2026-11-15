# Print Summary of RiskMap Spatial Cross-Validation Scores

This function prints the matrix of cross-validation scores produced by
`summary.RiskMap_cross_validation` in a readable format.

## Usage

``` r
# S3 method for class 'summary.RiskMap_cross_validation'
print(x, ...)
```

## Arguments

- x:

  An object of class `"summary.RiskMap_cross_validation"`, typically the
  output of `summary.RiskMap_cross_validation`.

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

Invisibly returns `x`.

## Details

This method prints the scores of each model averaged over all test
observations, with what they mean, and the scores of each test set.
