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

  Additional arguments passed to or from other methods.

## Value

Invisibly returns `x`.

## Details

This method is primarily used to format and display the summary score
matrix, printing it to the console. It provides a clear view of the
cross-validation performance metrics across different spatial models.
