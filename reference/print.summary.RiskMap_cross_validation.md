# Print the summary of a cross-validation

Prints the scores of each model averaged over all test observations,
with what they mean, from
[summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md).

## Usage

``` r
# S3 method for class 'summary.RiskMap_cross_validation'
print(x, ...)
```

## Arguments

- x:

  A summary of class `summary.RiskMap_cross_validation`.

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored. Some columns of the summary print as a data frame, and the
  arguments go to its print method.

## Value

Invisibly returns `x`.
