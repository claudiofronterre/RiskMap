# Print the assessment of models against simulated data

Prints the models, the number of simulated datasets and, for each model
and statistic, where the observed value lies among the simulated ones,
with what the joint envelopes show.

## Usage

``` r
# S3 method for class 'RiskMap_predictive_check'
print(x, ...)

# S3 method for class 'summary.RiskMap_predictive_check'
print(x, ...)
```

## Arguments

- x:

  The result of
  [`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md),
  or its summary.

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored. Some columns of the summary print as a data frame, and the
  arguments go to its print method.

## Value

Invisibly returns `x`.

## See also

[`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md)
