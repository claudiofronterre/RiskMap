# Summarise the assessment of models against simulated data

Gives, for each model and statistic, the observed value, the median and
the central 95% interval of the simulated values, and the proportion of
simulated values at or above the observed one.

## Usage

``` r
# S3 method for class 'RiskMap_predictive_check'
summary(object, ...)
```

## Arguments

- object:

  The result of
  [`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md).

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

A data frame of class `summary.RiskMap_predictive_check` with one row
for each statistic and model, and the columns `model`, `statistic`,
`observed`, `median`, `lower`, `upper`, `at_or_above`, the proportion of
simulated values at or above the observed one, and `n_not_finite`, the
number of simulated datasets in which the statistic isn't finite, which
are left out.

## Details

The proportion at or above the observed value is not a p-value. The
parameters are estimated from the same data, so that it is closer to 0.5
than a p-value would be, and a statistic that the model was fitted to
can hardly have an extreme one. See
[`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md).

The printed summary also says in how many distance classes the observed
variogram is outside its joint envelope, and where the observed
distribution is.

## See also

[`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md),
and
[`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)
for a table to render.
