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

The observed value is unusual for the model when it is outside the
interval.

The proportion at or above the observed value is not a p-value. The
parameters are estimated from the same data, so that it is closer to 0.5
than a p-value would be, and a statistic that the model was fitted to
can hardly have an extreme one. It counts the simulated values equal to
the observed one, and where there are many of them it can be 1, or near
1, for a value that isn't unusual. See
[`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md).

The printed summary also says in how many distance classes the observed
variogram is outside its joint envelope, where the observed distribution
is, and, for binomial models, whether the prevalences are those of
observations or of locations.

## See also

[`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md),
and
[`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)
for a table to render.
