# Summarise a parametric bootstrap

Computes, for each parameter, the bootstrap median, bias, standard error
and the percentile and bias-corrected intervals described in
[`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md).

## Usage

``` r
# S3 method for class 'RiskMap_bootstrap'
summary(object, level = object$settings$level, ...)
```

## Arguments

- object:

  A `RiskMap_bootstrap` object from
  [`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md).

- level:

  Probability covered by the intervals. Defaults to the level of the
  bootstrap.

- ...:

  Not used.

## Value

An object of class `summary.RiskMap_bootstrap`, a list with the table of
results (`table`), the `level`, the numbers of replicates and of failed
refits (`n_sim`, `n_failed`) and the `seed`.
