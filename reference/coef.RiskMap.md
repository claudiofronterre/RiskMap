# Extract Parameter Estimates from a "RiskMap" Model Fit

This `coef` method for the "RiskMap" class extracts the maximum
likelihood estimates from model fits obtained from the
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
functions.

## Usage

``` r
# S3 method for class 'RiskMap'
coef(object, ...)
```

## Arguments

- object:

  An object of class "RiskMap" obtained as a result of a call to
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

- ...:

  other parameters.

## Value

A named list of every parameter of the model on its natural scale, the
estimates and the values held by `fixed_pars`, in the order of
[riskmap_parameters](https://claudiofronterre.github.io/RiskMap/reference/riskmap_parameters.md):

- beta:

  The regression coefficients, named by the columns of the design
  matrix.

- sigma2:

  The variance of the spatial process.

- phi:

  The scale of the spatial correlation, in the model's distance units.

- tau2:

  The variance of the nugget, for a model with one.

- sigma2_me:

  The measurement error variance, for a Gaussian model.

- theta:

  The size, for a negative binomial model.

- sigma2_re:

  The variances of the random effects, named by their terms, for a model
  with
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  terms.

The list can be given as `start_pars` of
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

## See also

[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[riskmap_parameters](https://claudiofronterre.github.io/RiskMap/reference/riskmap_parameters.md)
