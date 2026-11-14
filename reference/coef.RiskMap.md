# Extract Parameter Estimates from a "RiskMap" Model Fit

This `coef` method for the "RiskMap" class extracts the maximum
likelihood estimates from model fits obtained from the
[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
functions.

## Usage

``` r
# S3 method for class 'RiskMap'
coef(object, ...)
```

## Arguments

- object:

  An object of class "RiskMap" obtained as a result of a call to
  [`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

- ...:

  other parameters.

## Value

A list containing the maximum likelihood estimates. For standard models:

- beta:

  A vector of coefficient estimates.

- sigma2:

  The estimate for the variance parameter \\\sigma^2\\.

- phi:

  The estimate for the scale parameter \\\phi\\ of the spatial
  correlation function.

- tau2:

  The estimate for the nugget effect parameter \\\tau^2\\, if
  applicable.

- sigma2_me:

  The estimate for the measurement error variance \\\sigma^2\_{me}\\, if
  applicable.

- sigma2_re:

  A named vector of variance estimates for the random effects, if
  applicable.

- beta:

  Coefficient estimates for log mean worm burden.

- k:

  Negative binomial overdispersion parameter.

- rho:

  Egg detection rate (fecundity).

- sigma2:

  Spatial process variance.

- phi:

  Spatial correlation scale.

## See also

[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
