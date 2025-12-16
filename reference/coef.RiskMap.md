# Extract Parameter Estimates from a "RiskMap" Model Fit

This `coef` method for the "RiskMap" class extracts the maximum
likelihood estimates from model fits obtained from the
[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
function.

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

A list containing the maximum likelihood estimates:

- beta:

  A vector of coefficient estimates.

- sigma2:

  The estimate for the variance parameter \\\sigma^2\\.

- phi:

  The estimate for the spatial range parameter \\\phi\\.

- tau2:

  The estimate for the nugget effect parameter \\\tau^2\\, if
  applicable.

- sigma2_me:

  The estimate for the measurement error variance \\\sigma^2\_{me}\\, if
  applicable.

- sigma2_re:

  A vector of variance estimates for the random effects, if applicable.

## Details

The function processes the `RiskMap` object to extract and name the
estimated parameters appropriately, transforming them if necessary.

## Note

This function handles both Gaussian and non-Gaussian families, and
accounts for fixed and random effects in the model.

## See also

[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
