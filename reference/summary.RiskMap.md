# Summarize Model Fits

Provides a `summary` method for the "RiskMap" class that computes the
standard errors and p-values for likelihood-based model fits.

## Usage

``` r
# S3 method for class 'RiskMap'
summary(object, ..., conf_level = 0.95)
```

## Arguments

- object:

  An object of class "RiskMap" obtained as a result of a call to
  [`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

- ...:

  other parameters.

- conf_level:

  The confidence level for the intervals (default is 0.95).

## Value

A list containing:

- reg_coef:

  A matrix with the estimates, standard errors, z-values, p-values, and
  confidence intervals for the regression coefficients.

- me:

  A matrix with the estimates and confidence intervals for the
  measurement error variance, if applicable.

- sp:

  A matrix with the estimates and confidence intervals for the spatial
  process parameters.

- tau2:

  The fixed nugget variance, if applicable.

- ranef:

  A matrix with the estimates and confidence intervals for the random
  effects variances, if applicable.

- conf_level:

  The confidence level used for the intervals.

- family:

  The family of the model (e.g., "gaussian").

- kappa:

  The kappa parameter of the model.

- log.lik:

  The log-likelihood of the model fit.

- cov_offset_used:

  A logical indicating if a covariance offset was used.

- aic:

  The Akaike Information Criterion (AIC) for the model, if applicable.

## Details

This function computes the standard errors and p-values for the
parameters of a "RiskMap" model, adjusting for the covariance structure
if needed.

## Note

Handles both Gaussian and non-Gaussian families, and accounts for fixed
and random effects in the model.

## See also

[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`coef.RiskMap`](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
