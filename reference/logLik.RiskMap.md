# Log-Likelihood, Covariance and Number of Observations of a Fit

The [`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) methods for fits from
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
with which [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work too.

## Usage

``` r
# S3 method for class 'RiskMap'
logLik(object, ...)

# S3 method for class 'RiskMap'
vcov(object, ...)

# S3 method for class 'RiskMap'
nobs(object, ...)
```

## Arguments

- object:

  A fit from
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

[`logLik()`](https://rdrr.io/r/stats/logLik.html) gives an object of
class `logLik`, [`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix
and [`nobs()`](https://rdrr.io/r/stats/nobs.html) an integer.

## Details

[`logLik()`](https://rdrr.io/r/stats/logLik.html) gives the
log-likelihood at the estimates. It is exact for a Gaussian model. For a
binomial, Poisson or negative binomial model it is a Monte Carlo
estimate, so that differences between such models that are not well
above its Monte Carlo error mean nothing, which is why
[summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)
gives no AIC for them. That of a negative binomial model leaves out the
terms that do not depend on the parameters, as that of a Poisson model
does (see
[`negative_binomial()`](https://claudiofronterre.github.io/RiskMap/reference/negative_binomial.md)).
Its degrees of freedom, `df`, are the number of estimated parameters,
without those that `fixed_pars` holds.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) gives the covariance
matrix of the estimates on the scale of
[coef()](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md),
by the delta method from their covariance on the scale on which the
model is fitted, where the variances and the scale are on the log scale.
Its rows and columns are named as `unlist(coef(object))` names them, for
the estimated parameters only. A parameter without a standard error,
such as a variance estimated at zero, has a row and a column of `NA`,
and a fit that holds every parameter has an empty matrix.

[`nobs()`](https://rdrr.io/r/stats/nobs.html) gives the number of
observations, the rows of the data, whatever their denominators.

## See also

[`summary.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)
for standard errors and confidence intervals.

## Examples

``` r
data(italy_sim)

fit <- glgpm(
  formula = y ~ gp(),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

logLik(fit)
#> 'log Lik.' 8.618636 (df=4)
AIC(fit)
#> [1] -9.237272
vcov(fit)
#>                  beta.(Intercept)        sigma2           phi     sigma2_me
#> beta.(Intercept)     6.853929e-01 -1.185935e-01 -5.836712e+00  3.579153e-07
#> sigma2              -1.185935e-01  5.882057e+00  5.605035e+01 -7.744215e-05
#> phi                 -5.836712e+00  5.605035e+01  2.766884e+03  5.171587e-04
#> sigma2_me            3.579153e-07 -7.744215e-05  5.171587e-04  7.087906e-04
nobs(fit)
#> [1] 100
```
