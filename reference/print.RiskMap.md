# Print a Fitted RiskMap Model

Print method for objects of class `"RiskMap"` from
[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md):
the call, the estimates of the regression coefficients and of the
covariance parameters on their natural scale, the log-likelihood and,
for binomial, Poisson and negative binomial models, whether the Monte
Carlo error of the estimate is certified. Standard errors and confidence
intervals are given by
[`summary.RiskMap`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md),
and the checks on the computation by
[`check_fit`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md).

## Usage

``` r
# S3 method for class 'RiskMap'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"RiskMap"`.

- digits:

  The number of significant digits of the estimates.

- ...:

  other parameters.

## Value

Invisibly returns `x`.

## Examples

``` r
data(italy_sim)

fit <- glgpm(
  formula = y ~ gp(),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

fit
#> ── Linear geostatistical model ─────────────────────────────────────────────────
#> glgpm(formula = y ~ gp(), data = italy_sim[1:100, ], family = "gaussian",
#>     messages = FALSE)
#> 100 observations at 10 locations.
#> 
#> Regression coefficients
#>             Estimate
#> (Intercept)    9.323
#> 
#> Covariance parameters (Matern, kappa = 0.5)
#>                          Estimate
#> Spatial process var.        4.849
#> Spatial corr. scale (km)    39.01
#> Measurement error var.     0.1786
#> Variance of the nugget effect fixed at 0
#> 
#> Log-likelihood  8.618636
#> 
#> ℹ `summary()` gives standard errors and confidence intervals.
```
