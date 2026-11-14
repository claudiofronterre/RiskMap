# Summarize Model Fits

Provides a `summary` method for the "RiskMap" class that computes
standard errors and confidence intervals for likelihood-based model fits

## Usage

``` r
# S3 method for class 'RiskMap'
summary(object, ..., level = 0.95)
```

## Arguments

- object:

  An object of class "RiskMap" from
  [`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)

- ...:

  other parameters.

- level:

  Confidence level for intervals (default 0.95).

## Value

A list of class `"summary.RiskMap"` with parameter estimates, standard
errors, and confidence intervals. For binomial, Poisson and negative
binomial models, its component `mcml` reports certified Monte Carlo
maximum likelihood: `certificate` (an approximate bound, with
probability `level`, on the Monte Carlo error of the estimate in
statistical standard errors, `NA` when there is none), `tolerance`,
`certified`, `outcome`, `correction_refused` (whether the Newton
correction of the last stage did not converge), `stages`,
`fitting_draws` and `validation_draws` (at the last stage),
`certification_bank` (the effective draws and Pareto shape of the
certification draws), `checks_to_review`, the number of checks of
[`check_fit`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md)
that need a look, `starting_values`, where the starting values came
from, either `"user"` (some or all given in `start_pars`, with RiskMap's
defaults for the others), `"current"` (RiskMap's defaults) or
`"transformed_gaussian"` (a Gaussian fit to the transformed data, which
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
uses when the Laplace approximation of the likelihood prefers it to the
defaults), and `starting_reason`, why the defaults were kept other than
by being preferred in that comparison, such as `"re_terms"` when it is
not made for a model with
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
terms, `"transformed_failed"` when the Gaussian fit failed, or
`"not_recorded"` for a fit saved by an earlier version that did not
record why (`NULL` when the comparison chose the start). `mcml` is
`NULL` for Gaussian models. The parameters of the family, such as the
measurement error variance of a Gaussian model, are in
`family_parameters`, a table of those that were estimated, and
`held_family_parameters`, a list of those held at a given value.

## Examples

``` r
data(italy_sim)

fit <- glgpm(
  formula = y ~ gp(),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

summary(fit)
#> ── Linear geostatistical model ─────────────────────────────────────────────────
#> glgpm(formula = y ~ gp(), data = italy_sim[1:100, ], family = "gaussian",
#>     messages = FALSE)
#> The identity link.
#> The limits are those of 95% confidence intervals.
#> 
#> Regression coefficients
#>             Estimate Lower limit Upper limit   StdErr z.value   p.value    
#> (Intercept)  9.32340     7.70078    10.94603  0.82788  11.262 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Measurement error
#>                        Estimate Lower limit Upper limit
#> Measurement error var.  0.17859     0.13335      0.2392
#> 
#> Spatial Gaussian process, Matern correlation with kappa = 0.5
#>                          Estimate Lower limit Upper limit
#> Spatial process var.       4.8494      1.8196      12.924
#> Spatial corr. scale (km)  39.0150      2.7773     548.084
#> Variance of the nugget effect fixed at 0
#> 
#> Log-likelihood   8.618636
#> AIC             -9.237272
summary(fit, level = 0.98)
#> ── Linear geostatistical model ─────────────────────────────────────────────────
#> glgpm(formula = y ~ gp(), data = italy_sim[1:100, ], family = "gaussian",
#>     messages = FALSE)
#> The identity link.
#> The limits are those of 98% confidence intervals.
#> 
#> Regression coefficients
#>             Estimate Lower limit Upper limit   StdErr z.value   p.value    
#> (Intercept)  9.32340     7.39745    11.24935  0.82788  11.262 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Measurement error
#>                        Estimate Lower limit Upper limit
#> Measurement error var.  0.17859     0.12626      0.2526
#> 
#> Spatial Gaussian process, Matern correlation with kappa = 0.5
#>                          Estimate Lower limit Upper limit
#> Spatial process var.       4.8494      1.5150      15.523
#> Spatial corr. scale (km)  39.0150      1.6947     898.204
#> Variance of the nugget effect fixed at 0
#> 
#> Log-likelihood   8.618636
#> AIC             -9.237272
```
