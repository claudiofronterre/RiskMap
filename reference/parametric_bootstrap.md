# Parametric bootstrap of a fitted geostatistical model

Simulates datasets from a model fitted with
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
at its estimated parameters, refits the model to each, and summarises
the sampling distribution of the estimates with standard errors, bias
and two kinds of bootstrap confidence interval.

## Usage

``` r
parametric_bootstrap(
  object,
  n_sim = 1000,
  seed = NULL,
  level = 0.95,
  cores = 1L,
  messages = getOption("RiskMap.messages", TRUE)
)
```

## Arguments

- object:

  A fitted `RiskMap` object returned by
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
  from a Gaussian, binomial, Poisson or negative binomial model, with or
  without a nugget and
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  terms.

- n_sim:

  Number of bootstrap replicates. Percentile intervals at the 95% level
  need about 1000 for stable limits.

- seed:

  Optional non-negative integer that makes the bootstrap reproducible.
  When `NULL` (the default), one integer is drawn from the session's
  random number generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  also reproduces it. The seed used is returned in `settings`. The
  caller's random number state is otherwise left unchanged.

- level:

  Probability covered by the intervals.

- cores:

  Number of processes. With more than one, the replicates run in
  parallel by forking, which is not available on Windows, where they run
  one after another.

- messages:

  Logical, whether to report on the bootstrap, with a progress bar as
  described in
  [riskmap_progress](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md).
  Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.

## Value

An object of class `RiskMap_bootstrap`, a list with

- estimates:

  A matrix of the bootstrap estimates on their natural scale, with one
  row per replicate (`NA` for failed refits) and one column per
  estimated parameter, named as `unlist(coef(object))` names them.
  Parameters held by `fixed_pars` are held in the refits and have no
  column.

- estimate:

  The estimate of `object`, on the same scale.

- labels:

  Readable names of the parameters.

- failures:

  A data frame with the `replicate`, `seed` and error `message` of each
  failed refit.

- warnings:

  A data frame with the `replicate`, `seed` and `message` of each
  warning from the refits.

- settings:

  A list with `n_sim`, `seed`, `level`, `cores`, the `seed_spacing`
  between replicates, the `family` and the `distance_units` of the
  model, and the `elapsed_seconds`.

## Details

Each replicate simulates a response at the original locations, with the
original covariates, denominators, offsets and random-effect groups,
from the fitted model (see
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)),
and refits the model to it by
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
starting from the estimate, without the final conditional samples and
quietly. The refits reuse the design matrix, the locations and the
settings of `control_mcml` stored in the fit, so the model formula and
the data are not evaluated again.

Replicate \\i\\ uses one seed, \\s_i = s + 10000 (i - 1)\\, where \\s\\
is `seed`. It simulates its data with `simulate_glgpm(seed =` \\s_i\\`)`
and refits with `glgpm(seed =` \\s_i\\`)`. The Monte Carlo streams of a
binomial, Poisson or negative binomial fit are derived from its seed by
offsets between 1 and about 6000, so consecutive integers would make
replicates share streams, and seeds 10000 apart do not. A replicate
depends only on its seed, so the results are identical for any number of
`cores` and the first replicates of a larger bootstrap are those of a
smaller one. Each replicate runs with one BLAS and OpenMP thread when
the RhpcBLASctl package is installed, and the user's thread settings are
restored afterwards.

A refit that fails is recorded in `failures` with its error message, and
its estimates are `NA`. Warnings from the refits, for example about
variances estimated at zero, are recorded in `warnings`.

[`summary()`](https://rdrr.io/r/base/summary.html) reports, for each
parameter, the estimate, the median of the bootstrap estimates, the bias
(their mean minus the estimate), their standard deviation as the
standard error, the proportion below the estimate, and two intervals.
The percentile interval takes the quantiles of the bootstrap estimates
at \\\alpha/2\\ and \\1 - \alpha/2\\, where \\\alpha = 1 -\\ `level`.
The bias-corrected (BC) interval takes them at \\\Phi(2 z_0 +
z\_{\alpha/2})\\ and \\\Phi(2 z_0 + z\_{1 - \alpha/2})\\, where \\z_0 =
\Phi^{-1}(p)\\ and \\p\\ is the proportion of bootstrap estimates below
the estimate (Efron and Tibshirani, 1993, Chapter 14). The BC interval
is not available when \\p\\ is 0 or 1.

The difference matters for the spatial variance \\\sigma^2\\ and the
scale \\\phi\\. Their maximum likelihood estimators are biased downward,
so most bootstrap estimates fall below the estimate, just as the
estimate tends to fall below the true value. The percentile interval
then shifts downward, the wrong way, because it reproduces the bias of
the estimator instead of correcting it. With \\p \> 0.5\\, \\z_0 \> 0\\
and the BC interval moves both limits upward. Both intervals are
unchanged by monotone transformations, so they are the same whether
computed for \\\phi\\ or \\\log \phi\\.

## References

Efron, B. and Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall, New York.

## See also

The methods for a bootstrap are
[print()](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_bootstrap.md),
[summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_bootstrap.md)
and
[`plot_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/plot_bootstrap.md),
which [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the
result also calls.
[`profile.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md)
gives profile intervals, and
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
simulates from a fit.

## Examples

``` r
data(italy_sim)
fit <- glgpm(y ~ gp(), data = italy_sim[1:50, ], family = "gaussian",
             messages = FALSE)
bootstrap <- parametric_bootstrap(fit, n_sim = 10, seed = 1,
                                  messages = FALSE)
bootstrap
#> ── Parametric bootstrap ────────────────────────────────────────────────────────
#> Linear geostatistical model, 10 replicates (seed 1). Variances and the scale
#> are on their natural scale.
#> 
#> Bootstrap distribution
#>                 Parameter Estimate Median     Bias Std. error Below
#>               (Intercept)    9.986  10.04  0.06578     0.6331  0.50
#>      Spatial process var.    2.554  1.443   -1.216     0.6673  1.00
#>  Spatial corr. scale (km)    44.47  8.297    8.477      105.8  0.80
#>    Measurement error var.   0.2188 0.2093 -0.01434     0.0374  0.50
#> ℹ Bias is the mean of the bootstrap estimates minus the estimate, and Below the
#>   proportion of them below the estimate.
#> 
#> 95% intervals, percentile and bias-corrected (BC)
#>                 Parameter Percentile lower Percentile upper BC lower BC upper
#>               (Intercept)            8.982             10.9    8.982     10.9
#>      Spatial process var.           0.3294             2.23       NA       NA
#>  Spatial corr. scale (km)            3.181            286.1    3.654    344.9
#>    Measurement error var.           0.1593           0.2499   0.1593   0.2499
#> ℹ NA means that every bootstrap estimate is on one side of the estimate, so the
#>   BC interval is not available.

# \donttest{
data(liberia)
liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
fit <- glgpm(npos ~ log(elevation) + gp(), denominator = ntest,
             data = liberia, distance_crs = 32629, family = "binomial",
             seed = 2026)
#> ℹ The model is fitted in EPSG:32629, with distances in km.
#> ✔ Certified. The Monte Carlo error of the estimate is within 0.0091 statistical
#>   standard errors (an approximate 95% bound, tolerance 0.1).
# Ten simulations keep the example short. An analysis needs many more.
bootstrap <- parametric_bootstrap(fit, n_sim = 10, seed = 2026)
#> ℹ Simulating and refitting 10 replicates.
#> ✔ Done in 18 seconds. 10 of 10 refits succeeded.
summary(bootstrap)
#> ── Parametric bootstrap ────────────────────────────────────────────────────────
#> Binomial geostatistical model, 10 replicates (seed 2026). Variances and the
#> scale are on their natural scale.
#> 
#> Bootstrap distribution
#>                 Parameter Estimate  Median     Bias Std. error Below
#>               (Intercept)   -3.196  -3.052  0.08698     0.2914  0.30
#>            log(elevation)   0.3375  0.3135 -0.01789    0.05728  0.80
#>      Spatial process var.  0.08825 0.09345 0.006961     0.0455  0.40
#>  Spatial corr. scale (km)    45.67   42.09    3.081      24.98  0.60
#> ℹ Bias is the mean of the bootstrap estimates minus the estimate, and Below the
#>   proportion of them below the estimate.
#> 
#> 95% intervals, percentile and bias-corrected (BC)
#>                 Parameter Percentile lower Percentile upper BC lower BC upper
#>               (Intercept)            -3.57            -2.69   -3.578   -2.901
#>            log(elevation)           0.2481           0.4226   0.3021    0.429
#>      Spatial process var.          0.03496           0.1661  0.03205   0.1632
#>  Spatial corr. scale (km)            20.32            93.91    22.38    95.55
plot(bootstrap)

# }
```
