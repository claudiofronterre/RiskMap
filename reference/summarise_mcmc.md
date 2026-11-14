# Effective sample sizes of the MCMC samples of the random effects

Computes the effective sample size of the Markov chain Monte Carlo
(MCMC) samples of the random effects stored in a fitted binomial,
Poisson or negative binomial model, for the mean of the spatial random
effects over the locations, for the spatial random effect at each
location and for each level of every unstructured random effect.

## Usage

``` r
summarise_mcmc(object)
```

## Arguments

- object:

  A fitted `RiskMap` object returned by
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  for a binomial, Poisson or negative binomial model with
  `return_samples = TRUE`.

## Value

A tibble with one row per chain and columns

- effect:

  `"mean of S"` for the mean of the spatial random effects over the
  locations, `"S"` for the spatial random effect at a location, or the
  name of the variable of an unstructured random effect.

- component:

  The location, numbered as the rows of `object$coords`, or the level of
  the random effect, numbered as in `object$re`. It is `NA` for the
  mean.

- ess:

  The effective sample size.

The number of samples is stored in the attribute `"n_samples"`.

## Details

The samples are those drawn by
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
at the fitted parameter values, which are used for the fitted values and
for spatial prediction. When the model has a nugget, the spatial random
effect at a location includes it.

The effective sample size of a chain is the number of independent draws
that would estimate the mean of its target distribution as precisely as
the correlated MCMC samples do. It is the number of samples divided by
the integrated autocorrelation time, \\1 + 2 \sum\_{k \ge 1} \rho_k\\,
where \\\rho_k\\ is the autocorrelation at lag \\k\\. The sum is
estimated with the initial monotone sequence estimator of Geyer (1992),
which adds pairs of consecutive sample autocorrelations while their sum
is positive and forces these sums to decrease. As in Stan, the
autocorrelation time is not allowed to fall below \\1 / \log\_{10} N\\
for \\N\\ samples, so the effective sample size of an antithetic chain
is at most \\N \log\_{10} N\\. A chain that never moves has an effective
sample size of zero.

## References

Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statistical
Science*, 7(4), 473-483.

## See also

[`plot_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/plot_mcmc.md)

## Examples

``` r
# \donttest{
data(liberia)
liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
fit <- glgpm(npos ~ log(elevation) + gp(), denominator = ntest,
             data = liberia, model_crs = 32629, family = "binomial",
             messages = FALSE)
ess <- summarise_mcmc(fit)
# The mean and the five locations with the smallest effective sample size
ess[1, ]
#> # A tibble: 1 × 3
#>   effect    component   ess
#>   <chr>         <int> <dbl>
#> 1 mean of S        NA 1294.
head(ess[order(ess$ess), ], 5)
#> # A tibble: 5 × 3
#>   effect    component   ess
#>   <chr>         <int> <dbl>
#> 1 mean of S        NA 1294.
#> 2 S                44 1325.
#> 3 S                28 1337.
#> 4 S                56 1338.
#> 5 S                53 1352.
# }
```
