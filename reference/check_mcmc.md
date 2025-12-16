# Check MCMC Convergence for Spatial Random Effects

This function checks the Markov Chain Monte Carlo (MCMC) convergence of
spatial random effects for either a `RiskMap` or `RiskMap.pred.re`
object. It plots the trace plot and autocorrelation function (ACF) for
the MCMC chain and calculates the effective sample size (ESS).

## Usage

``` r
check_mcmc(object, check_mean = TRUE, component = NULL, ...)
```

## Arguments

- object:

  An object of class `RiskMap` or `RiskMap.pred.re`. `RiskMap` is the
  output from
  [`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  function, and `RiskMap.pred.re` is obtained from the
  [`pred_over_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_over_grid.md)
  function.

- check_mean:

  Logical. If `TRUE`, checks the MCMC chain for the mean of the spatial
  random effects. If `FALSE`, checks the chain for a specific component
  of the random effects vector.

- component:

  Integer. The index of the spatial random effects component to check
  when `check_mean = FALSE`. Must be a positive integer corresponding to
  a location in the data. Ignored if `check_mean = TRUE`.

- ...:

  Additional arguments passed to the
  [`acf`](https://rdrr.io/r/stats/acf.html) function for customizing the
  ACF plot.

## Value

No return value, called for side effects (plots and warnings).

## Details

The function first checks that the input object is either of class
`RiskMap` or `RiskMap.pred.re`. Depending on the value of `check_mean`,
it either calculates the mean of the spatial random effects across all
locations for each iteration or uses the specified component. It then
generates two plots: - A trace plot of the selected spatial random
effect over iterations. - An autocorrelation plot (ACF) with the
effective sample size (ESS) displayed in the title.

The ESS is computed using the
[`ess`](https://rdrr.io/pkg/sns/man/ess.html) function, which provides a
measure of the effective number of independent samples in the MCMC
chain.

If `check_mean = TRUE`, the `component` argument is ignored, and a
warning is issued. To specify a particular component of the random
effects vector, set `check_mean = FALSE` and provide a valid `component`
value.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>
