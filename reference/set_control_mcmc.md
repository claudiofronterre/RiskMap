# Set Control Parameters for Simulation

The settings of the Langevin sampler that draws the latent effects of a
binomial, Poisson or negative binomial model given the data, in
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
and
[predict()](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md).

## Usage

``` r
set_control_mcmc(
  n_draws = 2000,
  burnin = 2000,
  thin = 10,
  h = NULL,
  c1_h = 1,
  c2_h = 0.6
)
```

## Arguments

- n_draws:

  Positive integer. The number of draws kept, after burn-in and
  thinning, on which fitted values and predictions are based. Default is
  2,000. The chain runs `burnin + n_draws * thin` iterations.

- burnin:

  Non-negative integer. The number of iterations discarded at the start
  of the chain, during which the step size is adapted. Default is 2,000.

- thin:

  Positive integer. One iteration in `thin` is kept after burn-in.
  Default is 10.

- h:

  Numeric. An optional positive step size for Langevin MCMC.

- c1_h:

  Numeric. Positive scale for the burn-in Robbins-Monro step-size
  adaptation. Default is 1.

- c2_h:

  Numeric. Exponent controlling the diminishing adaptation gain. Must be
  larger than 0.5 and no larger than 1. Default is 0.6.

## Value

A list of class `"RiskMap_control_mcmc"` with `n_draws`, `burnin`,
`thin`, `h`, `c1_h` and `c2_h`.

## Details

During burn-in the step size of the Langevin chain is adapted towards an
acceptance rate of 0.574, starting from `h` (by default \\1.65 /
d^{1/6}\\ for \\d\\ latent effects). The step kept after burn-in is the
geometric mean of the adapted step over the second half of burn-in.

The latent effects of a Gaussian model are drawn independently, so its
predictions take `n_draws` draws without burn-in or thinning.

Seeds are not part of the control settings. Pass `seed` to
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md),
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
or
[`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md).

## See also

[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)

## Examples

``` r
# Default parameters (MCMC)
control_mcmc <- set_control_mcmc()

# Custom MCMC parameters: 600 draws from 3,000 + 600 * 20 iterations
control_mcmc <- set_control_mcmc(n_draws = 600, burnin = 3000, thin = 20)
```
