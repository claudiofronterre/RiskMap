# Set Control Parameters for Simulation

This function sets control parameters for running simulations,
supporting MCMC methods for glgpm models.

## Usage

``` r
set_control_mcmc(
  n_sim = 22000,
  burnin = 2000,
  thin = 10,
  h = NULL,
  c1.h = 1,
  c2.h = 0.6,
  linear_model = FALSE
)
```

## Arguments

- n_sim:

  Integer. The total number of simulations to run. Default is 22000,
  which with the default burn-in and thinning keeps 2,000 draws. With
  `linear_model = TRUE` it is the number of independent draws, 2,000 by
  default.

- burnin:

  Integer. The number of initial simulations to discard (burn-in/warmup
  period). Default is 2000.

- thin:

  Integer. The interval at which simulations are recorded (thinning
  interval, MCMC only). Default is 10.

- h:

  Numeric. An optional positive step size for Langevin MCMC.

- c1.h:

  Numeric. Positive scale for the burn-in Robbins-Monro step-size
  adaptation. Default is 1.

- c2.h:

  Numeric. Exponent controlling the diminishing adaptation gain. Must be
  larger than 0.5 and no larger than 1. Default is 0.6.

- linear_model:

  Logical. If TRUE, sets up parameters for a linear model. Default is
  FALSE.

## Value

A list of control parameters with class "RiskMap_control_mcmc". Contents
depend on `sampler`:

- For "mcmc": n_sim, burnin, thin, h, c1.h, c2.h, linear_model

## Details

If `linear_model = TRUE`, only `n_sim` is required.

During burn-in the step size of the Langevin chain is adapted towards an
acceptance rate of 0.574, starting from `h` (by default \\1.65 /
d^{1/6}\\ for \\d\\ latent effects). The step kept after burn-in is the
geometric mean of the adapted step over the second half of burn-in.

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

# Custom MCMC parameters
control_mcmc <- set_control_mcmc(n_sim = 15000, burnin = 3000, thin = 20)
```
