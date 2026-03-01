# Simulation from Decay-adjusted Spatio-temporal (DAST) models

Simulates data from a fitted DAST model object (output from `dast`) or
from user-specified DAST parameters.

## Usage

``` r
dast_sim(
  n_sim,
  model_fit = NULL,
  formula = NULL,
  data = NULL,
  den = NULL,
  time = NULL,
  mda_times = NULL,
  int_mat = NULL,
  power_val = NULL,
  cov_offset = NULL,
  crs = NULL,
  convert_to_crs = NULL,
  scale_to_km = TRUE,
  sim_pars = list(beta = NULL, sigma2 = NULL, tau2 = NULL, phi = NULL, psi = NULL,
    sigma2_re = NULL, alpha = NULL, gamma = NULL),
  messages = TRUE
)
```

## Arguments

- n_sim:

  Number of simulations.

- model_fit:

  Optional fitted DAST model object of class `RiskMap`. If supplied, it
  overrides model specification arguments.

- formula:

  Model formula including a
  [`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md)
  term.

- data:

  Data frame or `sf` object used for simulation.

- den:

  Binomial denominator variable. If missing, it is assumed to be 1.

- time:

  Survey-time information. For simulations from scratch this can be a
  column in `data` (unquoted name or character string) or a numeric
  vector of length `nrow(data)`.

- mda_times:

  Vector of MDA times.

- int_mat:

  Intervention matrix (n x length(mda_times)) with coverage values.

- power_val:

  Power value for the MDA impact function.

- cov_offset:

  Optional offset for the linear predictor.

- crs:

  Coordinate reference system (CRS) code.

- convert_to_crs:

  Optional CRS to transform coordinates to before simulation.

- scale_to_km:

  Logical; if TRUE distances are computed in kilometers.

- sim_pars:

  List of simulation parameters. Used only when `model_fit` is `NULL`.
  Includes `beta`, `sigma2`, `tau2`, `phi`, `psi`, `sigma2_re`, `alpha`,
  and `gamma`.

- messages:

  Logical; if TRUE print progress messages.

## Value

A list with simulated outcomes in `data_sim` as an `n x n_sim` matrix
(rows are observations, columns are simulations), simulated latent
components and parameter values used for simulation.
