# Fitting of decay-adjusted spatio-temporal (DAST) model

The function fits a decay-adjusted spatio-temporal (DAST) model using
Monte Carlo maximum likelihood. The DAST model allows for the
incorporation of temporal decay in disease prevalence due to the impact
of mass drug administration (MDA). The function requires the full MDA
history as detailed in the arguments below.

Spatial and spatio-temporal dependence is specified through the
[`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md)
term in the model formula:

- `gp(x, y)` fits a purely spatial Gaussian process.

- `gp(x, y, t_gp)` fits a spatio-temporal Gaussian process, where `t_gp`
  is used as the GP temporal index.

In all cases, the `time` argument must be specified separately and
provides the observation-level survey times used in modelling MDA
impact. These survey times may differ from the GP temporal index.

## Usage

``` r
dast(
  formula,
  data,
  den = NULL,
  time,
  mda_times,
  int_mat,
  penalty = NULL,
  drop = NULL,
  power_val,
  crs = NULL,
  convert_to_crs = NULL,
  scale_to_km = TRUE,
  control_mcmc = set_control_sim(),
  par0 = NULL,
  S_samples = NULL,
  return_samples = TRUE,
  messages = TRUE,
  start_pars = list(beta = NULL, sigma2 = NULL, tau2 = NULL, phi = NULL, psi = NULL,
    sigma2_re = NULL, gamma = NULL, alpha = NULL),
  alpha_start = 0.5,
  gamma_start = 2.5
)
```

## Arguments

- formula:

  A model formula specifying the response variable, predictors, and the
  GP structure through
  [`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md).

- data:

  A `data.frame` or `sf` object containing the dataset.

- den:

  The denominator for binomial models.

- time:

  A variable in `data` giving the survey times of observations
  (required).

- mda_times:

  A vector specifying the mass drug administration (MDA) times.

- int_mat:

  Intervention matrix specifying the timing and coverage of MDA; the
  dimension of the matrix must be `n * n_mda`, where `n` is the number
  of rows of `data` and `n_mda` is the length of `mda_times`.

- penalty:

  Optional list specifying penalty functions for regularization, used in
  the estimation of the "drop" parameter `alpha`.

- drop:

  Optional value used for fixing the "drop" parameter of the MDA impact
  function.

- power_val:

  Value expressing the power of the MDA impact function.

- crs:

  Optional coordinate reference system (CRS) for spatial data.

- convert_to_crs:

  CRS to which spatial data should be converted.

- scale_to_km:

  Logical; whether to scale distances to kilometers (default: `TRUE`).

- control_mcmc:

  A list of MCMC control parameters, typically from
  [`set_control_sim()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_sim.md).

- par0:

  Optional list of initial parameter values.

- S_samples:

  Number of posterior samples to retain.

- return_samples:

  Logical; whether to return posterior samples (default: `TRUE`).

- messages:

  Logical; whether to print messages (default: `TRUE`).

- start_pars:

  List of starting values for parameters.

## Value

A list containing model estimates, posterior samples, and metadata,
including:

- `y`: Response variable values.

- `D`: Covariate matrix.

- `coords`: Unique spatial coordinates.

- `mda_times`: MDA time points.

- `survey_times_data`: Survey time data from the `time` argument.

- `time`: GP temporal index if specified in `gp(x,y,t_gp)`.

- `int_mat`: Intervention matrix.

- `ID_coords`: Indices of spatial locations (and time if spatio-temporal
  GP).

- `re`: Random effects levels (if applicable).

- `ID_re`: Indices of random effects (if applicable).

- `power_val`: Power of the MDA impact function.

- `fix_tau2`: Fixed tau-squared value (if applicable).

- `fix_alpha`: Fixed alpha value (if applicable).

- `formula`: Model formula.

- `crs`: Coordinate reference system.

- `scale_to_km`: Indicator of distance scaling.

- `data_sf`: Processed spatial dataset.

- `family`: Model family (e.g., "binomial").

- `sst`: Logical indicator of whether a spatio-temporal GP was used.

- `kappa`: Smoothness parameter.

- `units_m`: Denominator for binomial models.

- `cov_offset`: Offset for covariates.

- `call`: Function call.

- `penalty`: Penalty function details (if applicable).

- `posterior_samples`: Posterior samples if `return_samples = TRUE`.

## See also

[`set_control_sim`](https://claudiofronterre.github.io/RiskMap/reference/set_control_sim.md),
[`summary.RiskMap`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md),
[`to_table`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
