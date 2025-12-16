# Estimation of Generalized Linear Gaussian Process Models

Fits generalized linear Gaussian process models to spatial data,
incorporating spatial Gaussian processes with a Matern correlation
function. Supports Gaussian, binomial, and Poisson response families.

## Usage

``` r
glgpm(
  formula,
  data,
  family,
  invlink = NULL,
  den = NULL,
  crs = NULL,
  convert_to_crs = NULL,
  scale_to_km = TRUE,
  control_mcmc = set_control_sim(),
  par0 = NULL,
  S_samples = NULL,
  return_samples = TRUE,
  messages = TRUE,
  fix_var_me = NULL,
  start_pars = list(beta = NULL, sigma2 = NULL, tau2 = NULL, phi = NULL, sigma2_me =
    NULL, sigma2_re = NULL)
)
```

## Arguments

- formula:

  A formula object specifying the model to be fitted. The formula should
  include fixed effects, random effects (specified using
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)),
  and spatial effects (specified using
  [`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md)).

- data:

  A data frame or sf object containing the variables in the model.

- family:

  A character string specifying the distribution of the response
  variable. Must be one of "gaussian", "binomial", or "poisson".

- invlink:

  A function that defines the inverse of the link function for the
  distribution of the data given the random effects.

- den:

  Optional offset for binomial or Poisson distributions. If not
  provided, defaults to 1 for binomial.

- crs:

  Optional integer specifying the Coordinate Reference System (CRS) if
  data is not an sf object. Defaults to 4326 (long/lat).

- convert_to_crs:

  Optional integer specifying a CRS to convert the spatial coordinates.

- scale_to_km:

  Logical indicating whether to scale coordinates to kilometers.
  Defaults to TRUE.

- control_mcmc:

  Control parameters for MCMC sampling. Must be an object of class
  "mcmc.RiskMap" as returned by
  [`set_control_sim`](https://claudiofronterre.github.io/RiskMap/reference/set_control_sim.md).

- par0:

  Optional list of initial parameter values for the MCMC algorithm.

- S_samples:

  Optional matrix of pre-specified sample paths for the spatial random
  effect.

- return_samples:

  Logical indicating whether to return MCMC samples when fitting a
  Binomial or Poisson model. Defaults to FALSE.

- messages:

  Logical indicating whether to print progress messages. Defaults to
  TRUE.

- fix_var_me:

  Optional fixed value for the measurement error variance.

- start_pars:

  Optional list of starting values for model parameters: beta
  (regression coefficients), sigma2 (spatial process variance), tau2
  (nugget effect variance), phi (spatial correlation scale), sigma2_me
  (measurement error variance), and sigma2_re (random effects
  variances).

## Value

An object of class "RiskMap" containing the fitted model and relevant
information:

- y:

  Response variable.

- D:

  Covariate matrix.

- coords:

  Unique spatial coordinates.

- ID_coords:

  Index of coordinates.

- re:

  Random effects.

- ID_re:

  Index of random effects.

- fix_tau2:

  Fixed nugget effect variance.

- fix_var_me:

  Fixed measurement error variance.

- formula:

  Model formula.

- family:

  Response family.

- crs:

  Coordinate Reference System.

- scale_to_km:

  Indicator if coordinates are scaled to kilometers.

- data_sf:

  Original data as an sf object.

- kappa:

  Spatial correlation parameter.

- units_m:

  Distribution offset for binomial/Poisson.

- cov_offset:

  Covariate offset.

- call:

  Matched call.

## Details

Generalized linear Gaussian process models extend generalized linear
models (GLMs) by incorporating spatial Gaussian processes to account for
spatial correlation in the data. This function fits GLGPMs using maximum
likelihood methods, allowing for Gaussian, binomial, and Poisson
response families. In the case of the Binomial and Poisson families, a
Monte Carlo maximum likelihood algorithm is used.

The spatial Gaussian process is modeled with a Matern correlation
function, which is flexible and commonly used in geostatistical
modeling. The function supports both spatial covariates and unstructured
random effects, providing a comprehensive framework to analyze spatially
correlated data across different response distributions.

Additionally, the function allows for the inclusion of unstructured
random effects, specified through the
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
term in the model formula. These random effects can capture unexplained
variability at specific locations beyond the fixed and spatial covariate
effects, enhancing the model's flexibility in capturing complex spatial
patterns.

The `convert_to_crs` argument can be used to reproject the spatial
coordinates to a different CRS. The `scale_to_km` argument scales the
coordinates to kilometers if set to TRUE.

The `control_mcmc` argument specifies the control parameters for MCMC
sampling. This argument must be an object returned by
[`set_control_sim`](https://claudiofronterre.github.io/RiskMap/reference/set_control_sim.md).

The `start_pars` argument allows for specifying starting values for the
model parameters. If not provided, default starting values are used.

## See also

[`set_control_sim`](https://claudiofronterre.github.io/RiskMap/reference/set_control_sim.md),
[`summary.RiskMap`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md),
[`to_table`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
