# Simulation from Generalized Linear Gaussian Process Models

Simulates data from a fitted Generalized Linear Gaussian Process Model
(GLGPM) or a specified model formula and data.

## Usage

``` r
glgpm_sim(
  n_sim,
  model_fit = NULL,
  formula = NULL,
  data = NULL,
  family = NULL,
  den = NULL,
  cov_offset = NULL,
  crs = NULL,
  convert_to_crs = NULL,
  scale_to_km = TRUE,
  sim_pars = list(beta = NULL, sigma2 = NULL, tau2 = NULL, phi = NULL, sigma2_me = NULL,
    sigma2_re = NULL),
  messages = TRUE
)
```

## Arguments

- n_sim:

  Number of simulations to perform.

- model_fit:

  Fitted GLGPM model object of class 'RiskMap'. If provided, overrides
  'formula', 'data', 'family', 'crs', 'convert_to_crs', 'scale_to_km',
  and 'control_mcmc' arguments.

- formula:

  Model formula indicating the variables of the model to be simulated.

- data:

  Data frame or 'sf' object containing the variables in the model
  formula.

- family:

  Distribution family for the response variable. Must be one of
  'gaussian', 'binomial', or 'poisson'.

- den:

  Required for 'binomial' to denote the denominator (i.e. number of
  trials) of the Binomial distribution. For the 'poisson' family, the
  argument is optional and is used a multiplicative term to express the
  mean counts.

- cov_offset:

  Offset for the covariate part of the GLGPM.

- crs:

  Coordinate reference system (CRS) code for spatial data.

- convert_to_crs:

  CRS code to convert spatial data if different from 'crs'.

- scale_to_km:

  Logical; if TRUE, distances between locations are computed in
  kilometers; if FALSE, in meters.

- sim_pars:

  List of simulation parameters including 'beta', 'sigma2', 'tau2',
  'phi', 'sigma2_me', and 'sigma2_re'.

- messages:

  Logical; if TRUE, display progress and informative messages.

## Value

A list containing simulated data, simulated spatial random effects (if
applicable), and other simulation parameters.

## Details

Generalized Linear Gaussian Process Models (GLGPMs) extend generalized
linear models (GLMs) by incorporating spatial Gaussian processes to
model spatial correlation. This function simulates data from GLGPMs
using Markov Chain Monte Carlo (MCMC) methods. It supports Gaussian,
binomial, and Poisson response families, utilizing a Matern correlation
function to model spatial dependence.

The simulation process involves generating spatially correlated random
effects and simulating responses based on the fitted or specified model
parameters. For 'gaussian' family, the function simulates response
values by adding measurement error.

Additionally, GLGPMs can incorporate unstructured random effects
specified through the
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
term in the model formula, allowing for capturing additional variability
beyond fixed and spatial covariate effects.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
