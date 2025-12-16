# Laplace Sampling Markov Chain Monte Carlo (MCMC) for Generalized Linear Gaussian Process Models

Performs MCMC sampling using Laplace approximation for Generalized
Linear Gaussian Process Models (GLGPMs).

## Usage

``` r
Laplace_sampling_MCMC(
  y,
  units_m,
  mu,
  Sigma,
  ID_coords,
  ID_re = NULL,
  sigma2_re = NULL,
  family,
  control_mcmc,
  invlink = NULL,
  Sigma_pd = NULL,
  mean_pd = NULL,
  messages = TRUE
)
```

## Arguments

- y:

  Response variable vector.

- units_m:

  Units of measurement for the response variable.

- mu:

  Mean vector of the response variable.

- Sigma:

  Covariance matrix of the spatial process.

- ID_coords:

  Indices mapping response to locations.

- ID_re:

  Indices mapping response to unstructured random effects.

- sigma2_re:

  Variance of the unstructured random effects.

- family:

  Distribution family for the response variable. Must be one of
  'gaussian', 'binomial', or 'poisson'.

- control_mcmc:

  List with control parameters for the MCMC algorithm:

  n_sim

  :   Number of MCMC iterations.

  burnin

  :   Number of burn-in iterations.

  thin

  :   Thinning parameter for saving samples.

  h

  :   Step size for proposal distribution. Defaults to
      1.65/(n_tot^(1/6)).

  c1.h, c2.h

  :   Parameters for adaptive step size tuning.

- Sigma_pd:

  Precision matrix (optional) for Laplace approximation.

- mean_pd:

  Mean vector (optional) for Laplace approximation.

- messages:

  Logical; if TRUE, print progress messages.

## Value

An object of class "mcmc.RiskMap" containing:

- samples\$S:

  Samples of the spatial process.

- samples\$\<re_names\[i\]\>:

  Samples of each unstructured random effect, named according to columns
  of ID_re if provided.

- tuning_par:

  Vector of step size (h) values used during MCMC iterations.

- acceptance_prob:

  Vector of acceptance probabilities across MCMC iterations.

## Details

This function implements a Laplace sampling MCMC approach for GLGPMs. It
maximizes the integrand using \`maxim.integrand\` function for Laplace
approximation if \`Sigma_pd\` and \`mean_pd\` are not provided.

The MCMC procedure involves adaptive step size adjustment based on the
acceptance probability (\`acc_prob\`) and uses a Gaussian proposal
distribution centered on the current mean (\`mean_curr\`) with variance
\`h\`.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
