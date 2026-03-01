# Laplace-sampling MCMC for Generalized Linear Gaussian Process Models

Runs Markov chain Monte Carlo (MCMC) sampling using a Laplace
approximation for Generalized Linear Gaussian Process Models (GLGPMs).
The latent Gaussian field is integrated via a second-order Taylor
expansion around the mode, and a Gaussian proposal is used for
Metropolis–Hastings updates with adaptive step-size tuning.

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

  Numeric vector of responses of length \\n\\. For `family = "binomial"`
  this is the number of successes, for `family = "poisson"` counts, and
  for `family = "gaussian"` real values.

- units_m:

  Numeric vector giving the binomial totals (number of trials) when
  `family = "binomial"`; ignored for other families (can be `NULL`).

- mu:

  Numeric vector of length equal to the number of unique locations
  providing the mean of the latent spatial process on the link scale.

- Sigma:

  Numeric positive-definite covariance matrix for the latent spatial
  process \\S\\ at the unique locations referenced by `ID_coords`.

- ID_coords:

  Integer vector of length \\n\\ mapping each response in `y` to a
  row/column of `Sigma` (i.e., the index of the corresponding location).

- ID_re:

  Optional matrix or data.frame with one column per unstructured random
  effect (RE). Each column is an integer vector of length \\n\\ mapping
  observations in `y` to RE levels (e.g., cluster, survey, etc.). Use
  `NULL` to exclude REs.

- sigma2_re:

  Optional named numeric vector of RE variances. Names must match the
  column names of `ID_re`. Ignored if `ID_re = NULL`.

- family:

  Character string: one of `"gaussian"`, `"binomial"`, or `"poisson"`.

- control_mcmc:

  List of control parameters:

  n_sim

  :   Total number of MCMC iterations (including burn-in).

  burnin

  :   Number of initial iterations to discard.

  thin

  :   Thinning interval for saving samples.

  h

  :   Initial step size for the Gaussian proposal. Defaults to \\1.65 /
      n\_\mathrm{tot}^{1/6}\\ if not supplied.

  c1.h, c2.h

  :   Positive tuning constants for adaptive step-size updates.

- invlink:

  Optional inverse-link function. If `NULL`, defaults are used:
  `identity` (gaussian), `plogis` (binomial), and `exp` (poisson).

- Sigma_pd:

  Optional precision matrix used in the Laplace approximation. If
  `NULL`, it is obtained internally at the current mode.

- mean_pd:

  Optional mean vector used in the Laplace approximation. If `NULL`, it
  is obtained internally as the mode of the integrand.

- messages:

  Logical; if `TRUE`, prints progress and acceptance diagnostics.

## Value

An object of class `"mcmc.RiskMap"` with components:

- samples:

  A list containing posterior draws. Always includes `$S` (latent
  spatial field). If `ID_re` is supplied, each unstructured RE is
  returned under `$<re_name>`.

- tuning_par:

  Numeric vector of step sizes (`h`) used over iterations.

- acceptance_prob:

  Numeric vector of Metropolis–Hastings acceptance probabilities.

## Details

The algorithm alternates between:

1.  Locating the mode of the joint integrand for the latent variables
    (via `maxim.integrand`) when `Sigma_pd` and `mean_pd` are not
    provided, yielding a Gaussian approximation.

2.  Metropolis–Hastings updates using a Gaussian proposal centered at
    the current approximate mean with proposal variance governed by `h`.
    The step size is adapted based on empirical acceptance probability.

Dimensions must be consistent: `length(y) = n`,
`nrow(Sigma) = ncol(Sigma) = n_loc`, and `length(ID_coords) = n` with
entries in \\1,\dots,n\_\mathrm{loc}\\. If `ID_re` is provided, each
column must have length \\n\\; when `sigma2_re` is supplied, it must be
named and match `colnames(ID_re)`.

## Default links

The default inverse links are: identity (gaussian), logistic (binomial),
and exponential (poisson). Supply `invlink` to override.

## See also

[`maxim.integrand`](https://claudiofronterre.github.io/RiskMap/reference/maxim.integrand.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterre@lancaster.ac.uk>
