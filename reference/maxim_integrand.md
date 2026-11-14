# Maximization of the Integrand for Generalized Linear Gaussian Process Models

Finds the conditional mode of the latent effects in a binomial or
Poisson generalized linear Gaussian process model (GLGPM), that is, the
maximiser of the joint log-density of the responses and the latent
spatial and unstructured random effects, viewed as a function of the
latent effects.

## Usage

``` r
maxim_integrand(
  y,
  denominator,
  mu,
  Sigma,
  ID_coords,
  ID_re = NULL,
  family,
  sigma2_re = NULL,
  hessian = FALSE,
  gradient = FALSE,
  invlink = NULL
)
```

## Arguments

- y:

  Response variable vector: numbers of positives for
  `family = "binomial"`, or counts for `family = "poisson"`.

- denominator:

  Numeric vector with one positive value per response: the number of
  trials for `family = "binomial"`, or the exposure (the expected count
  is `denominator` times the inverse link of the linear predictor) for
  `family = "poisson"`. This is the same quantity as the `denominator`
  argument of
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

- mu:

  Vector of the fixed part of the linear predictor, one value per
  response.

- Sigma:

  Covariance matrix of the spatial process.

- ID_coords:

  Indices mapping response to locations.

- ID_re:

  Indices mapping response to unstructured random effects.

- family:

  The distribution of the response, `"binomial"` or `"poisson"`, a
  binomial or Poisson family object with one of the links that
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  accepts, such as `binomial("probit")`, or
  [`negative_binomial()`](https://claudiofronterre.github.io/RiskMap/reference/negative_binomial.md)
  with its size, such as `negative_binomial(theta = 2)`, with
  `denominator` the exposure.

- sigma2_re:

  Variance of the unstructured random effects.

- hessian:

  Logical; if TRUE, compute the Hessian matrix.

- gradient:

  Logical; if TRUE, compute the gradient vector.

- invlink:

  A function that maps the linear predictor to the mean of the response
  given the random effects, for a link that `family` does not name, as
  in
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

## Value

A list with the mode (`mode`), the name of the inverse link used
(`invlink_used`) and either the Hessian of the integrand at the mode
(`hessian`, if `hessian = TRUE`) or the inverse of its negative
(`Sigma.tilde`, otherwise). If `gradient = TRUE`, the gradient at the
mode is also returned (`gradient`).

## Details

The mode is found by Newton-type optimisation
([`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)) using the
analytic gradient and Hessian of the integrand. The integrand includes
terms for the spatial process (`Sigma`), the unstructured random effects
(`sigma2_re`), and the log-likelihood of the responses for the specified
family. With the canonical links, the terms of a response whose
probability rounds to 0 or 1, or whose mean underflows to 0, are
computed in their exact canonical form, so that the optimisation is not
stopped there.

## Examples

``` r
coords <- cbind(c(0, 1, 2), c(0, 0, 1))
Sigma <- exp(-as.matrix(dist(coords)))
mode <- maxim_integrand(
  y = c(3, 5, 1), denominator = c(10, 10, 10), mu = rep(-0.5, 3),
  Sigma = Sigma, ID_coords = 1:3, family = "binomial"
)
mode$mode
#> [1] -0.1997816  0.2532412 -0.8857314
```
