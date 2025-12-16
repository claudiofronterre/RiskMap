# Maximization of the Integrand for Generalized Linear Gaussian Process Models

Maximizes the integrand function for Generalized Linear Gaussian Process
Models (GLGPMs), which involves the evaluation of likelihood functions
with spatially correlated random effects.

## Usage

``` r
maxim.integrand(
  y,
  units_m,
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

- family:

  Distribution family for the response variable. Must be one of
  'gaussian', 'binomial', or 'poisson'.

- sigma2_re:

  Variance of the unstructured random effects.

- hessian:

  Logical; if TRUE, compute the Hessian matrix.

- gradient:

  Logical; if TRUE, compute the gradient vector.

## Value

A list containing the mode estimate, and optionally, the Hessian matrix
and gradient vector.

## Details

This function maximizes the integrand for GLGPMs using the Nelder-Mead
optimization algorithm. It computes the likelihood function
incorporating spatial covariance and unstructured random effects, if
provided.

The integrand includes terms for the spatial process (Sigma),
unstructured random effects (sigma2_re), and the likelihood function
(llik) based on the specified distribution family ('gaussian',
'binomial', or 'poisson').

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
