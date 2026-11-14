# Matern Correlation Function

Computes the Matern correlation function.

## Usage

``` r
matern_correlation(u, phi, kappa, return_sym_matrix = FALSE)
```

## Arguments

- u:

  A vector or matrix of distances between pairs of data locations.

- phi:

  The scale parameter \\\phi\\.

- kappa:

  The smoothness parameter \\\kappa\\.

- return_sym_matrix:

  A logical value indicating whether to return a symmetric correlation
  matrix. Defaults to `FALSE`.

## Value

A vector or matrix of the same shape as `u` with the values of the
Matern correlation function for the given distances, if
`return_sym_matrix=FALSE`. If `return_sym_matrix=TRUE`, a symmetric
correlation matrix is returned.

## Details

The Matern correlation function is defined as \$\$\rho(u; \phi; \kappa)
= (2^{\kappa-1})^{-1}(u/\phi)^\kappa K\_{\kappa}(u/\phi)\$\$ where
\\\phi\\ and \\\kappa\\ are the scale and smoothness parameters, and
\\K\_{\kappa}(\cdot)\\ denotes the modified Bessel function of the third
kind of order \\\kappa\\. The parameters \\\phi\\ and \\\kappa\\ must be
positive.

## Examples

``` r
u <- c(10, 20, 30)
matern_correlation(u, phi = 1, kappa = 0.8)
#> [1] 1.141202e-04 6.320912e-09 3.230773e-13
```
