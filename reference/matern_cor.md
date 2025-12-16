# Matern Correlation Function

Computes the Matern correlation function.

## Usage

``` r
matern_cor(u, phi, kappa, return_sym_matrix = FALSE)
```

## Arguments

- u:

  A vector of distances between pairs of data locations.

- phi:

  The scale parameter \\\phi\\.

- kappa:

  The smoothness parameter \\\kappa\\.

- return_sym_matrix:

  A logical value indicating whether to return a symmetric correlation
  matrix. Defaults to `FALSE`.

## Value

A vector of the same length as `u` with the values of the Matern
correlation function for the given distances, if
`return_sym_matrix=FALSE`. If `return_sym_matrix=TRUE`, a symmetric
correlation matrix is returned.

## Details

The Matern correlation function is defined as

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk> \$\$\rho(u; \phi; \kappa)
= (2^{\kappa-1})^{-1}(u/\phi)^\kappa K\_{\kappa}(u/\phi)\$\$ where
\\\phi\\ and \\\kappa\\ are the scale and smoothness parameters, and
\\K\_{\kappa}(\cdot)\\ denotes the modified Bessel function of the third
kind of order \\\kappa\\. The parameters \\\phi\\ and \\\kappa\\ must be
positive.
