# First Derivative with Respect to \\\phi\\

Computes the first derivative of the Matern correlation function with
respect to \\\phi\\.

## Usage

``` r
matern.grad.phi(U, phi, kappa)
```

## Arguments

- U:

  A vector of distances between pairs of data locations.

- phi:

  The scale parameter \\\phi\\.

- kappa:

  The smoothness parameter \\\kappa\\.

## Value

A matrix with the values of the first derivative of the Matern function
with respect to \\\phi\\ for the given distances.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
