# Second Derivative with Respect to \\\phi\\

Computes the second derivative of the Matern correlation function with
respect to \\\phi\\.

## Usage

``` r
matern_hessian_phi(u, phi, kappa)
```

## Arguments

- u:

  A distance matrix returned from
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html)

- phi:

  The scale parameter \\\phi\\.

- kappa:

  The smoothness parameter \\\kappa\\.

## Value

A matrix with the values of the second derivative of the Matern function
with respect to \\\phi\\ for the given distances.

## Examples

``` r
x <- matrix(rnorm(20), nrow = 5)
u <- dist(x)
matern_hessian_phi(u, phi = 1, kappa = 0.8)
#>               [,1]          [,2]          [,3]          [,4]          [,5]
#> [1,] -1.180330e-14  1.920059e-01  1.463741e-01  1.771557e-01 -3.818350e-01
#> [2,]  1.920059e-01 -1.180330e-14 -2.119411e-01  1.046032e-01  2.062844e-01
#> [3,]  1.463741e-01 -2.119411e-01 -1.180330e-14  1.441201e-01  1.299346e-01
#> [4,]  1.771557e-01  1.046032e-01  1.441201e-01 -1.180330e-14  2.151108e-01
#> [5,] -3.818350e-01  2.062844e-01  1.299346e-01  2.151108e-01 -1.180330e-14
```
