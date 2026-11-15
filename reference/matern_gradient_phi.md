# First Derivative with Respect to \\\phi\\

Computes the first derivative of the Matern correlation function with
respect to \\\phi\\.

## Usage

``` r
matern_gradient_phi(u, phi, kappa)
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

A matrix with the values of the first derivative of the Matern function
with respect to \\\phi\\ for the given distances.

## Examples

``` r
x <- matrix(rnorm(20), nrow = 5)
u <- dist(x)
matern_gradient_phi(u, phi = 1, kappa = 0.8)
#>            [,1]      [,2]      [,3]       [,4]      [,5]
#> [1,] 0.00000000 0.3541365 0.2366189 0.09682591 0.3168688
#> [2,] 0.35413651 0.0000000 0.4411938 0.21998710 0.4329899
#> [3,] 0.23661893 0.4411938 0.0000000 0.35643536 0.4272446
#> [4,] 0.09682591 0.2199871 0.3564354 0.00000000 0.2327614
#> [5,] 0.31686878 0.4329899 0.4272446 0.23276138 0.0000000
```
