# Empirical Logit Transformation for Binomial Data

Computes the empirical logit transformation for binomial counts.

## Usage

``` r
elogit(y, m)
```

## Arguments

- y:

  A numeric vector of observed positives.

- m:

  A numeric vector of denominators, i.e. the numbers tested.

## Value

A numeric vector containing the empirical logit transformation
\\\log((y + 0.5) / (m - y + 0.5))\\.

## Details

The empirical logit is often used as a finite transformation for
binomial data, including cases where \\y = 0\\ or \\y = m\\.

## Examples

``` r
y <- c(0, 3, 7, 10)
m <- c(10, 10, 10, 10)
elogit(y, m)
#> [1] -3.0445224 -0.7621401  0.7621401  3.0445224
```
