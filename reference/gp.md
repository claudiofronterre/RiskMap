# Gaussian Process Model Specification

Specifies the terms, smoothness, and nugget effect for a Gaussian
Process (GP) model.

## Usage

``` r
gp(..., kappa = 0.5, nugget = FALSE)
```

## Arguments

- ...:

  Variable representing the spatial coordinates for the GP model. If
  left blank the `geometry` column from the data is used automatically.

- kappa:

  The smoothness parameter \\\kappa\\. Default is `0.5`.

- nugget:

  The nugget effect, which represents the variance of the measurement
  error. Default is `FALSE` in which case it is not estimated. If `TRUE`
  the value will be estimated or a positive numeric value can be
  provided instead to fix the effect.

## Value

A list of class `RiskMap_gp_spec` containing the following elements:

- term:

  A character vector of the specified terms.

- kappa:

  The smoothness parameter \\\kappa\\.

- nugget:

  The nugget effect.

- dim:

  The number of specified terms.

- label:

  A character string representing the full call for the GP model.

## Details

The function constructs a list that includes the specified terms
(spatial coordinates or covariates), the smoothness parameter
\\\kappa\\, and the nugget effect. This list can be used as a
specification for a Gaussian Process model. A model has one spatial
Gaussian process, so a formula can contain only one `gp()` term.

## Examples

``` r
# this function is not intended for use by itself
gp()
#> $term
#> [1] "sf"
#> 
#> $kappa
#> [1] 0.5
#> 
#> $nugget
#> [1] 0
#> 
#> $dim
#> [1] 0
#> 
#> $label
#> [1] "gp()"
#> 
#> attr(,"class")
#> [1] "RiskMap_gp_spec"
gp(nugget = TRUE)
#> $term
#> [1] "sf"
#> 
#> $kappa
#> [1] 0.5
#> 
#> $nugget
#> [1] TRUE
#> 
#> $dim
#> [1] 0
#> 
#> $label
#> [1] "gp()"
#> 
#> attr(,"class")
#> [1] "RiskMap_gp_spec"

# use in a model formula
data(italy_sim)
fit <- glgpm(
  formula = y ~ gp(),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)
```
