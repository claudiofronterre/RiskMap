# Random Effect Model Specification

Specifies the terms for a random effect model.

## Usage

``` r
re(...)
```

## Arguments

- ...:

  Variables representing the random effects in the model. Each variable
  adds its own independent random effect.

## Value

A list of class `RiskMap_re_spec` containing the following elements:

- term:

  A character vector of the specified terms.

- dim:

  The number of specified terms.

- label:

  A character string representing the full call for the random effect
  model.

## Details

The function constructs a list that includes the specified terms for the
random effects. This list can be used as a specification for a random
effect model. Several `re()` calls in a formula are combined, so
`re(region) + re(province)` is the same model as `re(region, province)`.

## Note

At least one variable must be provided as input.

## Examples

``` r
# this function is not intended for use by itself
re(province)
#> $term
#> [1] "province"
#> 
#> $dim
#> [1] 1
#> 
#> $label
#> [1] "re(province)"
#> 
#> attr(,"class")
#> [1] "RiskMap_re_spec"

# use in a model formula
data(italy_sim)
fit <- glgpm(
  formula = y ~ gp() + re(province),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)
```
