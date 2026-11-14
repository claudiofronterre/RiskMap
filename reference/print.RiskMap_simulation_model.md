# Print a simulation model or simulations

Prints a short description of a model specification from
[`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
or of simulations from
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md),
instead of every stored value.

## Usage

``` r
# S3 method for class 'RiskMap_simulation_model'
print(x, ...)

# S3 method for class 'RiskMap_simulation'
print(x, ...)
```

## Arguments

- x:

  Output from
  [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
  or
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md).

- ...:

  Not used.

## Value

Invisibly returns `x`.
