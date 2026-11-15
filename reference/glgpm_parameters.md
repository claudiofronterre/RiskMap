# List the parameters of a model before fitting it

The parameters
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
estimates or holds for a model, with the names that `start_pars` and
[`coef.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md)
use
([riskmap_parameters](https://claudiofronterre.github.io/RiskMap/reference/riskmap_parameters.md)).
Only the formula, the data and the family are read, so this is quick.

## Usage

``` r
glgpm_parameters(
  formula,
  data,
  family,
  fix_var_me = NULL,
  distance_units = c("km", "m")
)
```

## Arguments

- formula, data, family, fix_var_me:

  As in
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

- distance_units:

  The units of distance, `"km"` (the default) or `"m"`, for the label of
  the scale.

## Value

A data frame of class `riskmap_parameters`, with a row per element of
each parameter and the columns `name`, `element` (the column of the
design matrix of a coefficient, the term of a random-effect variance),
`flat` (the name in `unlist(coef(fit))`), `label` (as in
[`summary.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)),
`start` (the rule of its automatic start, see
[riskmap_parameters](https://claudiofronterre.github.io/RiskMap/reference/riskmap_parameters.md))
and `held` (the value at which the model holds the parameter, `NA` when
it is estimated). Its print shows the template of a `start_pars` list.

## See also

[riskmap_parameters](https://claudiofronterre.github.io/RiskMap/reference/riskmap_parameters.md),
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)

## Examples

``` r
data(liberia)
liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
glgpm_parameters(npos ~ log(elevation) + gp(nugget = TRUE), data = liberia,
                 family = "binomial")
#>    name        element                flat                    label
#>    beta    (Intercept)    beta.(Intercept)              (Intercept)
#>    beta log(elevation) beta.log(elevation)           log(elevation)
#>  sigma2                             sigma2     Spatial process var.
#>     phi                                phi Spatial corr. scale (km)
#>    tau2                               tau2   Variance of the nugget
#>              start held
#>         regression     
#>         regression     
#>                one     
#>  distance_quantile     
#>                one     
#> 
#> A template for `start_pars`:
#> list(beta = c("(Intercept)" = , "log(elevation)" = ), sigma2 = , phi = , tau2 = )
```
