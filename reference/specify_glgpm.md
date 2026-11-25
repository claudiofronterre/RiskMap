# Specify a model for simulation

Defines a data-generating model without fitting it. Parameters are
supplied on their natural scale and remain fixed across simulations.

## Usage

``` r
specify_glgpm(
  formula,
  data,
  family,
  parameters,
  invlink = NULL,
  denominator = NULL,
  distance_crs = NULL,
  distance_units = c("km", "m")
)
```

## Arguments

- formula:

  A formula object specifying the model to be fitted. The formula should
  include fixed effects and spatial effects specified using
  [`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md)
  and optionally, random effects specified using
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md).
  It can also include one
  [`offset()`](https://rdrr.io/r/stats/offset.html) term that names a
  column of `data`, e.g. `offset(log_exposure)`.

- data:

  An sf object containing the variables in the model.

- family:

  The distribution of the response given the random effects, with its
  link. Either one of the strings `"gaussian"`, `"binomial"` and
  `"poisson"`, which take the identity, logit and log links, a family
  object as in [`stats::glm()`](https://rdrr.io/r/stats/glm.html), such
  as `binomial("probit")`, or
  [`negative_binomial()`](https://claudiofronterre.github.io/RiskMap/reference/negative_binomial.md)
  for overdispersed counts with the log link, whose size is estimated
  unless it is held, by `fixed_pars` of
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  or by the family object `MASS::negative.binomial(2)`. Where nothing is
  estimated, as in `specify_glgpm()`, the size must be given, as `theta`
  in its `parameters`. Binomial responses take the logit, probit,
  complementary log-log (`"cloglog"`) and Cauchy (`"cauchit"`) links,
  Poisson responses the log link and Gaussian responses the identity
  link. `binomial("log")` and `poisson("identity")` are not accepted,
  because they can give means outside the range of the family, nor is
  `poisson("sqrt")`, whose inverse gives the same mean to a positive and
  a negative linear predictor. Quasi-likelihood families, such as
  [`stats::quasibinomial()`](https://rdrr.io/r/stats/family.html), have
  no likelihood to maximise and are not accepted either. Any other link
  can be given in `invlink`. A named link gives the same estimates as
  its inverse given in `invlink`, and is faster when that inverse cannot
  be differentiated symbolically, as for `invlink = pnorm`.

- parameters:

  A named list containing `beta`, `sigma2` and `phi`. `beta` has one
  entry per design-matrix column, including the intercept when present.
  Supply `tau2` when `gp(nugget = TRUE)`, `sigma2_me` for Gaussian
  measurement error, `theta` for the size of a negative binomial family,
  and `sigma2_re` for each
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  term. Variances may be zero. `phi` is expressed in `distance_units`.

- invlink:

  A function that maps the linear predictor to the mean of the response
  given the random effects, for a binomial or Poisson model whose link
  is not one of those that `family` names, such as a link for an
  imperfect diagnostic test. Its first and second derivatives are found
  symbolically where possible and numerically otherwise. They only steer
  the importance sampler, so they affect the speed of a fit and not its
  estimates. Not applicable for Gaussian models, and not to be combined
  with a family object that names a link other than logit or log. Custom
  links are experimental, and how they are given and checked may change
  in a minor release.

- denominator:

  Optional denominator for binomial distributions, the number of trials,
  or exposure for Poisson and negative binomial ones, the mean of a
  count being the exposure times the inverse link of the linear
  predictor. A fit needs positive whole numbers of trials, and positive
  exposures. Counts and numbers of trials that are whole to
  floating-point error are rounded. Passed as a bare/unquoted column
  name present in `data`. If not provided, defaults to `1`.

- distance_crs:

  The projected coordinate reference system in which distances, and so
  `phi`, are measured, for example an EPSG code. `data` are transformed
  to it. The default, `NULL`, keeps projected data as they are and
  projects longitude/latitude data to the CRS from
  [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md),
  with a message naming it. A warning of class `riskmap_warning_crs` is
  raised when the CRS used distorts distances by more than 1% over the
  data.

- distance_units:

  Character string, either `"km"` or `"m"`, indicating the units used
  for internal coordinate distances and spatial parameters such as
  `phi`. The coordinates of the data keep the linear units of their CRS.
  Defaults to `"km"`.

## Value

A validated `RiskMap_simulation_model` for
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md).
It contains the formula, transformed model data, family, fixed
parameters, denominator name, distance units, inverse link, fixed-effect
design information, random-effect terms, offset, Matern smoothness and
response.

## Details

The response in `formula` need not exist in `data`.

## See also

[print()](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_simulation_model.md)
describes the model, and
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
simulates from it.

## Examples

``` r
library(sf)
locations <- st_as_sf(data.frame(x = c(0, 1000, 2000), y = 0),
                      coords = c("x", "y"), crs = 32629)

grid <- st_as_sf(data.frame(x = c(500, 1500), y = 500),
                 coords = c("x", "y"), crs = 32629)

model <- specify_glgpm(response ~ gp(),
                       data = locations,
                       family = "gaussian",
                       parameters = list(beta = 1,
                                         sigma2 = 1,
                                         phi = 2,
                                         sigma2_me = 0.1))
```
