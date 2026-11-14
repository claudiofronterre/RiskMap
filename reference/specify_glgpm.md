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
  denominator = NULL,
  model_crs = NULL,
  distance_units = c("km", "m"),
  invlink = NULL
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
  unless a value is given (`negative_binomial(theta = 2)` or
  `MASS::negative.binomial(2)`). Where nothing is estimated, as in
  `specify_glgpm()`, the size must be given. Binomial responses take the
  logit, probit, complementary log-log (`"cloglog"`) and Cauchy
  (`"cauchit"`) links, Poisson responses the log link and Gaussian
  responses the identity link. `binomial("log")` and
  `poisson("identity")` are not accepted, because they can give means
  outside the range of the family, nor is `poisson("sqrt")`, whose
  inverse gives the same mean to a positive and a negative linear
  predictor. Quasi-likelihood families, such as
  [`stats::quasibinomial()`](https://rdrr.io/r/stats/family.html), have
  no likelihood to maximise and are not accepted either. Any other link
  can be given in `invlink`. A named link gives the same estimates as
  its inverse given in `invlink`, and is faster when that inverse cannot
  be differentiated symbolically, as for `invlink = pnorm`.

- parameters:

  A named list containing `beta`, `sigma2` and `phi`. `beta` has one
  entry per design-matrix column, including the intercept when present.
  Supply `tau2` when `gp(nugget = TRUE)`, `sigma2_me` for Gaussian
  measurement error, and `sigma2_re` for each
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  term. Variances may be zero. `phi` is expressed in `distance_units`.

- denominator:

  Optional denominator for binomial distributions, or exposure for
  Poisson and negative binomial ones, the mean of a count being the
  exposure times the inverse link of the linear predictor. Passed as a
  bare/unquoted column name present in `data`. If not provided, defaults
  to `1`.

- model_crs:

  Optional CRS (e.g. an EPSG code) used internally for model fitting. If
  `data` is already in a projected CRS and `model_crs` is not provided,
  the data are used as-is. If `data` are in longitude/latitude and
  `model_crs` is not provided, the data are automatically reprojected to
  the CRS from
  [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md),
  a UTM zone or a projection centred on the data, and a message names
  the CRS and the units of distance. If provided, `model_crs` must be a
  projected (not longitude/latitude) CRS. A warning of class
  `riskmap_warning_crs` is raised when the CRS used distorts distances
  by more than 1% over the data.

- distance_units:

  Character string, either `"km"` or `"m"`, indicating the units used
  for internal coordinate distances and spatial parameters such as
  `phi`. The coordinates in the returned `data` retain the linear units
  declared by its CRS. Defaults to `"km"`.

- invlink:

  A function that maps the linear predictor to the mean of the response
  given the random effects, for a binomial or Poisson model whose link
  is not one of those that `family` names, such as a link for an
  imperfect diagnostic test. Its first and second derivatives are found
  symbolically where possible and numerically otherwise. They only steer
  the importance sampler, so they affect the speed of a fit and not its
  estimates. Not applicable for Gaussian models, and not to be combined
  with a family object that names a link other than logit or log.

## Value

A validated `RiskMap_simulation_model` for
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md).
It contains the formula, transformed model data, family, fixed
parameters, denominator name, distance units, inverse link, fixed-effect
design information, random-effect terms, offset, Matern smoothness and
response.

## Details

The response in `formula` need not exist in `data`. CRS handling follows
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md):
`model_crs` is used when supplied, otherwise longitude/latitude data are
transformed to the CRS from
[`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md).
A message reports automatic transformations, and a warning a CRS that
distorts distances by more than 1% over the data.

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
