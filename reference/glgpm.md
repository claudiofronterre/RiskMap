# Estimation of Generalized Linear Gaussian Process Models

Fits generalized linear Gaussian process models to spatial data,
incorporating spatial Gaussian processes with a Matern correlation
function. Supports Gaussian, binomial, Poisson and negative binomial
response families.

## Usage

``` r
glgpm(
  formula,
  data,
  family,
  invlink = NULL,
  denominator = NULL,
  distance_crs = NULL,
  distance_units = c("km", "m"),
  start_pars = NULL,
  fixed_pars = NULL,
  control_mcml = set_control_mcml(),
  control_mcmc = set_control_mcmc(),
  return_samples = TRUE,
  seed = NULL,
  messages = getOption("RiskMap.messages", TRUE)
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
  unless it is held, by `fixed_pars` of `glgpm()` or by the family
  object `MASS::negative.binomial(2)`. Where nothing is estimated, as in
  [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md),
  the size must be given, as `theta` in its `parameters`. Binomial
  responses take the logit, probit, complementary log-log (`"cloglog"`)
  and Cauchy (`"cauchit"`) links, Poisson responses the log link and
  Gaussian responses the identity link. `binomial("log")` and
  `poisson("identity")` are not accepted, because they can give means
  outside the range of the family, nor is `poisson("sqrt")`, whose
  inverse gives the same mean to a positive and a negative linear
  predictor. Quasi-likelihood families, such as
  [`stats::quasibinomial()`](https://rdrr.io/r/stats/family.html), have
  no likelihood to maximise and are not accepted either. Any other link
  can be given in `invlink`. A named link gives the same estimates as
  its inverse given in `invlink`, and is faster when that inverse cannot
  be differentiated symbolically, as for `invlink = pnorm`.

- invlink:

  A function that maps the linear predictor to the mean of the response
  given the random effects, for a binomial or Poisson model whose link
  is not one of those that `family` names, such as a link for an
  imperfect diagnostic test. Its first and second derivatives are found
  symbolically where possible and numerically otherwise. They only steer
  the importance sampler, so they affect the speed of a fit and not its
  estimates. Not applicable for Gaussian models, and not to be combined
  with a family object that names a link other than logit or log.

- denominator:

  Optional denominator for binomial distributions, or exposure for
  Poisson and negative binomial ones, the mean of a count being the
  exposure times the inverse link of the linear predictor. Passed as a
  bare/unquoted column name present in `data`. If not provided, defaults
  to `1`.

- distance_crs:

  The projected coordinate reference system in which distances, and so
  `phi`, are measured, for example an EPSG code. `data` are transformed
  to it, and the model is fitted in it. The default, `NULL`, keeps
  projected data as they are and projects data in longitude and latitude
  to the CRS from
  [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md),
  a UTM zone or a projection centred on the data, with a message naming
  the CRS and the units of distance. A CRS in longitude and latitude is
  refused. A warning of class `riskmap_warning_crs` is raised when the
  CRS used distorts distances by more than 1% over the data.
  [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md),
  [`summarise_distance()`](https://claudiofronterre.github.io/RiskMap/reference/summarise_distance.md),
  [`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md)
  and
  [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
  take the same argument, so that their distances can be measured as the
  fit's are.

- distance_units:

  Character string, either `"km"` or `"m"`, indicating the units used
  for internal coordinate distances and spatial parameters such as
  `phi`. The coordinates of the data keep the linear units of their CRS.
  Defaults to `"km"`.

- start_pars:

  `NULL` (the default), or a named list of starting values for any of
  the model's estimated parameters, one, several or all of them, named
  as in
  [riskmap_parameters](https://claudiofronterre.github.io/RiskMap/reference/riskmap_parameters.md).
  `beta` is matched to the columns of the design matrix by name, or by
  position if it is unnamed and has a value for each, and `sigma2_re` to
  the
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  terms in the same way. The flat names of `unlist(coef(fit))`, such as
  `"beta.x"`, are accepted too, in a list or as that named vector, so
  `start_pars = coef(fit)` starts a fit from the estimates of another.
  Every other parameter takes its automatic start.
  [`glgpm_parameters()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm_parameters.md)
  lists the parameters of a model before it is fitted.

- fixed_pars:

  `NULL` (the default), or a named list of values at which to hold any
  of the model's parameters instead of estimating them, one, several or
  all of them, with the names and structure of `start_pars`. Every value
  must be finite, and every one but the coefficients positive, except
  the measurement error variance of a Gaussian model, which may be zero.
  A held coefficient enters the linear predictor as a known term, as an
  offset does. The other parameters are estimated with the held ones at
  their values, and with every parameter held the fit evaluates the
  likelihood at the values (exactly for a Gaussian model, by Monte Carlo
  on one sample for the others) and draws the conditional samples there,
  so that [`predict()`](https://rdrr.io/r/stats/predict.html) works with
  given parameters. A parameter can't be both started and held, except
  that a start equal to its held value, as in `start_pars = coef(fit)`
  of a fit with the same holds, is left out. With one observation per
  location, a Gaussian model can estimate the nugget or the measurement
  error variance, not both, and a negative binomial model the nugget or
  the size, not both.

- control_mcml:

  Control parameters for adaptive Monte Carlo maximum likelihood
  estimation. Must be returned by
  [`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md).

- control_mcmc:

  Control parameters for MCMC sampling for binomial or Poisson models.
  Must be an object of class `RiskMap_control_mcmc` as returned by
  [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md).

- return_samples:

  Logical indicating whether to draw and retain conditional samples at
  the fitted parameters. These samples are used to compute fitted
  values. If `FALSE`, final sampling is skipped and fitted values are
  not computed. Defaults to `TRUE`. Ignored for Gaussian models.

- seed:

  Optional non-negative integer that seeds every random step of a
  binomial or Poisson fit: the Monte Carlo samples of the estimation and
  the final conditional samples. When `NULL` (the default), one integer
  is drawn from the session's random number generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  also reproduces the fit. The seed used is returned as `seed`. The
  caller's random number state is otherwise left unchanged. Ignored for
  Gaussian models, which have no random step.

- messages:

  Logical indicating whether to print messages, and a line for each
  stage of Monte Carlo maximum likelihood as described in
  [riskmap_progress](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md).
  Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.

## Value

A fitted model of class `RiskMap`, a list. These components are part of
its interface.

- data:

  The `sf` data of the fit, in the CRS in which the model was fitted, so
  that distances between them are those that `phi` describes.
  `sf::st_crs(fit$data)` gives that CRS, and
  `sf::st_transform(fit$data, fit$input_crs)` gives the data in the CRS
  they were given in.

- input_crs:

  The CRS of `data` as it was given.

- formula:

  The model formula.

- family:

  The name of the family, such as `"binomial"`.

- call:

  The matched call.

- seed:

  The seed used for the random steps of a binomial, Poisson or negative
  binomial fit, whether given or drawn, so that
  `glgpm(..., seed = fit$seed)` reproduces the fit. `NULL` for Gaussian
  models.

- distance_units:

  The units of distance, `"km"` or `"m"`, in which `phi` is given.

The other components are internal and may change in any version of
RiskMap. The estimates, their covariance and the likelihood are read
with
[coef()](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md),
[vcov()](https://claudiofronterre.github.io/RiskMap/reference/logLik.RiskMap.md),
[logLik()](https://claudiofronterre.github.io/RiskMap/reference/logLik.RiskMap.md)
and
[nobs()](https://claudiofronterre.github.io/RiskMap/reference/logLik.RiskMap.md),
fitted values with
[fitted()](https://claudiofronterre.github.io/RiskMap/reference/fitted.RiskMap.md),
standard errors and intervals with
[summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md),
and how the fit was computed with
[`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md).

Every result of RiskMap that can be saved and read again, such as a fit,
a prediction or a profile, has the attribute `riskmap_version`, a list
with the `version` of RiskMap that made it, as a
[`package_version()`](https://rdrr.io/r/base/numeric_version.html), and
`sha`, the commit it was installed from when it was installed from
GitHub and `NA` otherwise. An object without it was made before RiskMap
stamped its results, by an earlier version or an early development build
of 2.0.0.

## Details

Generalized linear Gaussian process models extend generalized linear
models (GLMs) by incorporating spatial Gaussian processes to account for
spatial correlation in the data. This function fits GLGPMs using maximum
likelihood methods, allowing for Gaussian, binomial, Poisson and
negative binomial response families. In the case of the binomial,
Poisson and negative binomial families, a Monte Carlo maximum likelihood
algorithm is used. A negative binomial model is fitted as the binomial
model of its counts at each value of its size (see
[`negative_binomial()`](https://claudiofronterre.github.io/RiskMap/reference/negative_binomial.md)).

The spatial Gaussian process is modeled with a Matern correlation
function, which is flexible and commonly used in geostatistical
modeling. The function supports both spatial covariates and unstructured
random effects, providing a comprehensive framework to analyze spatially
correlated data across different response distributions.

Additionally, the function allows for the inclusion of unstructured
random effects, specified through the
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
term in the model formula. These random effects can capture unexplained
variability at specific locations beyond the fixed and spatial covariate
effects, enhancing the model's flexibility in capturing complex spatial
patterns.

The model is fitted in the projected CRS of `distance_crs`, or of `data`
when it is `NULL`, and the data in longitude and latitude are projected
to the CRS from
[`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
with a message. The `distance_units` argument controls whether distances
between locations (and so parameters such as `phi`) are expressed in
kilometers (`"km"`, the default) or meters (`"m"`).

The `control_mcmc` argument specifies the control parameters for MCMC
sampling. This argument must be an object returned by
[`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md).
Non-Gaussian models are fitted by certified Monte Carlo maximum
likelihood (MCML), controlled by `control_mcml` (see
[`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md)).
The estimate is corrected with an independent validation sample and
certified with a second, fresh one, which bounds its Monte Carlo error
in statistical standard errors, with an allowance for that sample's own
error taken from it and from a replicate of it (see
[`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md)).
The samples are enlarged only when the certificate exceeds the
tolerance, and a warning is issued if the tolerance is not met within
the configured limits.
[`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md)
reports every stage. Warnings are kept for what changes how a result is
read, and
[`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md)
goes through all the checks on the computation, with what to do about
any that need a look. The estimation samples and the final conditional
samples are drawn from separate streams derived from `seed`, so
`glgpm(..., seed = fit$seed)` reproduces a fit exactly.

Parameters not given in `start_pars` take RiskMap's defaults, which are
the coefficients of the regression without latent effects, `sigma2 = 1`,
`phi` at the 10th percentile of the distances between locations, 1 for
an estimated nugget variance, a measurement error variance and the
variance of each
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
term, and for the size of a negative binomial model that of the
regression without latent effects. A binomial, Poisson or negative
binomial model with the canonical link and without
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
terms compares them with the estimates of a Gaussian model fitted to the
empirical link of the data (the empirical logit of a binomial model),
with one observation per location, the totals of the location where it
has several records, and starts from the one with the higher Laplace
approximation of the likelihood. The values given in `start_pars`
replace their own in both, and
`set_control_mcml(start_comparison = FALSE)` skips the comparison. A
start of `phi` far from the distances between the locations, below a
hundredth of the smallest or above ten times the largest, gives a
warning, since it may be in other units.
[`summary.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)
says which starting values a fit used and why.

[`stats::update()`](https://rdrr.io/r/stats/update.html) refits a model
with some of its arguments changed, as in
`update(fit, . ~ . + elevation)` or
`update(fit, fixed_pars = list(phi = 20))`. It evaluates the call again,
with the data it names in the environment where
[`update()`](https://rdrr.io/r/stats/update.html) is called.

## See also

[`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md),
[`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md),
[`summary.RiskMap`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md),
[`to_table`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)

## Examples

``` r

data(italy_sim)

# intercept only
fit <- glgpm(
  formula = y ~ gp(),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

summary(fit)
#> ── Linear geostatistical model ─────────────────────────────────────────────────
#> glgpm(formula = y ~ gp(), data = italy_sim[1:100, ], family = "gaussian",
#>     messages = FALSE)
#> The identity link.
#> The limits are those of 95% confidence intervals.
#> 
#> Regression coefficients
#>             Estimate Lower limit Upper limit   StdErr z.value   p.value    
#> (Intercept)  9.32340     7.70078    10.94603  0.82788  11.262 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Measurement error
#>                        Estimate Lower limit Upper limit
#> Measurement error var.  0.17859     0.13335      0.2392
#> 
#> Spatial Gaussian process, Matern correlation with kappa = 0.5
#>                          Estimate Lower limit Upper limit
#> Spatial process var.       4.8494      1.8196      12.924
#> Spatial corr. scale (km)  39.0150      2.7773     548.084
#> Variance of the nugget effect fixed at 0
#> 
#> Log-likelihood   8.618636
#> AIC             -9.237272

# add a covariate
fit <- update(fit, . ~ . + pop_dens)

summary(fit)
#> ── Linear geostatistical model ─────────────────────────────────────────────────
#> glgpm(formula = y ~ gp() + pop_dens, data = italy_sim[1:100,
#>     ], family = "gaussian", messages = FALSE)
#> The identity link.
#> The limits are those of 95% confidence intervals.
#> 
#> Regression coefficients
#>               Estimate Lower limit Upper limit     StdErr z.value   p.value    
#> (Intercept) 7.70187620  6.37937916  9.02437323 0.67475578 11.4143 < 2.2e-16 ***
#> pop_dens    0.00046825  0.00023917  0.00069733 0.00011688  4.0063 6.168e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Measurement error
#>                        Estimate Lower limit Upper limit
#> Measurement error var.  0.17853     0.13331      0.2391
#> 
#> Spatial Gaussian process, Matern correlation with kappa = 0.5
#>                          Estimate Lower limit Upper limit
#> Spatial process var.      1.77451     0.59997      5.2484
#> Spatial corr. scale (km) 75.04610     9.01027    625.0556
#> Variance of the nugget effect fixed at 0
#> 
#> Log-likelihood   14.34117
#> AIC             -18.68233

# add a random effect
fit <- glgpm(
  formula = y ~ gp() + re(province),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

summary(fit)
#> ── Linear geostatistical model ─────────────────────────────────────────────────
#> glgpm(formula = y ~ gp() + re(province), data = italy_sim[1:100,
#>     ], family = "gaussian", messages = FALSE)
#> The identity link.
#> The limits are those of 95% confidence intervals.
#> 
#> Regression coefficients
#>             Estimate Lower limit Upper limit   StdErr z.value   p.value    
#> (Intercept)  9.34077     7.92351    10.75804  0.72311  12.918 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Measurement error
#>                        Estimate Lower limit Upper limit
#> Measurement error var.  0.17859     0.13335      0.2392
#> 
#> Spatial Gaussian process, Matern correlation with kappa = 0.5
#>                          Estimate Lower limit Upper limit
#> Spatial process var.     0.438837    0.022556      8.5376
#> Spatial corr. scale (km) 0.917141    0.000000         Inf
#> Variance of the nugget effect fixed at 0
#> 
#> Unstructured random effects
#>                             Estimate Lower limit Upper limit
#> province (random eff. var.)   4.2742      1.4733        12.4
#> 
#> Log-likelihood   9.082286
#> AIC             -8.164572

# estimate the nugget
fit <- glgpm(
  formula = y ~ gp(nugget = TRUE),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

summary(fit)
#> ── Linear geostatistical model ─────────────────────────────────────────────────
#> glgpm(formula = y ~ gp(nugget = TRUE), data = italy_sim[1:100,
#>     ], family = "gaussian", messages = FALSE)
#> The identity link.
#> The limits are those of 95% confidence intervals.
#> 
#> Regression coefficients
#>             Estimate Lower limit Upper limit   StdErr z.value   p.value    
#> (Intercept)  9.32340     7.70078    10.94603  0.82788  11.262 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Measurement error
#>                        Estimate Lower limit Upper limit
#> Measurement error var.  0.17859     0.13335      0.2392
#> 
#> Spatial Gaussian process, Matern correlation with kappa = 0.5
#>                            Estimate Lower limit Upper limit
#> Spatial process var.     4.8494e+00  1.8196e+00      12.924
#> Spatial corr. scale (km) 3.9015e+01  2.7773e+00     548.084
#> Variance of the nugget   7.3615e-08  0.0000e+00         Inf
#> 
#> Log-likelihood   8.618636
#> AIC             -7.237272

# binomial model
data(loaloa)

# reduce iterations to speed up the fit, and set a
# seed for reproducibility
control_mcmc <- set_control_mcmc(n_draws = 8,
                                 burnin = 20)

fit <- glgpm(NO_INF ~ gp(),
             data = loaloa[1:50,],
             denominator = NO_EXAM,
             family = "binomial",
             control_mcmc = control_mcmc,
             seed = 42,
             messages = FALSE)
#> ℹ `data` are in longitude and latitude, so they are projected to a custom
#>   oblique Mercator projection ("RiskMap oblique Mercator 11.5E 5N azimuth -72.5
#>   k0 0.9998") for fitting, with distances in km. No UTM zone keeps distances
#>   within 0.5% over them. Set `distance_crs` to choose the projection.

summary(fit)
#> ── Binomial geostatistical model ───────────────────────────────────────────────
#> glgpm(formula = NO_INF ~ gp(), data = loaloa[1:50, ], family = "binomial",
#>     denominator = NO_EXAM, control_mcmc = control_mcmc, seed = 42,
#>     messages = FALSE)
#> The canonical logit link.
#> The limits are those of 95% confidence intervals.
#> 
#> Regression coefficients
#>             Estimate Lower limit Upper limit  StdErr z.value   p.value    
#> (Intercept)  -2.1233     -3.2385     -1.0081  0.5690 -3.7317 0.0001902 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Spatial Gaussian process, Matern correlation with kappa = 0.5
#>                          Estimate Lower limit Upper limit
#> Spatial process var.       2.1220      0.9347      4.8173
#> Spatial corr. scale (km) 102.1364     34.7203    300.4540
#> Variance of the nugget effect fixed at 0
#> 
#> Log-likelihood  -2018.33
#> 
#> Monte Carlo maximum likelihood
#> ✔ Certified. The Monte Carlo error of the estimate is within 0.0055 statistical
#>   standard errors (an approximate 95% bound, tolerance 0.1).
#> ℹ The starting values came from a Gaussian fit to the transformed data, which
#>   the Laplace approximation of the likelihood preferred to RiskMap's defaults.
#> ✔ No check needs a look. See `check_fit()` for the details.

# family = binomial("probit") fits the same model with a probit link
```
