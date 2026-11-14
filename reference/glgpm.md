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
  model_crs = NULL,
  distance_units = c("km", "m"),
  control_mcmc = set_control_mcmc(),
  control_mcml = set_control_mcml(),
  seed = NULL,
  return_samples = TRUE,
  messages = getOption("RiskMap.messages", TRUE),
  fix_var_me = NULL,
  start_pars = list(beta = NULL, sigma2 = NULL, tau2 = NULL, phi = NULL, sigma2_me =
    NULL, sigma2_re = NULL)
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
  [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md),
  the size must be given. Binomial responses take the logit, probit,
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

- control_mcmc:

  Control parameters for MCMC sampling for binomial or Poisson models.
  Must be an object of class `RiskMap_control_mcmc` as returned by
  [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md).

- control_mcml:

  Control parameters for adaptive Monte Carlo maximum likelihood
  estimation. Must be returned by
  [`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md).

- seed:

  Optional non-negative integer that seeds every random step of a
  binomial or Poisson fit: the Monte Carlo samples of the estimation and
  the final conditional samples. When `NULL` (the default), one integer
  is drawn from the session's random number generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  also reproduces the fit. The seed used is returned as `seed`. The
  caller's random number state is otherwise left unchanged. Ignored for
  Gaussian models, which have no random step.

- return_samples:

  Logical indicating whether to draw and retain conditional samples at
  the fitted parameters. These samples are used to compute fitted
  values. If `FALSE`, final sampling is skipped and fitted values are
  not computed. Defaults to `TRUE`. Ignored for Gaussian models.

- messages:

  Logical indicating whether to print messages, and a line for each
  stage of Monte Carlo maximum likelihood as described in
  [riskmap_progress](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md).
  Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.

- fix_var_me:

  Optional fixed value for the measurement error variance when fitting a
  Gaussian model. When not provided, the value will be estimated, but
  cannot be if each location only has one sample and the `nugget` term
  in
  [`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md)
  is also set to `TRUE`.

- start_pars:

  Optional list of starting values for model parameters:

  beta

  :   regression coefficients

  sigma2

  :   spatial process variance

  tau2

  :   nugget effect variance

  phi

  :   spatial correlation scale

  sigma2_me

  :   measurement error variance - only for Gaussian models

  sigma2_re

  :   random effects variances - only when random effects are included

  theta

  :   size of the negative binomial family - only when it is estimated
      (see
      [`negative_binomial()`](https://claudiofronterre.github.io/RiskMap/reference/negative_binomial.md))

## Value

An object of class `RiskMap` containing the fitted model and relevant
information:

- estimate:

  Estimated parameters as a named list: `beta` (named regression
  coefficients on the scale of the supplied covariates), `sigma2`,
  `phi`, `nu2` \\= \tau^2/\sigma^2\\ when the nugget is estimated,
  `sigma2_me` for Gaussian models with an estimated measurement error
  variance, `theta` for negative binomial models whose size is
  estimated, stored as \\\log(1/\theta)\\, the opposite sign to the
  logarithm of a variance ([`coef()`](https://rdrr.io/r/stats/coef.html)
  gives \\\theta\\), and `sigma2_re` (a named vector, one entry per
  random effect) when random effects are included. Keeping these in
  separate list elements - rather than one flat named vector - avoids a
  covariate whose name collides with a parameter name (e.g. a covariate
  literally called `"sigma2"`) corrupting name-based lookups.
  `covariance` is estimated as one joint matrix over every parameter;
  its rows/columns are named to match `unlist(estimate)`.

- grad_MLE:

  Gradient of the maximum likelihood function

- covariance:

  Covariance

- log_lik:

  Log likelihood

- y:

  Response variable

- D:

  Covariate matrix

- coords:

  Unique spatial coordinates

- ID_coords:

  Index of coordinates

- re:

  Random effects if
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  is included in `formula`

- ID_re:

  Index of random effects if
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  is included in `formula`

- fix_tau2:

  Fixed nugget effect variance

- fix_var_me:

  Fixed measurement error variance

- family_parameters:

  For a negative binomial model whose size is held, the size as
  `list(theta = )`. Absent otherwise.

- formula:

  Model formula

- family:

  Response family

- link_function:

  For binomial and Poisson models, the inverse link (`inv`), its first
  and second derivatives (`d1` and `d2`) and the name of the link
  (`link`), such as `"probit"`, which is `NULL` for an inverse given in
  `invlink`. `NULL` for Gaussian models.

- distance_units:

  Units (`"km"` or `"m"`) used for internal coordinate distances and
  spatial parameters

- data:

  The `sf` data used for model fitting, in the CRS used internally

- input_crs:

  The CRS of the input data

- kappa:

  Spatial correlation parameter

- units_m:

  The denominators or exposures of a binomial, Poisson or negative
  binomial model

- cov_offset:

  Covariate offset

- call:

  Matched call

- S_samples:

  For non-Gaussian models, fresh conditional MCMC samples evaluated at
  the fitted parameters when `return_samples` is `TRUE`. These are not
  the reference samples used for MCML estimation. The first block is the
  location-level Gaussian effect (including the nugget when present),
  followed by any unstructured random-effect blocks.

- linear_predictors:

  Conditional fitted linear predictors. For non-Gaussian models these
  are Monte Carlo means over the final-parameter latent samples and are
  `NULL` when `return_samples = FALSE`; for Gaussian models they use the
  conditional latent mean.

- fitted_values:

  Conditional fitted response means. Binomial models return
  probabilities; Poisson and negative binomial models return expected
  counts and therefore include any supplied exposure in `denominator`.
  These are `NULL` for non-Gaussian models when
  `return_samples = FALSE`.

- mcml_history:

  For non-Gaussian models, one entry per stage of certified Monte Carlo
  maximum likelihood (see
  [`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md)),
  containing the sizes of the fitting and validation samples, the
  importance distribution used, the optimum of the fitting sample
  (`fitted_estimate`) and the corrected estimate (`estimate`), the
  optimizer result, the size of the correction, the certificate and its
  two parts (`certificate_distance` and `certificate_radius`), effective
  sample sizes and Pareto tail diagnostics of each sample, the outcome
  of the stage and its timing. A stage where the global importance
  distribution failed also holds why, in `proposal_failure` or, at the
  start of the first stage, `start_failure`: the `reason`, the number of
  `refits` and the last `change` of its factors. Estimates use the same
  parameterisation as `estimate`, with regression coefficients and
  spatial scale restored to the supplied covariate and distance units.
  `NULL` for Gaussian models.

- mcml_certified:

  Whether the certificate of the returned estimate is within the
  tolerance. The certificate is an approximate 95% bound, in statistical
  standard errors, on the Monte Carlo error of the estimate, that is its
  distance from the maximum of the exact likelihood that the fit has
  located. It says nothing about whether that maximum is the highest.
  `NULL` for Gaussian models.

- mcml_selected_iteration:

  The stage whose estimate is returned. `NULL` for Gaussian models.

- seed:

  The seed used for the random steps of a binomial or Poisson fit,
  whether supplied or drawn. `NULL` for Gaussian models.

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

The `model_crs` argument can be used to reproject the spatial
coordinates to a different, projected CRS for model fitting. If `data`
are in longitude/latitude and `model_crs` is not supplied, the
coordinates are automatically reprojected to the CRS from
[`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
and a message details the conversion; automatic projection is avoided
for `model_crs` explicitly supplied in longitude/latitude, which raises
an error. The `distance_units` argument controls whether distances
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
tolerance. Every stage is recorded in `mcml_history`, and a warning is
issued if the tolerance is not met within the configured limits.
Warnings are kept for what changes how a result is read, and
[`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md)
goes through all the checks on the computation, with what to do about
any that need a look. The estimation samples and the final conditional
samples are drawn from separate streams derived from `seed`, so
`glgpm(..., seed = fit$seed)` reproduces a fit exactly.

The `start_pars` argument allows starting values to be supplied
explicitly. Otherwise a binomial, Poisson or negative binomial model
starts from RiskMap's defaults, which are the coefficients of the
regression without latent effects, `sigma2 = 1`, `phi` at the 10th
percentile of the distances between locations, and 1 for an estimated
nugget variance and for the variance of each
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
term. A model with the canonical link and without
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
terms compares them with the estimates of a Gaussian model fitted to the
empirical link of the data (the empirical logit of a binomial model),
with one observation per location, the totals of the location where it
has several records, and starts from the one with the higher Laplace
approximation of the likelihood.
[`summary.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)
says which starting values a fit used and why.

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
fit <- glgpm(
  formula = y ~ gp() + pop_dens,
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

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
control_mcmc <- set_control_mcmc(n_sim = 100,
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
#>   within 0.5% over them. Set `model_crs` to choose the projection.

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
#> Spatial corr. scale (km) 102.1364     34.7203    300.4537
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
