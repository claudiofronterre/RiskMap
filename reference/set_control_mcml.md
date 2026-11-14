# Set control parameters for certified MCML estimation

Non-Gaussian models are fitted by Monte Carlo maximum likelihood with
importance sampling. The estimated log-likelihood is maximised on a
fitting sample of draws, the estimate is corrected by one Newton step
with the score from an independent validation sample, followed by
further steps on that sample when the fitting sample predicts that one
step leaves a residual, and the corrected estimate is certified with a
second, fresh validation sample. The radius of the certificate, its
allowance for the Monte Carlo error of that sample's score, is the
larger of the radii of that sample and of a replicate of the same size
drawn at the same estimate. The radius of one sample is estimated from
its own draws, and is too small when the sample misses the rare draws of
very large weight that dominate its error, which the replicate misses
independently. The certificate is an approximate 95% bound, in
statistical standard errors, on the Monte Carlo error of the estimate,
that is its distance from the maximum of the exact likelihood that the
fit has located. It is not an exact finite-sample guarantee, and it says
nothing about whether that maximum is the highest (see `check_maxima`).
If the certificate exceeds `tolerance`, the samples are enlarged to the
size predicted to meet it and the fit is repeated.

## Usage

``` r
set_control_mcml(
  initial_samples = 1000L,
  validation_samples = 5000L,
  max_samples = 16000L,
  max_stages = 4L,
  tolerance = 0.1,
  proposal = c("global", "laplace"),
  defensive_weight = 0.05,
  laplace_start = NULL,
  integrate_noise = "auto",
  check_maxima = TRUE,
  pc_prior = FALSE,
  pc_sigma_threshold = 1,
  pc_sigma_probability = 0.05,
  pc_scale_threshold = NULL,
  pc_scale_probability = 0.05
)
```

## Arguments

- initial_samples:

  Positive integer giving the number of draws in the fitting sample at
  the first stage. Defaults to 1,000.

- validation_samples:

  Positive integer giving the number of draws in each validation sample,
  and in the replicate of the second, at the first stage. Later stages
  enlarge the validation samples in the same proportion as the fitting
  sample. Defaults to 5,000.

- max_samples:

  Positive integer giving the largest permitted fitting sample. Defaults
  to 16,000.

- max_stages:

  Positive integer giving the maximum number of stages. Defaults to 4.

- tolerance:

  Positive number giving the largest accepted certificate, in
  statistical standard errors. Defaults to 0.1.

- proposal:

  The importance distribution for the latent effects. `"global"` (the
  default) fits a Gaussian approximation to the whole conditional
  distribution of the latent effects, which keeps far more of the draws
  effective in sparse data than the Laplace approximation; `"laplace"`
  uses the Laplace approximation at the conditional mode. If the global
  approximation cannot be built at some parameter value, the stage is
  repeated with the Laplace approximation.

- defensive_weight:

  Number in \[0, 1) giving the share of importance draws taken from the
  prior distribution of the latent effects instead of the Gaussian
  approximation. This defensive mixture keeps the variance of the
  importance weights finite, which a Gaussian proposal alone does not in
  most spatial datasets, so that Monte Carlo errors are valid. Defaults
  to 0.05. Set to 0 to use the Gaussian approximation alone.

- laplace_start:

  Logical or `NULL`. If `TRUE`, the Laplace approximation to the
  likelihood is first maximised from the starting values, and Monte
  Carlo maximum likelihood starts from that estimate. If the global
  importance distribution cannot be built at that estimate, as can
  happen in sparse data when it has a very large nugget, the fit starts
  from the starting values instead. `NULL` (the default) means `TRUE`,
  except where `integrate_noise = TRUE` integrates a term out, where it
  means `FALSE`, because there the Laplace approximation's spurious
  maximum moves to the spatial effect and the start from it costs more
  without certifying more fits. With `integrate_noise = "auto"` and a
  nugget as the only term that can be integrated out, the Laplace
  estimate chooses the form of the first stage, which with the nugget
  integrated out starts from the starting values.

- integrate_noise:

  Logical, or `"auto"`. If `TRUE`, a binomial or Poisson model
  integrates its independent Gaussian terms out of the likelihood
  numerically wherever the data allow it, and Monte Carlo maximum
  likelihood samples the spatial effect and the terms that stay. The
  model, its parameters and its likelihood are the same, so the option
  changes only how the same maximum is found. Sparse surveys with a
  large nugget, or with household effects, then keep far more effective
  draws. The terms are the nugget, whose levels are the locations, and
  each
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  term. Records at a location with the same covariates, offset and
  random-effect levels are first combined into one unit. Then four rules
  apply.

  - A level of an
    [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
    term that holds a single unit, such as a household with one record
    or with records that do not differ, is integrated out of that unit's
    likelihood, and so is the nugget of a location with a single unit.

  - At a location with several units the nugget is integrated out as
    well, through an integral over the nugget of the units' own
    integrals, when the units share their covariates and offset and
    differ only in the levels of
    [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
    terms that each hold a single unit, as households of a village with
    covariates measured at the village do.

  - At a location with several units whose covariates or offsets differ,
    the units are integrated out in groups when every term but one has a
    single level shared by all of them and the levels of that one lie
    within the location. Its levels are the groups, such as the
    households of a cluster whose children differ in sex or birth order,
    each integrated over its own effect with its units' own covariates,
    and the nugget is integrated over the groups' integrals. Where every
    term is shared, the location is one group under its nugget.

  - At any other location the nugget stays sampled.

  Every other level stays sampled. That includes the levels of an
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  term that span several locations, such as villages or districts with
  several sampled points, and the levels at a location where two
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  terms vary, such as households and the individuals within them. A term
  can be integrated out at some of its levels and sampled at others.
  [`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md)
  reports, for each term, at how many levels it was integrated out and
  at how many it stayed sampled. It works with any link, including a
  custom `invlink`. Where the nugget variance is close to zero the
  likelihood is nearly flat along it, and there the Monte Carlo estimate
  of its effect is noisier with the nugget integrated out than sampled,
  so fewer such fits are certified. With `"auto"` each stage of the fit
  integrates out the nugget if its variance at the start of the stage is
  at least 0.1 and samples it otherwise. The first stage takes the
  variance from the Laplace estimate, or from the starting values if
  `laplace_start = FALSE`. Every validation bank is then drawn with the
  nugget sampled and integrated out, and the two estimates of the score
  are combined with weights inverse to their Monte Carlo covariances, so
  the certificate rests on the more precise form, or on both. A bank
  whose importance weights have a heavy tail, with a Pareto k above 0.8,
  can understate its covariance and is left out of the combination. A
  model with an
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  term that can be integrated out is fitted with `"auto"` as with
  `TRUE`, its
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  terms and its nugget integrated out at every stage wherever the rules
  above allow. On simulated surveys of households that certified every
  fit, and in a fraction of the time of sampling, also where the nugget
  variance was zero, where switching the nugget alone sometimes failed.
  Defaults to `"auto"`, which on simulated surveys with a nugget
  certified more fits than either form alone. A model without a nugget
  or
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
  terms is fitted in the same way whatever the value.

- check_maxima:

  Logical. The certificate establishes Monte Carlo precision at the
  maximum found, not that no higher maximum exists. In sparse surveys an
  optimiser can stop at a maximum with a spatial scale far below the
  spacing between locations while a higher one exists at a longer scale,
  or the reverse. If `TRUE` (the default), a fitted scale below half the
  median distance between neighbouring locations is refitted from four
  dispersed starting values on the same Monte Carlo draws. After a fit
  at a longer scale the likelihood is evaluated at four shorter scales,
  with the coefficients and variance of the model with independent
  location effects (the limit of the model as the scale goes to zero,
  computed exactly by quadrature), and the fit is repeated from the best
  of them if it comes within 0.25 log-likelihood units of the fit. Exact
  lower bounds on the likelihood at short scales force a further search
  if they exceed the fitted log-likelihood. Separate maxima are compared
  on independent draws and the highest is kept, with a warning if it
  differs. A fit whose spatial correlation is numerically negligible is
  reported as numerically independent, with the scale unresolved and
  without standard errors or certificate. The screen of short scales and
  the exact bounds are available for the canonical links and models
  without unstructured random effects.

- pc_prior:

  Logical. If `TRUE`, the spatial variance and scale are penalised by
  penalised-complexity priors (Fuglstad et al., 2019) and the estimate
  maximises the penalised likelihood on the log scale of the covariance
  parameters. This shrinks estimates away from degenerate solutions in
  small or sparse datasets. It is recommended for small or
  low-prevalence surveys, especially before prediction with
  `parameter_uncertainty = TRUE` in
  [`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md).
  When the spatial variance is large and well identified, the penalty
  shrinks it somewhat towards zero, so the default is `FALSE`.

- pc_sigma_threshold, pc_sigma_probability:

  The penalty sets P(sigma \> `pc_sigma_threshold`) =
  `pc_sigma_probability`, where sigma is the spatial standard deviation
  on the scale of the linear predictor. Default 1 and 0.05.

- pc_scale_threshold, pc_scale_probability:

  The penalty sets P(phi \< `pc_scale_threshold`) =
  `pc_scale_probability`, where phi is the spatial scale in the units of
  the coordinates. By default the threshold is the median distance
  between neighbouring locations and the probability is 0.05.

## Value

A control object of class `RiskMap_control_mcml`.

## See also

[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md)

## Examples

``` r
control_mcml <- set_control_mcml(max_samples = 8000)
```
