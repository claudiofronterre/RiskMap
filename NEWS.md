RiskMap 2.0.0
=============

- This version introduces many breaking changes from v1. The package should not
be considered stable, but further breaking changes will be handled gracefully.
- Test coverage has been increased from 0 to 92 %.
- All functions have examples
- Many functions have been renamed or removed:

| v1 | v2 | 
|----------|----------|
| `assess_pp` |  `assess_prediction` |
| `assess_sim` | `assess_simulation` |
| `check_mcmc` | `plot_mcmc` |
| `compute_ID_coords` | `create_ids` |
| `convex_hull_sf` | `create_convex_hull` |
| `dast()` | Removed; will be replaced in a separate package |
| `glgpm_sim` | `simulate_glgpm` |
| `Laplace_sampling_MCMC` | Removed; the replacement is internal |
| `matern.grad.phi` | `matern_gradient_phi ` |
| `matern.hessian.phi` | `matern_hessian_phi` |
| `matern_cor` | `matern_correlation` |
| `maxim.integrand` | `maxim_integrand` |
| `plot_score` / `plot_AnPIT` | `plot.RiskMap_cross_validation` (the two are merged; by default it draws every plot the `RiskMap_cross_validation` object has data for) |
| `plot_s_variogram` | `plot.RiskMap_variogram` |
| `pred_over_grid` | `predict` |
| `pred_target_grid` | `predict` |
| `pred_target_shp` | `predict` with `boundaries` |
| `set_control_sim` | `set_control_mcmc` |
| `surf_sim` | `simulate_glgpm` |
| `s_variogram` | `variogram` |

## Datasets
- Apart from `liberia`, all datasets are now in an simple features (sf) format.

## Coordinate reference systems

- Datasets must have a coordinate reference system. When coordinates are in
latitude/longitude and `model_crs` is not supplied, they are now automatically
reprojected to an appropriate UTM zone. `predict()`, `simulate_glgpm()` and
`assess_prediction()` have been updated to match; `simulate_glgpm()`'s own 
`convert_to_crs`/`scale_to_km` arguments have likewise been replaced with
`model_crs`/`distance_units`.
- providing `model_crs` in longitude/latitude now raises an error.
- In `glgpm()`, `convert_to_crs` has been replaced with `model_crs` and
`scale_to_km` has been replaced with `distance_units` (`"km"` or `"m"`).

## Exploratory analysis

- `plot()` for `variogram()` results gains `show_pairs` (default `FALSE`),
which labels each point with the number of pairs of locations in its distance
class, so that unstable semivariances from few pairs can be recognised.
- `plot()` for `variogram()` results draws the semivariances as points without
joining lines, omits distance classes without pairs (which previously caused
"Removed rows" warnings), starts both axes at zero, uses `theme_minimal()` as
other RiskMap plots do, and draws in red the points that fall outside the
envelope. `plot_envelope` is now checked.
- `variogram()` results gain a `print()` method that reports the distance
classes and the envelope settings before the table.
- In `summarise_distance()` and `variogram()`, `scale_to_km` has been
replaced with `distance_units` (`"km"` or `"m"`), for consistency with `glgpm()`.
Both now default to `"km"`, as `glgpm()` does. For `variogram()` this changes
the default from `"m"`, so code that gives `breaks` or `max_dist` in metres
should set `distance_units = "m"`.
- `variogram()` accepts `breaks` or `max_dist` beyond the largest distance
between locations, and warns only when whole distance classes lie beyond it.
When no pair falls in any class, the error reports the range of the distances
instead of always suggesting a change of units. `max_dist = NULL` can be passed
together with `breaks`.
- `variogram()` no longer fails when the envelope is built for a single
distance class, and runs about twice as fast with 2,000 or more locations,
with identical results. Rounding error no longer makes the envelope exclude one
curve too few at some levels, such as 0.9 with 999 permutations. The default
level of 0.95 with 999 permutations was not affected.
- In `summarise_distance()` and `variogram()`, `convert_to_utm` has been replaced
with `distance_crs`, matching `glgpm()`'s `model_crs` convention: an
already-projected CRS is retained as-is, longitude/latitude data is automatically
reprojected to an appropriate UTM zone with a message, and an explicitly supplied
`distance_crs` must be projected. Distances now correctly account for the CRS's
native linear units rather than assuming metres.
- `variogram()` now uses a simultaneous extreme-rank-length permutation envelope
  instead of separate pointwise intervals. It remains an exploratory diagnostic
  and does not report a hypothesis-test p-value. Its pairwise permutation engine
  is compiled, distances and lag bins are computed once, and the default is 999
  permutations plus the observed curve. The returned columns are now `distance`,
  `semivariance`, `n_pairs`, `lower_envelope` and `upper_envelope`. A `seed` can
  be supplied without changing the caller's random-number state.
  `plot_variogram()` (now `plot.RiskMap_variogram()`) shows the envelope by
  default and warns when it was omitted with `n_permutations = 0`.
- The `bins` parameter in `variogram()` has been removed and replaced with `breaks`.

## Model fitting

- Breaking change: the seed has moved out of `set_control_mcmc()`, whose
`seed` argument is removed. `glgpm()`, `predict()`, `assess_prediction()` and
`assess_simulation()` take a `seed` argument instead, as `simulate_glgpm()`,
`variogram()` and `profile()` already did. Replacing
`control_mcmc = set_control_mcmc(seed = s)` with `seed = s` gives identical
results. Without a seed, `glgpm()` draws one from the session's random number
generator, so `set.seed()` before the call reproduces the fit, and two calls
in a row differ. The seed used is returned as `seed` in the fitted object, and
`glgpm(..., seed = fit$seed)` reproduces any binomial or Poisson fit,
including its final conditional samples.
- The `den` argument in `glgpm()` has been renamed to `denominator`.
- `glpgm()` now only accepts data in an sf format. 
Consequently the locations do not need to be passed to `gp()` when fitting a model
as they are included automatically.
- The `nugget` parameter in `gp()` is now `FALSE` by default and should be set to
`TRUE` to estimate it.

- Canonical binomial and Poisson models now use compiled likelihood and MALA
kernels. The implementation retains the existing sampler, including support
for nugget effects and multiple unstructured random effects, while avoiding
the main interpreted-R bottlenecks. Custom inverse-link functions continue to
use the general R implementation.

- `plot_mcmc()` draws its diagnostics with ggplot2 and patchwork and returns
them as a `patchwork` object, instead of drawing with base graphics and
resetting `par()`. By default it shows the trace plot and the autocorrelation
function of the mean of the spatial random effects and of the location with
the smallest effective sample size, with a histogram of the effective sample
sizes at all the locations. The trace plots show the iterations at which the
samples were kept, so fits now record the `n_sim`, `burnin` and `thin` of
their samples in `attr(fit, "mcmc")`. `component` now defaults to the
location with the smallest effective sample size, can select several
locations and can be given together with `check_mean = TRUE`. The mean and
`component` previously included the levels of unstructured random effects,
which are stored after the locations in `S_samples`. They now refer to the
locations only.

- The new `summarise_mcmc()` returns a tibble of the effective sample sizes of
the mean of the spatial random effects, of each location and of each level of
an unstructured random effect. They are computed with the initial monotone
sequence estimator of Geyer (1992), so the package no longer depends on sns.
patchwork has moved from Suggests to Imports.

- Non-Gaussian optimisation now reuses the sample-wise Monte Carlo log
integrand when the objective, gradient, and Hessian are evaluated at identical
parameters, avoiding a repeated likelihood calculation.

- `profile()` computes profile likelihoods of binomial and Poisson models
fitted by certified Monte Carlo maximum likelihood, for single parameters and
for pairs, on one fixed bank of draws so that the profiles are smooth. It
reports profile and Wald intervals and says if a point on a profile is higher
than the estimate. Marginal grids of 11 points are extended with doubling steps
where the profile is flat, up to 11 points a side, and an end the profile does
not reach is reported as open. Interval ends are interpolated in the signed
root of the deviance, with one more point at a crossing between distant
points, and the points around each end are also solved from the estimate,
because constrained maxima can be multimodal. The constrained maximisations
work on the scale of the fit, with standardised covariates, which avoids the
false convergence that coefficients of covariates with large ranges caused in
their own units, and each starts from the extrapolation of the two previous
solutions. One that still stops with a non-zero code is solved again with a
central-difference gradient. A fixed coefficient whose constraint the start
cannot meet within the bounds is first moved to the nearest feasible point,
and each point records whether it was solved, is infeasible or failed. Joint
profiles are computed only when `pairs` asks for them, and the time left is estimated as the profiles
run. With `cores` above one, the two sides of each profile, the ends of the
intervals and the joint profiles run in parallel by forking (not on Windows),
with the same results as on one core with one BLAS thread. `plot()` draws the marginal profiles against their quadratic
approximation and the joint profiles as confidence regions with the marginal
ridges, as a corner plot for one pair.

- `profile()` also computes profile likelihoods of linear geostatistical
(Gaussian) models, from their exact likelihood, so `n_draws` and `seed` are
not used. Every estimated parameter can be profiled: the regression
coefficients, the spatial variance and scale, the variance of the nugget
(`tau2`, as `summary()` reports it, rather than its ratio to the spatial
variance), the measurement error variance unless it is fixed, and the
variances of the random effects, with the labels of `summary()`. For each
value of the covariance parameters the coefficients are found by generalised
least squares, and models with more than three observations per latent effect
use the Woodbury identity, as `glgpm()` does. The ends of the profile
intervals are refined until the deviance of a solved point is within 0.01 of
the cut-off, grid steps of the variances and the scale are at most 1 on the
log scale, so that a variance estimated near zero, with a huge standard
error, is still profiled, and each point of a joint profile is also solved
from the estimate, since a neighbour's solution with a variance near zero can
leave the optimiser on a plateau. A variance that the data cannot tell from
zero has a profile that stays below the cut-off down to its lower bound, and
its interval is reported as open there. The smoothness `kappa` stays fixed
and can be profiled by refitting over fixed values.

- In `plot()` of a profile, the subtitle of each marginal panel gives the
estimate, the profile interval and the Wald interval on separate lines, so
that it is no longer cut off in panels half or a third of a page wide, and
the two entries of the ridge legend of a joint profile are stacked.

- Non-Gaussian models are now fitted by certified Monte Carlo maximum
likelihood. The estimated likelihood is maximised on a fitting sample of 1,000
draws, the estimate is corrected by one Newton step with the score from an
independent validation sample of 5,000 draws, and the corrected estimate is
certified with a second, fresh validation sample. The certificate is an
approximate 95% bound, in statistical standard errors, on the Monte Carlo
error of the estimate at the maximum found. When it exceeds the tolerance
(0.1 by default), both samples are enlarged to the size predicted to meet it,
up to `max_samples`. The importance distribution is a Gaussian approximation fitted
to the whole conditional distribution of the latent effects (efficient
importance sampling), which keeps far more draws effective than the Laplace
approximation in sparse data. Binomial and Poisson models, with or without a
nugget and unstructured random effects, share one implementation, and the
Laplace approximation remains as a fallback and as an option. If the global
proposal fails while an estimate is being corrected or certified, the stage
is repeated with the Laplace proposal on fresh draws, so that every
certificate rests on validation draws taken after everything that led to it.
`set_control_mcml()` gains `validation_samples`, `tolerance` and `proposal`,
and loses `min_effective_samples`, `max_pareto_k`, `max_mcse_ratio` and
`score_probability`. The defaults of `max_samples` and `max_stages` are now
16,000 and 4. Fitted objects report `mcml_certified` in place of
`mcml_converged`.

- Certified MCML continues the Newton correction when one step leaves a
residual. Where the log-likelihood is far from quadratic between the fitted
estimate and the maximum, as along a spatial variance close to zero, one
Newton step stops short of the maximum, and the certificate, measured in the
curvature at the fitted estimate, sees only part of the shortfall. The fitting
sample now predicts the residual of the first step, at the cost of one
gradient on that sample, and when it exceeds a quarter of the validation
sample's radius, further Newton steps are taken with the score from the same
validation sample, each in the fitting sample's curvature at the current
point, until a step is below a tenth of the radius (at most four more steps).
The estimate is then certified in the curvature at the last step. If the steps
stop short by more than the certificate tolerates, they are discarded, the
stage is not certified, and the next stage starts from the estimate of the
first step. The check does not act where the validation sample's radius is
below 1e-6 statistical standard errors, as when the spatial variance is held
at zero and the likelihood is exact, because the predicted residual is then
the rounding error of the differences divided by the radius. Fits where the
check does not fire
are unchanged. Where it fires, the estimate and the certificate change, and so
do the returned standard errors, which now come from the curvature at the last
step. On a simulated survey with a log spatial variance of -4.56 at the
maximum, the standard error of the log variance moved from 14.5 to 29.5,
against 30.2 at the maximum. Each stage of `mcml_history` records the
predicted residual and the number of Newton steps.

- Certified MCML distinguishes Monte Carlo error from weak
identification. A variance component estimated as zero (and the spatial scale
when the spatial variance is zero) is held at the boundary, the first stage is
restarted once from an interior value, and precision is assessed on the
remaining parameters, which keep their standard errors. A spatial scale
estimated below the median distance between neighbouring locations is reported
as estimated, with a warning that gives its practical range, explains that
such a scale may be weakly identified, points to its profile likelihood and
suggests a nugget model, a fixed scale or the PC prior. If the information
matrix is still singular, fitting stops instead of drawing more samples. Each
case produces a plain-language warning.

- `maxim_integrand()` now takes `denominator` instead of `units_m`, matching
`glgpm()`. Its documentation no longer lists the Gaussian family, which it has
never supported, and describes the Newton-type optimisation it actually uses.

- Removed unused internal helpers (`log_mean_exp()`,
`importance_effective_sample_size()`, `restore_spatial_range()`) and the
compiled kernel `cpp_canonical_likelihood()`, which no longer had callers.

- `set_control_mcml()` gains `pc_prior` (default `FALSE`) and its thresholds
`pc_sigma_threshold`, `pc_sigma_probability`, `pc_scale_threshold` and
`pc_scale_probability`. When `pc_prior = TRUE`, the spatial variance and scale
are penalised by penalised-complexity priors (by default
P(sigma > 1) = 0.05 and P(phi < median nearest-neighbour distance) = 0.05),
which stabilises estimation in small or low-prevalence surveys. Warnings about
a spatial variance estimated as zero, a small spatial scale or a singular
information matrix now suggest this option.

- The final conditional sampling in `glgpm()` (used for fitted values and
stored prediction samples) now uses the compiled MALA kernel for canonical
links, as intended. It previously passed the inverse-link function itself,
which selected the general R implementation with numerical derivatives; a
binomial fit to the `liberia` data now takes about 2 s instead of 27 s.
Fresh conditional runs in `predict()` now honour a custom inverse link.

- `set_control_mcml()` gains `laplace_start` (default `TRUE`), which maximises
the Laplace approximation from the starting values before MCML.

- `set_control_mcml()` gains `check_maxima` (default `TRUE`). The
certificate establishes Monte Carlo precision at the maximum found, and in
sparse surveys an optimiser can stop at a maximum with a spatial scale far
below the spacing between locations while a higher one exists at a longer
scale. When the fitted scale is below half the median distance between
neighbouring locations, the fit is repeated from four dispersed starting
values on the same draws. Candidates within 0.1 standard errors of each other
are treated as one maximum, and separate maxima are compared on independent
draws. A message says when and why the check runs, and a warning says when it
moves the estimate to another maximum or finds two maxima that fit about
equally well. The check does not prove that no higher maximum exists.

- The check for other maxima also looks at short scales after a fit at a long
one. The likelihood with independent location effects, the limit of the model
as the spatial scale goes to zero, is computed exactly by adaptive
Gauss-Hermite quadrature and maximised once per survey. The estimated
likelihood is evaluated at four scales from an eighth of the spacing between
neighbouring locations to the spacing, with the coefficients and variance of
that limit and one three-point step in the log variance where the value is
close, and the fit is repeated from the best point if it comes within 0.25
log-likelihood units of the fit or if an evaluation fails. The limit and Jensen
lower bounds at those scales are exact, and if one of them exceeds the fitted
log-likelihood, estimated on independent draws, a further search is forced and
a warning is given if it does not resolve the alarm. A fit whose scale reaches
the point where the correlation between the closest locations is below 1e-6
ends with the outcome "numerically independent", without standard errors or
certificate, since the fit is then one of independent location effects and the
spatial scale is unresolved. A fit at a short positive scale whose likelihood is not clearly
above the limit at the current Monte Carlo precision is assessed as usual, with
a warning that points to `profile()` for the scale. Separate maxima whose
comparison on independent draws is ambiguous are compared again on four times
as many draws, and that comparison serves every later one. A stage whose
global proposal fails at its start goes straight to the Laplace proposal,
instead of being completed and then discarded, and the factorisation of the
spatial covariance is reused while only the regression coefficients change.
The screen
and the exact bounds need the canonical links and no unstructured random
effects, and other models keep the earlier check.

- Certified MCML computes the importance weights (for the canonical binomial
and Poisson links) and the dense linear algebra of the proposal in compiled
code, with the same LAPACK and BLAS routines R uses, and gives results
identical to the previous R code (up to the last bits on platforms whose
compilers fuse multiply and add). A fit no longer creates several large
temporary matrices at every evaluation of the likelihood, whose garbage
collection took about half of a fit.

- Fitted models have a `print()` method, which gives the call, the estimates
of the regression coefficients and covariance parameters on their natural
scale, the log-likelihood and, for binomial and Poisson models, whether the
Monte Carlo error of the estimate is certified. Printing a fit previously
listed every component of the object, including the data. `summary()` of a
binomial or Poisson fit now reports the certificate (its bound in statistical
standard errors, the level of the bound, the tolerance, the number of stages
and the draws at the last stage), with short notes on what the warnings of the
fit reported and the advice they gave, so that both are available for a saved
model. The certificate is in the `mcml` component of the summary.

- The messages and warnings of certified MCML quote the spatial scale of the
returned estimate, after the Newton correction. The message that announces
the refits from dispersed starting values, which are decided before the
correction, now says that it quotes the scale of the first fit, before its
correction. When the repeated Newton correction does not converge at the last
stage, the warning says so, instead of describing the bound as above the
tolerance.

## Predictions

- Grids created by `create_grid()` now have the class `RiskMap_grid` before
`data.frame`, as vctrs requires of data frame subclasses. dplyr verbs such as
`filter()` previously failed with an error on these grids and on the
predictions made from them.
- `create_grid()` no longer raises an error when `boundaries` is in
longitude/latitude; like `glgpm()`, it now automatically reprojects to an
appropriate UTM zone and reports the conversion with a message.
- The old three-stage prediction interface has been replaced by
  `predict.RiskMap()`. Prediction locations, predictors, offsets, grouping
  variables and weights are supplied together in one `sf` object. Supplying
  `boundaries` requests blockwise joint prediction and areal aggregation;
  grids created by `create_grid()` use exact fractional cell coverage.
  Location- and area-level predictions support exceedance probabilities and
  named draw-by-draw transformations for additional scientific targets.
  Response predictions can be standardized over a common empirical or
  counterfactual reference population, with optional population weights.
- Gaussian prediction uses a stable reduced covariance solve, avoiding cancellation
with small measurement-error variances while retaining the same conditional model.
- `predict()` gains `parameter_uncertainty` (default `FALSE`) and
`parameter_draws` (default 50). When requested, parameters are drawn from the
approximate sampling distribution of the estimate and, for each draw, the
latent effects are sampled at that value with a short run of the same Langevin
sampler used for plug-in prediction. The pooled draws give intervals,
exceedance probabilities and areal summaries that also reflect parameter
uncertainty, which matters for small or low-prevalence surveys.
- `plot()` maps the output of `predict()`. Default colours depend on what the
  map shows: a calm-to-alarm palette for predicted levels, a single neutral
  hue for uncertainty (`sd`, `cv`, `interval_width`), a blue-to-red palette
  centred at 50% for exceedance probabilities and a palette centred at zero
  for signed effects. The legend title names the quantity (for example
  "Mean prevalence" for binomial models) and probabilities are labelled as
  percentages. `spatial_scale = "both"` draws locations and areas side by
  side with a shared colour scale, several summaries can be drawn in one
  call, `breaks` gives a binned map, and boundaries supplied to `predict()`
  are outlined automatically. `print()` lists what a prediction contains and
  how to map it. `predict()` output now records `response_name`.
- `create_grid()` gains `keep`. `keep = "overlap"` keeps every cell that
  overlaps the boundary, including edge cells whose centre lies outside, so
  that totals such as population along a coastline are not lost. The grid
  metadata now also records the study boundary.
- Prediction maps draw grid cells as squares clipped to the grid's study
  boundary, so further `sf` layers (for example survey locations) can be added
  to a map in any coordinate reference system. Legend titles can be replaced
  with `labs(fill = ...)`. `plot()` also maps the fixed-effect contribution
  with `quantity = "fixed"`, several summaries at one scale are arranged in a
  grid, and side-by-side scales no longer carry panel titles.
- `plot()` now respects `limits` for maps of signed quantities such as the
  spatial effect, which were previously rescaled, and arranges several
  summaries in two columns.
- Maps of `create_grid()` grids outline the grid's study boundary, so the
  edge of the study area is visible where the clipped cells end.
- Several maps drawn in one call have wider side margins, so neighbouring
  legends no longer run into each other at narrow widths.
- `plot()` accepts several quantities, for example
  `quantity = c("fixed", "spatial")`, and draws them side by side as it
  does for several summaries.
- New `assess_extrapolation()` reports, for each prediction location,
  whether a prediction relies on extrapolation: covariates outside their
  observed range on the model's scale, fixed-effect leverage relative to the
  data, and the share of the spatial variance not explained by the sampled
  locations. Locations are classed as supported, spatial, covariate, or
  spatial and covariate extrapolation, and `plot()` maps each indicator.
  Locations with missing covariates are drawn in grey without a legend key.
  The legend shows the colour of every class, including classes absent from
  the map, whose keys ggplot2 otherwise leaves blank.
- `create_grid()` gains `raster`. The grid then consists of the cells of a
  projected `SpatRaster` with square cells: their centres are the grid points
  and the raster's layers become columns, so values such as population counts
  aggregated to the grid's resolution need no extraction. `boundaries` is
  optional in that case and, when given, selects cells according to `keep`.
- New `get_worldpop()` downloads WorldPop population counts for a country
  and year, for the whole population or for chosen age groups and sexes, at
  100 m or 1 km, from the unconstrained 2000-2020 or constrained 2015-2030
  estimates, and returns one raster of people per cell.
- `to_table()` is now generic and gains a method for `predict()` results with
  polygon summaries: one row per polygon, the mean and its interval in one
  column, probabilities as percentages and sums with thousands separators.
  All `to_table()` methods accept `col_names`. `predict()` output records the
  name of the `weights` column.
- `to_table()` for predictions formats each sum on its own, so that values in a
  column are no longer padded to a common number of decimals.
- Documentation and messages call phi the scale parameter of the spatial
  correlation function, rather than the spatial range, consistently with
  `summary()` and the MBG book. The practical range, where correlation falls
  to 0.05, is about 3 phi for the exponential correlation function.

## Cross validation

- `"user"` has been added as an option to `method` the parameter of
`assess_prediction()` replacing the previous behaviour where providing the
`user_split` parameter overrode any provided `method`
- In `assess_prediction()`, `n_size` has been renamed to `size` (for consistency
with `iter` and `fold`) and `which_metric` has been renamed to `metrics`.
Parameters have also been reordered: mandatory arguments first, then those
required per `method` (in the order `"cluster"`, `"regularized"`, `"user"`),
then the remaining optional arguments.
- `assess_prediction()` takes a `seed` that makes the random splits generated
for `method = "cluster"` or `"regularized"` reproducible, together with the
refits and predictions of each test set.
- `assess_prediction()`'s `plot_fold = TRUE` plots are now consistent across
`method`s: `"cluster"` is unchanged (`spatialsample::autoplot()`), while
`"user"` and `"regularized"` now plot the testing and training datasets.
- `assess_simulation()` now takes fitted objects returned by `glgpm()` or
  model specifications returned by `specify_glgpm()` as its `models`
  templates. `f_grid_target` and `f_area_target` have been renamed to
  `target_transform` and `area_summary`; their returned values are now checked
  before they are used. Grid- and area-level objectives can be requested
  together in one run with `spatial_scale = c("grid", "area")`.
- In `plot_metric()`; its `which_score`/`which_model`
parameters have been renamed to `metric`/`model`.

## Simulations

- `simulate_glgpm()` now accepts a fitted model or `specify_glgpm()` model,
  with `nsim`, `what`, `sample_locations`, `prediction_grid` and `seed`.
  `what = c("data", "surface")` simulates jointly at the exact union of both
  sets of locations. This replaces `simulate_surface()` and its nearest-grid
  approximation. The old simulation arguments and output structure are removed.
  Parameters remain fixed, offsets and fitted inverse links are retained, and
  Poisson responses use exposure times the inverse-link mean.
  Use `simulated_data()`, `simulated_surface()` and `simulated_values()` to
  extract results.

## Other

- Argument checks in the compiled MALA, likelihood and distance kernels have
been removed or moved to their R callers, following the contributing guide;
the R functions already validated everything the MALA and likelihood kernels
checked.
- `control_sim` parameters have been renamed to `control_mcmc`.
- `to_table()` now returns a directly renderable `knitr_kable` table, with an
  explicit `digits` argument that preserves trailing zeroes. This replaces the
  previous `xtable` return value. Values that are compared share one
  notation: when a value needs scientific notation, so do the estimate and
  limits of the same parameter, or the same metric for the other models.
- `summary()` of a model with a nugget gave the wrong confidence interval for
  the nugget variance: its standard error was that of the logarithm of
  tau2 / sigma2^2 instead of tau2. The interval now uses the delta method for
  log(tau2) = log(tau2 / sigma2) + log(sigma2).
- When one parameter has no standard error, as when a variance is estimated
  at zero, `summary()` and `to_table()` still give the confidence intervals of
  the other parameters. Previously all of them were missing.
- The `data_sf` parameter in `create_ids` is now `data`.
- The `conf_level` parameter in `summary.RiskMap` is now `level`
- `plot_metric()` and `plot_AnPIT()` have been merged into a single
  `plot.RiskMap_cross_validation()` method. `metric` and `model` now default
  to every metric/model available rather than being required, so
  `plot(cv_result)` alone produces every plot the object has data for.
  `metric` accepts `"PIT"` as well as `"CRPS"`/`"SCRPS"`/`"PIT_area"` to
  select the calibration curve. Each metric gets its own plot (or, with more
  than one model, its own small grid with one panel per model); with more
  than one metric, each one is drawn in turn and the full set is returned
  invisibly as a named list.
- `plot_variogram()` is now `plot.RiskMap_variogram()`, a method rather than
  a standalone function; call it as `plot(variogram_result)`.
- `plot_sim_surf()` is now `plot.RiskMap_simulation()`; call it as
  `plot(sim_result, simulation = ...)`. It now draws a `ggplot` map matching
  the prediction maps, with `palette` and `reverse_palette`
  arguments. `simulated_surface()` returns an object of class
  `RiskMap_simulated_surface` with a matching `plot()` method.
- In `assess_prediction()`, the `"AnPIT"` metric is now `"PIT"` and the
  `AnPIT_area` score is now `PIT_area`, for every family.
- In `plot.RiskMap_cross_validation()`, the calibration-curve arguments `mode`,
  `test_set` and `combine_panels` are now `pit_mode`, `pit_test_set` and
  `combine_pit`.
- Test coverage has been increased from 0 to 85 %.
