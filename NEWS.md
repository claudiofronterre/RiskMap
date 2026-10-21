RiskMap 2.0.0
=============

- This version introduces many breaking changes from v1. The package should not
be considered stable, but further breaking changes will be handled gracefully.
- Test coverage has been increased from 0 to 92 %.
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
| `plot_score` | `plot_metric` |
| `plot_s_variogram` | `plot_variogram` |
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

- In `summarise_distance()` and `variogram()`, `scale_to_km` has been
replaced with `distance_units` (`"km"` or `"m"`), for consistency with `glgpm()`.
Default behaviour is unchanged (`"km"` for `summarise_distance()`, `"m"` for
`variogram()`).
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
  `plot_variogram()` now shows the envelope by default and warns when it
  was omitted with `n_permutations = 0`.
- The `bins` parameter in `variogram()` has been removed and replaced with `breaks`.

## Model fitting

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

- Non-Gaussian optimisation now reuses the sample-wise Monte Carlo log
integrand when the objective, gradient, and Hessian are evaluated at identical
parameters, avoiding a repeated likelihood calculation.

- Non-Gaussian models now use automatically validated local-Laplace MCML.
Independent fitting and validation banks must contain at least 100 effective
draws, have acceptable importance-weight tails, agree on the score, and keep
Monte Carlo standard error below 10% of statistical standard error in every
parameter direction. Bank size increases only when those diagnostics require
it. The defaults were chosen in a 1,200-dataset simulation study and checked
against exact likelihoods on small datasets, Stan-based likelihood ratios and
1,200 further datasets simulated with fresh seeds; see #194.

- Adaptive MCML now distinguishes Monte Carlo error from weak
identification. A variance component estimated as zero (and the spatial range
when the spatial variance is zero) is held at the boundary, the first stage is
restarted once from an interior value, and precision is assessed on the
remaining parameters, which keep their standard errors. A spatial range shorter
than the median distance between neighbouring locations is held at that
distance (`set_control_mcml(hold_short_range = FALSE)` disables this). If the
information matrix is still singular, fitting stops instead of drawing more
samples. Each case produces a plain-language warning.

- `maxim_integrand()` now takes `denominator` instead of `units_m`, matching
`glgpm()`. Its documentation no longer lists the Gaussian family, which it has
never supported, and describes the Newton-type optimisation it actually uses.

- Removed unused internal helpers (`log_mean_exp()`,
`importance_effective_sample_size()`, `restore_spatial_range()`) and the
compiled kernel `cpp_canonical_likelihood()`, which no longer had callers.

- `set_control_mcml()` gains `pc_prior` (default `FALSE`) and its thresholds
`pc_sigma_threshold`, `pc_sigma_probability`, `pc_range_threshold` and
`pc_range_probability`. When `pc_prior = TRUE`, the spatial variance and range
are penalised by penalised-complexity priors (by default
P(sigma > 1) = 0.05 and P(range < median nearest-neighbour distance) = 0.05),
which stabilises estimation in small or low-prevalence surveys. Warnings about
a spatial variance estimated as zero, a held spatial range or a singular
information matrix now suggest this option.

- `set_control_mcml()` gains `laplace_start` (default `FALSE`) to maximise the
Laplace approximation from several starting values before MCML.

## Predictions

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

## Cross validation

- `"user"` has been added as an option to `method` the parameter of
`assess_prediction()` replacing the previous behaviour where providing the
`user_split` parameter overrode any provided `method`
- In `assess_prediction()`, `n_size` has been renamed to `size` (for consistency
with `iter` and `fold`) and `which_metric` has been renamed to `metrics`.
Parameters have also been reordered: mandatory arguments first, then those
required per `method` (in the order `"cluster"`, `"regularized"`, `"user"`),
then the remaining optional arguments.
- `assess_prediction()` now honours a `seed` set on its `control_mcmc` argument
to make the random splits generated for `method = "cluster"` or `"regularized"`
reproducible.
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

- A `seed` parameter can be passed to `set_control_mcmc()` to make non-Gaussian
outputs reproducible.
- `control_sim` parameters have been renamed to `control_mcmc`.
- `to_table()` now returns a directly renderable `knitr_kable` table, with an
  explicit `digits` argument that preserves trailing zeroes. This replaces the
  previous `xtable` return value.
- The `data_sf` parameter in `create_ids` is now `data`.
- The `conf_level` parameter in `summary.RiskMap` is now `level`
