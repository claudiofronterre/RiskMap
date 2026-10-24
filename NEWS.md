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
  `plot_variogram()` (now `plot.RiskMap_variogram()`) shows the envelope by
  default and warns when it was omitted with `n_permutations = 0`.
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

- Non-Gaussian models are now fitted by certified Monte Carlo maximum
likelihood. The estimated likelihood is maximised on a fitting sample of 1,000
draws, the estimate is corrected by one Newton step with the score from an
independent validation sample of 5,000 draws, and the corrected estimate is
certified with a second, fresh validation sample. The certificate bounds the
distance from the exact maximum likelihood estimate in statistical standard
errors, with probability about 0.95. When it exceeds the tolerance (0.1 by
default), both samples are enlarged to the size predicted to meet it, up to
`max_samples`. The importance distribution is a Gaussian approximation fitted
to the whole conditional distribution of the latent effects (efficient
importance sampling), which keeps far more draws effective than the Laplace
approximation in sparse data. Binomial and Poisson models, with or without a
nugget and unstructured random effects, share one implementation, and the
Laplace approximation remains as a fallback and as an option.
`set_control_mcml()` gains `validation_samples`, `tolerance` and `proposal`,
and loses `min_effective_samples`, `max_pareto_k`, `max_mcse_ratio` and
`score_probability`. The defaults of `max_samples` and `max_stages` are now
16,000 and 4. Fitted objects report `mcml_certified` in place of
`mcml_converged`.

- Certified MCML distinguishes Monte Carlo error from weak
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

- The final conditional sampling in `glgpm()` (used for fitted values and
stored prediction samples) now uses the compiled MALA kernel for canonical
links, as intended. It previously passed the inverse-link function itself,
which selected the general R implementation with numerical derivatives; a
binomial fit to the `liberia` data now takes about 2 s instead of 27 s.
Fresh conditional runs in `predict()` now honour a custom inverse link.

- `set_control_mcml()` gains `laplace_start` (default `TRUE`), which maximises
the Laplace approximation from the starting values before MCML.

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

- Argument checks in the compiled MALA, likelihood and distance kernels have
been removed or moved to their R callers, following the contributing guide;
the R functions already validated everything the MALA and likelihood kernels
checked.
- A `seed` parameter can be passed to `set_control_mcmc()` to make non-Gaussian
outputs reproducible.
- `control_sim` parameters have been renamed to `control_mcmc`.
- `to_table()` now returns a directly renderable `knitr_kable` table, with an
  explicit `digits` argument that preserves trailing zeroes. This replaces the
  previous `xtable` return value.
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
