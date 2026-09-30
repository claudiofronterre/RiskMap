RiskMap 2.0.0
=============

- This version introduces many breaking changes from v1. The package should not
be considered stable, but further breaking changes will be handled gracefully.
- `dast()` has been removed.
- Many functions have been renamed:

| v1 | v2 | 
|----------|----------|
| `assess_pp` |  `assess_prediction` |
| `assess_sim` | `assess_simulation` |
| `check_mcmc` | `plot_mcmc` |
| `compute_ID_coords` | `create_ids` |
| `convex_hull_sf` | `create_convex_hull` |
| `glgpm_sim` | `simulate_glgpm` |
| `Laplace_sampling_MCMC` | `laplace_sampling_mcmc` |
| `matern.grad.phi` | `matern_gradient_phi ` |
| `matern.hessian.phi` | `matern_hessian_phi` |
| `matern_cor` | `matern_correlation` |
| `maxim.integrand` | `maxim_integrand` |
| `plot_score` | `plot_metric` |
| `plot_s_variogram` | `plot_variogram` |
| `pred_over_grid` | `setup_prediction` |
| `pred_target_grid` | `predict_grid_target` |
| `pred_target_shp` | `predict_areal_target` |
| `set_control_sim` | `set_control_mcmc` |
| `surf_sim` | `simulate_glgpm` |
| `s_variogram` | `variogram` |

- Apart from `liberia`, all datasets are now in an sf format.
- The `den` argument in `glgpm()` has been renamed to `denominator`.
- `glpgm()` now only accepts data in an sf format. 
Consequently the locations do not need to be passed to `gp()` when fitting a model
as they are included automatically.
- In `glgpm()`, `convert_to_crs` has been replaced with `model_crs` and
`scale_to_km` has been replaced with `distance_units` (`"km"` or `"m"`).
If `data` are in longitude/latitude and `model_crs` is not supplied, the
coordinates are now automatically reprojected to an appropriate UTM zone;
providing `model_crs` in longitude/latitude now raises an error.
`setup_prediction()`, `simulate_glgpm()` and `assess_prediction()` have
been updated to match; `simulate_glgpm()`'s own `convert_to_crs`/`scale_to_km`
arguments have likewise been replaced with `model_crs`/`distance_units`.
- In `summarise_distance()` and `variogram()`, `scale_to_km` has similarly been
replaced with `distance_units` (`"km"` or `"m"`), for consistency with `glgpm()`.
Default behaviour is unchanged (`"km"` for `summarise_distance()`, `"m"` for
`variogram()`).
- In `summarise_distance()` and `variogram()`, `convert_to_utm` has been replaced
with `distance_crs`, matching `glgpm()`'s `model_crs` convention: an
already-projected CRS is retained as-is, longitude/latitude data is automatically
reprojected to an appropriate UTM zone with a message, and an explicitly supplied
`distance_crs` must be projected. Distances now correctly account for the CRS's
native linear units rather than assuming metres.
- `plot_variogram()`'s `plot_envelope` now defaults to `TRUE`; if `n_permutations`
was 0 or 1 in the `variogram()` call, a warning is raised and the envelope is
skipped instead of erroring.
- `create_grid()` no longer raises an error when `boundaries` is in
longitude/latitude; like `glgpm()`, it now automatically reprojects to an
appropriate UTM zone and reports the conversion with a message.
- The `shp` parameter has been renamed to `boundaries` in `create_grid()`,
`predict_areal_target()` and `assess_simulation()`; `shp_target` and `return_shp`
in `predict_areal_target()` have similarly been renamed to `areal_target` and
`return_boundaries`. `predict_areal_target()`'s returned `shp` component is now
named `boundaries`.
- The `bins` parameter in `variogram()` has been removed and replaced with `breaks`.
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
- In `plot.RiskMap_predict_grid_target()` and `plot.RiskMap_predict_areal_target()`,
`which_target`/`which_summary` have been renamed to `target`/`summary`, for
consistency elsewhere in the package. `target` now defaults to `NULL`, which
plots the first available target; an unrecognised `target` or `summary` now
raises an informative error instead of a blank plot.
- In `plot_metric()`; its `which_score`/`which_model`
parameters have been renamed to `metric`/`model`.
- The `nugget` parameter in `gp()` is now `FALSE` by default and should be set to
`TRUE` to estimate it.
- A `seed` parameter can be passed to `set_control_mcmc()` to make non-Gaussian
outputs reproducible.
- `control_sim` parameters have been renamed to `control_mcmc`.
- `to_table()` now returns a directly renderable `knitr_kable` table, with an
  explicit `digits` argument that preserves trailing zeroes. This replaces the
  previous `xtable` return value.
- `simulate_glgpm()` now accepts a fitted model or `specify_glgpm()` model,
  with `nsim`, `what`, `sample_locations`, `prediction_grid` and `seed`.
  `what = c("data", "surface")` simulates jointly at the exact union of both
  sets of locations. This replaces `simulate_surface()` and its nearest-grid
  approximation. The old simulation arguments and output structure are removed.
  Parameters remain fixed, offsets and fitted inverse links are retained, and
  Poisson responses use exposure times the inverse-link mean.
  Use `simulated_data()`, `simulated_surface()` and `simulated_values()` to
  extract results.
- Gaussian prediction uses a stable reduced covariance solve, avoiding cancellation
with small measurement-error variances while retaining the same conditional model.
- Test coverage has been increased from 0 to 70 %.
