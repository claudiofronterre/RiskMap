RiskMap 2.0.0
=============

- This version introduces many breaking changes from v1. The package should not be considered stable, but further breaking changes will be handled gracefully.
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
| `plot_s_variogram` | `plot_variogram` |
| `pred_over_grid` | `setup_prediction` |
| `pred_target_grid` | `predict_grid_target` |
| `pred_target_shp` | `predict_areal_target` |
| `set_control_sim` | `set_control_mcmc` |
| `surf_sim` | `simulate_surface` |
| `s_variogram` | `variogram` |

- Apart from `liberia`, all datasets are now in an sf format.
- `glpgm()` now only accepts data in an sf format. 
Consequently the locations do not need to be passed to `gp()` when fitting a model as they are included automatically.
- In `glgpm()`, `convert_to_crs` has been replaced with `model_crs` and `scale_to_km` has been replaced with `coordinate_units` (`"km"` or `"m"`).
If `data` are in longitude/latitude and `model_crs` is not supplied, the coordinates are now automatically reprojected to an appropriate UTM zone; providing `model_crs` in longitude/latitude now raises an error.
`setup_prediction()`, `simulate_glgpm()` and `assess_prediction()` have been updated to match; `simulate_glgpm()`'s own `convert_to_crs`/`scale_to_km` arguments have likewise been replaced with `model_crs`/`coordinate_units`.
- The `bins` parameter in `variogram()` has been removed and replaced with `breaks`.
- The `nugget` parameter in `gp()` is now `FALSE` by default.
- A `seed` parameter can be passed to `set_control_mcmc()` to make non-Gaussian outputs reproducible.
- Test coverage has been increased from 0 to 70 %.
