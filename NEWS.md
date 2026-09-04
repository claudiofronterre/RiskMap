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
- The `bins` parameter in `variogram()` has been removed and replaced with `breaks`.
- The `nugget` parameter in `gp()` is now `FALSE` by default.
- A `seed` parameter can be passed to `set_control_mcmc()` to make non-Gaussian outputs reproducible.
- Test coverage has been increased from 0 to 70 %.
