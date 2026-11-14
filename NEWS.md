# RiskMap 2.0.0

RiskMap 2.0.0 is a redesign of the package. Many functions of 1.0.0 have been
renamed, merged or removed, so code written for 1.0.0 needs changing. The table
gives the replacement of each function of 1.0.0, and the main changes follow.

## Functions of 1.0.0 and their replacements

| 1.0.0 | 2.0.0 |
|---|---|
| `assess_pp()` | `assess_prediction()` |
| `assess_sim()` | `assess_simulation()` |
| `check_mcmc()` | `plot_mcmc()` |
| `compute_ID_coords()` | `create_ids()` |
| `convex_hull_sf()` | `create_convex_hull()` |
| `dast()` | Removed, to move to a separate package |
| `dist_summaries()` | `summarise_distance()` |
| `glgpm_sim()` | `simulate_glgpm()` |
| `Laplace_sampling_MCMC()` | Removed, now internal |
| `matern.grad.phi()` | `matern_gradient_phi()` |
| `matern.hessian.phi()` | `matern_hessian_phi()` |
| `matern_cor()` | `matern_correlation()` |
| `maxim.integrand()` | `maxim_integrand()` |
| `plot_AnPIT()`, `plot_score()` | `plot()` of `assess_prediction()` results |
| `plot_mda()` | Removed, to move with `dast()` |
| `plot_s_variogram()` | `plot()` of `variogram()` results |
| `plot_sim_surf()` | `plot()` of `simulate_glgpm()` results |
| `pred_over_grid()`, `pred_target_grid()` | `predict()` |
| `pred_target_shp()` | `predict()` with `boundaries` |
| `propose_utm()` | `propose_crs()` |
| `s_variogram()` | `variogram()` |
| `set_control_sim()` | `set_control_mcmc()` |
| `surf_sim()` | `simulate_glgpm()` |
| `update_predictors()` | Removed, give the predictors in the `newdata` of `predict()` |

## Main changes

### Model fitting

- `glgpm()` fits Gaussian, binomial, Poisson and negative binomial models. A
  link can be named in a family object, such as `binomial("probit")`, and any
  other link given as an inverse link function.
- Non-Gaussian models are fitted by certified Monte Carlo maximum likelihood,
  which enlarges its samples until the Monte Carlo error of the estimates is
  small and reports a bound on that error.
- In binomial and Poisson fits, household effects whose children differ in
  covariates are now integrated out, so such fits are much faster.
- `check_fit()` reports how a fit was computed and flags anything that needs
  attention. `profile()` gives profile likelihoods and `parametric_bootstrap()`
  bootstrap intervals for the parameters.

### Prediction

- A single `predict()` method replaces the three-stage interface of 1.0.0. It
  predicts at new locations or over the areas given in `boundaries`, with
  exceedance probabilities, population weights, standardisation and, if asked,
  the uncertainty of the estimated parameters.
- `plot()` maps predictions, `to_table()` tabulates predictions over areas,
  and `assess_extrapolation()` shows where predictions go beyond the data.

### Data and projections

- Data are `sf` objects with a coordinate reference system. `propose_crs()`
  proposes a projection that keeps distances accurate, and `glgpm()` uses it
  for data in longitude and latitude.
- `create_grid()` builds prediction grids, also from the cells of a raster,
  `get_worldpop()` downloads WorldPop population counts, and `variogram()`
  gives the empirical variogram with a permutation envelope.

### Assessment and simulation

- `assess_prediction()` compares models by spatial cross-validation.
- `specify_glgpm()` and `simulate_glgpm()` simulate data and surfaces from a
  model, and `assess_simulation()` assesses prediction targets on them.

### Messages

- Errors, warnings and messages are written with cli and have classes, listed
  in `?riskmap_conditions`, and long computations show their progress
  (`?riskmap_progress`). `options(RiskMap.messages = FALSE)` silences the
  messages.
