# RiskMap 2.0.0

RiskMap 2.0.0 is a redesign of the package. Almost every function of 1.0.0 was
renamed, merged or removed, and no name of 1.0.0 remains as a stub, so code
written for 1.0.0 needs changing. The table gives the replacement of each
function of 1.0.0, and the main changes follow. The full record of changes is
in [`notes/changes-2.0.md`](https://github.com/claudiofronterre/RiskMap/blob/main/notes/changes-2.0.md)
in the source repository.

## Functions of 1.0.0 and their replacements

| 1.0.0 | 2.0.0 |
|---|---|
| `glgpm(den, crs, convert_to_crs, scale_to_km, par0, S_samples, fix_var_me)` | `glgpm(denominator, distance_crs, distance_units)`. `fixed_pars` replaces `par0` and `fix_var_me`, `start_pars` takes any parameters by name, and `S_samples` is gone |
| `gp(nugget = 0)` | `gp(nugget = FALSE)`, or `TRUE` for a nugget, which `fixed_pars = list(tau2 = )` can hold |
| `set_control_sim()` | `set_control_mcmc(n_draws, burnin, thin, h, c1_h, c2_h)` for sampling and `set_control_mcml()` for fitting |
| `summary(fit, conf_level)` | `summary(fit, level)` |
| `to_table()` | `to_table()`, which returns a `knitr_kable` table instead of an `xtable` |
| `pred_over_grid()`, `pred_target_grid()`, `update_predictors()` | `predict()`, with every predictor in `newdata` |
| `pred_target_shp()` | `predict(boundaries, aggregate)` |
| `plot()` of `pred_target_grid()` and `pred_target_shp()` results | `plot_prediction()`, or `plot()` of a `predict()` result |
| `assess_pp()` | `assess_prediction()`, which takes `models`, `size`, `metrics` and `control_mcmc` for `object`, `n_size`, `which_metric` and `control_sim`, and `"PIT"` for `"AnPIT"` |
| `plot_score()`, `plot_AnPIT()` | `plot_cross_validation()`, or `plot()` of an `assess_prediction()` result |
| `assess_sim()` | `assess_simulation()`, which takes `simulation`, `target`, `aggregate`, `boundaries`, `label`, `objective` and `breaks` for `obj_sim`, `f_grid_target`, `f_area_target`, `shp`, `col_names`, `pred_objective` and `categories`, and `spatial_scale = "location"` for `"grid"` |
| `summary()` of `assess_pp()` and `assess_sim()` results | `summary()` of `assess_prediction()` and `assess_simulation()` results |
| `glgpm_sim()`, `surf_sim()` | `simulate_glgpm()`, with `specify_glgpm()` for a model given by its parameters and `what = "surface"` for surfaces |
| `plot_sim_surf()` | `plot_simulation()`, or `plot()` of a `simulate_glgpm()` result |
| `s_variogram()`, `plot_s_variogram()` | `variogram()` and `plot_variogram()`, with `breaks` for `bins`, `n_permutations` for `n_permutation`, and `distance_crs` and `distance_units` for `convert_to_utm` and `scale_to_km` |
| `dist_summaries()` | `summarise_distance()`, with `distance_crs` and `distance_units` |
| `propose_utm()` | `propose_crs()`, which chooses the CRS by its distortion of distances |
| `convex_hull_sf()` | `create_convex_hull()` |
| `create_grid(shp, spat_res, grid_crs)` | `create_grid(boundaries, spacing, distance_crs)` |
| `check_mcmc()` | `plot_mcmc()`, with `locations` for `component`, and `summarise_mcmc()` |
| `matern_cor()` | `matern_correlation()` |
| `dast()`, `plot_mda()` | Removed. Models of mass drug administration move to a separate package |
| `Laplace_sampling_MCMC()`, `maxim.integrand()`, `matern.grad.phi()`, `matern.hessian.phi()`, `compute_ID_coords()` | Removed. Their work is internal |

## Main changes

### Fitting

- `glgpm()` fits Gaussian, binomial, Poisson and negative binomial models with
  a spatial Gaussian process, `gp()`, and random effects, `re()`, and takes
  any link. Binomial, Poisson and negative binomial models are fitted by
  certified Monte Carlo maximum likelihood, which enlarges its samples until
  the Monte Carlo error of the estimates is small and reports a bound on that
  error. In binomial and Poisson fits, household effects whose children differ
  in covariates are integrated out, which makes such fits much faster.
- Every parameter has one name, the same in `start_pars`, `fixed_pars`,
  `coef()`, `profile()` and `parametric_bootstrap()`, with `tau2` for the
  variance of the nugget. `glgpm_parameters()` lists the parameters of a model
  and `?riskmap_parameters` describes them. `start_pars` starts any of them,
  so `start_pars = coef(fit)` starts from another fit, and `fixed_pars` holds
  any of them, or all, in which case `predict()` works with the given values.
- Data are `sf` objects. `distance_crs` and `distance_units` set the
  projection and the units of distances in every function that measures them,
  and `propose_crs()` proposes a projection that keeps distances accurate.
- `?glgpm` lists the components of a fit that code may rely on. `logLik()`,
  `vcov()` and `nobs()` work on fits, and with them `AIC()`, `BIC()` and
  `update()`. Every result records the version of RiskMap that made it.
- `set_control_mcmc()` sets the number of draws it keeps in `n_draws`, and
  `set_control_mcml()` takes its sample sizes as `initial_draws`,
  `validation_draws` and `max_draws`.

### Checking a fit

- `check_fit()` reports how a fit was computed and flags anything that needs
  attention, `profile()` gives profile likelihoods and
  `parametric_bootstrap()` bootstrap intervals for the parameters.

### Prediction

- A single `predict()` method predicts at new locations or over the areas
  given in `boundaries`, with exceedance probabilities, population weights,
  standardisation and, if asked, the uncertainty of the estimated parameters.
  `plot_prediction()` maps the predictions, `to_table()` tabulates them, and
  `assess_extrapolation()` shows where they go beyond the data.
- The covariates of the prediction grid `tz_covariates` now match those of
  `tz_malaria`. Temperature and night-time light come from the same rasters,
  and ITN is a percentage with `NA` for missing values, where it was a
  proportion with -9999.

### Assessment and simulation

- `assess_prediction()` compares models by spatial cross-validation. Every
  score is a loss, so lower is better, and `PIT_area` is given once for all
  test observations, next to its value under perfect calibration.
  `specify_glgpm()` and `simulate_glgpm()` simulate data and surfaces from a
  model, and `assess_simulation()` assesses predictions on them.

### Plots, printing and messages

- Every plot has a function of its own, such as `plot_variogram()`,
  `plot_prediction()` and `plot_cross_validation()`, which `plot()` of the
  result also calls, so users can call either. Each returns one ggplot or
  patchwork object, which is drawn when printed. Methods and plot functions
  stop on an argument they don't use.
- Errors, warnings and messages are written with cli and have classes
  (`?riskmap_conditions`), long computations show their progress
  (`?riskmap_progress`), and `options(RiskMap.messages = FALSE)` silences the
  messages.
- Packages that only some functions need, such as spatialsample and spatialEco
  for some splits of `assess_prediction()`, ggspatial for north arrows and
  scale bars, and jsonlite for `get_worldpop()`, are suggested rather than
  imported, and those functions ask for them.
- Custom links, the negative binomial family, the penalised complexity priors
  of `set_control_mcml()`, `crs_distortion()` and the spatial indicator of
  `assess_extrapolation()` are experimental and may change in a minor release,
  as their help pages say.
