# Changelog

## RiskMap 2.0.0

RiskMap 2.0.0 is a redesign of the package. Almost every function of
1.0.0 was renamed, merged or removed, and no name of 1.0.0 remains as a
stub, so code written for 1.0.0 needs changing. The table gives the
replacement of each function of 1.0.0, and the main changes follow. The
full record of changes is in
[`notes/changes-2.0.md`](https://github.com/claudiofronterre/RiskMap/blob/main/notes/changes-2.0.md)
in the source repository.

### Functions of 1.0.0 and their replacements

| 1.0.0 | 2.0.0 |
|----|----|
| `glgpm(den, crs, convert_to_crs, scale_to_km, par0, S_samples, fix_var_me)` | `glgpm(denominator, distance_crs, distance_units)`. `fixed_pars` replaces `par0` and `fix_var_me`, `start_pars` takes any parameters by name, and `S_samples` is gone |
| `gp(nugget = 0)` | `gp(nugget = FALSE)`, or `TRUE` for a nugget, which `fixed_pars = list(tau2 = )` can hold |
| `set_control_sim()` | `set_control_mcmc(n_draws, burnin, thin, h, c1_h, c2_h)` for sampling and [`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md) for fitting |
| `summary(fit, conf_level)` | `summary(fit, level)` |
| [`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md) | [`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md), which returns a `knitr_kable` table instead of an `xtable` |
| `pred_over_grid()`, `pred_target_grid()`, `update_predictors()` | [`predict()`](https://rdrr.io/r/stats/predict.html), with every predictor in `newdata` |
| `pred_target_shp()` | `predict(boundaries, aggregate)` |
| [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of `pred_target_grid()` and `pred_target_shp()` results | [`plot_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot_prediction.md), or [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a [`predict()`](https://rdrr.io/r/stats/predict.html) result |
| `assess_pp()` | [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md), which takes `models`, `size`, `metrics` and `control_mcmc` for `object`, `n_size`, `which_metric` and `control_sim`, and `"PIT"` for `"AnPIT"` |
| `plot_score()`, `plot_AnPIT()` | [`plot_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_cross_validation.md), or [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of an [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md) result |
| `assess_sim()` | [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md), which takes `simulation`, `target`, `aggregate`, `boundaries`, `label`, `objective` and `breaks` for `obj_sim`, `f_grid_target`, `f_area_target`, `shp`, `col_names`, `pred_objective` and `categories`, and `spatial_scale = "location"` for `"grid"` |
| [`summary()`](https://rdrr.io/r/base/summary.html) of `assess_pp()` and `assess_sim()` results | [`summary()`](https://rdrr.io/r/base/summary.html) of [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md) and [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md) results |
| `glgpm_sim()`, `surf_sim()` | [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md), with [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md) for a model given by its parameters and `what = "surface"` for surfaces |
| `plot_sim_surf()` | [`plot_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_simulation.md), or [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md) result |
| `s_variogram()`, `plot_s_variogram()` | [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md) and [`plot_variogram()`](https://claudiofronterre.github.io/RiskMap/reference/plot_variogram.md), with `breaks` for `bins`, `n_permutations` for `n_permutation`, and `distance_crs` and `distance_units` for `convert_to_utm` and `scale_to_km` |
| `dist_summaries()` | [`summarise_distance()`](https://claudiofronterre.github.io/RiskMap/reference/summarise_distance.md), with `distance_crs` and `distance_units` |
| `propose_utm()` | [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md), which chooses the CRS by its distortion of distances |
| `convex_hull_sf()` | [`create_convex_hull()`](https://claudiofronterre.github.io/RiskMap/reference/create_convex_hull.md) |
| `create_grid(shp, spat_res, grid_crs)` | `create_grid(boundaries, spacing, distance_crs)` |
| `check_mcmc()` | [`plot_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/plot_mcmc.md), with `locations` for `component`, and [`summarise_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/summarise_mcmc.md) |
| `matern_cor()` | [`matern_correlation()`](https://claudiofronterre.github.io/RiskMap/reference/matern_correlation.md) |
| `dast()`, `plot_mda()` | Removed. Models of mass drug administration move to a separate package |
| `Laplace_sampling_MCMC()`, `maxim.integrand()`, `matern.grad.phi()`, `matern.hessian.phi()`, `compute_ID_coords()` | Removed. Their work is internal |

### Main changes

#### Fitting

- [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  fits Gaussian, binomial, Poisson and negative binomial models with a
  spatial Gaussian process,
  [`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md),
  and random effects,
  [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md),
  and takes any link. Binomial, Poisson and negative binomial models are
  fitted by certified Monte Carlo maximum likelihood, which enlarges its
  samples until the Monte Carlo error of the estimates is small and
  reports a bound on that error. In binomial and Poisson fits, household
  effects whose children differ in covariates are integrated out, which
  makes such fits much faster.
- Every parameter has one name, the same in `start_pars`, `fixed_pars`,
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`profile()`](https://rdrr.io/r/stats/profile.html) and
  [`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md),
  with `tau2` for the variance of the nugget.
  [`glgpm_parameters()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm_parameters.md)
  lists the parameters of a model and
  [`?riskmap_parameters`](https://claudiofronterre.github.io/RiskMap/reference/riskmap_parameters.md)
  describes them. `start_pars` starts any of them, so
  `start_pars = coef(fit)` starts from another fit, and `fixed_pars`
  holds any of them, or all, in which case
  [`predict()`](https://rdrr.io/r/stats/predict.html) works with the
  given values.
- Data are `sf` objects. `distance_crs` and `distance_units` set the
  projection and the units of distances in every function that measures
  them, and
  [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
  proposes a projection that keeps distances accurate.
- [`?glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  lists the components of a fit that code may rely on.
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) work on fits, and with
  them [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`update()`](https://rdrr.io/r/stats/update.html). Every result
  records the version of RiskMap that made it.
- [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md)
  sets the number of draws it keeps in `n_draws`, and
  [`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md)
  takes its sample sizes as `initial_draws`, `validation_draws` and
  `max_draws`.

#### Checking a fit

- [`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md)
  reports how a fit was computed and flags anything that needs
  attention, [`profile()`](https://rdrr.io/r/stats/profile.html) gives
  profile likelihoods and
  [`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)
  bootstrap intervals for the parameters.

#### Prediction

- A single [`predict()`](https://rdrr.io/r/stats/predict.html) method
  predicts at new locations or over the areas given in `boundaries`,
  with exceedance probabilities, population weights, standardisation
  and, if asked, the uncertainty of the estimated parameters.
  [`plot_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot_prediction.md)
  maps the predictions,
  [`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)
  tabulates them, and
  [`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md)
  shows where they go beyond the data.
- The covariates of the prediction grid `tz_covariates` now match those
  of `tz_malaria`. EVI, temperature, night-time light and precipitation
  come from the same rasters, so temperature is missing at 2 points and
  precipitation at 1, and ITN is a percentage with `NA` for missing
  values, where it was a proportion with -9999. The help pages give the
  source and units of each covariate, and population density is not in
  thousands, as they said.

#### Assessment and simulation

- [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
  compares models by spatial cross-validation. Every score is a loss, so
  lower is better, and `PIT_area` is given once for all test
  observations, next to its value under perfect calibration. Counts are
  scored with their exact distribution given the draws of the mean,
  where `assess_pp()` simulated them.
  [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
  and
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
  simulate data and surfaces from a model, and
  [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md)
  assesses predictions on them.
- In
  [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md)
  every value belongs to a class, those beyond the `breaks` to the end
  classes, covariates can have any name, and a prediction at a sampled
  location uses what the data say about its nugget.

#### Plots, printing and messages

- Every plot has a function of its own, such as
  [`plot_variogram()`](https://claudiofronterre.github.io/RiskMap/reference/plot_variogram.md),
  [`plot_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot_prediction.md)
  and
  [`plot_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_cross_validation.md),
  which [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the
  result also calls, so users can call either. Each returns one ggplot
  or patchwork object, which is drawn when printed. Methods and plot
  functions stop on an argument they don’t use.
- Errors, warnings and messages are written with cli and have classes
  ([`?riskmap_conditions`](https://claudiofronterre.github.io/RiskMap/reference/riskmap_conditions.md)),
  long computations show their progress
  ([`?riskmap_progress`](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md)),
  and `options(RiskMap.messages = FALSE)` silences the messages.
- Packages that only some functions need, such as spatialsample and
  spatialEco for some splits of
  [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md),
  ggspatial for north arrows and scale bars, and jsonlite for
  [`get_worldpop()`](https://claudiofronterre.github.io/RiskMap/reference/get_worldpop.md),
  are suggested rather than imported, and those functions ask for them.
- RiskMap needs R 4.5.0 or later.
- Custom links, the negative binomial family, the penalised complexity
  priors of
  [`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md),
  [`crs_distortion()`](https://claudiofronterre.github.io/RiskMap/reference/crs_distortion.md)
  and the spatial indicator of
  [`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md)
  are experimental and may change in a minor release, as their help
  pages say.
