# Changelog

## RiskMap 2.0.0

RiskMap 2.0.0 is a redesign of the package. Many functions of 1.0.0 have
been renamed, merged or removed, so code written for 1.0.0 needs
changing. The table gives the replacement of each function of 1.0.0, and
the main changes follow.

### Functions of 1.0.0 and their replacements

| 1.0.0 | 2.0.0 |
|----|----|
| `assess_pp()` | [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md) |
| `assess_sim()` | [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md) |
| `check_mcmc()` | [`plot_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/plot_mcmc.md) |
| `compute_ID_coords()` | [`create_ids()`](https://claudiofronterre.github.io/RiskMap/reference/create_ids.md) |
| `convex_hull_sf()` | [`create_convex_hull()`](https://claudiofronterre.github.io/RiskMap/reference/create_convex_hull.md) |
| `dast()` | Removed, to move to a separate package |
| `dist_summaries()` | [`summarise_distance()`](https://claudiofronterre.github.io/RiskMap/reference/summarise_distance.md) |
| `glgpm_sim()` | [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md) |
| `Laplace_sampling_MCMC()` | Removed, now internal |
| `matern.grad.phi()` | [`matern_gradient_phi()`](https://claudiofronterre.github.io/RiskMap/reference/matern_gradient_phi.md) |
| `matern.hessian.phi()` | [`matern_hessian_phi()`](https://claudiofronterre.github.io/RiskMap/reference/matern_hessian_phi.md) |
| `matern_cor()` | [`matern_correlation()`](https://claudiofronterre.github.io/RiskMap/reference/matern_correlation.md) |
| `maxim.integrand()` | [`maxim_integrand()`](https://claudiofronterre.github.io/RiskMap/reference/maxim_integrand.md) |
| `plot_AnPIT()`, `plot_score()` | [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md) results |
| `plot_mda()` | Removed, to move with `dast()` |
| `plot_s_variogram()` | [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md) results |
| `plot_sim_surf()` | [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md) results |
| `pred_over_grid()`, `pred_target_grid()` | [`predict()`](https://rdrr.io/r/stats/predict.html) |
| `pred_target_shp()` | [`predict()`](https://rdrr.io/r/stats/predict.html) with `boundaries` |
| `propose_utm()` | [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md) |
| `s_variogram()` | [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md) |
| `set_control_sim()` | [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md) |
| `surf_sim()` | [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md) |
| `update_predictors()` | Removed, give the predictors in the `newdata` of [`predict()`](https://rdrr.io/r/stats/predict.html) |

### Main changes

#### Model fitting

- [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  fits Gaussian, binomial, Poisson and negative binomial models. A link
  can be named in a family object, such as `binomial("probit")`, and any
  other link given as an inverse link function.
- Non-Gaussian models are fitted by certified Monte Carlo maximum
  likelihood, which enlarges its samples until the Monte Carlo error of
  the estimates is small and reports a bound on that error.
- In binomial and Poisson fits, household effects whose children differ
  in covariates are now integrated out, so such fits are much faster.
- [`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md)
  reports how a fit was computed and flags anything that needs
  attention. [`profile()`](https://rdrr.io/r/stats/profile.html) gives
  profile likelihoods and
  [`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)
  bootstrap intervals for the parameters.

#### Prediction

- A single [`predict()`](https://rdrr.io/r/stats/predict.html) method
  replaces the three-stage interface of 1.0.0. It predicts at new
  locations or over the areas given in `boundaries`, with exceedance
  probabilities, population weights, standardisation and, if asked, the
  uncertainty of the estimated parameters.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) maps
  predictions,
  [`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)
  tabulates predictions over areas, and
  [`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md)
  shows where predictions go beyond the data.

#### Data and projections

- Data are `sf` objects with a coordinate reference system.
  [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
  proposes a projection that keeps distances accurate, and
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  uses it for data in longitude and latitude.
- [`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md)
  builds prediction grids, also from the cells of a raster,
  [`get_worldpop()`](https://claudiofronterre.github.io/RiskMap/reference/get_worldpop.md)
  downloads WorldPop population counts, and
  [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md)
  gives the empirical variogram with a permutation envelope.

#### Assessment and simulation

- [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
  compares models by spatial cross-validation.
- [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
  and
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
  simulate data and surfaces from a model, and
  [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md)
  assesses prediction targets on them.

#### Messages

- Errors, warnings and messages are written with cli and have classes,
  listed in
  [`?riskmap_conditions`](https://claudiofronterre.github.io/RiskMap/reference/riskmap_conditions.md),
  and long computations show their progress
  ([`?riskmap_progress`](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md)).
  `options(RiskMap.messages = FALSE)` silences the messages.
