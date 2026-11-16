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
| `compute_ID_coords()` | Removed, now internal |
| `convex_hull_sf()` | `create_convex_hull()` |
| `dast()` | Removed, to move to a separate package |
| `dist_summaries()` | `summarise_distance()` |
| `glgpm_sim()` | `simulate_glgpm()` |
| `Laplace_sampling_MCMC()` | Removed, now internal |
| `matern.grad.phi()` | Removed, now internal |
| `matern.hessian.phi()` | Removed, now internal |
| `matern_cor()` | `matern_correlation()` |
| `maxim.integrand()` | Removed, now internal |
| `plot_AnPIT()`, `plot_score()` | `plot_cross_validation()`, or `plot()` of `assess_prediction()` results |
| `plot_mda()` | Removed, to move with `dast()` |
| `plot_s_variogram()` | `plot_variogram()`, or `plot()` of `variogram()` results |
| `plot_sim_surf()` | `plot_simulation()`, or `plot()` of `simulate_glgpm()` results |
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
- `start_pars` gives starting values for any of a model's parameters by
  name, one, several or all of them, and `start_pars = coef(fit)` starts from
  another fit. For binomial, Poisson and negative binomial models the automatic
  comparison of starts still runs for the parameters not given, which
  `set_control_mcml(start_comparison = FALSE)` switches off.
  `glgpm_parameters()` lists a model's parameters before it is fitted, and
  `?riskmap_parameters` describes them.
- `fixed_pars` in `glgpm()` holds any parameters at given values, such as a
  coefficient, the scale or the nugget variance, or all of them, in which case
  `predict()` works with the given parameters. It replaces
  `gp(nugget = <value>)`, `fix_var_me` and `negative_binomial(theta)`.
  `coef()` gives held parameters with the estimated ones, in the order of
  `?riskmap_parameters`, and `profile()` names the nugget `tau2` for every
  family.
- Gaussian fits with more than three records per latent effect start the
  measurement error and random-effect variances at the values given in
  `start_pars`. Such fits given a spatial variance other than 1, and parametric
  bootstraps of such models, can move very slightly.
- `glgpm()` takes its CRS in `distance_crs`, with `seed` and `messages` as
  its last arguments. `logLik()`, `vcov()` and `nobs()` work on fits, so
  `AIC()` and `BIC()` do too, and `?glgpm` lists the components of a fit that
  are part of its interface. `print()` and `summary()` of a fit stop on an
  argument they don't use, and every result records the version of RiskMap
  that made it.
- `set_control_mcmc()` takes the number of draws it keeps in `n_draws`,
  2,000 by default, instead of the number of iterations in `n_sim`, and its
  step-size settings are `c1_h` and `c2_h`. It no longer has `linear_model`,
  since the draws of a Gaussian model are always independent.
  `set_control_mcml()` takes its sample sizes as `initial_draws`,
  `validation_draws` and `max_draws`.
- `check_fit()` reports how a fit was computed and flags anything that needs
  attention. `profile()` gives profile likelihoods and `parametric_bootstrap()`
  bootstrap intervals for the parameters.
- `plot()` of a profile returns one plot instead of the list of its panels,
  and draws nothing, as every `plot()` method now does. `profile()`, `print()`
  and `plot()` of profiles and extrapolations, and `print()`, `summary()` and
  `plot()` of bootstraps and checks stop on an argument they don't use.
- Every `plot()` method has a function of its own, `plot_variogram()`,
  `plot_prediction()`, `plot_extrapolation()`, `plot_check()`,
  `plot_profile()`, `plot_bootstrap()`, `plot_cross_validation()`,
  `plot_simulation()` and `plot_simulated_surface()`, so users can call either.

### Prediction

- A single `predict()` method replaces the three-stage interface of 1.0.0. It
  predicts at new locations or over the areas given in `boundaries`, with
  exceedance probabilities, population weights, standardisation and, if asked,
  the uncertainty of the estimated parameters.
- `plot()` maps predictions, `to_table()` tabulates predictions over areas,
  and `assess_extrapolation()` shows where predictions go beyond the data.
- `assess_extrapolation()` counts the sampling noise of the data and flags a
  location when the data say less about it than about 95% of the surveyed
  locations, each left out in turn (#236). `plot()` of predictions takes `summaries` and
  `thresholds`, the names `predict()` uses.
- `predict()` takes `joint` in place of `type` and `standardise` in place of
  `standardize`, and stops on an argument it doesn't use. `plot()` of
  simulations maps the mean by default and takes `quantity`, `palette` and
  `limits` as `plot()` of predictions does.

### Data and projections

- Data are `sf` objects with a coordinate reference system. `propose_crs()`
  proposes a projection that keeps distances accurate, and `glgpm()` uses it
  for data in longitude and latitude.
- `create_grid()` builds prediction grids, also from the cells of a raster,
  `get_worldpop()` downloads WorldPop population counts, and `variogram()`
  gives the empirical variogram with a permutation envelope.
- `create_grid()` takes `distance_crs`, as `variogram()` and
  `summarise_distance()` do. A grid keeps its cells when its columns are
  selected, and `predict()` says when it counts points on a grid whole in the
  polygon that contains them (#214).
- ggspatial and jsonlite are suggested rather than imported. `get_worldpop()`
  asks for jsonlite, and maps ask for ggspatial only to draw a north arrow or
  a scale bar.

### Assessment and simulation

- `assess_prediction()` compares models by spatial cross-validation.
- `specify_glgpm()` and `simulate_glgpm()` simulate data and surfaces from a
  model, and `assess_simulation()` assesses prediction targets on them.
- `specify_glgpm()` takes its CRS in `distance_crs` and `plot_mcmc()` its
  locations in `locations`.
- `assess_prediction()` reports the SCRPS as a loss, so that lower is better
  for every score, and `PIT_area` once for all test observations, next to its
  value under perfect calibration (#230). It takes its fits in `models` and
  `min_dist` in their distance units, and draws the splits only with
  `plot_fold = TRUE`. `plot()` of its result returns one plot of one metric
  (#213), and `summary()` loses `view_all`. `assess_simulation()` takes
  `simulation`, `aggregate`, `label`, `objective` and `breaks`, and calls the
  scale of the grid `"location"`. `method = "cluster"` and `"regularized"`
  need spatialsample and spatialEco, which move to Suggests.

### Messages

- Errors, warnings and messages are written with cli and have classes, listed
  in `?riskmap_conditions`, and long computations show their progress
  (`?riskmap_progress`). `options(RiskMap.messages = FALSE)` silences the
  messages.
