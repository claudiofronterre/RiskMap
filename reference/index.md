# Package index

## Data and projections

Choose a projected coordinate reference system, summarise the distances
between locations, and build prediction grids and covariates.

- [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
  : Propose a projected CRS for the data
- [`crs_distortion()`](https://claudiofronterre.github.io/RiskMap/reference/crs_distortion.md)
  : Distortion of distances by a projected CRS
- [`summarise_distance()`](https://claudiofronterre.github.io/RiskMap/reference/summarise_distance.md)
  [`print(`*`<RiskMap_distance_summary>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/summarise_distance.md)
  : Summaries of the distances
- [`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md)
  : Create Grid of Points Within Boundaries
- [`create_convex_hull()`](https://claudiofronterre.github.io/RiskMap/reference/create_convex_hull.md)
  : Convex Hull of an sf Object
- [`get_worldpop()`](https://claudiofronterre.github.io/RiskMap/reference/get_worldpop.md)
  : Download gridded population counts from WorldPop

## Exploratory analysis

Look at the data and their spatial correlation before fitting a model.

- [`elogit()`](https://claudiofronterre.github.io/RiskMap/reference/elogit.md)
  : Empirical Logit Transformation for Binomial Data
- [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md)
  : Empirical variogram with a global permutation envelope
- [`print(`*`<RiskMap_variogram>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_variogram.md)
  : Print Method for RiskMap_variogram Objects

## Specifying and fitting

Write the model as a formula with a Gaussian process and random effects,
choose the family and link, and fit it by maximum likelihood.

- [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  : Estimation of Generalized Linear Gaussian Process Models

### Model terms

The spatial Gaussian process and unstructured random effects, as terms
of the model formula.

- [`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md) :
  Gaussian Process Model Specification
- [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md) :
  Random Effect Model Specification

### Families and links

Gaussian, binomial and Poisson responses use the families of the stats
package, such as `binomial("probit")`, and other links can be given to
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
as an inverse link function.

- [`negative_binomial()`](https://claudiofronterre.github.io/RiskMap/reference/negative_binomial.md)
  : The Negative Binomial Family

### Controls

Settings of Monte Carlo maximum likelihood and of the MCMC sampler of
the random effects.

- [`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md)
  : Set control parameters for certified MCML estimation
- [`set_control_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcmc.md)
  : Set Control Parameters for Simulation

### Fitted models

Print, summarise and tabulate a fit, and extract its estimates and
fitted values.

- [`print(`*`<RiskMap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap.md)
  : Print a Fitted RiskMap Model
- [`summary(`*`<RiskMap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)
  : Summarize Model Fits
- [`print(`*`<summary.RiskMap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.summary.RiskMap.md)
  : Print Summary of RiskMap Model
- [`coef(`*`<RiskMap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md)
  : Extract Parameter Estimates from a "RiskMap" Model Fit
- [`fitted(`*`<RiskMap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/fitted.RiskMap.md)
  : Extract fitted values from a RiskMap model
- [`to_table()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.md)
  : Format RiskMap Results as a Table

## Checking a fit

See how a fit was computed, check the Monte Carlo samples, and measure
the uncertainty of the estimates by profile likelihood or the parametric
bootstrap.

- [`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md)
  : Check How a Fitted Model Was Computed
- [`print(`*`<RiskMap_check>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_check.md)
  : Print the Checks on a Fitted Model
- [`summary(`*`<RiskMap_check>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_check.md)
  : Summarise the Checks on a Fitted Model
- [`print(`*`<summary.RiskMap_check>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.summary.RiskMap_check.md)
  : Print the Summary of the Checks on a Fitted Model
- [`plot_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/plot_mcmc.md)
  : Check the mixing of the MCMC samples of the random effects
- [`summarise_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/summarise_mcmc.md)
  : Effective sample sizes of the MCMC samples of the random effects
- [`profile(`*`<RiskMap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md)
  : Profile likelihood of a fitted geostatistical model
- [`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)
  : Parametric bootstrap of a fitted geostatistical model
- [`print(`*`<RiskMap_bootstrap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_bootstrap.md)
  [`print(`*`<summary.RiskMap_bootstrap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_bootstrap.md)
  : Print a parametric bootstrap
- [`summary(`*`<RiskMap_bootstrap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_bootstrap.md)
  : Summarise a parametric bootstrap

## Prediction

Predict the spatial process, the linear predictor and targets on the
scale of the response, at new locations or aggregated over areas.

- [`predict(`*`<RiskMap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
  : Spatial prediction from a fitted RiskMap model
- [`print(`*`<RiskMap_prediction>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_prediction.md)
  : Print a RiskMap prediction
- [`to_table(`*`<RiskMap_prediction>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/to_table.RiskMap_prediction.md)
  : Table of polygon predictions
- [`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md)
  [`print(`*`<RiskMap_extrapolation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md)
  : Assess where predictions extrapolate beyond the data

## Assessment

Assess the predictive performance of models by spatial cross-validation
and on simulated data.

- [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
  : Assess Predictive Performance via Spatial Cross-Validation
- [`print(`*`<RiskMap_cross_validation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_cross_validation.md)
  : Print the Results of a Spatial Cross-Validation
- [`summary(`*`<RiskMap_cross_validation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md)
  : Summarize Cross-Validation Scores for Spatial RiskMap Models
- [`print(`*`<summary.RiskMap_cross_validation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.summary.RiskMap_cross_validation.md)
  : Print Summary of RiskMap Spatial Cross-Validation Scores
- [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md)
  : Assess Simulations
- [`print(`*`<RiskMap_assess_simulation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_assess_simulation.md)
  : Print an Assessment of Simulations
- [`summary(`*`<RiskMap_assess_simulation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_assess_simulation.md)
  : Summarize Simulation Results
- [`print(`*`<summary.RiskMap_assess_simulation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.summary.RiskMap_assess_simulation.md)
  : Print Simulation Results

## Simulation

Specify a model with known parameters and simulate data and spatial
surfaces from it.

- [`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
  : Specify a model for simulation
- [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md)
  : Simulate replicated data and spatial surfaces
- [`simulated_data()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md)
  [`simulated_surface()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_data.md)
  : Extract a simulated dataset or surface
- [`simulated_values()`](https://claudiofronterre.github.io/RiskMap/reference/simulated_values.md)
  : Extract simulated values in tidy form
- [`print(`*`<RiskMap_simulation_model>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_simulation_model.md)
  [`print(`*`<RiskMap_simulation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap_simulation_model.md)
  : Print a simulation model or simulations

## Plots

Plot methods for the objects that RiskMap returns.

- [`plot(`*`<RiskMap_variogram>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_variogram.md)
  : Plot Method for RiskMap_variogram Objects
- [`plot(`*`<RiskMap_check>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_check.md)
  : Plot the Checks on a Fitted Model
- [`plot(`*`<RiskMap_profile>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_profile.md)
  : Plot profile likelihoods of a fitted geostatistical model
- [`plot(`*`<RiskMap_bootstrap>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_bootstrap.md)
  : Plot a parametric bootstrap
- [`plot(`*`<RiskMap_prediction>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_prediction.md)
  : Plot a RiskMap prediction
- [`plot(`*`<RiskMap_extrapolation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_extrapolation.md)
  : Map indicators of extrapolation
- [`plot(`*`<RiskMap_cross_validation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_cross_validation.md)
  : Plot Method for RiskMap_cross_validation Objects
- [`plot(`*`<RiskMap_simulation>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_simulation.md)
  : Plot Method for RiskMap_simulation Objects
- [`plot(`*`<RiskMap_simulated_surface>`*`)`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_simulated_surface.md)
  : Plot Method for RiskMap_simulated_surface Objects

## Datasets

Data used in the examples and in the book.

- [`abund_sma`](https://claudiofronterre.github.io/RiskMap/reference/abund_sma.md)
  :

  Female *Culex pipiens* abundance (collections) in the Sacramento
  Metropolitan Area

- [`anopheles`](https://claudiofronterre.github.io/RiskMap/reference/anopheles.md)
  :

  *Anopheles* mosquitoes in Southern Cameroon

- [`galicia`](https://claudiofronterre.github.io/RiskMap/reference/galicia.md)
  : Heavy metal biomonitoring in Galicia

- [`infect_sma`](https://claudiofronterre.github.io/RiskMap/reference/infect_sma.md)
  :

  West Nile virus pool tests for female *Culex pipiens* in the
  Sacramento Metropolitan Area

- [`italy_sim`](https://claudiofronterre.github.io/RiskMap/reference/italy_sim.md)
  : Simulated dataset on the Italian peninsula

- [`liberia`](https://claudiofronterre.github.io/RiskMap/reference/liberia.md)
  : River-blindness in Liberia

- [`loaloa`](https://claudiofronterre.github.io/RiskMap/reference/loaloa.md)
  : Loa loa prevalence data from 197 village surveys

- [`malkenya`](https://claudiofronterre.github.io/RiskMap/reference/malkenya.md)
  : Malaria Transmission in the Western Kenyan Highlands

- [`malnutrition`](https://claudiofronterre.github.io/RiskMap/reference/malnutrition.md)
  : Malnutrition in Ghana

- [`tz_covariates`](https://claudiofronterre.github.io/RiskMap/reference/tz_covariates.md)
  : Covariates Dataset for Malaria Prediction in Tanzania

- [`tz_malaria`](https://claudiofronterre.github.io/RiskMap/reference/tz_malaria.md)
  : Malaria Dataset from Tanzania Demographic Health Surveys 2015

## Messages, conditions and progress

Control what RiskMap reports while it works, and catch its errors and
warnings by class.

- [`riskmap_conditions`](https://claudiofronterre.github.io/RiskMap/reference/riskmap_conditions.md)
  : Conditions raised by RiskMap
- [`riskmap_progress`](https://claudiofronterre.github.io/RiskMap/reference/riskmap_progress.md)
  : Progress of long computations

## Building blocks

Lower-level functions used by the fitting engine.

- [`matern_correlation()`](https://claudiofronterre.github.io/RiskMap/reference/matern_correlation.md)
  : Matern Correlation Function
- [`matern_gradient_phi()`](https://claudiofronterre.github.io/RiskMap/reference/matern_gradient_phi.md)
  : First Derivative with Respect to \\\phi\\
- [`matern_hessian_phi()`](https://claudiofronterre.github.io/RiskMap/reference/matern_hessian_phi.md)
  : Second Derivative with Respect to \\\phi\\
- [`maxim_integrand()`](https://claudiofronterre.github.io/RiskMap/reference/maxim_integrand.md)
  : Maximization of the Integrand for Generalized Linear Gaussian
  Process Models
- [`create_ids()`](https://claudiofronterre.github.io/RiskMap/reference/create_ids.md)
  : Compute Unique Coordinate Identifiers
