# Simulation API migration (#115)

One public simulator now covers replicated data and spatial surfaces.
`simulate_surface()` and the old arguments of `simulate_glgpm()` are removed.
This is an intentional development API change. Existing stored simulation
objects must be regenerated before using the new accessors or assessment code.

## Chapter 3: parametric bootstrap

Use a fitted model as the data-generating mechanism. Its estimates stay fixed
while new latent effects and outcomes are generated. These draws are
unconditional on the observed responses and are not posterior predictions.

```r
liberia_boot <- simulate_glgpm(fit_liberia_no_nugget,
                              nsim = 100,
                              what = "data",
                              seed = 2026)

sim_data <- simulated_data(liberia_boot, simulation = 1)
fit_sim <- glgpm(npos ~ log(elevation) + gp(),
                data = sim_data,
                family = "binomial",
                den = ntest,
                distance_units = fit_liberia_no_nugget$distance_units,
                par0 = coef(fit_liberia_no_nugget))
```

`simulated_data()` replaces the response column with the simulated values.
There are no `npos_sim1` columns to select, and `liberia_boot[[1]]` is not a
dataset. Repeat extraction and refitting for each simulation.

Required book changes:

Tracked in [book issue #82](https://github.com/claudiofronterre/book_MBG/issues/82)
and [saved-object regeneration #80](https://github.com/claudiofronterre/book_MBG/issues/80).

- Update the introductory simulation call, output explanation, refitting loop
  and simulation exercise in Chapter 3.
- Update `R/CH3.R` to match the displayed workflow and regenerate bootstrap
  outputs with a recorded seed after the new API is merged.
- Keep expensive refits in the maintenance script, with the producing code
  visible and the saved results loaded silently in the chapter.

## Chapter 4: data and a known surface

```r
sim_study <- simulate_glgpm(true_model,
                           nsim = 200,
                           what = c("data", "surface"),
                           prediction_grid = pred_grid,
                           seed = 2026)

survey <- simulated_data(sim_study, simulation = 1)
surface <- simulated_surface(sim_study, simulation = 1)
values <- simulated_values(sim_study, component = "mean")
```

The survey defaults to the fitted model's original locations. For a different
fixed sampling design, provide `sample_locations`. Both location sets must
contain the required covariates and share the model's projected CRS. Survey
data also need the original denominator column. A surface alone does not need
denominators because its mean is a probability or rate, not a sampled count.

Both sets share one joint spatial realization. Repeated coordinates share the
spatial effect and nugget, and matching group labels share a group effect.
Response noise is drawn independently for each observation. The nugget enters
the linear predictor wherever it is included by the generating model.

Required book changes:

- Replace `simulate_surface()` and the `sampling_f_lib()` wrapper with this
  joint call. All covariates must be supplied at their actual locations.
- Replace direct accesses to the old simulation lists with the accessors.
- Regenerate stored simulation and assessment results in `R/CH4.R` after the
  implementation and the scientific targets have been reviewed.
- Keep the statistical redesign of assessment and classification summaries
  tracked separately under #108 and #109.

## Hypothetical generating models

`specify_glgpm()` keeps model parameters outside the simulation call.
It accepts an `sf` point dataset, a formula and parameters on their natural scale.
For example, using the package's test data:

```r
scenario <- specify_glgpm(y ~ cov + gp(),
                          data = gaussian_data,
                          family = "gaussian",
                          parameters = list(beta = c(1, 0.5),
                                            sigma2 = 1,
                                            phi = 2,
                                            sigma2_me = 0.1),
                          distance_units = "km")
sim <- simulate_glgpm(scenario, nsim = 10, seed = 1)
```

Use `denominator = "column_name"` for binomial totals or Poisson exposure.
Offsets remain in the formula. One generating model is accepted per call.
Support for a named collection of generating scenarios remains a separate
extension, as it was an open question in #115 rather than a prerequisite.

## Output and current limits

Geometry and covariates are stored once per location set. Simulation values
are stored as numeric arrays to avoid repeating geometry and strings for every
draw. The public `simulated_values()` accessor returns a tibble keyed by
`simulation`, `location_set`, `location_id` and `component`.

Components are `spatial_effect`, `nugget`, `group_effect`, `linear_predictor`,
`mean` and, for simulated data, `response`. The mean is the inverse-link target,
so expected counts additionally involve the denominator or exposure.

Exact joint simulation uses a dense covariance matrix over the unique
locations. It avoids approximation but remains expensive for very large grids.
Sampling locations are fixed within a call. A study with a different design
for each replicate can call the simulator separately for each design.

The current `assess_simulation()` interface accepts joint output for generating
models without grouped effects or custom inverse links. It explicitly rejects
these unsupported assessment cases even though the simulator supports them.
The candidate-model API and metric definitions remain the subject of #108.
