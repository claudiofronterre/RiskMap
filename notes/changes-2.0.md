# RiskMap 2.0.0, full record of changes

The record of every change made while RiskMap 2.0.0 was developed, for the
people working on the package. `NEWS.md` has the summary for users. A change
that users will notice gets one line in `NEWS.md` and its details here.

- This version introduces many breaking changes from 1.0.0, and keeps no
function of 1.0.0 as a stub. The table of the functions of 1.0.0 and their
replacements is in `NEWS.md`. Every exported function has an example.

## Errors, warnings and messages

- Every error, warning and message of RiskMap is raised through cli. Errors
name the function that the user called, even when the problem is found deep
inside it, and an input error describes the value supplied. Messages have
been reworded to one style, with arguments, values and functions marked up.
Code that matched the old wording of a message should match its class
instead.
- Conditions have classes, documented in `?riskmap_conditions`, so that code
can catch them by class rather than by their text. Every error RiskMap raises
on purpose has class `riskmap_error`. An invalid argument has
`riskmap_error_input`, whose field `arg` names the argument, a fit that fails
`riskmap_error_fit`, and a state RiskMap should never reach
`riskmap_error_internal`, which asks for a report. Warnings that a fit may not
be reliable, such as an uncertified Monte Carlo error, have
`riskmap_warning_fit`.
- Errors from the compiled Langevin sampler are raised again as
`riskmap_error_fit`, with the original error as the cause.
- The `messages` argument of `glgpm()`, `predict()`, `profile()`,
`parametric_bootstrap()`, `check_fit()`, `assess_prediction()` and
`assess_simulation()` defaults to the option `RiskMap.messages`, as
`get_worldpop()`'s does, so `options(RiskMap.messages = FALSE)` silences them
all.
- Some arguments are checked more strictly. `glgpm()` rejects a starting value
of 0 for a variance of `re()` terms, and `predict()` rejects points with Z or
M coordinates or empty points in `newdata`. `match.arg()` errors are replaced
by RiskMap input errors, and a unique abbreviation is now accepted wherever a
choice is made, as `match.arg()` allows.
- Leaving out an argument that has no default gives a RiskMap input error
that names the argument and the function, such as "`object` is missing, with
no default." from `check_fit()`, instead of R's error where the argument is
first used.
- Messages start with a symbol, as the printed output does. Most have an
information sign, a step done has a tick, such as a certified fit or a
finished bootstrap, and something to note has an exclamation mark. `glgpm()`
says in one message which coordinate reference system the model is fitted in
and the units of distance. For data in longitude and latitude this is the
message on their projection, which also says how to choose another.
- A joint prediction in blocks warns once when the covariances of some
blocks needed a diagonal jitter to be factorised, with how many blocks and
the largest jitter, instead of once for each block.
- `get_worldpop()` takes the default of `messages` from the option
`RiskMap.messages`, so `options(RiskMap.messages = FALSE)` silences it.
- A `max_dist` or `phi` that isn't positive is now reported as "must be a
single positive finite number" instead of "strictly between 0 and Inf".

## Progress

- `parametric_bootstrap()`, `profile()`, `assess_prediction()` and
`assess_simulation()` show a progress bar with the work done and an estimate
of the time left, and `predict()` shows one while it predicts at many
locations in blocks or propagates parameter uncertainty. `glgpm()` shows a
line for each stage of Monte Carlo maximum likelihood, with the time the stage
took. See `?riskmap_progress`.
- While a step of a fit runs in the console, its line has a spinner and the
time so far, which move each time the likelihood is evaluated, so a long stage
shows it is still running. This applies to the starting values from the
Laplace approximation, each stage of Monte Carlo maximum likelihood, and the
re-maximisation that `profile()` does before profiling. A step computed in one
call to compiled code keeps a still line. These are the conditional samples
drawn at the estimate, and the Langevin sampler that `predict()` runs when
given `control_mcmc`, which now has a line of its own. A line too wide for the
console is shortened so that its time still shows.
- `profile()` shows its re-maximisation as a step in the console, and as the
message it was before where progress is hidden.
- Progress is shown when `messages` is `TRUE` and the option
`RiskMap.progress` allows it. The option's default is
`rlang::is_interactive()`, so progress appears in the console but not in
scripts, in documents rendered by knitr or in tests, even when these are
started from the console. `options(RiskMap.progress = TRUE)` shows it in the log of a script,
where a bar writes a line for each tenth of the work.
- The bars replace the progress messages of the bootstrap ("k of n replicates
done"), the profiles ("Profiling sigma2 (1 of 3)"), cross-validation
("Re-estimating model for subset i") and simulation studies ("Simulation j of
n").
- With `cores` above 1, the bootstrap and the profiles start one process per
task, at most `cores` at a time, and collect each as it finishes, instead of
dividing the tasks between the processes in advance. The results don't
change. A profile on several cores now also switches the main process to one
BLAS thread while the workers run, as the bootstrap already did.

## Printed output

- The printed output of a fit, of its `summary()`, of `check_fit()` and of the
summary of the checks is formatted with cli. Each starts with a rule that
carries its title, the call is shown faint, and each verdict has a symbol. A
certified fit and an OK check have a green tick, a check that needs care a
yellow exclamation mark, and a problem a red cross. Hints have an information
sign, and the fix of a check follows an arrow.
- In a console that supports links, such as RStudio, the functions named in
the printed output link to their help pages.
- The output is still printed with `cat()`, not raised as messages, so it
appears in knitr chunks with `message: false` and in `capture.output()`.
Without colour support, as in knitr, it is plain text, and
`options(cli.unicode = FALSE)` gives ASCII symbols, which a PDF without
Unicode fonts needs.
- Some printed wording changes. The certificate starts "Certified." or "Not
certified." instead of "Certified:" or "Not certified:", `summary()` gives the
link in one sentence, such as "The canonical logit link.", the log-likelihood
and the AIC have no colon, and a fix follows an arrow instead of "Fix:".
Argument and function names in the findings and fixes of `check_fit()` are
quoted with backticks, as in RiskMap's messages. Code that matched the old
printed text needs updating.
- The printed checks give each status as its symbol instead of the words OK,
SUSPICIOUS and PROBLEM, which stay in the column `status` of the element
`checks`. Code that searched the printed checks for them should read that
column.
- RiskMap now needs cli 3.6.0 or later, for the links to help pages.
- The `summary()` of a binomial, Poisson or negative binomial fit says where
its starting values came from. Either they were given in `start_pars`, with
RiskMap's defaults for any not given, or they were all RiskMap's defaults, or
they came from a Gaussian fit to the transformed data, which `glgpm()` uses
when the Laplace approximation of the likelihood prefers them. When the
defaults were kept without that comparison, as for a model with `re()` terms
or a link other than the canonical one, or after the comparison failed, the
line says why. Before, a failure was recorded only in
`attr(fit, "starting_values")$fallback_reason`. The choice is in the new
element `starting_values` of `mcml`, and the reason in `starting_reason`.
`glgpm()` no longer says "Using transformed-Gaussian automatic starting
values." when it makes that choice.
- In a console with colours, the labels of rows that line up labels and
values, such as the log-likelihood, the settings of a prediction or the
evidence of `check_fit()`, are grey, so that the values stand out.
- Printing the result of `assess_prediction()` gives the models, the test
sets, how the parameters were set and the scores averaged over the test sets,
instead of the whole list of scores and test sets. Its `summary()` gives the
averages and the scores of each test set as tables, instead of a line for
each score.
- The summary of `parametric_bootstrap()`, the profiles of `profile()`, the
variogram of `variogram()` and the distances of `summarise_distance()` are
printed in the same style. The variogram also says in how many distance
classes it lies outside the envelope.
- Predictions, models for simulation, simulations and assessments of
simulations are printed in the same style. A model for simulation names its
link as a fit does, and gives its parameters one per line. The summary of an
assessment of simulations gives each spatial scale a rule of its own, and a
failed model fit becomes a warning bullet.
- The result of `predict()` records the link of its fit, when it is not the
canonical one, in the new element `link`, so that its print names the link as
the print of the fit does.
- The `summary()` of a fit gives the parameters of its family, such as the
measurement error variance of a Gaussian model, in the element
`family_parameters`, a table of those that were estimated, and in
`held_family_parameters`, a list of those held at a given value, such as
`fix_var_me`. They replace the element `me`, which held either. So
`to_table()` of a fit with `fix_var_me` no longer adds a row that repeats the
fixed value in every column, and the row of the intercept of a fit with no
covariates is named "(Intercept)", where it was blank.

## Plots

- Each `plot()` method of RiskMap has an exported function of its own, which
holds its arguments, checks and documentation. They are `plot_variogram()`,
`plot_prediction()`, `plot_extrapolation()`, `plot_check()`, `plot_profile()`,
`plot_bootstrap()`, `plot_cross_validation()`, `plot_simulation()` and
`plot_simulated_surface()`. The method is a one-line call of the function and
is documented on its page, so `plot(x, ...)` and the function give the same
object. A function given the result of another function stops with an input
error that names the function whose result it draws, as in "`plot_variogram()`
needs the result of `variogram()`". Their unused arguments stop as before, and
the maps of predictions and simulations still pass theirs to the fill scale.
The help pages `plot.RiskMap_*` become aliases of the new pages, and the
reference index lists the functions in their task groups instead of a group of
plots. `plot_mcmc()` is unchanged.
- `plot()` of a variogram, which already stopped on `plot_envelope` and
`color`, now says that they are `envelope` and `colour`, and its error is a
RiskMap input error like those of the other methods.
- `plot()` of predictions, and of the indicators of `assess_extrapolation()`,
at points that did not come from `create_grid()` drew an empty map, because
ggplot2 dropped every point for its border of colour `NA`. The border is now
transparent.
- The palette centred at zero, for signed effects in `plot_prediction()`,
`plot_simulation()` and `plot_simulated_surface()`, goes from teal to orange
instead of teal to magenta. The two sides of zero were nearly the same colour in
red-green colour blindness, and they now have the same lightness.

## Datasets
- Apart from `liberia`, all datasets are now in an simple features (sf) format.
- `anopheles` is now in WGS 84 / UTM zone 32N (EPSG:32632), where it was in
Pseudo-Mercator (EPSG:3857), and `italy_sim` in WGS 84 / UTM zone 33N
(EPSG:32633), where it was in UTM zone 34N. The old CRSs distorted distances by
up to 1.0% and 1.4% over the data, and the new ones by less than 0.1% and 0.5%.
Distances, and so fits, variograms and grids made from these datasets, change
slightly. Pseudo-Mercator is also not conformal, so it stretched distances
more north-south than east-west.
- `?malnutrition` gives `age` in years, where it said months. The ages run
from 0 to 4.96 in steps of one day.
- `tz_covariates` has its covariates on the scales of `tz_malaria`, so that a
model fitted to `tz_malaria` predicts over the grid. Temperature (4.7 degrees
lower than at the nearest cluster at the median) and NTL (mostly 0) came from
other products than those of `tz_malaria`, and EVI and Precipitation from other
processing of the same products, with differences of up to 0.3 and 761 mm.
All four are now the values at the grid points of the 2015 rasters of MODIS
EVI and daytime land surface temperature, VIIRS night-time light and CHIRPS
precipitation that give `tz_malaria` its values. Temperature is missing at 2
points and Precipitation at 1, where the rasters have no value. ITN was a
proportion with -9999 for 8 missing values and is now a percentage with `NA`.
Population was already on the scale of `tz_malaria` and is unchanged.
`?tz_malaria` and `?tz_covariates` give the source and units of every
covariate. Population was said to be in thousands, but it is about one
hundredth of the density in people per square kilometre, as people per
WorldPop cell of 3 arc-seconds. The script `data-raw/tz_covariates.R` records
the change.
- `?liberia` spells river blindness correctly.

## Coordinate reference systems

- Datasets must have a coordinate reference system. When coordinates are in
latitude/longitude and `model_crs` is not supplied, they are now automatically
reprojected to the CRS from `propose_crs()`. `predict()`, `simulate_glgpm()` and
`assess_prediction()` have been updated to match; the `convert_to_crs` and
`scale_to_km` arguments of `glgpm_sim()` have been replaced with the
`distance_crs` and `distance_units` of `specify_glgpm()`.
- providing `model_crs` in longitude/latitude now raises an error.
- The new `propose_crs()` chooses the projected CRS that distorts distances
least over the data. It compares the WGS 84 UTM zone at the centre of the data
and its two neighbours with four conformal projections centred on the data
(oblique stereographic, transverse Mercator, oblique Mercator and Lambert
conformal conic), whose parameters minimise the largest scale error over the
region. It returns the EPSG code of the best UTM zone when that zone distorts
distances by at most 0.5% (`epsg_tolerance`), and otherwise a `crs` object for
the best centred projection, named "RiskMap ..." in its WKT. The new
`crs_distortion()` measures the distortion of any projected CRS over the data.
- `glgpm()`, `summarise_distance()`, `variogram()`, `create_grid()` and
`specify_glgpm()` now project data in longitude and latitude with
`propose_crs()`. Earlier development versions used the UTM zone holding most
locations, which could distort distances by more than 3%, for example over the
Democratic Republic of the Congo, where the CRS now chosen distorts them by
less than 0.5%. Data that a UTM zone covers well get the zone at their centre,
which may differ from the zone used before. Fits, distances and grids of data
in longitude and latitude can therefore change. Set `model_crs` or
`distance_crs`, or project the data first, to keep a given CRS. The message
names the CRS and the units of distance, and says how to choose another CRS.
It doesn't give the distortion, which the warning below reports when it is
above 1%, and which `crs_distortion()` measures for any CRS.
- These functions warn, with class `riskmap_warning_crs`, when the CRS used,
whether given, carried by the data or chosen, distorts distances by more than
1% over the data. When no planar CRS can do better, as over sub-Saharan
Africa, the warning says so.
- `propose_utm()` has been removed. Use `propose_crs()`, which returns the
same UTM zone for data that a zone covers well. `propose_utm()` took the zone
holding most locations, which could be far from the centre of the data.
- In `glgpm()`, `convert_to_crs` has been replaced with `model_crs` and
`scale_to_km` has been replaced with `distance_units` (`"km"` or `"m"`).
- `create_grid()` gains `distance_crs`, which projects `boundaries` as it
projects the data in `variogram()` and `summarise_distance()`, so the CRS of
distances has one name in all three. With `raster` it must be `NULL` or the
CRS of the raster. It comes before `distance_units`, so a call that gave
`distance_units` by position must now name it. The message for boundaries in
longitude and latitude now says to set `distance_crs`.
- `propose_crs()`'s `epsg_tolerance` is now `utm_tolerance`, since it is the
largest distortion for which a UTM zone is proposed.

## Exploratory analysis

- `plot()` for `variogram()` results gains `show_pairs` (default `FALSE`),
which labels each point with the number of pairs of locations in its distance
class, so that unstable semivariances from few pairs can be recognised.
- `plot()` for `variogram()` results draws the semivariances as points without
joining lines, omits distance classes without pairs (which previously caused
"Removed rows" warnings), starts both axes at zero, uses `theme_minimal()` as
other RiskMap plots do, and draws in red the points that fall outside the
envelope. `plot_envelope` is now checked.
- `variogram()` results gain a `print()` method that reports the distance
classes and the envelope settings before the table.
- `summarise_distance()` results gain a `print()` method that shows the
distances to three significant digits (set by `digits`), with the number of
unique locations and the distance units. The values are still stored in full.
- In `summarise_distance()` and `variogram()`, `scale_to_km` has been
replaced with `distance_units` (`"km"` or `"m"`), for consistency with `glgpm()`.
Both now default to `"km"`, as `glgpm()` does. For `variogram()` this changes
the default from `"m"`, so code that gives `breaks` or `max_dist` in metres
should set `distance_units = "m"`.
- `variogram()` accepts `breaks` or `max_dist` beyond the largest distance
between locations, and warns only when whole distance classes lie beyond it.
When no pair falls in any class, the error reports the range of the distances
instead of always suggesting a change of units. `max_dist = NULL` can be passed
together with `breaks`.
- `variogram()` no longer fails when the envelope is built for a single
distance class, and runs about twice as fast with 2,000 or more locations,
with identical results. Rounding error no longer makes the envelope exclude one
curve too few at some levels, such as 0.9 with 999 permutations. The default
level of 0.95 with 999 permutations was not affected.
- In `summarise_distance()` and `variogram()`, `convert_to_utm` has been replaced
with `distance_crs`, matching `glgpm()`'s `model_crs` convention: an
already-projected CRS is retained as-is, longitude/latitude data is automatically
reprojected to an appropriate UTM zone with a message, and an explicitly supplied
`distance_crs` must be projected. Distances now correctly account for the CRS's
native linear units rather than assuming metres.
- `variogram()` now uses a simultaneous extreme-rank-length permutation envelope
  instead of separate pointwise intervals. It remains an exploratory diagnostic
  and does not report a hypothesis-test p-value. Its pairwise permutation engine
  is compiled, distances and lag bins are computed once, and the default is 999
  permutations plus the observed curve. The returned columns are now `distance`,
  `semivariance`, `n_pairs`, `lower_envelope` and `upper_envelope`. A `seed` can
  be supplied without changing the caller's random-number state.
  `plot_variogram()`, which `plot()` of a variogram also calls, shows the
  envelope by default and warns when it was omitted with `n_permutations = 0`.
- The `bins` parameter in `variogram()` has been removed and replaced with `breaks`.
- `plot()` of a variogram takes `envelope` and `colour` instead of
`plot_envelope` and `color`, in the British spelling of RiskMap's other
names. It now stops when `...` is given anything, and so does `print()` of
`summarise_distance()` results, so that an old name or a misspelt one is not
ignored.
- `variogram()` takes `seed` last, after `distance_crs` and `distance_units`.
- `create_convex_hull()` takes `data` instead of `sf_object`.

## Model fitting

- `glgpm()` fits negative binomial models of overdispersed counts, with
`family = negative_binomial()`, the log link and the exposure in
`denominator`. The size `theta` is estimated, starting from the negative
binomial regression without random effects or from `start_pars$theta`, unless
it is given, as in `negative_binomial(theta = 2)` or
`MASS::negative.binomial(2)`, at which it is held. A size that goes to the
Poisson limit is held there and reported by `check_fit()`. An estimated size
can't be fitted with a nugget or a random effect with one level per
observation when there is one observation per location. `predict()`,
`simulate_glgpm()`, `specify_glgpm()`, `assess_prediction()`, `profile()`,
`parametric_bootstrap()` and `check_fit()` take the size of a fit, and
`maxim_integrand()` takes the family with its size. Its log-likelihood leaves
out the same terms as that of a Poisson fit, so the two can be compared.
- Several `re()` calls in a formula are now combined, so `re(a) + re(b)` is the
same model as `re(a, b)`. Previously only the last call was fitted, with no
warning. A second, different `gp()` term, which also used to replace the first
with no warning, is now an error.
- A second `offset()` term, which used to be ignored with no warning, is now an
error. `offset()` must name a column of `data`, and an expression such as
`offset(log(n))` now stops with a message that says so. An offset is also read
correctly when the formula has several `re()` calls.
- Breaking change: the seed has moved out of `set_control_mcmc()`, whose
`seed` argument is removed. `glgpm()`, `predict()`, `assess_prediction()` and
`assess_simulation()` take a `seed` argument instead, as `simulate_glgpm()`,
`variogram()` and `profile()` already did. Replacing
`control_mcmc = set_control_mcmc(seed = s)` with `seed = s` gives identical
results. Without a seed, `glgpm()` draws one from the session's random number
generator, so `set.seed()` before the call reproduces the fit, and two calls
in a row differ. The seed used is returned as `seed` in the fitted object, and
`glgpm(..., seed = fit$seed)` reproduces any binomial or Poisson fit,
including its final conditional samples.
- The `den` argument in `glgpm()` has been renamed to `denominator`.
- `glpgm()` now only accepts data in an sf format. 
Consequently the locations do not need to be passed to `gp()` when fitting a model
as they are included automatically.
- The `nugget` parameter in `gp()` is now `FALSE` by default and should be set to
`TRUE` to estimate it.

- Canonical binomial and Poisson models now use compiled likelihood and MALA
kernels. The implementation retains the existing sampler, including support
for nugget effects and multiple unstructured random effects, while avoiding
the main interpreted-R bottlenecks. Custom inverse-link functions continue to
use the general R implementation.

- `plot_mcmc()` draws its diagnostics with ggplot2 and patchwork and returns
them as a `patchwork` object, instead of drawing with base graphics and
resetting `par()`. By default it shows the trace plot and the autocorrelation
function of the mean of the spatial random effects and of the location with
the smallest effective sample size, with a histogram of the effective sample
sizes at all the locations. The trace plots show the iterations at which the
samples were kept, so fits now record the `n_sim`, `burnin` and `thin` of
their samples in `attr(fit, "mcmc")`. The chains shown are chosen with
`locations`, which was `component` in `check_mcmc()`. It defaults to the
location with the smallest effective sample size, can select several
locations and can be given together with `check_mean = TRUE`. The mean and
the chosen chains previously included the levels of unstructured random
effects, which are stored after the locations in `S_samples`. They now refer
to the locations only. An argument in `...` that doesn't match one of
`acf()`'s `lag.max`, `type`, `na.action` and `demean`, partially as `acf()`
matches them, is an error, since `acf()` would ignore it. `component` is one
such argument.

- The new `summarise_mcmc()` returns a tibble of the effective sample sizes of
the mean of the spatial random effects, of each location and of each level of
an unstructured random effect. Its column `index` numbers the locations and
the levels. It was called `component`, which now names only the quantities
of `simulated_values()`. The sizes are computed with the initial monotone
sequence estimator of Geyer (1992), so the package no longer depends on sns.
patchwork has moved from Suggests to Imports.

- Non-Gaussian optimisation now reuses the sample-wise Monte Carlo log
integrand when the objective, gradient, and Hessian are evaluated at identical
parameters, avoiding a repeated likelihood calculation.

- `profile()` computes profile likelihoods of binomial and Poisson models
fitted by certified Monte Carlo maximum likelihood, for single parameters and
for pairs, on one fixed bank of draws so that the profiles are smooth. It
reports profile and Wald intervals and says if a point on a profile is higher
than the estimate. Marginal grids of 11 points are extended with doubling steps
where the profile is flat, up to 11 points a side, and an end the profile does
not reach is reported as open. Interval ends are interpolated in the signed
root of the deviance, with one more point at a crossing between distant
points, and the points around each end are also solved from the estimate,
because constrained maxima can be multimodal. The constrained maximisations
work on the scale of the fit, with standardised covariates, which avoids the
false convergence that coefficients of covariates with large ranges caused in
their own units, and each starts from the extrapolation of the two previous
solutions. One that still stops with a non-zero code is solved again with a
central-difference gradient. A fixed coefficient whose constraint the start
cannot meet within the bounds is first moved to the nearest feasible point,
and each point records whether it was solved, is infeasible or failed. Joint
profiles are computed only when `pairs` asks for them, and the time left is estimated as the profiles
run. With `cores` above one, the two sides of each profile, the ends of the
intervals and the joint profiles run in parallel by forking (not on Windows),
with the same results as on one core with one BLAS thread. `plot()` draws the marginal profiles against their quadratic
approximation and the joint profiles as confidence regions with the marginal
ridges, as a corner plot for one pair.

- `profile()` also computes profile likelihoods of linear geostatistical
(Gaussian) models, from their exact likelihood, so `n_draws` and `seed` are
not used. Every estimated parameter can be profiled: the regression
coefficients, the spatial variance and scale, the variance of the nugget
(`tau2`, as `summary()` reports it, rather than its ratio to the spatial
variance), the measurement error variance unless it is fixed, and the
variances of the random effects, with the labels of `summary()`. For each
value of the covariance parameters the coefficients are found by generalised
least squares, and models with more than three observations per latent effect
use the Woodbury identity, as `glgpm()` does. The ends of the profile
intervals are refined until the deviance of a solved point is within 0.01 of
the cut-off, grid steps of the variances and the scale are at most 1 on the
log scale, so that a variance estimated near zero, with a huge standard
error, is still profiled, and each point of a joint profile is also solved
from the estimate, since a neighbour's solution with a variance near zero can
leave the optimiser on a plateau. A variance that the data cannot tell from
zero has a profile that stays below the cut-off down to its lower bound, and
its interval is reported as open there. The smoothness `kappa` stays fixed
and can be profiled by refitting over fixed values.

- In `plot()` of a profile, the subtitle of each marginal panel gives the
estimate, the profile interval and the Wald interval on separate lines, so
that it is no longer cut off in panels half or a third of a page wide, and
the two entries of the ridge legend of a joint profile are stacked.

- `plot()` of a profile returns one ggplot for one panel and one patchwork
object otherwise, and draws nothing, so the figure can be saved or changed
before it is printed. It drew the panels with `gridExtra::grid.arrange()`, or
printed the corner plot, and returned the list of panels invisibly. The panels
and the corner plot are unchanged, and patchwork arranges any other set of
panels in up to three columns. Asking for no panel stops. One panel is
`plot(x, parameters = "phi", pairs = character())`, since `pairs` defaults to
every joint profile.

- `profile()` and `plot_profile()` take a pair as two names, as in
`list(c("sigma2", "phi"))`, or as one, as in `"sigma2:phi"`, the name of its
joint profile. `profile()` took only the first and `plot_profile()` only the
second. `plot_profile()` also takes the names of coefficients without the
prefix `beta.`, as `profile()` does. A name with more than one colon, from a
coefficient of an interaction, is split where both sides are parameters of the
fit. A pair of one parameter twice is an input error, and two names in a
character vector stop with a hint to give them in a list or as
`"first:second"`. A pair or parameter named twice, in any form, is profiled
and drawn once. `plot_profile()` names a pair in the order given to
`profile()`.

- `profile()`, the print of a profile and its `plot()`, `summary()`, `print()`
and `plot()` of a parametric bootstrap and the print of its summary, and
`print()`, `summary()` and `plot()` of the checks of `check_fit()` and the
print of their summary stop on an argument they don't use instead of ignoring
it (`check_dots_unused()`).

- Non-Gaussian models are now fitted by certified Monte Carlo maximum
likelihood. The estimated likelihood is maximised on a fitting sample of 1,000
draws, the estimate is corrected by one Newton step with the score from an
independent validation sample of 5,000 draws, and the corrected estimate is
certified with a second, fresh validation sample. The certificate is an
approximate 95% bound, in statistical standard errors, on the Monte Carlo
error of the estimate at the maximum found. When it exceeds the tolerance
(0.1 by default), both samples are enlarged to the size predicted to meet it,
up to `max_samples`. The importance distribution is a Gaussian approximation fitted
to the whole conditional distribution of the latent effects (efficient
importance sampling), which keeps far more draws effective than the Laplace
approximation in sparse data. Binomial and Poisson models, with or without a
nugget and unstructured random effects, share one implementation, and the
Laplace approximation remains as a fallback and as an option. If the global
proposal fails while an estimate is being corrected or certified, the stage
is repeated with the Laplace proposal on fresh draws, so that every
certificate rests on validation draws taken after everything that led to it.
`set_control_mcml()` gains `validation_samples`, `tolerance` and `proposal`,
and loses `min_effective_samples`, `max_pareto_k`, `max_mcse_ratio` and
`score_probability`. The defaults of `max_samples` and `max_stages` are now
16,000 and 4. Fitted objects report `mcml_certified` in place of
`mcml_converged`.

- Certified MCML continues the Newton correction when one step leaves a
residual. Where the log-likelihood is far from quadratic between the fitted
estimate and the maximum, as along a spatial variance close to zero, one
Newton step stops short of the maximum, and the certificate, measured in the
curvature at the fitted estimate, sees only part of the shortfall. The fitting
sample now predicts the residual of the first step, at the cost of one
gradient on that sample, and when it exceeds a quarter of the validation
sample's radius, further Newton steps are taken with the score from the same
validation sample, each in the fitting sample's curvature at the current
point, until a step is below a tenth of the radius (at most four more steps).
The estimate is then certified in the curvature at the last step. If the steps
stop short by more than the certificate tolerates, they are discarded, the
stage is not certified, and the next stage starts from the estimate of the
first step. The check does not act where the validation sample's radius is
below 1e-6 statistical standard errors, as when the spatial variance is held
at zero and the likelihood is exact, because the predicted residual is then
the rounding error of the differences divided by the radius. Fits where the
check does not fire
are unchanged. Where it fires, the estimate and the certificate change, and so
do the returned standard errors, which now come from the curvature at the last
step. On a simulated survey with a log spatial variance of -4.56 at the
maximum, the standard error of the log variance moved from 14.5 to 29.5,
against 30.2 at the maximum. Each stage of `mcml_history` records the
predicted residual and the number of Newton steps.

- Certified MCML distinguishes Monte Carlo error from weak
identification. A variance component estimated as zero (and the spatial scale
when the spatial variance is zero) is held at the boundary, the first stage is
restarted once from an interior value, and precision is assessed on the
remaining parameters, which keep their standard errors. A spatial scale
estimated below the median distance between neighbouring locations is reported
as estimated, with a warning that gives its practical range, explains that
such a scale may be weakly identified, points to its profile likelihood and
suggests a nugget model, a fixed scale or the PC prior. If the information
matrix is still singular, fitting stops instead of drawing more samples. Each
case produces a plain-language warning.

- Certified MCML recognises data whose likelihood has no maximum, which
previously could be certified with a regression coefficient at the bound of the
search. When no observation has a positive, when every individual tested is
positive, or more generally when the covariates separate the data (for
example, a category of a factor whose locations have no positives), the
likelihood increases without limit as some coefficients tend to infinity. An
exact check of the data finds such a separation before the fit, and a
coefficient that reaches the bound of the search is treated in the same way.
The fit then ends with the outcome "no maximum", without standard errors or
certificate, and a warning names the cause and the covariates involved. On a
simulated survey with no positives among 5,000 individuals tested, six of ten
runs had been certified with an intercept of -30 and a standard error of about
23,500. Surveys with at least one positive in each such category keep a finite
maximum and are fitted as before.

- A binomial fit whose covariate separates the data completely, such as one
where every location with a positive covariate has only positives and every
other location none, now ends with the outcome "no maximum" as intended. It
stopped with "the estimated log-likelihood is not finite at the optimum",
although the separation had been found. With the logit link, the linear
predictor of a location where every individual is positive passes 36.74 well
inside the bounds of the search, the inverse logit rounds to 1 there, and the
formulas of the conditional mode of the latent effects divided zero by zero.
The canonical links now use the exact form of each observation's
log-likelihood, score and curvature where the probability rounds to 0 or 1 or
the mean of a count underflows to 0, in the fit and in `maxim_integrand()`.
The probit and complementary log-log links already ended such fits with "no
maximum". Fits that never reach these limits are unchanged, since the new forms
apply only where the previous ones could not be evaluated.

- `maxim_integrand()` now takes `denominator` instead of `units_m`, matching
`glgpm()`. Its documentation no longer lists the Gaussian family, which it has
never supported, and describes the Newton-type optimisation it actually uses.

- Removed unused internal helpers (`log_mean_exp()`,
`importance_effective_sample_size()`, `restore_spatial_range()`) and the
compiled kernel `cpp_canonical_likelihood()`, which no longer had callers.

- `set_control_mcml()` gains `pc_prior` (default `FALSE`) and its thresholds
`pc_sigma_threshold`, `pc_sigma_probability`, `pc_scale_threshold` and
`pc_scale_probability`. When `pc_prior = TRUE`, the spatial variance and scale
are penalised by penalised-complexity priors (by default
P(sigma > 1) = 0.05 and P(phi < median nearest-neighbour distance) = 0.05),
which stabilises estimation in small or low-prevalence surveys. Warnings about
a spatial variance estimated as zero, a small spatial scale or a singular
information matrix now suggest this option.

- The final conditional sampling in `glgpm()` (used for fitted values and
stored prediction samples) now uses the compiled MALA kernel for canonical
links, as intended. It previously passed the inverse-link function itself,
which selected the general R implementation with numerical derivatives; a
binomial fit to the `liberia` data now takes about 2 s instead of 27 s.
Fresh conditional runs in `predict()` now honour a custom inverse link.

- `set_control_mcml()` gains `laplace_start` (default `TRUE`), which maximises
the Laplace approximation from the starting values before MCML.

- `set_control_mcml()` gains `check_maxima` (default `TRUE`). The
certificate establishes Monte Carlo precision at the maximum found, and in
sparse surveys an optimiser can stop at a maximum with a spatial scale far
below the spacing between locations while a higher one exists at a longer
scale. When the fitted scale is below half the median distance between
neighbouring locations, the fit is repeated from four dispersed starting
values on the same draws. Candidates within 0.1 standard errors of each other
are treated as one maximum, and separate maxima are compared on independent
draws. A message says when and why the check runs, and a warning says when it
moves the estimate to another maximum or finds two maxima that fit about
equally well. The check does not prove that no higher maximum exists.

- The check for other maxima also looks at short scales after a fit at a long
one. The likelihood with independent location effects, the limit of the model
as the spatial scale goes to zero, is computed exactly by adaptive
Gauss-Hermite quadrature and maximised once per survey. The estimated
likelihood is evaluated at four scales from an eighth of the spacing between
neighbouring locations to the spacing, with the coefficients and variance of
that limit and one three-point step in the log variance where the value is
close, and the fit is repeated from the best point if it comes within 0.25
log-likelihood units of the fit or if an evaluation fails. The limit and Jensen
lower bounds at those scales are exact, and if one of them exceeds the fitted
log-likelihood, estimated on independent draws, a further search is forced and
a warning is given if it does not resolve the alarm. A fit whose scale reaches
the point where the correlation between the closest locations is below 1e-6
ends with the outcome "numerically independent", without standard errors or
certificate, since the fit is then one of independent location effects and the
spatial scale is unresolved. A fit at a short positive scale whose likelihood is not clearly
above the limit at the current Monte Carlo precision is assessed as usual, with
a warning that points to `profile()` for the scale. Separate maxima whose
comparison on independent draws is ambiguous are compared again on four times
as many draws, and that comparison serves every later one. A stage whose
global proposal fails at its start goes straight to the Laplace proposal,
instead of being completed and then discarded, and the factorisation of the
spatial covariance is reused while only the regression coefficients change.
The screen
and the exact bounds need the canonical links and no unstructured random
effects, and other models keep the earlier check.

- Certified MCML computes the importance weights (for the canonical binomial
and Poisson links) and the dense linear algebra of the proposal in compiled
code, with the same LAPACK and BLAS routines R uses, and gives results
identical to the previous R code (up to the last bits on platforms whose
compilers fuse multiply and add). A fit no longer creates several large
temporary matrices at every evaluation of the likelihood, whose garbage
collection took about half of a fit.

- Fitted models have a `print()` method, which gives the call, the estimates
of the regression coefficients and covariance parameters on their natural
scale, the log-likelihood and, for binomial and Poisson models, whether the
Monte Carlo error of the estimate is certified. Printing a fit previously
listed every component of the object, including the data. `summary()` of a
binomial or Poisson fit now reports the certificate (its bound in statistical
standard errors, the level of the bound, the tolerance, the number of stages
and the draws at the last stage), with short notes on what the warnings of the
fit reported and the advice they gave, so that both are available for a saved
model. The certificate is in the `mcml` component of the summary.

- The messages and warnings of certified MCML quote the spatial scale of the
returned estimate, after the Newton correction. The message that announces
the refits from dispersed starting values, which are decided before the
correction, now says that it quotes the scale of the first fit, before its
correction. When the repeated Newton correction does not converge at the last
stage, the warning says so, instead of describing the bound as above the
tolerance.

- A failure of the global importance distribution no longer condemns the whole
fit. In sparse data the Laplace estimate can have a very large nugget, where
the fixed-point iteration stalls and whether it converges can depend on
rounding, while it converges quickly elsewhere. A failure at the Laplace
estimate now restarts the first stage from the starting values with the global
importance distribution, and a stage that does fall back to the Laplace
approximation is followed by one that tries the global distribution again. Fits
without a failure are unchanged. `mcml_history` records why the iteration
failed.

- New `check_fit()` goes through the checks on how a model was computed, says
which are OK and which need a look, and gives the fix where there is one. Its
`plot()` method shows the certificate, the effective draws, the Pareto shape and
the changes in the estimates by stage. It reports what used to crowd the output,
including the effective draws and Pareto shape of the certification draws, a
heavy tail of the importance weights, changes of importance distribution and a
search for the conditional mode that stopped at its iteration limit.

- Certified MCML is quieter. A fit warns only about what changes how a result is
read, which is a certificate that is missing or above the tolerance, a likelihood
without a maximum, weak identification, other maxima and a spatial variance at
zero, in a sentence or two that point to `check_fit()`. A nugget or random effect
at zero, a short spatial scale and a scale close to the limit of independent
effects no longer warn. With `messages = TRUE` a certified fit gives one line.
`summary()` ends with the verdict and the number of checks that need a look, in
place of its notes, and `summary()$mcml$notes` is replaced by `checks_to_review`.

- Certified MCML now tests a nugget or random-effect variance that it holds at
zero. A fit moves in the log of a variance, whose gradient vanishes near zero
whatever the slope in the variance itself, so reaching the boundary does not
show that zero is a maximum. The slope at zero is estimated on an independent
bank of draws. If it is clearly positive the stage is not certified, because the
held point is not a maximum, and the next stage starts the variance from an
interior value, even at the largest number of draws if stages remain. If it is
clearly negative, `check_fit()` reports that zero is the
maximum and the variance no longer counts as a check that needs a look. Fits that
hold no variance are unchanged, and those with a custom link are not tested.

- The radius of the certificate, its allowance for the Monte Carlo error of
the certification sample's own score, is now the larger of that sample's radius
and the radius of a replicate of it, with the same number of draws at the same
estimate on a stream of its own. The radius of one sample is estimated from its
own draws, and it is too small when the sample misses the rare draws of very
large weight that dominate its error. On 13 simulated binomial surveys with 500
samples of 5,000 draws each, how often the error that the radius allows for
exceeded it fell from at most 6.7% to at most 2.9%, against the 5% the radius
is meant to hold, and how often the Monte Carlo error of the estimate exceeded
the certificate fell from at most 0.62% to at most 0.20%. On 16 simulated
Poisson surveys the first fell from at most 10.0% to at most 4.6%. The
certificate changes for every fit, and the replicate costs about 9% more fit
time. In 1,270 paired fits of 127 simulated binomial surveys every fit stayed
certified, the median certificate grew from 0.0171 to 0.0175 statistical
standard errors, and 12 more fits needed a second stage. `check_fit()` says
which sample set the radius in the `radius_from` column of its `stages` table,
and each stage of `mcml_history` records both radii (`certificate_radii`) and
the replicate's weight diagnostics (`replicate_validation`).

- With a nugget and the spatial variance at zero, the exact alarm of the check for
other maxima now also tries scales up to 16 times the spacing, with 2% to 30% of
the total variance in the spatial effect, and fires when a lower bound on the
likelihood is above the maximum of the model with independent location effects by
more than half the squared tolerance. The screen stopped at the spacing and only
tried a share near zero, which could miss a spatial effect that a large nugget
was absorbing.

- Certified MCML now fits custom inverse links that reach the limits of their
range or are not log-concave, which previously stopped with "no valid fit". A
probability that rounds to exactly one, as cloglog does from a linear predictor
of 3.6 and probit from 8.3, set the likelihood of every draw in the bank to zero.
Such a draw now keeps a finite likelihood, with each failure charged at the bound
of 2^-53, and the other draws are unaffected. The conditional mode is found where
the likelihood is convex, as it is for the link of an imperfect test, and the
curvature no longer underflows far into a tail. Final conditional sampling no
longer rejects a proposal whose probability rounds to one, which cut off the
samples at sites where every trial was positive. A custom link that returns a
value that is not a number, or is outside the range of the family, anywhere from
-50 to 50 is now refused when it is supplied. Fits with the canonical links are
unchanged.

- A custom inverse link now gets the checks of the canonical links. A nugget or
random-effect variance held at zero is tested by its slope there, and the check
for other maxima runs its screen of short scales and its exact alarms, which
take the link's terms in the likelihood with independent location effects.
Neither test needs the link to be log-concave. The sampler of the final latent
draws and the sampler behind predictions with parameter uncertainty now receive
the derivatives that the fit used, in place of deriving them again from the
inverse link. Errors while an importance distribution is built at a trial
parameter value, which give that value no weight, are counted and reported by
`check_fit()`. Fits with the canonical links are unchanged.

- The likelihood with independent location effects, behind the screen of
short scales and the exact alarms, is now right at every location for a custom
link. It integrated each location's effect by Gauss-Hermite quadrature at a
mode found by Newton steps from zero. Where a probability rounds to one, as
for cloglog from a linear predictor of 3.6 and probit from 8.3, the likelihood
of a location with a failure is flat, the search did not move, and the value
could be wrong by several log units without the change from 30 to 60 nodes
showing it. With a link that is not log-concave, such as that of an imperfect
test, the search could stop at a lower mode. A custom link now takes the rule
of the site table, the trapezium rule over the whole range where the integrand
is not negligible, which assumes no mode, and its recorded error is the change
when the step is doubled. Fits with the canonical links are unchanged.

- A failure of the global proposal while the check for other maxima evaluates
the likelihood far from the fitted maximum, as at the short scales of its
screen, no longer makes the stage repeat with the Laplace proposal. That
evaluation already uses the Laplace proposal at its own parameter values, and
the failure says nothing about the proposal at the maximum. Since the screen
runs for custom links, such failures in it made the first stage of a cauchit
fit of `loaloa` repeat with the Laplace proposal, which suits that link far
worse, and the fit was no longer certified. It is now certified at its first
stage, with the global proposal. The failures still count when the check keeps
the maximum of one of its refits, which is then the fit of the stage. Fits
without such a failure are unchanged.

- A refit of the global proposal in which a location's quadratic fit has too
few weighted quadrature nodes is now a failed refit for the canonical links
too, as it already was for custom links, so the stage repeats with the Laplace
proposal. It stopped inside the proposal with "missing value where TRUE/FALSE
needed", which gave the trial parameter value no weight. The case arises when
the proposal for a location with many tests lies many widths of its likelihood
away from it, as in refits from an anchor where the spatial variance is near
zero. A fit started there could not leave that variance, reported it as
effectively zero and certified only the regression coefficients, at a
log-likelihood 240 units below the maximum it now reaches. No fit in the
package's tests reaches this case, nor do fits of surveys from Liberia, Loa
loa, Malawi and Tanzania, and those fits are unchanged.

- A custom inverse link that `Deriv` cannot differentiate, such as a bare
`invlink = pnorm`, now has its derivatives computed for all values at once.
They were computed one value at a time, which made such a fit much slower than
the same link written as `function(eta) pnorm(eta)`.

- The numerical derivatives of a custom inverse link that `Deriv` cannot
differentiate are now accurate near a linear predictor of zero. Their steps
were relative to the linear predictor down to 1.8e-5, so just above it they
fell to about 2e-9, and the second derivative, the numerical derivative of the
first, was wrong by up to about 2,000 for linear predictors between about 2e-5
and 1e-2 (-40 at 1e-4 for `function(eta) 0.01 + 0.89 * plogis(eta)`, where it
is -1.1e-5). No step is now below 1e-4, and the steps are unchanged where the
linear predictor is at least 1 in absolute value. The second derivative enters
the curvature of the mode search and of the proposals but not the importance
weights, so estimates were valid, but fits with such a link took a different
path.

- A custom inverse link built from `plogis`, `pt` or `pcauchy`, such as
`function(eta) 0.01 + 0.89 * plogis(eta)` for an imperfect test or
`function(eta) pt(eta, 4)`, and a bare `invlink = pnorm`, `plogis`, `pcauchy`
or `exp`, now has exact derivatives, as a link built from `pnorm` or `exp`
already had. `Deriv` has no rules for these functions and cannot read the
compiled body of a bare function, so their derivatives were numerical.
Numerical derivatives are accurate but not smooth, and the search for
starting values by the Laplace approximation, which differences the
approximation, could stop short of its maximum. A call of these functions
with a location, scale, tail or logarithm other than the default is still
differentiated numerically. A link built inside another function, such as
`make_link <- function(sens) function(eta) 0.01 + sens * plogis(eta)`, is now
differentiated exactly too, where its derivatives were numerical because
they could not find the link's own variables. In 31 paired fits every outcome was the same and
the starting values were better in 24 of 30 and worse in none, and `plogis`
given as a function now reproduces the fits of the canonical logit link to
within 1e-6 standard errors.

- With a link not known to make the likelihood log-concave, cauchit or an
inverse link given in `invlink`, the global proposal now floors the curvature
of its factors smoothly at the bound that keeps its precision positive
definite, so that the proposal depends smoothly on the curvatures, as the
certificate assumes of its dependence on the parameters. The floor had a kink
where it starts to act, which it does in nearly every refit under cauchit and
the imperfect test. The bound comes from the smallest eigenvalue of the prior
precision, which keeps a kink where the smallest eigenvalues of two of its
blocks cross, as with `re()` terms. The eigenvalue is now computed once for
each value of the covariance parameters, where cauchit computed it for every
proposal and a link given as a function did not need it. Fits with the logit,
probit, cloglog and log links are unchanged. In 56 paired fits every certified
fit kept its certificate and its stages, and its estimates moved by at most
2e-6 standard errors.

- With `integrate_noise`, the global proposal's iteration now converges where
the curvature of a factor is small. It floored the candidates of its
acceleration with the same smooth floor as its refits, which moves a value it
has already floored, so wherever a curvature was below about 0.0115 the
iteration could not reach its tolerance, stopped at its limit of 500 refits and
was accepted as a stalled anchor. This happened in 70 of 171 integrated fits
measured, among them every fit of the Tanzania survey with a nugget, and under
`integrate_noise = "auto"` it went unrecorded in stages fitted in the sampled
form. The candidates are now floored with `pmax()` at zero, as those of links
not known to be log-concave are floored at their bound, so the proposal is the
same fixed point, now reached in a median of 9 refits. In the paired fits no
iteration stalled, and every certified fit stayed certified in the same stages
with its certificate within 1.3%. The estimates of certified fits moved by at
most 2.5e-4 standard errors, apart from three fits of one survey whose spatial
variance was held at its limit and one under an imperfect test (up to 0.08).
Imperfect-test fits that certified in neither version moved more, by up to
0.25 standard errors on a simulated survey at 5% prevalence and up to 5.6 on
degenerate ones, as much as a new seed moves them, and one of them was
reported as weakly identified after one stage rather than as having no maximum
after two.
A global iteration that stops just short of its tolerance is still accepted.

- Fits with the probit, cloglog or cauchit link named in `family`, such as
`family = binomial("probit")`, now compute their importance weights in
compiled code, as the canonical logit and log links do. Each weight takes half
the time it did in R, and such fits are 1.3 to 1.8 times faster on surveys of
200 to 400 sites, from 90 to 65 seconds for probit on the Tanzania survey of
387 sites and from 29 to 16 minutes for a cauchit fit that needs two stages.
The compiled weights repeat the R computation operation for operation, so fits
are unchanged (up to the last bits on platforms whose compilers fuse multiply
and add). An inverse link given as a function
in `invlink` still computes its weights in R.

- `glgpm()`, `specify_glgpm()` and `maxim_integrand()` now take a family object
in `family`, as `glm()` does, so a link can be named as in
`family = binomial("probit")`. Binomial models take the logit, probit, cloglog
and cauchit links, Poisson models the log link and Gaussian models the identity
link, and `family = binomial` is read as `binomial()`. The logit and log links
fit exactly as `family = "binomial"` and `family = "poisson"` do. A named link
gives the estimates of its inverse supplied in `invlink`, with exact
derivatives, so it is never differentiated during the fit and is faster than
an inverse such as a bare `pnorm`. `binomial("log")` and `poisson("identity")`,
whose means can leave the range of the family, `poisson("sqrt")`, whose inverse
is not monotone, quasi families, which have no likelihood, and other families
are refused with a message pointing to `invlink`, which still takes any other
inverse link as a function. A link given both in `family` and in `invlink` is
an error, and `invlink` does not take the name of a link. `print()` and
`summary()` name the link, and simulations, cross-validation refits and the
parametric bootstrap keep it. The undocumented list form of `invlink` now takes
only `inv`, `d1` and `d2`, and `specify_glgpm()` no longer accepts an
abbreviated family such as `"bin"`, as `glgpm()` never did. `summary()` writes
an inverse link taken from a package, such as a bare `invlink = pnorm`, by its
name instead of its compiled body, `.Call(C_pnorm, ...)`.

- `check_fit()` adds a row when a link is not known to make the likelihood of
each location log-concave in its latent effect, as the cauchit link and links
given as a function in `invlink` are not, since the latent effects can then
have more than one conditional mode. The row is `OK`, because this belongs to
the model and not to how it was fitted. A fit with a named link takes its
starting regression coefficients from a binomial regression with that link
instead of a logit regression, and falls back on the logit regression if that
one cannot be fitted.

- The row "Link" of `check_fit()` also says when it has not been measured on
simulated surveys how often the radius of the bound on the Monte Carlo error is
too small, as it has for the logit and probit links of binomial data and the
log link of Poisson and negative binomial counts, and points to
`check_fit(extra = "radius")`, which predicts it for the fit. A fit with the
cloglog link, which makes the likelihood log-concave, now has the row for that
reason. For a link given as a function in `invlink`, the row says where the
log-likelihood of a record was found to be convex on a grid of linear
predictors from -30 to 30, so that a link such as that of an imperfect test is
known not to be log-concave.

- A binomial model fitted without an intercept, such as
`y ~ 0 + x1 + x2 + gp()`, now takes its automatic starting regression
coefficients from a regression without an intercept. The start regression
dropped the first column of the design matrix and added an intercept, so the
starting coefficients were shifted by one column, and the first covariate
started from the intercept of another model. With the logit link the
transformed-Gaussian start usually replaced them, but with a named link such as
probit, or a link given in `invlink`, MCML started from them. In simulated
surveys this made fits with a nugget up to ten times slower, and with
`laplace_start = FALSE` it led one fit to a lower maximum. Fits of models with
an intercept are unchanged.

- `glgpm()` no longer compares automatic starting values for a binomial,
Poisson or negative binomial model with `re()` terms, and starts it from
RiskMap's defaults, as it always did in effect. The other start came from a
Gaussian model of the transformed records whose covariance was over all the
records, or over all the levels of the terms, and both starts were scored with
a dense matrix over all the latent effects. On simulated surveys of eight
households at each of 50, 100 and 150 locations, a fit took 77, 403 and 736
seconds with the comparison and 24, 35 and 12 seconds without it. On larger
surveys the comparison could ask for gigabytes at once, so whether a fit
compared the starts or fell back to the defaults depended on the memory
available. It also could never choose the other start, because the Laplace
approximation at that start carried the name of the `re()` term and was read
as missing. Fits of models with `re()` terms are unchanged.

- Without `re()` terms, the Gaussian model that gives the other automatic
start is fitted to one observation per location, the totals of the location
where it has several records, with the covariates averaged over its records.
Its covariance was over all the records, so at 1,000 locations with three
records each the comparison took 407 seconds against 80 for the fit, and now
takes 32. A survey with one record per location is fitted exactly as before.
With several records at a location the other start changes, and the empirical
logit of a total follows the proportion at the location, where that of a
binary record is always log(3) or -log(3).

- `fixed_pars` holds any parameter, a coefficient, the spatial variance or
scale, a random-effect variance, or every parameter. These stay in the working
vector at their values, which the optimisers keep by equal bounds and the
Gaussian paths by maximising over the others, and which the observed
information, the certificate and the radius check leave out, as they leave out a
variance held at zero. Their standard errors are NA. A held coefficient enters
the linear predictor as a known term, and a fit that holds its intercept
standardises its covariates without centring, so that every held coefficient is
one coordinate. The search for separation keeps held coefficients fixed, so a
held coefficient can resolve a separation. The check for other maxima is off
when the spatial variance or scale is held, and its screen of short scales and
exact alarm, whose likelihood with independent location effects is maximised
over every coefficient, are off when a coefficient is held. Held parameters add
no variance to the nugget's, whose standard error is that of its ratio to the
spatial variance when that variance is held. Profiles hold the held parameters,
and `profile()` refuses one.
With every parameter held, a Gaussian fit gives its exact log-likelihood at the
values, and a binomial, Poisson or negative binomial fit its Monte Carlo
log-likelihood on one sample of `validation_samples` draws from the global
proposal at the values, with conditional samples there and no certificate.
`check_fit()` says that nothing is estimated and that this log-likelihood's
Monte Carlo error is not assessed. `glgpm_parameters()` shows the held values.
`summary()` and `print()` say which parameters are held, and the AIC of a
Gaussian fit counts the estimated ones.

- `glgpm(fixed_pars = )` holds parameters at given values, with the names and
structure of `start_pars`. It holds the nugget variance (`tau2`, with
`gp(nugget = TRUE)`, positive), the measurement error variance of a Gaussian
model (`sigma2_me`, zero or positive) and the size of a negative binomial model
(`theta`). `fix_var_me` is removed, `gp(nugget = )` takes only `TRUE` or
`FALSE`, and `negative_binomial()` takes no size and stops with a pointer to
`fixed_pars` when given one. `MASS::negative.binomial(theta)` is still read as a
held size, which must agree with one in `fixed_pars`, and `specify_glgpm()`
takes the size as `theta` in its `parameters`. A start for a held parameter
stops unless it equals the held value, so `start_pars = coef(fit)` works for a
refit with the same holds. A fit stores `nugget` (`"none"`, `"estimated"` or
`"fixed"`) and `fixed_pars`, the held values, in place of `fix_tau2`,
`fix_var_me` and `family_parameters`, and fits saved before are read through
those fields. `coef()` gives every parameter, the held ones at their values, in
the order of the parameter table (`beta`, `sigma2`, `phi`, `tau2`, the family's
parameters, `sigma2_re`), where the family's parameters came before `tau2` and
held ones were left out. The columns of `parametric_bootstrap()` and the
estimates that `assess_simulation()` records are the estimated parameters only.
`profile()` names the nugget `tau2` and holds its variance for every family. For
a binomial, Poisson or negative binomial fit it profiled `nu2`, the ratio of the
nugget to the spatial variance, so the profiles of those fits move, and
`check_fit(extra = "profile")` names the nugget variance. The fits themselves
are unchanged.

- `start_pars` of `glgpm()` defaults to `NULL` and is a partial named list over
the model's estimated parameters, with the names of `?riskmap_parameters`.
`beta` is matched to the columns of the design matrix by name, or by position
when it is unnamed with one value per column, and `sigma2_re` to the `re()`
terms in the same way. The flat names of `unlist(coef(fit))`, such as
`"beta.x"`, are accepted too. Every other parameter takes its automatic start.
Before, any one value given switched off the comparison of automatic starts for
every parameter and left the others at the plain defaults, `beta` and
`sigma2_re` were matched by position with only their lengths checked,
`sigma2_me` was accepted and ignored for every family but the Gaussian, and
`beta = c(NA, 1)` passed the checks. Now the comparison still runs, each
candidate taking the given values in place of its own, unless every estimated
parameter is given or `set_control_mcml(start_comparison = FALSE)` skips it.
Every value must be finite, and positive but for the coefficients. A name the
model lacks, or one it holds, stops with the names it has, and a start of
`phi` below a hundredth of the smallest distance between locations or above
ten times the largest warns that it may be in other units. `summary()` says
which values were given (`mcml$starting_given`) and where the others came from.
`glgpm_parameters(formula, data, family)` lists the parameters of a model,
with their elements, flat names, labels, start rules and held values, before it
is fitted, and prints a template of `start_pars`. Fits without `start_pars` are
unchanged by this.

- Gaussian fits no longer depend on the units of the response. Every variance
started at 1 and the coefficients were in the units of the response, so a
response whose residual variance was far from 1 started far from the maximum,
and the optimiser's first steps, taken to bring the total variance to that of
the data, also took the scale of correlation to a value far below or above the
distances, where the likelihood is flat. On 20 simulated surveys of 50 locations
(`sigma2` 1, `phi` 20, an error standard deviation of 0.2) multiplied by 1e-6,
1e-3, 0.1, 1, 10, 100, 1000 and 1e6, the log-likelihood was more than 0.1 below
the maximum, or the fit failed, in 100, 15, 5, 0, 0, 55, 90 and 90% of them, by
up to 55 in the small units and 17.8 in the large ones. The engines now fit the residuals of the
regression divided by their standard deviation (`gaussian_regression()` and
`scale_gaussian_response()`), and the estimates, their covariance, the score and
the log-likelihood are returned in the units of the response. Held values and
starts keep those units. Every variance starts at the mean square of those
residuals (the rule `residual_variance` of `glgpm_parameters()`). In those 160
fits, and in models with and without a nugget, `re()` terms, an offset,
covariates, held parameters and given starts, the fits in every unit agree to
about 1e-7, and a response far from zero is fitted as the same one around zero.
The bounds of the variances in the profiles of a Gaussian fit, which cut the
profiles of a response multiplied by 1e-6 or 1e6 short, are on that scale too.
The regression that starts the coefficients now takes away the offset, which it
ignored, and the held coefficients, and is fitted on the standardised design, so
that a design such as a cubic in the calendar year is no longer refused as rank
deficient. One that is rank deficient stops even when `start_pars` gives every
coefficient, where it was fitted with a warning. The trace that
`messages = TRUE` shows is of the standardised response. A value of `sigma2_me` or `sigma2_re` in
`start_pars` without one of `sigma2` no longer starts `sigma2` too, which it did
in models of every family. Other fits of count models are unchanged.
On the golden battery 68 of the 82 tasks are identical, with every fit of counts.
In the 14 Gaussian ones the estimates move by at most 1.2e-7 of their value and
the log-likelihoods by at most 6e-12, except three bootstrap estimates of a
variance at zero (below 1e-9). In the Gaussian fits of the book's Chapters 3 and
5 the log-likelihoods move by at most 2e-8 and the estimates by at most 6e-5 of
their value, except two measurement error variances estimated at zero in the
Galicia profile over kappa (from 1.2e-10 to 7.8e-11 and from 6.1e-10 to
3.7e-10).

- On the Woodbury path of Gaussian fits, which `glgpm()` takes above three
records per latent effect unless the nugget is held at a positive value or the
measurement error variance at zero, the starts of the measurement error and random-effect
variances were divided by the start of the spatial variance, while its
likelihood reads them on their own scale. Only the nugget, which every engine
fits as its ratio to the spatial variance, should be divided. Every engine now
takes its start from the parameter table (`working_start()`). A fit whose start
of the spatial variance is 1, the default then, is unchanged. Fits given `start_pars`
change, and so do the refits of `parametric_bootstrap()`, which start from
`coef()`. On 20 simulated surveys of 8 records at each of 60 locations, with a
cluster effect, with and without a nugget, and started from `coef()` of the fit,
from a spatial variance of 3 or from one of 0.3, every fit converged with and
without the change, with about as many evaluations, its estimates moved by at
most 3e-4 standard errors and its log-likelihood by at most 1e-8. In bootstraps
of ten replicates of five of them, 292 of the 300 replicate estimates moved by at
most 1e-6 of their value, and the other 8, cluster variances estimated below
1e-4, moved along the nearly flat likelihood there. On the golden battery only
the bootstrap of the Gaussian fit with an `re()` term changes, its replicate
estimates by at most 2e-7 of their value except a cluster variance estimated near
1e-10.

- Certified MCML now builds every global proposal of a stage, other than the one
at its start, from three plain refits of the factors found at the start, in place
of iterating with Anderson acceleration until no factor moves by more than 1e-8.
The stopping rule of the iteration made the estimated likelihood jump wherever
the number of refits changed, by up to 4.5e-10 on the Tanzania survey, whereas
three refits are a smooth function of the parameters, as the theory of the
certificate assumes, and jumped by at most 2e-12 there. The proposal at the start
of a stage is still iterated to the tolerance. In a paired test on 127 simulated
surveys with ten runs each, every run was certified with and without the change,
the median ratio of the certificates was 1.000, the largest rise in the error was
0.002 standard errors, and the fits took 0.81 of the time. A refit that fails
counts as a failure of the global proposal, as before. Each stage of
`mcml_history` records `refit_change`, the largest change of a factor at the last
of the three refits.

- `set_control_mcml()` gains `integrate_noise`. With `TRUE`, a binomial or
Poisson model with a nugget, with any link, integrates the nugget out of the
likelihood of each location numerically, and certified MCML samples only the
spatial effect and any unstructured random effects. The model, its parameters,
its likelihood and `predict()` are unchanged, so the option changes only how the
same maximum is found. The integral is tabled for each distinct pair of positives
and trials on a grid of the linear predictor, by adaptive Gauss-Hermite
quadrature up to a nugget variance of 0.2 and the trapezium rule above it for
the canonical links, and by the trapezium rule with a step set by the narrowest
width of the integrand for a custom link, and read through a quintic Hermite
interpolant. Records at a location with the same covariates, offset and
random-effect levels are combined, and at a location that still has several
records the nugget stays sampled. A nugget held at zero is held at exactly zero,
where every record takes the plain kernel. On simulated surveys the integrated
form kept 3 to 6 times as many effective draws at the same parameters. In this
form a global iteration that stops at its limit of 500 refits with a last change
below 0.01 is accepted, instead of sending the stage to the Laplace proposal.
This is safe because the certificate rests on validation banks drawn after the
proposal is fixed, so a nearly settled anchor changes how efficient the sampler
is, not whether the estimate or its certificate is valid. On simulated surveys
of 700 sites with a spatial variance of 5 the iteration often stopped just short
of its tolerance of 1e-8, and the Laplace proposal then left almost no effective
draws. The stages where this happened record `stalled_anchors` and
`stalled_change` in `mcml_history`, summarised as `mcml_stalled_anchors`. The
sampled form keeps its rule, so its fits are unchanged. The
default of `laplace_start` becomes `NULL`, meaning `TRUE` for the sampled form
and `FALSE` with `integrate_noise = TRUE`. `check_fit()` says at how many
locations the nugget was integrated out, and suggests the option for a fit with
a nugget variance of at least 0.1 that was not certified or whose weights have a
heavy tail. Where the nugget variance is close to zero the likelihood is nearly
flat along the nugget, and the integrated estimate of the nugget's score is far
noisier than the sampled one, so there the integrated form certified fewer
fits (on simulated surveys of 700 sites, 5 against 12 of 21 runs below a nugget
variance of 0.1, 68 against 28 of 75 above 0.3). With
`integrate_noise = "auto"` both forms are built and each stage of the fit
integrates the nugget out when its variance at the start of the stage is at
least 0.1 and samples it otherwise, the first stage taking the variance from the
Laplace estimate. Each stage records its `form` in `mcml_history`, the fit
records them as `mcml_noise_forms`, and `check_fit()` reports them. With
`"auto"` every validation bank is drawn in both forms, and the two estimates of
the score are combined with weights inverse to their Monte Carlo covariances,
so the certificate rests on the more precise form or on both. A bank whose
importance weights have a Pareto k above 0.8 is left out of the combination,
since a heavy tail can make its estimated covariance far too small. On 96
simulated surveys of 700 sites with nugget variances of 0.15 and 0.8, three
runs each, the sampled form certified 124 of the 288 fits and `"auto"` 224.
Among the fits with a Markov chain reference, `"auto"` had no error above the
tolerance of 0.1 in the reference's information (222 fits), against two for the
sampled form (122). In every cell of the design it gave at least as many
certificates with an error within the tolerance as the sampled form, and more
in 31 of the 48 cells. It took 31 core-hours against 38. In a build before the
replicate bank's radius the sampled form certified 142 of the fits and the
integrated form 213, and without the ceiling on k the combination certified
245, 9 of them with an error above the tolerance. The stages record the combination in `combined`
of each validation bank in `mcml_history`. The extra check `"radius"` of
`check_fit()` rebuilds its samples in the form of the nugget of the fit's last
stage, and where that stage combined the two forms it resamples banks in each
form and combines them in pairs as the fit does, which takes about twice as
long. `"auto"` is the default, so fits of binomial and Poisson models with a
nugget change, while a model without a nugget is fitted as before.

- `integrate_noise` now also integrates out `re()` terms whose levels lie within
a single location, such as households in a village cluster, alongside the
nugget. Records at a location with the same covariates, offset and
random-effect levels are combined into one unit, and three rules decide each
level and location. A level that holds a single unit is integrated out of that
unit's likelihood, with the variances of all such levels added, and at a
location with a single unit so is the nugget. At a location with several units
the nugget is integrated out as well when the units share their covariates and
offset and differ only in levels that each hold a single unit, through an
integral over the nugget of the units' own integrals, tabled for each pattern
of units and pair of variances in the same way as the records. At any other
location the nugget stays sampled. Every other level stays sampled, including
levels that span several locations and levels that hold units with different
covariates, and a term can be integrated out at some levels and sampled at
others. The integrated levels leave the latent vector, so the sampled dimension
falls by one for each of them. A variance held at zero is held at exactly zero
in this form, as the nugget is, and the test of a held variance takes the slope
of a nested location's likelihood in its inner variance directly. With
`"auto"` a model with such a term is fitted as with `TRUE`, its `re()` terms and
its nugget integrated out at every stage wherever the rules allow, while a
model whose only such term is the nugget switches it by stage as before.
`check_fit()` says for each term at how many levels it was integrated out and
at how many it stayed sampled, and `summary()` gives them as `integration` in
place of `integrated_locations`. The radius check and `profile()` rebuild the
integrated form for such a fit. On
48 simulated surveys of 200 locations with about 800 households each, fitted
three times in each form, `TRUE`, and so `"auto"`, certified all 144 fits
against 128 sampled, at a median bound of 0.022 against 0.063, in 9% of the
sampled form's time. Against MCMC references the median error was 0.011
against 0.017 sampled, and no fit was off by more than the tolerance. Switching
the nugget alone with the households integrated certified 142 of 143, and on a
survey whose nugget variance was zero two of its three runs ended without a
certificate. On a survey built from `malkenya`'s households the bounds
fell from about 0.08 sampled to 0.007 integrated. A model without such terms is
fitted exactly as before.

- `integrate_noise` now also integrates out households whose members differ in
covariates, such as children of a household who differ in sex or birth
order, which were sampled before. At a location whose units differ in their
covariates or offset, the units are integrated out in groups when every term
but one has a single level shared by all of them and the levels of that one
lie within the location. Each level, a household for instance, is integrated
over its own effect with each member's own linear predictor, and the nugget
over the households' integrals, so the latent effects are the spatial effects
at the locations. Where every term is shared, the location is one group under
its nugget, which covers a nugget with child covariates and no `re()` term.
The integrals depend on the regression coefficients, so they are filled again
as these change, except that a household of one member takes the record's
own table, a household's integral is kept while only the coefficients of
covariates its members share change, and households alike in their members'
outcomes and covariate differences share one. On 200 clusters of a
Demographic and Health Survey, a fit with 8 child and household covariates,
a nugget and household effects took 23 minutes and one stage instead of 6
hours and three, on the same node, at a bound of 0.008 with 4,830 effective
draws against 0.086 with 1,350, with every estimate within 0.03 standard
errors. On a
simulated survey of that shape, with 28,471 children in 16,614
households at 1,375 locations, a fit with a nugget and household effects
certified at the first stage in 9 minutes on one core. On simulated surveys
of 200 locations with about 800 households of children who differ in sex, a
fit with a nugget and household effects took 36 to 45 seconds and one stage
against about 180 seconds and two stages with the households sampled. Under a
link not known to be log-concave, such as that of an imperfect test, the
households take the general trapezium rule, and small surveys fit more
slowly than with the households sampled (72 against 13 seconds on 40
locations). A fit whose units share their covariates and offset at every
location is unchanged.

- The integral over the independent effects of a record or a location, which
`integrate_noise` tables, no longer slows down as a variance grows. Records
with no events, or with every trial an event, have a side where the
likelihood is flat, and the trapezium rule walked across it a distance
proportional to the standard deviation, so at the upper bound of the search
for a log variance, 20, one evaluation took seconds to minutes, and a fit whose
optimiser stepped there could run for hours. Such a stretch is now taken in
closed form wherever the likelihood is linear across it, the steps of a
location's integral over its nugget grow with the standard deviation of its
inner effects, and with a custom or named link the parts beyond the range on
which the link is tabled, where its likelihood is constant, are taken in
closed form as well. The search for the mode of the integrand now finds a mode
far from the linear predictor, and the integral of a Poisson record whose
linear predictor is in the hundreds is now right, where it overflowed. Up to a
variance of about 100 the values are unchanged, and above it they agree with
numerical integration to about 1e-10.

- A Newton correction of a certified Monte Carlo maximum likelihood estimate
that would take a variance or the spatial scale past the upper bound of the
search is no longer applied. Such a correction follows a direction along which
the likelihood is nearly flat. The fit keeps its estimate and is reported as
not certified, with a warning, and `check_fit()` says why and suggests a larger
`max_samples`. Before, the parameter was clamped at the bound and the
certificate was computed there.

- `profile()` of a binomial or Poisson fit now rebuilds the likelihood with the
link of the fit. It used the logit or log link whatever link the fit had, so
the profiles of a fit with a named link such as `binomial("probit")`, or with
an inverse link given in `invlink`, were those of another model, and nothing
warned of it. The fix does not change the profiles of fits with the logit
and log links.

- `check_fit()` gains `extra`, which names the checks that cost computing time,
with `"all"` for every one, so that `check_fit()` is the one place to check a
fit (#262). The first is `"profile"`, the marginal profile likelihood from
`profile()` of every covariance parameter of the fit, which are the spatial
variance and scale, the nugget, an estimated measurement error variance and the
variances of the random effects, but not the regression coefficients. Its row
names each parameter whose profile is flat, so that its 95% interval is open at
an end, or has a second maximum, and is then `SUSPICIOUS`, and is a `PROBLEM`
when a point of a profile is above the estimate, which is then not the maximum
likelihood estimate. `cores` runs the profiles in parallel, with the same
results. Before the extra checks start, a message says how long they should
take, from the time the fit took, the number of parameters and the cores. Each
draws from a random number stream of its own after the seed of the fit, so the
same call gives the same result, and the random number state of the session is
left as it was. A profile already computed can be passed to `check_fit()` in
`profile`, which moves there from `plot()` of the checks, and is used as it is,
whatever parameters it covers. The checks returned keep the result of every
check that ran with the numbers behind its finding, and the new `summary()` of
the checks gives them, with the certificate and draws of each stage and the full
results of the extra checks, while `print()` still gives the verdict. `plot()` of
the checks gains `which` to choose its panels, draws only what was computed, and
says which extra check would compute a panel that is missing. Each profiled
parameter has its own panel, such as `"profile_phi"`, and `which = "profile"`
draws them all. The print of a fit without extra checks is unchanged.

- `check_fit()` gains the extra check `"radius"` for binomial and Poisson fits,
which predicts how often the radius in the bound on the Monte Carlo error is too
small (#262). The radius is estimated from the certification draws themselves,
which can miss the rare draws of large importance weight. Sixteen short Markov
chains on the latent effects reach those draws, and samples of certification
draws resampled from them and a fresh sample of the fit's size give the rate at
which the error that the radius allows for exceeds it, with a 95% interval,
against the 5% of the bound. On 35 simulated binomial and Poisson surveys the
prediction was within a factor of 1.65 of the rate measured on independent
samples. Outside binomial fits with the logit link it is read with a margin of
15%, since it ran low on the Poisson surveys. The row is `SUSPICIOUS` when the
radius, with that margin, would have to widen for the rate to be 5%, and then
gives the bound with the radius so widened. The chains can miss a separate mode
of the latent effects, so `OK` is not a proof. The check should take about
twice as long as the fit, and draws from the streams 9500 to 9502 after the
seed of the fit. `plot()` of the checks gains the panel `"radius"`, the error of
resampled samples against the radius the fit uses for them, whose subtitle
gives the predicted rate first.

- `print()` of the checks from `check_fit()` now gives a check that needs a
look in full with its fix, and a check that is `OK` on one line with its
finding in a few words, and ends by pointing to `summary()`, which now gives
every finding in full before the numbers behind it (#262). A fix that a costly
check would settle names it, `check_fit(fit, extra = "profile")` for a scale
that may be weakly identified and for other maxima, and
`check_fit(fit, extra = "radius")` for heavy-tailed importance weights, or
points to the row of that check when it ran. `plot()` of the checks gains two
panels drawn from what the fit recorded. `"held_variances"` shows the slope of
the log-likelihood at each variance held at zero with its interval of 2.58
standard errors, the test of the fit, and `"screen"` shows the screen of short
scales of the check for other maxima, the log-likelihood estimated at each short
scale and the exact lower bound there, against the estimate and the limit with
independent location effects. Each says why when it has nothing to show. A fit
now records the exact bound at each scale of the screen and, with the spatial
variance at zero and a nugget, over the scales and spatial shares, where it
kept only the largest, and `summary()` of the checks gives them. Estimates and
the decisions of a fit are unchanged.

- `plot()` of the checks from `check_fit()` returns one object and draws
nothing until it is printed (#213), as the other `plot()` methods do. One
panel is a ggplot and several are a patchwork in two columns. It used to print
the panels itself and return the list of them invisibly, so a figure could not
be saved or changed before it was drawn.

- `check_fit()` measures its evidence on other maxima from the estimate. When the
check for other maxima had kept another candidate than the first fit, or had
compared them on four times as many draws, "Largest exact lower bound", "Limit
with independent location effects" and the short scales were above the first fit
on the first draws, although they said above the estimate. In one fit without a
maximum the limit was given as 4.3 above the estimate and is 1.3 above it. In the
panel `"screen"` the gate is 0.25 below the first fit where another replaced it,
and the subtitle says "first fit" where the fit was repeated, also when an exact
bound repeated it behind a closed gate. The finding of an exact alarm no longer
says "at another spatial scale", since the bound can be that of the limit with
independent location effects. The warning of a fit that is numerically one of
independent location effects no longer gives a variance, which was that of the
limit and not of the estimate. Estimates are unchanged.

- `summary()` of the checks says "not assessed" where the likelihood was not
compared with the limit of independent location effects, instead of "no", and
answers the question at any scale when the fit recorded both log-likelihoods. It was answered
only for a scale below the shortest of the screen. A fit near the limit at a
scale longer than the spacing stays `OK`, and its finding says that the fitted
likelihood is not clearly above the limit. The status of every other fit is
unchanged.

- The panel `"screen"` of the checks names the limit with independent location
effects in the legend, instead of beside its line, where the name could overprint
the estimate or the bounds.

- Certified MCML no longer stops with "the optimiser returned non-finite
parameter values" when the search for the conditional mode of the latent
effects goes wrong. Each search starts from the mode found at the parameter
values tried before it. On a simulated survey at 1% prevalence a search had
accepted steps to -507, where the curvature of the logit link underflows to
0 / 0, and the next search started there. No step was then accepted, and the
search returned its start as the mode, with a Hessian that was not a number, as
if it had converged. 15 of the 48 searches of the Laplace step were broken this
way, which moved the Laplace estimate, and the importance distribution at that
estimate was not a number, so that every log weight was missing. A step to a
point where the likelihood terms are not finite is no longer accepted. A search
fails when its terms, its Newton step or its Hessian are not finite, or when no
step increases the objective at a point that is not a mode, and a failed search
from the previous mode is repeated from zero, so the next search starts from a
mode that was found. An importance distribution that is not finite is an error,
which gives that parameter value no weight, as other failures to build one do.
`check_fit()` reports how many searches started again from zero. The three runs
on that survey that stopped now fit, within 0.05 standard errors of a Markov
chain reference, and two are certified. The third reaches the largest budget
of 16,000 draws with a bound of 0.18. Fits whose searches never went wrong are
unchanged, as on two other surveys of the same study.

- `parametric_bootstrap()` simulates datasets from a Gaussian, binomial or
Poisson fit at its estimate, refits the model to each from the estimate, and
reports for each parameter the bootstrap median, bias, standard error and
percentile and bias-corrected (BC) intervals. The BC interval corrects for the
downward bias of the estimators of the spatial variance and the scale, which
shifts percentile intervals the wrong way. Each replicate has its own seed,
10000 apart from the next so that no two replicates share a Monte Carlo
stream, and the results are the same for any number of `cores`. Failed refits
are recorded instead of stopping the bootstrap. `plot()` draws a histogram of
each parameter with the estimate and both intervals.

- `plot()` of a profile now shows the relative profile log-likelihood, the
profile log-likelihood minus its maximum, which is 0 at the estimate and
negative elsewhere, so that the cut-off of a 95% interval is at -1.92 instead
of 3.84. `scale = "deviance"` draws the profile deviance as before. The panel
of a joint profile is titled "Joint profile", and its legend calls the two
ridges "best y for each x" and "best x for each y", since the axes name the
parameters. `plot()` of a parametric bootstrap draws each histogram on the natural scale of the
parameter, and gains `log_scale = TRUE` for the logarithmic axis of the
variances and the scale that it used before. Logarithmic axes of both plots are
labelled with plain numbers, such as 0.00001 and 1,000 instead of 1e-05 and
1e+03. `profile()` of a binomial or Poisson fit labels the parameters as
`summary()`, `parametric_bootstrap()` and `profile()` of a Gaussian fit do,
for example "Spatial process var." instead of "Spatial variance". The printed
intervals of a profile have 4 significant digits, with the lower and upper ends
of each interval under a common heading, so that the table fits 80 characters.

- `glgpm()` takes the CRS in which distances are measured in `distance_crs`,
which was `model_crs`, the name that `variogram()`, `summarise_distance()`,
`create_grid()` and `specify_glgpm()` use. Its arguments are in the order
`formula`, `data`, `family`, `invlink`, `denominator`, `distance_crs`,
`distance_units`, `start_pars`, `fixed_pars`, `control_mcml`, `control_mcmc`,
`return_samples`, `seed` and `messages`, so that `seed` and `messages` come
last as in the other functions. This changes calls that give the arguments
after `denominator` by position, and the order of the arguments in
`fit$call`. `summary()` of a fit takes `level` before `...`.
- `?glgpm` documents seven components of a fit as its interface, `data`, in
the CRS in which the model was fitted, `input_crs`, the CRS the data were
given in, `formula`, `family`, `call`, `seed` and `distance_units`. The other
components are internal and may change. The new `logLik()`, `vcov()` and
`nobs()` methods read the rest. `logLik()` has the number of estimated
parameters, without those that `fixed_pars` holds, as its `df` and the number
of records as its `nobs`, so that `AIC()` and `BIC()` work. `vcov()` gives the
covariance of the estimates on the scale of `coef()`, by the delta method from
the working scale, named as `unlist(coef(fit))` and for the estimated
parameters only, with `NA` rows for a parameter without a standard error.
`?glgpm` says that `update()` refits a model with some of its arguments
changed, and its example uses it.
- `print()` and `summary()` of a fit, `print()` of its summary and `print()`
of `negative_binomial()` stop on an argument they do not use instead of
ignoring it, and `summary(fit, conf_level = )` says to use `level`.
- `negative_binomial()` returns an object of class `RiskMap_family`, which
was `riskmap_family_specification`, the one result class not named in the
`RiskMap_` style.
- `maxim_integrand()`, `matern_gradient_phi()` and `matern_hessian_phi()` are
no longer exported, since only the fitting engine uses them, and their help
pages are gone. `create_ids()` is removed. `assess_prediction()` with
`keep_par_fixed = TRUE` indexes the locations of each training set with
`index_coordinate_rows()`, as `glgpm()` does, which gives the same indices in
time linear in the records instead of quadratic.
- `set_control_mcmc(n_draws = 2000, burnin = 2000, thin = 10, h = NULL,
c1_h = 1, c2_h = 0.6)`. `n_draws` is the number of draws kept, on which fitted
values and predictions rest, as in `profile(n_draws = )` and
`predict(parameter_draws = )`, and the chain runs `burnin + n_draws * thin`
iterations. It replaces `n_sim`, the number of iterations, which had to exceed
`burnin` by a multiple of `thin`. The default keeps 2,000 draws from 22,000
iterations, as before. `c1.h` and `c2.h`, the only dotted argument names of the
package, are `c1_h` and `c2_h`. `linear_model` is gone. It said that the draws
of a Gaussian model are independent, which prediction knows from the family,
and the parameter draws of a Gaussian prediction now take
`set_control_mcmc(n_draws = , burnin = 0, thin = 1)`. The settings recorded
with a fit's conditional samples, which `plot_mcmc()` reads, are `n_draws`,
`burnin` and `thin`, and the message of a prediction still gives the number
of iterations. A control saved from an earlier build, with `n_sim`, stops with
a hint to make it again.
- `set_control_mcml()` takes `initial_draws`, `validation_draws` and
`max_draws` in place of `initial_samples`, `validation_samples` and
`max_samples`, the words of `summary(fit)$mcml` and of the messages, and the
hints of `summary()` and `check_fit()` name them. A control stored with a fit
made before, with the old names, still serves its refits.
- Binomial and Poisson models with a nugget and `re()` terms none of whose levels
can be integrated out, as levels that span locations, stopped with the internal
error "A switching engine needs a sampled and an integrated form of the same
model" unless `integrate_noise` was `TRUE` or `FALSE`. The check of the two
forms compared the dimension of their latent effects with its name, which only
one of them had. Such a model is now fitted in either form, stage by stage, as a
model without `re()` terms is.

## Predictions

- Grids created by `create_grid()` now have the class `RiskMap_grid` before
`data.frame`, as vctrs requires of data frame subclasses. dplyr verbs such as
`filter()` previously failed with an error on these grids and on the
predictions made from them.
- `create_grid()` no longer raises an error when `boundaries` is in
longitude/latitude; like `glgpm()`, it now automatically reprojects to an
appropriate UTM zone and reports the conversion with a message.
- The old three-stage prediction interface has been replaced by
  `predict.RiskMap()`. Prediction locations, predictors, offsets, grouping
  variables and weights are supplied together in one `sf` object. Supplying
  `boundaries` requests blockwise joint prediction and areal aggregation;
  grids created by `create_grid()` use exact fractional cell coverage.
  Location- and area-level predictions support exceedance probabilities and
  named draw-by-draw transformations for additional scientific targets.
  Response predictions can be standardized over a common empirical or
  counterfactual reference population, with optional population weights.
- Gaussian prediction uses a stable reduced covariance solve, avoiding cancellation
with small measurement-error variances while retaining the same conditional model.
- `predict()` gains `parameter_uncertainty` (default `FALSE`) and
`parameter_draws` (default 50). When requested, parameters are drawn from the
approximate sampling distribution of the estimate and, for each draw, the
latent effects are sampled at that value with a short run of the same Langevin
sampler used for plug-in prediction. The pooled draws give intervals,
exceedance probabilities and areal summaries that also reflect parameter
uncertainty, which matters for small or low-prevalence surveys.
- The Langevin sampler of the latent effects, used by `predict()` and by the
final draw of `glgpm()`, is faster, most of all with `re()` terms. Its Laplace
preconditioner is factored through the Schur complement of the `re()` term
with most levels, and each iteration applies that factor without the prior
precision, so that it costs two triangular solves of the size of the
locations and passes over the records, instead of dense algebra in every
location and `re()` level. The conditional mode is found by Newton's method
with the same factor instead of `nlminb()` with a dense Hessian, and links
other than the canonical ones, including named links such as
`binomial("cloglog")`, run in the compiled chain. For a chain of 12,000
iterations a call takes 0.16 s instead of 0.28 s for 200 sites with one
record each, 1.8 s instead of 18 s for 1,000 sites, 0.5 s instead of 19 s for
200 sites with about 800 households as `re()` levels, and 4 s instead of 13
minutes for a survey of 46 clusters with 4,868 households. The chain
now starts from a draw of the Laplace approximation rather than at the mode.
From the mode a chain in many dimensions accepted almost nothing until its
step had shrunk, while from a draw it is in its stationary spread at once. The
step size kept after burn-in is the geometric mean of the tuned step over the
second half of burn-in, which keeps the acceptance rate closer to its target,
and a long chain can be interrupted. With `parameter_uncertainty = TRUE` the
mode at each parameter draw is found by the same search, started from the mode
at the previous draw, instead of the dense Laplace approximation of the fitting
engine, and each chain keeps the step tuned at the estimate with a burn-in of
100 iterations instead of 500, so that `predict()` takes 10 s instead
of 80 s for the 200 sites with households. A seed still gives the same
results every time, but they differ from those the same seed gave before,
although they follow the same distribution.
- `set_control_mcmc()` keeps 2,000 draws by default instead of 1,000, with
`n_sim = 22000` and the same burn-in of 2,000 iterations and thinning by 10.
`glgpm()` stores these draws and `predict()` reuses them, so the Monte Carlo
error of predictions is about 1.4 times smaller. Predictions on a grid hold
twice as many draws in memory and take two to two and a half times as long,
for example 8 s instead of 4 s for a new chain on the 8,740 cells of
`tz_covariates`.
Gaussian models also take 2,000 independent draws by default, including with
`parameter_uncertainty = TRUE`, which took 12,000, and with
`set_control_mcmc(linear_model = TRUE)`, which also took 12,000.
- `plot()` maps the output of `predict()`. Default colours depend on what the
  map shows: a calm-to-alarm palette for predicted levels, a single neutral
  hue for uncertainty (`sd`, `cv`, `interval_width`), a blue-to-red palette
  centred at 50% for exceedance probabilities and a palette centred at zero
  for signed effects. The legend title names the quantity (for example
  "Mean prevalence" for binomial models) and probabilities are labelled as
  percentages. `spatial_scale = "both"` draws locations and areas side by
  side with a shared colour scale, several summaries can be drawn in one
  call, `breaks` gives a binned map, and boundaries supplied to `predict()`
  are outlined automatically. `print()` lists what a prediction contains and
  how to map it. `predict()` output now records `response_name`.
- `create_grid()` gains `keep`. `keep = "overlap"` keeps every cell that
  overlaps the boundary, including edge cells whose centre lies outside, so
  that totals such as population along a coastline are not lost. The grid
  metadata now also records the study boundary.
- Prediction maps draw grid cells as squares clipped to the grid's study
  boundary, so further `sf` layers (for example survey locations) can be added
  to a map in any coordinate reference system. Legend titles can be replaced
  with `labs(fill = ...)`. `plot()` also maps the fixed-effect contribution
  with `quantity = "fixed"`, several summaries at one scale are arranged in a
  grid, and side-by-side scales no longer carry panel titles.
- `plot()` now respects `limits` for maps of signed quantities such as the
  spatial effect, which were previously rescaled, and arranges several
  summaries in two columns.
- Maps of `create_grid()` grids outline the grid's study boundary, so the
  edge of the study area is visible where the clipped cells end.
- Several maps drawn in one call have wider side margins, so neighbouring
  legends no longer run into each other at narrow widths.
- `plot()` accepts several quantities, for example
  `quantity = c("fixed", "spatial")`, and draws them side by side as it
  does for several summaries.
- New `assess_extrapolation()` reports, for each prediction location,
  whether a prediction relies on extrapolation: covariates outside their
  observed range on the model's scale, fixed-effect leverage relative to the
  data, and the share of the spatial variance not explained by the sampled
  locations. Locations are classed as supported, spatial, covariate, or
  spatial and covariate extrapolation, and `plot()` maps each indicator.
  Locations with missing covariates are drawn in grey without a legend key.
  The legend shows the colour of every class, including classes absent from
  the map, whose keys ggplot2 otherwise leaves blank.
- `assess_extrapolation()` includes the sampling noise of the outcomes in
  `spatial_unexplained` (#236). The share is
  v(x) = 1 - rho(x)' (R + W / sigma2)^-1 rho(x), with W diagonal and
  w_j = tau2 + c_j. For binomial, Poisson and negative binomial models c_j is
  the inverse of the summed curvature -d^2 log f(y | eta) / d eta^2 of the
  records at location j, at the conditional mode of the latent effects given
  the data at the estimates, the Laplace approximation that Emanuele Giorgi
  proposed. For Gaussian models it is the measurement-error variance divided
  by the number of records at the location. The share used to treat the
  spatial effect of binomial, Poisson and negative binomial models as observed
  without noise at every sampled location, and to take the measurement-error
  variance of a Gaussian model whatever the number of records, so Gaussian
  values change only at locations with several records. In the MBG book's
  Liberia example the mean share on the grid of Chapter 4 rises from 0.38 to
  0.63, close to the share computed from the draws of `predict()`.
  `sampling_noise = FALSE` gives the share without the noise, c_j = 0 for
  every family. The variances of `re()` terms are not included yet.
- The default `spatial_threshold` of `assess_extrapolation()` is `NULL`
  instead of 0.8. A location is then classed as spatial extrapolation when
  its share is above the 95th percentile of the shares left at the sampled
  locations by the data of the others,
  v_-j = 1 / [(R + W / sigma2)^-1]_jj - w_j / sigma2, which come from the
  same factor. A first version took their largest value, which a single
  remote sampled location sets, and which is 1 as soon as one location shares
  no correlation with the others. A number strictly between 0 and 1 gives a
  fixed threshold, where 0 and 1 were accepted before. The percentile is 1
  only when more than 5% of the sampled locations get no spatial information
  from the others, and then no location can be flagged and
  `assess_extrapolation()` warns with class `riskmap_warning`. The result
  keeps the threshold in the attribute `spatial_threshold`, with
  `spatial_rule` and `sampling_noise`, and has a `print()` method that counts
  the locations in each class and states the thresholds. In Liberia the
  threshold is 0.818 and 3.3% of the grid is flagged, where 0.8 without the
  noise flagged 0.2%. The example uses 150 locations of `italy_sim` instead
  of 50, whose locations share almost no correlation.
- `plot()` of predictions takes `summaries` and `thresholds` instead of
  `summary` and `threshold`, the names that `predict()` and `to_table()` use.
  Several thresholds give one exceedance map each, arranged as several
  summaries are, and `print()` of a prediction names the new argument. A call
  that names `summary` stops with an error naming `summaries`, while
  `threshold` still reaches `thresholds` by partial matching.
- `create_grid()` gains `raster`. The grid then consists of the cells of a
  projected `SpatRaster` with square cells: their centres are the grid points
  and the raster's layers become columns, so values such as population counts
  aggregated to the grid's resolution need no extraction. `boundaries` is
  optional in that case and, when given, selects cells according to `keep`.
- New `get_worldpop()` downloads WorldPop population counts for a country
  and year, for the whole population or for chosen age groups and sexes, at
  100 m or 1 km, from the unconstrained 2000-2020 or constrained 2015-2030
  estimates, and returns one raster of people per cell.
- `to_table()` is now generic and gains a method for `predict()` results with
  polygon summaries: one row per polygon, the mean and its interval in one
  column, probabilities as percentages and sums with thousands separators.
  All `to_table()` methods accept `col_names`. `predict()` output records the
  name of the `weights` column.
- `to_table()` for predictions formats each sum on its own, so that values in a
  column are no longer padded to a common number of decimals.
- Documentation and messages call phi the scale parameter of the spatial
  correlation function, rather than the spatial range, consistently with
  `summary()` and the MBG book. The practical range, where correlation falls
  to 0.05, is about 3 phi for the exponential correlation function.
- `predict()` calculates its location summaries for all locations at once
  instead of one location at a time. The median and the interval limits come
  from a single partial sort of each location's draws and are unchanged, and
  means and standard deviations agree with the previous values up to rounding
  in the last binary digit. On the 8,740 cells of `tz_covariates`, a
  `predict()` with the default summaries is about twice as fast, and one with
  every summary about three times as fast.
- A grid from `create_grid()` keeps its class and the metadata of its cells
when its columns are selected with `[`, `dplyr::select()`, `mutate()` or a
dplyr join (#214). Before, `predict()` with `boundaries` then assigned each
point to the polygon containing it without a word, instead of weighting cells
by their coverage, which changed areal estimates. `cbind()`, `merge()` and
`st_as_sf()` still drop the class, and so does keeping only the geometry
column. `rbind()` keeps the cells of its first grid. `predict()` now says
when it assigns points to polygons one by one although they lie on a square
lattice, as the centres of a grid's cells do.
- `predict()` takes `joint = NULL` in place of `type = c("auto", "marginal",
"joint")`, since `type` names the scale of the prediction in `stats::predict()`
methods and in `fitted()`. `TRUE` draws the locations jointly, `FALSE` draws
each location on its own, and `NULL`, the default, draws jointly with
`boundaries` and marginally without them, as `"auto"` did. `joint = FALSE` with
`boundaries` stops, as `type = "marginal"` did, and so does `include =
"random_effects"` without `joint = TRUE`. The result records the choice in
`joint`, `TRUE` or `FALSE`, instead of `type`, and its print reads it there.
`standardize` is spelt `standardise`, as RiskMap's own names are British
(`summarise_distance()`, `summarise_mcmc()`), and the result's `standardized`
is `standardised`. Predictions are unchanged.
- `predict()` and the print of a prediction stop on an argument they don't
use instead of ignoring it (`check_dots_unused()` in `checks.R`, which other
methods whose `...` is unused can share). The error names the argument, says
to use `joint` for `type`, which `predict(type = "response")` would otherwise
have ignored silently, and `standardise` for `standardize`, and otherwise asks
to check the spelling.
- `print()` and `plot()` of the indicators of `assess_extrapolation()` stop on
an argument they don't use. A subset that lost the classes or the threshold
still prints as an sf object and passes its arguments, such as `n`, to sf.
- The indicators of `assess_extrapolation()` and the table of
`glgpm_parameters()`, which are data frames, print as they do in the console
in R Markdown and Quarto documents whose `df-print` option shows data frames
as tables, without a call to `print()`. They have `knit_print()` methods. A
subset of rows of the indicators, or one that lost the classes, the threshold
or the geometry, is a table like any other data frame. Without the geometry
the indicators now also print as a data frame in the console, where a count of
the classes printed as if its rows were locations. A table of parameters
without the columns of the template prints as a data frame, where it stopped
with an error.
- Predictions of Gaussian models with `re()` terms read the codes of the groups
where they compared them with their labels, so the groups of labels that are
letters, factors or numbers other than 1, 2, ... were lost. The spatial
prediction, which conditions on all the data, then moved too, by up to 1.36
with a spatial variance of 2.26 in one example, and with numbers in another
order than the data's the effects went to the wrong groups. Labels 1, 2, ... in
the order in which they first appear gave the right result and still give the
same one. Cross-validation of these models changes with it.
- Joint Gaussian predictions with `include = "random_effects"` and without
`boundaries` had the sign of the covariance between the spatial and the random
effects reversed, stopped with three or more `re()` terms, and stopped with
`sigma2_me` held at 0, where the spatial prediction ignored the groups. Areal
predictions had the sign right and their draws are the same.
- Binomial, Poisson and negative binomial fits with `phi` or coefficients held in
`fixed_pars` drew their stored samples with the held values in the user's units
and everything else on the scale of the fit, with standardised covariates and
rescaled distances. The stored samples, `fitted()` and the default `predict()`
were then wrong, by up to 0.27 for a fitted prevalence in one example. The
estimates, the log-likelihood and the certificate, profiles, bootstraps,
predictions with `control_mcmc` and those with parameter uncertainty were not
affected. A fit saved with such a hold has to be fitted again, or predicted from
with `control_mcmc`.
- `predict()` without `newdata` gives each record the spatial effect of its own
location. With several records at some locations the effects of the locations
were recycled over the records, with a warning only when their numbers didn't
divide. The samples it keeps have one row for each record.
- The rows of `newdata` at one location share the spatial effect and the nugget
in joint predictions. Their nuggets were independent, and their spatial effects
were drawn from a singular covariance with a diagonal jitter, which warned. The
draws at the other locations of such a prediction change with the seed held,
and predictions at distinct locations don't.
- A prediction with `parameter_uncertainty = TRUE` keeps the estimate of a
variance, a scale or a size whose logarithm has a standard error above 3, and
draws the others given it. A Gaussian fit lets a variance estimated at zero drift
on the log scale, with a standard error of thousands, so that half of the draws
of a nugget estimated at zero were infinite and the prediction warned or stopped.
With standard errors between 3 and 10 the draws reached millions without a
warning. With a spatial variance at zero, the nugget is drawn by its own
variance. A warning names the parameters that keep their estimates, unless they
are variances below a hundredth of the total, or the scale of such a spatial
variance, and a fit with none to draw stops. The limit of 3 is provisional.
- A level of a random effect that isn't in the data of the fit stops a
prediction with `include = "random_effects"`, as before, and the help says so.
The engine behind it draws the effect of such a level from its distribution,
for the assessment of predictions (`allow_new_levels` of
`validate_prediction_data()`).
- `predict()` takes `"exceedance"` or `"cv"` as the only summary. The
components, such as the spatial effect, then have no summaries, and `print()`
doesn't list them. Both stopped. The call that `print()` suggests for the map of
an exceedance probability gives a threshold when there are several.
- A row of `newdata` that gives the model a covariate that isn't finite has `NA`
predictions, with a warning, as a row with missing values has. `log(0)` gave a
prevalence of exactly 0 and `log(-1)` an error about missing values in an
object of `model.frame()`. The indicators of `assess_extrapolation()` are `NA`
there too. A row of `standardise` that gives such a covariate stops.
- An offset that is zero at every record of the data is added in `newdata`. It
was ignored with a warning.
- The column of an exceedance probability is named by its own threshold, as the
help says, so `thresholds = c(0.2, 0.25)` gives `response_exceedance_0.2` and
not `response_exceedance_0.20`. Two thresholds never share a column.
- `plot_prediction()` finds a threshold that differs from one of the prediction's
by rounding, as `seq(0.1, 0.5, 0.1)[3]` does from 0.3.
- The help of `predict()` and the error for `include = "nugget"` at the observed
locations say that a Gaussian prediction there is of the spatial process without
the nugget, as at new locations, and that `fitted()` includes the nugget. They
said that the nugget was included, which holds for the conditional samples of
the other families. The predictions are unchanged.
- `messages = FALSE` silences the messages of `predict()` on the projection of
`newdata` and `boundaries`. The message on grid points that are counted whole
(#214) still shows, since it is about the result.
- Polygons that overlap once projected, because one has a vertex on the edge of
its neighbour in longitude and latitude, are no longer refused. They are
projected again with vertices added along their edges, so that the edges stay
together and no cell counts in two polygons. An overlap of up to a millionth of
the smaller polygon is rounding. Real overlaps and repeated polygons still stop,
and polygons that don't overlap when projected are projected as before.
- `assess_extrapolation()` computes the leverage from the QR decomposition of the
design, so `el + I(el^2) + I(el^3)` is no longer reported as singular. Other
leverages change in the last digits.

## Assessing a fit

- `assess_fit()` takes `models`, a fit or a list of fits to the same data,
`statistics`, `breaks`, `nsim`, `seed`, `return_simulations` and `messages`. It
simulates datasets from each model with `simulate_glgpm()`, with new surfaces and
the parameters at their estimates, and compares three features of them with those
of the data.
- The statistics are a function of a dataset, or a named list of functions, that
gives named numbers. The default ones are the mean and the SD of the prevalence
or rate and the proportion of zeros, or the mean, SD, minimum and maximum of a
Gaussian response. They are of the observations, and of the locations for
binomial models. Simulated datasets in which a statistic isn't finite are left
out of its summary, with a warning.
- The distribution is the empirical one of the same values, and the variogram
that of the residuals from the covariates. The response is transformed (empirical
logit, log rate with 0.5 added to the count, or the Gaussian response as it is),
with one value for each location, and the offset and the least squares fit on the
covariates are taken out, again in each simulated dataset. An offset under a link
other than the canonical one is one of the covariates of that fit. Both have the
95% joint envelope of `variogram()`, whose distance classes and pairs are
now in internal functions that both use, with no change to `variogram()`.
- The models share the seed, so their surfaces share random numbers, and the
simulations are kept only with `return_simulations = TRUE`. `summary()` gives the
observed value, the median and the 95% interval of the simulated ones and the
proportion at or above the observed one, `to_table()` tabulates them, and
`plot_predictive_check()`, or `plot()`, draws the statistics, the distribution or
the variogram. There is no test and no mark that a model passes.
- `?assessment` says which function answers which question about a model.

## Cross validation

- `"user"` has been added as an option to `method` the parameter of
`assess_prediction()` replacing the previous behaviour where providing the
`user_split` parameter overrode any provided `method`
- In `assess_prediction()`, `n_size` has been renamed to `size` (for consistency
with `iter` and `fold`) and `which_metric` has been renamed to `metrics`.
Parameters have also been reordered: mandatory arguments first, then those
required per `method` (in the order `"cluster"`, `"regularized"`, `"user"`),
then the remaining optional arguments.
- `assess_prediction()` takes a `seed` that makes the random splits generated
for `method = "cluster"` or `"regularized"` reproducible, together with the
refits and predictions of each test set.
- `assess_prediction()`'s `plot_fold = TRUE` plots are now consistent across
`method`s: `"cluster"` is unchanged (`spatialsample::autoplot()`), while
`"user"` and `"regularized"` now plot the testing and training datasets.
- `assess_simulation()` now takes fitted objects returned by `glgpm()` or
  model specifications returned by `specify_glgpm()` as its `models`
  templates. `f_grid_target` has been replaced by `target`, a function of the
  mean that defaults to the identity (see Simulations), and `f_area_target`
  has been renamed to `area_summary`. Their returned values are now checked
  before they are used. Grid- and area-level objectives can be requested
  together in one run with `spatial_scale = c("grid", "area")`.
- In `plot_metric()`; its `which_score`/`which_model`
parameters have been renamed to `metric`/`model`.
- `assess_prediction()` reports the SCRPS as a loss, ½[1 + CRPS/c + log(2c)]
with c = E|X − X′|/2, the negative of the scaled CRPS of Bolin and Wallin
(2023), so that lower is better for the CRPS and the SCRPS alike (#230). Every
SCRPS changes sign and nothing else moves. The help page states the
orientation and that the CRPS of binomial, Poisson and negative binomial
models is in counts, and so do the print of the results, their `summary()` and
the legends of their maps.
- `PIT_area` is computed once, from the calibration curve of all test
observations together, which is the average of the test sets' curves weighted
by their sizes, instead of once per test set and averaged by `summary()`
(#230). It is stored in `model$<name>$PIT_area`, no longer in `metric`, with
`PIT_area_expected`, its mean over 200 sets of test observations simulated
from their own predictive distributions. An observation in several test sets
is simulated once, from its first predictive distribution, and scored under
each. The simulation uses a seed drawn after the scores, the same for every
model, so the CRPS and SCRPS are those of before and identical models get
identical values. Under perfect calibration the area is about 0.31/√n for n
test observations, so a per-test-set area mostly measured the size of the test
sets.
- `assess_prediction()` scores binomial, Poisson and negative binomial counts
with F(k), the average over the draws of the mean of the family's distribution
function, instead of with one count simulated for each draw. The simulated
counts had no seed, so the same fit got different scores in two positions of
`models` (mean CRPS 3.477 and 3.469 in one example). CRPS, SCRPS, the
calibration curve, `PIT_area` and `PIT_area_expected` of these families change
by the size of that noise. The draws of the linear predictor and the seeds are
those of before, and Gaussian results are identical bit for bit. F is
evaluated at every count up to 100 and on an adaptive grid beyond, where the
sum over a cell is in closed form, with a relative error of about 10^-5 in the
CRPS. With 2,000 draws the score of a count takes about 5 ms with 20 trials,
25 ms with 100, and 0.1 to 0.5 s for predicted counts in the thousands or
more, most where the draws of the mean are widely dispersed
(`count_distribution()`). A rare event whose simulated counts were all zero
has a SCRPS. A test set of counts is recorded as a failure if a draw of the
mean is not finite, is outside the range of the family or gives counts above
2^53, or if a binomial denominator is not a whole number. A binomial fit was
scored without such a draw, and its scores with such a denominator were wrong.
- The families have the trait `distribution_function()`, which `check_family()`
compares with the density, and `score_limit()`, the largest count that is
scored, instead of `score_support()`, which needed the simulated counts.
- `summary()` of a cross-validation gives `CRPS` and `SCRPS` averaged over every
test observation with a score, `PIT_area` and `PIT_area_expected`. Before, the
means of the test sets were weighted by their sizes, which counted the missing
SCRPS of a predictive distribution without spread. Its print gives the scores
of each test set and a note on what they mean. `view_all`
and the attribute `overall_averages`, the average across models, are gone.
`print()` of the result gives the averages, as before.
- `plot()` of a cross-validation returns one ggplot and draws nothing (#213).
It draws one metric, by default the calibration curve when it was computed
and the first score otherwise, and stops with the list of metrics when given
several. Several models are panels of one plot with one colour scale. The
average calibration curve is the curve of all test observations, from which
`PIT_area` comes, its axis is labelled "Probability level u", and `PIT_area`
is no longer mapped. `...` comes last and stops on any argument, so ggplot2
components are added with `+`. `combine_pit = TRUE` stops unless `pit_mode` is
`"average"`. dplyr and the cross-validation's use of gridExtra are gone.
- `assess_prediction()` takes its list of fits in `models` instead of `object`,
and `plot_fold = FALSE` by default, so that the splits are drawn only when
asked for. With `method = "user"` the number of test sets is that of
`user_split`, so `iter` need not repeat it, and a different `iter` stops.
`min_dist` is in the distance units of the fits, converted to metres, in which
`spatialEco::subsample.distance()` compares distances whatever the CRS,
instead of always in kilometres. Fits in kilometres, as in every book example,
give the same splits as before, fits in metres take `min_dist` in metres, and
models with different distance units stop. `user_split` given as one list with
`in_id` and `out_id`, as a data frame or with no test sets stops. spatialEco and spatialsample move
to Suggests, and `method = "regularized"` and `"cluster"` ask for them.
- `assess_simulation()` takes `simulation` instead of `obj_sim`, `aggregate`
instead of `area_summary`, `label` instead of `col_names`, `objective`
instead of `pred_objective` and `breaks` instead of `categories`, and calls the
scale of the grid `"location"` instead of `"grid"`, as `plot()` of predictions
does. Its arguments are `simulation`, `models`, `spatial_scale`, `target`,
`aggregate`, `boundaries`, `label`, `objective`, `breaks`, `control_mcmc`,
`seed` and `messages`. `aggregate` takes the names of `predict()`
(`"mean"`, `"median"`, `"sum"`, `"min"`, `"max"`) or a function of a
location-by-draw matrix and the weights, as `predict()` does, with equal
weights. The result holds `objective$location` and `objective$area` instead of
`pred_objective$grid` and `pred_objective$area`, and so does its `summary()`.
- The print and summary methods of both assessments stop on an argument they
don't use.
- `create_test_sets()` makes the test sets of a cross-validation, for an `sf`
object or a fit, and `plot_test_sets()` maps them or draws the distances of the
test locations from their training data, with those of a prediction grid. The
methods are `"leave_one_out"`, one location in each test set, `"random"`,
`"cluster"`, k-means on the projected coordinates with 20 starts, and
`"spaced"`, the former `"regularized"`, locations at least `min_dist` apart, in
the units of the model, drawn by a sampler that stops only when no location is
left. Records at one location are always together, `"random"` and `"cluster"`
take `repeats`, and data in longitude and latitude are projected as in
`glgpm()`. The clusters are those of k-means on the coordinates, where
spatialsample clustered a matrix of distances, so they differ from the former
ones. With `"regularized"` a test set could stop short of `size` locations
although more were available, and `size = 1` gave more than one.
- `assess_prediction()` takes `models`, `test_sets`, `refit`, `control_mcmc`,
`seed` and `messages`. `models` is a fit or a list of fits, and `test_sets` comes
from `create_test_sets()` or is a list of row numbers or of lists with `train`
and `test`, so the splits of other packages can be used. `refit = TRUE` replaces
`keep_par_fixed = FALSE`, and all metrics are always computed. The former
arguments `method`, `fold`, `min_dist`, `size`, `iter`, `plot_fold`,
`user_split`, `metrics` and `keep_par_fixed` stop and say what replaces them.
Models with different responses or denominators are refused, since their scores
can't be compared, and were scored before.
- `assess_prediction()` predicts each test set jointly, with the nugget and the
random effects, which were left out of the predictions and of the scores of
models with `re()` terms. A level that the training data don't have is drawn
from the distribution of the random effect. In a Gaussian example the mean CRPS
is 0.6535, where the exact one is 0.6544 and the former one was 0.7394. The
scores of models without `re()` terms are those of before within Monte Carlo
error, since the draws differ. The seeds of the test sets are those of before.
With 1,000 training locations and 2,000 draws, a test set of 2,000 locations takes
23 seconds for a binomial model, as before, and 1.6 for a Gaussian one (1.1), and
4,000 locations 4.8 (2.2), with the same memory.
- The reference for calibration is simulated from draws of the joint predictive
distribution of each test set, 999 times, or as many as there are predictive draws.
`PIT_area_expected` is the mean of the areas of these curves, and the new
`PIT_area_limit` their order statistic that a correct model is below with a
probability of 95%. The band is the joint envelope of the curves and of that of the
data (`global_rank_envelope()`, as in `variogram()`), and so is that of the
histograms. The observations were simulated independently before, which made the
reference too small with clusters, where the errors of a test set go together. With
one test set and the parameters at their true values (100 locations in a square of
side 100, range 20, 60 for training, 500 or 1,000 datasets, all four families), a
correct model is above the limit in 4.4 to 6.0% of the datasets, its curve leaves the
band in 4.5 to 5.8% and its bars in 2.4 to 5.8%. The former band, for each level
separately, is left somewhere by 15 to 28% of the curves (5% for one binomial trial).
Test sets are independent of each other in the reference, so it is cautious with
several of them. With five or ten clusters a correct model is above the limit in 0 to
2.0% of the datasets and leaves the band in 0 to 3.8%, and with random or leave-one-out
test sets in none and at most 0.3%. In simulations
with ranges of 5 and 20, the area of a correct model is 0.76 to 0.98 times the
expected one with clusters, for Gaussian, binomial, Poisson and negative binomial
data, 0.80 and 0.90 when the clusters are repeated three times, 0.74 with random test
sets, 0.67 with leave-one-out and 0.60 with leave-one-out and one binomial trial. An
observation in several test sets has one uniform, its randomised PIT value in the
first of them.
- The result of `assess_prediction()` is a list with `scores`, a data frame with
a row for each model and test observation (`model`, `test_set`, `row`, `CRPS`,
`SCRPS`, `pit_lower`, `pit_upper`), `calibration`, with the areas, the
calibration curves with the 95% band of a correct model and the simulated curves
on 101 levels, `test_sets`, `refit`, `family` and `failures`. It held nested lists
for each model, metric and test set.
- `summary()` of a cross-validation is a data frame with a row for each model.
Observations in several test sets count once, with the average of their scores.
With `reference` it gives the differences of the scores from those of a model,
with standard errors that treat the observations as independent. The scores of
each test set are no longer printed. `to_table()` takes `reference` too.
- `plot_cross_validation()` takes `type`, `"calibration"` or `"score"`, `style`,
`"curve"` or `"histogram"`, `score`, `reference`, `model`, `bins` and `colour`,
instead of `metric`, `pit_mode`, `pit_test_set`, `combine_pit` and `model`. Each model
has a calibration panel with the 95% joint band of a correct model and the area with
its limit, and the histogram is the nonrandomised one of Czado, Gneiting and Held
(2009). With `reference`, the maps show differences of scores on a diverging scale,
green where they are better and orange where they are worse, which readers with red
and green deficiencies can tell apart. The scales of the maps end at the 98th
percentile.
- The plots of assessments are green (`#00857F` and paler versions of it) by default,
and their `colour` argument changes it.
- spatialEco and spatialsample are no longer suggested.

## Simulations

- `simulate_glgpm()` now accepts a fitted model or `specify_glgpm()` model,
  with `nsim`, `what`, `sample_locations`, `prediction_grid` and `seed`.
  `what = c("data", "surface")` simulates jointly at the exact union of both
  sets of locations. This replaces `simulate_surface()` and its nearest-grid
  approximation. The old simulation arguments and output structure are removed.
  Parameters remain fixed, offsets and fitted inverse links are retained, and
  Poisson responses use exposure times the inverse-link mean.
  Use `simulated_data()`, `simulated_surface()` and `simulated_values()` to
  extract results.

- `specify_glgpm()` takes the projected CRS in which distances are measured
  in `distance_crs`, as `summarise_distance()` and `variogram()` do, instead
  of `model_crs`. Its arguments are now `formula`, `data`, `family`,
  `parameters`, `invlink`, `denominator`, `distance_crs` and
  `distance_units`, with `invlink` moved from last, so calls that give the
  later arguments by position must name them. Simulations with a given seed
  are unchanged.

- `assess_simulation()` scores every model on the scale of the mean, each
  through its own inverse link. The truth is the simulated mean, which the
  generating model's inverse link gives, and the predictive draws of each
  refitted model are mapped to the mean through that model's inverse link.
  Generating models can therefore have any link, named or custom, and models
  whose links differ from the generating one, such as a probit model on
  simulations from a logit model, are scored on the same prevalence as the
  truth. A single function of the linear predictor was applied to every
  model before, which scored such a model through the wrong inverse link.
  The optional `target` is applied to the mean (the probability for a
  binomial model, the mean per unit exposure for a Poisson model and the
  value for a Gaussian model), so prevalence needs no argument, and other
  targets such as `qlogis` for the logit scale remain available.

- `plot()` of a simulation and of a `simulated_surface()` maps the simulated
  mean by default, the prevalence for a binomial model, instead of the linear
  predictor, and `quantity` maps any other component of `simulated_values()`
  that the surface has (`"linear_predictor"`, `"spatial_effect"`, `"nugget"`,
  `"group_effect"`). Both take the colours, legend and `limits` of `plot()` of
  predictions. `palette = NULL` gives the sequential palette of predicted
  levels for the mean and the palette centred at zero, with symmetric limits,
  for the linear predictor and the effects. A `"-"` before a palette's name
  reverses it, and `reverse_palette` is gone. The legend names the component,
  for example "Prevalence" with percentages, "Rate", "Mean" or "Spatial
  effect", and the maps lose their axes, as prediction maps have none.
  `simulated_surface()` keeps the family of the model as an attribute for the
  legend. The internal map helpers that only these methods used, with their
  unreachable north arrow and scale bar, are gone.

- The print of a negative binomial model from `specify_glgpm()`, and of
  simulations from it or from a fit, gives the size `theta` with the other
  parameters. It was left out.

- `assess_simulation()` puts the simulated response and the denominators in
  columns of its own, `.riskmap_response` and `.riskmap_denominator`, so a
  candidate keeps its covariates whatever their names. Columns called `y` and
  `units_m` were overwritten, and a covariate `y` was not asked for on the
  prediction grid, without a message. A covariate called `units_m` gave mean
  squared errors of 0.049 and 0.030 in one example, against 0.0069 and 0.0082
  under another name. A candidate with a transformed response, which was
  refitted without the transformation, stops.
- In the classification of `assess_simulation()` the lowest class includes the
  first break and values beyond the breaks belong to the end classes, for the
  truth and for the draws, as in the maps of `plot_prediction()`. A true value
  at the first break or beyond the breaks was left out of every metric, draws
  beyond the breaks belonged to no class, and draws all above the last break
  were classed in the lowest class. Classifications without such values are
  unchanged. The help states the rule, and that the predicted class is the most
  probable one, with ties to the lowest.
- When the prediction grid of `assess_simulation()` has sampled locations and
  the generating model has a nugget, the truth there has the nugget of the data,
  and the prediction had a new one. The draws are now of the spatial effect and
  that nugget given the data. They are the conditional samples of the latent
  effect of the site for a binomial, Poisson or negative binomial candidate, and
  draws from its closed-form distribution for a Gaussian one, joint with the
  other locations of the same area (`assessment_response_draws()`). The mean
  squared error at 60 sampled locations with equal variances of the spatial
  effect and the nugget is 0.048, which is that of the exact conditional mean,
  instead of 0.69. Draws at other locations are unchanged, apart from Gaussian
  candidates in areas that contain a sampled location, whose spatial effect is
  drawn again. When areas are assessed, rows of the prediction grid with the same
  coordinates share their new nugget, as the truth does, where each had its own.
  The other results are identical bit for bit.
- `assess_simulation()` stops before any fit when no area contains a prediction
  location, where the mean squared errors were `NaN`, and when a model to
  classify is called `"Class"`, whose results were lost. A covariate missing from
  the prediction grid is named, where the message failed with an error of cli.
  `summary()` counts the fits with warnings, such as those that found no maximum,
  in `fits`, and `print()` of the assessment and of its summary report them.
- `plot()` of a simulation and of a `simulated_surface()` draws points when the
  prediction locations are not on a regular grid, as for a grid in another CRS
  than the model's. They were drawn as tiles too small to see.

## Other

- Argument checks in the compiled MALA, likelihood and distance kernels have
been removed or moved to their R callers, following the contributing guide;
the R functions already validated everything the MALA and likelihood kernels
checked.
- `control_sim` parameters have been renamed to `control_mcmc`.
- `to_table()` now returns a directly renderable `knitr_kable` table, with an
  explicit `digits` argument that preserves trailing zeroes. This replaces the
  previous `xtable` return value. The estimate and limits of a parameter, or
  one metric for all the models, share one notation and one number of decimal
  places. When one of them is below `10^(-digits)`, they get as many decimals
  as the smallest needs for two significant digits. Scientific notation is for
  more than six decimals and for a value of a million or more.
- `summary()` of a model with a nugget gave the wrong confidence interval for
  the nugget variance: its standard error was that of the logarithm of
  tau2 / sigma2^2 instead of tau2. The interval now uses the delta method for
  log(tau2) = log(tau2 / sigma2) + log(sigma2).
- When one parameter has no standard error, as when a variance is estimated
  at zero, `summary()` and `to_table()` still give the confidence intervals of
  the other parameters. Previously all of them were missing.
- The `data_sf` parameter in `create_ids` is now `data`.
- The `conf_level` parameter in `summary.RiskMap` is now `level`
- `plot_metric()` and `plot_AnPIT()` have been merged into
  `plot_cross_validation()`, which `plot()` of a cross-validation also calls.
  `metric` takes `"PIT"` for the calibration curve as well as `"CRPS"` and
  `"SCRPS"`, and one plot draws one metric (see Cross validation).
- `plot_variogram()` is the function that `plot()` of a variogram calls (see
  Plots).
- `plot_sim_surf()` is now `plot_simulation()`, which `plot()` of a
  simulation also calls, as in `plot(sim_result, simulation = ...)`. It draws
  a `ggplot` map matching the prediction maps, with a `palette` argument.
  `simulated_surface()` returns an object of class
  `RiskMap_simulated_surface`, which `plot_simulated_surface()` or `plot()`
  maps.
- In `assess_prediction()`, the `"AnPIT"` metric is now `"PIT"` and the
  `AnPIT_area` score is now `PIT_area`, for every family.
- In `plot.RiskMap_cross_validation()`, the calibration-curve arguments `mode`,
  `test_set` and `combine_panels` are now `pit_mode`, `pit_test_set` and
  `combine_pit`.
- ggspatial and jsonlite move from Imports to Suggests, which leaves 20
packages in Imports. `plot()` of predictions calls
`rlang::check_installed("ggspatial")` when `north_arrow` or `scale_bar` is
`TRUE`, and `get_worldpop()` calls `rlang::check_installed("jsonlite")` before
it reads the WorldPop catalogue, which the 1 km age and sex layers of
2000-2020 do not need. In an interactive session both offer to install the
package. The example of `get_worldpop()` runs only in interactive sessions,
since `R CMD check --as-cran` runs `\donttest{}` examples and it downloads
from WorldPop. RiskMap needs rlang 1.0.0 or later, whose
`check_installed()` and `caller_arg()` it uses.
- Every result that may be saved and read again has the attribute
`riskmap_version`, a list with the `version` of RiskMap that made it, as a
`package_version`, and `sha`, the `RemoteSha` of its DESCRIPTION when RiskMap
was installed from GitHub and `NA` otherwise. Fits, predictions, profiles,
bootstraps, checks, test sets, cross-validations, simulation assessments, simulations,
simulation models, variograms and extrapolations have it. It is not printed,
and an object without it was made before RiskMap stamped its results, by an
earlier version or an early development build of 2.0.0, so that code which
reads saved objects can tell when one is too old.
- `to_table()` has methods for fits, their summaries, cross-validations and
their summaries, in place of a default method that recognised them by their
class, and its default method stops, naming the objects it accepts. The tables
are unchanged.
- The unused `importFrom()` of `Deriv::Deriv()`, `graphics::points()`,
`numDeriv::grad()` and `numDeriv::hessian()`, `rlang::sym()` and
`utils::flush.console()` are removed. The code calls Deriv and numDeriv with
`::`.
- dplyr is no longer imported, and spatialEco and spatialsample move to
Suggests, which leaves 17 packages in Imports, parallel among them. dplyr stays
in Suggests for the tests of grids with dplyr verbs.
- gridExtra is no longer imported, since patchwork arranges the panels of
`plot()` of a profile, which leaves 16 packages in Imports, parallel among
them.
- `profile()` takes `seed = NULL` by default instead of 1, as every other
function with a seed does. Without a seed, a profile of a binomial, Poisson
or negative binomial fit draws one integer from the session's random number
generator, so `set.seed()` before the call reproduces it, and returns it as
`settings$seed`. The draw is the only change to the session's stream. A
Gaussian profile, which has no bank, draws none. `seed` must be a
non-negative integer. A profile that fails leaves the session's stream as it
was. The profiles of `check_fit(extra = "profile")` keep their own seed,
derived from the fit's.
- `check_fit()` takes `cores` before `messages`, as `profile()` and
`parametric_bootstrap()` do, so that `cores` and `messages` come last in that
order everywhere, after `seed` and `level` where a function has them.
`parametric_bootstrap()` defaults to `cores = 1L`.
- `print()` of a simulation model, of simulations and of `glgpm_parameters()`
stops on an argument it doesn't use, and `print()` of a distance summary
raises a RiskMap input error for one, as the other print methods do.
- The `importFrom()` of `grid::unit()` and `stats::profile()` are removed,
since ggplot2 and stats are imported whole. The rules that `Deriv` uses for
the distribution functions of links are built with `quote()` instead of
`alist()`, which gives the same lists without the NOTE of `R CMD check` about
an undefined global variable `x`.
- The documentation is written in roxygen's markdown throughout (#54). Rd
macros such as `\code{}`, `\link{}`, `\emph{}` and `\itemize{}` became
markdown, links to functions read as `glgpm()`, and links to generics name
the RiskMap method or the package. The argument `method` of
`assess_prediction()` no longer wraps its list in `\describe{}`, which broke
the PDF manual of `R CMD check`.
- The help page of each function whose result has methods lists them under
See also (#178), and `print()` of profiles, of `glgpm_parameters()` and of
`negative_binomial()` is documented on the page of the function that makes the
object. Custom links, the negative binomial family, the penalties of
`set_control_mcml(pc_prior = TRUE)`, `crs_distortion()` and the spatial
indicator of `assess_extrapolation()` are marked experimental in their help
pages, as parts that may change in a minor release.
- The check that CRAN runs (`R CMD check --as-cran` without `NOT_CRAN`) is
kept short. The 70 tests that took more than a second, apart from the outputs
of binomial and Poisson fits, skip on CRAN, so its tests take about 80 seconds
instead of six to eight minutes. They all still run when `NOT_CRAN=true`, as on
CI, but not in R-hub's checks, which run as CRAN does. The examples of
`check_fit()`, `summarise_mcmc()` and `plot_mcmc()` leave `\donttest{}`, since
they run in under three seconds. The Gaussian example of
`parametric_bootstrap()` uses 10 simulations instead of 20. The examples that
stay in `\donttest{}`, which `--as-cran` runs too, use 10 simulations in the
binomial `parametric_bootstrap()` instead of 100, `plot_check()` checks the
radius with a short profile of the scale instead of every extra check, the
binomial example of `profile()` profiles the scale alone, and the joint
profiles of `profile()` and `plot_profile()` use a grid of 7 by 7. Each example
outside `\donttest{}` runs in under three seconds.
- RiskMap asks for R 4.5.0 or later instead of 3.5.0, because the current
version of Deriv, which it imports, needs it. When 2.0.0 was prepared, CI checked
R 4.5.3, which was oldrel-1, on pushes to main and weekly. The full test suite
last passed on R 4.1.3 with Deriv 4.1.6, and nothing below 4.5 is checked any
longer. Two changes made for older R stay. A fit that holds every parameter
doesn't fail on R before 4.4.0, where `rcond()` of its empty information raised a
LAPACK error, and two tests of Gaussian prediction compare values without the
dimnames that older versions of Matrix add.
