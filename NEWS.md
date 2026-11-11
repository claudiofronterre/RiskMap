RiskMap 2.0.0
=============

- This version introduces many breaking changes from v1. The package should not
be considered stable, but further breaking changes will be handled gracefully.
- Test coverage has been increased from 0 to 92 %.
- All functions have examples
- Many functions have been renamed or removed:

| v1 | v2 | 
|----------|----------|
| `assess_pp` |  `assess_prediction` |
| `assess_sim` | `assess_simulation` |
| `check_mcmc` | `plot_mcmc` |
| `compute_ID_coords` | `create_ids` |
| `convex_hull_sf` | `create_convex_hull` |
| `dast()` | Removed; will be replaced in a separate package |
| `glgpm_sim` | `simulate_glgpm` |
| `Laplace_sampling_MCMC` | Removed; the replacement is internal |
| `matern.grad.phi` | `matern_gradient_phi ` |
| `matern.hessian.phi` | `matern_hessian_phi` |
| `matern_cor` | `matern_correlation` |
| `maxim.integrand` | `maxim_integrand` |
| `plot_score` / `plot_AnPIT` | `plot.RiskMap_cross_validation` (the two are merged; by default it draws every plot the `RiskMap_cross_validation` object has data for) |
| `plot_s_variogram` | `plot.RiskMap_variogram` |
| `pred_over_grid` | `predict` |
| `pred_target_grid` | `predict` |
| `pred_target_shp` | `predict` with `boundaries` |
| `set_control_sim` | `set_control_mcmc` |
| `surf_sim` | `simulate_glgpm` |
| `s_variogram` | `variogram` |

## Datasets
- Apart from `liberia`, all datasets are now in an simple features (sf) format.

## Coordinate reference systems

- Datasets must have a coordinate reference system. When coordinates are in
latitude/longitude and `model_crs` is not supplied, they are now automatically
reprojected to an appropriate UTM zone. `predict()`, `simulate_glgpm()` and
`assess_prediction()` have been updated to match; `simulate_glgpm()`'s own 
`convert_to_crs`/`scale_to_km` arguments have likewise been replaced with
`model_crs`/`distance_units`.
- providing `model_crs` in longitude/latitude now raises an error.
- In `glgpm()`, `convert_to_crs` has been replaced with `model_crs` and
`scale_to_km` has been replaced with `distance_units` (`"km"` or `"m"`).

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
  `plot_variogram()` (now `plot.RiskMap_variogram()`) shows the envelope by
  default and warns when it was omitted with `n_permutations = 0`.
- The `bins` parameter in `variogram()` has been removed and replaced with `breaks`.

## Model fitting

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
their samples in `attr(fit, "mcmc")`. `component` now defaults to the
location with the smallest effective sample size, can select several
locations and can be given together with `check_mean = TRUE`. The mean and
`component` previously included the levels of unstructured random effects,
which are stored after the locations in `S_samples`. They now refer to the
locations only.

- The new `summarise_mcmc()` returns a tibble of the effective sample sizes of
the mean of the spatial random effects, of each location and of each level of
an unstructured random effect. They are computed with the initial monotone
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
`invlink = pnorm`, now has its derivatives computed for all values at once, with
the same values as before. They were computed one value at a time, which made
such a fit much slower than the same link written as `function(eta) pnorm(eta)`.

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

## Other

- Argument checks in the compiled MALA, likelihood and distance kernels have
been removed or moved to their R callers, following the contributing guide;
the R functions already validated everything the MALA and likelihood kernels
checked.
- `control_sim` parameters have been renamed to `control_mcmc`.
- `to_table()` now returns a directly renderable `knitr_kable` table, with an
  explicit `digits` argument that preserves trailing zeroes. This replaces the
  previous `xtable` return value. Values that are compared share one
  notation: when a value needs scientific notation, so do the estimate and
  limits of the same parameter, or the same metric for the other models.
- `summary()` of a model with a nugget gave the wrong confidence interval for
  the nugget variance: its standard error was that of the logarithm of
  tau2 / sigma2^2 instead of tau2. The interval now uses the delta method for
  log(tau2) = log(tau2 / sigma2) + log(sigma2).
- When one parameter has no standard error, as when a variance is estimated
  at zero, `summary()` and `to_table()` still give the confidence intervals of
  the other parameters. Previously all of them were missing.
- The `data_sf` parameter in `create_ids` is now `data`.
- The `conf_level` parameter in `summary.RiskMap` is now `level`
- `plot_metric()` and `plot_AnPIT()` have been merged into a single
  `plot.RiskMap_cross_validation()` method. `metric` and `model` now default
  to every metric/model available rather than being required, so
  `plot(cv_result)` alone produces every plot the object has data for.
  `metric` accepts `"PIT"` as well as `"CRPS"`/`"SCRPS"`/`"PIT_area"` to
  select the calibration curve. Each metric gets its own plot (or, with more
  than one model, its own small grid with one panel per model); with more
  than one metric, each one is drawn in turn and the full set is returned
  invisibly as a named list.
- `plot_variogram()` is now `plot.RiskMap_variogram()`, a method rather than
  a standalone function; call it as `plot(variogram_result)`.
- `plot_sim_surf()` is now `plot.RiskMap_simulation()`; call it as
  `plot(sim_result, simulation = ...)`. It now draws a `ggplot` map matching
  the prediction maps, with `palette` and `reverse_palette`
  arguments. `simulated_surface()` returns an object of class
  `RiskMap_simulated_surface` with a matching `plot()` method.
- In `assess_prediction()`, the `"AnPIT"` metric is now `"PIT"` and the
  `AnPIT_area` score is now `PIT_area`, for every family.
- In `plot.RiskMap_cross_validation()`, the calibration-curve arguments `mode`,
  `test_set` and `combine_panels` are now `pit_mode`, `pit_test_set` and
  `combine_pit`.
- Test coverage has been increased from 0 to 85 %.
