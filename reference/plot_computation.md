# Plot the Checks on a Fitted Model

Plots the checks of a fit from
[`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md),
as [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the
result does. It makes the figures of the checks that were computed, and
never computes anything itself, so an extra check runs once in
[`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md).
By stage of certified Monte Carlo maximum likelihood, the panels show
the bound on the Monte Carlo error against the tolerance
(`"certificate"`), the effective certification draws
(`"effective_draws"`), the tail shape of the importance weights against
its limit of 0.7 (`"pareto_k"`) and how far the estimates moved
(`"moves"`). Two more panels come from what the fit recorded. For a
variance held at zero whose slope was tested, `"held_variances"` shows
the slope of the log-likelihood in that variance at zero with its
interval of 2.58 standard errors each side, the test of the fit, so that
the likelihood clearly rises as the variance leaves zero when the whole
interval is above zero, and clearly falls when it is below. The panel
`"screen"` shows the screen of short scales of the check for other
maxima. On a log scale of distance it gives, at each short scale, the
log-likelihood estimated on the fitting draws above that of the
estimate, which repeats the fit from the best such scale when one comes
within 0.25 of it (the dashed line), and the exact lower bound on the
likelihood there, with the likelihood with independent location effects,
the limit as the scale goes to zero, as a dotted line. When another fit
replaced the one that was screened, the values are above that fit and
the dashed line is 0.25 below the replaced one. An exact bound clearly
above the estimate shows that the estimate is not the maximum likelihood
estimate. For a fit with the spatial variance at zero and a nugget, the
panel adds the exact bound over longer scales for each share of the
variance in a spatial effect. With the extra check `"radius"`, the panel
`"radius"` shows the error of each resampled sample of certification
draws against the radius the fit would use for it, the larger of its own
and a replicate's, in standard errors, with the line where they are
equal. Its subtitle gives the predicted rate at which the radius is
exceeded, as in the row of the check, and the share of samples whose
error exceeds the radius, as in
[`summary.RiskMap_check()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_check.md),
from many more pairs of samples than the panel shows. With the extra
check `"profile"`, each profiled parameter has a panel of its profile
likelihood, named `"profile_"` followed by the parameter, such as
`"profile_phi"`.

## Usage

``` r
plot_computation(x, which = NULL, ...)

# S3 method for class 'RiskMap_check'
plot(x, ...)
```

## Arguments

- x:

  The result of
  [`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md).

- which:

  Names of the panels to draw, among `"certificate"`,
  `"effective_draws"`, `"pareto_k"`, `"moves"`, `"held_variances"`,
  `"screen"`, `"radius"`, the panels of the profiles, such as
  `"profile_phi"`, and `"profile"` for every profile. The default is
  every panel that has something to show. A chosen panel with nothing to
  show gives a message saying why, and which extra check of
  [`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md)
  would compute it. If no panel has anything to show,
  `plot_computation()` stops with those messages.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  arguments of `plot_computation()`. For `plot_computation()`, not used,
  and an error if given, so that a misspelt argument is not ignored.

## Value

A `ggplot` object when one panel is drawn, and otherwise a `patchwork`
object that arranges the panels in two columns. Nothing is drawn until
the object is printed, so it can be saved or changed first.

## See also

[`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md),
[`summary.RiskMap_check()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_check.md)

## Examples

``` r
# \donttest{
data(liberia)
liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
fit <- glgpm(
  npos ~ log_elevation + gp(kappa = 0.5), denominator = ntest,
  data = liberia, distance_crs = 32629, family = "binomial", messages = FALSE
)
# A short profile of the scale keeps the example quick
profiles <- profile(fit, parameters = "phi", n_points = 7, n_draws = 1000,
                    seed = 1, messages = FALSE)
checked <- check_computation(fit, extra = "radius", profile = profiles)
#> ℹ Computing the predictions of how often the radius of the bound is too small,
#>   which should take about 4 seconds, from the 2 seconds the fit took.
plot_computation(checked)

plot_computation(checked,
                 which = c("certificate", "screen", "radius", "profile_phi"))

# The same panels
plot(checked)

# }
```
