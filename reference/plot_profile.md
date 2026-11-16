# Plot profile likelihoods of a fitted geostatistical model

Draws the profiles computed by
[`profile.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md),
as [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the
result does. Each marginal profile shows the relative profile
log-likelihood, the profile log-likelihood minus its maximum, against
the parameter. It is 0 at the maximum, which is the estimate unless a
point on a profile is higher, and negative elsewhere. The panel also
shows the quadratic approximation implied by the standard error
(dashed), the cut-off for the requested level and the profile interval
shaded, with the estimate, the profile interval and the Wald interval in
the subtitle, one per line. Each joint profile shows the relative
log-likelihood over the grid as filled regions at the 50%, 80% and 95%
levels of a chi-squared distribution on two degrees of freedom, the
estimate, and the ridges of the two marginal profiles (the value of one
parameter that maximises the likelihood for each value of the other).
With one joint profile and both its marginal profiles, the panels are
arranged as a corner plot. Variances and the scale are shown on a
logarithmic axis.

## Usage

``` r
plot_profile(
  x,
  parameters = NULL,
  pairs = NULL,
  scale = c("loglik", "deviance"),
  ...
)

# S3 method for class 'RiskMap_profile'
plot(x, ...)
```

## Arguments

- x:

  The result of
  [`profile.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md).

- parameters:

  Names of the marginal profiles to draw. The default is all of them.

- pairs:

  Names of the joint profiles to draw, as `"first:second"`. The default
  is all of them.

- scale:

  The vertical scale of the profiles: `"loglik"` (the default) for the
  relative log-likelihood, or `"deviance"` for the profile deviance.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  arguments of `plot_profile()`. For `plot_profile()`, not used, and an
  error if given, so that a misspelt argument is not ignored.

## Value

A `ggplot` object for one panel, and otherwise a `patchwork` object that
arranges the panels. Nothing is drawn until the object is printed, so it
can be saved or changed first. One marginal panel is
`plot_profile(x, parameters = "phi", pairs = character())`, since
`pairs` defaults to every joint profile.

## Details

With `scale = "deviance"`, the panels show the profile deviance instead,
twice the drop in the log-likelihood from its maximum, which is minus
twice the relative log-likelihood. The cut-off of a profile interval at
level \\1 - \alpha\\ is then the quantile \\\chi^2\_{1, 1 - \alpha}\\,
and on the relative log-likelihood it is \\-\chi^2\_{1, 1 - \alpha}/2\\
(about -1.92 at the 95% level).

## See also

[`profile.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md)

## Examples

``` r
# \donttest{
data(galicia)
fit <- glgpm(log(lead) ~ gp(kappa = 1.5), data = galicia,
             family = "gaussian", distance_units = "km", messages = FALSE)
profiles <- profile(fit, pairs = list(c("sigma2", "phi")), n_grid = 7,
                    messages = FALSE)
plot_profile(profiles)

plot_profile(profiles, parameters = "phi", pairs = character(),
             scale = "deviance")

# The same panels
plot(profiles)

# }
```
