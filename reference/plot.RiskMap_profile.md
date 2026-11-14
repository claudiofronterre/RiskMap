# Plot profile likelihoods of a fitted geostatistical model

Draws the profiles computed by
[`profile.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md).
Each marginal profile shows the relative profile log-likelihood, the
profile log-likelihood minus its maximum, against the parameter. It is 0
at the maximum, which is the estimate unless a point on a profile is
higher, and negative elsewhere. The panel also shows the quadratic
approximation implied by the standard error (dashed), the cut-off for
the requested level and the profile interval shaded, with the estimate,
the profile interval and the Wald interval in the subtitle, one per
line. Each joint profile shows the relative log-likelihood over the grid
as filled regions at the 50%, 80% and 95% levels of a chi-squared
distribution on two degrees of freedom, the estimate, and the ridges of
the two marginal profiles (the value of one parameter that maximises the
likelihood for each value of the other). With one joint profile and both
its marginal profiles, the panels are arranged as a corner plot.
Variances and the scale are shown on a logarithmic axis.

## Usage

``` r
# S3 method for class 'RiskMap_profile'
plot(x, parameters = NULL, pairs = NULL, scale = c("loglik", "deviance"), ...)
```

## Arguments

- x:

  A `RiskMap_profile` object from
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

  Not used.

## Value

Invisibly, a list of `ggplot` objects, one per panel. The combined
figure is drawn.

## Details

With `scale = "deviance"`, the panels show the profile deviance instead,
twice the drop in the log-likelihood from its maximum, which is minus
twice the relative log-likelihood. The cut-off of a profile interval at
level \\1 - \alpha\\ is then the quantile \\\chi^2\_{1, 1 - \alpha}\\,
and on the relative log-likelihood it is \\-\chi^2\_{1, 1 - \alpha}/2\\
(about -1.92 at the 95% level).

## See also

[`profile.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md)
