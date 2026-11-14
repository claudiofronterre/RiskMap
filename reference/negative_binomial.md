# The Negative Binomial Family

Specifies a negative binomial response for
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
and
[`maxim_integrand()`](https://claudiofronterre.github.io/RiskMap/reference/maxim_integrand.md),
for counts that vary more than a Poisson model allows. A count \\y\\
with exposure \\m\\ given the random effects has mean \\\mu = m
\exp(\eta)\\, with the log link, and variance \\\mu + \mu^2 / \theta\\,
as in `MASS::negative.binomial(theta)`. The size \\\theta\\ is estimated
unless a value is given, at which it is held.

## Usage

``` r
negative_binomial(theta = NULL)
```

## Arguments

- theta:

  The size, a positive number at which it is held, or `NULL` (the
  default) to estimate it.

## Value

An object of class `riskmap_family_specification`.

## Details

The exposure is given in the `denominator` of
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
and is 1 for every count when it is not given. As \\\theta\\ grows the
model tends to the Poisson model with the same mean, and a fit whose
size goes to that limit holds it there, as a variance estimated as zero
is held, and says so in
[`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md).
A fit that estimates the size reports it in
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) as `theta`. Its start
is the size of the negative binomial regression without the random
effects, or `start_pars$theta` in
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).
With one observation per location, an estimated size can't be told apart
from a nugget or from a random effect with one level per observation, so
those are refused. A fit with an estimated size uses the sampled form of
Monte Carlo maximum likelihood (see `integrate_noise` in
[`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md)),
and its check for other maxima has no screen of short scales.
`MASS::negative.binomial(theta)` is read as `negative_binomial(theta)`.

Given the size, the log-likelihood of a count is that of a binomial
count of \\y\\ successes in \\y + \theta\\ trials with the logit link
and the offset \\\log(m / \theta)\\, plus terms in \\y\\ and \\\theta\\,
which is how RiskMap fits it. The log-likelihood that a fit reports
leaves out \\y \log m - \log y!\\, as that of a Poisson fit does, so the
two can be compared.

## See also

[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`stats::poisson()`](https://rdrr.io/r/stats/family.html)

## Examples

``` r
negative_binomial()
#> The negative binomial family with the log link and its size estimated.
negative_binomial(theta = 2)
#> The negative binomial family with the log link and its size held at 2.
```
