# The Negative Binomial Family

Specifies a negative binomial response for
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
and
[`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md),
for counts that vary more than a Poisson model allows. A count \\y\\
with exposure \\m\\ given the random effects has mean \\\mu = m
\exp(\eta)\\, with the log link, and variance \\\mu + \mu^2 / \theta\\,
as in `MASS::negative.binomial(theta)`. The size \\\theta\\ is estimated
unless `fixed_pars = list(theta = )` of
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
holds it, and
[`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
takes it as `theta` in its `parameters`. The negative binomial family is
experimental and may change in a minor release.

## Usage

``` r
negative_binomial(...)

# S3 method for class 'RiskMap_family'
print(x, ...)
```

## Arguments

- ...:

  For `negative_binomial()`, not used, and a size given here stops,
  since it is held with `fixed_pars` instead. For
  [`print()`](https://rdrr.io/r/base/print.html), not used, and an error
  if given.

- x:

  A family from `negative_binomial()`, whose
  [`print()`](https://rdrr.io/r/base/print.html) method describes it and
  returns `x` invisibly.

## Value

An object of class `RiskMap_family`.

## Details

The exposure is given in the `denominator` of
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
and is 1 for every count when it is not given. As \\\theta\\ grows the
model tends to the Poisson model with the same mean, and a fit whose
size goes to that limit holds it there, as a variance estimated as zero
is held, and says so in
[`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md).
A fit that estimates the size reports it in
[coef()](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md),
[print()](https://claudiofronterre.github.io/RiskMap/reference/print.RiskMap.md)
and
[summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)
as `theta`. Its start is the size of the negative binomial regression
without the random effects, or `start_pars$theta` in
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).
With one observation per location, an estimated size can't be told apart
from a nugget or from a random effect with one level per observation, so
those are refused. A fit with an estimated size uses the sampled form of
Monte Carlo maximum likelihood (see `integrate_noise` in
[`set_control_mcml()`](https://claudiofronterre.github.io/RiskMap/reference/set_control_mcml.md)),
and its check for other maxima has no screen of short scales.
`MASS::negative.binomial(theta)` is read as the negative binomial family
with its size held at `theta`.

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
```
