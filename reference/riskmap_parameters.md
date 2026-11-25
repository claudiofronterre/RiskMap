# The parameters of a geostatistical model

Every parameter of a model fitted by
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
has one name, which `start_pars` of
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
the result of
[`coef.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md),
`parameters` of
[`specify_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/specify_glgpm.md)
(for the parameters a simulation takes) and the columns of
[`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)
share. The names a model has follow from its formula and family, so they
can be listed before it is fitted, with
[`glgpm_parameters()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm_parameters.md).

|  |  |  |
|----|----|----|
| Name | When the model has it | Elements and scale |
| `beta` | Always | One coefficient per column of the design matrix, named as the columns, such as `"(Intercept)"`, `"log(elevation)"` or `"typeurban"` for a level of a factor |
| `sigma2` | Always | The variance of the spatial process |
| `phi` | Always | The scale of the spatial correlation, in the model's `distance_units` |
| `tau2` | With `gp(nugget = TRUE)` | The variance of the nugget |
| `sigma2_me` | Gaussian family | The variance of the measurement error |
| `theta` | Negative binomial family | The size, see [`negative_binomial()`](https://claudiofronterre.github.io/RiskMap/reference/negative_binomial.md) |
| `sigma2_re` | With [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md) terms | One variance per grouping variable, named by the variable as written in [`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md) |

The flat names of `unlist(coef(fit))`, such as `"beta.log(elevation)"`
and `"sigma2_re.household"`, name a single element. `start_pars` accepts
them, in a list or as the named vector `unlist(coef(fit))`, and they
name the columns of
[`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md).

## Starting values

`start_pars` takes any of a model's estimated parameters, one, several
or all of them. `beta` can be a named subset of the coefficients, or
unnamed with one value per column. `sigma2_re` can be a named subset of
the variances, or unnamed with one value per term in the order of the
formula. Every value must be finite, and every one but the coefficients
positive. The others take their automatic start, given as `start` by
[`glgpm_parameters()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm_parameters.md).
`regression` is the coefficient of the regression without the latent
effects, `family_regression` the size of the negative binomial
regression without them, `residual_variance` the variance of the
residuals of a Gaussian model's regression, `one` the value 1 and
`distance_quantile` the 10th percentile of the distances between the
locations. A name the model lacks stops with the names it has.

## Held values

`fixed_pars` of
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
takes the same names and structure, and holds any of the parameters it
gives at their values instead of estimating them, such as
`list(phi = 20)`, `list(beta = c(x = 0.5))` or
`list(sigma2_re = c(household = 0.2))`. The nugget variance is held only
in a model with `gp(nugget = TRUE)`, and a nugget or random-effect
variance of zero is the model without that term.
[`coef.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md)
gives the held parameters with the estimated ones, so
`start_pars = coef(fit)` works for a refit with the same holds, and
[`summary.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)
and [`print()`](https://rdrr.io/r/base/print.html) show them as fixed,
without standard errors.
[`profile.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md)
has no profile of a held parameter and holds them in the profiles of the
others, and
[`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md)
refits with the same holds.

## See also

[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`glgpm_parameters()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm_parameters.md),
[`coef.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md)
