# Summarise the Checks on a Fitted Model

Summary method for objects of class `"RiskMap_check"` from
[`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md).
Where [`print()`](https://rdrr.io/r/base/print.html) gives the verdict,
[`summary()`](https://rdrr.io/r/base/summary.html) gives the evidence,
the finding of each check in full with the numbers behind it, the
certificate and the certification draws of each stage with the largest
change in the estimates from the stage before, and the full results of
the extra checks that ran. For the check of other maxima the numbers
include, at each short scale of the screen, the log-likelihood estimated
on the fitting draws and its exact lower bound, both above the estimate,
and for a fit with the spatial variance at zero and a nugget the exact
bound by scale and spatial share.

## Usage

``` r
# S3 method for class 'RiskMap_check'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"RiskMap_check"`.

- ...:

  other parameters.

## Value

An object of class `"summary.RiskMap_check"`, a list with the `checks`
table, the `results` of the default and extra checks in the order of the
table, the `stages` table with the column `largest_move`, and the
`draws`, `fit_seconds` and `not_computed` of `object`.

## See also

[`check_fit()`](https://claudiofronterre.github.io/RiskMap/reference/check_fit.md),
[`plot.RiskMap_check()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_check.md)

## Examples

``` r
data(italy_sim)

fit <- glgpm(
  formula = y ~ gp(),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)

summary(check_fit(fit))
#> ── Evidence for the checks on this fit ─────────────────────────────────────────
#> 
#> ✔ Monte Carlo accuracy
#>   A Gaussian fit has no Monte Carlo error.
#>   Nothing recorded beyond the finding.
#> 
#> ℹ These checks concern the computation, not whether the model suits the data.
```
