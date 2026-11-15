# Assess where predictions extrapolate beyond the data

Computes, for each prediction location, indicators of whether a
prediction from a fitted geostatistical model relies on extrapolation,
either in space (where the data say little about the spatial effect) or
in the covariates (outside the conditions observed in the data).

## Usage

``` r
assess_extrapolation(
  object,
  newdata,
  spatial_threshold = NULL,
  sampling_noise = TRUE
)

# S3 method for class 'RiskMap_extrapolation'
print(x, ...)
```

## Arguments

- object:

  A fitted `RiskMap` model.

- newdata:

  An `sf` object of prediction points with the covariates used by the
  model, typically a grid from
  [`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md).

- spatial_threshold:

  Threshold of the share of the spatial variance left unexplained above
  which a location is classed as spatial extrapolation. `NULL` (the
  default) uses the largest share left at a sampled location when its
  own data are left out (see Details). A number strictly between 0 and 1
  gives a fixed threshold.

- sampling_noise:

  Whether the share of the spatial variance left unexplained includes
  the sampling noise of the outcomes (the default). With `FALSE`, it
  depends only on where the data were collected and on the nugget, as if
  the spatial effect had been observed without noise at every sampled
  location.

- x:

  A `RiskMap_extrapolation` object, whose
  [`print()`](https://rdrr.io/r/base/print.html) method counts the
  locations in each class and states the thresholds of the classes,
  returning `x` invisibly.

- ...:

  Not used.

## Value

An `sf` object of class `RiskMap_extrapolation` with the rows and
geometry of `newdata` and the columns

- covariates_outside:

  Number of covariates outside their observed range.

- covariate_novelty:

  Largest distance outside an observed range, as a share of that range
  (0 when all are within range).

- most_novel_covariate:

  Name of the design-matrix column with the largest distance outside its
  range, or `NA`.

- covariate_leverage:

  Leverage relative to its maximum over the data.

- spatial_unexplained:

  Share of the spatial variance left unexplained by the data, \\v(x)\\.

- extrapolation:

  Factor with levels `"Supported"`, `"Spatial"`, `"Covariate"` and
  `"Spatial and covariate"`.

Rows with missing covariates have missing indicators. The attribute
`spatial_threshold` holds the threshold used, and `covariate_ranges` the
range of each covariate in the data. Map the result with
[`plot.RiskMap_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_extrapolation.md).

## Details

Covariates are assessed on the model's scale, that is, as the columns of
the fixed-effect design matrix. A model with `log(elevation)` is
therefore assessed on log elevation, and a spline in age through each of
its basis columns. For each column, the distance outside the range
observed in the data is expressed as a share of that range.

The leverage \\h(x) = d(x)^\top (D^\top D)^{-1} d(x)\\, where \\D\\ is
the fitted design matrix, is reported relative to its maximum over the
data. Values above 1 indicate a combination of covariates that is more
unusual than any in the data, even when each covariate is within its
range.

The share of the spatial variance left unexplained is \$\$v(x) = 1 -
\rho(x)^\top (R + \sigma^{-2} W)^{-1} \rho(x),\$\$ where \\\rho(x)\\
holds the estimated correlations between \\x\\ and the sampled
locations, \\R\\ their correlation matrix, \\\sigma^2\\ the spatial
variance, and \\W\\ a diagonal matrix with \\w_j = \tau^2 + c_j\\, the
nugget variance plus the sampling variance of the linear predictor at
the sampled location \\j\\. For a Gaussian model, \\c_j\\ is the
measurement-error variance divided by the number of records at the
location, and \\v(x)\\ is \\\mathrm{Var}\[S(x) \mid y\] / \sigma^2\\ at
the estimated parameters. For binomial, Poisson and negative binomial
models, \\c_j\\ is the inverse of the information \\-\partial^2 \log f(y
\mid \eta) / \partial \eta^2\\ of the records at the location, summed
over them and evaluated at the conditional mode of the spatial effect
given the data at the estimated parameters, so that \\v(x)\\ is the
Laplace approximation of \\\mathrm{Var}\[S(x) \mid y\] / \sigma^2\\.
With `sampling_noise = FALSE`, \\c_j = 0\\.

\\v(x)\\ is close to 0 where the data determine the spatial effect and
close to 1 where they carry no information about it. Unlike the distance
to the nearest sampled location, it accounts for the estimated scale of
spatial correlation, for clustering of the sampled locations and for how
much the data at each location say about the spatial effect, which is
little where few individuals were sampled or where a prevalence is close
to 0 or 1. The variances of
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
terms are not included in \\W\\, so for models with
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md)
terms \\v(x)\\ is lower than the share they leave.

A location is classed as covariate extrapolation when any covariate lies
outside its observed range or its relative leverage exceeds 1, and as
spatial extrapolation when \\v(x)\\ exceeds the threshold. By default
the threshold is the largest over the sampled locations \\s_j\\ of
\$\$v\_{-j} = \mathrm{Var}\[S(s_j) \mid y\_{-j}\] / \sigma^2 = 1 /
\[(R + \sigma^{-2} W)^{-1}\]\_{jj} - w_j / \sigma^2,\$\$ the share left
at \\s_j\\ by the data of the other locations. Leave-one-out
cross-validation predicts at each sampled location from the others, so a
location above this threshold has less spatial information than any
prediction that held-out data could check. A sampled location too far
from the others to share any correlation with them makes the threshold
1, and then no location is classed as spatial extrapolation, so a fixed
`spatial_threshold` is then more informative. This is the rule of the
covariate indicators, which compare each covariate with its observed
range and the leverage with its largest value in the data. The two
classes are kept separate because they call for different responses:
more sampling where data are sparse, and caution about the form of the
covariate effects outside the observed conditions.

## See also

[`plot.RiskMap_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_extrapolation.md),
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)

## Examples

``` r
data(italy_sim)
example_data <- italy_sim[1:50, ]
fit <- glgpm(y ~ log(pop_dens) + gp(), data = example_data,
             family = "gaussian", messages = FALSE)
extrapolation <- assess_extrapolation(fit, italy_sim[51:100, ])
extrapolation
#> ── Extrapolation ───────────────────────────────────────────────────────────────
#> Indicators at 50 locations.
#> 
#> Supported              40
#> Spatial                 0
#> Covariate              10
#> Spatial and covariate   0
#> 
#> ℹ Spatial: more of the spatial variance left unexplained than at any sampled
#>   location when its own data are left out (100%).
#> ! A sampled location shares no spatial correlation with the others, so no
#>   location is classed as spatial extrapolation. A number in `spatial_threshold`
#>   gives a fixed threshold.
#> ℹ Covariate: a covariate outside its range in the data, or a leverage above its
#>   largest value there.
#> ℹ `plot(x)` maps the classes, and `indicator = "spatial_unexplained"` the
#>   share.
table(extrapolation$extrapolation)
#> 
#>             Supported               Spatial             Covariate 
#>                    40                     0                    10 
#> Spatial and covariate 
#>                     0 
```
