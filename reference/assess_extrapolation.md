# Assess where predictions extrapolate beyond the data

Computes, for each prediction location, indicators of whether a
prediction from a fitted geostatistical model relies on extrapolation,
either in space (where the data say little about the spatial effect) or
in the covariates (outside the conditions observed in the data). The
spatial indicator, with its sampling noise and its default threshold, is
experimental and may change in a minor release.

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
  default) uses the 95th percentile of the shares left at the sampled
  locations when their own data are left out (see Details). A number
  strictly between 0 and 1 gives a fixed threshold.

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

  Not used, and an error if given, so that a misspelt argument is not
  ignored. A subset that lost the classes or the threshold prints as an
  `sf` object, or as a data frame if it lost the geometry, and the
  arguments go to its print method.

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

Rows with missing covariates, or with ones that give the model a
covariate that isn't finite, such as `log(0)`, have missing indicators.
The attribute `spatial_threshold` holds the threshold used, and
`covariate_ranges` the range of each covariate in the data. Map the
result with
[`plot_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_extrapolation.md),
or [`plot()`](https://rdrr.io/r/graphics/plot.default.html), which does
the same.

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
the threshold is the 95th percentile over the sampled locations \\s_j\\
of \$\$v\_{-j} = \mathrm{Var}\[S(s_j) \mid y\_{-j}\] / \sigma^2 = 1 /
\[(R + \sigma^{-2} W)^{-1}\]\_{jj} - w_j / \sigma^2,\$\$ the share left
at \\s_j\\ by the data of the other locations. Leave-one-out
cross-validation predicts at each sampled location from the others, so
the \\v\_{-j}\\ are the amounts of spatial information under which
predictions can be checked against held-out data, and a location above
their 95th percentile has less than nearly all of them. A percentile
rather than the largest value keeps a single remote sampled location
from setting the threshold. The percentile is 1 only when the data of
the other locations say nothing about the spatial effect at more than 5%
of the sampled locations, as when the locations share almost no spatial
correlation. No location can then be classed as spatial extrapolation,
so `assess_extrapolation()` warns, and a fixed `spatial_threshold` is
more informative. The two classes are kept separate because they call
for different responses: more sampling where data are sparse, and
caution about the form of the covariate effects outside the observed
conditions.

## See also

[`plot_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/plot_extrapolation.md),
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md),
and
[assessment](https://claudiofronterre.github.io/RiskMap/reference/assessment.md)
for the other assessments of a model.

## Examples

``` r
data(italy_sim)
example_data <- italy_sim[1:150, ]
fit <- glgpm(y ~ log(pop_dens) + gp(), data = example_data,
             family = "gaussian", messages = FALSE)
extrapolation <- assess_extrapolation(fit, italy_sim[151:200, ])
extrapolation
#> ── Extrapolation ───────────────────────────────────────────────────────────────
#> Indicators at 50 locations.
#> 
#> Supported              30
#> Spatial                20
#> Covariate               0
#> Spatial and covariate   0
#> 
#> ℹ Spatial: more than 82.7% of the spatial variance left unexplained, the 95th
#>   percentile of the shares at the sampled locations when their own data are
#>   left out.
#> ℹ Covariate: a covariate outside its range in the data, or a leverage above its
#>   largest value there.
#> ℹ `plot(x)` maps the classes, and `indicator = "spatial_unexplained"` the
#>   share.
table(extrapolation$extrapolation)
#> 
#>             Supported               Spatial             Covariate 
#>                    30                    20                     0 
#> Spatial and covariate 
#>                     0 
```
