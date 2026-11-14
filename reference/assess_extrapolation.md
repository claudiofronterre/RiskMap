# Assess where predictions extrapolate beyond the data

Computes, for each prediction location, indicators of whether a
prediction from a fitted geostatistical model relies on extrapolation,
either in space (far from the data relative to the estimated spatial
correlation) or in the covariates (outside the conditions observed in
the data).

## Usage

``` r
assess_extrapolation(object, newdata, spatial_threshold = 0.8)
```

## Arguments

- object:

  A fitted `RiskMap` model.

- newdata:

  An `sf` object of prediction points with the covariates used by the
  model, typically a grid from
  [`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md).

- spatial_threshold:

  Share of the spatial variance left unexplained by the data above which
  a location is classed as spatial extrapolation. Defaults to 0.8.

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

  Share of the spatial variance not explained by the data.

- extrapolation:

  Factor with levels `"Supported"`, `"Spatial"`, `"Covariate"` and
  `"Spatial and covariate"`.

Rows with missing covariates have missing indicators. Map the result
with
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

The unexplained spatial variance is \\1 - \rho(x)^\top (R + \nu^2
I)^{-1} \rho(x)\\, where \\\rho(x)\\ holds the estimated correlations
between \\x\\ and the sampled locations, \\R\\ their correlation matrix,
and \\\nu^2\\ the ratio of the nugget (and, for Gaussian models, the
measurement-error variance) to the spatial variance. It is close to 0
next to well-sampled locations and close to 1 where the data carry no
information about the spatial effect. Unlike the distance to the nearest
sampled location, it accounts for the estimated scale of spatial
correlation and for clustering of the sampled locations. It depends only
on where the data were collected, not on the outcomes or on the number
of individuals at each location. For binomial, Poisson and negative
binomial outcomes, which carry less information than a noise-free
observation of the linear predictor, it is therefore an optimistic
(lower) value.

A location is classed as covariate extrapolation when any covariate lies
outside its observed range or its relative leverage exceeds 1, and as
spatial extrapolation when its unexplained spatial variance exceeds
`spatial_threshold`. The two are kept separate because they call for
different responses: more sampling where data are sparse, and caution
about the form of the covariate effects outside the observed conditions.

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
table(extrapolation$extrapolation)
#> 
#>             Supported               Spatial             Covariate 
#>                     0                    40                     0 
#> Spatial and covariate 
#>                    10 
```
