# RiskMap

Model-based geostatistics in R. Fit spatial models to survey data,
predict disease risk at locations and over areas, and quantify its
uncertainty.

[Get started](#get-started) [ Read the book](https://www.mbgr.org/)
[Reference](https://claudiofronterre.github.io/RiskMap/reference/index.md)

[![R-CMD-check](https://github.com/claudiofronterre/RiskMap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/claudiofronterre/RiskMap/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/RiskMap)](https://CRAN.R-project.org/package=RiskMap)

RiskMap implements generalised linear geostatistical models for
Gaussian, binomial, Poisson and negative binomial outcomes, fitted by
Monte Carlo maximum likelihood, or by exact maximum likelihood for
Gaussian outcomes. It was developed for public health, for example to
map the prevalence of a disease from surveys, and applies to any outcome
measured at geo-referenced locations.

------------------------------------------------------------------------

## Installation

RiskMap 2.0 is not yet on CRAN, so install it from GitHub. It contains
C++ code, which needs
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows and
the Xcode command line tools on macOS.

``` r

# install.packages("remotes")
remotes::install_github("claudiofronterre/RiskMap")
```

`install.packages("RiskMap")` installs RiskMap 1.0.0 from CRAN, whose
interface is different.

## Get started

The `liberia` data give, for 90 villages, the number of people examined
for the nodules of river blindness and the number found with them. This
example fits a binomial geostatistical model to the counts and maps the
predicted prevalence across Liberia.

``` r

library(RiskMap)
library(sf)
library(geobounds)

# Load the survey data as points
data(liberia)
liberia <- st_as_sf(liberia,
                    coords = c("long", "lat"),
                    crs = 4326)

# Choose a projection and project the data and the border
crs <- propose_crs(liberia)
liberia <- st_transform(liberia, crs)
border <- st_transform(gb_get_adm0("Liberia"), crs)

# Fit a binomial model with a spatial Gaussian process
fit <- glgpm(npos ~ gp(),
             data = liberia,
             family = "binomial",
             denominator = ntest)

# Predict prevalence on a 2.5 km grid and map it
grid <- create_grid(border, spacing = 2.5)
prediction <- predict(fit, newdata = grid)
plot(prediction)
```

![Map of the predicted prevalence of river blindness nodules across
Liberia, from about 5% on the south coast to over 25% in the
north-west.](reference/figures/README-example-1.png)

Predicted prevalence of nodules across Liberia. Boundary from
[geoBoundaries](https://www.geoboundaries.org/) (UNMIL and OCHA ROWCA),
CC BY 3.0 IGO.

## What RiskMap does

### Explore

Choose a projection that keeps distances accurate with
[`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md),
look for spatial correlation with
[`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md),
and build prediction grids and population weights with
[`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md)
and
[`get_worldpop()`](https://claudiofronterre.github.io/RiskMap/reference/get_worldpop.md).

[Data and
projections](https://claudiofronterre.github.io/RiskMap/reference/index.html#data-and-projections)

### Fit

Write the model as a formula with a spatial Gaussian process,
[`gp()`](https://claudiofronterre.github.io/RiskMap/reference/gp.md),
and random effects,
[`re()`](https://claudiofronterre.github.io/RiskMap/reference/re.md),
and fit it with
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).
[`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md)
flags anything in the fit that needs attention.

[Specifying and
fitting](https://claudiofronterre.github.io/RiskMap/reference/index.html#specifying-and-fitting)

### Predict

Predict prevalence at new locations or over areas with
[`predict()`](https://rdrr.io/r/stats/predict.html), with exceedance
probabilities and averages weighted by population from
[`get_worldpop()`](https://claudiofronterre.github.io/RiskMap/reference/get_worldpop.md).

[Prediction](https://claudiofronterre.github.io/RiskMap/reference/index.html#prediction)

### Assess

Compare models by spatial cross-validation with
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md),
and test prediction targets on simulated data with
[`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md).

[Assessment](https://claudiofronterre.github.io/RiskMap/reference/index.html#assessment)

## Learn with the book

**[Model-based geostatistics for global public health using
R](https://www.mbgr.org/)**, by Emanuele Giorgi and Claudio Fronterre,
teaches model-based geostatistics with RiskMap through analyses of real
data from public health. It is free to read online.

## RiskMap 2.0

RiskMap 2.0 is a redesign of version 1.0.0. Many functions have been
renamed, merged or removed, and prediction at locations and over areas
now goes through a single
[`predict()`](https://rdrr.io/r/stats/predict.html) method, so code
written for 1.0.0 needs changing. The
[NEWS](https://claudiofronterre.github.io/RiskMap/news/index.md) lists
the main changes, with a table of the functions of 1.0.0 and their
replacements.

## Citation

To cite RiskMap, please use

Giorgi, E., Fronterre, C. and Smart, S. (2026). *RiskMap: Geostatistical
Modeling of Spatially Referenced Data*. R package version 2.0.0.
<https://claudiofronterre.github.io/RiskMap/>

`citation("RiskMap")` in R, or [Citing
RiskMap](https://claudiofronterre.github.io/RiskMap/authors.html#citation),
gives the same reference with a BibTeX entry. Please report bugs and
suggest features in the [issue
tracker](https://github.com/claudiofronterre/RiskMap/issues).
