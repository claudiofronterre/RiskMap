<!-- The home page of the website. README.md is the page shown on GitHub and
CRAN, so keep the two in step. The map comes from knitting README.Rmd. -->

# RiskMap

<p class="rm-pitch">Model-based geostatistics in R. Fit spatial models to
survey data, predict disease risk at locations and over areas, and quantify
its uncertainty.</p>

<div class="rm-actions">
<a class="btn btn-primary btn-lg" href="#get-started">Get started</a>
<a class="btn btn-outline-primary btn-lg" href="https://www.mbgr.org/"><i class="fa-solid fa-book-open" aria-hidden="true"></i> Read the book</a>
<a class="btn btn-outline-secondary btn-lg" href="reference/index.html">Reference</a>
</div>

::: {.rm-badges}
<!-- badges: start -->
[![R-CMD-check](https://github.com/claudiofronterre/RiskMap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/claudiofronterre/RiskMap/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/RiskMap)](https://CRAN.R-project.org/package=RiskMap)
<!-- badges: end -->
:::

RiskMap implements generalised linear geostatistical models for Gaussian,
binomial, Poisson and negative binomial outcomes. It was developed for public
health, for example to map the prevalence of a disease from surveys, and
applies to any outcome measured at geo-referenced locations.

<hr class="rm-rule">

## Installation

RiskMap 2.0 is not yet on CRAN, so install it from GitHub. It contains C++
code, which needs [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on
Windows and the Xcode command line tools on macOS.

``` r
# install.packages("remotes")
remotes::install_github("claudiofronterre/RiskMap")
```

`install.packages("RiskMap")` installs RiskMap 1.0.0 from CRAN, whose
interface is different.

## Get started {#get-started}

The `liberia` data give, for 90 villages, the number of people examined for
the nodules of river blindness and the number found with them. This example
fits a binomial geostatistical model to the counts and maps the predicted
prevalence across Liberia.

:::: {.row .g-4 .align-items-center}

::: {.col-xl-7}
``` r
library(RiskMap)
library(sf)
library(geobounds)

data(liberia)
liberia <- st_as_sf(liberia,
                    coords = c("long", "lat"),
                    crs = 4326)
crs <- propose_crs(liberia)
liberia <- st_transform(liberia, crs)
border <- st_transform(gb_get_adm0("Liberia"), crs)

fit <- glgpm(npos ~ gp(),
             data = liberia,
             family = "binomial",
             denominator = ntest)

grid <- create_grid(border, spacing = 2.5)
prediction <- predict(fit, newdata = grid)
plot(prediction)
```
:::

::: {.col-xl-5}
<img src="man/figures/README-example-1.png" class="rm-figure" alt="Map of the predicted prevalence of river blindness nodules across Liberia, from about 5% on the south coast to over 25% in the north-west.">
<p class="rm-caption">Predicted prevalence of nodules across Liberia.
Boundary from <a href="https://www.geoboundaries.org/">geoBoundaries</a>
(UNMIL and OCHA ROWCA), CC BY 3.0 IGO.</p>
:::

::::

## What RiskMap does

::::: {.row .row-cols-1 .row-cols-md-2 .g-3 .rm-cards}

:::: {.col}
::: {.card}
<i class="fa-solid fa-magnifying-glass-location rm-icon" aria-hidden="true"></i>

### Explore

Choose a projection that keeps distances accurate with `propose_crs()`, look
for spatial correlation with `variogram()`, and build prediction grids and
population weights with `create_grid()` and `get_worldpop()`.

[Data and projections](reference/index.html#data-and-projections)
:::
::::

:::: {.col}
::: {.card}
<i class="fa-solid fa-chart-line rm-icon" aria-hidden="true"></i>

### Fit

Write the model as a formula with a spatial Gaussian process, `gp()`, and
random effects, `re()`, and fit it with `glgpm()`. `check_fit()` flags
anything in the fit that needs attention.

[Specifying and fitting](reference/index.html#specifying-and-fitting)
:::
::::

:::: {.col}
::: {.card}
<i class="fa-solid fa-map-location-dot rm-icon" aria-hidden="true"></i>

### Predict

Predict prevalence at new locations or over areas with `predict()`, with
exceedance probabilities and averages weighted by population from
`get_worldpop()`.

[Prediction](reference/index.html#prediction)
:::
::::

:::: {.col}
::: {.card}
<i class="fa-solid fa-bullseye rm-icon" aria-hidden="true"></i>

### Assess

Compare models by spatial cross-validation with `assess_prediction()`, and
test prediction targets on simulated data with `assess_simulation()`.

[Assessment](reference/index.html#assessment)
:::
::::

:::::

## Learn with the book

::: {.rm-callout}
**[Model-based geostatistics for global public health using
R](https://www.mbgr.org/)**, by Emanuele Giorgi and Claudio Fronterre, teaches
model-based geostatistics with RiskMap through analyses of real data from
public health. It is free to read online.
:::

## RiskMap 2.0

RiskMap 2.0 is a redesign of version 1.0.0. Many functions have been renamed,
merged or removed, and prediction at locations and over areas now goes through
a single `predict()` method, so code written for 1.0.0 needs changing. The
[NEWS](news/index.html) lists the main changes, with a table of the functions
of 1.0.0 and their replacements.

## Citation

To cite RiskMap, please use

Giorgi, E., Fronterre, C. and Smart, S. (2026). *RiskMap: Geostatistical
Modeling of Spatially Referenced Data*. R package version 2.0.0.
<https://claudiofronterre.github.io/RiskMap/>

`citation("RiskMap")` in R, or [Citing RiskMap](authors.html#citation), gives
the same reference with a BibTeX entry. Please report bugs and suggest
features in the [issue tracker](https://github.com/claudiofronterre/RiskMap/issues).
