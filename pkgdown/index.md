<!-- The home page of the website. README.md is the page shown on GitHub and
CRAN, so keep the two in step. The map comes from knitting README.Rmd. -->

# RiskMap

<p class="rm-pitch">Turn surveys at scattered locations into maps of disease
risk, with geostatistical models that tell you how certain each map is.</p>

<div class="rm-actions">
<a class="btn btn-primary btn-lg" href="#from-data-to-a-map">Get started</a>
<a class="btn btn-outline-primary btn-lg" href="https://www.mbgr.org/"><i class="fa-solid fa-book-open" aria-hidden="true"></i> Read the book</a>
<a class="btn btn-outline-secondary btn-lg" href="reference/index.html">Reference</a>
</div>

::: {.rm-badges}
<!-- badges: start -->
[![R-CMD-check](https://github.com/claudiofronterre/RiskMap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/claudiofronterre/RiskMap/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/RiskMap)](https://CRAN.R-project.org/package=RiskMap)
<!-- badges: end -->
:::

RiskMap fits generalised linear models with a spatial Gaussian process to
Gaussian, binomial, Poisson and negative binomial data, and predicts risk at
single locations and over areas. It was written for public health, such as
mapping the prevalence of a disease from village surveys, but it applies to
any continuous or count outcome measured at points in space.

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

## From data to a map {#from-data-to-a-map}

The `liberia` data give, for 90 villages, the number of people examined for
the nodules of river blindness and the number found with them. A binomial
model with a Gaussian process, fitted to the counts, gives the predicted
prevalence everywhere in the area the villages span.

:::: {.row .g-4 .align-items-center}

::: {.col-xl-7}
``` r
library(RiskMap)
library(sf)

data(liberia)
liberia <- st_as_sf(liberia,
                    coords = c("long", "lat"),
                    crs = 4326)
liberia <- st_transform(liberia, propose_crs(liberia))

fit <- glgpm(npos ~ gp(),
             data = liberia,
             family = "binomial",
             denominator = ntest)

grid <- create_grid(create_convex_hull(liberia), spacing = 2.5)
prediction <- predict(fit, newdata = grid)
plot(prediction)
```
:::

::: {.col-xl-5}
<img src="man/figures/README-example-1.png" class="rm-figure" alt="Map of the predicted prevalence of river blindness nodules over the area of the Liberian villages, from under 10% in the south to over 30% in the north-west.">
<p class="rm-caption">The predicted prevalence of nodules, as drawn by <code>plot(prediction)</code>.</p>
:::

::::

## What RiskMap does

::::: {.row .row-cols-1 .row-cols-md-2 .g-3 .rm-cards}

:::: {.col}
::: {.card}
<i class="fa-solid fa-chart-line rm-icon" aria-hidden="true"></i>

### Fit

Write the model as a formula with a Gaussian process, `gp()`, and random
effects, `re()`, and fit it with `glgpm()`. Monte Carlo maximum likelihood
enlarges its samples until the Monte Carlo error is small.

[Specifying and fitting](reference/index.html#specifying-and-fitting)
:::
::::

:::: {.col}
::: {.card}
<i class="fa-solid fa-magnifying-glass-chart rm-icon" aria-hidden="true"></i>

### Check

See how a fit was computed with `check_fit()`, and measure the uncertainty of
the estimates with `profile()` likelihoods and `parametric_bootstrap()`.

[Checking a fit](reference/index.html#checking-a-fit)
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
public health. It is free to read online. The methods are set out in Diggle
and Giorgi (2019).
:::

## RiskMap 2.0

RiskMap 2.0 is a redesign of version 1.0.0. Many functions have been renamed,
merged or removed, and prediction at locations and over areas now goes through
a single `predict()` method, so code written for 1.0.0 needs changing. The
[NEWS](news/index.html) describes every change and has a table of the
functions of 1.0.0 with their replacements.

## Citation

Run `citation("RiskMap")` in R, or see [Citing RiskMap](authors.html#citation),
for the reference to the package. If you use the methods, please also cite

Diggle, P. J. and Giorgi, E. (2019). *Model-based Geostatistics for Global
Public Health: Methods and Applications*. Chapman and Hall/CRC.

Please report bugs and suggest features in the [issue
tracker](https://github.com/claudiofronterre/RiskMap/issues).
