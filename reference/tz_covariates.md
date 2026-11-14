# Covariates Dataset for Malaria Prediction in Tanzania

This dataset provides covariates over a 10 by 10 km regular grid
covering Tanzania. It is intended to be used together with the
`tz_malaria` dataset for spatial prediction of malaria prevalence.

## Usage

``` r
data(tz_covariates)
```

## Format

An sf object with 8740 observations of 7 variables:

- Population:

  Population density in the area (in thousands).

- ITN:

  Percentage of households with at least one insecticide-treated net
  (ITN).

- EVI:

  Enhanced Vegetation Index, indicating vegetation density.

- Temperature:

  Average temperature in degrees Celsius.

- NTL:

  Nighttime light intensity, indicating urbanization and infrastructure.

- Precipitation:

  Total precipitation in millimeters.

- geometry:

  Simple feature geometry (POINT) giving the sampling coordinates in
  metres

## Source

Giorgi E, Fronterrè C, Macharia PM, Alegana VA, Snow RW, Diggle PJ.
(2021) Model building and assessment of the impact of covariates for
disease prevalence mapping in low-resource settings: to explain and to
predict. *J. R. Soc. Interface* 18: 20210104.
[doi:10.1098/rsif.2021.0104](https://doi.org/10.1098/rsif.2021.0104)

## Details

The coordinate reference system is 32736.
