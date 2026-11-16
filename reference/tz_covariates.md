# Covariates Dataset for Malaria Prediction in Tanzania

Covariates at the centres of the cells of a regular 10 by 10 km grid
covering Tanzania, for spatial prediction of malaria prevalence from the
`tz_malaria` dataset. Each covariate has the units it has in
`tz_malaria`, and EVI, Temperature, NTL and Precipitation come from the
same rasters.

## Usage

``` r
data(tz_covariates)
```

## Format

An sf object with 8740 observations of 7 variables:

- Population:

  Population density from WorldPop at a resolution of 3 arc-seconds
  (about 100 m), in people per grid cell of about a hectare, which is
  about one hundredth of the density in people per square kilometre.

- ITN:

  Percentage of households with at least one insecticide-treated net
  (ITN), from a modelled surface whose source is not recorded. Missing
  at 8 points.

- EVI:

  Enhanced vegetation index from MODIS (about 500 m), the mean over
  2015, without units, which indicates the density of vegetation.

- Temperature:

  Daytime land surface temperature from MODIS (band `LST_Day_1km`, about
  1 km), the mean over 2015, in degrees Celsius. Missing at 2 points,
  where the raster has no value.

- NTL:

  Night-time light from the VIIRS Day/Night Band monthly composites
  (about 500 m), the mean radiance over 2015, in nanowatts per square
  centimetre per steradian. Small negative values at 41 points are noise
  of the sensor.

- Precipitation:

  Total precipitation over 2015 from CHIRPS (about 5 km), in
  millimetres. Missing at 1 point, where the raster has no value.

- geometry:

  Simple feature geometry (POINT) giving the coordinates of the centres
  of the cells in metres

## Source

Giorgi E, Fronterrè C, Macharia PM, Alegana VA, Snow RW, Diggle PJ.
(2021) Model building and assessment of the impact of covariates for
disease prevalence mapping in low-resource settings: to explain and to
predict. *J. R. Soc. Interface* 18: 20210104.
[doi:10.1098/rsif.2021.0104](https://doi.org/10.1098/rsif.2021.0104)

The covariates come from MODIS (EVI and land surface temperature), the
VIIRS Day/Night Band (night-time light), CHIRPS (precipitation) and
WorldPop (population density).

## Details

The coordinate reference system is 32736. EVI, Temperature, NTL and
Precipitation are the values at the grid points of the rasters that give
these covariates in `tz_malaria`, and ITN is a percentage, as in
`tz_malaria`. Before RiskMap 2.0.0 they came from other products or
processing, ITN was a proportion and -9999 marked a missing value, so
that predictions from models with Temperature, NTL or ITN were wrong.
The script `data-raw/tz_covariates.R` in the source repository of
RiskMap made the change.
