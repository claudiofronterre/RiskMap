# Malaria Dataset from Tanzania Demographic Health Surveys 2015

Malaria prevalence and associated variables from the 2015 Tanzania
Demographic Health Surveys. The data includes geographical, demographic,
environmental, and health-related variables.

## Usage

``` r
data(tz_malaria)
```

## Format

An sf object with 387 rows and 18 columns, containing the following
variables:

- cluster.number:

  Cluster number, identifying the survey cluster.

- MM:

  Month of the survey (in two-digit format).

- YY:

  Year of the survey.

- UpAge:

  Upper age limit of the surveyed individuals in years.

- LoAge:

  Lower age limit of the surveyed individuals in years.

- Ex:

  Number of individuals examined for malaria.

- Pf:

  Number of individuals tested positive for Plasmodium falciparum
  (malaria parasite).

- PfPR2.10:

  Plasmodium falciparum parasite rate in the population (aged 2-10
  years).

- Method:

  Method used for malaria diagnosis (e.g., Rapid Diagnostic Test (RDT)).

- EVI:

  Enhanced vegetation index from MODIS (about 500 m), the mean over
  2015, without units, which indicates the density of vegetation.

- Temperature:

  Daytime land surface temperature from MODIS (band `LST_Day_1km`, about
  1 km), the mean over 2015, in degrees Celsius.

- Precipitation:

  Total precipitation over 2015 from CHIRPS (about 5 km), in
  millimetres. Missing for 2 clusters.

- Population:

  Population density from WorldPop at a resolution of 3 arc-seconds
  (about 100 m), in people per grid cell of about a hectare, which is
  about one hundredth of the density in people per square kilometre.

- ITN:

  Percentage of households with at least one insecticide-treated net
  (ITN), from a modelled surface whose source is not recorded.

- NTL:

  Night-time light from the VIIRS Day/Night Band monthly composites
  (about 500 m), the mean radiance over 2015, in nanowatts per square
  centimetre per steradian.

- Urban.Rural:

  Indicator of whether the area is urban ('U') or rural ('R').

- elogit:

  Empirical logit transformation of the positive cases -
  `elogit(Pf, Ex)`

- geometry:

  Simple feature geometry (POINT) giving the sampling coordinates in
  metres

## Source

[Tanzania Demographic Health Surveys 2015](https://dhsprogram.com),
Giorgi E, Fronterrè C, Macharia PM, Alegana VA, Snow RW, Diggle PJ.
(2021) Model building and assessment of the impact of covariates for
disease prevalence mapping in low-resource settings: to explain and to
predict. J. R. Soc. Interface 18: 20210104.
[doi:10.1098/rsif.2021.0104](https://doi.org/10.1098/rsif.2021.0104)

The covariates come from MODIS (EVI and land surface temperature), the
VIIRS Day/Night Band (night-time light), CHIRPS (precipitation) and
WorldPop (population density), at the coordinates of the clusters.

## Details

The coordinate reference system is 32736.
