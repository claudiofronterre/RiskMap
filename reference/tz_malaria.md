# Malaria Dataset from Tanzania Demographic Health Surveys 2015

This dataset contains information on malaria prevalence and associated
variables from the 2015 Tanzania Demographic Health Surveys. The data
includes geographical, demographic, environmental, and health-related
variables.

- cluster.number:

  Cluster number, identifying the survey cluster.

- Lat:

  Latitude of the survey cluster.

- Long:

  Longitude of the survey cluster.

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

  Enhanced Vegetation Index, indicating vegetation density.

- Temperature:

  Average temperature in degrees Celsius.

- Precipitation:

  Total precipitation in millimeters.

- Population:

  Population density in the area (in thousands).

- ITN:

  Percentage of households with at least one insecticide-treated net
  (ITN).

- NTL:

  Nighttime light intensity, indicating urbanization and infrastructure.

- Urban.Rural:

  Indicator of whether the area is urban ('U') or rural ('R').

- utm_x:

  UTM (Universal Transverse Mercator) x-coordinate of the survey
  cluster.

- utm_y:

  UTM (Universal Transverse Mercator) y-coordinate of the survey
  cluster.

## Usage

``` r
data(tz_malaria)
```

## Format

A data frame with 387 rows and 20 columns, containing the following
variables: The CRS of the UTM coordinates is 32736.

## Source

[Tanzania Demographic Health Surveys 2015](https://dhsprogram.com),
Giorgi E, Fronterrè C, Macharia PM, Alegana VA, Snow RW, Diggle PJ. 2021
Model building and assessment of the impact of covariates for disease
prevalence mapping in low-resource settings: to explain and to predict.
J. R. Soc. Interface 18: 20210104.
<https://doi.org/10.1098/rsif.2021.0104>
