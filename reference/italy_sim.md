# Simulated data-set on the Italian peninsula

This is a simulated data-set over Italy for a continuous outcome. The
data-set contains 10 repeated observations for each of the 200
geo-referenced locations. The variables are as follows:

- x1 ordinate of the spatial locations.

- x2 abscissa of the spatial locations.

- y simulated continuous outcome.

- region the name of the region within which a given observation falls.

- province the name of the province within which a given observation
  falls.

- pop_dens the population density at the location of the observation.

- ID_loc an ID identifying the location to which the observation belong.

The coordinate reference system of the data is 32634.

## Usage

``` r
data(italy_sim)
```

## Format

A data frame with 2000 rows and 7 variables
