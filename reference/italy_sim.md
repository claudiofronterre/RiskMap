# Simulated dataset on the Italian peninsula

Simulated dataset over Italy for a continuous outcome.

## Usage

``` r
data(italy_sim)
```

## Format

An sf object with 2000 rows (10 repeated observations for each of the
200 geo-referenced locations) and 6 variables:

- y:

  simulated continuous outcome.

- region:

  the name of the region within which a given observation falls.

- province:

  the name of the province within which a given observation falls.

- pop_dens:

  the population density at the location of the observation.

- ID_loc:

  an ID identifying the location to which the observation belong.

- geometry:

  Simple feature geometry (POINT) giving the sampling coordinates in
  metres

## Details

The coordinate reference system is WGS 84 / UTM zone 33N (EPSG:32633),
which distorts distances by less than 0.5% over the data. The outcome
was simulated in UTM zone 34N (EPSG:32634), so its spatial correlation
is isotropic in those coordinates, and its range in the coordinates of
the data varies by up to about 1.4% over Italy.
