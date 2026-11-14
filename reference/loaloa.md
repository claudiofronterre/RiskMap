# Loa loa prevalence data from 197 village surveys

Prevalence of Loa loa (eyeworm) in a series of surveys undertaken in 197
villages in west Africa (Cameroon and southern Nigeria).

## Usage

``` r
data(loaloa)
```

## Format

An sf object with 197 rows and 10 variables:

- ROW:

  row id: 1 to 197.

- VILLCODE:

  village id.

- NO_EXAM:

  Number of people tested.

- NO_INF:

  Number of positive test results.

- ELEVATION:

  Height above sea-level in metres.

- MEAN9901:

  Mean of all NDVI values recorded at village location, 1999-2001

- MAX9901:

  Maximum of all NDVI values recorded at village location, 1999-2001

- MIN9901:

  Minimum of all NDVI values recorded at village location, 1999-2001

- STDEV9901:

  standard deviation of all NDVI values recorded at village location,
  1999-2001

- geometry:

  Simple feature geometry (POINT) giving the sampling coordinates in
  degrees latitude and longitude

## Source

Diggle, P.J., Thomson, M.C., Christensen, O.F., Rowlingson, B., Obsomer,
V., Gardon, J., Wanji, S., Takougang, I., Enyong, P., Kamgno, J., Remme,
H., Boussinesq, M. and Molyneux, D.H. (2007). Spatial modelling and
prediction of Loa loa risk: decision making under uncertainty. *Annals
of Tropical Medicine and Parasitology*, 101, 499-509.

## Details

The coordinate reference system is 4326.
