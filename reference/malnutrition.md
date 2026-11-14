# Malnutrition in Ghana

Malnutrition in Ghana from the Demographic and Health Survey 2014.

## Usage

``` r
data(malnutrition)
```

## Format

An sf object with 2671 rows and 9 variables:

- age:

  age in months of the child.

- sex:

  sex of the child.

- HAZ:

  height-for-age Z-score.

- WAZ:

  weight-for-age Z-score.

- urb:

  binary indicator: urban area=1; rural area=0.

- etn:

  ethnic group.

- edu:

  level of education of the mother, which takes integer values from
  1="Poorly educated" to 3="Highly educated".

- wealth:

  wealth score of the household, which takes integer values from
  1="Poor" to 3="Rich".

- geometry:

  Simple feature geometry (POINT) giving the sampling coordinates in
  degrees latitude and longitude

## Source

Demographic and Health Survey, dhsprogram.com

## Details

The coordinate reference system is 4326.
