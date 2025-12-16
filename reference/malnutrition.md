# Malnutrition in Ghana

This geostatistical dataset was extracted from the Demographic and
Health Survey 2014 conducted in Ghana.

- lng Longitude of the sampling cluster.

- lat Latitude of the sampling cluster.

- age age in months of the child.

- sex sex of the child.

- HAZ height-for-age Z-score.

- WAZ weight-for-age Z-score

- urb binary indicator: urban area=1; rural area=0.

- etn ethnic group.

- edu level of education of the mother, which takes integer values from
  1="Poorly educated" to 3="Highly educated".

- wealth wealth score of the household, which takes integer values from
  1="Poor" to 3="Rich".

The coordinate reference system is 3857.

## Usage

``` r
data(malnutrition)
```

## Format

A data frame with 2671 rows and 10 variables

## Source

Demographic and Health Survey, dhsprogram.com
