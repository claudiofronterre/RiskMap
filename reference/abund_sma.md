# Female Culex pipiens abundance (collections) in the Sacramento Metropolitan Area

A dataset of mosquito collection events used to quantify abundance of
female Culex pipiens in the Sacramento Metropolitan Area (Sacramento,
Placer, El Dorado counties, California, USA). Each row represents a
single collection event at a specific location and date.

## Usage

``` r
data(abund_sma)
```

## Format

A data frame with one row per collection event (for a total of 1552
collection events) and 6 variables:

- lon:

  Numeric. Longitude in decimal degrees (WGS84).

- lat:

  Numeric. Latitude in decimal degrees (WGS84).

- total_females:

  Integer. Total number of female Cx. pipiens captured in the event.

- date:

  Date. Exact collection date (local).

- trap_nights:

  Integer. Number of trap-nights for the event.

- trap_type:

  Character. Acronym for trap type used to collect mosquitoes. These
  include: `"NJLT"` (New Jersey Light Trap), `"GRVD"` (Gravid Trap), and
  `"MMT"` (Mosquito Magnet Trap).

## Source

Constructed from `sample_collections` in the `vectorsurvR` R package.

## Details

Derived from the vectorsurvR sample datasets by:

1.  Filtering to female *Culex pipiens* collection records within the
    Sacramento Metropolitan Area (SMA) using point-in-polygon against
    the union of Sacramento, Placer, and El Dorado county boundaries.

2.  Summing counts and trap-nights per unique *location*–*date*–*trap
    type*.

3.  Coordinates are kept as WGS84 (EPSG:4326).

This sample dataset is intended for examples, mapping, and vector-index
workflows. It is not official surveillance output. Dates span the same
period as the vectorsurvR sample (e.g., 2015–2021).
