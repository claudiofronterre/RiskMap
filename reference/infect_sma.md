# West Nile virus pool tests for female \*Culex pipiens\* in the Sacramento Metropolitan Area

A dataset of PCR-tested mosquito pools used to summarize infection for
female \*Culex pipiens\* in the Sacramento Metropolitan Area
(Sacramento, Placer, El Dorado counties, California, USA). Each row
represents a tested pool at a specific location and date, with an
estimated pool size.

## Usage

``` r
data(infect_sma)
```

## Format

A data frame with one row per tested pool (for a total of 596 pools) and
5 variables:

- lon:

  Numeric. Longitude in decimal degrees (WGS84).

- lat:

  Numeric. Latitude in decimal degrees (WGS84).

- est_pool_n:

  Integer. Estimated number of mosquitoes in the pool (see Details).

- wnv_pos:

  Logical. Whether the pool tested positive for West Nile virus
  (`TRUE`/`FALSE`).

- date:

  Date. Pool collection date (local).

## Source

Constructed from `sample_collections` in the `vectorsurvR` R package.

## Details

Derived from the vectorsurvR sample datasets by:

1.  Filtering to female *Culex pipiens* pools with WNV testing within
    the SMA using point-in-polygon against the union of Sacramento,
    Placer, and El Dorado counties.

2.  Estimating pool size `est_pool_n` for each pool by:

    - Summing total female counts from nearby collection points within
      the *same week* and within a spatial radius (e.g., 2 km) to obtain
      \\T\_{\mathrm{near}}\\.

    - Counting nearby pools within the same week/radius to obtain
      \\m\_{\mathrm{near}}\\.

    - Setting `est_pool_n = round(max(1, T_near / m_near))`; if no
      nearby collections are found, a conservative fallback (default 25)
      is used.

3.  Coordinates are kept as WGS84 (EPSG:4326).

**Important:** `est_pool_n` is an estimate for demonstration and
pooled-modelling examples; it is not the recorded laboratory pool size.
