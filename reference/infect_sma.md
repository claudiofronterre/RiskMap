# West Nile virus pool tests for female *Culex pipiens* in the Sacramento Metropolitan Area

A dataset of PCR-tested mosquito pools used to summarize infection for
female *Culex pipiens* in the Sacramento Metropolitan Area (Sacramento,
Placer, El Dorado counties, California, USA). Each row represents a
tested pool at a specific location and date, with an estimated pool
size.

## Usage

``` r
data(infect_sma)
```

## Format

An sf object with one row per tested pool (for a total of 596 pools) and
4 variables:

- est_pool_n:

  Integer. Estimated number of mosquitoes in the pool (see Details).

- wnv_pos:

  Logical. Whether the pool tested positive for West Nile virus
  (`TRUE`/`FALSE`).

- date:

  Date. Pool collection date (local).

- geometry:

  Simple feature geometry (POINT) giving the sampling coordinates in
  degrees latitude and longitude

## Source

Constructed from `sample_collections` in the `vectorsurvR` R package.

## Details

Derived from the vectorsurvR sample datasets in two steps.

- Filtering to female *Culex pipiens* pools with WNV testing within the
  SMA, using point-in-polygon against the union of the Sacramento,
  Placer and El Dorado counties.

- Estimating the size `est_pool_n` of each pool.

  - The total female counts from nearby collection points within the
    *same week* and within a spatial radius (for example 2 km) are
    summed to obtain \\T\_{\mathrm{near}}\\.

  - The nearby pools within the same week and radius are counted to
    obtain \\m\_{\mathrm{near}}\\.

  - `est_pool_n` is `round(max(1, T_near / m_near))`, or a conservative
    fallback of 25 if no nearby collections are found.

The coordinate reference system is 4326. `est_pool_n` is an estimate for
demonstration and pooled-modelling examples, not the recorded laboratory
pool size.
