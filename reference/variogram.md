# Empirical variogram with a global permutation envelope

Computes the empirical semivariogram using lag-distance breakpoints and,
optionally, an extreme-rank-length global envelope for spatial
independence.

## Usage

``` r
variogram(
  data,
  variable,
  breaks = NULL,
  n_bins = 14L,
  max_dist = NULL,
  n_permutations = 999L,
  level = 0.95,
  seed = NULL,
  distance_crs = NULL,
  distance_units = c("km", "m")
)
```

## Arguments

- data:

  an object of class `sf` containing the variable for which the
  variogram is to be computed and the coordinates

- variable:

  a character indicating the name of variable for which the variogram is
  to be computed.

- breaks:

  an optional numeric vector of lag-distance breakpoints. If supplied,
  these breakpoints are used directly to define the distance classes. A
  warning is raised when some classes lie entirely beyond the largest
  distance between locations, since they contain no pairs.

- n_bins:

  the number of lag-distance classes to generate when `breaks = NULL`.
  By default `n_bins = 14`.

- max_dist:

  an optional maximum lag distance. When `breaks = NULL`, the
  breakpoints are generated as
  `seq(0, max_dist, length.out = n_bins + 1)`. By default
  `max_dist = NULL` and the upper lag distance is set to `d_max/3`,
  where `d_max` is the maximum observed distance in the data.

- n_permutations:

  a non-negative integer indicating the number of random permutations
  used with the observed curve to construct the global envelope. By
  default `n_permutations = 999`; set it to zero to omit the envelope.

- level:

  the simultaneous coverage probability of the global envelope. By
  default `level = 0.95`.

- seed:

  an optional non-negative integer seed for the permutations. The
  caller's random-number state is restored before the function returns.

- distance_crs:

  `NULL` to retain an existing projected CRS, or automatically reproject
  longitude/latitude data to the CRS from
  [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
  (with a message naming the CRS and the units of distance).
  Alternatively, a CRS to reproject `data` to, which must itself be
  projected. A warning of class `riskmap_warning_crs` is raised when the
  CRS used distorts distances by more than 1% over the data.

- distance_units:

  Character string, either `"km"` or `"m"`, indicating whether the
  distances used in the variogram are expressed in kilometers or meters.
  By default `distance_units = "km"`, as in
  [`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

## Value

an object of class `RiskMap_variogram` which is a list containing the
following components:

- variogram:

  a data-frame containing the following columns:

  distance

  :   the mean pair distance within each lag class

  semivariance

  :   the observed empirical semivariance

  n_pairs

  :   the number of pairs in the lag class

  If `n_permutations > 0`, the data-frame also contains the following
  columns:

  lower_envelope

  :   the lower bound of the simultaneous envelope

  upper_envelope

  :   the upper bound of the simultaneous envelope

- distance_units:

  the value passed to `distance_units`

- n_permutations:

  the number of permutations

- breaks:

  the calculated breaks

- level:

  the simultaneous coverage probability

- envelope_method:

  `"global_extreme_rank_length"`, or `NULL` if no permutations were
  requested

## Details

The observed variogram is included among the permuted curves so that
they are exchangeable under the null hypothesis of spatial independence.
Curves are ordered by their extreme-rank-length ordering. Unlike
separate pointwise intervals, the resulting envelope controls the
probability that the variogram leaves the envelope anywhere across the
lag distances. The envelope is intended as an exploratory diagnostic
rather than a formal hypothesis test.

## Examples

``` r
data(italy_sim)

italy_variogram <- variogram(
                     data = italy_sim[1:200,],
                     variable = "y",
                     n_bins = 10,
                     n_permutations = 199,
                     seed = 123)

plot(italy_variogram)

```
