# Empirical variogram

Computes the empirical variogram using “bins” of distance provided by
the user.

## Usage

``` r
s_variogram(
  data,
  variable,
  bins = NULL,
  n_permutation = 0,
  convert_to_utm = TRUE,
  scale_to_km = FALSE
)
```

## Arguments

- data:

  an object of class `sf` containing the variable for which the
  variogram is to be computed and the coordinates

- variable:

  a character indicating the name of variable for which the variogram is
  to be computed.

- bins:

  a vector indicating the \`bins\` to be used to define the classes of
  distance used in the computation of the variogram. By default
  `bins=NULL` and bins are then computed as `seq(0, d_max/2, length=15)`
  where `d_max` is the maximum distance observed in the data.

- n_permutation:

  a non-negative integer indicating the number of permutation used to
  compute the 95 level envelope under the assumption of spatial
  independence. By default `n_permutation=0`, and no envelope is
  generated.

- convert_to_utm:

  a logical value, indicating if the conversion to UTM shuold be
  performed (`convert_to_utm = TRUE`) or the coordinate reference system
  of the data must be used without any conversion
  (`convert_to_utm = FALSE`). By default `convert_to_utm = TRUE`. Note:
  if `convert_to_utm = TRUE` the conversion to UTM is performed using
  the epsg provided by
  [`propose_utm`](https://claudiofronterre.github.io/RiskMap/reference/propose_utm.md).

- scale_to_km:

  a logical value, indicating if the distances used in the variogram
  must be scaled to kilometers (`scale_to_km = TRUE`) or left in meters
  (`scale_to_km = FALSE`). By default `scale_to_km = FALSE`

## Value

an object of class 'variogram' which is a list containing the following
components

`variogram` a data-frame containing the following columns: `mid_points`,
the middle points of the classes of distance provided by `bins`;
`obs_vari` the values of the observed variogram; `obs_vari` the number
of pairs. If `n_permutation > 0`, the data-frame also contains
`lower_bound` and `upper_bound` corresponding to the lower and upper
bounds of the 95 used to assess the departure of the observed variogram
from the assumption of spatial independence.

`scale_to_km` the value passed to `scale_to_km`

`n_permutation` the number of permutations

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
