# Create test sets for cross-validation

Splits the locations of a dataset into training and test sets for
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md).
Each test set is predicted from the rest of the data, its training set.
Records at the same location are always in the same set.

## Usage

``` r
create_test_sets(
  data,
  method = c("leave_one_out", "random", "cluster", "spaced"),
  n_sets = NULL,
  size = NULL,
  min_dist = NULL,
  repeats = 1,
  seed = NULL
)
```

## Arguments

- data:

  An `sf` object with point geometries, the data the models are fitted
  to, or a fit from
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).

- method:

  How the test sets are made.

  - `"leave_one_out"` (the default) has one location in each test set.

  - `"random"` divides the locations at random into `n_sets` sets of
    about equal size.

  - `"cluster"` divides the locations into `n_sets` spatial clusters, by
    k-means on their coordinates.

  - `"spaced"` draws `n_sets` test sets of `size` locations that are at
    least `min_dist` apart.

- n_sets:

  The number of test sets, 10 by default. With `"leave_one_out"` it is
  the number of locations, and it can't be given.

- size:

  For `"spaced"` only, the number of locations in each test set.

- min_dist:

  For `"spaced"` only, the smallest distance between two locations of a
  test set, in the units of distance of the model, which are kilometres
  for an `sf` object and `distance_units` for a fit.

- repeats:

  The number of times that a `"random"` or `"cluster"` division is
  repeated, which gives `repeats` times `n_sets` test sets.

- seed:

  Optional non-negative integer that makes the test sets reproducible.
  The caller's random number state is restored on exit.

## Value

An object of class `RiskMap_test_sets`, a list with

- sets:

  A list with one element for each test set, each a list of the row
  numbers of the training data, `train`, and of the test data, `test`.

- repetition:

  The repeat each test set belongs to.

- method:

  The method, or `"user"` for test sets given to
  [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
  as a list.

- geometry:

  The locations of the data, in their own CRS.

- coordinates, distance_units:

  The projected coordinates of the rows of the data, in the units of
  distance.

- crs:

  The projected CRS.

- location:

  The distinct location of each row of the data.

- min_dist, size:

  As given.

## Details

Choose the test sets so that the test locations are about as far from
their training data as the locations to be predicted are from all the
data (Milà *et al.* 2022).
[`plot_test_sets()`](https://claudiofronterre.github.io/RiskMap/reference/plot_test_sets.md)
with `type = "distance"` compares these distances. Leave-one-out and
random test sets suit predictions between the locations of the data, and
clusters suit predictions in areas without data. The distance `min_dist`
of `"spaced"` is between the locations of a test set, which keeps them
spread over the region. It isn't a distance from the training data.

Distances and clusters are computed on projected coordinates, as in
[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md).
Data in longitude and latitude are projected to the CRS that
[`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
gives, and a fit has the CRS and the units of distance of its model.

The clusters of k-means are the best of 20 random starts. For
`"spaced"`, the locations are visited in random order and one is kept if
it is at least `min_dist` from those kept so far, until `size` locations
are kept or every location has been visited. A test set can therefore
have fewer than `size` locations, and a warning says so.

## Test sets made by another package

[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
also takes a list with one element for each test set, the row numbers of
the test observations or a list of the training rows, `train`, and the
test rows, `test`. Observations in neither are left out, for example
those in a buffer around the test locations. For the folds of the
rsample or the spatialsample package, this list is

    lapply(folds$splits, function(split) {
      list(train = as.integer(split),
           test = as.integer(split, data = "assessment"))
    })

and for `folds_list` of the blockCV package it is

    lapply(folds$folds_list, function(fold) {
      list(train = fold[[1]], test = fold[[2]])
    })

## References

Milà, C., Mateu, J., Pebesma, E. & Meyer, H. (2022). Nearest neighbour
distance matching Leave-One-Out Cross-Validation for map validation.
*Methods in Ecology and Evolution*, 13(6), 1304-1316.
[doi:10.1111/2041-210X.13851](https://doi.org/10.1111/2041-210X.13851) .

## See also

[`plot_test_sets()`](https://claudiofronterre.github.io/RiskMap/reference/plot_test_sets.md),
which [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of the
test sets also calls, and
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md),
which scores the predictions.

## Examples

``` r
data(italy_sim)
example_data <- italy_sim[1:100, ]

# One location at a time
create_test_sets(example_data)
#> ── Test sets for cross-validation ──────────────────────────────────────────────
#> 10 test sets of 10 observations each from 100 observations at 10 locations.
#> Leave-one-out, with one location in each test set.
#> The distance from a test location to its nearest training location is 51.4 km
#> in the median and at most 397.9 km.
#> 
#> ℹ `plot()` maps the test sets, and `assess_prediction()` scores the predictions
#>   of models on them.

# Five spatial clusters
clusters <- create_test_sets(example_data, method = "cluster", n_sets = 5,
                             seed = 1)
clusters
#> ── Test sets for cross-validation ──────────────────────────────────────────────
#> 5 test sets of 10 to 40 observations from 100 observations at 10 locations.
#> Spatial clusters of locations, by k-means.
#> The distance from a test location to its nearest training location is 153.6 km
#> in the median and at most 397.9 km.
#> 
#> ℹ `plot()` maps the test sets, and `assess_prediction()` scores the predictions
#>   of models on them.
plot(clusters)


# Ten test sets of five locations at least 50 km apart
create_test_sets(example_data, method = "spaced", n_sets = 10, size = 5,
                 min_dist = 50, seed = 1)
#> ── Test sets for cross-validation ──────────────────────────────────────────────
#> 10 test sets of 50 observations each from 100 observations at 10 locations.
#> Locations at least 50 km apart.
#> The distance from a test location to its nearest training location is 120.2 km
#> in the median and at most 739.1 km.
#> 
#> ℹ `plot()` maps the test sets, and `assess_prediction()` scores the predictions
#>   of models on them.
```
