# Test sets for cross-validation, create_test_sets() and those given by rows.

# Thirty locations, the last ten of them with two records.
test_set_data <- function(crs = 32637) {
  set.seed(30)
  coordinates <- cbind(runif(30, 0, 1e5), runif(30, 0, 1e5))
  rows <- c(1:30, 21:30)
  # A field of scale 30 km and a measurement error, without changing the
  # session's random numbers after the locations.
  value <- withr::with_seed(31, {
    correlation <- exp(-as.matrix(dist(coordinates)) / 30000)
    field <- drop(t(chol(correlation)) %*% rnorm(30))
    field[rows] + rnorm(length(rows), sd = sqrt(0.3))
  })
  data <- data.frame(x = coordinates[rows, 1], y = coordinates[rows, 2],
                     value = value)
  st_transform(st_as_sf(data, coords = c("x", "y"), crs = 32637), crs)
}

# Every test set has its own rows of test data, all the others for training,
# and no location in both.
expect_complementary <- function(test_sets) {
  n <- length(test_sets$location)
  for (set in test_sets$sets) {
    expect_identical(sort(c(set$train, set$test)), seq_len(n))
    expect_false(any(test_sets$location[set$test] %in%
                       test_sets$location[set$train]))
    expect_gt(length(set$train), 0)
    expect_gt(length(set$test), 0)
  }
}

test_that("create_test_sets() divides the locations", {
  data <- test_set_data()
  loo <- create_test_sets(data)
  expect_s3_class(loo, "RiskMap_test_sets")
  expect_complementary(loo)
  expect_length(loo$sets, 30L)
  expect_identical(lengths(lapply(loo$sets, `[[`, "test")),
                   rep(1:2, c(20, 10)))
  for (method in c("random", "cluster")) {
    sets <- create_test_sets(data, method = method, n_sets = 4, seed = 1)
    expect_complementary(sets)
    expect_length(sets$sets, 4L)
    expect_identical(sort(unlist(lapply(sets$sets, `[[`, "test"))), 1:40)
    expect_length(create_test_sets(data, method = method)$sets, 10L)
  }
  # The numbers of locations of random sets differ by one at most.
  sets <- create_test_sets(data, "random", n_sets = 4, seed = 1)
  n_locations <- vapply(sets$sets, function(set) {
    length(unique(sets$location[set$test]))
  }, integer(1))
  expect_identical(sort(n_locations), c(7L, 7L, 8L, 8L))
})

test_that("a division is repeated, each time another one", {
  data <- test_set_data()
  sets <- create_test_sets(data, "random", n_sets = 3, repeats = 2, seed = 1)
  expect_complementary(sets)
  expect_length(sets$sets, 6L)
  expect_identical(sets$repetition, rep(1:2, each = 3L))
  tested <- lapply(sets$sets, `[[`, "test")
  expect_identical(sort(unlist(tested[1:3])), 1:40)
  expect_identical(sort(unlist(tested[4:6])), 1:40)
  expect_false(identical(tested[1:3], tested[4:6]))
})

test_that("clusters are compact, numbered from north to south", {
  # Three groups of locations far apart.
  set.seed(31)
  centres <- cbind(c(0, 5e4, 0), c(0, 5e4, 1e5))
  group <- rep(1:3, each = 5)
  coordinates <- centres[group, ] + matrix(rnorm(30, 0, 500), 15)
  data <- st_as_sf(data.frame(x = coordinates[, 1], y = coordinates[, 2]),
                   coords = c("x", "y"), crs = 32637)
  for (seed in 1:3) {
    sets <- create_test_sets(data, "cluster", n_sets = 3, seed = seed)
    expect_identical(lapply(sets$sets, `[[`, "test"),
                     list(11:15, 6:10, 1:5))
  }
})

test_that("test sets have one location, or as many sets as locations", {
  data <- test_set_data()
  for (method in c("random", "cluster")) {
    sets <- create_test_sets(data, method, n_sets = 30, seed = 1)
    expect_complementary(sets)
    expect_identical(sort(vapply(sets$sets, function(set) {
      unique(sets$location[set$test])
    }, integer(1))), 1:30)
  }
  sets <- create_test_sets(data, "spaced", n_sets = 2, size = 1, min_dist = 1,
                           seed = 1)
  expect_complementary(sets)
  expect_length(sets$sets, 2L)
  two <- data[1:2, ]
  expect_complementary(create_test_sets(two, "cluster", n_sets = 2))
  expect_complementary(create_test_sets(two))
})

test_that("the locations of spaced test sets are at least min_dist apart", {
  data <- test_set_data()
  sets <- create_test_sets(data, "spaced", n_sets = 5, size = 6, min_dist = 15,
                           seed = 1)
  expect_complementary(sets)
  for (set in sets$sets) {
    tested <- set$test[!duplicated(sets$location[set$test])]
    expect_length(tested, 6L)
    # In km.
    expect_gte(min(dist(st_coordinates(data)[tested, ])), 15e3)
  }
  # All the locations that fit, and not fewer.
  expect_warning(
    short <- create_test_sets(data, "spaced", n_sets = 2, size = 25,
                              min_dist = 15, seed = 1),
    "2 of 2 test sets have fewer than `size` = 25 locations"
  )
  coordinates <- st_coordinates(data)
  for (set in short$sets) {
    to_tested <- vapply(set$train, function(row) {
      min(sqrt(rowSums((coordinates[set$test, , drop = FALSE] -
                          rep(coordinates[row, ],
                              each = length(set$test)))^2)))
    }, numeric(1))
    expect_lt(max(to_tested), 15e3)
  }
  # A minimum distance beyond the region leaves one location.
  expect_warning(
    one <- create_test_sets(data, "spaced", n_sets = 1, size = 2,
                            min_dist = 1e4, seed = 1),
    "1 of 1 test set has fewer"
  )
  expect_length(unique(one$location[one$sets[[1]]$test]), 1L)
})

test_that("min_dist is in the units of distance of a fit", {
  expect_identical(gaussian_intercept_model$distance_units, "km")
  in_metres <- glgpm(y ~ gp(), data = gaussian_intercept_data, family = "gaussian",
                     distance_units = "m", messages = FALSE)
  metres <- create_test_sets(in_metres, "spaced", n_sets = 3, size = 3,
                             min_dist = 3000, seed = 1)
  kilometres <- create_test_sets(gaussian_intercept_model, "spaced",
                                 n_sets = 3, size = 3, min_dist = 3, seed = 1)
  expect_identical(metres$sets, kilometres$sets)
  expect_identical(metres$distance_units, "m")
  expect_equal(metres$coordinates, 1000 * kilometres$coordinates,
               ignore_attr = TRUE)
  # An sf object has kilometres, whatever the units of its CRS.
  in_km <- st_transform(gaussian_data, sub("units=m", "units=km",
                                           st_crs(gaussian_data)$proj4string))
  from_km <- create_test_sets(in_km, "spaced", n_sets = 3, size = 3,
                              min_dist = 3, seed = 1)
  expect_identical(from_km$sets, kilometres$sets)
})

test_that("data in longitude and latitude are projected as a fit's are", {
  data <- test_set_data()
  lonlat <- test_set_data(4326)
  expect_message(
    sets <- create_test_sets(lonlat, "cluster", n_sets = 4, seed = 1),
    "in longitude and latitude"
  )
  expect_complementary(sets)
  expect_false(isTRUE(st_is_longlat(sets$crs)))
  expect_identical(st_crs(sets$geometry), st_crs(lonlat))
  fit <- suppressMessages(glgpm(value ~ gp(), data = lonlat,
                                family = "gaussian", messages = FALSE))
  expect_identical(st_crs(fit$data), sets$crs)
  expect_equal(sets$coordinates, coordinates_in_units(fit$data, "km"),
               ignore_attr = TRUE)
  # The distances don't matter to the other methods, which say nothing.
  expect_no_message(create_test_sets(lonlat, "random", seed = 1))
  expect_no_message(loo <- create_test_sets(lonlat))
  # A projection doesn't merge or split locations.
  expect_length(loo$sets, 30L)
  # The test sets of the data serve a fit to them.
  result <- assess_prediction(fit, sets, seed = 1, messages = FALSE)
  expect_identical(nrow(result$scores), 40L)
  expect_identical(nrow(result$failures), 0L)
})

test_that("seed reproduces test sets and leaves the session's random numbers", {
  data <- test_set_data()
  set.seed(1)
  state <- .Random.seed
  first <- create_test_sets(data, "cluster", seed = 5)
  expect_identical(.Random.seed, state)
  expect_identical(create_test_sets(data, "cluster", seed = 5), first)
  # Without a seed they come from the session.
  one <- create_test_sets(data, "random")
  expect_false(identical(.Random.seed, state))
  set.seed(1)
  expect_identical(create_test_sets(data, "random"), one)
})

test_that("create_test_sets() checks its arguments", {
  data <- test_set_data()
  expect_input_error(create_test_sets(), "data")
  expect_input_error(create_test_sets(data.frame(x = 1)), "data", "sf")
  expect_input_error(create_test_sets(st_buffer(data, 10)), "data")
  expect_input_error(create_test_sets(data[c(21, 31), ]), "data",
                     "at least two locations")
  expect_error(create_test_sets(data, "regularized"),
               class = "riskmap_error_input")
  expect_input_error(create_test_sets(data, n_sets = 3), "n_sets",
                     "one test set for each location")
  for (method in c("random", "cluster")) {
    expect_input_error(create_test_sets(data, method, n_sets = 1), "n_sets",
                       "at least 2")
    expect_input_error(create_test_sets(data, method, n_sets = 31), "n_sets",
                       "number of locations of `data`, 30")
    expect_input_error(create_test_sets(data, method, n_sets = 2.5), "n_sets")
    expect_input_error(create_test_sets(data, method, size = 2), "size",
                       "only when `method` is")
    expect_input_error(create_test_sets(data, method, min_dist = 2),
                       "min_dist", "only when `method` is")
  }
  # Fewer locations than the ten sets of the default.
  expect_input_error(create_test_sets(data[1:5, ], "cluster"), "n_sets",
                     "number of locations of `data`, 5")
  spaced <- function(...) create_test_sets(data, "spaced", ...)
  expect_input_error(spaced(size = 2), "min_dist", "must be given")
  expect_input_error(spaced(min_dist = 2), "size", "must be given")
  expect_input_error(spaced(size = 30, min_dist = 2), "size",
                     "below the number of locations of `data`, 30")
  expect_input_error(spaced(size = 2, min_dist = -1), "min_dist")
  expect_input_error(spaced(size = 2, min_dist = NA), "min_dist")
  expect_input_error(spaced(size = 1.5, min_dist = 1), "size")
  expect_input_error(spaced(size = 2, min_dist = 1, repeats = 2), "repeats",
                     "Give the number of test sets in `n_sets`")
  expect_input_error(create_test_sets(data, repeats = 2), "repeats",
                     "the same every time")
  expect_input_error(create_test_sets(data, "random", seed = -1), "seed")
})

test_that("print() describes test sets", {
  data <- test_set_data()
  printed <- squish(capture.output(
    result <- print(create_test_sets(data, "cluster", n_sets = 4, repeats = 2,
                                     seed = 1))
  ))
  expect_s3_class(result, "RiskMap_test_sets")
  expect_match(printed, "8 test sets of [0-9]+ to [0-9]+ observations from 40 observations at 30 locations\\.")
  expect_match(printed, "Spatial clusters of locations, by k-means. The division is repeated 2 times.", fixed = TRUE)
  expect_match(printed, "nearest training location is [0-9.]+ km in the median and at most [0-9.]+ km")
  expect_output(print(create_test_sets(data)),
                "30 test sets of 1 to 2 observations")
  expect_output(print(create_test_sets(data, "spaced", n_sets = 1, size = 1,
                                       min_dist = 2.5, seed = 1)),
                "Locations at\\s+least 2.5 km apart")
  expect_output(print(create_test_sets(data, "random", n_sets = 2, seed = 1)),
                "2 test sets of [0-9]+ and [0-9]+ observations")
  user <- read_test_sets(list(1:3, 4:6), gaussian_intercept_model)
  expect_output(print(user), paste("2 test sets of 3 observations each from", fixture_n))
  expect_input_error(print(user, digits = 2), "digits")
})

test_that("the distances of test locations are those to their training data", {
  data <- test_set_data()
  sets <- create_test_sets(data, "cluster", n_sets = 4, seed = 1)
  distances <- test_set_distances(sets)
  expect_identical(nrow(distances), 30L)
  all <- as.matrix(dist(st_coordinates(data))) / 1000
  direct <- vapply(seq_len(nrow(distances)), function(i) {
    min(all[distances$row[i], sets$sets[[distances$test_set[i]]]$train])
  }, numeric(1))
  expect_equal(distances$distance, direct)
  # In blocks, as for a large grid.
  from <- cbind(runif(23), runif(23))
  to <- cbind(runif(7), runif(7))
  expect_equal(nearest_distances(from, to, block_size = 5L),
               nearest_distances(from, to))
  expect_equal(nearest_distances(from, to[1, , drop = FALSE]),
               sqrt(rowSums((from - rep(to[1, ], each = 23))^2)))
})

test_that("test sets can be given as rows of the data", {
  read <- function(x) read_test_sets(x, gaussian_intercept_model)
  sets <- read(list(1:3, c(fixture_n, fixture_n - 1L)))
  expect_s3_class(sets, "RiskMap_test_sets")
  expect_identical(sets$method, "user")
  expect_identical(sets$sets, list(list(train = 4:fixture_n, test = 1:3),
                                   list(train = 1:(fixture_n - 2L), test = c(fixture_n, fixture_n - 1L))))
  # A training set of its own leaves out the other rows, as for a buffer.
  expect_identical(read(list(list(train = 1:5, test = (fixture_n - 1L):fixture_n)))$sets,
                   list(list(train = 1:5, test = (fixture_n - 1L):fixture_n)))
  expect_identical(read(list(list(test = (fixture_n - 1L):fixture_n)))$sets,
                   list(list(train = 1:(fixture_n - 2L), test = (fixture_n - 1L):fixture_n)))

  refused <- function(x, ...) expect_input_error(read(x), "test_sets", ...)
  refused("cluster", "create_test_sets(data, method = \"cluster\")")
  refused(1:3, "list(rows)")
  refused(matrix(c(1, rep(0, 9)), ncol = 1), "apply(x == 1, 2, which")
  refused(data.frame(a = 1:3))
  refused(list(), "at least one test set")
  refused(list(train = 1:5, test = 6:10), "list(list(train = , test = ))")
  refused(list(list(in_id = 1:5, out_id = 6:10)), "now called `train` and `test`")
  refused(list(list(train = 1:5)), "test rows")
  refused(list("a"))
  refused(list(seq_len(fixture_n)), "at least one observation for training")
  refused(list(integer(0)), "at least one row number")
  between <- paste("between 1 and", fixture_n)
  refused(list(1:3, 0),
          paste("test rows of test_sets[[2]] must be row numbers", between))
  refused(list(fixture_n + 1), between)
  refused(list(1.5), between)
  refused(list(c(1, NA)), "without missing values")
  refused(list(c(2, 2)), "without repeats")
  refused(list(list(train = 1:5, test = 5:6)), "must not overlap")
  refused(list(list(train = fixture_n + 1, test = 5:6)), "training rows")
  refused(list(list(train = integer(0), test = 5:6)), "training rows")
})

test_that("test sets must belong to the data of the fits", {
  sets <- create_test_sets(gaussian_data, "random", n_sets = 2, seed = 1)
  read <- function(x) read_test_sets(x, gaussian_intercept_model)
  expect_identical(read(sets), sets)
  # The same locations in another CRS are the same data.
  elsewhere <- create_test_sets(latlon_data, "random", n_sets = 2, seed = 1)
  expect_identical(read(elsewhere), elsewhere)
  expect_input_error(
    read(create_test_sets(gaussian_data[1:9, ], "random", n_sets = 2)),
    "test_sets", sprintf("%d rows, not 9", fixture_n)
  )
  expect_input_error(
    read(create_test_sets(gaussian_data[fixture_n:1, ], "random", n_sets = 2)),
    "test_sets", "same rows in the same order"
  )
})
