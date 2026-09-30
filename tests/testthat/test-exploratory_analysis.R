test_that("summarise_distance produces errors", {

  test_that("data must be an object of class 'sf'", {
    df <- data.frame(x = 1:5, y = 1:5, y = rnorm(5))
    expect_error(
      summarise_distance(df),
      "'data' must be of class 'sf'"
    )
  })

  test_that("distance_crs must be a valid CRS", {
    expect_error(
      summarise_distance(gaussian_data, distance_crs = "a"),
      "The 'distance_crs' provided is not a valid CRS"
    )
  })

  test_that("distance_crs must be projected", {
    expect_error(
      summarise_distance(gaussian_data, distance_crs = 4326),
      "'distance_crs' must be a projected CRS, not longitude/latitude"
    )
  })

  test_that("distance_units must be 'km' or 'm'", {
    expect_error(
      summarise_distance(gaussian_data, distance_units = "a"),
      "'distance_units' must be either 'km' or 'm'"
    )
  })
})

test_that("summarise_distance produces correct output", {

  square <- data.frame(
    x = c(1, 2, 1, 2),
    y = c(1, 1, 2, 2)
  )

  square_sf <- sf::st_as_sf(square, coords = c("x", "y"), crs = 32630)

  # already projected: distance_crs = NULL retains the existing CRS as-is
  result <- summarise_distance(square_sf, distance_units = "m")
  expect_length(result, 4)
  expect_setequal(names(result), c("min", "max", "mean", "median"))
  expect_equal(result[["min"]], 1)
  expect_equal(result[["max"]], sqrt(2))
  expect_equal(result[["mean"]], (4 + (2*sqrt(2))) / 6)
  expect_equal(result[["median"]], 1)

  # reprojected to an explicit distance_crs (tolerance to account for round trip)
  result <- summarise_distance(square_sf, distance_crs = 3857, distance_units = "m")
  expect_equal(result[["min"]], 1, tolerance = 0.1)
  expect_equal(result[["max"]], sqrt(2), tolerance = 0.1)
  expect_equal(result[["mean"]], (4 + (2*sqrt(2))) / 6, tolerance = 0.1)
  expect_equal(result[["median"]], 1, tolerance = 0.1)

  # scaled (default distance_units = "km")
  result <- summarise_distance(square_sf)
  expect_equal(result[["min"]], 0.001)
  expect_equal(result[["max"]], sqrt(2) / 1000)
  expect_equal(result[["mean"]], (4 + (2*sqrt(2))) / 6000)
  expect_equal(result[["median"]], 0.001)

})

test_that("summarise_distance follows the glgpm() CRS convention", {

  latlon <- st_transform(gaussian_data, 4326)
  suggested_crs <- propose_utm(latlon)

  expect_message(
    result <- summarise_distance(latlon, distance_units = "m"),
    "automatically reprojecting to EPSG"
  )
  expect_equal(result, summarise_distance(st_transform(latlon, suggested_crs),
                                          distance_units = "m"))

  # an already-projected CRS is retained as-is, with no message
  expect_no_message(summarise_distance(gaussian_data))

  # an explicit distance_crs is honoured
  result <- summarise_distance(gaussian_data, distance_crs = 32637, distance_units = "m")
  expect_equal(result, summarise_distance(gaussian_data, distance_units = "m"))
})

test_that("summarise_distance derives distance units from a non-metre CRS", {

  foot_data <- st_transform(gaussian_data, 2263)
  coords <- coordinates_in_units(foot_data, "m")

  result <- summarise_distance(foot_data, distance_units = "m")
  expect_equal(result[["max"]], max(dist(unique(coords))))
})


test_that("variogram produces errors", {

  test_that("data must be an object of class 'sf'", {
    df <- data.frame(x = 1:5, y = 1:5, y = rnorm(5))
    expect_error(
      variogram(df, variable = "y"),
      "'data' must be of class 'sf'"
    )
  })

  test_that("variable must be a character vector", {
    expect_error(
      variogram(gaussian_data, variable = 1),
      "'variable' must be a single object of class 'character'"
    )
  })

  test_that("variable must have length 1", {
    expect_error(
      variogram(gaussian_data, variable = c("y", "y")),
      "'variable' must be a single object of class 'character'"
    )
  })

  test_that("variable must exist in the columns of data", {
    expect_error(
      variogram(gaussian_data, variable = "not_a_column"),
      "'variable' must be one of the columns in 'data'"
    )
  })

  test_that("breaks and n_bins cannot both be explicitly supplied", {
    expect_error(
      variogram(gaussian_data, variable = "y",
                breaks = c(0, 1, 2), n_bins = 5),
      "'breaks' and 'n_bins' cannot both be supplied"
    )
  })

  test_that("supplying breaks alone does not trigger the n_bins conflict", {
    expect_no_error(
      variogram(gaussian_data,
                variable = "y",
                breaks = seq(0, 9000, 1000)),
    )
  })

  test_that("breaks and max_dist cannot both be explicitly supplied", {
    expect_error(
      variogram(gaussian_data, variable = "y",
                breaks = c(0, 1, 2), max_dist = 5),
      "'breaks' and 'max_dist' cannot both be supplied"
    )
  })

  test_that("n_bins must be numeric", {
    expect_error(
      variogram(gaussian_data, variable = "y", n_bins = "10"),
      "'n_bins' must be a single positive integer"
    )
  })

  test_that("n_bins must have length 1", {
    expect_error(
      variogram(gaussian_data, variable = "y", n_bins = c(5, 10)),
      "'n_bins' must be a single positive integer"
    )
  })

  test_that("n_bins must be >= 1", {
    expect_error(
      variogram(gaussian_data, variable = "y", n_bins = 0),
      "'n_bins' must be a single positive integer"
    )
  })

  test_that("n_bins must be a whole number", {
    expect_error(
      variogram(gaussian_data, variable = "y", n_bins = 4.5),
      "'n_bins' must be a single positive integer"
    )
  })

  test_that("max_dist must be numeric", {
    expect_error(
      variogram(gaussian_data, variable = "y", max_dist = "100"),
      "'max_dist' must be a positive numeric value"
    )
  })

  test_that("max_dist must have length 1", {
    expect_error(
      variogram(gaussian_data, variable = "y", max_dist = c(10, 20)),
      "'max_dist' must be a positive numeric value"
    )
  })

  test_that("max_dist must be strictly positive", {
    expect_error(
      variogram(gaussian_data, variable = "y", max_dist = 0),
      "'max_dist' must be a positive numeric value"
    )
    expect_error(
      variogram(gaussian_data, variable = "y", max_dist = -5),
      "'max_dist' must be a positive numeric value"
    )
  })

  test_that("n_permutations must be non-negative", {
    expect_error(
      variogram(gaussian_data, variable = "y", n_permutations = -1),
      "'n_permutations' must be a single non-negative integer"
    )
  })

  test_that("n_permutations must be a whole number", {
    expect_error(
      variogram(gaussian_data, variable = "y", n_permutations = 10.5),
      "'n_permutations' must be a single non-negative integer"
    )
  })

  test_that("too few permutations for the requested level are rejected", {
    expect_error(
      variogram(gaussian_data, variable = "y", n_permutations = 2),
      "'n_permutations' is too small to construct a global envelope"
    )
  })

  test_that("the requested level controls the minimum permutation count", {
    expect_no_warning(
      variogram(gaussian_data, variable = "y", n_permutations = 50,
                level = 0.95)
    )
    expect_error(
      variogram(gaussian_data, variable = "y", n_permutations = 50,
                level = 0.99),
      "too small to construct a global envelope"
    )
  })

  test_that("n_permutations of 0 is allowed and does not warn about being low", {
    expect_no_warning(
        variogram(gaussian_data, variable = "y", n_permutations = 0)
    )
  })

  test_that("n_permutations of exactly 100 does not trigger the low-permutation warning", {
    expect_no_warning(
        variogram(gaussian_data, variable = "y", n_permutations = 100)
    )
  })

  test_that("an explicit distance_crs that is longitude/latitude produces an error", {
    expect_error(
        variogram(latlon_data, variable = "y", n_permutations = 100, distance_crs = 4326),
      "'distance_crs' must be a projected CRS, not longitude/latitude"
    )
  })

  test_that("longitude/latitude data is automatically reprojected with a message", {
    expect_message(
        variogram(latlon_data, variable = "y", n_permutations = 100),
      "automatically reprojecting to EPSG"
    )
  })

  test_that("distance_crs = NULL (default) on already-projected data does not emit that message", {
    expect_no_message(
        variogram(gaussian_data, variable = "y", n_permutations = 100)
    )
  })

  test_that("breaks must be numeric", {
    expect_error(
      variogram(gaussian_data, variable = "y", breaks = c("a", "b", "c")),
      "'breaks' must be a numeric vector with at least two values"
    )
  })

  test_that("breaks must have at least two values", {
    expect_error(
      variogram(gaussian_data, variable = "y", breaks = 5),
      "'breaks' must be a numeric vector with at least two values"
    )
  })

  test_that("breaks must be strictly increasing", {
    expect_error(
      variogram(gaussian_data, variable = "y", breaks = c(0, 5, 5, 10)),
      "'breaks' must be strictly increasing"
    )
    expect_error(
      variogram(gaussian_data, variable = "y", breaks = c(0, 10, 5)),
      "'breaks' must be strictly increasing"
    )
  })

  test_that("breaks must be non-negative", {
    expect_error(
      variogram(gaussian_data, variable = "y", breaks = c(-5, 0, 5)),
      "'breaks' must be non-negative"
    )
  })

  test_that("breaks exceeding the maximum observed distance raise an error", {
    expect_error(
      variogram(gaussian_data, variable = "y", max_dist = 1e9),
      "the provided lag distances go beyond the maximum observed distance"
    )
  })

  test_that("distance_crs must be a valid, projected CRS and distance_units must be 'km' or 'm'", {
    expect_error(
      variogram(gaussian_data, variable = "y", distance_crs = "a"),
      "The 'distance_crs' provided is not a valid CRS")

    expect_error(
      variogram(gaussian_data, variable = "y", distance_units = "a"),
      "'distance_units' must be either 'km' or 'm'")
  })

})


test_that("variogram produces expected output", {

  expected_output <- c("variogram", "distance_units", "n_permutations", "breaks",
                       "level", "envelope_method")
  expected_columns <- c("distance", "semivariance", "n_pairs",
                        "lower_envelope", "upper_envelope")
  expected_columns_zero <- c("distance", "semivariance", "n_pairs")

  result <- variogram(gaussian_data, variable = "y", n_bins = 10,
                      n_permutations = 99, seed = 123)

  expect_s3_class(result, "RiskMap_variogram")
  expect_setequal(names(result), expected_output)
  expect_setequal(names(result$variogram), expected_columns)
  expect_equal(nrow(result$variogram), 10)
  expect_equal(result$distance_units, "m")
  expect_equal(result$n_permutations, 99)
  expect_length(result$breaks, 10 + 1)
  expect_equal(result$level, 0.95)
  expect_equal(result$envelope_method, "global_extreme_rank_length")

  breaks <- seq(0, 9000, 1000)
  result <- variogram(gaussian_data, variable = "y", breaks = breaks,
                      n_permutations = 99, seed = 123)
  expect_s3_class(result, "RiskMap_variogram")
  expect_setequal(names(result), expected_output)
  expect_setequal(names(result$variogram), expected_columns)
  expect_equal(result$distance_units, "m")
  expect_equal(result$n_permutations, 99)
  expect_equal(nrow(result$variogram), length(breaks) - 1)
  expect_equal(result$breaks, breaks)

  result <- variogram(gaussian_data, variable = "y", n_bins = 10, n_permutations = 0)
  expect_s3_class(result, "RiskMap_variogram")
  expect_setequal(names(result), expected_output)
  expect_setequal(names(result$variogram), expected_columns_zero)
  expect_equal(nrow(result$variogram), 10)
  expect_equal(result$distance_units, "m")
  expect_equal(result$n_permutations, 0)
  expect_length(result$breaks, 10 + 1)
  expect_null(result$envelope_method)
})

test_that("variogram values are correct for a simple dataset", {

  square <- data.frame(
    x = c(1, 2, 1, 2, 10),
    y = c(1, 1, 2, 2, 10),
    z = c( 1, 2, 3, 4, 10)
  )

  # 10 points will be removed, so only the square distances are considered
  square_sf <- sf::st_as_sf(square, coords = c("x", "y"), crs = 32630)

  manual_observed <- c(
    mean(c((1 - 2)^2, (1 - 3)^2, (2 - 4)^2, (3 - 4)^2) / 2),
    mean(c((1 - 4)^2, (2 - 3)^2) / 2)
  )

  result <- variogram(square_sf,
                      variable = "z",
                      breaks = seq(0.1, 1.5, 0.2),
                      n_permutations = 119,
                      seed = 12)

  expect_equal(sum(result$variogram$n_pairs), 6)

  bin1 <- result$variogram$n_pairs == 4
  bin2 <- result$variogram$n_pairs == 2

  expect_equal(result$variogram$distance[bin1], 1)
  expect_equal(result$variogram$distance[bin2], sqrt(2))

  expect_equal(result$variogram$semivariance[bin1], manual_observed[1])
  expect_equal(result$variogram$semivariance[bin2], manual_observed[2])
})

test_that("compiled binned semivariances agree with a direct calculation", {
  values <- c(2, 5, 11, 13)
  permutations <- cbind(1:4, c(4, 2, 1, 3))
  first <- c(1L, 2L, 3L, 2L, 3L, 3L)
  second <- c(0L, 0L, 0L, 1L, 1L, 2L)
  bins <- c(0L, 0L, 1L, 1L, 0L, 1L)

  permuted_values <- matrix(values[permutations], nrow = length(values))
  result <- cpp_binned_semivariances(permuted_values, first, second, bins, 2L)
  direct <- vapply(seq_len(ncol(permutations)), function(column) {
    permuted <- values[permutations[, column]]
    pair_values <- (permuted[first + 1L] - permuted[second + 1L])^2 / 2
    c(mean(pair_values[bins == 0L]), mean(pair_values[bins == 1L]))
  }, numeric(2))

  expect_equal(result, direct)
})

test_that("global rank envelope uses the whole curve and handles ties conservatively", {
  curves <- rbind(
    c(10, 0, 4, 4, 4),
    c(0, 10, 4, 4, 4)
  )

  envelope <- global_rank_envelope(curves, level = 0.8)

  expect_equal(envelope$lower, c(0, 0))
  expect_equal(envelope$upper, c(10, 10))
})

test_that("a variogram seed is reproducible without changing the caller RNG", {
  set.seed(91)
  state_before <- .Random.seed
  first <- variogram(gaussian_data, "y", n_permutations = 99, seed = 5)
  expect_identical(.Random.seed, state_before)
  second <- variogram(gaussian_data, "y", n_permutations = 99, seed = 5)

  expect_equal(first, second)
})

test_that("plot_variogram defaults to plotting the envelope when available", {

  result <- variogram(gaussian_data, variable = "y", n_bins = 10,
                      n_permutations = 99, seed = 123)

  expect_no_warning(p <- plot_variogram(result))
  expect_true("GeomRibbon" %in% vapply(p$layers, function(l) class(l$geom)[1], character(1)))
})

test_that("plot_variogram warns and skips the envelope when n_permutations <= 1", {

  result <- variogram(gaussian_data, variable = "y", n_bins = 10, n_permutations = 0)

  expect_warning(
    p <- plot_variogram(result),
    "'n_permutations' was 0"
  )
  expect_false("GeomRibbon" %in% vapply(p$layers, function(l) class(l$geom)[1], character(1)))

  expect_no_warning(plot_variogram(result, plot_envelope = FALSE))
})
