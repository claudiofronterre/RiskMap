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
      variogram(sf_data, variable = 1),
      "'variable' must be a single object of class 'character'"
    )
  })

  test_that("variable must have length 1", {
    expect_error(
      variogram(sf_data, variable = c("y", "y")),
      "'variable' must be a single object of class 'character'"
    )
  })

  test_that("variable must exist in the columns of data", {
    expect_error(
      variogram(sf_data, variable = "not_a_column"),
      "'variable' must be one of the columns in 'data'"
    )
  })

  test_that("breaks and n_bins cannot both be explicitly supplied", {
    expect_error(
      variogram(sf_data, variable = "y",
                breaks = c(0, 1, 2), n_bins = 5),
      "'breaks' and 'n_bins' cannot both be supplied"
    )
  })

  test_that("supplying breaks alone does not trigger the n_bins conflict", {
    expect_no_error(
      variogram(sf_data, variable = "y", breaks = seq(0, 100, 10), scale_to_km = TRUE),
    )
  })

  test_that("breaks and max_dist cannot both be explicitly supplied", {
    expect_error(
      variogram(sf_data, variable = "y",
                breaks = c(0, 1, 2), max_dist = 5),
      "'breaks' and 'max_dist' cannot both be supplied"
    )
  })

  test_that("n_bins must be numeric", {
    expect_error(
      variogram(sf_data, variable = "y", n_bins = "10"),
      "'n_bins' must be a positive integer"
    )
  })

  test_that("n_bins must have length 1", {
    expect_error(
      variogram(sf_data, variable = "y", n_bins = c(5, 10)),
      "'n_bins' must be a positive integer"
    )
  })

  test_that("n_bins must be >= 1", {
    expect_error(
      variogram(sf_data, variable = "y", n_bins = 0),
      "'n_bins' must be a positive integer"
    )
  })

  test_that("n_bins must be a whole number", {
    expect_error(
      variogram(sf_data, variable = "y", n_bins = 4.5),
      "'n_bins' must be a positive integer"
    )
  })

  test_that("max_dist must be numeric", {
    expect_error(
      variogram(sf_data, variable = "y", max_dist = "100"),
      "'max_dist' must be a positive numeric value"
    )
  })

  test_that("max_dist must have length 1", {
    expect_error(
      variogram(sf_data, variable = "y", max_dist = c(10, 20)),
      "'max_dist' must be a positive numeric value"
    )
  })

  test_that("max_dist must be strictly positive", {
    expect_error(
      variogram(sf_data, variable = "y", max_dist = 0),
      "'max_dist' must be a positive numeric value"
    )
    expect_error(
      variogram(sf_data, variable = "y", max_dist = -5),
      "'max_dist' must be a positive numeric value"
    )
  })

  test_that("n_permutations must be non-negative", {
    expect_error(
      variogram(sf_data, variable = "y", n_permutations = -1),
      "'n_permutations' must be a positive integer number"
    )
  })

  test_that("n_permutations must be a whole number", {
    expect_error(
      variogram(sf_data, variable = "y", n_permutations = 10.5),
      "'n_permutations' must be a positive integer number"
    )
  })

  test_that("n_permutations of exactly 2 is explicitly rejected", {
    expect_error(
      variogram(sf_data, variable = "y", n_permutations = 2),
      "'n_permutations' must be greater than 2"
    )
  })

  test_that("n_permutations below 100 raises a warning (but is not an error)", {
    expect_warning(
        variogram(sf_data, variable = "y", n_permutations = 50),
      "'n_permutations' is set very low - consider increasing it"
    )
  })

  test_that("n_permutations of 0 is allowed and does not warn about being low", {
    expect_no_warning(
        variogram(sf_data, variable = "y", n_permutations = 0)
    )
  })

  test_that("n_permutations of exactly 100 does not trigger the low-permutation warning", {
    expect_no_warning(
        variogram(sf_data, variable = "y", n_permutations = 100, convert_to_utm = TRUE)
    )
  })

  test_that("convert_to_utm = FALSE produces an error if coords are lon lat", {
    expect_error(
        variogram(sf_data, variable = "y", n_permutations = 100, convert_to_utm = FALSE),
      "The dataset coordinates are in longitude and latitude"
    )
  })

  test_that("convert_to_utm = TRUE (default) does not emit that message", {
    expect_no_message(
        variogram(sf_data, variable = "y", n_permutations = 100)
    )
  })

  test_that("breaks must be numeric", {
    expect_error(
      variogram(sf_data, variable = "y", breaks = c("a", "b", "c")),
      "'breaks' must be a numeric vector with at least two values"
    )
  })

  test_that("breaks must have at least two values", {
    expect_error(
      variogram(sf_data, variable = "y", breaks = 5),
      "'breaks' must be a numeric vector with at least two values"
    )
  })

  test_that("breaks must be strictly increasing", {
    expect_error(
      variogram(sf_data, variable = "y", breaks = c(0, 5, 5, 10)),
      "'breaks' must be strictly increasing"
    )
    expect_error(
      variogram(sf_data, variable = "y", breaks = c(0, 10, 5)),
      "'breaks' must be strictly increasing"
    )
  })

  test_that("breaks must be non-negative", {
    expect_error(
      variogram(sf_data, variable = "y", breaks = c(-5, 0, 5)),
      "'breaks' must be non-negative"
    )
  })

  test_that("breaks exceeding the maximum observed distance raise an error", {
    expect_error(
      variogram(sf_data, variable = "y", max_dist = 1e9),
      "the provided lag distances go beyond the maximum observed distance"
    )
  })

  test_that("convert_to_utm and scale_to_km must be logical", {
    expect_error(
      variogram(sf_data, variable = "y", convert_to_utm = "a"),
      "'convert_to_utm' must be either TRUE or FALSE")

    expect_error(
      variogram(sf_data, variable = "y", scale_to_km = "a"),
      "'scale_to_km' must be either TRUE or FALSE")
  })

})


test_that("variogram produces expected output", {

  expected_output <- c("variogram", "scale_to_km", "n_permutations", "breaks")
  expected_columns <- c("mid_points", "obs_vari", "n_obs", "lower_bound", "upper_bound")
  expected_columns_zero <- c("mid_points", "obs_vari", "n_obs")

  result <- variogram(sf_data, variable = "y", n_bins = 10, n_permutations = 100)

  expect_s3_class(result, "RiskMap_variogram")
  expect_setequal(names(result), expected_output)
  expect_setequal(names(result$variogram), expected_columns)
  expect_equal(nrow(result$variogram), 10)
  expect_false(result$scale_to_km)
  expect_equal(result$n_permutations, 100)
  expect_length(result$breaks, 10 + 1)
  expect_equal(mean(result$variogram$mid_points), mean(result$breaks))

  breaks <- seq(0, 35000, 2500)
  result <- variogram(sf_data, variable = "y", breaks = breaks, n_permutations = 100)
  expect_s3_class(result, "RiskMap_variogram")
  expect_setequal(names(result), expected_output)
  expect_setequal(names(result$variogram), expected_columns)
  expect_false(result$scale_to_km)
  expect_equal(result$n_permutations, 100)
  expect_equal(nrow(result$variogram), length(breaks) - 1)
  expect_equal(result$breaks, breaks)
  expect_equal(mean(result$variogram$mid_points), mean(result$breaks))

  result <- variogram(sf_data, variable = "y", n_bins = 10, n_permutations = 0)
  expect_s3_class(result, "RiskMap_variogram")
  expect_setequal(names(result), expected_output)
  expect_setequal(names(result$variogram), expected_columns_zero)
  expect_equal(nrow(result$variogram), 10)
  expect_false(result$scale_to_km)
  expect_equal(result$n_permutations, 0)
  expect_length(result$breaks, 10 + 1)
  expect_equal(mean(result$variogram$mid_points), mean(result$breaks))
})

test_that("variogram values are correct for a simple dataset", {

  square <- data.frame(
    x = c(1, 2, 1, 2, 10),
    y = c(1, 1, 2, 2, 10),
    z = c( 1, 2, 3, 4, 10)
  )

  # 10 points will be removed, so only the square distances are considered
  square_sf <- sf::st_as_sf(square, coords = c("x", "y"), crs = 32630)

  # independently calculate the variogram for all possible permutations
  # get all possible combinations of 1:5
  permutations <- expand.grid(1:5, 1:5, 1:5, 1:5, 1:5)
  permutations <- as.matrix(permutations[apply(permutations, 1, function(x) length(unique(x)) == 5), ])

  # compute semivariances for a permutation
  get_bin_averages <- function(perm) {
    vals <- c(1, 2, 3, 4, 10)[perm]

    s <- c(
      (vals[1]-vals[2])^2/2,  # 1-2: bin1
      (vals[1]-vals[3])^2/2,  # 1-3: bin1
      (vals[1]-vals[4])^2/2,  # 1-4: bin2
      (vals[2]-vals[3])^2/2,  # 2-3: bin2
      (vals[2]-vals[4])^2/2,  # 2-4: bin1
      (vals[3]-vals[4])^2/2   # 3-4: bin1
    )

    c(bin1 = mean(s[c(1,2,5,6)]),
      bin2 = mean(s[c(3,4)]))
  }

  manual_results <- t(apply(permutations, 1, get_bin_averages))
  manual_envelope <- apply(manual_results, 2, quantile, probs = c(0.025, 0.975))
  manual_observed <- get_bin_averages(1:5)

  result <- variogram(square_sf,
                      variable = "z",
                      breaks = seq(0.1, 1.5, 0.2),
                      convert_to_utm = FALSE)

  expect_equal(sum(result$variogram$n_obs), 6)

  bin1 <- result$variogram$mid_points == 1
  bin2 <- result$variogram$mid_points + 0.1 > sqrt(2)

  expect_equal(result$variogram$n_obs[bin1], 4)
  expect_equal(result$variogram$n_obs[bin2], 2)

  expect_equal(result$variogram$obs_vari[bin1], unname(manual_observed[1]))
  expect_equal(result$variogram$obs_vari[bin2], unname(manual_observed[2]))

  expect_equal(result$variogram$lower_bound[bin1], manual_envelope[1,1])
  expect_equal(result$variogram$lower_bound[bin2], manual_envelope[1,2])

  expect_equal(result$variogram$upper_bound[bin1], manual_envelope[2,1])
  expect_equal(result$variogram$upper_bound[bin2], manual_envelope[2,2])
})
