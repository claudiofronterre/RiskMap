test_that("check_formula functions correctly", {

  df <- data.frame(
    y = c(1, 2, 3),
    x = c(0, 1, 2),
    z = c(0, 1, 0),
    c = c(1, 2, 3),
    gp = c(1, 2, 3)
  )

  data <- sf::st_as_sf(df, coords = c("x", "z"), crs = 4326)

  expect_no_error(check_formula(y ~ gp(), data))
  expect_error(check_formula("not formula", data), "'formula' must be a 'formula'")
  expect_error(check_formula(y ~ c, data), "The 'formula' must contain a Gaussian Process term")
  expect_error(check_formula(y ~ gp, data), "The 'formula' must contain a Gaussian Process term")
  expect_error(check_formula(y ~ gp(xx, c), data), "The 'formula' term 'xx'")
  expect_error(check_formula(y ~ gp(xx, zz), data), "The 'formula' terms 'xx', 'zz'")
  expect_error(check_formula(y ~ gp(c) + re(xx, zz), data), "The 'formula' terms 'xx', 'zz'")
 })
  
test_that("check_binomial functions correctly", {

  expect_no_error(check_binomial(0:3, NULL))
  expect_no_error(check_binomial(0:3, 0:3))
  expect_no_error(check_binomial(0:3, 1:4))
  expect_no_error(check_binomial(c(0, 1.00000001, 2), NULL))
  expect_error(check_binomial(-1:2, NULL), "'y' must only consist of zero or positive integers")
  expect_error(check_binomial(c(0, 1.1, 2), NULL), "'y' must only consist of zero or positive integers")
  expect_error(check_binomial(1:4, 0:3), "Values of 'den' must be greater")
})

test_that("check_data functions correctly", {

  data <- data.frame(
    x = c(1, 2, 3),
    y = c(0, 1, 2),
    z = c(0, 1, 0)
  )

  sf_data <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(4326))
  sf_no_crs <- sf::st_as_sf(data, coords = c("x", "y"))

  data$y[2] <- 100
  sf_wrong_coord <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(4326))

  polygon <- sf::st_polygon(list(matrix(c(4,4, 5,4, 5,5, 4,5, 4,4), ncol = 2, byrow = TRUE)))
  sf_polygon <- sf::st_sf(z = 3, geometry = sf::st_sfc(polygon), crs = sf::st_crs(4326))
  sf_merged <- rbind(sf_data, sf_polygon)

  expect_no_error(check_data(sf_data))
  expect_error(check_data(data), "'data' must be of class 'sf'")
  expect_error(check_data(sf_no_crs), "'data' must contain a coordinate reference system")
  expect_error(check_data(sf_merged), "'data' can only contain point geometry")
  expect_error(check_data(sf_wrong_coord), "'data' contains impossible latitude or longitude values")
})

