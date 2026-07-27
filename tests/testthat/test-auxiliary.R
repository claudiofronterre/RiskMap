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
