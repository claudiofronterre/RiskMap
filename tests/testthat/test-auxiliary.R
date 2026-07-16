test_that("check_data functions correctly", {

  data <- data.frame(
    x = c(1, 2, 3),
    y = c(0, 1, 2),
    z = c(0, 1, 0)
  )

  sf_data <- sf::st_as_sf(data, coords = c("x", "y"))

  polygon <- sf::st_polygon(list(matrix(c(4,4, 5,4, 5,5, 4,5, 4,4), ncol = 2, byrow = TRUE)))
  sf_polygon <- sf::st_sf(z = 3, geometry = sf::st_sfc(polygon))
  sf_merged <- rbind(sf_data, sf_polygon)

  expect_no_error(check_data(sf_data))
  expect_error(check_data(data), "'data' must be of class 'sf'")
  expect_error(check_data(sf_merged), "'data' can only contain point geometry")

})
