polygon <- st_polygon(list(matrix(c(4,4, 5,4, 5,5, 4,5, 4,4), ncol = 2, byrow = TRUE)))
sf_polygon <- st_sf(geometry = st_sfc(polygon), crs = st_crs(4326))

square_coords <- c(40,40,
                   50,40,
                   50,50,
                   40,50,
                   40,40)

test_that("create_grid produces errors", {

  expect_error(create_grid("sf_data", 1), "'shp' must be of class 'sf'")
  expect_error(create_grid(sf_data, 1), "'shp' can only contain 'POLYGON' or 'MULTIPOLYGON' geometry")
  expect_error(create_grid(sf_polygon, 1), "The coordinates of 'shp' are in longitude and latitude")
  expect_error(create_grid(sf_polygon, 1, 12121.2), "'grid_crs' must be a single positive integer")
  expect_error(create_grid(sf_polygon, 1, 121212), "The 'grid_crs' provided is not a valid CRS")

  square <- st_polygon(list(matrix(square_coords, ncol = 2, byrow = TRUE)))
  sf_square <- st_sf(geometry = st_sfc(square), crs = st_crs(32638))
  expect_error(create_grid(sf_square, 1), "No points intersect with the 'shp'")
})


test_that("create_grid functions correctly", {

  result <- create_grid(sf_polygon, 10, 32638)
  expect_s3_class(result, "sfc_POINT")

  square <- st_polygon(list(matrix(square_coords * 1000, ncol = 2, byrow = TRUE)))
  sf_square <- st_sf(geometry = st_sfc(square), crs = st_crs(32638))
  result <- create_grid(sf_square, 10)
  expect_length(result, 1)

  result <- create_grid(sf_square, 5)
  expect_length(result, 4)

  result <- create_grid(sf_square, 1)
  expect_length(result, 100)

  # confirm for sf with multiple polygons
  square2 <- st_polygon(list(matrix((square_coords + 10) * 1000, ncol = 2, byrow = TRUE)))
  sf_square2 <- st_sf(geometry = st_sfc(square2), crs = st_crs(32638))
  sf_combined <- rbind(sf_square, sf_square2)
  result <- create_grid(sf_combined, 10)
  expect_length(result, 2)

})
