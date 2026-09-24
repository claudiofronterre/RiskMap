polygon <- st_polygon(list(matrix(c(4,4, 5,4, 5,5, 4,5, 4,4), ncol = 2, byrow = TRUE)))
sf_polygon <- st_sf(geometry = st_sfc(polygon), crs = st_crs(4326))

square_coords <- c(40,40,
                   50,40,
                   50,50,
                   40,50,
                   40,40)

test_that("create_grid produces errors", {

  expect_error(create_grid("gaussian_data", 1), "'boundaries' must be of class 'sf'")
  expect_error(create_grid(gaussian_data, 1), "'boundaries' can only contain 'POLYGON' or 'MULTIPOLYGON' geometry")
  expect_error(create_grid(sf_polygon, "not numeric"), "The value for 'spacing' must be a single positive number")
  expect_error(create_grid(sf_polygon, 1, distance_units = "feet"),
               "'distance_units' must be either 'km' or 'm'")


  square <- st_polygon(list(matrix(square_coords, ncol = 2, byrow = TRUE)))
  sf_square <- st_sf(geometry = st_sfc(square), crs = st_crs(32638))
  expect_error(create_grid(sf_square, 1), "No grid-point centres fall within 'boundaries'")
})

test_that("create_grid functions correctly", {

  result <- create_grid(st_transform(sf_polygon, 32638), 10)
  expect_s3_class(result, "sf")
  expect_equal(as.character(unique(st_geometry_type(result))), "POINT")

  square <- st_polygon(list(matrix(square_coords * 1000, ncol = 2, byrow = TRUE)))
  sf_square <- st_sf(geometry = st_sfc(square), crs = st_crs(32638))
  result <- create_grid(sf_square, 10)
  expect_equal(nrow(result), 1)

  result <- create_grid(sf_square, 5)
  expect_equal(nrow(result), 4)

  result <- create_grid(sf_square, 1)
  expect_equal(nrow(result), 100)

  result <- create_grid(sf_square, 2.5)
  expect_equal(nrow(result), 16)

  # confirm for sf with multiple polygons
  square2 <- st_polygon(list(matrix((square_coords + 10) * 1000, ncol = 2, byrow = TRUE)))
  sf_square2 <- st_sf(geometry = st_sfc(square2), crs = st_crs(32638))
  sf_combined <- rbind(sf_square, sf_square2)
  result <- create_grid(sf_combined, 10)
  expect_equal(nrow(result), 2)

})

test_that("create_grid converts longitude/latitude data automatically, like glgpm()", {
  expect_message(
    result <- create_grid(sf_polygon, 1),
    "longitude/latitude"
  )
  expect_s3_class(result, "sf")
  expect_false(st_is_longlat(result))
  expect_equal(as.character(unique(st_geometry_type(result))), "POINT")
  expect_gt(nrow(result), 0)
})

test_that("create_grid converts spacing to the projected CRS units", {
  square <- st_polygon(list(matrix(square_coords * 1000,
                                   ncol = 2,
                                   byrow = TRUE)))
  square_m <- st_sf(geometry = st_sfc(square), crs = st_crs(32638))
  square_ft <- st_transform(
    square_m,
    "+proj=utm +zone=38 +datum=WGS84 +units=us-ft +no_defs"
  )

  grid_m <- create_grid(square_m, spacing = 2.5)
  grid_ft <- create_grid(square_ft, spacing = 2.5)

  expect_length(grid_m$geometry, 16)
  expect_length(grid_ft$geometry, 16)
  expect_equal(
    as.numeric(st_distance(grid_m[1, ], grid_m[2, ])),
    2500,
    tolerance = 1e-6
  )
  expect_equal(
    as.numeric(st_distance(grid_ft[1, ], grid_ft[2, ])),
    2500 / 0.304800609601219,
    tolerance = 1e-4
  )
})
