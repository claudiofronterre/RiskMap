test_that("check_formula functions correctly", {

  df <- data.frame(
    y = c(1, 2, 3),
    x = c(0, 1, 2),
    z = c(0, 1, 0),
    c = c(1, 2, 3),
    gp = c(1, 2, 3)
  )

  data <- sf::st_as_sf(df, coords = c("x", "z"), crs = 4326)

  test_kappa <- 1
  expect_no_error(check_formula(y ~ gp(), data))
  expect_no_error(check_formula(y ~ gp(kappa = test_kappa), data))
  expect_no_error(check_formula(y ~ log(c) + gp(), data))

  expect_error(check_formula("not formula", data), "'formula' must be a 'formula'")
  expect_error(check_formula(~ gp(), data), "The 'formula' must specify a response")
  expect_error(check_formula(~ gp(), data, response_required = FALSE),
               "The 'formula' must specify a response")
  expect_error(check_formula(y ~ c, data), "The 'formula' must contain a Gaussian Process term")
  expect_error(check_formula(y ~ gp, data), "The 'formula' must contain a Gaussian Process term")
  expect_error(check_formula(y ~ gp(xx, c), data), "The 'formula' term 'xx'")
  expect_error(check_formula(y ~ gp(xx, zz), data), "The 'formula' terms 'xx', 'zz'")
  expect_error(check_formula(y ~ gp(c) + re(xx, zz), data), "The 'formula' terms 'xx', 'zz'")

  # c not included in formula so should not error
  data$c[1] <- NA
  expect_no_error(check_formula(y ~ gp(), data))

  data$y[1] <- NA
  expect_error(check_formula(y ~ gp(), data), "'data' contains rows with missing data")

  data_without_response <- data[, "c"]
  data_without_response$c[1] <- 1
  expect_no_error(check_formula(y ~ c + gp(), data_without_response,
                                response_required = FALSE))

  data_without_response$c[1] <- NA
  expect_error(check_formula(y ~ c + gp(), data_without_response,
                             response_required = FALSE),
               "'data' contains rows with missing data")
 })

test_that("check_complete_data only checks the specified columns and names the caller's argument", {

  df <- data.frame(a = c(1, 2, 3), b = c(1, NA, 3))

  # "b" has a missing value, but it is ignored when only checking "a"
  expect_no_error(check_complete_data(df, "a"))
  expect_error(check_complete_data(df, "b"), "'df' contains rows with missing data")
  expect_error(check_complete_data(df, c("a", "b")), "'df' contains rows with missing data")
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

  gaussian_data <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(4326))
  sf_no_crs <- sf::st_as_sf(data, coords = c("x", "y"))

  data$y[2] <- 100
  sf_wrong_coord <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(4326))

  polygon <- sf::st_polygon(list(matrix(c(4,4, 5,4, 5,5, 4,5, 4,4), ncol = 2, byrow = TRUE)))
  multipolygon <- sf::st_multipolygon(list(polygon))
  sf_polygon <- sf::st_sf(z = 3, geometry = sf::st_sfc(polygon), crs = sf::st_crs(4326))
  sf_multipolygon <- sf::st_sf(z = 3, geometry = sf::st_sfc(multipolygon), crs = sf::st_crs(4326))
  sf_merged <- rbind(gaussian_data, sf_polygon)

  expect_error(check_data(gaussian_data, geometry = "n"), "'geometry' must be either 'point', 'polygon' or 'any'")
  expect_error(check_data(gaussian_data, type = "n"), "'type' must be either 'sf' or 'sfc'")

  expect_no_error(check_data(gaussian_data))
  expect_error(check_data(data), "'data' must be of class 'sf'")

  expect_error(check_data(sf_no_crs), "'sf_no_crs' must contain a coordinate reference system")
  expect_error(check_data(sf_merged), "'sf_merged' can only contain 'POINT' geometry")
  expect_error(check_data(sf_wrong_coord), "'sf_wrong_coord' contains impossible latitude or longitude values")

  expect_no_error(check_data(sf_polygon, "polygon"))
  expect_error(check_data(gaussian_data, "polygon"), "'gaussian_data' can only contain 'POLYGON' or 'MULTIPOLYGON' geometry")
  expect_error(check_data(sf_merged, "polygon"), "'sf_merged' can only contain 'POLYGON' or 'MULTIPOLYGON' geometry")

  expect_no_error(check_data(sf_multipolygon, "polygon"))

  expect_no_error(check_data(st_geometry(gaussian_data), "point", "sfc"))
  expect_no_error(check_data(st_geometry(sf_multipolygon), "polygon", "sfc"))

})

test_that("check_crs functions correctly", {
  crs <- "invalid"
  expect_error(check_crs(crs), "The 'crs' provided is not a valid CRS")
  crs <- 121212
  expect_error(check_crs(crs), "The 'crs' provided is not a valid CRS")
  dif_crs <- 12121.2
  expect_error(check_crs(dif_crs), "The 'dif_crs' provided is not a valid CRS")

  crs <- 2648
  expect_no_error(check_crs(crs))
  expect_no_error(check_crs(2648))

  custom_crs <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
  expect_no_error(check_crs(custom_crs))
})

test_that("gp functions correctly", {

  expected_output <- c("term", "kappa", "nugget", "dim", "label")

  expect_error(gp(nugget = 0), "'nugget' must be either 'TRUE'")
  expect_error(gp(nugget = -1), "'nugget' must be either 'TRUE'")
  expect_error(gp(nugget = "TRUE"), "'nugget' must be either 'TRUE'")

  expect_error(gp(kappa = 0), "'kappa' must be positive")
  expect_error(gp(kappa = "TRUE"), "'kappa' must be positive")

  default_result <- gp()

  expect_setequal(names(default_result), expected_output)
  expect_identical(default_result$term, "sf")
  expect_identical(default_result$kappa, 0.5)
  expect_identical(default_result$nugget, 0)
  expect_identical(default_result$dim, 0L)
  expect_identical(default_result$label, "gp()")

  true_result <- gp(nugget = TRUE)
  expect_true(true_result$nugget)

  custom_result <- gp(a, b, kappa = 2, nugget = 3)
  expect_identical(custom_result$term, c("a", "b"))
  expect_identical(custom_result$kappa, 2)
  expect_identical(custom_result$nugget, 3)
  expect_identical(custom_result$dim, 2L)
  expect_identical(custom_result$label, "gp(a,b)")

})

test_that("check_positive_integer functions correctly", {
  expect_no_error(check_positive_integer(1, "a"))
  expect_no_error(check_positive_integer(999, "a"))
  expect_no_error(check_positive_integer(NULL, "a", allow_null = TRUE))
  expect_no_error(check_positive_integer(0, "a", allow_zero = TRUE))

  expect_error(check_positive_integer(0.1, "a"), "'a' must be a single positive integer")
  expect_error(check_positive_integer(0, "a"), "'a' must be a single positive integer")
  expect_error(check_positive_integer(c(0, 1), "a"), "'a' must be a single positive integer")
  expect_error(check_positive_integer("not", "a"), "'a' must be a single positive integer")
  expect_error(check_positive_integer(NULL, "a"), "'a' must be a single positive integer")
  expect_error(check_positive_integer(NA, "a"), "'a' must be a single positive integer")
  expect_error(check_positive_integer(Inf, "a"), "'a' must be a single positive integer")
})

test_that("check_logical functions correctly and names the caller's argument", {
  my_flag <- TRUE
  expect_no_error(check_logical(my_flag))

  my_flag <- FALSE
  expect_no_error(check_logical(my_flag))

  my_flag <- NA
  expect_error(check_logical(my_flag), "'my_flag' must be either TRUE or FALSE")

  my_flag <- c(TRUE, FALSE)
  expect_error(check_logical(my_flag), "'my_flag' must be either TRUE or FALSE")

  my_flag <- "TRUE"
  expect_error(check_logical(my_flag), "'my_flag' must be either TRUE or FALSE")

  my_flag <- 1
  expect_error(check_logical(my_flag), "'my_flag' must be either TRUE or FALSE")
})

test_that("check_positive_number functions correctly", {
  expect_no_error(check_positive_number(1, ""))
  expect_no_error(check_positive_number(0.1, ""))

  a <- 0
  expect_error(check_positive_number(a, ""), "The value for 'a' must be a single positive number")
  a <- c(0, 1)
  expect_error(check_positive_number(a, ""), "The value for 'a' must be a single positive number")
  a <- "not"
  expect_error(check_positive_number(a, ""), "The value for 'a' must be a single positive number")
  a <- NULL
  expect_error(check_positive_number(a, ""), "The value for 'a' must be a single positive number")
  a <- NA
  expect_error(check_positive_number(a, ""), "The value for 'a' must be a single positive number")
  a <- Inf
  expect_error(check_positive_number(a, ""), "The value for 'a' must be a single positive number")
})

test_that("check_zero_one validates the closed unit interval", {
  expect_no_error(check_zero_one(0, "a"))
  expect_no_error(check_zero_one(0.5, "a"))
  expect_no_error(check_zero_one(1, "a"))

  expect_error(check_zero_one(-0.1, "a"), "between 0 and 1")
  expect_error(check_zero_one(1.1, "a"), "between 0 and 1")
  expect_error(check_zero_one(c(0, 1), "a"), "between 0 and 1")
  expect_error(check_zero_one(NA, "a"), "between 0 and 1")
  expect_error(check_zero_one(Inf, "a"), "between 0 and 1")
})

test_that("plot_metric validates model before metric (#156)", {
  ## an invalid model should be reported even when the metric is also invalid -
  ## previously this incorrectly complained about the metric first
  expect_error(
    plot_metric(cross_validation, metric = "not_a_metric", model = "not_a_model"),
    "'model' 'not_a_model' was not found"
  )

  expect_error(
    plot_metric(cross_validation, metric = "not_a_metric", model = "model_a"),
    "'metric' 'not_a_metric' was not computed for model 'model_a'"
  )
})

test_that("plot_metric plots the requested metric for the requested model", {
  p <- plot_metric(cross_validation, metric = "SCRPS", model = "model_b")

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_metric applies extra ggplot components passed via ...", {
  p_plain <- plot_metric(cross_validation, metric = "CRPS", model = "model_a")
  expect_length(p_plain$scales$scales, 0)

  p_custom <- plot_metric(cross_validation, metric = "CRPS", model = "model_a",
                          ggplot2::scale_color_gradient(low = "yellow", high = "red"))
  expect_length(p_custom$scales$scales, 1)
  expect_no_error(ggplot2::ggplot_build(p_custom))
})

test_that("create_convex_hull generates errors", {
  expect_error(create_convex_hull("a"), "'sf_object' must be an sf object")
})

test_that("create_convex_hull functions correctly", {
  result <- create_convex_hull(gaussian_data)
  expect_s3_class(result, "sf")
  expect_identical(nrow(result), 1L)

  intersect <- unlist(st_intersects(gaussian_data, result))
  expect_identical(intersect, rep(1L, 10))
})

test_that("propose_utm generates errors", {
  expect_error(propose_utm("a"), "'data' must be of class 'sf'")
})

test_that("propose_utm handles both polygons and points", {
  expect_no_error(propose_utm(areal))
  expect_no_error(propose_utm(gaussian_data))
})

test_that("propose_utm warns when locations are split equally", {
  lonlat_points <- function(lon, lat) {
    st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  }

  # longitudes 33 and 39 fall in UTM zones 36 and 37
  expect_warning(
    epsg <- propose_utm(lonlat_points(lon = c(33, 39), lat = c(5, 5))),
    "An equal amount of locations falls in different UTM zones"
  )
  # default to min
  expect_identical(epsg, 32636)

  expect_warning(
    epsg <- propose_utm(lonlat_points(lon = c(33, 33), lat = c(5, -5))),
    "An equal amount of locations falls north and south of the Equator"
  )
  # default to south
  expect_identical(epsg, 32736)

  expect_no_warning(
    epsg <- propose_utm(lonlat_points(lon = c(33, 33, 39), lat = c(5, 5, 5)))
  )
  expect_identical(epsg, 32636)

  expect_no_warning(
    epsg <- propose_utm(lonlat_points(lon = c(33, 33, 33, 33), lat = c(5, 6, 7, -5)))
  )
  expect_identical(epsg, 32636)
})

test_that("re errors when no variables are provided",{
  expect_error(re(), "You need to provide at least one variable to re")
})

test_that("re produces expected output", {
  result <- re(variable)
  expect_s3_class(result, "RiskMap_re_spec")
  expect_identical(result$term, "variable")
  expect_identical(result$dim, 1L)
  expect_identical(result$label, "re(variable)")

  result <- re(variable, region)
  expect_s3_class(result, "RiskMap_re_spec")
  expect_identical(result$term, c("variable", "region"))
  expect_identical(result$dim, 2L)
  expect_identical(result$label, "re(variable,region)")
})

