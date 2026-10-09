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

  expect_error(check_data(gaussian_data, geometry = "n"), "'geometry' must be either 'point' or 'polygon'")
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
  expect_equal(default_result$term, "sf")
  expect_equal(default_result$kappa, 0.5)
  expect_equal(default_result$nugget, 0)
  expect_equal(default_result$dim, 0)
  expect_equal(default_result$label, "gp()")

  true_result <- gp(nugget = TRUE)
  expect_equal(true_result$nugget, TRUE)

  custom_result <- gp(a, b, kappa = 2, nugget = 3)
  expect_equal(custom_result$term, c("a", "b"))
  expect_equal(custom_result$kappa, 2)
  expect_equal(custom_result$nugget, 3)
  expect_equal(custom_result$dim, 2)
  expect_equal(custom_result$label, "gp(a,b)")

})

test_that("estimates are consistent between coef and summary", {

  sum <- summary(gaussian_model)
  cof <- coef(gaussian_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_me, sum$me[1,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)

  sum <- summary(binomial_model)
  cof <- coef(binomial_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)

  sum <- summary(poisson_model)
  cof <- coef(poisson_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)
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

test_that("plot.RiskMap_cross_validation validates model before metric (#156)", {
  ## an invalid model should be reported even when the metric is also invalid -
  ## previously this incorrectly complained about the metric first
  expect_error(
    plot(cross_validation, metric = "not_a_metric", model = "not_a_model"),
    "'model' 'not_a_model' was not found"
  )

  expect_error(
    plot(cross_validation, metric = "not_a_metric", model = "model_a"),
    "'metric' 'not_a_metric' was not computed for model 'model_a'"
  )
})

test_that("plot.RiskMap_cross_validation plots the requested metric for the requested model", {
  p <- plot(cross_validation, metric = "SCRPS", model = "model_b")

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot.RiskMap_cross_validation applies extra ggplot components passed via ...", {
  p_plain <- plot(cross_validation, metric = "CRPS", model = "model_a")
  expect_length(p_plain$scales$scales, 0)

  p_custom <- plot(cross_validation, metric = "CRPS", model = "model_a",
                   ggplot2::scale_color_gradient(low = "yellow", high = "red"))
  expect_length(p_custom$scales$scales, 1)
  expect_no_error(ggplot2::ggplot_build(p_custom))
})

test_that("plot.RiskMap_cross_validation plots each metric separately rather than one combined grid (#87)", {
  ## the fixture was fit with the default metrics = c("PIT", "CRPS", "SCRPS"),
  ## and "PIT" also computes the scalar "PIT_area" score, so there are 4
  ## metrics here. With more than one metric, each gets its own plot/grid
  ## (printed in turn as a side effect) instead of everything being squeezed
  ## into one combined grid; the list returned invisibly has one entry per
  ## metric.
  groups <- plot(cross_validation)

  expect_type(groups, "list")
  expect_named(groups, c("PIT", "CRPS", "SCRPS", "PIT_area"))
  for (g in groups) {
    expect_s3_class(g, "gtable")
    expect_equal(length(g$grobs), 2) # one panel per model, no padding
  }
})

test_that("plot.RiskMap_cross_validation sizes each metric's grid from its own models only (#87)", {
  ## 3 models, 2 metrics - each metric's grid should be sized from its own
  ## 3 models, independent of the other metric (no cross-metric padding)
  cv3 <- cross_validation
  cv3$model$model_c <- cv3$model$model_a

  groups <- plot(cv3, metric = c("PIT", "CRPS"))

  expect_named(groups, c("PIT", "CRPS"))
  for (g in groups) {
    expect_s3_class(g, "gtable")
    expect_equal(length(g$grobs), 3)
  }
})

test_that("plot.RiskMap_cross_validation can select only the calibration curve", {
  p <- plot(cross_validation, metric = "PIT", model = "model_a")

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot.RiskMap_cross_validation's calibration curve actually filters to the requested test_set", {

  p1 <- plot(cross_validation, metric = "PIT", model = "model_a",
            pit_mode = "single", pit_test_set = 1)
  p2 <- plot(cross_validation, metric = "PIT", model = "model_a",
            pit_mode = "single", pit_test_set = 2)

  expect_false(identical(p1$data$value, p2$data$value))
  expect_equal(p1$labels$title, "Model model_a: PIT (test set 1)")
  expect_error(
    plot(cross_validation, metric = "PIT", model = "model_a",
        pit_mode = "single", pit_test_set = 99),
    "No data for test set 99"
  )
})

test_that("plot.RiskMap_cross_validation combines calibration curves with combine_pit (#183)", {
  p <- plot(cross_validation, metric = "PIT", combine_pit = TRUE)

  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$model), c("model_a", "model_b"))
})

test_that("plot.RiskMap_cross_validation labels metric maps by metric and model (#183)", {
  p <- plot(cross_validation, metric = "CRPS", model = "model_a")

  expect_equal(p$labels$title, "Model model_a: CRPS")
  expect_equal(p$labels$colour, "CRPS")
})

test_that("plot.RiskMap_cross_validation labels PIT as AnPIT for discrete families (#183)", {
  p <- plot(cross_validation, metric = "PIT_area", model = "model_a")
  expect_equal(p$labels$title, "Model model_a: PIT_area")
  expect_equal(p$labels$colour, "PIT_area")

  ## discrete families store AnPIT curves rather than PIT values
  discrete_cv <- cross_validation
  discrete_cv$model$model_a$AnPIT <- list(seq(0, 1, length.out = 11),
                                          seq(0, 1, length.out = 11)^2)
  discrete_cv$model$model_a$PIT <- NULL

  p_area <- plot(discrete_cv, metric = "PIT_area", model = "model_a")
  expect_equal(p_area$labels$title, "Model model_a: AnPIT_area")
  expect_equal(p_area$labels$colour, "AnPIT_area")

  p_curve <- plot(discrete_cv, metric = "PIT", model = "model_a")
  expect_equal(p_curve$labels$title, "Model model_a: AnPIT (average)")
  expect_equal(p_curve$labels$y, "AnPIT")
})
