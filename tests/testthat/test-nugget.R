test_that("nugget and fix_var_me cannot both be estimated when family is gaussian", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(0, 1, 2),
    z = c(0, 1, 2),
    den = c(4, 4, 4)
  )

  sf_data <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(4326))

  expect_error(glgpm(z ~ gp(nugget = TRUE), sf_data, "gaussian", messages = FALSE), "When there is only one observation per location")
  expect_no_error(glgpm(z ~ gp(nugget = TRUE), sf_data, "gaussian", fix_var_me = 1, messages = FALSE))
  expect_no_error(glgpm(z ~ gp(nugget = TRUE), sf_data, "binomial", den = den, messages = FALSE))

  two_location_data <- rbind(data,
                             data.frame(
                              x = 1,
                              y = 0,
                              z = 3,
                              den = 4
                            ))

  sf_data <- sf::st_as_sf(two_location_data, coords = c("x", "y"), crs = sf::st_crs(4326))

  expect_no_error(glgpm(z ~ gp(nugget = TRUE), sf_data, "gaussian", messages = FALSE))

  result <- glgpm(z ~ gp(nugget = TRUE), sf_data, "gaussian", fix_var_me = 1, messages = FALSE)
  expect_true("tau2" %in% names(coef(result)))

  result <- glgpm(z ~ gp(nugget = 1), sf_data, "gaussian", messages = FALSE)
  expect_equal(summary(result)$tau2, 1)

})
