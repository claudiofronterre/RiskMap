test_that("conditions carry the RiskMap classes", {
  raise <- function() riskmap_abort("Nothing to fit.", class = "riskmap_error_fit")
  error <- expect_error(raise(), class = "riskmap_error_fit")
  expect_s3_class(error, c("riskmap_error_fit", "riskmap_error", "rlang_error"))
  expect_identical(conditionCall(error), quote(raise()))

  expect_warning(riskmap_warn("Check the fit."), class = "riskmap_warning")
  expect_message(riskmap_inform("Fitting."), class = "riskmap_message")

  bug <- function() riskmap_bug("The cache is empty.")
  error <- expect_error(bug(), class = "riskmap_error_internal")
  expect_s3_class(error, "riskmap_error")
  expect_match(conditionMessage(error), "report it", fixed = TRUE)
})

test_that("check_positive_integer functions correctly", {
  expect_no_error(check_positive_integer(1, "a"))
  expect_no_error(check_positive_integer(999, "a"))
  expect_no_error(check_positive_integer(NULL, "a", allow_null = TRUE))
  expect_no_error(check_positive_integer(0, "a", allow_zero = TRUE))
  expect_identical(check_positive_integer(5, "a"), 5)

  for (value in list(0.1, 0, c(0, 1), "not", NULL, NA, Inf)) {
    expect_input_error(check_positive_integer(value, "a"), "a")
  }
  expect_error(check_positive_integer(-1, "a", allow_zero = TRUE),
               "non-negative", class = "riskmap_error_input")
})

test_that("check_logical functions correctly and names the caller's argument", {
  my_flag <- TRUE
  expect_no_error(check_logical(my_flag))

  my_flag <- FALSE
  expect_no_error(check_logical(my_flag))

  for (my_flag in list(NA, c(TRUE, FALSE), "TRUE", 1)) {
    expect_input_error(check_logical(my_flag), "my_flag")
  }
})

test_that("check_positive_number functions correctly", {
  expect_no_error(check_positive_number(1, ""))
  expect_no_error(check_positive_number(0.1, ""))

  for (a in list(0, c(0, 1), "not", NULL, NA, Inf)) {
    expect_input_error(check_positive_number(a, ""), "a")
  }
  start_pars <- list(sigma2 = -1)
  error <- expect_input_error(check_positive_number(start_pars[["sigma2"]]),
                              'start_pars[["sigma2"]]')
  expect_match(conditionMessage(error), "starting value for `sigma2`",
               fixed = TRUE)
})

test_that("check_zero_one validates the closed unit interval", {
  expect_no_error(check_zero_one(0, "a"))
  expect_no_error(check_zero_one(0.5, "a"))
  expect_no_error(check_zero_one(1, "a"))

  for (value in list(-0.1, 1.1, c(0, 1), NA, Inf)) {
    expect_input_error(check_zero_one(value, "a"), "a")
  }
})

test_that("check_class accepts any of the classes", {
  fit <- structure(list(), class = c("RiskMap", "list"))
  expect_identical(check_class(fit, "RiskMap"), fit)
  expect_no_error(check_class(fit, c("RiskMap_prediction", "RiskMap")))
  expect_no_error(check_class(NULL, "RiskMap", allow_null = TRUE))
  expect_input_error(check_class(data.frame(), "RiskMap", arg = "object"),
                     "object")
  expect_input_error(check_class(NULL, "RiskMap", arg = "object"), "object")
})

test_that("match_choice works as match.arg does", {
  choose <- function(family = c("gaussian", "binomial", "poisson")) {
    match_choice(family)
  }
  expect_identical(choose(), "gaussian")
  expect_identical(choose("binomial"), "binomial")
  expect_identical(choose("pois"), "poisson")
  expect_input_error(choose("gamma"), "family")
  expect_input_error(choose(c("gaussian", "binomial")), "family")
  expect_input_error(choose(NA_character_), "family")
  expect_identical(match_choice("km", c("km", "m"), arg = "distance_units"),
                   "km")
  expect_input_error(match_choice("cm", c("km", "m"), arg = "distance_units"),
                     "distance_units")
  no_default <- function(units) match_choice(units)
  expect_error(no_default("km"), class = "riskmap_error_internal")
})

test_that("check_crs functions correctly", {
  crs <- "invalid"
  error <- expect_input_error(check_crs(crs), "crs")
  # the error from sf is kept as the cause
  expect_s3_class(error$parent, "condition")
  crs <- 121212
  expect_input_error(check_crs(crs), "crs")
  dif_crs <- 12121.2
  expect_input_error(check_crs(dif_crs), "dif_crs")

  crs <- 2648
  expect_no_error(check_crs(crs))
  expect_no_error(check_crs(2648))

  custom_crs <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
  expect_no_error(check_crs(custom_crs))
})

test_that("input errors name the calling function and describe the value", {
  fit_example <- function(n_iter, verbose, level = 0.95) {
    check_positive_integer(n_iter, "n_iter")
    check_logical(verbose)
    check_range(level, min = 0, max = 1, allow_equal = FALSE)
    invisible(TRUE)
  }
  expect_snapshot(error = TRUE, {
    fit_example(2.5, TRUE)
    fit_example(c(1, 2), TRUE)
    fit_example(10, "yes")
    fit_example(10, TRUE, level = 1)
  })
})

test_that("check_range leaves infinite bounds out of the message", {
  range_example <- function(max_dist, weight = 1, shift = -1) {
    check_range(max_dist, min = 0, max = Inf, allow_equal = FALSE)
    check_range(weight, min = 0)
    check_range(shift, min = -Inf, max = 2, allow_equal = FALSE)
    invisible(TRUE)
  }
  expect_snapshot(error = TRUE, {
    range_example(-5)
    range_example(Inf)
    range_example(10, weight = -1)
    range_example(10, shift = 2)
  })
  expect_error(check_range(1, min = 2, arg = "a"),
               "greater than or equal to 2", class = "riskmap_error_input")
  expect_error(check_range(Inf, arg = "a"), "single finite number, not Inf",
               class = "riskmap_error_input")
})

test_that("class and choice errors read well", {
  summarise_example <- function(object, units = c("km", "m")) {
    check_class(object, c("RiskMap", "RiskMap_prediction"))
    match_choice(units)
  }
  fit <- structure(list(), class = "RiskMap")
  expect_snapshot(error = TRUE, {
    summarise_example(data.frame())
    summarise_example(fit, units = "miles")
  })
})
