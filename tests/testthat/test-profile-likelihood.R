test_that("profiles of a binomial fit have their grids, deviances and intervals", {
  profiles <- profile(binomial_model, parameters = c("cov", "sigma2"),
                      pairs = list(c("sigma2", "phi")), n_points = 5,
                      n_grid = 3, n_draws = 200, messages = FALSE)
  expect_s3_class(profiles, "RiskMap_profile")
  expect_named(profiles$marginal, c("beta.cov", "sigma2"))
  expect_named(profiles$joint, "sigma2:phi")
  expect_equal(nrow(profiles$joint[[1]]), 9L)
  for (marginal in profiles$marginal) {
    expect_gte(nrow(marginal), 5L)
    expect_gte(min(marginal$deviance, na.rm = TRUE), -1e-8)
  }
  intervals <- profiles$intervals
  expect_equal(intervals$parameter, c("beta.cov", "sigma2"))
  inside <- with(intervals, (is.na(profile_lower) | profile_lower <= estimate) &
                   (is.na(profile_upper) | estimate <= profile_upper))
  expect_true(all(inside))
  # Variances are reported on their natural scale.
  expect_equal(intervals$estimate[2], exp(binomial_model$estimate$sigma2))
})

test_that("profiling restores the random number state", {
  set.seed(42)
  expected <- stats::runif(1)
  set.seed(42)
  invisible(profile(binomial_model, parameters = "cov", pairs = list(),
                    n_points = 5, n_draws = 200, messages = FALSE))
  expect_equal(stats::runif(1), expected)
})

test_that("profile rejects fits and settings it cannot handle", {
  expect_error(profile(gaussian_model), "binomial and Poisson")
  expect_error(profile(binomial_model, parameters = "kappa", messages = FALSE),
               "Unknown parameter")
  expect_error(profile(binomial_model, n_points = 2), "n_points")
  expect_error(profile(binomial_model, level = 1.5), "level")
})

test_that("the profile plot draws its panels and returns them", {
  profiles <- profile(binomial_model, parameters = c("sigma2", "phi"),
                      pairs = list(c("sigma2", "phi")), n_points = 5,
                      n_grid = 3, n_draws = 200, messages = FALSE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  panels <- plot(profiles)
  expect_named(panels, c("sigma2", "phi", "sigma2:phi"))
  expect_true(all(vapply(panels, ggplot2::is_ggplot, logical(1))))
  expect_output(print(profiles), "Profile likelihood")
})
