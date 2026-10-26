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

# A profile problem with an exact log-likelihood in two parameters.
make_toy_profile <- function(loglik, se = c(1, 1)) {
  names <- c("beta.a", "b")
  problem <- list(
    value = loglik, estimate = stats::setNames(c(0, 0), names), names = names,
    labels = list(beta.a = "a", b = "b"),
    standard_errors = stats::setNames(se, names),
    lower = c(-30, -30), upper = c(30, 30)
  )
  problem$control <- make_optimiser_control(1)
  best <- maximise_profile(problem, integer(), numeric(), c(0, 0))
  problem$bank_estimate <- best$par
  problem$max_loglik <- best$value
  problem
}

test_that("profile intervals are exact for a quadratic log-likelihood", {
  problem <- make_toy_profile(function(w) -0.5 * w[1]^2 - 0.5 * (w[2] - w[1])^2)
  profile <- profile_one_parameter(problem, "beta.a", 5L, 3)
  intervals <- profile_intervals(
    problem, list(beta.a = add_profile_deviance(profile, problem$max_loglik)),
    0.95)
  expect_equal(c(intervals$profile_lower, intervals$profile_upper),
               c(-1, 1) * stats::qnorm(0.975), tolerance = 1e-5)
  expect_type(profile$second_start, "logical")
})

test_that("a flat profile is extended with doubling steps and left open", {
  problem <- make_toy_profile(function(w) -0.001 * w[1]^2 - 0.5 * w[2]^2)
  profile <- profile_one_parameter(problem, "beta.a", 5L, 3)
  # Two regular steps of 1.5, then steps of 3, 6 and 12, five points a side.
  expect_equal(sort(profile$value),
               c(-24, -12, -6, -3, -1.5, 0, 1.5, 3, 6, 12, 24), tolerance = 1e-8)
  intervals <- profile_intervals(
    problem, list(beta.a = add_profile_deviance(profile, problem$max_loglik)),
    0.95)
  expect_true(is.na(intervals$profile_lower) && is.na(intervals$profile_upper))
})

test_that("a crossing between distant points is refined", {
  # The deviance is quadratic near the estimate and much steeper beyond 2,
  # so the interval ends are at 1.96 but the grid jumps over them.
  loglik <- function(w) {
    a <- abs(w[1])
    -0.5 * a^2 - 5 * max(a - 2, 0)^3 - 0.5 * w[2]^2
  }
  problem <- make_toy_profile(loglik, se = c(0.5, 1))
  profile <- profile_one_parameter(problem, "beta.a", 5L, 1)
  intervals <- profile_intervals(
    problem, list(beta.a = add_profile_deviance(profile, problem$max_loglik)),
    0.95)
  expect_equal(intervals$profile_upper, stats::qnorm(0.975), tolerance = 1e-4)
})
