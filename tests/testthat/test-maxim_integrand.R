coords_mi <- cbind(c(0, 1, 2, 3), c(0, 0, 1, 1))
sigma_mi <- exp(-as.matrix(dist(coords_mi)))

test_that("maxim_integrand accepts the denominator positionally and by name", {
  named <- maxim_integrand(
    y = c(3, 5, 1, 0), denominator = rep(10, 4), mu = rep(-0.5, 4),
    Sigma = sigma_mi, ID_coords = 1:4, family = "binomial"
  )
  positional <- maxim_integrand(
    c(3, 5, 1, 0), rep(10, 4), rep(-0.5, 4), sigma_mi, 1:4,
    family = "binomial"
  )
  expect_identical(named, positional)
})

test_that("maxim_integrand rejects families other than binomial and Poisson", {
  expect_error(
    maxim_integrand(
      y = c(0.3, 1.2, -0.4, 0.8), denominator = rep(1, 4), mu = rep(0, 4),
      Sigma = sigma_mi, ID_coords = 1:4, family = "gaussian"
    ),
    "must be either \"binomial\" or \"poisson\""
  )
  expect_error(
    maxim_integrand(
      y = c(3, 5, 1, 0), mu = rep(-0.5, 4), Sigma = sigma_mi,
      ID_coords = 1:4, family = "binomial"
    ),
    "'denominator' must be provided"
  )
})

test_that("a Poisson denominator acts as an exposure", {
  # With a log link, an exposure e is equivalent to adding log(e) to mu.
  exposure <- c(0.5, 1, 2, 4)
  with_exposure <- maxim_integrand(
    y = c(1, 3, 2, 6), denominator = exposure, mu = rep(0.2, 4),
    Sigma = sigma_mi, ID_coords = 1:4, family = "poisson"
  )
  with_offset <- maxim_integrand(
    y = c(1, 3, 2, 6), denominator = rep(1, 4), mu = 0.2 + log(exposure),
    Sigma = sigma_mi, ID_coords = 1:4, family = "poisson"
  )
  expect_equal(with_exposure$mode, with_offset$mode, tolerance = 1e-6)
})
