test_that("compiled log-likelihood agrees with stable R calculations", {
  eta <- matrix(c(
    -30, -3, 0, 3, 30,
    -12, -1, 0.5, 8, 20
  ), nrow = 2, byrow = TRUE)

  binomial_y <- c(0, 1, 3, 7, 10)
  binomial_units <- c(1, 4, 6, 9, 10)
  expected_binomial <- rowSums(
    sweep(eta, 2, binomial_y, `*`) -
      sweep(log1p(exp(-abs(eta))) + pmax(eta, 0), 2, binomial_units, `*`)
  )
  expect_equal(
    cpp_canonical_log_likelihood(eta, binomial_y, binomial_units,
                                 poisson = FALSE),
    expected_binomial, tolerance = 1e-12
  )

  poisson_y <- c(0, 1, 2, 5, 9)
  exposure <- c(0.5, 1, 1.5, 2, 3)
  poisson_eta <- eta / 3
  expected_poisson <- rowSums(
    sweep(poisson_eta, 2, poisson_y, `*`) -
      sweep(exp(poisson_eta), 2, exposure, `*`)
  )
  expect_equal(
    cpp_canonical_log_likelihood(poisson_eta, poisson_y, exposure,
                                 poisson = TRUE),
    expected_poisson, tolerance = 1e-12
  )
})

test_that("compiled MALA supports multiple random-effect blocks", {
  control <- set_control_mcmc(
    n_sim = 100, burnin = 40, thin = 3, h = 0.5, seed = 4821
  )
  fit <- laplace_sampling_mcmc(
    y = c(1, 2, 3, 4), units_m = rep(6, 4),
    mu = c(-1, -1, -0.5, -0.5),
    Sigma = matrix(c(0.8, 0.2, 0.2, 0.7), 2),
    ID_coords = c(1L, 1L, 2L, 2L),
    ID_re = cbind(
      cluster = c(1L, 1L, 2L, 2L),
      period = c(1L, 2L, 1L, 2L)
    ),
    sigma2_re = c(cluster = 0.3, period = 0.2),
    family = "binomial", control_mcmc = control, messages = FALSE
  )

  expect_named(fit$samples, c("S", "cluster", "period"))
  expect_identical(dim(fit$samples$S), c(20L, 2L))
  expect_identical(dim(fit$samples$cluster), c(20L, 2L))
  expect_identical(dim(fit$samples$period), c(20L, 2L))
  expect_true(all(is.finite(unlist(fit$samples))))
  expect_true(all(fit$tuning_par[41:100] == fit$tuning_par[40]))
})

test_that("compiled MALA reproduces the general R transition", {
  control <- set_control_mcmc(
    n_sim = 100, burnin = 40, thin = 3, h = 0.5, seed = 99
  )
  arguments <- list(
    y = c(1, 3), units_m = c(5, 6), mu = c(-0.5, 0.2),
    Sigma = matrix(c(0.8, 0.2, 0.2, 0.7), 2), ID_coords = 1:2,
    family = "binomial", control_mcmc = control, messages = FALSE
  )
  compiled <- do.call(laplace_sampling_mcmc, arguments)
  # Supplying the same link explicitly selects the general R implementation.
  r_reference <- do.call(
    laplace_sampling_mcmc,
    c(arguments, list(invlink = function(x) plogis(x)))
  )

  expect_identical(compiled$acceptance, r_reference$acceptance)
  # Equivalent floating-point update orders differ at machine precision on
  # some BLAS/platform combinations.
  expect_equal(compiled$tuning_par, r_reference$tuning_par,
               tolerance = 1e-14)
  # The two paths use algebraically equivalent derivative calculations.
  expect_equal(compiled$samples$S, r_reference$samples$S, tolerance = 5e-5)
})

test_that("compiled log-likelihood guards Poisson overflow", {
  overflow <- matrix(c(0, 800), nrow = 2)
  expect_identical(
    cpp_canonical_log_likelihood(overflow, 1, 1, poisson = TRUE)[2], -Inf
  )
})
