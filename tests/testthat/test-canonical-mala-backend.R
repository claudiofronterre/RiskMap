test_that("compiled canonical likelihood agrees with stable R calculations", {
  eta <- matrix(c(
    -30, -3, 0, 3, 30,
    -12, -1, 0.5, 8, 20
  ), nrow = 2, byrow = TRUE)

  binomial_y <- c(0, 1, 3, 7, 10)
  binomial_units <- c(1, 4, 6, 9, 10)
  binomial <- cpp_canonical_likelihood(
    eta, binomial_y, binomial_units, "binomial"
  )
  probability <- plogis(eta)
  expected_binomial_log_likelihood <- rowSums(
    sweep(eta, 2, binomial_y, `*`) -
      sweep(log1p(exp(-abs(eta))) + pmax(eta, 0), 2,
            binomial_units, `*`)
  )
  expected_binomial_score <- sweep(probability, 2, binomial_units, `*`)
  expected_binomial_score <- sweep(
    -expected_binomial_score, 2, binomial_y, `+`
  )
  expected_binomial_curvature <- sweep(
    probability * (1 - probability), 2, binomial_units, `*`
  )

  expect_equal(
    binomial$log_likelihood, expected_binomial_log_likelihood,
    tolerance = 1e-12
  )
  expect_equal(binomial$score, expected_binomial_score, tolerance = 1e-12)
  expect_equal(
    binomial$curvature, expected_binomial_curvature, tolerance = 1e-12
  )

  poisson_y <- c(0, 1, 2, 5, 9)
  exposure <- c(0.5, 1, 1.5, 2, 3)
  poisson_eta <- eta / 3
  poisson <- cpp_canonical_likelihood(
    poisson_eta, poisson_y, exposure, "poisson"
  )
  poisson_mean <- exp(poisson_eta)
  expected_poisson_log_likelihood <- rowSums(
    sweep(poisson_eta, 2, poisson_y, `*`) -
      sweep(poisson_mean, 2, exposure, `*`)
  )
  expected_poisson_score <- sweep(poisson_mean, 2, exposure, `*`)
  expected_poisson_score <- sweep(-expected_poisson_score, 2, poisson_y, `+`)
  expected_poisson_curvature <- sweep(poisson_mean, 2, exposure, `*`)

  expect_equal(
    poisson$log_likelihood, expected_poisson_log_likelihood,
    tolerance = 1e-12
  )
  expect_equal(poisson$score, expected_poisson_score, tolerance = 1e-12)
  expect_equal(
    poisson$curvature, expected_poisson_curvature, tolerance = 1e-12
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
  expect_identical(compiled$tuning_par, r_reference$tuning_par)
  # The two paths use algebraically equivalent derivative calculations.
  expect_equal(compiled$samples$S, r_reference$samples$S, tolerance = 5e-5)
})

test_that("compiled likelihood validates its observation contract", {
  expect_error(
    cpp_canonical_likelihood(matrix(0, 2, 2), 1, c(1, 1), "binomial"),
    "inconsistent lengths"
  )
  expect_error(
    cpp_canonical_likelihood(matrix(0, 1, 1), 0, 1, "gaussian"),
    "Unsupported family"
  )
})
