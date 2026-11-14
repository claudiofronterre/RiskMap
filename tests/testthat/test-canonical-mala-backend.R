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
    cpp_link_log_likelihood(eta, binomial_y, binomial_units, 0L),
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
    cpp_link_log_likelihood(poisson_eta, poisson_y, exposure, 1L),
    expected_poisson, tolerance = 1e-12
  )
})

test_that("compiled MALA supports multiple random-effect blocks", {
  control <- set_control_mcmc(
    n_sim = 100, burnin = 40, thin = 3, h = 0.5
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
    family = "binomial", control_mcmc = control, seed = 4821
  )

  expect_named(fit$samples, c("S", "cluster", "period"))
  expect_identical(dim(fit$samples$S), c(20L, 2L))
  expect_identical(dim(fit$samples$cluster), c(20L, 2L))
  expect_identical(dim(fit$samples$period), c(20L, 2L))
  expect_true(all(is.finite(unlist(fit$samples))))
  expect_true(all(fit$tuning_par[41:100] == fit$tuning_par[40]))
})

test_that("a canonical link given as a function gives the same chain", {
  control <- set_control_mcmc(
    n_sim = 100, burnin = 40, thin = 3, h = 0.5
  )
  arguments <- list(
    y = c(1, 3), units_m = c(5, 6), mu = c(-0.5, 0.2),
    Sigma = matrix(c(0.8, 0.2, 0.2, 0.7), 2), ID_coords = 1:2,
    family = "binomial", control_mcmc = control, seed = 99
  )
  compiled <- do.call(laplace_sampling_mcmc, arguments)
  # The same link given as a function is evaluated in R once per iteration.
  r_reference <- do.call(
    laplace_sampling_mcmc,
    c(arguments, list(invlink = function(x) plogis(x)))
  )

  expect_identical(compiled$acceptance, r_reference$acceptance)
  # Equivalent floating-point update orders differ at machine precision on
  # some BLAS/platform combinations.
  expect_equal(compiled$tuning_par, r_reference$tuning_par,
               tolerance = 1e-14)
  # The two evaluations of the link are algebraically equivalent.
  expect_equal(compiled$samples$S, r_reference$samples$S, tolerance = 5e-5)
})

test_that("compiled log-likelihood guards Poisson overflow", {
  overflow <- matrix(c(0, 800), nrow = 2)
  expect_identical(
    cpp_link_log_likelihood(overflow, 1, 1, 1L)[2], -Inf
  )
})

test_that("final conditional sampling uses the compiled kernel for canonical links", {
  expect_identical(attr(binomial_model, "mcmc")$invlink_used, "canonical")
  expect_identical(attr(poisson_model, "mcmc")$invlink_used, "canonical")
})

test_that("the compiled rows of a named link are R's", {
  # Linear predictors where the probability rounds to 0 or 1 and the bounds of
  # bound_log_probability() act, a missing one that gives -Inf, and records
  # with no events and with every trial an event.
  eta <- rbind(seq(-40, 40, length.out = 6), c(-8, -1, 0, 1, 8, 39),
               c(-0.3, 0.2, NaN, 1.5, -2, 0))
  y <- c(0, 1, 3, 7, 10, 10)
  units <- c(1, 4, 6, 9, 10, 10)
  for (link in c("probit", "cloglog", "cauchit")) {
    linkf <- make_inverse_link_functions(
      "binomial", c(named_inverse_link(link), link = link), 5L
    )
    expected <- compute_link_log_likelihood_rows(eta, y, units, "binomial",
                                                 linkf)
    code <- make_family("binomial", linkf)$code
    expect_identical(code, match(link, c("probit", "cloglog", "cauchit")) + 1L)
    # Equal to R's to the last bit where the compiler does not fuse multiply
    # and add, as on x86-64.
    expect_equal(cpp_link_log_likelihood(eta, y, units, code), expected,
                 tolerance = 1e-14)
    expect_identical(cpp_link_log_likelihood(eta, y, units, code)[3], -Inf)
  }
})
