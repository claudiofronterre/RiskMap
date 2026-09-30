test_that("MALA adaptation stops after burn-in", {
  control <- set_control_mcmc(
    n_sim = 80, burnin = 40, thin = 2, h = 0.5,
    seed = 918
  )

  fit <- laplace_sampling_mcmc(
    y = 3, units_m = 8, mu = 0, Sigma = matrix(0.7),
    ID_coords = 1L, family = "binomial", control_mcmc = control,
    messages = FALSE
  )

  expect_equal(nrow(fit$samples$S), 20)
  expect_true(all(fit$tuning_par[41:80] == fit$tuning_par[40]))
  expect_length(fit$acceptance, 80)
  expect_named(fit$acceptance_rate, c("burnin", "sampling"))
})

test_that("one-dimensional binomial samples reproduce exact posterior moments", {
  control <- set_control_mcmc(
    n_sim = 7000, burnin = 1000, thin = 1, h = 0.8,
    seed = 2718
  )
  sigma2 <- 0.9
  y <- 4
  m <- 10

  fit <- laplace_sampling_mcmc(
    y = y, units_m = m, mu = 0, Sigma = matrix(sigma2),
    ID_coords = 1L, family = "binomial", control_mcmc = control,
    messages = FALSE
  )

  log_kernel <- function(s) {
    dbinom(y, size = m, prob = plogis(s), log = TRUE) +
      dnorm(s, sd = sqrt(sigma2), log = TRUE)
  }
  normalizer <- integrate(function(s) exp(log_kernel(s)), -Inf, Inf)$value
  exact_mean <- integrate(
    function(s) s * exp(log_kernel(s)), -Inf, Inf
  )$value / normalizer
  exact_second <- integrate(
    function(s) s^2 * exp(log_kernel(s)), -Inf, Inf
  )$value / normalizer

  expect_equal(mean(fit$samples$S), exact_mean, tolerance = 0.06)
  expect_equal(var(as.numeric(fit$samples$S)),
               exact_second - exact_mean^2, tolerance = 0.06)
})

test_that("sampler controls and location indices are validated", {
  expect_error(
    set_control_mcmc(n_sim = 100, burnin = 100),
    "larger than burnin"
  )
  expect_error(set_control_mcmc(h = 0), "positive finite")
  expect_error(set_control_mcmc(c2.h = 0.5), "larger than 0.5")

  control <- set_control_mcmc(n_sim = 20, burnin = 10, thin = 1)
  expect_error(
    laplace_sampling_mcmc(
      y = c(1, 2), units_m = c(4, 4), mu = 0,
      Sigma = matrix(1), ID_coords = c(1, 2), family = "binomial",
      control_mcmc = control, messages = FALSE
    ),
    "valid integer location index"
  )
})
