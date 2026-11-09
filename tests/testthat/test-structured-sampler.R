# A survey of 12 locations with households and a smaller second re() term,
# with the household term given first so that the sampler reorders the latent
# effects (its factor puts the term with most levels last).
structured_survey <- function(seed = 7) {
  set.seed(seed)
  sites <- 12L
  coords <- cbind(runif(sites), runif(sites)) * 30
  Sigma <- exp(-as.matrix(dist(coords)) / 8) + diag(0.2, sites)
  location <- rep(seq_len(sites), 1L + rpois(sites, 2))
  records <- length(location)
  ID_re <- cbind(household = seq_len(records),
                 village = (location - 1L) %/% 4L + 1L)
  units <- 1 + rpois(records, 3)
  mu <- rep(-0.8, records)
  y <- rbinom(records, units, plogis(mu + rnorm(records)))
  list(Sigma = Sigma, location = location, ID_re = ID_re, units = units,
       mu = mu, y = y, sigma2_re = c(household = 0.6, village = 0.4))
}

# The Laplace precision Q = P + A' W A in the sampler's order.
explicit_precision <- function(structure, weights) {
  design <- matrix(0, nrow(structure$columns), structure$dimension)
  for (j in seq_len(ncol(structure$columns))) {
    design[cbind(seq_len(nrow(design)), structure$columns[, j] + 1L)] <- 1
  }
  prior <- matrix(0, structure$dimension, structure$dimension)
  prior[seq_len(structure$n_u), seq_len(structure$n_u)] <- structure$prior_u
  diagonal <- structure$n_u + seq_len(structure$n_d)
  prior[cbind(diagonal, diagonal)] <- structure$prior_d
  prior + crossprod(design, weights * design)
}

test_that("the structured factor is a square root of the Laplace covariance", {
  survey <- structured_survey()
  structure <- sampler_latent_structure(
    survey$Sigma, survey$location, survey$ID_re, survey$sigma2_re
  )
  # The household term has most levels and goes last.
  expect_identical(structure$n_d, nrow(survey$ID_re))
  expect_identical(structure$n_u, nrow(survey$Sigma) + max(survey$ID_re[, 2]))
  set.seed(1)
  # Some negative weights, as a link that is not log-concave can give, with
  # Q still positive definite.
  weights <- runif(length(survey$y), -0.05, 2)
  Q <- explicit_precision(structure, weights)
  factor <- cpp_structured_factor_matrix(
    structure$columns, structure$n_u, structure$n_d, structure$prior_u,
    structure$prior_d, weights
  )
  expect_equal(tcrossprod(factor), solve(Q), tolerance = 1e-10)
  right_hand_side <- rnorm(structure$dimension)
  expect_equal(
    cpp_structured_solve(structure$columns, structure$n_u, structure$n_d,
                         structure$prior_u, structure$prior_d, weights,
                         right_hand_side),
    solve(Q, right_hand_side), tolerance = 1e-10
  )
  # Without re() terms the factor is that of the locations alone.
  plain <- sampler_latent_structure(survey$Sigma, seq_len(12L), NULL, NULL)
  plain_weights <- runif(12L, 0.1, 1)
  expect_equal(
    tcrossprod(cpp_structured_factor_matrix(
      plain$columns, plain$n_u, plain$n_d, plain$prior_u, plain$prior_d,
      plain_weights
    )),
    solve(explicit_precision(plain, plain_weights)), tolerance = 1e-10
  )
})

test_that("a Laplace precision that is not positive definite is refused", {
  survey <- structured_survey()
  structure <- sampler_latent_structure(
    survey$Sigma, survey$location, survey$ID_re, survey$sigma2_re
  )
  weights <- rep(-10, length(survey$y))
  expect_null(cpp_structured_solve(
    structure$columns, structure$n_u, structure$n_d, structure$prior_u,
    structure$prior_d, weights, numeric(structure$dimension)
  ))
})

test_that("the sampler finds the mode that maxim_integrand() finds", {
  survey <- structured_survey()
  for (family in c("binomial", "poisson")) {
    y <- if (family == "binomial") survey$y else survey$y + 1
    fit <- laplace_sampling_mcmc(
      y = y, units_m = survey$units, mu = survey$mu, Sigma = survey$Sigma,
      ID_coords = survey$location, ID_re = survey$ID_re,
      sigma2_re = survey$sigma2_re, family = family,
      control_mcmc = set_control_mcmc(n_sim = 30, burnin = 10, thin = 1),
      seed = 3
    )
    reference <- maxim_integrand(
      y = y, denominator = survey$units, mu = survey$mu,
      Sigma = survey$Sigma, ID_coords = survey$location,
      ID_re = survey$ID_re, sigma2_re = survey$sigma2_re, family = family
    )
    expect_equal(fit$mode, reference$mode, tolerance = 1e-6)
    expect_named(fit$samples, c("S", "household", "village"))
    expect_identical(dim(fit$samples$household), c(20L, nrow(survey$ID_re)))
    expect_identical(dim(fit$samples$village),
                     c(20L, max(survey$ID_re[, 2])))
  }
})

test_that("random-effect samples reproduce grid-integrated moments", {
  # One location with two households: latent effects S, d1 and d2.
  y <- c(3, 1)
  m <- c(6, 5)
  sigma2 <- 0.8
  sigma2_household <- 0.5
  fit <- laplace_sampling_mcmc(
    y = y, units_m = m, mu = c(0, 0), Sigma = matrix(sigma2),
    ID_coords = c(1L, 1L), ID_re = cbind(household = 1:2),
    sigma2_re = c(household = sigma2_household), family = "binomial",
    control_mcmc = set_control_mcmc(n_sim = 42000, burnin = 2000, thin = 2),
    seed = 2024
  )
  draws <- cbind(fit$samples$S, fit$samples$household)

  axis <- seq(-4.5, 4.5, length.out = 71)
  grid <- as.matrix(expand.grid(s = axis, d1 = axis, d2 = axis))
  log_weights <- dbinom(y[1], m[1], plogis(grid[, 1] + grid[, 2]), log = TRUE) +
    dbinom(y[2], m[2], plogis(grid[, 1] + grid[, 3]), log = TRUE) -
    grid[, 1]^2 / (2 * sigma2) -
    (grid[, 2]^2 + grid[, 3]^2) / (2 * sigma2_household)
  weights <- exp(log_weights - max(log_weights))
  weights <- weights / sum(weights)
  exact_mean <- colSums(grid * weights)
  centred <- sweep(grid, 2, exact_mean)
  exact_cov <- crossprod(centred, centred * weights)

  expect_equal(unname(colMeans(draws)), unname(exact_mean), tolerance = 0.06)
  expect_equal(unname(cov(draws)), unname(exact_cov), tolerance = 0.07)
})

test_that("a link that is not log-concave is sampled correctly", {
  # An imperfect test, p = 0.01 + 0.89 expit(eta), whose log-likelihood is
  # convex in places.
  inverse <- function(eta) 0.01 + 0.89 * plogis(eta)
  first <- function(eta) 0.89 * dlogis(eta)
  second <- function(eta) 0.89 * dlogis(eta) * (1 - 2 * plogis(eta))
  y <- 8
  m <- 20
  sigma2 <- 4
  fit <- laplace_sampling_mcmc(
    y = y, units_m = m, mu = 0, Sigma = matrix(sigma2), ID_coords = 1L,
    family = "binomial",
    invlink = list(inv = inverse, d1 = first, d2 = second),
    control_mcmc = set_control_mcmc(n_sim = 22000, burnin = 2000, thin = 1),
    seed = 99
  )
  axis <- seq(-15, 15, length.out = 6001)
  log_weights <- dbinom(y, m, inverse(axis), log = TRUE) - axis^2 / (2 * sigma2)
  weights <- exp(log_weights - max(log_weights))
  weights <- weights / sum(weights)
  exact_mean <- sum(axis * weights)
  exact_variance <- sum((axis - exact_mean)^2 * weights)

  draws <- as.numeric(fit$samples$S)
  expect_equal(mean(draws), exact_mean, tolerance = 0.06)
  expect_equal(var(draws), exact_variance, tolerance = 0.08)
})

test_that("averaging the step changes the chain only from the end of burn-in", {
  run <- function(average_step) {
    laplace_sampling_mcmc(
      y = c(2, 7), units_m = c(8, 9), mu = c(0, 0),
      Sigma = matrix(c(0.9, 0.45, 0.45, 0.8), 2), ID_coords = 1:2,
      family = "binomial",
      control_mcmc = set_control_mcmc(n_sim = 400, burnin = 200, thin = 1),
      seed = 5, average_step = average_step
    )
  }
  averaged <- run(TRUE)
  last <- run(FALSE)
  expect_identical(averaged$tuning_par[1:199], last$tuning_par[1:199])
  expect_identical(averaged$acceptance[1:200], last$acceptance[1:200])
  # The kept step lies within the range of the second half of burn-in.
  second_half <- averaged$tuning_par[101:199]
  kept <- averaged$tuning_par[200]
  expect_true(kept >= min(second_half) / 1.02 && kept <= max(second_half) * 1.02)
  expect_true(all(averaged$tuning_par[201:400] == kept))
})

test_that("a chain can start its mode search from a given point", {
  survey <- structured_survey()
  arguments <- list(
    y = survey$y, units_m = survey$units, mu = survey$mu,
    Sigma = survey$Sigma, ID_coords = survey$location, ID_re = survey$ID_re,
    sigma2_re = survey$sigma2_re, family = "binomial",
    control_mcmc = set_control_mcmc(n_sim = 30, burnin = 10, thin = 1),
    seed = 3
  )
  first <- do.call(laplace_sampling_mcmc, arguments)
  restarted <- do.call(
    laplace_sampling_mcmc, c(arguments, list(start = first$mode))
  )
  expect_equal(restarted$mode, first$mode, tolerance = 1e-8)
  expect_error(
    do.call(laplace_sampling_mcmc, c(arguments, list(start = 1:3))),
    "'start' must contain"
  )
})

test_that("a failed mode search falls back to the search of maxim_integrand()", {
  survey <- structured_survey()
  arguments <- list(
    y = survey$y, units_m = survey$units, mu = survey$mu,
    Sigma = survey$Sigma, ID_coords = survey$location, ID_re = survey$ID_re,
    sigma2_re = survey$sigma2_re, family = "binomial",
    control_mcmc = set_control_mcmc(n_sim = 30, burnin = 10, thin = 1),
    seed = 3
  )
  usual <- do.call(laplace_sampling_mcmc, arguments)
  local_mocked_bindings(
    find_sampler_mode = function(...) stop("no mode")
  )
  fallback <- do.call(laplace_sampling_mcmc, arguments)
  expect_equal(fallback$mode, usual$mode, tolerance = 1e-6)
  expect_true(all(is.finite(fallback$samples$S)))
})

test_that("a chain started from a draw of the Laplace approximation needs no burn-in", {
  # About 300 latent effects. From the mode a chain at its default step
  # accepts almost nothing at first, and from a draw of the Laplace
  # approximation it accepts at about the target rate at once.
  set.seed(11)
  sites <- 60L
  coords <- cbind(runif(sites), runif(sites)) * 60
  Sigma <- exp(-as.matrix(dist(coords)) / 10) + diag(0.1, sites)
  location <- rep(seq_len(sites), 1L + rpois(sites, 3))
  records <- length(location)
  units <- 1 + rpois(records, 3)
  y <- rbinom(records, units, 0.2)
  run <- function(start_from_draw) {
    laplace_sampling_mcmc(
      y = y, units_m = units, mu = rep(qlogis(0.2), records), Sigma = Sigma,
      ID_coords = location, ID_re = cbind(household = seq_len(records)),
      sigma2_re = c(household = 0.5), family = "binomial",
      control_mcmc = set_control_mcmc(n_sim = 300, burnin = 0, thin = 1),
      seed = 21, start_from_draw = start_from_draw
    )
  }
  from_draw <- run(TRUE)$acceptance[1:50]
  from_mode <- run(FALSE)$acceptance[1:50]
  expect_gt(mean(from_draw), 0.35)
  expect_lt(mean(from_draw), 0.85)
  expect_lt(mean(from_mode), 0.15)
})

test_that("without adaptation the step stays as given", {
  fit <- laplace_sampling_mcmc(
    y = c(2, 7), units_m = c(8, 9), mu = c(0, 0),
    Sigma = matrix(c(0.9, 0.45, 0.45, 0.8), 2), ID_coords = 1:2,
    family = "binomial",
    control_mcmc = set_control_mcmc(n_sim = 200, burnin = 100, thin = 1,
                                    h = 0.7),
    seed = 5, adapt_step = FALSE
  )
  # exp(log(0.7)) can differ from 0.7 in the last binary digit.
  expect_equal(fit$tuning_par, rep(0.7, 200), tolerance = 1e-14)
  expect_length(fit$samples$S[, 1], 100L)
})
