test_that("Gaussian fitting reports the exact marginal likelihood", {
  fit <- glgpm(
    y ~ cov + gp(), data = gaussian_data, family = "gaussian",
    fix_var_me = 0.2, messages = FALSE
  )

  sigma2 <- exp(fit$estimate$sigma2)
  phi <- exp(fit$estimate$phi)
  covariance <- sigma2 * matern_correlation(
    pairwise_distances(fit$coords), phi = phi, kappa = fit$kappa,
    return_sym_matrix = TRUE
  ) + diag(length(fit$y)) * 0.2
  residual <- fit$y - as.numeric(fit$D %*% fit$estimate$beta)
  root <- chol(covariance)
  exact_log_likelihood <- -0.5 * (
    2 * sum(log(diag(root))) +
      sum(residual * backsolve(root, forwardsolve(t(root), residual)))
  )

  expect_equal(fit$log_lik, exact_log_likelihood, tolerance = 1e-9)
})

test_that("independent MCMC seeds agree within estimated Monte Carlo error", {
  cases <- list(
    binomial = list(y = 4, units_m = 10, sigma2 = 0.9, h = 0.8),
    poisson = list(y = 5, units_m = 1.4, sigma2 = 0.8, h = 0.75)
  )

  for (family in names(cases)) {
    case <- cases[[family]]
    draws <- lapply(c(2718, 3141, 1618), function(seed) {
      fit <- laplace_sampling_mcmc(
        y = case$y, units_m = case$units_m, mu = 0,
        Sigma = matrix(case$sigma2), ID_coords = 1L, family = family,
        control_mcmc = set_control_mcmc(
          n_sim = 3500, burnin = 500, thin = 1, h = case$h, seed = seed
        ),
        messages = FALSE
      )
      as.numeric(fit$samples$S)
    })
    estimates <- vapply(draws, mean, numeric(1))
    mcse <- vapply(draws, function(x) {
      sqrt(var(x) / max(1, as.numeric(sns::ess(x))))
    }, numeric(1))

    differences <- abs(outer(estimates, estimates, "-"))
    joint_mcse <- sqrt(outer(mcse^2, mcse^2, "+"))
    expect_lt(max(differences[upper.tri(differences)] /
                    joint_mcse[upper.tri(joint_mcse)]), 4.5)
  }
})

test_that("joint Gaussian prediction draws reproduce target moments", {
  mean <- c(-0.4, 0.7, 1.2)
  covariance <- matrix(
    c(1.0, 0.3, 0.1,
      0.3, 0.8, 0.25,
      0.1, 0.25, 0.6),
    3, 3, byrow = TRUE
  )
  n_samples <- 20000
  set.seed(8128)
  draws <- sample_correlated_gaussian(
    mean = matrix(mean, nrow = 3, ncol = n_samples),
    lower_root = t(chol(covariance)),
    n_samples = n_samples
  )

  mean_mcse <- sqrt(diag(covariance) / n_samples)
  expect_true(all(abs(rowMeans(draws) - mean) <= 4 * mean_mcse))

  empirical_covariance <- cov(t(draws))
  covariance_mcse <- outer(seq_len(3), seq_len(3), Vectorize(function(i, j) {
    sqrt((covariance[i, i] * covariance[j, j] + covariance[i, j]^2) /
           n_samples)
  }))
  expect_true(all(abs(empirical_covariance - covariance) <=
                    4 * covariance_mcse))
})
