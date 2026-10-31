test_that("MCMC diagnostics reject fits without samples", {
  expect_error(plot_mcmc("not risk"), "'object' must be a 'RiskMap'")
  expect_error(summarise_mcmc("not risk"), "'object' must be a 'RiskMap'")
  expect_error(plot_mcmc(gaussian_model), "'object' is a Gaussian model")

  binomial_without_samples <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    control_mcml = fixture_control_mcml,
    return_samples = FALSE, messages = FALSE
  ))
  expect_null(binomial_without_samples$S_samples)
  expect_null(binomial_without_samples$linear_predictors)
  expect_null(fitted(binomial_without_samples))
  expect_error(plot_mcmc(binomial_without_samples),
               "'object' does not contain any MCMC chains")
  expect_error(summarise_mcmc(binomial_without_samples),
               "'object' does not contain any MCMC chains")
})

test_that("plot_mcmc checks its arguments", {
  # The model has 10 locations followed by 5 levels of re(i), whose samples
  # are not locations.
  for (component in list(11, 20, 10.5, 0, "no", c(1, 1), Inf, integer())) {
    expect_error(plot_mcmc(poisson_model, check_mean = FALSE,
                           component = component),
                 "'component' must contain distinct whole numbers between 1 and the number of locations, 10")
  }
  expect_error(plot_mcmc(poisson_model, check_mean = NA),
               "'check_mean' must be either TRUE or FALSE")
})

test_that("plot_mcmc shows the mean and the location with the smallest ESS", {
  diagnostics <- plot_mcmc(poisson_model)
  traces <- diagnostics[[1L]][[1L]]
  correlations <- diagnostics[[1L]][[2L]]
  histogram <- diagnostics[[2L]]
  spatial <- poisson_model$S_samples[, 1:10]
  smallest <- which.min(mcmc_effective_size(spatial))

  expect_s3_class(diagnostics, "patchwork")
  expect_identical(levels(traces$data$chain), c(
    "Mean of S over 10 locations",
    sprintf("S at location %d, the smallest ESS", smallest)
  ))
  expect_equal(traces$data$value,
               c(rowMeans(spatial), spatial[, smallest]))
  expect_identical(nlevels(correlations$data$chain), 2L)
  expect_equal(histogram$data$ess, mcmc_effective_size(spatial))
  expect_no_error(print(diagnostics))

  correlations <- plot_mcmc(poisson_model, lag.max = 5)[[1L]][[2L]]
  expect_identical(unique(correlations$data$lag), as.numeric(0:5))
})

test_that("plot_mcmc keeps the book's call for one chosen location", {
  diagnostics <- plot_mcmc(poisson_model, check_mean = FALSE, component = 10)
  traces <- diagnostics[[1L]][[1L]]

  expect_identical(levels(traces$data$chain), "S at location 10")
  expect_equal(traces$data$value, poisson_model$S_samples[, 10])
  expect_no_error(print(diagnostics))
  expect_identical(
    levels(plot_mcmc(poisson_model, component = c(3, 1))[[1L]][[1L]]$data$chain),
    c("Mean of S over 10 locations", "S at location 3", "S at location 1")
  )
})

test_that("traces are drawn against the iterations at which samples were kept", {
  expect_identical(attr(poisson_model, "mcmc")[c("n_sim", "burnin", "thin")],
                   list(n_sim = 1000, burnin = 200, thin = 10))
  traces <- plot_mcmc(poisson_model)[[1L]][[1L]]
  expect_identical(unique(traces$data$iteration), 200 + 10 * seq_len(80))
  expect_identical(ggplot2::get_labs(traces)$x, "Iteration")

  # Fits saved before the settings were recorded number the samples instead.
  older_fit <- poisson_model
  attr(older_fit, "mcmc")[c("n_sim", "burnin", "thin")] <- NULL
  traces <- plot_mcmc(older_fit)[[1L]][[1L]]
  expect_identical(unique(traces$data$iteration), seq_len(80))
  expect_identical(ggplot2::get_labs(traces)$x, "Sample")
})

test_that("summarise_mcmc returns the ESS of the mean and of every component", {
  summary <- summarise_mcmc(poisson_model)
  samples <- poisson_model$S_samples

  expect_s3_class(summary, "data.frame")
  expect_named(summary, c("effect", "component", "ess"))
  expect_identical(summary$effect,
                   c("mean of S", rep("S", 10), rep("i", 5)))
  expect_identical(summary$component, c(NA, 1:10, 1:5))
  expect_identical(attr(summary, "n_samples"), 80L)
  expect_equal(summary$ess, c(
    mcmc_effective_size(rowMeans(samples[, 1:10])),
    mcmc_effective_size(samples)
  ))
})

test_that("the effective sample size matches known autocorrelation times", {
  set.seed(2026)
  independent <- rnorm(20000)
  autoregressive <- as.numeric(stats::arima.sim(list(ar = 0.5), 20000))
  slow <- as.numeric(stats::arima.sim(list(ar = 0.9), 20000))

  # An AR(1) chain with coefficient rho has ESS n (1 - rho) / (1 + rho). The
  # tolerances are about three standard deviations of the estimate, from 300
  # simulated chains of each kind.
  expect_equal(mcmc_effective_size(independent), 20000, tolerance = 0.07)
  expect_equal(mcmc_effective_size(autoregressive), 20000 / 3,
               tolerance = 0.12)
  expect_equal(mcmc_effective_size(slow), 20000 / 19, tolerance = 0.25)
  expect_equal(
    mcmc_effective_size(cbind(independent, autoregressive, slow),
                        block_size = 2L),
    mcmc_effective_size(cbind(independent, autoregressive, slow))
  )
})

test_that("the effective sample size uses Geyer's estimator", {
  set.seed(1)
  chain <- as.numeric(stats::arima.sim(list(ar = 0.7), 300))
  # The Fourier transform gives the same autocovariances as acf()
  autocovariance <- acf(chain, type = "covariance", lag.max = 299,
                        plot = FALSE)$acf
  expect_equal(mcmc_effective_size(chain),
               initial_monotone_ess(as.numeric(autocovariance)))

  # Autocorrelations 0.5, 0.3, 0.2, 0.1 and then zero give pair sums 1.5,
  # 0.5 and 0.1, which are already monotone, so the autocorrelation time is
  # 2 * 2.1 - 1 = 3.2.
  expect_equal(initial_monotone_ess(c(1, 0.5, 0.3, 0.2, 0.1, 0, 0, 0)),
               8 / 3.2)
  # A rise in the pair sums is cut back to the previous sum.
  expect_equal(initial_monotone_ess(c(1, 0.2, 0.1, 0.1, 0.3, 0.1, 0, 0)),
               8 / (2 * (1.2 + 0.2 + 0.2) - 1))

  # A chain that never moves has no effective samples, and the
  # autocorrelation time of an antithetic chain is at least 1 / log10(n).
  expect_identical(mcmc_effective_size(rep(1.5, 100)), 0)
  expect_equal(mcmc_effective_size(rep(c(-1, 1), 500)), 3000)
})
