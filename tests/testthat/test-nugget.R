test_that("nugget and fix_var_me cannot both be estimated when family is gaussian", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(0, 1, 2),
    z = c(0, 1, 2),
    denominator = c(4, 4, 4)
  )

  gaussian_data <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(4326))

  expect_error(glgpm(z ~ gp(nugget = TRUE), gaussian_data, "gaussian", messages = FALSE), "When there is only one observation per location")
  expect_no_error(glgpm(z ~ gp(nugget = TRUE), gaussian_data, "gaussian", fix_var_me = 1, messages = FALSE))
  zero_error_fit <- glgpm(
    z ~ gp(), gaussian_data, "gaussian", fix_var_me = 0, messages = FALSE
  )
  expect_identical(zero_error_fit$fix_var_me, 0)
  expect_lt(max(abs(zero_error_fit$grad_MLE)), 0.001)
  expect_warning(
    expect_no_error(glgpm(
      z ~ gp(nugget = TRUE), gaussian_data, "binomial",
      denominator = denominator, control_mcmc = control_mcmc,
      control_mcml = fixture_control_mcml, messages = FALSE
    )),
    "estimated as effectively zero"
  )

  two_location_data <- rbind(data,
                             data.frame(
                              x = 1,
                              y = 0,
                              z = 3,
                              denominator = 4
                            ))

  gaussian_data <- sf::st_as_sf(two_location_data, coords = c("x", "y"), crs = sf::st_crs(4326))

  expect_error(
    glgpm(z ~ gp(), gaussian_data, "gaussian", fix_var_me = 0,
          messages = FALSE),
    "structurally singular"
  )

  expect_no_error(glgpm(z ~ gp(nugget = TRUE), gaussian_data, "gaussian", messages = FALSE))

  result <- glgpm(z ~ gp(nugget = TRUE), gaussian_data, "gaussian", fix_var_me = 1, messages = FALSE)
  expect_true("tau2" %in% names(coef(result)))

  result <- glgpm(z ~ gp(nugget = 1), gaussian_data, "gaussian", messages = FALSE)
  expect_identical(summary(result)$tau2, 1)
  expect_lt(max(abs(result$grad_MLE)), 0.001)

})

test_that("nugget_ratio distinguishes estimated and fixed nuggets", {
  expect_identical(nugget_ratio(FALSE, sigma2 = 2), 0)
  expect_identical(nugget_ratio(0.6, sigma2 = 2), 0.3)
  expect_equal(nugget_ratio(TRUE, sigma2 = 2, log_nu2 = log(0.3)), 0.3)
  expect_error(
    nugget_ratio(TRUE, sigma2 = 2),
    "'log_nu2' is required when the nugget is estimated"
  )
})

test_that("non-Gaussian estimated nugget uses consistent derivatives", {
  set.seed(2029)
  n <- 80L
  denominator <- sample(10:50, n, replace = TRUE)
  covariate <- rnorm(n)
  x <- runif(n, 0, 100000)
  z <- runif(n, 0, 100000)
  # real spatial and nugget components, so sigma2 and tau2 are identifiable
  sigma2 <- 1; phi <- 20000; kappa <- 1
  Sigma <- sigma2 * matern_correlation(dist(cbind(x, z)), phi = phi, kappa = kappa,
                                       return_sym_matrix = TRUE)
  S <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
  nugget_true <- rnorm(n, sd = 0.6)
  data <- data.frame(
    y = rbinom(n, denominator, plogis(-0.5 + 0.4 * covariate + S + nugget_true)),
    covariate = covariate,
    denominator = denominator,
    x = x,
    z = z
  )
  spatial_data <- sf::st_as_sf(data, coords = c("x", "z"), crs = 32630)
  control <- set_control_mcmc(
    n_sim = 500,
    burnin = 100,
    thin = 4,
    seed = 2029
  )

  fit <- glgpm(
    y ~ covariate + gp(nugget = TRUE),
    data = spatial_data,
    denominator = denominator,
    family = "binomial",
    control_mcmc = control,
    control_mcml = fixture_control_mcml,
    messages = FALSE
  )

  # A large score at the reported optimum catches disagreement between the
  # Monte Carlo likelihood and its analytical derivatives. The 0.02 threshold
  # allows for the finite 500-draw fixture and platform-specific BLAS rounding
  # while remaining small relative to a unit score contribution.
  expect_lt(max(abs(fit$grad_MLE)), 0.02)
  expect_true(is.finite(coef(fit)[["tau2"]]))
  expect_identical(attr(fit, "optimizer")$convergence, 0L)
  stage <- fit$mcml_history[[1L]]
  expect_gte(stage$fit_effective_samples, 1)
  expect_lte(stage$fit_diagnostics$relative_ess, 1)
})
