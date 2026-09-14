## Validation for issue #92: `estimate` is named directly by the fitting
## engines (glgpm_lm()/glgpm_nong(), via name_estimates()), and coef.RiskMap()/
## summary.RiskMap() subset it (and `covariance`, which shares the same
## dimnames) by these names rather than recomputing positions themselves.
## These tests check that the names correctly identify each raw (working-
## scale) parameter, by manually combining them with the appropriate
## exp()/addition and confirming the result matches what coef.RiskMap()
## reports.

test_that("estimate is named for every model configuration", {
  fits <- list(gaussian_model, gaussian_offset_model, gaussian_intercept_model,
               binomial_model, poisson_model)

  for (fit in fits) {
    expect_false(is.null(names(fit$estimate)))
    expect_equal(length(names(fit$estimate)), length(fit$estimate))
  }
})

test_that("covariance shares estimate's names on both dimensions", {
  fits <- list(gaussian_model, gaussian_offset_model, gaussian_intercept_model,
               binomial_model, poisson_model)

  for (fit in fits) {
    expect_equal(dimnames(fit$covariance), list(names(fit$estimate), names(fit$estimate)))
  }
})

test_that("estimate reproduces coef.RiskMap() for a model with random effects (no nugget)", {
  co <- coef(gaussian_model)
  est <- gaussian_model$estimate

  expect_setequal(names(est), c("(Intercept)", "cov", "sigma2", "phi", "sigma2_me", "sigma2_re_i"))
  expect_equal(est[c("(Intercept)", "cov")], co$beta, ignore_attr = TRUE)
  expect_equal(exp(est["sigma2"]), co$sigma2, ignore_attr = TRUE)
  expect_equal(exp(est["phi"]), co$phi, ignore_attr = TRUE)
  expect_equal(exp(est["sigma2_me"]), co$sigma2_me, ignore_attr = TRUE)
  expect_equal(exp(est["sigma2_re_i"]), co$sigma2_re, ignore_attr = TRUE)
  expect_equal(names(co$sigma2_re), "i")
})

test_that("estimate reproduces coef.RiskMap() for an intercept-only model", {
  co <- coef(gaussian_intercept_model)
  est <- gaussian_intercept_model$estimate

  expect_setequal(names(est), c("(Intercept)", "sigma2", "phi", "sigma2_me"))
  expect_equal(est["(Intercept)"], co$beta, ignore_attr = TRUE)
  expect_equal(exp(est["sigma2"]), co$sigma2, ignore_attr = TRUE)
  expect_equal(exp(est["phi"]), co$phi, ignore_attr = TRUE)
  expect_equal(exp(est["sigma2_me"]), co$sigma2_me, ignore_attr = TRUE)
})

test_that("estimate reproduces coef.RiskMap() for a model with an estimated nugget and fix_var_me", {
  ## gaussian_offset_model: gp(nugget = TRUE), fix_var_me = 0, no random effects
  co <- coef(gaussian_offset_model)
  est <- gaussian_offset_model$estimate

  expect_setequal(names(est), c("(Intercept)", "cov", "sigma2", "phi", "nu2"))
  ## tau2 isn't a raw parameter: the engines fit nu2 = tau2 / sigma2 on the log
  ## scale, so reconstructing tau2 needs both named entries added before exp()
  expect_equal(exp(est["nu2"] + est["sigma2"]), co$tau2, ignore_attr = TRUE)
})

test_that("estimate reproduces coef.RiskMap() for binomial and poisson models with random effects", {
  for (fit in list(binomial_model, poisson_model)) {
    co <- coef(fit)
    est <- fit$estimate

    expect_setequal(names(est), c("(Intercept)", "cov", "sigma2", "phi", "sigma2_re_i"))
    expect_equal(est[c("(Intercept)", "cov")], co$beta, ignore_attr = TRUE)
    expect_equal(exp(est["sigma2"]), co$sigma2, ignore_attr = TRUE)
    expect_equal(exp(est["phi"]), co$phi, ignore_attr = TRUE)
    expect_equal(exp(est["sigma2_re_i"]), co$sigma2_re, ignore_attr = TRUE)
    expect_equal(names(co$sigma2_re), "i")
  }
})

test_that("estimate correctly names sigma2_me even when the nugget is also estimated", {
  ## Before the #92 swap-over, coef.RiskMap() never assigned the name
  ## "sigma2_me" internally in this combination (nugget estimated +
  ## sigma2_me estimated) - see the #92 issue comment. Its *value* was still
  ## extracted correctly by position, so output was numerically unaffected;
  ## this guards against that gap reappearing now that lookup is name-based.
  set.seed(2)
  n_loc <- 6
  coords_u <- cbind(runif(n_loc, 0, 10000), runif(n_loc, 0, 10000))
  coords_rep <- coords_u[rep(seq_len(n_loc), each = 2), ]
  d <- data.frame(x = coords_rep[, 1], z = coords_rep[, 2], cov = rnorm(2 * n_loc))
  d$y <- 1 + 0.5 * d$cov + rnorm(2 * n_loc)
  dup_data <- sf::st_as_sf(d, coords = c("x", "z"), crs = 32637)

  fit <- glgpm(y ~ cov + gp(nugget = TRUE), data = dup_data, family = "gaussian",
              messages = FALSE)

  co <- coef(fit)
  est <- fit$estimate

  expect_setequal(names(est), c("(Intercept)", "cov", "sigma2", "phi", "nu2", "sigma2_me"))
  expect_equal(exp(est["sigma2_me"]), co$sigma2_me, ignore_attr = TRUE)
  expect_equal(exp(est["nu2"] + est["sigma2"]), co$tau2, ignore_attr = TRUE)
})
