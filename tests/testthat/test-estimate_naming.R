## Validation for issue #92: `estimate` is a named list (beta, sigma2, phi,
## [nu2], [sigma2_me], [sigma2_re]) built directly by the fitting engines
## (glgpm_lm()/glgpm_nong(), via structure_estimate()), and coef.RiskMap()/
## summary.RiskMap() read from it directly rather than recomputing positions
## themselves.
##
## A flat named vector can't be keyed safely by name: a covariate literally
## named e.g. "sigma2" collides with the spatial variance parameter's own
## name, silently corrupting name-based lookups (reported against the
## flat-vector version of this fix - see the last test below). Nesting
## `estimate` keeps those namespaces separate; `unlist(estimate)` (used
## internally by summary.RiskMap() to line up against `covariance`, which is
## estimated as one joint matrix over every parameter) disambiguates the same
## way, e.g. "beta.sigma2" vs "sigma2".

test_that("estimate is a list with numeric beta/sigma2/phi for every model configuration", {
  fits <- list(gaussian_model, gaussian_offset_model, gaussian_intercept_model,
               binomial_model, poisson_model)

  for (fit in fits) {
    expect_type(fit$estimate, "list")
    expect_true(is.numeric(fit$estimate$beta))
    expect_true(is.numeric(fit$estimate$sigma2))
    expect_true(is.numeric(fit$estimate$phi))
  }
})

test_that("covariance is named to match unlist(estimate) on both dimensions", {
  fits <- list(gaussian_model, gaussian_offset_model, gaussian_intercept_model,
               binomial_model, poisson_model)

  for (fit in fits) {
    flat_names <- names(unlist(fit$estimate))
    expect_equal(dimnames(fit$covariance), list(flat_names, flat_names))
  }
})

test_that("estimate reproduces coef.RiskMap() for a model with random effects (no nugget)", {
  co <- coef(gaussian_model)
  est <- gaussian_model$estimate

  expect_setequal(names(est), c("beta", "sigma2", "phi", "sigma2_me", "sigma2_re"))
  expect_equal(est$beta, co$beta, ignore_attr = TRUE)
  expect_equal(exp(est$sigma2), co$sigma2, ignore_attr = TRUE)
  expect_equal(exp(est$phi), co$phi, ignore_attr = TRUE)
  expect_equal(exp(est$sigma2_me), co$sigma2_me, ignore_attr = TRUE)
  expect_equal(exp(est$sigma2_re), co$sigma2_re, ignore_attr = TRUE)
  expect_equal(names(co$sigma2_re), "i")
})

test_that("estimate reproduces coef.RiskMap() for an intercept-only model", {
  co <- coef(gaussian_intercept_model)
  est <- gaussian_intercept_model$estimate

  expect_setequal(names(est), c("beta", "sigma2", "phi", "sigma2_me"))
  expect_equal(est$beta, co$beta, ignore_attr = TRUE)
  expect_equal(exp(est$sigma2), co$sigma2, ignore_attr = TRUE)
  expect_equal(exp(est$phi), co$phi, ignore_attr = TRUE)
  expect_equal(exp(est$sigma2_me), co$sigma2_me, ignore_attr = TRUE)
})

test_that("coef.RiskMap() does not mislabel a single covariate as the intercept when the model has none (#92)", {
  ## Bug: coef.RiskMap() used to relabel beta to "Intercept" whenever there
  ## was exactly one coefficient, regardless of whether it actually was the
  ## intercept - so `y ~ 0 + cov + gp()` (no intercept, one covariate)
  ## reported "Intercept" instead of "cov".
  fit <- glgpm(y ~ 0 + cov + gp(), data = gaussian_data, family = "gaussian",
              messages = FALSE)

  expect_equal(colnames(fit$D), "cov")
  expect_equal(names(fit$estimate$beta), "cov")
  expect_equal(names(coef(fit)$beta), "cov")
  expect_equal(rownames(summary(fit)$reg_coef), "cov")
})

test_that("estimate reproduces coef.RiskMap() for a model with an estimated nugget and fix_var_me", {
  ## gaussian_offset_model: gp(nugget = TRUE), fix_var_me = 0, no random effects
  co <- coef(gaussian_offset_model)
  est <- gaussian_offset_model$estimate

  expect_setequal(names(est), c("beta", "sigma2", "phi", "nu2"))
  ## tau2 isn't a raw parameter: the engines fit nu2 = tau2 / sigma2 on the log
  ## scale, so reconstructing tau2 needs both entries added before exp()
  expect_equal(exp(est$nu2 + est$sigma2), co$tau2, ignore_attr = TRUE)
})

test_that("estimate reproduces coef.RiskMap() for binomial and poisson models with random effects", {
  for (fit in list(binomial_model, poisson_model)) {
    co <- coef(fit)
    est <- fit$estimate

    expect_setequal(names(est), c("beta", "sigma2", "phi", "sigma2_re"))
    expect_equal(est$beta, co$beta, ignore_attr = TRUE)
    expect_equal(exp(est$sigma2), co$sigma2, ignore_attr = TRUE)
    expect_equal(exp(est$phi), co$phi, ignore_attr = TRUE)
    expect_equal(exp(est$sigma2_re), co$sigma2_re, ignore_attr = TRUE)
    expect_equal(names(co$sigma2_re), "i")
  }
})

test_that("estimate correctly names sigma2_me even when the nugget is also estimated", {
  ## Before the #92 name-based swap-over, coef.RiskMap() never assigned the
  ## name "sigma2_me" internally in this combination (nugget estimated +
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

  expect_setequal(names(est), c("beta", "sigma2", "phi", "nu2", "sigma2_me"))
  expect_equal(exp(est$sigma2_me), co$sigma2_me, ignore_attr = TRUE)
  expect_equal(exp(est$nu2 + est$sigma2), co$tau2, ignore_attr = TRUE)
})

test_that("a covariate named after a parameter no longer corrupts estimates (#92)", {
  ## Regression test for the bug a colleague reported: with `estimate` as one
  ## flat named vector, a covariate literally called "sigma2" collided with
  ## the spatial variance parameter's own name, and name-based lookups (e.g.
  ## estimate["sigma2"] in coef.RiskMap()) silently returned the wrong value
  ## - the covariate's beta coefficient, exponentiated, instead of the actual
  ## spatial variance. Structuring `estimate` as a list (beta$sigma2 vs
  ## top-level $sigma2) keeps those namespaces separate.
  clashing_data <- gaussian_data
  clashing_data$sigma2 <- clashing_data$cov

  fit_clash <- glgpm(y ~ sigma2 + gp() + re(i), data = clashing_data,
                     family = "gaussian", messages = FALSE)
  fit_ref   <- glgpm(y ~ cov + gp() + re(i), data = gaussian_data,
                     family = "gaussian", messages = FALSE)

  co_clash <- coef(fit_clash)
  co_ref   <- coef(fit_ref)

  expect_equal(co_clash$sigma2, co_ref$sigma2)
  expect_equal(co_clash$phi, co_ref$phi)
  expect_equal(co_clash$beta, co_ref$beta, ignore_attr = TRUE)
  expect_equal(co_clash$sigma2_re, co_ref$sigma2_re, ignore_attr = TRUE)

  ## summary()'s standard errors go through the separate unlist()-based
  ## delta-method path, so check that's unaffected by the name clash too
  expect_equal(summary(fit_clash)$sp, summary(fit_ref)$sp, ignore_attr = TRUE)
})
