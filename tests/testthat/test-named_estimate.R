## Validation for issue #92: `named_estimate` is a named copy of `estimate`,
## added alongside the untouched, unnamed `estimate` vector. These tests
## establish that the names attached in glgpm_lm()/glgpm_nong() correctly
## identify each raw (working-scale) parameter, by checking that combining
## them with the appropriate exp()/addition reproduces exactly what
## coef.RiskMap()/summary.RiskMap() already report today via their
## independent, position-based logic. Nothing about coef.RiskMap() or
## summary.RiskMap() has changed yet; this is groundwork ahead of switching
## them over to look estimates up by name.

test_that("named_estimate has the same values as the untouched estimate vector", {
  fits <- list(gaussian_model, gaussian_offset_model, gaussian_intercept_model,
               binomial_model, poisson_model)

  for (fit in fits) {
    expect_false(is.null(names(fit$named_estimate)))
    expect_equal(length(fit$named_estimate), length(fit$estimate))
    ## `estimate` itself is untouched (still whatever incidental/partial names
    ## nlminb happened to carry through, e.g. from quantile()'s "10%"); only
    ## the values are guaranteed to match
    expect_equal(unname(fit$named_estimate), unname(fit$estimate))
  }
})

test_that("named_estimate reproduces coef.RiskMap() for a model with random effects (no nugget)", {
  co <- coef(gaussian_model)
  ne <- gaussian_model$named_estimate

  expect_setequal(names(ne), c("(Intercept)", "cov", "sigma2", "phi", "sigma2_me", "i_sigma2_re"))
  expect_equal(unname(ne[c("(Intercept)", "cov")]), unname(co$beta))
  expect_equal(exp(unname(ne["sigma2"])), co$sigma2)
  expect_equal(exp(unname(ne["phi"])), co$phi)
  expect_equal(exp(unname(ne["sigma2_me"])), co$sigma2_me)
  expect_equal(exp(unname(ne["i_sigma2_re"])), unname(co$sigma2_re))
})

test_that("named_estimate reproduces coef.RiskMap() for an intercept-only model", {
  co <- coef(gaussian_intercept_model)
  ne <- gaussian_intercept_model$named_estimate

  expect_setequal(names(ne), c("(Intercept)", "sigma2", "phi", "sigma2_me"))
  expect_equal(unname(ne["(Intercept)"]), unname(co$beta))
  expect_equal(exp(unname(ne["sigma2"])), co$sigma2)
  expect_equal(exp(unname(ne["phi"])), co$phi)
  expect_equal(exp(unname(ne["sigma2_me"])), co$sigma2_me)
})

test_that("named_estimate reproduces coef.RiskMap() for a model with an estimated nugget and fix_var_me", {
  ## gaussian_offset_model: gp(nugget = TRUE), fix_var_me = 0, no random effects
  co <- coef(gaussian_offset_model)
  ne <- gaussian_offset_model$named_estimate

  expect_setequal(names(ne), c("(Intercept)", "cov", "sigma2", "phi", "nu2"))
  ## tau2 isn't a raw parameter: the engines fit nu2 = tau2 / sigma2 on the log
  ## scale, so reconstructing tau2 needs both named entries added before exp()
  expect_equal(exp(unname(ne["nu2"]) + unname(ne["sigma2"])), unname(co$tau2))
})

test_that("named_estimate reproduces coef.RiskMap() for binomial and poisson models with random effects", {
  for (fit in list(binomial_model, poisson_model)) {
    co <- coef(fit)
    ne <- fit$named_estimate

    expect_setequal(names(ne), c("(Intercept)", "cov", "sigma2", "phi", "i_sigma2_re"))
    expect_equal(unname(ne[c("(Intercept)", "cov")]), unname(co$beta))
    expect_equal(exp(unname(ne["sigma2"])), co$sigma2)
    expect_equal(exp(unname(ne["phi"])), co$phi)
    expect_equal(exp(unname(ne["i_sigma2_re"])), unname(co$sigma2_re))
  }
})

test_that("named_estimate correctly names sigma2_me even when the nugget is also estimated", {
  ## coef.RiskMap() never assigns the name "sigma2_me" internally in this
  ## combination (nugget estimated + sigma2_me estimated) - see the #92 issue
  ## comment. Its *value* is still extracted correctly by position, so
  ## coef()/summary() output is unaffected; named_estimate names this slot
  ## correctly, unlike coef.RiskMap()'s own internal bookkeeping.
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
  ne <- fit$named_estimate

  expect_setequal(names(ne), c("(Intercept)", "cov", "sigma2", "phi", "nu2", "sigma2_me"))
  expect_equal(exp(unname(ne["sigma2_me"])), co$sigma2_me)
  expect_equal(exp(unname(ne["nu2"]) + unname(ne["sigma2"])), unname(co$tau2))
})
