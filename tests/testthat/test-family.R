## The family object (#253, step 1a). Its methods are the arithmetic the engine
## used before, and its properties the engine's tests of the canonical link.

probit_family_link <- function() {
  make_inverse_link_functions(
    "binomial", c(named_inverse_link("probit"), link = "probit"), 5L
  )
}

test_that("the properties of a family object follow the link", {
  logit <- make_family("binomial")
  expect_true(logit$canonical)
  expect_true(logit$compiled)
  expect_true(logit$assumes_concave)
  expect_false(logit$moves_start)
  expect_false(logit$needs_finite_weight)
  expect_false(logit$poisson)
  expect_identical(logit$linkf$link, "logit")
  expect_true(make_family("poisson")$poisson)

  probit <- make_family("binomial", probit_family_link())
  expect_false(probit$canonical)
  expect_false(probit$compiled)
  expect_false(probit$assumes_concave)
  expect_true(probit$moves_start)
  expect_true(probit$needs_finite_weight)

  sites <- make_site_table(c(0, 20, 3), c(20, 20, 20), "binomial")
  integrated <- make_family("binomial", make_integrated_kernel(sites, 0.5))
  expect_true(integrated$integrated)
  expect_false(integrated$canonical)
  expect_false(integrated$moves_start)
  expect_null(integrated$slope_parts)
})

test_that("a family object is defined only for the binomial and Poisson families", {
  expect_error(make_family("gaussian"), "binomial and Poisson")
})

test_that("the methods of the canonical object agree with those of a custom logit", {
  set.seed(7)
  eta <- matrix(rnorm(12, -0.5, 1.5), 3)
  y <- c(2, 0, 5, 1)
  units_m <- c(10, 4, 12, 3)
  canonical <- make_family("binomial")
  custom <- make_family("binomial", make_inverse_link_functions(
    "binomial",
    list(inv = stats::plogis, d1 = stats::dlogis,
         d2 = function(x) stats::dlogis(x) * (1 - 2 * stats::plogis(x))),
    5L
  ))
  a <- canonical$terms(eta[1, ], y, units_m)
  b <- custom$terms(eta[1, ], y, units_m)
  expect_equal(b$log_likelihood, a$log_likelihood, tolerance = 1e-12)
  expect_equal(b$score, a$score, tolerance = 1e-12)
  expect_equal(b$curvature, a$curvature, tolerance = 1e-12)
  expect_equal(custom$log_likelihood(eta, y, units_m),
               canonical$log_likelihood(eta, y, units_m), tolerance = 1e-12)
  expect_equal(custom$slope_parts(eta, y, units_m),
               canonical$slope_parts(eta, y, units_m), tolerance = 1e-12)
})

test_that("the old entry points give what the object gives", {
  set.seed(8)
  eta <- matrix(rnorm(8), 2)
  y <- c(1, 3, 0, 2)
  units_m <- rep(5, 4)
  for (family in c("binomial", "poisson")) {
    object <- make_family(family)
    expect_identical(
      laplace_observation_terms(eta[1, ], y, units_m, family, object$linkf),
      object$terms(eta[1, ], y, units_m)
    )
    expect_identical(
      compute_observation_log_likelihood(eta, y, units_m, family,
                                         object$linkf),
      object$log_likelihood(eta, y, units_m)
    )
    expect_identical(
      compute_variance_slope_terms(eta, y, units_m, family, c(1, 1, 2, 2)),
      compute_slope_terms(eta, y, units_m, object, c(1, 1, 2, 2))
    )
  }
})

test_that("a fit's family object is rebuilt from what the fit stores", {
  object <- fit_family(binomial_model)
  expect_s3_class(object, "riskmap_family")
  expect_identical(object$family, "binomial")
  expect_true(object$canonical)
  expect_identical(fit_family(poisson_model)$family, "poisson")
  probit_fit <- list(family = "binomial", link_function = probit_family_link())
  expect_false(fit_family(probit_fit)$canonical)
  expect_identical(fit_family(probit_fit)$linkf$link, "probit")
})
