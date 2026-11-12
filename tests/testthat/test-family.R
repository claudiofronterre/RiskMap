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
})

test_that("an object with integrated records keeps the properties of its link", {
  y <- c(0, 20, 3)
  units_m <- c(20, 20, 20)
  kernel <- make_integrated_kernel(make_site_table(y, units_m, "binomial"), 0.5)
  for (linkf in list(NULL, probit_family_link())) {
    plain <- make_family("binomial", linkf)
    integrated <- make_family("binomial", linkf, kernel = kernel)
    expect_true(integrated$integrated)
    expect_identical(integrated$canonical, plain$canonical)
    expect_identical(integrated$compiled, plain$compiled)
    expect_identical(integrated$needs_finite_weight, plain$needs_finite_weight)
    expect_false(integrated$assumes_concave)
    expect_false(integrated$moves_start)
    expect_identical(integrated$terms(c(-1, 2, 0), y, units_m),
                     kernel$terms(c(-1, 2, 0)))
    eta <- matrix(c(-1, 2, 0, 0.5, 1, -2), 2, byrow = TRUE)
    expect_identical(integrated$slope_parts(eta, y, units_m),
                     plain$slope_parts(eta, y, units_m))
  }
  expect_error(make_family("binomial", kernel), "goes in",
               class = "riskmap_error_internal")
})

test_that("a family object is defined only for the binomial and Poisson families", {
  expect_error(make_family("gaussian"), "binomial and Poisson",
               class = "riskmap_error_internal")
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

test_that("the canonical methods give the closed forms of their links", {
  set.seed(8)
  eta <- matrix(rnorm(8), 2)
  y <- c(1, 3, 0, 2)
  units_m <- c(5, 4, 6, 5)
  binomial <- make_family("binomial")
  probability <- stats::plogis(eta[1, ])
  terms <- binomial$terms(eta[1, ], y, units_m)
  expect_equal(terms$log_likelihood,
               sum(y * log(probability) + (units_m - y) * log1p(-probability)),
               tolerance = 1e-12)
  expect_equal(terms$score, y - units_m * probability, tolerance = 1e-12)
  expect_equal(terms$curvature, units_m * probability * (1 - probability),
               tolerance = 1e-12)
  rows <- apply(eta, 1L, function(row) {
    sum(stats::dbinom(y, units_m, stats::plogis(row), log = TRUE) -
          lchoose(units_m, y))
  })
  expect_equal(binomial$log_likelihood(eta, y, units_m), rows,
               tolerance = 1e-12)
  slope <- binomial$slope_parts(eta, y, units_m)
  rate <- stats::plogis(eta)
  expect_equal(slope$score, matrix(y, 2, 4, byrow = TRUE) -
                 sweep(rate, 2L, units_m, "*"), tolerance = 1e-12)
  expect_equal(slope$curvature, sweep(rate * (1 - rate), 2L, units_m, "*"),
               tolerance = 1e-12)

  poisson <- make_family("poisson")
  mean <- units_m * exp(eta[1, ])
  terms <- poisson$terms(eta[1, ], y, units_m)
  expect_equal(terms$log_likelihood, sum(y * log(mean) - mean),
               tolerance = 1e-12)
  expect_equal(terms$score, y - mean, tolerance = 1e-12)
  expect_equal(terms$curvature, mean, tolerance = 1e-12)
  rows <- apply(eta, 1L, function(row) sum(y * row - units_m * exp(row)))
  expect_equal(poisson$log_likelihood(eta, y, units_m), rows,
               tolerance = 1e-12)
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

test_that("the nonlinear parts of the canonical links agree with the direct difference", {
  set.seed(9)
  centre <- rnorm(4, -0.5)
  step <- matrix(rnorm(12, 0, 0.7), 4)
  y <- c(2, 0, 5, 1)
  units_m <- c(10, 4, 12, 3)
  custom_logit <- make_inverse_link_functions(
    "binomial",
    list(inv = stats::plogis, d1 = stats::dlogis,
         d2 = function(x) stats::dlogis(x) * (1 - 2 * stats::plogis(x))),
    5L
  )
  canonical <- make_family("binomial")$nonlinear_parts(centre, step, y, units_m)
  direct <- make_family("binomial", custom_logit)$nonlinear_parts(centre, step,
                                                                  y, units_m)
  expect_equal(direct$nonlinear, canonical$nonlinear, tolerance = 1e-10)
  expect_equal(direct$slope, canonical$slope, tolerance = 1e-12)
  rate <- units_m * exp(centre)
  poisson <- make_family("poisson")$nonlinear_parts(centre, step, y, units_m)
  expect_equal(poisson$nonlinear,
               (y * (centre + step) - units_m * exp(centre + step)) -
                 (y * centre - rate) - step * (y - rate),
               tolerance = 1e-10)
})

test_that("an object with integrated records takes its nonlinear parts from the kernel", {
  y <- c(0, 20, 3)
  units_m <- c(20, 20, 20)
  kernel <- make_integrated_kernel(make_site_table(y, units_m, "binomial"), 0.5)
  integrated <- make_family("binomial", kernel = kernel)
  step <- matrix(c(-0.5, 0.3, 0.1, 0.4, -0.2, 0.6), 3)
  expect_identical(integrated$nonlinear_parts(c(-1, 2, 0), step, y, units_m),
                   kernel$nonlinear_parts(c(-1, 2, 0), step))
})

test_that("only a custom link has kernel arrays for the site table", {
  expect_null(make_family("binomial")$kernel_arrays())
  arrays <- make_family("binomial", probit_family_link())$kernel_arrays()
  expect_identical(arrays,
                   compute_custom_kernel_arrays(probit_family_link(),
                                                "binomial"))
})

test_that("a family string with a flag that disagrees with its link is refused", {
  expect_identical(as_family_object("binomial", NULL, TRUE)$canonical, TRUE)
  expect_error(as_family_object("binomial", NULL, FALSE), "disagrees",
               class = "riskmap_error_internal")
  expect_error(as_family_object("binomial", probit_family_link(), TRUE),
               "disagrees", class = "riskmap_error_internal")
  object <- make_family("poisson")
  expect_identical(as_family_object(object), object)
})

test_that("an object comes alone, and a kernel is refused where the link is read", {
  object <- make_family("binomial")
  expect_error(as_family_object(object, probit_family_link(), FALSE),
               "takes no", class = "riskmap_error_internal")
  expect_error(as_family_object(object, NULL, TRUE), "takes no",
               class = "riskmap_error_internal")
  y <- c(0, 20, 3)
  units_m <- c(20, 20, 20)
  kernel <- make_integrated_kernel(make_site_table(y, units_m, "binomial"), 0.5)
  expect_error(make_site_table(y, units_m, "binomial", linkf = kernel),
               "model's link", class = "riskmap_error_internal")
  expect_error(make_site_table(y, units_m, make_family("binomial",
                                                        kernel = kernel)),
               "model's link", class = "riskmap_error_internal")
  expect_error(make_independent_site_likelihood(cbind(c(1, 1, 1)), y, units_m,
                                                1:3, numeric(3), "binomial",
                                                linkf = kernel),
               "model's link", class = "riskmap_error_internal")
})
