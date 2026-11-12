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
  kernel <- make_integrated_kernel(
    make_site_table(y, units_m, make_family("binomial")), 0.5
  )
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

test_that("the Gaussian object is exact, and the engine's readers refuse it", {
  gaussian <- make_family("gaussian")
  expect_s3_class(gaussian, "riskmap_family")
  expect_true(gaussian$exact)
  expect_true(gaussian$measurement_error)
  expect_null(gaussian$linkf)
  expect_null(gaussian$canonical)
  expect_null(gaussian$terms)
  expect_null(gaussian$log_likelihood)
  expect_error(check_family_object(gaussian), "fitted exactly",
               class = "riskmap_error_internal")
  expect_error(make_independent_site_likelihood(cbind(1), 8, 20, 1L, 0,
                                                make_family("gaussian")),
               "fitted exactly", class = "riskmap_error_internal")
  expect_false(make_family("binomial")$exact)
  expect_false(make_family("poisson")$exact)
  expect_error(make_family("gaussian", probit_family_link()), "exact",
               class = "riskmap_error_internal")
  expect_error(make_family("gamma"), "must be one of",
               class = "riskmap_error_internal")
  expect_identical(fit_family(gaussian_model)$family, "gaussian")
  expect_true(fit_family(gaussian_model)$exact)
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
  kernel <- make_integrated_kernel(
    make_site_table(y, units_m, make_family("binomial")), 0.5
  )
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

test_that("the engine takes a family object, not a family string", {
  object <- make_family("poisson")
  expect_identical(check_family_object(object), object)
  latent_map <- make_latent_map(cbind(1L), 1L)
  expect_error(find_latent_mode(0, 1, 5, latent_map, diag(1),
                                make_family("gaussian"), 0),
               "fitted exactly", class = "riskmap_error_internal")
  expect_error(find_latent_mode(0, 1, 5, latent_map, diag(1), "binomial", 0),
               "family object", class = "riskmap_error_internal")
  expect_error(compute_slope_terms(matrix(0), 1, 5, make_family("gaussian"),
                                   1L),
               "fitted exactly", class = "riskmap_error_internal")
  expect_error(check_family_object("binomial"), "family object",
               class = "riskmap_error_internal")
  expect_error(make_site_table(c(0, 2), c(5, 5), "binomial"), "family object",
               class = "riskmap_error_internal")
  expect_error(fit_global_proposal(0, 0, 1, 5, make_latent_map(cbind(1L), 1L),
                                   diag(1), "binomial"),
               "family object", class = "riskmap_error_internal")
})

test_that("a kernel is refused where the link is read", {
  y <- c(0, 20, 3)
  units_m <- c(20, 20, 20)
  kernel <- make_integrated_kernel(
    make_site_table(y, units_m, make_family("binomial")), 0.5
  )
  expect_error(make_site_table(y, units_m, family_object("binomial", kernel)),
               "model's link", class = "riskmap_error_internal")
  expect_error(make_site_table(y, units_m,
                               make_family("binomial", kernel = kernel)),
               "model's link", class = "riskmap_error_internal")
  expect_error(make_independent_site_likelihood(
    cbind(c(1, 1, 1)), y, units_m, 1:3, numeric(3),
    family_object("binomial", kernel)
  ), "model's link", class = "riskmap_error_internal")
})

test_that("the traits of each family give their closed forms", {
  y <- c(0, 3, 5)
  units_m <- c(5, 5, 6)
  eta <- c(-1, 0.5, 2)
  binomial <- family_traits("binomial")
  poisson <- family_traits("poisson")
  expect_identical(binomial$empirical_link(y, units_m),
                   log((y + 0.5) / (units_m - y + 0.5)))
  expect_identical(poisson$empirical_link(y, units_m),
                   log((y + 0.5) / units_m))
  expect_identical(binomial$transform_variance(y, units_m),
                   median(1 / (y + 0.5) + 1 / (units_m - y + 0.5)))
  expect_identical(poisson$transform_variance(y, units_m),
                   median(1 / (y + 0.5)))
  expect_identical(binomial$start_log_likelihood(eta, y, units_m),
                   sum(y * eta - units_m * softplus(eta)))
  expect_identical(poisson$start_log_likelihood(eta, y, units_m),
                   sum(y * eta - units_m * exp(eta)))
  expect_null(poisson$start_log_likelihood(c(eta, 800), c(y, 1),
                                           c(units_m, 1)))
  expect_identical(poisson$fitted_mean(c(0.1, 0.2), c(2, 3)), c(2, 3) * c(0.1, 0.2))
  expect_identical(binomial$fitted_mean(c(0.1, 0.2), c(2, 3)), c(0.1, 0.2))
  expect_error(binomial$validate(c(1, 7), c(5, 5)), "exceed")
  expect_true(poisson$validate(c(1, 7), c(5, 5)))
  expect_true(binomial$upper_limit)
  expect_false(poisson$upper_limit)
  expect_true(family_traits("gaussian")$exact)
})

test_that("only the families the engine fits pass its check", {
  expect_identical(check_engine_family("poisson")$family, "poisson")
  expect_error(check_engine_family("gaussian"),
               class = "riskmap_error_internal")
  expect_error(check_engine_family(c("binomial", "poisson")),
               class = "riskmap_error_internal")
  expect_null(fit_traits(list(family = "gamma")))
  expect_true(fit_traits(gaussian_model)$exact)
  expect_false(fit_traits(binomial_model)$exact)
})

test_that("the start of each family is the regression its caller ran", {
  y <- poisson_data$y
  units_m <- poisson_data$denominator
  D <- cbind(1, poisson_data$cov)
  expected <- coef(glm(y ~ cov + offset(log(denominator)) + offset(offset),
                       data = poisson_data, family = poisson))
  start <- family_traits("poisson")$start_coefficients(
    y, units_m, D, poisson_data$offset, NULL, y ~ cov, poisson_data
  )
  expect_equal(start, expected, tolerance = 1e-10)

  y <- c(0, 3, 5, 2, 4)
  units_m <- c(5, 5, 6, 4, 7)
  D <- cbind(1, c(-1, 0, 1, 2, 0.5))
  offset <- c(0.1, 0, -0.2, 0.3, 0)
  binomial <- family_traits("binomial")
  for (link in list(NULL, "probit")) {
    glm_family <- if (is.null(link)) stats::binomial() else stats::binomial(link)
    expected <- stats::glm.fit(D, cbind(y, units_m - y), offset = offset,
                               family = glm_family)$coefficients
    expect_identical(
      binomial$start_coefficients(y, units_m, D, offset, link, NULL, NULL),
      expected
    )
  }
  response <- c(1.2, 0.4, -0.3, 2.2, 0.9)
  expect_identical(
    family_traits("gaussian")$start_coefficients(response, NULL, D, NULL,
                                                 NULL, NULL, NULL),
    linear_start_values(D, response)
  )
})

test_that("the family's simulator and scoring support take the draws of their callers", {
  means <- c(0.2, 0.7, 0.5)
  units_m <- c(4, 10, 6)
  set.seed(3)
  expected <- rbinom(3, size = units_m, prob = means)
  set.seed(3)
  expect_identical(family_traits("binomial")$simulate(means, units_m, NULL),
                   expected)
  set.seed(4)
  expected <- rpois(3, lambda = means * units_m)
  set.seed(4)
  expect_identical(family_traits("poisson")$simulate(means, units_m, NULL),
                   expected)
  set.seed(5)
  expected <- means + rnorm(3, sd = sqrt(0.3))
  set.seed(5)
  expect_identical(
    family_traits("gaussian")$simulate(means, NULL, list(sigma2_me = 0.3)),
    expected
  )
  draws <- c(0, 3, 9)
  expect_identical(family_traits("binomial")$score_support(draws, 2, 10, means),
                   0:10)
  # Each of the three terms of the Poisson support is the largest once: the
  # draws, the response, and the 0.999 quantile of the mean rate.
  poisson <- family_traits("poisson")
  expect_identical(poisson$score_support(draws, 4, 2, means), 0:9)
  expect_identical(poisson$score_support(draws, 12, 2, means), 0:12)
  rates <- c(6, 7, 8)
  expect_gt(qpois(0.999, mean(4 * rates)), 12)
  expect_identical(poisson$score_support(draws, 4, 4, rates),
                   0:qpois(0.999, mean(4 * rates)))
  # A matrix of means with a denominator for each row, as draw_simulation()
  # passes them.
  matrix_means <- matrix(c(0.2, 0.7, 0.5, 0.1, 0.9, 0.4), 3)
  set.seed(6)
  expected <- rbinom(6, size = units_m, prob = as.numeric(matrix_means))
  set.seed(6)
  expect_identical(
    family_traits("binomial")$simulate(matrix_means, units_m, NULL), expected
  )
  set.seed(7)
  expected <- rpois(6, lambda = as.numeric(matrix_means) * units_m)
  set.seed(7)
  expect_identical(
    family_traits("poisson")$simulate(matrix_means, units_m, NULL), expected
  )
  expect_null(family_traits("gaussian")$score_support(draws, 1, 1, means))
  expect_identical(family_traits("binomial")$inverse_link, stats::plogis)
  expect_identical(family_traits("poisson")$inverse_link, exp)
  expect_identical(family_traits("gaussian")$inverse_link, identity)
  expect_false(family_traits("gaussian")$discrete)
  expect_true(family_traits("poisson")$discrete)
  expect_true(family_traits("binomial")$whole_denominator)
  expect_false(family_traits("poisson")$whole_denominator)
  # A fit or a template without a family is refused, not taken as inexact.
  no_family <- binomial_model
  no_family["family"] <- list(NULL)
  expect_error(assessment_refit_args(no_family, no_family$formula,
                                     no_family$data, control_mcmc),
               class = "riskmap_error_internal")
  expect_error(simulation_model(no_family), class = "riskmap_error_internal")
})

test_that("a fit's family says whether it is exact and how its mean is named", {
  expect_true(fit_is_exact(gaussian_model))
  expect_false(fit_is_exact(binomial_model))
  expect_false(fit_is_exact(list(family = "gamma")))
  expect_identical(family_traits("binomial")$response_quantity, "prevalence")
  expect_identical(family_traits("poisson")$response_quantity, "rate")
  expect_null(family_traits("gaussian")$response_quantity)
  expect_true(family_traits("binomial")$probability)
  expect_false(family_traits("poisson")$probability)
})
