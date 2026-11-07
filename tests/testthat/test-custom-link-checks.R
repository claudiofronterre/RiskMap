## Custom links get the checks of the canonical links (#254): the test of a
## held variance, the screen of short scales and the exact alarms, and the
## derivatives the fit used in the later samplers.

logit_as_custom <- function() {
  make_inverse_link_functions("binomial", list(
    inv = function(eta) stats::plogis(eta),
    d1 = function(eta) stats::dlogis(eta),
    d2 = function(eta) stats::dlogis(eta) * (1 - 2 * stats::plogis(eta))
  ), 5L)
}

probit_link <- function() {
  make_inverse_link_functions("binomial", list(
    inv = function(eta) stats::pnorm(eta),
    d1 = function(eta) stats::dnorm(eta),
    d2 = function(eta) -eta * stats::dnorm(eta)
  ), 5L)
}

test_that("the slope terms of a custom link match those of the canonical link", {
  set.seed(3)
  eta <- matrix(rnorm(12, -1, 1.5), 3)
  y <- c(2, 0, 5, 1)
  units_m <- c(10, 4, 12, 3)
  groups <- list(1:4, c(1L, 1L, 2L, 2L))
  for (group in groups) {
    expect_equal(
      compute_variance_slope_terms(eta, y, units_m, "binomial", group,
                                   logit_as_custom()),
      compute_variance_slope_terms(eta, y, units_m, "binomial", group),
      tolerance = 1e-12
    )
  }
})

test_that("the slope terms of a probit link use its own derivatives", {
  set.seed(4)
  eta <- matrix(rnorm(8, 0, 1.2), 2)
  y <- c(3, 0, 6, 2)
  units_m <- c(8, 5, 9, 2)
  k <- function(value, i) {
    y[i] * stats::pnorm(value, log.p = TRUE) +
      (units_m[i] - y[i]) * stats::pnorm(value, lower.tail = FALSE, log.p = TRUE)
  }
  expected <- apply(eta, 1L, function(row) {
    first <- vapply(seq_along(row), function(i) {
      numDeriv::grad(function(x) k(x, i), row[i])
    }, numeric(1))
    second <- vapply(seq_along(row), function(i) {
      numDeriv::hessian(function(x) k(x, i), row[i])
    }, numeric(1))
    0.5 * (sum(first^2) + sum(second))
  })
  expect_equal(
    compute_variance_slope_terms(eta, y, units_m, "binomial", 1:4,
                                 probit_link()),
    expected, tolerance = 1e-6
  )
})

test_that("a held nugget is tested under a custom link", {
  coordinates <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, 4L, 1L), y = c(9, 0, 0, 9), units_m = rep(10, 4L),
    distances = as.matrix(dist(coordinates)), family = "binomial",
    fix_tau2 = TRUE, invlink = function(eta) stats::pnorm(eta)
  )
  index <- engine$layout$index
  expect_identical(engine$testable_variances(index$nu2), index$nu2)
  parameters <- numeric(4L)
  parameters[index$sigma2] <- log(0.5)
  parameters[index$phi] <- log(0.5)
  parameters[index$nu2] <- -18
  held <- assess_held_variances(engine, parameters, index$nu2, 2000L, 11L)
  expect_identical(held$parameter, index$nu2)
  expect_true(is.finite(held$slope))
  expect_gt(held$se, 0)
})

site_data <- function() {
  set.seed(8)
  list(D = cbind(1, seq(-1, 1, length.out = 12)), y = c(0, 1, 3, 0, 2, 5, 1, 0, 4, 6, 2, 1),
       units_m = rep(8, 12), location = rep(1:6, each = 2), offset = numeric(12))
}

test_that("the independent-site likelihood of a custom logit matches the canonical one", {
  s <- site_data()
  canonical <- make_independent_site_likelihood(s$D, s$y, s$units_m, s$location,
                                                s$offset, "binomial")
  custom <- make_independent_site_likelihood(s$D, s$y, s$units_m, s$location,
                                             s$offset, "binomial",
                                             linkf = logit_as_custom())
  a <- canonical$log_likelihood(c(-0.8, 0.4), log(0.7))
  b <- custom$log_likelihood(c(-0.8, 0.4), log(0.7))
  expect_equal(b$value, a$value, tolerance = 1e-10)
  expect_equal(b$gradient, a$gradient, tolerance = 1e-8)
  expect_equal(b$mean, a$mean, tolerance = 1e-8)
  expect_equal(b$variance, a$variance, tolerance = 1e-8)
})

test_that("the independent-site likelihood of a probit link matches integration", {
  s <- site_data()
  sites <- make_independent_site_likelihood(s$D, s$y, s$units_m, s$location,
                                            s$offset, "binomial",
                                            linkf = probit_link())
  beta <- c(-0.6, 0.5)
  variance <- 0.8
  exact <- sum(vapply(unique(s$location), function(site) {
    rows <- which(s$location == site)
    integrand <- function(effects) vapply(effects, function(effect) {
      eta <- drop(s$D[rows, , drop = FALSE] %*% beta) + effect
      exp(sum(s$y[rows] * stats::pnorm(eta, log.p = TRUE) +
                (s$units_m[rows] - s$y[rows]) *
                stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE))) *
        stats::dnorm(effect, 0, sqrt(variance))
    }, numeric(1))
    log(stats::integrate(integrand, -Inf, Inf, rel.tol = 1e-12)$value)
  }, numeric(1)))
  expect_equal(sites$log_likelihood(beta, log(variance))$value, exact,
               tolerance = 1e-8)
})

test_that("a fit with a custom link runs the screen and the exact alarms", {
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, invlink = function(eta) stats::pnorm(eta),
    control_mcmc = control_mcmc, seed = fixture_seed,
    control_mcml = fixture_control_mcml, messages = FALSE
  ))
  screen <- attr(fit, "optimizer")$mcml_maxima_screen
  expect_false(identical(screen$gate, "unavailable"))
  expect_null(screen$reason)
  expect_false(is.null(screen$independent))
})

test_that("errors while a proposal is built are counted", {
  ## A link that fails only beyond the range RiskMap checks, so it is accepted
  ## and then fails at a large linear predictor.
  failing <- list(
    inv = function(eta) {
      if (any(abs(eta) > 60)) stop("beyond the range of this link")
      stats::plogis(eta)
    },
    d1 = function(eta) stats::dlogis(eta),
    d2 = function(eta) stats::dlogis(eta) * (1 - 2 * stats::plogis(eta))
  )
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, 4L, 1L), y = c(1, 0, 0, 1), units_m = rep(10, 4L),
    distances = as.matrix(dist(1:4)), family = "binomial", invlink = failing
  )
  expect_identical(engine$proposal_errors()$count, 0L)
  draws <- matrix(rnorm(8), 2L)
  expect_identical(engine$evaluate(c(200, 0, 0), draws), -1e100)
  errors <- engine$proposal_errors()
  expect_identical(errors$count, 1L)
  expect_match(errors$first, "beyond the range of this link")

  row <- fit_check_importance_distribution(list(
    mcml_proposal = "global",
    mcml_proposal_errors = list(count = 3L, first = "boom")
  ))
  expect_identical(row$status, "OK")
  expect_match(row$finding, "could not be built at 3 parameter values")
  expect_match(row$finding, "boom")
  quiet <- fit_check_importance_distribution(list(
    mcml_proposal = "global", mcml_proposal_errors = list(count = 0L)
  ))
  expect_no_match(quiet$finding, "could not be built")
})

test_that("the later samplers get the derivatives the fit used", {
  link_function <- probit_link()
  captured <- NULL
  local_mocked_bindings(
    laplace_sampling_mcmc = function(..., invlink) {
      captured <<- invlink
      stop("captured")
    }
  )
  expect_error(sample_at_nongaussian_fit(
    estimate = binomial_model$estimate, y = binomial_model$y,
    units_m = binomial_model$units_m, D = binomial_model$D,
    coords = binomial_model$coords, cov_offset = numeric(length(binomial_model$y)),
    ID_coords = binomial_model$ID_coords, ID_re = NULL,
    fix_tau2 = binomial_model$fix_tau2, family = "binomial",
    link_function = link_function, kappa = 0.5,
    control_mcmc = control_mcmc, messages = FALSE
  ), "captured")
  expect_named(captured, c("inv", "d1", "d2"))
  expect_identical(captured$d2, link_function$d2)
})
