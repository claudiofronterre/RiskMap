test_that("the radii of several banks are those of compute_validation_radius()", {
  set.seed(3)
  covariances <- lapply(1:3, function(bank) {
    factor <- matrix(stats::rnorm(9), 3)
    crossprod(factor) / 100
  })
  expect_identical(
    compute_validation_radii(covariances),
    vapply(covariances, function(covariance) {
      compute_validation_radius(diag(3), covariance)
    }, numeric(1))
  )
})

test_that("the bridge equation and the split R-hat behave as they should", {
  # With pi~ / q constant the normalising constant is that constant.
  expect_equal(solve_bridge_equation(rep(2, 10), 6L, 4L), 2, tolerance = 1e-8)
  set.seed(4)
  expect_equal(compute_split_rhat(matrix(stats::rnorm(4000), 1000)), 1,
               tolerance = 0.02)
  shifted <- cbind(stats::rnorm(500), stats::rnorm(500, mean = 3))
  expect_gt(compute_split_rhat(shifted), 1.5)
})

test_that("elliptical slice sampling targets the prior times the likelihood", {
  # A N(0, 1) prior and a likelihood exp(-3 z^2 / 2) give N(0, 1 / 4).
  core <- function(draws) -1.5 * rowSums(draws^2)
  set.seed(5)
  state <- matrix(stats::rnorm(4000), 4000)
  log_likelihood <- core(state)
  for (round in 1:30) {
    moved <- move_elliptical_slice(state, log_likelihood, core)
    state <- moved$state
    log_likelihood <- moved$log_likelihood
  }
  expect_equal(mean(state), 0, tolerance = 0.03)
  expect_equal(stats::var(as.numeric(state)), 0.25, tolerance = 0.03)
})

test_that("the rates of the radius are measured on triples of banks", {
  # Errors from N(0, I) with the radius at their exact 95% point.
  set.seed(6)
  errors <- matrix(stats::rnorm(8000), 4000)
  radius <- sqrt(stats::qchisq(0.95, 2))
  rates <- measure_radius_rule(errors, rep(radius, 4000))
  expect_lt(abs(rates$own - 0.05), 0.006)
  # The error the radius allows for is N(0, I / 2), which exceeds the radius
  # with probability exp(-radius^2).
  expect_lt(abs(rates$level - exp(-radius^2)), 0.0015)
  expect_lt(abs(rates$factor - 1 / sqrt(2)), 0.02)
  expect_lt(rates$exceedance, rates$level)
  # The radius of a triple is the larger of the certification bank's and its
  # replicate's, so it is half as large only when both radii are, in a
  # quarter of triples, where ||G|| exceeds it with probability
  # exp(-radius^2 / 8).
  halved <- measure_radius_rule(errors, rep(c(radius, radius / 2), 2000))
  expect_lt(abs(halved$own - (0.75 * 0.05 + 0.25 * exp(-radius^2 / 8))),
            0.01)
})

test_that("the rebuilt problem reproduces the fit's certification bank", {
  # With the Laplace proposal the importance distribution does not depend on
  # an anchor, so the bank drawn on the fit's certification stream is the
  # fit's.
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, messages = FALSE,
    control_mcml = set_control_mcml(
      initial_samples = 500, validation_samples = 500, max_samples = 500,
      max_stages = 1, proposal = "laplace"
    )
  ))
  last <- fit$mcml_history[[1L]]
  skip_if(!is.finite(last$certificate_radius %||% NA_real_))
  problem <- rebuild_certification_problem(fit, last)
  expect_identical(problem$draws, 500L)
  bank <- draw_radius_bank(problem, fixture_seed + 3L)
  score <- compute_weighted_score(bank$log_weights, bank$gradients)
  expect_equal(compute_validation_radii(list(score$covariance)),
               last$certificate_radii[["certification"]], tolerance = 1e-6)
  expect_equal(sqrt(sum(score$score^2)), last$certificate_distance,
               tolerance = 1e-6)
})

test_that("the radius of a fit is measured reproducibly", {
  set.seed(7)
  caller_seed <- .Random.seed
  measured <- measure_certificate_radius(binomial_model, chains = 4,
                                         rounds = 40, banks = 20)
  expect_identical(.Random.seed, caller_seed)
  for (rate in c("level", "own", "exceedance")) {
    expect_gte(measured[[rate]], 0)
    expect_lte(measured[[rate]], 100)
  }
  expect_gt(measured$factor, 0)
  expect_identical(measured$margin, 0)
  expect_gte(measured$inflated, measured$certificate)
  expect_identical(measured$settings$draws, 500L)
  expect_identical(measured$settings$seed,
                   as.integer(binomial_model$seed + 9500L))
  expect_gte(measured$seconds, 0)
  again <- measure_certificate_radius(binomial_model, chains = 4, rounds = 40,
                                      banks = 20)
  expect_identical(again[names(again) != "seconds"],
                   measured[names(measured) != "seconds"])
  other <- measure_certificate_radius(binomial_model, chains = 4, rounds = 40,
                                      banks = 20, seed = 1)
  expect_identical(other$settings$seed, 1L)
  poisson <- measure_certificate_radius(poisson_model, chains = 4,
                                        rounds = 40, banks = 20)
  expect_identical(poisson$margin, 0.15)
})

test_that("the radius is measured only where there is a certificate", {
  measure <- function(object, ...) measure_certificate_radius(object, ...)
  expect_error(measure(list()), "'RiskMap' object")
  expect_error(measure(gaussian_model), "binomial and Poisson")
  expect_error(measure(binomial_model, rounds = 5), "at least 10")
  expect_error(measure(binomial_model, banks = 2), "at least 3")
  expect_error(measure(binomial_model, chains = 0))
  uncertified <- binomial_model
  uncertified$mcml_history[[1L]]$certificate_radius <- NA_real_
  expect_error(measure(uncertified), "no certificate")
})
