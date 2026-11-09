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
  # The interval of the level comes from resampling the banks.
  expect_lte(rates$level_interval[1], rates$level)
  expect_gte(rates$level_interval[2], rates$level)
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

test_that("triples of banks are distinct and uniform", {
  set.seed(8)
  triples <- draw_distinct_triples(5L, 60000L)
  expect_true(all(triples[, 1] != triples[, 2] & triples[, 1] != triples[, 3] &
                    triples[, 2] != triples[, 3]))
  expect_true(all(triples >= 1L & triples <= 5L))
  # Each of the 60 ordered triples has probability 1 / 60.
  counts <- table(paste(triples[, 1], triples[, 2], triples[, 3]))
  expect_length(counts, 60L)
  expect_lt(max(abs(counts / 60000 - 1 / 60)), 0.003)
})

test_that("the row of the radius says whether the radius holds", {
  details <- list(level = 3.2, level_interval = c(2.5, 4.1), margin = 0,
                  widening = 1, widened = 0.05, certificate = 0.05,
                  tolerance = 0.1)
  holds <- describe_radius_check(details)
  expect_identical(holds$status, "OK")
  expect_match(holds$finding, "3.2% of samples (2.5% to 4.1%)", fixed = TRUE)
  expect_identical(holds$fix, "")
  details[c("level", "margin", "widening", "widened")] <-
    list(6.5, 0.15, 1.3, 0.08)
  wide <- describe_radius_check(details)
  expect_identical(wide$status, "SUSPICIOUS")
  expect_match(wide$finding, "margin of 15%", fixed = TRUE)
  expect_match(wide$finding, "factor of 1.3", fixed = TRUE)
  expect_match(wide$fix, "max_samples", fixed = TRUE)
  expect_true(is_binomial_logit(binomial_model))
  expect_false(is_binomial_logit(poisson_model))
})

test_that("the radius of a fit is checked reproducibly", {
  set.seed(7)
  caller_seed <- .Random.seed
  checked <- fit_check_radius(binomial_model, chains = 4, rounds = 40,
                              banks = 20)
  expect_identical(.Random.seed, caller_seed)
  expect_named(checked, c("row", "data"))
  expect_identical(names(checked$row), c("check", "status", "finding", "fix"))
  expect_identical(checked$row$check, "Radius of the bound")
  expect_s3_class(checked$data, "RiskMap_radius_check")
  # One row for each resampled sample.
  expect_identical(dim(checked$data$samples), c(20L, 2L))
  expect_named(checked$data$samples, c("error", "radius"))
  expect_true(all(checked$data$samples$radius > 0))
  details <- checked$data
  expect_gte(details$own_alone, 0)
  expect_lte(details$own_alone, 100)
  for (rate in c("level", "own", "exceedance")) {
    expect_gte(details[[rate]], 0)
    expect_lte(details[[rate]], 100)
  }
  expect_length(details$level_interval, 2L)
  expect_lte(details$level_interval[1], details$level_interval[2])
  expect_gt(details$factor, 0)
  expect_identical(details$margin, 0)
  expect_gte(details$widening, 1)
  expect_gte(details$widened, details$certificate)
  expect_true(is.finite(details$pareto_k))
  expect_identical(details$settings$draws, 500L)
  expect_identical(details$settings$seed,
                   as.integer(binomial_model$seed + 9500L))
  expect_identical(details$settings$seed,
                   make_extra_check_seed(binomial_model$seed, 9500L))
  expect_gte(details$seconds, 0)
  without_time <- function(result) {
    result$data$seconds <- NULL
    result
  }
  again <- fit_check_radius(binomial_model, chains = 4, rounds = 40,
                            banks = 20)
  expect_identical(without_time(again), without_time(checked))
  other <- fit_check_radius(binomial_model, chains = 4, rounds = 40,
                            banks = 20, seed = 1)
  expect_identical(other$data$settings$seed, 1L)
  poisson <- fit_check_radius(poisson_model, chains = 4, rounds = 40,
                              banks = 20)
  expect_identical(poisson$data$margin, 0.15)
  # The streams wrap as those of a fit, and a fit without a seed starts at 1.
  expect_identical(make_radius_check_streams(NULL, .Machine$integer.max),
                   9500:9502)
  expect_identical(make_radius_check_streams(NULL, NULL), 1:3)
  expect_identical(make_radius_check_streams(7L, 12345L), 7:9)
})

test_that("the radius is checked only where there is a certificate", {
  check <- function(object, ...) fit_check_radius(object, ...)
  expect_error(check(list()), "'RiskMap' object")
  expect_error(check(gaussian_model),
               "cannot be checked, because a Gaussian fit has no Monte Carlo error")
  no_stages <- binomial_model
  no_stages$mcml_history <- NULL
  attr(no_stages, "optimizer")$mcml_history <- NULL
  expect_identical(find_radius_check_obstacle(no_stages),
                   "the fit has no stages of certified Monte Carlo maximum likelihood")
  expect_null(find_radius_check_obstacle(binomial_model))
  expect_error(check(binomial_model, rounds = 5), "at least 10")
  expect_error(check(binomial_model, banks = 2), "at least 3")
  expect_error(check(binomial_model, chains = 0))
  uncertified <- binomial_model
  uncertified$mcml_history[[1L]]$certificate_radius <- NA_real_
  expect_error(check(uncertified), "no certificate")
})

test_that("the radius check rebuilds the form of the fit's last stage", {
  fit <- list()
  attr(fit, "optimizer") <- list(mcml_integrate_noise = TRUE)
  expect_identical(find_radius_check_forms(fit, list()), "integrated")
  expect_identical(find_radius_check_forms(fit, list(form = "sampled")),
                   "sampled")
  attr(fit, "optimizer") <- list(mcml_integrate_noise = FALSE)
  expect_identical(find_radius_check_forms(fit, list()), "sampled")
  expect_identical(find_radius_check_forms(fit, list(form = "integrated")),
                   "integrated")
  combined <- list(form = "integrated",
                   second_validation = list(combined = list(used = "both")))
  expect_identical(find_radius_check_forms(fit, combined),
                   c("sampled", "integrated"))
  expect_identical(make_radius_check_streams(NULL, .Machine$integer.max, 2L),
                   9500:9505)
})

# A fit whose importance distribution does not depend on an anchor, so that
# the banks drawn on its streams are the fit's.
laplace_nugget_fit <- function(integrate_noise) {
  suppressWarnings(glgpm(
    y ~ cov + gp(nugget = TRUE), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, messages = FALSE,
    control_mcml = set_control_mcml(
      initial_samples = 500, validation_samples = 500, max_samples = 500,
      max_stages = 1, proposal = "laplace", integrate_noise = integrate_noise
    )
  ))
}

test_that("the rebuilt problem reproduces an integrated fit's bank", {
  skip_on_cran()
  fit <- laplace_nugget_fit(TRUE)
  last <- fit$mcml_history[[1L]]
  skip_if(!is.finite(last$certificate_radius %||% NA_real_))
  expect_identical(find_radius_check_forms(fit, last), "integrated")
  problem <- rebuild_certification_problem(fit, last, "integrated")
  bank <- draw_radius_bank(problem, fixture_seed + 3L)
  score <- compute_weighted_score(bank$log_weights, bank$gradients)
  expect_equal(compute_validation_radii(list(score$covariance)),
               last$certificate_radii[["certification"]], tolerance = 1e-6)
  expect_equal(sqrt(sum(score$score^2)), last$certificate_distance,
               tolerance = 1e-6)
})

test_that("combined banks are rebuilt and combined as the fit combined them", {
  skip_on_cran()
  fit <- laplace_nugget_fit("auto")
  last <- fit$mcml_history[[1L]]
  skip_if(!is.finite(last$certificate_radius %||% NA_real_))
  record <- last$second_validation$combined
  expect_false(is.null(record))
  expect_identical(record$k_max, 0.8)
  forms <- find_radius_check_forms(fit, last)
  expect_identical(forms, c("sampled", "integrated"))
  # The integrated bank is drawn on the stream 250 after the sampled one.
  streams <- c(sampled = fixture_seed + 3L, integrated = fixture_seed + 253L)
  banks <- lapply(stats::setNames(forms, forms), function(form) {
    problem <- rebuild_certification_problem(fit, last, form)
    draw_radius_bank(problem, streams[[form]])
  })
  scores <- lapply(banks, function(bank) {
    compute_weighted_score(bank$log_weights, bank$gradients)
  })
  pareto_k <- vapply(banks, `[[`, numeric(1), "pareto_k")
  expect_equal(pareto_k, record$pareto_k, tolerance = 1e-6)
  combined <- combine_two_scores(lapply(scores, `[[`, "score"),
                                 lapply(scores, `[[`, "covariance"),
                                 pareto_k, last$form, record$k_max)
  expect_identical(combined$used, record$used)
  expect_equal(compute_validation_radii(list(combined$covariance)),
               last$certificate_radii[["certification"]], tolerance = 1e-6)
  expect_equal(sqrt(sum(combined$score^2)), last$certificate_distance,
               tolerance = 1e-6)
  # The check combines its resampled banks in the same way.
  checked <- fit_check_radius(fit, chains = 4, rounds = 40, banks = 20)
  details <- checked$data
  expect_identical(details$forms, c("sampled", "integrated"))
  expect_named(details$used, c("both", "sampled", "integrated"))
  expect_equal(sum(details$used), 1)
  expect_identical(dim(details$samples), c(20L, 2L))
  expect_true(all(details$samples$radius > 0))
  expect_match(checked$row$finding, "combines a bank with the nugget sampled")
  expect_identical(details$settings$seed,
                   make_extra_check_seed(fit$seed, 9500L))
  expect_match(describe_radius_forms(details), "^sampled and integrated out")
})
