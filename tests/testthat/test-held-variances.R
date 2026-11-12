# An exact small model on a tensor grid of Gauss-Hermite nodes: `covariance` is
# the prior covariance of the latent effects at the observations.
make_exact_model <- function(y, units_m, offset, family, nodes = 40L) {
  rule <- gauss_hermite_rule(nodes)
  grid <- as.matrix(expand.grid(rep(list(sqrt(2) * rule$x), length(y))))
  weight <- Reduce(`*`, expand.grid(rep(list(rule$w / sqrt(pi)), length(y))))
  log_likelihood <- function(eta) {
    rowSums(sweep(eta, 2L, y, "*") -
              if (family == "poisson") {
                sweep(exp(eta), 2L, units_m, "*")
              } else {
                sweep(log1p(exp(eta)), 2L, units_m, "*")
              })
  }
  list(
    value = function(covariance) {
      eta <- sweep(grid %*% chol(covariance), 2L, offset, "+")
      log(sum(weight * exp(log_likelihood(eta))))
    },
    mean_of = function(covariance, statistic) {
      eta <- sweep(grid %*% chol(covariance), 2L, offset, "+")
      density <- weight * exp(log_likelihood(eta))
      sum(density * statistic(eta)) / sum(density)
    }
  )
}

expect_slope_matches_exact_model <- function(family, group, effect, y, units_m) {
  offset <- c(-1, -2, -0.5)
  covariance <- matrix(c(1, 0.5, 0.2, 0.5, 1, 0.4, 0.2, 0.4, 1), 3L)
  model <- make_exact_model(y, units_m, offset, family, nodes = 60L)
  step <- 1e-4
  values <- vapply(0:2 * step, function(variance) {
    model$value(covariance + variance * tcrossprod(effect))
  }, numeric(1))
  # A one-sided difference, because the variance cannot be negative.
  numerical <- (-3 * values[1L] + 4 * values[2L] - values[3L]) / (2 * step)
  terms <- model$mean_of(covariance, function(eta) {
    compute_slope_terms(eta, y, units_m, family_object(family), group)
  })
  expect_equal(terms, numerical, tolerance = 1e-5)
}

test_that("the slope terms give the derivative in a variance at zero", {
  y <- c(3, 0, 5)
  units_m <- c(10, 10, 12)
  # One effect for each observation, a nugget, and one shared by two of them.
  for (family in c("binomial", "poisson")) {
    expect_slope_matches_exact_model(family, 1:3, diag(3), y, units_m)
    expect_slope_matches_exact_model(family, c(1L, 1L, 2L),
                                     cbind(c(1, 1, 0), c(0, 0, 1)), y, units_m)
  }
})

# A small engine with a nugget, on four locations.
make_nugget_engine <- function(y, units_m, ID_re = NULL) {
  coordinates <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
  make_adaptive_mcml_engine(
    D = matrix(1, length(y), 1L), y = y, units_m = units_m,
    distances = as.matrix(dist(coordinates)), family = "binomial",
    ID_coords = seq_along(y), ID_re = ID_re, fix_tau2 = TRUE
  )
}

test_that("a held nugget is tested against the exact slope", {
  y <- c(9, 0, 0, 9)
  units_m <- rep(10, 4L)
  engine <- make_nugget_engine(y, units_m)
  index <- engine$layout$index
  parameters <- numeric(length(index$beta) + 3L)
  parameters[index$sigma2] <- log(0.5)
  parameters[index$phi] <- log(0.5)
  parameters[index$nu2] <- -18
  expect_identical(engine$testable_variances(index$nu2), index$nu2)
  expect_identical(engine$testable_variances(c(index$sigma2, index$nu2)),
                   integer())

  held <- assess_held_variances(engine, parameters, index$nu2, 40000L, 11L)
  expect_named(held, c("parameter", "slope", "se", "rising", "falling"))
  expect_identical(held$parameter, index$nu2)
  covariance <- 0.5 * exp(-as.matrix(dist(cbind(c(0, 1, 0, 1), c(0, 0, 1, 1)))) / 0.5)
  model <- make_exact_model(y, units_m, rep(0, 4L), "binomial", nodes = 24L)
  exact <- model$mean_of(covariance, function(eta) {
    compute_slope_terms(eta, y, units_m, make_family("binomial"), 1:4)
  })
  expect_lt(abs(held$slope - exact), 4 * held$se)
  expect_true(held$rising)
  expect_false(held$falling)
})

test_that("the test of a held variance is skipped only when it cannot be made", {
  # An engine from another package, without the test.
  expect_null(assess_held_variances(list(), c(0, 0), 2L, 10L, 1L))
  engine <- make_nugget_engine(c(1, 0, 0, 1), rep(10, 4L))
  expect_null(assess_held_variances(engine, numeric(4L), integer(), 10L, 1L))
  # A custom link is tested through its own derivatives (#254).
  custom <- make_adaptive_mcml_engine(
    D = matrix(1, 4L, 1L), y = c(1, 0, 0, 1), units_m = rep(10, 4L),
    distances = as.matrix(dist(1:4)), family = "binomial", fix_tau2 = TRUE,
    invlink = function(eta) stats::pnorm(eta)
  )
  expect_identical(custom$testable_variances(custom$layout$index$nu2),
                   custom$layout$index$nu2)
})

test_that("a random-effect variance is tested by its groups", {
  y <- c(9, 8, 0, 1)
  ID_re <- matrix(c(1L, 1L, 2L, 2L), ncol = 1L, dimnames = list(NULL, "cluster"))
  engine <- make_nugget_engine(y, rep(10, 4L), ID_re)
  index <- engine$layout$index
  expect_identical(engine$testable_variances(index$sigma2_re), index$sigma2_re)
  parameters <- numeric(max(unlist(index)))
  parameters[index$sigma2] <- log(0.5)
  parameters[index$phi] <- log(0.5)
  parameters[index$nu2] <- -2
  parameters[index$sigma2_re] <- -18
  held <- assess_held_variances(engine, parameters, index$sigma2_re, 5000L, 11L)
  expect_identical(held$parameter, index$sigma2_re)
  expect_true(is.finite(held$slope))
  expect_gt(held$se, 0)
})

test_that("a stage with a rising held variance is not certified", {
  engine <- list(latent_dimension = 1L, proposal_type = function() "laplace",
                 proposal_failures = function() 0L,
                 evaluate = function(parameters, draws, log_weights_only = FALSE,
                                     ...) {
                   if (log_weights_only) rep(0, nrow(draws)) else 0
                 })
  limits <- list(variance_index = 2L, linked_index = list(),
                 boundary_log_variance = -9, near_zero_log_variance = -4,
                 tolerance = 0.1)
  run_with <- function(rising) {
    local_mocked_bindings(
      optimise_on_bank = function(...) {
        list(fit = list(par = c(0, -12), objective = 1, convergence = 0L),
             boundary = 2L, restarted = FALSE, restart_improved = FALSE)
      },
      compute_observed_information = function(...) {
        list(information = diag(2), root = matrix(1), free = 1L, boundary = 2L,
             gradient = c(0, 0))
      },
      assess_held_variances = function(engine, parameters, boundary, n_draws,
                                       seed) {
        expect_identical(seed, 6L)
        data.frame(parameter = 2L, slope = 3, se = 0.5, rising = rising,
                   falling = FALSE)
      },
      correct_and_certify = function(...) {
        list(estimate = c(0, -12), outcome = "certified", certified = TRUE,
             certificate = 0.05)
      }
    )
    run_certified_stage(engine, function(...) 0, NULL, c(0, 0), 10L, 10L, 1:6,
                        c(-10, -20), c(10, 20), FALSE, limits)
  }
  rising <- run_with(TRUE)
  expect_identical(rising$outcome, "not certified")
  expect_true(rising$held_variance_rising)
  expect_identical(rising$estimate, c(0, -12))
  expect_identical(rising$next_start, c(0, log(0.25)))
  expect_match(rising$reason, "rises as a variance held at zero")

  undecided <- run_with(FALSE)
  expect_null(undecided$next_start)
  expect_null(undecided$held_variance_rising)
})

test_that("the next stage starts a rising held variance from the interior", {
  stages <- list(
    list(outcome = "not certified", estimate = c(0, -12), certificate = NA_real_,
         held_variance_rising = TRUE, next_start = c(0, log(0.25)),
         held_variances = data.frame(parameter = 2L, slope = 3, se = 0.5,
                                     rising = TRUE, falling = FALSE)),
    list(outcome = "certified", certified = TRUE, estimate = c(0, -1.5),
         certificate = 0.05)
  )
  starts <- list()
  budgets <- integer()
  seeds <- list()
  local_mocked_bindings(
    run_certified_stage = function(engine, objective, penalty, start, budget,
                                   validation_budget, seeds_of_stage, ...) {
      starts[[length(starts) + 1L]] <<- start
      budgets[length(budgets) + 1L] <<- budget
      seeds[[length(seeds) + 1L]] <<- seeds_of_stage
      stages[[length(starts)]]
    }
  )
  engine <- list(proposal_failures = function() 0L,
                 reset_proposal_failures = function() NULL,
                 proposal_type = function() "laplace",
                 latent_dimension = 1L)
  fit <- fit_adaptive_mcml(engine, start = c(0, 0), initial_budget = 1000L,
                           validation_budget = 5000L, seed = 1L)
  expect_identical(starts[[2L]], c(0, log(0.25)))
  # Five times the draws, as for a certificate of twice the tolerance.
  expect_identical(budgets, c(1000L, 5000L))
  expect_length(seeds[[1L]], 7L)
  expect_identical(seeds[[1L]][[6L]], 5002L)
  expect_identical(seeds[[1L]][[7L]], 8002L)
  expect_identical(fit$outcome, "certified")
  expect_false(fit$held_variance_rising)

  # The last stage keeps the estimate it fitted, not the start of a next one.
  local_mocked_bindings(
    run_certified_stage = function(...) stages[[1L]]
  )
  last <- fit_adaptive_mcml(engine, start = c(0, 0), maximum_stages = 1L)
  expect_identical(last$estimate, c(0, -12))
  expect_true(last$held_variance_rising)
  expect_identical(last$outcome, "not certified")
  expect_warning(
    report_mcml_outcome(last, list(index = list(sigma2 = 3L)), 0.1, FALSE),
    "rises as a variance held at zero"
  )

  # At the largest budget a rising held variance still gets a stage from the
  # interior, but more draws alone would not help any other uncertified stage.
  count <- 0L
  local_mocked_bindings(
    run_certified_stage = function(engine, objective, penalty, start, budget,
                                   ...) {
      count <<- count + 1L
      starts[[count]] <<- start
      budgets[count] <<- budget
      stages[[count]]
    }
  )
  capped <- fit_adaptive_mcml(engine, start = c(0, 0), initial_budget = 16000L,
                              maximum_budget = 16000L)
  expect_identical(count, 2L)
  expect_identical(budgets[1:2], c(16000L, 16000L))
  expect_identical(starts[[2L]], c(0, log(0.25)))
  expect_identical(capped$outcome, "certified")

  count <- 0L
  stages[[1L]]$held_variance_rising <- NULL
  stages[[1L]]$next_start <- NULL
  plain <- fit_adaptive_mcml(engine, start = c(0, 0), initial_budget = 16000L,
                             maximum_budget = 16000L)
  expect_identical(count, 1L)
  expect_identical(plain$outcome, "not certified")
})

# Sites on a grid whose effects are spatial, or the same effects shuffled.
make_grid_survey <- function(shuffle) {
  coordinates <- as.matrix(expand.grid(1:12, 1:12))
  distances <- as.matrix(dist(coordinates))
  set.seed(5L)
  effects <- drop(rnorm(144L) %*% chol(2 * exp(-distances / 4)))
  if (shuffle) effects <- sample(effects)
  units_m <- rep(30, 144L)
  y <- rbinom(144L, 30L, stats::plogis(-1 + effects))
  list(D = matrix(1, 144L, 1L), y = y, units_m = units_m, distances = distances)
}

judge_survey <- function(survey, tolerance = 0.1) {
  engine <- make_adaptive_mcml_engine(survey$D, survey$y, survey$units_m,
                                      survey$distances, family = "binomial",
                                      fix_tau2 = TRUE)
  sites <- make_independent_site_likelihood(
    survey$D, survey$y, survey$units_m, seq_along(survey$y),
    numeric(length(survey$y)), make_family("binomial")
  )
  screen <- make_short_scale_screen(engine, sites, survey$distances)
  fit <- screen$independent(0)
  index <- engine$layout$index
  estimate <- numeric(max(unlist(index)))
  estimate[index$beta] <- fit$beta
  estimate[index$sigma2] <- -12
  estimate[index$phi] <- 0
  estimate[index$nu2] <- fit$log_variance + 12
  list(alarm = judge_spatial_share(screen, fit, estimate, tolerance), fit = fit,
       screen = screen, engine = engine, estimate = estimate)
}

test_that("a spatial variance held at zero is judged over scales and shares", {
  judged <- judge_survey(make_grid_survey(FALSE))
  alarm <- judged$alarm
  expect_true(alarm$fired)
  expect_contains(judged$screen$spatial_shares, alarm$share)
  expect_contains(judged$screen$wide_log_phi, alarm$log_phi)
  # The bound is exact, against the maximum of the independent model.
  expect_gt(alarm$statistic, judged$fit$log_likelihood + 0.1^2 / 2)
  expect_identical(
    alarm$start[judged$engine$layout$index$phi], alarm$log_phi
  )
  # The start keeps the total variance of the independent model.
  values <- judged$engine$unpack(alarm$start)
  expect_equal(values$sigma2 * (1 + values$nugget_ratio),
               exp(judged$fit$log_variance), tolerance = 1e-12)

  # The bounds over the grid are kept, whether or not the alarm fires.
  grid <- expand.grid(log_phi = judged$screen$wide_log_phi,
                      share = judged$screen$spatial_shares)
  grid <- data.frame(log_phi = grid$log_phi, share = grid$share)
  expect_identical(alarm$shares[c("log_phi", "share")], grid)
  expect_identical(max(alarm$shares$bound), alarm$statistic)

  # The same effects without spatial structure raise no alarm.
  quiet <- judge_survey(make_grid_survey(TRUE))$alarm
  expect_false(quiet$fired)
  expect_named(quiet, c("fired", "shares"))
  expect_identical(quiet$shares[c("log_phi", "share")], grid)
  # Nor a gain below the tolerance's own scale.
  expect_false(judge_survey(make_grid_survey(FALSE), tolerance = 100)$alarm$fired)
})

test_that("a model without a nugget has no spatial share to judge", {
  screen <- list(spatial_shares = NULL)
  expect_null(judge_spatial_share(screen, list(), numeric(3L), 0.1))
})

test_that("the variances check reads the test of a held variance", {
  name <- "nugget variance (tau2)"
  check_with <- function(rising, falling) {
    fit_check_variances(list(
      mcml_boundary = name,
      mcml_held_variances = data.frame(parameter = 5L, slope = 1, se = 1,
                                       rising = rising, falling = falling,
                                       name = name)
    ))
  }
  rising <- check_with(TRUE, FALSE)
  expect_identical(rising$status, "PROBLEM")
  expect_match(rising$finding, "rises as the nugget variance")
  expect_match(rising$fix, "max_samples")
  falling <- check_with(FALSE, TRUE)
  expect_identical(falling$status, "OK")
  expect_match(falling$finding, "falls as the variance leaves zero")
  undecided <- check_with(FALSE, FALSE)
  expect_identical(undecided$status, "SUSPICIOUS")
  expect_match(undecided$finding, "effectively zero")
  expect_match(undecided$finding, "could not tell whether the likelihood rises")
  # Without a test, as with a custom link, the check does not claim one.
  untested <- fit_check_variances(list(mcml_boundary = name))
  expect_no_match(untested$finding, "could not tell")
})

test_that("the alarm's finding does not name a short scale", {
  alarm <- fit_check_maxima(list(mcml_maxima_alarm = TRUE))
  expect_identical(alarm$status, "PROBLEM")
  expect_match(alarm$finding, "at another spatial scale")
})

test_that("a rising held variance is reported in the certificate's sentence", {
  mcml <- list(outcome = "not certified", certificate = NA_real_, level = 0.95,
               tolerance = 0.1, correction_refused = FALSE,
               variance_rising = TRUE)
  expect_match(describe_mcml_certificate(mcml), "is not a maximum")
  expect_match(fit_check_monte_carlo(mcml)$finding, "not a maximum")
})
