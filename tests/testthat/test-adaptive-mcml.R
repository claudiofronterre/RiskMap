test_that("adaptive importance likelihood agrees with one-dimensional quadrature", {
  y <- 2
  trials <- 5
  beta <- -0.4
  sigma2 <- 0.8
  parameters <- c(beta, log(sigma2), log(1))
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 1), y = y, units_m = trials,
    distances = matrix(0), family = "binomial"
  )
  exact <- log(integrate(
    function(spatial_effect) {
      linear_predictor <- beta + spatial_effect
      exp(y * linear_predictor - trials * softplus(linear_predictor)) *
        dnorm(spatial_effect, sd = sqrt(sigma2))
    },
    lower = -Inf, upper = Inf, rel.tol = 1e-12
  )$value)

  set.seed(8128)
  bank <- matrix(rnorm(20000), ncol = 1)
  observed <- engine$evaluate(parameters, bank)

  expect_equal(observed, exact, tolerance = 0.005)
})

test_that("a certified stage corrects the fit and certifies the correction", {
  coordinate <- seq(0, 1, length.out = 30)
  covariate <- seq(-1, 1, length.out = 30)
  trials <- rep(20, 30)
  y <- round(trials * plogis(-0.4 + 0.5 * covariate +
                               0.4 * sin(2 * pi * coordinate)))
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y,
    units_m = trials, distances = as.matrix(dist(coordinate)),
    family = "binomial"
  )
  fit <- fit_adaptive_mcml(
    engine = engine, start = c(-0.4, 0.5, log(0.8), log(0.4)),
    initial_budget = 1000, validation_budget = 2000, maximum_budget = 1000,
    maximum_stages = 1, seed = 2718,
    lower = c(-10, -10, -8, -8), upper = c(10, 10, 4, 4)
  )
  stage <- fit$history[[1L]]

  expect_identical(stage$proposal, "global")
  expect_identical(stage$proposal_failures, 0L)
  expect_true(is.finite(stage$certificate))
  expect_equal(stage$certificate,
               stage$certificate_distance + stage$certificate_radius)
  expect_identical(fit$certified, stage$certificate <= 0.1)
  expect_identical(fit$outcome,
                   if (fit$certified) "certified" else "not certified")
  # The returned estimate is the corrected one, not the optimum of the bank.
  expect_identical(fit$estimate, stage$estimate)
  expect_false(identical(stage$estimate, stage$fitted_estimate))
  expect_equal(fit$score, stage$second_validation$score)
  expect_identical(fit$validation_budget, 2000L)
})

test_that("boundary variance components and their linked range are detected", {
  estimate <- c(0.1, 0.2, -12, 1.5, -0.5)
  expect_identical(
    adaptive_boundary_parameters(estimate, c(3L, 5L), -9),
    3L
  )
  expect_identical(
    adaptive_boundary_parameters(estimate, c(3L, 5L), -9,
                                 linked_index = list(`3` = 4L)),
    c(3L, 4L)
  )
  expect_identical(
    adaptive_boundary_parameters(c(0.1, 0.2, -1, 1.5), 3L, -9,
                                 linked_index = list(`3` = 4L)),
    integer()
  )
})

test_that("variance parameters are described in plain language", {
  layout <- cov_par_layout(2L, fix_tau2 = TRUE, n_re = 1L)
  expect_identical(
    describe_variance_parameters(c(5L, 6L), layout, re_names = "village"),
    c("nugget variance (tau2)", "random-effect variance for 'village'")
  )
  expect_identical(describe_variance_parameters(integer(), layout),
                   character())
})

test_that("the Laplace start returns the Laplace value at its estimate", {
  coordinate <- seq(0, 1, length.out = 25)
  covariate <- seq(-1, 1, length.out = 25)
  trials <- rep(15, 25)
  y <- round(trials * plogis(-0.3 + 0.4 * covariate +
                               0.5 * sin(2 * pi * coordinate)))
  distances <- as.matrix(dist(coordinate))
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y, units_m = trials,
    distances = distances, family = "binomial", proposal = "laplace"
  )
  start <- c(0, 0, log(0.5), log(0.3))
  laplace <- maximise_laplace(engine$evaluate, start, 25,
                              lower = c(-10, -10, -8, -8),
                              upper = c(10, 10, 4, 4))
  zero <- matrix(0, nrow = 1, ncol = 25)

  expect_equal(laplace$log_likelihood,
               engine$evaluate(laplace$estimate, zero))
  expect_gte(laplace$log_likelihood, engine$evaluate(start, zero))
})

test_that("a zero spatial variance is held at the boundary without escalation", {
  set.seed(31)
  coordinate <- seq(0, 1, length.out = 60)
  covariate <- rnorm(60)
  trials <- rep(25, 60)
  # No spatial effect: the maximum likelihood spatial variance is zero.
  y <- rbinom(60, trials, plogis(-0.5 + 0.6 * covariate))
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y, units_m = trials,
    distances = as.matrix(dist(coordinate)), family = "binomial"
  )
  fit <- fit_adaptive_mcml(
    engine = engine, start = c(-0.5, 0.6, -14, log(0.2)),
    initial_budget = 1000, validation_budget = 2000, maximum_budget = 8000,
    maximum_stages = 4, seed = 99,
    lower = c(-10, -10, -20, -8), upper = c(10, 10, 4, 4),
    variance_index = 3L, linked_index = list(`3` = 4L)
  )

  expect_identical(fit$boundary_parameters, c(3L, 4L))
  expect_identical(fit$outcome, "certified")
  expect_identical(fit$budget, 1000L)
  # Held parameters are not corrected.
  expect_identical(fit$estimate[3:4], fit$history[[1L]]$fitted_estimate[3:4])
  expect_true(all(is.finite(diag(fit$covariance)[1:2])))
  expect_true(all(is.na(diag(fit$covariance)[3:4])))
})

test_that("set_control_mcml validates laplace_start", {
  expect_true(set_control_mcml()$laplace_start)
  expect_false(set_control_mcml(laplace_start = FALSE)$laplace_start)
  expect_error(set_control_mcml(laplace_start = NA), "laplace_start")
  expect_error(set_control_mcml(laplace_start = "yes"), "laplace_start")
})

test_that("median nearest-neighbour distance ignores duplicate locations", {
  distances <- as.matrix(dist(c(0, 0, 1, 3)))
  expect_identical(median_nearest_neighbour_distance(distances), 1)
})

test_that("a scale smaller than the location spacing is reported, not held", {
  fit <- list(estimate = c(0, 0, log(0.5), log(0.3)),
              boundary_parameters = integer())
  layout <- cov_par_layout(2L, fix_tau2 = FALSE)
  distances <- as.matrix(dist(seq(0, 10, length.out = 11)))
  small <- detect_small_scale(fit, layout, distances, kappa = 0.5,
                               distance_scale = 2)
  expect_equal(small, c(scale = 0.6, practical_range = 0.6 * -log(0.05),
                        spacing = 2), tolerance = 1e-8)
  # With the spatial variance at zero the scale is meaningless and is not
  # reported as small.
  fit$boundary_parameters <- 3:4
  expect_null(detect_small_scale(fit, layout, distances, 0.5, 2))

  outcome <- list(boundary_parameters = integer(), budget = 1000L,
                  validation_budget = 5000L, outcome = "certified",
                  certified = TRUE, certificate = 0.05, history = list(1))
  expect_warning(
    report_mcml_outcome(outcome, character(), layout, pc_prior = FALSE,
                        tolerance = 0.1, messages = FALSE,
                        small_scale = small),
    "smaller than the median distance between neighbouring locations"
  )
})

test_that("the practical range factor matches the exponential correlation", {
  expect_equal(compute_practical_range_factor(0.5), -log(0.05),
               tolerance = 1e-8)
  expect_gt(compute_practical_range_factor(1.5),
            compute_practical_range_factor(0.5))
})

test_that("the PC penalty is a proper density on the working scale", {
  layout <- cov_par_layout(2L, fix_tau2 = FALSE)
  penalty <- make_pc_penalty(layout, scale_threshold = 2,
                             sigma_threshold = 1, sigma_probability = 0.05,
                             scale_probability = 0.05)
  at <- function(log_sigma2, log_phi) penalty(c(0, 0, log_sigma2, log_phi))
  density_sigma <- function(a) exp(vapply(a, at, 0, log_phi = 0))
  density_range <- function(b) exp(vapply(b, function(x) at(0, x), 0))
  sigma_total <- integrate(density_sigma, -40, 15)$value
  range_total <- integrate(density_range, -30, 30)$value
  # The penalty is separable, log pi = g(log sigma^2) + h(log phi), so it
  # integrates to one exactly when sigma_total * range_total = exp(g(0) + h(0)).
  expect_equal(sigma_total * range_total, exp(at(0, 0)), tolerance = 1e-6)
  # Tail probabilities match the specification.
  expect_equal(integrate(density_sigma, log(1), 15)$value / sigma_total, 0.05,
               tolerance = 1e-6)
  expect_equal(integrate(density_range, -30, log(2))$value / range_total, 0.05,
               tolerance = 1e-6)
})

test_that("penalising the engine leaves importance weights unchanged", {
  coordinate <- seq(0, 1, length.out = 20)
  covariate <- seq(-1, 1, length.out = 20)
  trials <- rep(10, 20)
  y <- round(trials * plogis(-0.2 + 0.3 * covariate))
  engine <- make_adaptive_mcml_engine(
    D = cbind(1, covariate), y = y, units_m = trials,
    distances = as.matrix(dist(coordinate)), family = "binomial"
  )
  penalty <- function(parameters) -sum(parameters^2)
  penalised <- penalise_engine(engine$evaluate, penalty)
  theta <- c(-0.2, 0.3, log(0.5), log(0.3))
  set.seed(3)
  z <- matrix(rnorm(200 * 20), 200)
  expect_identical(penalised(theta, z, log_weights_only = TRUE),
                   engine$evaluate(theta, z, log_weights_only = TRUE))
  expect_equal(penalised(theta, z), engine$evaluate(theta, z) + penalty(theta))
})

test_that("set_control_mcml validates the PC prior options", {
  expect_false(set_control_mcml()$pc_prior)
  control <- set_control_mcml(pc_prior = TRUE, pc_scale_threshold = 5)
  expect_true(control$pc_prior)
  expect_identical(control$pc_scale_threshold, 5)
  expect_error(set_control_mcml(pc_prior = NA), "pc_prior")
  expect_error(set_control_mcml(pc_sigma_probability = 1),
               "pc_sigma_probability")
})

test_that("glgpm fits with the PC prior and reports the unpenalised likelihood", {
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    control_mcml = set_control_mcml(
      initial_samples = 500, validation_samples = 500, max_samples = 500,
      max_stages = 1, pc_prior = TRUE
    ),
    messages = FALSE
  ))
  expect_true(attr(fit, "optimizer")$mcml_pc_prior)
  expect_true(all(is.finite(unlist(fit$estimate))))
  expect_true(is.finite(fit$log_lik))
})

test_that("identification warnings point to the PC penalty only when it is off", {
  fit_warnings <- function(pc_prior) {
    warnings <- character()
    withCallingHandlers(
      glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial",
            denominator = denominator, control_mcmc = control_mcmc,
            control_mcml = set_control_mcml(
              initial_samples = 500, validation_samples = 500,
              max_samples = 500, max_stages = 1, pc_prior = pc_prior
            ),
            messages = FALSE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    warnings
  }
  unpenalised <- fit_warnings(FALSE)
  identification <- grepl("effectively zero|smaller than the median|singular",
                          unpenalised)
  skip_if_not(any(identification), "fixture fit is well identified")
  expect_true(all(grepl("pc_prior = TRUE", unpenalised[identification])))
  expect_false(any(grepl("pc_prior = TRUE", fit_warnings(TRUE))))
})

test_that("defensive rows are spread evenly and absent from a single row", {
  expect_identical(sum(defensive_rows(1, 0.05)), 0L)
  expect_identical(sum(defensive_rows(19, 0.05)), 0L)
  expect_identical(which(defensive_rows(60, 0.05)), c(20L, 40L, 60L))
  expect_identical(sum(defensive_rows(1000, 0.05)), 50L)
  expect_false(any(defensive_rows(1000, 0)))
})

test_that("defensive mixture weights use the full mixture density", {
  coordinate <- seq(0, 1, length.out = 6)
  trials <- rep(10, 6)
  y <- c(1, 3, 2, 6, 4, 5)
  distances <- as.matrix(dist(coordinate))
  parameters <- c(-0.3, log(0.7), log(0.3))
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 6), y = y, units_m = trials,
    distances = distances, family = "binomial", defensive_weight = 0.25
  )
  set.seed(31)
  bank <- matrix(rnorm(8 * 6), nrow = 8)
  log_weights <- engine$evaluate(parameters, bank, log_weights_only = TRUE)

  proposal <- engine$proposal_at(parameters)
  gaussian_covariance <- crossprod(proposal$root)
  prior_covariance <- proposal$covariance
  defensive <- defensive_rows(8, 0.25)
  draws <- bank %*% proposal$root
  draws[defensive, ] <- bank[defensive, ] %*% chol(prior_covariance)
  draws <- sweep(draws, 2L, proposal$centre, "+")
  log_normal <- function(x, mean, covariance) {
    root <- chol(covariance)
    residual <- backsolve(root, x - mean, transpose = TRUE)
    -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(root))) +
              sum(residual^2))
  }
  expected <- vapply(seq_len(8), function(row) {
    u <- draws[row, ]
    # The engine omits the binomial coefficients, which do not depend on the
    # parameters.
    linear_predictor <- -0.3 + u
    log_target <- sum(y * linear_predictor -
                        trials * softplus(linear_predictor)) +
      log_normal(u, numeric(6), prior_covariance)
    log_proposal <- log(
      0.75 * exp(log_normal(u, proposal$centre, gaussian_covariance)) +
        0.25 * exp(log_normal(u, proposal$centre, prior_covariance))
    )
    log_target - log_proposal
  }, numeric(1))

  expect_equal(log_weights, expected, tolerance = 1e-8)
})

test_that("defensive mixture keeps the Laplace fit and the likelihood", {
  y <- 2
  trials <- 5
  beta <- -0.4
  sigma2 <- 0.8
  parameters <- c(beta, log(sigma2), log(1))
  engine_for <- function(weight) {
    make_adaptive_mcml_engine(
      D = matrix(1, nrow = 1), y = y, units_m = trials,
      distances = matrix(0), family = "binomial", defensive_weight = weight,
      proposal = "laplace"
    )
  }
  laplace_only <- engine_for(0)
  defensive <- engine_for(0.05)
  exact <- log(integrate(
    function(spatial_effect) {
      linear_predictor <- beta + spatial_effect
      exp(y * linear_predictor - trials * softplus(linear_predictor)) *
        dnorm(spatial_effect, sd = sqrt(sigma2))
    },
    lower = -Inf, upper = Inf, rel.tol = 1e-12
  )$value)

  zero <- matrix(0, nrow = 1, ncol = 1)
  expect_identical(defensive$evaluate(parameters, zero),
                   laplace_only$evaluate(parameters, zero))

  set.seed(8128)
  bank <- matrix(rnorm(20000), ncol = 1)
  expect_equal(defensive$evaluate(parameters, bank), exact, tolerance = 0.005)
})

test_that("set_control_mcml validates the defensive weight", {
  expect_identical(set_control_mcml()$defensive_weight, 0.05)
  expect_identical(set_control_mcml(defensive_weight = 0)$defensive_weight, 0)
  expect_error(set_control_mcml(defensive_weight = 1), "defensive_weight")
  expect_error(set_control_mcml(defensive_weight = -0.1), "defensive_weight")
  expect_error(set_control_mcml(defensive_weight = c(0.1, 0.2)),
               "defensive_weight")
})

test_that("stable nonlinear parts equal direct likelihood differences", {
  centre <- c(-3, 0.2, 2)
  step <- outer(c(0.5, 1, 2), seq(-3, 3, length.out = 7))
  y <- c(0, 4, 9)
  units_m <- c(20, 10, 10)
  direct <- function(family) {
    linkf <- make_inverse_link_functions(family, NULL, 3)
    log_likelihood <- function(eta) {
      if (family == "poisson") {
        y * eta - units_m * exp(eta)
      } else {
        y * eta - units_m * log1p(exp(eta))
      }
    }
    terms <- laplace_observation_terms(centre, y, units_m, family, linkf)
    matrix(log_likelihood(as.vector(centre + step)), 3) -
      log_likelihood(centre) - step * terms$score
  }
  for (family in c("binomial", "poisson")) {
    stable <- compute_nonlinear_parts(
      centre, step, y, units_m, family,
      make_inverse_link_functions(family, NULL, 3), canonical = TRUE
    )
    expect_equal(stable$nonlinear, direct(family), tolerance = 1e-10)
  }
  # Very large steps switch to the log-sum-exp form without overflow.
  extreme <- compute_nonlinear_parts(
    0, matrix(c(-800, 800), 1), 1, 10, "binomial",
    make_inverse_link_functions("binomial", NULL, 1), canonical = TRUE
  )
  expect_true(all(is.finite(extreme$nonlinear)))
})

test_that("binomial nonlinear parts stay exact for extreme linear predictors", {
  # At these centres p rounds to 0 or 1, so log(p) and log(1 - p) must come
  # from the linear predictor itself.
  centre <- c(40, -750, 30, -30, 40)
  step <- matrix(c(-50, 800, -35, 35, -10), ncol = 1)
  softplus_stable <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  exact <- -softplus_stable(centre + step) + softplus_stable(centre) +
    step * stats::plogis(centre)
  stable <- compute_nonlinear_parts(
    centre, step, rep(0, 5), rep(1, 5), "binomial",
    make_inverse_link_functions("binomial", NULL, 5), canonical = TRUE
  )
  expect_equal(stable$nonlinear, exact, tolerance = 1e-10)
  expect_equal(stable$nonlinear[1:2], c(-10.0000454, -50), tolerance = 1e-7)
})

test_that("weighted quadratic fits are exact and use QR when needed", {
  quadrature <- build_normal_quadrature()
  nodes <- quadrature$nodes
  expect_equal(sum(quadrature$weights), 1)
  expect_equal(sum(quadrature$weights * nodes^2), 1)
  coefficients <- rbind(c(0.3, -1, 2), c(1, 0.5, 0.1))
  basis <- cbind(1, nodes, -0.5 * nodes^2)
  values <- coefficients %*% t(basis)
  root_weights <- rbind(sqrt(quadrature$weights),
                        sqrt(quadrature$weights * exp(-nodes^2)))
  fitted <- fit_observation_quadratics(values, root_weights, nodes)
  expect_equal(unname(fitted), coefficients, tolerance = 1e-8)
  # The QR refit, used for badly conditioned rows, agrees with the closed
  # form; a minimum condition of one sends every row to QR.
  expect_equal(fit_observation_quadratics(values, root_weights, nodes,
                                          minimum_condition = 1),
               fitted, tolerance = 1e-10)
  # Weight on two nodes only: the fit is not identified.
  concentrated <- matrix(0, 1, length(nodes))
  concentrated[1, 10:11] <- 1
  expect_false(all(is.finite(suppressWarnings(
    fit_observation_quadratics(values[1, , drop = FALSE], concentrated,
                               nodes)[1, ]
  ))))
})

test_that("the global proposal keeps more effective draws than Laplace", {
  # A sparse site: twenty tested and none positive, so the conditional
  # distribution of its latent effect is skewed.
  parameters <- c(-1, log(4), log(1))
  engine_for <- function(proposal) {
    make_adaptive_mcml_engine(
      D = matrix(1, nrow = 1), y = 0, units_m = 20, distances = matrix(0),
      family = "binomial", defensive_weight = 0, proposal = proposal
    )
  }
  set.seed(5)
  bank <- matrix(rnorm(20000), ncol = 1)
  relative_ess <- function(proposal) {
    summarise_importance_weights(engine_for(proposal)$evaluate(
      parameters, bank, log_weights_only = TRUE
    ))$relative_ess
  }
  exact <- log(integrate(function(u) {
    dbinom(0, 20, plogis(-1 + u)) * dnorm(u, sd = 2)
  }, -Inf, Inf, rel.tol = 1e-12)$value)

  # Over 20 banks of this size the relative effective sample size averages
  # 0.94 for the global proposal and 0.70 for the Laplace one.
  expect_gt(relative_ess("global"), relative_ess("laplace") + 0.1)
  expect_gt(relative_ess("global"), 0.9)
  expect_equal(engine_for("global")$evaluate(parameters, bank), exact,
               tolerance = 0.005)
})

test_that("the global proposal works through a random-effects latent map", {
  set.seed(12)
  n <- 40
  coordinate <- runif(n)
  village <- rep(1:8, each = 5)
  trials <- rep(15, n)
  y <- rbinom(n, trials, plogis(-1 + rnorm(8, sd = 0.5)[village]))
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = n), y = y, units_m = trials,
    distances = as.matrix(dist(coordinate)), family = "binomial",
    ID_re = cbind(village = village)
  )
  proposal <- engine$proposal_at(c(-1, log(0.5), log(0.3), log(0.25)))
  expect_identical(proposal$type, "global")
  expect_identical(engine$proposal_failures(), 0L)
  expect_length(proposal$centre, n + 8)
})

test_that("certificate algebra matches its definition", {
  information <- matrix(c(4, 1, 1, 3), 2)
  score <- c(0.2, -0.1)
  covariance <- matrix(c(0.01, 0.002, 0.002, 0.02), 2)
  certificate <- certify_estimate(chol(information), score, covariance)
  expect_equal(certificate$correction, solve(information, score))
  expect_equal(certificate$distance,
               sqrt(sum(score * solve(information, score))))
  expect_equal(certificate$bound, certificate$distance + certificate$radius)
  # With C = c H the radius is sqrt(c times the 95% point of chi-square(2)).
  radius <- compute_validation_radius(chol(information), 0.01 * information)
  expect_equal(radius, sqrt(0.01 * qchisq(0.95, 2)), tolerance = 0.01)
})

test_that("the predicted residual of the correction is that of a Newton step", {
  information <- matrix(c(4, 1, 1, 3), 2)
  centre <- c(0.5, -0.2)
  quadratic <- function(parameters) {
    -0.5 * drop(t(parameters - centre) %*% information %*%
                  (parameters - centre))
  }
  step <- c(0.3, 0.1)
  expect_equal(predict_correction_residual(quadratic, c(0, 0), step, 1:2,
                                           chol(information), step, 0.1),
               0, tolerance = 1e-6)
  # Along exp(x) - x the score at the corrected point is left over.
  curved <- function(parameters) -(exp(parameters[1]) - parameters[1]) -
    parameters[2]^2
  curvature <- diag(c(exp(1), 2))
  step <- c(-0.5, 0)
  left <- exp(1) / 2 - exp(0.5)
  expect_equal(predict_correction_residual(curved, c(1, 0), c(0.5, 0), 1:2,
                                           chol(curvature), step, 0.1),
               abs(left) / sqrt(exp(1)) / 0.1, tolerance = 1e-6)
})

test_that("repeated corrections reach the maximum of the correction bank", {
  log_likelihood <- function(parameters) {
    -(exp(parameters[1]) - 2 * parameters[1]) - (parameters[2] - 1)^2 / 2
  }
  local_mocked_bindings(
    score_validation_bank = function(engine, parameters, free, ...) {
      list(score = c(2 - exp(parameters[1]), 1 - parameters[2])[free])
    }
  )
  repeated <- repeat_correction(NULL, NULL, c(1.2, 1), chol(diag(2)), 1:2,
                                10L, 1L, log_likelihood, c(-10, -10),
                                c(10, 10), radius = 0.01)
  expect_equal(repeated$estimate, c(log(2), 1), tolerance = 1e-8)
  expect_true(repeated$converged)
  expect_true(all(diff(repeated$lengths) < 0))
  # The curvature is the fitting bank's at the last point it was computed.
  expect_equal(repeated$information[2, 2], 1, tolerance = 1e-5)
})

test_that("the repeated correction applies only past its gate", {
  # The correction bank's score is exact, so one Newton step reaches the
  # maximum of a quadratic and falls short along exp(x) - 2x.
  correct <- function(log_likelihood, score, start, correction) {
    local_mocked_bindings(
      score_validation_bank = function(engine, parameters, free, ...) {
        list(score = score(parameters)[free],
             covariance = diag(1e-6, length(free)), log_likelihood = 0,
             seconds = 0)
      }
    )
    information <- compute_bank_curvature(log_likelihood, start, 1:2)
    correct_and_certify(NULL, NULL, start,
                        c(information, list(free = 1:2)), 10L, 1:3,
                        c(-10, -10), c(10, 10), 0.1,
                        bank_objective = log_likelihood,
                        correction = correction)
  }
  quadratic <- function(parameters) -sum((parameters - c(0.5, 1))^2) / 2
  quadratic_score <- function(parameters) c(0.5, 1) - parameters
  single <- correct(quadratic, quadratic_score, c(1.5, 0), "single")
  repeated <- correct(quadratic, quadratic_score, c(1.5, 0), "repeat")
  expect_equal(repeated$estimate, single$estimate)
  expect_identical(repeated$correction_steps, 1L)
  expect_lt(repeated$predicted_residual, 0.25)
  curved <- function(parameters) {
    -(exp(parameters[1]) - 2 * parameters[1]) - (parameters[2] - 1)^2 / 2
  }
  curved_score <- function(parameters) {
    c(2 - exp(parameters[1]), 1 - parameters[2])
  }
  single <- correct(curved, curved_score, c(1.5, 0), "single")
  repeated <- correct(curved, curved_score, c(1.5, 0), "repeat")
  refused <- correct(curved, curved_score, c(1.5, 0), "refuse")
  expect_gt(abs(single$estimate[1] - log(2)), 0.1)
  expect_equal(repeated$estimate, c(log(2), 1), tolerance = 1e-6)
  expect_gt(repeated$correction_steps, 2L)
  expect_true(refused$residual_refused)
  expect_false(refused$certified)
  expect_equal(refused$estimate, single$estimate)
})

test_that("the validation radius leaves the caller's random numbers", {
  set.seed(77)
  before <- .Random.seed
  compute_validation_radius(chol(diag(2)), diag(2))
  expect_identical(.Random.seed, before)
})

test_that("a bank evaluated in blocks of rows gives the same weights", {
  coordinate <- seq(0, 1, length.out = 15)
  trials <- rep(10, 15)
  y <- c(0, 1, 0, 2, 3, 1, 0, 4, 5, 2, 1, 0, 3, 2, 1)
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 15), y = y, units_m = trials,
    distances = as.matrix(dist(coordinate)), family = "binomial"
  )
  parameters <- c(-1.5, log(0.6), log(0.2))
  set.seed(4)
  bank <- matrix(rnorm(1000 * 15), 1000)
  whole <- engine$evaluate(parameters, bank, log_weights_only = TRUE)
  # Blocks of ten rows, most of them without a defensive row, use the share
  # of defensive rows in the whole bank.
  defensive <- defensive_rows(1000, engine$defensive_weight)
  blocks <- unlist(lapply(split(1:1000, rep(1:100, each = 10)), function(rows) {
    engine$evaluate(parameters, bank[rows, , drop = FALSE],
                    log_weights_only = TRUE, defensive = defensive[rows],
                    defensive_share = mean(defensive))
  }), use.names = FALSE)
  expect_equal(blocks, whole)
  # A single block is the bank drawn as one matrix from the same seed.
  score <- score_validation_bank(engine, parameters, 1:3, 1000, seed = 4,
                            block_rows = 1000)
  expect_equal(score$log_likelihood, engine$evaluate(parameters, bank))
  blocked <- score_validation_bank(engine, parameters, 1:3, 1000, seed = 4,
                              block_rows = 100)
  expect_true(all(is.finite(blocked$score)))
  expect_identical(dim(blocked$covariance), c(3L, 3L))
})

test_that("the budget rule targets the tolerance and respects the limit", {
  budgets <- plan_next_budgets(1000L, 5000L, certificate = 0.2,
                               tolerance = 0.1, maximum_budget = 16000L)
  expect_identical(budgets$budget, 5000L)
  expect_identical(budgets$validation_budget, 25000L)
  capped <- plan_next_budgets(1000L, 5000L, certificate = 1,
                              tolerance = 0.1, maximum_budget = 16000L)
  expect_identical(capped$budget, 16000L)
  small <- plan_next_budgets(1000L, 5000L, certificate = 0.101,
                             tolerance = 0.1, maximum_budget = 16000L)
  expect_identical(small$budget, 1500L)
  # Certificates far above the tolerance, or not finite, go to the largest
  # budget instead of overflowing the integer range.
  for (certificate in c(200, Inf, NaN)) {
    largest <- plan_next_budgets(1000L, 5000L, certificate = certificate,
                                 tolerance = 0.1, maximum_budget = 16000L)
    expect_identical(largest$budget, 16000L)
    expect_identical(largest$validation_budget, 80000L)
  }
})

test_that("a failed global proposal is counted and the engine can switch", {
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 3), y = c(0, 1, 2), units_m = c(5, 5, 5),
    distances = as.matrix(dist(1:3)), family = "binomial"
  )
  expect_identical(engine$proposal_type(), "global")
  local_mocked_bindings(fit_global_proposal = function(...) NULL)
  proposal <- engine$proposal_at(c(0, 0, 0))
  expect_identical(proposal$type, "laplace")
  expect_identical(engine$proposal_failures(), 1L)
  engine$reset_proposal_failures()
  engine$use_proposal("laplace")
  expect_identical(engine$proposal_type(), "laplace")
  engine$proposal_at(c(0.1, 0, 0))
  expect_identical(engine$proposal_failures(), 0L)
})

test_that("set_control_mcml validates the certification settings", {
  control <- set_control_mcml()
  expect_identical(control$initial_samples, 1000L)
  expect_identical(control$validation_samples, 5000L)
  expect_identical(control$max_samples, 16000L)
  expect_identical(control$max_stages, 4L)
  expect_identical(control$tolerance, 0.1)
  expect_identical(control$proposal, "global")
  expect_identical(set_control_mcml(proposal = "laplace")$proposal, "laplace")
  expect_error(set_control_mcml(proposal = "mean"), "should be one of")
  expect_error(set_control_mcml(tolerance = 0), "positive number")
  expect_error(set_control_mcml(validation_samples = 0), "positive integer")
})

test_that("the accelerated global iteration reaches the plain fixed point", {
  set.seed(21)
  coordinate <- runif(40)
  trials <- rep(15, 40)
  y <- rbinom(40, trials, plogis(-2 + rnorm(40, sd = 0.8)))
  precision <- solve(0.8 * exp(-as.matrix(dist(coordinate)) / 0.2))
  latent_map <- make_latent_map(cbind(seq_len(40)), 40)
  linkf <- make_inverse_link_functions("binomial", NULL, 40)
  mu <- rep(-2, 40)
  mode <- adaptive_laplace_mode(mu, y, trials, latent_map, precision,
                                "binomial", linkf, numeric(40))$mode
  fit <- function(memory) {
    fit_global_proposal(mode, mu, y, trials, latent_map, precision,
                        "binomial", linkf, canonical = TRUE, memory = memory)
  }
  plain <- fit(0L)
  accelerated <- fit(3L)
  expect_true(plain$converged && accelerated$converged)
  expect_lte(accelerated$iterations, plain$iterations)
  expect_equal(accelerated$centre, plain$centre, tolerance = 1e-6)
  expect_equal(accelerated$precision, plain$precision, tolerance = 1e-6)
})

test_that("dispersed starts span the site spacing and half the extent", {
  layout <- cov_par_layout(2L, fix_tau2 = FALSE)
  distances <- as.matrix(dist(c(0, 1, 2, 3, 10)))
  dispersed <- make_dispersed_starts(layout, distances)
  expect_equal(exp(dispersed$log_phi), c(1, 1, 5, 5))
  expect_equal(exp(dispersed$log_sigma2), c(0.25, 4, 0.25, 4))
  expect_equal(dispersed$trigger, log(0.5))
  expect_identical(dispersed$phi_index, 4L)
})

test_that("the check for maxima finds a higher one at a longer scale", {
  # Two maxima in the log scale: a local one at -3 and a higher one at 2.
  objective <- function(parameters, z) {
    log(exp(-(parameters[2] + 3)^2) + 3 * exp(-(parameters[2] - 2)^2)) -
      parameters[1]^2
  }
  # A stand-in engine whose estimated log-likelihood on any bank is exactly
  # the objective, so the independent comparison has no Monte Carlo error.
  engine <- list(latent_dimension = 1L,
                 evaluate = function(parameters, z, log_weights_only) {
                   rep(objective(parameters, z), nrow(z))
                 })
  announced <- NULL
  limits <- list(variance_index = integer(), linked_index = list(),
                 boundary_log_variance = -9, tolerance = 0.1,
                 dispersed = list(sigma2_index = 1L, phi_index = 2L,
                                  log_sigma2 = log(c(0.25, 4, 0.25, 4)),
                                  log_phi = log(c(1, 1, 50, 50)),
                                  labels = c("a", "b", "c", "d"),
                                  trigger = log(0.5),
                                  announce = function(log_phi) {
                                    announced <<- log_phi
                                  }))
  lower <- c(-10, -10)
  upper <- c(10, 10)
  local <- optimise_on_bank(objective, c(0, -3.2), NULL, lower, upper,
                            TRUE, limits)
  expect_equal(local$fit$par[2], -3, tolerance = 1e-3)
  checked <- check_other_maxima(engine, objective, NULL, local, NULL,
                                    lower, upper, limits, 100L, 1L)
  expect_equal(announced, local$fit$par[2])
  expect_true(checked$separate)
  expect_false(checked$ambiguous)
  expect_false(checked$check$chosen[1])
  expect_equal(checked$optimised$fit$par[2], 2, tolerance = 1e-3)
  expect_gt(checked$check$distance_from_fitted[checked$check$chosen], 0.1)

  announced <- NULL
  higher <- optimise_on_bank(objective, c(0, 1.5), NULL, lower, upper, TRUE,
                             limits)
  unchecked <- check_other_maxima(engine, objective, NULL, higher, NULL,
                                      lower, upper, limits, 100L, 1L)
  expect_null(unchecked$check)
  expect_null(announced)
})

test_that("a stage whose optimiser returns non-finite values has no valid fit", {
  # The Laplace fallback on a ridge returned NaN for every parameter with a
  # convergence code of zero, and the check for other maxima then failed on
  # the comparison of the scale with its trigger.
  engine <- list(latent_dimension = 1L,
                 proposal_type = function() "laplace",
                 proposal_failures = function() 0L,
                 evaluate = function(...) stop("evaluated at an invalid fit"))
  limits <- list(dispersed = list(sigma2_index = 1L, phi_index = 2L,
                                  trigger = log(0.5)))
  local_mocked_bindings(
    optimise_on_bank = function(...) {
      list(fit = list(par = c(NaN, NaN), objective = 1e100, convergence = 0L),
           boundary = integer(), restarted = FALSE, restart_improved = FALSE)
    },
    check_other_maxima = function(...) stop("checked an invalid fit")
  )
  stage <- run_certified_stage(engine, function(...) 0, NULL, c(0, 0), 10L,
                               10L, 1:4, c(-10, -10), c(10, 10), TRUE, limits)
  expect_identical(stage$outcome, "no valid fit")
  expect_match(stage$reason, "non-finite")
  expect_null(stage$maxima_check)
  expect_null(stage$fit_diagnostics)
})

test_that("the check for other maxima leaves out refits that fail", {
  objective <- function(parameters, z) {
    log(exp(-(parameters[2] + 3)^2) + 3 * exp(-(parameters[2] - 2)^2)) -
      parameters[1]^2
  }
  engine <- list(latent_dimension = 1L,
                 evaluate = function(parameters, z, log_weights_only) {
                   rep(objective(parameters, z), nrow(z))
                 })
  limits <- list(variance_index = integer(), linked_index = list(),
                 boundary_log_variance = -9, tolerance = 0.1,
                 dispersed = list(sigma2_index = 1L, phi_index = 2L,
                                  log_sigma2 = log(c(0.25, 4, 0.25, 4)),
                                  log_phi = log(c(1, 1, 50, 50)),
                                  labels = c("a", "b", "c", "d"),
                                  trigger = log(0.5)))
  lower <- c(-10, -10)
  upper <- c(10, 10)
  local <- optimise_on_bank(objective, c(0, -3.2), NULL, lower, upper,
                            TRUE, limits)
  original <- optimise_on_bank
  # Refits from the long-scale starts fail as the Laplace fallback did.
  local_mocked_bindings(optimise_on_bank = function(objective, start, ...) {
    if (start[2] > 3) {
      return(list(fit = list(par = c(NaN, NaN), objective = 1e100),
                  boundary = integer()))
    }
    original(objective, start, ...)
  })
  checked <- check_other_maxima(engine, objective, NULL, local, NULL,
                                    lower, upper, limits, 100L, 1L)
  expect_identical(checked$check$start, c("fitted", "a", "b"))
  expect_true(all(is.finite(checked$check$log_likelihood)))
  expect_true(all(is.finite(checked$optimised$fit$par)))
})

test_that("candidate maxima are tied, separate or ambiguous", {
  engine <- list(latent_dimension = 1L,
                 evaluate = function(parameters, z, log_weights_only) {
                   # Equal likelihoods at the two separated points, with
                   # Monte Carlo noise from the bank.
                   -abs(parameters[1]) * 0 + z[, 1] * 0.01
                 })
  information <- diag(100, 2)
  # Within 0.1 standard errors (0.01 apart with information 100): one maximum.
  tied <- choose_maximum(engine, NULL, list(c(0, 0), c(0.005, 0)),
                         c(-10, -9.999), information, 50L, 1L, 0.1)
  expect_false(tied$separate)
  expect_false(tied$ambiguous)
  expect_identical(tied$chosen, 2L)
  expect_true(all(is.na(tied$table$independent_log_likelihood)))
  # Far apart with equal likelihoods on independent draws: ambiguous.
  separated <- choose_maximum(engine, NULL, list(c(0, 0), c(1, 0)),
                              c(-10, -9.999), information, 50L, 1L, 0.1)
  expect_true(separated$ambiguous)
  expect_identical(sum(is.finite(separated$table$independent_log_likelihood)),
                   2L)
})

test_that("set_control_mcml validates check_maxima", {
  expect_true(set_control_mcml()$check_maxima)
  expect_false(set_control_mcml(check_maxima = FALSE)$check_maxima)
  expect_error(set_control_mcml(check_maxima = NA), "check_maxima")
})

test_that("Newton steps with a known curvature reach the maximum", {
  curvature <- matrix(c(4, 1, 1, 3), 2)
  centre <- c(0.7, -1.2)
  negative <- function(parameters) {
    0.5 * sum((parameters - centre) * (curvature %*% (parameters - centre)))
  }
  exact <- run_newton_steps(negative, c(0, 0), chol(curvature), c(-5, -5),
                          c(5, 5))
  expect_true(exact$converged)
  expect_equal(exact$par, centre, tolerance = 1e-6)
  # A badly scaled curvature can stall the steps; the optimisation as a
  # whole still reaches the maximum because nlminb() finishes it.
  limits <- list(variance_index = integer(), linked_index = list(),
                 boundary_log_variance = -9)
  rough <- optimise_on_bank(function(parameters, z) -negative(parameters),
                            c(0, 0), NULL, c(-5, -5), c(5, 5), FALSE, limits,
                            curvature_root = chol(diag(c(1, 50))))
  expect_equal(rough$fit$par, centre, tolerance = 1e-5)
  # Steps are kept within the bounds.
  bounded <- run_newton_steps(negative, c(0, 0), chol(curvature), c(-5, -1),
                            c(5, 5))
  expect_gte(min(bounded$par[2]), -1)
  expect_false(bounded$converged)
})

test_that("the optimisation finishes with nlminb after Newton steps", {
  negative_objective <- function(parameters, z) {
    -sum((parameters - c(1, 2))^2) - 3
  }
  limits <- list(variance_index = integer(), linked_index = list(),
                 boundary_log_variance = -9)
  plain <- optimise_on_bank(negative_objective, c(0, 0), NULL, c(-5, -5),
                            c(5, 5), FALSE, limits)
  seeded <- optimise_on_bank(negative_objective, c(0, 0), NULL, c(-5, -5),
                             c(5, 5), FALSE, limits,
                             curvature_root = chol(diag(2, 2)))
  expect_identical(plain$fit$method, "nlminb")
  expect_identical(seeded$fit$method, "newton, then nlminb")
  expect_equal(seeded$fit$par, c(1, 2), tolerance = 1e-6)
  expect_equal(plain$fit$par, seeded$fit$par, tolerance = 1e-6)
})

test_that("curvatures are reused only when well conditioned", {
  expect_null(factor_curvature(NULL))
  expect_null(factor_curvature(diag(c(1, 1e-12))))
  expect_null(factor_curvature(matrix(c(1, 2, 2, 1), 2)))
  expect_equal(factor_curvature(diag(c(4, 9))), diag(c(2, 3)))
})

test_that("the stopping tolerance is fixed in log-likelihood units", {
  expect_identical(make_optimiser_control()$rel.tol, 1e-10)
  expect_equal(make_optimiser_control(1000)$rel.tol, 1e-7)
  expect_equal(make_optimiser_control(10)$rel.tol, 1e-6)
  expect_equal(make_optimiser_control(1e9)$rel.tol, 1e-12)
  expect_identical(make_optimiser_control(Inf)$rel.tol, 1e-10)
})

test_that("the one-round central Hessian matches an exact quadratic", {
  curvature <- matrix(c(3, -1, 0.5, -1, 2, 0.2, 0.5, 0.2, 1), 3)
  objective <- function(point) {
    0.5 * sum(point * (curvature %*% point)) + sum(point^3) / 6
  }
  point <- c(0.3, -0.4, 1.2)
  expected <- curvature + diag(point)
  expect_equal(compute_central_hessian(objective, point), expected,
               tolerance = 1e-6)
  both <- compute_central_hessian(objective, point, with_gradient = TRUE)
  expect_identical(both$hessian, compute_central_hessian(objective, point))
  expect_equal(both$gradient, drop(curvature %*% point) + point^2 / 2,
               tolerance = 1e-6)
})

test_that("weight diagnostics survive degenerate banks", {
  none <- summarise_importance_weights(rep(-Inf, 10))
  expect_identical(none$log_likelihood, -Inf)
  expect_identical(none$relative_ess, 0)
  single <- summarise_importance_weights(c(-Inf, 0))
  expect_equal(single$log_likelihood, log(0.5))
  expect_true(is.na(single$pareto_k) || is.finite(single$pareto_k))
})

test_that("a failed evaluation gives every draw zero weight", {
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 3), y = c(0, 1, 2), units_m = c(5, 5, 5),
    distances = as.matrix(dist(1:3)), family = "binomial"
  )
  local_mocked_bindings(
    cpp_spatial_prior_factors = function(...) stop("not positive")
  )
  bank <- matrix(rnorm(12), 4)
  expect_identical(engine$evaluate(c(0, 0, 0), bank, log_weights_only = TRUE),
                   rep(-Inf, 4))
  expect_identical(engine$evaluate(c(0, 0, 0), bank), -1e100)
})

test_that("the anchored start reaches the proposal of the start at the mode", {
  set.seed(14)
  coordinate <- runif(30)
  trials <- rep(20, 30)
  y <- rbinom(30, trials, plogis(-1.5 + rnorm(30, sd = 0.7)))
  engine_for <- function() {
    make_adaptive_mcml_engine(
      D = matrix(1, nrow = 30), y = y, units_m = trials,
      distances = as.matrix(dist(coordinate)), family = "binomial"
    )
  }
  anchored <- engine_for()
  expect_true(anchored$set_anchor(c(-1.4, log(0.5), log(0.2))))
  plain <- engine_for()
  parameters <- c(-1.6, log(0.6), log(0.25))
  expect_equal(anchored$proposal_at(parameters)$centre,
               plain$proposal_at(parameters)$centre, tolerance = 1e-6)
  expect_equal(anchored$proposal_at(parameters)$root,
               plain$proposal_at(parameters)$root, tolerance = 1e-6)
})

test_that("the check for other maxima reports a move and an ambiguity", {
  layout <- cov_par_layout(1L, fix_tau2 = FALSE)
  check <- data.frame(
    log_likelihood = c(-12, -10), distance_from_fitted = c(0, 3),
    maximum = 1:2, independent_log_likelihood = c(-12.1, -10.2),
    independent_se = c(0.05, 0), chosen = c(FALSE, TRUE),
    start = c("fitted", "a"), log_phi = log(c(0.1, 20))
  )
  moved <- list(maxima_check = check, better_maximum = TRUE,
                ambiguous_maxima = FALSE)
  expect_warning(report_maxima_check(moved, layout, 1),
                 "found a separate maximum.*scale 20 instead of 0.1")
  ambiguous <- moved
  ambiguous$better_maximum <- FALSE
  ambiguous$ambiguous_maxima <- TRUE
  expect_warning(report_maxima_check(ambiguous, layout, 1),
                 "fit about equally well")
  quiet <- moved
  quiet$better_maximum <- FALSE
  expect_silent(report_maxima_check(quiet, layout, 1))
})

test_that("the objective carries the PC penalty only when asked", {
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 4), y = c(0, 1, 2, 1), units_m = rep(5, 4),
    distances = as.matrix(dist(1:4)), family = "binomial"
  )
  plain <- make_mcml_objective(engine, set_control_mcml(),
                               as.matrix(dist(1:4)), 1)
  expect_null(plain$penalty)
  penalised <- make_mcml_objective(engine, set_control_mcml(pc_prior = TRUE),
                                   as.matrix(dist(1:4)), 1)
  zero <- matrix(0, 1, 4)
  parameters <- c(-1, 0, 0)
  expect_equal(penalised$objective(parameters, zero),
               engine$evaluate(parameters, zero) +
                 penalised$penalty(parameters))
})

test_that("the table of the check for other maxima is restored to user units", {
  history <- list(list(estimate = c(0, 0, 1, log(2)),
                       maxima_check = data.frame(log_phi = log(c(0.1, 3)))))
  restored <- restore_mcml_history_scale(history, diag(2), distance_scale = 10)
  expect_equal(restored[[1]]$maxima_check$log_phi, log(c(1, 30)))
  expect_equal(restored[[1]]$estimate[4], log(20))
})

# The R computations that the compiled routines of the engine replaced. The
# compiled versions must reproduce them exactly.
reference_inverse_cholesky_factor <- function(hessian) {
  dimension <- nrow(hessian)
  reversed <- dimension:1
  root <- chol(hessian[reversed, reversed, drop = FALSE])
  inverse_root <- backsolve(root, diag(dimension))
  t(inverse_root)[reversed, reversed, drop = FALSE]
}

make_compiled_test_problem <- function(with_random_effect) {
  set.seed(31)
  locations <- 12
  coordinates <- matrix(runif(2 * locations), ncol = 2)
  precision <- chol2inv(chol(exp(-as.matrix(dist(coordinates)) / 0.3)))
  levels <- 4
  dimension <- locations + if (with_random_effect) levels else 0
  if (with_random_effect) {
    full <- matrix(0, dimension, dimension)
    full[seq_len(locations), seq_len(locations)] <- precision
    full[cbind(locations + 1:levels, locations + 1:levels)] <- 2
    precision <- full
  }
  observations <- 20
  columns <- cbind(sample(locations, observations, replace = TRUE))
  if (with_random_effect) {
    columns <- cbind(columns, locations + sample(levels, observations,
                                                 replace = TRUE))
  }
  list(precision = precision, latent_map = make_latent_map(columns, dimension),
       linear = rnorm(observations), curvature = runif(observations, 0.1, 2))
}

test_that("compiled marginals of the global iteration equal the R code", {
  for (with_random_effect in c(FALSE, TRUE)) {
    problem <- make_compiled_test_problem(with_random_effect)
    latent_map <- problem$latent_map
    proposal_precision <- problem$precision +
      latent_data_hessian(latent_map, problem$curvature)
    inverse <- chol2inv(chol(proposal_precision))
    centre <- as.numeric(inverse %*%
                           latent_gradient(latent_map, problem$linear))
    compiled <- cpp_global_marginals(problem$precision, problem$linear,
                                     problem$curvature, latent_map$columns)
    expect_identical(compiled$marginal_mean,
                     latent_linear_predictor(latent_map, centre))
    expect_identical(compiled$marginal_variance, as.numeric(rowsum(
      inverse[latent_map$pair_linear], latent_map$pair_observation,
      reorder = FALSE
    )))
  }
  expect_null(cpp_global_marginals(-diag(3), numeric(3), numeric(3),
                                   cbind(1:3)))
})

test_that("the compiled converged proposal equals the R code", {
  for (with_random_effect in c(FALSE, TRUE)) {
    problem <- make_compiled_test_problem(with_random_effect)
    latent_map <- problem$latent_map
    proposal_precision <- problem$precision +
      latent_data_hessian(latent_map, problem$curvature)
    proposal_precision <- (proposal_precision + t(proposal_precision)) / 2
    root <- chol(proposal_precision)
    compiled <- cpp_global_proposal_factors(
      problem$precision, problem$linear, problem$curvature, latent_map$columns
    )
    expect_identical(compiled$precision, proposal_precision)
    expect_identical(compiled$centre, as.numeric(backsolve(root, forwardsolve(
      t(root), latent_gradient(latent_map, problem$linear)
    ))))
    expect_identical(compiled$root,
                     reference_inverse_cholesky_factor(proposal_precision))
    expect_identical(inverse_cholesky_factor(proposal_precision),
                     compiled$root)
  }
  expect_null(cpp_global_proposal_factors(-diag(3), numeric(3), numeric(3),
                                          cbind(1:3)))
  expect_error(inverse_cholesky_factor(-diag(3)), "not positive definite")
})

test_that("the compiled spatial prior equals the R code and defers jitter", {
  set.seed(32)
  distances <- as.matrix(dist(matrix(runif(20), ncol = 2)))
  for (kappa in c(0.5, 1.5, 2.5)) {
    compiled <- cpp_spatial_prior_factors(as.numeric(distances), 10L, 0.4,
                                          kappa, 1.7, 0.2)
    correlation <- matrix(cpp_half_integer_matern(as.numeric(distances), 0.4,
                                                  kappa), 10, 10)
    diag(correlation) <- diag(correlation) + 0.2
    covariance <- 1.7 * correlation
    root <- factor_covariance(covariance)
    expect_identical(compiled$covariance, covariance)
    expect_identical(compiled$root, root)
    expect_identical(compiled$precision, chol2inv(root))
  }
  # Two sites at the same place without a nugget give a singular covariance,
  # which the engine factorises with jitter as before.
  repeated <- as.matrix(dist(c(0, 0, 1)))
  singular <- cpp_spatial_prior_factors(as.numeric(repeated), 3L, 1, 0.5, 1, 0)
  expect_null(singular$root)
  expect_identical(singular$covariance, unname(exp(-repeated)))
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 3), y = c(1, 2, 1), units_m = rep(4, 3),
    distances = repeated, family = "binomial", proposal = "laplace"
  )
  expect_warning(proposal <- engine$proposal_at(c(0, 0, 0)), "jitter")
  expect_gt(attr(proposal$spatial_root, "jitter"), 0)
})

test_that("compiled importance weights equal the R computation", {
  set.seed(33)
  locations <- 15
  distances <- as.matrix(dist(matrix(runif(2 * locations), ncol = 2)))
  cluster <- sample(4, locations, replace = TRUE)
  design <- cbind(1, rnorm(locations))
  linear_predictor <- -0.5 + 0.3 * design[, 2] + rnorm(locations, sd = 0.5)
  cases <- list(
    list(family = "binomial", y = rbinom(locations, 10, plogis(linear_predictor)),
         units = rep(10, locations), random = NULL, defensive = 0.05),
    list(family = "binomial", y = rbinom(locations, 10, plogis(linear_predictor)),
         units = rep(10, locations), random = cbind(cluster), defensive = 0.05),
    list(family = "poisson", y = rpois(locations, exp(linear_predictor + 1)),
         units = rep(1, locations), random = cbind(cluster), defensive = 0.05),
    list(family = "binomial", y = rbinom(locations, 10, plogis(linear_predictor)),
         units = rep(10, locations), random = NULL, defensive = 0)
  )
  for (case in cases) {
    engine <- make_adaptive_mcml_engine(
      D = design, y = case$y, units_m = case$units, distances = distances,
      family = case$family, ID_re = case$random,
      defensive_weight = case$defensive
    )
    scope <- environment(engine$evaluate)
    parameters <- c(-0.5, 0.3, log(0.5), log(0.3),
                    if (!is.null(case$random)) log(0.4))
    for (shift in c(0, 0.05)) for (rows in c(1L, 40L)) {
      target <- parameters + c(0, 0, shift, -shift,
                               if (!is.null(case$random)) shift)
      draws <- matrix(rnorm(rows * engine$latent_dimension), rows)
      proposal <- engine$proposal_at(parameters)
      same <- identical(target, parameters)
      values <- if (same) proposal$values else scope$unpack(target)
      spatial_root <- if (same) proposal$spatial_root else
        factor_covariance(scope$covariance_at(values))
      defensive <- defensive_rows(rows, case$defensive)
      in_r <- compute_log_weights_in_r(
        draws, proposal, defensive, mean(defensive), values, spatial_root,
        same, as.numeric(design %*% values$beta), scope$latent,
        scope$latent_map$columns,
        function(eta) {
          compute_observation_log_likelihood(eta, case$y, case$units,
                                             case$family, scope$linkf)
        }
      )
      # The R code can carry a random effect's name on a single-row result.
      # The weights are identical where the compiler does not fuse multiply
      # and add, as on x86-64. Where it does, as on arm64, the compiled sums
      # differ from R's in the last bits.
      expect_equal(
        engine$evaluate(target, draws, log_weights_only = TRUE,
                        proposal_parameters = parameters),
        unname(in_r), tolerance = 1e-12
      )
    }
  }
})

test_that("a late proposal failure repeats the stage on fresh streams", {
  failures <- 0L
  type <- "global"
  engine <- list(
    reset_proposal_failures = function() failures <<- 0L,
    proposal_failures = function() failures,
    proposal_type = function() type,
    use_proposal = function(new) type <<- new
  )
  calls <- list()
  scenario <- "late"
  local_mocked_bindings(
    run_certified_stage = function(engine, objective, penalty, start, budget,
                                   validation_budget, seeds, ...) {
      calls[[length(calls) + 1L]] <<- seeds
      global <- identical(type, "global")
      before <- if (global && scenario == "early") 1L else 0L
      failures <<- before + if (global && scenario == "late") 1L else 0L
      list(estimate = start, outcome = "certified", certified = TRUE,
           certificate = 0.01, free = seq_along(start),
           information = diag(length(start)),
           information_root = diag(length(start)),
           failures_before_validation = before)
    }
  )
  run <- function(which) {
    scenario <<- which
    calls <<- list()
    type <<- "global"
    fit_adaptive_mcml(engine, start = c(0, 0, 0), seed = 7L)
  }
  late <- run("late")
  expect_length(calls, 2L)
  expect_true(all(calls[[1]] != calls[[2]]))
  expect_true(late$history[[1]]$late_fallback)
  expect_identical(late$proposal, "laplace")
  early <- run("early")
  expect_length(calls, 2L)
  expect_identical(calls[[1]], calls[[2]])
  expect_false(early$history[[1]]$late_fallback)
  none <- run("none")
  expect_length(calls, 1L)
  expect_null(none$history[[1]]$proposal_fallback)
})

test_that("a stage whose global proposal fails at the start ends at once", {
  failures <- 0L
  engine <- list(latent_dimension = 1L,
                 proposal_type = function() "global",
                 set_anchor = function(parameters) {
                   failures <<- failures + 1L
                   invisible(FALSE)
                 },
                 proposal_failures = function() failures)
  local_mocked_bindings(
    optimise_on_bank = function(...) stop("the condemned stage was fitted")
  )
  stage <- run_certified_stage(engine, function(...) 0, NULL, c(0, 0), 10L,
                               10L, 1:5, c(-10, -10), c(10, 10), TRUE, list())
  expect_identical(stage$outcome, "no valid fit")
  expect_match(stage$reason, "failed at the start")
  expect_null(stage$failures_before_validation)
})

test_that("reusing the spatial factorisation leaves the likelihood unchanged", {
  make <- function() {
    make_adaptive_mcml_engine(
      D = cbind(1, c(-1, 0.5, 2, -0.3, 1.1)), y = c(0, 3, 5, 1, 2),
      units_m = rep(6, 5), distances = as.matrix(dist(cbind(1:5, c(2, 1, 4, 3, 5)))),
      family = "binomial"
    )
  }
  set.seed(4)
  bank <- matrix(stats::rnorm(50 * 5), nrow = 50)
  first <- c(-0.5, 0.2, 0, log(2))
  second <- replace(first, 2, 0.4)
  reused <- make()
  invisible(reused$evaluate(first, bank))
  # The engine also starts each mode search from the previous mode, so the
  # two paths agree to rounding, not to the last bit.
  expect_equal(reused$evaluate(second, bank), make()$evaluate(second, bank),
               tolerance = 1e-10)
})
