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
  small <- small_spatial_scale(fit, layout, distances, kappa = 0.5,
                               distance_scale = 2)
  expect_equal(small, c(scale = 0.6, practical_range = 0.6 * -log(0.05),
                        spacing = 2), tolerance = 1e-8)
  # With the spatial variance at zero the scale is meaningless and is not
  # reported as small.
  fit$boundary_parameters <- 3:4
  expect_null(small_spatial_scale(fit, layout, distances, 0.5, 2))

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
  expect_equal(practical_range_factor(0.5), -log(0.05), tolerance = 1e-8)
  expect_gt(practical_range_factor(1.5), practical_range_factor(0.5))
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

test_that("defensive mixture keeps the Laplace approximation and the likelihood", {
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

test_that("stable nonlinear parts equal the direct log-likelihood differences", {
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
    stable <- observation_nonlinear_part(
      centre, step, y, units_m, family,
      make_inverse_link_functions(family, NULL, 3), canonical = TRUE
    )
    expect_equal(stable$nonlinear, direct(family), tolerance = 1e-10)
  }
  # Very large steps switch to the log-sum-exp form without overflow.
  extreme <- observation_nonlinear_part(
    0, matrix(c(-800, 800), 1), 1, 10, "binomial",
    make_inverse_link_functions("binomial", NULL, 1), canonical = TRUE
  )
  expect_true(all(is.finite(extreme$nonlinear)))
})

test_that("weighted quadratic fits recover exact quadratics and use QR when needed", {
  quadrature <- standard_normal_quadrature()
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

test_that("the global proposal keeps more effective draws than the Laplace one", {
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
    importance_diagnostics(engine_for(proposal)$evaluate(
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

test_that("the global proposal works through a latent map with random effects", {
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
  certificate <- mcml_certificate(chol(information), score, covariance)
  expect_equal(certificate$correction, solve(information, score))
  expect_equal(certificate$distance,
               sqrt(sum(score * solve(information, score))))
  expect_equal(certificate$bound, certificate$distance + certificate$radius)
  # With C = c H the radius is sqrt(c times the 95% point of chi-square(2)).
  radius <- validation_radius(chol(information), 0.01 * information)
  expect_equal(radius, sqrt(0.01 * qchisq(0.95, 2)), tolerance = 0.01)
})

test_that("the validation radius leaves the caller's random numbers unchanged", {
  set.seed(77)
  before <- .Random.seed
  validation_radius(chol(diag(2)), diag(2))
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
  score <- validation_score(engine, parameters, 1:3, 1000, seed = 4,
                            block_rows = 1000)
  expect_equal(score$log_likelihood, engine$evaluate(parameters, bank))
  blocked <- validation_score(engine, parameters, 1:3, 1000, seed = 4,
                              block_rows = 100)
  expect_true(all(is.finite(blocked$score)))
  expect_identical(dim(blocked$covariance), c(3L, 3L))
})

test_that("the budget rule targets the tolerance and respects the limit", {
  budgets <- next_mcml_budgets(1000L, 5000L, certificate = 0.2,
                               tolerance = 0.1, maximum_budget = 16000L)
  expect_identical(budgets$budget, 5000L)
  expect_identical(budgets$validation_budget, 25000L)
  capped <- next_mcml_budgets(1000L, 5000L, certificate = 1,
                              tolerance = 0.1, maximum_budget = 16000L)
  expect_identical(capped$budget, 16000L)
  small <- next_mcml_budgets(1000L, 5000L, certificate = 0.101,
                             tolerance = 0.1, maximum_budget = 16000L)
  expect_identical(small$budget, 1500L)
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
  dispersed <- dispersed_starts(layout, distances)
  expect_equal(exp(dispersed$log_phi), c(1, 1, 5, 5))
  expect_equal(exp(dispersed$log_sigma2), c(0.25, 4, 0.25, 4))
  expect_equal(dispersed$trigger, log(0.5))
  expect_identical(dispersed$phi_index, 4L)
})

test_that("the optimisation check finds a higher maximum at a longer scale", {
  # Two maxima in the log scale: a local one at -3 and a higher one at 2.
  objective <- function(parameters, z) {
    log(exp(-(parameters[2] + 3)^2) + 3 * exp(-(parameters[2] - 2)^2)) -
      parameters[1]^2
  }
  announced <- NULL
  limits <- list(variance_index = integer(), linked_index = list(),
                 boundary_log_variance = -9,
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
  checked <- check_dispersed_maxima(objective, local, NULL, lower, upper,
                                    limits)
  expect_equal(announced, local$fit$par[2])
  expect_false(checked$check$chosen[1])
  expect_equal(checked$optimised$fit$par[2], 2, tolerance = 1e-3)

  announced <- NULL
  higher <- optimise_on_bank(objective, c(0, 1.5), NULL, lower, upper, TRUE,
                             limits)
  unchecked <- check_dispersed_maxima(objective, higher, NULL, lower, upper,
                                      limits)
  expect_null(unchecked$check)
  expect_null(announced)
})

test_that("set_control_mcml validates optimisation_check", {
  expect_true(set_control_mcml()$optimisation_check)
  expect_false(set_control_mcml(optimisation_check = FALSE)$optimisation_check)
  expect_error(set_control_mcml(optimisation_check = NA), "optimisation_check")
})
