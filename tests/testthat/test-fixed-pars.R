# fixed_pars holds the nugget, the measurement error variance and the size
# of a negative binomial family (step 2d of the parameter table, P3, P4 and
# P6 of the interface review), in place of gp(nugget = <value>), fix_var_me
# and negative_binomial(theta).

test_that("fixed_pars holds the nugget and the measurement error variance", {
  fit <- glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data,
               family = "gaussian", fixed_pars = list(tau2 = 0.2),
               messages = FALSE)
  expect_identical(fit$nugget, "fixed")
  expect_identical(fit$fixed_pars, list(tau2 = 0.2))
  expect_null(fit$estimate$nu2)
  expect_identical(fit_fix_tau2(fit), 0.2)
  expect_named(coef(fit), c("beta", "sigma2", "phi", "tau2", "sigma2_me"))
  expect_identical(coef(fit)$tau2, 0.2)
  expect_false("tau2" %in% names(estimated_coef(fit)))
  expect_identical(summary(fit)$tau2, 0.2)
  expect_output(print(fit), "Variance of the nugget effect fixed at 0.2")
  # A model without a nugget, or with one estimated.
  expect_identical(gaussian_model$nugget, "none")
  expect_null(gaussian_model$fixed_pars)
  expect_identical(fit_fix_tau2(gaussian_model), 0)
  expect_identical(gaussian_offset_model$nugget, "estimated")
  expect_identical(gaussian_offset_model$fixed_pars, list(sigma2_me = 0))
  expect_identical(fit_fix_var_me(gaussian_offset_model), 0)
  # The refits of a fit hold what it held.
  expect_identical(fit_fixed_pars(fit), list(tau2 = 0.2))
  expect_identical(
    assessment_refit_args(fit, fit$formula, gaussian_data, control_mcmc)$fixed_pars,
    list(tau2 = 0.2)
  )
})

test_that("holds that don't fit the model stop", {
  holds <- list(
    list(tau2 = 0), list(tau2 = -1), list(sigma2_me = -0.1), list(phi = -1),
    list(beta = c(age = 1)), list(theta = 2), list(tau2 = c(1, 2)),
    list(sigma2_me = NA), list(sigma2 = 0), list(beta = c(cov = Inf))
  )
  for (fixed_pars in holds) {
    expect_input_error(
      glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data,
            family = "gaussian", fixed_pars = fixed_pars, messages = FALSE),
      "fixed_pars"
    )
  }
  # A nugget the model doesn't have can't be held.
  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          fixed_pars = list(tau2 = 0.2), messages = FALSE),
    "fixed_pars"
  )
  # A parameter can't be started and held, unless the start is its value.
  expect_input_error(
    glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data,
          family = "gaussian", fixed_pars = list(tau2 = 0.2),
          start_pars = list(tau2 = 0.3), messages = FALSE),
    "start_pars"
  )
  expect_no_error(
    glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data,
          family = "gaussian", fixed_pars = list(tau2 = 0.2),
          start_pars = list(tau2 = 0.2), messages = FALSE)
  )
  # The old forms of a hold stop with a pointer to fixed_pars.
  error <- expect_input_error(gp(nugget = 0.5), "nugget")
  expect_match(conditionMessage(error), "fixed_pars = list(tau2 = 0.5)",
               fixed = TRUE)
  expect_input_error(negative_binomial(theta = 2), "theta")
  expect_error(glgpm(y ~ gp(), data = gaussian_data, family = "gaussian",
                     fix_var_me = 0.1),
               "unused argument")
})

test_that("fixed_pars holds the size of a negative binomial family", {
  counts <- function(...) {
    suppressWarnings(glgpm(
      y ~ cov + gp(), data = poisson_data, denominator = denominator,
      control_mcmc = control_mcmc, seed = fixture_seed,
      control_mcml = fixture_control_mcml, return_samples = FALSE,
      messages = FALSE, ...
    ))
  }
  held <- counts(family = negative_binomial(), fixed_pars = list(theta = 2))
  expect_identical(held$fixed_pars, list(theta = 2))
  expect_identical(fit_family_held(held), list(theta = 2))
  expect_identical(held_family_parameters(held), list(theta = 2))
  expect_identical(estimated_family_parameters(held), character())
  expect_identical(coef(held)$theta, 2)
  # As MASS::negative.binomial(2) holds it.
  mass <- counts(family = held_negative_binomial(2))
  expect_identical(mass$estimate, held$estimate)
  expect_identical(mass$fixed_pars, list(theta = 2))
  expect_no_error(counts(family = held_negative_binomial(2),
                         fixed_pars = list(theta = 2)))
  expect_input_error(counts(family = held_negative_binomial(2),
                            fixed_pars = list(theta = 3)),
                     "fixed_pars")
  # A model that simulates from the fit, and specify_glgpm(), take the size
  # in their parameters.
  expect_identical(simulation_model(held)$family_parameters, list(theta = 2))
  model <- specify_glgpm(
    y ~ cov + gp(), data = poisson_data, family = negative_binomial(),
    parameters = list(beta = c(0.1, 0.2), sigma2 = 0.5, phi = 2, theta = 3),
    denominator = denominator
  )
  expect_identical(model$family_parameters, list(theta = 3))
  expect_input_error(
    specify_glgpm(y ~ cov + gp(), data = poisson_data,
                  family = negative_binomial(),
                  parameters = list(beta = c(0.1, 0.2), sigma2 = 0.5, phi = 2),
                  denominator = denominator),
    "family"
  )
  expect_input_error(
    specify_glgpm(y ~ cov + gp(), data = poisson_data,
                  family = held_negative_binomial(2),
                  parameters = list(beta = c(0.1, 0.2), sigma2 = 0.5, phi = 2,
                                    theta = 3),
                  denominator = denominator),
    "parameters"
  )
})

test_that("fits saved before fixed_pars are read through their old fields", {
  saved <- gaussian_offset_model
  saved$nugget <- NULL
  saved$fixed_pars <- NULL
  saved$fix_tau2 <- TRUE
  saved$fix_var_me <- 0
  expect_identical(fit_fix_tau2(saved), TRUE)
  expect_identical(fit_fix_var_me(saved), 0)
  expect_identical(fit_fixed_pars(saved), list(sigma2_me = 0))
  expect_identical(coef(saved), coef(gaussian_offset_model))
  held <- saved
  held$fix_tau2 <- 0.3
  held$estimate$nu2 <- NULL
  expect_identical(fit_fixed_pars(held), list(tau2 = 0.3, sigma2_me = 0))
  expect_identical(coef(held)$tau2, 0.3)
  size <- list(family = "negative_binomial", fix_tau2 = 0,
               family_parameters = list(theta = 2))
  expect_identical(fit_family_held(size), list(theta = 2))
  expect_identical(fit_fixed_pars(size), list(theta = 2))
})

test_that("a profile of an MCML fit holds the nugget variance", {
  set.seed(2)
  data <- binomial_data[rep(seq_len(nrow(binomial_data)), each = 3), ]
  data$y <- rbinom(nrow(data), data$denominator, 0.4)
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(nugget = TRUE), data = data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    return_samples = FALSE, messages = FALSE
  ))
  problem <- make_profile_problem(fit, 500L, 1L)
  expect_true("tau2" %in% problem$names)
  expect_false("nu2" %in% problem$names)
  nugget <- match("tau2", problem$names)
  expect_equal(problem$estimate[["tau2"]],
               fit$estimate$nu2 + fit$estimate$sigma2, tolerance = 1e-12)
  expect_identical(problem$labels$tau2, "Variance of the nugget")
  # At a held value of log tau2, the spatial variance moves and tau2 stays.
  target <- problem$estimate[[nugget]] + 0.5
  point <- maximise_profile(problem, nugget, target, problem$bank_estimate)
  expect_equal(point$par[[nugget]], target, tolerance = 1e-12)
  expect_identical(list_covariance_parameters(fit), c("sigma2", "phi", "tau2"))
})

test_that("a fit saved with gp(nugget = <value>) still predicts and refits", {
  # A fit of the earlier interface kept its formula as the user gave it and
  # held the nugget in fix_tau2.
  current <- glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data,
                   family = "gaussian", fixed_pars = list(tau2 = 0.2),
                   messages = FALSE)
  saved <- current
  saved$formula <- y ~ cov + gp(nugget = 0.2)
  saved$nugget <- NULL
  saved$fixed_pars <- NULL
  saved$fix_tau2 <- 0.2
  saved["fix_var_me"] <- list(NULL)
  expect_identical(fit_formula(saved), y ~ cov + gp(nugget = TRUE))
  expect_identical(fit_formula(current), current$formula)
  expect_identical(fit_fixed_pars(saved), list(tau2 = 0.2))
  expect_identical(coef(saved), coef(current))
  grid <- gaussian_data[1:3, ]
  expect_identical(
    predict(saved, grid, seed = 1L, messages = FALSE)$predictions,
    predict(current, grid, seed = 1L, messages = FALSE)$predictions
  )
  expect_identical(simulation_model(saved)$parameters,
                   simulation_model(current)$parameters)
})

test_that("gp() points a value of the nugget to fixed_pars or to FALSE", {
  expect_match(conditionMessage(expect_input_error(gp(nugget = 0.4), "nugget")),
               "fixed_pars = list(tau2 = 0.4)", fixed = TRUE)
  zero <- conditionMessage(expect_input_error(gp(nugget = 0), "nugget"))
  expect_match(zero, "nugget = FALSE", fixed = TRUE)
  expect_no_match(zero, "fixed_pars")
  for (value in list(-1, NA_real_, Inf)) {
    expect_no_match(conditionMessage(expect_input_error(gp(nugget = value),
                                                        "nugget")),
                    "fixed_pars")
  }
})

replicated_gaussian <- function() {
  set.seed(3)
  data <- do.call(rbind, replicate(4, gaussian_data, simplify = FALSE))
  data$y <- data$y + rnorm(nrow(data), sd = 0.3)
  data
}

test_that("a held slope is the model with its term as an offset", {
  data <- replicated_gaussian()
  held <- glgpm(y ~ cov + gp(), data = data, family = "gaussian",
                fixed_pars = list(beta = c(cov = 0.5)), messages = FALSE)
  data$known <- 0.5 * data$cov
  offset <- glgpm(y ~ 1 + offset(known) + gp(), data = data,
                  family = "gaussian", messages = FALSE)
  expect_identical(coef(held)$beta[["cov"]], 0.5)
  expect_equal(held$log_lik, offset$log_lik, tolerance = 1e-8)
  expect_equal(unlist(coef(held)[c("sigma2", "phi", "sigma2_me")]),
               unlist(coef(offset)[c("sigma2", "phi", "sigma2_me")]),
               tolerance = 1e-5)
  expect_equal(coef(held)$beta[["(Intercept)"]],
               coef(offset)$beta[["(Intercept)"]], tolerance = 1e-5)
  # The held coefficient has no standard error, and the others keep theirs.
  se <- sqrt(diag(held$covariance))
  expect_true(is.na(se[["beta.cov"]]))
  expect_true(all(is.finite(se[names(se) != "beta.cov"])))
  expect_identical(held$fixed_pars, list(beta = c(cov = 0.5)))
  expect_false("beta.cov" %in% names(estimated_coef(held)))
  expect_match(capture.output(print(summary(held))), "cov fixed at 0.5",
               fixed = TRUE, all = FALSE)
})

test_that("a held intercept is held on an uncentred design", {
  data <- replicated_gaussian()
  fit <- glgpm(y ~ cov + gp(), data = data, family = "gaussian",
               fixed_pars = list(beta = c("(Intercept)" = 1)),
               messages = FALSE)
  expect_identical(coef(fit)$beta[["(Intercept)"]], 1)
  expect_equal(fit$estimate$beta[["(Intercept)"]], 1, tolerance = 1e-12)
  scaling <- model_design_scaling(fit$D, fit_parameter_table(fit))
  expect_identical(unname(scaling$center), c(0, 0))
  expect_true(is.finite(sqrt(fit$covariance["beta.cov", "beta.cov"])))
})

test_that("a held scale is the profile of the scale at that value", {
  data <- replicated_gaussian()
  fit <- glgpm(y ~ cov + gp(), data = data, family = "gaussian",
               messages = FALSE)
  held <- glgpm(y ~ cov + gp(), data = data, family = "gaussian",
                fixed_pars = list(phi = 2), messages = FALSE)
  expect_identical(coef(held)$phi, 2)
  expect_true(is.na(held$covariance["phi", "phi"]))
  problem <- make_gaussian_profile_problem(fit)
  index <- match("phi", problem$names)
  point <- maximise_profile(problem, index, log(2), problem$bank_estimate)
  expect_equal(point$value, held$log_lik, tolerance = 1e-6)
  # The held parameter has no profile, and the profiles of the others hold
  # it.
  expect_input_error(profile(held, parameters = "phi", messages = FALSE),
                     "parameters")
  expect_false("phi" %in% list_covariance_parameters(held))
  profiles <- profile(held, parameters = "sigma2", n_points = 5,
                      messages = FALSE)
  expect_true(all(profiles$marginal$sigma2$phi == held$estimate$phi))
  expect_equal(held$estimate$phi, log(2), tolerance = 1e-12)
})

test_that("held variances of a binomial fit stay at their values", {
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    return_samples = FALSE, messages = FALSE,
    fixed_pars = list(sigma2_re = c(i = 0.3), sigma2 = 0.5)
  ))
  expect_identical(coef(fit)$sigma2, 0.5)
  expect_identical(coef(fit)$sigma2_re, c(i = 0.3))
  expect_true(all(is.na(diag(fit$covariance)[c("sigma2", "sigma2_re.i")])))
  expect_true(all(is.finite(diag(fit$covariance)[c("beta.(Intercept)",
                                                    "beta.cov", "phi")])))
  expect_true(is.finite(summary(fit)$mcml$certificate))
  # The bootstrap refits with the holds, with columns for the others only.
  boot <- suppressWarnings(parametric_bootstrap(fit, n_sim = 2, seed = 3L,
                                                messages = FALSE))
  expect_identical(colnames(boot$estimates),
                   c("beta.(Intercept)", "beta.cov", "phi"))
  # start_pars = coef(fit) works for a refit with the same holds.
  refit <- suppressWarnings(glgpm(
    y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    return_samples = FALSE, messages = FALSE, start_pars = coef(fit),
    fixed_pars = list(sigma2_re = c(i = 0.3), sigma2 = 0.5)
  ))
  expect_identical(attr(refit, "starting_values")$selected, "user")
  expect_input_error(
    glgpm(y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
          denominator = denominator, messages = FALSE,
          start_pars = list(sigma2 = 1), fixed_pars = list(sigma2 = 0.5)),
    "start_pars"
  )
  expect_input_error(
    glgpm(y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
          denominator = denominator, messages = FALSE,
          fixed_pars = list(sigma2_re = c(i = 0))),
    "fixed_pars"
  )
})

test_that("a model with every parameter held is evaluated at the values", {
  values <- list(beta = c(1, 0.3), sigma2 = 0.5, phi = 1.2, sigma2_me = 0.1)
  fit <- glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
               fixed_pars = values, messages = FALSE)
  expect_identical(coef(fit)[names(values)],
                   c(list(beta = c("(Intercept)" = 1, cov = 0.3)),
                     values[-1]))
  expect_true(all(is.na(fit$covariance)))
  # The exact log-likelihood at the values, as the profile computes it with
  # its constant.
  free <- glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
                messages = FALSE)
  problem <- make_gaussian_profile_problem(free)
  constant <- free$log_lik - problem$value(problem$estimate)
  working <- c(1, 0.3, log(0.5), log(1.2), log(0.1))
  expect_equal(fit$log_lik, problem$value(working) + constant,
               tolerance = 1e-8)
  binomial <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    messages = FALSE,
    fixed_pars = list(beta = c(0.6, 0), sigma2 = 0.5, phi = 0.3)
  ))
  expect_true(is.finite(binomial$log_lik))
  expect_null(summary(binomial)$mcml)
  expect_false(is.null(binomial$S_samples))
  prediction <- predict(binomial, binomial_data[1:3, ], seed = 1L,
                        messages = FALSE)
  expect_true(all(is.finite(prediction$predictions$response_mean)))
})

test_that("predictions with parameter uncertainty keep held parameters", {
  data <- replicated_gaussian()
  held <- glgpm(y ~ cov + gp(), data = data, family = "gaussian",
                fixed_pars = list(phi = 2), messages = FALSE)
  set.seed(5)
  draws <- draw_parameter_values(held, 50L)
  expect_true(all(draws[, "phi"] == held$estimate$phi))
  expect_gt(stats::sd(draws[, "sigma2"]), 0)
  all_held <- glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
                    fixed_pars = list(beta = c(1, 0.3), sigma2 = 0.5,
                                      phi = 1.2, sigma2_me = 0.1),
                    messages = FALSE)
  expect_input_error(draw_parameter_values(all_held, 5L),
                     "parameter_uncertainty")
})

test_that("check_fit() says that a fit holding every parameter estimates nothing", {
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    messages = FALSE,
    fixed_pars = list(beta = c(0.6, 0), sigma2 = 0.5, phi = 0.3)
  ))
  checks <- check_fit(fit, messages = FALSE)$checks
  row <- checks[checks$check == "Monte Carlo accuracy", ]
  expect_identical(row$status, "OK")
  expect_match(row$finding, "Every parameter is held by `fixed_pars`",
               fixed = TRUE)
  # The log-likelihood is still a Monte Carlo estimate.
  expect_match(row$finding, "Monte Carlo estimate from one bank of", fixed = TRUE)
  gaussian <- check_fit(gaussian_model, messages = FALSE)$checks
  expect_identical(gaussian$finding[gaussian$check == "Monte Carlo accuracy"],
                   "A Gaussian fit has no Monte Carlo error.")
})

test_that("held coefficients are fixed directions in the search for separation", {
  x <- rep(c(0, 1), each = 20)
  design <- cbind("(Intercept)" = 1, x = x)
  y <- ifelse(x == 1, 0, rep(c(3, 5, 4, 6), 10))
  trials <- rep(10, 40)
  expect_identical(find_separation(design, y, trials, "binomial")$coefficients,
                   2L)
  expect_null(find_separation(design, y, trials, "binomial", held = 2L))
  # With no positives, a held intercept and a covariate of both signs leave
  # a finite maximum.
  centred <- cbind("(Intercept)" = 1, z = seq(-1, 1, length.out = 40))
  expect_identical(find_separation(centred, 0 * y, trials,
                                   "binomial")$everywhere, "no positives")
  expect_null(find_separation(centred, 0 * y, trials, "binomial", held = 1L))
})

test_that("the screen of short scales doesn't run with held coefficients", {
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    return_samples = FALSE, messages = FALSE,
    fixed_pars = list(beta = c(cov = 0))
  ))
  record <- attr(fit, "optimizer")$mcml_maxima_screen
  expect_identical(record$gate, "unavailable")
  expect_identical(record$reason, "held coefficients")
  expect_null(record$alarm)
  expect_true(is.finite(summary(fit)$mcml$certificate))
})

test_that("glgpm_parameters() shows the parameters that fixed_pars holds", {
  parameters <- glgpm_parameters(y ~ cov + gp(), data = binomial_data,
                                 family = "binomial",
                                 fixed_pars = list(phi = 20,
                                                   beta = c(cov = 0)))
  expect_identical(parameters$held[parameters$name == "phi"], 20)
  expect_identical(parameters$held[parameters$element %in% "cov"], 0)
  expect_output(print(parameters),
                "list(beta = c(\"(Intercept)\" = ), sigma2 = )", fixed = TRUE)
})

test_that("a held spatial variance leaves the nugget its standard error", {
  data <- replicated_gaussian()
  fit <- glgpm(y ~ cov + gp(nugget = TRUE), data = data, family = "gaussian",
               fixed_pars = list(sigma2 = 0.5), messages = FALSE)
  # log tau2 = nu2 + log sigma2, and sigma2 is known.
  limits <- summary(fit)$sp["Variance of the nugget", ]
  expect_true(all(is.finite(limits)))
  profiled <- profile_nugget(fit, fit_parameter_table(fit))
  expect_equal(profiled$standard_errors[["tau2"]],
               sqrt(fit$covariance["nu2", "nu2"]), tolerance = 1e-12)
})
