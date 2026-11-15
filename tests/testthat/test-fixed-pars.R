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
