# start_pars is a partial named list over the model's estimated parameters
# (step 2c of the parameter table, P1 and P2 of the interface review).

two_terms <- function() {
  parameter_table(c("(Intercept)", "x"), fix_tau2 = TRUE,
                  re_names = c("village", "household"), distance_units = "km")
}

test_that("starts are read by name, flat name and position", {
  parameters <- two_terms()
  expect_identical(read_start_pars(parameters, NULL), list())
  given <- read_start_pars(parameters, list(
    phi = 3, beta = c(x = 0.5), sigma2_re = c(household = 0.2)
  ))
  expect_identical(given, list(phi = 3, beta = c(x = 0.5),
                               sigma2_re = c(household = 0.2)))
  # The flat names of unlist(coef(fit)), and NULL elements left out.
  given <- read_start_pars(parameters, list(
    "beta.(Intercept)" = -1, "sigma2_re.village" = 0.3, tau2 = 0.1,
    sigma2 = NULL
  ))
  expect_identical(given, list(beta = c("(Intercept)" = -1),
                               sigma2_re = c(village = 0.3), tau2 = 0.1))
  # Unnamed vectors at full length, read by position.
  given <- read_start_pars(parameters, list(beta = c(1, 2),
                                            sigma2_re = c(0.3, 0.4)))
  expect_identical(given$beta, c("(Intercept)" = 1, x = 2))
  expect_identical(given$sigma2_re, c(village = 0.3, household = 0.4))
  # The coefficients in the order of the design, whatever their order given.
  given <- read_start_pars(parameters, list(beta = c(x = 2, "(Intercept)" = 1)))
  expect_identical(given$beta, c("(Intercept)" = 1, x = 2))
  one_term <- parameter_table("(Intercept)", fix_tau2 = 0, re_names = "g")
  expect_identical(read_start_pars(one_term, list(sigma2_re = 2))$sigma2_re,
                   c(g = 2))
})

test_that("starts that don't fit the model stop", {
  parameters <- two_terms()
  bad <- list(
    list(sigma2_me = 1), list(theta = 2), list(beta = c(z = 1)),
    list(beta = 1), list(sigma2_re = 0.2), list(sigma2_re = c(school = 1)),
    list(beta = c(x = NA)), list(sigma2 = -1), list(phi = Inf),
    list(sigma2 = c(1, 2)), list(beta = c(x = 1), beta.x = 2),
    list(1), "sigma2", list(sigma2 = "1")
  )
  for (start in bad) {
    expect_error(read_start_pars(parameters, start),
                 class = "riskmap_error_input")
  }
  error <- tryCatch(read_start_pars(parameters, list(sigma2_me = 1)),
                    error = identity)
  expect_match(conditionMessage(error), "glgpm_parameters", fixed = TRUE)
  expect_identical(error$arg, "start_pars")
  # A held parameter takes no start.
  held <- model_parameter_table("(Intercept)", fix_tau2 = 0.5, exact = TRUE,
                                fix_var_me = 0.1, estimated_family = character(),
                                family_parameters = NULL, re_names = NULL,
                                distance_units = "km")
  error <- tryCatch(read_start_pars(held, list(tau2 = 1)), error = identity)
  expect_match(conditionMessage(error), "holds at 0.5", fixed = TRUE)
  expect_error(read_start_pars(held, list(sigma2_me = 1)),
               class = "riskmap_error_input")
})

test_that("a start of phi far from the distances warns", {
  parameters <- two_terms()
  coords <- cbind(c(0, 1, 4), c(0, 0, 3))
  expect_no_warning(read_start_pars(parameters, list(phi = 2), coords, "km"))
  expect_warning(read_start_pars(parameters, list(phi = 0.001), coords, "km"),
                 "km", class = "riskmap_warning")
  expect_warning(read_start_pars(parameters, list(phi = 500), coords, "km"),
                 class = "riskmap_warning")
})

test_that("given starts replace their elements of a complete start", {
  parameters <- two_terms()
  start <- list(beta = c(0, 0), sigma2 = 1, phi = 2, tau2 = 1,
                sigma2_re = c(1, 1))
  given <- list(beta = c(x = 0.5), sigma2_re = c(household = 0.2), phi = 3)
  expect_identical(apply_given_start(parameters, start, given),
                   list(beta = c(0, 0.5), sigma2 = 1, phi = 3, tau2 = 1,
                        sigma2_re = c(1, 0.2)))
  expect_false(all_starts_given(parameters, given))
  expect_true(all_starts_given(parameters, list(
    beta = c("(Intercept)" = 1, x = 2), sigma2 = 1, phi = 2, tau2 = 1,
    sigma2_re = c(village = 1, household = 1)
  )))
  expect_identical(given_start_names(given),
                   c("beta.x", "sigma2_re.household", "phi"))
  expect_null(given_start_names(list()))
})

fit_binomial_with <- function(start_pars = NULL,
                              control_mcml = fixture_control_mcml) {
  suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = control_mcml,
    return_samples = FALSE, messages = FALSE, start_pars = start_pars
  ))
}

test_that("a partial start keeps the comparison of starts", {
  skip_on_cran()
  default <- fit_binomial_with()
  partial <- fit_binomial_with(list(phi = 3))
  starting <- attr(partial, "starting_values")
  expect_identical(starting$given, "phi")
  expect_true(starting$selected %in% c("current", "transformed_gaussian"))
  expect_true(is.finite(starting$current_laplace))
  expect_null(attr(default, "starting_values")$given)
  expect_match(describe_starting_values(starting$selected, starting$reason,
                                        starting$given),
               "^The starting values of `phi` were given in `start_pars`\\. The other starting values")
  expect_identical(summary(partial)$mcml$starting_given, "phi")
  expect_null(summary(default)$mcml$starting_given)
  # set_control_mcml(start_comparison = FALSE) skips it.
  skipped <- fit_binomial_with(
    list(phi = 3),
    control_mcml = set_control_mcml(initial_draws = 500,
                                    validation_draws = 500,
                                    max_draws = 500, max_stages = 1,
                                    start_comparison = FALSE)
  )
  starting <- attr(skipped, "starting_values")
  expect_identical(starting[c("selected", "reason")],
                   list(selected = "current", reason = "not_compared"))
  expect_true(is.na(starting$current_laplace))
  expect_match(describe_starting_values(starting$selected, starting$reason,
                                        starting$given),
               "start_comparison = FALSE", fixed = TRUE)
  # A complete start needs no comparison.
  complete <- fit_binomial_with(coef(default))
  expect_identical(attr(complete, "starting_values")$selected, "user")
  expect_identical(describe_starting_values("user", NULL, "phi"),
                   "The starting values were all given in `start_pars`.")
  expect_error(set_control_mcml(start_comparison = NA),
               class = "riskmap_error_input")
})

test_that("start_pars = coef(fit) starts the engine at the fit's estimate", {
  skip_on_cran()
  # The engine's start, captured, against the fit's estimate on the scale of
  # the engine, with standardised covariates, distances rescaled by the
  # largest distance and, for a Gaussian model, the standardised response.
  engine_start <- function(fit, refit) {
    captured <- NULL
    real <- working_start
    try(with_mocked_bindings(
      refit(),
      working_start = function(table, start) {
        captured <<- real(table, start)
        stop("captured the start")
      }
    ), silent = TRUE)
    maps <- profile_scale_maps(
      standardize_design_matrix(fit$D)$coefficient_transform,
      log(scale_spatial_coordinates(fit$coords)$distance_scale),
      ncol(fit$D)
    )
    expected <- unname(maps$to(unlist(fit$estimate)))
    if (fit$family == "gaussian") {
      # The Gaussian engines fit the standardised residuals of the regression.
      regression <- lm.fit(fit$D, fit$y)
      expected <- unlist(fit$estimate)
      expected[seq_len(ncol(fit$D))] <- unlist(fit$estimate$beta) -
        regression$coefficients
      expected <- unname(maps$to(expected))
      log_scale <- log(mean(regression$residuals^2)) / 2
      coefficients <- seq_len(ncol(fit$D))
      expected[coefficients] <- expected[coefficients] / exp(log_scale)
      variances <- grepl("^sigma2", names(unlist(fit$estimate)))
      expected[variances] <- expected[variances] - 2 * log_scale
    }
    list(captured = captured, expected = expected)
  }
  count_fit <- function(family, data, start_pars = NULL) {
    suppressWarnings(glgpm(
      y ~ cov + gp(), data = data, family = family, denominator = denominator,
      control_mcmc = control_mcmc, seed = fixture_seed,
      control_mcml = fixture_control_mcml, return_samples = FALSE,
      messages = FALSE, start_pars = start_pars
    ))
  }
  set.seed(2)
  replicated <- do.call(rbind, replicate(6, gaussian_data, simplify = FALSE))
  replicated$y <- replicated$y + rnorm(nrow(replicated), sd = 0.3)
  gaussian_fit <- function(data, start_pars = NULL) {
    glgpm(y ~ cov + gp(nugget = TRUE) + re(i), data = data,
          family = "gaussian", messages = FALSE, start_pars = start_pars)
  }
  fits <- list(
    binomial = list(fit = binomial_model, refit = function() {
      suppressWarnings(glgpm(
        y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
        denominator = denominator, control_mcmc = control_mcmc,
        seed = fixture_seed, control_mcml = fixture_control_mcml,
        messages = FALSE, start_pars = coef(binomial_model)
      ))
    }),
    poisson = local({
      fit <- count_fit("poisson", poisson_data)
      list(fit = fit, refit = function() {
        count_fit("poisson", poisson_data, coef(fit))
      })
    }),
    negative_binomial = local({
      fit <- count_fit(negative_binomial(), poisson_data)
      list(fit = fit, refit = function() {
        count_fit(negative_binomial(), poisson_data, coef(fit))
      })
    }),
    gaussian_direct = list(fit = gaussian_model, refit = function() {
      glgpm(y ~ cov + gp() + re(i), data = gaussian_data, family = "gaussian",
            messages = FALSE, start_pars = coef(gaussian_model))
    }),
    gaussian_woodbury = local({
      fit <- gaussian_fit(replicated)
      list(fit = fit, refit = function() gaussian_fit(replicated, coef(fit)))
    })
  )
  for (name in names(fits)) {
    starts <- engine_start(fits[[name]]$fit, fits[[name]]$refit)
    expect_length(starts$captured, length(starts$expected))
    expect_equal(starts$captured, starts$expected, tolerance = 1e-10,
                 label = name)
  }
})

test_that("glgpm_parameters() lists the parameters before fitting", {
  listed <- glgpm_parameters(y ~ cov + gp(nugget = TRUE) + re(i, j),
                             data = binomial_data, family = "binomial")
  expect_s3_class(listed, "riskmap_parameters")
  expect_identical(listed$flat, c("beta.(Intercept)", "beta.cov", "sigma2",
                                  "phi", "tau2", "sigma2_re.i",
                                  "sigma2_re.j"))
  expect_identical(listed$start[c(1, 4)], c("regression", "distance_quantile"))
  expect_true(all(is.na(listed$held)))
  # A model without a nugget lists none.
  expect_false("tau2" %in% glgpm_parameters(y ~ gp(), data = binomial_data,
                                            family = "binomial")$name)
  # The names are those of coef() of the fit.
  expect_identical(
    glgpm_parameters(y ~ cov + gp() + re(i), data = binomial_data,
                     family = "binomial")$flat,
    names(unlist(coef(binomial_model)))
  )
  # Held parameters have a value.
  held <- glgpm_parameters(y ~ cov + gp(nugget = TRUE), data = gaussian_data,
                           family = "gaussian",
                           fixed_pars = list(tau2 = 0.2, sigma2_me = 0.1))
  expect_identical(held$held[held$name %in% c("tau2", "sigma2_me")],
                   c(0.2, 0.1))
  counts <- glgpm_parameters(y ~ gp(), data = poisson_data,
                             family = negative_binomial())
  expect_identical(counts$start[counts$name == "theta"], "family_regression")
  output <- capture.output(print(listed))
  expect_match(output[length(output)],
               "list(beta = c(\"(Intercept)\" = , \"cov\" = ), sigma2 = , phi = , tau2 = , sigma2_re = c(\"i\" = , \"j\" = ))",
               fixed = TRUE)
  expect_error(glgpm_parameters(y ~ gp(), data = binomial_data,
                                family = "binomial",
                                fixed_pars = list(sigma2_me = 0.1)),
               class = "riskmap_error_input")
  # A size held by fixed_pars.
  held_size <- glgpm_parameters(y ~ gp(), data = poisson_data,
                                family = negative_binomial(),
                                fixed_pars = list(theta = 2))
  expect_identical(held_size$held[held_size$name == "theta"], 2)
})

test_that("a partial start reaches the engine whichever candidate wins", {
  # The transformed candidate is made to win, so the value given must have
  # replaced the Gaussian fit's own on the way to the engine.
  scores <- 0
  captured <- NULL
  real <- working_start
  try(with_mocked_bindings(
    fit_binomial_with(list(phi = 3)),
    starting_laplace_log_marginal = function(...) {
      scores <<- scores + 1
      -100 + scores
    },
    working_start = function(table, start) {
      # The transformed Gaussian model, which estimates a measurement error
      # variance, starts as usual.
      if ("sigma2_me" %in% table$rows$name) return(real(table, start))
      captured <<- start
      stop("captured the start")
    }
  ), silent = TRUE)
  expect_identical(scores, 2)
  expect_identical(captured$phi,
                   3 / scale_spatial_coordinates(binomial_model$coords)$distance_scale)
})

test_that("a start that leaves nothing to compare skips the comparison", {
  poisson_fit <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = poisson_data, family = "poisson",
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    return_samples = FALSE, messages = FALSE
  ))
  fitted <- 0
  counts <- with_mocked_bindings(
    suppressWarnings(glgpm(
      y ~ cov + gp(), data = poisson_data, family = negative_binomial(),
      denominator = denominator, control_mcmc = control_mcmc,
      seed = fixture_seed, control_mcml = fixture_control_mcml,
      return_samples = FALSE, messages = FALSE,
      start_pars = coef(poisson_fit)
    )),
    glgpm_lm = function(...) {
      fitted <<- fitted + 1
      stop("the transformed Gaussian model was fitted")
    }
  )
  expect_identical(fitted, 0)
  starting <- attr(counts, "starting_values")
  expect_identical(starting[c("selected", "reason")],
                   list(selected = "current", reason = "candidates_given"))
  expect_match(describe_starting_values(starting$selected, starting$reason,
                                        starting$given),
               "could change.", fixed = TRUE)
})

test_that("the named vector unlist(coef(fit)) is a start", {
  default <- fit_binomial_with()
  expect_identical(
    read_start_pars(fit_parameter_table(default), unlist(coef(default))),
    read_start_pars(fit_parameter_table(default), coef(default))
  )
  expect_error(glgpm_parameters(y ~ gp(), data = gaussian_data,
                                family = "gaussian",
                                fixed_pars = list(sigma2_me = -1)),
               class = "riskmap_error_input")
})

test_that("print() of glgpm_parameters() stops on an argument it doesn't use (N9)", {
  listed <- glgpm_parameters(y ~ gp(), data = gaussian_data,
                             family = "gaussian")
  expect_input_error(print(listed, digits = 3), "digits")
  expect_no_error(capture.output(print(listed)))
})
