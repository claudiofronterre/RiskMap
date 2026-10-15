test_that("check_dast produces errors", {

  expect_error(
    check_dast(y ~ gp(a, b), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = NULL, start_pars = list()),
    "Specify gp() for a purely spatial model, or gp(t) for a spatio-temporal model",
    fixed = TRUE
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = NULL,
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = NULL, start_pars = list()),
    "You must supply a time column via the 'time' argument"
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(not_present),
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = NULL, start_pars = list()),
    "'time' column not found in 'data'"
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = NULL,
               penalty = NULL, start_pars = list()),
    "'int_mat' must be provided"
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = dast_int_mat[-1, ],
               penalty = NULL, start_pars = list()),
    sprintf("'int_mat' must have %d rows", nrow(dast_data))
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = c(dast_mda_times, 5), int_mat = dast_int_mat,
               penalty = NULL, start_pars = list()),
    "'int_mat' must have 3 columns"
  )

  no_penalty <- function(x) 0

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = no_penalty, start_pars = list()),
    "'penalty' must be a list of 6 functions"
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = rep(list(no_penalty), 5), start_pars = list()),
    "'penalty' must be a list of 6 functions"
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = c(rep(list(no_penalty), 5), 0), start_pars = list()),
    "'penalty' must be a list of 6 functions"
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = NULL, start_pars = list(psi = -1)),
    "The starting value for 'psi' must be a single positive number"
  )

  expect_error(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = NULL, start_pars = list(psi = "a")),
    "The starting value for 'psi' must be a single positive number"
  )

  expect_true(
    check_dast(y ~ gp(), data = dast_data, time = quote(survey_time),
               mda_times = dast_mda_times, int_mat = dast_int_mat,
               penalty = rep(list(no_penalty), 6),
               start_pars = list(psi = 1))
  )
})

expected_output <- c("estimate", "grad_MLE", "covariance", "log_lik",
                     "S_samples", "y", "D", "coords", "sst", "mda_times",
                     "survey_times_data", "int_mat", "ID_coords", "power_val",
                     "fix_tau2", "formula", "crs", "input_crs",
                     "distance_units", "data_sf", "family", "kappa",
                     "units_m", "cov_offset", "call", "linkf")

test_that("dast produces expected output for spatial models", {
  fit <- dast(y ~ gp(), data = dast_data, denominator = denominator, time = survey_time,
             mda_times = dast_mda_times, int_mat = dast_int_mat, power_val = 1,
             control_mcmc = control_mcmc, messages = FALSE)

  expect_s3_class(fit, c("RiskMap_dast", "RiskMap"), exact = TRUE)
  expect_setequal(names(fit), expected_output)
  expect_identical(fit$family, "binomial")
  expect_false(fit$sst)
  expect_identical(fit$y, as.double(dast_data$y))
  expect_identical(fit$units_m, as.double(dast_data$denominator))
  expect_identical(fit$survey_times_data, dast_data$survey_time)
  expect_identical(fit$mda_times, dast_mda_times)
  expect_identical(fit$int_mat, dast_int_mat)
  expect_identical(fit$power_val, 1)
  expect_named(fit$estimate, c("beta", "sigma2", "phi", "alpha", "gamma"))
  expect_true(all(is.finite(unlist(fit$estimate))))
  expect_identical(dim(fit$covariance), rep(length(unlist(fit$estimate)), 2))
  expect_identical(ncol(fit$S_samples), nrow(fit$coords))

  # coef() back-transforms the decay parameters
  estimates <- coef(fit)
  expect_identical(estimates$alpha, plogis(fit$estimate$alpha))
  expect_identical(estimates$gamma, exp(fit$estimate$gamma))
  expect_gt(estimates$alpha, 0)
  expect_lt(estimates$alpha, 1)
  expect_null(estimates$psi)

  fit_summary <- summary(fit)
  expect_s3_class(fit_summary, c("summary.RiskMap_dast", "summary.RiskMap"),
                  exact = TRUE)
  expect_identical(rownames(fit_summary$mda),
               c("Drop (alpha)", "Scale of the decay (gamma)"))
  expect_identical(colnames(fit_summary$mda),
               c("Estimate", "Lower limit", "Upper limit"))
  expect_identical(unname(fit_summary$mda[, "Estimate"]),
               c(estimates$alpha, estimates$gamma))
  expect_true(all(fit_summary$mda[, "Lower limit"] <= fit_summary$mda[, "Estimate"]))
  expect_true(all(fit_summary$mda[, "Upper limit"] >= fit_summary$mda[, "Estimate"]))
  expect_output(print(fit_summary), "MDA impact function")
})

test_that("dast produces expected output for spatio-temporal models with a fixed drop", {
  st_data <- dast_data
  st_data$t_gp <- rep(1:2, length.out = nrow(st_data))

  fit <- dast(y ~ gp(t_gp), data = st_data, denominator = denominator, time = survey_time,
             mda_times = dast_mda_times, int_mat = dast_int_mat, power_val = 1,
             drop = 0.5, control_mcmc = control_mcmc, messages = FALSE)

  expect_s3_class(fit, "RiskMap_dast")
  expect_setequal(names(fit), c(expected_output, "time", "fix_alpha"))
  expect_true(fit$sst)
  expect_identical(fit$fix_alpha, 0.5)
  expect_named(fit$estimate, c("beta", "sigma2", "phi", "gamma", "psi"))
  expect_true(all(is.finite(unlist(fit$estimate))))

  estimates <- coef(fit)
  expect_null(estimates$alpha)
  expect_identical(estimates$psi, exp(fit$estimate$psi))

  fit_summary <- summary(fit)
  expect_identical(rownames(fit_summary$mda), "Scale of the decay (gamma)")
  expect_true("Temporal corr. scale (psi)" %in% rownames(fit_summary$sp))
  expect_output(print(fit_summary), "Drop (alpha) fixed at 0.5", fixed = TRUE)
})
