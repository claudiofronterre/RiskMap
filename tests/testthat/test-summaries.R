# Tests for the summary() methods and their print() methods.
#
# The print methods only format summary objects, so their tests build the
# summaries by hand with fixed values. Snapshots of summaries from real fits
# would change with platform, BLAS or optimiser details rather than with the
# formatting being tested.

# summary.RiskMap ------------------------------------------------------------

test_that("estimates are consistent between coef and summary", {

  sum <- summary(gaussian_model)
  cof <- coef(gaussian_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_me, sum$family_parameters[1,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)

  sum <- summary(binomial_model)
  cof <- coef(binomial_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)

  sum <- summary(poisson_model)
  cof <- coef(poisson_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)
})

test_that("summary.RiskMap returns the expected components for each family", {
  gaussian_summary <- summary(gaussian_model)
  expect_s3_class(gaussian_summary, "summary.RiskMap")
  expect_setequal(names(gaussian_summary),
                  c("reg_coef", "family_parameters", "held_family_parameters",
                    "sp", "tau2", "ranef", "level",
                    "family", "kappa", "log_lik", "cov_offset_used", "aic",
                    "call", "link_name", "invlink_expression"))
  expect_identical(gaussian_summary$family, "gaussian")
  expect_identical(gaussian_summary$link_name, "identity")
  expect_identical(gaussian_summary$call, gaussian_model$call)
  expect_identical(gaussian_summary$log_lik, gaussian_model$log_lik)
  expect_identical(gaussian_summary$kappa, gaussian_model$kappa)
  expect_false(gaussian_summary$cov_offset_used)
  expect_equal(gaussian_summary$aic,
               2 * length(unlist(gaussian_model$estimate)) - 2 * gaussian_model$log_lik)
  expect_identical(rownames(gaussian_summary$ranef), "i (random eff. var.)")
  expect_identical(rownames(gaussian_summary$sp),
                   c("Spatial process var.", "Spatial corr. scale (km)"))

  binomial_summary <- summary(binomial_model)
  expect_identical(binomial_summary$link_name, "canonical (logit)")
  expect_null(binomial_summary$family_parameters)
  expect_identical(binomial_summary$held_family_parameters, list())
  expect_null(binomial_summary$aic)

  poisson_summary <- summary(poisson_model)
  expect_identical(poisson_summary$link_name, "canonical (log)")
  expect_null(poisson_summary$family_parameters)
  expect_null(poisson_summary$aic)
})

test_that("summary.RiskMap reports a custom inverse link", {
  custom_link <- list(
    inv = function(x) 1 / (1 + exp(-x)),
    d1 = function(x) plogis(x) * (1 - plogis(x)),
    d2 = function(x) {
      probability <- plogis(x)
      probability * (1 - probability) * (1 - 2 * probability)
    }
  )
  fit <- glgpm(y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
               denominator = denominator, control_mcmc = control_mcmc,
               seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
               invlink = custom_link, messages = FALSE)

  fit_summary <- summary(fit)
  expect_identical(fit_summary$link_name, "custom")
  expect_identical(fit_summary$invlink_expression,
                   "Inverse link function = 1/(1 + exp(-x))")
})

test_that("summary.RiskMap names a package or primitive inverse link", {
  # The fit is used only for its structure, so its link is replaced.
  fit <- binomial_model
  fit$link_function <- list(inv = pnorm, name = "custom")
  fit$call$invlink <- quote(pnorm)
  expect_identical(summary(fit)$invlink_expression,
                   "Inverse link function = pnorm(x)")
  fit$call$invlink <- quote(stats::pnorm)
  expect_identical(summary(fit)$invlink_expression,
                   "Inverse link function = stats::pnorm(x)")
  fit$link_function$inv <- exp
  fit$call$invlink <- quote(exp)
  expect_identical(summary(fit)$invlink_expression,
                   "Inverse link function = exp(x)")
  # A refit's call holds the function itself, not its name.
  fit$call$invlink <- exp
  expect_identical(summary(fit)$invlink_expression,
                   "Inverse link function = <user-supplied function>")
})

test_that("summary.RiskMap computes standard errors and confidence intervals", {
  fit_summary <- summary(binomial_model)
  estimate <- unlist(binomial_model$estimate)
  standard_error <- sqrt(diag(binomial_model$covariance))
  is_beta <- startsWith(names(estimate), "beta.")
  z_crit <- qnorm(0.975)

  reg_coef <- fit_summary$reg_coef
  expect_identical(colnames(reg_coef),
                   c("Estimate", "Lower limit", "Upper limit", "StdErr",
                     "z.value", "p.value"))
  expect_identical(rownames(reg_coef), colnames(binomial_model$D))
  expect_equal(reg_coef[, "Estimate"], estimate[is_beta], ignore_attr = TRUE)
  expect_equal(reg_coef[, "StdErr"], standard_error[is_beta], ignore_attr = TRUE)
  expect_equal(reg_coef[, "Lower limit"],
               estimate[is_beta] - z_crit * standard_error[is_beta],
               ignore_attr = TRUE)
  expect_equal(reg_coef[, "z.value"],
               estimate[is_beta] / standard_error[is_beta], ignore_attr = TRUE)
  expect_equal(reg_coef[, "p.value"],
               2 * pnorm(-abs(reg_coef[, "z.value"])))

  # covariance parameters are estimated on the log scale, so their intervals
  # are log-normal
  sp_index <- match(c("sigma2", "phi"), names(estimate))
  expect_equal(fit_summary$sp[, "Estimate"], exp(estimate[sp_index]),
               ignore_attr = TRUE)
  expect_equal(fit_summary$sp[, "Lower limit"],
               exp(estimate[sp_index] - z_crit * standard_error[sp_index]),
               ignore_attr = TRUE)
  expect_equal(fit_summary$sp[, "Upper limit"],
               exp(estimate[sp_index] + z_crit * standard_error[sp_index]),
               ignore_attr = TRUE)

  re_index <- match("sigma2_re.i", names(estimate))
  expect_equal(fit_summary$ranef[, "Upper limit"],
               exp(estimate[re_index] + z_crit * standard_error[re_index]),
               ignore_attr = TRUE)
})

test_that("summary.RiskMap respects level", {
  default_summary <- summary(gaussian_model)
  narrow_summary <- summary(gaussian_model, level = 0.9)

  expect_identical(narrow_summary$level, 0.9)
  expect_equal(narrow_summary$reg_coef[, "Estimate"],
               default_summary$reg_coef[, "Estimate"])
  default_width <- default_summary$reg_coef[, "Upper limit"] -
    default_summary$reg_coef[, "Lower limit"]
  narrow_width <- narrow_summary$reg_coef[, "Upper limit"] -
    narrow_summary$reg_coef[, "Lower limit"]
  expect_equal(narrow_width / default_width,
               rep(qnorm(0.95) / qnorm(0.975), 2), ignore_attr = TRUE)
})

test_that("summary.RiskMap reports fixed and estimated nugget and measurement error", {
  # gaussian_offset_model: gp(nugget = TRUE), sigma2_me held at 0, offset
  offset_summary <- summary(gaussian_offset_model)
  expect_true(offset_summary$tau2)
  expect_true(offset_summary$cov_offset_used)
  expect_identical(offset_summary$held_family_parameters, list(sigma2_me = 0))
  expect_identical(rownames(offset_summary$sp)[3], "Variance of the nugget")
  # the nugget is estimated as a ratio to sigma2 but reported as a variance
  expect_equal(offset_summary$sp["Variance of the nugget", "Estimate"],
               coef(gaussian_offset_model)$tau2)

  # log(tau2) = nu2 + sigma2 on the working scale, so its standard error
  # comes from the delta method with both parameters
  covariance <- gaussian_offset_model$covariance
  gradient <- setNames(numeric(nrow(covariance)), rownames(covariance))
  gradient[c("sigma2", "nu2")] <- 1
  log_tau2_se <- sqrt(drop(gradient %*% covariance %*% gradient))
  tau2 <- coef(gaussian_offset_model)$tau2
  expect_equal(
    offset_summary$sp["Variance of the nugget",
                      c("Lower limit", "Upper limit")],
    exp(log(tau2) + c(-1, 1) * qnorm(0.975) * log_tau2_se),
    ignore_attr = TRUE
  )

  # A nugget without a standard error, as when it is estimated at zero,
  # leaves the intervals of the other parameters unchanged
  no_nugget_se <- gaussian_offset_model
  no_nugget_se$covariance["nu2", ] <- NA
  no_nugget_se$covariance[, "nu2"] <- NA
  missing_summary <- summary(no_nugget_se)
  expect_equal(missing_summary$reg_coef, offset_summary$reg_coef)
  expect_equal(missing_summary$sp[1:2, ], offset_summary$sp[1:2, ])
  expect_true(all(is.na(
    missing_summary$sp["Variance of the nugget",
                       c("Lower limit", "Upper limit")]
  )))

  # gaussian_model: no nugget, estimated measurement error
  gaussian_summary <- summary(gaussian_model)
  expect_identical(gaussian_summary$tau2, 0)
  expect_identical(rownames(gaussian_summary$family_parameters),
                   "Measurement error var.")
  expect_identical(gaussian_summary$held_family_parameters, list())
  expect_equal(gaussian_summary$family_parameters[, "Estimate"],
               coef(gaussian_model)$sigma2_me, ignore_attr = TRUE)
})

# A summary.RiskMap for a gaussian model; 'changes' replaces fields, and
# setting a field to NULL removes it.
build_summary_riskmap <- function(changes = list()) {
  reg_coef <- cbind(
    Estimate      = c(1.2, -0.4),
    "Lower limit" = c(0.8, -0.9),
    "Upper limit" = c(1.6, 0.1),
    StdErr        = c(0.2, 0.25),
    z.value       = c(6, -1.6),
    p.value       = c(1.97e-09, 0.1096)
  )
  rownames(reg_coef) <- c("(Intercept)", "cov")

  me <- cbind(Estimate = 0.3, "Lower limit" = 0.1, "Upper limit" = 0.9)
  rownames(me) <- "Measurement error var."

  sp <- cbind(
    Estimate      = c(0.5, 0.2, 0.1),
    "Lower limit" = c(0.25, 0.05, 0.02),
    "Upper limit" = c(1, 0.8, 0.5)
  )
  rownames(sp) <- c("Spatial process var.", "Spatial corr. scale (km)",
                    "Variance of the nugget")

  ranef <- cbind(Estimate = 0.4, "Lower limit" = 0.2, "Upper limit" = 0.8)
  rownames(ranef) <- "i (random eff. var.)"

  summary_object <- list(
    reg_coef = reg_coef,
    family_parameters = me,
    held_family_parameters = list(),
    sp = sp,
    tau2 = TRUE,
    ranef = ranef,
    level = 0.95,
    family = "gaussian",
    kappa = 0.5,
    log_lik = -12.3456,
    cov_offset_used = TRUE,
    aic = 38.6912,
    call = quote(glgpm(formula = y ~ cov + gp(nugget = TRUE) + re(i),
                       data = gaussian_data, family = "gaussian")),
    link_name = "identity",
    invlink_expression = "Inverse link function = x"
  )
  summary_object <- modifyList(summary_object, changes)
  class(summary_object) <- "summary.RiskMap"
  summary_object
}

test_that("print.summary.RiskMap formats a gaussian model with every optional section", {
  summary_object <- build_summary_riskmap()

  expect_snapshot(print(summary_object))
  expect_output(expect_invisible(print(summary_object)))
})

test_that("print.summary.RiskMap reports a fixed measurement error variance", {
  summary_object <- build_summary_riskmap(list(
    family_parameters = NULL, held_family_parameters = list(sigma2_me = 0.1),
    call = NULL
  ))

  expect_snapshot(print(summary_object))
})

test_that("print.summary.RiskMap formats binomial and poisson models", {
  sp_without_nugget <- build_summary_riskmap()$sp[1:2, ]

  binomial_summary <- build_summary_riskmap(list(
    family = "binomial", family_parameters = NULL, sp = sp_without_nugget, tau2 = 0,
    ranef = NULL, cov_offset_used = FALSE, aic = NULL, call = NULL,
    link_name = "canonical (logit)",
    invlink_expression = "Inverse link function = 1 / (1 + exp(-x))"
  ))
  expect_snapshot(print(binomial_summary))

  poisson_summary <- build_summary_riskmap(list(
    family = "poisson", family_parameters = NULL, sp = sp_without_nugget, tau2 = 0.25,
    ranef = NULL, cov_offset_used = FALSE, aic = NULL, call = NULL,
    link_name = "canonical (log)",
    invlink_expression = "Inverse link function = exp(x)"
  ))
  expect_snapshot(print(poisson_summary))
})

test_that("print.summary.RiskMap reports the certificate and the checks to review", {
  sp_without_nugget <- build_summary_riskmap()$sp[1:2, ]
  certified <- list(certificate = 0.0594, level = 0.95, tolerance = 0.1,
                    certified = TRUE, outcome = "certified",
                    correction_refused = FALSE, stages = 1L,
                    fitting_draws = 1000L, validation_draws = 5000L,
                    checks_to_review = 0L)
  binomial_summary <- build_summary_riskmap(list(
    family = "binomial", family_parameters = NULL, sp = sp_without_nugget, tau2 = 0,
    ranef = NULL, cov_offset_used = FALSE, aic = NULL, call = NULL,
    link_name = "canonical (logit)",
    invlink_expression = "Inverse link function = 1 / (1 + exp(-x))",
    mcml = certified
  ))
  expect_snapshot(print(binomial_summary))

  refused <- modifyList(certified, list(
    certificate = 0.0412, certified = FALSE, outcome = "not certified",
    correction_refused = TRUE, stages = 4L, fitting_draws = 2500L,
    validation_draws = 12500L, checks_to_review = 2L
  ))
  binomial_summary$mcml <- refused
  expect_snapshot(print(binomial_summary))
})

test_that("summary.RiskMap reads the certificate of a non-Gaussian fit", {
  expect_null(summary(gaussian_model)$mcml)
  optimizer <- attr(binomial_model, "optimizer")
  mcml <- summary(binomial_model)$mcml
  expect_identical(mcml$certificate, optimizer$mcml_certificate)
  expect_identical(mcml$certified, isTRUE(optimizer$mcml_certified))
  expect_identical(mcml$outcome, optimizer$mcml_outcome)
  expect_identical(mcml$stages, optimizer$mcml_stages)
  expect_identical(mcml$tolerance,
                   attr(binomial_model, "control_mcml")$tolerance)
  expect_identical(mcml$fitting_draws, optimizer$mcml_budget)
  expect_identical(mcml$validation_draws, optimizer$mcml_validation_budget)
  expect_output(print(summary(binomial_model)),
                "Monte Carlo maximum likelihood")
})

test_that("print.RiskMap is short and gives the certificate", {
  printed <- capture.output(result <- print(binomial_model))
  expect_identical(result, binomial_model)
  expect_lt(length(printed), 30L)
  expect_match(printed[1L], "Binomial geostatistical model")
  expect_true(any(grepl(paste0("^", n, " observations at "), printed)))
  expect_true(any(grepl("Regression coefficients", printed, fixed = TRUE)))
  expect_true(any(grepl("certified", printed, ignore.case = TRUE)))
  expect_false(any(grepl("concerns the maximum found", printed)))
  estimates <- coef(binomial_model)
  variance <- grep("^Spatial process var.", printed, value = TRUE)
  expect_equal(as.numeric(sub("^Spatial process var. +", "", variance)),
               estimates$sigma2, tolerance = 1e-3)
  expect_true(any(grepl("`summary()` gives standard errors", printed,
                        fixed = TRUE)))

  gaussian_printed <- capture.output(print(gaussian_model))
  expect_lt(length(gaussian_printed), 30L)
  expect_true(any(grepl("i (random eff. var.)", gaussian_printed,
                        fixed = TRUE)))
  expect_false(any(grepl("Monte Carlo", gaussian_printed)))
  expect_output(print(poisson_model), "glgpm(formula = y ~ cov", fixed = TRUE)
})

test_that("print.summary.RiskMap labels a summary of a real fit", {
  expect_output(print(summary(gaussian_model)), "Linear geostatistical model")
  expect_output(print(summary(gaussian_model)), "Unstructured random effects")
  expect_output(print(summary(binomial_model)), "Binomial geostatistical model")
  expect_output(print(summary(poisson_model)), "Poisson geostatistical model")
})

# summary.RiskMap_cross_validation -------------------------------------------

# A cross-validation of two models on two test sets, of rows 1 to 3 and 3 to 4,
# so that row 3 is scored twice.
build_cross_validation <- function(...) {
  scores <- data.frame(
    model = rep(c("model_a", "model_b"), each = 5),
    test_set = rep(c(1L, 1L, 1L, 2L, 2L), 2), row = rep(c(1L, 2L, 3L, 3L, 4L), 2),
    CRPS = c(0.2, 0.3, 0.4, 0.6, 0.5, 0.6, 0.1, 0.2, 0.4, 0.3),
    SCRPS = c(-0.5, -0.6, -0.7, -0.7, -0.8, -1, -0.2, -0.3, -0.5, -0.4),
    pit_lower = 0.5, pit_upper = 0.5
  )
  object <- list(
    scores = scores,
    calibration = list(area = data.frame(
      model = c("model_a", "model_b"), PIT_area = c(0.04, 0.05),
      PIT_area_expected = c(0.03, 0.031), PIT_area_limit = c(0.06, 0.07)
    )),
    test_sets = list(sets = list(list(train = 4:10, test = 1:3),
                                 list(train = 5:10, test = 3:4))),
    refit = FALSE, family = c(model_a = "gaussian", model_b = "gaussian"),
    failures = data.frame(model = character(), test_set = integer(),
                          message = character())
  )
  changed <- list(...)
  object[names(changed)] <- changed
  structure(object, class = "RiskMap_cross_validation")
}

test_that("the summary of a cross-validation averages over the observations", {
  scores <- summary(build_cross_validation())
  expect_s3_class(scores, c("summary.RiskMap_cross_validation", "data.frame"),
                  exact = TRUE)
  expect_identical(row.names(scores), c("model_a", "model_b"))
  expect_identical(names(scores),
                   c("CRPS", "SCRPS", "PIT_area", "PIT_area_expected",
                     "PIT_area_limit"))
  # Row 3 counts once, with the average of its scores.
  expect_equal(scores["model_a", "CRPS"], mean(c(0.2, 0.3, 0.5, 0.5)))
  expect_equal(scores["model_b", "SCRPS"], mean(c(-1, -0.2, -0.4, -0.4)))
  expect_identical(scores$PIT_area, c(0.04, 0.05))
  expect_identical(scores$PIT_area_expected, c(0.03, 0.031))
  expect_identical(scores$PIT_area_limit, c(0.06, 0.07))
  expect_equal(attr(scores, "n_observations"), c(model_a = 4, model_b = 4))
  expect_identical(attr(scores, "n_failed"), c(model_a = 0L, model_b = 0L))
  expect_false(attr(scores, "refit"))
  expect_null(attr(scores, "reference"))
})

test_that("the summary of a cross-validation compares with a reference", {
  scores <- summary(build_cross_validation(), reference = "model_a")
  expect_identical(names(scores),
                   c("CRPS", "CRPS_difference", "CRPS_se", "SCRPS",
                     "SCRPS_difference", "SCRPS_se", "PIT_area",
                     "PIT_area_expected", "PIT_area_limit"))
  difference <- c(0.4, -0.2, -0.2, -0.2)
  expect_equal(scores["model_b", "CRPS_difference"], mean(difference))
  expect_equal(scores["model_b", "CRPS_se"], sd(difference) / 2)
  expect_equal(unlist(scores["model_a", c("CRPS_difference", "CRPS_se")]),
               c(CRPS_difference = 0, CRPS_se = 0))
  expect_identical(attr(scores, "reference"), "model_a")
  # One observation has no standard error.
  one <- build_cross_validation()
  one$scores <- one$scores[one$scores$row == 1L, ]
  expect_true(is.na(summary(one, reference = "model_a")["model_b", "CRPS_se"]))

  expect_input_error(summary(build_cross_validation(), reference = "model_c"),
                     "reference", "\"model_a\" and \"model_b\"")
  expect_input_error(summary(build_cross_validation(), reference = 1),
                     "reference")
  expect_input_error(summary(build_cross_validation(),
                             reference = c("model_a", "model_b")), "reference")
  expect_input_error(summary(build_cross_validation(), view_all = FALSE),
                     "view_all", "The scores of each test set")
})

test_that("the summary of a cross-validation leaves out failed test sets", {
  failed <- build_cross_validation(
    refit = TRUE, failures = data.frame(model = "model_a", test_set = 2L,
                                        message = "failed")
  )
  missing <- failed$scores$model == "model_a" & failed$scores$test_set == 2L
  failed$scores[missing, c("CRPS", "SCRPS")] <- NA
  scores <- summary(failed, reference = "model_b")
  expect_equal(scores["model_a", "CRPS"], mean(c(0.2, 0.3, 0.4)))
  # The differences are those of the observations both models have scores for.
  difference <- c(0.2 - 0.6, 0.3 - 0.1, 0.4 - 0.3)
  expect_equal(scores["model_a", "CRPS_difference"], mean(difference))
  expect_equal(scores["model_a", "CRPS_se"], sd(difference) / sqrt(3))
  expect_identical(attr(scores, "n_failed"), c(model_a = 1L, model_b = 0L))
  # The models have scores for different numbers of observations.
  expect_equal(attr(scores, "n_observations"), c(model_a = 3, model_b = 4))
  printed <- squish(capture.output(print(scores)))
  expect_match(printed, "Scores over 3 to 4 test observations", fixed = TRUE)
  expect_match(printed, "Model Scored CRPS", fixed = TRUE)
  expect_match(printed, "model_a 3 ", fixed = TRUE)
  expect_no_match(capture.output(print(scores["model_b", ])), "Scored")
  expect_match(printed, "Parameters are re-estimated in each training set.",
               fixed = TRUE)
  expect_match(printed,
               "1 test-set assessment of model_a failed, and its scores are missing.",
               fixed = TRUE)

  # Nothing was scored, printed as NA and not NaN.
  failed$scores$CRPS <- NA_real_
  scores <- summary(failed, reference = "model_b")
  expect_identical(scores$CRPS, c(NA_real_, NA_real_))
  expect_identical(cv_summary_table(scores, 4L)$CRPS, c("NA", "NA"))
  expect_output(print(scores), "Scores over 0 test observations")
  expect_identical(cv_summary_table(scores, 4L)[["CRPS difference"]],
                   c("NA", ""))
})

test_that("the note on the limit of the area needs a limit", {
  without <- build_cross_validation()
  expect_output(print(summary(without)), "PIT area of a correct model")
  without$calibration$area$PIT_area_limit <- NA_real_
  printed <- capture.output(print(summary(without)))
  expect_no_match(printed, "PIT area of a correct model")
})

test_that("the summary of a cross-validation by assess_prediction() is finite", {
  scores <- summary(cross_validation, reference = "model_a")
  expect_identical(row.names(scores), c("model_a", "model_b"))
  expect_true(all(is.finite(as.matrix(as.data.frame(unclass(scores))))))
  expect_equal(attr(scores, "n_observations"), c(model_a = 10, model_b = 10))
})

test_that("print() of a cross-validation gives the averages, not the raw list", {
  printed <- capture.output(result <- print(cross_validation))
  expect_identical(result, cross_validation)
  expect_match(printed[1L], "Cross-validation", fixed = TRUE)
  expect_match(squish(printed), "2 models assessed on 2 test sets of 4 and 6 observations. Parameters are fixed at the estimates from all the data.", fixed = TRUE)
  expect_true(any(startsWith(trimws(printed), "model_a ")))
  expect_true(any(startsWith(trimws(printed), "model_b ")))
  expect_match(squish(printed), "`summary()` with `reference` compares the models with one of them, and `plot()` draws the calibration curves or maps the scores.", fixed = TRUE)
  expect_match(squish(printed), "Lower is better for the CRPS and the SCRPS.", fixed = TRUE)
  expect_lt(length(printed), 24L)
  # One model isn't compared.
  alone <- build_cross_validation()
  alone$scores <- alone$scores[alone$scores$model == "model_a", ]
  alone$calibration$area <- alone$calibration$area[1, ]
  alone$family <- c(model_a = "poisson")
  printed <- squish(capture.output(print(alone)))
  expect_match(printed, "1 model assessed on 2 test sets of 3 and 2 observations")
  expect_no_match(printed, "reference", fixed = TRUE)
  expect_match(printed, "the response, which are counts", fixed = TRUE)
  expect_input_error(print(cross_validation, digits = 2), "digits")
})

test_that("the print of a cross-validation describes its test sets", {
  expect_identical(describe_test_sets(c(80L, 20L)),
                   "2 test sets of 80 and 20 observations")
  expect_identical(describe_test_sets(rep(50L, 10L)),
                   "10 test sets of 50 observations each")
  expect_identical(describe_test_sets(1L), "1 test set of 1 observation")
  expect_identical(describe_test_sets(rep(1L, 20L)),
                   "20 test sets of 1 observation each")
  expect_identical(describe_test_sets(c(8L, 12L, 9L, 10L)),
                   "4 test sets of 8 to 12 observations")
})

test_that("print() of the summary of a cross-validation formats its scores", {
  scores <- summary(build_cross_validation())
  expect_snapshot(print(scores))
  compared <- summary(build_cross_validation(), reference = "model_a")
  expect_snapshot(print(compared))
  expect_output(expect_invisible(print(scores)))
  expect_input_error(print(scores, digits = 2), "digits")
  # A document that prints data frames as tables prints it as the console does.
  expect_output(knitr::knit_print(scores), "Cross-validation scores")
  # Some of the models are a summary too, and some of the columns a data frame.
  printed <- squish(capture.output(print(compared["model_b", ])))
  expect_match(printed, "model_b 0.3250 -0.0500 (0.1500)", fixed = TRUE)
  expect_no_match(printed, "model_a 0.3750", fixed = TRUE)
  expect_match(printed, "Differences are from model_a", fixed = TRUE)
  expect_no_match(capture.output(print(scores[c("CRPS", "SCRPS")], digits = 2)),
                  "Cross-validation")
  expect_input_error(to_table(scores["CRPS"]), "object", "whole summary")
})

# summary.RiskMap_assess_simulation -------------------------------------------

# A RiskMap_assess_simulation with two simulations of one model; summary()
# of it is deterministic, so it is used to build the object being printed.
build_assess_simulation <- function(spatial_scale) {
  classify_result <- function(sensitivity, cc) {
    data.frame(Class = c("low", "high"), Sensitivity = sensitivity,
               Specificity = c(0.8, 0.7), PPV = c(0.6, 0.5),
               NPV = c(0.9, 0.85), CC = cc)
  }
  scale_results <- list(
    mse = matrix(c(0.12, 0.18), nrow = 1,
                 dimnames = list("model", c("sim_1", "sim_2"))),
    classify = list(
      model = list(by_cat = list(classify_result(c(0.5, 0.6), c(0.7, 0.7)),
                                 classify_result(c(0.7, 0.4), c(0.8, 0.8))),
                   CC = c(0.7, 0.8)),
      Class = factor(c("low", "high"), levels = c("low", "high"))
    )
  )
  objective <- rep(list(scale_results), length(spatial_scale))
  names(objective) <- spatial_scale

  structure(
    list(objective = objective, n_sim = 2L,
         spatial_scale = spatial_scale),
    class = "RiskMap_assess_simulation"
  )
}

test_that("classification summaries handle one simulation and metric-wise missingness #108", {
  one_result <- data.frame(
    Class = c("low", "high"),
    Sensitivity = c(0.5, NA),
    Specificity = c(NA, 0.75),
    PPV = c(0.4, NA),
    NPV = c(NA, 0.8),
    CC = c(0.6, 0.7)
  )
  one_simulation <- structure(
    list(
      objective = list(
        location = list(
          classify = list(
            model = list(by_cat = list(one_result), CC = 0.65),
            Class = factor(c("low", "high"), levels = c("low", "high"))
          )
        )
      ),
      n_sim = 1L,
      spatial_scale = "location"
    ),
    class = "RiskMap_assess_simulation"
  )

  expect_warning(summary_one <- summary(one_simulation), "one simulation")
  model_one <- summary_one$location$classify$model
  expect_identical(model_one$classify_res$Sensitivity, c(0.5, NA))
  expect_identical(model_one$n_valid$Sensitivity, c(1L, 0L))
  expect_identical(model_one$cc_summary$mean, 0.65)
  expect_true(is.na(model_one$cc_summary$sd))
  expect_true(is.na(model_one$cc_summary$lower))
  expect_output(print(summary_one), "uncertainty unavailable")

  second_result <- one_result
  second_result$Sensitivity <- c(NA, 0.25)
  second_result$Specificity <- c(0.5, NA)
  two_simulations <- one_simulation
  two_simulations$n_sim <- 2L
  two_simulations$objective$location$classify$model$by_cat <-
    list(one_result, second_result)
  two_simulations$objective$location$classify$model$CC <- c(0.65, 0.75)

  summary_two <- summary(two_simulations)$location$classify$model
  expect_identical(summary_two$classify_res$Sensitivity, c(0.5, 0.25))
  expect_identical(summary_two$n_valid$Sensitivity, c(1L, 1L))
  expect_identical(summary_two$n_valid$Specificity, c(1L, 1L))
  expect_identical(summary_two$cc_summary$n_valid, 2L)
  expect_identical(summary_two$cc_summary$n_sim, 2L)
})

test_that("summary.RiskMap_assess_simulation summarises MSE and classification", {
  assessment <- build_assess_simulation(c("location", "area"))
  assessment$objective$location$mse[1, 2] <- NA

  expect_error(summary.RiskMap_assess_simulation(list()), "RiskMap_assess_simulation")

  simulation_summary <- summary(assessment)
  expect_s3_class(simulation_summary, "summary.RiskMap_assess_simulation")
  expect_named(simulation_summary, c("location", "area"))

  # non-finite MSE values are dropped from the mean and counted as invalid
  expect_identical(simulation_summary$location$mse$n_sim, 2L)
  expect_identical(simulation_summary$location$mse$n_valid, 1)
  expect_identical(simulation_summary$location$mse$MSE_mean, 0.12)
  expect_true(is.na(simulation_summary$location$mse$MSE_sd))
  expect_equal(simulation_summary$area$mse$MSE_mean, 0.15)
  expect_equal(simulation_summary$area$mse$MSE_sd, sd(c(0.12, 0.18)))

  model_summary <- simulation_summary$location$classify$model
  expect_identical(model_summary$classify_res$Class, c("low", "high"))
  expect_equal(model_summary$classify_res$Sensitivity, c(0.6, 0.5))
  expect_equal(model_summary$cc_summary$mean, 0.75)
  expect_equal(model_summary$cc_summary$lower,
               unname(quantile(c(0.7, 0.8), 0.025)))

  assessment$objective$location$mse <- c(0.12, 0.18)
  expect_error(summary(assessment), "must be a matrix")
})

test_that("print.summary.RiskMap_assess_simulation formats one spatial scale", {
  summary_object <- summary(build_assess_simulation("location"))

  expect_snapshot(print(summary_object))
  expect_output(expect_invisible(print(summary_object)))
})

test_that("print.summary.RiskMap_assess_simulation labels each spatial scale", {
  summary_object <- summary(build_assess_simulation(c("location", "area")))

  expect_snapshot(print(summary_object))
})

test_that("the summary of a fit says where its starting values came from", {
  expect_match(describe_starting_values("user"),
               "given in `start_pars`, with RiskMap's defaults for any others.",
               fixed = TRUE)
  expect_match(describe_starting_values("transformed_gaussian"),
               "Gaussian fit to the transformed data", fixed = TRUE)
  expect_match(describe_starting_values("current"), "defaults", fixed = TRUE)
  # A fit saved by an earlier version doesn't record it.
  expect_null(describe_starting_values(NULL))
  fitted <- summary(binomial_model)
  selected <- attr(binomial_model, "starting_values")$selected
  expect_false(is.null(selected))
  expect_identical(fitted$mcml$starting_values, selected)
  expect_match(squish(capture.output(print(fitted))),
               describe_starting_values(selected, fitted$mcml$starting_reason),
               fixed = TRUE)
  expect_null(summary(gaussian_model)$mcml)
  # A start given in part keeps the comparison for the rest, and the summary
  # names the values given.
  given <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, start_pars = list(beta = c(0, 0.3)),
    control_mcml = fixture_control_mcml, seed = fixture_seed, messages = FALSE
  ))
  expect_identical(summary(given)$mcml$starting_given,
                   c("beta.(Intercept)", "beta.cov"))
  expect_true(summary(given)$mcml$starting_values %in%
                c("current", "transformed_gaussian"))
  expect_match(squish(capture.output(print(summary(given)))),
               "The starting values of `beta.(Intercept)` and `beta.cov` were given in `start_pars`.",
               fixed = TRUE)
  # A fit saved before starts could be given in part says so.
  expect_match(describe_starting_values("user"),
               "with RiskMap's defaults for any others", fixed = TRUE)
})

test_that("the summary of a fit says why it kept the default starting values", {
  expect_match(describe_starting_values("transformed_gaussian"),
               "came from a Gaussian fit to the transformed data, which the Laplace approximation of the likelihood preferred",
               fixed = TRUE)
  expect_match(describe_starting_values("current"),
               "defaults, which the Laplace approximation of the likelihood preferred",
               fixed = TRUE)
  expect_match(describe_starting_values("current", "re_terms"),
               "a model with `re()` terms takes without comparing them",
               fixed = TRUE)
  expect_match(describe_starting_values("current", "link"),
               "a link other than the canonical one", fixed = TRUE)
  expect_match(describe_starting_values("current", "transformed_failed"),
               "because a Gaussian fit to the transformed data failed.",
               fixed = TRUE)
  # Every reason that select_nongaussian_start() records has its own words.
  reasons <- c("re_terms", "link", "transformed_failed",
               "transformed_invalid", "comparison_failed",
               "transformed_not_finite", "neither_finite")
  sentences <- vapply(reasons, describe_starting_values, character(1),
                      selected = "current")
  expect_identical(anyDuplicated(sentences), 0L)
  expect_false(any(sentences == "The starting values were RiskMap's defaults."))
  # A fit saved by an earlier version doesn't record them.
  expect_null(describe_starting_values(NULL))
  expect_identical(describe_starting_values("current", "not_recorded"),
                   "The starting values were RiskMap's defaults.")
  # One saved after the comparison was added but before reasons were, with
  # re() terms, kept the defaults with a fallback reason and no other.
  old <- binomial_model
  attr(old, "starting_values") <- list(
    selected = "current", current_laplace = -100, transformed_laplace = NA,
    fallback_reason = "transformed Gaussian candidate had a non-finite Laplace value"
  )
  expect_identical(summary(old)$mcml$starting_reason, "not_recorded")
  expect_match(squish(capture.output(print(summary(old)))),
               "The starting values were RiskMap's defaults.", fixed = TRUE)
  # One where the comparison chose the defaults records neither.
  attr(old, "starting_values")$fallback_reason <- NULL
  expect_null(summary(old)$mcml$starting_reason)

  # The fixture has an re() term.
  fitted <- summary(binomial_model)
  expect_identical(fitted$mcml$starting_values, "current")
  expect_identical(fitted$mcml$starting_reason, "re_terms")
  expect_match(squish(capture.output(print(fitted))),
               describe_starting_values("current", "re_terms"), fixed = TRUE)
  expect_null(summary(gaussian_model)$mcml)
})
