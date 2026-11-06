# A fixture with some of the flags that glgpm() writes to the optimizer attribute.
with_flags <- function(fit = binomial_model, ...) {
  optimizer <- attr(fit, "optimizer")
  flags <- list(...)
  optimizer[names(flags)] <- flags
  attr(fit, "optimizer") <- optimizer
  fit
}

find_check <- function(fit, name) {
  checks <- check_fit(fit)$checks
  checks[checks$check == name, ]
}

certified_mcml <- list(certificate = 0.0594, tolerance = 0.1,
                       outcome = "certified", correction_refused = FALSE)

test_that("check_fit lists the checks of a fit with a status for each", {
  expect_error(check_fit(1), "RiskMap")
  checked <- check_fit(binomial_model)
  expect_s3_class(checked, "RiskMap_check")
  expect_named(checked$checks, c("check", "status", "finding", "fix"))
  expect_contains(checked$checks$check,
                  c("Monte Carlo accuracy", "Other maxima", "Variances",
                    "Importance distribution"))
  expect_contains(c("PROBLEM", "SUSPICIOUS", "OK"), checked$checks$status)
  expect_identical(checked$tolerance,
                   attr(binomial_model, "control_mcml")$tolerance)
  expect_identical(summary(binomial_model)$mcml$checks_to_review,
                   sum(checked$checks$status != "OK"))
})

test_that("a Gaussian fit has one check and nothing to plot", {
  checked <- check_fit(gaussian_model)
  expect_identical(nrow(checked$checks), 1L)
  expect_identical(checked$checks$status, "OK")
  expect_match(checked$checks$finding, "no Monte Carlo error")
  expect_null(checked$stages)
  expect_null(checked$draws)
  expect_output(print(checked), "Checks on this fit")
  expect_error(plot(checked), "no stages to plot")
})

test_that("those that need a look come first", {
  checks <- check_fit(with_flags(
    mcml_maxima_alarm = TRUE, mcml_boundary = "nugget variance (tau2)"
  ))$checks
  expect_identical(checks$status[1:2], c("PROBLEM", "SUSPICIOUS"))
  expect_identical(checks$check[1:2], c("Other maxima", "Variances"))
  expect_false(is.unsorted(match(checks$status,
                                 c("PROBLEM", "SUSPICIOUS", "OK"))))
  expect_identical(rownames(checks), as.character(seq_len(nrow(checks))))
})

test_that("the Monte Carlo check follows the outcome of the fit", {
  certified <- fit_check_monte_carlo(certified_mcml)
  expect_identical(certified$status, "OK")
  expect_identical(certified$finding,
                   "The error of the estimate is within 0.059 standard errors (limit 0.1).")
  expect_identical(certified$fix, "")

  above <- fit_check_monte_carlo(modifyList(certified_mcml, list(
    certificate = 1.66, outcome = "not certified"
  )))
  expect_identical(above$status, "PROBLEM")
  expect_match(above$finding, "1.7 standard errors, above the limit of 0.1")
  expect_match(above$fix, "'max_samples' or 'max_stages'")

  refused <- fit_check_monte_carlo(modifyList(certified_mcml, list(
    outcome = "not certified", correction_refused = TRUE
  )))
  expect_identical(refused$status, "PROBLEM")
  expect_match(refused$finding, "did not converge")
  expect_no_match(refused$finding, "above the limit")

  for (outcome in c("weak identification", "numerically independent")) {
    none <- fit_check_monte_carlo(modifyList(certified_mcml, list(
      certificate = NA_real_, outcome = outcome
    )))
    expect_identical(none$status, "PROBLEM")
    expect_match(none$finding, "no bound")
    expect_match(none$fix, "spatial scale")
  }
})

test_that("the weights check flags a heavy tail", {
  light <- fit_check_weights(list(effective_draws = 1081.4, pareto_k = 0.46))
  expect_identical(light$status, "OK")
  expect_identical(
    light$finding,
    "The certification draws had 1,081 effective draws and a tail shape of 0.46 (limit 0.7)."
  )
  heavy <- fit_check_weights(list(effective_draws = 94, pareto_k = 1.005))
  expect_identical(heavy$status, "SUSPICIOUS")
  expect_match(heavy$finding, "tail shape of 1 (limit 0.7), which", fixed = TRUE)
  expect_match(heavy$fix, "'max_samples'")
  expect_identical(fit_check_weights(list(effective_draws = 1, pareto_k = NA))$status,
                   "SUSPICIOUS")
  expect_null(fit_check_weights(NULL))
})

test_that("the maxima check says what the check for other maxima found", {
  none <- fit_check_maxima(list())
  expect_identical(none$status, "OK")
  expect_match(none$finding, "cannot rule one out")

  alarm <- fit_check_maxima(list(mcml_maxima_alarm = TRUE))
  expect_identical(alarm$status, "PROBLEM")
  expect_match(alarm$finding, "not the maximum likelihood estimate")
  expect_match(alarm$fix, "profile(fit, parameters = \"phi\")", fixed = TRUE)

  check <- data.frame(
    maximum = 1:3, independent_log_likelihood = c(-12.1, -10.2, -10.3),
    independent_se = c(0.05, 0, 0.04), chosen = c(FALSE, TRUE, FALSE),
    log_phi = log(c(0.1, 20, 3))
  )
  better <- fit_check_maxima(list(mcml_better_maximum = TRUE,
                                  mcml_maxima_check = check), "km")
  expect_identical(better$status, "SUSPICIOUS")
  expect_identical(
    better$finding,
    paste0("The check found a maximum higher than the one reached from the ",
           "starting values, and the estimate is from it. The two differ by 1.9 ",
           "log-likelihood units (standard error 0.05) on independent draws, ",
           "and the other has a spatial scale of 0.1 km.")
  )
  # An ambiguity is with the closest rival, not the maximum first fitted.
  ambiguous <- fit_check_maxima(list(mcml_ambiguous_maxima = TRUE,
                                     mcml_maxima_check = check), "km")
  expect_identical(ambiguous$status, "SUSPICIOUS")
  expect_match(ambiguous$finding, "about equally well")
  expect_match(ambiguous$finding, "differ by 0.1 .* scale of 3 km")
  expect_match(ambiguous$fix, "profile(fit)", fixed = TRUE)
  expect_identical(describe_other_maximum(NULL, FALSE), "")
})

test_that("the variances check separates the spatial variance from the others", {
  expect_identical(fit_check_variances(list())$status, "OK")

  nugget <- fit_check_variances(list(mcml_boundary = "nugget variance (tau2)"))
  expect_identical(nugget$status, "SUSPICIOUS")
  expect_match(nugget$finding, "was estimated as effectively zero, so its standard error")
  expect_match(nugget$fix, "gp(nugget = FALSE)", fixed = TRUE)

  spatial <- fit_check_variances(list(mcml_boundary = "spatial variance (sigma2)",
                                      mcml_pc_prior = FALSE))
  expect_identical(spatial$status, "PROBLEM")
  expect_match(spatial$finding, "spatial scale cannot be estimated")
  expect_match(spatial$fix, "pc_prior = TRUE")

  two <- fit_check_variances(list(
    mcml_boundary = c("spatial variance (sigma2)", "nugget variance (tau2)"),
    mcml_pc_prior = TRUE
  ))
  expect_match(two$finding, "were estimated as effectively zero.*their standard errors")
  expect_identical(two$fix, "Consider refitting without them.")
})

test_that("the scale check says whether the data determine the scale", {
  phi <- signif(coef(binomial_model)$phi, 3)
  long <- fit_check_spatial_scale(binomial_model, list())
  expect_identical(long$status, "OK")
  expect_match(long$finding, paste0("The scale of ", phi), fixed = TRUE)

  short <- fit_check_spatial_scale(
    binomial_model,
    list(mcml_small_scale = TRUE, mcml_spacing = 12.34, mcml_pc_prior = FALSE)
  )
  expect_identical(short$status, "SUSPICIOUS")
  expect_match(short$finding, "shorter than the median distance")
  expect_match(short$finding, "(12.3 ", fixed = TRUE)
  expect_match(short$fix, "profile(fit, parameters = \"phi\")", fixed = TRUE)
  expect_match(short$fix, "pc_prior = TRUE")
  penalised <- fit_check_spatial_scale(
    binomial_model, list(mcml_small_scale = TRUE, mcml_pc_prior = TRUE)
  )
  expect_no_match(penalised$fix, "pc_prior")

  near <- fit_check_spatial_scale(binomial_model, list(mcml_near_limit = TRUE))
  expect_identical(near$status, "SUSPICIOUS")
  expect_match(near$finding, "not clearly above the limit")

  independent <- fit_check_spatial_scale(
    binomial_model, list(mcml_numerically_independent = TRUE)
  )
  expect_identical(independent$status, "PROBLEM")
  expect_match(independent$finding, "unresolved")
  weak <- fit_check_spatial_scale(binomial_model,
                                  list(mcml_identifiability_failure = TRUE))
  expect_identical(weak$status, "PROBLEM")
  expect_match(weak$finding, "not identifiable")
  expect_null(fit_check_spatial_scale(
    binomial_model, list(mcml_boundary = "spatial variance (sigma2)")
  ))
})

test_that("changes of importance distribution and the mode limit stay OK", {
  plain <- fit_check_importance_distribution(list(mcml_proposal = "global"))
  expect_identical(plain$status, "OK")
  expect_identical(plain$finding, "The global distribution was used at every stage.")
  expect_identical(
    fit_check_importance_distribution(list(mcml_proposal = "laplace"))$finding,
    "The Laplace approximation was used at every stage."
  )
  restarted <- fit_check_importance_distribution(
    list(mcml_proposal = "global", mcml_start_replaced = TRUE)
  )
  expect_identical(restarted$status, "OK")
  expect_match(restarted$finding,
               "restarted from the starting values and used it at every stage")
  later <- fit_check_importance_distribution(
    list(mcml_start_replaced = TRUE, mcml_proposal_fallbacks = 2L)
  )
  expect_match(later$finding, "starting values. The global distribution could not")
  fallback <- fit_check_importance_distribution(
    list(mcml_proposal = "global", mcml_proposal_fallbacks = c(1L, 3L))
  )
  expect_identical(fallback$status, "OK")
  expect_match(fallback$finding, "at stages 1, 3, which used the Laplace")

  expect_null(fit_check_mode_search(NULL))
  expect_null(fit_check_mode_search(list(searches = 40L, at_limit = 0L)))
  limit <- fit_check_mode_search(list(searches = 4000L, at_limit = 3L))
  expect_identical(limit$status, "OK")
  expect_match(limit$finding, "limit in 3 of 4,000 searches", fixed = TRUE)
})

test_that("a fit without a maximum gets that check alone, with what to do", {
  everywhere <- check_fit(with_flags(
    mcml_no_maximum = TRUE, mcml_outcome = "no maximum",
    mcml_separation = list(everywhere = "no positives")
  ))$checks
  expect_identical(nrow(everywhere), 1L)
  expect_identical(everywhere$status, "PROBLEM")
  expect_match(everywhere$finding, "no individual tested is positive")
  expect_match(everywhere$fix, "exact upper confidence limit for the prevalence")

  counts <- fit_check_no_maximum(list(mcml_separation = list(everywhere = "no positives")),
                                 "poisson")
  expect_match(counts$finding, "every count is zero")
  expect_match(counts$fix, "for the rate")

  only <- fit_check_no_maximum(list(mcml_separation = list(everywhere = "only positives")),
                               "binomial")
  expect_match(only$finding, "every individual tested is positive")
  expect_match(only$fix, "exact lower confidence limit")

  separated <- fit_check_no_maximum(
    list(mcml_separation = list(coefficients = "urban", no_positives = 3L,
                                only_positives = 0L)),
    "binomial"
  )
  expect_match(separated$finding, "covariates separate the data, through 'urban'")
  expect_match(separated$fix, "that has no positives with other categories")

  bound <- fit_check_no_maximum(list(), "binomial")
  expect_match(bound$finding, "reached the bound of the search")
  expect_match(bound$fix, "Check whether a covariate")
})

test_that("the summary counts the checks that need a look", {
  expect_null(describe_checks_to_review(NULL))
  expect_identical(describe_checks_to_review(0L),
                   "No check needs a look, see check_fit() for the details.")
  expect_identical(describe_checks_to_review(1L),
                   "1 check needs a look, see check_fit().")
  expect_identical(describe_checks_to_review(2L),
                   "2 checks need a look, see check_fit().")
  flagged <- with_flags(mcml_maxima_alarm = FALSE,
                        mcml_near_limit = TRUE,
                        mcml_boundary = "nugget variance (tau2)")
  count <- sum(check_fit(flagged)$checks$status != "OK")
  expect_gte(count, 2L)
  expect_identical(summary(flagged)$mcml$checks_to_review, count)
  expect_output(print(summary(flagged)), paste0(count, " checks need a look"))
  expect_output(print(flagged), paste0(count, " checks need a look"))
  expect_no_match(capture.output(print(gaussian_model)), "check_fit")
})

test_that("the stages and the moves of the estimates are summarised", {
  fit <- binomial_model
  estimate <- as.numeric(unlist(fit$estimate))
  errors <- seq_along(estimate) / 10
  fit$covariance <- diag(errors^2, length(estimate))
  rownames(fit$covariance) <- names(unlist(fit$estimate))
  later <- estimate + 2 * errors
  fit$mcml_history <- list(
    list(proposal = "laplace", certificate = 0.5, estimate = estimate,
         validation_budget = 100L,
         second_validation = list(relative_ess = 0.1, pareto_k = 1.1)),
    list(proposal = "global", certificate = 0.05, estimate = later,
         validation_budget = 1000L,
         second_validation = list(relative_ess = 0.5, pareto_k = 0.3)),
    list(proposal = "global", certificate = NA_real_, estimate = later)
  )
  checked <- check_fit(fit)
  expect_identical(checked$stages$stage, 1:3)
  expect_identical(checked$stages$proposal, c("laplace", "global", "global"))
  expect_identical(checked$stages$certificate, c(0.5, 0.05, NA))
  expect_equal(checked$stages$effective_draws, c(10, 500, NA), tolerance = 1e-12)
  expect_identical(checked$stages$pareto_k, c(1.1, 0.3, NA))

  expect_identical(unique(checked$moves$stage), 2:3)
  expect_identical(checked$moves$parameter[seq_along(estimate)],
                   names(unlist(fit$estimate)))
  expect_equal(checked$moves$move[checked$moves$stage == 2L],
               rep(2, length(estimate)), tolerance = 1e-12)
  expect_equal(checked$moves$move[checked$moves$stage == 3L],
               rep(0, length(estimate)), tolerance = 1e-12)

  # A parameter without a standard error is left out, and a fit with none has
  # no moves to show, as has one with a single stage.
  held <- fit
  held$covariance[1L, ] <- NA_real_
  held$covariance[, 1L] <- NA_real_
  expect_disjoint(check_fit(held)$moves$parameter, rownames(fit$covariance)[1L])
  held$covariance[] <- NA_real_
  expect_null(check_fit(held)$moves)
  fit$mcml_history <- fit$mcml_history[1L]
  expect_null(check_fit(fit)$moves)
})

test_that("printing gives a line for each check and the fix where there is one", {
  flagged <- with_flags(mcml_boundary = "nugget variance (tau2)")
  checked <- check_fit(flagged)
  printed <- capture.output(result <- print(checked))
  expect_identical(result, checked)
  expect_identical(printed[1L], "Checks on this fit")
  expect_match(printed[2L], "^Fitted in ")
  expect_true(any(grepl("^SUSPICIOUS +Variances +The nugget", printed)))
  expect_true(any(grepl("^OK +Monte Carlo accuracy", printed)))
  expect_true(any(grepl("^ +Fix: Consider refitting without it", printed)))
  expect_match(printed[length(printed)], "not whether the model suits the data")
  expect_false(any(grepl("^ +$", printed)))
  withr::local_options(width = 60L)
  narrow <- capture.output(print(checked))
  expect_lte(max(nchar(narrow)), 80L)
})

test_that("the plot shows the stages, and a profile if given", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  checked <- check_fit(binomial_model)
  panels <- plot(checked)
  expect_contains(names(panels), "certificate")
  for (panel in panels) expect_s3_class(panel, "ggplot")
  expect_error(plot(checked, profile = 1), "'profile' must be a profile")

  scale_profile <- profile(binomial_model, parameters = "phi", n_points = 5,
                           n_draws = 200, messages = FALSE)
  with_profile <- plot(checked, profile = scale_profile)
  expect_identical(names(with_profile), c(names(panels), "profile"))
  expect_s3_class(with_profile$profile, "ggplot")
  variance_profile <- profile(binomial_model, parameters = "sigma2", n_points = 5,
                              n_draws = 200, messages = FALSE)
  expect_error(plot(checked, profile = variance_profile),
               "'profile' must be a profile")

  stages <- data.frame(stage = 1:2, certificate = c(NA, 0.05))
  expect_null(plot_stage_panel(stages[1L, ], "certificate", "t", "y"))
  expect_s3_class(plot_stage_panel(stages, "certificate", "t", "y", limit = 0.1,
                                   log_scale = TRUE), "ggplot")
  expect_null(plot_estimate_moves(NULL))
  moves <- data.frame(parameter = c("a", "b", "a", "b"), stage = c(2L, 2L, 3L, 3L),
                      move = c(1, 2, 0.1, 0.2))
  expect_s3_class(plot_estimate_moves(moves), "ggplot")
})
