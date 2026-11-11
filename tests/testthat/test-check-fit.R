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
  expect_named(checked$checks, c("check", "status", "finding", "fix", "brief"))
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
  expect_null(checked$fit_seconds)
  expect_output(print(checked), "Checks on this fit")
  expect_error(plot(checked), "no stages to plot.*extra = \"profile\"")
  expect_error(plot(checked, which = "certificate"), "no stages to plot")
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
  expect_identical(alarm$fix, paste0("Inspect the profile likelihood, ",
                                     "check_fit(fit, extra = \"profile\")."))

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
  expect_match(ambiguous$fix,
               "see the profile likelihood, check_fit(fit, extra = \"profile\").",
               fixed = TRUE)
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
  expect_match(short$fix, paste0("^See the profile likelihood, ",
                                 "check_fit\\(fit, extra = \"profile\"\\)\\. "))
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
  expect_null(fit_check_mode_search(
    list(searches = 40L, at_limit = 0L, restarted = 0L)
  ))
  restarted <- fit_check_mode_search(
    list(searches = 4000L, at_limit = 0L, restarted = 2L)
  )
  expect_identical(restarted$status, "OK")
  expect_match(restarted$finding,
               "started again from zero in 2 of 4,000 searches", fixed = TRUE)
  both <- fit_check_mode_search(
    list(searches = 4000L, at_limit = 3L, restarted = 2L)
  )
  expect_match(both$finding, "limit in 3, and failed from the mode found",
               fixed = TRUE)
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
         validation_budget = 100L, radius_from = "replicate",
         second_validation = list(relative_ess = 0.1, pareto_k = 1.1)),
    list(proposal = "global", certificate = 0.05, estimate = later,
         validation_budget = 1000L, radius_from = "certification",
         second_validation = list(relative_ess = 0.5, pareto_k = 0.3)),
    list(proposal = "global", certificate = NA_real_, estimate = later)
  )
  checked <- check_fit(fit)
  expect_identical(checked$stages$stage, 1:3)
  expect_identical(checked$stages$proposal, c("laplace", "global", "global"))
  expect_identical(checked$stages$certificate, c(0.5, 0.05, NA))
  expect_identical(checked$stages$radius_from,
                   c("replicate", "certification", NA))
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

test_that("printing gives a check that needs a look in full and an OK one on a line", {
  flagged <- with_flags(mcml_boundary = "nugget variance (tau2)")
  checked <- check_fit(flagged)
  withr::local_options(width = 80L)
  printed <- capture.output(result <- print(checked))
  expect_identical(result, checked)
  expect_identical(printed[1L], "Checks on this fit")
  expect_match(printed[2L], "^Fitted in ")
  expect_true(any(grepl("^SUSPICIOUS +Variances +The nugget", printed)))
  expect_true(any(grepl("^ +Fix: Consider refitting without it", printed)))
  # Each OK check is one line with its brief finding, the next line starting
  # another check or the closing lines.
  checks <- checked$checks
  ok <- checks[checks$status == "OK", ]
  expect_gte(nrow(ok), 3L)
  for (row in seq_len(nrow(ok))) {
    line <- which(startsWith(printed, paste0(formatC("OK", width = -10L), "  ",
                                             ok$check[[row]], " ")))
    expect_length(line, 1L)
    expect_match(printed[[line]], paste0("  ", ok$brief[[row]]), fixed = TRUE)
    expect_false(startsWith(printed[[line + 1L]], " "))
  }
  expect_true(any(grepl("^OK +Monte Carlo accuracy +Within [0-9.]+ standard errors\\.$",
                        printed)))
  expect_false(any(grepl("though it cannot rule one out", printed)))
  expect_identical(
    utils::tail(printed, 2L),
    c("summary() of these checks gives the evidence behind each finding.",
      "These checks concern the computation, not whether the model suits the data.")
  )
  expect_false(any(grepl("^ +$", printed)))
  withr::local_options(width = 60L)
  narrow <- capture.output(print(checked))
  expect_lte(max(nchar(narrow)), 80L)
  # A table made before the brief findings gives the findings.
  old <- checked
  old$checks$brief <- NULL
  expect_true(any(grepl("though it cannot",
                        capture.output(print(old)))))
})

test_that("each OK check has a brief finding for its line", {
  expect_identical(fit_check_monte_carlo(certified_mcml)$brief,
                   "Within 0.059 standard errors.")
  expect_identical(
    fit_check_weights(list(effective_draws = 1081.4, pareto_k = 0.46))$brief,
    "1,081 effective draws, tail shape 0.46."
  )
  expect_identical(fit_check_maxima(list())$brief,
                   "The check found no higher maximum.")
  expect_identical(fit_check_variances(list())$brief,
                   "No variance was estimated as zero.")
  falling <- fit_check_variances(list(
    mcml_boundary = c("nugget variance (tau2)", "random-effect variance (sigma2_re.i)"),
    mcml_held_variances = data.frame(
      name = c("nugget variance (tau2)", "random-effect variance (sigma2_re.i)"),
      slope = c(-3, -2), se = c(0.5, 0.5), rising = FALSE, falling = TRUE
    )
  ))
  expect_identical(falling$status, "OK")
  expect_identical(falling$brief, "Zero is the maximum for tau2 and sigma2_re.i.")
  expect_identical(
    fit_check_spatial_scale(binomial_model, list())$brief,
    paste0(signif(coef(binomial_model)$phi, 3), " ",
           binomial_model$distance_units %||% "distance units",
           ", longer than the median spacing.")
  )
  expect_identical(fit_check_importance_distribution(list())$brief,
                   "The global distribution at every stage.")
  expect_identical(
    fit_check_importance_distribution(list(
      mcml_start_replaced = TRUE, mcml_proposal_fallbacks = c(1L, 3L),
      mcml_proposal_errors = list(count = 2L, first = "not positive definite")
    ))$brief,
    paste0("Restarted from the starting values, the Laplace approximation at ",
           "stages 1, 3, not built at 2 parameter values.")
  )
  expect_identical(
    fit_check_mode_search(list(searches = 4000L, at_limit = 3L,
                               restarted = 2L))$brief,
    "At its limit in 3, restarted in 2 of 4,000 searches."
  )
  expect_identical(fit_check_link(list(name = "custom"), "binomial")$brief,
                   "A function, not known to be log-concave.")
  expect_identical(
    fit_check_link(list(name = "custom", link = "cauchit"), "binomial")$brief,
    "cauchit, not known to be log-concave."
  )
  # The brief of a check that needs a look is its finding, which print()
  # gives in full.
  heavy <- fit_check_weights(list(effective_draws = 94, pareto_k = 1.005))
  expect_identical(heavy$brief, heavy$finding)
})

test_that("a fix names the extra check that would settle it", {
  expect_identical(describe_profile_pointer(),
                   ", check_fit(fit, extra = \"profile\")")
  expect_match(describe_profile_pointer(TRUE), "the check \"Profile likelihood\"",
               fixed = TRUE)
  heavy <- list(effective_draws = 94, pareto_k = 1.005)
  expect_match(fit_check_weights(heavy, radius = describe_radius_pointer())$fix,
               paste0("read the bound with caution. The check of the radius, ",
                      "check_fit(fit, extra = \"radius\"), predicts whether it ",
                      "is too small."), fixed = TRUE)
  expect_no_match(fit_check_weights(heavy)$fix, "extra")

  # In check_fit(), the fixes follow the extras that ran.
  weak <- with_flags(mcml_small_scale = TRUE, mcml_ambiguous_maxima = TRUE)
  fixes <- function(checked) {
    stats::setNames(checked$checks$fix, checked$checks$check)
  }
  plain <- fixes(check_fit(weak))
  expect_match(plain[["Spatial scale"]], "check_fit(fit, extra = \"profile\")",
               fixed = TRUE)
  expect_match(plain[["Other maxima"]], "check_fit(fit, extra = \"profile\")",
               fixed = TRUE)
  scale_profile <- profile(weak, parameters = "phi", n_points = 5,
                           n_draws = 200, messages = FALSE)
  profiled <- fixes(check_fit(weak, profile = scale_profile))
  for (check in c("Spatial scale", "Other maxima")) {
    expect_no_match(profiled[[check]], "extra = ", fixed = TRUE)
    expect_match(profiled[[check]],
                 "in the check \"Profile likelihood\" and its panels in plot()",
                 fixed = TRUE)
  }
  # The radius is named for a fit with a certificate, and not for a Gaussian
  # fit, where the check does not apply.
  results <- run_fit_checks(
    binomial_model,
    modifyList(summarise_mcml_fit(binomial_model),
               list(certification_bank = heavy))
  )
  expect_match(results[["Importance weights"]]$row$fix,
               "check_fit(fit, extra = \"radius\")", fixed = TRUE)
  ran <- run_fit_checks(
    binomial_model,
    modifyList(summarise_mcml_fit(binomial_model),
               list(certification_bank = heavy)),
    extras = "radius"
  )
  expect_match(ran[["Importance weights"]]$row$fix,
               "See the check \"Radius of the bound\".", fixed = TRUE)
})

test_that("the plot shows the stages, and a profile once it is computed", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  checked <- check_fit(binomial_model)
  panels <- plot(checked)
  expect_identical(names(panels), c("certificate", "effective_draws", "pareto_k"))
  for (panel in panels) expect_s3_class(panel, "ggplot")

  # A panel that was not computed says which extra to run, and plot() draws
  # the others, or stops when there are none.
  expect_message(chosen <- plot(checked, which = c("certificate", "profile")),
                 "needs the extra check \"profile\", check_fit(fit, extra = \"profile\")",
                 fixed = TRUE)
  expect_identical(names(chosen), "certificate")
  expect_error(plot(checked, which = "profile"),
               "needs the extra check \"profile\"")
  expect_message(plot(checked, which = c("certificate", "moves")),
                 "no changes in the estimates to plot")
  expect_error(plot(checked, which = "radius"),
               "needs the extra check \"radius\", check_fit(fit, extra = \"radius\")",
               fixed = TRUE)
  expect_error(plot(checked, which = "radii"),
               "Unknown panel \"radii\". The panels are \"certificate\"")
  expect_error(plot(checked, which = "radii"),
               "\"held_variances\", \"screen\", \"radius\", \"profile_\"",
               fixed = TRUE)
  expect_error(plot(checked, which = 1), "'which' must be a character vector")

  scale_profile <- profile(binomial_model, parameters = "phi", n_points = 5,
                           n_draws = 200, messages = FALSE)
  given <- check_fit(binomial_model, profile = scale_profile)
  with_profile <- plot(given)
  expect_identical(names(with_profile), c(names(panels), "profile_phi"))
  expect_s3_class(with_profile$profile_phi, "ggplot")
  expect_identical(names(plot(given, which = "profile")), "profile_phi")
  expect_identical(names(plot(given, which = c("profile_phi", "pareto_k"))),
                   c("profile_phi", "pareto_k"))
  expect_message(plot(given, which = c("profile_sigma2", "profile_phi")),
                 "\"profile_sigma2\" has no profile, since the profiles are of phi",
                 fixed = TRUE)

  stages <- data.frame(stage = 1:2, certificate = c(NA, 0.05))
  expect_null(plot_stage_panel(stages[1L, ], "certificate", "t", "y"))
  expect_s3_class(plot_stage_panel(stages, "certificate", "t", "y", limit = 0.1,
                                   log_scale = TRUE), "ggplot")
  expect_null(plot_estimate_moves(NULL))
  moves <- data.frame(parameter = c("a", "b", "a", "b"), stage = c(2L, 2L, 3L, 3L),
                      move = c(1, 2, 0.1, 0.2))
  expect_s3_class(plot_estimate_moves(moves), "ggplot")
})

test_that("the panels of a held variance and of the screen say why they are empty", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  # The binomial fixture holds no variance, and its random effects keep the
  # screen of short scales from running.
  checked <- check_fit(binomial_model)
  expect_error(plot(checked, which = "held_variances"),
               "No variance was held at zero, so there is no slope at zero to plot.",
               fixed = TRUE)
  expect_error(plot(checked, which = "screen"),
               paste0("The screen of short scales did not run on this fit, ",
                      "because of its unstructured random effects."),
               fixed = TRUE)
  expect_message(chosen <- plot(checked, which = c("certificate", "screen")),
                 "unstructured random effects", fixed = TRUE)
  expect_identical(names(chosen), "certificate")

  spatial <- check_fit(with_flags(mcml_boundary = "spatial variance (sigma2)"))
  expect_error(plot(spatial, which = "held_variances"),
               paste0("The spatial variance (sigma2) held at zero was not ",
                      "tested, since with the spatial variance at zero the ",
                      "slope depends on the scale"), fixed = TRUE)
  untested <- check_fit(with_flags(mcml_boundary = "nugget variance (tau2)"))
  expect_error(plot(untested, which = "held_variances"),
               paste0("The fit recorded no test of the slope at zero for the ",
                      "nugget variance (tau2)."), fixed = TRUE)
  failed <- check_fit(with_flags(mcml_maxima_screen = list(gate = "failed")))
  expect_error(plot(failed, which = "screen"),
               "likelihood with independent location effects could not be computed",
               fixed = TRUE)
  unrecorded <- check_fit(with_flags(mcml_maxima_screen = NULL))
  expect_error(plot(unrecorded, which = "screen"),
               "The fit recorded no screen of short scales.", fixed = TRUE)

  gaussian <- check_fit(gaussian_model)
  expect_error(plot(gaussian, which = "screen"),
               "A Gaussian fit is fitted exactly, so it has no screen of short scales to plot.",
               fixed = TRUE)
  expect_error(plot(gaussian, which = "held_variances"),
               "so it has no test of a variance held at zero to plot.", fixed = TRUE)
  no_maximum <- check_fit(with_flags(
    mcml_no_maximum = TRUE, mcml_outcome = "no maximum",
    mcml_separation = list(everywhere = "no positives")
  ))
  expect_error(plot(no_maximum, which = c("held_variances", "screen")),
               "The likelihood has no maximum, so there is no test of a held variance to plot.",
               fixed = TRUE)
})

test_that("the panel of the held variances draws the test the fit recorded", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  slopes <- data.frame(
    parameter = c(5L, 6L, 7L), slope = c(0.5, 3, -2), se = c(0.4, 0.5, 0.5),
    rising = c(FALSE, TRUE, FALSE), falling = c(FALSE, FALSE, TRUE),
    name = c("nugget variance (tau2)", "random-effect variance (sigma2_re.a)",
             "random-effect variance (sigma2_re.b)")
  )
  held <- with_flags(mcml_boundary = slopes$name, mcml_held_variances = slopes)
  checked <- check_fit(held)
  expect_identical(
    checked$results$Variances$data[["Slope of the likelihood at zero"]],
    slopes[c("name", "slope", "se", "rising", "falling")]
  )
  panels <- plot(checked, which = "held_variances")
  expect_named(panels, "held_variances")
  expect_s3_class(panels$held_variances, "ggplot")
  # The interval is that of the test of the fit, 2.58 standard errors each
  # side, so that a slope is rising or falling when its interval excludes zero.
  expect_identical(formals(assess_held_variances)$multiplier, 2.58)
  ranges <- ggplot2::layer_data(panels$held_variances, 2L)
  expect_equal(sort(ranges$xmin), sort(slopes$slope - 2.58 * slopes$se),
               tolerance = 1e-12)
  expect_equal(sort(ranges$xmax), sort(slopes$slope + 2.58 * slopes$se),
               tolerance = 1e-12)
  matched <- match(round(ranges$xmin, 10),
                   round(slopes$slope - 2.58 * slopes$se, 10))
  expect_identical(ranges$xmin > 0, slopes$rising[matched])
  expect_identical(ranges$xmax < 0, slopes$falling[matched])
  labels <- ggplot2::layer_data(panels$held_variances, 4L)$label
  expect_setequal(labels, c("undecided", "rises from zero", "falls from zero"))
  # The default plot includes it.
  expect_contains(names(plot(checked)), "held_variances")
  expect_null(plot_held_variance_panel(NULL))
  expect_null(plot_held_variance_panel(slopes[0L, ]))
})

test_that("the screen of short scales is tabulated from the fit", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  screened <- suppressWarnings(glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
    control_mcml = fixture_control_mcml, return_samples = FALSE,
    messages = FALSE
  ))
  record <- attr(screened, "optimizer")$mcml_maxima_screen
  checked <- check_fit(screened)
  data <- checked$results[["Other maxima"]]$data
  units <- screened$distance_units %||% "distance units"
  expect_identical(data[["Largest exact lower bound, above the estimate"]],
                   record$alarm$statistic - record$alarm$fitted_comparison)
  expect_identical(
    data[["Limit with independent location effects, above the estimate"]],
    record$independent$log_likelihood - record$alarm$fitted_comparison
  )
  expect_identical(data[[paste0("Scale of the estimate (", units, ")")]],
                   coef(screened)$phi)
  # The estimate of this fit is at a short scale, so the screen did not run.
  expect_identical(record$gate, "not run")
  expect_null(data[["Short scales, log-likelihood above the estimate"]])
  expect_error(plot(checked, which = "screen"),
               paste0("The screen of short scales did not run on this fit, ",
                      "because its estimate is at a short scale or its ",
                      "spatial variance is at zero."), fixed = TRUE)
})

test_that("the screen table keeps what the fit recorded", {
  screen <- list(
    gate = "closed",
    screen = data.frame(log_phi = log(c(1, 2)), gap = c(-3, -0.5)),
    independent = list(log_likelihood = -105),
    alarm = list(statistic = -104, fitted_comparison = -100)
  )
  # The scales are exponentiated from the log scale of the record.
  expect_equal(tabulate_short_scales(screen, "km"),
               data.frame(scale_km = c(1, 2), screened = c(-3, -0.5)),
               tolerance = 1e-12)
  expect_null(tabulate_short_scales(list(gate = "not run"), "km"))
  checked <- check_fit(with_flags(mcml_maxima_screen = screen))
  data <- checked$results[["Other maxima"]]$data
  expect_identical(data[["Largest exact lower bound, above the estimate"]], -4)
  expect_identical(
    data[["Limit with independent location effects, above the estimate"]], -5
  )
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  panel <- plot(checked, which = "screen")$screen
  expect_s3_class(panel, "ggplot")
  built <- ggplot2::ggplot_build(panel)
  labels <- unlist(lapply(built$data, function(layer) as.character(layer$label)))
  expect_contains(labels, c("gate", "independent location effects"))
  expect_identical(
    gsub("\n", " ", panel$labels$subtitle),
    paste0("No short scale came within 0.25 of the estimate on the fitting ",
           "draws, so the fit was not repeated. No exact lower bound is ",
           "clearly above the estimate. The largest exact lower bound is 4 ",
           "below the estimate.")
  )
  expect_lte(max(nchar(strsplit(panel$labels$subtitle, "\n")[[1L]])), 62L)
  expect_contains(capture.output(print(summary(checked))),
                  "  Short scales, log-likelihood above the estimate")
  expect_match(describe_screen_outcome(list("Screen of short scales" = "opened",
                                            "Exact alarm" = TRUE), 0.25),
               "fit was repeated from it. An exact lower bound is above")
  expect_null(plot_screen_panel(list()))
})
