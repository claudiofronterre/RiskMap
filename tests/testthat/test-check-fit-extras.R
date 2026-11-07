# Marginal profiles built by hand, with the columns and elements that the
# profile check reads. Each is given by `profile_part()` with the values of
# its parameter on the natural scale, and `higher` is how far the highest
# point is above the maximum on the bank at the estimate.
profile_part <- function(value, loglik, lower, upper) {
  list(value = value, loglik = loglik, lower = lower, upper = upper)
}
make_profiles <- function(..., higher = 0) {
  parts <- list(...)
  top <- max(unlist(lapply(parts, `[[`, "loglik")))
  marginal <- lapply(names(parts), function(name) {
    part <- parts[[name]]
    data.frame(parameter = name, value = log(part$value), loglik = part$loglik,
               status = "ok", deviance = 2 * (top - part$loglik))
  })
  names(marginal) <- names(parts)
  ends <- function(end) vapply(parts, function(part) as.numeric(part[[end]]), 0)
  structure(list(
    marginal = marginal,
    intervals = data.frame(parameter = names(parts), estimate = 2,
                           profile_lower = ends("lower"),
                           profile_upper = ends("upper")),
    level = 0.95, max_loglik = top, higher_than_estimate = higher
  ), class = "RiskMap_profile")
}

no_maximum_fit <- function() {
  optimizer <- attr(binomial_model, "optimizer")
  optimizer[c("mcml_no_maximum", "mcml_outcome", "mcml_separation")] <- list(
    TRUE, "no maximum", list(everywhere = "no positives")
  )
  fit <- binomial_model
  attr(fit, "optimizer") <- optimizer
  fit
}

test_that("the extra checks are named, with all for every one", {
  expect_identical(match_extra_checks(character()), character())
  expect_identical(match_extra_checks("profile"), "profile")
  expect_identical(match_extra_checks(c("profile", "profile")), "profile")
  expect_identical(match_extra_checks("all"), names(list_extra_checks()))
  expect_error(check_fit(binomial_model, extra = "radius"),
               "Unknown extra check \"radius\". The extra checks are \"profile\", or \"all\"",
               fixed = TRUE)
  expect_error(check_fit(binomial_model, extra = c("profile", "a", "b")),
               "Unknown extra checks \"a\", \"b\"", fixed = TRUE)
  expect_error(check_fit(binomial_model, extra = TRUE),
               "'extra' must be a character vector")
  expect_error(check_fit(binomial_model, extra = NA_character_),
               "'extra' must be a character vector")
  expect_error(check_fit(binomial_model, messages = NA), "'messages' must be")
})

test_that("a fit without extras keeps every default check and no extra", {
  checked <- check_fit(binomial_model)
  expect_named(checked, c("checks", "draws", "stages", "moves", "tolerance",
                          "results", "extras", "not_computed", "fit_seconds"))
  expect_identical(checked$extras, list())
  expect_identical(checked$not_computed, character())
  expect_identical(checked$fit_seconds,
                   attr(binomial_model, "optimizer")$elapsed_seconds)
  expect_setequal(names(checked$results), checked$checks$check)
  for (name in names(checked$results)) {
    result <- checked$results[[name]]
    expect_named(result, c("row", "data"))
    expect_identical(result$row$check, name)
    row <- checked$checks[checked$checks$check == name, ]
    rownames(row) <- NULL
    expect_identical(result$row, row)
  }
  expect_identical(checked$checks,
                   list_fit_checks(binomial_model, summarise_mcml_fit(binomial_model)))
  expect_identical(names(check_fit(gaussian_model)$results), "Monte Carlo accuracy")
  expect_null(check_fit(gaussian_model)$results[[1L]]$data)
})

test_that("the print of a fit without extras is unchanged", {
  # The print reads only what check_fit() gave before the extras, as on dev.
  checked <- check_fit(binomial_model)
  before <- structure(unclass(checked)[c("checks", "draws", "stages", "moves",
                                         "tolerance")],
                      class = "RiskMap_check")
  expect_identical(capture.output(print(checked)), capture.output(print(before)))
  # The layout, on rows that do not depend on a fit.
  rows <- list(
    fit_check_monte_carlo(list(certificate = 0.0594, tolerance = 0.1,
                               outcome = "certified", correction_refused = FALSE)),
    fit_check_weights(list(effective_draws = 1081.4, pareto_k = 0.46)),
    fit_check_maxima(list()),
    fit_check_variances(list(mcml_boundary = "nugget variance (tau2)")),
    fit_check_importance_distribution(list(mcml_proposal = "global"))
  )
  fixed <- structure(list(
    checks = tabulate_fit_checks(lapply(rows, check_result)),
    draws = describe_mcml_draws(list(stages = 2L, fitting_draws = 2000L,
                                     validation_draws = 5000L))
  ), class = "RiskMap_check")
  withr::local_options(width = 80L)
  expect_snapshot(print(fixed))
})

test_that("the covariance parameters are those the fit estimated", {
  expect_identical(list_covariance_parameters(binomial_model),
                   c("sigma2", "phi", "sigma2_re.i"))
  expect_identical(list_covariance_parameters(gaussian_model),
                   c("sigma2", "phi", "sigma2_me", "sigma2_re.i"))
  # A Gaussian nugget is profiled as its variance, and a fixed measurement
  # error variance is not profiled.
  expect_identical(list_covariance_parameters(gaussian_offset_model),
                   c("sigma2", "phi", "tau2"))
  expect_identical(list_covariance_parameters(gaussian_intercept_model),
                   c("sigma2", "phi", "sigma2_me"))
  coefficients_only <- binomial_model
  coefficients_only$estimate <- binomial_model$estimate["beta"]
  expect_identical(list_covariance_parameters(coefficients_only), character())
  expect_identical(find_extra_check_obstacle("profile", coefficients_only),
                   "the model has no covariance parameters")
  expect_null(find_extra_check_obstacle("profile", binomial_model))
})

test_that("the profile check runs on request, from a stream of the fit's seed", {
  withr::local_preserve_seed()
  set.seed(11)
  before <- .Random.seed
  expect_message(
    checked <- check_fit(binomial_model, extra = "profile"),
    paste0("Computing the profiles of the 3 covariance parameters, which should ",
           "take about .*, from the .* the fit took")
  )
  expect_identical(.Random.seed, before)
  expect_named(checked$extras, "profile")
  expect_s3_class(checked$extras$profile$data, "RiskMap_profile")
  expect_named(checked$extras$profile$data$marginal,
               c("sigma2", "phi", "sigma2_re.i"))
  expect_identical(checked$extras$profile$data$settings$seed,
                   make_extra_check_seed(binomial_model$seed, 9600L))
  expect_contains(checked$checks$check, "Profile likelihood")
  expect_identical(nrow(checked$checks), nrow(check_fit(binomial_model)$checks) + 1L)
  # Reproducible, with or without a random number state in the session.
  rm(".Random.seed", envir = .GlobalEnv)
  expect_silent(again <- check_fit(binomial_model, extra = "all", messages = FALSE))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(again$extras, checked$extras)
  expect_identical(again$checks, checked$checks)
  direct <- profile(binomial_model, parameters = c("sigma2", "phi", "sigma2_re.i"),
                    seed = make_extra_check_seed(binomial_model$seed, 9600L),
                    messages = FALSE)
  expect_identical(checked$extras$profile$data, direct)
  printed <- capture.output(print(checked))
  expect_true(any(grepl("^SUSPICIOUS +Profile likelihood +The 95% profile", printed)))
})

test_that("a profile computed by the user is used as it is", {
  scale_profile <- profile(binomial_model, parameters = "phi", n_points = 5,
                           n_draws = 200, messages = FALSE)
  expect_silent(checked <- check_fit(binomial_model, profile = scale_profile))
  expect_identical(checked$extras$profile$data, scale_profile)
  expect_identical(
    check_fit(binomial_model, extra = "profile", profile = scale_profile)$extras,
    checked$extras
  )
  expect_error(check_fit(binomial_model, profile = 1), "'profile' must be a profile")
  # Whatever parameters it covers, coefficients included.
  other_profile <- profile(binomial_model, parameters = c("cov", "sigma2"),
                           n_points = 5, n_draws = 200, messages = FALSE)
  other <- check_fit(binomial_model, profile = other_profile)
  expect_identical(other$extras$profile$data, other_profile)
  expect_named(other$extras$profile$data$marginal, c("beta.cov", "sigma2"))
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  expect_named(plot(other, which = "profile"),
               c("profile_beta.cov", "profile_sigma2"))
  expect_error(check_fit(poisson_model, profile = scale_profile),
               "'profile' must be a profile of this fit")
  expect_error(check_fit(binomial_model, cores = 0),
               "'cores' must be a positive whole number")
  expect_error(check_fit(binomial_model, cores = 1.5),
               "'cores' must be a positive whole number")
})

test_that("an extra check that does not apply or fails says why", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  expect_silent(checked <- check_fit(no_maximum_fit(), extra = "profile"))
  expect_identical(checked$extras, list())
  reason <- paste0("The profiles of the covariance parameters were not computed, ",
                   "because the likelihood has no maximum.")
  expect_identical(checked$not_computed, c(profile = reason))
  expect_identical(checked$checks$check, "Maximum of the likelihood")
  printed <- capture.output(print(checked))
  expect_match(paste(printed, collapse = " "), reason, fixed = TRUE)
  expect_message(plot(checked, which = c("certificate", "profile")), reason,
                 fixed = TRUE)
  expect_error(plot(checked, which = "profile"), reason, fixed = TRUE)

  # The number of cores goes to profile().
  cores <- NULL
  local_mocked_bindings(profile.RiskMap = function(fitted, ..., cores) {
    cores <<- cores
    stop("no bank")
  })
  failed <- check_fit(binomial_model, extra = "profile", messages = FALSE,
                      cores = 3)
  expect_identical(cores, 3L)
  expect_identical(failed$extras, list())
  expect_identical(
    failed$not_computed[["profile"]],
    paste0("The profiles of the covariance parameters were not computed, ",
           "because their computation stopped with the error \"no bank\".")
  )
  expect_identical(failed$checks, check_fit(binomial_model)$checks)
})

test_that("the profiles of a fit with a named link use that link", {
  probit <- glgpm(y ~ cov + gp(), data = binomial_data,
                  family = binomial("probit"), denominator = denominator,
                  control_mcmc = control_mcmc, seed = fixture_seed,
                  control_mcml = fixture_control_mcml, messages = FALSE)
  captured <- NULL
  local_mocked_bindings(
    make_adaptive_mcml_engine = function(..., invlink = list(missing = TRUE)) {
      captured <<- invlink
      stop("engine captured")
    }
  )
  checked <- check_fit(probit, extra = "profile", messages = FALSE)
  expect_match(checked$not_computed[["profile"]], "engine captured")
  expect_identical(captured$link, "probit")
  expect_identical(captured$inv, probit$link_function$inv)
  expect_identical(captured$d1, probit$link_function$d1)
})

test_that("a Gaussian profile needs no seed and can be plotted alone", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  expect_message(checked <- check_fit(gaussian_model, extra = "profile"),
                 "^Computing the profiles of the 4 covariance parameters.\n$")
  expect_null(make_extra_check_seed(gaussian_model$seed, 9600L))
  expect_identical(checked$checks$check,
                   c("Profile likelihood", "Monte Carlo accuracy"))
  expect_identical(names(plot(checked)),
                   paste0("profile_", c("sigma2", "phi", "sigma2_me", "sigma2_re.i")))
  expect_identical(names(plot(checked, which = "profile_sigma2_me")),
                   "profile_sigma2_me")
  expect_error(plot(checked, which = "moves"), "no stages to plot")
})

test_that("the seed of an extra check wraps like the streams of a fit", {
  expect_identical(make_extra_check_seed(12345L, 9600L), 21945L)
  expect_identical(make_extra_check_seed(.Machine$integer.max, 9600L), 9600L)
  expect_null(make_extra_check_seed(NULL, 9600L))
})

test_that("the time of the extra checks comes from the time of the fit", {
  fit <- binomial_model
  optimizer <- attr(fit, "optimizer")
  optimizer$elapsed_seconds <- 12
  attr(fit, "optimizer") <- optimizer
  # Ten times the fit for each of the three parameters, less on several
  # cores, of which at most two for each parameter help.
  expect_identical(
    describe_extra_check_time(fit, "profile"),
    paste0("Computing the profiles of the 3 covariance parameters, which ",
           "should take about 6 minutes, from the 12 seconds the fit took.")
  )
  skip_on_os("windows")
  expect_match(describe_extra_check_time(fit, "profile", cores = 2L),
               "about 3 minutes", fixed = TRUE)
  expect_match(describe_extra_check_time(fit, "profile", cores = 16L),
               "about 60 seconds", fixed = TRUE)
  optimizer$elapsed_seconds <- 0.01
  attr(fit, "optimizer") <- optimizer
  expect_match(describe_extra_check_time(fit, "profile"),
               "about 1 second, from the 1 second the fit took.", fixed = TRUE)
  expect_identical(describe_extra_check_time(gaussian_model, "profile"),
                   "Computing the profiles of the 4 covariance parameters.")
  single <- gaussian_model
  single$estimate <- gaussian_model$estimate[c("beta", "phi")]
  expect_identical(describe_extra_check_time(single, "profile"),
                   "Computing the profile of the 1 covariance parameter.")
})

test_that("separate maxima of a profile need a fall between them", {
  expect_identical(find_profile_maxima(numeric()), integer())
  expect_identical(find_profile_maxima(c(-3, -1, 0, -1, -3)), 3L)
  expect_identical(find_profile_maxima(c(-3, -1, 0, -1, -0.2, -2)), c(3L, 5L))
  # A dip of less than half a unit does not separate two maxima.
  expect_identical(find_profile_maxima(c(-3, 0, -0.4, -0.1, -2)), 2L)
  # A profile rising to the end of its grid has a maximum there.
  expect_identical(find_profile_maxima(c(0, -1, -2, -1.2)), c(1L, 4L))
  expect_identical(find_profile_maxima(rep(0, 4)), 1L)
})

test_that("a parameter of a profile is named in words", {
  expect_identical(describe_profile_parameter("phi"), "the spatial scale (phi)")
  expect_identical(describe_profile_parameter("sigma2"),
                   "the spatial variance (sigma2)")
  expect_identical(describe_profile_parameter("nu2"),
                   "the ratio of the nugget to the spatial variance (nu2)")
  expect_identical(describe_profile_parameter("tau2"), "the nugget variance (tau2)")
  expect_identical(describe_profile_parameter("sigma2_me"),
                   "the measurement error variance (sigma2_me)")
  expect_identical(describe_profile_parameter("sigma2_re.village"),
                   "the random-effect variance for 'village' (sigma2_re.village)")
  expect_identical(describe_profile_parameter("beta.cov"),
                   "the coefficient of 'cov' (beta.cov)")
  expect_identical(join_phrases(character()), "")
  expect_identical(join_phrases("a"), "a")
  expect_identical(join_phrases(c("a", "b", "c")), "a, b and c")
})

test_that("the profile row names each parameter that needs a look", {
  grid <- c(0.5, 1, 2, 4, 8, 16, 32)
  peaked <- c(-6, -2, 0, -2, -6, -9, -12)
  one <- fit_check_profile(
    make_profiles(phi = profile_part(grid, peaked, 1.2, 3.4),
                  sigma2 = profile_part(grid, peaked, 1.1, 3.9)), "km"
  )
  expect_identical(one$status, "OK")
  expect_identical(one$fix, "")
  expect_identical(one$finding, paste0(
    "The profiles of the spatial scale (phi) and the spatial variance (sigma2) ",
    "each have one maximum and a closed 95% interval."
  ))
  expect_identical(
    fit_check_profile(make_profiles(phi = profile_part(grid, peaked, 1, 3)))$finding,
    "The profile of the spatial scale (phi) has one maximum and a closed 95% interval."
  )

  second <- fit_check_profile(make_profiles(
    phi = profile_part(grid, c(-6, -2, 0, -2, -1, -5, -9), 1.2, 3.4),
    sigma2 = profile_part(grid, peaked, 1.1, 3.9)
  ))
  expect_identical(second$status, "SUSPICIOUS")
  expect_identical(second$finding, paste0(
    "The profile of the spatial scale (phi) has a second maximum at 8 distance ",
    "units, 1 log-likelihood units below the highest. Wald intervals may ",
    "describe the uncertainty poorly."
  ))
  expect_identical(second$fix, paste0(
    "Use the profile interval of phi rather than the Wald interval, and see ",
    "its panel in plot()."
  ))

  flat <- fit_check_profile(make_profiles(
    phi = profile_part(grid, c(0, 0, 0, 0, 0, 0, 0), NA, NA),
    sigma2 = profile_part(grid, c(0, 0, 0, -0.1, -2, -6, -9), NA, 3.4),
    nu2 = profile_part(grid, c(0, 0, 0, -0.5, -3, -6, -9), NA, 2.1),
    sigma2_re.i = profile_part(grid, peaked, 1, 3)
  ), "km", advice = " Or the penalty.")
  expect_identical(flat$status, "SUSPICIOUS")
  expect_identical(flat$finding, paste0(
    "The 95% profile interval is open at both ends for the spatial scale (phi), ",
    "and below for the spatial variance (sigma2) and the ratio of the nugget to ",
    "the spatial variance (nu2), so the profile is flat there and the data do ",
    "not bound these parameters on that side. Wald intervals may describe the ",
    "uncertainty poorly."
  ))
  expect_identical(flat$fix, paste0(
    "Use the profile intervals of phi, sigma2 and nu2 rather than the Wald ",
    "intervals, and see their panels in plot(). Consider a scale fixed from ",
    "external knowledge. Or the penalty. A variance whose interval is open ",
    "below may be zero, so consider refitting without nu2."
  ))
  above <- fit_check_profile(make_profiles(
    sigma2 = profile_part(grid, -(0:6), 0.5, NA)
  ))
  expect_match(above$finding, paste0(
    "open above for the spatial variance (sigma2), so the profile is flat there ",
    "and the data do not bound the parameter on that side."
  ), fixed = TRUE)

  both <- fit_check_profile(make_profiles(
    phi = profile_part(grid, c(-1, 0, -1, -0.2, 0, 0, 0), NA, NA),
    sigma2_re.i = profile_part(grid, c(-2, 0, -2, -0.5, -4, -6, -9), NA, 3)
  ), "m")
  expect_match(both$finding, paste0(
    "^The profile of the spatial scale \\(phi\\) has a second maximum at 8 m, ",
    "as high as the highest. The profile of the random-effect variance for 'i' ",
    "\\(sigma2_re.i\\) has a second maximum at 4, 0.5 log-likelihood units ",
    "below the highest. The ",
    "95% profile interval is open at both ends for the spatial scale \\(phi\\), ",
    "and below for the random-effect variance for 'i' \\(sigma2_re.i\\), so"
  ))
  expect_match(both$fix, "intervals of phi and sigma2_re.i rather", fixed = TRUE)

  higher <- fit_check_profile(make_profiles(
    phi = profile_part(grid, c(-6, -2, 0, -2, 0.4, -5, -9), 1.2, 3.4),
    sigma2 = profile_part(grid, peaked, NA, 3.9), higher = 0.4
  ), "km")
  expect_identical(higher$status, "PROBLEM")
  expect_identical(higher$finding, paste0(
    "The profile of the spatial scale (phi) rises above the likelihood of the ",
    "estimate at 8 km, by 0.4 log-likelihood units, so the estimate is not the ",
    "maximum likelihood estimate."
  ))
  expect_identical(higher$fix, paste0(
    "Refit with starting values at the highest point of that profile, through ",
    "'start_pars' in glgpm(), and keep the fit with the higher log-likelihood."
  ))
  two_higher <- fit_check_profile(make_profiles(
    phi = profile_part(grid, c(-6, -2, 0, -2, 0.4, -5, -9), 1.2, 3.4),
    sigma2 = profile_part(grid, c(-6, -2, 0, 0.2, -6, -9, -12), 1.1, 3.9),
    higher = 0.4
  ), "km")
  expect_identical(two_higher$finding, paste0(
    "The profile of the spatial scale (phi) rises above the likelihood of the ",
    "estimate at 8 km, by 0.4 log-likelihood units. The profile of the spatial ",
    "variance (sigma2) rises above the likelihood of the estimate at 4, by 0.2 ",
    "log-likelihood units. The estimate is therefore not the maximum ",
    "likelihood estimate."
  ))
  expect_match(two_higher$fix, "these profiles", fixed = TRUE)
})

test_that("the summary gives the numbers behind each finding", {
  flagged <- binomial_model
  optimizer <- attr(flagged, "optimizer")
  optimizer$mcml_boundary <- "nugget variance (tau2)"
  attr(flagged, "optimizer") <- optimizer
  checked <- check_fit(flagged)
  summarised <- summary(checked)
  expect_s3_class(summarised, "summary.RiskMap_check")
  expect_identical(unname(vapply(summarised$results,
                                  function(result) result$row$check, "")),
                   checked$checks$check)
  expect_identical(summarised$stages$largest_move, NA_real_)
  expect_error(summary.RiskMap_check(1), "RiskMap_check")

  data <- checked$results[["Monte Carlo accuracy"]]$data
  last <- binomial_model$mcml_history[[1L]]
  expect_identical(data[["Certificate, in standard errors"]],
                   attr(binomial_model, "optimizer")$mcml_certificate)
  expect_identical(data[["95% radius of its Monte Carlo error"]],
                   last$certificate_radius)
  expect_identical(checked$results$Variances$data[["Held at zero"]],
                   "nugget variance (tau2)")
  expect_s3_class(checked$results[["Other maxima"]]$data[[
    "Maxima from the dispersed starts"]], "data.frame")

  printed <- capture.output(print(summarised))
  expect_identical(printed[1L], "Evidence for the checks on this fit")
  expect_match(printed[2L], "^Fitted in 1 stage")
  expect_true(any(grepl("^SUSPICIOUS +Variances$", printed)))
  expect_true(any(grepl("^ +Held at zero +nugget variance \\(tau2\\)$", printed)))
  expect_true(any(grepl("^ +Pareto shape +0\\.", printed)))
  expect_contains(printed, "Certificate and certification draws by stage")
  expect_match(printed[length(printed)], "not whether the model suits the data")
  expect_output(result <- print(summarised))
  expect_identical(result, summarised)

  gaussian <- capture.output(print(summary(check_fit(gaussian_model))))
  expect_identical(gaussian[2:3], c("", "OK  Monte Carlo accuracy"))
  expect_contains(gaussian, "  Nothing recorded beyond the finding.")
})

test_that("the summary gives the full results of the extra checks", {
  scale_profile <- profile(binomial_model, parameters = "phi", n_points = 5,
                           n_draws = 200, messages = FALSE)
  printed <- capture.output(print(summary(check_fit(binomial_model,
                                                    profile = scale_profile))))
  expect_true(any(grepl("^SUSPICIOUS +Profile likelihood$", printed)))
  expect_contains(printed, paste0("  ", capture.output(print(scale_profile))))
  points <- which(printed == "  Points of the profile of the spatial scale (phi)")
  expect_length(points, 1L)
  # A heading and a line for each point.
  expect_gte(length(printed) - points, nrow(scale_profile$marginal$phi) + 1L)

  missing <- capture.output(print(summary(check_fit(no_maximum_fit(),
                                                    extra = "profile"))))
  expect_match(paste(missing, collapse = " "),
               paste0("The profiles of the covariance parameters were not ",
                      "computed, because the likelihood has no maximum."),
               fixed = TRUE)
})

test_that("values of the evidence are written in words where they are logical", {
  expect_identical(format_check_value(TRUE), "yes")
  expect_identical(format_check_value(c(1L, 3L)), "1, 3")
  expect_identical(format_check_value(5000L), "5,000")
  expect_identical(format_check_value(0.0123456), "0.01235")
  expect_identical(format_check_value("global"), "global")
  expect_null(describe_screen_gate(NULL))
  expect_identical(describe_screen_gate(list(gate = "closed")), "closed")
  expect_null(tabulate_maxima_check(NULL, "km"))
  check <- data.frame(start = c("fitted", "short"), log_phi = log(c(1, 2)),
                      log_likelihood = c(-10, -11),
                      independent_log_likelihood = c(NA, -11.5),
                      independent_se = c(NA, 0.1), chosen = c(TRUE, FALSE))
  expect_named(tabulate_maxima_check(check, "km"),
               c("start", "scale_km", "loglik", "independent_loglik",
                 "independent_se", "chosen"))
  check$independent_log_likelihood <- NA_real_
  expect_named(tabulate_maxima_check(check, "distance units"),
               c("start", "scale_distance_units", "loglik", "chosen"))
})
