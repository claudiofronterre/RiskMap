##' The extra checks of `check_fit()`, each with what it computes
##'
##' The names are those that `extra` takes, in the order the checks run. The
##' description completes a sentence such as "Computing ...".
##' @noRd
list_extra_checks <- function() {
  c(profile = "the profiles of the covariance parameters",
    radius = "the predictions of how often the radius of the bound is too small")
}

##' The extra checks asked for, in the order they run
##'
##' `"all"` gives every one. An unknown name stops with the valid ones.
##' @noRd
match_extra_checks <- function(extra) {
  if (!is.character(extra) || anyNA(extra)) {
    abort_input("{.arg extra} must be a character vector without missing values.", arg = "extra")
  }
  known <- names(list_extra_checks())
  unknown <- setdiff(extra, c(known, "all"))
  if (length(unknown)) {
    abort_input(c("{.arg extra} has unknown extra check{?s} {.val {unknown}}.", "i" = "The extra checks are {.val {known}}, or {.val all} for every one."), arg = "extra")
  }
  if ("all" %in% extra) known else known[known %in% extra]
}

##' Stop unless `profile` is NULL or a profile of this fit
##'
##' A profile records the estimate it was centred on, which must be that of
##' the fit for the parameters the two share. The nugget of a Gaussian
##' profile is its variance, `tau2`, where the fit has its ratio, `nu2`.
##' @noRd
check_given_profile <- function(profile, object) {
  if (is.null(profile)) return(invisible(TRUE))
  wrong_profile <- function() {
    abort_input("{.arg profile} must be a profile of this fit with at least one marginal profile.", arg = "profile")
  }
  if (!inherits(profile, "RiskMap_profile")) wrong_profile()
  estimate <- unlist(object$estimate)
  shared <- intersect(names(profile$estimate), names(estimate))
  if (!length(profile$marginal) || !length(shared) ||
      !isTRUE(all.equal(unname(profile$estimate[shared]),
                        unname(estimate[shared])))) {
    wrong_profile()
  }
  invisible(TRUE)
}

##' Run the extra checks named in `extra`
##'
##' Gives the `results` of those that ran, each with its row and data, and
##' `not_computed`, a sentence for each of the others saying why. A profile
##' given by the user is used as it is. `cores` goes to the checks that can
##' use several processes. The random number state of the session is
##' restored on exit.
##' @noRd
run_extra_checks <- function(object, extra, profile = NULL, messages = TRUE,
                             cores = 1L) {
  results <- list()
  not_computed <- character()
  if (!length(extra)) {
    return(list(results = results, not_computed = not_computed))
  }
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  descriptions <- list_extra_checks()
  given <- extra == "profile" & !is.null(profile)
  obstacles <- vapply(extra, function(name) {
    find_extra_check_obstacle(name, object) %||% NA_character_
  }, character(1L))
  obstacles[given] <- NA_character_
  computed <- extra[is.na(obstacles) & !given]
  if (messages && length(computed)) {
    riskmap_inform("{describe_extra_check_time(object, computed, cores)}")
  }
  for (name in extra) {
    result <- if (is.na(obstacles[[name]])) {
      tryCatch(run_extra_check(name, object, profile, cores),
               error = function(error) error)
    }
    reason <- if (!is.na(obstacles[[name]])) {
      obstacles[[name]]
    } else if (inherits(result, "error")) {
      paste0("their computation stopped with the error \"",
             condition_text(result), "\"")
    }
    if (is.null(reason)) {
      results[[name]] <- result
    } else {
      not_computed[[name]] <- paste0(
        toupper(substr(descriptions[[name]], 1L, 1L)),
        substring(descriptions[[name]], 2L), " were not computed, because ",
        reason, "."
      )
    }
  }
  list(results = results, not_computed = not_computed)
}

##' Why an extra check does not apply to a fit, NULL when it does
##'
##' Both checks describe a maximum of the likelihood. The profiles need a
##' covariance parameter to profile, and the check of the radius a certificate
##' (`find_radius_check_obstacle()`).
##' @noRd
find_extra_check_obstacle <- function(name, object) {
  if (isTRUE(attr(object, "optimizer")$mcml_no_maximum)) {
    return("the likelihood has no maximum")
  }
  if (identical(name, "profile") && !length(list_covariance_parameters(object))) {
    return("the model has no covariance parameters")
  }
  if (identical(name, "radius")) return(find_radius_check_obstacle(object))
  NULL
}

##' Run one extra check, giving its row and data
##' @noRd
run_extra_check <- function(name, object, profile = NULL, cores = 1L) {
  switch(name, profile = run_profile_check(object, profile, cores),
         radius = fit_check_radius(object))
}

##' The seed of an extra check, a stream of its own after the seed of the fit
##'
##' Wrapped as `fit_adaptive_mcml()` wraps the streams of a fit, which stay
##' below 9,500 after its seed for fewer than 500 stages and so do not reach
##' the offsets used here, 9500 to 9502 for the radius and 9600 for the
##' profiles. NULL for a fit without a seed, such as a Gaussian fit.
##' @noRd
make_extra_check_seed <- function(fit_seed, offset) {
  if (is.null(fit_seed)) return(NULL)
  as.integer(((as.numeric(fit_seed) - 1 + offset) %% .Machine$integer.max) + 1)
}

##' The covariance parameters of a fit, as `profile.RiskMap()` names them
##'
##' Every estimated parameter except the regression coefficients, so a fit
##' without a nugget or random effects has none to profile. The nugget of a
##' Gaussian fit is profiled as its variance, `tau2`.
##' @noRd
list_covariance_parameters <- function(object) {
  names <- names(unlist(object$estimate))
  names <- names[!startsWith(names, "beta.")]
  if (identical(object$family, "gaussian")) names[names == "nu2"] <- "tau2"
  names
}

##' How long the extra checks should take, from the time the fit took
##'
##' With the default settings of `profile.RiskMap()`, a profile of the scale
##' of a binomial or Poisson fit took between 1.8 and 14 times as long as the
##' fit, with a median of 8.5, on 11 fits of one stage to the liberia,
##' anopheles and loaloa data with and without a nugget, and 4.1 times on a
##' fit of two stages, whose profile costs no more. Profiles of every
##' covariance parameter took 13.4, 9.5 and 11.3 times the fit for each
##' parameter on liberia, liberia with a nugget and anopheles with a nugget,
##' so ten times for each parameter is about right. With several cores the
##' two sides of each profile run at once, so at most twice as many cores as
##' parameters help, and a Windows session runs on one. A Gaussian fit
##' records no time, and its exact profiles are quick. The check of the
##' radius took between 0.3 and 2.1 times as long as the fit on twelve
##' simulated Poisson surveys, with a median of 1.7, so twice the fit, and it
##' runs on one core.
##' @noRd
describe_extra_check_time <- function(object, extras, cores = 1L) {
  descriptions <- list_extra_checks()[extras]
  parameters <- length(list_covariance_parameters(object))
  if ("profile" %in% extras) {
    descriptions[["profile"]] <- paste0(
      "the profile", if (parameters != 1L) "s", " of the ", parameters,
      " covariance parameter", if (parameters != 1L) "s"
    )
  }
  listed <- paste(descriptions, collapse = " and ")
  fit_seconds <- attr(object, "optimizer")$elapsed_seconds
  if (is.null(fit_seconds)) return(paste0("Computing ", listed, "."))
  if (.Platform$OS.type == "windows") cores <- 1L
  multiples <- c(profile = 10 * parameters / min(cores, 2L * parameters),
                 radius = 2)
  paste0("Computing ", listed, ", which should take about ",
         format_check_time(sum(multiples[extras]) * fit_seconds), ", from the ",
         format_check_time(fit_seconds), " the fit took.")
}

##' A duration in seconds as readable text, with one second in the singular
##' @noRd
format_check_time <- function(seconds) {
  if (seconds < 1.5) return("1 second")
  format_profile_time(seconds)
}

##' The profile check, from a given profile or one computed for it
##'
##' The profiles are of every covariance parameter of the fit, with the
##' default settings of `profile.RiskMap()` and its bank drawn from the
##' stream 9600 after the seed of the fit.
##' @noRd
run_profile_check <- function(object, profile = NULL, cores = 1L) {
  if (is.null(profile)) {
    profile <- profile.RiskMap(
      object, parameters = list_covariance_parameters(object),
      seed = make_extra_check_seed(object$seed, 9600L) %||% 1L,
      cores = cores, messages = FALSE
    )
  }
  advice <- if (!fit_is_exact(object)) {
    describe_pc_prior_advice(attr(object, "optimizer"))
  } else {
    ""
  }
  check_result(fit_check_profile(profile, object$distance_units, advice),
               profile)
}

##' A parameter of a profile in words, with its name
##' @noRd
describe_profile_parameter <- function(name) {
  words <- switch(
    name,
    sigma2 = "spatial variance",
    phi = "spatial scale",
    nu2 = "ratio of the nugget to the spatial variance",
    tau2 = "nugget variance",
    sigma2_me = "measurement error variance",
    if (startsWith(name, "sigma2_re.")) {
      paste0("random-effect variance for '", sub("^sigma2_re\\.", "", name), "'")
    } else if (startsWith(name, "beta.")) {
      paste0("coefficient of '", sub("^beta\\.", "", name), "'")
    } else {
      name
    }
  )
  paste0("the ", words, " (", name, ")")
}

##' What each marginal profile shows
##'
##' One row per parameter of `profile`, sorted by its value, with the
##' parameter in words, the `top` of the profile above the maximum on the
##' bank at the estimate, its separate `maxima` (`find_profile_maxima()`) and
##' the ends where its interval is open. Values are on the natural scale,
##' with the distance `units` for the scale.
##' @noRd
read_marginal_profiles <- function(profile, units) {
  bank_maximum <- profile$max_loglik - max(profile$higher_than_estimate, 0)
  lapply(stats::setNames(nm = names(profile$marginal)), function(name) {
    points <- profile$marginal[[name]]
    points <- points[is.finite(points$loglik), ]
    points <- points[order(points$value), ]
    interval <- profile$intervals[profile$intervals$parameter == name, ]
    value_at <- function(index) {
      paste0(signif(profile_natural(name, points$value[[index]]), 3),
             if (identical(name, "phi")) paste0(" ", units))
    }
    maxima <- find_profile_maxima(points$loglik)
    top <- which.max(points$loglik)
    list(
      label = describe_profile_parameter(name),
      above = points$loglik[[top]] - bank_maximum, top_value = value_at(top),
      second_value = if (length(maxima) >= 2L) value_at(maxima[2L]),
      second_below = if (length(maxima) >= 2L) {
        max(points$loglik) - points$loglik[[maxima[2L]]]
      },
      open = c(below = is.na(interval$profile_lower),
               above = is.na(interval$profile_upper))
    )
  })
}

##' Phrases joined by commas and a final "and"
##' @noRd
join_phrases <- function(phrases) {
  if (length(phrases) <= 1L) return(paste(phrases, collapse = ""))
  paste(paste(phrases[-length(phrases)], collapse = ", "),
        "and", phrases[length(phrases)])
}

##' What the profiles of the parameters show
##'
##' A point of a profile above the likelihood of the estimate by more than
##' 0.01, as in `profile.RiskMap()`, shows that the estimate is not the
##' maximum likelihood estimate. A second maximum (`find_profile_maxima()`)
##' or an interval open at an end, where the profile stays within the cut-off
##' to the end of its grid or the bound of the parameter, means that Wald
##' intervals describe the uncertainty poorly. The row names each parameter
##' concerned. `advice` adds the PC penalty to the fix of a flat profile of
##' the spatial variance or scale.
##' @noRd
fit_check_profile <- function(profile, units = NULL, advice = "") {
  name <- "Profile likelihood"
  read <- read_marginal_profiles(profile, units %||% "distance units")
  level <- paste0(round(100 * profile$level), "%")
  higher <- Filter(function(item) item$above > 0.01, read)
  second <- Filter(function(item) !is.null(item$second_value), read)
  flat <- Filter(function(item) any(item$open), read)
  if (length(higher)) {
    rises <- vapply(higher, function(item) {
      paste0("The profile of ", item$label, " rises above the likelihood of ",
             "the estimate at ", item$top_value, ", by ", signif(item$above, 2),
             " log-likelihood units")
    }, "")
    return(fit_check_row(
      name, "PROBLEM",
      if (length(rises) == 1L) {
        paste0(rises, ", so the estimate is not the maximum likelihood estimate.")
      } else {
        paste0(paste0(rises, ".", collapse = " "), " The estimate is therefore ",
               "not the maximum likelihood estimate.")
      },
      paste0("Refit with starting values at the highest point of ",
             if (length(higher) > 1L) "these profiles" else "that profile",
             ", through 'start_pars' in glgpm(), and keep the fit with the ",
             "higher log-likelihood.")
    ))
  }
  if (!length(second) && !length(flat)) {
    labels <- vapply(read, `[[`, "", "label")
    return(fit_check_row(name, "OK", paste0(
      "The profile", if (length(read) > 1L) "s", " of ", join_phrases(labels),
      if (length(read) > 1L) " each have" else " has", " one maximum and a ",
      "closed ", level, " interval."
    ), brief = if (length(read) > 1L) {
      paste0("One maximum each, closed ", level, " intervals.")
    } else {
      paste0("One maximum and a closed ", level, " interval.")
    }))
  }
  fit_check_row(
    name, "SUSPICIOUS",
    paste0(
      paste0(vapply(second, function(item) {
        paste0("The profile of ", item$label, " has a second maximum at ",
               item$second_value, ", ",
               if (item$second_below < 0.01) "as high as the highest" else
                 paste0(signif(item$second_below, 2),
                        " log-likelihood units below the highest"), ". ")
      }, ""), collapse = ""),
      if (length(flat)) paste0(
        "The ", level, " profile interval is open ", describe_open_ends(flat),
        ", so the profile is flat there and the data do not bound ",
        if (length(flat) > 1L) "these parameters" else "the parameter",
        " on that side. "
      ),
      "Wald intervals may describe the uncertainty poorly."
    ),
    describe_profile_fix(names(second), flat, advice)
  )
}

##' Where the intervals of the flat profiles are open, grouped by end
##'
##' For example "at both ends for the spatial scale (phi), and below for the
##' spatial variance (sigma2)".
##' @noRd
describe_open_ends <- function(flat) {
  ends <- vapply(flat, function(item) {
    if (all(item$open)) "at both ends" else
      if (item$open[["below"]]) "below" else "above"
  }, "")
  groups <- vapply(intersect(c("at both ends", "below", "above"), ends),
                   function(end) {
    paste0(end, " for ", join_phrases(vapply(flat[ends == end], `[[`, "",
                                             "label")))
  }, "")
  if (length(groups) == 1L) return(groups)
  paste0(paste(groups[-length(groups)], collapse = ", "), ", and ",
         groups[length(groups)])
}

##' The fix of a profile with a second maximum or an open interval
##'
##' `flat` is the list of profiles with an open end, from
##' `read_marginal_profiles()`.
##' @noRd
describe_profile_fix <- function(second, flat, advice) {
  concerned <- union(second, names(flat))
  spatial <- intersect(names(flat), c("sigma2", "phi"))
  towards_zero <- vapply(flat, function(item) item$open[["below"]], logical(1L))
  variances <- setdiff(names(flat)[towards_zero], c("phi", spatial))
  variances <- variances[!startsWith(variances, "beta.")]
  paste0(
    "Use the profile interval", if (length(concerned) > 1L) "s", " of ",
    join_phrases(concerned), " rather than the Wald interval",
    if (length(concerned) > 1L) "s", ", and see ",
    if (length(concerned) > 1L) "their panels" else "its panel", " in plot().",
    if (length(spatial)) paste0(" Consider a scale fixed from external ",
                                "knowledge.", advice),
    if (length(variances)) paste0(
      " A variance whose interval is open below may be zero, so consider ",
      "refitting without ", join_phrases(variances), "."
    )
  )
}

##' The separate maxima of a profile, highest first
##'
##' `loglik` is ordered by the profiled value. A local maximum, an end of the
##' grid included, counts when the profile falls by more than `drop`
##' log-likelihood units between it and every higher maximum. Half a unit is
##' a change of deviance of one, ten times the largest standard deviation of
##' the deviance between banks of 5,000 draws in the profile study of
##' certified MCML, so that a maximum is not the noise of one bank.
##' @noRd
find_profile_maxima <- function(loglik, drop = 0.5) {
  count <- length(loglik)
  if (!count) return(integer())
  before <- c(-Inf, loglik[-count])
  after <- c(loglik[-1L], -Inf)
  candidates <- which(loglik >= before & loglik >= after)
  candidates <- candidates[order(loglik[candidates], decreasing = TRUE)]
  kept <- integer()
  for (index in candidates) {
    separate <- vapply(kept, function(higher) {
      between <- seq(min(index, higher), max(index, higher))
      loglik[[index]] - min(loglik[between]) > drop
    }, logical(1L))
    if (all(separate)) kept <- c(kept, index)
  }
  kept
}
