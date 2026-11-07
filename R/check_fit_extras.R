##' The extra checks of `check_fit()`, each with what it computes
##'
##' The names are those that `extra` takes, in the order the checks run. The
##' description completes a sentence such as "Computing ...".
##' @noRd
list_extra_checks <- function() {
  c(profile = "the profile of the spatial scale")
}

##' The extra checks asked for, in the order they run
##'
##' `"all"` gives every one. An unknown name stops with the valid ones.
##' @noRd
match_extra_checks <- function(extra) {
  stopifnot("'extra' must be a character vector" =
              is.character(extra) && !anyNA(extra))
  known <- names(list_extra_checks())
  unknown <- setdiff(extra, c(known, "all"))
  if (length(unknown)) {
    stop("Unknown extra check", if (length(unknown) > 1L) "s", " ",
         paste0("\"", unknown, "\"", collapse = ", "), ". The extra checks ",
         "are ", paste0("\"", known, "\"", collapse = ", "),
         ", or \"all\" for every one.", call. = FALSE)
  }
  if ("all" %in% extra) known else known[known %in% extra]
}

##' Stop unless `profile` is NULL or a profile of `phi` for this fit
##'
##' A profile records the estimate it was centred on, which must be that of
##' the fit.
##' @noRd
check_given_profile <- function(profile, object) {
  if (is.null(profile)) return(invisible(TRUE))
  stopifnot(
    "'profile' must be a profile of this fit with a marginal profile of 'phi'" =
      inherits(profile, "RiskMap_profile") &&
      "phi" %in% names(profile$marginal) &&
      isTRUE(all.equal(unname(profile$estimate[["phi"]]),
                       unname(object$estimate$phi)))
  )
  invisible(TRUE)
}

##' Run the extra checks named in `extra`
##'
##' Gives the `results` of those that ran, each with its row and data, and
##' `not_computed`, a sentence for each of the others saying why. A profile
##' given by the user is used as it is. The random number state of the
##' session is restored on exit.
##' @noRd
run_extra_checks <- function(object, extra, profile = NULL, messages = TRUE) {
  results <- list()
  not_computed <- character()
  if (!length(extra)) {
    return(list(results = results, not_computed = not_computed))
  }
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  descriptions <- list_extra_checks()
  obstacles <- vapply(extra, function(name) {
    if (identical(name, "profile") && !is.null(profile)) return(NA_character_)
    find_extra_check_obstacle(name, object) %||% NA_character_
  }, character(1L))
  computed <- extra[is.na(obstacles) &
                      !(extra == "profile" & !is.null(profile))]
  if (messages && length(computed)) {
    message(describe_extra_check_time(object, computed))
  }
  for (name in extra) {
    result <- if (is.na(obstacles[[name]])) {
      tryCatch(run_extra_check(name, object, profile),
               error = function(error) error)
    }
    reason <- if (!is.na(obstacles[[name]])) {
      obstacles[[name]]
    } else if (inherits(result, "error")) {
      paste0("its computation stopped with the error \"",
             conditionMessage(result), "\"")
    }
    if (is.null(reason)) {
      results[[name]] <- result
    } else {
      not_computed[[name]] <- paste0(
        toupper(substr(descriptions[[name]], 1L, 1L)),
        substring(descriptions[[name]], 2L), " was not computed, because ",
        reason, "."
      )
    }
  }
  list(results = results, not_computed = not_computed)
}

##' Why an extra check does not apply to a fit, NULL when it does
##'
##' The profile of the scale describes a maximum of the likelihood.
##' @noRd
find_extra_check_obstacle <- function(name, object) {
  if (isTRUE(attr(object, "optimizer")$mcml_no_maximum)) {
    return("the likelihood has no maximum")
  }
  NULL
}

##' Run one extra check, giving its row and data
##' @noRd
run_extra_check <- function(name, object, profile = NULL) {
  switch(name, profile = run_profile_check(object, profile))
}

##' The seed of an extra check, a stream of its own after the seed of the fit
##'
##' Wrapped as `fit_adaptive_mcml()` wraps the streams of a fit, which stay
##' below 10,000 after its seed and do not reach the offsets used here. NULL
##' for a fit without a seed, such as a Gaussian fit.
##' @noRd
make_extra_check_seed <- function(fit_seed, offset) {
  if (is.null(fit_seed)) return(NULL)
  as.integer(((as.numeric(fit_seed) - 1 + offset) %% .Machine$integer.max) + 1)
}

##' How long the extra checks should take, from the time the fit took
##'
##' With the default settings of `profile.RiskMap()`, a profile of the scale
##' of a binomial or Poisson fit took between 1.8 and 14 times as long as the
##' fit, with a median of 8.5, on 11 fits of one stage to the liberia,
##' anopheles and loaloa data with and without a nugget, and 4.1 times on a
##' fit of two stages, whose profile costs no more. Ten times errs on the long
##' side. A Gaussian fit records no time, and its exact profile is quick.
##' @noRd
describe_extra_check_time <- function(object, extras) {
  descriptions <- list_extra_checks()[extras]
  listed <- paste(descriptions, collapse = " and ")
  fit_seconds <- attr(object, "optimizer")$elapsed_seconds
  if (is.null(fit_seconds)) return(paste0("Computing ", listed, "."))
  multiples <- c(profile = 10)[extras]
  paste0("Computing ", listed, ", which should take about ",
         format_check_time(sum(multiples) * fit_seconds), ", from the ",
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
##' The profile has the default settings of `profile.RiskMap()`, with its bank
##' drawn from the stream 9600 after the seed of the fit.
##' @noRd
run_profile_check <- function(object, profile = NULL) {
  if (is.null(profile)) {
    profile <- profile.RiskMap(
      object, parameters = "phi",
      seed = make_extra_check_seed(object$seed, 9600L) %||% 1L,
      messages = FALSE
    )
  }
  advice <- if (!identical(object$family, "gaussian")) {
    describe_pc_prior_advice(attr(object, "optimizer"))
  } else {
    ""
  }
  check_result(fit_check_profile(profile, object$distance_units, advice),
               profile)
}

##' What the profile of the spatial scale shows
##'
##' A point of the profile above the estimate shows that the estimate is not
##' the maximum likelihood estimate (`higher_than_estimate`, above 0.01 as in
##' `profile.RiskMap()`). A second maximum (`find_profile_maxima()`) or a 95%
##' interval open at an end, where the profile stays within the cut-off to
##' the end of its grid or the bound of the scale, means that Wald intervals
##' describe the uncertainty poorly. `advice` adds the PC penalty to the fix
##' of a flat profile.
##' @noRd
fit_check_profile <- function(profile, units = NULL, advice = "") {
  name <- "Profile of the scale"
  units <- paste0(" ", units %||% "distance units")
  marginal <- profile$marginal$phi
  marginal <- marginal[is.finite(marginal$loglik), ]
  marginal <- marginal[order(marginal$value), ]
  scale_at <- function(index) signif(exp(marginal$value[[index]]), 3)
  interval <- profile$intervals[profile$intervals$parameter == "phi", ]
  level <- paste0(round(100 * profile$level), "%")
  higher <- profile$higher_than_estimate
  if (isTRUE(higher > 0.01)) {
    top <- which.max(marginal$loglik)
    return(fit_check_row(
      name, "PROBLEM",
      paste0("The profile at a scale of ", scale_at(top), units, " is above ",
             "the likelihood of the estimate by ", signif(higher, 2),
             " log-likelihood units, so the estimate is not the maximum ",
             "likelihood estimate."),
      paste0("Refit from that scale, with 'start_pars = list(phi = ",
             scale_at(top), ")' in glgpm(), and keep the fit with the higher ",
             "log-likelihood.")
    ))
  }
  maxima <- find_profile_maxima(marginal$loglik)
  open <- c(below = is.na(interval$profile_lower),
            above = is.na(interval$profile_upper))
  if (length(maxima) < 2L && !any(open)) {
    return(fit_check_row(name, "OK", paste0(
      "The profile of the spatial scale has one maximum, and its ", level,
      " interval is ", signif(interval$profile_lower, 3), " to ",
      signif(interval$profile_upper, 3), units, "."
    )))
  }
  below <- if (length(maxima) >= 2L) {
    max(marginal$loglik) - marginal$loglik[[maxima[2L]]]
  }
  fit_check_row(
    name, "SUSPICIOUS",
    paste0(
      if (length(maxima) >= 2L) paste0(
        "The profile of the spatial scale has a second maximum at ",
        scale_at(maxima[2L]), units, ", ", if (below < 0.01) {
          "as high as the highest"
        } else {
          paste0(signif(below, 2), " log-likelihood units below the highest")
        }, ". "
      ),
      if (any(open)) paste0(
        if (length(maxima) >= 2L) "It" else "The profile of the spatial scale",
        " is flat ", if (all(open)) {
          "towards short and long scales, so its "
        } else if (open[["below"]]) {
          "towards short scales, so its "
        } else {
          "towards long scales, so its "
        }, level, " interval is open ",
        if (all(open)) "at both ends" else if (open[["below"]]) "below" else "above",
        " and the data do not bound the scale there. "
      ),
      "Wald intervals may describe the uncertainty poorly."
    ),
    paste0(
      "Use the profile interval of the scale rather than its Wald interval, ",
      "and see the profile panel of plot().",
      if (any(open)) paste0(" Consider a scale fixed from external ",
                            "knowledge.", advice)
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
