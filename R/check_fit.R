##' @title Check How a Fitted Model Was Computed
##' @description Goes through the checks on the computation of a model fitted
##'   by [glgpm()], says which are OK and which need a look, and gives the fix
##'   where there is one. Checks that cost computing time run only when `extra`
##'   names them. `print()` gives the verdict, `summary()` the numbers behind
##'   each finding and `plot()` the figures. The model's own `summary()` gives
##'   only the number of checks that need a look.
##'
##' @details The checks concern the computation, not whether the model suits the
##'   data, so `OK` does not mean that the model is adequate. They read what the
##'   fit recorded, and add no statistical judgement of their own. The status
##'   `PROBLEM` means that the estimate or its certificate should not be relied
##'   on as it stands. `SUSPICIOUS` means that the computation is valid but the
##'   result needs care. The limit of 0.7 on the tail shape of the importance
##'   weights, the Pareto shape, is that of Vehtari et al. (2024).
##'
##'   A Gaussian fit has no Monte Carlo error, so its only default check says so.
##'
##'   A link other than logit, probit or cloglog for binomial data, or log for
##'   counts, adds a row saying that it is not known to make the likelihood of
##'   each location log-concave in its latent effect, so that the latent
##'   effects can have more than one conditional mode. That belongs to the
##'   model and not to how it was fitted, so the row is `OK`.
##'
##'   The default checks use only what the fit stored, so they take no time.
##'   The extra checks compute something new and add a row each. Before they
##'   start, a message says how long they should take, from the time the fit
##'   took. The extra checks are
##'   * `"profile"`, the marginal profile likelihood of each covariance
##'     parameter of the fit from [profile.RiskMap()], with its default
##'     settings. These are the spatial variance and scale, the nugget (its
##'     ratio to the spatial variance, `nu2`, or for a Gaussian fit its
##'     variance, `tau2`), an estimated measurement error variance and the
##'     variances of the random effects, but not the regression coefficients.
##'     Its row names each parameter whose profile is flat, so that its 95%
##'     interval is open at an end, or has a second maximum, where it falls by
##'     more than half a log-likelihood unit between two maxima, and is then
##'     `SUSPICIOUS`. It is a `PROBLEM` when a point of a profile is above the
##'     estimate, which is then not the maximum likelihood estimate. The
##'     profiles of a binomial or Poisson fit should take about 10 times as
##'     long as the fit for each parameter, less with `cores`.
##'   * `"radius"`, for a binomial or Poisson fit, a prediction of how often
##'     the radius in the bound on the Monte Carlo error is too small. The
##'     radius is estimated from the certification draws themselves, which
##'     can miss the rare draws of large importance weight. Short Markov
##'     chains on the latent effects, 16 chains of 4,000 rounds, reach those
##'     draws, and samples of certification draws resampled from them and a
##'     fresh sample of the fit's size show how often the error that the
##'     radius allows for exceeds it, against the 5% of the bound. Outside
##'     binomial fits with the logit link the prediction is read with a margin
##'     of 15%, since it ran low on simulated Poisson surveys. Its row is
##'     `SUSPICIOUS` when the radius, with that margin, would have to widen for
##'     the rate to be 5%, and then gives the bound with the radius so widened.
##'     The chains can miss a separate mode of the latent effects, so `OK` is
##'     not a proof. The check should take about twice as long as the fit.
##'
##'   Each extra check draws from random number streams of its own, derived
##'   from the seed of the fit (`fit$seed`) as the streams of the fit are, so
##'   the same call gives the same result. The profiles use the stream 9600
##'   after the seed and the check of the radius the streams 9500 to 9502,
##'   which no fit uses. The random number state of the session is left as it
##'   was. An extra check that does not apply to the fit, or that
##'   fails, is not computed, and the printed checks say why.
##'
##' @param object An object of class `"RiskMap"` from [glgpm()].
##' @param extra Names of the extra checks to run, `"all"` for every one. See
##'   Details. The default is none.
##' @param profile An optional object of class `"RiskMap_profile"` from
##'   [profile.RiskMap()] with marginal profiles for this fit. The profile
##'   check then uses it as it is, whatever parameters it covers, instead of
##'   computing one, whether or not `extra` names it.
##' @param messages Logical, whether to say how long the extra checks should
##'   take before they start.
##' @param cores Number of processes for the extra checks that can use
##'   several, passed to [profile.RiskMap()], whose results are the same on
##'   any number of cores.
##' @return An object of class `"RiskMap_check"`, a list with:
##' \item{checks}{A data frame with one row per check, those that need a look
##' first, and the columns `check`, `status` (`"PROBLEM"`, `"SUSPICIOUS"` or
##' `"OK"`), `finding`, `fix`, which is empty where there is nothing to do,
##' and `brief`, the finding in a few words, which `print()` gives for a check
##' that is `OK`.}
##' \item{draws}{One sentence on the stages and draws of the fit, `NULL` for
##' Gaussian models.}
##' \item{stages}{For non-Gaussian models, a data frame with one row per stage
##' and the columns `stage`, `proposal`, `certificate`, `radius_from`,
##' `effective_draws` and `pareto_k`. The radius of the certificate is the
##' larger of the radii of the certification draws and of a replicate of them,
##' and `radius_from` says which it was (`"certification"` or `"replicate"`).
##' The last two columns describe the certification draws. All three are `NA`
##' for a stage that did not reach them.}
##' \item{moves}{A data frame with the absolute change in each `parameter`
##' between one `stage` and the one before it, in standard errors (`move`).
##' `NULL` for a fit with one stage or without standard errors.}
##' \item{tolerance}{The tolerance on the certificate.}
##' \item{results}{A list with an element for each default check that ran,
##' named by the check. Each holds its `row` of `checks` and its `data`, a list
##' of the numbers behind the finding.}
##' \item{extras}{A list with an element for each extra check that ran, named
##' as in `extra`, each with its `row` and its `data`. The `data` of
##' `"profile"` is the `"RiskMap_profile"` object. That of `"radius"` is a list
##' with the predicted percentages of samples whose error exceeds the radius
##' (`level`, with its `level_interval`), in which the error of the
##' certification draws alone exceeds the radius (`own`), the share above the
##' line in its panel, and whose error exceeds the bound (`exceedance`), the `factor` by which the radius would have to grow for
##' `level` to be 5%, the `margin`, the `widening` of the radius with that
##' margin, at least one, the bound so widened (`widened`), the `settings` of
##' the check and `samples`, a data frame with the `error` of each resampled
##' sample and the `radius` the fit would use for it, the larger of its own
##' and a replicate's, in standard errors.}
##' \item{not_computed}{A named character vector giving, for each extra check
##' that was asked for but not computed, the reason.}
##' \item{fit_seconds}{The time the fit took in seconds, `NULL` where the fit
##' did not record it.}
##' @seealso [summary.RiskMap_check()], [plot.RiskMap_check()],
##'   [summary.RiskMap()], [profile.RiskMap()]
##' @examples
##' data(italy_sim)
##'
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' check_fit(fit)
##' \donttest{
##' checked <- check_fit(fit, extra = "profile")
##' checked
##' summary(checked)
##' plot(checked, which = "profile")
##' }
##' @export
check_fit <- function(object, extra = character(), profile = NULL,
                      messages = TRUE, cores = 1L) {
  check_class(object, "RiskMap")
  extra <- match_extra_checks(extra)
  check_given_profile(profile, object)
  check_logical(messages)
  check_positive_integer(cores)
  if (!is.null(profile)) extra <- union(extra, "profile")
  mcml <- summarise_mcml_fit(object)
  history <- object$mcml_history %||% attr(object, "optimizer")$mcml_history
  extras <- run_extra_checks(object, extra, profile, messages,
                             as.integer(cores))
  results <- run_fit_checks(object, mcml, names(extras$results))
  structure(
    list(
      checks = tabulate_fit_checks(c(results, extras$results)),
      draws = if (!is.null(mcml)) describe_mcml_draws(mcml),
      stages = summarise_mcml_stages(history),
      moves = summarise_estimate_moves(object, history),
      tolerance = mcml$tolerance,
      results = results,
      extras = extras$results,
      not_computed = extras$not_computed,
      fit_seconds = attr(object, "optimizer")$elapsed_seconds
    ),
    class = "RiskMap_check"
  )
}

##' The table of checks, those that need a look first
##'
##' `results` is a list of check results, each with its `row`.
##' @noRd
tabulate_fit_checks <- function(results) {
  checks <- do.call(rbind, lapply(unname(results), `[[`, "row"))
  severity <- match(checks$status, c("PROBLEM", "SUSPICIOUS", "OK"))
  checks <- checks[order(severity), ]
  rownames(checks) <- NULL
  checks
}

##' The table of the default checks on a fit, those that need a look first
##'
##' `mcml` is the list of `summarise_mcml_fit()`, NULL for a Gaussian fit.
##' @noRd
list_fit_checks <- function(object, mcml) {
  tabulate_fit_checks(run_fit_checks(object, mcml))
}

##' The default checks on a fit, each with the numbers behind its finding
##'
##' Gives a list with an element for each check that applies, named by the
##' check, holding its `row` and its `data` (`collect_check_data()`). `mcml`
##' is the list of `summarise_mcml_fit()`, NULL for a Gaussian fit. A fit whose
##' likelihood has no maximum gets that check alone, since the others describe
##' a maximum. `extras` names the extra checks that ran, so that a fix points
##' to the row of an extra check that ran and otherwise to the call that runs
##' it.
##' @noRd
run_fit_checks <- function(object, mcml, extras = character()) {
  optimizer <- attr(object, "optimizer")
  radius <- if (is.null(find_radius_check_obstacle(object))) {
    describe_radius_pointer("radius" %in% extras)
  } else {
    ""
  }
  profile <- describe_profile_pointer("profile" %in% extras)
  rows <- if (is.null(mcml)) {
    list(fit_check_row("Monte Carlo accuracy", "OK",
                       "A Gaussian fit has no Monte Carlo error."))
  } else if (isTRUE(optimizer$mcml_no_maximum)) {
    list(fit_check_no_maximum(optimizer, object$family))
  } else {
    list(fit_check_monte_carlo(mcml),
         fit_check_weights(mcml$certification_bank,
                           isTRUE(mcml$suggest_integration), radius),
         fit_check_maxima(optimizer, object$distance_units, profile),
         fit_check_variances(optimizer),
         fit_check_spatial_scale(object, optimizer, profile),
         fit_check_importance_distribution(optimizer),
         fit_check_mode_search(optimizer$mcml_mode_searches),
         fit_check_link(object$link_function, object$family))
  }
  rows <- Filter(Negate(is.null), rows)
  names(rows) <- vapply(rows, `[[`, "", "check")
  lapply(rows, function(row) {
    check_result(row, if (!is.null(mcml)) {
      collect_check_data(row$check, object, mcml)
    })
  })
}

##' The result of a check, its row of the table and the data behind it
##'
##' `data` holds the numbers that the finding rests on, which `summary()`
##' prints and `plot()` draws from.
##' @noRd
check_result <- function(row, data = NULL) {
  list(row = row, data = data)
}

##' One row of the table of checks
##'
##' `brief` is the finding in a few words, which `print()` gives on one line
##' for a check that is OK.
##' @noRd
fit_check_row <- function(check, status, finding, fix = "", brief = finding) {
  data.frame(check = check, status = status, finding = finding, fix = fix,
             brief = brief)
}

##' Text with its first letter in capitals
##' @noRd
capitalise_first <- function(text) {
  paste0(toupper(substr(text, 1L, 1L)), substring(text, 2L))
}

##' Where a fix sends the user for the profile likelihood
##'
##' The phrase completes "see the profile likelihood", with the call of
##' `check_fit()` that computes the profiles or, when they ran, their row and
##' panels.
##' @noRd
describe_profile_pointer <- function(ran = FALSE) {
  if (ran) {
    return(" in the check \"Profile likelihood\" and its panels in plot()")
  }
  ", check_fit(fit, extra = \"profile\")"
}

##' Where a fix sends the user for the check of the radius of the bound
##'
##' A sentence for the fix of heavy-tailed weights, which can make the radius
##' too small.
##' @noRd
describe_radius_pointer <- function(ran = FALSE) {
  if (ran) return(" See the check \"Radius of the bound\".")
  paste0(" The check of the radius, check_fit(fit, extra = \"radius\"), ",
         "predicts whether it is too small.")
}

##' The advice to use the PC penalty, unless the fit already does
##'
##' Weak identification of the spatial variance or scale is what the penalty
##' addresses.
##' @noRd
describe_pc_prior_advice <- function(optimizer) {
  if (isTRUE(optimizer$mcml_pc_prior)) return("")
  paste0(" Or refit with 'control_mcml = set_control_mcml(pc_prior = TRUE)', a ",
         "mild penalty that stabilises the spatial variance and scale.")
}

##' Whether the certificate bounds the Monte Carlo error within the tolerance
##' @noRd
fit_check_monte_carlo <- function(mcml) {
  name <- "Monte Carlo accuracy"
  bound <- signif(mcml$certificate, 2)
  see_scale <- "See the check of the spatial scale."
  switch(
    mcml$outcome,
    "certified" = fit_check_row(name, "OK", paste0(
      "The error of the estimate is within ", bound, " standard errors (limit ",
      mcml$tolerance, ")."
    ), brief = paste0("Within ", bound, " standard errors.")),
    "weak identification" = fit_check_row(
      name, "PROBLEM",
      "There is no bound on the error, because the observed information is singular.",
      see_scale
    ),
    "numerically independent" = fit_check_row(
      name, "PROBLEM",
      paste0("There is no bound on the error, because the fit is numerically ",
             "one of independent location effects."),
      see_scale
    ),
    if (isTRUE(mcml$variance_rising)) {
      fit_check_row(
        name, "PROBLEM",
        "There is no bound on the error, because the estimate is not a maximum.",
        "See the check of the variances."
      )
    } else if (isTRUE(mcml$correction_left_bounds)) {
      fit_check_row(
        name, "PROBLEM",
        paste0("The bound of ", bound, " standard errors may not hold, because ",
               "the Newton correction of the estimate would have taken a ",
               "variance or scale past the upper bound of the search, along a ",
               "direction where the likelihood is nearly flat."),
        "Refit with a larger 'max_samples' in set_control_mcml()."
      )
    } else if (mcml$correction_refused) {
      fit_check_row(
        name, "PROBLEM",
        paste0("The bound of ", bound, " standard errors may not hold, because ",
               "the Newton correction of the estimate did not converge."),
        "Refit with a larger 'max_samples' in set_control_mcml()."
      )
    } else {
      fit_check_row(
        name, "PROBLEM",
        paste0("The bound on the error is ", bound, " standard errors, above ",
               "the limit of ", mcml$tolerance, "."),
        paste0("Refit with a larger 'max_samples' or 'max_stages' in ",
               "set_control_mcml(), or relax 'tolerance'.",
               describe_integration_advice(isTRUE(mcml$suggest_integration)))
      )
    }
  )
}

##' Whether the importance weights of the certification draws have a heavy tail
##'
##' Above a Pareto shape of 0.7 the Gaussian approximation behind the radius of
##' the certificate is unreliable at the usual numbers of draws (Vehtari et
##' al., 2024). NULL for a fit that did not reach those draws. With
##' `suggest_integration`, the fix also suggests integrating the nugget out.
##' `radius` is the sentence of the fix on the check of the radius
##' (`describe_radius_pointer()`), empty where that check does not apply.
##' @noRd
fit_check_weights <- function(bank, suggest_integration = FALSE, radius = "") {
  if (is.null(bank)) return(NULL)
  name <- "Importance weights"
  draws <- format(round(bank$effective_draws), big.mark = ",", trim = TRUE)
  found <- paste0(
    "The certification draws had ", draws, " effective draws and a tail ",
    "shape of ", signif(bank$pareto_k, 2), " (limit 0.7)"
  )
  if (isTRUE(bank$pareto_k <= 0.7)) {
    return(fit_check_row(name, "OK", paste0(found, "."), brief = paste0(
      draws, " effective draws, tail shape ", signif(bank$pareto_k, 2), "."
    )))
  }
  fit_check_row(
    name, "SUSPICIOUS",
    paste0(found, ", which makes the bound on the error less reliable."),
    paste0("More draws may help, through 'max_samples' in set_control_mcml(). ",
           "If the shape stays above 0.7, read the bound with caution.",
           radius, describe_integration_advice(suggest_integration))
  )
}

##' The advice to integrate the nugget or random effects out, for a fit that
##' sampled them
##' @noRd
describe_integration_advice <- function(suggest) {
  if (!suggest) return("")
  paste0(" Or refit with the default 'integrate_noise = \"auto\"' in ",
         "set_control_mcml(), which integrates the nugget and random effects ",
         "within locations out where their variances are large and often ",
         "keeps far more effective draws.")
}

##' What the check for other maxima found
##'
##' An exact alarm shows that the estimate is not the global maximum. The check
##' does not prove that no higher maximum exists, which the certificate says
##' nothing about either. `profile` completes the fix
##' (`describe_profile_pointer()`).
##' @noRd
fit_check_maxima <- function(optimizer, units = NULL,
                             profile = describe_profile_pointer()) {
  name <- "Other maxima"
  ambiguous <- isTRUE(optimizer$mcml_ambiguous_maxima)
  better <- isTRUE(optimizer$mcml_better_maximum)
  if (isTRUE(optimizer$mcml_maxima_alarm)) {
    return(fit_check_row(
      name, "PROBLEM",
      paste0("An exact lower bound on the likelihood at another spatial scale ",
             "is above the likelihood of the estimate, so the estimate is not ",
             "the maximum likelihood estimate."),
      paste0("Inspect the profile likelihood", profile, ".")
    ))
  }
  if (!ambiguous && !better) {
    return(fit_check_row(
      name, "OK",
      "The check found no higher maximum, though it cannot rule one out.",
      brief = "The check found no higher maximum."
    ))
  }
  fit_check_row(
    name, "SUSPICIOUS",
    paste0(
      if (ambiguous) {
        "Separate maxima fit about equally well, and the estimate is from the higher."
      } else {
        paste0("The check found a maximum higher than the one reached from the ",
               "starting values, and the estimate is from it.")
      },
      describe_other_maximum(optimizer$mcml_maxima_check, ambiguous, units)
    ),
    paste0("Wald intervals may describe the uncertainty poorly, so see the ",
           "profile likelihood", profile, ".")
  )
}

##' How far the estimate's maximum is from the other, on independent draws
##'
##' The other maximum is the fitted one when the estimate moved, otherwise the
##' closest rival. `check` is on the user's distance scale.
##' @noRd
describe_other_maximum <- function(check, ambiguous, units = NULL) {
  if (is.null(check)) return("")
  judged <- check[is.finite(check$independent_log_likelihood), ]
  kept <- judged[judged$chosen, ]
  rest <- judged[!judged$chosen, ]
  other <- if (ambiguous) {
    rest[which.max(rest$independent_log_likelihood), ]
  } else {
    rest[rest$maximum == check$maximum[[1L]], ][1L, ]
  }
  paste0(
    " The two differ by ",
    signif(kept$independent_log_likelihood - other$independent_log_likelihood, 2),
    " log-likelihood units (standard error ", signif(other$independent_se, 2),
    ") on independent draws, and the other has a spatial scale of ",
    signif(exp(other$log_phi), 3), if (!is.null(units)) paste0(" ", units), "."
  )
}

##' Whether a variance component was estimated as effectively zero
##'
##' A nugget or random-effect variance at zero has been tested by the slope of
##' the likelihood there (`assess_held_variances()`). If the likelihood clearly
##' rises, zero is not a maximum, and if it clearly falls, it is one and there
##' is nothing to review. A spatial variance at zero has no such test, because
##' its slope depends on the scale, which is unidentified there.
##' @noRd
fit_check_variances <- function(optimizer) {
  name <- "Variances"
  boundary <- optimizer$mcml_boundary
  if (!length(boundary)) {
    return(fit_check_row(name, "OK",
                         "No variance was estimated as effectively zero.",
                         brief = "No variance was estimated as zero."))
  }
  held <- optimizer$mcml_held_variances
  listed <- function(names) paste(names, collapse = " and the ")
  if (any(held$rising)) {
    return(fit_check_row(
      name, "PROBLEM",
      paste0("The likelihood rises as the ", listed(held$name[held$rising]),
             " leaves zero, so holding it at zero does not give a maximum."),
      paste0("Refit with a larger 'max_samples' or 'max_stages' in ",
             "set_control_mcml().")
    ))
  }
  if (setequal(boundary, held$name) && all(held$falling)) {
    symbols <- sub("^.*\\((.*)\\)$", "\\1", boundary)
    return(fit_check_row(name, "OK", paste0(
      "The ", listed(boundary), " was estimated as zero, where the likelihood ",
      "is at its maximum, because it falls as the variance leaves zero."
    ), brief = paste0("Zero is the maximum for ", join_phrases(symbols), ".")))
  }
  single <- length(boundary) == 1L
  spatial <- any(grepl("sigma2", boundary, fixed = TRUE))
  nugget_alone <- single && grepl("tau2", boundary, fixed = TRUE)
  fit_check_row(
    name, if (spatial) "PROBLEM" else "SUSPICIOUS",
    paste0(
      "The ", paste(boundary, collapse = " and the "), " ",
      if (single) "was" else "were", " estimated as effectively zero, so ",
      if (single) "its standard error is" else "their standard errors are",
      " not reported", if (spatial) ", and the spatial scale cannot be estimated",
      ".",
      if (!is.null(held) && !all(held$falling)) paste0(
        " The draws could not tell whether the likelihood rises as the ",
        listed(held$name[!held$falling]), " leaves zero."
      )
    ),
    paste0(
      "Consider refitting without ",
      if (nugget_alone) "it, gp(nugget = FALSE)." else if (single) "it." else "them.",
      if (spatial) describe_pc_prior_advice(optimizer)
    )
  )
}

##' Whether the data determine the spatial scale
##'
##' NULL with the spatial variance at zero, where there is no scale to check.
##' `profile` completes the fix (`describe_profile_pointer()`).
##' @noRd
fit_check_spatial_scale <- function(object, optimizer,
                                    profile = describe_profile_pointer()) {
  name <- "Spatial scale"
  advice <- describe_pc_prior_advice(optimizer)
  see_profile <- paste0("See the profile likelihood", profile, ".")
  if (isTRUE(optimizer$mcml_identifiability_failure)) {
    return(fit_check_row(
      name, "PROBLEM",
      "Some parameters, often the spatial scale, are not identifiable from these data.",
      paste0("Consider a simpler model or fixing the scale.", advice)
    ))
  }
  if (isTRUE(optimizer$mcml_numerically_independent)) {
    return(fit_check_row(
      name, "PROBLEM",
      paste0("The scale is so short that the correlation between the closest ",
             "locations is below 1e-6, so it is unresolved."),
      paste0(see_profile, " If spatial correlation is not expected, consider an ",
             "unstructured effect for each location instead of gp().", advice)
    ))
  }
  if (any(grepl("sigma2", optimizer$mcml_boundary, fixed = TRUE))) return(NULL)
  units <- object$distance_units %||% "distance units"
  value <- paste0(signif(stats::coef(object)$phi, 3), " ", units)
  scale <- paste0("The scale of ", value)
  small <- isTRUE(optimizer$mcml_small_scale)
  near <- isTRUE(optimizer$mcml_near_limit)
  if (!small && !near) {
    return(fit_check_row(name, "OK", paste0(
      scale, " is longer than the median distance between neighbouring locations."
    ), brief = paste0(value, ", longer than the median spacing.")))
  }
  fit_check_row(
    name, "SUSPICIOUS",
    paste0(
      scale,
      if (small) {
        paste0(" is shorter than the median distance between neighbouring ",
               "locations",
               if (!is.null(optimizer$mcml_spacing)) paste0(
                 " (", signif(optimizer$mcml_spacing, 3), " ", units, ")"
               ),
               ", so it may be weakly identified and its interval very wide.")
      } else {
        " is short."
      },
      if (near) paste0(" The fitted likelihood is not clearly above the limit ",
                       "with independent location effects.")
    ),
    paste0(
      see_profile, " Consider ",
      if (is.numeric(object$fix_tau2)) "a nugget, gp(nugget = TRUE), or ",
      "a scale fixed from external knowledge.", advice
    )
  )
}

##' Which importance distribution the fit used
##'
##' A change of distribution does not change how a fit is read, because the
##' certificate is computed from the draws of the distribution used, so this
##' check is always OK. It also says how often the distribution could not be
##' built at a parameter value the search tried, which gives that value no
##' weight.
##' @noRd
fit_check_importance_distribution <- function(optimizer) {
  fallbacks <- optimizer$mcml_proposal_fallbacks
  errors <- optimizer$mcml_proposal_errors
  errors_found <- if (isTRUE(errors$count > 0L)) paste0(
    " It could not be built at ", format(errors$count, big.mark = ","),
    if (errors$count == 1L) " parameter value" else " parameter values",
    " that the search tried, the first time with the error: ", errors$first
  )
  restart <- if (isTRUE(optimizer$mcml_start_replaced)) paste0(
    "The global distribution could not be built at the Laplace estimate, so ",
    "the fit restarted from the starting values"
  )
  used <- if (identical(optimizer$mcml_proposal, "laplace")) {
    "The Laplace approximation"
  } else {
    "The global distribution"
  }
  brief <- c(
    if (!is.null(restart)) "restarted from the starting values",
    if (length(fallbacks)) paste0(
      "the Laplace approximation at stage", if (length(fallbacks) > 1L) "s",
      " ", paste(fallbacks, collapse = ", ")
    ),
    if (!is.null(errors_found)) paste0(
      "not built at ", format(errors$count, big.mark = ","),
      if (errors$count == 1L) " parameter value" else " parameter values"
    )
  )
  brief <- if (length(brief)) {
    paste0(capitalise_first(paste(brief, collapse = ", ")), ".")
  } else {
    paste0(used, " at every stage.")
  }
  finding <- paste0(
    if (length(fallbacks)) {
      paste0(if (!is.null(restart)) paste0(restart, ". "),
             "The global distribution could not be built at stage",
             if (length(fallbacks) == 1L) " " else "s ",
             paste(fallbacks, collapse = ", "),
             ", which used the Laplace approximation instead.")
    } else if (!is.null(restart)) {
      paste0(restart, " and used it at every stage.")
    } else {
      paste0(used, " was used at every stage.")
    },
    errors_found
  )
  fit_check_row("Importance distribution", "OK", finding, brief = brief)
}

##' Whether the search for the conditional mode reached its iteration limit
##'
##' The search stops at the limit with its last iterate. A search that fails
##' from the mode found before starts again from zero
##' (`adaptive_laplace_mode()`). The mode only shapes the importance
##' distribution, so the cost is speed and not validity. NULL when no search
##' reached the limit or started again. A fit made before the restarts were
##' counted has no count of them, which is read as none.
##' @noRd
fit_check_mode_search <- function(counts) {
  restarted <- counts$restarted %||% 0L
  if (is.null(counts) || (counts$at_limit == 0L && restarted == 0L)) {
    return(NULL)
  }
  outcomes <- c(
    if (counts$at_limit > 0L) paste0(
      "stopped at its iteration limit in ",
      format(counts$at_limit, big.mark = ",")
    ),
    if (restarted > 0L) paste0(
      "failed from the mode found before and started again from zero in ",
      format(restarted, big.mark = ",")
    )
  )
  brief <- c(
    if (counts$at_limit > 0L) paste0(
      "at its limit in ", format(counts$at_limit, big.mark = ",")
    ),
    if (restarted > 0L) {
      paste0("restarted in ", format(restarted, big.mark = ","))
    }
  )
  searches <- paste0(" of ", format(counts$searches, big.mark = ","),
                     " searches")
  fit_check_row("Latent mode search", "OK", paste0(
    "The search for the mode of the latent effects ",
    paste(outcomes, collapse = ", and "), searches,
    ", which slows the fit but does not affect the certificate."
  ), brief = paste0(capitalise_first(paste(brief, collapse = ", ")), searches,
                    "."))
}

##' Whether the link is one known to make the likelihood log-concave
##'
##' With the logit, probit and complementary log-log links for binomial data,
##' and the log link for counts, the likelihood of each location is
##' log-concave in its latent effect, so the latent effects have one
##' conditional mode. Other links can give more than one, which belongs to the
##' model and the data, not to how the fit was computed, so the check is OK.
##' Concavity cannot be established numerically, so the links come from a
##' list. NULL for the links on it.
##' @noRd
fit_check_link <- function(link_function, family) {
  if (!identical(link_function$name, "custom") ||
      isTRUE(link_function$link %in% c("logit", "probit", "cloglog", "log"))) {
    return(NULL)
  }
  consequence <- paste0(
    "the latent effects can have more than one conditional mode. Every check ",
    "still runs, but the importance distribution can fit less well, so the ",
    "fit may take longer to certify."
  )
  finding <- if (is.null(link_function$link)) {
    paste0("The link was given as a function, so it is not known whether it ",
           "makes the likelihood of each location log-concave in its latent ",
           "effect, as ", if (identical(family, "poisson")) {
             "the log link does"
           } else {
             "the logit, probit and cloglog links do"
           }, ". If it does not, ",
           consequence)
  } else {
    paste0("The ", link_function$link, " link is not among those known to ",
           "make the likelihood of each location log-concave in its latent ",
           "effect, so ", consequence)
  }
  fit_check_row("Link", "OK", finding, brief = paste0(
    link_function$link %||% "A function", ", not known to be log-concave."
  ))
}

##' Why the likelihood has no maximum, and what can be done instead
##' @noRd
fit_check_no_maximum <- function(optimizer, family) {
  name <- "Maximum of the likelihood"
  binomial <- identical(family, "binomial")
  closing <- paste0(" The estimate is where the optimiser stopped, so it has no ",
                    "standard errors.")
  separation <- optimizer$mcml_separation
  everywhere <- separation$everywhere
  if (!is.null(everywhere)) {
    none <- identical(everywhere, "no positives")
    return(fit_check_row(
      name, "PROBLEM",
      paste0(
        "The likelihood has no maximum, because ",
        if (!binomial) {
          "every count is zero"
        } else if (none) {
          "no individual tested is positive"
        } else {
          "every individual tested is positive"
        },
        ".", closing
      ),
      paste0(
        "Maximum likelihood cannot fit the model to such data. An exact ",
        if (none) "upper" else "lower", " confidence limit for the ",
        if (binomial) "prevalence" else "rate", ", or a Bayesian analysis with ",
        "an informative prior on the intercept, can be used instead."
      )
    ))
  }
  if (!is.null(separation)) {
    return(fit_check_row(
      name, "PROBLEM",
      paste0("The likelihood has no maximum, because the covariates separate ",
             "the data, through ",
             paste0("'", separation$coefficients, "'", collapse = " and "), ".",
             closing),
      paste0(
        "Consider merging a category of a factor that has ",
        if (separation$only_positives && !separation$no_positives) {
          "only positives"
        } else {
          "no positives"
        },
        " with other categories, or removing the covariate from the model."
      )
    ))
  }
  fit_check_row(
    name, "PROBLEM",
    paste0("The likelihood appears to have no maximum, because a regression ",
           "coefficient reached the bound of the search.", closing),
    paste0("Check whether a covariate, or a category of a factor, has only ",
           "observations with no positives or only observations where every ",
           "individual is positive.")
  )
}

##' The certificate and the certification draws of each stage
##' @noRd
summarise_mcml_stages <- function(history) {
  if (!length(history)) return(NULL)
  do.call(rbind, lapply(seq_along(history), function(stage) {
    record <- history[[stage]]
    bank <- summarise_certification_bank(record)
    data.frame(stage = stage, proposal = record$proposal %||% NA_character_,
               certificate = record$certificate %||% NA_real_,
               radius_from = record$radius_from %||% NA_character_,
               effective_draws = bank$effective_draws %||% NA_real_,
               pareto_k = bank$pareto_k %||% NA_real_)
  }))
}

##' How far each estimate moved between stages, in standard errors
##'
##' NULL for a fit with one stage or without standard errors. A parameter held
##' at a bound has no standard error, so it is left out.
##' @noRd
summarise_estimate_moves <- function(object, history) {
  if (length(history) < 2L) return(NULL)
  standard_errors <- sqrt(diag(object$covariance))
  estimates <- do.call(cbind, lapply(history, function(stage) {
    as.numeric(unlist(stage$estimate))
  }))
  if (nrow(estimates) != length(standard_errors)) return(NULL)
  moves <- abs(estimates[, -1L, drop = FALSE] -
                 estimates[, -ncol(estimates), drop = FALSE]) / standard_errors
  table <- data.frame(
    parameter = rep(rownames(object$covariance), times = ncol(moves)),
    stage = rep(seq(2L, ncol(estimates)), each = nrow(moves)),
    move = as.vector(moves)
  )
  table <- table[is.finite(table$move), ]
  if (nrow(table)) table else NULL
}

##' @title Print the Checks on a Fitted Model
##' @description Print method for objects of class `"RiskMap_check"` from
##'   [check_fit()], which gives the verdict. Each check has an entry, those
##'   that need a look first. A check that needs a look is given in full,
##'   with its fix where there is one, and a check that is `OK` takes one line,
##'   with its finding in a few words. The extra checks that were asked for
##'   but not computed are listed with the reason. [summary.RiskMap_check()]
##'   gives each finding in full with the numbers behind it.
##' @param x An object of class `"RiskMap_check"`.
##' @param ... other parameters.
##' @return Invisibly returns `x`.
##' @method print RiskMap_check
##' @export
print.RiskMap_check <- function(x, ...) {
  check_class(x, "RiskMap_check")
  cat("Checks on this fit\n")
  if (!is.null(x$draws)) write_wrapped(x$draws)
  cat("\n")
  checks <- x$checks
  status_width <- max(nchar(checks$status))
  name_width <- max(nchar(checks$check))
  finding_start <- status_width + name_width + 4L
  width <- max(finding_start + 30L, getOption("width", 80L) - 2L)
  # A table made before the column of brief findings gives the findings.
  briefs <- checks$brief %||% checks$finding
  for (row in seq_len(nrow(checks))) {
    check <- checks[row, ]
    ok <- identical(check$status, "OK")
    finding <- strwrap(if (ok) briefs[[row]] else check$finding,
                       width = width - finding_start)
    cat(formatC(check$status, width = -status_width), "  ",
        formatC(check$check, width = -name_width), "  ", finding[1L], "\n",
        sep = "")
    if (length(finding) > 1L) {
      cat(paste0(strrep(" ", finding_start), finding[-1L]), sep = "\n")
    }
    if (!ok && nzchar(check$fix)) {
      cat(strwrap(paste("Fix:", check$fix), width = width,
                  indent = status_width + 2L, exdent = status_width + 7L),
          sep = "\n")
    }
  }
  if (length(x$not_computed)) {
    cat("\n")
    write_wrapped(x$not_computed)
  }
  cat("\n")
  write_wrapped(paste0("summary() of these checks gives the evidence behind ",
                       "each finding."))
  write_wrapped(paste0("These checks concern the computation, not whether the ",
                       "model suits the data."))
  invisible(x)
}

##' @title Plot the Checks on a Fitted Model
##' @description Plot method for objects of class `"RiskMap_check"` from
##'   [check_fit()]. It draws the figures of the checks that were computed,
##'   and never computes anything itself, so an extra check runs once in
##'   [check_fit()]. By stage of certified Monte Carlo maximum likelihood, the
##'   panels show the bound on the Monte Carlo error against the tolerance
##'   (`"certificate"`), the effective certification draws
##'   (`"effective_draws"`), the tail shape of the importance weights against
##'   its limit of 0.7 (`"pareto_k"`) and how far the estimates moved
##'   (`"moves"`). Two more panels come from what the fit recorded. For a
##'   variance held at zero whose slope was tested, `"held_variances"` shows
##'   the slope of the log-likelihood in that variance at zero with its
##'   interval of 2.58 standard errors each side, the test of the fit, so that
##'   the likelihood clearly rises as the variance leaves zero when the whole
##'   interval is above zero, and clearly falls when it is below. The panel
##'   `"screen"` shows the screen of short scales of the check for other
##'   maxima. On a log scale of distance it gives, at each short scale, the
##'   log-likelihood estimated on the fitting draws above that of the
##'   estimate, which repeats the fit from the best such scale when one comes
##'   within 0.25 of it (the dashed line), and the exact lower bound on the
##'   likelihood there, with the likelihood with independent location effects,
##'   the limit as the scale goes to zero, as a dotted line. An exact bound
##'   clearly above the estimate shows that the estimate is not the maximum
##'   likelihood estimate. For a fit with the spatial variance at zero and a
##'   nugget, the panel adds the exact bound over longer scales for each share
##'   of the variance in a spatial effect. With the extra check `"radius"`,
##'   the panel `"radius"` shows the error of each resampled sample of
##'   certification draws against the
##'   radius the fit would use for it, the larger of its own and a
##'   replicate's, in standard errors, with the line where they are equal. Its
##'   subtitle gives the predicted rate at which the radius is exceeded, as in
##'   the row of the check, and the share of samples whose error exceeds the
##'   radius, as in [summary.RiskMap_check()], from many more pairs of samples
##'   than the panel shows. With
##'   the extra check `"profile"`, each profiled parameter has a panel of its
##'   profile likelihood, named `"profile_"` followed by the parameter, such
##'   as `"profile_phi"`.
##' @param x An object of class `"RiskMap_check"`.
##' @param which Names of the panels to draw, among `"certificate"`,
##'   `"effective_draws"`, `"pareto_k"`, `"moves"`, `"held_variances"`,
##'   `"screen"`, `"radius"`, the panels of the profiles, such as
##'   `"profile_phi"`, and `"profile"` for every profile.
##'   The default is every panel that has something to show. A chosen panel
##'   with nothing to show gives a message saying why, and which extra check
##'   of [check_fit()] would compute it. If no panel has anything to show,
##'   `plot()` stops with those messages.
##' @param ... other parameters.
##' @return Invisibly returns the list of panels drawn, each a `ggplot`
##'   object.
##' @seealso [check_fit()], [summary.RiskMap_check()]
##' @examples
##' \donttest{
##' data(liberia)
##' liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
##' fit <- glgpm(
##'   npos ~ log_elevation + gp(kappa = 0.5), denominator = ntest,
##'   data = liberia, model_crs = 32629, family = "binomial", messages = FALSE
##' )
##' checked <- check_fit(fit, extra = "all")
##' plot(checked)
##' plot(checked, which = c("certificate", "screen", "radius", "profile_phi"))
##' }
##' @method plot RiskMap_check
##' @export
plot.RiskMap_check <- function(x, which = NULL, ...) {
  check_class(x, "RiskMap_check")
  panel_names <- c("certificate", "effective_draws", "pareto_k", "moves",
                   "held_variances", "screen", "radius")
  if (!is.null(which)) {
    if (!is.character(which) || !length(which) || anyNA(which)) {
      abort_input("{.arg which} must be a character vector, not {.obj_type_friendly {which}}.", arg = "which")
    }
    unknown <- which[!which %in% c(panel_names, "profile") &
                       !startsWith(which, "profile_")]
    if (length(unknown)) {
      abort_input(c("{.arg which} has unknown panel{?s} {.val {unknown}}.", "i" = "The panels are {.val {panel_names}}, {.val profile_} followed by a parameter, such as {.val profile_phi}, and {.val profile} for every profile."), arg = "which")
    }
  }
  panels <- draw_check_panels(x)
  profiles <- grep("^profile_", names(panels), value = TRUE)
  requested <- unique(unlist(lapply(which %||% c(panel_names, "profile"),
                                    function(name) {
    if (identical(name, "profile") && length(profiles)) profiles else name
  })))
  empty <- requested[vapply(panels[requested], is.null, logical(1L))]
  reasons <- explain_missing_panels(x, empty)
  panels <- Filter(Negate(is.null), panels[requested])
  if (!length(panels)) {
    explanation <- paste(reasons, collapse = " ")
    abort_input("{explanation}", arg = "which")
  }
  # Without a choice of panels, those with nothing to show are left out quietly.
  if (!is.null(which)) for (reason in reasons) riskmap_inform("{reason}")
  print(patchwork::wrap_plots(panels, ncol = min(2L, length(panels))))
  invisible(panels)
}

##' Every panel of the checks, NULL where there is nothing to show
##' @noRd
draw_check_panels <- function(x) {
  stages <- if (!is.null(x$stages)) list(
    certificate = plot_stage_panel(
      x$stages, "certificate", "Bound on the Monte Carlo error",
      "Standard errors (log scale)", limit = x$tolerance, log_scale = TRUE
    ),
    effective_draws = plot_stage_panel(
      x$stages, "effective_draws", "Effective certification draws",
      "Draws (log scale)", log_scale = TRUE
    ),
    pareto_k = plot_stage_panel(
      x$stages, "pareto_k", "Tail shape of the importance weights",
      "Pareto shape", limit = 0.7
    ),
    moves = plot_estimate_moves(x$moves)
  )
  profile <- x$extras$profile$data
  profiled <- names(profile$marginal)
  profiles <- lapply(profiled, function(name) plot_marginal_profile(profile, name))
  names(profiles) <- if (length(profiled)) paste0("profile_", profiled)
  c(list(certificate = stages$certificate,
         effective_draws = stages$effective_draws,
         pareto_k = stages$pareto_k, moves = stages$moves,
         held_variances = plot_held_variance_panel(
           x$results$Variances$data[["Slope of the likelihood at zero"]]
         ),
         screen = plot_screen_panel(x$results[["Other maxima"]]$data),
         radius = plot_radius_panel(x$extras$radius$data)),
    profiles)
}

##' Why each of the named panels has nothing to show
##'
##' A panel of an extra check that did not run names the extra to run.
##' @noRd
explain_missing_panels <- function(x, names) {
  gaussian <- is.null(x$stages)
  profiled <- names(x$extras$profile$data$marginal)
  reasons <- vapply(names, function(name) {
    if (identical(name, "profile") || startsWith(name, "profile_")) {
      if ("profile" %in% names(x$not_computed)) {
        return(x$not_computed[["profile"]])
      }
      if (!length(profiled)) {
        return(paste0("The panel \"", name, "\" needs the extra check ",
                      "\"profile\", check_fit(fit, extra = \"profile\")."))
      }
      return(paste0("The panel \"", name, "\" has no profile, since the ",
                    "profiles are of ", join_phrases(profiled), "."))
    }
    if (identical(name, "radius")) {
      if ("radius" %in% names(x$not_computed)) {
        return(x$not_computed[["radius"]])
      }
      if (!gaussian) {
        return(paste0("The panel \"radius\" needs the extra check ",
                      "\"radius\", check_fit(fit, extra = \"radius\")."))
      }
    }
    if (name %in% c("held_variances", "screen")) {
      return(explain_missing_fit_panel(x, name))
    }
    if (gaussian) {
      return(paste0("A Gaussian fit has no Monte Carlo error, so it has no ",
                    "stages to plot."))
    }
    switch(
      name,
      certificate = "No stage of the fit has a certificate to plot.",
      effective_draws = ,
      pareto_k = "No stage of the fit reached the certification draws.",
      moves = paste0("The fit has one stage, or no standard errors, so there ",
                     "are no changes in the estimates to plot.")
    )
  }, character(1L))
  unique(unname(reasons))
}

##' Why the panel of the held variances or of the screen has nothing to show
##'
##' Both draw only what the fit recorded, so the reason is about the fit, not
##' an extra check.
##' @noRd
explain_missing_fit_panel <- function(x, name) {
  held_panel <- identical(name, "held_variances")
  if (is.null(x$stages)) {
    return(paste0("A Gaussian fit is fitted exactly, so it has no ",
                  if (held_panel) "test of a variance held at zero" else
                    "screen of short scales", " to plot."))
  }
  if ("Maximum of the likelihood" %in% x$checks$check) {
    return(paste0("The likelihood has no maximum, so there is no ",
                  if (held_panel) "test of a held variance" else
                    "screen of short scales", " to plot."))
  }
  if (held_panel) {
    held <- x$results$Variances$data[["Held at zero"]]
    if (is.null(held) || identical(held, "none")) {
      return(paste0("No variance was held at zero, so there is no slope at ",
                    "zero to plot."))
    }
    if (any(grepl("sigma2", held, fixed = TRUE))) {
      return(paste0("The ", join_phrases(held), " held at zero ",
                    if (length(held) > 1L) "were" else "was", " not tested, ",
                    "since with the spatial variance at zero the slope ",
                    "depends on the scale, which is then unidentified."))
    }
    return(paste0("The fit recorded no test of the slope at zero for the ",
                  join_phrases(held), "."))
  }
  data <- x$results[["Other maxima"]]$data
  gate <- data[["Screen of short scales"]]
  if (is.null(gate)) {
    return("The fit recorded no screen of short scales.")
  }
  if (startsWith(gate, "unavailable")) {
    return(paste0("The screen of short scales did not run on this fit, ",
                  "because of its ",
                  sub("^unavailable \\((.*)\\)$", "\\1", gate), "."))
  }
  if (identical(gate, "failed")) {
    return(paste0("The screen of short scales did not run, because the ",
                  "likelihood with independent location effects could not be ",
                  "computed."))
  }
  largest <- data[["Largest exact lower bound, above the estimate"]]
  if (length(largest)) {
    return(paste0("The fit kept only the largest exact lower bound of the ",
                  "screen of short scales (", signif(largest, 2), " above the ",
                  "estimate), as fits did before the bound at each scale was ",
                  "recorded, so there is nothing to plot. A refit records it."))
  }
  "The fit recorded no values of the screen of short scales to plot."
}

##' The slope of the log-likelihood at each variance held at zero
##'
##' From `slopes`, the test of `assess_held_variances()` as the check of the
##' variances keeps it, with an interval of `multiplier` standard errors each
##' side, the multiplier of that test. The verdict of each slope is the one the
##' fit recorded, written next to it so that it is not told by colour alone.
##' NULL without a tested variance.
##' @noRd
plot_held_variance_panel <- function(slopes, multiplier = 2.58) {
  if (is.null(slopes) || !nrow(slopes)) return(NULL)
  colours <- profile_palette()
  slopes$lower <- slopes$slope - multiplier * slopes$se
  slopes$upper <- slopes$slope + multiplier * slopes$se
  slopes$verdict <- ifelse(slopes$rising, "rises from zero",
                           ifelse(slopes$falling, "falls from zero",
                                  "undecided"))
  slopes$name <- factor(slopes$name, rev(unique(slopes$name)))
  verdict_colours <- c("rises from zero" = colours$cutoff,
                       "falls from zero" = colours$profile,
                       "undecided" = colours$wald)
  ggplot2::ggplot(slopes, ggplot2::aes(.data$slope, .data$name,
                                       colour = .data$verdict)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = colours$ridge) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lower,
                                         xmax = .data$upper), linewidth = 0.9) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_text(ggplot2::aes(.data$upper, label = .data$verdict),
                       hjust = -0.12, size = 3.2, colour = "#3b4350") +
    ggplot2::scale_colour_manual(values = verdict_colours, guide = "none") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.3))
    ) +
    ggplot2::labs(
      title = "Slope of the likelihood at a variance held at zero",
      subtitle = paste(strwrap(paste0(
        "The interval is ", multiplier, " standard errors each side. Above ",
        "zero, the likelihood rises as the variance leaves zero, so zero is ",
        "not a maximum. Below zero, it falls, so zero is one."
      ), width = 62L), collapse = "\n"),
      x = "Slope of the log-likelihood in the variance at zero", y = NULL
    ) +
    profile_theme()
}

##' The screen of short scales and its exact bounds, against the estimate
##'
##' From the data of the check of other maxima (`tabulate_short_scales()`,
##' `tabulate_spatial_shares()`), on a log scale of distance, with the
##' log-likelihood above that of the estimate. The dashed line at `margin`
##' below the estimate is the gate of the screen, the margin of
##' `make_short_scale_screen()`, and the dotted line the limit with
##' independent location effects. The bounds by spatial share are labelled at
##' their right end. The vertical axis is cut a unit below the lowest of the
##' screened values, the gate and the limit, since a bound far below them,
##' which the bound at a scale near the estimate often is, says nothing
##' about other maxima. NULL with nothing to show.
##' @noRd
plot_screen_panel <- function(data, margin = 0.25) {
  short <- data[["Short scales, log-likelihood above the estimate"]]
  shares <- data[["Exact bound by scale and spatial share, above the estimate"]]
  if (is.null(short) && is.null(shares)) return(NULL)
  colours <- profile_palette()
  table <- if (!is.null(short)) short else shares
  units <- gsub("_", " ", sub("^scale_", "", names(table)[1L]))
  series <- c("Screened, on the fitting draws" = "screened",
              "Exact lower bound" = "exact_bound")
  series_colours <- stats::setNames(c(colours$profile, colours$ridge),
                                    names(series))
  lines <- reshape_screen_columns(short, series)
  share_lines <- if (!is.null(shares)) {
    reshape_screen_columns(shares, stats::setNames(
      names(shares)[-1L], paste(names(shares)[-1L], "spatial")
    ))
  }
  estimate <- unlist(data[startsWith(names(data), "Scale of the estimate")])
  limit <- data[["Limit with independent location effects, above the estimate"]]
  panel <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "#b7bec8")
  left <- min(table[[1L]], shares[[1L]])
  if (!is.null(short$screened)) {
    panel <- panel +
      ggplot2::geom_hline(yintercept = -margin, linetype = "dashed",
                          colour = colours$cutoff) +
      ggplot2::annotate("text", x = left, y = -margin, label = "gate",
                        hjust = 0, vjust = 1.5, size = 2.9,
                        colour = colours$cutoff)
  }
  if (length(limit) && is.finite(limit)) {
    panel <- panel +
      ggplot2::geom_hline(yintercept = limit, linetype = "dotted",
                          colour = colours$wald, linewidth = 0.8)
  }
  if (!is.null(share_lines)) {
    levels <- levels(share_lines$series)
    ends <- share_lines[!duplicated(share_lines$series, fromLast = TRUE), ]
    ramp <- grDevices::colorRampPalette(c("#7fb2db", "#08306b"))(length(levels))
    panel <- panel +
      ggplot2::geom_line(data = share_lines,
                         ggplot2::aes(.data$scale, .data$value,
                                      colour = .data$series),
                         linewidth = 0.7) +
      ggplot2::geom_text(data = ends,
                         ggplot2::aes(.data$scale, .data$value,
                                      label = .data$series),
                         hjust = -0.1, size = 2.9, colour = "#3b4350")
    series_colours <- c(series_colours, stats::setNames(ramp, levels))
  }
  if (!is.null(lines)) {
    panel <- panel +
      ggplot2::geom_line(data = lines,
                         ggplot2::aes(.data$scale, .data$value,
                                      colour = .data$series),
                         linewidth = 0.8) +
      ggplot2::geom_point(data = lines,
                          ggplot2::aes(.data$scale, .data$value,
                                       colour = .data$series), size = 2.2)
  }
  if (length(estimate) && is.finite(estimate)) {
    panel <- panel +
      ggplot2::annotate("point", x = estimate, y = 0, size = 3, shape = 21,
                        fill = "white", colour = colours$ridge, stroke = 1) +
      ggplot2::annotate("text", x = estimate, y = 0, label = "estimate",
                        vjust = -1, size = 3.2, colour = "#3b4350")
  }
  if (length(limit) && is.finite(limit)) {
    panel <- panel +
      ggplot2::annotate("text", x = max(table[[1L]], shares[[1L]], estimate),
                        y = limit, label = "independent location effects",
                        hjust = 1, vjust = 1.5, size = 2.9, colour = "#3b4350")
  }
  shown <- c(lines$value, share_lines$value, limit, 0)
  focus <- c(0, -margin, limit, short$screened,
             if (is.null(short$screened)) share_lines$value)
  focus <- focus[is.finite(focus)]
  shown <- shown[is.finite(shown)]
  panel +
    ggplot2::scale_x_log10(labels = plain_number_labels(),
                           expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
    ggplot2::scale_colour_manual(
      values = series_colours,
      breaks = intersect(names(series), as.character(lines$series))
    ) +
    ggplot2::labs(
      title = "Screen of short scales",
      subtitle = paste(strwrap(describe_screen_outcome(data, margin),
                               width = 62L), collapse = "\n"),
      x = paste0("Spatial scale (", units, ", log scale)"),
      y = "Log-likelihood above the estimate", colour = NULL
    ) +
    ggplot2::coord_cartesian(ylim = c(min(focus) - 1, max(shown) + 0.12 *
                                        (max(shown) - min(focus) + 1))) +
    profile_theme() +
    ggplot2::theme(legend.position = "bottom")
}

##' Columns of a table of the screen as lines to draw
##'
##' `columns` names the columns to take, by the label of their series. One row
##' per finite value, with the `scale` from the first column of `table` and
##' the `series` as a factor in the order of `columns`. NULL without values.
##' @noRd
reshape_screen_columns <- function(table, columns) {
  columns <- columns[columns %in% names(table)]
  if (is.null(table) || !length(columns)) return(NULL)
  lines <- do.call(rbind, lapply(names(columns), function(label) {
    data.frame(scale = table[[1L]], value = table[[columns[[label]]]],
               series = label)
  }))
  lines <- lines[is.finite(lines$value), ]
  if (!nrow(lines)) return(NULL)
  lines$series <- factor(lines$series, names(columns))
  lines
}

##' What the screen of short scales and its exact alarm found, for the panel
##'
##' A fit made before the bound at each scale was recorded has only the
##' largest, which the subtitle gives.
##' @noRd
describe_screen_outcome <- function(data, margin) {
  gate <- data[["Screen of short scales"]] %||% ""
  largest <- data[["Largest exact lower bound, above the estimate"]]
  short <- data[["Short scales, log-likelihood above the estimate"]]
  paste0(
    switch(gate,
           closed = paste0("No short scale came within ", margin, " of the ",
                           "estimate on the fitting draws, so the fit was not ",
                           "repeated. "),
           opened = paste0("A short scale came within ", margin, " of the ",
                           "estimate on the fitting draws, so the fit was ",
                           "repeated from it. "),
           ""),
    if (isTRUE(data[["Exact alarm"]])) {
      paste0("An exact lower bound is above the estimate, which is therefore ",
             "not the maximum likelihood estimate.")
    } else {
      "No exact lower bound is clearly above the estimate."
    },
    if (length(largest) && is.null(short$exact_bound)) paste0(
      " The fit kept only the largest exact lower bound (", signif(largest, 2),
      "), so the bound at each scale is not shown."
    )
  )
}

##' One quantity by stage, with a dashed line at its limit if it has one
##'
##' NULL when no stage has a finite value.
##' @noRd
plot_stage_panel <- function(stages, column, title, axis_title, limit = NULL,
                             log_scale = FALSE) {
  shown <- stages[is.finite(stages[[column]]), ]
  if (!nrow(shown)) return(NULL)
  colours <- profile_palette()
  panel <- ggplot2::ggplot(shown, ggplot2::aes(.data$stage, .data[[column]]))
  if (!is.null(limit)) {
    panel <- panel + ggplot2::geom_hline(yintercept = limit, linetype = "dashed",
                                         colour = colours$cutoff)
  }
  if (nrow(shown) > 1L) {
    panel <- panel + ggplot2::geom_line(colour = colours$profile)
  }
  panel <- panel +
    ggplot2::geom_point(colour = colours$profile, size = 2) +
    ggplot2::scale_x_continuous(breaks = stages$stage) +
    ggplot2::labs(title = title, x = "Stage", y = axis_title) +
    profile_theme()
  if (log_scale) {
    panel <- panel + ggplot2::scale_y_log10(breaks = scales::breaks_log(n = 5L),
                                            labels = plain_number_labels())
  }
  panel
}

##' The change in each estimate between stages, with a line through the largest
##'
##' NULL without any changes to show.
##' @noRd
plot_estimate_moves <- function(moves) {
  if (is.null(moves)) return(NULL)
  colours <- profile_palette()
  largest <- stats::aggregate(move ~ stage, moves, max)
  panel <- ggplot2::ggplot(moves, ggplot2::aes(.data$stage, .data$move)) +
    ggplot2::geom_point(colour = colours$wald, size = 1.2)
  if (nrow(largest) > 1L) {
    panel <- panel + ggplot2::geom_line(data = largest, colour = colours$profile)
  }
  panel +
    ggplot2::geom_point(data = largest, colour = colours$profile, size = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(moves$stage))) +
    ggplot2::labs(title = "Change in the estimates from the stage before",
                  subtitle = "The line follows the largest change.",
                  x = "Stage", y = "Standard errors") +
    profile_theme()
}

##' The error of resampled samples against the radius the fit uses
##'
##' From the data of the extra check `"radius"`, on log scales. The radius is
##' the larger of the sample's own and a replicate's. Samples above the
##' dashed line, where the two are equal, have an error larger than the
##' radius. The subtitle starts with the predicted rate at which the radius
##' is exceeded, the first number of the row, and gives the share above the
##' line as `own`, the figure of `summary()` from 40,000 pairs of samples,
##' rather than the share of the points shown. NULL without samples to show.
##' @noRd
plot_radius_panel <- function(check) {
  samples <- check$samples
  if (is.null(samples)) return(NULL)
  samples <- samples[is.finite(samples$error) & is.finite(samples$radius) &
                       samples$error > 0 & samples$radius > 0, , drop = FALSE]
  if (!nrow(samples)) return(NULL)
  colours <- profile_palette()
  percent <- format_radius_percent
  subtitle <- paste0(
    "Predicted to be exceeded in ", percent(check$level), " of samples (",
    percent(check$level_interval[1L]), " to ",
    percent(check$level_interval[2L]), "), against 5%. The error of the ",
    "certification draws alone exceeds it in ", percent(check$own),
    ", above the dashed line."
  )
  ggplot2::ggplot(samples, ggplot2::aes(.data$radius, .data$error)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = colours$cutoff) +
    ggplot2::geom_point(colour = colours$profile, size = 1, alpha = 0.35) +
    ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 5L),
                           labels = plain_number_labels()) +
    ggplot2::scale_y_log10(breaks = scales::breaks_log(n = 5L),
                           labels = plain_number_labels()) +
    ggplot2::labs(
      title = "Error and radius of resampled certification draws",
      subtitle = paste(strwrap(subtitle, width = 62L), collapse = "\n"),
      x = "Radius of the bound (standard errors, log scale)",
      y = "Error (standard errors, log scale)"
    ) +
    profile_theme()
}
