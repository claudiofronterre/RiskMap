##' The numbers behind the finding of a default check
##'
##' A named list for `summary()` to print, whose elements are single values or
##' tables. Elements that the fit did not record are left out. `check` is the
##' name of the check, as in its row, and `mcml` the list of
##' `summarise_mcml_fit()`.
##' @noRd
collect_check_data <- function(check, object, mcml) {
  optimizer <- attr(object, "optimizer")
  history <- object$mcml_history %||% optimizer$mcml_history
  last <- if (length(history)) history[[length(history)]]
  units <- object$distance_units %||% "distance units"
  screen <- optimizer$mcml_maxima_screen
  spatial_held <- any(grepl("sigma2", optimizer$mcml_boundary, fixed = TRUE))
  data <- switch(
    check,
    "Monte Carlo accuracy" = list(
      "Certificate, in standard errors" = mcml$certificate,
      "Newton correction on the certification draws" = last$certificate_distance,
      "95% radius of its Monte Carlo error" = last$certificate_radius,
      "Tolerance" = mcml$tolerance,
      "Outcome" = mcml$outcome
    ),
    "Importance weights" = list(
      "Certification draws" = last$validation_budget,
      "Effective draws" = mcml$certification_bank$effective_draws,
      "Largest weight, as a share of the total" =
        last$second_validation$maximum_weight,
      "Pareto shape" = mcml$certification_bank$pareto_k,
      "Limit of the Pareto shape" = 0.7
    ),
    "Other maxima" = list(
      "Exact alarm" = isTRUE(optimizer$mcml_maxima_alarm),
      "Higher maximum found by the check" = isTRUE(optimizer$mcml_better_maximum),
      "Separate maxima about equally high" =
        isTRUE(optimizer$mcml_ambiguous_maxima),
      "Screen of short scales" = describe_screen_gate(screen),
      "Scale of the estimate" = if (!spatial_held) stats::coef(object)$phi,
      "Short scales, log-likelihood above the estimate" =
        tabulate_short_scales(screen, units),
      "Largest exact lower bound, above the estimate" =
        screen$alarm$statistic - screen$alarm$fitted_comparison,
      "Limit with independent location effects, above the estimate" =
        if (!isTRUE(optimizer$mcml_pc_prior)) {
          screen$independent$log_likelihood - screen$alarm$fitted_comparison
        },
      "Exact bound by scale and spatial share, above the estimate" =
        tabulate_spatial_shares(screen, units),
      "Maxima from the dispersed starts" =
        tabulate_maxima_check(optimizer$mcml_maxima_check, units)
    ),
    "Variances" = list(
      "Held at zero" = if (length(optimizer$mcml_boundary)) {
        optimizer$mcml_boundary
      } else {
        "none"
      },
      "Slope of the likelihood at zero" = optimizer$mcml_held_variances[
        intersect(c("name", "slope", "se", "rising", "falling"),
                  names(optimizer$mcml_held_variances))
      ]
    ),
    "Spatial scale" = stats::setNames(list(
      stats::coef(object)$phi, optimizer$mcml_spacing,
      isTRUE(optimizer$mcml_small_scale), isTRUE(optimizer$mcml_near_limit),
      isTRUE(optimizer$mcml_numerically_independent),
      isTRUE(optimizer$mcml_identifiability_failure)
    ), c(paste0("Scale (", units, ")"),
         paste0("Median distance between neighbouring locations (", units, ")"),
         "Shorter than that distance",
         "Likelihood near the limit of independent location effects",
         "Numerically independent location effects",
         "Parameters not identifiable")),
    "Importance distribution" = list(
      "Distribution" = optimizer$mcml_proposal,
      "Restarted from the starting values" = isTRUE(optimizer$mcml_start_replaced),
      "Stages that used the Laplace approximation" =
        if (length(optimizer$mcml_proposal_fallbacks)) {
          optimizer$mcml_proposal_fallbacks
        } else {
          "none"
        },
      "Parameter values where it could not be built" =
        optimizer$mcml_proposal_errors$count %||% 0L,
      "First error" = optimizer$mcml_proposal_errors$first
    ),
    "Latent mode search" = list(
      "Searches" = optimizer$mcml_mode_searches$searches,
      "Stopped at the iteration limit" = optimizer$mcml_mode_searches$at_limit,
      "Started again from zero after failing from the mode found before" =
        optimizer$mcml_mode_searches$restarted
    ),
    "Link" = list("Link" = object$link_function$link %||% "a function"),
    "Maximum of the likelihood" = list(
      "Every observation" = optimizer$mcml_separation$everywhere,
      "Covariates that separate the data" = optimizer$mcml_separation$coefficients,
      "Coefficients at the bound of the search" =
        optimizer$mcml_coefficients_at_bound
    )
  )
  scale_name <- match("Scale of the estimate", names(data))
  if (!is.na(scale_name)) {
    names(data)[scale_name] <- paste0("Scale of the estimate (", units, ")")
  }
  Filter(function(value) length(value) > 0L, data)
}

##' What the screen of short scales did, in a few words
##' @noRd
describe_screen_gate <- function(screen) {
  if (is.null(screen$gate)) return(NULL)
  paste0(screen$gate, if (!is.null(screen$reason)) paste0(" (", screen$reason, ")"))
}

##' The name of a column of spatial scales in the distance `units`
##' @noRd
name_scale_column <- function(units) {
  paste0("scale_", gsub("[^A-Za-z]+", "_", units))
}

##' The log-likelihood at the short scales of the screen, above the estimate
##'
##' At each scale of the screen of short scales (`screen_short_scales()`),
##' `screened` is the log-likelihood estimated on the fitting draws at the
##' point from the maximum with independent location effects, or after a step
##' in the spatial variance where that was higher, less the fitted
##' log-likelihood on the same draws. `exact_bound` is the exact lower bound
##' at that point (`judge_exact_alarm()`) less the log-likelihood of the
##' estimate on independent draws. A column is left out when the fit did not
##' record it, as when the screen did not run because the estimate is at a
##' short scale, and NULL when neither was recorded. The scale is in `units`.
##' @noRd
tabulate_short_scales <- function(screen, units) {
  screened <- screen$screen
  bounds <- screen$alarm$bounds
  reference <- screen$alarm$fitted_comparison
  if (!is.finite(reference %||% NA_real_)) bounds <- NULL
  log_phi <- if (!is.null(screened)) screened$log_phi else bounds$log_phi
  if (!length(log_phi) || (is.null(screened) && is.null(bounds))) return(NULL)
  table <- data.frame(scale = exp(log_phi))
  if (!is.null(screened)) table$screened <- screened$gap
  if (!is.null(bounds)) {
    table$exact_bound <- bounds$bound[match(log_phi, bounds$log_phi)] -
      reference
  }
  names(table)[1L] <- name_scale_column(units)
  table
}

##' The exact bound over scales and spatial shares, above the estimate
##'
##' For a fit with the spatial variance at zero and a nugget, the exact lower
##' bound on the likelihood with a share of the total variance in a spatial
##' effect (`judge_spatial_share()`), less the maximum with independent
##' location effects, which is then the estimate. One row per scale in
##' `units` and one column per share. NULL where the fit did not record them.
##' @noRd
tabulate_spatial_shares <- function(screen, units) {
  shares <- screen$alarm$shares
  reference <- screen$independent$log_likelihood
  if (is.null(shares) || !is.finite(reference %||% NA_real_)) return(NULL)
  log_phi <- sort(unique(shares$log_phi))
  table <- data.frame(scale = exp(log_phi))
  for (share in sort(unique(shares$share))) {
    rows <- shares[shares$share == share, ]
    table[[paste0(100 * share, "%")]] <-
      rows$bound[match(log_phi, rows$log_phi)] - reference
  }
  names(table)[1L] <- name_scale_column(units)
  table
}

##' The maxima of the check for other maxima, with the scale in `units`
##'
##' The log-likelihood on independent draws and its standard error are left
##' out when no maximum was judged on them.
##' @noRd
tabulate_maxima_check <- function(check, units) {
  if (is.null(check)) return(NULL)
  table <- data.frame(start = check$start, scale = exp(check$log_phi),
                      loglik = check$log_likelihood,
                      independent_loglik = check$independent_log_likelihood,
                      independent_se = check$independent_se,
                      chosen = check$chosen)
  names(table)[2L] <- name_scale_column(units)
  judged <- is.finite(table$independent_loglik)
  if (!any(judged)) table[c("independent_loglik", "independent_se")] <- NULL
  table
}

##' @title Summarise the Checks on a Fitted Model
##' @description Summary method for objects of class `"RiskMap_check"` from
##'   [check_fit()]. Where `print()` gives the verdict, `summary()` gives the
##'   evidence, the finding of each check in full with the numbers behind it,
##'   the certificate and the certification draws of each stage with the
##'   largest change in the estimates from the stage before, and the full
##'   results of the extra checks that ran. For the check of other maxima the
##'   numbers include, at each short scale of the screen, the log-likelihood
##'   estimated on the fitting draws and its exact lower bound, both above the
##'   estimate, and for a fit with the spatial variance at zero and a nugget
##'   the exact bound by scale and spatial share.
##' @param object An object of class `"RiskMap_check"`.
##' @param ... other parameters.
##' @return An object of class `"summary.RiskMap_check"`, a list with the
##'   `checks` table, the `results` of the default and extra checks in the
##'   order of the table, the `stages` table with the column `largest_move`,
##'   and the `draws`, `fit_seconds` and `not_computed` of `object`.
##' @seealso [check_fit()], [plot.RiskMap_check()]
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
##' summary(check_fit(fit))
##' @method summary RiskMap_check
##' @export
summary.RiskMap_check <- function(object, ...) {
  stopifnot(inherits(object, "RiskMap_check"))
  results <- c(object$results, object$extras)
  checks <- vapply(results, function(result) result$row$check, character(1L))
  results <- results[match(object$checks$check, checks)]
  structure(
    list(checks = object$checks, results = results,
         extras = names(object$extras),
         stages = add_largest_moves(object$stages, object$moves),
         draws = object$draws, fit_seconds = object$fit_seconds,
         not_computed = object$not_computed),
    class = "summary.RiskMap_check"
  )
}

##' The stages table with the largest change in an estimate at each stage
##'
##' `largest_move` is in standard errors, `NA` at the first stage and where no
##' estimate has a standard error.
##' @noRd
add_largest_moves <- function(stages, moves) {
  if (is.null(stages)) return(NULL)
  stages$largest_move <- NA_real_
  if (!is.null(moves)) {
    largest <- stats::aggregate(move ~ stage, moves, max)
    stages$largest_move[match(largest$stage, stages$stage)] <- largest$move
  }
  stages
}

##' @title Print the Summary of the Checks on a Fitted Model
##' @description Print method for objects of class `"summary.RiskMap_check"`
##'   from [summary.RiskMap_check()].
##' @param x An object of class `"summary.RiskMap_check"`.
##' @param ... other parameters.
##' @return Invisibly returns `x`.
##' @method print summary.RiskMap_check
##' @export
print.summary.RiskMap_check <- function(x, ...) {
  stopifnot(inherits(x, "summary.RiskMap_check"))
  cat("Evidence for the checks on this fit\n")
  timing <- if (!is.null(x$fit_seconds)) {
    paste0("The fit took ", format_check_time(x$fit_seconds), ".")
  }
  if (!is.null(x$draws) || !is.null(timing)) {
    write_wrapped(paste(c(x$draws, timing), collapse = " "))
  }
  status_width <- max(nchar(x$checks$status))
  for (result in x$results) {
    cat("\n", formatC(result$row$status, width = -status_width), "  ",
        result$row$check, "\n", sep = "")
    write_indented(strwrap(result$row$finding,
                           width = getOption("width", 80L) - 4L))
    if (inherits(result$data, "RiskMap_profile")) {
      write_profile_evidence(result$data)
    } else if (inherits(result$data, "RiskMap_radius_check")) {
      write_radius_evidence(result$data)
    } else {
      write_check_data(result$data)
    }
  }
  if (!is.null(x$stages)) {
    cat("\nCertificate and certification draws by stage\n")
    write_indented(utils::capture.output(
      print(x$stages, digits = 4L, row.names = FALSE)
    ))
  }
  if (length(x$not_computed)) {
    cat("\n")
    write_wrapped(x$not_computed)
  }
  cat("\n")
  write_wrapped(paste0("These checks concern the computation, not whether the ",
                       "model suits the data."))
  invisible(x)
}

##' Lines of text indented by two spaces
##' @noRd
write_indented <- function(lines) {
  cat(paste0("  ", lines), sep = "\n")
  invisible(NULL)
}

##' The data of a check, single values in two columns and tables under a title
##' @noRd
write_check_data <- function(data) {
  if (!length(data)) {
    write_indented("Nothing recorded beyond the finding.")
    return(invisible(NULL))
  }
  tables <- vapply(data, is.data.frame, logical(1L))
  values <- vapply(data[!tables], format_check_value, character(1L))
  if (length(values)) {
    label_width <- max(nchar(names(values)))
    write_indented(paste0(formatC(names(values), width = -label_width), "  ",
                          values))
  }
  for (name in names(data)[tables]) {
    write_indented(name)
    write_indented(paste0("  ", utils::capture.output(
      print(data[[name]], digits = 4L, row.names = FALSE)
    )))
  }
  invisible(NULL)
}

##' One value of the data of a check as text
##' @noRd
format_check_value <- function(value) {
  if (is.logical(value)) return(paste(ifelse(value, "yes", "no"), collapse = ", "))
  if (is.numeric(value)) {
    return(paste(format(value, digits = 4L, big.mark = ",", trim = TRUE),
                 collapse = ", "))
  }
  paste(value, collapse = ", ")
}

##' The full results of the profile check
##'
##' The intervals as `print()` of the profile gives them, then each point of
##' each marginal profile, with the parameter on its natural scale.
##' @noRd
write_profile_evidence <- function(profile) {
  write_indented(utils::capture.output(print(profile)))
  for (name in names(profile$marginal)) {
    points <- profile$marginal[[name]]
    points <- points[order(points$value), ]
    write_indented(paste0("Points of the profile of ",
                          describe_profile_parameter(name)))
    write_indented(paste0("  ", utils::capture.output(print(
      data.frame(value = profile_natural(name, points$value),
                 log_likelihood = points$loglik, deviance = points$deviance,
                 status = points$status),
      digits = 4L, row.names = FALSE
    ))))
  }
  invisible(NULL)
}

##' The full results of the check of the radius
##'
##' The predicted rates, the bound so widened and the settings, then the
##' error and radius of the resampled samples at a few quantiles. Rates are
##' percentages of samples of certification draws, to two figures as in the
##' row and the panel. The rate with the radius of a sample's own draws
##' alone is a diagnostic, since the fit uses the larger of two radii.
##' @noRd
write_radius_evidence <- function(check) {
  percent <- format_radius_percent
  settings <- check$settings
  probabilities <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  quantiles <- function(values) {
    stats::quantile(values, probabilities, names = FALSE)
  }
  write_check_data(list(
    "Samples whose error exceeds the radius, predicted" = percent(check$level),
    "Its 95% interval, without the error of the chains" =
      paste(percent(check$level_interval), collapse = " to "),
    "Samples whose certification draws alone exceed the radius" =
      percent(check$own),
    "Diagnostic, the same with the radius of their own draws alone" =
      percent(check$own_alone),
    "Samples whose error exceeds the bound" = percent(check$exceedance),
    "Growth of the radius that brings the first rate to 5%" = check$factor,
    "Margin on that growth" = paste0(100 * check$margin, "%"),
    "Widening of the radius" = check$widening,
    "Bound, in standard errors" = check$certificate,
    "Bound with the radius so widened" = check$widened,
    "Tolerance" = check$tolerance,
    "Nugget in the samples" = describe_radius_forms(check),
    "Pareto shape of the fresh sample" = check$bank_pareto_k,
    "Median Pareto shape of 100 resampled samples" = check$pareto_k,
    "Split R-hat of the chains' log weight" = check$rhat,
    "Chains" = paste(settings$chains, "of", format(settings$rounds,
                                                   big.mark = ","), "rounds"),
    "Resampled samples" = paste(format(settings$banks, big.mark = ","), "of",
                                format(settings$draws, big.mark = ","),
                                "draws"),
    "Seed of the fresh sample" = as.character(settings$seed),
    "Time of the check" = format_check_time(check$seconds),
    "Error and radius of the resampled samples, in standard errors" =
      data.frame(quantile = paste0(100 * probabilities, "%"),
                 error = quantiles(check$samples$error),
                 radius = quantiles(check$samples$radius),
                 error_over_radius = quantiles(check$samples$error /
                                                 check$samples$radius))
  ))
}

##' The form of the nugget in the samples of the check of the radius
##' @noRd
describe_radius_forms <- function(check) {
  forms <- check$forms %||% "sampled"
  if (length(forms) > 1L) {
    return(paste0("sampled and integrated out, combined in pairs (both ",
                  "banks in ", format_radius_percent(100 * check$used[["both"]]),
                  " of pairs)"))
  }
  if (identical(forms, "integrated")) "integrated out" else "sampled"
}
