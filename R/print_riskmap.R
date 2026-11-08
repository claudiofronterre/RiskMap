##' Certified MCML part of a model summary
##'
##' Reads the certificate of a non-Gaussian fit from the attributes that
##' `glgpm()` writes, so that it is available after fitting, for example from a
##' saved model. NULL for Gaussian fits, which have no Monte Carlo error.
##' @return A list with the `certificate` (an approximate bound at `level`, in
##' statistical standard errors, on the Monte Carlo error of the estimate,
##' NA when there is none), the `tolerance`, whether the fit is `certified`,
##' its `outcome`, whether the repeated Newton correction of the last stage
##' was refused (`correction_refused`), whether its stage ended because the
##' likelihood rises as a held variance leaves zero (`variance_rising`), the
##' number of `stages`, the
##' `fitting_draws` and `validation_draws` of the last stage, the
##' `certification_bank`, the number of locations whose nugget was integrated
##' out of the likelihood at the last stage (`integrated_locations`, of
##' `locations`), the forms of the stages of an automatic fit
##' (`noise_forms`), the largest
##' change of a factor at the last fixed refit of the last stage
##' (`refit_change`), whether to suggest `integrate_noise = TRUE`
##' (`suggest_integration`, `suggests_integration()`) and the number of checks
##' of `check_fit()` that need a look (`checks_to_review`).
##' @noRd
summarise_mcml_fit <- function(object) {
  optimizer <- attr(object, "optimizer")
  if (is.null(optimizer$mcml_outcome)) return(NULL)
  history <- object$mcml_history %||% optimizer$mcml_history
  last <- if (length(history)) history[[length(history)]] else NULL
  mcml <- list(
    certificate = optimizer$mcml_certificate %||% NA_real_,
    level = 0.95,
    tolerance = attr(object, "control_mcml")$tolerance %||% NA_real_,
    certified = isTRUE(optimizer$mcml_certified),
    outcome = optimizer$mcml_outcome,
    correction_refused = isTRUE(last$residual_refused),
    variance_rising = isTRUE(last$held_variance_rising),
    stages = optimizer$mcml_stages %||% length(history),
    fitting_draws = optimizer$mcml_budget,
    validation_draws = optimizer$mcml_validation_budget,
    certification_bank = if (!is.null(last)) summarise_certification_bank(last),
    integrated_locations = optimizer$mcml_integrated_locations %||% 0L,
    locations = optimizer$mcml_locations,
    noise_forms = optimizer$mcml_noise_forms,
    refit_change = last$refit_change %||% NA_real_,
    suggest_integration = suggests_integration(object, optimizer)
  )
  mcml$checks_to_review <- sum(list_fit_checks(object, mcml)$status != "OK")
  mcml
}

##' Whether `check_fit()` suggests integrating the nugget out
##'
##' For a fit that sampled a nugget whose variance is estimated, or fixed, at
##' `threshold` or more. Where the nugget variance at the maximum is small the
##' likelihood is nearly flat along the nugget, and there the integrated
##' estimate of the nugget's score is the noisier of the two forms. In the
##' larger test of the integrated form it certified 5 runs against 12 of 21 for
##' the sampled form where that variance was below 0.1, 9 and 9 of 18 between
##' 0.1 and 0.3, and 68 against 28 of 75 above 0.3. A nugget estimated as zero
##' gets no suggestion, nor does a fit that chose the form itself
##' (`integrate_noise = "auto"`).
##' @noRd
suggests_integration <- function(object, optimizer, threshold = 0.1) {
  if (isTRUE(optimizer$mcml_integrate_noise) ||
      !is.null(optimizer$mcml_noise_forms)) {
    return(FALSE)
  }
  variance <- if (is.numeric(object$fix_tau2)) {
    object$fix_tau2
  } else if (!is.null(object$estimate$nu2)) {
    exp(object$estimate$nu2 + object$estimate$sigma2)
  }
  isTRUE(variance >= threshold)
}

##' One sentence on the Monte Carlo certificate of a fit
##'
##' `mcml` is the list of `summarise_mcml_fit()`.
##' @noRd
describe_mcml_certificate <- function(mcml) {
  bound <- signif(mcml$certificate, 2)
  level <- paste0(100 * mcml$level, "%")
  switch(
    mcml$outcome,
    "certified" = paste0(
      "Certified: the Monte Carlo error of the estimate is within ", bound,
      " statistical standard errors (an approximate ", level, " bound, ",
      "tolerance ", mcml$tolerance, ")."
    ),
    "weak identification" = paste0(
      "No certificate: the observed information is singular, so the Monte ",
      "Carlo precision of the estimate is not established."
    ),
    "numerically independent" = paste0(
      "No certificate: the fit is numerically one of independent location ",
      "effects, and the certificate needs an interior maximum."
    ),
    "no maximum" = paste0(
      "No certificate: the likelihood has no maximum, so there is no ",
      "maximum likelihood estimate whose Monte Carlo error could be bounded."
    ),
    if (isTRUE(mcml$variance_rising)) {
      paste0(
        "Not certified: the likelihood rises as a variance held at zero leaves ",
        "it, so the estimate is not a maximum. Refit with a larger ",
        "'max_samples' or 'max_stages' in set_control_mcml()."
      )
    } else if (mcml$correction_refused) {
      paste0(
        "Not certified: the Newton correction of the estimate did not ",
        "converge on the validation draws, so the approximate ", level,
        " bound of ", bound, " statistical standard errors may not hold. ",
        "Refit with a larger 'max_samples' in set_control_mcml()."
      )
    } else {
      paste0(
        "Not certified: the approximate ", level, " bound on the Monte Carlo ",
        "error of the estimate is ", bound, " statistical standard errors, ",
        "above the tolerance of ", mcml$tolerance, ". Refit with a larger ",
        "'max_samples' or 'max_stages' in set_control_mcml()."
      )
    }
  )
}

##' One sentence on the stages and draws of a certified MCML fit
##' @noRd
describe_mcml_draws <- function(mcml) {
  if (is.null(mcml$fitting_draws)) return(NULL)
  paste0(
    "Fitted in ", mcml$stages, " stage",
    if (identical(as.integer(mcml$stages), 1L)) "" else "s", ", with ",
    format(mcml$fitting_draws, big.mark = ","), " fitting and ",
    format(mcml$validation_draws, big.mark = ","),
    " validation draws at the last stage.",
    describe_noise_forms(mcml$noise_forms),
    if (isTRUE(mcml$integrated_locations > 0)) paste0(
      " The nugget was integrated out at ",
      format(mcml$integrated_locations, big.mark = ","), " of ",
      format(mcml$locations, big.mark = ","), " locations",
      if (length(mcml$noise_forms) > 1L) " in the last stage", "."
    )
  )
}

##' Which stages of an automatic fit sampled the nugget and which integrated it
##'
##' `forms` are the forms of the stages, NULL for a fit that did not choose
##' them (`integrate_noise = "auto"` in `set_control_mcml()`).
##' @noRd
describe_noise_forms <- function(forms) {
  if (is.null(forms)) return(NULL)
  stages <- function(form) {
    which <- which(forms == form)
    if (length(which) == 1L) return(paste("stage", which))
    paste0("stages ", paste(utils::head(which, -1L), collapse = ", "), " and ",
           utils::tail(which, 1L))
  }
  action <- c(sampled = "sampled the nugget",
              integrated = "integrated the nugget out")
  if (length(unique(forms)) == 1L) {
    return(paste0(" The automatic choice ", action[[forms[1L]]],
                  if (length(forms) > 1L) " in every stage", "."))
  }
  paste0(" The automatic choice sampled the nugget in ", stages("sampled"),
         " and integrated it out in ", stages("integrated"), ".")
}

##' One sentence on how many checks of a fit need a look
##'
##' NULL when the number is not known.
##' @noRd
describe_checks_to_review <- function(count) {
  if (is.null(count)) return(NULL)
  if (count == 0L) return("No check needs a look, see check_fit() for the details.")
  paste0(count, if (count == 1L) " check needs" else " checks need",
         " a look, see check_fit().")
}

##' Write sentences wrapped to the console width
##'
##' Each element of `text` starts a new line.
##' @noRd
write_wrapped <- function(text) {
  width <- max(40L, getOption("width", 80L) - 2L)
  for (sentence in text) {
    cat(strwrap(sentence, width = width), sep = "\n")
  }
  invisible(NULL)
}

##' @title Print a Fitted RiskMap Model
##' @description Print method for objects of class \code{"RiskMap"} from
##'   \code{\link{glgpm}}: the call, the estimates of the regression
##'   coefficients and of the covariance parameters on their natural scale,
##'   the log-likelihood and, for binomial and Poisson models, whether the
##'   Monte Carlo error of the estimate is certified. Standard errors and
##'   confidence intervals are given by \code{\link{summary.RiskMap}}, and the
##'   checks on the computation by \code{\link{check_fit}}.
##' @param x An object of class \code{"RiskMap"}.
##' @param digits The number of significant digits of the estimates.
##' @param ... other parameters.
##' @return Invisibly returns \code{x}.
##' @method print RiskMap
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
##' fit
##' @export
print.RiskMap <- function(x, digits = max(3L, getOption("digits") - 3L),
                          ...) {
  stopifnot(inherits(x, "RiskMap"))
  check_positive_integer(digits, "digits")
  if (!is.null(x$call)) {
    cat("Call:\n")
    cat(paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  }
  model <- c(gaussian = "Linear", binomial = "Binomial",
             poisson = "Poisson")[x$family]
  # Only a link other than the canonical one is named, since the
  # coefficients are then on another scale.
  link <- if (identical(x$link_function$name, "custom")) {
    paste0(" with ", x$link_function$link %||% "a custom", " link")
  }
  cat(if (is.na(model)) "Geostatistical" else model, " geostatistical model",
      link, ", ", length(x$y), " observations at ", NROW(x$coords),
      " locations\n", sep = "")

  estimates <- coef(x)
  cat("\nRegression coefficients\n")
  print_estimates(estimates$beta, digits)
  covariance <- c(
    "Spatial process var." = estimates$sigma2,
    stats::setNames(estimates$phi, paste0("Spatial corr. scale (",
                                          x$distance_units, ")")),
    "Variance of the nugget" = estimates$tau2,
    "Measurement error var." = estimates$sigma2_me,
    if (!is.null(estimates$sigma2_re)) stats::setNames(
      estimates$sigma2_re, paste0(names(x$re) %||% names(estimates$sigma2_re),
                                  " (random eff. var.)")
    )
  )
  cat("\nCovariance parameters (Matern, kappa = ", x$kappa, ")\n", sep = "")
  print_estimates(covariance, digits)
  if (is.numeric(x$fix_tau2)) {
    cat("Variance of the nugget effect fixed at ", x$fix_tau2, "\n", sep = "")
  }
  if (identical(x$family, "gaussian") && is.numeric(x$fix_var_me)) {
    cat("Measurement error var. fixed at ", x$fix_var_me, "\n", sep = "")
  }

  cat("\nLog-likelihood: ", format(x$log_lik, digits = digits + 3L), "\n",
      sep = "")
  mcml <- summarise_mcml_fit(x)
  if (!is.null(mcml)) {
    write_wrapped(c(describe_mcml_certificate(mcml),
                    describe_checks_to_review(mcml$checks_to_review)))
  }
  cat("\nsummary() gives standard errors and confidence intervals.\n")
  invisible(x)
}

##' Print named estimates as a column, each with its own notation
##' @noRd
print_estimates <- function(values, digits) {
  column <- matrix(vapply(values, format, character(1), digits = digits),
                   ncol = 1L, dimnames = list(names(values), "Estimate"))
  print.default(column, quote = FALSE, right = TRUE)
  invisible(NULL)
}
