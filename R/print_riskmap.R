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
##' `certification_bank` and the number of checks of `check_fit()` that need a
##' look (`checks_to_review`).
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
    certification_bank = if (!is.null(last)) summarise_certification_bank(last)
  )
  mcml$checks_to_review <- sum(list_fit_checks(object, mcml)$status != "OK")
  mcml
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
    " validation draws at the last stage."
  )
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
  cat(if (is.na(model)) "Geostatistical" else model, " geostatistical model, ",
      length(x$y), " observations at ", NROW(x$coords), " locations\n",
      sep = "")

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
