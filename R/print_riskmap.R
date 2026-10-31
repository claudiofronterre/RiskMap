##' Certified MCML part of a model summary
##'
##' Reads the certificate of a non-Gaussian fit from the attributes that
##' `glgpm()` writes, with short notes on what the warnings of the fit
##' reported, so that both are available after fitting, for example from a
##' saved model. NULL for Gaussian fits, which have no Monte Carlo error.
##' @return A list with the `certificate` (an approximate bound at `level`, in
##' statistical standard errors, on the Monte Carlo error of the estimate,
##' NA when there is none), the `tolerance`, whether the fit is `certified`,
##' its `outcome`, whether the repeated Newton correction of the last stage
##' was refused (`correction_refused`), the number of `stages`, the
##' `fitting_draws` and `validation_draws` of the last stage, and the
##' `notes`.
##' @noRd
summarise_mcml_fit <- function(object) {
  optimizer <- attr(object, "optimizer")
  if (is.null(optimizer$mcml_outcome)) return(NULL)
  history <- object$mcml_history %||% optimizer$mcml_history
  last <- if (length(history)) history[[length(history)]] else NULL
  list(
    certificate = optimizer$mcml_certificate %||% NA_real_,
    level = 0.95,
    tolerance = attr(object, "control_mcml")$tolerance %||% NA_real_,
    certified = isTRUE(optimizer$mcml_certified),
    outcome = optimizer$mcml_outcome,
    correction_refused = isTRUE(last$residual_refused),
    stages = optimizer$mcml_stages %||% length(history),
    fitting_draws = optimizer$mcml_budget,
    validation_draws = optimizer$mcml_validation_budget,
    notes = describe_mcml_notes(optimizer)
  )
}

##' One sentence on the Monte Carlo certificate of a fit
##'
##' `mcml` is the list of `summarise_mcml_fit()`. With `caveat`, the sentence
##' on a certified fit adds that the bound says nothing about other maxima.
##' @noRd
describe_mcml_certificate <- function(mcml, caveat = TRUE) {
  bound <- signif(mcml$certificate, 2)
  level <- paste0(100 * mcml$level, "%")
  switch(
    mcml$outcome,
    "certified" = paste0(
      "Certified: the Monte Carlo error of the estimate is within ", bound,
      " statistical standard errors (an approximate ", level, " bound, ",
      "tolerance ", mcml$tolerance, ").",
      if (caveat) paste0(" The bound concerns the maximum found, not whether ",
                         "it is the highest.")
    ),
    "weak identification" = paste0(
      "No certificate: the observed information is singular, so the Monte ",
      "Carlo precision of the estimate is not established."
    ),
    "numerically independent" = paste0(
      "No certificate: the fit is numerically one of independent location ",
      "effects, and the certificate needs an interior maximum."
    ),
    if (mcml$correction_refused) {
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

##' Short notes on the warnings of a certified MCML fit
##'
##' One sentence for each condition that a warning reported during the fit,
##' read from the flags of the `optimizer` attribute, with the advice that
##' the warning gave.
##' @noRd
describe_mcml_notes <- function(optimizer) {
  boundary <- optimizer$mcml_boundary
  single <- length(boundary) == 1L
  notes <- c(
    if (length(boundary)) paste0(
      "The ", paste(boundary, collapse = " and the "), " ",
      if (single) "was" else "were", " estimated as effectively zero, so ",
      if (single) "its standard error is" else "their standard errors are",
      " not reported. Consider refitting without ",
      if (single) "it." else "them."
    ),
    if (isTRUE(optimizer$mcml_numerically_independent)) paste0(
      "The spatial correlation at the fitted scale is numerically negligible, ",
      "so the scale is unresolved. See profile(fit, parameters = \"phi\"), or ",
      "consider a model without spatial correlation."
    ),
    if (isTRUE(optimizer$mcml_identifiability_failure)) paste0(
      "Some parameters, often the spatial scale, are not identifiable from ",
      "these data. Consider a simpler model or fixing the scale."
    ),
    if (isTRUE(optimizer$mcml_small_scale)) paste0(
      "The spatial scale is shorter than the median distance between ",
      "neighbouring locations, so it may be weakly identified and its ",
      "interval very wide. See profile(fit, parameters = \"phi\")."
    ),
    if (isTRUE(optimizer$mcml_near_limit)) paste0(
      "The fitted likelihood is not clearly above the limit with independent ",
      "location effects at this Monte Carlo precision."
    ),
    if (isTRUE(optimizer$mcml_maxima_alarm)) paste0(
      "An exact bound shows that the estimate is not the maximum likelihood ",
      "estimate. Inspect the profile likelihood of the spatial scale."
    ),
    if (isTRUE(optimizer$mcml_ambiguous_maxima)) paste0(
      "Separate maxima of the likelihood fit about equally well, and the ",
      "estimate is from the higher. Wald intervals may describe the ",
      "uncertainty poorly, so consider the profile likelihood (profile())."
    ) else if (isTRUE(optimizer$mcml_better_maximum)) paste0(
      "The check for other maxima found a maximum higher than the one ",
      "reached from the starting values, and the estimate is taken from it."
    )
  )
  weak_spatial <- any(grepl("sigma2", boundary, fixed = TRUE)) ||
    isTRUE(optimizer$mcml_small_scale) ||
    isTRUE(optimizer$mcml_near_limit) ||
    isTRUE(optimizer$mcml_numerically_independent)
  if (weak_spatial && !isTRUE(optimizer$mcml_pc_prior)) {
    notes <- c(notes, paste0(
      "Refitting with 'control_mcml = set_control_mcml(pc_prior = TRUE)' ",
      "stabilises the spatial variance and scale with a mild penalty."
    ))
  }
  notes
}

##' Write sentences wrapped to the console width
##'
##' Each element of `text` starts a new line, after `bullet` if given, with
##' continuation lines indented to match.
##' @noRd
write_wrapped <- function(text, bullet = "") {
  width <- max(40L, getOption("width", 80L) - 2L)
  for (sentence in text) {
    cat(strwrap(sentence, width = width, initial = bullet,
                exdent = nchar(bullet)), sep = "\n")
  }
  invisible(NULL)
}

##' @title Print a Fitted RiskMap Model
##' @description Print method for objects of class \code{"RiskMap"} from
##'   \code{\link{glgpm}}: the call, the estimates of the regression
##'   coefficients and of the covariance parameters on their natural scale,
##'   the log-likelihood and, for binomial and Poisson models, whether the
##'   Monte Carlo error of the estimate is certified. Standard errors,
##'   confidence intervals and notes on the warnings of the fit are given by
##'   \code{\link{summary.RiskMap}}.
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
    write_wrapped(describe_mcml_certificate(mcml, caveat = FALSE))
  }
  notes <- length(mcml$notes)
  cat("\nsummary() gives standard errors",
      if (notes) paste0(", confidence intervals and ", notes, " note",
                        if (notes > 1L) "s", " on this fit") else
        " and confidence intervals", ".\n", sep = "")
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
