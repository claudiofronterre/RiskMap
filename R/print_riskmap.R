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
##' `certification_bank`, how many levels of each term were integrated out of
##' the likelihood at the last stage (`integration`, a data frame with the
##' `term`, the levels `integrated` and the `levels` in all, the nugget's
##' levels being the locations), the forms of the stages of an automatic fit
##' (`noise_forms`), the largest
##' change of a factor at the last fixed refit of the last stage
##' (`refit_change`), whether to suggest integrating the nugget out
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
    correction_left_bounds = isTRUE(last$correction_left_bounds),
    variance_rising = isTRUE(last$held_variance_rising),
    stages = optimizer$mcml_stages %||% length(history),
    fitting_draws = optimizer$mcml_budget,
    validation_draws = optimizer$mcml_validation_budget,
    certification_bank = if (!is.null(last)) summarise_certification_bank(last),
    integration = describe_fit_integration(optimizer),
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
##' `threshold` or more, or a random effect with such a variance and some
##' levels within one location. Where the variance at the maximum is small
##' the likelihood is nearly flat along it, and there the integrated estimate
##' of its score is the noisier of the two forms. In the larger test of the
##' integrated nugget it certified 5 runs against 12 of 21 for the sampled
##' form where that variance was below 0.1, 9 and 9 of 18 between 0.1 and 0.3,
##' and 68 against 28 of 75 above 0.3. A variance estimated as zero gets no
##' suggestion, nor does a fit that integrated terms out or chose the form
##' itself (`integrate_noise = "auto"`).
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
  if (isTRUE(variance >= threshold)) return(TRUE)
  ID_re <- object$ID_re
  re_variances <- exp(unlist(object$estimate$sigma2_re %||% numeric()))
  if (is.null(ID_re) || !length(re_variances) || is.null(object$ID_coords)) {
    return(FALSE)
  }
  ID_re <- as.matrix(ID_re)
  any(vapply(seq_len(ncol(ID_re)), function(j) {
    within <- tapply(object$ID_coords, ID_re[, j], function(locations) {
      length(unique(locations)) == 1L
    })
    isTRUE(re_variances[j] >= threshold) && any(within)
  }, logical(1)))
}

##' How many levels of each term a fit integrated out at its last stage
##'
##' The data frame of `make_adaptive_mcml_engine()`'s `integration_counts`
##' from the optimizer's record, or for a fit recorded before random effects
##' could be integrated out, the nugget's row from its integrated locations.
##' NULL for a fit that integrated nothing out.
##' @noRd
describe_fit_integration <- function(optimizer) {
  counts <- optimizer$mcml_integration_counts
  if (!is.null(counts)) {
    return(if (isTRUE(optimizer$mcml_integrate_noise)) counts)
  }
  located <- optimizer$mcml_integrated_locations %||% 0L
  if (located > 0) {
    data.frame(term = "nugget", integrated = located,
               levels = optimizer$mcml_locations %||% NA_integer_,
               stringsAsFactors = FALSE)
  }
}

##' One sentence on the Monte Carlo certificate of a fit
##'
##' `mcml` is the list of `summarise_mcml_fit()`. The sentence is formatted
##' with cli, so argument and function names look as in RiskMap's messages.
##' `certificate_symbol()` gives the bullet it is printed with.
##' @noRd
describe_mcml_certificate <- function(mcml) {
  bound <- signif(mcml$certificate, 2)
  level <- paste0(100 * mcml$level, "%")
  tolerance <- mcml$tolerance
  refit <- function(arguments) {
    paste0("Refit with a larger ", arguments, " in {.fun set_control_mcml}.")
  }
  text <- switch(
    mcml$outcome,
    "certified" = paste0(
      "Certified. The Monte Carlo error of the estimate is within {bound} ",
      "statistical standard errors (an approximate {level} bound, tolerance ",
      "{tolerance})."
    ),
    "weak identification" = paste0(
      "No certificate. The observed information is singular, so the Monte ",
      "Carlo precision of the estimate is not established."
    ),
    "numerically independent" = paste0(
      "No certificate. The fit is numerically one of independent location ",
      "effects, and the certificate needs an interior maximum."
    ),
    "no maximum" = paste0(
      "No certificate. The likelihood has no maximum, so there is no ",
      "maximum likelihood estimate whose Monte Carlo error could be bounded."
    ),
    if (isTRUE(mcml$variance_rising)) {
      paste0(
        "Not certified. The likelihood rises as a variance held at zero ",
        "leaves it, so the estimate is not a maximum. ",
        refit("{.arg max_samples} or {.arg max_stages}")
      )
    } else if (mcml$correction_refused) {
      paste0(
        "Not certified. The Newton correction of the estimate did not ",
        "converge on the validation draws, so the approximate {level} bound ",
        "of {bound} statistical standard errors may not hold. ",
        refit("{.arg max_samples}")
      )
    } else {
      paste0(
        "Not certified. The approximate {level} bound on the Monte Carlo ",
        "error of the estimate is {bound} statistical standard errors, above ",
        "the tolerance of {tolerance}. ",
        refit("{.arg max_samples} or {.arg max_stages}")
      )
    }
  )
  cli::format_inline(text)
}

##' The bullet of the certificate of a fit, a tick only when it is certified
##'
##' Every other outcome is a problem in the check of the Monte Carlo accuracy
##' of `check_fit()`, so it gets a cross.
##' @noRd
certificate_symbol <- function(mcml) {
  if (identical(mcml$outcome, "certified")) "v" else "x"
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
    describe_integration_counts(mcml$integration,
                                length(mcml$noise_forms) > 1L)
  )
}

##' Sentences on the levels of each term integrated out at the last stage
##'
##' `counts` is the data frame of `describe_fit_integration()`, and
##' `last_stage` adds that the counts are those of the last stage. A term
##' integrated out nowhere is said to have stayed sampled. NULL when nothing
##' was integrated out.
##' @noRd
describe_integration_counts <- function(counts, last_stage = FALSE) {
  if (is.null(counts) || !any(counts$integrated > 0)) return(NULL)
  number <- function(value) format(value, big.mark = ",", trim = TRUE)
  sentences <- vapply(seq_len(nrow(counts)), function(k) {
    nugget <- identical(counts$term[k], "nugget")
    name <- if (nugget) "The nugget" else paste("The random effect",
                                                 counts$term[k])
    unit <- if (nugget) "locations" else "levels"
    integrated <- counts$integrated[k]
    levels <- counts$levels[k]
    if (integrated == 0) return(paste0(" ", name, " stayed sampled."))
    paste0(" ", name, " was integrated out at ", number(integrated), " of ",
           number(levels), " ", unit,
           if (integrated < levels) {
             paste0(" and sampled at the other ", number(levels - integrated))
           },
           if (last_stage) " in the last stage", ".")
  }, character(1))
  paste(sentences, collapse = "")
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
##' With its bullet as its name, a tick when no check needs a look. NULL when
##' the number is not known.
##' @noRd
describe_checks_to_review <- function(count) {
  if (is.null(count)) return(NULL)
  check_fit <- out_fun("check_fit")
  if (count == 0L) {
    return(c(v = paste0("No check needs a look. See ", check_fit,
                        " for the details.")))
  }
  c("!" = paste0(count, if (count == 1L) " check needs" else " checks need",
                 " a look. See ", check_fit, "."))
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
  check_class(x, "RiskMap")
  check_positive_integer(digits, "digits")
  out_print(format_fit(x, digits))
  invisible(x)
}

##' The lines that print a fitted model
##' @noRd
format_fit <- function(x, digits) {
  traits <- fit_traits(x)
  # Only a link other than the canonical one is named, since the
  # coefficients are then on another scale.
  title <- paste0(model_title(traits), describe_link(fit_link(x)))
  estimates <- coef(x)
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
  fixed <- c(
    if (is.numeric(x$fix_tau2)) {
      paste0("Variance of the nugget effect fixed at ", x$fix_tau2)
    },
    if (isTRUE(traits$measurement_error) && is.numeric(x$fix_var_me)) {
      paste0("Measurement error var. fixed at ", x$fix_var_me)
    }
  )
  mcml <- summarise_mcml_fit(x)
  certificate <- if (!is.null(mcml)) {
    c(stats::setNames(describe_mcml_certificate(mcml),
                      certificate_symbol(mcml)),
      describe_checks_to_review(mcml$checks_to_review))
  }
  summary_hint <- paste0(out_fun("summary", "summary.RiskMap"),
                         " gives standard errors and confidence intervals.")
  c(out_rule(title),
    format_call(x$call),
    out_text(paste0(format_count(length(x$y)), " observations at ",
                    format_count(NROW(x$coords)), " locations.")),
    "",
    out_heading("Regression coefficients"),
    format_estimates(estimates$beta, digits),
    "",
    out_heading(paste0("Covariance parameters (Matern, kappa = ", x$kappa,
                       ")")),
    format_estimates(covariance, digits),
    out_text(fixed),
    "",
    out_pairs("Log-likelihood", format(x$log_lik, digits = digits + 3L)),
    "",
    out_bullets(c(certificate, i = summary_hint)))
}

##' The title of a fitted model
##'
##' "Linear geostatistical model" for a family whose likelihood is exact, and
##' otherwise the name of the family, as in "Binomial geostatistical model".
##' `traits` are those of `fit_traits()`, NULL for a family RiskMap doesn't
##' know.
##' @noRd
model_title <- function(traits) {
  if (is.null(traits)) return("Geostatistical model")
  if (isTRUE(traits$exact)) return("Linear geostatistical model")
  family <- traits$family
  paste0(toupper(substr(family, 1L, 1L)), substring(family, 2L),
         " geostatistical model")
}

##' The link of a fit when it is not the canonical one
##'
##' The name of a named link, such as "probit", or "custom" for an inverse
##' link given as a function. NULL for the canonical link, whose coefficients
##' are on the usual scale.
##' @noRd
fit_link <- function(fit) {
  if (!identical(fit$link_function$name, "custom")) return(NULL)
  fit$link_function$link %||% "custom"
}

##' The words that add a link other than the canonical one to a title, as in
##' " with probit link" or " with a custom link", and nothing for NULL
##' @noRd
describe_link <- function(link) {
  if (is.null(link)) return("")
  if (identical(link, "custom")) return(" with a custom link")
  paste0(" with ", link, " link")
}

##' The call of a fit, faint, as R deparses it
##' @noRd
format_call <- function(call) {
  if (is.null(call)) return(character())
  out_faint(sub("\\s+$", "", deparse(call)))
}

##' A count with a comma between thousands
##' @noRd
format_count <- function(count) {
  format(count, big.mark = ",", scientific = FALSE, trim = TRUE)
}

##' Named estimates as a column, each with its own notation
##' @noRd
format_estimates <- function(values, digits) {
  column <- matrix(vapply(values, format, character(1), digits = digits),
                   ncol = 1L, dimnames = list(names(values), "Estimate"))
  out_capture(print.default(column, quote = FALSE, right = TRUE))
}
