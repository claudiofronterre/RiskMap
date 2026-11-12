##' Logarithms of a custom link's probability that stay finite at its limits
##'
##' A custom inverse link is known only through its value. A probability within
##' 2^-53 of one rounds to exactly one, for cloglog from a linear predictor of 3.6
##' and for probit from 8.3, although the model reaches one only in the limit. The
##' failure probability is then known only to be at most 2^-53, and its logarithm
##' is taken at that bound. A probability near zero keeps its relative precision,
##' so its logarithm is bounded only where it underflows. Each bound matters only
##' for an observation with a failure, or a success, in that position, and it
##' overstates the likelihood of a draw that is already at least 2^-53 times less
##' likely per failure than the same draw below the rounding point. The share of
##' such a draw in an importance-sampling mean is negligible near a maximum, and
##' being finite it cannot make a sum or a central difference undefined.
##' @noRd
bound_log_probability <- function(probability) {
  list(success = log(pmax(probability, .Machine$double.xmin)),
       failure = log1p(-pmin(probability, 1 - .Machine$double.eps / 2)))
}

##' Multiply a term by a count that may be zero
##'
##' A zero count contributes nothing to the log-likelihood or its derivatives,
##' whatever the term, which can overflow at a bounded probability.
##' @noRd
multiply_by_count <- function(count, term) {
  product <- count * term
  product[count == 0] <- 0
  product
}

##' Log-likelihood, score and curvature of an observation with a custom link
##'
##' `curvature` is the negative second derivative of the log-likelihood, written
##' through d log p = p' / p and d log(1 - p) = -p' / (1 - p), so that no product
##' of tail probabilities underflows, as (p (1 - p))^2 does below a probit linear
##' predictor of -27. Probabilities are bounded as in `bound_log_probability()`.
##' The curvature is negative where the log-likelihood is convex, as it is for a
##' link that is not log-concave, such as that of an imperfect test. Returns NULL
##' if the mean is outside the support of the family, or if a derivative is not
##' finite. `with_curvature = FALSE` omits the curvature for callers that have no
##' second derivative. `log_likelihoods` holds each observation's term of the
##' log-likelihood.
##' @noRd
compute_custom_link_terms <- function(eta, y, units_m, family, linkf,
                                      with_curvature = TRUE) {
  mean <- linkf$inv(eta)
  if (any(!is.finite(mean)) || any(mean < 0)) return(NULL)
  first <- linkf$d1(eta)
  second <- if (with_curvature) linkf$d2(eta)
  if (family == "poisson") {
    positive <- pmax(mean, .Machine$double.xmin)
    rate <- first / positive
    log_likelihoods <- y * log(pmax(units_m * mean, .Machine$double.xmin)) -
      units_m * mean
    score <- multiply_by_count(y, rate) - units_m * first
    curvature <- if (with_curvature) {
      multiply_by_count(y, rate^2 - second / positive) + units_m * second
    }
  } else {
    if (any(mean > 1)) return(NULL)
    success <- pmax(mean, .Machine$double.xmin)
    failure <- 1 - pmin(mean, 1 - .Machine$double.eps / 2)
    log_probability <- bound_log_probability(mean)
    failures <- units_m - y
    success_rate <- first / success
    failure_rate <- first / failure
    log_likelihoods <- y * log_probability$success +
      failures * log_probability$failure
    score <- multiply_by_count(y, success_rate) -
      multiply_by_count(failures, failure_rate)
    curvature <- if (with_curvature) {
      multiply_by_count(y, success_rate^2 - second / success) +
        multiply_by_count(failures, failure_rate^2 + second / failure)
    }
  }
  if (any(!is.finite(score)) ||
      (with_curvature && any(!is.finite(curvature)))) {
    return(NULL)
  }
  list(log_likelihood = sum(log_likelihoods), score = score,
       curvature = curvature, log_likelihoods = log_likelihoods)
}

##' Move the start of a search off the limits of a custom link
##'
##' Past the point where a probability rounds to one, the likelihood of an
##' observation with a failure is flat (`bound_log_probability()`), so a search
##' started there has no slope to follow and stops at once. Each such observation,
##' and each with a success whose probability has rounded to zero, or a positive
##' count whose mean has, moves by whole units of its first latent effect until its
##' mean has left the limit. `columns` gives that effect for each observation and
##' `linear_predictor` maps the latent effects to the linear predictors.
##' @noRd
move_start_off_limits <- function(latent, linear_predictor, y, units_m, family,
                                  linkf, columns, maximum_moves = 100L) {
  for (move in seq_len(maximum_moves)) {
    mean <- linkf$inv(linear_predictor(latent))
    too_low <- which(y > 0 & mean <= 0)
    too_high <- if (family == "binomial") {
      which(y < units_m & mean >= 1 - .Machine$double.eps / 2)
    } else {
      integer()
    }
    if (!length(too_low) && !length(too_high)) break
    latent[columns[too_low]] <- latent[columns[too_low]] + 1
    latent[columns[too_high]] <- latent[columns[too_high]] - 1
  }
  latent
}

##' Exact inverse and derivatives of a named binomial link
##'
##' The inverses are written without the clamping of [stats::binomial()],
##' because `bound_log_probability()` already handles probabilities of 0 and
##' 1, and the derivatives are exact so that no numerical differences are
##' needed. The second derivative of the complementary log-log link is
##' written as a difference of two exponentials, which stays finite where
##' `exp(eta)` overflows.
##' @noRd
named_inverse_link <- function(link) {
  switch(
    link,
    probit = list(
      inv = function(eta) pnorm(eta),
      d1 = function(eta) dnorm(eta),
      d2 = function(eta) -eta * dnorm(eta)
    ),
    cloglog = list(
      inv = function(eta) -expm1(-exp(eta)),
      d1 = function(eta) exp(eta - exp(eta)),
      d2 = function(eta) exp(eta - exp(eta)) - exp(2 * eta - exp(eta))
    ),
    cauchit = list(
      inv = function(eta) pcauchy(eta),
      d1 = function(eta) dcauchy(eta),
      d2 = function(eta) -2 * eta / (pi * (1 + eta^2)^2)
    ),
    riskmap_bug("No inverse is defined for the {link} link.")
  )
}

##' Validate and construct inverse-link functions
##'
##' This is the single internal representation used by the conditional-mode,
##' sampling and non-Gaussian likelihood code. Keeping the validation and
##' derivative construction here prevents those paths from interpreting a
##' user-supplied link differently.
##'
##' @noRd
make_inverse_link_functions <- function(family, invlink, ncheck,
                                        second_derivative = TRUE) {
  # `what` names the function checked, with cli markup.
  validate <- function(f, what) {
    what <- cli::format_inline(what)
    if (!is.function(f)) {
      abort_input("{what} must be a function, not {.obj_type_friendly {f}}.", arg = "invlink")
    }
    x <- rep(0, ncheck)
    out <- tryCatch(f(x), error = function(e) e)
    if (inherits(out, "error")) {
      abort_input("{what} failed on a numeric vector.", arg = "invlink", parent = out)
    }
    if (!is.numeric(out)) {
      abort_input("{what} must return a numeric vector, not {.obj_type_friendly {out}}.", arg = "invlink")
    }
    if (length(out) != length(x)) {
      abort_input("{what} must return a vector as long as its input.", arg = "invlink")
    }
    if (!all(is.finite(out))) {
      abort_input("{what} returns values that aren't finite.", arg = "invlink")
    }
    invisible(TRUE)
  }

  if (is.null(invlink)) {
    if (family == "poisson") {
      inv <- function(x) exp(x)
      d1 <- function(x) exp(x)
      d2 <- if (second_derivative) function(x) exp(x) else NULL
    } else {
      inv <- function(x) plogis(x)
      d1 <- function(x) {
        probability <- inv(x)
        probability * (1 - probability)
      }
      d2 <- if (second_derivative) function(x) {
        probability <- inv(x)
        derivative <- probability * (1 - probability)
        derivative * (1 - 2 * probability)
      } else NULL
    }
    validate(inv, "The canonical inverse link")
    validate(d1, "The derivative of the canonical inverse link")
    if (second_derivative) {
      validate(d2, "The second derivative of the canonical inverse link")
    }
    return(list(inv = inv, d1 = d1, d2 = d2, name = "canonical",
                link = if (family == "poisson") "log" else "logit"))
  }

  # The list form is internal. It carries a named link from read_family(),
  # and the stored link of a fit into its refits (refit_invlink()).
  link <- NULL
  if (is.function(invlink)) {
    inv <- invlink
    d1 <- NULL
    d2 <- NULL
  } else if (is.list(invlink)) {
    unknown <- setdiff(names(invlink), c("inv", "d1", "d2", "link"))
    if (is.null(names(invlink)) || length(unknown)) {
      abort_input("A list in {.arg invlink} takes only the elements {.code inv}, {.code d1}, {.code d2} and {.code link}.", arg = "invlink")
    }
    inv <- invlink$inv
    d1 <- invlink$d1
    d2 <- invlink$d2
    link <- invlink$link
  } else {
    abort_input("{.arg invlink} must be {.code NULL} or a function, not {.obj_type_friendly {invlink}}.", arg = "invlink")
  }

  validate(inv, "{.arg invlink}")
  check_invlink(inv, family)
  check_invlink_tails(inv, family)
  if (is.null(d1)) d1 <- differentiate(inv, ncheck)
  validate(d1, "The derivative of {.arg invlink}")
  if (second_derivative) {
    if (is.null(d2)) d2 <- differentiate(d1, ncheck)
    validate(d2, "The second derivative of {.arg invlink}")
  } else {
    d2 <- NULL
  }

  list(inv = inv, d1 = d1, d2 = d2, name = "custom", link = link)
}


##' The family of a model and the likelihood of its records
##'
##' One object through which the MCML engine sees the family and the link of a
##' model (#253). `family` is "binomial" or "poisson", and `linkf` the inverse
##' link functions of `make_inverse_link_functions()`, NULL for the canonical
##' link. Where records have effects integrated out, `kernel` is their kernel
##' at the variances of the classes (`make_integrated_kernel()`), which gives
##' the likelihood of the records in the place of the link, and the object
##' keeps every other property and method of the model's link. Each method is
##' the arithmetic the engine used before, moved without change, so that fits
##' stay identical to the last bit, and each property is a test of the
##' canonical link that the engine made, under the name of what it decides.
##'
##' Properties.
##' - `canonical`, `poisson`, of the model's link, and `integrated`, whether
##'   `kernel` is given.
##' - `compiled`, the weights and rows of the log-likelihood are compiled
##'   (`cpp_importance_log_weights()`, `cpp_canonical_log_likelihood()`), with
##'   `poisson` as their flag. The compiled weights evaluate an integrated
##'   kernel through its table themselves.
##' - `assumes_concave`, the mode search factorises its Hessian without a
##'   guard and returns it unfloored, for the canonical link without a kernel.
##' - `moves_start`, the mode search starts off the limits of a custom link
##'   without a kernel (`move_start_off_limits()`).
##' - `needs_finite_weight`, a bank with no finite weight is a failed
##'   evaluation.
##' - `legacy_link`, what the code that step 1 has not moved yet takes in the
##'   place of the link, the kernel where there is one and the link otherwise.
##'
##' Methods, for records with `y` events of `units_m` trials, or counts with
##' exposure `units_m`, at linear predictors `eta`.
##' - `terms(eta, y, units_m)`, the log-likelihood, score and curvature, or NULL
##'   where they are not defined (`laplace_observation_terms()`).
##' - `log_likelihood(eta, y, units_m)`, the log-likelihood of each row of a
##'   matrix of linear predictors, one column per record
##'   (`compute_observation_log_likelihood()`).
##' - `slope_parts(eta, y, units_m)`, the score and curvature of each record
##'   at each row of `eta` under the model's link, for the slope at a held
##'   variance, or NULL where a custom link cannot give them
##'   (`compute_slope_terms()`). Records integrated out take theirs from their
##'   table instead.
##' - `move_start(latent, linear_predictor, y, units_m, columns)`, as
##'   `move_start_off_limits()`.
##' @noRd
make_family <- function(family, linkf = NULL, kernel = NULL) {
  if (!family %in% c("binomial", "poisson")) {
    riskmap_bug("The family object is defined for the binomial and Poisson families only.")
  }
  if (is.null(linkf)) linkf <- make_inverse_link_functions(family, NULL, 1L)
  if (is_integrated_link(linkf)) {
    riskmap_bug("The kernel of integrated records goes in {.arg kernel}, not in {.arg linkf}.")
  }
  canonical <- identical(linkf$name, "canonical")
  integrated <- !is.null(kernel)
  poisson <- family == "poisson"
  terms <- if (integrated) {
    function(eta, y, units_m) kernel$terms(eta)
  } else if (!canonical) {
    function(eta, y, units_m) {
      compute_custom_link_terms(eta, y, units_m, family, linkf)
    }
  } else if (poisson) {
    function(eta, y, units_m) canonical_poisson_terms(eta, y, units_m, linkf)
  } else {
    function(eta, y, units_m) canonical_binomial_terms(eta, y, units_m, linkf)
  }
  log_likelihood <- if (integrated) {
    function(eta, y, units_m) kernel$log_likelihood(eta)
  } else if (canonical) {
    function(eta, y, units_m) {
      cpp_canonical_log_likelihood(eta, y, units_m, poisson = poisson)
    }
  } else {
    function(eta, y, units_m) {
      compute_link_log_likelihood_rows(eta, y, units_m, family, linkf)
    }
  }
  slope_parts <- if (canonical) {
    function(eta, y, units_m) {
      canonical_slope_parts(eta, y, units_m, family)
    }
  } else {
    function(eta, y, units_m) custom_slope_parts(eta, y, units_m, family, linkf)
  }
  structure(
    list(
      family = family, linkf = linkf, kernel = kernel, canonical = canonical,
      integrated = integrated, poisson = poisson, compiled = canonical,
      assumes_concave = canonical && !integrated,
      moves_start = !canonical && !integrated,
      needs_finite_weight = !canonical,
      legacy_link = if (integrated) kernel else linkf,
      terms = terms, log_likelihood = log_likelihood,
      slope_parts = slope_parts,
      move_start = function(latent, linear_predictor, y, units_m, columns) {
        move_start_off_limits(latent, linear_predictor, y, units_m, family,
                              linkf, columns)
      }
    ),
    class = "riskmap_family"
  )
}

##' The family object for a family string and a link or an integrated kernel
##'
##' For the old entry points, which take in `linkf` either a link or the kernel
##' of integrated records (`make_integrated_kernel()`). With a kernel the
##' model's link is not known, and the object is built on the canonical link,
##' whose properties those entry points do not read. Removed with them in step
##' 1d of #253.
##' @noRd
family_of_link <- function(family, linkf) {
  if (is_integrated_link(linkf)) return(make_family(family, kernel = linkf))
  make_family(family, linkf)
}

##' The family object of a fitted model
##'
##' Rebuilt from what the fit stores, the family string and the
##' `link_function` list, for the binomial and Poisson families. It gives the
##' likelihood of the fit's records under its link. Where the fit integrated
##' effects out, the likelihood of those records also needs the fit's plan of
##' integrated records and the variances of its classes (`kernel`), which
##' this does not rebuild.
##' @noRd
fit_family <- function(fit) {
  make_family(fit$family, fit$link_function)
}

##' Observation terms of the canonical log link of the Poisson family
##'
##' The arm of `laplace_observation_terms()` for the canonical Poisson link,
##' with the formulas of `maxim_integrand()`, and the exact canonical terms
##' where those formulas fail at the limits of the link
##' (`replace_at_canonical_limits()`). NULL if a linear predictor or a mean is
##' not finite.
##' @noRd
canonical_poisson_terms <- function(eta, y, units_m, linkf) {
  if (any(!is.finite(eta))) return(NULL)
  mean <- linkf$inv(eta)
  if (any(!is.finite(mean))) return(NULL)
  first <- linkf$d1(eta)
  second <- linkf$d2(eta)
  residual <- y / mean - units_m
  terms <- list(
    log_likelihoods = y * log(pmax(units_m * mean, .Machine$double.eps)) -
      units_m * mean,
    score = residual * first,
    curvature = y * first^2 / mean^2 - residual * second
  )
  terms <- replace_at_canonical_limits(terms, mean <= 0, eta, y, units_m,
                                       "poisson")
  list(log_likelihood = sum(terms$log_likelihoods), score = terms$score,
       curvature = terms$curvature)
}

##' Observation terms of the canonical logit link of the binomial family
##'
##' The arm of `laplace_observation_terms()` for the canonical binomial link,
##' with the formulas of `maxim_integrand()`, and the exact canonical terms
##' where those formulas fail at the limits of the link
##' (`replace_at_canonical_limits()`). NULL if a linear predictor is not
##' finite.
##' @noRd
canonical_binomial_terms <- function(eta, y, units_m, linkf) {
  if (any(!is.finite(eta))) return(NULL)
  probability <- linkf$inv(eta)
  first <- linkf$d1(eta)
  second <- linkf$d2(eta)
  variance <- probability * (1 - probability)
  residual <- y - units_m * probability
  terms <- list(
    log_likelihoods = y * log(pmax(probability, .Machine$double.eps)) +
      (units_m - y) * log(pmax(1 - probability, .Machine$double.eps)),
    score = residual * first / variance,
    curvature = units_m * first^2 / variance -
      residual * (second / variance - first^2 * (1 - 2 * probability) /
                    variance^2)
  )
  terms <- replace_at_canonical_limits(
    terms, probability <= 0 | probability >= 1, eta, y, units_m, "binomial"
  )
  list(log_likelihood = sum(terms$log_likelihoods), score = terms$score,
       curvature = terms$curvature)
}

##' Exact observation terms of the canonical links at the limits of the link
##'
##' The inverse logit rounds to 0 below a linear predictor of about -709.8
##' and to 1 above about 36.74, and the mean of a count underflows to 0 below
##' about -745. There the formulas of `laplace_observation_terms()` and
##' `maxim_integrand()` divide zero by zero, although the likelihood is smooth
##' and finite. Complete separation by a covariate takes the linear predictor
##' of locations where every individual is positive past 36.74 well inside the
##' bounds of the search, so the latent mode must be found there too. The
##' canonical forms need no division. With p = plogis(eta) they are
##' log p = -softplus(-eta), log(1 - p) = -softplus(eta), the score
##' y (1 - p) - (m - y) p and the curvature m p (1 - p). For the log link
##' they are y log(m) + y eta - m exp(eta), y - m exp(eta) and m exp(eta),
##' where y log(m) is the constant of the usual formula. They replace the
##' per-observation `log_likelihoods`, `score` and `curvature` of `terms` only
##' where `at` is TRUE, so that terms away from the limits are unchanged bit
##' for bit. Short of the limits the usual formulas keep two weaknesses. Their
##' curvature is not finite between about -709.8 and -372 for the logit, or
##' below -372 for a positive count, and `adaptive_laplace_mode()` rejects a
##' step to such a point. They also bound log p, log(1 - p) and log(m mu)
##' below at log(.Machine$double.eps). The log-likelihood of an observation
##' with any positives is then flat below a linear predictor of about -36,
##' that of a binomial observation with any negatives is flat above about
##' 36.04, and both fall where the exact forms take over. Neither is reached
##' under complete separation, where the locations with positives have only
##' positives and lie in the upper tail, and the others lie in the lower one.
##' Exact forms throughout would remove both but would change every canonical
##' fit by rounding.
##' @noRd
replace_at_canonical_limits <- function(terms, at, eta, y, units_m, family) {
  if (!any(at)) return(terms)
  eta <- eta[at]
  y <- y[at]
  units_m <- units_m[at]
  if (family == "poisson") {
    mean <- exp(eta)
    log_likelihoods <- -units_m * mean
    # An observation with no count adds no y log(m), even when m is zero.
    positive <- y > 0
    log_likelihoods[positive] <- log_likelihoods[positive] +
      y[positive] * (log(units_m[positive]) + eta[positive])
    score <- y - units_m * mean
    curvature <- units_m * mean
  } else {
    success <- stats::plogis(eta)
    failure <- stats::plogis(-eta)
    log_likelihoods <- -y * softplus(-eta) - (units_m - y) * softplus(eta)
    score <- y * failure - (units_m - y) * success
    curvature <- units_m * success * failure
  }
  terms$log_likelihoods[at] <- log_likelihoods
  terms$score[at] <- score
  terms$curvature[at] <- curvature
  terms
}

##' Log-likelihood of every row of a matrix of linear predictors, custom link
##'
##' The arm of `compute_observation_log_likelihood()` for a custom link. A draw
##' whose mean is not finite or is outside the support has log-likelihood -Inf,
##' and a mean that has rounded to a limit of its range is inside the support
##' (`bound_log_probability()`).
##' @noRd
compute_link_log_likelihood_rows <- function(eta, y, units_m, family, linkf) {
  mean <- matrix(linkf$inv(as.vector(eta)), nrow = nrow(eta))
  outside <- !is.finite(mean) | mean < 0
  log_likelihood <- if (family == "poisson") {
    rowSums(sweep(log(pmax(mean, .Machine$double.xmin)), 2L, y, "*") -
              sweep(mean, 2L, units_m, "*"))
  } else {
    outside <- outside | mean > 1
    log_probability <- bound_log_probability(mean)
    rowSums(sweep(log_probability$success, 2L, y, "*") +
              sweep(log_probability$failure, 2L, units_m - y, "*"))
  }
  log_likelihood[rowSums(outside) > 0] <- -Inf
  log_likelihood
}

##' Score and curvature of each record at each row, canonical link
##'
##' The closed forms of `compute_variance_slope_terms()` for the canonical
##' links, one row per draw and one column per record.
##' @noRd
canonical_slope_parts <- function(eta, y, units_m, family) {
  rate <- if (identical(family, "poisson")) exp(eta) else stats::plogis(eta)
  expected <- sweep(rate, 2L, units_m, "*")
  curvature <- if (identical(family, "poisson")) expected else
    expected * (1 - rate)
  score <- rep(y, each = nrow(eta)) - expected
  list(score = score, curvature = curvature)
}

##' Score and curvature of each record at each row, custom link
##'
##' The arm of `compute_variance_slope_terms()` for a custom link, through
##' `compute_custom_link_terms()`, or NULL where they cannot be computed.
##' @noRd
custom_slope_parts <- function(eta, y, units_m, family, linkf) {
  terms <- compute_custom_link_terms(
    as.vector(eta), rep(y, each = nrow(eta)), rep(units_m, each = nrow(eta)),
    family, linkf
  )
  if (is.null(terms)) return(NULL)
  list(score = matrix(terms$score, nrow(eta)),
       curvature = matrix(terms$curvature, nrow(eta)))
}
