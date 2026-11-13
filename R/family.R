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
    # The log link of the Poisson and negative binomial families, and
    # otherwise the logit. A negative binomial fit keeps the log link of its
    # mean, although the engine fits its records with the logit
    # (engine_records()).
    log_link <- identical(family_traits(family)$default_link, "log")
    if (log_link) {
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
                link = if (log_link) "log" else "logit"))
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

##' The families of RiskMap and their links
##'
##' One entry for each family a user can give (#253, step 1c), and the
##' negative binomial family, internal until step 6b. `links` are the
##' links `read_family_link()` accepts, `default_link` the canonical one, which
##' a family string implies and which is fitted on the canonical path, and
##' `refused` the links of a stats family that it refuses, with the reason.
##' The name is not `canonical_link`, because `$canonical` would then match it
##' partially on an object without the property `canonical`. The other fields
##' do not depend on the link, and `family_traits()` gives them with the
##' family's methods.
##' - `exact`, the likelihood integrates the latent field out in closed form,
##'   so the fit is exact and has no Monte Carlo (`glgpm_lm()`).
##' - `measurement_error`, the model has a measurement error variance
##'   `sigma2_me`, which an observation-level random effect or the nugget can
##'   be confounded with.
##' - `warn_missing_denominator`, a missing `denominator`, taken to be 1, is
##'   worth a warning.
##' - `upper_limit`, the responses are bounded above by their denominators, so
##'   a record can sit at the upper limit of its mean, with every trial
##'   positive.
##' - `mean_range`, the range of the mean that an inverse link must respect,
##'   with `range_message`, the error when it does not.
##' - `quantity`, `no_events`, `none_everywhere` and `every_everywhere`, the
##'   words of `report_no_maximum()` for the fitted mean, a record without
##'   events, and data with no events or only events.
##' - `link_hint`, the advice when a stats family names a link that is not
##'   accepted, NULL for the advice to give an inverse link in `invlink`.
##' - `inverse_link`, the inverse of the canonical link as a function, which a
##'   simulation model and the scores of a fit without a stored link use.
##' - `discrete`, the responses are counts, so a prediction is scored by the
##'   nonrandomised PIT (AnPIT) and the CRPS of a discrete distribution.
##' - `whole_denominator`, the denominators are numbers of trials and so whole.
##' - `response_quantity`, the name of the predicted mean of the response in
##'   plots, NULL for the name of the response itself, and `probability`, the
##'   mean is a probability, shown as a percentage.
##' - `title`, the family's name at the start of a model's title, and
##'   `inverse_link_formula`, its canonical inverse link as `summary()` shows
##'   it.
##' - `parameters`, the names of the family's own parameters, held at the
##'   values given with the family (`negative_binomial()`), and
##'   `kernel_family`, for a family whose records the MCML engine fits as
##'   those of another family, that family (`engine_records()`), NULL for a
##'   family the engine fits as itself.
##' @noRd
riskmap_families <- list(
  gaussian = list(
    links = "identity", default_link = "identity", refused = character(),
    exact = TRUE, measurement_error = TRUE,
    warn_missing_denominator = FALSE,
    upper_limit = FALSE, mean_range = c(-Inf, Inf), range_message = NULL,
    quantity = "mean", no_events = NULL, none_everywhere = NULL,
    every_everywhere = NULL,
    link_hint = "Gaussian responses are fitted with the identity link only.",
    inverse_link = identity, discrete = FALSE, whole_denominator = FALSE,
    response_quantity = NULL, probability = FALSE,
    title = "Gaussian", inverse_link_formula = "x",
    parameters = character(), kernel_family = NULL
  ),
  binomial = list(
    links = c("logit", "probit", "cloglog", "cauchit"),
    default_link = "logit",
    refused = c(log = "the log link can give probabilities above 1"),
    exact = FALSE, measurement_error = FALSE,
    warn_missing_denominator = TRUE,
    upper_limit = TRUE, mean_range = c(0, 1),
    range_message = "For a binomial model, {.arg invlink} must return probabilities between 0 and 1.",
    quantity = "prevalence", no_events = "with no positives",
    none_everywhere = "None of the individuals tested is positive",
    every_everywhere = "Every individual tested is positive",
    link_hint = NULL,
    inverse_link = stats::plogis, discrete = TRUE, whole_denominator = TRUE,
    response_quantity = "prevalence", probability = TRUE,
    title = "Binomial", inverse_link_formula = "1 / (1 + exp(-x))",
    parameters = character(), kernel_family = NULL
  ),
  poisson = list(
    links = "log", default_link = "log",
    refused = c(
      identity = "the identity link can give negative means",
      sqrt = paste("its inverse, the square of the linear",
                   "predictor, gives the same mean to a positive",
                   "and a negative linear predictor")
    ),
    exact = FALSE, measurement_error = FALSE,
    warn_missing_denominator = FALSE,
    upper_limit = FALSE, mean_range = c(0, Inf),
    range_message = "For a Poisson model, {.arg invlink} must return non-negative means.",
    quantity = "mean", no_events = "with a zero count",
    none_everywhere = "Every count is zero", every_everywhere = NULL,
    link_hint = NULL,
    inverse_link = exp, discrete = TRUE, whole_denominator = FALSE,
    response_quantity = "rate", probability = FALSE,
    title = "Poisson", inverse_link_formula = "exp(x)",
    parameters = character(), kernel_family = NULL
  ),
  # Internal until its size can be estimated (#253, step 6b). Its records go
  # to the engine as binomial ones (engine_records()).
  negative_binomial = list(
    links = "log", default_link = "log", refused = character(),
    exact = FALSE, measurement_error = FALSE,
    warn_missing_denominator = FALSE,
    upper_limit = FALSE, mean_range = c(0, Inf),
    range_message = "For a negative binomial model, {.arg invlink} must return non-negative means.",
    quantity = "mean", no_events = "with a zero count",
    none_everywhere = "Every count is zero", every_everywhere = NULL,
    link_hint = "Negative binomial responses are fitted with the log link only.",
    inverse_link = exp, discrete = TRUE, whole_denominator = FALSE,
    response_quantity = "rate", probability = FALSE,
    title = "Negative binomial", inverse_link_formula = "exp(x)",
    parameters = "theta", kernel_family = "binomial"
  )
)

##' The properties and methods of a family that do not depend on its link
##'
##' The entry of `riskmap_families` for `family`, with `family` itself and
##' the family's methods for its data, each the arithmetic its caller used
##' before (#253, step 1c). For records with `y` events of `units_m` trials,
##' or counts with exposure `units_m`, or Gaussian responses `y`,
##' - `validate(y, units_m)`, stops if the responses are not of the family
##'   (`check_binomial()`), and otherwise returns TRUE invisibly.
##' - `start_coefficients(y, units_m, D, cov_offset, link, formula, data)`,
##'   the starting regression coefficients of `glgpm()`, from least squares,
##'   `binomial_start_regression()` with the link named `link`, or a Poisson
##'   regression of `formula` on `data` with the log link.
##' - `empirical_link(y, units_m)` and `transform_variance(y, units_m)`, the
##'   empirical link of the transformed-Gaussian start and the median of its
##'   delta-method variance (`select_nongaussian_start()`), NULL for an exact
##'   family.
##' - `start_log_likelihood(eta, y, units_m)`, the log-likelihood under the
##'   canonical link that compares starting values
##'   (`starting_laplace_log_marginal()`), up to a constant, or NULL where a
##'   Poisson mean is not finite, which rules the candidate out. NULL for an
##'   exact family.
##' - `fitted_mean(mean, units_m)`, the fitted value of a record from the mean
##'   of its inverse link, the expected count for the Poisson and negative
##'   binomial families.
##' - `simulate(means, units_m, parameters)`, one response for each mean, from
##'   `rbinom()` with `units_m` trials, `rpois()` with exposure `units_m`,
##'   `rnbinom()` with exposure `units_m` and size `parameters$theta`, or
##'   `rnorm()` with the measurement error variance `parameters$sigma2_me`
##'   (`draw_simulation()`, `assess_prediction()`).
##' - `score_support(draws, y, units_m, means, parameters)`, the counts over
##'   which a discrete prediction of a record is scored, from the draws of its
##'   response and their means, NULL for a family that is not discrete.
##'
##' The negative binomial family hands the engine its records as binomial ones
##' (`engine_records()`), so the methods that start a fit, `start_coefficients`,
##' `empirical_link`, `transform_variance` and `start_log_likelihood`, are the
##' binomial family's on those records, and its own are NULL. It has two more.
##' - `engine_records(y, units_m, offset, parameters)`, the records of the
##'   binomial kernel, with `y` events of `y + theta` trials and the offset
##'   `offset + log(units_m / theta)`, for exposures `units_m`.
##' - `log_likelihood_constant(y, parameters)`, the terms of the
##'   log-likelihood that the binomial kernel leaves out,
##'   `sum(lgamma(y + theta) - lgamma(theta))`. With them the log-likelihood
##'   tends to the Poisson family's as `theta` grows, which also leaves out
##'   `-lgamma(y + 1)`.
##' @noRd
family_traits <- function(family) {
  if (!is.character(family) || length(family) != 1L ||
      !family %in% names(riskmap_families)) {
    riskmap_bug("{.arg family} must be one of {.or {.val {names(riskmap_families)}}}.")
  }
  methods <- switch(
    family,
    gaussian = list(
      validate = function(y, units_m) invisible(TRUE),
      start_coefficients = function(y, units_m, D, cov_offset, link, formula,
                                    data) {
        linear_start_values(D, y)
      },
      empirical_link = NULL, transform_variance = NULL,
      start_log_likelihood = NULL,
      fitted_mean = function(mean, units_m) mean,
      simulate = function(means, units_m, parameters) {
        as.numeric(means) + rnorm(length(means),
                                  sd = sqrt(parameters$sigma2_me))
      },
      score_support = function(draws, y, units_m, means, parameters) NULL
    ),
    binomial = list(
      validate = function(y, units_m) check_binomial(y, units_m),
      start_coefficients = function(y, units_m, D, cov_offset, link, formula,
                                    data) {
        glm_fitted <- binomial_start_regression(y, units_m, D, cov_offset,
                                                link = link)
        coef(glm_fitted)
      },
      empirical_link = function(y, units_m) {
        log((y + 0.5) / (units_m - y + 0.5))
      },
      transform_variance = function(y, units_m) {
        median(1 / (y + 0.5) + 1 / (units_m - y + 0.5))
      },
      start_log_likelihood = function(eta, y, units_m) {
        sum(y * eta - units_m * softplus(eta))
      },
      fitted_mean = function(mean, units_m) mean,
      simulate = function(means, units_m, parameters) {
        rbinom(length(means), size = units_m, prob = as.numeric(means))
      },
      score_support = function(draws, y, units_m, means, parameters) {
        0:units_m
      }
    ),
    poisson = list(
      validate = function(y, units_m) invisible(TRUE),
      start_coefficients = function(y, units_m, D, cov_offset, link, formula,
                                    data) {
        pf_aux <- update(formula, . ~ . + offset(log(units_m)) +
                           offset(cov_offset))
        data_aux <- data
        data_aux$units_m <- units_m; data_aux$cov_offset <- cov_offset
        glm_fitted <- glm(pf_aux, data = data_aux, family = poisson)
        coef(glm_fitted)
      },
      empirical_link = function(y, units_m) log((y + 0.5) / units_m),
      transform_variance = function(y, units_m) median(1 / (y + 0.5)),
      start_log_likelihood = function(eta, y, units_m) {
        means <- exp(eta)
        if (any(!is.finite(means))) {
          return(NULL)
        }
        sum(y * eta - units_m * means)
      },
      fitted_mean = function(mean, units_m) units_m * mean,
      simulate = function(means, units_m, parameters) {
        rpois(length(means), lambda = as.numeric(means) * units_m)
      },
      score_support = function(draws, y, units_m, means, parameters) {
        lambda <- units_m * means
        0:max(max(draws), y, qpois(0.999, mean(lambda)))
      }
    ),
    negative_binomial = list(
      validate = function(y, units_m) check_counts(y, units_m),
      start_coefficients = NULL, empirical_link = NULL,
      transform_variance = NULL, start_log_likelihood = NULL,
      fitted_mean = function(mean, units_m) units_m * mean,
      simulate = function(means, units_m, parameters) {
        rnbinom(length(means), size = family_size(parameters),
                mu = as.numeric(means) * units_m)
      },
      score_support = function(draws, y, units_m, means, parameters) {
        lambda <- units_m * means
        0:max(max(draws), y, qnbinom(0.999, size = family_size(parameters),
                                     mu = mean(lambda)))
      },
      engine_records = function(y, units_m, offset, parameters) {
        theta <- family_size(parameters)
        units_m <- units_m %||% rep(1, length(y))
        offset <- offset %||% numeric(length(y))
        list(y = y, units_m = y + theta, offset = offset + log(units_m / theta))
      },
      log_likelihood_constant = function(y, parameters) {
        theta <- family_size(parameters)
        sum(lgamma(y + theta) - lgamma(theta))
      }
    )
  )
  c(list(family = family), riskmap_families[[family]], methods)
}

##' The size theta of a negative binomial family, from its `parameters`
##'
##' Stops if `parameters` has no single positive finite `theta`, which
##' `negative_binomial()` checks when a user gives it, so that a caller that
##' lost the family's parameters fails rather than fitting another model.
##' @noRd
family_size <- function(parameters) {
  theta <- parameters$theta
  if (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta) ||
      theta <= 0) {
    riskmap_bug("The negative binomial family needs its size {.code theta} in {.arg parameters}.")
  }
  theta
}

##' Check counts and their exposures
##'
##' Stops unless the responses `y` are non-negative whole numbers and the
##' exposures `units_m` positive finite numbers, as the negative binomial
##' family needs, whose records have the offset `log(units_m / theta)`.
##' @noRd
check_counts <- function(y, units_m) {
  tolerance <- sqrt(.Machine$double.eps)
  if (!all(y >= 0) || !all(abs(y - round(y)) < tolerance)) {
    abort_input("The responses of a negative binomial model must be non-negative whole numbers.", arg = "data")
  }
  if (!all(is.finite(units_m)) || !all(units_m > 0)) {
    abort_input("The denominators of a negative binomial model must be positive finite numbers.", arg = "denominator")
  }
  invisible(TRUE)
}

##' The records that the MCML engine fits, and the family it fits them as
##'
##' For a family with a `kernel_family`, its records as that family's
##' (`family_traits()`), and the kernel family. For the negative binomial
##' family with log link and size theta, the log-likelihood of a count y with
##' mean m exp(eta) is
##' y eta' - (y + theta) log(1 + exp(eta')) + lgamma(y + theta) -
##' lgamma(theta) - lgamma(y + 1), with eta' = eta + log(m / theta),
##' which is the binomial logit kernel of y events of y + theta trials with
##' the offset log(m / theta) (#253, step 6a). So the engine fits these
##' records with every part of the canonical binomial path, compiled weights,
##' tables of integrated records and checks, and no code of its own. The
##' constant is added to the reported log-likelihood
##' (`log_likelihood_constant()`). `offset` is the offset of the records, or
##' their linear predictors where a caller passes them. A family the engine
##' fits as itself returns its records unchanged. `invlink` must be NULL for a
##' family with a kernel family, whose link is the log, and `parameters` are
##' the family's own (`negative_binomial()`).
##' @return A list with `family`, `y`, `units_m` and `offset`.
##' @noRd
engine_records <- function(family, y, units_m, offset, parameters = NULL,
                           invlink = NULL) {
  traits <- family_traits(family)
  if (is.null(traits$kernel_family)) {
    return(list(family = family, y = y, units_m = units_m, offset = offset))
  }
  if (!is.null(invlink)) {
    riskmap_bug("The {family} family has the log link only, and takes no {.arg invlink}.")
  }
  c(list(family = traits$kernel_family),
    traits$engine_records(y, units_m, offset, parameters))
}

##' The terms of a model's log-likelihood that its engine records leave out
##'
##' Zero for a family the engine fits as itself, and otherwise the family's
##' `log_likelihood_constant()` (`engine_records()`).
##' @noRd
log_likelihood_constant <- function(family, y, parameters = NULL) {
  traits <- family_traits(family)
  if (is.null(traits$kernel_family)) return(0)
  traits$log_likelihood_constant(y, parameters)
}

##' The families that the MCML engine fits, those whose likelihood is not exact
##' @noRd
engine_families <- function() {
  exact <- vapply(riskmap_families, function(entry) entry$exact, logical(1))
  names(riskmap_families)[!exact]
}

##' The traits of a family that the MCML engine fits
##'
##' Stops unless `family` is one string naming a family whose likelihood is not
##' exact, binomial or Poisson, and otherwise returns its traits
##' (`family_traits()`). For the engine's internal functions, whose callers
##' have read the family already, so a failure is RiskMap's.
##' @noRd
check_engine_family <- function(family) {
  engine_families <- engine_families()
  if (!is.character(family) || length(family) != 1L ||
      !family %in% engine_families) {
    riskmap_bug("{.arg family} must be {.or {.val {engine_families}}}.")
  }
  family_traits(family)
}

##' Whether a fitted model's family is exact, fitted in closed form
##'
##' TRUE for a Gaussian fit, whose prediction conditions in closed form and
##' which has no Monte Carlo error, and FALSE for a binomial or Poisson fit. A
##' fit whose family is not one of RiskMap's, such as one whose `family` was
##' removed, stops with an input error, so that it is not predicted or checked
##' as if it were a binomial or Poisson fit.
##' @noRd
fit_is_exact <- function(fit) {
  traits <- fit_traits(fit)
  if (is.null(traits)) {
    families <- names(riskmap_families)
    abort_input("{.arg object} must be a fit of the family {.or {.val {families}}}, but its {.field family} is {.obj_type_friendly {fit$family}}.", arg = "object")
  }
  traits$exact
}

##' The traits of the family of a fitted model
##'
##' `family_traits()` of the family string the fit stores, or NULL if it is
##' not one of RiskMap's families.
##' @noRd
fit_traits <- function(fit) {
  family <- fit$family
  if (!is.character(family) || length(family) != 1L ||
      !family %in% names(riskmap_families)) {
    return(NULL)
  }
  family_traits(family)
}

##' Read the response family and link that a user gives
##'
##' The one reader of `family` for every exported function that takes it.
##' `family` is one of the strings `"gaussian"`, `"binomial"` and
##' `"poisson"`, a family function such as `binomial`, or a family object
##' such as `binomial("probit")`. A family object becomes its family string
##' and, for a link other than the canonical one, the list form of `invlink`
##' that `make_inverse_link_functions()` reads, with the link's name and its
##' exact inverse and derivatives (`named_inverse_link()`). The canonical
##' logit and log links keep `invlink = NULL`, so they take the canonical path
##' and fit exactly as the strings do, and the identity link of the Gaussian
##' family leaves nothing to carry. A string leaves `invlink` as it is, which
##' keeps the internal list form open to refits. A family with parameters of
##' its own is given only by its function, `negative_binomial()`, whose
##' specification is accepted by every caller and adds the family's
##' `parameters`. It is internal until #253, step 6b, so it is not among the
##' families that errors name.
##' @param family What the user gave as `family`.
##' @param invlink What the user gave as `invlink`.
##' @param families The family strings the calling function accepts.
##' @return A list with `family`, one of `families` or a family given by its
##'   specification, and `invlink`, with `parameters` for a family that has
##'   them.
##' @noRd
read_family <- function(family, invlink = NULL,
                        families = c("gaussian", "binomial", "poisson")) {
  if (inherits(invlink, "link-glm")) {
    link_name <- invlink$name
    abort_input(c("{.arg invlink} must be an inverse-link function, not a {.cls link-glm} object.", "i" = "Give a named link in {.arg family}, as in {.code binomial(\"{link_name}\")}, or the inverse alone, as in {.code make.link(\"{link_name}\")$linkinv}."), arg = "invlink")
  }
  if (is.character(invlink)) {
    abort_input(c("{.arg invlink} must be a function, not the name of a link.", "i" = "Give a named link in {.arg family}, as in {.code binomial(\"probit\")}."), arg = "invlink")
  }
  if (inherits(family, "riskmap_family_specification")) {
    if (!is.null(invlink)) {
      abort_input(family_traits(family$family)$link_hint, arg = "invlink")
    }
    return(list(family = family$family, invlink = NULL,
                parameters = family$parameters))
  }
  # A family with parameters of its own is not given as a string.
  families <- families[vapply(families, function(name) {
    !length(family_traits(name)$parameters)
  }, logical(1))]
  if (is.function(family)) {
    family <- tryCatch(family(), error = function(e) NULL)
  }
  if (inherits(family, "family")) {
    link <- read_family_link(family$family %||% "", family$link %||% "")
    if (!is.null(link$invlink) && !is.null(invlink)) {
      family_call <- paste0(family$family, "(\"", family$link, "\")")
      abort_input("Give the link either in {.arg family}, as in {.code {family_call}}, or as a function in {.arg invlink}, not both.", arg = "invlink")
    }
    family <- link$family
    if (is.null(invlink)) invlink <- link$invlink
  }
  if (!is.character(family) || length(family) != 1L ||
      !family %in% families) {
    abort_input("{.arg family} must be {.or {.val {families}}}, or a family object such as {.code binomial(\"probit\")}.", arg = "family")
  }
  list(family = family, invlink = invlink)
}

##' The negative binomial family at a known size
##'
##' Counts y with mean m exp(eta), for exposure m, and variance
##' mu + mu^2 / theta, as in `MASS::negative.binomial(theta)`, with the size
##' `theta` held at the value given. Internal until the size can be estimated
##' (#253, step 6b), when it is to be exported with `theta` estimated by
##' default. `read_family()` reads it.
##' @param theta The size, a positive number.
##' @return A specification of class `riskmap_family_specification`.
##' @noRd
negative_binomial <- function(theta) {
  check_supplied(theta)
  check_positive_number(theta, type = "")
  structure(list(family = "negative_binomial", link = "log",
                 parameters = list(theta = theta)),
            class = "riskmap_family_specification")
}

##' The `family` argument that refits a fitted model
##'
##' The family string the fit stores, or for a family with parameters of its
##' own its specification at the fit's values (`negative_binomial()`), for the
##' exported functions that refit a fit, such as `assess_prediction()`.
##' @noRd
fit_family_argument <- function(fit) {
  traits <- fit_traits(fit)
  if (is.null(traits) || !length(traits$parameters)) return(fit$family)
  negative_binomial(theta = family_size(fit$family_parameters))
}

##' Family string and link of a family object
##'
##' Takes the `family` and `link` names of a stats family object and returns
##' its family string and `invlink`, or stops with the reason the link is not
##' accepted and what to give instead, from the links of `riskmap_families`.
##' @noRd
read_family_link <- function(family, link) {
  object <- paste0(family, "(\"", link, "\")")
  # A template, interpolated by the error that uses it.
  use_invlink <- "Other links can be given as an inverse-link function in {.arg invlink}, with {.code family = \"{family}\"}."
  if (startsWith(family, "quasi")) {
    abort_input(c("{.fn {family}} is a quasi-likelihood family, which has no likelihood for RiskMap to maximise.", "i" = "Use {.fn binomial} or {.fn poisson}, and give a link they don't name as an inverse-link function in {.arg invlink}."), arg = "family")
  }
  if (!family %in% names(riskmap_families)) {
    abort_input(c("{.arg family} must be a Gaussian, binomial or Poisson family, not {.val {family}}.", "i" = "For a binomial or Poisson response, any link can be given as an inverse-link function in {.arg invlink}."), arg = "family")
  }
  traits <- family_traits(family)
  if (length(traits$parameters)) {
    riskmap_bug("The {family} family has parameters of its own, so it is given by its own function, not as a stats family.")
  }
  refused <- if (link %in% names(traits$refused)) traits$refused[[link]]
  if (!is.null(refused)) {
    abort_input(c("{.code {object}} can't be used, because {refused}.", "i" = use_invlink), arg = "family")
  }
  accepted <- traits$links
  if (!link %in% accepted) {
    hint <- traits$link_hint %||% use_invlink
    abort_input(c("The {link} link isn't available for the {family} family, which takes {.or {.val {accepted}}}.", "i" = hint), arg = "family")
  }
  canonical <- identical(link, traits$default_link)
  list(family = family,
       invlink = if (!canonical) c(named_inverse_link(link), link = link))
}

##' The family of a model and the likelihood of its records
##'
##' One object through which the MCML engine sees the family and the link of a
##' model (#253). `family` is "gaussian", "binomial" or "poisson", and `linkf`
##' the inverse link functions of `make_inverse_link_functions()`, NULL for the
##' canonical link. The object has every field and method of
##' `family_traits()`. The Gaussian family is exact (`exact`), fitted without
##' the engine, so its object has only its traits, with no link, properties
##' of the engine or likelihood methods, and takes no `linkf` or `kernel`.
##' `check_family_object()`, through which the engine reads a family, refuses
##' it. Where records have effects
##' integrated out, `kernel` is their kernel at the variances of the classes
##' (`make_integrated_kernel()`), which gives
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
##'
##' Methods, for records with `y` events of `units_m` trials, or counts with
##' exposure `units_m`, at linear predictors `eta`.
##' - `terms(eta, y, units_m)`, the log-likelihood, score and curvature, or NULL
##'   where they are not defined, for the mode search (`find_latent_mode()`).
##' - `log_likelihood(eta, y, units_m)`, the log-likelihood of each row of a
##'   matrix of linear predictors, one column per record
##'   (`compute_log_weights_in_r()`).
##' - `slope_parts(eta, y, units_m)`, the score and curvature of each record
##'   at each row of `eta` under the model's link, for the slope at a held
##'   variance, or NULL where a custom link cannot give them
##'   (`compute_slope_terms()`). Records integrated out take theirs from their
##'   table instead.
##' - `move_start(latent, linear_predictor, y, units_m, columns)`, as
##'   `move_start_off_limits()`.
##' - `nonlinear_parts(centre, step, y, units_m)`, the part of each record's
##'   log-likelihood that is not linear in a step from `centre`, with the slope
##'   at `centre`, for the global proposal (`update_observation_factors()`).
##' - `kernel_arrays()`, the arrays of a custom link's kernel for the site
##'   table (`compute_custom_kernel_arrays()`), NULL for the canonical link.
##' @noRd
make_family <- function(family, linkf = NULL, kernel = NULL) {
  traits <- family_traits(family)
  if (!is.null(traits$kernel_family)) {
    riskmap_bug("The engine fits the records of the {family} family as {traits$kernel_family} ones ({.fn engine_records}), so it has no family object of its own.")
  }
  if (traits$exact) {
    if (!is.null(linkf) || !is.null(kernel)) {
      riskmap_bug("The {family} family is exact and takes no {.arg linkf} or {.arg kernel}.")
    }
    return(structure(c(traits, list(linkf = NULL, kernel = NULL)),
                     class = "riskmap_family"))
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
  nonlinear_parts <- if (integrated) {
    function(centre, step, y, units_m) kernel$nonlinear_parts(centre, step)
  } else if (canonical && poisson) {
    function(centre, step, y, units_m) {
      canonical_poisson_nonlinear_parts(centre, step, y, units_m)
    }
  } else if (canonical) {
    function(centre, step, y, units_m) {
      canonical_binomial_nonlinear_parts(centre, step, y, units_m)
    }
  } else {
    function(centre, step, y, units_m) {
      link_nonlinear_parts(centre, step, y, units_m, family, linkf, terms)
    }
  }
  kernel_arrays <- if (canonical) {
    function() NULL
  } else {
    function() compute_custom_kernel_arrays(linkf, family)
  }
  structure(
    c(traits, list(
      linkf = linkf, kernel = kernel, canonical = canonical,
      integrated = integrated, poisson = poisson, compiled = canonical,
      assumes_concave = canonical && !integrated,
      moves_start = !canonical && !integrated,
      needs_finite_weight = !canonical,
      terms = terms, log_likelihood = log_likelihood,
      slope_parts = slope_parts, nonlinear_parts = nonlinear_parts,
      kernel_arrays = kernel_arrays,
      move_start = function(latent, linear_predictor, y, units_m, columns) {
        move_start_off_limits(latent, linear_predictor, y, units_m, family,
                              linkf, columns)
      }
    )),
    class = "riskmap_family"
  )
}

##' Stops if a family object is of an exact family, which the engine does not
##' fit, and otherwise returns it.
##' @noRd
refuse_exact_family <- function(object) {
  if (isTRUE(object$exact)) {
    riskmap_bug("The {object$family} family is fitted exactly, not by the engine, which can't read its family object.")
  }
  object
}

##' The family object that an engine function takes in its `family` argument
##'
##' Stops unless `family` is a family object (`make_family()`) that the engine
##' can read, which an exact family's is not (`refuse_exact_family()`), and
##' otherwise returns it. Since step 1d of #253 the engine's functions take the
##' object only, not a family string with a link.
##' @noRd
check_family_object <- function(family) {
  if (!inherits(family, "riskmap_family")) {
    riskmap_bug("{.arg family} must be a family object from {.fn make_family}, not {.obj_type_friendly {family}}.")
  }
  refuse_exact_family(family)
}

##' The family object of a fitted model
##'
##' Rebuilt from what the fit stores, the family string and the
##' `link_function` list, which is NULL for a Gaussian fit. It gives the
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
##' The `terms()` of the family object for the canonical Poisson link,
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
##' The `terms()` of the family object for the canonical binomial link,
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
##' about -745. There the formulas of `canonical_binomial_terms()`,
##' `canonical_poisson_terms()` and `maxim_integrand()` divide zero by zero,
##' although the likelihood is smooth and finite. Complete separation by a
##' covariate takes the linear predictor
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
##' below -372 for a positive count, and `find_latent_mode()` rejects a
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
##' The `log_likelihood()` of the family object for a custom link. A draw
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

##' Nonlinear part of each record's log-likelihood, canonical binomial link
##'
##' The `nonlinear_parts()` of the family object for the canonical binomial link. The
##' difference is written with log1p and expm1 so that nearly equal
##' log-likelihood values are never subtracted, and for large steps the
##' logarithms of p and 1 - p are taken directly from the linear predictor,
##' because p itself rounds to 0 or 1 when the predictor is extreme.
##' @noRd
canonical_binomial_nonlinear_parts <- function(centre, step, y, units_m) {
  probability <- stats::plogis(centre)
  slope <- y - units_m * probability
  # log(1 - p + p exp(step)), switching to log-sum-exp for large steps.
  log_success <- stats::plogis(centre, log.p = TRUE) + step
  log_failure <- matrix(
    stats::plogis(centre, lower.tail = FALSE, log.p = TRUE),
    nrow(step), ncol(step)
  )
  log_ratio <- ifelse(
    abs(step) < 20,
    log1p(probability * expm1(step)),
    pmax(log_failure, log_success) +
      log1p(exp(-abs(log_failure - log_success)))
  )
  nonlinear <- units_m * (probability * step - log_ratio)
  list(nonlinear = nonlinear, slope = slope)
}

##' Nonlinear part of each record's log-likelihood, canonical Poisson link
##'
##' The `nonlinear_parts()` of the family object for the canonical Poisson link.
##' @noRd
canonical_poisson_nonlinear_parts <- function(centre, step, y, units_m) {
  rate <- units_m * exp(centre)
  slope <- y - rate
  nonlinear <- -rate * (expm1(step) - step)
  list(nonlinear = nonlinear, slope = slope)
}

##' Nonlinear part of each record's log-likelihood, other links
##'
##' The `nonlinear_parts()` of the family object for a custom link, the direct
##' difference, with the slope from the observation terms
##' `observation_terms` of the family object.
##' @noRd
link_nonlinear_parts <- function(centre, step, y, units_m, family, linkf,
                                 observation_terms) {
  log_likelihood <- function(eta) {
    mean <- linkf$inv(eta)
    if (family == "poisson") {
      y * log(pmax(units_m * mean, .Machine$double.xmin)) - units_m * mean
    } else {
      log_probability <- bound_log_probability(mean)
      y * log_probability$success + (units_m - y) * log_probability$failure
    }
  }
  terms <- observation_terms(centre, y, units_m)
  if (is.null(terms)) {
    riskmap_abort("The linear predictor at the centre of the global proposal is invalid.", class = "riskmap_error_fit")
  }
  per_observation_centre <- log_likelihood(centre)
  nonlinear <- matrix(log_likelihood(as.vector(centre + step)), nrow(step)) -
    per_observation_centre - step * terms$score
  list(nonlinear = nonlinear, slope = terms$score)
}

##' Score and curvature of each record at each row, canonical link
##'
##' The closed forms of the `slope_parts()` of the family object for the
##' canonical links, one row per draw and one column per record.
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
##' The `slope_parts()` of the family object for a custom link, through
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

##' The kernel of a custom link on a lattice of linear predictors
##'
##' Each record's log-likelihood is y A(eta) + m B(eta) + c, with A and B the
##' same for every record, so one set of arrays serves all of them. For the
##' binomial family A = log p - log(1 - p) and B = log(1 - p), and for the
##' Poisson family A = log(mean) and B = -mean (with c = y log m). The
##' arrays hold A, B and their first two derivatives from the link's own
##' derivatives, on the lattice from `low` to `high` with step 0.25 / 128, so
##' the coarsest steps of the table's rule fall on lattice points. As in the
##' sampled form (`bound_log_probability()`), a probability that rounds to a
##' limit is charged at the bound of 2^-53, where that logarithm is flat. Where
##' the link is not finite, beyond the range on which it was checked, the
##' kernel is continued from the nearest point where it is.
##' @noRd
compute_custom_kernel_arrays <- function(linkf, family, low = -100,
                                         high = 80, step = 0.25 / 128) {
  eta <- low + step * (0:round((high - low) / step))
  mean <- suppressWarnings(linkf$inv(eta))
  first <- suppressWarnings(linkf$d1(eta))
  second <- suppressWarnings(linkf$d2(eta))
  invalid <- !is.finite(mean) | !is.finite(first) | !is.finite(second) |
    mean < 0 | (family == "binomial" & mean > 1)
  if (all(invalid)) {
    abort_input("{.arg invlink} gives no finite mean and derivatives at any of the points tried.", arg = "invlink")
  }
  if (any(invalid)) {
    valid <- which(!invalid)
    nearest <- valid[pmax(1L, findInterval(seq_along(eta), valid))]
    below <- seq_along(eta) < valid[1L]
    nearest[below] <- valid[1L]
    mean[invalid] <- mean[nearest[invalid]]
    first[invalid] <- 0
    second[invalid] <- 0
  }
  if (family == "binomial") {
    bound <- 2^-53
    success <- mean > bound
    failure <- 1 - mean > bound
    log_success <- log(pmax(mean, bound))
    log_failure <- log(pmax(1 - mean, bound))
    success_first <- ifelse(success, first / mean, 0)
    success_second <- ifelse(success, second / mean - (first / mean)^2, 0)
    failure_first <- ifelse(failure, -first / (1 - mean), 0)
    failure_second <- ifelse(failure,
                             -second / (1 - mean) - (first / (1 - mean))^2, 0)
    list(low = low, step = step, A = log_success - log_failure,
         A1 = success_first - failure_first,
         A2 = success_second - failure_second, B = log_failure,
         B1 = failure_first, B2 = failure_second)
  } else {
    positive <- mean > 0
    list(low = low, step = step, A = log(pmax(mean, .Machine$double.xmin)),
         A1 = ifelse(positive, first / mean, 0),
         A2 = ifelse(positive, second / mean - (first / mean)^2, 0),
         B = -mean, B1 = -first, B2 = -second)
  }
}

##' The largest term of each record of the likelihood over its linear predictor
##'
##' y log(y / m) + (m - y) log(1 - y / m) for a binomial record and
##' y log(y) - y for a Poisson one, with 0 log 0 = 0, which bound the terms of
##' the MCML engine for any link. The bounds of `bound_log_probability()` only
##' move a probability within [0, 1], or a mean above zero, so they hold there
##' too.
##' @noRd
compute_largest_record_terms <- function(y, units_m, family) {
  if (identical(family, "poisson")) {
    return(ifelse(y > 0, y * log(pmax(y, 1)) - y, 0))
  }
  proportion <- y / units_m
  ifelse(y > 0, y * log(pmax(proportion, .Machine$double.xmin)), 0) +
    ifelse(units_m > y,
           (units_m - y) * log1p(-pmin(proportion, 1 - .Machine$double.eps)),
           0)
}
