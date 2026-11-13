##' Check the likelihood, mean and simulator of a family on a grid
##'
##' The checks of step 5 of #253. The tests run them on every family of
##' `riskmap_families` with each of its links, so that a family or a link
##' added to the table is checked as soon as it is there, and `check_fit()`
##' runs them on a link given as a function to learn whether its likelihood
##' is convex anywhere. On records of the kinds the family meets
##' (`check_family_records()`), at the linear predictors `eta`, it checks that
##' - the inverse link gives a mean in the range of the family,
##' - the log-likelihood by which the engine weights draws, the
##'   `log_likelihood()` of the family object, is the family's log density
##'   (`log_density()`) plus a term free of the linear predictor, so that it is
##'   bounded above where the density is,
##' - the score and curvature of the mode search, the `terms()` of the family
##'   object, agree with central differences of the log-likelihood and of the
##'   score, with a relative `tolerance`,
##' - and the draws of `simulate()`, on a fixed seed, have the mean and
##'   variance of the family's density, to within five standard errors.
##'
##' A second difference of a record's log-likelihood over three neighbouring
##' points of the grid that is positive shows that it is convex somewhere
##' between them, from the log-likelihood alone, whereas a curvature of a link
##' given as a function rests on numerical derivatives that can lose its sign
##' where two terms cancel. The point of the largest such difference is
##' returned, for the first record that has one. A link given as a function
##' may compute its mean with a cancellation, as 1 - exp(-exp(eta)) does in
##' its lower tail, where the relative error grows as the mean falls, so
##' convexity is looked for only where the mean is at least 10^-6 inside its
##' range, and the difference must exceed 10^-8 times the size of the
##' log-likelihoods and counts involved, about ten times what such rounding
##' can give there. No grid can show that a likelihood is concave, since a t link
##' with 10,000 degrees of freedom looks concave from -50 to 30 while its
##' log F is convex below -101, so the links known to be concave are listed
##' in the table (`concave_links`). A probability given by its value loses the
##' precision of 1 - p as p nears one, which central differences would
##' multiply, and is bounded where it rounds (`bound_log_probability()`), so a
##' record with failures is compared only where 1 - p is above 10^-6, and one
##' with events only where its mean is above 10^-200. Negative binomial
##' records are those of the binomial kernel (`engine_records()`) with the
##' terms they leave out (`log_likelihood_constant()`), at the size in
##' `parameters`. The Gaussian family is exact, with no kernel, so only its
##' mean and simulator are checked, at the measurement error variance in
##' `parameters`. `linkf` is a link of `make_inverse_link_functions()`, NULL
##' for the family's own inverse link (`inverse_link`). The derivatives of the
##' built-in links are exact, which the `tolerance` of 10^-6 suits. Those of a
##' link given as a function come from `Deriv` or numerical differences,
##' accurate to about 10^-5 and much worse near a linear predictor of zero,
##' so `check_fit()` reads only `convex_at` for such a link.
##'
##' Returns a list with `problems`, a sentence for each check that fails,
##' empty if none does, and `convex_at`, NULL or the linear predictor at
##' which a record's log-likelihood is most convex, with the record's `y` and
##' `units_m`. The random number state of the session is left as it was.
##' @noRd
check_family <- function(family, linkf = NULL, parameters = NULL,
                         eta = seq(-30, 30, by = 0.25), step = 1e-4,
                         tolerance = 1e-6, draws = 20000L) {
  traits <- family_traits(family)
  inverse <- if (is.null(linkf)) traits$inverse_link else linkf$inv
  problems <- character()
  mean <- inverse(eta)
  range <- traits$mean_range
  outside <- is.na(mean) | mean < range[1L] | mean > range[2L]
  if (any(outside)) {
    problems <- c(problems, paste0(
      "The mean is outside the range of the family at a linear predictor of ",
      eta[which(outside)[1L]], "."
    ))
  }
  convex_at <- NULL
  if (!traits$exact) {
    records <- check_family_records(traits)
    for (record in split(records, seq_len(nrow(records)))) {
      checked <- check_family_record(traits, linkf, inverse, parameters,
                                     eta, record$y, record$units_m, step,
                                     tolerance)
      problems <- c(problems, checked$problems)
      if (is.null(convex_at) && !is.null(checked$convex_at)) {
        convex_at <- list(eta = checked$convex_at, y = record$y,
                          units_m = record$units_m)
      }
    }
  }
  problems <- c(problems, check_family_simulator(traits, inverse,
                                                 parameters, draws))
  list(problems = problems, convex_at = convex_at)
}

##' The records on which `check_family()` checks a family
##'
##' For binomial records, one trial with an event and one without, which give
##' log F and log(1 - F) alone, so that their concavity is that of every
##' record, and records of 20 trials. For counts, a zero count, which gives
##' minus the mean, a count of one with an exposure so small that it gives
##' nearly the log of the mean alone, and larger counts.
##' @noRd
check_family_records <- function(traits) {
  if (traits$upper_limit) {
    return(data.frame(y = c(1, 0, 3, 20, 0), units_m = c(1, 1, 20, 20, 20)))
  }
  data.frame(y = c(0, 1, 1, 4, 30), units_m = c(1, 1e-6, 1, 2, 0.5))
}

##' The checks of `check_family()` on one record
##'
##' Each linear predictor of `eta` is a record of its own with `y` events of
##' `units_m` trials, or a count with exposure `units_m`. `inverse` is the
##' inverse link. Returns its `problems` and `convex_at`, the linear predictor
##' where the second difference of its log-likelihood shows it most convex, or
##' NULL.
##' @noRd
check_family_record <- function(traits, linkf, inverse, parameters, eta, y,
                                units_m, step, tolerance) {
  # The engine sees a negative binomial record as a binomial one.
  kernel_family <- traits$kernel_family %||% traits$family
  object <- make_family(kernel_family,
                        if (is.null(traits$kernel_family)) linkf)
  describe <- paste0("for a record of ", y, " with ", units_m,
                     if (traits$upper_limit) " trials" else " exposure")
  kernel_at <- function(eta) {
    count <- length(eta)
    records <- if (is.null(traits$engine_records)) {
      list(y = y, units_m = units_m, offset = eta)
    } else {
      traits$engine_records(rep(y, count), rep(units_m, count), eta,
                            parameters)
    }
    constant <- if (is.null(traits$log_likelihood_constant)) {
      0
    } else {
      traits$log_likelihood_constant(y, units_m, parameters)
    }
    log_likelihood <- object$log_likelihood(
      matrix(records$offset, ncol = 1L), records$y[1L], records$units_m[1L]
    ) + constant
    terms <- object$terms(records$offset, rep(records$y[1L], count),
                          rep(records$units_m[1L], count))
    list(log_likelihood = as.numeric(log_likelihood),
         score = terms$score, curvature = terms$curvature)
  }
  centre <- kernel_at(eta)
  if (is.null(centre$score)) {
    undefined <- Filter(function(at) is.null(kernel_at(at)$score), eta)
    return(list(problems = paste0(
      "The score or curvature is not defined at a linear predictor of ",
      undefined[1L], ", ", describe, "."
    ), convex_at = NULL))
  }
  lower <- kernel_at(eta - step)
  upper <- kernel_at(eta + step)
  usable <- function(eta) {
    mean <- inverse(eta)
    (y == 0 | mean > 1e-200) &
      (!traits$upper_limit | y == units_m | 1 - mean > 1e-6)
  }
  used <- usable(eta - step) & usable(eta) & usable(eta + step) &
    is.finite(centre$log_likelihood) & is.finite(lower$score) &
    is.finite(upper$score)
  problems <- character()
  first_failure <- function(failed, what) {
    if (!any(failed)) return(character())
    paste0(what, " at a linear predictor of ", eta[which(failed)[1L]], ", ",
           describe, ".")
  }
  if (any(is.nan(centre$log_likelihood) | centre$log_likelihood == Inf)) {
    problems <- c(problems, first_failure(
      is.nan(centre$log_likelihood) | centre$log_likelihood == Inf,
      "The log-likelihood is not a number or is infinite"
    ))
  }
  if (!any(used)) {
    return(list(problems = c(problems, paste0(
      "No linear predictor of the grid can be checked, ", describe, "."
    )), convex_at = NULL))
  }
  mean <- inverse(eta[used])
  density <- traits$log_density(y, mean, units_m, parameters)
  difference <- centre$log_likelihood[used] - density
  scale <- 1 + abs(centre$log_likelihood[used])
  shifted <- abs(difference - difference[1L]) > 1e-8 * scale
  problems <- c(problems, first_failure(
    replace(logical(length(eta)), which(used), shifted),
    "The log-likelihood is not the log density plus a term free of the linear predictor"
  ))
  difference_score <- (upper$log_likelihood - lower$log_likelihood) /
    (2 * step)
  difference_curvature <- -(upper$score - lower$score) / (2 * step)
  score_off <- used & abs(centre$score - difference_score) >
    tolerance * (1 + abs(centre$score))
  curvature_off <- used & abs(centre$curvature - difference_curvature) >
    tolerance * (1 + abs(centre$curvature))
  problems <- c(problems,
                first_failure(score_off, "The score differs from a central difference of the log-likelihood"),
                first_failure(curvature_off, "The curvature differs from a central difference of the score"))
  list(problems = problems,
       convex_at = find_convex_point(traits, eta, centre$log_likelihood,
                                     used, inverse, y, units_m))
}

##' Where the second difference of a record's log-likelihood shows it convex
##'
##' Over each three neighbouring points of `eta` that are `used` and whose
##' mean is at least 10^-6 inside the range of the family, the second divided
##' difference of `log_likelihood`, positive beyond 10^-8 times the size of
##' the log-likelihoods and counts involved. Returns the middle point of the
##' largest, or NULL if there is none.
##' @noRd
find_convex_point <- function(traits, eta, log_likelihood, used, inverse, y,
                              units_m) {
  count <- length(eta)
  if (count < 3L) return(NULL)
  mean <- inverse(eta)
  inside <- used & mean >= traits$mean_range[1L] + 1e-6 &
    mean <= traits$mean_range[2L] - 1e-6
  before <- seq_len(count - 2L)
  middle <- before + 1L
  after <- before + 2L
  rise <- (log_likelihood[after] - log_likelihood[middle]) /
    (eta[after] - eta[middle]) -
    (log_likelihood[middle] - log_likelihood[before]) /
    (eta[middle] - eta[before])
  second <- 2 * rise / (eta[after] - eta[before])
  size <- abs(log_likelihood[before]) + abs(log_likelihood[middle]) +
    abs(log_likelihood[after]) + 1 + y + units_m
  spacing <- (eta[after] - eta[before]) / 2
  convex <- inside[before] & inside[middle] & inside[after] &
    second * spacing^2 > 1e-8 * size
  convex[is.na(convex)] <- FALSE
  if (!any(convex)) return(NULL)
  eta[middle][convex][which.max(second[convex])]
}

##' The check of a family's simulator in `check_family()`
##'
##' At linear predictors of -1, 0 and 1, with 20 trials or an exposure of 2,
##' `draws` responses from `simulate()` on a fixed seed are compared with the
##' mean and variance of the family's density, found by summing it over the
##' support of a discrete family (`score_support()`) and by integrating it for
##' a continuous one. Returns a sentence for each linear predictor where
##' either is more than five standard errors away. The random number state of
##' the session is left as it was.
##' @noRd
check_family_simulator <- function(traits, inverse, parameters, draws) {
  restore <- preserve_random_seed()
  on.exit(restore(), add = TRUE)
  set.seed(1L)
  units_m <- if (traits$whole_denominator) 20 else 2
  problems <- character()
  for (eta in c(-1, 0, 1)) {
    mean <- inverse(eta)
    simulated <- traits$simulate(rep(mean, draws), rep(units_m, draws),
                                 parameters)
    moments <- if (traits$discrete) {
      support <- traits$score_support(simulated, 0, units_m, mean, parameters)
      weights <- exp(traits$log_density(support, mean, units_m, parameters))
      centred <- support - sum(support * weights)
      c(sum(support * weights), sum(centred^2 * weights),
        sum(centred^4 * weights))
    } else {
      integral <- function(f) {
        stats::integrate(function(y) {
          f(y) * exp(traits$log_density(y, mean, units_m, parameters))
        }, -Inf, Inf, rel.tol = 1e-10)$value
      }
      first <- integral(identity)
      c(first, integral(function(y) (y - first)^2),
        integral(function(y) (y - first)^4))
    }
    mean_error <- (base::mean(simulated) - moments[1L]) /
      sqrt(moments[2L] / draws)
    variance_error <- (stats::var(simulated) - moments[2L]) /
      sqrt((moments[3L] - moments[2L]^2) / draws)
    if (abs(mean_error) > 5 || abs(variance_error) > 5) {
      problems <- c(problems, paste0(
        "The draws of the simulator do not have the mean and variance of the ",
        "density at a linear predictor of ", eta, "."
      ))
    }
  }
  problems
}
