##' Control settings for the optimiser in Monte Carlo maximum likelihood
##'
##' `nlminb()` stops when the predicted relative decrease of the objective is
##' below `rel.tol`. Near a maximum a loss of D in log-likelihood corresponds
##' to a distance of sqrt(2 D) statistical standard errors whatever the size
##' of the likelihood, so the tolerance is set from `objective_scale`, the
##' size of the objective at the start, to stop within about 1e-4
##' log-likelihood units, or 0.014 standard errors, of the maximum. That is
##' far below the certificate's tolerance, and the Newton correction with the
##' validation score removes most of it. Without a scale the tolerance is
##' 1e-10.
##' @noRd
make_optimiser_control <- function(objective_scale = NULL) {
  relative <- if (is.null(objective_scale) || !is.finite(objective_scale) ||
                  objective_scale <= 0) {
    1e-10
  } else {
    min(1e-6, max(1e-12, 1e-4 / objective_scale))
  }
  list(eval.max = 600, iter.max = 300, rel.tol = relative, x.tol = 1e-8)
}

##' Effective sample size and tail diagnostics of importance weights
##' @noRd
summarise_importance_weights <- function(log_weights) {
  valid <- is.finite(log_weights)
  if (!any(valid)) {
    return(list(log_likelihood = -Inf, relative_ess = 0, maximum_weight = 1,
                pareto_k = Inf))
  }
  largest <- max(log_weights[valid])
  relative <- exp(log_weights - largest)
  relative[!valid] <- 0
  weights <- relative / sum(relative)
  # The Pareto tail diagnostic is only a diagnostic, so a bank too small or
  # degenerate for it records NA instead of stopping the fit.
  pareto_k <- tryCatch(
    as.numeric(loo::pareto_k_values(
      suppressWarnings(loo::psis(log_weights[valid]))
    )),
    error = function(error) NA_real_
  )
  list(
    log_likelihood = largest + log(mean(relative)),
    relative_ess = 1 / sum(weights^2) / length(weights),
    maximum_weight = max(weights),
    pareto_k = pareto_k
  )
}

##' Score of the estimated log-likelihood on an independent validation bank
##'
##' The score is the importance-weighted mean of the pathwise derivatives of
##' the log weights, by central differences along the moving proposal, and its
##' Monte Carlo covariance is estimated from the same draws, for the
##' parameters indexed by `free`. The bank is drawn from `seed` in blocks of
##' rows, so that a large bank is never held in memory at once, and each
##' proposal is built once and reused for every block. A penalty on the
##' parameters does not change the weights, so its gradient is added to the
##' score. Returns NULL if a proposal or a derivative cannot be computed.
##' @noRd
score_validation_bank <- function(engine, parameters, free, n_draws, seed,
                             penalty = NULL, block_rows = NULL) {
  started <- proc.time()[["elapsed"]]
  dimension <- engine$latent_dimension
  if (is.null(block_rows)) block_rows <- max(250L, floor(2.5e6 / dimension))
  step <- .Machine$double.eps^(1 / 3) * pmax(1, abs(parameters[free]))
  shifted <- list(parameters)
  for (parameter in seq_along(free)) {
    upper <- lower <- parameters
    upper[free[parameter]] <- upper[free[parameter]] + step[parameter]
    lower[free[parameter]] <- lower[free[parameter]] - step[parameter]
    shifted <- c(shifted, list(upper, lower))
  }
  proposals <- lapply(shifted, function(values) {
    tryCatch(engine$proposal_at(values), error = function(error) NULL)
  })
  if (any(vapply(proposals, is.null, logical(1)))) return(NULL)
  defensive <- defensive_rows(n_draws, engine$defensive_weight)
  share <- mean(defensive)
  log_weights <- numeric(n_draws)
  gradients <- matrix(0, n_draws, length(free))
  set.seed(seed)
  for (first in seq.int(1L, n_draws, by = block_rows)) {
    rows <- first:min(n_draws, first + block_rows - 1L)
    block <- matrix(stats::rnorm(length(rows) * dimension),
                    nrow = length(rows))
    weights_at <- function(variant) {
      engine$evaluate(shifted[[variant]], block, log_weights_only = TRUE,
                      built_proposal = proposals[[variant]],
                      defensive = defensive[rows], defensive_share = share)
    }
    log_weights[rows] <- weights_at(1L)
    for (parameter in seq_along(free)) {
      gradients[rows, parameter] <- (
        weights_at(2L * parameter) - weights_at(2L * parameter + 1L)
      ) / (2 * step[parameter])
    }
  }
  # Draws outside the support of the target carry no weight.
  valid <- is.finite(log_weights)
  if (!any(valid)) return(NULL)
  gradients[!valid, ] <- 0
  if (any(!is.finite(gradients))) return(NULL)
  diagnostics <- summarise_importance_weights(log_weights)
  relative <- exp(log_weights - max(log_weights[valid]))
  relative[!valid] <- 0
  weights <- relative / sum(relative)
  score <- colSums(gradients * weights)
  centred <- sweep(gradients, 2L, score) * weights
  if (!is.null(penalty)) {
    score <- score + numDeriv::grad(penalty, parameters)[free]
  }
  c(list(score = score, covariance = crossprod(centred)), diagnostics,
    list(seconds = proc.time()[["elapsed"]] - started))
}

##' Radius of the validation error in the metric of the observed information
##'
##' For e ~ N(0, C), ||H^{-1/2} e||^2 is a sum of independent chi-square
##' variables with one degree of freedom, weighted by the eigenvalues of
##' H^{-1/2} C H^{-1/2}, and the radius is the square root of its
##' `probability` point. The quantile is computed from fixed standard normal
##' draws made with a private seed, so that the radius is reproducible and the
##' user's random number stream is left unchanged. Numerical inversion of the
##' characteristic function (Imhof, 1961) is unreliable here, because the
##' integrand decays too slowly when one weight dominates.
##' @param information_root Upper Cholesky factor R of H = R'R.
##' @noRd
compute_validation_radius <- function(information_root, covariance,
                              probability = 0.95, n_draws = 100000L) {
  lower_root <- t(information_root)
  standardised <- t(forwardsolve(lower_root,
                                 t(forwardsolve(lower_root, covariance))))
  weights <- pmax(eigen((standardised + t(standardised)) / 2, symmetric = TRUE,
                        only.values = TRUE)$values, 0)
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  set.seed(424242L, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  squares <- matrix(stats::rnorm(n_draws * length(weights)), n_draws)^2
  sqrt(stats::quantile(as.numeric(squares %*% weights), probability,
                       names = FALSE))
}

##' Monte Carlo certificate of an estimate
##'
##' With H the observed information over the free parameters, S the score on
##' an independent validation bank at the estimate and C its Monte Carlo
##' covariance, the Newton correction is d = H^-1 S and the certificate is
##' U = ||d||_H + q, where q is the 95% radius of the validation error in the
##' same metric. Under the conditions of Fronterre (2026, Proposition 4), U
##' is an approximate 95% bound, in statistical standard errors, on the
##' distance from the estimate to the maximum of the exact likelihood near
##' it. It applies at any estimate, not only at a maximiser, and says nothing
##' about other maxima.
##' @param information_root Upper Cholesky factor R of H = R'R.
##' @noRd
certify_estimate <- function(information_root, score, covariance) {
  standardised_score <- forwardsolve(t(information_root), score)
  distance <- sqrt(sum(standardised_score^2))
  radius <- compute_validation_radius(information_root, covariance)
  list(correction = as.numeric(backsolve(information_root,
                                         standardised_score)),
       distance = distance, radius = radius, bound = distance + radius)
}

##' Newton steps with a known curvature
##'
##' Minimises `negative`, the negative estimated log-likelihood on the current
##' bank with the moving proposal, by Newton steps that use a fixed
##' positive-definite curvature (an observed information) through its upper
##' Cholesky factor `curvature_root`. Gradients are central differences of the
##' objective itself, so the objective is not approximated. Each step is
##' halved until it lowers the objective, and kept within the bounds. The
##' iteration stops when the information-scaled gradient
##' sqrt(g' C^-1 g), the length of the next step in standard errors, is below
##' `tolerance`.
##' @return A list in the form of an `nlminb()` result, with `converged` FALSE
##' and the best point reached if a step cannot lower the objective or the
##' iterations run out, so that the caller can continue with `nlminb()`.
##' @noRd
run_newton_steps <- function(negative, start, curvature_root, lower, upper,
                           tolerance = 1e-3, maximum_iterations = 20L,
                           relative_step = 1e-5) {
  evaluations <- 0L
  evaluate <- function(parameters) {
    evaluations <<- evaluations + 1L
    negative(parameters)
  }
  gradient_at <- function(parameters) {
    step <- relative_step * pmax(1, abs(parameters))
    vapply(seq_along(parameters), function(parameter) {
      shift <- replace(numeric(length(parameters)), parameter, step[parameter])
      (evaluate(parameters + shift) - evaluate(parameters - shift)) /
        (2 * step[parameter])
    }, numeric(1))
  }
  current <- pmin(pmax(start, lower), upper)
  value <- evaluate(current)
  result <- function(converged, iterations) {
    list(par = current, objective = value,
         convergence = if (converged) 0L else 1L, converged = converged,
         iterations = iterations,
         evaluations = c("function" = evaluations, gradient = 0L),
         message = if (converged) {
           "information-scaled gradient below tolerance"
         } else {
           "Newton steps with a fixed curvature stopped"
         },
         method = "newton")
  }
  if (!is.finite(value) || value >= 1e99) return(result(FALSE, 0L))
  for (iteration in seq_len(maximum_iterations)) {
    gradient <- gradient_at(current)
    if (any(!is.finite(gradient))) return(result(FALSE, iteration))
    standardised <- forwardsolve(t(curvature_root), gradient)
    if (sqrt(sum(standardised^2)) < tolerance) {
      return(result(TRUE, iteration - 1L))
    }
    direction <- -as.numeric(backsolve(curvature_root, standardised))
    accepted <- FALSE
    for (halving in 0:10) {
      candidate <- pmin(pmax(current + direction / 2^halving, lower), upper)
      candidate_value <- evaluate(candidate)
      if (is.finite(candidate_value) && candidate_value < value) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) return(result(FALSE, iteration))
    current <- candidate
    value <- candidate_value
  }
  result(FALSE, maximum_iterations)
}

##' Maximise the estimated log-likelihood on one fitting bank
##'
##' With `curvature_root`, the upper Cholesky factor of a known information
##' matrix, Newton steps (`run_newton_steps()`) move the start close to the
##' maximum first. `nlminb()` then always finishes the optimisation, from the
##' Newton end point, so that convergence is judged by the same criterion as
##' without the curvature; on a flat or badly scaled surface an old curvature
##' can make a step look shorter than it is. At the first stage a fit that
##' ends with a variance component at its boundary is restarted once from an
##' interior value.
##' @noRd
optimise_on_bank <- function(objective, start, fit_bank, lower, upper,
                             first_stage, limits, curvature_root = NULL) {
  negative <- function(parameters) -objective(parameters, fit_bank)
  newton <- if (is.null(curvature_root)) NULL else
    run_newton_steps(negative, start, curvature_root, lower, upper)
  from <- start
  if (!is.null(newton) && is.finite(newton$objective) &&
      newton$objective < 1e99) {
    from <- newton$par
  }
  fit <- nlminb(from, negative, lower = lower, upper = upper,
                control = make_optimiser_control(abs(negative(from))))
  fit$method <- if (is.null(newton)) "nlminb" else "newton, then nlminb"
  if (!is.null(newton)) {
    fit$evaluations[["function"]] <- fit$evaluations[["function"]] +
      newton$evaluations[["function"]]
  }
  restart <- if (first_stage && length(limits$variance_index)) {
    restart_from_interior(negative, fit, limits$variance_index,
                          limits$boundary_log_variance, lower, upper)
  } else {
    list(fit = fit, restarted = FALSE)
  }
  fit <- restart$fit
  boundary <- adaptive_boundary_parameters(fit$par, limits$variance_index,
                                           limits$boundary_log_variance,
                                           limits$linked_index)
  list(fit = fit, boundary = boundary, restarted = restart$restarted,
       restart_improved = isTRUE(restart$restart_improved))
}

##' Hessian by central differences with one round
##'
##' The formulas of `numDeriv::genD()` without Richardson extrapolation, with
##' 2k + k(k - 1) + 1 evaluations of `objective` for k parameters, against
##' about twice as many with two rounds, and an error of order step^2.
##' @noRd
compute_central_hessian <- function(objective, point, relative_step = 1e-4) {
  count <- length(point)
  step <- relative_step * pmax(1, abs(point))
  centre <- objective(point)
  hessian <- matrix(0, count, count)
  for (row in seq_len(count)) {
    shift <- replace(numeric(count), row, step[row])
    hessian[row, row] <- (objective(point + shift) - 2 * centre +
                            objective(point - shift)) / step[row]^2
  }
  for (row in seq_len(count)[-1L]) for (column in seq_len(row - 1L)) {
    pair <- c(row, column)
    shift <- replace(numeric(count), pair, step[pair])
    hessian[row, column] <- hessian[column, row] <- (
      objective(point + shift) - 2 * centre + objective(point - shift) -
        hessian[row, row] * step[row]^2 -
        hessian[column, column] * step[column]^2
    ) / (2 * step[row] * step[column])
  }
  hessian
}

##' Observed information on the fitting bank over the free parameters
##'
##' If it is singular, variance components that are small but not exactly
##' zero are the usual cause, so they are held as well before the model is
##' judged weakly identified. `root` is NULL when the information is singular
##' over the free parameters.
##' @noRd
compute_observed_information <- function(objective, estimate, fit_bank,
                                         boundary, limits) {
  information <- tryCatch(
    # One round of central differences suffices for this smooth,
    # deterministic objective: on binomial, Poisson, nugget and cluster
    # models it agreed with four Richardson rounds to about 1e-5 in the
    # metric of the information, closer than two rounds, at half their cost.
    compute_central_hessian(function(parameters) {
      -objective(parameters, fit_bank)
    }, estimate),
    error = function(error) NULL
  )
  if (!is.null(information) && any(!is.finite(information))) {
    information <- NULL
  }
  factor_free <- function(free) {
    if (is.null(information) || !length(free)) return(NULL)
    tryCatch(
      factor_covariance(information[free, free, drop = FALSE],
                        "adaptive MCML observed information",
                        allow_jitter = FALSE),
      error = function(error) NULL
    )
  }
  free <- setdiff(seq_along(estimate), boundary)
  root <- factor_free(free)
  variance_index <- limits$variance_index
  near_zero <- setdiff(
    variance_index[estimate[variance_index] <= limits$near_zero_log_variance],
    boundary
  )
  if (is.null(root) && !is.null(information) && length(near_zero)) {
    boundary <- sort(unique(c(
      boundary, near_zero,
      unlist(limits$linked_index[as.character(near_zero)], use.names = FALSE)
    )))
    free <- setdiff(seq_along(estimate), boundary)
    root <- factor_free(free)
  }
  list(information = information, root = root, free = free,
       boundary = boundary)
}

##' One stage of certified Monte Carlo maximum likelihood
##'
##' Maximises the estimated log-likelihood on a fitting bank of `budget` draws,
##' corrects the estimate by one Newton step with the score on an independent
##' validation bank, and certifies the corrected estimate with a second, fresh
##' validation bank drawn at that point. `seeds` gives the streams of the
##' three banks and of the independent bank of the check for other maxima. At
##' the first stage that check (`check_dispersed_maxima()`) may replace the
##' maximum. `curvature_root`, if given, seeds the optimisation with Newton
##' steps (`optimise_on_bank()`).
##' @noRd
run_certified_stage <- function(engine, objective, penalty, start, budget,
                            validation_budget, seeds, lower, upper,
                            first_stage, limits, curvature_root = NULL) {
  stage_started <- proc.time()[["elapsed"]]
  # Every global proposal of the stage starts from the factors at its start,
  # so each is a deterministic function of the parameters.
  if (identical(engine$proposal_type(), "global")) engine$set_anchor(start)
  set.seed(seeds[[1L]])
  fit_bank <- matrix(stats::rnorm(budget * engine$latent_dimension),
                     nrow = budget)
  optimised <- optimise_on_bank(objective, start, fit_bank, lower, upper,
                                first_stage, limits, curvature_root)
  check <- NULL
  separate <- ambiguous <- FALSE
  if (first_stage && !is.null(limits$dispersed)) {
    checked <- check_dispersed_maxima(engine, objective, penalty, optimised,
                                      fit_bank, lower, upper, limits,
                                      validation_budget, seeds[[4L]])
    optimised <- checked$optimised
    check <- checked$check
    separate <- checked$separate
    ambiguous <- checked$ambiguous
  }
  fit <- optimised$fit
  record <- list(
    budget = budget, validation_budget = validation_budget,
    proposal = engine$proposal_type(), fitted_estimate = fit$par,
    estimate = fit$par, convergence = fit$convergence,
    optimiser = fit$method,
    objective = fit$objective, iterations = fit$iterations,
    evaluations = fit$evaluations, message = fit$message,
    boundary_parameters = optimised$boundary,
    restarted = optimised$restarted,
    restart_improved = optimised$restart_improved,
    maxima_check = check, separate_maximum = separate,
    ambiguous_maxima = ambiguous,
    fit_diagnostics = summarise_importance_weights(
      engine$evaluate(fit$par, fit_bank, log_weights_only = TRUE)
    ),
    outcome = "no valid fit", certified = FALSE, certificate = NA_real_,
    timing = list(fit = proc.time()[["elapsed"]] - stage_started)
  )
  finish <- function(record, outcome, reason = NULL) {
    record$outcome <- outcome
    record$reason <- reason
    record$elapsed_seconds <- proc.time()[["elapsed"]] - stage_started
    record
  }
  if (!is.finite(fit$objective) || fit$objective >= 1e99) {
    return(finish(record, "no valid fit",
                  "the estimated log-likelihood is not finite at the optimum"))
  }
  started <- proc.time()[["elapsed"]]
  information <- compute_observed_information(objective, fit$par, fit_bank,
                                      optimised$boundary, limits)
  record$timing$information <- proc.time()[["elapsed"]] - started
  record$boundary_parameters <- information$boundary
  record$information <- information$information
  if (is.null(information$information)) {
    return(finish(record, "no valid fit",
                  "the observed information could not be computed"))
  }
  if (is.null(information$root)) {
    return(finish(record, "weak identification",
                  "singular observed information"))
  }
  record$free <- information$free
  record$information_root <- information$root
  certified <- correct_and_certify(engine, penalty, fit$par, information,
                                   validation_budget, seeds, lower, upper,
                                   limits$tolerance)
  if (is.null(certified)) {
    return(finish(record, "no valid fit",
                  "the validation score could not be computed"))
  }
  record[names(certified)] <- certified
  record$timing$validation <- certified$first_validation$seconds +
    certified$second_validation$seconds
  finish(record, if (record$certified) "certified" else "not certified")
}

##' Correct an estimate with one validation bank and certify it with another
##'
##' The estimate is corrected by one Newton step, d = H^-1 S, with the score
##' S on the first validation bank, and the corrected estimate is certified on
##' a second, fresh bank drawn at that point (`certify_estimate()`). The
##' second bank differentiates every parameter, including held ones, so that
##' the reported gradient is complete. Returns NULL if a validation score
##' cannot be computed.
##' @noRd
correct_and_certify <- function(engine, penalty, estimate, information,
                                validation_budget, seeds, lower, upper,
                                tolerance) {
  free <- information$free
  first <- score_validation_bank(engine, estimate, free, validation_budget,
                                 seeds[[2L]], penalty)
  if (is.null(first)) return(NULL)
  correction <- certify_estimate(information$root, first$score,
                                 first$covariance)
  corrected <- estimate
  corrected[free] <- pmin(pmax(corrected[free] + correction$correction,
                               lower[free]), upper[free])
  second <- score_validation_bank(engine, corrected, seq_along(corrected),
                                  validation_budget, seeds[[3L]], penalty)
  if (is.null(second)) return(NULL)
  certificate <- certify_estimate(information$root, second$score[free],
                                  second$covariance[free, free, drop = FALSE])
  list(estimate = corrected, correction_size = correction$distance,
       first_validation = first, second_validation = second,
       certificate = certificate$bound,
       certificate_distance = certificate$distance,
       certificate_radius = certificate$radius,
       certified = certificate$bound <= tolerance)
}

##' Check for other maxima on the fitting bank
##'
##' If the fitted spatial scale is below the trigger of `limits$dispersed`
##' (see `make_dispersed_starts()`) and the spatial variance is not at zero, the
##' estimated log-likelihood is also maximised from the dispersed starts on
##' the same bank, with the other parameters at their fitted values, and
##' `choose_maximum()` decides which maximum to keep. `n_draws` and `seed`
##' define the independent bank on which separate maxima are compared.
##' @return The kept `optimised` result, a `check` table (NULL if the check
##' was not triggered), whether the kept maximum is `separate` from the fitted
##' one, and whether the choice is `ambiguous`.
##' @noRd
check_dispersed_maxima <- function(engine, objective, penalty, optimised,
                                   fit_bank, lower, upper, limits, n_draws,
                                   seed) {
  dispersed <- limits$dispersed
  estimate <- optimised$fit$par
  if (dispersed$sigma2_index %in% optimised$boundary ||
      !(estimate[dispersed$phi_index] < dispersed$trigger)) {
    return(list(optimised = optimised, check = NULL, separate = FALSE,
                ambiguous = FALSE))
  }
  if (!is.null(dispersed$announce)) {
    dispersed$announce(estimate[[dispersed$phi_index]])
  }
  refits <- lapply(seq_along(dispersed$labels), function(dispersed_start) {
    start <- estimate
    start[dispersed$sigma2_index] <- dispersed$log_sigma2[dispersed_start]
    start[dispersed$phi_index] <- dispersed$log_phi[dispersed_start]
    optimise_on_bank(objective, pmin(pmax(start, lower), upper), fit_bank,
                     lower, upper, TRUE, limits)
  })
  candidates <- c(list(optimised), refits)
  information <- tryCatch(
    compute_central_hessian(function(parameters) {
      -objective(parameters, fit_bank)
    }, estimate),
    error = function(error) NULL
  )
  choice <- choose_maximum(
    engine, penalty, lapply(candidates, function(candidate) {
      candidate$fit$par
    }),
    vapply(candidates, function(candidate) -candidate$fit$objective,
           numeric(1)),
    information, n_draws, seed, limits$tolerance
  )
  choice$table$start <- c("fitted", dispersed$labels)
  choice$table$log_phi <- vapply(candidates, function(candidate) {
    candidate$fit$par[[dispersed$phi_index]]
  }, numeric(1))
  list(optimised = candidates[[choice$chosen]], check = choice$table,
       separate = choice$separate, ambiguous = choice$ambiguous)
}

##' Choose among candidate maxima of the estimated log-likelihood
##'
##' The first candidate is the fitted maximum and `values` are the
##' log-likelihoods on the fitting bank. Candidates within `tolerance`
##' standard errors of each other, in the metric of the observed `information`
##' at the fitted maximum (directions of zero or negative curvature carry no
##' weight), are the same maximum at the resolution of the certificate, and
##' the best of them on the fitting bank represents it. Separate maxima are
##' compared on an independent bank (`adjudicate_maxima()`), because the
##' fitting bank favours the points optimised on it, and the highest there is
##' kept. The choice is ambiguous if another separate maximum is within 2.58
##' standard errors of it on that bank; a small difference in likelihood
##' alone does not make two separated points the same maximum.
##' @noRd
choose_maximum <- function(engine, penalty, parameters, values, information,
                           n_draws, seed, tolerance) {
  metric_root <- NULL
  if (!is.null(information) && all(is.finite(information))) {
    decomposition <- eigen((information + t(information)) / 2,
                           symmetric = TRUE)
    metric_root <- sqrt(pmax(decomposition$values, 0)) *
      t(decomposition$vectors)
  }
  distance <- function(first, second) {
    if (is.null(metric_root)) {
      return(if (max(abs(first - second)) > 1e-6) Inf else 0)
    }
    sqrt(sum((metric_root %*% (first - second))^2))
  }
  representatives <- integer()
  maximum <- integer(length(values))
  for (candidate in order(values, decreasing = TRUE)) {
    near <- vapply(representatives, function(representative) {
      distance(parameters[[candidate]], parameters[[representative]]) <=
        tolerance
    }, logical(1))
    if (any(near)) {
      maximum[candidate] <- maximum[representatives[which(near)[1L]]]
    } else {
      representatives <- c(representatives, candidate)
      maximum[candidate] <- length(representatives)
    }
  }
  table <- data.frame(
    log_likelihood = values,
    distance_from_fitted = vapply(parameters, distance, numeric(1),
                                  second = parameters[[1L]]),
    maximum = maximum, independent_log_likelihood = NA_real_,
    independent_se = NA_real_, chosen = FALSE
  )
  chosen <- representatives[1L]
  ambiguous <- FALSE
  if (length(representatives) > 1L) {
    judged <- adjudicate_maxima(engine, penalty, parameters[representatives],
                                n_draws, seed)
    table$independent_log_likelihood[representatives] <- judged$values
    table$independent_se[representatives] <- judged$se
    chosen <- representatives[judged$top]
    others <- setdiff(seq_along(representatives), judged$top)
    gaps <- judged$values[judged$top] - judged$values[others]
    ambiguous <- isTRUE(any(gaps <= 2.58 * judged$se[others], na.rm = TRUE))
  }
  table$chosen[chosen] <- TRUE
  list(chosen = chosen, table = table, ambiguous = ambiguous,
       separate = maximum[chosen] != maximum[1L])
}

##' Log-likelihoods of candidate maxima on an independent bank
##'
##' Each candidate's estimated log-likelihood (with the penalty, if any) on a
##' fresh bank of `n_draws` draws from `seed`, and the standard error of its
##' difference from the highest, paired over the common draws by the delta
##' method.
##' @noRd
adjudicate_maxima <- function(engine, penalty, parameters, n_draws, seed) {
  set.seed(seed)
  bank <- matrix(stats::rnorm(n_draws * engine$latent_dimension),
                 nrow = n_draws)
  log_weights <- lapply(parameters, function(values) {
    engine$evaluate(values, bank, log_weights_only = TRUE)
  })
  # A candidate whose weights are all zero has log-likelihood -Inf and an
  # infinite standard error.
  scaled <- lapply(log_weights, function(weights) {
    if (!any(is.finite(weights))) return(rep(NA_real_, length(weights)))
    relative <- exp(weights - max(weights))
    relative / mean(relative)
  })
  values <- vapply(seq_along(parameters), function(candidate) {
    weights <- log_weights[[candidate]]
    if (!any(is.finite(weights))) return(-Inf)
    largest <- max(weights)
    largest + log(mean(exp(weights - largest))) +
      if (is.null(penalty)) 0 else penalty(parameters[[candidate]])
  }, numeric(1))
  top <- which.max(values)
  se <- vapply(seq_along(parameters), function(candidate) {
    difference <- scaled[[top]] - scaled[[candidate]]
    if (any(!is.finite(difference))) return(Inf)
    stats::sd(difference) / sqrt(n_draws)
  }, numeric(1))
  list(values = values, se = se, top = top)
}

##' Upper Cholesky factor of a curvature for Newton steps
##'
##' A curvature is used only if it is finite, positive definite and not
##' nearly singular (reciprocal condition number at least 1e-8). Returns
##' NULL otherwise.
##' @noRd
factor_curvature <- function(curvature) {
  if (is.null(curvature) || any(!is.finite(curvature))) return(NULL)
  curvature <- (curvature + t(curvature)) / 2
  if (rcond(curvature) < 1e-8) return(NULL)
  tryCatch(chol(curvature), error = function(error) NULL)
}

##' Bank sizes for the next stage
##'
##' The certificate decreases approximately with the square root of the
##' number of draws, so the next fitting bank is the one predicted to bring it
##' to the tolerance, with a safety factor of 1.25 against the variability of
##' the current certificate. The validation banks grow in the same proportion.
##' The arithmetic is in double precision and the cap is applied before
##' conversion to integers, and a certificate that is not finite moves
##' straight to the largest budget.
##' @noRd
plan_next_budgets <- function(budget, validation_budget, certificate,
                              tolerance, maximum_budget) {
  wanted <- if (is.finite(certificate)) {
    budget * (certificate / tolerance)^2 * 1.25
  } else {
    Inf
  }
  next_budget <- min(max(500 * ceiling(wanted / 500), budget + 500),
                     maximum_budget)
  next_validation <- 500 * ceiling(validation_budget * next_budget / budget /
                                     500)
  list(budget = as.integer(next_budget),
       validation_budget = as.integer(next_validation))
}

##' Fit and certify Monte Carlo maximum likelihood estimates
##'
##' Runs certified stages with growing banks until the certificate is within
##' `tolerance`, the model is found to be weakly identified, no valid fit can
##' be obtained, or the largest budget is reached. If the global proposal
##' fails anywhere within a stage, the stage is repeated with the Laplace
##' proposal, which is then kept. Stage s uses the random number streams
##' 3s - 2 (fitting bank), 3s - 1 and 3s (validation banks) after `seed`, and
##' stream 1000 + s for the independent bank of the check for other maxima.
##' `dispersed`, from `make_dispersed_starts()`, turns on the check for other
##' maxima at the first stage. Each later stage starts its optimisation with
##' Newton steps on the previous stage's information when every parameter was
##' free.
##' @noRd
fit_adaptive_mcml <- function(engine, start, objective = engine$evaluate,
                              penalty = NULL, dispersed = NULL,
                              initial_budget = 1000L,
                              validation_budget = 5000L,
                              maximum_budget = 16000L, maximum_stages = 4L,
                              tolerance = 0.1, seed = 1L,
                              lower = NULL, upper = NULL,
                              variance_index = integer(),
                              linked_index = list(),
                              boundary_log_variance = -9,
                              near_zero_log_variance = -4) {
  if (is.null(lower)) lower <- rep(-10, length(start))
  if (is.null(upper)) upper <- rep(10, length(start))
  validation_budget <- as.integer(validation_budget)
  limits <- list(
    variance_index = variance_index, linked_index = linked_index,
    boundary_log_variance = boundary_log_variance,
    near_zero_log_variance = near_zero_log_variance, tolerance = tolerance,
    dispersed = dispersed
  )
  stream <- function(offset) {
    as.integer(((seed - 1 + offset) %% .Machine$integer.max) + 1)
  }
  total_started <- proc.time()[["elapsed"]]
  history <- list()
  current <- start
  budget <- as.integer(initial_budget)
  curvature_root <- NULL
  for (stage in seq_len(maximum_stages)) {
    run_stage <- function() {
      engine$reset_proposal_failures()
      run_certified_stage(engine, objective, penalty, current, budget,
                      validation_budget,
                      c(stream(3L * stage - 2:0), stream(1000L + stage)),
                      lower, upper, stage == 1L, limits, curvature_root)
    }
    record <- run_stage()
    record$proposal_failures <- engine$proposal_failures()
    if (record$proposal_failures > 0L &&
        identical(engine$proposal_type(), "global")) {
      engine$use_proposal("laplace")
      failures <- record$proposal_failures
      record <- run_stage()
      record$proposal_failures <- failures
      record$proposal_fallback <- TRUE
    }
    record$stage <- stage
    history[[stage]] <- record
    current <- record$estimate
    curvature_root <- if (length(record$free) == length(current)) {
      factor_curvature(record$information)
    }
    if (!identical(record$outcome, "not certified") ||
        budget >= maximum_budget || stage == maximum_stages) break
    budgets <- plan_next_budgets(budget, validation_budget,
                                 record$certificate, tolerance,
                                 maximum_budget)
    history[[stage]]$next_budget <- budgets$budget
    budget <- budgets$budget
    validation_budget <- budgets$validation_budget
  }
  last <- history[[length(history)]]
  covariance <- NULL
  if (!is.null(last$information_root)) {
    covariance <- matrix(NA_real_, length(current), length(current))
    covariance[last$free, last$free] <- chol2inv(last$information_root)
  }
  score <- if (is.null(last$second_validation)) {
    rep(NA_real_, length(current))
  } else {
    last$second_validation$score
  }
  list(
    estimate = current, outcome = last$outcome,
    maxima_check = history[[1L]]$maxima_check,
    better_maximum = isTRUE(history[[1L]]$separate_maximum),
    ambiguous_maxima = isTRUE(history[[1L]]$ambiguous_maxima),
    certified = isTRUE(last$certified), certificate = last$certificate,
    reason = last$reason, budget = last$budget,
    validation_budget = last$validation_budget,
    proposal = engine$proposal_type(),
    boundary_parameters = last$boundary_parameters,
    score = score,
    log_likelihood = if (is.null(last$second_validation)) {
      last$fit_diagnostics$log_likelihood
    } else {
      last$second_validation$log_likelihood
    },
    covariance = covariance, information = last$information,
    history = history,
    elapsed_seconds = proc.time()[["elapsed"]] - total_started
  )
}
