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
##' about twice as many with two rounds, and an error of order step^2. With
##' `with_gradient`, it returns a list with the Hessian and the gradient by
##' central differences from the same evaluations.
##' @noRd
compute_central_hessian <- function(objective, point, relative_step = 1e-4,
                                    with_gradient = FALSE) {
  count <- length(point)
  step <- relative_step * pmax(1, abs(point))
  centre <- objective(point)
  hessian <- matrix(0, count, count)
  gradient <- numeric(count)
  for (row in seq_len(count)) {
    shift <- replace(numeric(count), row, step[row])
    above <- objective(point + shift)
    below <- objective(point - shift)
    hessian[row, row] <- (above - 2 * centre + below) / step[row]^2
    gradient[row] <- (above - below) / (2 * step[row])
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
  if (with_gradient) list(hessian = hessian, gradient = gradient) else hessian
}

##' Observed information on the fitting bank over the free parameters
##'
##' If it is singular, variance components that are small but not exactly
##' zero are the usual cause, so they are held as well before the model is
##' judged weakly identified. `root` is NULL when the information is singular
##' over the free parameters. `gradient` is the fitting bank's score at
##' `estimate`, from the same evaluations.
##' @noRd
compute_observed_information <- function(objective, estimate, fit_bank,
                                         boundary, limits) {
  differences <- tryCatch(
    # One round of central differences suffices for this smooth,
    # deterministic objective: on binomial, Poisson, nugget and cluster
    # models it agreed with four Richardson rounds to about 1e-5 in the
    # metric of the information, closer than two rounds, at half their cost.
    compute_central_hessian(function(parameters) {
      -objective(parameters, fit_bank)
    }, estimate, with_gradient = TRUE),
    error = function(error) NULL
  )
  information <- differences$hessian
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
       boundary = boundary, gradient = -differences$gradient)
}

##' One stage of certified Monte Carlo maximum likelihood
##'
##' Maximises the estimated log-likelihood on a fitting bank of `budget` draws,
##' corrects the estimate by one Newton step with the score on an independent
##' validation bank, and certifies the corrected estimate with a second, fresh
##' validation bank drawn at that point. `seeds` gives the streams of the
##' three banks and of the independent bank of the check for other maxima. At
##' the first stage that check (`check_other_maxima()`) may replace the
##' maximum, and if the fitted spatial correlation is numerically negligible
##' the stage ends with the outcome "numerically independent", without
##' information or certificate. The fit is then one of independent location
##' effects and the spatial scale is unresolved. `curvature_root`, if
##' given, seeds the optimisation with Newton steps (`optimise_on_bank()`).
##' @noRd
run_certified_stage <- function(engine, objective, penalty, start, budget,
                            validation_budget, seeds, lower, upper,
                            first_stage, limits, curvature_root = NULL) {
  stage_started <- proc.time()[["elapsed"]]
  # Every global proposal of the stage starts from the factors at its start,
  # so each is a deterministic function of the parameters.
  if (identical(engine$proposal_type(), "global")) {
    engine$set_anchor(start)
    # A global proposal that fails at the start already condemns the stage,
    # which fit_adaptive_mcml() repeats with the Laplace proposal on the same
    # streams, so the rest of it is not computed.
    if (engine$proposal_failures() > 0L) {
      return(list(
        budget = budget, validation_budget = validation_budget,
        proposal = "global", fitted_estimate = start, estimate = start,
        outcome = "no valid fit", certified = FALSE, certificate = NA_real_,
        reason = "the global proposal failed at the start of the stage",
        elapsed_seconds = proc.time()[["elapsed"]] - stage_started
      ))
    }
  }
  set.seed(seeds[[1L]])
  fit_bank <- matrix(stats::rnorm(budget * engine$latent_dimension),
                     nrow = budget)
  optimised <- optimise_on_bank(objective, start, fit_bank, lower, upper,
                                first_stage, limits, curvature_root)
  check <- screen <- NULL
  separate <- ambiguous <- independent <- FALSE
  if (first_stage && !is.null(limits$dispersed) &&
      is_valid_fit(optimised$fit)) {
    checked <- check_other_maxima(engine, objective, penalty, optimised,
                                  fit_bank, lower, upper, limits,
                                  validation_budget, seeds[[4L]],
                                  if (length(seeds) >= 5L) seeds[[5L]])
    optimised <- checked$optimised
    check <- checked$check
    screen <- checked$screen
    separate <- checked$separate
    ambiguous <- checked$ambiguous
    independent <- checked$numerically_independent
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
    maxima_check = check, maxima_screen = screen,
    separate_maximum = separate, ambiguous_maxima = ambiguous,
    fit_diagnostics = if (all(is.finite(fit$par))) {
      summarise_importance_weights(
        engine$evaluate(fit$par, fit_bank, log_weights_only = TRUE)
      )
    },
    outcome = "no valid fit", certified = FALSE, certificate = NA_real_,
    timing = list(fit = proc.time()[["elapsed"]] - stage_started)
  )
  finish <- function(record, outcome, reason = NULL) {
    record$outcome <- outcome
    record$reason <- reason
    record$elapsed_seconds <- proc.time()[["elapsed"]] - stage_started
    record
  }
  if (!all(is.finite(fit$par))) {
    return(finish(record, "no valid fit",
                  "the optimiser returned non-finite parameter values"))
  }
  if (!is_valid_fit(fit)) {
    return(finish(record, "no valid fit",
                  "the estimated log-likelihood is not finite at the optimum"))
  }
  if (independent) {
    return(finish(record, "numerically independent", paste(
      "the fitted spatial correlation is numerically negligible, so the fit",
      "is one of independent location effects and the spatial scale is",
      "unresolved"
    )))
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
  # Proposal failures so far depend only on the fitting draws. Any later one
  # depends on the correction bank (fit_adaptive_mcml()).
  record$failures_before_validation <- engine$proposal_failures()
  certified <- correct_and_certify(
    engine, penalty, fit$par, information, validation_budget, seeds, lower,
    upper, limits$tolerance,
    bank_objective = function(parameters) objective(parameters, fit_bank),
    correction = limits$correction %||% "single"
  )
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
##'
##' One step leaves a residual where the log-likelihood is far from quadratic
##' between the fitted estimate and the maximum, as along a log variance near
##' zero, and the certificate, measured in the curvature at the fitted
##' estimate, then sees only part of it. `bank_objective`, the log-likelihood
##' on the fitting bank as a function of the parameters, allows the
##' `correction` methods other than "single":
##' - "recurve" certifies the corrected estimate in the fitting bank's
##'   curvature at that estimate;
##' - "repeat" predicts the residual from the fitting bank
##'   (`predict_correction_residual()`) and, if it exceeds `residual_gate`
##'   radii of the correction bank, continues the correction on that bank
##'   (`repeat_correction()`), and certifies in the last curvature, unless
##'   those steps stop short of the bank's maximum by more than 2 - sqrt(2)
##'   radii, the most residual that Lemma 4.2 of Fronterre (2026) tolerates;
##' - "refuse" does not certify an estimate whose predicted residual exceeds
##'   2 - sqrt(2) radii.
##' A stage that is not certified for its residual goes on to the next stage
##' (`fit_adaptive_mcml()`), which starts from the corrected estimate.
##' This choice of methods is a prototype for comparing them.
##' @noRd
correct_and_certify <- function(engine, penalty, estimate, information,
                                validation_budget, seeds, lower, upper,
                                tolerance, bank_objective = NULL,
                                correction = "single", residual_gate = 0.25) {
  free <- information$free
  root <- information$root
  first <- score_validation_bank(engine, estimate, free, validation_budget,
                                 seeds[[2L]], penalty)
  if (is.null(first)) return(NULL)
  correction_step <- certify_estimate(root, first$score, first$covariance)
  corrected <- estimate
  corrected[free] <- pmin(pmax(corrected[free] + correction_step$correction,
                               lower[free]), upper[free])
  if (is.null(bank_objective)) correction <- "single"
  record <- list(correction_method = correction, correction_steps = 1L)
  if (identical(correction, "recurve")) {
    curvature <- compute_bank_curvature(bank_objective, corrected, free)
    if (!is.null(curvature$root)) {
      root <- curvature$root
      record$information <- curvature$information
      record$information_root <- root
    }
  } else if (correction %in% c("repeat", "refuse")) {
    record$predicted_residual <- predict_correction_residual(
      bank_objective, estimate, corrected, free, root,
      corrected[free] - estimate[free], correction_step$radius,
      information$gradient[free]
    )
    if (identical(correction, "repeat") &&
        record$predicted_residual > residual_gate) {
      repeated <- repeat_correction(engine, penalty, corrected, root, free,
                                    validation_budget, seeds[[2L]],
                                    bank_objective, lower, upper,
                                    correction_step$radius)
      corrected <- repeated$estimate
      root <- repeated$root
      record$correction_steps <- 1L + length(repeated$lengths)
      record$correction_lengths <- repeated$lengths
      record$correction_converged <- repeated$converged
      record$correction_remaining <- repeated$remaining
      if (!is.null(repeated$information)) {
        record$information <- repeated$information
        record$information_root <- root
      }
    }
  }
  second <- score_validation_bank(engine, corrected, seq_along(corrected),
                                  validation_budget, seeds[[3L]], penalty)
  if (is.null(second)) return(NULL)
  certificate <- certify_estimate(root, second$score[free],
                                  second$covariance[free, free, drop = FALSE])
  refused <- switch(
    correction,
    refuse = record$predicted_residual > 2 - sqrt(2),
    "repeat" = isFALSE(record$correction_converged) &&
      record$correction_remaining > 2 - sqrt(2),
    FALSE
  )
  c(list(estimate = corrected, correction_size = correction_step$distance,
         first_validation = first, second_validation = second,
         certificate = certificate$bound,
         certificate_distance = certificate$distance,
         certificate_radius = certificate$radius,
         certified = certificate$bound <= tolerance && !refused,
         residual_refused = refused),
    record)
}

##' Gradient by central differences over the parameters indexed by `index`
##' @noRd
compute_central_gradient <- function(objective, point,
                                     index = seq_along(point),
                                     relative_step = 1e-4) {
  step <- relative_step * pmax(1, abs(point))
  vapply(index, function(parameter) {
    shift <- replace(numeric(length(point)), parameter, step[parameter])
    (objective(point + shift) - objective(point - shift)) /
      (2 * step[parameter])
  }, numeric(1))
}

##' Observed information on the fitting bank at a point
##'
##' `bank_objective` is the log-likelihood on the fitting bank. Returns the
##' information over all parameters and its upper Cholesky factor over
##' `free`, which is NULL if the information is not positive definite there.
##' @noRd
compute_bank_curvature <- function(bank_objective, point, free) {
  information <- tryCatch(
    compute_central_hessian(function(parameters) {
      -bank_objective(parameters)
    }, point),
    error = function(error) NULL
  )
  if (is.null(information) || any(!is.finite(information))) {
    return(list(information = NULL, root = NULL))
  }
  root <- tryCatch(
    factor_covariance(information[free, free, drop = FALSE],
                      "adaptive MCML observed information",
                      allow_jitter = FALSE),
    error = function(error) NULL
  )
  list(information = information, root = root)
}

##' Residual of the Newton correction predicted on the fitting bank
##'
##' If the log-likelihood were quadratic between the fitted estimate and the
##' corrected one, the fitting bank's score would change by -H d over the
##' correction d. What it changes by in addition is the score left at the
##' corrected estimate, and its length in the metric of H, in `radius` units
##' (the radius of the correction bank), predicts the step that a further
##' Newton correction would take. It costs 2k evaluations on the fitting bank
##' for k free parameters, and 2k more without `fitted_gradient`, the
##' fitting bank's score at the fitted estimate over the free parameters,
##' which the evaluations for H already give.
##' @param information_root Upper Cholesky factor R of H = R'R.
##' @noRd
predict_correction_residual <- function(bank_objective, fitted, corrected,
                                        free, information_root, step,
                                        radius, fitted_gradient = NULL) {
  if (is.null(fitted_gradient)) {
    fitted_gradient <- compute_central_gradient(bank_objective, fitted, free)
  }
  change <- compute_central_gradient(bank_objective, corrected, free) -
    fitted_gradient
  left <- change + drop(crossprod(information_root) %*% step)
  sqrt(sum(forwardsolve(t(information_root), left)^2)) / radius
}

##' Continue the Newton correction on the correction bank
##'
##' Newton steps from `corrected` on the score of the correction bank (the
##' bank of `n_draws` draws from `seed`), each in the fitting bank's curvature
##' at the current point, or the last positive-definite one, until a step is
##' shorter than `step_tolerance` radii of that bank, a step would be longer
##' than the one before, or `maximum_steps` steps have been taken. The
##' estimate then approaches the correction bank's own maximum, whose error is
##' the noise of that bank alone. Each step costs one score on the
##' correction bank and one Hessian on the fitting bank. Returns the
##' estimate, the curvature and its factor at the last point where it was
##' computed, the lengths of the steps taken in radii, whether the last was
##' below the tolerance, and the length of the last step computed, taken or
##' not, as a bound on the residual left.
##' @param root Upper Cholesky factor of the curvature used for the first
##' correction, kept if the fitting bank's is never positive definite.
##' @noRd
repeat_correction <- function(engine, penalty, corrected, root, free,
                              n_draws, seed, bank_objective, lower, upper,
                              radius, step_tolerance = 0.1,
                              maximum_steps = 4L) {
  lengths <- numeric()
  information <- NULL
  converged <- FALSE
  remaining <- Inf
  for (step in seq_len(maximum_steps)) {
    curvature <- compute_bank_curvature(bank_objective, corrected, free)
    rescored <- score_validation_bank(engine, corrected, free, n_draws, seed,
                                      penalty)
    if (is.null(rescored)) break
    step_root <- curvature$root %||% root
    standardised <- forwardsolve(t(step_root), rescored$score)
    length <- sqrt(sum(standardised^2)) / radius
    remaining <- length
    if (length(lengths) && length > lengths[length(lengths)]) break
    root <- step_root
    if (!is.null(curvature$root)) information <- curvature$information
    lengths <- c(lengths, length)
    corrected[free] <- pmin(pmax(corrected[free] +
                                   backsolve(root, standardised),
                                 lower[free]), upper[free])
    if (length < step_tolerance) {
      converged <- TRUE
      break
    }
  }
  list(estimate = corrected, root = root, information = information,
       lengths = lengths, converged = converged, remaining = remaining)
}

##' Whether an optimiser result can be used as a maximum
##'
##' The estimates and the objective must be finite, and the objective must
##' not be the penalty value that the estimated log-likelihood takes where it
##' cannot be computed.
##' @noRd
is_valid_fit <- function(fit) {
  all(is.finite(fit$par)) && is.finite(fit$objective) && fit$objective < 1e99
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
##' alone does not make two separated points the same maximum. An ambiguous
##' choice is made again on a bank of four times `n_draws` from
##' `enlarged_seed`, when given, and that bank decides.
##' @noRd
choose_maximum <- function(engine, penalty, parameters, values, information,
                           n_draws, seed, tolerance, enlarged_seed = NULL) {
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
    independent_se = NA_real_, independent_draws = NA_integer_,
    chosen = FALSE
  )
  chosen <- representatives[1L]
  ambiguous <- FALSE
  judge <- function(draws, stream) {
    judged <- adjudicate_maxima(engine, penalty, parameters[representatives],
                                draws, stream)
    others <- setdiff(seq_along(representatives), judged$top)
    gaps <- judged$values[judged$top] - judged$values[others]
    judged$ambiguous <- isTRUE(any(gaps <= 2.58 * judged$se[others],
                                   na.rm = TRUE))
    judged$draws <- as.integer(draws)
    judged
  }
  if (length(representatives) > 1L) {
    judged <- judge(n_draws, seed)
    # An ambiguous choice is judged again on a bank four times larger, which
    # resolves differences that the first bank's error hides.
    if (judged$ambiguous && !is.null(enlarged_seed)) {
      judged <- judge(4L * n_draws, enlarged_seed)
    }
    table$independent_log_likelihood[representatives] <- judged$values
    table$independent_se[representatives] <- judged$se
    table$independent_draws[representatives] <- judged$draws
    chosen <- representatives[judged$top]
    ambiguous <- judged$ambiguous
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
##' be obtained, or the largest budget is reached. Stage s uses the random
##' number streams 3s - 2 (fitting bank), 3s - 1 and 3s (validation banks)
##' after `seed`, stream 1000 + s for the independent bank of the check for
##' other maxima and stream 4000 + s for its enlarged bank.
##'
##' If the global proposal fails within a stage, the stage is repeated with
##' the Laplace proposal, which is then kept. A failure before the validation
##' banks (in the fit, the check for other maxima or the information) depends
##' only on the fitting draws, and the stage is repeated on the same streams. A
##' failure at the start of the stage ends it at once (`run_certified_stage()`),
##' before any bank is used.
##' A failure while correcting or certifying depends on the correction bank,
##' so the repeated stage draws every bank afresh, from streams
##' 2000 + 3s - 2 to 2000 + 3s and 3000 + s. The certificate of the repeated
##' stage then rests on a pair of validation banks drawn after everything that
##' led to it, as the conditional guarantee of the procedure requires.
##' `dispersed`, from `make_dispersed_starts()`, turns on the check for other
##' maxima at the first stage, with the screen of short scales and the exact
##' alarms when it holds their settings (`make_short_scale_screen()`). Each later stage starts its optimisation with
##' Newton steps on the previous stage's information when every parameter was
##' free. `correction` chooses how the fitted estimate is corrected and
##' certified (`correct_and_certify()`). The option `RiskMap.mcml_correction`
##' that sets its default exists only to compare the methods in this
##' prototype.
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
                              near_zero_log_variance = -4,
                              correction = getOption("RiskMap.mcml_correction",
                                                     "repeat")) {
  if (is.null(lower)) lower <- rep(-10, length(start))
  if (is.null(upper)) upper <- rep(10, length(start))
  validation_budget <- as.integer(validation_budget)
  correction <- match.arg(correction,
                          c("repeat", "single", "recurve", "refuse"))
  limits <- list(
    variance_index = variance_index, linked_index = linked_index,
    boundary_log_variance = boundary_log_variance,
    near_zero_log_variance = near_zero_log_variance, tolerance = tolerance,
    dispersed = dispersed, correction = correction
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
    run_stage <- function(offset) {
      engine$reset_proposal_failures()
      run_certified_stage(engine, objective, penalty, current, budget,
                      validation_budget,
                      c(stream(offset + 3L * stage - 2:0),
                        stream(1000L + offset + stage),
                        stream(4000L + offset + stage)),
                      lower, upper, stage == 1L, limits, curvature_root)
    }
    record <- run_stage(0L)
    record$proposal_failures <- engine$proposal_failures()
    if (record$proposal_failures > 0L &&
        identical(engine$proposal_type(), "global")) {
      engine$use_proposal("laplace")
      failures <- record$proposal_failures
      late <- failures > (record$failures_before_validation %||% failures)
      record <- run_stage(if (late) 2000L else 0L)
      record$proposal_failures <- failures
      record$proposal_fallback <- TRUE
      record$late_fallback <- late
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
    maxima_screen = history[[1L]]$maxima_screen,
    maxima_alarm = isTRUE(history[[1L]]$maxima_screen$alarm$persists),
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
