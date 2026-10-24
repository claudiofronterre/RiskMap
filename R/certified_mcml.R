##' Control settings for the optimiser in Monte Carlo maximum likelihood
##'
##' The common-random-number objective is smooth, so tight tolerances are
##' reached without difficulty.
##' @noRd
mcml_optimiser_control <- function() {
  list(eval.max = 600, iter.max = 300, rel.tol = 1e-10, x.tol = 1e-8)
}

##' Effective sample size and tail diagnostics of importance weights
##' @noRd
importance_diagnostics <- function(log_weights) {
  valid <- is.finite(log_weights)
  if (!any(valid)) {
    return(list(log_likelihood = -Inf, relative_ess = 0, maximum_weight = 1,
                pareto_k = Inf))
  }
  largest <- max(log_weights[valid])
  relative <- exp(log_weights - largest)
  relative[!valid] <- 0
  weights <- relative / sum(relative)
  psis_fit <- suppressWarnings(loo::psis(log_weights[valid]))
  list(
    log_likelihood = largest + log(mean(relative)),
    relative_ess = 1 / sum(weights^2) / length(weights),
    maximum_weight = max(weights),
    pareto_k = as.numeric(loo::pareto_k_values(psis_fit))
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
validation_score <- function(engine, parameters, free, n_draws, seed,
                             penalty = NULL, block_rows = NULL) {
  started <- proc.time()[["elapsed"]]
  dimension <- engine$latent_dimension
  if (is.null(block_rows)) block_rows <- max(250L, floor(2.5e6 / dimension))
  step <- .Machine$double.eps^(1 / 3) * pmax(1, abs(parameters[free]))
  shifted <- list(parameters)
  for (j in seq_along(free)) {
    upper <- lower <- parameters
    upper[free[j]] <- upper[free[j]] + step[j]
    lower[free[j]] <- lower[free[j]] - step[j]
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
    z <- matrix(stats::rnorm(length(rows) * dimension), nrow = length(rows))
    weights_at <- function(k) {
      engine$evaluate(shifted[[k]], z, log_weights_only = TRUE,
                      built_proposal = proposals[[k]],
                      defensive = defensive[rows], defensive_share = share)
    }
    log_weights[rows] <- weights_at(1L)
    for (j in seq_along(free)) {
      gradients[rows, j] <-
        (weights_at(2L * j) - weights_at(2L * j + 1L)) / (2 * step[j])
    }
  }
  # Draws outside the support of the target carry no weight.
  valid <- is.finite(log_weights)
  if (!any(valid)) return(NULL)
  gradients[!valid, ] <- 0
  if (any(!is.finite(gradients))) return(NULL)
  diagnostics <- importance_diagnostics(log_weights)
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
validation_radius <- function(information_root, covariance,
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
##' same metric. U bounds the distance from the estimate to the exact maximum
##' likelihood estimate, in statistical standard errors, with probability
##' about 0.95 (Fronterre, 2026, Proposition 4). It holds at any estimate, not
##' only at a maximiser.
##' @param information_root Upper Cholesky factor R of H = R'R.
##' @noRd
mcml_certificate <- function(information_root, score, covariance) {
  standardised_score <- forwardsolve(t(information_root), score)
  distance <- sqrt(sum(standardised_score^2))
  radius <- validation_radius(information_root, covariance)
  list(correction = as.numeric(backsolve(information_root,
                                         standardised_score)),
       distance = distance, radius = radius, bound = distance + radius)
}

##' Maximise the estimated log-likelihood on one fitting bank
##'
##' At the first stage a fit that ends with a variance component at its
##' boundary is restarted once from an interior value. A range shorter than
##' `minimum_log_range` is held there and the other parameters re-optimised,
##' unless the spatial variance is already zero.
##' @noRd
optimise_on_bank <- function(objective, start, fit_bank, lower, upper,
                             first_stage, limits) {
  negative <- function(parameters) -objective(parameters, fit_bank)
  fit <- nlminb(start, negative, lower = lower, upper = upper,
                control = mcml_optimiser_control())
  restart <- if (first_stage && length(limits$variance_index)) {
    restart_from_interior(negative, fit, limits$variance_index,
                          limits$boundary_log_variance, lower, upper)
  } else {
    list(fit = fit, restarted = FALSE)
  }
  fit <- restart$fit
  boundary_of <- function(estimate) {
    adaptive_boundary_parameters(estimate, limits$variance_index,
                                 limits$boundary_log_variance,
                                 limits$linked_index)
  }
  boundary <- boundary_of(fit$par)
  range_index <- limits$range_index
  # A range shorter than the spacing between locations cannot be
  # distinguished from unstructured variation. Unless the spatial variance is
  # already zero (when the range is meaningless anyway), the range is held at
  # that spacing and the remaining parameters are re-optimised.
  if (!is.null(range_index) && !(range_index %in% boundary) &&
      fit$par[range_index] < limits$minimum_log_range) {
    lower[range_index] <- limits$minimum_log_range
    held <- tryCatch(
      nlminb(pmax(fit$par, lower), negative, lower = lower, upper = upper,
             control = mcml_optimiser_control()),
      error = function(error) NULL
    )
    if (!is.null(held)) fit <- held
    boundary <- boundary_of(fit$par)
  }
  range_held <- !is.null(range_index) && !(range_index %in% boundary) &&
    fit$par[range_index] <= limits$minimum_log_range + 1e-6
  if (range_held) boundary <- sort(c(boundary, range_index))
  list(fit = fit, boundary = boundary, range_held = range_held, lower = lower,
       restarted = restart$restarted,
       restart_improved = isTRUE(restart$restart_improved))
}

##' Observed information on the fitting bank over the free parameters
##'
##' If it is singular, variance components that are small but not exactly
##' zero are the usual cause, so they are held as well before the model is
##' judged weakly identified. `root` is NULL when the information is singular
##' over the free parameters.
##' @noRd
observed_information <- function(objective, estimate, fit_bank, boundary,
                                 limits) {
  information <- tryCatch(
    # Two Richardson rounds suffice for this smooth, deterministic
    # objective: standard errors agree with the default four rounds to about
    # 1e-4 relative error at half the cost.
    numDeriv::hessian(function(parameters) -objective(parameters, fit_bank),
                      estimate, method.args = list(r = 2)),
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
##' three banks.
##' @noRd
certified_stage <- function(engine, objective, penalty, start, budget,
                            validation_budget, seeds, lower, upper,
                            first_stage, limits) {
  stage_started <- proc.time()[["elapsed"]]
  set.seed(seeds[[1L]])
  fit_bank <- matrix(stats::rnorm(budget * engine$latent_dimension),
                     nrow = budget)
  optimised <- optimise_on_bank(objective, start, fit_bank, lower, upper,
                                first_stage, limits)
  fit <- optimised$fit
  record <- list(
    budget = budget, validation_budget = validation_budget,
    proposal = engine$proposal_type(), fitted_estimate = fit$par,
    estimate = fit$par, convergence = fit$convergence,
    objective = fit$objective, iterations = fit$iterations,
    evaluations = fit$evaluations, message = fit$message,
    boundary_parameters = optimised$boundary,
    range_held = optimised$range_held, restarted = optimised$restarted,
    restart_improved = optimised$restart_improved, lower = optimised$lower,
    fit_diagnostics = importance_diagnostics(
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
  information <- observed_information(objective, fit$par, fit_bank,
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
  free <- information$free
  record$free <- free
  record$information_root <- information$root
  first <- validation_score(engine, fit$par, free, validation_budget,
                            seeds[[2L]], penalty)
  if (is.null(first)) {
    return(finish(record, "no valid fit",
                  "the validation score could not be computed"))
  }
  correction <- mcml_certificate(information$root, first$score,
                                 first$covariance)
  corrected <- fit$par
  corrected[free] <- pmin(pmax(corrected[free] + correction$correction,
                               lower[free]), upper[free])
  # The reported gradient covers every parameter, so the second bank
  # differentiates the held parameters as well.
  second <- validation_score(engine, corrected, seq_along(corrected),
                             validation_budget, seeds[[3L]], penalty)
  if (is.null(second)) {
    return(finish(record, "no valid fit",
                  "the validation score could not be computed"))
  }
  certificate <- mcml_certificate(information$root, second$score[free],
                                  second$covariance[free, free, drop = FALSE])
  record$estimate <- corrected
  record$correction_size <- correction$distance
  record$first_validation <- first
  record$second_validation <- second
  record$certificate <- certificate$bound
  record$certificate_distance <- certificate$distance
  record$certificate_radius <- certificate$radius
  record$certified <- certificate$bound <= limits$tolerance
  record$timing$validation <- first$seconds + second$seconds
  finish(record, if (record$certified) "certified" else "not certified")
}

##' Bank sizes for the next stage
##'
##' The certificate decreases approximately with the square root of the
##' number of draws, so the next fitting bank is the one predicted to bring it
##' to the tolerance, with a safety factor of 1.25 against the variability of
##' the current certificate. The validation banks grow in the same proportion.
##' @noRd
next_mcml_budgets <- function(budget, validation_budget, certificate,
                              tolerance, maximum_budget) {
  wanted <- budget * (certificate / tolerance)^2 * 1.25
  next_budget <- as.integer(500L * ceiling(wanted / 500L))
  next_budget <- min(max(next_budget, budget + 500L),
                     as.integer(maximum_budget))
  next_validation <- as.integer(
    500L * ceiling(validation_budget * next_budget / budget / 500L)
  )
  list(budget = next_budget, validation_budget = next_validation)
}

##' Fit and certify Monte Carlo maximum likelihood estimates
##'
##' Runs certified stages with growing banks until the certificate is within
##' `tolerance`, the model is found to be weakly identified, no valid fit can
##' be obtained, or the largest budget is reached. If the global proposal
##' fails anywhere within a stage, the stage is repeated with the Laplace
##' proposal, which is then kept. Stage s uses the random number streams
##' 3s - 2 (fitting bank), 3s - 1 and 3s (validation banks) after `seed`.
##' @noRd
fit_adaptive_mcml <- function(engine, start, objective = engine$evaluate,
                              penalty = NULL, initial_budget = 1000L,
                              validation_budget = 5000L,
                              maximum_budget = 16000L, maximum_stages = 4L,
                              tolerance = 0.1, seed = 1L,
                              lower = NULL, upper = NULL,
                              variance_index = integer(),
                              linked_index = list(),
                              boundary_log_variance = -9,
                              range_index = NULL,
                              minimum_log_range = -Inf,
                              near_zero_log_variance = -4) {
  if (is.null(lower)) lower <- rep(-10, length(start))
  if (is.null(upper)) upper <- rep(10, length(start))
  validation_budget <- as.integer(validation_budget)
  limits <- list(
    variance_index = variance_index, linked_index = linked_index,
    boundary_log_variance = boundary_log_variance, range_index = range_index,
    minimum_log_range = minimum_log_range,
    near_zero_log_variance = near_zero_log_variance, tolerance = tolerance
  )
  stream <- function(offset) {
    as.integer(((seed - 1 + offset) %% .Machine$integer.max) + 1)
  }
  total_started <- proc.time()[["elapsed"]]
  history <- list()
  current <- start
  budget <- as.integer(initial_budget)
  for (stage in seq_len(maximum_stages)) {
    run_stage <- function() {
      engine$reset_proposal_failures()
      certified_stage(engine, objective, penalty, current, budget,
                      validation_budget, stream(3L * stage - 2:0), lower,
                      upper, stage == 1L, limits)
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
    lower <- record$lower
    current <- record$estimate
    if (!identical(record$outcome, "not certified") ||
        budget >= maximum_budget || stage == maximum_stages) break
    budgets <- next_mcml_budgets(budget, validation_budget,
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
    certified = isTRUE(last$certified), certificate = last$certificate,
    reason = last$reason, budget = last$budget,
    validation_budget = last$validation_budget,
    proposal = engine$proposal_type(),
    boundary_parameters = last$boundary_parameters,
    range_held = isTRUE(last$range_held), score = score,
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
