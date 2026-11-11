##' Observation log-likelihood, score and curvature in the linear predictor
##'
##' The `terms()` of the family object (`make_family()`). Kept for the callers
##' that steps 1b and 1c of #253 move to the object, and removed in 1d.
##' `curvature` is the negative second derivative of the log-likelihood.
##' @noRd
laplace_observation_terms <- function(eta, y, units_m, family, linkf) {
  family_of_link(family, linkf)$terms(eta, y, units_m)
}

##' Index structure mapping latent variables to observations
##'
##' `columns` holds, for each observation, the latent columns that enter its
##' linear predictor (its location and one level per random effect). The
##' data part of the Hessian, sum_i w_i a_i a_i', touches only the pairs of
##' those columns, which are precomputed so that it can be accumulated with a
##' single `rowsum()` instead of sparse matrix products.
##' @noRd
make_latent_map <- function(columns, dimension) {
  columns <- as.matrix(columns)
  storage.mode(columns) <- "integer"
  width <- ncol(columns)
  pair_rows <- columns[, rep(seq_len(width), each = width), drop = FALSE]
  pair_cols <- columns[, rep(seq_len(width), times = width), drop = FALSE]
  linear <- as.vector(pair_rows) + (as.vector(pair_cols) - 1) * dimension
  observation <- rep(seq_len(nrow(columns)), times = width^2)
  # A random-effect level integrated out of an observation has no latent
  # column, marked 0 (plan_integrated_records()), and enters no pair.
  if (any(columns == 0L)) {
    live <- as.vector(pair_rows) > 0L & as.vector(pair_cols) > 0L
    linear <- linear[live]
    observation <- observation[live]
  }
  list(columns = columns, dimension = dimension,
       pair_linear = linear, pair_observation = observation)
}

latent_linear_predictor <- function(latent_map, latent) {
  columns <- latent_map$columns
  eta <- latent[columns[, 1L]]
  if (ncol(columns) > 1L) {
    for (j in 2:ncol(columns)) {
      index <- columns[, j]
      eta <- eta + if (any(index == 0L)) c(0, latent)[index + 1L] else
        latent[index]
    }
  }
  eta
}

latent_gradient <- function(latent_map, score) {
  out <- numeric(latent_map$dimension)
  for (j in seq_len(ncol(latent_map$columns))) {
    sums <- rowsum(score, latent_map$columns[, j], reorder = FALSE)
    index <- as.integer(rownames(sums))
    live <- index > 0L
    out[index[live]] <- out[index[live]] + sums[live, 1L]
  }
  out
}

##' The latent samples of each observation's random-effect level
##'
##' Columns `index` of the matrix of latent `samples`, with zeros where an
##' observation's level is integrated out (column 0).
##' @noRd
latent_level_samples <- function(samples, index) {
  if (!any(index == 0L)) return(samples[, index, drop = FALSE])
  cbind(0, samples)[, index + 1L, drop = FALSE]
}

latent_data_hessian <- function(latent_map, curvature) {
  sums <- rowsum(curvature[latent_map$pair_observation],
                 latent_map$pair_linear, reorder = FALSE)
  out <- numeric(latent_map$dimension^2)
  out[as.numeric(rownames(sums))] <- sums[, 1L]
  dim(out) <- c(latent_map$dimension, latent_map$dimension)
  out
}

##' Conditional mode of the latent field by Newton's method
##'
##' Maximises log f(y | u) - u' P u / 2, where `precision` is the prior
##' precision P and `latent_map` maps latent variables to observations. Every
##' linear solve uses a Cholesky factor, and the search starts from `start`
##' (the previous mode), so successive calls along an optimisation path need
##' very few iterations. The objective is log-concave for canonical links. For
##' other links the Hessian can be indefinite away from a mode, wherever the
##' likelihood is convex, and a Cholesky factor does not exist. The negative
##' curvatures are then set to zero, which leaves the prior precision and a
##' positive semi-definite sum and so a direction of ascent, and step halving
##' keeps every iteration an ascent step. The Hessian returned is positive
##' definite for the same reason. A custom link's search starts off the limits of
##' the link (`move_start_off_limits()`). `reached_limit` says that the search
##' stopped at `maximum_iterations` without meeting its stopping rule. `kernel`
##' is the family object (`make_family()`), whose properties say which of these
##' a link takes.
##'
##' Far in the lower tail the curvature of the logit link underflows to 0 / 0,
##' below a linear predictor of -372, so a step to a point where the terms
##' are not finite is not accepted. A search fails when the terms, the Newton
##' step or the Hessian are not finite, or when no step increases the objective
##' at a point that is not a mode, as from a previous mode of -507 on a survey
##' at 1% prevalence. A failed search from a nonzero `start` is repeated from
##' zero, and `restarted` says so. If the search from zero fails too there is
##' no mode to return, and the function stops with an error.
##' @noRd
find_latent_mode <- function(mu, y, units_m, latent_map, precision, kernel,
                             start, tolerance = 1e-10,
                             maximum_iterations = 100L) {
  objective <- function(latent, terms) {
    terms$log_likelihood - 0.5 * sum(latent * as.numeric(precision %*% latent))
  }
  # The terms, or NULL where they are not finite and so cannot guide a step.
  terms_at <- function(latent) {
    terms <- kernel$terms(
      mu + latent_linear_predictor(latent_map, latent), y, units_m
    )
    finite <- !is.null(terms) && is.finite(terms$log_likelihood) &&
      all(is.finite(terms$score)) && all(is.finite(terms$curvature))
    if (finite) terms
  }
  hessian_at <- function(curvature) {
    precision + latent_data_hessian(latent_map, curvature)
  }
  # The Hessian and its upper factor, for a custom link with any negative
  # curvatures set to zero if the Hessian is not positive definite.
  factorise_hessian <- function(curvature) {
    hessian <- hessian_at(curvature)
    root <- if (kernel$assumes_concave) chol(hessian) else {
      tryCatch(chol(hessian), error = function(error) NULL)
    }
    if (is.null(root)) {
      hessian <- hessian_at(pmax(curvature, 0))
      root <- chol(hessian)
    }
    list(hessian = hessian, root = root)
  }
  start_at <- function(latent) {
    if (!kernel$moves_start) return(latent)
    kernel$move_start(
      latent,
      function(point) mu + latent_linear_predictor(latent_map, point),
      y, units_m, latent_map$columns[, 1L]
    )
  }
  # The search from one start, or NULL if it fails.
  search_from <- function(start) {
    latent <- start_at(start)
    terms <- terms_at(latent)
    if (is.null(terms)) return(NULL)
    current <- objective(latent, terms)
    if (!is.finite(current)) return(NULL)
    reached_limit <- TRUE
    for (iteration in seq_len(maximum_iterations)) {
      gradient <- latent_gradient(latent_map, terms$score) -
        as.numeric(precision %*% latent)
      root <- factorise_hessian(terms$curvature)$root
      step <- backsolve(root, forwardsolve(t(root), gradient))
      if (!all(is.finite(step))) return(NULL)
      accepted <- FALSE
      for (halving in 0:30) {
        candidate <- latent + step / 2^halving
        candidate_terms <- terms_at(candidate)
        if (is.null(candidate_terms)) next
        value <- objective(candidate, candidate_terms)
        if (is.finite(value) && value >= current - 1e-12 * abs(current)) {
          accepted <- TRUE
          break
        }
      }
      if (!accepted) {
        # The full step promises an increase of half of gradient' step. A
        # point where that is within the rounding allowed above is a mode
        # that rounding keeps the search from confirming, and any other point
        # is not one.
        if (sum(gradient * step) / 2 > 1e-12 * abs(current)) return(NULL)
        reached_limit <- FALSE
        break
      }
      change <- max(abs(candidate - latent))
      latent <- candidate
      terms <- candidate_terms
      current <- value
      if (change <= tolerance * max(1, max(abs(latent)))) {
        reached_limit <- FALSE
        break
      }
    }
    hessian <- if (kernel$assumes_concave) {
      hessian_at(terms$curvature)
    } else {
      factorise_hessian(terms$curvature)$hessian
    }
    if (!all(is.finite(hessian))) return(NULL)
    list(mode = latent, hessian = (hessian + t(hessian)) / 2,
         iterations = iteration, reached_limit = reached_limit)
  }
  found <- NULL
  restarted <- FALSE
  if (!isTRUE(all(start == 0))) {
    found <- tryCatch(search_from(start), error = function(error) NULL)
    restarted <- is.null(found)
  }
  if (is.null(found)) found <- search_from(numeric(length(start)))
  if (is.null(found)) {
    stop("The search for the mode of the latent effects failed.",
         call. = FALSE)
  }
  c(found, list(restarted = restarted))
}

##' Conditional mode of the latent field for a family string and a link
##'
##' `find_latent_mode()` with the family object of `family` and `linkf`. Kept
##' for the callers that step 1 of #253 has not moved yet, and removed in 1d.
##' @noRd
adaptive_laplace_mode <- function(mu, y, units_m, latent_map, precision,
                                  family, linkf, start,
                                  tolerance = 1e-10, maximum_iterations = 100L) {
  find_latent_mode(mu, y, units_m, latent_map, precision,
                   family_of_link(family, linkf), start, tolerance,
                   maximum_iterations)
}

##' Upper Cholesky factor of the inverse of a positive-definite matrix
##'
##' Returns the same factor as `chol(solve(hessian))` without forming the
##' inverse: with J the exchange matrix, chol(J H J) = R gives
##' chol(H^{-1}) = J R^{-T} J. The fixed standard-normal draws are mapped by
##' this factor, so it must agree with the established construction for
##' results to be reproduced. The compiled version uses the same LAPACK and
##' BLAS calls as chol() and backsolve().
##' @noRd
inverse_cholesky_factor <- function(hessian) {
  root <- cpp_inverse_cholesky_factor(as.matrix(hessian))
  if (is.null(root)) {
    stop("The matrix is not positive definite.", call. = FALSE)
  }
  root
}

##' Rows of a standard-normal bank drawn from the defensive component
##'
##' A share `defensive_weight` of the rows, spread evenly through the bank, is
##' mapped through the prior covariance instead of the Laplace proposal. The
##' rows are fixed by their position, so the draws remain smooth functions of
##' the parameters within a stage. A bank with a single row (the Laplace
##' approximation) has no defensive rows.
##' @noRd
defensive_rows <- function(n_draws, defensive_weight) {
  draw <- seq_len(n_draws)
  floor(draw * defensive_weight) > floor((draw - 1) * defensive_weight)
}

##' Centre and precision of the Gaussian proposal at one parameter value
##'
##' With `type = "global"` and an `anchor`, the converged factors at another
##' parameter value, the proposal is three plain refits from the anchor
##' (`refit_from_anchor()`), so given the anchor it is a smooth function of the
##' parameters and the conditional mode is not needed. If a refit fails, the
##' Laplace approximation at the mode found by `find_mode()` is used. Without
##' an anchor, the global iteration (`fit_global_proposal()`) starts from the
##' conditional mode and runs to its tolerance, and if it fails or does not
##' converge the Laplace approximation is used. Returns the centre and
##' precision, the upper factor of the proposal covariance when the global fit
##' computed it (otherwise NULL), the global factors (NULL for the Laplace
##' approximation), whether the global fit `failed` and, if it did, the
##' `failure`: the `reason`, the number of `refits` and the last `change` of
##' the factors. A proposal refitted from the anchor also returns that last
##' change as `refit_change`.
##'
##' With `stall_tolerance`, which the engine sets for the integrated form only,
##' an iteration that reaches its limit of refits with a last change below that
##' tolerance is accepted as it stands, not counted as a failure, and its last
##' change is returned as `stalled_change`. Such an iteration has settled to
##' within a small fraction of its factors, closer than the fixed refits from
##' an anchor routinely leave them. This is safe because the certificate rests
##' on validation banks drawn after the proposal is fixed. A nearly settled
##' anchor therefore changes how efficient the importance sampler is, not
##' whether its estimate or its certificate is valid. Without it the stage would
##' repeat with the Laplace proposal, which on surveys with a large spatial
##' variance leaves almost no effective draws.
##' @noRd
build_gaussian_proposal <- function(type, mu, y, units_m, latent_map,
                                    precision, family, linkf, canonical,
                                    anchor, find_mode, stall_tolerance = NULL) {
  approximation <- NULL
  global <- NULL
  failed <- FALSE
  stalled <- FALSE
  if (identical(type, "global") && !is.null(anchor)) {
    refitted <- refit_from_anchor(anchor, mu, y, units_m, latent_map,
                                  precision, family, linkf, canonical)
    if (!refitted$failed) {
      return(list(centre = refitted$centre, precision = refitted$precision,
                  root = refitted$root, factors = refitted$factors,
                  failed = FALSE, failure = NULL,
                  refit_change = refitted$change))
    }
    approximation <- find_mode()
    return(list(centre = approximation$mode,
                precision = approximation$hessian, root = NULL,
                factors = NULL, failed = TRUE,
                failure = list(reason = "a refit failed",
                               refits = refitted$refits,
                               change = refitted$change)))
  }
  if (identical(type, "global")) {
    approximation <- find_mode()
    global <- fit_global_proposal(
      mode = approximation$mode, mu = mu, y = y, units_m = units_m,
      latent_map = latent_map, precision = precision, family = family,
      linkf = linkf, canonical = canonical
    )
    stalled <- !is.null(global) && !global$converged &&
      !is.null(stall_tolerance) && is.finite(global$change) &&
      global$change < stall_tolerance
    failed <- is.null(global) || (!global$converged && !stalled)
  }
  failure <- if (failed) {
    list(reason = if (is.null(global)) "a refit failed" else "no convergence",
         refits = global$iterations %||% NA_integer_,
         change = global$change %||% NA_real_)
  }
  if (is.null(global)) {
    if (is.null(approximation)) approximation <- find_mode()
    return(list(centre = approximation$mode,
                precision = approximation$hessian, root = NULL,
                factors = NULL, failed = failed, failure = failure))
  }
  list(centre = global$centre, precision = global$precision,
       root = global$root, factors = global$factors, failed = failed,
       failure = failure, stalled_change = if (stalled) global$change)
}

##' Log-likelihood of every row of a matrix of linear predictors
##'
##' The `log_likelihood()` of the family object (`make_family()`). `eta` has one
##' row per draw and one column per observation. A draw whose mean is not
##' finite or is outside the support has log-likelihood -Inf, and other draws
##' are not affected. Kept for its callers outside the engine, and removed in
##' step 1d of #253.
##' @noRd
compute_observation_log_likelihood <- function(eta, y, units_m, family,
                                               linkf) {
  family_of_link(family, linkf)$log_likelihood(eta, y, units_m)
}

##' Deviations of draws from the centre under the prior of the latent effects
##'
##' The spatial block through the prior Cholesky factor, random-effect levels
##' scaled by their standard deviations. `latent` gives the latent columns of
##' the spatial block and of each random effect.
##' @noRd
compute_prior_deviations <- function(draws, proposal, latent) {
  deviations <- matrix(0, nrow(draws), ncol(draws))
  deviations[, latent$spatial] <-
    draws[, latent$spatial, drop = FALSE] %*% proposal$spatial_root
  for (j in seq_along(latent$random)) {
    index <- latent$random[[j]]
    deviations[, index] <- draws[, index, drop = FALSE] *
      sqrt(proposal$values$sigma2_re[j])
  }
  deviations
}

##' Log density of draws under the defensive mixture proposal
##'
##' Squared norms of the draws standardised under each mixture component. A
##' draw's residual under its own component is its standard-normal vector, so
##' only the other component needs a triangular solve. For the Gaussian
##' proposal's draws the spatial residuals under the prior can be taken from
##' the samples already standardised for the prior density
##' (`standardised_spatial`). The realised share of defensive draws is used as
##' the mixture weight, so the estimate is unbiased for any bank size. A bank
##' evaluated in row blocks passes the share of the whole bank.
##' @noRd
compute_log_mixture_proposal <- function(draws, deviations, defensive,
                                         proposal, latent,
                                         standardised_spatial = NULL,
                                         share = mean(defensive)) {
  draw_norms <- rowSums(draws^2)
  log_normalising <- latent$dimension * log(2 * pi)
  log_laplace_determinant <- 2 * sum(log(diag(proposal$root)))
  if (share == 0) {
    # The latent values are centre + draws %*% root, so their standardised
    # residuals are the draws themselves.
    return(-0.5 * (log_normalising + log_laplace_determinant + draw_norms))
  }
  laplace_norms <- draw_norms
  laplace_norms[defensive] <- colSums(forwardsolve(
    t(proposal$root), t(deviations[defensive, , drop = FALSE])
  )^2)
  laplace_draws <- !defensive
  spatial_norms <- if (is.null(standardised_spatial)) {
    colSums(forwardsolve(
      t(proposal$spatial_root),
      t(deviations[laplace_draws, latent$spatial, drop = FALSE])
    )^2)
  } else {
    colSums((standardised_spatial[, laplace_draws, drop = FALSE] -
               as.numeric(forwardsolve(t(proposal$spatial_root),
                                       proposal$centre[latent$spatial])))^2)
  }
  for (j in seq_along(latent$random)) {
    index <- latent$random[[j]]
    spatial_norms <- spatial_norms +
      rowSums(deviations[laplace_draws, index, drop = FALSE]^2) /
      proposal$values$sigma2_re[j]
  }
  prior_norms <- draw_norms
  prior_norms[laplace_draws] <- spatial_norms
  log_prior_determinant <- 2 * sum(log(diag(proposal$spatial_root))) +
    if (length(latent$random)) {
      sum(lengths(latent$random) * log(proposal$values$sigma2_re))
    } else {
      0
    }
  log_laplace <- log1p(-share) - 0.5 * (
    log_normalising + log_laplace_determinant + laplace_norms
  )
  log_defensive <- log(share) - 0.5 * (
    log_normalising + log_prior_determinant + prior_norms
  )
  largest <- pmax(log_laplace, log_defensive)
  largest + log(exp(log_laplace - largest) + exp(log_defensive - largest))
}

##' Log importance weights of a bank of draws, computed in R
##'
##' Used for links without the compiled version, `cpp_importance_log_weights()`,
##' which reproduces this computation exactly (to the last bits where the
##' compiler fuses multiply and add). `columns` holds each
##' observation's latent columns (location, then random-effect levels) and
##' `log_likelihood` maps a matrix of linear predictors to the log-likelihood
##' of each row.
##' @noRd
compute_log_weights_in_r <- function(draws, proposal, defensive,
                                     defensive_share, values, spatial_root,
                                     same_parameters, linear_offset, latent,
                                     columns, log_likelihood) {
  deviations <- draws %*% proposal$root
  if (any(defensive)) {
    deviations[defensive, ] <- compute_prior_deviations(
      draws[defensive, , drop = FALSE], proposal, latent
    )
  }
  samples <- deviations + rep(proposal$centre, each = nrow(deviations))
  spatial_samples <- samples[, latent$spatial, drop = FALSE]
  # Samples standardised by the prior factor. Their squared norms give the
  # prior quadratic form and, when the target and proposal parameters
  # coincide, they also give the residuals under the defensive component.
  standardised_samples <- forwardsolve(t(spatial_root), t(spatial_samples))
  log_prior <- -0.5 * (
    length(latent$spatial) * log(2 * pi) + 2 * sum(log(diag(spatial_root))) +
      colSums(standardised_samples^2)
  )
  eta <- spatial_samples[, columns[, 1L], drop = FALSE] +
    rep(linear_offset, each = nrow(samples))
  for (j in seq_along(latent$random)) {
    random_samples <- samples[, latent$random[[j]], drop = FALSE]
    eta <- eta + latent_level_samples(samples, columns[, j + 1L])
    log_prior <- log_prior - 0.5 * (
      length(latent$random[[j]]) * log(2 * pi * values$sigma2_re[j]) +
        rowSums(random_samples^2) / values$sigma2_re[j]
    )
  }
  log_target <- log_prior + log_likelihood(eta)
  log_proposal <- compute_log_mixture_proposal(
    draws, deviations, defensive, proposal, latent,
    if (same_parameters) standardised_samples, share = defensive_share
  )
  log_target - log_proposal
}

##' Terms of the slope of the log-likelihood in a variance at zero
##'
##' For an effect of variance v shared by the observations of each `group`, the
##' derivative of the log-likelihood at v = 0 is the mean, under the conditional
##' distribution of the other latent effects, of half the sum over groups of
##' (sum k')^2 + sum k'', with k the log-likelihood of an observation in its
##' linear predictor. Returns that quantity for each row of the matrix `eta` of
##' linear predictors. The family object `kernel` (`make_family()`) gives k' and
##' k'' (`slope_parts()`), and NA marks a draw where a custom link cannot give
##' them. Concavity is not needed. Records with the nugget integrated out take
##' them from `derivatives`, the first and second derivatives `h1` and `h2` of
##' `cpp_site_evaluate()`, one column per record, whatever the link.
##' @noRd
compute_slope_terms <- function(eta, y, units_m, kernel, group,
                                derivatives = NULL) {
  if (!is.null(derivatives)) {
    score <- derivatives$h1
    curvature <- -derivatives$h2
  } else {
    parts <- kernel$slope_parts(eta, y, units_m)
    if (is.null(parts)) return(rep(NA_real_, nrow(eta)))
    score <- parts$score
    curvature <- parts$curvature
  }
  grouped <- rowsum(t(score), group, reorder = FALSE)
  0.5 * (colSums(grouped^2) - rowSums(curvature))
}

##' Terms of the slope at a held variance for a family string and a link
##'
##' `compute_slope_terms()` with the family object of `family` and `linkf`,
##' NULL for the canonical link. Kept for its callers outside the engine, and
##' removed in step 1d of #253.
##' @noRd
compute_variance_slope_terms <- function(eta, y, units_m, family, group,
                                         linkf = NULL, derivatives = NULL) {
  compute_slope_terms(eta, y, units_m, make_family(family, linkf), group,
                      derivatives)
}

##' Stages whose global iteration was accepted although it stopped at its limit
##'
##' From the stage records of `fit_adaptive_mcml()`, the `stages` with
##' `stalled_anchors`, how many iterations were accepted in all (`count`) and the
##' largest last change (`change`). NULL when no stage has one, as in every fit
##' of the sampled form.
##' @noRd
summarise_stalled_anchors <- function(history) {
  counts <- vapply(history, function(stage) {
    as.numeric(stage$stalled_anchors %||% 0)
  }, numeric(1))
  if (!any(counts > 0)) return(NULL)
  changes <- vapply(history, function(stage) {
    as.numeric(stage$stalled_change %||% NA_real_)
  }, numeric(1))
  list(stages = which(counts > 0), count = sum(counts),
       change = max(changes, na.rm = TRUE))
}

##' Construct the adaptive importance-sampling engine
##'
##' The standard-normal bank is held fixed within a stage, while the Gaussian
##' map from the bank to the latent effects changes smoothly with the trial
##' parameters. This preserves a deterministic common-random-number objective
##' without tying the proposal to the user's starting values.
##'
##' The importance distribution is the defensive mixture
##' (1 - a) N(centre, Q^-1) + a N(centre, Sigma), where Sigma is the prior
##' covariance of the latent effects and a = `defensive_weight`. With
##' `proposal = "global"` (the default) the centre and precision Q come from
##' the global Gaussian approximation of `fit_global_proposal()`, which is
##' fitted to the whole conditional distribution of the latent effects. With
##' `proposal = "laplace"`, or where the global fit fails, they are the
##' conditional mode and the curvature there (the Laplace approximation).
##' Failures of the global fit are counted, so that a fit can be repeated with
##' the Laplace proposal throughout. The count can be marked and restored to
##' the mark, which sets aside the failures of the check for other maxima
##' (`set_aside_check_failures()`). With a Gaussian proposal alone the
##' importance weights have infinite variance whenever the largest eigenvalue of
##' Sigma^{1/2} W Sigma^{1/2} exceeds one (W the likelihood curvature), which
##' holds in practically every spatial dataset (Pitt et al., 2013). The mixture
##' bounds every moment of the weights, so that the Monte Carlo error estimates
##' used for stopping are valid, at a small cost in efficiency. Weights use the
##' full mixture density. With `defensive_weight = 0` the proposal is the
##' Gaussian component alone.
##'
##' With `integrate_noise`, the nugget and the levels of random effects that
##' lie within one location are integrated out of the likelihood of the records
##' by the rules of `plan_integrated_records()` (`make_integrated_kernel()`).
##' `integrate_noise` is `TRUE` for every term the rules allow and `FALSE` for
##' none. The latent effects are then the
##' spatial effect at the locations and the sampled random-effect levels, with
##' the nugget left in the prior only at the locations where it is sampled. The
##' model and its parameters are the same. An integrated variance whose log
##' parameter is at `hold_log_ratio` is zero, and the records it entered then
##' take the plain kernel or the rest of their integrated variance.
##' `spatial_covariance()` gives the covariance of the location effects with
##' the whole nugget in either form, `integrated_terms` the terms integrated
##' out at one level or more, `integration_counts` how many levels of each
##' term are integrated out, and `hold_indices` the parameters held at exactly
##' zero (`hold_noise_at_zero()`). In this form a
##' global iteration that stops at its limit of refits with a last change below
##' 0.01 is accepted rather than counted as a failure (`build_gaussian_proposal()`),
##' and `stalled_anchors()` gives how many were accepted since the failures were
##' last reset and the largest such change. The sampled form keeps the rule of
##' the release, so its fits are unchanged.
##'
##' The prior factors, the global proposal's dense algebra and, for canonical
##' links, the importance weights are computed in compiled code, which
##' reproduces the R computations (`compute_log_weights_in_r()` for the
##' weights) exactly, or to the last bits where the compiler fuses multiply
##' and add, but keeps their large working matrices off the R heap.
##' @noRd
make_adaptive_mcml_engine <- function(
    D, y, units_m, distances, family = c("binomial", "poisson"),
    kappa = 0.5, cov_offset = NULL, ID_coords = NULL,
    fix_tau2 = FALSE, ID_re = NULL, invlink = NULL,
    defensive_weight = 0.05, proposal = c("global", "laplace"),
    integrate_noise = FALSE) {
  family <- match.arg(family)
  default_proposal_type <- match.arg(proposal)
  proposal_failures <- 0L
  # The largest change of a factor at the last of the fixed refits of the
  # global proposals since the failures were last reset.
  largest_refit_change <- -Inf
  # Global iterations accepted although they stopped at their limit of refits,
  # in the integrated form only, and the largest last change among them, since
  # the failures were last reset.
  stalled_anchors <- 0L
  largest_stalled_change <- -Inf
  first_failure <- NULL
  # Proposals built so far, which date the cached one for
  # restore_proposal_failures().
  builds <- 0L
  # Errors while a proposal is built, which make an evaluation fail without
  # failing the global iteration, are counted apart from those failures.
  proposal_errors <- 0L
  first_error <- NULL
  mode_searches <- 0L
  mode_searches_at_limit <- 0L
  # Searches that failed from the previous mode and were repeated from zero.
  mode_searches_restarted <- 0L
  n_loc <- nrow(distances)
  p <- ncol(D)
  if (is.null(cov_offset)) cov_offset <- numeric(length(y))
  if (is.null(ID_coords)) ID_coords <- seq_len(length(y))
  if (!is.null(ID_re)) ID_re <- as.matrix(ID_re)
  # The terms that can be integrated out: the nugget, then each random effect.
  has_nugget <- isTRUE(fix_tau2) ||
    (is.numeric(fix_tau2) && isTRUE(fix_tau2 > 0))
  re_terms <- if (!is.null(ID_re)) {
    paste0("re(", colnames(ID_re) %||% paste0("effect", seq_len(ncol(ID_re))),
           ")")
  }
  term_names <- c("nugget", re_terms)
  allowed <- stats::setNames(rep(isTRUE(integrate_noise), length(term_names)),
                             term_names)
  allowed[["nugget"]] <- allowed[["nugget"]] && has_nugget
  plan <- if (any(allowed)) {
    plan_integrated_records(y, units_m, D, cov_offset, ID_coords, ID_re,
                            n_loc, allowed)
  }
  integrated <- isTRUE(plan$any)
  # The last change below which a global iteration that stops at its limit
  # of refits is accepted, in the integrated form (build_gaussian_proposal()).
  stall_tolerance <- 1e-2
  # Whether the nugget is integrated out at any location, and where it stays
  # sampled in the prior.
  nugget_integrated <- integrated && isTRUE(allowed[["nugget"]]) &&
    any(plan$nugget_integrated)
  nugget_sampled <- if (nugget_integrated) !plan$nugget_integrated
  n_re <- if (is.null(ID_re)) 0L else ncol(ID_re)
  if (integrated) {
    y <- plan$y
    units_m <- plan$units_m
    D <- plan$D
    cov_offset <- plan$cov_offset
    ID_coords <- plan$ID_coords
    ID_re <- plan$ID_re
  }
  n <- length(y)
  n_dim_re <- if (!n_re) integer() else if (integrated) {
    as.integer(plan$sampled_levels)
  } else {
    apply(ID_re, 2L, max)
  }
  layout <- cov_par_layout(p, fix_tau2 = fix_tau2,
                           sigma2_me = FALSE, n_re = n_re)

  linkf <- make_inverse_link_functions(family, invlink, n)
  model <- make_family(family, linkf)
  # The proposal (build_gaussian_proposal()) still takes the test of the
  # canonical link, until step 1b of #253 moves it to the object.
  canonical <- model$canonical
  sites <- if (integrated) {
    make_site_table(y, units_m, family, linkf = linkf, kind = plan$kind,
                    a = plan$a, b = plan$b, patterns = plan$patterns)
  }
  # A variance whose log parameter is at this value is held at exactly zero,
  # and the records it entered lose it from their integrated variance.
  hold_log_ratio <- -20
  # The integrated terms in the order of the plan's classes, with the index of
  # each one's log parameter.
  term_index <- c(if (isTRUE(fix_tau2)) layout$index$nu2 else NA_integer_,
                  if (n_re) layout$index$sigma2_re)
  hold_indices <- if (integrated) {
    held <- term_index[term_names %in% plan$integrated_terms]
    held[!is.na(held)]
  }

  unpack <- function(parameters) {
    sigma2 <- exp(parameters[layout$index$sigma2])
    nugget <- nugget_ratio(
      fix_tau2, sigma2,
      if (!is.null(layout$index$nu2)) parameters[layout$index$nu2] else NULL
    )
    values <- list(
      beta = parameters[layout$index$beta], sigma2 = sigma2,
      phi = exp(parameters[layout$index$phi]), nugget_ratio = nugget,
      sigma2_re = if (n_re) exp(parameters[layout$index$sigma2_re]) else NULL
    )
    if (integrated) {
      # The variance of each term, zero where it is held, the variances of the
      # classes of integrated records, and the nugget ratio left in the prior
      # at the locations where the nugget stays sampled.
      held <- !is.na(term_index) &
        parameters[pmax(term_index, 1L)] <= hold_log_ratio
      variances <- c(if (has_nugget) sigma2 * nugget else 0,
                     values$sigma2_re)
      variances[held] <- 0
      values$term_variances <- stats::setNames(variances, term_names)
      values$class_variances <- if (nrow(plan$classes)) {
        as.numeric((plan$classes + 0) %*% variances)
      } else {
        numeric()
      }
      if (nugget_integrated) {
        values$prior_ratio <- if (held[1L] || !any(nugget_sampled)) 0 else
          nugget * nugget_sampled
      }
    }
    values
  }
  # The nugget ratio in the prior of the latent effects.
  latent_ratio <- function(values) {
    if (nugget_integrated) values$prior_ratio else values$nugget_ratio
  }

  half_integer_kappa <- kappa %in% c(0.5, 1.5, 2.5)
  distance_values <- as.numeric(distances)
  # The prior covariance of the latent location effects, or with `marginal`
  # that of the location effects with the whole nugget, which is the same in
  # both forms.
  covariance_at <- function(values, marginal = FALSE) {
    # Evaluate the closed-form half-integer Matern correlation directly on the
    # stored distance matrix rather than rebuilding it from a dist object.
    correlation <- if (half_integer_kappa) {
      matrix(cpp_half_integer_matern(distance_values, values$phi, kappa),
             n_loc, n_loc)
    } else {
      matern_correlation(stats::as.dist(distances), phi = values$phi,
                         kappa = kappa, return_sym_matrix = TRUE)
    }
    diag(correlation) <- diag(correlation) +
      if (marginal) values$nugget_ratio else latent_ratio(values)
    values$sigma2 * correlation
  }

  # A random effect integrated out at every level has no latent columns.
  latent_layout <- function() {
    indices <- vector("list", n_re)
    first <- integer(n_re)
    next_index <- n_loc + 1L
    if (n_re) for (j in seq_len(n_re)) {
      indices[[j]] <- next_index - 1L + seq_len(n_dim_re[j])
      first[j] <- next_index
      next_index <- next_index + n_dim_re[j]
    }
    list(spatial = seq_len(n_loc), random = indices, first = first,
         dimension = next_index - 1L)
  }
  latent <- latent_layout()
  random_first <- latent$first
  random_size <- as.integer(n_dim_re)

  # The spatial covariance at one parameter value, its upper Cholesky factor
  # and the precision. For the closed-form Matern correlations they are
  # computed in compiled code, without temporary copies on the R heap. A
  # covariance that needs jitter is factorised by factor_covariance(), as for
  # other values of kappa.
  # The last factorisation is kept, because an optimiser's finite
  # differences in the regression coefficients leave the covariance
  # unchanged. One entry bounds the memory used for large surveys.
  prior_cache <- NULL
  spatial_prior_at <- function(values) {
    key <- c(values$phi, values$sigma2, latent_ratio(values))
    if (!is.null(prior_cache) && identical(key, prior_cache$key)) {
      return(prior_cache$prior)
    }
    prior <- compute_spatial_prior(values)
    prior_cache <<- list(key = key, prior = prior)
    prior
  }
  compute_spatial_prior <- function(values) {
    if (half_integer_kappa) {
      prior <- cpp_spatial_prior_factors(
        distance_values, n_loc, values$phi, kappa, values$sigma2,
        latent_ratio(values)
      )
      if (!is.null(prior$root)) return(prior)
      covariance <- prior$covariance
    } else {
      covariance <- covariance_at(values)
    }
    root <- factor_covariance(covariance, "adaptive MCML target")
    list(covariance = covariance, root = root, precision = chol2inv(root))
  }

  # Index map from latent variables (spatial locations, then random-effect
  # levels) to observations, built once.
  latent_map <- make_latent_map(
    if (integrated) plan$latent_columns else
      cbind(ID_coords, if (n_re) vapply(seq_len(n_re), function(j) {
        latent$random[[j]][ID_re[, j]]
      }, integer(n))),
    latent$dimension
  )
  # The previous conditional mode is an excellent Newton start for nearby
  # parameter values; the most recent factorisations are reused when the
  # proposal and target parameters coincide.
  last_mode <- numeric(latent$dimension)
  # The likelihood of the records at one parameter value: the link, or with
  # the nugget integrated out the records' table at the nugget variance.
  kernel_at <- function(values) {
    if (!integrated) return(model)
    make_family(family, linkf,
                kernel = make_integrated_kernel(sites, values$class_variances))
  }
  # Converged global factors at a fixed parameter value (set by
  # set_anchor()), used to start the global iteration at other values.
  anchor <- NULL
  cached_parameters <- NULL
  cached_type <- NULL
  cached_proposal <- NULL
  cached_failed <- FALSE
  cached_build <- 0L

  proposal_at <- function(parameters, type = default_proposal_type) {
    if (!is.null(cached_parameters) &&
        identical(parameters, cached_parameters) &&
        identical(type, cached_type)) {
      return(cached_proposal)
    }
    values <- unpack(parameters)
    prior <- spatial_prior_at(values)
    precision <- prior$precision
    if (n_re) {
      precision <- matrix(0, latent$dimension, latent$dimension)
      precision[latent$spatial, latent$spatial] <- prior$precision
      for (j in seq_len(n_re)) {
        index <- latent$random[[j]]
        precision[cbind(index, index)] <- 1 / values$sigma2_re[j]
      }
    }
    mu <- as.numeric(D %*% values$beta + cov_offset)
    kernel <- kernel_at(values)
    # A search that fails stops with an error, so the next one starts from
    # the last mode that was found.
    find_mode <- function() {
      approximation <- find_latent_mode(
        mu = mu, y = y, units_m = units_m, latent_map = latent_map,
        precision = precision, kernel = kernel, start = last_mode
      )
      last_mode <<- approximation$mode
      mode_searches <<- mode_searches + 1L
      if (approximation$reached_limit) {
        mode_searches_at_limit <<- mode_searches_at_limit + 1L
      }
      if (approximation$restarted) {
        mode_searches_restarted <<- mode_searches_restarted + 1L
      }
      approximation
    }
    gaussian <- build_gaussian_proposal(
      type, mu, y, units_m, latent_map, precision, family, kernel$legacy_link,
      canonical && !integrated, anchor, find_mode,
      if (integrated) stall_tolerance
    )
    builds <<- builds + 1L
    if (!is.null(gaussian$stalled_change)) {
      stalled_anchors <<- stalled_anchors + 1L
      largest_stalled_change <<- max(largest_stalled_change,
                                     gaussian$stalled_change)
    }
    if (gaussian$failed) {
      proposal_failures <<- proposal_failures + 1L
      if (is.null(first_failure)) first_failure <<- gaussian$failure
    }
    if (!is.null(gaussian$refit_change)) {
      largest_refit_change <<- max(largest_refit_change, gaussian$refit_change)
    }
    root <- gaussian$root
    if (is.null(root)) {
      root <- tryCatch(
        inverse_cholesky_factor(gaussian$precision),
        error = function(error) factor_covariance(
          solve(gaussian$precision), "adaptive MCML proposal"
        )
      )
    }
    # The compiled factorisation does not stop at entries that are not
    # numbers, so a proposal that is not finite is an error here, which gives
    # the parameter value no weight in evaluate().
    if (!all(is.finite(gaussian$centre)) || !all(is.finite(root))) {
      stop("The importance distribution is not finite.", call. = FALSE)
    }
    proposal <- list(centre = gaussian$centre, root = root, values = values,
                     covariance = prior$covariance,
                     spatial_root = prior$root,
                     type = if (is.null(gaussian$factors)) "laplace" else
                       "global",
                     mu = mu, factors = gaussian$factors)
    cached_parameters <<- parameters
    cached_type <<- type
    cached_proposal <<- proposal
    cached_failed <<- gaussian$failed
    cached_build <<- builds
    proposal
  }

  # `built_proposal`, if given, is the proposal at `proposal_parameters`, so
  # that a bank evaluated in row blocks builds each proposal once. Such a bank
  # also passes its `defensive` rows and their share of the whole bank.
  evaluate <- function(parameters, draws, log_weights_only = FALSE,
                       proposal_parameters = parameters,
                       proposal_type = default_proposal_type,
                       built_proposal = NULL, defensive = NULL,
                       defensive_share = NULL) {
    # Keeping proposal_parameters fixed lets the calibration experiments
    # separate ordinary importance-sampling error from variability caused by
    # rebuilding the proposal as the optimizer moves.
    proposal <- if (is.null(built_proposal)) {
      tryCatch(proposal_at(proposal_parameters, proposal_type),
               error = function(error) {
                 proposal_errors <<- proposal_errors + 1L
                 if (is.null(first_error)) first_error <<- conditionMessage(error)
                 NULL
               })
    } else {
      built_proposal
    }
    # A failed evaluation gives zero weight to every draw, or the penalty
    # value of the objective.
    failed <- if (log_weights_only) rep(-Inf, nrow(draws)) else -1e100
    if (is.null(proposal)) return(failed)
    if (is.null(defensive)) {
      defensive <- defensive_rows(nrow(draws), defensive_weight)
    }
    if (is.null(defensive_share)) defensive_share <- mean(defensive)
    same_parameters <- identical(parameters, proposal_parameters)
    if (same_parameters) {
      values <- proposal$values
      spatial_root <- proposal$spatial_root
    } else {
      values <- unpack(parameters)
      spatial_root <- tryCatch(
        factor_covariance(covariance_at(values), "adaptive MCML target"),
        error = function(error) NULL
      )
    }
    if (is.null(spatial_root)) return(failed)
    linear_offset <- as.numeric(D %*% values$beta + cov_offset)
    log_weights <- if (model$compiled) {
      cpp_importance_log_weights(
        draws, proposal$root, proposal$centre, defensive, defensive_share,
        proposal$spatial_root, as.numeric(proposal$values$sigma2_re),
        spatial_root, as.numeric(values$sigma2_re), same_parameters,
        linear_offset, latent_map$columns, random_first, random_size, y,
        units_m, model$poisson, sites$pointer,
        values$class_variances %||% numeric(),
        sites$units %||% matrix(0L, 0L, 4L)
      )
    } else {
      compute_log_weights_in_r(
        draws, proposal, defensive, defensive_share, values, spatial_root,
        same_parameters, linear_offset, latent, latent_map$columns,
        function(eta) kernel_at(values)$log_likelihood(eta, y, units_m)
      )
    }
    if (log_weights_only) return(log_weights)
    maximum <- max(log_weights)
    if (model$needs_finite_weight && identical(maximum, -Inf)) return(failed)
    maximum + log(mean(exp(log_weights - maximum)))
  }

  # The groups of the effects whose variance can be tested at zero, by
  # parameter. A nugget is one effect per location. In a nested record a
  # random effect integrated out as an inner term has no group (NA), and its
  # slope comes from the record's own derivative in the inner variance.
  slope_groups <- c(if (!is.null(layout$index$nu2)) {
    stats::setNames(list(ID_coords), layout$index$nu2)
  }, if (n_re) {
    stats::setNames(lapply(seq_len(n_re), function(j) ID_re[, j]),
                    layout$index$sigma2_re)
  })
  inner_records <- if (integrated && n_re) {
    stats::setNames(plan$inner, layout$index$sigma2_re)
  }
  # With the spatial variance at zero the nugget ratio says nothing about the
  # nugget itself, so only the other effects can be tested.
  testable_variances <- function(boundary) {
    tested <- intersect(boundary, as.integer(names(slope_groups)))
    if (layout$index$sigma2 %in% boundary) {
      tested <- setdiff(tested, layout$index$nu2)
    }
    tested
  }
  # One row per draw and one column per parameter in `tested`
  # (compute_slope_terms()), from the same draws as the weights.
  variance_slope_terms <- function(draws, proposal, defensive, tested) {
    deviations <- draws %*% proposal$root
    if (any(defensive)) {
      deviations[defensive, ] <- compute_prior_deviations(
        draws[defensive, , drop = FALSE], proposal, latent
      )
    }
    samples <- deviations + rep(proposal$centre, each = nrow(deviations))
    eta <- samples[, latent_map$columns[, 1L], drop = FALSE] +
      rep(proposal$mu, each = nrow(samples))
    for (j in seq_along(latent$random)) {
      eta <- eta + latent_level_samples(samples, latent_map$columns[, j + 1L])
    }
    classes <- proposal$values$class_variances
    derivatives <- if (integrated && any(classes > 0)) {
      cpp_site_evaluate(sites$pointer, classes, sites$units, eta, FALSE, 2L)
    }
    vapply(tested, function(parameter) {
      group <- slope_groups[[as.character(parameter)]]
      inner <- inner_records[[as.character(parameter)]]
      if (!length(inner)) {
        return(compute_slope_terms(eta, y, units_m, model, group,
                                   derivatives))
      }
      # Records whose term is an inner one take their own slope, and the
      # others the grouped form.
      keep <- !is.na(group)
      grouped <- compute_slope_terms(
        eta[, keep, drop = FALSE], y[keep], units_m[keep], model,
        group[keep],
        if (!is.null(derivatives)) {
          list(h1 = derivatives$h1[, keep, drop = FALSE],
               h2 = derivatives$h2[, keep, drop = FALSE])
        }
      )
      slopes <- cpp_site_inner_slopes(sites$pointer, classes,
                                      sites$units[inner, , drop = FALSE],
                                      eta[, inner, drop = FALSE])
      grouped + rowSums(slopes)
    }, numeric(nrow(draws)))
  }

  list(evaluate = evaluate, proposal_at = proposal_at, layout = layout,
       latent_dimension = latent$dimension, unpack = unpack,
       testable_variances = testable_variances,
       variance_slope_terms = variance_slope_terms,
       spatial_covariance = function(parameters) {
         covariance_at(unpack(parameters), marginal = TRUE)
       },
       integrated = integrated,
       integrated_terms = if (integrated) plan$integrated_terms else
         character(),
       term_names = term_names,
       # Without a nugget only the random effects are counted.
       integration_counts = if (integrated) {
         plan$counts[c(has_nugget, rep(TRUE, n_re)), , drop = FALSE]
       },
       hold_log_ratio = if (length(hold_indices)) hold_log_ratio,
       hold_indices = hold_indices,
       integrated_locations = if (nugget_integrated) {
         sum(plan$nugget_integrated)
       } else {
         0L
       },
       site_table_info = function() {
         if (integrated) cpp_site_table_info(sites$pointer)
       },
       defensive_weight = defensive_weight,
       proposal_type = function() default_proposal_type,
       use_proposal = function(type) {
         default_proposal_type <<- match.arg(type, c("global", "laplace"))
         invisible(NULL)
       },
       proposal_failures = function() proposal_failures,
       first_proposal_failure = function() first_failure,
       largest_refit_change = function() {
         if (is.finite(largest_refit_change)) largest_refit_change else NA_real_
       },
       proposal_errors = function() {
         list(count = proposal_errors, first = first_error)
       },
       linkf = linkf,
       mode_search_counts = function() {
         list(searches = mode_searches, at_limit = mode_searches_at_limit,
              restarted = mode_searches_restarted)
       },
       set_anchor = function(parameters) {
         anchor <<- NULL
         cached_parameters <<- NULL
         proposal <- tryCatch(proposal_at(parameters, "global"),
                              error = function(error) NULL)
         if (!is.null(proposal$factors)) {
           anchor <<- c(list(mu = proposal$mu), proposal$factors)
         }
         cached_parameters <<- NULL
         invisible(!is.null(anchor))
       },
       reset_proposal_failures = function() {
         proposal_failures <<- 0L
         largest_refit_change <<- -Inf
         stalled_anchors <<- 0L
         largest_stalled_change <<- -Inf
         first_failure <<- NULL
         invisible(NULL)
       },
       stalled_anchors = function() {
         list(count = stalled_anchors,
              change = if (is.finite(largest_stalled_change)) {
                largest_stalled_change
              } else {
                NA_real_
              })
       },
       # The failures so far, the largest refit change and the stalled
       # iterations, to be restored after evaluations whose failures are set
       # aside (set_aside_check_failures()).
       mark_proposal_failures = function() {
         list(failures = proposal_failures, first = first_failure,
              builds = builds, refit_change = largest_refit_change,
              stalled = stalled_anchors, stalled_change = largest_stalled_change)
       },
       # A failed proposal built since the mark that is still cached would be
       # reused without being counted, so it is built again when next needed.
       restore_proposal_failures = function(mark) {
         set_aside <- proposal_failures - mark$failures
         proposal_failures <<- mark$failures
         first_failure <<- mark$first
         largest_refit_change <<- mark$refit_change %||% largest_refit_change
         stalled_anchors <<- mark$stalled %||% stalled_anchors
         largest_stalled_change <<- mark$stalled_change %||%
           largest_stalled_change
         if (cached_failed && cached_build > mark$builds) {
           cached_parameters <<- NULL
         }
         invisible(set_aside)
       })
}

##' An engine that switches between the sampled and the integrated nugget
##'
##' With `integrate_noise = "auto"` both forms of the engine are built for a
##' model with a nugget (`sampled` and `integrated`, from
##' `make_adaptive_mcml_engine()`), and each stage of the fit uses one of them.
##' `choose_form()` takes the integrated form when the nugget variance sigma2
##' nu2 at the given parameters is at least `threshold`, and the sampled form
##' otherwise. `fit_adaptive_mcml()` calls it at the start of every stage after
##' the first, whose form `glgpm_nong_adaptive()` sets. Near a zero nugget the
##' likelihood is nearly flat along the nugget, and there the sampled form's
##' estimate of the nugget's score is far less noisy, because its moving
##' proposal follows the nugget in the prior of the latent effects, while in
##' the integrated form the nugget changes the likelihood of every record. With
##' a larger nugget the integrated form keeps far more effective draws. In the
##' larger test of the integrated form it certified 5 runs against 12 of 21 for
##' the sampled form where the nugget variance at the maximum was below 0.1, 9
##' and 9 of 18 between 0.1 and 0.3, and 68 against 28 of 75 above 0.3.
##'
##' A model with an `re()` term that can be integrated out is fitted with all
##' its terms integrated out (`glgpm_nong_adaptive()`), so the engine serves
##' models whose only such term is the nugget.
##'
##' The two forms have the same parameters, layout and latent dimension, so
##' the returned engine has the interface of `make_adaptive_mcml_engine()`.
##' It is an environment, whose `integrated`, `hold_log_ratio`,
##' `hold_indices`, `integrated_terms`, `integration_counts` and
##' `integrated_locations` are those of the form in use. `form()` gives that
##' form, `set_form()` sets it, the proposal type is set in both forms, the
##' failures of both are reset together, and the counts of mode searches and
##' of proposal errors are summed over the two. Every other function is that
##' of the form in use, and the failures marked and restored within a stage
##' are those of its form.
##' @noRd
make_switching_mcml_engine <- function(sampled, integrated, threshold = 0.1) {
  stopifnot(isTRUE(integrated$integrated), !isTRUE(sampled$integrated),
            identical(sampled$latent_dimension, integrated$latent_dimension),
            identical(sampled$layout, integrated$layout))
  engines <- list(sampled = sampled, integrated = integrated)
  current <- "sampled"
  engine <- new.env(parent = emptyenv())
  delegate <- function(name) {
    force(name)
    function(...) engines[[current]][[name]](...)
  }
  for (name in c("evaluate", "proposal_at", "unpack", "testable_variances",
                 "variance_slope_terms", "spatial_covariance",
                 "site_table_info", "proposal_type", "proposal_failures",
                 "first_proposal_failure", "largest_refit_change",
                 "set_anchor", "stalled_anchors", "mark_proposal_failures",
                 "restore_proposal_failures")) {
    assign(name, delegate(name), envir = engine)
  }
  engine$layout <- sampled$layout
  engine$latent_dimension <- sampled$latent_dimension
  engine$defensive_weight <- sampled$defensive_weight
  engine$linkf <- sampled$linkf
  engine$threshold <- threshold
  # Every validation bank is drawn in both forms and combined
  # (combine_validation_banks()), so both forms keep the anchor of the stage.
  engine$combine_validation <- TRUE
  engine$set_anchor <- function(parameters) {
    engines[[setdiff(names(engines), current)]]$set_anchor(parameters)
    engines[[current]]$set_anchor(parameters)
  }
  makeActiveBinding("integrated", function() identical(current, "integrated"),
                    engine)
  for (name in c("hold_log_ratio", "hold_indices", "integrated_terms",
                 "integration_counts", "integrated_locations")) {
    local({
      field <- name
      makeActiveBinding(field, function() engines[[current]][[field]], engine)
    })
  }
  engine$form <- function() current
  engine$set_form <- function(form) {
    current <<- match.arg(form, c("sampled", "integrated"))
    invisible(current)
  }
  engine$choose_form <- function(parameters) {
    values <- sampled$unpack(parameters)
    variance <- values$sigma2 * values$nugget_ratio
    engine$set_form(if (isTRUE(variance >= threshold)) "integrated" else
      "sampled")
  }
  engine$use_proposal <- function(type) {
    for (form in engines) form$use_proposal(type)
    invisible(NULL)
  }
  engine$reset_proposal_failures <- function() {
    for (form in engines) form$reset_proposal_failures()
    invisible(NULL)
  }
  engine$mode_search_counts <- function() {
    counts <- lapply(engines, function(form) form$mode_search_counts())
    list(searches = sum(vapply(counts, `[[`, integer(1), "searches")),
         at_limit = sum(vapply(counts, `[[`, integer(1), "at_limit")),
         restarted = sum(vapply(counts, `[[`, integer(1), "restarted")))
  }
  engine$proposal_errors <- function() {
    errors <- lapply(engines, function(form) form$proposal_errors())
    firsts <- Filter(Negate(is.null), lapply(errors, `[[`, "first"))
    list(count = sum(vapply(errors, `[[`, integer(1), "count")),
         first = if (length(firsts)) firsts[[1L]])
  }
  engine
}

##' Identify variance parameters estimated as effectively zero
##'
##' `linked_index` is a list named by variance index giving parameters that
##' are not identifiable when that variance is zero.
##'
##' Variance components are optimised on the log scale, so a component whose
##' maximum likelihood estimate is zero drifts towards the lower bound along an
##' almost flat direction. Such components are held at the boundary when Monte
##' Carlo precision is assessed, instead of being treated as Monte Carlo error.
##' @noRd
adaptive_boundary_parameters <- function(estimate, variance_index,
                                         boundary_log_variance,
                                         linked_index = list()) {
  boundary <- variance_index[estimate[variance_index] <= boundary_log_variance]
  # A parameter that only acts through a zero variance (the spatial scale when
  # the spatial variance is zero) is not identifiable and is held with it.
  linked <- unlist(linked_index[as.character(boundary)], use.names = FALSE)
  sort(unique(c(boundary, linked)))
}

##' Re-optimise once from an interior value of boundary variance components
##'
##' A local optimiser that reaches a variance boundary from one start may have
##' missed an interior maximum. One restart from a moderate value (a quarter of
##' the reference scale) on the same common-random-number bank resolves this;
##' the better objective is retained, so the restart can never worsen a fit.
##' @noRd
restart_from_interior <- function(objective, fit, variance_index,
                                  boundary_log_variance, lower, upper) {
  boundary <- adaptive_boundary_parameters(
    fit$par, variance_index, boundary_log_variance
  )
  if (!length(boundary)) return(list(fit = fit, restarted = FALSE))
  interior <- fit$par
  interior[boundary] <- log(0.25)
  restart <- tryCatch(
    nlminb(interior, objective, lower = lower, upper = upper,
           control = make_optimiser_control(abs(fit$objective))),
    error = function(error) NULL
  )
  improved <- !is.null(restart) && is.finite(restart$objective) &&
    restart$objective < fit$objective - 1e-6
  list(fit = if (improved) restart else fit, restarted = TRUE,
       restart_improved = improved)
}

##' Median distance from each location to its nearest distinct neighbour
##' @noRd
median_nearest_neighbour_distance <- function(distances) {
  distances <- as.matrix(distances)
  distances[distances <= 0] <- Inf
  nearest <- apply(distances, 1L, min)
  nearest <- nearest[is.finite(nearest)]
  if (!length(nearest)) return(0)
  stats::median(nearest)
}

##' Describe variance parameters for user-facing messages
##' @noRd
describe_variance_parameters <- function(index, layout, re_names = NULL) {
  if (!length(index)) return(character())
  labels <- character()
  labels[as.character(layout$index$sigma2)] <- "spatial variance (sigma2)"
  if (!is.null(layout$index$nu2)) {
    labels[as.character(layout$index$nu2)] <- "nugget variance (tau2)"
  }
  if (!is.null(layout$index$sigma2_re)) {
    labels[as.character(layout$index$sigma2_re)] <- paste0(
      "random-effect variance for '", re_names, "'"
    )
  }
  # The scale is held only because the spatial variance is zero; it is
  # explained separately rather than listed as a variance component.
  index <- setdiff(index, layout$index$phi)
  unname(labels[as.character(index)])
}

##' Retain user-facing diagnostics of each certified MCML stage
##'
##' Drops the observed information and the per-draw workspaces, keeping the
##' estimates, budgets, the certificate and its parts, the importance-weight
##' diagnostics of each bank and the timing.
##' @noRd
compact_adaptive_mcml_history <- function(history) {
  bank_summary <- function(validation) {
    if (is.null(validation)) return(NULL)
    # A combined bank (combine_validation_banks()) keeps its record.
    validation[intersect(c("log_likelihood", "relative_ess", "maximum_weight",
                           "pareto_k", "combined"), names(validation))]
  }
  lapply(history, function(stage) {
    stage$first_validation <- bank_summary(stage$first_validation)
    stage$second_validation <- bank_summary(stage$second_validation)
    stage$replicate_validation <- bank_summary(stage$replicate_validation)
    stage$information <- NULL
    stage$information_root <- NULL
    stage
  })
}

##' Maximise the Laplace approximation to the log-likelihood
##'
##' With a single standard-normal draw at zero and the Laplace proposal, the
##' engine evaluates the Laplace approximation to the marginal log-likelihood
##' exactly. Its maximum is cheap to compute and close to the Monte Carlo
##' maximum, so Monte Carlo maximum likelihood starts from it. `objective`
##' takes the parameters and a bank.
##' @noRd
maximise_laplace <- function(objective, start, latent_dimension, lower,
                             upper) {
  zero <- matrix(0, nrow = 1L, ncol = latent_dimension)
  negative <- function(parameters) {
    value <- objective(parameters, zero)
    if (is.finite(value)) -value else 1e100
  }
  fit <- tryCatch(
    nlminb(pmin(pmax(start, lower), upper), negative, lower = lower,
           upper = upper, control = make_optimiser_control()),
    error = function(error) NULL
  )
  if (is.null(fit) || !is.finite(fit$objective) || fit$objective >= 1e99) {
    return(NULL)
  }
  list(estimate = fit$par, log_likelihood = -fit$objective,
       convergence = fit$convergence)
}

##' Dispersed starting values for the check for other maxima
##'
##' The certificate establishes Monte Carlo precision at the maximum found,
##' not that no higher maximum exists. In sparse surveys the Laplace
##' approximation can have its only maximum at a spatial scale far below the
##' site spacing while the likelihood has a higher maximum at a longer scale,
##' so an optimiser started from the Laplace estimate stays at the lower one.
##' A fitted scale below `trigger_fraction` times the median nearest-neighbour
##' distance is the warning sign, and the fit is then repeated from four
##' dispersed starts: spatial variance 0.25 or 4, crossed with a scale equal
##' to that distance or half the largest distance.
##' @noRd
make_dispersed_starts <- function(layout, distances, trigger_fraction = 0.5) {
  spacing <- median_nearest_neighbour_distance(distances)
  grid <- expand.grid(sigma2 = c(0.25, 4), phi = c(spacing, max(distances) / 2))
  list(sigma2_index = layout$index$sigma2, phi_index = layout$index$phi,
       log_sigma2 = log(grid$sigma2), log_phi = log(grid$phi),
       labels = paste0(c("small", "large"), " variance, ",
                       rep(c("short", "long"), each = 2), " scale"),
       trigger = log(trigger_fraction * spacing))
}

##' Penalised-complexity (PC) prior penalty for the spatial covariance
##'
##' Returns log pi(sigma, phi) for the spatial standard deviation sigma and the
##' scale phi of the Matern correlation in two dimensions (Fuglstad et al.,
##' 2019), expressed on the working parameters (log sigma^2, log phi) and
##' including the Jacobian of that transformation:
##' P(sigma > sigma_threshold) = sigma_probability and
##' P(phi < scale_threshold) = scale_probability. Fixed effects, nugget and
##' unstructured random-effect variances are not penalised.
##' @noRd
make_pc_penalty <- function(layout, scale_threshold,
                            sigma_threshold = 1, sigma_probability = 0.05,
                            scale_probability = 0.05) {
  rate_sigma <- -log(sigma_probability) / sigma_threshold
  rate_scale <- -log(scale_probability) * scale_threshold
  function(parameters) {
    log_sigma2 <- parameters[layout$index$sigma2]
    log_phi <- parameters[layout$index$phi]
    sigma <- exp(log_sigma2 / 2)
    phi <- exp(log_phi)
    # pi(sigma) = rate exp(-rate sigma), with d sigma / d log sigma^2 = sigma/2;
    # pi(phi) = rate phi^-2 exp(-rate / phi), with d phi / d log phi = phi.
    (log(rate_sigma) - rate_sigma * sigma + log(sigma / 2)) +
      (log(rate_scale) - log_phi - rate_scale / phi)
  }
}

##' Add a penalty to an adaptive MCML engine
##'
##' The penalty depends only on the parameters, so it is the same for every
##' draw: importance weights, effective sample sizes and Monte Carlo error are
##' unchanged, while the objective and its curvature include the penalty.
##' @noRd
penalise_engine <- function(evaluate, penalty) {
  function(parameters, draws, log_weights_only = FALSE, ...) {
    value <- evaluate(parameters, draws, log_weights_only = log_weights_only,
                      ...)
    if (log_weights_only) return(value)
    value + penalty(parameters)
  }
}

##' Starting values on the working scale of the MCML engine
##'
##' Regression coefficients, then the log spatial variance and log scale,
##' the log nugget ratio when the nugget is estimated, and the log
##' random-effect variances.
##' @noRd
make_working_start <- function(reference_pars, fix_tau2) {
  start <- c(reference_pars$beta, log(reference_pars$sigma2),
             log(reference_pars$phi))
  if (isTRUE(fix_tau2)) {
    start <- c(start, log(reference_pars$tau2 / reference_pars$sigma2))
  }
  if (length(reference_pars$sigma2_re)) {
    start <- c(start, log(reference_pars$sigma2_re))
  }
  start
}

##' Fit a non-Gaussian geostatistical model by certified MCML
##'
##' Maximises the Laplace approximation from the supplied start, then runs
##' certified Monte Carlo maximum likelihood (`fit_adaptive_mcml()`) from that
##' estimate and returns the corrected, certified estimate. Every Monte Carlo
##' stream of the fit is derived from `seed`, and the caller's random number
##' state is restored on exit.
##' @noRd
glgpm_nong_adaptive <- function(
    y, D, coords, units_m, kappa, reference_pars, cov_offset,
    ID_coords, ID_re, fix_tau2, family, invlink, control_mcml,
    seed, messages = TRUE, distance_scale = 1) {
  stopifnot(family %in% c("poisson", "binomial"))
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  fit_seed <- as.integer(seed)

  distances <- as.matrix(pairwise_distances(coords))
  engine_with <- function(integrate_noise) {
    make_adaptive_mcml_engine(
      D = D, y = y, units_m = units_m, distances = distances,
      family = family, kappa = kappa, cov_offset = cov_offset,
      ID_coords = ID_coords, fix_tau2 = fix_tau2, ID_re = ID_re,
      invlink = invlink,
      defensive_weight = control_mcml$defensive_weight %||% 0.05,
      proposal = control_mcml$proposal %||% "global",
      integrate_noise = integrate_noise
    )
  }
  # With integrate_noise = "auto", a model with an re() term that can be
  # integrated out is fitted as with TRUE, with that term and the nugget
  # integrated out wherever the rules allow, at every stage. A model whose only
  # such term is the nugget builds the sampled and the integrated form, and
  # each stage uses one of them (make_switching_mcml_engine()). The threshold
  # can be changed through control_mcml$integrate_threshold, for experiments.
  automatic <- identical(control_mcml$integrate_noise, "auto")
  engine <- engine_with(isTRUE(control_mcml$integrate_noise) || automatic)
  if (automatic && !isTRUE(engine$integrated)) {
    engine <- engine_with(FALSE)
  } else if (automatic && setequal(engine$integrated_terms, "nugget")) {
    engine <- make_switching_mcml_engine(
      engine_with(FALSE), engine,
      threshold = control_mcml$integrate_threshold %||% 0.1
    )
  }
  start <- make_working_start(reference_pars, fix_tau2)
  # Coordinates and the initial scale are rescaled before this function is
  # called. These broad bounds prevent overflow without constraining any
  # scientifically plausible fit on the internal scale.
  lower <- c(rep(-30, ncol(D)), -20, -20,
             rep(-20, length(start) - ncol(D) - 2L))
  upper <- c(rep(30, ncol(D)), 20, 10,
             rep(20, length(start) - ncol(D) - 2L))
  pc_prior <- isTRUE(control_mcml$pc_prior)
  penalised <- make_mcml_objective(engine, control_mcml, distances,
                                   distance_scale)
  objective <- penalised$objective
  penalty <- penalised$penalty

  laplace <- NULL
  start <- unname(start)
  restart_start <- NULL
  if (isTRUE(control_mcml$laplace_start %||% !engine$integrated)) {
    laplace <- maximise_laplace(
      function(parameters, draws) {
        objective(parameters, draws, proposal_type = "laplace")
      },
      start, engine$latent_dimension, lower, upper
    )
    if (!is.null(laplace)) {
      restart_start <- start
      start <- laplace$estimate
    }
  }
  # The first stage of an automatic fit takes its form from the Laplace
  # estimate, or without it from the starting values. In the integrated form
  # it starts from the starting values, as integrate_noise = TRUE does, unless
  # laplace_start asks for the Laplace estimate.
  if (is.function(engine$choose_form)) {
    engine$choose_form(start)
    if (engine$integrated && is.null(control_mcml$laplace_start) &&
        !is.null(restart_start)) {
      start <- restart_start
      restart_start <- NULL
    }
  }
  dispersed <- NULL
  if (isTRUE(control_mcml$check_maxima %||% TRUE)) {
    dispersed <- make_dispersed_starts(engine$layout, distances)
    # The screen of short scales and the exact alarms need the likelihood
    # with independent location effects, computed exactly without
    # unstructured random effects. A custom link supplies its own terms.
    if (is.null(ID_re)) {
      dispersed$screen <- make_short_scale_screen(
        engine,
        make_independent_site_likelihood(
          D, y, units_m, ID_coords %||% seq_along(y),
          cov_offset %||% numeric(length(y)), family,
          linkf = if (!identical(engine$linkf$name, "canonical")) engine$linkf
        ),
        distances, kappa
      )
    } else {
      dispersed$screen_unavailable <- "unstructured random effects"
    }
  }
  fit <- fit_adaptive_mcml(
    engine = engine, start = start, restart_start = restart_start,
    objective = objective, penalty = penalty, dispersed = dispersed,
    initial_budget = control_mcml$initial_samples,
    validation_budget = control_mcml$validation_samples,
    maximum_budget = control_mcml$max_samples,
    maximum_stages = control_mcml$max_stages,
    tolerance = control_mcml$tolerance,
    seed = fit_seed, lower = lower, upper = upper,
    variance_index = unlist(engine$layout$index[c("sigma2", "nu2",
                                                  "sigma2_re")],
                            use.names = FALSE),
    linked_index = stats::setNames(
      list(engine$layout$index$phi), engine$layout$index$sigma2
    ),
    coefficient_index = engine$layout$index$beta,
    separation = find_separation(D, y, units_m %||% rep(1, length(y)),
                                 family)
  )
  if (identical(fit$outcome, "no valid fit")) {
    stop("Monte Carlo maximum likelihood could not obtain a valid fit: ",
         fit$reason, ".", call. = FALSE)
  }
  re_names <- if (is.null(ID_re)) NULL else colnames(ID_re)
  boundary_names <- describe_variance_parameters(
    fit$boundary_parameters, engine$layout, re_names = re_names
  )
  if (!is.null(fit$held_variances)) {
    fit$held_variances$name <- describe_variance_parameters(
      fit$held_variances$parameter, engine$layout, re_names = re_names
    )
  }
  fit$proposal_errors <- engine$proposal_errors()
  small_scale <- if (!fit$outcome %in% c("numerically independent",
                                         "no maximum")) {
    detect_small_scale(fit, engine$layout, distances, kappa, distance_scale)
  }
  report_mcml_outcome(fit, engine$layout, control_mcml$tolerance, messages,
                      colnames(D), family)
  assemble_mcml_result(fit, laplace, boundary_names, D, fix_tau2, re_names,
                       family, invlink, length(y), pc_prior, fit_seed,
                       small_scale,
                       integration = list(
                         used = engine$integrated,
                         locations = engine$integrated_locations,
                         of = nrow(distances),
                         terms = engine$integrated_terms,
                         counts = engine$integration_counts,
                         forms = fit$forms
                       ))
}

##' The objective of certified MCML, with the PC penalty if requested
##'
##' Returns the engine's estimated log-likelihood, or with
##' `control_mcml$pc_prior` the penalised version (`make_pc_penalty()`), whose
##' scale threshold defaults to the median nearest-neighbour distance, and
##' the penalty itself (NULL without the prior).
##' @noRd
make_mcml_objective <- function(engine, control_mcml, distances,
                                distance_scale) {
  if (!isTRUE(control_mcml$pc_prior)) {
    return(list(objective = engine$evaluate, penalty = NULL))
  }
  scale_threshold <- if (is.null(control_mcml$pc_scale_threshold)) {
    median_nearest_neighbour_distance(distances)
  } else {
    control_mcml$pc_scale_threshold / distance_scale
  }
  penalty <- make_pc_penalty(
    engine$layout, scale_threshold = scale_threshold,
    sigma_threshold = control_mcml$pc_sigma_threshold,
    sigma_probability = control_mcml$pc_sigma_probability,
    scale_probability = control_mcml$pc_scale_probability
  )
  list(objective = penalise_engine(engine$evaluate, penalty),
       penalty = penalty)
}

##' A fitted scale smaller than the spacing between locations
##'
##' Returns the fitted scale phi, its practical range and the median
##' nearest-neighbour distance, in the model's distance units, when phi is
##' below that distance and the spatial variance is not at zero, and NULL
##' otherwise.
##' @noRd
detect_small_scale <- function(fit, layout, distances, kappa,
                                distance_scale) {
  if (layout$index$sigma2 %in% fit$boundary_parameters) return(NULL)
  spacing <- median_nearest_neighbour_distance(distances)
  scale <- exp(fit$estimate[layout$index$phi])
  if (!(scale < spacing)) return(NULL)
  practical_range <- scale * compute_practical_range_factor(kappa)
  c(scale = scale, practical_range = practical_range, spacing = spacing) *
    distance_scale
}

##' Practical range of a Matern correlation in units of its scale
##'
##' The distance, as a multiple of the scale phi, at which the Matern
##' correlation with smoothness `kappa` falls to `correlation` (about three
##' scales for the exponential correlation).
##' @noRd
compute_practical_range_factor <- function(kappa, correlation = 0.05) {
  matern <- function(scaled_distance) {
    2^(1 - kappa) / gamma(kappa) * scaled_distance^kappa *
      besselK(scaled_distance, kappa)
  }
  stats::uniroot(function(scaled_distance) {
    matern(scaled_distance) - correlation
  }, c(1e-6, 100), tol = 1e-10)$root
}

##' Build the fitted object returned by `glgpm_nong_adaptive()`
##'
##' `integration` says whether any term was integrated out (`used`) at the
##' last stage, at how many `locations` the nugget was and `of` how many, which
##' `terms` were integrated out with the `counts` of their levels, and for an
##' automatic choice of the form the `forms` of the stages.
##' @noRd
assemble_mcml_result <- function(fit, laplace, boundary_names, D, fix_tau2,
                                 re_names, family, invlink, n_observations,
                                 pc_prior, fit_seed, small_scale = NULL,
                                 integration = NULL) {
  last <- fit$history[[length(fit$history)]]
  estimate <- structure_estimate(
    fit$estimate, beta_names = colnames(D), fix_tau2 = fix_tau2,
    sigma2_me = FALSE, re_names = re_names
  )
  covariance <- fit$covariance
  if (is.null(covariance)) {
    covariance <- matrix(
      NA_real_, nrow = length(fit$estimate), ncol = length(fit$estimate)
    )
  }
  flat_names <- names(unlist(estimate))
  dimnames(covariance) <- list(flat_names, flat_names)
  out <- list(
    estimate = estimate,
    grad_MLE = fit$score,
    covariance = covariance,
    log_lik = fit$log_likelihood,
    S_samples = NULL,
    link_function = make_inverse_link_functions(family, invlink,
                                                n_observations)
  )
  attr(out, "optimizer") <- list(
    convergence = last$convergence,
    message = last$message,
    evaluations = last$evaluations,
    invalid_evaluations = NA_integer_,
    max_abs_gradient = if (all(is.finite(fit$score))) {
      max(abs(fit$score))
    } else {
      NA_real_
    },
    information_rcond = if (is.null(fit$information)) NA_real_ else
      rcond(fit$information),
    information_jitter = 0,
    mcml_certified = fit$certified,
    mcml_outcome = fit$outcome,
    mcml_certificate = fit$certificate,
    mcml_proposal = fit$proposal,
    mcml_start_replaced = fit$start_replaced,
    mcml_proposal_fallbacks = fit$proposal_fallbacks,
    mcml_proposal_errors = fit$proposal_errors,
    mcml_mode_searches = fit$mode_searches,
    mcml_integrate_noise = isTRUE(integration$used),
    mcml_noise_forms = integration$forms,
    mcml_integrated_locations = integration$locations %||% 0L,
    mcml_locations = integration$of,
    mcml_integrated_terms = integration$terms,
    mcml_integration_counts = integration$counts,
    mcml_pc_prior = pc_prior,
    laplace_start = laplace,
    mcml_boundary = boundary_names,
    mcml_held_variances = fit$held_variances,
    mcml_identifiability_failure = identical(fit$outcome,
                                             "weak identification"),
    mcml_budget = fit$budget,
    mcml_validation_budget = fit$validation_budget,
    mcml_stages = length(fit$history),
    mcml_history = compact_adaptive_mcml_history(fit$history),
    mcml_small_scale = !is.null(small_scale),
    mcml_spacing = small_scale[["spacing"]],
    mcml_maxima_check = fit$maxima_check,
    mcml_maxima_screen = fit$maxima_screen,
    mcml_maxima_alarm = fit$maxima_alarm,
    mcml_numerically_independent = identical(fit$outcome,
                                             "numerically independent"),
    mcml_no_maximum = identical(fit$outcome, "no maximum"),
    mcml_separation = if (!is.null(fit$separation)) {
      c(fit$separation[c("observations", "no_positives", "only_positives",
                         "everywhere")],
        list(coefficients = colnames(D)[fit$separation$coefficients]))
    },
    mcml_coefficients_at_bound = colnames(D)[fit$coefficients_at_bound],
    mcml_near_limit = !is.null(fit$maxima_screen$near_limit),
    mcml_better_maximum = fit$better_maximum,
    mcml_ambiguous_maxima = fit$ambiguous_maxima,
    seed = fit_seed,
    elapsed_seconds = fit$elapsed_seconds
  )
  class(out) <- "RiskMap"
  # The stages where the integrated form accepted a global iteration that
  # stopped at its limit of refits, with the largest last change. Absent when
  # there are none, as in every fit of the sampled form.
  stalled <- summarise_stalled_anchors(fit$history)
  if (!is.null(stalled)) {
    attr(out, "optimizer")$mcml_stalled_anchors <- stalled
  }
  out
}

##' Warn about what the check for other maxima found
##'
##' An exact alarm shows that the estimate is not the global maximum. A separate
##' maximum, higher or about equal, means that Wald intervals may describe the
##' uncertainty poorly. `check_fit()` gives the details and what to do.
##' @noRd
report_maxima_check <- function(fit) {
  if (isTRUE(fit$maxima_screen$alarm$persists)) {
    warning(
      "An exact lower bound on the likelihood shows that the estimate is not ",
      "the maximum likelihood estimate. See check_fit().",
      call. = FALSE
    )
  }
  if (isTRUE(fit$ambiguous_maxima)) {
    warning(
      "Separate maxima of the likelihood fit about equally well, so Wald ",
      "intervals may describe the uncertainty poorly. See check_fit().",
      call. = FALSE
    )
  } else if (isTRUE(fit$better_maximum)) {
    warning(
      "The check for other maxima found a maximum higher than the one reached ",
      "from the starting values, and the estimate is taken from it. See ",
      "check_fit().",
      call. = FALSE
    )
  }
  invisible(NULL)
}

##' Warn that the fit is numerically independent
##'
##' The spatial correlation at the fitted scale is numerically negligible, so
##' the fit is one of independent location effects and the scale is
##' unresolved. That is not evidence that the maximum lies at zero scale.
##' `record` is the record of the screen (`check_other_maxima()`).
##' @noRd
report_numerically_independent <- function(record) {
  variance <- record$independent$log_variance
  warning(
    "The fitted spatial scale is so short that the fit is numerically one of ",
    "independent location effects",
    if (!is.null(variance)) paste0(" with variance ", signif(exp(variance), 3)),
    ", and the spatial scale is unresolved. Standard errors and the Monte ",
    "Carlo certificate are not reported. See check_fit().",
    call. = FALSE
  )
}

##' Warn that the likelihood has no maximum
##'
##' Gives the cause, from the separation of the data (`fit$separation`, from
##' `find_separation()`) or else the regression coefficients that reached the
##' bound of the search (`fit$coefficients_at_bound`). `check_fit()` says what
##' can be done instead. `coefficient_names` names the regression coefficients.
##' @noRd
report_no_maximum <- function(fit, coefficient_names, family) {
  binomial <- identical(family, "binomial")
  listed <- function(index) {
    names <- paste0("'", coefficient_names[index], "'")
    if (length(names) <= 2L) return(paste(names, collapse = " and "))
    paste(paste(names[-length(names)], collapse = ", "), "and",
          names[length(names)])
  }
  returned <- paste0(
    " The estimate is where the optimiser stopped, and it is returned ",
    "without standard errors or a Monte Carlo certificate. See check_fit()."
  )
  separation <- fit$separation
  everywhere <- separation$everywhere
  if (!is.null(everywhere)) {
    none <- identical(everywhere, "no positives")
    quantity <- if (binomial) "prevalence" else "mean"
    warning(
      if (!binomial) "Every count is zero" else if (none) {
        "None of the individuals tested is positive"
      } else {
        "Every individual tested is positive"
      },
      ", so the likelihood increases as the fitted ", quantity, " tends to ",
      if (none) "zero" else "one", " and has no maximum.", returned,
      call. = FALSE
    )
  } else if (!is.null(separation)) {
    describe <- function(count, kind) {
      paste(count, if (count == 1L) "observation" else "observations", kind)
    }
    observations <- c(
      if (separation$no_positives) describe(
        separation$no_positives,
        if (binomial) "with no positives" else "with a zero count"
      ),
      if (separation$only_positives) describe(
        separation$only_positives, "where every individual is positive"
      )
    )
    single <- length(separation$coefficients) == 1L
    warning(
      "The covariates separate the data. The fit to ",
      paste(observations, collapse = " and "), " improves without limit as ",
      "the regression coefficient", if (!single) "s", " of ",
      listed(separation$coefficients), " ", if (single) "grows" else "grow",
      " in absolute value, so the likelihood has no maximum.", returned,
      call. = FALSE
    )
  } else {
    at_bound <- fit$coefficients_at_bound
    single <- length(at_bound) == 1L
    warning(
      "The estimate", if (!single) "s", " of the regression coefficient",
      if (!single) "s", " of ", listed(at_bound), " reached the bound of the ",
      "search, so the likelihood appears to have no maximum.", returned,
      call. = FALSE
    )
  }
  invisible(NULL)
}

##' Tell the user how a certified MCML fit ended
##'
##' Warns only about what changes how a result is read, in a sentence or two
##' that point to `check_fit()`: a certificate that is missing or above the
##' tolerance, a likelihood without a maximum, weak identification, other
##' maxima and a spatial variance at zero. With `messages`, a certified fit gets
##' one line. A fit whose likelihood has no maximum gets only the warning of
##' `report_no_maximum()`, since the others describe a maximum.
##' `coefficient_names` names the regression coefficients.
##' @noRd
report_mcml_outcome <- function(fit, layout, tolerance, messages,
                                coefficient_names = NULL, family = "binomial") {
  if (identical(fit$outcome, "no maximum")) {
    return(report_no_maximum(fit, coefficient_names, family))
  }
  if (layout$index$sigma2 %in% fit$boundary_parameters) {
    warning(
      "The spatial variance (sigma2) was estimated as effectively zero, so the ",
      "spatial scale cannot be estimated. See check_fit().",
      call. = FALSE
    )
  }
  report_maxima_check(fit)
  if (identical(fit$outcome, "numerically independent")) {
    report_numerically_independent(fit$maxima_screen)
  } else if (identical(fit$outcome, "weak identification")) {
    warning(
      "The observed information matrix is singular, so some parameters (often ",
      "the spatial scale) are not identifiable from these data, and the ",
      "estimate has no standard errors or Monte Carlo certificate. See ",
      "check_fit().",
      call. = FALSE
    )
  } else if (isTRUE(fit$held_variance_rising)) {
    warning(
      "The Monte Carlo error of the estimate is not certified, because the ",
      "likelihood rises as a variance held at zero leaves it, so the estimate ",
      "is not a maximum. See check_fit().",
      call. = FALSE
    )
  } else if (!fit$certified) {
    last <- fit$history[[length(fit$history)]]
    refused <- isTRUE(last$residual_refused)
    warning(
      "The Monte Carlo error of the estimate is not certified. ",
      if (isTRUE(last$correction_left_bounds)) paste0(
        "The Newton correction of the estimate would have taken a variance or ",
        "scale past the upper bound of the search, along a direction where ",
        "the likelihood is nearly flat, so it was not applied, and the ",
        "approximate 95% bound of ",
        signif(fit$certificate, 2), " statistical standard errors may not hold."
      ) else if (refused) paste0(
        "The Newton correction of the estimate did not converge on the ",
        "validation draws, so the approximate 95% bound of ",
        signif(fit$certificate, 2), " statistical standard errors may not hold."
      ) else paste0(
        "Its approximate 95% bound is ", signif(fit$certificate, 2),
        " statistical standard errors, above the tolerance of ", tolerance, "."
      ),
      " See check_fit().",
      call. = FALSE
    )
  } else if (messages) {
    message(
      "The Monte Carlo error of the estimate is within ",
      signif(fit$certificate, 2), " statistical standard errors (an ",
      "approximate 95% bound, tolerance ", tolerance, ")."
    )
  }
  invisible(NULL)
}

##' Effective draws and Pareto shape of the certification bank of a stage
##'
##' NULL for a stage that did not reach certification.
##' @noRd
summarise_certification_bank <- function(stage) {
  bank <- stage$second_validation
  if (is.null(bank)) return(NULL)
  list(effective_draws = bank$relative_ess * stage$validation_budget,
       pareto_k = bank$pareto_k)
}
