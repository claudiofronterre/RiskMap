##' Observation log-likelihood, score and curvature in the linear predictor
##'
##' Uses the same formulas as `maxim_integrand()`, including custom inverse
##' links. `w` is the negative second derivative of the log-likelihood.
##' @noRd
laplace_observation_terms <- function(eta, y, units_m, family, linkf) {
  if (family == "poisson") {
    mean <- linkf$inv(eta)
    if (any(!is.finite(mean)) || any(mean <= 0)) return(NULL)
    first <- linkf$d1(eta)
    second <- linkf$d2(eta)
    residual <- y / mean - units_m
    list(
      log_likelihood = sum(y * log(pmax(units_m * mean, .Machine$double.eps)) -
                             units_m * mean),
      score = residual * first,
      curvature = y * first^2 / mean^2 - residual * second
    )
  } else {
    probability <- linkf$inv(eta)
    if (any(!is.finite(probability)) ||
        any(probability <= 0 | probability >= 1)) return(NULL)
    first <- linkf$d1(eta)
    second <- linkf$d2(eta)
    variance <- probability * (1 - probability)
    residual <- y - units_m * probability
    list(
      log_likelihood = sum(
        y * log(pmax(probability, .Machine$double.eps)) +
          (units_m - y) * log(pmax(1 - probability, .Machine$double.eps))
      ),
      score = residual * first / variance,
      curvature = units_m * first^2 / variance -
        residual * (second / variance - first^2 * (1 - 2 * probability) /
                      variance^2)
    )
  }
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
  list(columns = columns, dimension = dimension,
       pair_linear = linear,
       pair_observation = rep(seq_len(nrow(columns)), times = width^2))
}

latent_linear_predictor <- function(latent_map, latent) {
  columns <- latent_map$columns
  eta <- latent[columns[, 1L]]
  if (ncol(columns) > 1L) {
    for (j in 2:ncol(columns)) eta <- eta + latent[columns[, j]]
  }
  eta
}

latent_gradient <- function(latent_map, score) {
  out <- numeric(latent_map$dimension)
  for (j in seq_len(ncol(latent_map$columns))) {
    sums <- rowsum(score, latent_map$columns[, j], reorder = FALSE)
    index <- as.integer(rownames(sums))
    out[index] <- out[index] + sums[, 1L]
  }
  out
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
##' very few iterations. The objective is log-concave for canonical links;
##' step halving keeps every iteration an ascent step for other links.
##' @noRd
adaptive_laplace_mode <- function(mu, y, units_m, latent_map, precision,
                                  family, linkf, start,
                                  tolerance = 1e-10, maximum_iterations = 100L) {
  objective <- function(latent, terms) {
    terms$log_likelihood - 0.5 * sum(latent * as.numeric(precision %*% latent))
  }
  terms_at <- function(latent) {
    laplace_observation_terms(mu + latent_linear_predictor(latent_map, latent),
                              y, units_m, family, linkf)
  }
  latent <- start
  terms <- terms_at(latent)
  if (is.null(terms)) {
    latent <- numeric(length(start))
    terms <- terms_at(latent)
    if (is.null(terms)) stop("Invalid starting linear predictor.", call. = FALSE)
  }
  current <- objective(latent, terms)
  for (iteration in seq_len(maximum_iterations)) {
    gradient <- latent_gradient(latent_map, terms$score) -
      as.numeric(precision %*% latent)
    hessian <- precision + latent_data_hessian(latent_map, terms$curvature)
    root <- chol(hessian)
    step <- backsolve(root, forwardsolve(t(root), gradient))
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
    if (!accepted) break
    change <- max(abs(candidate - latent))
    latent <- candidate
    terms <- candidate_terms
    current <- value
    if (change <= tolerance * max(1, max(abs(latent)))) break
  }
  hessian <- precision + latent_data_hessian(latent_map, terms$curvature)
  list(mode = latent, hessian = (hessian + t(hessian)) / 2,
       iterations = iteration)
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
##' With `type = "global"` the global approximation (`fit_global_proposal()`)
##' starts from `anchor`, the converged factors at another parameter value,
##' if one is given. The factors of the likelihood in the latent linear
##' predictor s carry over exactly to a new regression part `mu`, as
##' b - C (mu - mu_anchor) and C, so the iteration starts close to its fixed
##' point and the conditional mode is not needed. Without an anchor, or if
##' that start does not converge, the iteration starts from the conditional
##' mode found by `find_mode()`, and if it fails there the Laplace
##' approximation at the mode is used. Returns the centre and precision, the
##' upper factor of the proposal covariance when the global fit computed it
##' (otherwise NULL), the global factors (NULL for the Laplace approximation)
##' and whether the global fit `failed`.
##' @noRd
build_gaussian_proposal <- function(type, mu, y, units_m, latent_map,
                                    precision, family, linkf, canonical,
                                    anchor, find_mode) {
  global_from <- function(mode = NULL, initial = NULL) {
    fit_global_proposal(
      mode = mode, mu = mu, y = y, units_m = units_m,
      latent_map = latent_map, precision = precision, family = family,
      linkf = linkf, canonical = canonical, initial = initial
    )
  }
  approximation <- NULL
  global <- NULL
  failed <- FALSE
  if (identical(type, "global")) {
    if (!is.null(anchor)) {
      global <- global_from(initial = list(
        linear = anchor$linear - anchor$curvature * (mu - anchor$mu),
        curvature = anchor$curvature
      ))
      if (!is.null(global) && !global$converged) global <- NULL
    }
    if (is.null(global)) {
      approximation <- find_mode()
      global <- global_from(mode = approximation$mode)
      failed <- is.null(global) || !global$converged
    }
  }
  if (is.null(global)) {
    if (is.null(approximation)) approximation <- find_mode()
    return(list(centre = approximation$mode,
                precision = approximation$hessian, root = NULL,
                factors = NULL, failed = failed))
  }
  list(centre = global$centre, precision = global$precision,
       root = global$root, factors = global$factors, failed = failed)
}

##' Log-likelihood of every row of a matrix of linear predictors
##'
##' `eta` has one row per draw and one column per observation. A draw whose
##' mean is outside the support has log-likelihood -Inf.
##' @noRd
compute_observation_log_likelihood <- function(eta, y, units_m, family,
                                               linkf) {
  if (identical(linkf$name, "canonical")) {
    return(cpp_canonical_log_likelihood(eta, y, units_m,
                                        poisson = family == "poisson"))
  }
  mean <- matrix(linkf$inv(as.vector(eta)), nrow = nrow(eta))
  if (family == "poisson") {
    if (any(!is.finite(mean)) || any(mean <= 0))
      return(rep(-Inf, nrow(eta)))
    return(rowSums(sweep(log(mean), 2L, y, "*") -
                     sweep(mean, 2L, units_m, "*")))
  }
  if (any(!is.finite(mean)) || any(mean <= 0 | mean >= 1))
    return(rep(-Inf, nrow(eta)))
  rowSums(sweep(log(mean), 2L, y, "*") +
            sweep(log1p(-mean), 2L, units_m - y, "*"))
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
    eta <- eta + samples[, columns[, j + 1L], drop = FALSE]
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
##' the Laplace proposal throughout. With a Gaussian proposal alone the
##' importance weights have infinite variance whenever the largest eigenvalue of
##' Sigma^{1/2} W Sigma^{1/2} exceeds one (W the likelihood curvature), which
##' holds in practically every spatial dataset (Pitt et al., 2013). The mixture
##' bounds every moment of the weights, so that the Monte Carlo error estimates
##' used for stopping are valid, at a small cost in efficiency. Weights use the
##' full mixture density. With `defensive_weight = 0` the proposal is the
##' Gaussian component alone.
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
    defensive_weight = 0.05, proposal = c("global", "laplace")) {
  family <- match.arg(family)
  default_proposal_type <- match.arg(proposal)
  proposal_failures <- 0L
  n <- length(y)
  n_loc <- nrow(distances)
  p <- ncol(D)
  if (is.null(cov_offset)) cov_offset <- numeric(n)
  if (is.null(ID_coords)) ID_coords <- seq_len(n)
  if (!is.null(ID_re)) ID_re <- as.matrix(ID_re)
  n_re <- if (is.null(ID_re)) 0L else ncol(ID_re)
  n_dim_re <- if (n_re) apply(ID_re, 2L, max) else integer()
  layout <- cov_par_layout(p, fix_tau2 = fix_tau2,
                           sigma2_me = FALSE, n_re = n_re)

  linkf <- make_inverse_link_functions(family, invlink, n)
  canonical <- identical(linkf$name, "canonical")

  unpack <- function(parameters) {
    sigma2 <- exp(parameters[layout$index$sigma2])
    nugget <- nugget_ratio(
      fix_tau2, sigma2,
      if (!is.null(layout$index$nu2)) parameters[layout$index$nu2] else NULL
    )
    list(
      beta = parameters[layout$index$beta], sigma2 = sigma2,
      phi = exp(parameters[layout$index$phi]), nugget_ratio = nugget,
      sigma2_re = if (n_re) exp(parameters[layout$index$sigma2_re]) else NULL
    )
  }

  half_integer_kappa <- kappa %in% c(0.5, 1.5, 2.5)
  distance_values <- as.numeric(distances)
  covariance_at <- function(values) {
    # Evaluate the closed-form half-integer Matern correlation directly on the
    # stored distance matrix rather than rebuilding it from a dist object.
    correlation <- if (half_integer_kappa) {
      matrix(cpp_half_integer_matern(distance_values, values$phi, kappa),
             n_loc, n_loc)
    } else {
      matern_correlation(stats::as.dist(distances), phi = values$phi,
                         kappa = kappa, return_sym_matrix = TRUE)
    }
    diag(correlation) <- diag(correlation) + values$nugget_ratio
    values$sigma2 * correlation
  }

  latent_layout <- function() {
    indices <- vector("list", n_re)
    next_index <- n_loc + 1L
    if (n_re) for (j in seq_len(n_re)) {
      indices[[j]] <- next_index:(next_index + n_dim_re[j] - 1L)
      next_index <- next_index + n_dim_re[j]
    }
    list(spatial = seq_len(n_loc), random = indices,
         dimension = next_index - 1L)
  }
  latent <- latent_layout()
  random_first <- vapply(latent$random, function(index) index[1L], integer(1))
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
    key <- c(values$phi, values$sigma2, values$nugget_ratio)
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
        values$nugget_ratio
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
    cbind(ID_coords, if (n_re) vapply(seq_len(n_re), function(j) {
      latent$random[[j]][ID_re[, j]]
    }, integer(n))),
    latent$dimension
  )
  # The previous conditional mode is an excellent Newton start for nearby
  # parameter values; the most recent factorisations are reused when the
  # proposal and target parameters coincide.
  last_mode <- numeric(latent$dimension)
  # Converged global factors at a fixed parameter value (set by
  # set_anchor()), used to start the global iteration at other values.
  anchor <- NULL
  cached_parameters <- NULL
  cached_type <- NULL
  cached_proposal <- NULL

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
    find_mode <- function() {
      approximation <- adaptive_laplace_mode(
        mu = mu, y = y, units_m = units_m, latent_map = latent_map,
        precision = precision, family = family, linkf = linkf,
        start = last_mode
      )
      last_mode <<- approximation$mode
      approximation
    }
    gaussian <- build_gaussian_proposal(
      type, mu, y, units_m, latent_map, precision, family, linkf, canonical,
      anchor, find_mode
    )
    if (gaussian$failed) proposal_failures <<- proposal_failures + 1L
    root <- gaussian$root
    if (is.null(root)) {
      root <- tryCatch(
        inverse_cholesky_factor(gaussian$precision),
        error = function(error) factor_covariance(
          solve(gaussian$precision), "adaptive MCML proposal"
        )
      )
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
               error = function(error) NULL)
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
    log_weights <- if (canonical) {
      cpp_importance_log_weights(
        draws, proposal$root, proposal$centre, defensive, defensive_share,
        proposal$spatial_root, as.numeric(proposal$values$sigma2_re),
        spatial_root, as.numeric(values$sigma2_re), same_parameters,
        linear_offset, latent_map$columns, random_first, random_size, y,
        units_m, family == "poisson"
      )
    } else {
      compute_log_weights_in_r(
        draws, proposal, defensive, defensive_share, values, spatial_root,
        same_parameters, linear_offset, latent, latent_map$columns,
        function(eta) {
          compute_observation_log_likelihood(eta, y, units_m, family, linkf)
        }
      )
    }
    if (log_weights_only) return(log_weights)
    maximum <- max(log_weights)
    maximum + log(mean(exp(log_weights - maximum)))
  }

  list(evaluate = evaluate, proposal_at = proposal_at, layout = layout,
       latent_dimension = latent$dimension, unpack = unpack,
       spatial_covariance = function(parameters) {
         covariance_at(unpack(parameters))
       },
       defensive_weight = defensive_weight,
       proposal_type = function() default_proposal_type,
       use_proposal = function(type) {
         default_proposal_type <<- match.arg(type, c("global", "laplace"))
         invisible(NULL)
       },
       proposal_failures = function() proposal_failures,
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
         invisible(NULL)
       })
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
    validation[c("log_likelihood", "relative_ess", "maximum_weight",
                 "pareto_k")]
  }
  lapply(history, function(stage) {
    stage$first_validation <- bank_summary(stage$first_validation)
    stage$second_validation <- bank_summary(stage$second_validation)
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
##' to that distance or half the largest distance. With `messages = TRUE` the
##' user is told when and why the refit runs, with distances in the model's
##' units. The scale quoted is that of the fit on the fitting draws, before
##' the Newton correction that gives the returned estimate.
##' @noRd
make_dispersed_starts <- function(layout, distances, distance_scale = 1,
                             messages = FALSE, trigger_fraction = 0.5) {
  spacing <- median_nearest_neighbour_distance(distances)
  grid <- expand.grid(sigma2 = c(0.25, 4), phi = c(spacing, max(distances) / 2))
  announce <- function(log_phi) {
    if (!messages) return(invisible(NULL))
    message(
      "The spatial scale of the first fit (phi = ",
      signif(exp(log_phi) * distance_scale, 3), ", before its Monte Carlo ",
      "correction) is below half the median ",
      "distance between neighbouring locations (",
      signif(spacing * distance_scale, 3), "). In sparse data a maximum at so ",
      "short a scale can be a local one, with a higher maximum at a longer ",
      "scale, so the fit is repeated from four dispersed starting values to ",
      "check."
    )
  }
  list(sigma2_index = layout$index$sigma2, phi_index = layout$index$phi,
       log_sigma2 = log(grid$sigma2), log_phi = log(grid$phi),
       labels = paste0(c("small", "large"), " variance, ",
                       rep(c("short", "long"), each = 2), " scale"),
       trigger = log(trigger_fraction * spacing), announce = announce)
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
  engine <- make_adaptive_mcml_engine(
    D = D, y = y, units_m = units_m, distances = distances,
    family = family, kappa = kappa, cov_offset = cov_offset,
    ID_coords = ID_coords, fix_tau2 = fix_tau2, ID_re = ID_re,
    invlink = invlink,
    defensive_weight = control_mcml$defensive_weight %||% 0.05,
    proposal = control_mcml$proposal %||% "global"
  )
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
  if (isTRUE(control_mcml$laplace_start)) {
    laplace <- maximise_laplace(
      function(parameters, draws) {
        objective(parameters, draws, proposal_type = "laplace")
      },
      start, engine$latent_dimension, lower, upper
    )
    if (!is.null(laplace)) start <- laplace$estimate
  }
  dispersed <- NULL
  if (isTRUE(control_mcml$check_maxima %||% TRUE)) {
    dispersed <- make_dispersed_starts(engine$layout, distances,
                                       distance_scale, messages)
    # The screen of short scales and the exact alarms need the likelihood
    # with independent location effects, computed exactly for the canonical
    # links and without unstructured random effects.
    if (is.null(ID_re) && is.null(invlink)) {
      dispersed$screen <- make_short_scale_screen(
        engine,
        make_independent_site_likelihood(
          D, y, units_m, ID_coords %||% seq_along(y),
          cov_offset %||% numeric(length(y)), family
        ),
        distances, kappa, messages, distance_scale
      )
    } else {
      dispersed$screen_unavailable <- if (is.null(ID_re)) {
        "custom inverse link"
      } else {
        "unstructured random effects"
      }
    }
  }
  fit <- fit_adaptive_mcml(
    engine = engine, start = start, objective = objective,
    penalty = penalty, dispersed = dispersed,
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
    )
  )
  if (identical(fit$outcome, "no valid fit")) {
    stop("Monte Carlo maximum likelihood could not obtain a valid fit: ",
         fit$reason, ".", call. = FALSE)
  }
  re_names <- if (is.null(ID_re)) NULL else colnames(ID_re)
  boundary_names <- describe_variance_parameters(
    fit$boundary_parameters, engine$layout, re_names = re_names
  )
  small_scale <- if (!identical(fit$outcome, "numerically independent")) {
    detect_small_scale(fit, engine$layout, distances, kappa, distance_scale)
  }
  report_mcml_outcome(fit, boundary_names, engine$layout, pc_prior,
                      control_mcml$tolerance, messages, small_scale,
                      distance_scale)
  assemble_mcml_result(fit, laplace, boundary_names, D, fix_tau2, re_names,
                       family, invlink, length(y), pc_prior, fit_seed,
                       small_scale)
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
##' @noRd
assemble_mcml_result <- function(fit, laplace, boundary_names, D, fix_tau2,
                                 re_names, family, invlink, n_observations,
                                 pc_prior, fit_seed, small_scale = NULL) {
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
    mcml_pc_prior = pc_prior,
    laplace_start = laplace,
    mcml_boundary = boundary_names,
    mcml_identifiability_failure = identical(fit$outcome,
                                             "weak identification"),
    mcml_budget = fit$budget,
    mcml_validation_budget = fit$validation_budget,
    mcml_stages = length(fit$history),
    mcml_history = compact_adaptive_mcml_history(fit$history),
    mcml_small_scale = !is.null(small_scale),
    mcml_maxima_check = fit$maxima_check,
    mcml_maxima_screen = fit$maxima_screen,
    mcml_maxima_alarm = fit$maxima_alarm,
    mcml_numerically_independent = identical(fit$outcome,
                                             "numerically independent"),
    mcml_near_limit = !is.null(fit$maxima_screen$near_limit),
    mcml_better_maximum = fit$better_maximum,
    mcml_ambiguous_maxima = fit$ambiguous_maxima,
    seed = fit_seed,
    elapsed_seconds = fit$elapsed_seconds
  )
  class(out) <- "RiskMap"
  out
}

##' Tell the user what the check for other maxima found
##'
##' Warns when an exact alarm shows that the estimate is not the global
##' maximum, when the estimate comes from a maximum separate from the one
##' first fitted, and when two separate maxima fit about equally well on
##' independent draws. Neither of the last two proves that no higher maximum
##' exists elsewhere. The scale quoted for the kept maximum is that of the
##' returned estimate, after the Newton correction, and the other maximum's
##' is that of its fit on the fitting draws.
##' @noRd
report_maxima_check <- function(fit, layout, distance_scale) {
  alarm <- fit$maxima_screen$alarm
  if (isTRUE(alarm$persists)) {
    warning(
      "An exact lower bound on the likelihood at a short spatial scale, from ",
      "the limit in which the model has independent location effects, is ",
      "above the log-likelihood of the estimate. The estimate is therefore ",
      "not the maximum likelihood estimate, and the check for other maxima ",
      "did not reach a higher maximum. Inspect the profile likelihood of the ",
      "spatial scale.",
      call. = FALSE
    )
  }
  check <- fit$maxima_check
  if (is.null(check) || !(fit$better_maximum || fit$ambiguous_maxima)) {
    return(invisible(NULL))
  }
  judged <- check[is.finite(check$independent_log_likelihood), ]
  kept <- judged[judged$chosen, ]
  rest <- judged[!judged$chosen, ]
  # The maximum compared with the kept one: the fitted maximum when the
  # estimate moved, otherwise the closest rival on the independent draws.
  other <- if (fit$ambiguous_maxima) {
    rest[which.max(rest$independent_log_likelihood), ]
  } else {
    rest[rest$maximum == check$maximum[[1L]], ][1L, ]
  }
  scale <- function(row) signif(exp(row$log_phi) * distance_scale, 3)
  estimated_scale <- signif(exp(fit$estimate[[layout$index$phi]]) *
                              distance_scale, 3)
  difference <- paste0(
    signif(kept$independent_log_likelihood - other$independent_log_likelihood,
           2),
    " log-likelihood units (standard error ", signif(other$independent_se, 2),
    ") on independent Monte Carlo draws"
  )
  profile_advice <- paste0(
    " The likelihood has more than one maximum, so Wald intervals may ",
    "describe the uncertainty poorly. Consider inspecting the profile ",
    "likelihood of the spatial scale and variance. The check does not prove ",
    "that no higher maximum exists elsewhere."
  )
  if (fit$ambiguous_maxima) {
    warning(
      "The check for other maxima found separate maxima of the likelihood ",
      "that fit about equally well. They differ by ", difference, ". The ",
      "estimate is taken from the higher one, with spatial scale ",
      estimated_scale, ". The other has scale ", scale(other), ".",
      profile_advice,
      call. = FALSE
    )
  } else {
    warning(
      "The check for other maxima found a separate maximum of the ",
      "likelihood, higher than the one reached from the starting values by ",
      difference, ", and the estimate is taken from it (spatial scale ",
      estimated_scale, " instead of ", scale(other), ").",
      profile_advice,
      call. = FALSE
    )
  }
  invisible(NULL)
}

##' Tell the user that the fit is numerically independent
##'
##' The spatial correlation at the fitted scale is numerically negligible, so
##' the fit is one of independent location effects and the scale is
##' unresolved. That is not evidence that the maximum lies at zero scale.
##' Points to the profile likelihood of the scale, to a model without spatial
##' correlation and to the PC penalty. `record` is the record of the screen
##' (`check_other_maxima()`).
##' @noRd
report_numerically_independent <- function(record, pc_advice = "") {
  variance <- record$independent$log_variance
  warning(
    "The fitted spatial scale is so short that the correlation between the ",
    "closest locations is below 1e-6, so the fit is numerically one of ",
    "independent location effects",
    if (!is.null(variance)) paste0(" with variance ", signif(exp(variance), 3)),
    ", and the spatial scale is unresolved. Standard errors and the Monte ",
    "Carlo certificate, which need an interior maximum, are not reported. ",
    "Use profile(fit, parameters = \"phi\") to see how the likelihood changes ",
    "with the scale. If spatial correlation is not expected, consider a model ",
    "with an unstructured effect for each location instead of gp().",
    pc_advice,
    call. = FALSE
  )
}

##' Tell the user that a short fitted scale is close to the limit of zero
##'
##' The comparison with the limit has not shown that the fitted likelihood is
##' above it at the current Monte Carlo precision. That is not evidence that
##' the maximum is at zero, so the estimate is assessed as usual and only this
##' warning is given. `scale` is the spatial scale of the returned estimate,
##' in the model's distance units, which the warning quotes rather than the
##' scale of the fit on the fitting draws compared with the limit.
##' @noRd
report_near_limit <- function(record, scale) {
  near <- record$near_limit
  if (is.null(near)) return(invisible(NULL))
  warning(
    "The estimated spatial scale (phi = ", signif(scale, 3),
    ") is short, and the fitted likelihood is not clearly distinguishable from ",
    "the limit with independent location effects at the current Monte Carlo ",
    "precision (difference ", signif(near$gap, 2), ", standard error ",
    signif(near$se, 2), " on ", format(near$draws, big.mark = ","),
    " independent draws). Use profile(fit, parameters = \"phi\") to see how ",
    "well the data determine the scale.",
    call. = FALSE
  )
}

##' Tell the user how a certified MCML fit ended
##'
##' Warnings name variance components at zero, a scale smaller than the
##' location spacing, a singular information matrix, or a certificate above
##' the tolerance. Weak identification of the spatial variance or scale is what
##' the PC penalty addresses, so these warnings point to it when it is off.
##' @noRd
report_mcml_outcome <- function(fit, boundary_names, layout, pc_prior,
                                tolerance, messages, small_scale = NULL,
                                distance_scale = 1) {
  pc_advice <- if (pc_prior) "" else paste0(
    " Alternatively, refit with 'control_mcml = set_control_mcml(pc_prior = ",
    "TRUE)', which stabilises the spatial variance and scale with a mild ",
    "penalty. This is recommended for small or low-prevalence surveys and ",
    "makes prediction with 'parameter_uncertainty = TRUE' possible."
  )
  spatial_at_zero <- layout$index$sigma2 %in% fit$boundary_parameters
  if (length(boundary_names)) {
    single <- length(boundary_names) == 1L
    warning(
      "The ", if (single) boundary_names else paste(
        paste(boundary_names[-length(boundary_names)], collapse = ", "),
        "and the", boundary_names[length(boundary_names)]
      ), " ", if (single) "was" else "were",
      " estimated as effectively zero. ",
      if (single) "This component" else "These components",
      " cannot be distinguished from zero with these data, so ",
      if (single) "its standard error is" else "their standard errors are",
      " not reported",
      if (spatial_at_zero) ", and the spatial scale cannot be estimated",
      ". Consider refitting without ", if (single) "it." else "them.",
      if (spatial_at_zero) pc_advice,
      call. = FALSE
    )
  }
  if (!is.null(small_scale)) {
    warning(
      "The estimated spatial scale (phi = ", signif(small_scale[["scale"]], 3),
      ", practical range ", signif(small_scale[["practical_range"]], 3),
      ") is smaller than the median distance between neighbouring locations (",
      signif(small_scale[["spacing"]], 3), "), all in the distance units of ",
      "the model. The practical range is the distance at which the ",
      "correlation falls to about 0.05. Spatial correlation on so short a ",
      "scale is hard to distinguish from unstructured variation between ",
      "locations, so the scale may be weakly identified by these data and its ",
      "standard error may be large. Its profile likelihood shows how much the ",
      "data say about it. The estimate is the maximum likelihood estimate. ",
      "Consider a model with a nugget (unstructured) effect, ",
      "'gp(nugget = TRUE)', or fixing the scale from external knowledge.",
      pc_advice,
      call. = FALSE
    )
  }
  report_maxima_check(fit, layout, distance_scale)
  report_near_limit(fit$maxima_screen,
                    exp(fit$estimate[[layout$index$phi]]) * distance_scale)
  draws <- paste0(format(fit$budget, big.mark = ","), " fitting and ",
                  format(fit$validation_budget, big.mark = ","),
                  " validation draws")
  if (identical(fit$outcome, "numerically independent")) {
    report_numerically_independent(fit$maxima_screen, pc_advice)
  } else if (identical(fit$outcome, "weak identification")) {
    warning(
      "The observed information matrix is singular, so some parameters ",
      "(often the spatial scale when spatial correlation is weak) are not ",
      "identifiable from these data. The estimate is returned without ",
      "standard errors and its Monte Carlo precision is not established. ",
      "Consider a simpler model or fixing the scale.",
      pc_advice,
      call. = FALSE
    )
  } else if (!fit$certified) {
    refused <- isTRUE(fit$history[[length(fit$history)]]$residual_refused)
    warning(
      "The Monte Carlo precision of the estimate could not be established ",
      "within the permitted stages and draws (", draws, " at the last ",
      "stage). ",
      if (refused) paste0(
        "At the last stage the Newton correction of the estimate did not ",
        "converge on the validation draws, so the approximate 95% bound of ",
        signif(fit$certificate, 2), " statistical standard errors on its ",
        "Monte Carlo error may not hold. Increase 'max_samples' in ",
        "set_control_mcml()."
      ) else paste0(
        "The approximate 95% bound on its Monte Carlo error is ",
        signif(fit$certificate, 2), " statistical standard errors, above the ",
        "tolerance of ", tolerance, ". Increase 'max_samples' or 'max_stages' ",
        "in set_control_mcml(), or relax 'tolerance'."
      ),
      call. = FALSE
    )
  } else if (messages) {
    message(
      "The Monte Carlo error of the estimate is within ",
      signif(fit$certificate, 2), " statistical standard errors (an ",
      "approximate 95% bound, tolerance ", tolerance, "), after ",
      length(fit$history), " stage", if (length(fit$history) == 1L) "" else
        "s", " with ", draws, ". This concerns the maximum found, not ",
      "whether it is the highest."
    )
  }
  invisible(NULL)
}
