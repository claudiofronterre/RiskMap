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
##' results to be reproduced.
##' @noRd
inverse_cholesky_factor <- function(hessian) {
  dimension <- nrow(hessian)
  reversed <- dimension:1
  root <- chol(hessian[reversed, reversed, drop = FALSE])
  inverse_root <- backsolve(root, diag(dimension))
  t(inverse_root)[reversed, reversed, drop = FALSE]
}

##' Construct the adaptive local-Laplace importance engine
##'
##' The standard-normal bank is held fixed within a stage, while the local
##' Gaussian map changes smoothly with the trial parameters. This preserves a
##' deterministic common-random-number objective without tying the proposal to
##' the user's starting values.
##' @noRd
make_adaptive_mcml_engine <- function(
    D, y, units_m, distances, family = c("binomial", "poisson"),
    kappa = 0.5, cov_offset = NULL, ID_coords = NULL,
    fix_tau2 = FALSE, ID_re = NULL, invlink = NULL) {
  family <- match.arg(family)
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

  linkf <- make_inverse_link_functions(family, invlink)
  inverse_link <- linkf$inv
  inverse_derivative <- linkf$d1
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

  log_observation_likelihood <- function(eta) {
    if (canonical) {
      return(cpp_canonical_log_likelihood(eta, y, units_m,
                                          poisson = family == "poisson"))
    }
    mean <- matrix(inverse_link(as.vector(eta)), nrow = nrow(eta))
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

  observation_score <- function(eta) {
    if (canonical && family == "binomial") {
      probability <- plogis(eta)
      expected <- sweep(probability, 2L, units_m, "*")
      return(sweep(expected, 2L, y,
                   function(x, observed) observed - x))
    }
    if (canonical) {
      mean <- exp(eta)
      expected <- sweep(mean, 2L, units_m, "*")
      return(sweep(expected, 2L, y,
                   function(x, observed) observed - x))
    }
    mean <- matrix(inverse_link(as.vector(eta)), nrow = nrow(eta))
    first <- matrix(inverse_derivative(as.vector(eta)), nrow = nrow(eta))
    if (family == "poisson") {
      expected <- sweep(mean, 2L, units_m, "*")
      residual <- sweep(expected, 2L, y,
                        function(x, observed) observed - x)
      return(residual * first / mean)
    }
    expected <- sweep(mean, 2L, units_m, "*")
    residual <- sweep(expected, 2L, y,
                      function(x, observed) observed - x)
    residual * first / (mean * (1 - mean))
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
  cached_parameters <- NULL
  cached_proposal <- NULL

  proposal_at <- function(parameters) {
    if (!is.null(cached_parameters) && identical(parameters, cached_parameters)) {
      return(cached_proposal)
    }
    values <- unpack(parameters)
    covariance <- covariance_at(values)
    spatial_root <- factor_covariance(covariance, "adaptive MCML target")
    precision <- matrix(0, latent$dimension, latent$dimension)
    precision[latent$spatial, latent$spatial] <- chol2inv(spatial_root)
    if (n_re) for (j in seq_len(n_re)) {
      index <- latent$random[[j]]
      precision[cbind(index, index)] <- 1 / values$sigma2_re[j]
    }
    mu <- as.numeric(D %*% values$beta + cov_offset)
    approximation <- adaptive_laplace_mode(
      mu = mu, y = y, units_m = units_m, latent_map = latent_map,
      precision = precision, family = family, linkf = linkf,
      start = last_mode
    )
    last_mode <<- approximation$mode
    root <- tryCatch(
      inverse_cholesky_factor(approximation$hessian),
      error = function(error) factor_covariance(
        solve(approximation$hessian), "adaptive MCML proposal"
      )
    )
    proposal <- list(mode = approximation$mode, root = root, values = values,
                     covariance = covariance, spatial_root = spatial_root)
    cached_parameters <<- parameters
    cached_proposal <<- proposal
    proposal
  }

  evaluate <- function(parameters, z, diagnostics = FALSE,
                       log_weights_only = FALSE,
                       proposal_parameters = parameters) {
    # Keeping proposal_parameters fixed lets the calibration experiments
    # separate ordinary importance-sampling error from variability caused by
    # rebuilding the local Laplace proposal as the optimizer moves.
    proposal <- tryCatch(
      proposal_at(proposal_parameters), error = function(error) NULL
    )
    if (is.null(proposal)) return(if (diagnostics) NULL else -1e100)
    samples <- z %*% proposal$root
    samples <- samples + rep(proposal$mode, each = nrow(samples))
    if (identical(parameters, proposal_parameters)) {
      values <- proposal$values
      covariance <- proposal$covariance
      spatial_root <- proposal$spatial_root
    } else {
      values <- unpack(parameters)
      covariance <- covariance_at(values)
      spatial_root <- tryCatch(
        factor_covariance(covariance, "adaptive MCML target"),
        error = function(error) NULL
      )
    }
    if (is.null(spatial_root)) return(if (diagnostics) NULL else -1e100)

    spatial_samples <- samples[, latent$spatial, drop = FALSE]
    solved <- solve_from_cholesky(spatial_root, t(spatial_samples))
    log_prior <- -0.5 * (
      n_loc * log(2 * pi) + 2 * sum(log(diag(spatial_root))) +
        colSums(t(spatial_samples) * solved)
    )
    eta <- spatial_samples[, ID_coords, drop = FALSE] +
      rep(as.numeric(D %*% values$beta + cov_offset), each = nrow(samples))
    if (n_re) for (j in seq_len(n_re)) {
      random_samples <- samples[, latent$random[[j]], drop = FALSE]
      eta <- eta + random_samples[, ID_re[, j], drop = FALSE]
      log_prior <- log_prior - 0.5 * (
        n_dim_re[j] * log(2 * pi * values$sigma2_re[j]) +
          rowSums(random_samples^2) / values$sigma2_re[j]
      )
    }
    log_target <- log_prior + log_observation_likelihood(eta)

    # The draws are mode + z %*% root, so the proposal's standardised
    # residuals are exactly z and no triangular solve is needed.
    log_proposal <- -0.5 * (
      latent$dimension * log(2 * pi) + 2 * sum(log(diag(proposal$root))) +
        rowSums(z^2)
    )
    log_weights <- log_target - log_proposal
    if (log_weights_only) return(log_weights)
    maximum <- max(log_weights)
    shifted <- exp(log_weights - maximum)
    weights <- shifted / sum(shifted)
    value <- maximum + log(mean(shifted))
    if (!diagnostics) return(value)
    psis_fit <- suppressWarnings(loo::psis(log_weights))
    observation_scores <- observation_score(eta) %*% D
    inverse_covariance <- chol2inv(spatial_root)
    base_correlation <- matern_correlation(
      stats::as.dist(distances), phi = values$phi, kappa = kappa,
      return_sym_matrix = TRUE
    )
    covariance_derivatives <- list()
    covariance_derivatives[[1L]] <- if (isTRUE(fix_tau2)) {
      covariance
    } else {
      values$sigma2 * base_correlation
    }
    covariance_derivatives[[2L]] <- values$sigma2 * values$phi *
      matern_gradient_phi(stats::as.dist(distances), values$phi, kappa)
    if (isTRUE(fix_tau2)) {
      covariance_derivatives[[3L]] <-
        diag(values$sigma2 * values$nugget_ratio, n_loc)
    }
    covariance_scores <- vapply(
      covariance_derivatives,
      function(derivative) {
        quadratic_matrix <- inverse_covariance %*% derivative %*%
          inverse_covariance
        -0.5 * sum(diag(inverse_covariance %*% derivative)) +
          0.5 * rowSums(
            (spatial_samples %*% quadratic_matrix) * spatial_samples
          )
      },
      numeric(nrow(samples))
    )
    if (is.null(dim(covariance_scores))) {
      covariance_scores <- matrix(covariance_scores, ncol = 1L)
    }
    sample_scores <- cbind(observation_scores, covariance_scores)
    if (n_re) for (j in seq_len(n_re)) {
      random_samples <- samples[, latent$random[[j]], drop = FALSE]
      sample_scores <- cbind(
        sample_scores,
        -0.5 * n_dim_re[j] +
          0.5 * rowSums(random_samples^2) / values$sigma2_re[j]
      )
    }
    mean_score <- colSums(sweep(sample_scores, 1L, weights, "*"))
    centered_scores <- sweep(sample_scores, 2L, mean_score, "-")
    weighted_scores <- sweep(centered_scores, 1L, weights, "*")
    list(
      value = value,
      relative_ess = 1 / sum(weights^2) / length(weights),
      maximum_weight = max(weights),
      pareto_k = as.numeric(loo::pareto_k_values(psis_fit)),
      relative_mcse = sd(shifted) / sqrt(length(shifted)) / mean(shifted),
      score = mean_score,
      score_mc_covariance = crossprod(weighted_scores)
    )
  }

  list(evaluate = evaluate, proposal_at = proposal_at, layout = layout,
       latent_dimension = latent$dimension)
}

##' Estimate score Monte Carlo covariance for a moving proposal
##' @noRd
adaptive_mcml_pathwise_score <- function(
    engine, parameters, z, proposal_parameters = NULL) {
  proposal_for <- function(candidate) {
    if (is.null(proposal_parameters)) candidate else proposal_parameters
  }
  base_log_weights <- engine(
    parameters, z, log_weights_only = TRUE,
    proposal_parameters = proposal_for(parameters)
  )
  weights <- normalise_log_weights(base_log_weights)
  step <- .Machine$double.eps^(1 / 3) * pmax(1, abs(parameters))
  sample_scores <- matrix(
    NA_real_, nrow = length(base_log_weights), ncol = length(parameters)
  )
  for (j in seq_along(parameters)) {
    upper <- lower <- parameters
    upper[j] <- upper[j] + step[j]
    lower[j] <- lower[j] - step[j]
    sample_scores[, j] <- (
      engine(upper, z, log_weights_only = TRUE,
             proposal_parameters = proposal_for(upper)) -
        engine(lower, z, log_weights_only = TRUE,
               proposal_parameters = proposal_for(lower))
    ) / (2 * step[j])
  }
  mean_score <- colSums(sweep(sample_scores, 1L, weights, "*"))
  centered_scores <- sweep(sample_scores, 2L, mean_score, "-")
  weighted_scores <- sweep(centered_scores, 1L, weights, "*")
  list(
    score = mean_score,
    covariance = crossprod(weighted_scores)
  )
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
  # A parameter that only acts through a zero variance (the spatial range when
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
           control = list(eval.max = 150, iter.max = 100,
                          rel.tol = 1e-8, x.tol = 1e-6)),
    error = function(error) NULL
  )
  improved <- !is.null(restart) && is.finite(restart$objective) &&
    restart$objective < fit$objective - 1e-6
  list(fit = if (improved) restart else fit, restarted = TRUE,
       restart_improved = improved)
}

##' Monte Carlo precision of the free parameters
##'
##' Parameters on a variance boundary are held fixed, so the observed
##' information, the worst-direction MCSE ratio and the independent score
##' comparison are computed over the remaining (free) parameters only.
##' @noRd
adaptive_precision_diagnostics <- function(engine, estimate, fit_bank,
                                           validation_bank, free,
                                           score_agreement_probability) {
  out <- list(
    information = NULL, covariance = NULL, information_singular = TRUE,
    fit_pathwise_score = NULL, validation_pathwise_score = NULL,
    parameter_mcse_ratio = Inf, score_agreement_statistic = Inf,
    score_agreement_limit = stats::qchisq(score_agreement_probability,
                                          df = length(free)),
    timing = list(information = 0, fit_pathwise_score = 0,
                  validation_pathwise_score = 0)
  )
  started <- proc.time()[["elapsed"]]
  out$information <- tryCatch(
    # Two Richardson rounds suffice for this smooth, deterministic
    # objective: standard errors agree with the default four rounds to about
    # 1e-4 relative error at half the cost.
    numDeriv::hessian(function(parameters) -engine(parameters, fit_bank),
                      estimate, method.args = list(r = 2)),
    error = function(error) NULL
  )
  out$timing$information <- proc.time()[["elapsed"]] - started
  if (is.null(out$information)) return(out)
  information_root <- tryCatch(
    factor_covariance(out$information[free, free, drop = FALSE],
                      "adaptive MCML observed information",
                      allow_jitter = FALSE),
    error = function(error) NULL
  )
  if (is.null(information_root)) return(out)
  out$information_singular <- FALSE
  free_covariance <- chol2inv(information_root)
  out$covariance <- matrix(NA_real_, length(estimate), length(estimate))
  out$covariance[free, free] <- free_covariance

  started <- proc.time()[["elapsed"]]
  out$fit_pathwise_score <- adaptive_mcml_pathwise_score(
    engine, estimate, fit_bank
  )
  out$timing$fit_pathwise_score <- proc.time()[["elapsed"]] - started
  started <- proc.time()[["elapsed"]]
  out$validation_pathwise_score <- adaptive_mcml_pathwise_score(
    engine, estimate, validation_bank
  )
  out$timing$validation_pathwise_score <- proc.time()[["elapsed"]] - started

  parameter_ratio <- function(score_mc_covariance) {
    parameter_mc_covariance <- free_covariance %*%
      score_mc_covariance[free, free, drop = FALSE] %*% free_covariance
    relative_covariance <- information_root %*%
      parameter_mc_covariance %*% t(information_root)
    relative_covariance <- (relative_covariance + t(relative_covariance)) / 2
    sqrt(max(eigen(relative_covariance, symmetric = TRUE,
                   only.values = TRUE)$values))
  }
  out$parameter_mcse_ratio <- max(
    parameter_ratio(out$fit_pathwise_score$covariance),
    parameter_ratio(out$validation_pathwise_score$covariance)
  )
  score_difference <- (out$fit_pathwise_score$score -
                         out$validation_pathwise_score$score)[free]
  score_root <- tryCatch(
    factor_covariance(
      (out$fit_pathwise_score$covariance +
         out$validation_pathwise_score$covariance)[free, free, drop = FALSE],
      "adaptive MCML score-difference covariance"
    ),
    error = function(error) NULL
  )
  if (!is.null(score_root)) {
    out$score_agreement_statistic <- sum(
      forwardsolve(t(score_root), score_difference)^2
    )
  }
  out
}

##' Fit and independently validate adaptive local-Laplace MCML
##' @noRd
fit_adaptive_mcml <- function(engine, start, latent_dimension,
                            budgets = NULL, initial_budget = 1000L,
                            maximum_budget = 20000L, maximum_stages = 4L,
                            seed = 1L, min_relative_ess = 0.1,
                            min_absolute_ess = NULL,
                            max_pareto_k = 0.7,
                            max_parameter_mcse_ratio = 0.1,
                            score_agreement_probability = 0.99,
                            lower = NULL, upper = NULL,
                            variance_index = integer(),
                            linked_index = list(),
                            boundary_log_variance = -9,
                            range_index = NULL,
                            minimum_log_range = -Inf,
                            near_zero_log_variance = -4) {
  if (is.null(lower)) lower <- rep(-10, length(start))
  if (is.null(upper)) upper <- rep(10, length(start))
  if (!is.null(min_absolute_ess)) {
    stopifnot(length(min_absolute_ess) == 1L, is.finite(min_absolute_ess),
              min_absolute_ess > 0)
  }
  history <- list()
  current <- start
  total_started <- proc.time()[["elapsed"]]
  dynamic_budget <- is.null(budgets)
  if (dynamic_budget) budgets <- as.integer(initial_budget)
  maximum_stages <- if (dynamic_budget) maximum_stages else length(budgets)
  stage_seed <- function(offset) {
    as.integer(((seed - 1 + offset) %% .Machine$integer.max) + 1)
  }

  for (stage in seq_len(maximum_stages)) {
    budget <- if (dynamic_budget) budgets[length(budgets)] else budgets[stage]
    set.seed(stage_seed(2L * stage - 1L))
    fit_bank <- matrix(rnorm(budget * latent_dimension), nrow = budget)
    stage_started <- proc.time()[["elapsed"]]
    optimizer_started <- proc.time()[["elapsed"]]
    fit <- nlminb(
      current,
      function(parameters) -engine(parameters, fit_bank),
      lower = lower, upper = upper,
      control = list(eval.max = 150, iter.max = 100,
                     rel.tol = 1e-8, x.tol = 1e-6)
    )
    restart <- if (stage == 1L && length(variance_index)) {
      restart_from_interior(
        function(parameters) -engine(parameters, fit_bank), fit,
        variance_index, boundary_log_variance, lower, upper
      )
    } else {
      list(fit = fit, restarted = FALSE)
    }
    fit <- restart$fit
    boundary <- adaptive_boundary_parameters(
      fit$par, variance_index, boundary_log_variance, linked_index
    )
    # A range shorter than the spacing between locations cannot be
    # distinguished from unstructured variation. Unless the spatial variance is
    # already zero (when the range is meaningless anyway), the range is held at
    # that spacing and the remaining parameters are re-optimised.
    if (!is.null(range_index) && !(range_index %in% boundary) &&
        fit$par[range_index] < minimum_log_range) {
      lower[range_index] <- minimum_log_range
      held <- tryCatch(
        nlminb(pmax(fit$par, lower),
               function(parameters) -engine(parameters, fit_bank),
               lower = lower, upper = upper,
               control = list(eval.max = 150, iter.max = 100,
                              rel.tol = 1e-8, x.tol = 1e-6)),
        error = function(error) NULL
      )
      if (!is.null(held)) fit <- held
      boundary <- adaptive_boundary_parameters(
        fit$par, variance_index, boundary_log_variance, linked_index
      )
    }
    range_held <- !is.null(range_index) && !(range_index %in% boundary) &&
      fit$par[range_index] <= minimum_log_range + 1e-6
    if (range_held) boundary <- sort(c(boundary, range_index))
    free <- setdiff(seq_along(fit$par), boundary)
    optimizer_seconds <- proc.time()[["elapsed"]] - optimizer_started
    fit_diagnostics_started <- proc.time()[["elapsed"]]
    fit_diagnostics <- engine(fit$par, fit_bank, diagnostics = TRUE)
    fit_diagnostics_seconds <-
      proc.time()[["elapsed"]] - fit_diagnostics_started

    set.seed(stage_seed(2L * stage))
    validation_bank <- matrix(rnorm(budget * latent_dimension), nrow = budget)
    validation_diagnostics_started <- proc.time()[["elapsed"]]
    validation_diagnostics <- engine(
      fit$par, validation_bank, diagnostics = TRUE
    )
    validation_diagnostics_seconds <-
      proc.time()[["elapsed"]] - validation_diagnostics_started
    if (is.null(fit_diagnostics) || is.null(validation_diagnostics)) {
      history[[stage]] <- list(
        stage = stage, budget = budget, estimate = fit$par,
        convergence = fit$convergence,
        objective = fit$objective, iterations = fit$iterations,
        evaluations = fit$evaluations, message = fit$message,
        fit_diagnostics = fit_diagnostics,
        validation_diagnostics = validation_diagnostics,
        likelihood_difference = Inf, difference_limit = NA_real_,
        accepted = FALSE,
        timing = list(
          optimizer = optimizer_seconds,
          fit_diagnostics = fit_diagnostics_seconds,
          validation_diagnostics = validation_diagnostics_seconds
        ),
        elapsed_seconds = proc.time()[["elapsed"]] - stage_started
      )
      current <- fit$par
      next
    }
    likelihood_difference <- abs(
      fit_diagnostics$value - validation_diagnostics$value
    )
    difference_limit <- max(
      0.01,
      3 * sqrt(fit_diagnostics$relative_mcse^2 +
                 validation_diagnostics$relative_mcse^2)
    )
    fit_effective_samples <- fit_diagnostics$relative_ess * budget
    validation_effective_samples <-
      validation_diagnostics$relative_ess * budget
    sufficient_ess <- if (is.null(min_absolute_ess)) {
      fit_diagnostics$relative_ess >= min_relative_ess &&
        validation_diagnostics$relative_ess >= min_relative_ess
    } else {
      fit_effective_samples >= min_absolute_ess &&
        validation_effective_samples >= min_absolute_ess
    }
    preliminary_pass <- isTRUE(
      fit$convergence == 0L &&
      sufficient_ess &&
      fit_diagnostics$pareto_k <= max_pareto_k &&
      validation_diagnostics$pareto_k <= max_pareto_k
    )
    precision <- if (preliminary_pass) {
      adaptive_precision_diagnostics(
        engine, fit$par, fit_bank, validation_bank, free,
        score_agreement_probability
      )
    } else {
      list(information = NULL, covariance = NULL,
           information_singular = NA,
           fit_pathwise_score = NULL, validation_pathwise_score = NULL,
           parameter_mcse_ratio = Inf, score_agreement_statistic = Inf,
           score_agreement_limit = stats::qchisq(
             score_agreement_probability, df = length(free)
           ),
           timing = list(information = 0, fit_pathwise_score = 0,
                         validation_pathwise_score = 0))
    }
    # If the information is singular over the free parameters, variance
    # components that are small but not exactly zero are the usual cause. Hold
    # them as well before concluding that the model is not identifiable.
    near_zero <- setdiff(
      variance_index[fit$par[variance_index] <= near_zero_log_variance],
      boundary
    )
    if (preliminary_pass && isTRUE(precision$information_singular) &&
        length(near_zero)) {
      boundary <- sort(unique(c(
        boundary, near_zero,
        unlist(linked_index[as.character(near_zero)], use.names = FALSE)
      )))
      free <- setdiff(seq_along(fit$par), boundary)
      precision <- adaptive_precision_diagnostics(
        engine, fit$par, fit_bank, validation_bank, free,
        score_agreement_probability
      )
    }
    information <- precision$information
    covariance <- precision$covariance
    fit_pathwise_score <- precision$fit_pathwise_score
    validation_pathwise_score <- precision$validation_pathwise_score
    parameter_mcse_ratio <- precision$parameter_mcse_ratio
    score_agreement_statistic <- precision$score_agreement_statistic
    score_agreement_limit <- precision$score_agreement_limit
    information_seconds <- precision$timing$information
    fit_score_seconds <- precision$timing$fit_pathwise_score
    validation_score_seconds <- precision$timing$validation_pathwise_score
    # A singular information matrix over the free parameters is a property of
    # the data and model (weak identification), not of the Monte Carlo sample.
    identifiability_failure <- isTRUE(precision$information_singular)
    accepted <- preliminary_pass && is.finite(parameter_mcse_ratio) &&
      parameter_mcse_ratio <= max_parameter_mcse_ratio &&
      score_agreement_statistic <= score_agreement_limit

    history[[stage]] <- list(
      stage = stage, budget = budget, estimate = fit$par,
      convergence = fit$convergence,
      objective = fit$objective, iterations = fit$iterations,
      evaluations = fit$evaluations, message = fit$message,
      fit_diagnostics = fit_diagnostics,
      validation_diagnostics = validation_diagnostics,
      likelihood_difference = likelihood_difference,
      difference_limit = difference_limit,
      fit_effective_samples = fit_effective_samples,
      validation_effective_samples = validation_effective_samples,
      information = information,
      covariance = covariance,
      fit_pathwise_score = fit_pathwise_score,
      validation_pathwise_score = validation_pathwise_score,
      max_parameter_mcse_ratio = parameter_mcse_ratio,
      score_agreement_statistic = score_agreement_statistic,
      score_agreement_limit = score_agreement_limit,
      boundary_parameters = boundary,
      range_held = range_held,
      restarted = restart$restarted,
      restart_improved = isTRUE(restart$restart_improved),
      identifiability_failure = identifiability_failure,
      accepted = accepted,
      timing = list(
        optimizer = optimizer_seconds,
        fit_diagnostics = fit_diagnostics_seconds,
        validation_diagnostics = validation_diagnostics_seconds,
        information = information_seconds,
        fit_pathwise_score = fit_score_seconds,
        validation_pathwise_score = validation_score_seconds
      ),
      elapsed_seconds = proc.time()[["elapsed"]] - stage_started
    )
    current <- fit$par
    if (accepted) break
    # More draws cannot repair a non-identifiable information matrix.
    if (identifiability_failure) {
      history[[stage]]$budget_reason <- "singular observed information"
      break
    }

    if (dynamic_budget) {
      # Independent-importance MCSE decreases approximately with sqrt(N).
      # Estimate the next useful budget instead of following a fixed ladder;
      # the safety factor protects against variability in the pilot estimate.
      parameter_budget <- if (is.finite(parameter_mcse_ratio)) {
        ceiling(
          budget *
            (parameter_mcse_ratio / max_parameter_mcse_ratio)^2 * 1.25
        )
      } else {
        budget
      }
      disagreement <- likelihood_difference > difference_limit
      score_disagreement <-
        score_agreement_statistic > score_agreement_limit
      proposal_failure <-
        !sufficient_ess ||
        fit_diagnostics$pareto_k > max_pareto_k ||
        validation_diagnostics$pareto_k > max_pareto_k

      # More draws reduce Monte Carlo error but do not repair a systematically
      # poor proposal. Permit one larger confirmation stage; a persistent
      # proposal failure should trigger a different proposal or MCMC fallback.
      next_budget <- max(
        parameter_budget,
        if (score_disagreement || proposal_failure || fit$convergence != 0L) {
          2L * budget
        } else {
          budget + 500L
        }
      )
      next_budget <- as.integer(500L * ceiling(next_budget / 500L))
      next_budget <- min(next_budget, as.integer(maximum_budget))
      history[[stage]]$next_budget <- next_budget
      history[[stage]]$budget_reason <- paste(
        c(
          if (disagreement) "likelihood-bank disagreement (diagnostic only)",
          if (score_disagreement) "independent score disagreement",
          if (proposal_failure) "proposal diagnostics",
          if (parameter_mcse_ratio > max_parameter_mcse_ratio) {
            "parameter Monte Carlo error"
          },
          if (fit$convergence != 0L) "optimizer convergence"
        ),
        collapse = "; "
      )
      if (next_budget <= budget) break
      budgets <- c(budgets, next_budget)
    }
  }

  last <- history[[length(history)]]
  covariance <- last$covariance
  information <- last$information
  if (is.null(covariance)) {
    information <- tryCatch(
      numDeriv::hessian(function(parameters) -engine(parameters, fit_bank),
                        current, method.args = list(r = 2)),
      error = function(error) NULL
    )
    free <- setdiff(seq_along(current), last$boundary_parameters)
    information_root <- if (is.null(information)) NULL else tryCatch(
      factor_covariance(
        information[free, free, drop = FALSE],
        "adaptive MCML observed information",
        allow_jitter = FALSE
      ),
      error = function(error) NULL
    )
    if (!is.null(information_root)) {
      covariance <- matrix(NA_real_, length(current), length(current))
      covariance[free, free] <- chol2inv(information_root)
    }
  }

  list(
    estimate = current,
    accepted = isTRUE(last$accepted),
    budget = last$budget,
    boundary_parameters = last$boundary_parameters,
    range_held = isTRUE(last$range_held),
    identifiability_failure = isTRUE(last$identifiability_failure),
    history = history,
    covariance = covariance,
    information = information,
    elapsed_seconds = proc.time()[["elapsed"]] - total_started
  )
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
  # The range is held only because the spatial variance is zero; it is
  # explained separately rather than listed as a variance component.
  index <- setdiff(index, layout$index$phi)
  unname(labels[as.character(index)])
}

##' Retain user-facing adaptive diagnostics without derivative workspaces
##' @noRd
compact_adaptive_mcml_history <- function(history) {
  compact_bank <- function(diagnostics) {
    if (is.null(diagnostics)) return(NULL)
    diagnostics[c(
      "value", "relative_ess", "maximum_weight", "pareto_k", "relative_mcse"
    )]
  }
  lapply(history, function(stage) {
    stage$fit_diagnostics <- compact_bank(stage$fit_diagnostics)
    stage$validation_diagnostics <- compact_bank(stage$validation_diagnostics)
    stage$information <- NULL
    stage$covariance <- NULL
    stage$fit_pathwise_score <- NULL
    stage$validation_pathwise_score <- NULL
    stage
  })
}

##' Maximise the Laplace approximation from several starting values
##'
##' With a single standard-normal draw at zero the adaptive engine evaluates
##' the Laplace approximation to the marginal log-likelihood exactly. The
##' Laplace optimum is cheap to compute, is usually close to the MCML optimum,
##' and optimising it from several candidates guards against local optima and
##' spurious variance boundaries before any Monte Carlo work is done.
##' @noRd
laplace_multistart <- function(engine, candidates, latent_dimension,
                               lower, upper) {
  zero <- matrix(0, nrow = 1L, ncol = latent_dimension)
  objective <- function(parameters) {
    value <- engine(parameters, zero)
    if (is.finite(value)) -value else 1e100
  }
  fits <- lapply(candidates, function(candidate) {
    tryCatch(
      nlminb(pmin(pmax(candidate, lower), upper), objective,
             lower = lower, upper = upper,
             control = list(eval.max = 200, iter.max = 150,
                            rel.tol = 1e-8, x.tol = 1e-6)),
      error = function(error) NULL
    )
  })
  values <- vapply(fits, function(fit) {
    if (is.null(fit)) Inf else fit$objective
  }, numeric(1))
  if (!any(is.finite(values))) return(NULL)
  best <- which.min(values)
  list(estimate = fits[[best]]$par, log_likelihood = -values[[best]],
       candidate = names(candidates)[best],
       candidate_log_likelihoods = -values)
}

##' Candidate starting values for the Laplace multistart
##'
##' The supplied start is always included. Alternatives replace the range by a
##' data-based value (a third of the median inter-location distance) and move
##' every variance component to a moderate interior value, so that a start in
##' the basin of a zero-variance boundary does not determine the fit.
##' @noRd
laplace_start_candidates <- function(start, layout, distances) {
  candidates <- list(supplied = start)
  median_distance <- stats::median(distances[upper.tri(distances)])
  if (is.finite(median_distance) && median_distance > 0) {
    data_range <- start
    data_range[layout$index$phi] <- log(median_distance / 3)
    candidates$data_range <- data_range
  }
  interior <- start
  interior[layout$index$sigma2] <- max(start[layout$index$sigma2], log(0.25))
  if (!is.null(layout$index$nu2)) interior[layout$index$nu2] <- log(0.25)
  if (!is.null(layout$index$sigma2_re)) {
    interior[layout$index$sigma2_re] <- log(0.25)
  }
  candidates$interior_variances <- interior
  candidates
}

##' Fit a non-Gaussian geostatistical model by adaptive validated MCML
##' @noRd
glgpm_nong_adaptive <- function(
    y, D, coords, units_m, kappa, reference_pars, cov_offset,
    ID_coords, ID_re, fix_tau2, family, invlink, control_mcml,
    seed = NULL, messages = TRUE) {
  stopifnot(family %in% c("poisson", "binomial"))
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  fit_seed <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, 1L)
  } else {
    as.integer(seed)
  }

  distances <- as.matrix(pairwise_distances(coords))
  engine <- make_adaptive_mcml_engine(
    D = D, y = y, units_m = units_m, distances = distances,
    family = family, kappa = kappa, cov_offset = cov_offset,
    ID_coords = ID_coords, fix_tau2 = fix_tau2, ID_re = ID_re,
    invlink = invlink
  )
  start <- c(
    reference_pars$beta,
    log(reference_pars$sigma2),
    log(reference_pars$phi)
  )
  if (isTRUE(fix_tau2)) {
    start <- c(start, log(reference_pars$tau2 / reference_pars$sigma2))
  }
  if (length(reference_pars$sigma2_re)) {
    start <- c(start, log(reference_pars$sigma2_re))
  }

  # Coordinates and the initial range are scaled before this function is
  # called. These broad bounds prevent overflow without constraining any
  # scientifically plausible fit on the internal scale.
  lower <- c(rep(-30, ncol(D)), -20, -20,
             rep(-20, length(start) - ncol(D) - 2L))
  upper <- c(rep(30, ncol(D)), 20, 10,
             rep(20, length(start) - ncol(D) - 2L))
  laplace <- NULL
  if (isTRUE(control_mcml$laplace_start)) {
    laplace <- laplace_multistart(
      engine$evaluate,
      laplace_start_candidates(unname(start), engine$layout, distances),
      engine$latent_dimension, lower, upper
    )
    if (!is.null(laplace)) start <- laplace$estimate
  }
  fit <- fit_adaptive_mcml(
    engine = engine$evaluate, start = unname(start),
    latent_dimension = engine$latent_dimension,
    initial_budget = control_mcml$initial_samples,
    maximum_budget = control_mcml$max_samples,
    maximum_stages = control_mcml$max_stages,
    seed = fit_seed,
    min_absolute_ess = control_mcml$min_effective_samples,
    max_pareto_k = control_mcml$max_pareto_k,
    max_parameter_mcse_ratio = control_mcml$max_mcse_ratio,
    score_agreement_probability = control_mcml$score_probability,
    lower = lower, upper = upper,
    variance_index = unlist(engine$layout$index[c("sigma2", "nu2",
                                                  "sigma2_re")],
                            use.names = FALSE),
    linked_index = stats::setNames(
      list(engine$layout$index$phi), engine$layout$index$sigma2
    ),
    range_index = engine$layout$index$phi,
    minimum_log_range = if (isFALSE(control_mcml$hold_short_range)) -Inf else
      log(median_nearest_neighbour_distance(distances))
  )
  last <- fit$history[[length(fit$history)]]
  public_history <- compact_adaptive_mcml_history(fit$history)
  estimate <- structure_estimate(
    fit$estimate, beta_names = colnames(D), fix_tau2 = fix_tau2,
    sigma2_me = FALSE,
    re_names = if (is.null(ID_re)) NULL else colnames(ID_re)
  )
  covariance <- fit$covariance
  if (is.null(covariance)) {
    covariance <- matrix(
      NA_real_, nrow = length(fit$estimate), ncol = length(fit$estimate)
    )
  }
  flat_names <- names(unlist(estimate))
  dimnames(covariance) <- list(flat_names, flat_names)
  # The last stage already estimated the gradient at the estimate on its fitting
  # bank (the pathwise score), so it is reported directly. Otherwise the bank
  # is recreated and a central difference is used.
  last_score <- fit$history[[length(fit$history)]]$fit_pathwise_score$score
  gradient <- if (!is.null(last_score) && length(last_score) ==
                  length(fit$estimate)) {
    last_score
  } else {
    final_seed <- as.integer((
      (fit_seed - 1 + 2L * length(fit$history) - 1L) %% .Machine$integer.max
    ) + 1)
    set.seed(final_seed)
    final_bank <- matrix(
      rnorm(fit$budget * engine$latent_dimension), nrow = fit$budget
    )
    tryCatch(
      numDeriv::grad(
        function(parameters) engine$evaluate(parameters, final_bank),
        fit$estimate, method = "simple", method.args = list(eps = 1e-6)
      ),
      error = function(error) rep(NA_real_, length(fit$estimate))
    )
  }

  boundary_names <- describe_variance_parameters(
    fit$boundary_parameters, engine$layout,
    re_names = if (is.null(ID_re)) NULL else colnames(ID_re)
  )
  if (length(boundary_names)) {
    warning(
      "The ", if (length(boundary_names) == 1L) boundary_names else paste(
        paste(boundary_names[-length(boundary_names)], collapse = ", "),
        "and the", boundary_names[length(boundary_names)]
      ), " ",
      if (length(boundary_names) == 1L) "was" else "were",
      " estimated as effectively zero. ",
      if (length(boundary_names) == 1L) "This component" else "These components",
      " cannot be distinguished from zero with these data, so ",
      if (length(boundary_names) == 1L) "its standard error is" else
        "their standard errors are",
      " not reported",
      if (engine$layout$index$sigma2 %in% fit$boundary_parameters) {
        ", and the spatial range cannot be estimated"
      },
      ". Consider refitting without ",
      if (length(boundary_names) == 1L) "it." else "them.",
      call. = FALSE
    )
  }
  if (fit$range_held) {
    warning(
      "The estimated spatial range is shorter than the typical distance ",
      "between neighbouring locations, so spatial correlation cannot be ",
      "distinguished from unstructured variation with these data. The range ",
      "was held at the median nearest-neighbour distance and its standard ",
      "error is not reported. Consider a model with a nugget (unstructured) ",
      "effect instead, or fixing the range from external knowledge.",
      call. = FALSE
    )
  }
  if (fit$identifiability_failure) {
    warning(
      "The observed information matrix is singular, so some parameters ",
      "(often the spatial range when spatial correlation is weak) are not ",
      "identifiable from these data. The estimate is returned without ",
      "standard errors; consider a simpler model or fixing the range.",
      call. = FALSE
    )
  } else if (!fit$accepted) {
    warning(
      "Adaptive MCML reached its stopping limit without meeting all ",
      "Monte Carlo precision and independent-validation criteria. The best ",
      "available estimate is returned; inspect 'mcml_history' before using ",
      "its uncertainty estimates.",
      call. = FALSE
    )
  } else if (messages) {
    message("Adaptive MCML met its validation criteria after ",
            length(fit$history), " stage",
            if (length(fit$history) == 1L) "" else "s", " and ",
            format(fit$budget, big.mark = ","), " draws per bank.")
  }

  out <- list(
    estimate = estimate,
    grad_MLE = gradient,
    covariance = covariance,
    log_lik = last$fit_diagnostics$value,
    S_samples = NULL,
    link_function = make_inverse_link_functions(family, invlink)
  )
  optimizer <- list(
    convergence = last$convergence,
    message = last$message,
    evaluations = last$evaluations,
    invalid_evaluations = NA_integer_,
    max_abs_gradient = if (all(is.finite(gradient))) max(abs(gradient)) else NA_real_,
    stationary = all(is.finite(gradient)) && max(abs(gradient)) <= 1e-3,
    information_rcond = if (is.null(fit$information)) NA_real_ else
      rcond(fit$information),
    information_jitter = 0,
    mcml_supported = fit$accepted,
    mcml_range_held = fit$range_held,
    laplace_start = if (is.null(laplace)) NULL else list(
      estimate = laplace$estimate,
      log_likelihood = laplace$log_likelihood,
      candidate = laplace$candidate,
      candidate_log_likelihoods = laplace$candidate_log_likelihoods
    ),
    mcml_boundary = boundary_names,
    mcml_identifiability_failure = fit$identifiability_failure,
    mcml_budget = fit$budget,
    mcml_stages = length(fit$history),
    mcml_history = public_history,
    seed = fit_seed,
    elapsed_seconds = fit$elapsed_seconds
  )
  attr(out, "optimizer") <- optimizer
  class(out) <- "RiskMap"
  out
}

