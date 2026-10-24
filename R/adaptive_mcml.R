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
##' importance
##' weights have infinite variance whenever the largest eigenvalue of
##' Sigma^{1/2} W Sigma^{1/2} exceeds one (W the likelihood curvature), which
##' holds in practically every spatial dataset (Pitt et al., 2013). The mixture
##' bounds every moment of the weights, so that the Monte Carlo error estimates
##' used for stopping are valid, at a small cost in efficiency. Weights use the
##' full mixture density. With `defensive_weight = 0` the proposal is the
##' Gaussian component alone.
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
  cached_type <- NULL
  cached_proposal <- NULL

  # Deviations from the centre with the prior covariance of the latent effects:
  # the spatial block through its Cholesky factor, random-effect levels
  # scaled by their standard deviations.
  prior_deviations <- function(z, proposal) {
    deviations <- matrix(0, nrow(z), ncol(z))
    deviations[, latent$spatial] <-
      z[, latent$spatial, drop = FALSE] %*% proposal$spatial_root
    if (n_re) for (j in seq_len(n_re)) {
      index <- latent$random[[j]]
      deviations[, index] <- z[, index, drop = FALSE] *
        sqrt(proposal$values$sigma2_re[j])
    }
    deviations
  }

  # Squared norms of the draws standardised under each mixture component. A
  # draw's residual under its own component is its standard-normal vector, so
  # only the other component needs a triangular solve. For the Laplace draws
  # the spatial residuals under the prior can be taken from the samples
  # already standardised for the prior density (standardised_spatial). The
  # realised share of defensive draws is used as the mixture weight, so the
  # estimate is unbiased for any bank size. A bank evaluated in row blocks
  # passes the share of the whole bank.
  log_mixture_proposal <- function(z, deviations, defensive, proposal,
                                   standardised_spatial = NULL,
                                   share = mean(defensive)) {
    z_norms <- rowSums(z^2)
    log_normalising <- latent$dimension * log(2 * pi)
    log_laplace_determinant <- 2 * sum(log(diag(proposal$root)))
    if (share == 0) {
      # The draws are centre + z %*% root, so the residuals are exactly z.
      return(-0.5 * (log_normalising + log_laplace_determinant + z_norms))
    }
    laplace_norms <- z_norms
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
    if (n_re) for (j in seq_len(n_re)) {
      index <- latent$random[[j]]
      spatial_norms <- spatial_norms +
        rowSums(deviations[laplace_draws, index, drop = FALSE]^2) /
        proposal$values$sigma2_re[j]
    }
    prior_norms <- z_norms
    prior_norms[laplace_draws] <- spatial_norms
    log_prior_determinant <- 2 * sum(log(diag(proposal$spatial_root))) +
      if (n_re) sum(n_dim_re * log(proposal$values$sigma2_re)) else 0
    log_laplace <- log1p(-share) - 0.5 * (
      log_normalising + log_laplace_determinant + laplace_norms
    )
    log_defensive <- log(share) - 0.5 * (
      log_normalising + log_prior_determinant + prior_norms
    )
    largest <- pmax(log_laplace, log_defensive)
    largest + log(exp(log_laplace - largest) + exp(log_defensive - largest))
  }

  proposal_at <- function(parameters, type = default_proposal_type) {
    if (!is.null(cached_parameters) &&
        identical(parameters, cached_parameters) &&
        identical(type, cached_type)) {
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
    centre <- approximation$mode
    proposal_precision <- approximation$hessian
    global <- NULL
    if (identical(type, "global")) {
      global <- fit_global_proposal(
        mode = approximation$mode, mu = mu, y = y, units_m = units_m,
        latent_map = latent_map, precision = precision, family = family,
        linkf = linkf, canonical = canonical
      )
      if (is.null(global) || !global$converged) {
        proposal_failures <<- proposal_failures + 1L
      }
      if (!is.null(global)) {
        centre <- global$centre
        proposal_precision <- global$precision
      }
    }
    root <- tryCatch(
      inverse_cholesky_factor(proposal_precision),
      error = function(error) factor_covariance(
        solve(proposal_precision), "adaptive MCML proposal"
      )
    )
    proposal <- list(centre = centre, root = root, values = values,
                     covariance = covariance, spatial_root = spatial_root,
                     type = if (is.null(global)) "laplace" else "global",
                     mode = approximation$mode)
    cached_parameters <<- parameters
    cached_type <<- type
    cached_proposal <<- proposal
    proposal
  }

  # `built_proposal`, if given, is the proposal at `proposal_parameters`, so
  # that a bank evaluated in row blocks builds each proposal once. Such a bank
  # also passes its `defensive` rows and their share of the whole bank.
  evaluate <- function(parameters, z, log_weights_only = FALSE,
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
    if (is.null(proposal)) return(-1e100)
    if (is.null(defensive)) {
      defensive <- defensive_rows(nrow(z), defensive_weight)
    }
    if (is.null(defensive_share)) defensive_share <- mean(defensive)
    deviations <- z %*% proposal$root
    if (any(defensive)) {
      deviations[defensive, ] <- prior_deviations(
        z[defensive, , drop = FALSE], proposal
      )
    }
    samples <- deviations + rep(proposal$centre, each = nrow(deviations))
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
    if (is.null(spatial_root)) return(-1e100)

    spatial_samples <- samples[, latent$spatial, drop = FALSE]
    # Samples standardised by the prior factor. Their squared norms give the
    # prior quadratic form and, when the target and proposal parameters
    # coincide, they also give the residuals under the defensive component.
    standardised_samples <- forwardsolve(t(spatial_root), t(spatial_samples))
    log_prior <- -0.5 * (
      n_loc * log(2 * pi) + 2 * sum(log(diag(spatial_root))) +
        colSums(standardised_samples^2)
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

    log_proposal <- log_mixture_proposal(
      z, deviations, defensive, proposal,
      if (identical(parameters, proposal_parameters)) standardised_samples,
      share = defensive_share
    )
    log_weights <- log_target - log_proposal
    if (log_weights_only) return(log_weights)
    maximum <- max(log_weights)
    maximum + log(mean(exp(log_weights - maximum)))
  }

  list(evaluate = evaluate, proposal_at = proposal_at, layout = layout,
       latent_dimension = latent$dimension,
       defensive_weight = defensive_weight,
       proposal_type = function() default_proposal_type,
       use_proposal = function(type) {
         default_proposal_type <<- match.arg(type, c("global", "laplace"))
         invisible(NULL)
       },
       proposal_failures = function() proposal_failures,
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
           control = mcml_optimiser_control()),
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
           upper = upper, control = mcml_optimiser_control()),
    error = function(error) NULL
  )
  if (is.null(fit) || !is.finite(fit$objective) || fit$objective >= 1e99) {
    return(NULL)
  }
  list(estimate = fit$par, log_likelihood = -fit$objective,
       convergence = fit$convergence)
}

##' Dispersed starting values for the optimisation check
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
##' units.
##' @noRd
dispersed_starts <- function(layout, distances, distance_scale = 1,
                             messages = FALSE, trigger_fraction = 0.5) {
  spacing <- median_nearest_neighbour_distance(distances)
  grid <- expand.grid(sigma2 = c(0.25, 4), phi = c(spacing, max(distances) / 2))
  announce <- function(log_phi) {
    if (!messages) return(invisible(NULL))
    message(
      "The fitted spatial scale (phi = ",
      signif(exp(log_phi) * distance_scale, 3), ") is below half the median ",
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
##' P(phi < range_threshold) = range_probability. Fixed effects, nugget and
##' unstructured random-effect variances are not penalised.
##' @noRd
make_pc_penalty <- function(layout, range_threshold,
                            sigma_threshold = 1, sigma_probability = 0.05,
                            range_probability = 0.05) {
  rate_sigma <- -log(sigma_probability) / sigma_threshold
  rate_range <- -log(range_probability) * range_threshold
  function(parameters) {
    log_sigma2 <- parameters[layout$index$sigma2]
    log_phi <- parameters[layout$index$phi]
    sigma <- exp(log_sigma2 / 2)
    phi <- exp(log_phi)
    # pi(sigma) = rate exp(-rate sigma), with d sigma / d log sigma^2 = sigma/2;
    # pi(phi) = rate phi^-2 exp(-rate / phi), with d phi / d log phi = phi.
    (log(rate_sigma) - rate_sigma * sigma + log(sigma / 2)) +
      (log(rate_range) - log_phi - rate_range / phi)
  }
}

##' Add a penalty to an adaptive MCML engine
##'
##' The penalty depends only on the parameters, so it is the same for every
##' draw: importance weights, effective sample sizes and Monte Carlo error are
##' unchanged, while the objective and its curvature include the penalty.
##' @noRd
penalise_engine <- function(evaluate, penalty) {
  function(parameters, z, log_weights_only = FALSE, ...) {
    out <- evaluate(parameters, z, log_weights_only = log_weights_only, ...)
    if (log_weights_only) return(out)
    out + penalty(parameters)
  }
}

##' Starting values on the working scale of the MCML engine
##'
##' Regression coefficients, then the log spatial variance and log scale,
##' the log nugget ratio when the nugget is estimated, and the log
##' random-effect variances.
##' @noRd
mcml_working_start <- function(reference_pars, fix_tau2) {
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
##' estimate and returns the corrected, certified estimate.
##' @noRd
glgpm_nong_adaptive <- function(
    y, D, coords, units_m, kappa, reference_pars, cov_offset,
    ID_coords, ID_re, fix_tau2, family, invlink, control_mcml,
    seed = NULL, messages = TRUE, distance_scale = 1) {
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
    invlink = invlink,
    defensive_weight = control_mcml$defensive_weight %||% 0.05,
    proposal = control_mcml$proposal %||% "global"
  )
  start <- mcml_working_start(reference_pars, fix_tau2)
  # Coordinates and the initial scale are rescaled before this function is
  # called. These broad bounds prevent overflow without constraining any
  # scientifically plausible fit on the internal scale.
  lower <- c(rep(-30, ncol(D)), -20, -20,
             rep(-20, length(start) - ncol(D) - 2L))
  upper <- c(rep(30, ncol(D)), 20, 10,
             rep(20, length(start) - ncol(D) - 2L))
  pc_prior <- isTRUE(control_mcml$pc_prior)
  penalty <- NULL
  objective <- engine$evaluate
  if (pc_prior) {
    range_threshold <- if (is.null(control_mcml$pc_range_threshold)) {
      median_nearest_neighbour_distance(distances)
    } else {
      control_mcml$pc_range_threshold / distance_scale
    }
    penalty <- make_pc_penalty(
      engine$layout, range_threshold = range_threshold,
      sigma_threshold = control_mcml$pc_sigma_threshold,
      sigma_probability = control_mcml$pc_sigma_probability,
      range_probability = control_mcml$pc_range_probability
    )
    objective <- penalise_engine(engine$evaluate, penalty)
  }

  laplace <- NULL
  start <- unname(start)
  if (isTRUE(control_mcml$laplace_start)) {
    laplace <- maximise_laplace(
      function(parameters, z) {
        objective(parameters, z, proposal_type = "laplace")
      },
      start, engine$latent_dimension, lower, upper
    )
    if (!is.null(laplace)) start <- laplace$estimate
  }
  fit <- fit_adaptive_mcml(
    engine = engine, start = start, objective = objective,
    penalty = penalty,
    dispersed = if (isTRUE(control_mcml$optimisation_check %||% TRUE)) {
      dispersed_starts(engine$layout, distances, distance_scale, messages)
    },
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
  small_scale <- small_spatial_scale(fit, engine$layout, distances, kappa,
                                     distance_scale)
  report_mcml_outcome(fit, boundary_names, engine$layout, pc_prior,
                      control_mcml$tolerance, messages, small_scale)
  out <- assemble_mcml_result(fit, laplace, boundary_names, D, fix_tau2,
                              re_names, family, invlink, length(y), pc_prior,
                              fit_seed)
  optimizer <- attr(out, "optimizer")
  optimizer$mcml_small_scale <- !is.null(small_scale)
  optimizer$mcml_optimisation_check <- fit$optimisation_check
  optimizer$mcml_better_maximum <- fit$better_maximum
  attr(out, "optimizer") <- optimizer
  out
}

##' A fitted scale smaller than the spacing between locations
##'
##' Returns the fitted scale phi, its practical range and the median
##' nearest-neighbour distance, in the model's distance units, when phi is
##' below that distance and the spatial variance is not at zero, and NULL
##' otherwise.
##' @noRd
small_spatial_scale <- function(fit, layout, distances, kappa,
                                distance_scale) {
  if (layout$index$sigma2 %in% fit$boundary_parameters) return(NULL)
  spacing <- median_nearest_neighbour_distance(distances)
  scale <- exp(fit$estimate[layout$index$phi])
  if (!(scale < spacing)) return(NULL)
  c(scale = scale, practical_range = scale * practical_range_factor(kappa),
    spacing = spacing) * distance_scale
}

##' Practical range of a Matern correlation in units of its scale
##'
##' The distance, as a multiple of the scale phi, at which the Matern
##' correlation with smoothness `kappa` falls to `correlation` (about three
##' scales for the exponential correlation).
##' @noRd
practical_range_factor <- function(kappa, correlation = 0.05) {
  matern <- function(x) {
    2^(1 - kappa) / gamma(kappa) * x^kappa * besselK(x, kappa)
  }
  stats::uniroot(function(x) matern(x) - correlation, c(1e-6, 100),
                 tol = 1e-10)$root
}

##' Build the fitted object returned by `glgpm_nong_adaptive()`
##' @noRd
assemble_mcml_result <- function(fit, laplace, boundary_names, D, fix_tau2,
                                 re_names, family, invlink, n, pc_prior,
                                 fit_seed) {
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
    link_function = make_inverse_link_functions(family, invlink, n)
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
    seed = fit_seed,
    elapsed_seconds = fit$elapsed_seconds
  )
  class(out) <- "RiskMap"
  out
}

##' Tell the user how a certified MCML fit ended
##'
##' Warnings name variance components at zero, a scale smaller than the
##' location spacing, a singular information matrix, or a certificate above
##' the tolerance. Weak identification of the spatial variance or scale is what
##' the PC penalty addresses, so these warnings point to it when it is off.
##' @noRd
report_mcml_outcome <- function(fit, boundary_names, layout, pc_prior,
                                tolerance, messages, small_scale = NULL) {
  pc_advice <- if (pc_prior) "" else paste0(
    " Alternatively, refit with 'control_mcml = set_control_mcml(pc_prior = ",
    "TRUE)', which stabilises the spatial variance and scale with a mild ",
    "penalty; this is recommended for small or low-prevalence surveys and ",
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
      "scale cannot be distinguished from unstructured variation between ",
      "locations, so these data say little about the scale and its standard ",
      "error may be very large. The estimate is the maximum likelihood ",
      "estimate. Consider a model with a nugget (unstructured) effect, ",
      "'gp(nugget = TRUE)', or fixing the scale from external knowledge.",
      pc_advice,
      call. = FALSE
    )
  }
  if (isTRUE(fit$better_maximum)) {
    check <- fit$optimisation_check
    gain <- max(check$log_likelihood) - check$log_likelihood[[1L]]
    warning(
      "The optimisation check found a higher maximum of the likelihood than ",
      "the one reached from the starting values, by ", signif(gain, 2),
      " log-likelihood units on the same Monte Carlo draws, and the estimate ",
      "is taken from it. The likelihood has more than one maximum, so Wald ",
      "intervals may describe the uncertainty poorly. Consider inspecting ",
      "the profile likelihood of the spatial scale and variance.",
      call. = FALSE
    )
  }
  draws <- paste0(format(fit$budget, big.mark = ","), " fitting and ",
                  format(fit$validation_budget, big.mark = ","),
                  " validation draws")
  if (identical(fit$outcome, "weak identification")) {
    warning(
      "The observed information matrix is singular, so some parameters ",
      "(often the spatial scale when spatial correlation is weak) are not ",
      "identifiable from these data. The estimate is returned without ",
      "standard errors and its Monte Carlo precision is not established; ",
      "consider a simpler model or fixing the scale.",
      pc_advice,
      call. = FALSE
    )
  } else if (!fit$certified) {
    warning(
      "The Monte Carlo precision of the estimate could not be established ",
      "with the largest permitted number of draws (", draws, "): the ",
      "estimate may be up to ", signif(fit$certificate, 2), " standard ",
      "errors from the exact maximum likelihood estimate, above the ",
      "tolerance of ", tolerance, ". Increase 'max_samples' or 'max_stages' ",
      "in set_control_mcml(), or relax 'tolerance'.",
      call. = FALSE
    )
  } else if (messages) {
    message(
      "Monte Carlo maximum likelihood certified the estimate to within ",
      signif(fit$certificate, 2), " standard errors of the exact maximum ",
      "likelihood estimate (tolerance ", tolerance, "), after ",
      length(fit$history), " stage", if (length(fit$history) == 1L) "" else
        "s", " with ", draws, "."
    )
  }
  invisible(NULL)
}
