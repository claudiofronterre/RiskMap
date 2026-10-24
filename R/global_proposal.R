##' Gauss-Hermite nodes and weights for a standard normal distribution
##'
##' Golub-Welsch: the nodes are the eigenvalues of the Jacobi matrix of the
##' probabilists' Hermite polynomials and the weights the squared first
##' components of its eigenvectors, so the weights sum to one.
##' @noRd
standard_normal_quadrature <- function(number_of_nodes = 20L) {
  index <- seq_len(number_of_nodes - 1L)
  jacobi <- matrix(0, number_of_nodes, number_of_nodes)
  jacobi[cbind(index, index + 1L)] <- sqrt(index)
  jacobi[cbind(index + 1L, index)] <- sqrt(index)
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = decomposition$values, weights = decomposition$vectors[1L, ]^2)
}

##' Stable nonlinear part of each observation's log-likelihood
##'
##' For an observation with linear predictor `centre + step`, returns
##' l(centre + step) - l(centre) - step * l'(centre), the part of the
##' log-likelihood that is not linear in the step, together with the slope
##' l'(centre). For the canonical binomial and Poisson links the difference is
##' written with log1p and expm1 so that nearly equal log-likelihood values are
##' never subtracted; other links use the direct difference.
##' @param centre Linear predictor at the expansion point, one per observation.
##' @param step Matrix of steps, one row per observation.
##' @noRd
observation_nonlinear_part <- function(centre, step, y, units_m, family,
                                       linkf, canonical) {
  if (canonical && family == "binomial") {
    probability <- stats::plogis(centre)
    slope <- y - units_m * probability
    # log(1 - p + p exp(step)), switching to log-sum-exp for large steps.
    log_ratio <- ifelse(
      abs(step) < 20,
      log1p(probability * expm1(step)),
      pmax(log1p(-probability), log(probability) + step) +
        log1p(exp(-abs(log1p(-probability) - log(probability) - step)))
    )
    nonlinear <- units_m * (probability * step - log_ratio)
    return(list(nonlinear = nonlinear, slope = slope))
  }
  if (canonical && family == "poisson") {
    rate <- units_m * exp(centre)
    slope <- y - rate
    nonlinear <- -rate * (expm1(step) - step)
    return(list(nonlinear = nonlinear, slope = slope))
  }
  log_likelihood <- function(eta) {
    mean <- linkf$inv(eta)
    if (family == "poisson") {
      y * log(pmax(units_m * mean, .Machine$double.xmin)) - units_m * mean
    } else {
      y * log(pmax(mean, .Machine$double.xmin)) +
        (units_m - y) * log(pmax(1 - mean, .Machine$double.xmin))
    }
  }
  terms <- laplace_observation_terms(centre, y, units_m, family, linkf)
  if (is.null(terms)) {
    stop("Invalid linear predictor in the global proposal.", call. = FALSE)
  }
  per_observation_centre <- log_likelihood(centre)
  nonlinear <- matrix(log_likelihood(as.vector(centre + step)), nrow(step)) -
    per_observation_centre - step * terms$score
  list(nonlinear = nonlinear, slope = terms$score)
}

##' Weighted quadratic fit of every observation's log-likelihood
##'
##' For each observation, fits the nonlinear part of its log-likelihood at
##' standardised quadrature nodes z by weighted least squares on the basis
##' 1, z, -z^2/2. The fits are solved together in closed form; any observation
##' whose normal equations are poorly conditioned is refitted by QR, so the
##' closed form is never trusted without a check.
##' @param nonlinear Matrix of residuals, one row per observation.
##' @param root_weights Matrix of square-root weights of the same shape.
##' @param nodes Standardised nodes.
##' @return A matrix with columns `intercept`, `linear` and `curvature`.
##' @noRd
fit_observation_quadratics <- function(nonlinear, root_weights, nodes,
                                       minimum_condition = 1e-8) {
  weights <- root_weights^2
  moment <- function(power) as.numeric(weights %*% nodes^power)
  rhs <- function(power) rowSums(weights * nonlinear *
                                   matrix(nodes^power, nrow(weights),
                                          length(nodes), byrow = TRUE))
  m00 <- moment(0); m01 <- moment(1); m02 <- -0.5 * moment(2)
  m11 <- moment(2); m12 <- -0.5 * moment(3); m22 <- 0.25 * moment(4)
  t0 <- rhs(0); t1 <- rhs(1); t2 <- -0.5 * rhs(2)
  # Cofactors of the symmetric 3 x 3 normal matrix.
  c00 <- m11 * m22 - m12^2
  c01 <- -(m01 * m22 - m12 * m02)
  c02 <- m01 * m12 - m11 * m02
  c11 <- m00 * m22 - m02^2
  c12 <- -(m00 * m12 - m01 * m02)
  c22 <- m00 * m11 - m01^2
  determinant <- m00 * c00 + m01 * c01 + m02 * c02
  # Determinant of the matrix scaled to unit diagonal, a cheap conditioning
  # measure: one for orthogonal columns, zero for collinear ones.
  scaled_determinant <- determinant / (m00 * m11 * m22)
  coefficients <- cbind(
    intercept = (c00 * t0 + c01 * t1 + c02 * t2) / determinant,
    linear = (c01 * t0 + c11 * t1 + c12 * t2) / determinant,
    curvature = (c02 * t0 + c12 * t1 + c22 * t2) / determinant
  )
  poorly_conditioned <- which(!is.finite(scaled_determinant) |
                                scaled_determinant < minimum_condition |
                                rowSums(!is.finite(coefficients)) > 0)
  basis <- cbind(1, nodes, -0.5 * nodes^2)
  for (i in poorly_conditioned) {
    coefficients[i, ] <- qr.coef(qr(root_weights[i, ] * basis),
                                 root_weights[i, ] * nonlinear[i, ])
  }
  coefficients
}

##' Update the per-observation Gaussian factors of the global proposal
##'
##' Each observation's likelihood, as a function of the latent part s of its
##' linear predictor, is approximated by exp(b s - C s^2 / 2), fitted where the
##' current proposal puts its mass: at quadrature nodes of the proposal's
##' marginal distribution of s, weighted by the ratio of the observation's
##' likelihood to its current approximation. In standardised coordinates
##' s = m + sd z the fitted coefficients (a0, a1, a2) give C = a2 / sd^2 and
##' b = l'(m) + a1 / sd + C m.
##' @noRd
update_observation_factors <- function(marginal_mean, marginal_sd, offset,
                                       y, units_m, family, linkf, canonical,
                                       linear_factor, curvature_factor,
                                       quadrature) {
  nodes <- quadrature$nodes
  step <- outer(marginal_sd, nodes)
  parts <- observation_nonlinear_part(offset + marginal_mean, step, y,
                                      units_m, family, linkf, canonical)
  # Log ratio of the likelihood to its current approximation at each node,
  # up to a constant per observation.
  log_ratio <- parts$nonlinear +
    (marginal_sd * (parts$slope - linear_factor +
                      curvature_factor * marginal_mean)) %o% nodes +
    (0.5 * curvature_factor * marginal_sd^2) %o% nodes^2
  log_ratio[!is.finite(log_ratio)] <- -Inf
  largest <- log_ratio[, 1L]
  for (k in seq_along(nodes)[-1L]) largest <- pmax(largest, log_ratio[, k])
  relative <- exp(log_ratio - largest)
  root_weights <- sqrt(sweep(relative, 2L, quadrature$weights, "*"))
  nonlinear <- parts$nonlinear
  nonlinear[root_weights == 0 | !is.finite(nonlinear)] <- 0
  coefficients <- fit_observation_quadratics(nonlinear, root_weights, nodes)
  curvature <- coefficients[, "curvature"] / marginal_sd^2
  linear <- parts$slope + coefficients[, "linear"] / marginal_sd +
    curvature * marginal_mean
  list(linear = linear, curvature = curvature)
}

##' One refit of the factors of the global proposal
##'
##' The fixed-point map of `fit_global_proposal()`. It builds the proposal
##' from the current factors, stacked as c(linear, curvature), and refits each
##' observation's factor at the quadrature nodes of its marginal distribution
##' under that proposal. Returns the refitted factors, or NULL if the
##' precision is not positive definite or a fit fails.
##' @noRd
refit_global_factors <- function(factors, mu, y, units_m, latent_map,
                                 precision, family, linkf, canonical,
                                 quadrature, floor_curvature) {
  n <- length(factors) / 2
  linear <- factors[seq_len(n)]
  curvature <- factors[n + seq_len(n)]
  proposal_precision <- precision + latent_data_hessian(latent_map, curvature)
  root <- tryCatch(chol(proposal_precision), error = function(error) NULL)
  if (is.null(root)) return(NULL)
  inverse <- chol2inv(root)
  centre <- as.numeric(inverse %*% latent_gradient(latent_map, linear))
  marginal_mean <- latent_linear_predictor(latent_map, centre)
  marginal_variance <- as.numeric(rowsum(
    inverse[latent_map$pair_linear], latent_map$pair_observation,
    reorder = FALSE
  ))
  if (any(!is.finite(marginal_variance)) || any(marginal_variance <= 0)) {
    return(NULL)
  }
  updated <- tryCatch(
    update_observation_factors(marginal_mean, sqrt(marginal_variance), mu,
                               y, units_m, family, linkf, canonical,
                               linear, curvature, quadrature),
    error = function(error) NULL
  )
  if (is.null(updated)) return(NULL)
  refitted <- c(updated$linear, floor_curvature(updated$curvature))
  if (any(!is.finite(refitted))) return(NULL)
  refitted
}

##' Anderson acceleration step for a fixed-point iteration
##'
##' Given the residuals f_j = G(x_j) - x_j and refits G(x_j) of the recent
##' iterates, one column each with the oldest first, returns the combination
##' of refits whose combined residual is smallest in least squares (Walker and
##' Ni, 2011). With a single column, or collinear residual differences, it is
##' the plain step G(x_k).
##' @noRd
anderson_step <- function(residuals, refits) {
  last <- ncol(residuals)
  if (last < 2L) return(refits[, last])
  residual_differences <- residuals[, -1L, drop = FALSE] -
    residuals[, -last, drop = FALSE]
  refit_differences <- refits[, -1L, drop = FALSE] -
    refits[, -last, drop = FALSE]
  decomposition <- qr(residual_differences)
  if (decomposition$rank < ncol(residual_differences)) return(refits[, last])
  weights <- qr.coef(decomposition, residuals[, last])
  as.numeric(refits[, last] - refit_differences %*% weights)
}

##' Fit the global Gaussian proposal for the latent effects
##'
##' Starting from the Laplace approximation at the conditional mode, the
##' per-observation factors exp(b s - C s^2 / 2) are refitted until they stop
##' changing (efficient importance sampling with numerical integration,
##' Richard and Zhang, 2007; Koopman, Lucas and Scharth, 2015). The proposal
##' is N(centre, Q^-1) with Q = P + A' diag(C) A and centre = Q^-1 A' b, where
##' P is the prior precision and A maps latent effects to observations. The
##' factors are fitted in standardised coordinates by least squares on the
##' stable nonlinear part of each log-likelihood, which keeps the fit well
##' conditioned when a variance is near zero.
##'
##' The refits are combined by Anderson acceleration with `memory` previous
##' iterates, which reaches the fixed point in far fewer refits than plain
##' iteration. A combined step that fails or enlarges the residual is replaced
##' by the plain step. The iteration always starts from the conditional mode
##' and stops when one refit changes no factor by more than `tolerance`, as
##' plain iteration does, so the proposal is a deterministic function of the
##' parameters. Returns NULL if a refit fails or the precision is not positive
##' definite. An iteration that has not converged still defines a valid
##' proposal; it is returned with `converged = FALSE` so that the caller can
##' decide whether to fall back to the Laplace proposal. With `initial`, a
##' list of starting factors `linear` and `curvature`, the iteration starts
##' from them instead of the conditional mode, which is then not needed.
##' @noRd
fit_global_proposal <- function(mode, mu, y, units_m, latent_map, precision,
                                family, linkf, canonical,
                                tolerance = 1e-8, maximum_iterations = 500L,
                                quadrature = standard_normal_quadrature(),
                                memory = 3L, initial = NULL) {
  if (is.null(initial)) {
    mode_terms <- laplace_observation_terms(
      mu + latent_linear_predictor(latent_map, mode), y, units_m, family,
      linkf
    )
    if (is.null(mode_terms)) return(NULL)
    curvature <- mode_terms$curvature
    linear <- mode_terms$score +
      curvature * latent_linear_predictor(latent_map, mode)
  } else {
    linear <- initial$linear
    curvature <- initial$curvature
  }
  n <- length(linear)
  # A lower bound on the factors' curvature that keeps Q positive definite.
  # For C >= c < 0, A' diag(C) A >= c A'A, and by Gershgorin the largest
  # eigenvalue of A'A is at most the number of latent columns per observation
  # times the largest number of observations sharing a latent column. It is
  # computed only if a fitted curvature is negative.
  curvature_floor <- NULL
  floor_curvature <- function(curvature) {
    if (all(curvature >= 0)) return(curvature)
    if (is.null(curvature_floor)) {
      columns <- latent_map$columns
      largest_column_sum <- ncol(columns) *
        max(tabulate(columns, nbins = latent_map$dimension))
      smallest_prior <- min(eigen(precision, symmetric = TRUE,
                                  only.values = TRUE)$values)
      curvature_floor <<- -0.5 * smallest_prior / largest_column_sum
    }
    pmax(curvature, curvature_floor)
  }
  refit <- function(factors) {
    refit_global_factors(factors, mu, y, units_m, latent_map, precision,
                         family, linkf, canonical, quadrature,
                         floor_curvature)
  }
  factors <- c(linear, curvature)
  refitted <- refit(factors)
  if (is.null(refitted)) return(NULL)
  refits <- 1L
  residuals <- history <- NULL
  change <- Inf
  for (iteration in seq_len(maximum_iterations)) {
    change <- max(abs(refitted - factors))
    if (change < tolerance) break
    stored <- if (is.null(residuals)) 1L else ncol(residuals) + 1L
    keep <- max(1L, stored - memory):stored
    residuals <- cbind(residuals, refitted - factors)[, keep, drop = FALSE]
    history <- cbind(history, refitted)[, keep, drop = FALSE]
    candidate <- anderson_step(residuals, history)
    candidate[n + seq_len(n)] <- floor_curvature(candidate[n + seq_len(n)])
    next_refit <- refit(candidate)
    refits <- refits + 1L
    if (is.null(next_refit) ||
        max(abs(next_refit - candidate)) > 2 * change) {
      residuals <- history <- NULL
      candidate <- refitted
      next_refit <- refit(candidate)
      refits <- refits + 1L
      if (is.null(next_refit)) return(NULL)
    }
    factors <- candidate
    refitted <- next_refit
  }
  linear <- refitted[seq_len(n)]
  curvature <- refitted[n + seq_len(n)]
  proposal_precision <- precision + latent_data_hessian(latent_map, curvature)
  proposal_precision <- (proposal_precision + t(proposal_precision)) / 2
  root <- tryCatch(chol(proposal_precision), error = function(error) NULL)
  if (is.null(root)) return(NULL)
  centre <- as.numeric(backsolve(root, forwardsolve(
    t(root), latent_gradient(latent_map, linear)
  )))
  list(centre = centre, precision = proposal_precision, iterations = refits,
       change = change, converged = change < tolerance,
       factors = list(linear = linear, curvature = curvature))
}
