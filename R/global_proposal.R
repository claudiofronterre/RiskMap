##' Gauss-Hermite nodes and weights for a standard normal distribution
##'
##' Golub-Welsch: the nodes are the eigenvalues of the Jacobi matrix of the
##' probabilists' Hermite polynomials and the weights the squared first
##' components of its eigenvectors, so the weights sum to one.
##' @noRd
build_normal_quadrature <- function(number_of_nodes = 20L) {
  index <- seq_len(number_of_nodes - 1L)
  jacobi <- matrix(0, number_of_nodes, number_of_nodes)
  jacobi[cbind(index, index + 1L)] <- sqrt(index)
  jacobi[cbind(index + 1L, index)] <- sqrt(index)
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = decomposition$values, weights = decomposition$vectors[1L, ]^2)
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
  right_side <- function(power) {
    rowSums(weights * nonlinear *
              matrix(nodes^power, nrow(weights), length(nodes), byrow = TRUE))
  }
  m00 <- moment(0); m01 <- moment(1); m02 <- -0.5 * moment(2)
  m11 <- moment(2); m12 <- -0.5 * moment(3); m22 <- 0.25 * moment(4)
  t0 <- right_side(0); t1 <- right_side(1); t2 <- -0.5 * right_side(2)
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
  for (observation in poorly_conditioned) {
    coefficients[observation, ] <- qr.coef(
      qr(root_weights[observation, ] * basis),
      root_weights[observation, ] * nonlinear[observation, ]
    )
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
##' b = l'(m) + a1 / sd + C m. `family` is the family object, whose
##' `nonlinear_parts()` gives the part of each likelihood that is not linear in
##' the step, with its slope.
##' @noRd
update_observation_factors <- function(marginal_mean, marginal_sd, offset,
                                       y, units_m, family,
                                       current_linear, current_curvature,
                                       quadrature) {
  kernel <- check_family_object(family)
  nodes <- quadrature$nodes
  step <- outer(marginal_sd, nodes)
  parts <- kernel$nonlinear_parts(offset + marginal_mean, step, y, units_m)
  # Log ratio of the likelihood to its current approximation at each node,
  # up to a constant per observation.
  log_ratio <- parts$nonlinear +
    (marginal_sd * (parts$slope - current_linear +
                      current_curvature * marginal_mean)) %o% nodes +
    (0.5 * current_curvature * marginal_sd^2) %o% nodes^2
  log_ratio[!is.finite(log_ratio)] <- -Inf
  largest <- log_ratio[, 1L]
  for (node in seq_along(nodes)[-1L]) {
    largest <- pmax(largest, log_ratio[, node])
  }
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
##' under that proposal. The marginal means and variances come from compiled
##' code (`cpp_global_marginals()`), which keeps the dense factorisation and
##' inverse off the R heap. Returns the refitted factors, or NULL if the
##' precision is not positive definite or a fit fails, which includes a
##' quadratic fit with too few weighted nodes.
##' @noRd
refit_global_factors <- function(factors, mu, y, units_m, latent_map,
                                 precision, family, quadrature,
                                 floor_curvature) {
  kernel <- check_family_object(family)
  observations <- length(factors) / 2
  linear <- factors[seq_len(observations)]
  curvature <- factors[observations + seq_len(observations)]
  marginals <- cpp_global_marginals(precision, linear, curvature,
                                    latent_map$columns)
  if (is.null(marginals)) return(NULL)
  marginal_mean <- marginals$marginal_mean
  marginal_variance <- marginals$marginal_variance
  if (any(!is.finite(marginal_variance)) || any(marginal_variance <= 0)) {
    return(NULL)
  }
  updated <- tryCatch(
    update_observation_factors(marginal_mean, sqrt(marginal_variance), mu,
                               y, units_m, kernel, linear, curvature,
                               quadrature),
    error = function(error) NULL
  )
  if (is.null(updated)) return(NULL)
  # A quadratic fit with too few weighted nodes gives an NA curvature, which
  # fails the refit for every link. An infinite curvature fails it unless the
  # kernel is concave (the canonical link without integrated records), whose
  # -Inf is left to the floor.
  if (anyNA(updated$curvature) ||
      (!kernel$assumes_concave && any(!is.finite(updated$curvature)))) {
    return(NULL)
  }
  refitted <- c(updated$linear, floor_curvature(updated$curvature))
  if (any(!is.finite(refitted))) return(NULL)
  refitted
}

##' A lower bound on the factors' curvature that keeps Q positive definite
##'
##' For C >= c < 0, A' diag(C) A >= c A'A, and by Gershgorin the largest
##' eigenvalue of A'A is at most the number of latent columns per observation
##' times the largest number of observations sharing a latent column. Returns a
##' function that floors a vector of curvatures at that bound, which it computes
##' only if a curvature is negative.
##' @noRd
make_curvature_floor <- function(latent_map, precision) {
  curvature_floor <- NULL
  function(curvature) {
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
}

##' The floor of the factors' curvature in the global proposal
##'
##' The smooth floor at zero (`smooth_curvature_floor()`) for integrated
##' records and for a link whose curvature can be negative, which keeps the
##' proposal a smooth function of the parameters, as the theory of the
##' certificate assumes, where a floor at a bound has a kink. Otherwise the
##' floor at the Gershgorin bound (`make_curvature_floor()`), as before step 5b
##' of #253, since under a link of the family's `concave_links` a curvature is
##' negative only through rounding (`make_family()`, property `smooth_floor`).
##' @noRd
choose_curvature_floor <- function(kernel, latent_map, precision) {
  if (kernel$smooth_floor) {
    smooth_curvature_floor
  } else {
    make_curvature_floor(latent_map, precision)
  }
}

##' The global proposal from a fixed number of refits of an anchor
##'
##' At a parameter value other than the anchor's, the anchor's factors, carried
##' over to the new regression part `mu` as b - C (mu - mu_anchor) and C, are
##' refitted `count` times by plain refits (`refit_global_factors()`), with no
##' acceleration, tolerance or restart. Given the anchor, the proposal is then
##' a smooth function of the parameters, as the theory of the certificate
##' assumes, whereas an iteration stopped by a tolerance jumps wherever the
##' number of refits it takes changes. Returns the centre, precision, root and
##' factors as `fit_global_proposal()` does, with the number of `refits`, the
##' largest `change` of a factor at the last refit, and `failed = TRUE` (with
##' no proposal) if a refit fails or the precision is not positive definite.
##' @noRd
refit_from_anchor <- function(anchor, mu, y, units_m, latent_map, precision,
                              family, count = 3L,
                              quadrature = build_normal_quadrature()) {
  kernel <- check_family_object(family)
  floor_curvature <- choose_curvature_floor(kernel, latent_map, precision)
  factors <- c(anchor$linear - anchor$curvature * (mu - anchor$mu),
               anchor$curvature)
  change <- NA_real_
  for (refit in seq_len(count)) {
    refitted <- refit_global_factors(factors, mu, y, units_m, latent_map,
                                     precision, kernel, quadrature,
                                     floor_curvature)
    if (is.null(refitted)) {
      return(list(failed = TRUE, refits = refit, change = change))
    }
    change <- max(abs(refitted - factors))
    factors <- refitted
  }
  observations <- length(factors) / 2
  linear <- factors[seq_len(observations)]
  curvature <- factors[observations + seq_len(observations)]
  final <- cpp_global_proposal_factors(precision, linear, curvature,
                                       latent_map$columns)
  if (is.null(final)) {
    return(list(failed = TRUE, refits = count, change = change))
  }
  list(centre = final$centre, precision = final$precision, root = final$root,
       factors = list(linear = linear, curvature = curvature),
       refits = count, change = change, failed = FALSE)
}

##' Anderson acceleration step for a fixed-point iteration
##'
##' Given the residuals f_j = G(x_j) - x_j and refits G(x_j) of the recent
##' iterates, one column each with the oldest first, returns the combination
##' of refits whose combined residual is smallest in least squares (Walker and
##' Ni, 2011). With a single column, or collinear residual differences, it is
##' the plain step G(x_k).
##' @noRd
take_anderson_step <- function(residuals, refits) {
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
##' by the plain step. The iteration starts from the conditional mode and stops
##' when one refit changes no factor by more than `tolerance`, as plain
##' iteration does. The engine uses it for the anchor of a stage, and the
##' other proposals of the stage refit the anchor a fixed number of times
##' (`refit_from_anchor()`). Returns NULL if a refit fails or the precision is
##' not positive definite. An iteration that has not converged still defines a
##' valid proposal. It is returned with
##' `converged = FALSE` so that the caller can decide whether to fall back to
##' the Laplace proposal. The returned `root` is the upper factor of Q^-1 that
##' maps standard-normal draws (`inverse_cholesky_factor()`), or NULL if that
##' factorisation fails.
##' @noRd
fit_global_proposal <- function(mode, mu, y, units_m, latent_map, precision,
                                family,
                                tolerance = 1e-8, maximum_iterations = 500L,
                                quadrature = build_normal_quadrature(),
                                memory = 3L) {
  kernel <- check_family_object(family)
  floor_curvature <- choose_curvature_floor(kernel, latent_map, precision)
  mode_terms <- kernel$terms(
    mu + latent_linear_predictor(latent_map, mode), y, units_m
  )
  if (is.null(mode_terms)) return(NULL)
  # A link that is not log-concave can have negative curvature at the mode.
  curvature <- if (kernel$assumes_concave) {
    mode_terms$curvature
  } else {
    floor_curvature(mode_terms$curvature)
  }
  linear <- mode_terms$score +
    curvature * latent_linear_predictor(latent_map, mode)
  observations <- length(linear)
  curvature_part <- observations + seq_len(observations)
  refit <- function(factors) {
    refit_global_factors(factors, mu, y, units_m, latent_map, precision,
                         kernel, quadrature, floor_curvature)
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
    candidate <- take_anderson_step(residuals, history)
    candidate[curvature_part] <- floor_curvature(candidate[curvature_part])
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
  linear <- refitted[seq_len(observations)]
  curvature <- refitted[curvature_part]
  final <- cpp_global_proposal_factors(precision, linear, curvature,
                                       latent_map$columns)
  if (is.null(final)) return(NULL)
  list(centre = final$centre, precision = final$precision, root = final$root,
       iterations = refits, change = change, converged = change < tolerance,
       factors = list(linear = linear, curvature = curvature))
}
