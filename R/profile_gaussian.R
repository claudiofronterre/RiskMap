##' The exact likelihood of a Gaussian fit, rebuilt for profiling
##'
##' The Gaussian counterpart of `make_profile_problem()`, rebuilt from the data
##' and settings kept in the fitted object. The likelihood is exact, so there
##' is no bank of draws. It is the log-likelihood that [glgpm()] maximises,
##' without the constant. The parameters are those of the fit on the user's
##' working scale (coefficients, then the logarithms of the variances and the
##' scale), except that an estimated nugget enters as the logarithm of its
##' variance, `tau2`, as `summary()` reports it, and not of its ratio to the
##' spatial variance. The optimiser works on the scale of the fit, with
##' standardised covariates and distances rescaled by the largest distance.
##'
##' For given covariance parameters, the coefficients that maximise the
##' likelihood, subject to any fixed coefficients, are found by generalised
##' least squares (`internal$concentrate`), so the optimiser searches only
##' over the covariance parameters.
##' @noRd
make_gaussian_profile_problem <- function(fitted) {
  design_scaling <- standardize_design_matrix(fitted$D)
  spatial_scaling <- scale_spatial_coordinates(fitted$coords)
  p <- ncol(fitted$D)
  likelihood <- make_gaussian_likelihood(
    y = fitted$y - (fitted$cov_offset %||% 0), D = design_scaling$design,
    distances = pairwise_distances(spatial_scaling$coordinates),
    kappa = fitted$kappa, ID_coords = fitted$ID_coords,
    ID_re = if (is.null(fitted$ID_re)) NULL else as.matrix(fitted$ID_re),
    fix_tau2 = fitted$fix_tau2, fix_var_me = fitted$fix_var_me
  )
  maps <- profile_scale_maps(design_scaling$coefficient_transform,
                             log(spatial_scaling$distance_scale), p)
  value_internal <- function(internal) {
    result <- tryCatch(likelihood$value(internal),
                       error = function(error) NA_real_)
    if (is.finite(result)) result else NA_real_
  }
  value <- function(working) value_internal(maps$to(working))

  estimate <- unlist(fitted$estimate)
  variances <- diag(fitted$covariance)
  nugget <- match("nu2", names(estimate))
  if (!is.na(nugget)) {
    # log tau2 = log nu2 + log sigma2, with its variance by the delta method,
    # as in summary().
    estimate[nugget] <- estimate[nugget] + estimate[["sigma2"]]
    variances[nugget] <- variances[nugget] + variances[["sigma2"]] +
      2 * fitted$covariance["nu2", "sigma2"]
    names(estimate)[nugget] <- "tau2"
  }
  log_scale <- log(spatial_scaling$distance_scale)
  others <- length(estimate) - p - 2L
  problem <- list(
    value = value, estimate = estimate, names = names(estimate),
    labels = profile_labels(names(estimate), fitted$distance_units),
    standard_errors = stats::setNames(sqrt(pmax(variances, 0)),
                                      names(estimate)),
    # The coefficients are solved in closed form and need no bounds.
    lower = c(rep(-Inf, p), -20, -20 + log_scale, rep(-20, others)),
    upper = c(rep(Inf, p), 20, 10 + log_scale, rep(20, others)),
    internal = list(value = value_internal, to = maps$to, from = maps$from,
                    transform = design_scaling$coefficient_transform, p = p,
                    lower = c(rep(-Inf, p), -20, -20, rep(-20, others)),
                    upper = c(rep(Inf, p), 20, 10, rep(20, others)),
                    concentrate = likelihood$coefficients),
    predict_starts = TRUE, exact = TRUE, max_step = 1
  )
  problem$control <- make_optimiser_control(abs(value(estimate)))
  best <- maximise_profile(problem, integer(), numeric(), estimate)
  problem$bank_estimate <- best$par
  problem$max_loglik <- best$value
  problem
}

##' The exact log-likelihood of a linear geostatistical model
##'
##' The observations have covariance `V = Z G Z' + omega2 I`, where `Z` maps
##' observations to locations and to the levels of the random effects, and
##' `G` is block diagonal: `sigma2 R + tau2 I` over the locations, with `R`
##' the Matern correlation, and `sigma2_re I` for each random effect. For
##' given covariance parameters, the log-likelihood depends on the data
##' through `log det V` and the cross-products of the design and the response
##' weighted by the inverse of `V`, so these are computed once and give both
##' the log-likelihood at any coefficients and the coefficients that maximise
##' it (`coefficients()`). The internal parameters are the standardised
##' coefficients followed by the covariance parameters in the layout of
##' `cov_par_layout()`, with the logarithm of `tau2` in the place of the
##' nugget ratio.
##'
##' As in [glgpm()], `V` is factorised directly unless there are more than
##' three observations per latent effect, when the Woodbury identity works in
##' the dimension of `G` instead (`latent_cross_products()`).
##' @noRd
make_gaussian_likelihood <- function(y, D, distances, kappa, ID_coords, ID_re,
                                     fix_tau2, fix_var_me, woodbury = NULL) {
  p <- ncol(D)
  n_re <- if (is.null(ID_re)) 0L else ncol(ID_re)
  n_locations <- max(ID_coords)
  levels <- if (n_re > 0L) apply(ID_re, 2L, max) else integer()
  layout <- cov_par_layout(p, fix_tau2 = fix_tau2,
                           sigma2_me = is.null(fix_var_me), n_re = n_re)$index
  # Each evaluation fills the covariance matrix of the locations from the
  # distances of its lower triangle, with as few copies as possible: most of
  # the time of a small model goes to collecting discarded matrices.
  distances <- as.numeric(distances)
  lower <- which(lower.tri(diag(n_locations)), arr.ind = TRUE)
  mirror <- (lower[, "row"] - 1L) * n_locations + lower[, "col"]
  lower <- (lower[, "col"] - 1L) * n_locations + lower[, "row"]
  columns <- cbind(D, y)
  coefficients_index <- seq_len(p)
  if (is.null(woodbury)) {
    woodbury <- !use_direct_gaussian_covariance(
      length(y), n_locations + sum(levels), fix_var_me, fix_tau2
    )
  }
  direct <- direct_cross_products(columns, ID_coords, ID_re)
  through_latent <- if (woodbury) {
    latent_cross_products(columns, ID_coords, ID_re, n_locations, levels)
  }

  covariance_parameters <- function(internal) {
    list(sigma2 = exp(internal[layout$sigma2]),
         phi = exp(internal[layout$phi]),
         tau2 = if (isTRUE(fix_tau2)) exp(internal[layout$nu2]) else
           as.numeric(fix_tau2 %||% 0),
         omega2 = if (is.null(fix_var_me)) exp(internal[layout$sigma2_me]) else
           fix_var_me,
         sigma2_re = if (n_re > 0L) exp(internal[layout$sigma2_re]) else
           numeric())
  }

  # The last covariance parameters and their cross-products, since a
  # profile point asks for the coefficients and then the log-likelihood at
  # the same covariance parameters.
  cache <- new.env(parent = emptyenv())
  cross_products <- function(internal) {
    covariance <- internal[-coefficients_index]
    if (identical(covariance, cache$covariance)) return(cache$products)
    parameters <- covariance_parameters(internal)
    covariances <- parameters$sigma2 *
      matern_correlation(distances, phi = parameters$phi, kappa = kappa)
    spatial <- diag(parameters$sigma2 + parameters$tau2, n_locations)
    spatial[lower] <- covariances
    spatial[mirror] <- covariances
    products <- if (woodbury) through_latent(spatial, parameters)
    if (is.null(products)) products <- direct(spatial, parameters)
    cache$covariance <- covariance
    cache$products <- products
    products
  }

  value <- function(internal) {
    products <- cross_products(internal)
    beta <- internal[coefficients_index]
    cross <- products$cross
    quadratic <- cross[p + 1L, p + 1L] -
      2 * sum(beta * cross[coefficients_index, p + 1L]) +
      sum(beta * (cross[coefficients_index, coefficients_index] %*% beta))
    unname(-0.5 * (products$log_determinant + quadratic))
  }

  # Covariance parameters at which V cannot be factorised give missing
  # coefficients, and so a missing log-likelihood.
  coefficients <- function(internal, constraint = NULL, targets = NULL) {
    cross <- tryCatch(cross_products(internal)$cross,
                      error = function(error) NULL)
    if (is.null(cross)) return(rep(NA_real_, p))
    constrained_least_squares(cross, constraint, targets)
  }

  list(value = value, coefficients = coefficients)
}

##' Weighted cross-products from the observed covariance, factorised directly
##'
##' Returns a function of the covariance matrix of the locations (with the
##' nugget) and the covariance parameters, giving `log det V` and the
##' cross-products of `columns` weighted by the inverse of `V`.
##' @noRd
direct_cross_products <- function(columns, ID_coords, ID_re) {
  one_per_location <- identical(as.integer(ID_coords),
                                seq_len(nrow(columns)))
  same_level <- NULL
  function(spatial, parameters) {
    # The indicators of observations at the same level of each random
    # effect, made when first needed, since they are as large as V.
    if (is.null(same_level) && !is.null(ID_re)) {
      same_level <<- lapply(seq_len(ncol(ID_re)), function(effect) {
        1 * outer(ID_re[, effect], ID_re[, effect], "==")
      })
    }
    covariance <- if (one_per_location) spatial else
      spatial[ID_coords, ID_coords]
    for (effect in seq_along(same_level)) {
      covariance <- covariance +
        parameters$sigma2_re[effect] * same_level[[effect]]
    }
    diag(covariance) <- diag(covariance) + parameters$omega2
    # Symmetric by construction, so factorised as it is.
    root <- chol(covariance)
    whitened <- backsolve(root, columns, transpose = TRUE)
    list(log_determinant = log_determinant_from_cholesky(root),
         cross = crossprod(whitened))
  }
}

##' Weighted cross-products through the latent effects, by the Woodbury
##' identity
##'
##' With `G = F'F` and `B = Z F'`, `V = omega2 I + B B'`, so
##' `V^-1 = (I - B M^-1 B') / omega2` and `det V = omega2^(n - q) det M`,
##' where `M = omega2 I + B'B` and `q` is the number of latent effects. This
##' does not need the inverse of `G`. `F` is block diagonal, the Cholesky
##' factor of the covariance of the locations and the standard deviations of
##' the random effects over their levels, so `B'B = F Z'Z F'` and the
##' cross-products `B' columns = F Z' columns` are built block by block.
##' Returns a function like that of `direct_cross_products()`, which gives
##' `NULL` when the covariance of the locations cannot be factorised.
##' @noRd
latent_cross_products <- function(columns, ID_coords, ID_re, n_locations,
                                  levels) {
  n <- nrow(columns)
  latent <- cbind(ID_coords, if (length(levels)) ID_re +
                    rep(n_locations + cumsum(c(0, levels[-length(levels)])),
                        each = n))
  latent_design <- matrix(0, n, n_locations + sum(levels))
  for (effect in seq_len(ncol(latent))) {
    latent_design[cbind(seq_len(n), latent[, effect])] <- 1
  }
  locations <- seq_len(n_locations)
  # Z'Z is diagonal over the locations, since each observation has one.
  counts <- tabulate(ID_coords, n_locations)
  latent_cross <- crossprod(latent_design)
  location_levels <- latent_cross[locations, -locations, drop = FALSE]
  level_levels <- latent_cross[-locations, -locations, drop = FALSE]
  latent_columns <- crossprod(latent_design, columns)
  columns_cross <- crossprod(columns)
  function(spatial, parameters) {
    spatial_root <- tryCatch(chol(spatial), error = function(error) NULL)
    if (is.null(spatial_root)) return(NULL)
    middle <- tcrossprod(spatial_root *
                           rep(sqrt(counts), each = n_locations))
    projected <- spatial_root %*% latent_columns[locations, , drop = FALSE]
    if (length(levels)) {
      deviations <- rep(sqrt(parameters$sigma2_re), levels)
      cross <- (spatial_root %*% location_levels) *
        rep(deviations, each = n_locations)
      middle <- rbind(cbind(middle, cross),
                      cbind(t(cross), level_levels * tcrossprod(deviations)))
      projected <- rbind(projected,
                         latent_columns[-locations, , drop = FALSE] *
                           deviations)
    }
    diag(middle) <- diag(middle) + parameters$omega2
    middle_root <- chol(middle)
    projected <- backsolve(middle_root, projected, transpose = TRUE)
    list(log_determinant = (n - nrow(middle)) * log(parameters$omega2) +
           log_determinant_from_cholesky(middle_root),
         cross = (columns_cross - crossprod(projected)) / parameters$omega2)
  }
}

##' Generalised least squares from weighted cross-products, with some
##' coefficients fixed
##'
##' `cross` holds the cross-products of the design and the response weighted
##' by the inverse covariance, the response last. The coefficients minimise
##' the weighted residual sum of squares subject to
##' `constraint %*% coefficients = targets`: they are the constrained
##' solution of least norm plus the best combination of a basis of the null
##' space of the constraint.
##' @noRd
constrained_least_squares <- function(cross, constraint = NULL,
                                      targets = NULL) {
  p <- nrow(cross) - 1L
  weighted_design <- cross[seq_len(p), seq_len(p), drop = FALSE]
  weighted_response <- cross[seq_len(p), p + 1L]
  if (is.null(constraint)) {
    return(unname(drop(solve(weighted_design, weighted_response))))
  }
  particular <- drop(crossprod(constraint,
                               solve(tcrossprod(constraint), targets)))
  if (nrow(constraint) == p) return(unname(particular))
  null_space <- qr.Q(qr(t(constraint)), complete = TRUE)[
    , -seq_len(nrow(constraint)), drop = FALSE]
  free <- solve(crossprod(null_space, weighted_design %*% null_space),
                crossprod(null_space,
                          weighted_response - weighted_design %*% particular))
  unname(particular + drop(null_space %*% free))
}
