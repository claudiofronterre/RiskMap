##' Fit a Gaussian model through its observed covariance
##'
##' This path is used where the Woodbury implementation is either undefined
##' (zero measurement error) or has incorrect fixed-nugget derivatives.
##'
##' @noRd
use_direct_gaussian_covariance <- function(n_observations, latent_dimension,
                                           fix_var_me, fix_tau2) {
  correctness_requires_direct <- identical(fix_var_me, 0) ||
    (is.numeric(fix_tau2) && length(fix_tau2) == 1L && fix_tau2 > 0)

  # Benchmarks show that direct factorisation is faster with light replication,
  # while Woodbury is substantially faster once observations greatly outnumber
  # latent effects.
  correctness_requires_direct || n_observations <= 3 * latent_dimension
}

##' @noRd
glgpm_lm_direct <- function(y, D, coords, kappa, ID_coords, ID_re,
                            fix_var_me, fix_tau2, start_beta,
                            start_cov_pars, messages) {
  n <- length(y)
  p <- ncol(D)
  n_loc <- nrow(coords)
  n_re <- if (is.null(ID_re)) 0L else ncol(ID_re)

  spatial_design <- matrix(0, n, n_loc)
  spatial_design[cbind(seq_len(n), ID_coords)] <- 1
  random_designs <- vector("list", n_re)
  if (n_re > 0L) {
    for (j in seq_len(n_re)) {
      n_levels <- max(ID_re[, j])
      random_designs[[j]] <- matrix(0, n, n_levels)
      random_designs[[j]][cbind(seq_len(n), ID_re[, j])] <- 1
    }
  }

  # With no independent measurement error, repeated identical latent-design
  # rows make the observed covariance structurally singular for every
  # parameter value.
  if (identical(fix_var_me, 0)) {
    latent_design <- do.call(cbind, c(list(spatial_design), random_designs))
    if (qr(latent_design)$rank < n) {
      stop(
        "The model with 'fix_var_me = 0' has a structurally singular ",
        "observed covariance. Repeated observations require a positive ",
        "measurement-error variance or additional independent effects.",
        call. = FALSE
      )
    }
  }

  ind_beta <- seq_len(p)
  ind_sigma2 <- p + 1L
  ind_phi <- p + 2L
  next_index <- ind_phi
  if (isTRUE(fix_tau2)) {
    ind_nu2 <- next_index <- next_index + 1L
  }
  if (is.null(fix_var_me)) {
    ind_omega2 <- next_index <- next_index + 1L
  }
  if (n_re > 0L) {
    ind_sigma2_re <- next_index + seq_len(n_re)
  }

  covariance_start <- start_cov_pars
  start_index <- 2L
  start <- c(start_beta, log(covariance_start[1:2]))
  if (isTRUE(fix_tau2)) {
    start_index <- start_index + 1L
    start <- c(start, log(covariance_start[start_index] /
                            covariance_start[1L]))
  }
  if (is.null(fix_var_me)) {
    start_index <- start_index + 1L
    start <- c(start, log(covariance_start[start_index]))
  }
  if (n_re > 0L) {
    re_start <- covariance_start[start_index + seq_len(n_re)]
    start <- c(start, log(re_start))
  }

  distances <- pairwise_distances(coords)
  identity_n <- diag(n)
  spatial_cross <- tcrossprod(spatial_design)
  random_cross <- lapply(random_designs, tcrossprod)
  cache <- new.env(parent = emptyenv())

  covariance_state <- function(par) {
    if (exists("par", cache, inherits = FALSE) && identical(par, cache$par)) {
      return(cache$state)
    }

    sigma2 <- exp(par[ind_sigma2])
    phi <- exp(par[ind_phi])
    correlation <- matern_correlation(
      distances, phi = phi, kappa = kappa, return_sym_matrix = TRUE
    )
    correlation_phi <- matern_gradient_phi(distances, phi, kappa)
    correlation_phi2 <- matern_hessian_phi(distances, phi, kappa)

    spatial_covariance <- spatial_design %*% correlation %*%
      t(spatial_design)
    d_sigma2 <- sigma2 * spatial_covariance
    d_phi <- sigma2 * phi *
      (spatial_design %*% correlation_phi %*% t(spatial_design))
    d2_phi <- sigma2 * (
      spatial_design %*%
        (phi^2 * correlation_phi2 + phi * correlation_phi) %*%
        t(spatial_design)
    )
    covariance <- d_sigma2

    first <- vector("list", length(par))
    second <- vector("list", length(par))
    first[[ind_sigma2]] <- d_sigma2
    first[[ind_phi]] <- d_phi
    second[[ind_sigma2]] <- d_sigma2
    second[[ind_phi]] <- d2_phi

    cross_second <- list()
    cross_second[[paste(ind_sigma2, ind_phi, sep = ":")]] <- d_phi

    if (isTRUE(fix_tau2)) {
      nu2 <- exp(par[ind_nu2])
      d_nu2 <- sigma2 * nu2 * spatial_cross
      covariance <- covariance + d_nu2
      first[[ind_nu2]] <- d_nu2
      second[[ind_nu2]] <- d_nu2
      first[[ind_sigma2]] <- first[[ind_sigma2]] + d_nu2
      second[[ind_sigma2]] <- second[[ind_sigma2]] + d_nu2
      cross_second[[paste(ind_sigma2, ind_nu2, sep = ":")]] <- d_nu2
    } else if (is.numeric(fix_tau2) && fix_tau2 > 0) {
      covariance <- covariance + fix_tau2 * spatial_cross
    }

    if (is.null(fix_var_me)) {
      omega2 <- exp(par[ind_omega2])
      d_omega2 <- omega2 * identity_n
      covariance <- covariance + d_omega2
      first[[ind_omega2]] <- d_omega2
      second[[ind_omega2]] <- d_omega2
    } else if (fix_var_me > 0) {
      covariance <- covariance + fix_var_me * identity_n
    }

    if (n_re > 0L) {
      for (j in seq_len(n_re)) {
        variance <- exp(par[ind_sigma2_re[j]])
        derivative <- variance * random_cross[[j]]
        covariance <- covariance + derivative
        first[[ind_sigma2_re[j]]] <- derivative
        second[[ind_sigma2_re[j]]] <- derivative
      }
    }

    root <- factor_covariance(
      covariance, "Gaussian observed covariance", allow_jitter = FALSE
    )
    state <- list(
      root = root,
      precision = chol2inv(root),
      log_determinant = log_determinant_from_cholesky(root),
      first = first,
      second = second,
      cross_second = cross_second
    )
    cache$par <- par
    cache$state <- state
    state
  }

  log_likelihood <- function(par) {
    state <- covariance_state(par)
    residual <- y - as.numeric(D %*% par[ind_beta])
    alpha <- state$precision %*% residual
    -0.5 * (state$log_determinant + crossprod(residual, alpha))
  }

  gradient <- function(par) {
    state <- covariance_state(par)
    residual <- y - as.numeric(D %*% par[ind_beta])
    alpha <- as.numeric(state$precision %*% residual)
    out <- numeric(length(par))
    out[ind_beta] <- crossprod(D, alpha)
    covariance_indices <- setdiff(seq_along(par), ind_beta)
    for (j in covariance_indices) {
      derivative <- state$first[[j]]
      out[j] <- 0.5 * (crossprod(alpha, derivative %*% alpha) -
                         sum(state$precision * derivative))
    }
    out
  }

  hessian <- function(par) {
    state <- covariance_state(par)
    precision <- state$precision
    residual <- y - as.numeric(D %*% par[ind_beta])
    alpha <- as.numeric(precision %*% residual)
    out <- matrix(0, length(par), length(par))
    out[ind_beta, ind_beta] <- -crossprod(D, precision %*% D)
    covariance_indices <- setdiff(seq_along(par), ind_beta)

    for (j in covariance_indices) {
      derivative_j <- state$first[[j]]
      beta_cross <- -crossprod(D, precision %*% derivative_j %*% alpha)
      out[ind_beta, j] <- beta_cross
      out[j, ind_beta] <- beta_cross
      for (k in covariance_indices[covariance_indices >= j]) {
        derivative_k <- state$first[[k]]
        key <- paste(min(j, k), max(j, k), sep = ":")
        derivative_jk <- if (j == k) state$second[[j]] else
          state$cross_second[[key]]
        if (is.null(derivative_jk)) derivative_jk <- matrix(0, n, n)
        value <- 0.5 * (
          sum((precision %*% derivative_k %*% precision) * derivative_j) -
            sum(precision * derivative_jk) -
            2 * crossprod(alpha, derivative_j %*% precision %*%
                            derivative_k %*% alpha) +
            crossprod(alpha, derivative_jk %*% alpha)
        )
        out[j, k] <- out[k, j] <- value
      }
    }
    out
  }

  objective <- safe_optimizer_objective(function(par) -log_likelihood(par))
  estimate <- nlminb(
    start, objective,
    function(par) -gradient(par),
    function(par) -hessian(par),
    control = list(trace = 1 * messages)
  )
  score <- gradient(estimate$par)
  information <- -hessian(estimate$par)
  information_root <- factor_covariance(
    information, "observed information matrix"
  )

  out <- list(
    estimate = structure_estimate(
      estimate$par, colnames(D), fix_tau2,
      sigma2_me = is.null(fix_var_me),
      re_names = if (n_re > 0L) colnames(ID_re) else NULL
    ),
    grad_MLE = score,
    covariance = chol2inv(information_root),
    log_lik = -estimate$objective,
    link_function = NULL,
    units_m = NULL,
    S_samples = NULL
  )
  parameter_names <- names(unlist(out$estimate))
  dimnames(out$covariance) <- list(parameter_names, parameter_names)
  diagnostics <- collect_optimizer_diagnostics(
    estimate, objective, score, information, information_root
  )
  warn_unconverged_optimizer(diagnostics)
  attr(out, "optimizer") <- diagnostics
  class(out) <- "RiskMap"
  out
}
