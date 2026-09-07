##' Set Control Parameters for Simulation
##'
##' This function sets control parameters for running simulations, supporting MCMC methods
##' for glgpm models.
##'
##' @param n_sim Integer. The total number of simulations to run. Default is 12000.
##' @param burnin Integer. The number of initial simulations to discard (burn-in/warmup period). Default is 2000.
##' @param thin Integer. The interval at which simulations are recorded (thinning interval, MCMC only). Default is 10.
##' @param h Numeric. An optional parameter for Langevin MCMC. Must be non-negative if specified.
##' @param c1.h Numeric. A control parameter for Langevin MCMC. Must be positive. Default is 0.01.
##' @param c2.h Numeric. Another control parameter for Langevin MCMC. Must be between 0 and 1. Default is 1e-04.
##' @param seed Integer. Optional value passed to `set.seed` to control random number generation for
##' generating chains and make results reproducible. Defaults to `NULL`.
##' @param linear_model Logical. If TRUE, sets up parameters for a linear model. Default is FALSE.
##'
##' @details
##' If \code{linear_model = TRUE}, only \code{n_sim} is required
##'
##' @return A list of control parameters with class "RiskMap_control_mcmc". Contents depend on \code{sampler}:
##' \itemize{
##'   \item For "mcmc": n_sim, burnin, thin, h, c1.h, c2.h, linear_model
##' }
##'
##' @examples
##' # Default parameters (MCMC)
##' control_mcmc <- set_control_mcmc()
##'
##' # Custom MCMC parameters
##' control_mcmc <- set_control_mcmc(n_sim = 15000, burnin = 3000, thin = 20)
##'
##' @seealso \code{\link{glgpm}}
##' @importFrom Matrix Matrix forceSymmetric
##' @export
set_control_mcmc <- function(n_sim = 12000,
                             burnin = 2000,
                             thin = 10,
                             h = NULL,
                             c1.h = 0.01,
                             c2.h = 1e-04,
                             seed = NULL,
                             linear_model = FALSE){

  if (!is.null(seed))
    check_positive_integer(seed, "seed")

  # =============================================================================
  # LINEAR MODEL (simple case for both samplers)
  # =============================================================================

  if (linear_model) {
    res <- list(
      n_sim = n_sim,
      linear_model = linear_model
    )
    class(res) <- "RiskMap_control_mcmc"
    return(res)
  }

  # =============================================================================
  # MCMC SAMPLER (Langevin for glgpm)
  # =============================================================================

  # Validate MCMC parameters
  if (n_sim < burnin) {
    stop("n_sim cannot be smaller than burnin.")
  }

  if (thin <= 0) {
    stop("thin must be positive")
  }

  if ((n_sim - burnin) %% thin != 0) {
    stop("thin must be a divisor of (n_sim - burnin)")
  }

  if (!is.null(h) && h < 0) {
    stop("h must be non-negative.")
  }

  if (c1.h <= 0) {
    stop("c1.h must be positive.")
  }

  if (c2.h < 0 | c2.h > 1) {
    stop("c2.h must be between 0 and 1.")
  }

  res <- list(
    n_sim = n_sim,
    burnin = burnin,
    thin = thin,
    h = h,
    c1.h = c1.h,
    c2.h = c2.h,
    seed = seed,
    linear_model = FALSE
  )

  class(res) <- "RiskMap_control_mcmc"
  return(res)
}

##' @title Laplace-sampling MCMC for Generalized Linear Gaussian Process Models
##'
##' @description
##' Runs Markov chain Monte Carlo (MCMC) sampling using a Laplace
##' approximation for Generalized Linear Gaussian Process Models (GLGPMs).
##' The latent Gaussian field is integrated via a second-order Taylor
##' expansion around the mode, and a Gaussian proposal is used for
##' Metropolis–Hastings updates with adaptive step-size tuning.
##'
##' @param y Numeric vector of responses of length \eqn{n}.
##'   For \code{family = "binomial"} this is the number of successes,
##'   for \code{family = "poisson"} counts, and for \code{family = "gaussian"} real values.
##' @param units_m Numeric vector giving the binomial totals (number of trials)
##'   when \code{family = "binomial"}; ignored for other families (can be \code{NULL}).
##' @param mu Numeric vector of length equal to the number of unique locations
##'   providing the mean of the latent spatial process on the link scale.
##' @param Sigma Numeric positive-definite covariance matrix for the latent spatial
##'   process \eqn{S} at the unique locations referenced by \code{ID_coords}.
##' @param ID_coords Integer vector of length \eqn{n} mapping each response in
##'   \code{y} to a row/column of \code{Sigma} (i.e., the index of the corresponding location).
##' @param ID_re Optional matrix or data.frame with one column per unstructured
##'   random effect (RE). Each column is an integer vector of length \eqn{n}
##'   mapping observations in \code{y} to RE levels (e.g., cluster, survey, etc.).
##'   Use \code{NULL} to exclude REs.
##' @param sigma2_re Optional named numeric vector of RE variances. Names must
##'   match the column names of \code{ID_re}. Ignored if \code{ID_re = NULL}.
##' @param family Character string: one of \code{"gaussian"}, \code{"binomial"},
##'   or \code{"poisson"}.
##' @param control_mcmc List of control parameters:
##'   \describe{
##'     \item{n_sim}{Total number of MCMC iterations (including burn-in).}
##'     \item{burnin}{Number of initial iterations to discard.}
##'     \item{thin}{Thinning interval for saving samples.}
##'     \item{h}{Initial step size for the Gaussian proposal. Defaults to \eqn{1.65 / n_\mathrm{tot}^{1/6}} if not supplied.}
##'     \item{c1.h, c2.h}{Positive tuning constants for adaptive step-size updates.}
##'   }
##' @param invlink Optional inverse-link function. If \code{NULL}, defaults are used:
##'   \code{identity} (gaussian), \code{plogis} (binomial), and \code{exp} (poisson).
##' @param Sigma_pd Optional precision matrix used in the Laplace approximation.
##'   If \code{NULL}, it is obtained internally at the current mode.
##' @param mean_pd Optional mean vector used in the Laplace approximation.
##'   If \code{NULL}, it is obtained internally as the mode of the integrand.
##' @param messages Logical; if \code{TRUE}, prints progress and acceptance diagnostics.
##'
##' @details
##' The algorithm alternates between:
##' \enumerate{
##'   \item Locating the mode of the joint integrand for the latent variables
##'         (via \code{maxim_integrand}) when \code{Sigma_pd} and \code{mean_pd}
##'         are not provided, yielding a Gaussian approximation.
##'   \item Metropolis–Hastings updates using a Gaussian proposal centered at
##'         the current approximate mean with proposal variance governed by \code{h}.
##'         The step size is adapted based on empirical acceptance probability.
##' }
##'
##' Dimensions must be consistent:
##' \code{length(y) = n}, \code{nrow(Sigma) = ncol(Sigma) = n_loc},
##' and \code{length(ID_coords) = n} with entries in \eqn{1,\dots,n_\mathrm{loc}}.
##' If \code{ID_re} is provided, each column must have length \eqn{n}; when
##' \code{sigma2_re} is supplied, it must be named and match \code{colnames(ID_re)}.
##'
##' @return An object of class \code{"RiskMap_mcmc"} with components:
##' \describe{
##'   \item{samples}{A list containing posterior draws. Always includes
##'                 \code{$S} (latent spatial field). If \code{ID_re} is supplied,
##'                 each unstructured RE is returned under \code{$<re_name>}.}
##'   \item{tuning_par}{Numeric vector of step sizes (\code{h}) used over iterations.}
##'   \item{acceptance_prob}{Numeric vector of Metropolis–Hastings acceptance probabilities.}
##' }
##'
##' @section Default links:
##' The default inverse links are: identity (gaussian), logistic (binomial),
##' and exponential (poisson). Supply \code{invlink} to override.
##'
##' @seealso \code{\link{maxim_integrand}}
##'
##' @export
laplace_sampling_mcmc <- function(y,
                                  units_m,
                                  mu,
                                  Sigma,
                                  ID_coords,
                                  ID_re = NULL,
                                  sigma2_re = NULL,
                                  family,
                                  control_mcmc,
                                  invlink = NULL,
                                  Sigma_pd = NULL,
                                  mean_pd = NULL,
                                  messages = TRUE
){

  stopifnot(family %in% c("poisson", "binomial"))

  # set seed if it exists and reset on exit
  if (!is.null(control_mcmc$seed)){
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(control_mcmc$seed)
  }

  # ---------- utilities ----------
  check_vec_fun <- function(f, n, name) {
    if (!is.function(f)) stop(sprintf("`%s` must be a function.", name))
    x <- rep(0, n)
    out <- tryCatch(f(x), error = function(e) e)
    if (inherits(out, "error")) stop(sprintf("`%s` failed on a numeric vector: %s", name, out$message))
    if (!is.numeric(out)) stop(sprintf("`%s` must return a numeric vector.", name))
    if (length(out) != length(x)) stop(sprintf("`%s` must return a vector of the same length as its input.", name))
    if (!all(is.finite(out))) stop(sprintf("`%s` returns non-finite values.", name))
    invisible(TRUE)
  }

  sum_by_group <- function(v, grp, nlev) {
    f <- factor(grp, levels = seq_len(nlev))
    s <- tapply(v, f, sum)
    ans <- rep(0, nlev)
    if (!is.null(s)) ans[seq_along(s)] <- replace(s, is.na(s), 0)
    as.numeric(ans)
  }

  # ---------- dimensions ----------
  Sigma.inv <- solve(Sigma)
  n_loc <- nrow(Sigma)
  n <- length(y)

  if (( !is.null(ID_re) && is.null(sigma2_re)) ||
      (  is.null(ID_re) && !is.null(sigma2_re))) {
    stop("To introduce unstructured random effects both `ID_re` and `sigma2_re` must be provided.")
  }

  if (is.null(ID_re)) {
    n_re <- 0L
    n_dim_re <- integer(0)
    ind_re <- list()
  } else {
    n_re <- length(sigma2_re)
    n_dim_re <- vapply(seq_len(n_re), function(i) length(unique(ID_re[, i])), integer(1))
    ind_re <- vector("list", n_re)
    add_i <- 0L
    for (i in seq_len(n_re)) {
      ind_re[[i]] <- (add_i + n_loc + 1):(add_i + n_loc + n_dim_re[i])
      if (i < n_re) add_i <- sum(n_dim_re[1:i])
    }
  }
  n_tot <- n_loc + if (n_re > 0) sum(n_dim_re) else 0L

  # ---------- inverse link handling (inv, d1) ----------
  make_invlink_funs <- function(family, invlink, ncheck) {
    have_Deriv <- requireNamespace("Deriv", quietly = TRUE)
    have_numDeriv <- requireNamespace("numDeriv", quietly = TRUE)

    if (is.null(invlink)) {
      if (family == "poisson") {
        inv <- function(x) exp(x)
        d1  <- function(x) exp(x)
      } else {
        inv <- function(x) plogis(x)
        d1  <- function(x) { p <- inv(x); p * (1 - p) }
      }
      check_vec_fun(inv, ncheck, "canonical invlink")
      check_vec_fun(d1,  ncheck, "canonical invlink_prime")
      return(list(inv = inv, d1 = d1, name = "canonical"))
    }

    if (is.function(invlink)) {
      inv_user <- invlink
      d1_user  <- NULL
    } else if (is.list(invlink)) {
      inv_user <- invlink$inv %||% invlink$inv_link %||% invlink$invlink
      d1_user  <- invlink$d1  %||% invlink$inv_link_prime %||% invlink$mu_eta
    } else {
      stop("'invlink' must be NULL, a function, or a list with components inv and d1.")
    }

    check_vec_fun(inv_user, ncheck, "invlink")

    if (is.null(d1_user)) {
      if (have_Deriv) {
        inv_wrapped <- function(eta) inv_user(eta)
        d1_user <- Deriv::Deriv(inv_wrapped, "eta")
      } else if (have_numDeriv) {
        d1_user <- function(eta) vapply(eta, function(z)
          numDeriv::grad(function(x) inv_user(x), z), numeric(1))
      } else {
        stop("Cannot auto-derive first derivative. Install `Deriv` or `numDeriv`, or provide `d1`.")
      }
    }
    check_vec_fun(d1_user, ncheck, "invlink_prime")

    list(inv = inv_user, d1 = d1_user, name = "custom")
  }

  linkf <- make_invlink_funs(family, invlink, n)
  inv_link <- linkf$inv
  inv1     <- linkf$d1

  # ---------- default Laplace proposal if missing ----------
  if (is.null(Sigma_pd) || is.null(mean_pd)) {
    out_maxim <- maxim_integrand(y = y, units_m = units_m, Sigma = Sigma, mu = mu,
                                 ID_coords = ID_coords, ID_re = ID_re,
                                 sigma2_re = sigma2_re,
                                 family = family, invlink = invlink,
                                 hessian = FALSE, gradient = TRUE)
    if (is.null(Sigma_pd)) Sigma_pd <- out_maxim$Sigma.tilde
    if (is.null(mean_pd))  mean_pd  <- out_maxim$mode
  }

  # ---------- affine reparameterisation ----------
  n_sim   <- control_mcmc$n_sim
  Sigma_pd_sroot <- t(chol(Sigma_pd))
  A <- solve(Sigma_pd_sroot)

  if (n_re == 0) {
    Sigma_tot <- Sigma
  } else {
    Sigma_tot <- matrix(0, n_tot, n_tot)
    Sigma_tot[1:n_loc, 1:n_loc] <- Sigma
    for (i in seq_len(n_re)) diag(Sigma_tot)[ind_re[[i]]] <- sigma2_re[i]
  }
  Sigma_w_inv <- solve(A %*% Sigma_tot %*% t(A))
  mu_w <- -as.numeric(A %*% mean_pd)

  cond.dens.W <- function(W, S_tot) {
    S <- S_tot[1:n_loc]
    S_re_list <- if (n_re > 0) lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]]) else NULL

    eta <- mu + S[ID_coords]
    if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]

    if (family == "poisson") {
      mu_vec <- inv_link(eta)
      if (any(!is.finite(mu_vec)) || any(mu_vec <= 0)) stop("invlink must return positive means for Poisson.")
      llik <- sum(y * log(pmax(mu_vec, .Machine$double.eps)) - units_m * mu_vec)
    } else {
      p <- inv_link(eta)
      if (any(!is.finite(p)) || any(p <= 0 | p >= 1)) stop("invlink must return values in (0,1) for Binomial.")
      llik <- sum(y * log(pmax(p, .Machine$double.eps)) +
                    (units_m - y) * log(pmax(1 - p, .Machine$double.eps)))
    }
    diff_w <- W - mu_w
    as.numeric(-0.5 * crossprod(diff_w, Sigma_w_inv %*% diff_w) + llik)
  }

  lang.grad <- function(W, S_tot) {
    S <- S_tot[1:n_loc]
    S_re_list <- if (n_re > 0) lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]]) else NULL

    eta <- mu + S[ID_coords]
    if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]

    if (family == "poisson") {
      mu_vec <- inv_link(eta)
      if (any(!is.finite(mu_vec)) || any(mu_vec <= 0)) stop("invlink must return positive means for Poisson.")
      mu1 <- inv1(eta)
      g_eta <- (y - units_m * mu_vec) * (mu1 / mu_vec)
    } else {
      p <- inv_link(eta)
      if (any(!is.finite(p)) || any(p <= 0 | p >= 1)) stop("invlink must return values in (0,1) for Binomial.")
      p1 <- inv1(eta)
      den <- p * (1 - p)
      g_eta <- (y - units_m * p) * (p1 / den)
    }

    grad_S_tot <- numeric(n_tot)
    grad_S_tot[1:n_loc] <- sum_by_group(g_eta, ID_coords, n_loc)
    if (n_re > 0) {
      for (j in seq_len(n_re)) {
        grad_S_tot[ind_re[[j]]] <- sum_by_group(g_eta, ID_re[, j], n_dim_re[j])
      }
    }

    as.numeric(-Sigma_w_inv %*% (W - mu_w) + t(Sigma_pd_sroot) %*% grad_S_tot)
  }

  # ---------- MALA tuning ----------
  h      <- control_mcmc$h; if (is.null(h)) h <- 1.65/(n_tot^(1/6))
  burnin <- control_mcmc$burnin
  thin   <- control_mcmc$thin
  c1.h   <- control_mcmc$c1.h
  c2.h   <- control_mcmc$c2.h

  W_curr <- rep(0, n_tot)
  S_tot_curr <- as.numeric(Sigma_pd_sroot %*% W_curr + mean_pd)
  mean_curr <- as.numeric(W_curr + (h^2/2) * lang.grad(W_curr, S_tot_curr))
  lp_curr <- cond.dens.W(W_curr, S_tot_curr)
  acc <- 0L
  n_samples <- floor((n_sim - burnin) / thin)   # was: (n_sim - burnin) / thin
  sim <- matrix(NA_real_, nrow = n_samples, ncol = n_tot)

  # ---- safe progress output (no stderr, no flush.console) ----
  if (messages) message("\n - Conditional simulation (burnin=", burnin, ", thin=", thin, "):")
  pb <- NULL
  if (messages && interactive()) {
    pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  h.vec <- rep(NA_real_, n_sim)
  acc_prob <- rep(NA_real_, n_sim)

  for (i in seq_len(n_sim)) {
    W_prop <- mean_curr + h * rnorm(n_tot)
    S_tot_prop <- as.numeric(Sigma_pd_sroot %*% W_prop + mean_pd)
    mean_prop <- as.numeric(W_prop + (h^2/2) * lang.grad(W_prop, S_tot_prop))
    lp_prop <- cond.dens.W(W_prop, S_tot_prop)

    dprop_curr <- -sum((W_prop - mean_curr)^2) / (2 * h^2)
    dprop_prop <- -sum((W_curr - mean_prop)^2) / (2 * h^2)

    log_prob <- lp_prop + dprop_prop - lp_curr - dprop_curr

    if (log(runif(1)) < log_prob) {
      acc <- acc + 1L
      W_curr <- W_prop
      S_tot_curr <- S_tot_prop
      lp_curr <- lp_prop
      mean_curr <- mean_prop
    }

    if (i > burnin && (i - burnin) %% thin == 0) {
      cnt <- (i - burnin) %/% thin
      sim[cnt, ] <- S_tot_curr
    }

    acc_prob[i] <- acc / i
    h <- max(1e-19, h + c1.h * i^(-c2.h) * (acc / i - 0.57))
    h.vec[i] <- h

    if (messages) {
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, i)
      } else if (i %% max(1, floor(n_sim/20)) == 0) {
        # non-interactive fallback: update every ~5%
        message(sprintf("   %3d%%", round(100 * i / n_sim)))
      }
    }
  }

  if (!is.null(pb)) close(pb)
  if (messages) message(" done.\n")

  out_sim <- list()
  out_sim$samples <- list()
  out_sim$samples$S <- sim[, 1:n_loc, drop = FALSE]
  if (n_re > 0) {
    re_names <- if (!is.null(colnames(ID_re))) colnames(ID_re) else paste0("re", seq_len(n_re))
    for (i in seq_len(n_re)) {
      out_sim$samples[[re_names[i]]] <- sim[, ind_re[[i]], drop = FALSE]
    }
  }
  out_sim$tuning_par <- h.vec
  out_sim$acceptance_prob <- acc_prob
  out_sim$invlink_used <- linkf$name
  class(out_sim) <- "RiskMap_mcmc"
  out_sim
}




##' Check MCMC Convergence for Spatial Random Effects
##'
##' This function checks the Markov Chain Monte Carlo (MCMC) convergence of spatial random effects
##' for either a \code{RiskMap} or \code{RiskMap_pred} object.
##' It plots the trace plot and autocorrelation function (ACF) for the MCMC chain
##' and calculates the effective sample size (ESS).
##'
##' @param object An object of class \code{RiskMap} or \code{RiskMap_pred}.
##'  \code{RiskMap} is the output from \code{\link{glgpm}} function, and
##'  \code{RiskMap_pred} is obtained from the \code{\link{setup_prediction}} function.
##' @param check_mean Logical. If \code{TRUE}, checks the MCMC chain for the mean of the spatial random effects.
##'  If \code{FALSE}, checks the chain for a specific component of the random effects vector.
##' @param component Integer. The index of the spatial random effects component to check when \code{check_mean = FALSE}.
##'  Must be a positive integer corresponding to a location in the data. Ignored if \code{check_mean = TRUE}.
##' @param ... Additional arguments passed to the \code{\link[stats]{acf}} function for customizing the ACF plot.
##'
##' @details
##' The function first checks that the input object is either of class \code{RiskMap} or \code{RiskMap_pred}.
##' Depending on the value of \code{check_mean}, it either calculates the mean of the spatial random effects
##' across all locations for each iteration or uses the specified component.
##' It then generates two plots:
##' - A trace plot of the selected spatial random effect over iterations.
##' - An autocorrelation plot (ACF) with the effective sample size (ESS) displayed in the title.
##'
##' The ESS is computed using the \code{\link[sns]{ess}} function, which provides a measure of the effective number
##' of independent samples in the MCMC chain.
##'
##' If \code{check_mean = TRUE}, the \code{component} argument is ignored, and a warning is issued.
##' To specify a particular component of the random effects vector, set \code{check_mean = FALSE} and provide
##' a valid \code{component} value.
##'
##' @return
##' No return value, called for side effects (plots and warnings).
##' @importFrom sns ess
##' @importFrom graphics par
##' @export
plot_mcmc <- function(object, check_mean = TRUE,
                      component = NULL, ...) {
  if(!inherits(object, "RiskMap") &
     !inherits(object, "RiskMap_pred")) {
    stop("'object' must be either:
           a 'RiskMap' object obtained as an output from 'glgpm';
           a 'RiskMap_pred' object obtained as an output from 'setup_prediction'")
  }

  if (object$family == "gaussian")
    stop("'object' is a gaussian model which cannot contain MCMC chains")

  if (is.null(object$S_samples))
    stop("'object' does not contain any MCMC chains - rerun 'glgpm' with 'return_samples' = TRUE")

  if(inherits(object, "RiskMap")) {
    S_samples <- object$S_samples
  } else if (inherits(object, "RiskMap_pred")) {
    S_samples <- t(object$S_samples)
  }

  if(check_mean & !is.null(component)) {
    warning("if check_mean = TRUE, the value passed to 'component' is ignored;
            set check_mean = FALSE when specifying a value for 'component'")
  }
  n_samples <- nrow(S_samples)
  n_loc <- ncol(S_samples)
  if(check_mean) {
    S_chain <- apply(S_samples, 1, mean)
  } else {
    if(is.null(component)) stop("When 'check_mean' = FALSE a component of the
                                random effects vector must be specified through 'component'
                                by providing a positive integer")
    check_positive_integer(component, "component")
    if(component > n_loc) stop("'component' must be a single positive integer
                               between 1 and the number of locations in the data")
    S_chain <- S_samples[,component]
  }

  par(mfrow = c(1,2))
  plot(S_chain, type = "l",
       ylab = "", xlab = "Iteration",
       main = "Spatial random effect")

  S_chain_ess <- round(ess(S_chain),3)
  acf(S_chain, main = paste("Effective sample size:",
                            S_chain_ess), ...)
  par(mfrow = c(1,1))
}
