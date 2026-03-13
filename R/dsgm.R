##' @title Detect available Stan backend
##' @description Prefers cmdstanr because its compiled executable persists
##'   across R sessions without DSO issues. Falls back to rstan only when
##'   cmdstanr is unavailable. With rstan the model is recompiled every new
##'   R session (the DSO is process-local and cannot be cached on disk).
##' @keywords internal
detect_stan_backend <- function(messages = TRUE) {
  # Try cmdstanr first — executable survives across R sessions
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    ok <- tryCatch({ cmdstanr::cmdstan_version(); TRUE }, error = function(e) FALSE)
    if (ok) {
      if (messages) message("Using cmdstanr backend")
      return("cmdstanr")
    }
  }
  # Fall back to rstan
  if (requireNamespace("rstan", quietly = TRUE)) {
    if (messages) message("Using rstan backend (cmdstanr unavailable)")
    return("rstan")
  }
  stop("Neither rstan nor cmdstanr is available. Install cmdstanr for persistent caching across sessions.")
}


##' @title Compile and cache a Stan model
##' @description Returns a compiled Stan model for either the STH or LF model.
##'   Uses separate cache slots so both models can coexist in a session.
##'
##'   For \code{cmdstanr} (recommended), the compiled executable is written to
##'   a persistent user cache directory (\code{tools::R_user_dir("RiskMap", "cache")})
##'   and reused across R sessions without recompilation. For \code{rstan},
##'   the model DSO is process-local and must be relinked every new R session;
##'   use cmdstanr to avoid this.
##'
##' @param model \code{"sth"} uses \code{inst/stan/dsgm_spatial.stan};
##'   \code{"lf"} uses \code{inst/stan/dsgm_mdiag.stan}.
##' @param backend \code{"rstan"}, \code{"cmdstanr"}, or \code{NULL}
##'   (auto-detect).
##' @keywords internal
get_stan_model <- function(model = c("sth", "lf"), backend = NULL,
                           messages = TRUE) {
  model <- match.arg(model)

  cache_model   <- paste0("stan_model_",   model)
  cache_backend <- paste0("stan_backend_", model)

  # ------------------------------------------------------------------
  # 1. In-session cache: resolve backend FIRST so the key always matches
  #    what was stored. Passing backend=NULL from pred_over_grid would
  #    otherwise always miss the cache (NULL != "rstan").
  # ------------------------------------------------------------------
  if (is.null(backend)) {
    # If the model is already cached, re-use the backend it was compiled with
    # (avoids running detect_stan_backend() unnecessarily)
    if (!is.null(.dsgm_cache[[cache_backend]])) {
      backend <- .dsgm_cache[[cache_backend]]
    } else {
      backend <- detect_stan_backend(messages = messages)
    }
  }

  if (!is.null(.dsgm_cache[[cache_model]]) &&
      identical(.dsgm_cache[[cache_backend]], backend)) {
    if (messages) message("Using cached Stan model (", model, ")")
    return(list(model = .dsgm_cache[[cache_model]], backend = backend))
  }

  stan_file <- switch(model,
                      sth = system.file("stan/dsgm_spatial.stan", package = "RiskMap"),
                      lf  = system.file("stan/dsgm_mdiag.stan",   package = "RiskMap")
  )
  if (!file.exists(stan_file))
    stop("Stan model file not found: ", basename(stan_file))

  # ------------------------------------------------------------------
  # 2. Persistent on-disk cache directory (user-writable, survives sessions)
  # ------------------------------------------------------------------
  cache_dir <- tools::R_user_dir("RiskMap", "cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  if (backend == "rstan") {
    # rstan stan_model objects contain a live reference to a compiled DSO
    # (shared library). saveRDS/readRDS does NOT work across sessions because
    # the DSO is no longer loaded when the RDS is read back in a new session,
    # causing "NULL value passed for DllInfo".
    #
    # The correct approach is to let rstan manage its own caching via
    # auto_write = TRUE, which saves a .rds alongside the .stan file and
    # reloads the DSO correctly. The package library is read-only, so we copy
    # the .stan file into our writable cache dir once and always work from
    # there. rstan then reads/writes its .rds cache next to that copy.
    cached_stan <- file.path(cache_dir, basename(stan_file))

    # Copy .stan to cache dir if missing or source is newer
    if (!file.exists(cached_stan) ||
        file.mtime(stan_file) > file.mtime(cached_stan)) {
      file.copy(stan_file, cached_stan, overwrite = TRUE)
    }

    # stan_model with auto_write = TRUE: compiles once, reloads DSO on
    # subsequent calls without recompiling the C++ code
    rds_path <- sub("\\.stan$", ".rds", cached_stan)

    if (file.exists(rds_path)) {
      if (messages) message("Loading pre-compiled rstan model (", model, ")")
      # readRDS gives us the R-side object; we then need to ensure the DSO is
      # loaded. rstan::stan_model() with auto_write=TRUE does this correctly
      # (it reloads the .so without recompiling the C++). Calling it with the
      # cached .stan copy is safe and much faster than a full recompile.
      compiled <- rstan::stan_model(
        file       = cached_stan,
        model_name = paste0("dsgm_", model),
        verbose    = FALSE,
        auto_write = TRUE
      )
    } else {
      if (messages) message("Compiling Stan model '", basename(stan_file), "'...")
      compiled <- rstan::stan_model(
        file       = cached_stan,
        model_name = paste0("dsgm_", model),
        verbose    = FALSE,
        auto_write = TRUE
      )
      if (messages) message("  Saved to: ", rds_path)
    }

  } else {
    # cmdstanr: we manage staleness ourselves rather than relying on
    # cmdstan_model(compile = TRUE), because installed package files often have
    # their mtimes reset (e.g. on HPC systems), making the .stan file appear
    # newer than the cached exe even when nothing has changed.
    exe_path <- file.path(cache_dir,
                          paste0("dsgm_", model,
                                 if (.Platform$OS.type == "windows") ".exe" else ""))

    exe_is_fresh <- file.exists(exe_path) &&
      file.mtime(exe_path) >= file.mtime(stan_file)

    if (exe_is_fresh) {
      if (messages) message("Loading pre-compiled cmdstanr model (", model, ")")
      # Load without recompiling: pass the existing exe directly
      compiled <- cmdstanr::cmdstan_model(
        stan_file = stan_file,
        exe_file  = exe_path,
        compile   = FALSE
      )
    } else {
      if (messages) message("Compiling Stan model '", basename(stan_file), "'...")
      compiled <- cmdstanr::cmdstan_model(
        stan_file = stan_file,
        exe_file  = exe_path,
        compile   = TRUE
      )
      if (messages) message("  Saved to: ", exe_path)
    }
  }

  # Store in in-session cache
  .dsgm_cache[[cache_model]]   <- compiled
  .dsgm_cache[[cache_backend]] <- backend

  if (messages) message("Stan model ready (", model, ")")
  return(list(model = compiled, backend = backend))
}


# =============================================================================
# Stan samplers
# =============================================================================

##' @title Sample spatial process for the STH model
##' @description Samples S(x) from its posterior at fixed theta_0 using
##'   \code{inst/stan/dsgm_spatial.stan}.
##' @param par Named list with fields: \code{beta}, \code{k}, \code{rho},
##'   \code{alpha_W}, \code{gamma_W}, \code{sigma2}, \code{phi}.
##' @return An object of class \code{dsgm_spatial_samples}.
##' @keywords internal
sample_spatial_process_stan <- function(y_prev,
                                        intensity_data,
                                        D,
                                        coords,
                                        ID_coords,
                                        int_mat,
                                        survey_times_data,
                                        mda_times,
                                        par,
                                        n_samples     = 1000,
                                        n_warmup      = 1000,
                                        n_chains      = 4,
                                        n_cores       = 4,
                                        adapt_delta   = 0.8,
                                        max_treedepth = 10,
                                        backend       = NULL,
                                        messages      = TRUE) {

  n       <- length(y_prev)
  n_loc   <- nrow(coords)
  p       <- ncol(D)
  pos_idx <- which(y_prev == 1)
  n_pos   <- length(pos_idx)

  mda_impact <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                                   par$alpha_W, par$gamma_W, kappa = 1)

  stan_data <- list(
    n         = n,
    n_loc     = n_loc,
    n_pos     = n_pos,
    p         = p,
    y         = y_prev,
    C_pos     = intensity_data,
    pos_idx   = pos_idx,
    ID_coords = ID_coords,
    D_mat     = as.matrix(dist(coords)),
    eta_fixed = as.numeric(D %*% par$beta),
    mda_impact = mda_impact,
    k          = par$k,
    rho        = par$rho,
    sigma2     = par$sigma2,
    phi        = par$phi
  )

  mod <- get_stan_model(model = "sth", backend = backend, messages = messages)
  backend <- mod$backend

  if (messages)
    message(sprintf("Sampling %d iter (%d warmup), %d chain(s) [sth]...",
                    n_samples + n_warmup, n_warmup, n_chains))

  fit <- .run_stan(mod$model, backend, stan_data,
                   n_samples, n_warmup, n_chains, n_cores,
                   adapt_delta, max_treedepth, messages)
  S_samples <- .extract_S(fit, backend, messages)

  result <- list(S_samples = S_samples, stan_fit = fit,
                 n_samples = nrow(S_samples), n_loc = n_loc,
                 coords = coords, par = par, backend = backend)
  class(result) <- "dsgm_spatial_samples"
  result
}


##' @title Sample spatial process for the LF multi-diagnostic model
##' @description Samples S(x) from its posterior at fixed theta_0 using
##'   \code{inst/stan/dsgm_mdiag.stan}. All model parameters are fixed and
##'   passed as Stan data; only S_raw is sampled.
##'
##'   \strong{Parameter name translation:} internally the package uses
##'   \code{k}/\code{rho} for the NB aggregation and detection-rate parameters.
##'   The Stan file (\code{dsgm_mdiag.stan}) and the TMB template
##'   (\code{dsgm_mdiag.cpp}) use \code{omega}/\code{alpha} for the same
##'   quantities.  The translation is done here at the R/Stan boundary.
##'
##' @param par Named list: \code{beta}, \code{k}, \code{rho},
##'   \code{gamma_sens}, \code{sigma2}, \code{phi}, \code{tau2}.
##' @param mda_impact Pre-computed MDA decay vector (length n).
##'   Pass \code{rep(1, n)} when there is no MDA.
##' @return An object of class \code{dsgm_spatial_samples}.
##' @keywords internal
sample_spatial_process_stan_lf <- function(y_counts,
                                           units_m,
                                           which_diag,
                                           D,
                                           coords,
                                           ID_coords,
                                           par,
                                           mda_impact    = NULL,
                                           n_samples     = 1000,
                                           n_warmup      = 1000,
                                           n_chains      = 4,
                                           n_cores       = 4,
                                           adapt_delta   = 0.8,
                                           max_treedepth = 10,
                                           backend       = NULL,
                                           messages      = TRUE) {

  n     <- length(y_counts)
  n_loc <- nrow(coords)
  p     <- ncol(D)

  if (is.null(mda_impact)) mda_impact <- rep(1.0, n)
  use_mda_flag <- as.integer(!all(mda_impact == 1.0))

  # Translate R naming (k/rho) -> Stan naming (omega/alpha)
  stan_data <- list(
    n          = n,
    n_loc      = n_loc,
    p          = p,
    y          = as.integer(y_counts),
    units_m    = as.integer(units_m),
    is_mf      = as.integer(which_diag),
    ID_coords  = as.integer(ID_coords),
    D_mat      = as.matrix(dist(coords)),
    eta_fixed  = as.numeric(D %*% par$beta),
    mda_impact = as.numeric(mda_impact),
    use_mda    = use_mda_flag,
    # Fixed parameters (R k/rho -> Stan omega/alpha)
    omega      = par$k,
    alpha      = par$rho,
    gamma_sens = par$gamma_sens,
    sigma2     = par$sigma2,
    phi        = par$phi
    # NOTE: tau2 (nugget) is NOT passed to Stan; it is a TMB-only parameter
  )

  mod <- get_stan_model(model = "lf", backend = backend, messages = messages)
  backend <- mod$backend

  if (messages)
    message(sprintf("Sampling %d iter (%d warmup), %d chain(s) [lf_mdiag]...",
                    n_samples + n_warmup, n_warmup, n_chains))

  fit <- .run_stan(mod$model, backend, stan_data,
                   n_samples, n_warmup, n_chains, n_cores,
                   adapt_delta, max_treedepth, messages)
  S_samples <- .extract_S(fit, backend, messages)

  result <- list(S_samples = S_samples, stan_fit = fit,
                 n_samples = nrow(S_samples), n_loc = n_loc,
                 coords = coords, par = par, backend = backend)
  class(result) <- "dsgm_spatial_samples"
  result
}


# Internal helpers to avoid duplicating backend dispatch logic
.run_stan <- function(stan_model, backend, stan_data,
                      n_samples, n_warmup, n_chains, n_cores,
                      adapt_delta, max_treedepth, messages) {
  if (backend == "rstan") {
    rstan::sampling(
      stan_model,
      data    = stan_data,
      iter    = n_samples + n_warmup,
      warmup  = n_warmup,
      chains  = n_chains,
      cores   = n_cores,
      control = list(adapt_delta = adapt_delta,
                     max_treedepth = max_treedepth),
      refresh = ifelse(messages, max(1, (n_samples + n_warmup) %/% 10), 0),
      show_messages = messages,
      verbose       = messages
    )
  } else {
    stan_model$sample(
      data            = stan_data,
      iter_warmup     = n_warmup,
      iter_sampling   = n_samples,
      chains          = n_chains,
      parallel_chains = n_cores,
      adapt_delta     = adapt_delta,
      max_treedepth   = max_treedepth,
      refresh = ifelse(messages, max(1, (n_samples + n_warmup) %/% 10), 0),
      show_messages   = messages,
      show_exceptions = messages
    )
  }
}

.extract_S <- function(fit, backend, messages) {
  if (backend == "rstan") {
    S <- rstan::extract(fit, pars = "S")$S
    if (messages) {
      n_div <- sum(rstan::get_num_divergent(fit))
      n_mtr <- sum(rstan::get_num_max_treedepth(fit))
      message(sprintf("Extracted %d samples (%d locations) | divergent: %d | max-treedepth: %d",
                      nrow(S), ncol(S), n_div, n_mtr))
      if (n_div > 0) warning("Divergent transitions detected. Consider increasing adapt_delta.")
    }
  } else {
    S <- fit$draws("S", format = "matrix")
    if (messages) {
      diag  <- fit$diagnostic_summary()
      n_div <- sum(diag$num_divergent)
      n_mtr <- sum(diag$num_max_treedepth)
      message(sprintf("Extracted %d samples (%d locations) | divergent: %d | max-treedepth: %d",
                      nrow(S), ncol(S), n_div, n_mtr))
      if (n_div > 0) warning("Divergent transitions detected. Consider increasing adapt_delta.")
    }
  }
  S
}


# =============================================================================
# Utility functions for spatial samples
# =============================================================================

##' @title Thin spatial process samples
##' @description Reduces autocorrelation by keeping every \code{thin}-th row.
##' @param spatial_samples Output from \code{sample_spatial_process_stan} or
##'   \code{sample_spatial_process_stan_lf}.
##' @param thin Thinning interval (default 10).
##' @return A thinned \code{dsgm_spatial_samples} object.
##' @keywords internal
thin_spatial_samples <- function(spatial_samples, thin = 10) {
  keep <- seq(1, nrow(spatial_samples$S_samples), by = thin)
  spatial_samples$S_samples <- spatial_samples$S_samples[keep, , drop = FALSE]
  spatial_samples$n_samples  <- nrow(spatial_samples$S_samples)
  spatial_samples
}


##' @title Effective sample size for spatial process
##' @description Computes ESS for each spatial location using the \pkg{coda}
##'   package.
##' @param spatial_samples Output from a Stan sampler.
##' @return Numeric vector of ESS values, one per location.
##' @keywords internal
compute_spatial_ess <- function(spatial_samples) {
  if (!requireNamespace("coda", quietly = TRUE))
    stop("Package 'coda' is required for ESS calculation")
  sapply(seq_len(spatial_samples$n_loc), function(i)
    coda::effectiveSize(coda::as.mcmc(spatial_samples$S_samples[, i])))
}


##' @title Print method for dsgm_spatial_samples
##' @export
print.dsgm_spatial_samples <- function(x, ...) {
  cat("DSGM Spatial Process Samples\n")
  cat(sprintf("  Samples   : %d\n", x$n_samples))
  cat(sprintf("  Locations : %d\n", x$n_loc))
  cat(sprintf("  Matrix    : %d x %d\n", nrow(x$S_samples), ncol(x$S_samples)))
  if (requireNamespace("coda", quietly = TRUE)) {
    ess <- compute_spatial_ess(x)
    cat(sprintf("  ESS       : %.0f - %.0f  (median %.0f)\n",
                min(ess), max(ess), median(ess)))
  }
  invisible(x)
}


# =============================================================================
# Initial value functions
# =============================================================================

##' @title Compute initial parameter values for the STH model
##' @description Fits a simplified non-spatial model (no S(x)) to obtain
##'   starting values for \code{beta}, \code{k}, \code{rho}, \code{alpha_W},
##'   \code{gamma_W}, \code{sigma2}, and \code{phi}.
##' @keywords internal
dsgm_initial_value <- function(y_prev, intensity_data, D, coords, ID_coords,
                               int_mat, survey_times_data, mda_times,
                               penalty, fix_alpha_W, fix_gamma_W,
                               start_pars, messages) {

  n <- length(y_prev)
  p <- ncol(D)

  # Penalty helper
  .pen <- function(alpha_W, gamma_W) {
    pen <- 0
    if (!is.null(penalty)) {
      if (is.null(fix_alpha_W)) {
        if (is.function(penalty$alpha))
          pen <- pen + penalty$alpha(alpha_W)
        else if (!is.null(penalty$alpha_param1))
          pen <- pen - (penalty$alpha_param1-1)*log(alpha_W) -
            (penalty$alpha_param2-1)*log(1-alpha_W)
      }
      if (is.null(fix_gamma_W)) {
        if (is.function(penalty$gamma))
          pen <- pen + penalty$gamma(gamma_W)
        else if (!is.null(penalty$gamma_type)) {
          if (penalty$gamma_type == "gamma")
            pen <- pen - (penalty$gamma_shape-1)*log(gamma_W) + penalty$gamma_rate*gamma_W
          else {
            d <- gamma_W - penalty$gamma_mean
            pen <- pen + 0.5*d^2/penalty$gamma_sd^2
          }
        } else if (!is.null(penalty$gamma_shape)) {
          pen <- pen - (penalty$gamma_shape-1)*log(gamma_W) + penalty$gamma_rate*gamma_W
        } else if (!is.null(penalty$gamma_mean)) {
          d <- gamma_W - penalty$gamma_mean
          pen <- pen + 0.5*d^2/penalty$gamma_sd^2
        }
      }
    }
    pen
  }

  nll <- function(par) {
    beta    <- par[1:p]
    k       <- exp(par[p+1])
    rho     <- exp(par[p+2])
    idx     <- p+3
    alpha_W <- if (is.null(fix_alpha_W)) { v <- plogis(par[idx]); idx <<- idx+1; v } else fix_alpha_W
    gamma_W <- if (is.null(fix_gamma_W))   exp(par[idx])                          else fix_gamma_W

    mu_W <- exp(D %*% beta) *
      compute_mda_effect(survey_times_data, mda_times, int_mat,
                         alpha_W, gamma_W, kappa = 1)

    pr <- pmax(pmin(1 - (k/(k + mu_W*(1-exp(-rho))))^k, 1-1e-10), 1e-10)
    ll <- sum(log(1-pr[y_prev==0])) + sum(log(pr[y_prev==1]))

    pos   <- which(y_prev == 1)
    mupos <- mu_W[pos];  prpos <- pr[pos]
    mu_C  <- rho*mupos / prpos
    s2_C  <- pmax((rho*mupos*(1+rho))/prpos +
                    (rho^2*mupos^2/prpos)*(1/k + 1 - 1/prpos), 1e-10)
    kC    <- pmax((mu_C-1)^2/s2_C, 1e-10)
    tC    <- pmax(s2_C/(mu_C-1), 1e-10)
    li    <- dgamma(intensity_data-1, shape=kC, scale=tC, log=TRUE)
    li[!is.finite(li)] <- -1e10
    ll <- ll + sum(li)

    nll <- -(ll - .pen(alpha_W, gamma_W))
    if (!is.finite(nll) || nll > 1e10) 1e10 else nll
  }

  # Starting values
  mc1 <- mean(intensity_data-1, na.rm=TRUE)
  vc1 <- var(intensity_data-1,  na.rm=TRUE)
  k0  <- if (!is.null(start_pars$k)) start_pars$k else
    if (!is.na(vc1) && vc1>0 && mc1>0) max(mc1^2/vc1, 0.1) else 0.5
  rho0 <- if (!is.null(start_pars$rho)) start_pars$rho else 1
  aW0  <- if (!is.null(start_pars$alpha_W)) start_pars$alpha_W else 0.5
  gW0  <- if (!is.null(start_pars$gamma_W)) start_pars$gamma_W else 2.0
  b0   <- if (!is.null(start_pars$beta)) start_pars$beta else {
    b <- rep(0, p)
    b[1] <- log(max(mean(y_prev)*mean(intensity_data, na.rm=TRUE), 1)); b }

  par0 <- c(b0, log(k0), log(rho0))
  if (is.null(fix_alpha_W)) par0 <- c(par0, qlogis(aW0))
  if (is.null(fix_gamma_W)) par0 <- c(par0, log(gW0))

  if (messages) message("Optimising simplified model for initial values (STH)...")
  fit <- nlminb(par0, nll,
                control=list(eval.max=2000, iter.max=1000,
                             trace=ifelse(messages,1,0)))
  pe <- if (fit$convergence != 0) {
    warning("STH initial value optimisation did not converge (code=",
            fit$convergence, "). Using starting values.")
    par0
  } else {
    if (messages) message(sprintf("  Converged: nll = %.2f", fit$objective))
    fit$par
  }

  beta_e <- pe[1:p]; k_e <- exp(pe[p+1]); rho_e <- exp(pe[p+2])
  idx <- p+3
  aW_e <- if (is.null(fix_alpha_W)) { v <- plogis(pe[idx]); idx <- idx+1; v } else fix_alpha_W
  gW_e <- if (is.null(fix_gamma_W))   exp(pe[idx])                          else fix_gamma_W

  s2_0  <- if (!is.null(start_pars$sigma2)) start_pars$sigma2 else 1.0
  phi_0 <- if (!is.null(start_pars$phi)) start_pars$phi else {
    d <- as.matrix(dist(coords)); median(d[upper.tri(d)])/3 }

  if (messages)
    message(sprintf("  k=%.3f  rho=%.3f  alpha_W=%.3f  gamma_W=%.3f  sigma2=%.3f  phi=%.3f",
                    k_e, rho_e, aW_e, gW_e, s2_0, phi_0))

  list(beta=beta_e, k=k_e, rho=rho_e, alpha_W=aW_e, gamma_W=gW_e,
       sigma2=s2_0, phi=phi_0)
}


##' @title Compute initial parameter values for the LF multi-diagnostic model
##' @description Fits a simplified non-spatial binomial model using both MF
##'   (parasitological) and CFA/ICT (serological) rows to obtain starting
##'   values for \code{beta}, \code{k}, \code{rho}, \code{sigma2}, \code{phi},
##'   \code{tau2}, and optionally \code{alpha_W}/\code{gamma_W}.
##' @keywords internal
dsgm_initial_value_lf <- function(y_counts, units_m, which_diag, D, coords,
                                  int_mat, survey_times_data, mda_times,
                                  gamma_sens, penalty,
                                  fix_k, fix_alpha_W, fix_gamma_W,
                                  use_mda, start_pars, messages) {

  n <- length(y_counts)
  p <- ncol(D)

  nll <- function(par) {
    beta <- par[1:p]
    k    <- if (!is.null(fix_k)) fix_k else exp(par[p+1])
    rho  <- exp(par[p+2])
    idx  <- p+3

    alpha_W <- if (use_mda && is.null(fix_alpha_W)) { v <- plogis(par[idx]); idx <<- idx+1; v } else fix_alpha_W
    gamma_W <- if (use_mda && is.null(fix_gamma_W))   exp(par[idx])                          else fix_gamma_W

    mu <- exp(D %*% beta)
    if (use_mda)
      mu <- mu * compute_mda_effect(survey_times_data, mda_times, int_mat,
                                    alpha_W, gamma_W, kappa=1)

    prob <- numeric(n)
    mf  <- which_diag == 1
    cfa <- which_diag == 0
    prob[mf]  <- pmax(pmin(1 - (k/(k + mu[mf] *(1-exp(-rho))))^k, 1-1e-10), 1e-10)
    prob[cfa] <- pmax(pmin(gamma_sens*(1-(k/(k + mu[cfa]))^k),     1-1e-10), 1e-10)

    ll <- sum(y_counts*log(prob) + (units_m-y_counts)*log(1-prob))
    if (!is.finite(ll)) 1e10 else -ll
  }

  # Starting values
  obs_prev <- mean(y_counts / pmax(units_m, 1))
  b0  <- if (!is.null(start_pars$beta)) start_pars$beta else {
    b <- rep(0, p); b[1] <- log(max(obs_prev, 0.01)); b }
  k0  <- if (!is.null(fix_k)) fix_k else
    if (!is.null(start_pars$k)) start_pars$k else 0.5
  rho0 <- if (!is.null(start_pars$rho)) start_pars$rho else 0.5
  aW0  <- if (!is.null(fix_alpha_W)) fix_alpha_W else
    if (!is.null(start_pars$alpha_W)) start_pars$alpha_W else 0.5
  gW0  <- if (!is.null(fix_gamma_W)) fix_gamma_W else
    if (!is.null(start_pars$gamma_W)) start_pars$gamma_W else 2.0

  par0 <- c(b0, log(k0), log(rho0))
  if (use_mda) {
    if (is.null(fix_alpha_W)) par0 <- c(par0, qlogis(aW0))
    if (is.null(fix_gamma_W)) par0 <- c(par0, log(gW0))
  }

  if (messages) message("Optimising simplified model for initial values (LF)...")
  fit <- nlminb(par0, nll,
                control=list(eval.max=2000, iter.max=1000,
                             trace=ifelse(messages,1,0)))
  pe <- if (fit$convergence != 0) {
    warning("LF initial value optimisation did not converge (code=",
            fit$convergence, "). Using starting values.")
    par0
  } else {
    if (messages) message(sprintf("  Converged: nll = %.2f", fit$objective))
    fit$par
  }

  beta_e <- pe[1:p]
  k_e    <- if (!is.null(fix_k)) fix_k else exp(pe[p+1])
  rho_e  <- exp(pe[p+2])
  idx    <- p+3
  aW_e <- if (use_mda && is.null(fix_alpha_W)) { v <- plogis(pe[idx]); idx <- idx+1; v } else aW0
  gW_e <- if (use_mda && is.null(fix_gamma_W))   exp(pe[idx])                          else gW0

  s2_0   <- if (!is.null(start_pars$sigma2)) start_pars$sigma2 else 1.0
  tau2_0 <- if (!is.null(start_pars$tau2))   start_pars$tau2   else 0.1
  phi_0  <- if (!is.null(start_pars$phi)) start_pars$phi else {
    d <- as.matrix(dist(coords)); median(d[upper.tri(d)])/3 }

  if (messages)
    message(sprintf("  k=%.3f  rho=%.3f  sigma2=%.3f  phi=%.3f  tau2=%.3f",
                    k_e, rho_e, s2_0, phi_0, tau2_0))

  list(beta=beta_e, k=k_e, rho=rho_e, gamma_sens=gamma_sens,
       sigma2=s2_0, phi=phi_0, tau2=tau2_0, alpha_W=aW_e, gamma_W=gW_e)
}


# =============================================================================
# Main user-facing function
# =============================================================================

##' @title Fit a Doubly Stochastic Geostatistical Model (DSGM)
##'
##' @description
##' Fits a doubly stochastic geostatistical model using Monte Carlo Maximum
##' Likelihood (MCML).  Two model families are supported via \code{model}:
##'
##' \describe{
##'   \item{\code{"sth"}}{Joint prevalence-intensity model for
##'     soil-transmitted helminths or intestinal schistosomiasis.  The response
##'     variable is an egg count (EPG); zeros indicate uninfected individuals.
##'     Latent worm burden \eqn{W \sim \mathrm{NegBin}(\mu(x),\,k)} drives
##'     EPG via \eqn{Y \mid W \sim \mathrm{Poisson}(\rho W)}.}
##'   \item{\code{"lf_mdiag"}}{Joint model for two diagnostic outcomes for
##'     lymphatic filariasis.  A parasitological count diagnostic (MF) and a
##'     binary serological diagnostic (CFA/ICT) are modelled simultaneously
##'     through the same NB latent worm burden.}
##' }
##'
##' Both models use Stan (\code{dsgm_spatial.stan} for STH,
##' \code{dsgm_mdiag.stan} for LF) to sample S(x) at fixed theta_0, then
##' TMB to optimise the MCML objective.
##'
##' @param formula Model formula.  For \code{"sth"}, the response is the egg
##'   count (0 = uninfected).  For \code{"lf_mdiag"}, the response is ignored;
##'   supply diagnostic data via the formula response and \code{diagnostic} column.  Both
##'   the \code{diagnostic} column.  Both models require a spatial GP term: \code{gp(x,y)} or
##'   \code{gp(sf)}.
##' @param data A \code{data.frame} or \code{sf} object.
##' @param model \code{"sth"} (default) or \code{"lf_mdiag"}.
##'
##' @section STH-specific arguments:
##' \describe{
##'   \item{time}{Column in \code{data} giving survey time per observation.
##'     Required for \code{"sth"}.}
##'   \item{mda_times}{Numeric vector of MDA round times.}
##'   \item{int_mat}{Coverage matrix (n \eqn{\times} n_mda).}
##' }
##'
##' @param den Name of a column in \code{data} containing the binomial
##'   denominator (number of individuals examined) for each observation.
##'   For \code{model = "lf_mdiag"} this is the number tested per row.
##'   If not supplied, defaults to 1 for every observation (individual-level
##'   data, e.g. standard STH surveys).
##'
##' @section LF-specific arguments:
##' \describe{
##'   \item{which_diag}{Integer 0/1 vector derived from the \code{diagnostic} column: 1 = parasitological (\code{"par"}), 0 = serological (\code{"ser"}).}
##'   \item{gamma_sens}{Fixed serological sensitivity in (0,1]; default 0.97.}
##'   \item{fix_k}{Fix the NB aggregation parameter k; \code{NULL} to estimate.}
##'   \item{use_mda}{Logical; include MDA decay for LF (default \code{FALSE}).
##'     When \code{TRUE} also supply \code{time}, \code{mda_times}, \code{int_mat}.}
##' }
##'
##' @param penalty Optional list of penalty functions for MDA parameters.
##' @param drop_W Fixed value for the immediate worm reduction alpha_W
##'   (\code{NULL} to estimate).
##' @param decay_W Fixed value for the worm recovery rate gamma_W
##'   (\code{NULL} to estimate).
##' @param crs Optional CRS for spatial data.
##' @param convert_to_crs CRS to project data to before fitting.
##' @param scale_to_km Scale coordinates to km (default \code{TRUE}).
##' @param par0 Optional named list of initial parameter values.  If
##'   \code{NULL}, computed automatically.
##' @param n_samples Number of Stan MCMC samples (default 1000).
##' @param n_warmup Number of Stan warmup iterations (default 1000).
##' @param n_chains Number of Stan chains (default 1).
##' @param adapt_delta Target Stan acceptance rate (default 0.8).
##' @param max_treedepth Stan tree depth limit (default 10).
##' @param return_samples Include spatial samples in result (default \code{TRUE}).
##' @param backend Stan backend: \code{"rstan"}, \code{"cmdstanr"}, or
##'   \code{NULL} (auto-detect).
##' @param messages Print progress messages (default \code{TRUE}).
##' @param start_pars Named list of starting values for auto-initialisation
##'   (\code{beta}, \code{k}, \code{rho}, \code{sigma2}, \code{phi},
##'   \code{tau2}, \code{alpha_W}, \code{gamma_W}).
##'
##' @return An object of class \code{"RiskMap"} with fields \code{family}
##'   (\code{"intprev"} or \code{"lf_mdiag"}), \code{model_params},
##'   \code{params_se}, \code{tmb_sdr}, and optionally
##'   \code{spatial_samples}.
##'
##' @seealso \code{\link{dast}}
##' @export
dsgm <- function(formula,
                 data,
                 model          = c("sth", "lf_mdiag"),
                 # STH / shared MDA
                 time           = NULL,
                 mda_times      = NULL,
                 int_mat        = NULL,
                 # LF-specific
                 den            = NULL,
                 gamma_sens     = 0.97,
                 fix_k          = NULL,
                 use_mda        = NULL,
                 # Shared
                 penalty        = NULL,
                 drop_W         = NULL,
                 decay_W        = NULL,
                 crs            = NULL,
                 convert_to_crs = NULL,
                 scale_to_km    = TRUE,
                 par0           = NULL,
                 n_samples      = 1000,
                 n_warmup       = 1000,
                 n_chains       = 1,
                 adapt_delta    = 0.8,
                 max_treedepth  = 10,
                 return_samples = TRUE,
                 backend        = NULL,
                 messages       = TRUE,
                 start_pars     = list(beta    = NULL,
                                       k       = NULL,
                                       rho     = NULL,
                                       sigma2  = NULL,
                                       phi     = NULL,
                                       tau2    = NULL,
                                       alpha_W = NULL,
                                       gamma_W = NULL)) {

  model <- match.arg(model)

  # ---------------------------------------------------------------------------
  # Common validation
  # ---------------------------------------------------------------------------
  if (!inherits(formula, "formula"))
    stop("'formula' must be a formula object")
  if (!inherits(data, c("data.frame", "sf")))
    stop("'data' must be a data.frame or sf object")
  if (n_samples <= 0 || n_warmup < 0)
    stop("'n_samples' must be positive and 'n_warmup' non-negative")
  if (n_chains < 1)
    stop("'n_chains' must be at least 1")
  if (adapt_delta <= 0 || adapt_delta >= 1)
    stop("'adapt_delta' must be in (0, 1)")

  # Default use_mda: enable only when mda_times have been supplied
  if (is.null(use_mda)) use_mda <- !is.null(mda_times)

  fix_alpha_W <- drop_W
  fix_gamma_W <- decay_W
  no_penalty  <- is.null(penalty)
  n           <- nrow(data)

  # ---------------------------------------------------------------------------
  # Formula / design matrix
  # ---------------------------------------------------------------------------
  inter_f  <- interpret.formula(formula)
  gp_terms <- inter_f$gp.spec$term
  gp_dim   <- inter_f$gp.spec$dim

  if (!(length(gp_terms) == 1 && gp_terms[1] == "sf") && gp_dim != 2)
    stop("Specify a 2-D spatial GP: gp(x, y) or gp(sf)")

  # ---------------------------------------------------------------------------
  # tau2 (nugget) policy driven by the gp() specification:
  #   gp(...)               nugget = 0   (default) -> fix tau2 = 0, not estimated
  #   gp(..., nugget = v)   nugget = v > 0          -> fix tau2 = v, not estimated
  #   gp(..., nugget = NULL) nugget = NULL           -> tau2 is estimated freely
  # ---------------------------------------------------------------------------
  gp_nugget <- inter_f$gp.spec$nugget
  fix_tau2  <- if (is.null(gp_nugget)) NULL else as.numeric(gp_nugget)

  mf         <- model.frame(inter_f$pf, data = data, na.action = na.fail)
  D          <- as.matrix(model.matrix(attr(mf, "terms"), data = data))
  p          <- ncol(D)
  cov_offset <- if (is.null(inter_f$offset)) rep(0, n) else data[[inter_f$offset]]

  # ---------------------------------------------------------------------------
  # Denominator (number examined per observation)
  # ---------------------------------------------------------------------------
  den_name <- deparse(substitute(den))
  if (den_name == "NULL") {
    den_vals <- rep(1L, n)
  } else {
    if (!den_name %in% names(data))
      stop(sprintf("'den' column '%s' not found in 'data'", den_name))
    den_vals <- as.integer(data[[den_name]])
    if (any(den_vals < 1, na.rm = TRUE))
      stop("'den' values must all be >= 1")
    if (any(is.na(den_vals)))
      stop("Missing values detected in 'den' column")
  }

  # ---------------------------------------------------------------------------
  # Coordinates
  # ---------------------------------------------------------------------------
  if (inherits(data, "sf")) {
    if (!is.null(convert_to_crs)) {
      data <- sf::st_transform(data, convert_to_crs); crs <- convert_to_crs
    } else {
      crs <- sf::st_crs(data)$input
    }
    coords_all <- sf::st_coordinates(data)
  } else {
    cn <- gp_terms[1:2]
    if (!all(cn %in% names(data)))
      stop("Coordinate columns specified in gp() not found in 'data'")
    coords_all <- as.matrix(data[, cn])
  }
  if (scale_to_km) {
    coords_all <- coords_all / 1000
    if (messages) message("Coordinates scaled to kilometres")
  }
  coords_u  <- unique(coords_all)
  n_loc     <- nrow(coords_u)
  ID_coords <- apply(coords_all, 1, function(x)
    which(apply(coords_u, 1, function(y) all(abs(x-y) < 1e-10)))[1])
  if (messages) message(sprintf("Identified %d unique spatial locations", n_loc))

  # ---------------------------------------------------------------------------
  # MDA / time setup
  # ---------------------------------------------------------------------------
  if (use_mda) {
    time_name <- deparse(substitute(time))
    if (time_name == "NULL" || !time_name %in% names(data))
      stop("'time' column not found in 'data' (required when use_mda = TRUE)")
    survey_times_data <- data[[time_name]]
    if (any(is.na(survey_times_data))) stop("Missing values in 'time' variable")
    if (is.null(mda_times)) stop("'mda_times' required when use_mda = TRUE")
    if (is.null(int_mat))   stop("'int_mat' required when use_mda = TRUE")
    int_mat <- as.matrix(int_mat)
    if (nrow(int_mat) != n)
      stop(sprintf("'int_mat' must have %d rows", n))
    if (ncol(int_mat) != length(mda_times))
      stop(sprintf("'int_mat' must have %d columns", length(mda_times)))
  } else {
    survey_times_data <- rep(0, n)
    mda_times         <- 0
    int_mat           <- matrix(0, nrow=n, ncol=1)
  }

  # ===========================================================================
  # STH branch
  # ===========================================================================
  if (model == "sth") {

    egg_counts <- as.numeric(model.response(mf))
    if (any(egg_counts < 0, na.rm=TRUE)) stop("Egg counts cannot be negative")
    if (any(is.na(egg_counts))) stop("Missing values in egg count response")

    y_prev         <- as.integer(egg_counts > 0)
    intensity_data <- egg_counts[egg_counts > 0]
    n_pos          <- sum(y_prev)
    if (n_pos == 0) stop("No positive egg counts; model cannot be fitted")

    if (messages)
      message(sprintf("STH data: %d observations, %d (%.1f%%) egg-positive",
                      n, n_pos, 100*n_pos/n))

    if (is.null(par0)) {
      if (messages) message("\n=== Computing initial parameter values (STH) ===")
      par0 <- dsgm_initial_value(
        y_prev=y_prev, intensity_data=intensity_data, D=D,
        coords=coords_u, ID_coords=ID_coords,
        int_mat=int_mat, survey_times_data=survey_times_data,
        mda_times=mda_times, penalty=penalty,
        fix_alpha_W=fix_alpha_W, fix_gamma_W=fix_gamma_W,
        start_pars=start_pars, messages=messages)
    }
    req  <- c("beta","k","rho","sigma2","phi")
    if (is.null(fix_alpha_W)) req <- c(req, "alpha_W")
    if (is.null(fix_gamma_W)) req <- c(req, "gamma_W")
    miss <- setdiff(req, names(par0))
    if (length(miss) > 0)
      stop("Missing initial parameters: ", paste(miss, collapse=", "))

    if (messages) {
      message("\n=== Sampling spatial process (STH, Stan) ===")
      message(sprintf("  n_samples=%d  n_warmup=%d  n_chains=%d  adapt_delta=%.2f",
                      n_samples, n_warmup, n_chains, adapt_delta))
    }
    sp <- sample_spatial_process_stan(
      y_prev=y_prev, intensity_data=intensity_data, D=D,
      coords=coords_u, ID_coords=ID_coords,
      int_mat=int_mat, survey_times_data=survey_times_data,
      mda_times=mda_times, par=par0,
      n_samples=n_samples, n_warmup=n_warmup,
      n_chains=n_chains, n_cores=1,
      adapt_delta=adapt_delta, max_treedepth=max_treedepth,
      backend=backend, messages=messages)

    if (messages) message("\n=== Fitting via TMB-MCML (STH) ===")
    fit <- dsgm_fit_tmb(
      model="sth",
      y_prev=y_prev, intensity_data=intensity_data, D=D,
      coords=coords_u, ID_coords=ID_coords,
      int_mat=int_mat, survey_times_data=survey_times_data,
      mda_times=mda_times, par0=par0, cov_offset=cov_offset,
      fix_alpha_W=fix_alpha_W, fix_gamma_W=fix_gamma_W,
      fix_tau2=fix_tau2, use_mda=use_mda,
      penalty=penalty, S_samples_obj=sp,
      use_hessian_refinement=TRUE, messages=messages)

    res <- list(
      family            = "intprev",
      prevalence_data   = y_prev,
      intensity_data    = intensity_data,
      egg_counts        = egg_counts,
      n_positive        = n_pos,
      D                 = D,
      coords            = coords_u,
      mda_times         = mda_times,
      survey_times_data = survey_times_data,
      int_mat           = int_mat,
      ID_coords         = ID_coords,
      fix_alpha_W       = fix_alpha_W,
      fix_gamma_W       = fix_gamma_W,
      use_mda           = use_mda,
      formula           = formula,
      crs               = crs,
      scale_to_km       = scale_to_km,
      data_sf           = data,
      kappa             = 0.5,
      cov_offset        = cov_offset,
      penalty           = if (!no_penalty) penalty else NULL,
      model_params      = fit$params,
      params_se         = fit$params_se,
      tmb_sdr           = fit$tmb_sdr,
      tmb_obj           = fit$tmb_obj,
      convergence       = fit$convergence,
      log_likelihood    = fit$log_likelihood,
      n_locations       = n_loc,
      n_observations    = n,
      n_covariates      = p,
      n_mda_rounds      = if (use_mda) length(mda_times) else 0L,
      stan_settings     = list(n_samples=n_samples, n_warmup=n_warmup,
                               n_chains=n_chains, adapt_delta=adapt_delta,
                               max_treedepth=max_treedepth),
      call              = match.call()
    )
    if (return_samples) res$spatial_samples <- sp
    class(res) <- "RiskMap"
    return(res)
  }

  # ===========================================================================
  # LF multi-diagnostic branch
  # ===========================================================================
  if (model == "lf_mdiag") {

    # Derive which_diag from the 'diagnostic' column in data.
    # If the column is absent, assume all observations are parasitological
    # (single-diagnostic MF model).
    if ("diagnostic" %in% names(data)) {
      diag_raw <- data[["diagnostic"]]
      if (!all(diag_raw %in% c("par", "ser")))
        stop("Column 'diagnostic' must contain only 'par' (parasitological) or 'ser' (serological).")
      which_diag <- ifelse(diag_raw == "par", 1L, 0L)
    } else {
      if (messages)
        message("No 'diagnostic' column found: assuming all observations are parasitological (MF).")
      which_diag <- rep(1L, n)
    }
    if (gamma_sens <= 0 || gamma_sens > 1)
      stop("'gamma_sens' must be in (0, 1]")

    # Extract positive counts from formula response
    y_counts <- as.integer(model.response(mf))
    if (any(y_counts < 0, na.rm = TRUE)) stop("Response (positive counts) cannot be negative")
    if (any(is.na(y_counts)))            stop("Missing values in response variable")

    # Denominator: number examined per row (from 'den', or 1 if not supplied)
    units_m <- den_vals

    if (messages)
      message(sprintf("LF data: %d observations (%d MF, %d serological)",
                      n, sum(which_diag==1), sum(which_diag==0)))

    if (is.null(par0)) {
      if (messages) message("\n=== Computing initial parameter values (LF) ===")
      par0 <- dsgm_initial_value_lf(
        y_counts=y_counts, units_m=units_m, which_diag=which_diag, D=D,
        coords=coords_u,
        int_mat=int_mat, survey_times_data=survey_times_data,
        mda_times=mda_times, gamma_sens=gamma_sens, penalty=penalty,
        fix_k=fix_k, fix_alpha_W=fix_alpha_W, fix_gamma_W=fix_gamma_W,
        use_mda=use_mda, start_pars=start_pars, messages=messages)
    }
    # Ensure gamma_sens is set; tau2 follows gp() nugget specification
    par0$gamma_sens <- gamma_sens
    if (!is.null(fix_tau2)) {
      # nugget fixed (default is 0): override whatever initial value gave us
      par0$tau2 <- fix_tau2
    } else if (is.null(par0$tau2)) {
      # nugget estimated: use start_pars or a small default
      par0$tau2 <- if (!is.null(start_pars$tau2)) start_pars$tau2 else 0.1
    }

    # Pre-compute MDA impact at theta_0 for Stan
    mda_impact <- if (use_mda)
      compute_mda_effect(survey_times_data, mda_times, int_mat,
                         par0$alpha_W, par0$gamma_W, kappa=1)
    else
      rep(1.0, n)

    if (messages) {
      message("\n=== Sampling spatial process (LF, Stan) ===")
      message(sprintf("  n_samples=%d  n_warmup=%d  n_chains=%d  adapt_delta=%.2f",
                      n_samples, n_warmup, n_chains, adapt_delta))
    }
    sp <- sample_spatial_process_stan_lf(
      y_counts=y_counts, units_m=units_m, which_diag=which_diag, D=D,
      coords=coords_u, ID_coords=ID_coords, par=par0,
      mda_impact=mda_impact,
      n_samples=n_samples, n_warmup=n_warmup,
      n_chains=n_chains, n_cores=1,
      adapt_delta=adapt_delta, max_treedepth=max_treedepth,
      backend=backend, messages=messages)

    if (messages) message("\n=== Fitting via TMB-MCML (LF) ===")
    fit <- dsgm_fit_tmb(
      model="lf_mdiag",
      y_counts=y_counts, units_m=units_m, which_diag=which_diag, D=D,
      coords=coords_u, ID_coords=ID_coords,
      int_mat=int_mat, survey_times_data=survey_times_data,
      mda_times=mda_times, par0=par0, cov_offset=cov_offset,
      gamma_sens=gamma_sens, fix_k=fix_k, use_mda=use_mda,
      fix_alpha_W=fix_alpha_W, fix_gamma_W=fix_gamma_W,
      fix_tau2=fix_tau2,
      penalty=penalty, S_samples_obj=sp,
      use_hessian_refinement=TRUE, messages=messages)

    res <- list(
      family            = "lf_mdiag",
      y_counts          = y_counts,
      units_m           = units_m,
      which_diag        = which_diag,
      D                 = D,
      coords            = coords_u,
      mda_times         = mda_times,
      survey_times_data = survey_times_data,
      int_mat           = int_mat,
      ID_coords         = ID_coords,
      fix_alpha_W       = fix_alpha_W,
      fix_gamma_W       = fix_gamma_W,
      fix_k             = fix_k,
      gamma_sens        = gamma_sens,
      use_mda           = use_mda,
      formula           = formula,
      crs               = crs,
      scale_to_km       = scale_to_km,
      data_sf           = data,
      kappa             = 0.5,
      cov_offset        = cov_offset,
      penalty           = if (!no_penalty) penalty else NULL,
      model_params      = fit$params,
      params_se         = fit$params_se,
      tmb_sdr           = fit$tmb_sdr,
      tmb_obj           = fit$tmb_obj,
      convergence       = fit$convergence,
      log_likelihood    = fit$log_likelihood,
      n_locations       = n_loc,
      n_observations    = n,
      n_covariates      = p,
      n_mda_rounds      = length(mda_times),
      stan_settings     = list(n_samples=n_samples, n_warmup=n_warmup,
                               n_chains=n_chains, adapt_delta=adapt_delta,
                               max_treedepth=max_treedepth),
      call              = match.call()
    )
    if (return_samples) res$spatial_samples <- sp
    class(res) <- "RiskMap"
    return(res)
  }
}
