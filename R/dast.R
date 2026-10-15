compute_mda_effect <- function(survey_times_data, mda_times, intervention,
                               alpha, gamma, kappa) {
  n <- length(survey_times_data)
  effect <- rep(NA, n)

  mda_effect_f <- function(v, alpha, gamma, kappa) {
    alpha*exp(-(v/gamma)^kappa)
  }
  mda_effect_f <- Vectorize(mda_effect_f, "v")

  f <- function(t, mda_times, int, alpha, gamma, kappa) {
    ind_t <- which(t > mda_times)
    u_j <- mda_times[ind_t]
    if(length(u_j) > 0) {
      out <- prod((1-mda_effect_f(t-u_j, alpha, gamma, kappa))^int[ind_t])
    } else {
      out <- 1
    }
    return(out)
  }

  f <- Vectorize(f, "t")

  for(i in 1:n) {
    effect[i] <- f(survey_times_data[i], mda_times, intervention[i,],
                   alpha, gamma, kappa)
  }
  return(effect)
}


compute_mda_effect_derivatives <- function(survey_times_data, mda_times, intervention,
                                           alpha, gamma, kappa) {
  n <- length(survey_times_data)
  effect <- numeric(n)
  d_alpha <- numeric(n)
  d_gamma <- numeric(n)
  d2_alpha <- numeric(n)
  d2_gamma <- numeric(n)
  d2_alpha_gamma <- numeric(n)

  for (i in 1:n) {
    t <- survey_times_data[i]
    int <- intervention[i,]

    ind_t <- which(t > mda_times)
    u_j <- mda_times[ind_t]
    int_j <- int[ind_t]

    if (length(u_j) > 0) {
      v_j <- t - u_j
      z_j <- (v_j/gamma)^kappa
      e_j <- exp(-z_j)
      denom <- 1 - alpha * e_j
      effect[i] <- exp(sum(int_j * log(denom)))

      # First derivatives (CORRECT - matches numeric)
      term_alpha <- int_j * e_j / denom
      dlogP_dalpha <- -sum(term_alpha)

      gamma_term <- kappa * z_j / gamma
      term_gamma <- alpha * term_alpha * gamma_term
      dlogP_dgamma <- -sum(term_gamma)

      # Second derivatives
      # ∂²(log P)/∂α² (CORRECT - matches numeric)
      d2logP_dalpha2 <- -sum(int_j * e_j^2 / denom^2)

      #######################################################
      # CORRECTED ∂²(log P)/∂γ² (FIXED SIGN)
      # f = αe^{-z}, z = (v/γ)^κ
      # ∂f/∂γ = f * (κz/γ)
      # ∂²f/∂γ² = f * [κz(κ+1 - κz)/γ²]
      term1 <- (alpha * e_j * kappa * z_j * (kappa + 1 - kappa*z_j)) /
        (gamma^2 * denom)
      term2 <- (alpha * e_j * kappa * z_j / (gamma * denom))^2
      d2logP_dgamma2 <- sum(int_j * (term1 - term2))
      #######################################################

      # ∂²(log P)/∂α∂γ (CORRECT - matches numeric)
      term_cross <- e_j * gamma_term / denom
      d2logP_dalphadgamma <- sum(int_j * (-term_cross - alpha * e_j * term_cross / denom))

      # Convert to derivatives of P
      d_alpha[i] <- effect[i] * dlogP_dalpha
      d_gamma[i] <- effect[i] * dlogP_dgamma
      d2_alpha[i] <- effect[i] * (dlogP_dalpha^2 + d2logP_dalpha2)
      d2_gamma[i] <- effect[i] * (dlogP_dgamma^2 + d2logP_dgamma2)
      d2_alpha_gamma[i] <- effect[i] * (dlogP_dalpha * dlogP_dgamma + d2logP_dalphadgamma)

    } else {
      effect[i] <- 1
      d_alpha[i] <- 0
      d_gamma[i] <- 0
      d2_alpha[i] <- 0
      d2_gamma[i] <- 0
      d2_alpha_gamma[i] <- 0
    }
  }

  return(list(effect = effect,
              d_alpha = d_alpha,
              d_gamma = d_gamma,
              d2_alpha = d2_alpha,
              d2_gamma = d2_gamma,
              d2_alpha_gamma = d2_alpha_gamma))
}

dast_initial_value <- function(y, D, units_m, int_mat, survey_times_data,
                               penalty,
                               mda_times, fix_alpha, power_val) {

  p <- ncol(D)
  n <- nrow(D)

  llik <- function(par) {
    beta <- par[1:p]
    if(is.null(fix_alpha)) {
      alpha <- exp(par[p+1])/(1+exp(par[p+1]))
      gamma <- exp(par[p+2])
    } else {
      alpha <- fix_alpha
      gamma <- exp(par[p+1])
    }

    fact <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                               alpha, gamma, kappa = power_val)
    eta <- as.numeric(D%*%beta)
    prob_star <- 1/(1+exp(-eta))
    prob <- fact*prob_star

    # Apply penalty for both alpha and gamma (penalty[[4]] is the gamma penalty)
    out <- -(sum(y*log(prob/(1-prob)) + units_m*log(1-prob)) -
               penalty[[1]](alpha) - penalty[[4]](gamma))

    return(out)
  }

  start <- c(runif(p, -1, 1),
             0,
             log(mean(dist(survey_times_data))))
  est <- nlminb(start, llik)

  est$beta <- est$par[1:p]
  if(is.null(fix_alpha)) {
    est$alpha <- exp(est$par[p+1])/(1+exp(est$par[p+1]))
    est$gamma <- exp(est$par[p+2])
  } else {
    est$gamma <- exp(est$par[p+1])
  }

  return(est)
}

##' @title Check that dast()-specific arguments are valid
##' @description Validates the arguments to \code{dast()} that are not
##' already covered by \code{prepare_data()}/\code{check_formula()}: the
##' GP structure (purely spatial \code{gp()}, or spatio-temporal
##' \code{gp(t)}), the \code{time} column, the dimensions of
##' \code{int_mat}, the shape of \code{penalty}, and the starting value
##' for \code{psi}.
##' @param formula,data,mda_times,int_mat,penalty,start_pars Identical to
##' the corresponding arguments of \code{dast()}.
##' @param time The caller's substituted \code{time} expression, e.g.
##' \code{substitute(time)} from \code{dast()}.
##' @return \code{TRUE} if the dast()-specific arguments are valid, or
##' raises an error if not.
##' @noRd
check_dast <- function(formula, data, time, mda_times, int_mat, penalty, start_pars) {

  inter_f <- interpret_formula(formula)
  gp_terms <- inter_f$gp_spec$term
  gp_dim   <- inter_f$gp_spec$dim

  if (!(length(gp_terms) == 1 && gp_terms[1] == "sf") && gp_dim != 1) {
    stop("Specify gp() for a purely spatial model, or gp(t) for a spatio-temporal model, where 't' is the GP time index column.")
  }

  time_name <- deparse(time)
  if (time_name == "NULL") stop("You must supply a time column via the 'time' argument")
  if (!time_name %in% names(data)) stop("'time' column not found in 'data'")

  if (is.null(int_mat)) stop("'int_mat' must be provided")
  n <- nrow(data)
  if (nrow(as.matrix(int_mat)) != n)
    stop(sprintf("'int_mat' must have %d rows", n))
  if (ncol(as.matrix(int_mat)) != length(mda_times))
    stop(sprintf("'int_mat' must have %d columns", length(mda_times)))

  if (!is.null(penalty)) {
    if (!is.list(penalty) || length(penalty) != 6 ||
        !all(vapply(penalty, is.function, logical(1))))
      stop("'penalty' must be a list of 6 functions: penalty/d1/d2 for alpha (elements 1-3) ",
           "and penalty/d1/d2 for gamma (elements 4-6).")
  }

  if (!is.null(start_pars[["psi"]])) check_positive_number(start_pars[["psi"]])

  invisible(TRUE)
}

##' @title Fitting of decay-adjusted spatio-temporal (DAST) model
##'
##' @description
##' The function fits a decay-adjusted spatio-temporal (DAST) model using Monte Carlo maximum likelihood.
##' The DAST model allows for the incorporation of temporal decay in disease prevalence due
##' to the impact of mass drug administration (MDA). The function requires the full MDA history as detailed in the arguments below.
##'
##' Spatial and spatio-temporal dependence is specified through the \code{gp()} term in the model formula:
##' \itemize{
##'   \item \code{gp()} fits a purely spatial Gaussian process, using the coordinates of \code{data}.
##'   \item \code{gp(t_gp)} fits a spatio-temporal Gaussian process, where \code{t_gp} is a column of
##'   \code{data} used as the GP temporal index.
##' }
##' In all cases, the \code{time} argument must be specified separately and provides the observation-level survey times used
##' in modelling MDA impact. These survey times may differ from the GP temporal index.
##'
##' @param formula A formula object specifying the model to be fitted.
##' The formula should include fixed effects and spatial effects specified using
##' [gp()] and optionally, random effects specified using [re()]. To include a
##' spatiotemporal process, add the time variable to `gp()`
##' @param data An sf object containing the variables in the model.
##' @param denominator The denominator for binomial models.
##' @param time A variable in \code{data} giving the survey times of observations (required).
##' @param mda_times A vector specifying the mass drug administration (MDA) times.
##' @param int_mat Intervention matrix specifying the timing and coverage of MDA; the dimension of the matrix
##' must be \code{n * n_mda}, where \code{n} is the number of rows of \code{data} and \code{n_mda} is the length of \code{mda_times}.
##' @param power_val Value expressing the power of the MDA impact function.
##' @param penalty Optional list of length 6 specifying penalty functions for regularization.
##'   Elements \code{[[1]]}, \code{[[2]]}, \code{[[3]]} are the penalty, its first and second
##'   derivative with respect to \code{alpha}. Elements \code{[[4]]}, \code{[[5]]}, \code{[[6]]}
##'   are the penalty, its first and second derivative with respect to \code{gamma}.
##'   Any of these six elements may be set to \code{function(x) 0} to suppress penalisation
##'   for that parameter.
##' @param drop Optional value used for fixing the "drop" parameter of the MDA impact function.
##' @param model_crs Optional CRS (e.g. an EPSG code) used internally for model fitting.
##' If `data` is already in a projected CRS and `model_crs` is not provided, the data are used as-is.
##' If `data` are in longitude/latitude and `model_crs` is not provided, the data are automatically
##' reprojected to an appropriate UTM zone (see [propose_utm()]) and a message reports the conversion used.
##' If provided, `model_crs` must be a projected (not longitude/latitude) CRS.
##' @param distance_units Units in which to compute distances; either \code{"km"} (default) or \code{"m"}.
##' @param control_mcmc A list of MCMC control parameters, typically from \code{set_control_mcmc()}.
##' @param S_samples Number of posterior samples to retain.
##' @param return_samples Logical; whether to return posterior samples (default: \code{TRUE}).
##' @param messages Logical; whether to print messages (default: \code{TRUE}).
##' @param start_pars List of starting values for parameters.
##'
##' @return A list containing model estimates, posterior samples, and metadata, including:
##' \itemize{
##'   \item \code{y}: Response variable values.
##'   \item \code{D}: Covariate matrix.
##'   \item \code{coords}: Unique spatial coordinates.
##'   \item \code{mda_times}: MDA time points.
##'   \item \code{survey_times_data}: Survey time data from the \code{time} argument.
##'   \item \code{time}: GP temporal index if specified in \code{gp(t_gp)}.
##'   \item \code{int_mat}: Intervention matrix.
##'   \item \code{ID_coords}: Indices of spatial locations (and time if spatio-temporal GP).
##'   \item \code{re}: Random effects levels (if applicable).
##'   \item \code{ID_re}: Indices of random effects (if applicable).
##'   \item \code{power_val}: Power of the MDA impact function.
##'   \item \code{fix_tau2}: Fixed tau-squared value (if applicable).
##'   \item \code{fix_alpha}: Fixed alpha value (if applicable).
##'   \item \code{formula}: Model formula.
##'   \item \code{crs}: Coordinate reference system.
##'   \item \code{distance_units}: Units used for distance computations.
##'   \item \code{data_sf}: Processed spatial dataset.
##'   \item \code{family}: Model family (e.g., "binomial").
##'   \item \code{sst}: Logical indicator of whether a spatio-temporal GP was used.
##'   \item \code{kappa}: Smoothness parameter.
##'   \item \code{units_m}: Denominator for binomial models.
##'   \item \code{cov_offset}: Offset for covariates.
##'   \item \code{call}: Function call.
##'   \item \code{penalty}: Penalty function details (if applicable).
##'   \item \code{posterior_samples}: Posterior samples if \code{return_samples = TRUE}.
##' }
##'
##' @seealso \code{\link{set_control_mcmc}}, \code{\link{summary.RiskMap}}, \code{\link{to_table}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterre@@lancaster.ac.uk}
##' @export
dast <- function(formula,
                 data,
                 denominator = NULL,
                 time,
                 mda_times,
                 int_mat,
                 power_val,
                 penalty = NULL,
                 drop = NULL,
                 model_crs = NULL,
                 distance_units = c("km", "m"),
                 control_mcmc = set_control_mcmc(),
                 S_samples = NULL,
                 return_samples = TRUE,
                 messages = TRUE,
                 start_pars = list(beta = NULL,
                                   sigma2 = NULL,
                                   tau2 = NULL,
                                   phi = NULL,
                                   psi = NULL,
                                   sigma2_re = NULL,
                                   gamma = NULL,
                                   alpha = NULL)) {

  # only pass generic start_pars
  dast_start_pars <- start_pars[intersect(names(start_pars), c("gamma", "alpha", "psi"))]
  standard_start_pars <- start_pars[intersect(names(start_pars), c("beta", "sigma2", "tau2", "phi", "sigma2_re"))]

  prepared <- prepare_data(
    formula = formula, data = data, family = "binomial",
    denominator = substitute(denominator), model_crs = model_crs,
    distance_units = distance_units, control_mcmc = control_mcmc,
    return_samples = return_samples, messages = messages,
    start_pars = standard_start_pars
  )
  list2env(prepared, environment())
  start_pars[names(dast_start_pars)] <- dast_start_pars

  check_dast(formula = formula, data = data, time = substitute(time),
             mda_times = mda_times, int_mat = int_mat, penalty = penalty,
             start_pars = start_pars)

  # --- dast-specific preparation: not handled by prepare_data() ---

  gp_terms <- inter_f$gp_spec$term
  gp_dim   <- inter_f$gp_spec$dim

  # time= (required, checked by check_dast() above) -> survey_times_data
  time_name <- deparse(substitute(time))
  survey_times_data <- data[[time_name]]

  gp_time_obs <- NULL
  if (gp_dim == 1) {
    sst <- TRUE
    gp_time_obs <- data[[gp_terms[1]]]   # may differ from survey_times_data
  } else {
    sst <- FALSE
  }

  # Build default no-penalty list.
  # [[1]]-[[3]]: alpha penalty and its first/second derivatives.
  # [[4]]-[[6]]: gamma penalty and its first/second derivatives.
  if(is.null(penalty)) {
    no_penalty <- TRUE
    penalty <- list(pn        = function(x) return(0),
                    pn_d1     = function(x) return(0),
                    pn_d2     = function(x) return(0),
                    pn_gamma    = function(x) return(0),
                    pn_gamma_d1 = function(x) return(0),
                    pn_gamma_d2 = function(x) return(0))
  } else {
    no_penalty <- FALSE
  }

  n <- length(y)
  int_mat <- as.matrix(int_mat)

  if(!is.null(drop)) {
    fix_alpha <- drop
  } else {
    fix_alpha <- NULL
  }

  # deduplicate on x,y,t
  if (sst) {
    coords_o <- coordinates_in_units(data, distance_units)
    coords_time <- unique(cbind(coords_o, gp_time_obs))
    coords <- coords_time[, 1:2, drop = FALSE]
    time_gp <- coords_time[, 3]
    m <- nrow(coords_o)
    ID_coords <- sapply(seq_len(m), function(i)
      which(coords_o[i, 1] == coords[, 1] &
              coords_o[i, 2] == coords[, 2] &
              gp_time_obs[i] == time_gp))
    s_unique <- unique(ID_coords)
  } else {
    time_gp <- NULL
  }

  if(is.null(start_pars$gamma) | (is.null(drop) & is.null(start_pars$alpha)) ) {
    dast_i <- dast_initial_value(y, D, units_m, int_mat = int_mat, survey_times_data,
                                 fix_alpha = fix_alpha,penalty=penalty,
                                 mda_times, power_val = power_val)
    start_pars$beta <- dast_i$beta
    if(is.null(drop)) {
      start_pars$alpha <- dast_i$alpha
    }
    start_pars$gamma <- dast_i$gamma
  }

  if(sst) {
    if(is.null(start_pars$psi)) {
      # base distances from GP time index (not survey time)
      psi_base <- as.numeric(unique(gp_time_obs))
      if (length(psi_base) > 1) {
        start_pars$psi <- quantile(dist(psi_base), 0.1)
      } else {
        start_pars$psi <- 1
      }
    }
  }

  res <- dast_fit(y = y, D = D, coords = coords, time = time_gp,
                  units_m = units_m,
                  mda_times = mda_times, survey_times_data = survey_times_data,
                  sst = sst,
                  int_mat = int_mat,
                  kappa = kappa,
                  penalty = penalty,
                  ID_coords = ID_coords, ID_re = ID_re, s_unique = s_unique, re_unique = re_unique,
                  fix_tau2 = fix_tau2, fix_alpha = fix_alpha,
                  return_samples = return_samples,
                  par0 = start_pars, cov_offset = cov_offset,
                  power_val = power_val,
                  start_beta = start_pars$beta,
                  start_alpha = start_pars$alpha,
                  start_gamma = start_pars$gamma,
                  start_cov_pars = c(start_pars$sigma2,
                                     start_pars$phi,
                                     start_pars$tau2,
                                     start_pars$sigma2_re),
                  start_psi = start_pars$psi,
                  control_mcmc = control_mcmc,
                  messages = messages)



  res$y <- y
  res$D <- original_D
  res$coords <- coords
  res$sst <- sst
  if(sst) res$time <- time_gp
  res$mda_times <- mda_times
  res$survey_times_data <- survey_times_data
  res$int_mat <- int_mat
  res$ID_coords <- ID_coords
  if(n_re>0) {
    res$re <- re_unique_f
    res$ID_re <- as.data.frame(ID_re)
    colnames(res$ID_re) <- names_re
  }
  res$power_val <- power_val
  res$fix_tau2 <- fix_tau2
  res$fix_alpha <- fix_alpha
  res$formula <- formula
  if(no_penalty) {
    res$penalty <- NULL
  } else {
    res$penalty <- penalty
  }
  res$crs <- st_crs(data)$input
  res$input_crs <- input_crs
  res$distance_units <- distance_units
  res$data_sf <- data
  res$family <- "binomial"
  res$kappa <- kappa
  res$units_m <- units_m
  res$cov_offset <- cov_offset
  res$call <- match.call()


  invlink <- NULL
  if (is.null(invlink)) {
    inv <- function(eta) plogis(eta)
    d1  <- function(eta) { p <- inv(eta); p * (1 - p) }
    d2  <- function(eta) { p <- inv(eta); d <- p * (1 - 2 * p) }
    invlink <- list(inv = inv, d1 = d1, d2 = d2, name = "canonical")
  }
  res$linkf <- invlink

  return(res)
}


##' Monte Carlo maximum likelihood fitting engine for dast()
##'
##' Analogous to glgpm_nong() for the standard GLGPM models: performs the
##' Laplace approximation and MCMC sampling via the shared
##' maxim_integrand()/laplace_sampling_mcmc() backend (dast's modified
##' binomial probability, mda_effect * plogis(eta), is expressed as a
##' custom invlink), then the Monte Carlo maximum likelihood optimisation
##' for a DAST model.
##' @noRd
dast_fit <-
  function(y, D, coords, time, units_m, kappa, penalty,
           mda_times, survey_times_data,
           sst,
           int_mat,
           par0, cov_offset, power_val,
           ID_coords, ID_re, s_unique, re_unique,
           fix_tau2, fix_alpha, return_samples,
           start_beta,
           start_alpha, start_gamma,
           start_cov_pars,
           start_psi,
           control_mcmc,
           messages = TRUE) {

    beta0 <- par0$beta
    mu0 <- D%*%beta0+cov_offset

    sigma2_0 <- par0$sigma2

    phi0 <- par0$phi

    tau2_0 <- par0$tau2

    if(is.null(tau2_0)) tau2_0 <- fix_tau2

    sigma2_re_0 <- par0$sigma2_re

    n_loc <- nrow(coords)
    n_re <- length(sigma2_re_0)
    n_samples <- (control_mcmc$n_sim-control_mcmc$burnin)/control_mcmc$thin

    u = dist(coords)
    if(sst) v = dist(time)

    if(sst) {
      psi0 <- par0$psi
      Sigma0 <- sigma2_0*matern_correlation(u = u, phi = phi0, kappa = kappa,
                                    return_sym_matrix = TRUE) *
        matern_correlation(u = v, phi = psi0, kappa = 0.5,
                   return_sym_matrix = TRUE)
    } else {
      Sigma0 <- sigma2_0*matern_correlation(u = u, phi = phi0, kappa = kappa,
                                    return_sym_matrix = TRUE)
    }

    diag(Sigma0) <- diag(Sigma0) + tau2_0

    sigma2_re_0 <- par0$sigma2_re

    alpha0 <- par0$alpha
    if(is.null(alpha0)) alpha0<- fix_alpha

    gamma0 <- par0$gamma

    mda_effect0 <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                                      alpha0, gamma0, kappa = power_val)

    # dast's binomial probability is mda_effect * plogis(eta) - a fixed
    # (at this trial parameter vector) multiple of the canonical link, so it
    # can be expressed as a custom invlink for the shared glgpm backend
    # (maxim_integrand()/laplace_sampling_mcmc()) rather than dast's own
    # hand-rolled equivalents.
    # clamped like dast's own cond.dens.W()/lang.grad(): mda_effect0 can be
    # exactly 0 once alpha saturates near 1, which would otherwise make the
    # generic backend reject every state as having an invalid probability.
    mda_eps <- .Machine$double.eps
    mda_invlink0 <- list(
      inv = function(eta) pmin(pmax(mda_effect0 * plogis(eta), mda_eps), 1 - mda_eps),
      d1  = function(eta) { ps <- plogis(eta); mda_effect0 * ps * (1 - ps) },
      d2  = function(eta) { ps <- plogis(eta); d <- ps * (1 - ps); mda_effect0 * d * (1 - 2 * ps) }
    )

    if(messages) message("\n - Obtaining covariance matrix and mean for the proposal distribution of the MCMC \n")
    out_maxim <-
      maxim_integrand(y = y, units_m = units_m, Sigma = Sigma0, mu = mu0,
                      ID_coords = ID_coords, ID_re = ID_re,
                      family = "binomial", sigma2_re = sigma2_re_0,
                      hessian = FALSE, gradient = TRUE, invlink = mda_invlink0)

    Sigma_pd <- out_maxim$Sigma.tilde
    mean_pd <- out_maxim$mode

    simulation <-
      laplace_sampling_mcmc(y = y, units_m = units_m, mu = mu0,
                            Sigma = Sigma0, ID_coords = ID_coords, ID_re = ID_re,
                            sigma2_re = sigma2_re_0, family = "binomial",
                            control_mcmc = control_mcmc, invlink = mda_invlink0,
                            Sigma_pd = Sigma_pd, mean_pd = mean_pd,
                            messages = messages)

    S_tot_samples <- simulation$samples$S

    # Get the number of columns in D
    p <- ncol(D)

    # Define index for beta, sigma2, and phi
    ind_beta <- 1:p
    ind_sigma2 <- p + 1
    ind_phi <- p + 2

    # Conditional indexing based on fix_tau2 and n_re
    if (!is.null(fix_tau2)) {
      if (n_re > 0) {
        ind_sigma2_re <- (p + 2 + 1):(p + 2 + n_re)
        n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[, i])))
      } else {
        n_dim_re <- NULL
        ind_re <- NULL
      }

      if (is.null(fix_alpha)) {
        ind_alpha <- p + n_re + 3
        ind_gamma <- p + n_re + 4
        ind_psi <- p + n_re + 5  # psi when alpha is not fixed
      } else {
        ind_gamma <- p + n_re + 3
        ind_psi <- p + n_re + 4  # psi when alpha is fixed
      }
    } else {
      ind_nu2 <- p + 3
      if (n_re > 0) {
        ind_sigma2_re <- (p + 3 + 1):(p + 3 + n_re)
        n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[, i])))
      } else {
        n_dim_re <- NULL
        ind_re <- NULL
      }

      if (is.null(fix_alpha)) {
        ind_alpha <- p + n_re + 4
        ind_gamma <- p + n_re + 5
        ind_psi <- p + n_re + 6  # psi when alpha is not fixed
      } else {
        ind_gamma <- p + n_re + 4
        ind_psi <- p + n_re + 5  # psi when alpha is fixed
      }
    }


    if(n_re> 0) {
      for(i in 1:n_re) {
        S_tot_samples <- cbind(S_tot_samples, simulation$samples[[i+1]])
      }
      ind_re <- list()
      add_i <- 0
      for(i in 1:n_re) {
        ind_re[[i]] <- (add_i+n_loc+1):(add_i+n_loc+n_dim_re[i])
        if(i < n_re) add_i <- sum(n_dim_re[1:i])
      }
    }


    log.integrand <- function(S_tot, val) {
      n <- length(y)
      S <- S_tot[1:n_loc]

      q.f_re <- 0
      if(n_re > 0) {
        S_re <- NULL
        S_re_list <- list()
        for(i in 1:n_re) {
          S_re_list[[i]] <- S_tot[ind_re[[i]]]
          q.f_re <- q.f_re + n_dim_re[i]*log(val$sigma2_re[i])+
            sum(S_re_list[[i]]^2)/val$sigma2_re[i]
        }
      } else {
        q.f_re <- 0
      }

      eta <- val$mu + S[ID_coords]
      if(n_re > 0) {
        for(i in 1:n_re) {
          eta <- eta + S_re_list[[i]][ID_re[,i]]
        }
      }

      prob_star <- 1/(1+exp(-eta))
      prob <- val$mda_effect*prob_star

      llik <- sum(y*log(prob)+(units_m-y)*log(1-prob))

      q.f_S <- n_loc*log(val$sigma2)+val$ldetR+t(S)%*%val$R.inv%*%S/val$sigma2
      # Subtract both alpha and gamma penalties
      out <- -0.5*(q.f_S+q.f_re)+llik
      return(out)
    }

    compute.log.f <- function(par,ldetR=NA,R.inv=NA) {
      beta <- par[ind_beta]
      sigma2 <- exp(par[ind_sigma2])

      if(is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha])/(1+exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }

      gamma <- exp(par[ind_gamma])

      if(length(fix_tau2)>0) {
        nu2 <- fix_tau2/sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }
      phi <- exp(par[ind_phi])

      if(sst) psi <- exp(par[ind_psi])

      val <- list()
      val$sigma2 <- sigma2
      val$mu <- as.numeric(D%*%beta)+cov_offset
      val$mda_effect <- compute_mda_effect(survey_times_data, mda_times,
                                           intervention = int_mat,
                                           alpha, gamma, kappa = power_val)


      if(n_re > 0) {
        val$sigma2_re <- exp(par[ind_sigma2_re])
      }
      if(is.na(ldetR) & is.na(as.numeric(R.inv)[1])) {
        if(sst) {
          R <- matern_correlation(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)*
            matern_correlation(v, phi = psi, kappa=kappa,return_sym_matrix = TRUE)
        } else {
          R <- matern_correlation(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
        }
        diag(R) <- diag(R)+nu2
        val$ldetR <- determinant(R)$modulus
        val$R.inv <- solve(R)
      } else {
        val$ldetR <- ldetR
        val$R.inv <- R.inv
      }
      sapply(1:n_samples,
             function(i) log.integrand(S_tot_samples[i,],val))
    }


    if(!is.null(fix_tau2)) {
      par0_vec <- c(par0$beta, log(c(par0$sigma2, par0$phi)))
    } else {
      par0_vec <- c(par0$beta, log(c(par0$sigma2, par0$phi, par0$tau2/par0$sigma2)))
    }

    if(n_re > 0) {
      if(is.null(fix_alpha)) {
        par0_vec <- c(par0_vec,log(par0$sigma2_re),
                      log(par0$alpha/(1-par0$alpha)), log(par0$gamma))
      } else {
        par0_vec <- c(par0_vec, log(par0$sigma2_re),
                      log(par0$gamma))
      }
    } else {
      if(is.null(fix_alpha)) {
        par0_vec <- c(par0_vec,
                      log(par0$alpha/(1-par0$alpha)), log(par0$gamma))
      } else {
        par0_vec <- c(par0_vec,
                      log(par0$gamma))
      }
    }
    if(sst) par0_vec <- c(par0_vec, log(par0$psi))

    log.f.tilde <- compute.log.f(par0_vec)

    MC.log.lik <- function(par) {
      alpha <- if (is.null(fix_alpha)) exp(par[ind_alpha])/(1+exp(par[ind_alpha])) else fix_alpha
      gamma <- exp(par[ind_gamma])
      log(mean(exp(compute.log.f(par)-log.f.tilde))) -
        penalty[[1]](alpha) - penalty[[4]](gamma)
    }

    grad.MC.log.lik <- function(par) {
      beta <- par[ind_beta]
      mu <- as.numeric(D %*% beta) + cov_offset
      sigma2 <- exp(par[ind_sigma2])
      if (is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha]) / (1 + exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }
      gamma <- exp(par[ind_gamma])
      mda_effect_all <- compute_mda_effect_derivatives(survey_times_data, mda_times,
                                                       intervention = int_mat,
                                                       alpha, gamma, kappa = power_val)
      mda_effect <- mda_effect_all$effect
      mda_der_alpha <- mda_effect_all$d_alpha
      mda_der_gamma <- mda_effect_all$d_gamma

      if (length(fix_tau2) > 0) {
        nu2 <- fix_tau2 / sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }
      phi <- exp(par[ind_phi])

      if(sst) psi <- exp(par[ind_psi])
      if (n_re > 0) {
        sigma2_re <- exp(par[ind_sigma2_re])
      }

      if (sst) {
        R_u <- matern_correlation(u, phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R_v <- matern_correlation(v, phi = psi, kappa = 0.5, return_sym_matrix = TRUE)
        R <- R_u * R_v
      } else {
        R <- matern_correlation(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
      }
      diag(R) <- diag(R) + nu2
      R.inv <- solve(R)
      ldetR <- determinant(R)$modulus

      if(sst) {
        psi <- exp(par[ind_psi])
        R1.psi <- R_u * matern_gradient_phi(v, psi, 0.5)
        m1.psi <- R.inv %*% R1.psi
        t1.psi <- -0.5 * sum(diag(m1.psi))
        m2.psi <- m1.psi %*% R.inv
      }

      exp.fact <- exp(compute.log.f(par, ldetR, R.inv) - log.f.tilde)
      L.m <- sum(exp.fact)
      exp.fact <- exp.fact / L.m

      if (sst) {
        R1.phi <- matern_gradient_phi(u, phi, kappa) * R_v
      } else {
        R1.phi <- matern_gradient_phi(u, phi, kappa)
      }
      m1.phi <- R.inv %*% R1.phi
      t1.phi <- -0.5 * sum(diag(m1.phi))
      m2.phi <- m1.phi %*% R.inv

      if (is.null(fix_tau2)) {
        t1.nu2 <- -0.5 * sum(diag(R.inv))
        m2.nu2 <- R.inv %*% R.inv
      }

      gradient.S <- function(S_tot) {
        S <- S_tot[1:n_loc]
        if (n_re > 0) {
          S_re_list <- list()
          for (i in 1:n_re) {
            S_re_list[[i]] <- S_tot[ind_re[[i]]]
          }
        }

        eta <- mu + S[ID_coords]
        if (n_re > 0) {
          for (i in 1:n_re) {
            eta <- eta + S_re_list[[i]][ID_re[, i]]
          }
        }

        prob_star <- 1 / (1 + exp(-eta))
        prob <- mda_effect * prob_star
        d_S <- (y / prob - (units_m - y) / (1 - prob)) * mda_effect * prob_star / (1 + exp(eta))
        d_S <- as.numeric(d_S)

        q.f_S <- t(S) %*% R.inv %*% S
        grad.beta <- t(D) %*% d_S
        grad.log.sigma2 <- (-n_loc / (2 * sigma2) + 0.5 * q.f_S / (sigma2^2)) * sigma2
        grad.log.phi <- (t1.phi + 0.5 * as.numeric(t(S) %*% m2.phi %*% S) / sigma2) * phi

        if (is.null(fix_alpha)) {
          der.alpha <- exp(par[ind_alpha]) / ((1 + exp(par[ind_alpha]))^2)
          grad.alpha.t <- der.alpha * (
            sum((y / mda_effect - (units_m - y) * prob_star / (1 - prob)) * mda_der_alpha) -
              penalty[[2]](alpha)
          )
        }

        # Gradient wrt log(gamma): chain rule gives factor gamma.
        # Also subtract the derivative of the gamma penalty wrt log(gamma),
        # which is gamma * penalty[[5]](gamma) by the chain rule.
        grad.log.gamma <- gamma * (
          sum((y / mda_effect - (units_m - y) * prob_star / (1 - prob)) * mda_der_gamma) -
            penalty[[5]](gamma)
        )

        out <- c(grad.beta, grad.log.sigma2, grad.log.phi)

        if (is.null(fix_tau2)) {
          grad.log.nu2 <- (t1.nu2 + 0.5 * as.numeric(t(S) %*% m2.nu2 %*% S) / sigma2) * nu2
          out <- c(out, grad.log.nu2)
        }

        if (n_re > 0) {
          grad.log.sigma2_re <- rep(NA, n_re)
          for (i in 1:n_re) {
            grad.log.sigma2_re[i] <- (-n_dim_re[i] / (2 * sigma2_re[i]) +
                                        0.5 * sum(S_re_list[[i]]^2) / (sigma2_re[i]^2)) * sigma2_re[i]
          }
          out <- c(out, grad.log.sigma2_re)
        }

        if (is.null(fix_alpha)) {
          out <- c(out, grad.alpha.t, grad.log.gamma)
        } else {
          out <- c(out, grad.log.gamma)
        }

        if (sst) {
          grad.log.psi <- (t1.psi + 0.5 * as.numeric(t(S) %*% m2.psi %*% S) / sigma2) * psi
          out <- c(out, grad.log.psi)
        }

        return(out)
      }

      out <- rep(0, length(par))
      for (i in 1:n_samples) {
        out <- out + exp.fact[i] * gradient.S(S_tot_samples[i, ])
      }
      out
    }

    hess.MC.log.lik <- function(par) {
      # Unpack parameters
      beta   <- par[ind_beta]
      mu     <- as.numeric(D %*% beta) + cov_offset
      sigma2 <- exp(par[ind_sigma2])
      if (is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha]) / (1 + exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }
      gamma <- exp(par[ind_gamma])

      # MDA effect + derivatives
      mda_all           <- compute_mda_effect_derivatives(
        survey_times_data, mda_times,
        intervention = int_mat,
        alpha, gamma, kappa = power_val
      )
      mda_eff        <- mda_all$effect
      mda_der_alpha  <- mda_all$d_alpha
      mda_der_gamma  <- mda_all$d_gamma
      mda_der2_alpha <- mda_all$d2_alpha
      mda_der2_gamma <- mda_all$d2_gamma
      mda_der2_alpha_gamma <- mda_all$d2_alpha_gamma

      # nu2
      if (!is.null(fix_tau2)) {
        nu2 <- fix_tau2 / sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }

      # Build R and its derivatives
      phi <- exp(par[ind_phi])
      if (sst) {
        psi    <- exp(par[ind_psi])
        R_u    <- matern_correlation(u,   phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R_v    <- matern_correlation(v,   phi = psi, kappa = 0.5,    return_sym_matrix = TRUE)
        R      <- R_u * R_v
        R1.phi <- matern_gradient_phi(u, phi, kappa) * R_v
        R2.phi <- matern_hessian_phi(u, phi, kappa) * R_v
        R1.psi <- R_u * matern_gradient_phi(v, psi, 0.5)
        R2.psi <- R_u * matern_hessian_phi(v, psi, 0.5)
      } else {
        R      <- matern_correlation(u, phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R1.phi <- matern_gradient_phi(u, phi, kappa)
        R2.phi <- matern_hessian_phi(u, phi, kappa)
      }
      diag(R) <- diag(R) + nu2
      R.inv   <- solve(R)

      # Precompute phi derivatives
      m1.phi <- R.inv %*% R1.phi
      t1.phi <- -0.5 * sum(diag(m1.phi))
      m2.phi <- m1.phi %*% R.inv
      t2.phi <- -0.5 * sum(diag(R.inv %*% R2.phi - R.inv %*% R1.phi %*% R.inv %*% R1.phi))
      n2.phi <- R.inv %*% (2 * R1.phi %*% R.inv %*% R1.phi - R2.phi) %*% R.inv

      if (sst) {
        R2.psi.phi <- matern_gradient_phi(u, phi, 0.5)*matern_gradient_phi(v, psi, 0.5)
        m1.psi   <- R.inv %*% R1.psi
        t1.psi   <- -0.5 * sum(diag(m1.psi))
        m2.psi   <- m1.psi %*% R.inv
        t2.psi   <- -0.5 * sum(diag(R.inv %*% R2.psi - R.inv %*% R1.psi %*% R.inv %*% R1.psi))
        n2.psi   <- R.inv %*% (2 * R1.psi %*% R.inv %*% R1.psi - R2.psi) %*% R.inv

        t2.psi.phi <- -0.5*(sum(R.inv*R2.psi.phi)-
                              sum(m1.phi*t(m1.psi)))
        n2.psi.phi <- R.inv%*%(R1.phi%*%m1.psi+
                                 R1.psi%*%m1.phi-
                                 R2.psi.phi)%*%R.inv

      }

      # Precompute nu2 derivatives
      if (is.null(fix_tau2)) {
        m2.nu2     <- R.inv %*% R.inv
        t1.nu2     <- -0.5 * sum(diag(R.inv))
        t2.nu2     <-  0.5 * sum(diag(m2.nu2))
        n2.nu2     <-  2 * R.inv %*% m2.nu2
        t2.nu2.phi <-  0.5 * sum(diag(R.inv %*% R1.phi %*% R.inv))
        n2.nu2.phi <-  R.inv %*% (R.inv %*% R1.phi + R1.phi %*% R.inv) %*% R.inv
        if (sst) {
          t2.nu2.psi <-  0.5 * sum(diag(R.inv %*% R1.psi %*% R.inv))
          n2.nu2.psi <-  R.inv %*% (R.inv %*% R1.psi + R1.psi %*% R.inv) %*% R.inv
        }
      }

      # Monte Carlo weights
      ldetR <- determinant(R)$modulus
      expf  <- exp(compute.log.f(par, ldetR, R.inv) - log.f.tilde)
      w     <- expf / sum(expf)

      # Accumulators
      H_acc <- matrix(0, nrow = length(par), ncol = length(par))
      g_acc <- rep(0, length(par))

      hessian.S <- function(S_tot, ef) {
        S <- S_tot[1:n_loc]
        if (n_re > 0) {
          S_re_list <- lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]])
        }

        # Linear predictor
        eta <- mu + S[ID_coords]
        if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]
        prob_star <- 1 / (1 + exp(-eta))
        prob      <- mda_eff * prob_star

        # Derivatives wrt S
        d_S  <- (y/prob - (units_m-y)/(1-prob)) * mda_eff * prob_star / (1+exp(eta))
        d2_S <- (-y/prob^2 - (units_m-y)/(1-prob)^2) * (mda_eff*prob_star/(1+exp(eta)))^2 +
          (y/prob - (units_m-y)/(1-prob)) * mda_eff * exp(eta)*(1-exp(eta))/(1+exp(eta))^3
        d_S  <- as.numeric(d_S)

        # d_S_alpha
        if (is.null(fix_alpha)) {
          d_S_alpha <- (-y/prob^2 - (units_m-y)/(1-prob)^2) *
            mda_der_alpha * prob_star * mda_eff * prob_star/(1+exp(eta)) +
            (y/prob - (units_m-y)/(1-prob)) * mda_der_alpha * prob_star/(1+exp(eta))
        }

        # d_S_gamma
        d_S_gamma <- (-y/prob^2 - (units_m-y)/(1-prob)^2) *
          mda_der_gamma * prob_star * mda_eff * prob_star/(1+exp(eta)) +
          (y/prob - (units_m-y)/(1-prob)) * mda_der_gamma * prob_star/(1+exp(eta))

        qSS <- as.numeric(t(S) %*% R.inv %*% S)

        # Gradient components
        grad_beta  <- t(D) %*% d_S
        grad_ls2   <- (-n_loc/(2*sigma2) + 0.5*qSS/sigma2^2) * sigma2
        grad_lphi  <- (t1.phi + 0.5*as.numeric(t(S)%*%m2.phi%*%S)/sigma2) * phi

        g_vec <- c(grad_beta, grad_ls2, grad_lphi)

        if (is.null(fix_tau2)) {
          grad_lnu2 <- (t1.nu2 + 0.5*as.numeric(t(S)%*%m2.nu2%*%S)/sigma2)*nu2
          g_vec     <- c(g_vec, grad_lnu2)
        }

        if (is.null(fix_alpha)) {
          der_alpha  <- exp(par[ind_alpha])/(1+exp(par[ind_alpha]))^2
          grad_alpha <- der_alpha * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_alpha) -
              penalty[[2]](alpha)
          )
          g_vec <- c(g_vec, grad_alpha)
        }

        # Gradient wrt log(gamma), including gamma penalty derivative.
        # d/d(log gamma) = gamma * d/d(gamma), so subtract gamma * penalty[[5]](gamma).
        grad_lgamma <- gamma * (
          sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_gamma) -
            penalty[[5]](gamma)
        )
        g_vec <- c(g_vec, grad_lgamma)

        if (n_re > 0) {
          grad_ls2re <- sapply(seq_len(n_re), function(i, sigma2_re)
            (-n_dim_re[i]/(2*sigma2_re[i]) + 0.5*sum(S_re_list[[i]]^2)/sigma2_re[i]^2)*sigma2_re[i]
          )
          g_vec <- c(g_vec, grad_ls2re)
        }

        if (sst) {
          grad_lpsi <- (t1.psi + 0.5*as.numeric(t(S)%*%m2.psi%*%S)/sigma2)*psi
          g_vec     <- c(g_vec, grad_lpsi)
        }

        # Build H_loc
        H_loc <- matrix(0, nrow=length(par), ncol=length(par))

        # beta-beta
        H_loc[ind_beta, ind_beta] <- t(D) %*% (D * d2_S)

        # beta-alpha
        if (is.null(fix_alpha)) {
          H_loc[ind_beta, ind_alpha] <- t(D) %*% d_S_alpha * der_alpha
          H_loc[ind_alpha, ind_beta] <- H_loc[ind_beta, ind_alpha]
        }

        # beta-gamma
        H_loc[ind_beta, ind_gamma] <- t(D) %*% d_S_gamma * gamma
        H_loc[ind_gamma, ind_beta] <- H_loc[ind_beta, ind_gamma]

        # alpha-alpha
        if (is.null(fix_alpha)) {
          der2_alpha <- exp(par[ind_alpha])*(1-exp(par[ind_alpha]))/(1+exp(par[ind_alpha]))^3
          H_loc[ind_alpha, ind_alpha] <- der2_alpha * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_alpha) -
              penalty[[2]](alpha)
          ) + der_alpha^2 * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_alpha +
                  (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_alpha^2
            ) - penalty[[3]](alpha)
          )
        }

        # alpha-gamma
        if (is.null(fix_alpha)) {
          H_loc[ind_alpha, ind_gamma] <- (gamma*der_alpha)*sum(
            (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_alpha_gamma +
              (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_alpha*mda_der_gamma
          )
          H_loc[ind_gamma, ind_alpha] <- H_loc[ind_alpha, ind_gamma]
        }

        # gamma-gamma.
        # The Hessian of the penalized log-likelihood wrt log(gamma) is:
        #   H_lik[gamma,gamma]
        #   - gamma^2 * penalty[[6]](gamma)   [from d/d(log g) of gamma * pen'(g)]
        #   - gamma   * penalty[[5]](gamma)   [the grad_lgamma term itself]
        H_loc[ind_gamma, ind_gamma] <- gamma*sum(
          (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_gamma
        ) + gamma^2*sum(
          (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_gamma +
            (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_gamma^2
        ) - gamma^2 * penalty[[6]](gamma) - gamma * penalty[[5]](gamma)

        # sigma2-sigma2
        H_loc[ind_sigma2, ind_sigma2] <-
          (n_loc/(2*sigma2^2) - qSS/sigma2^3)*sigma2^2 + grad_ls2

        # sigma2-phi
        H_loc[ind_sigma2, ind_phi] <- (grad_lphi/phi - t1.phi)*(-phi)
        H_loc[ind_phi, ind_sigma2] <- H_loc[ind_sigma2, ind_phi]

        # phi-phi
        H_loc[ind_phi, ind_phi] <- (t2.phi - 0.5*as.numeric(t(S)%*%n2.phi%*%S)/sigma2)*phi^2 + grad_lphi

        # nu2-block
        if (is.null(fix_tau2)) {
          H_loc[ind_sigma2, ind_nu2] <- (grad_lnu2/nu2 - t1.nu2)*(-nu2)
          H_loc[ind_nu2, ind_sigma2] <- H_loc[ind_sigma2, ind_nu2]
          H_loc[ind_nu2, ind_nu2]   <- (t2.nu2 - 0.5*as.numeric(t(S)%*%n2.nu2%*%S)/sigma2)*nu2^2 + grad_lnu2
          H_loc[ind_phi, ind_nu2]   <- (t2.nu2.phi - 0.5*as.numeric(t(S)%*%n2.nu2.phi%*%S)/sigma2)*phi*nu2
          H_loc[ind_nu2, ind_phi]   <- H_loc[ind_phi, ind_nu2]
        }

        # psi-block
        if (sst) {
          H_loc[ind_psi, ind_psi]     <- (t2.psi - 0.5*as.numeric(t(S)%*%n2.psi%*%S)/sigma2)*psi^2 + grad_lpsi
          H_loc[ind_sigma2, ind_psi]    <- (grad_lpsi/psi - t1.psi)*(-psi)
          H_loc[ind_psi, ind_sigma2]    <- H_loc[ind_sigma2, ind_psi]
          H_loc[ind_phi, ind_psi]       <- (t2.psi.phi-0.5*t(S)%*%n2.psi.phi%*%S/sigma2)*phi*psi
          H_loc[ind_psi, ind_phi]       <- H_loc[ind_phi, ind_psi]
          if (is.null(fix_tau2)) {
            val <- (t2.nu2.psi - 0.5*as.numeric(t(S)%*%n2.nu2.psi%*%S)/sigma2)*nu2*psi
            H_loc[ind_nu2, ind_psi] <- val
            H_loc[ind_psi, ind_nu2] <- val
          }
        }

        list(
          mat1 = ef * (g_vec %*% t(g_vec) + H_loc),
          g    = ef * g_vec
        )
      }

      # Monte Carlo aggregation
      for (i in seq_len(n_samples)) {
        tmp   <- hessian.S(S_tot_samples[i, ], w[i])
        g_acc <- g_acc + tmp$g
        H_acc <- H_acc + tmp$mat1
      }

      H_acc - g_acc %*% t(g_acc)
    }

    start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)]/start_cov_pars[1]
    if(is.null(fix_alpha)) {
      start_par <- c(start_beta, log(start_cov_pars),
                     log(start_alpha/(1-start_alpha)), log(start_gamma))
    } else {
      start_par <- c(start_beta, log(start_cov_pars), log(start_gamma))
    }
    if(sst) start_par <- c(start_par, log(start_psi))


    out <- list()
    estim <- nlminb(
      start     = start_par,
      objective = function(x)   -MC.log.lik(x),
      gradient  = function(x)   -grad.MC.log.lik(x),
      hessian   = function(x)   -hess.MC.log.lik(x),
      control   = list(trace=1*messages)
    )

    out$estimate <- list(beta   = setNames(estim$par[ind_beta], colnames(D)),
                         sigma2 = unname(estim$par[ind_sigma2]),
                         phi    = unname(estim$par[ind_phi]))
    if (is.null(fix_tau2)) out$estimate$nu2 <- unname(estim$par[ind_nu2])
    if (n_re > 0) {
      sigma2_re <- unname(estim$par[ind_sigma2_re])
      names(sigma2_re) <- colnames(ID_re)
      out$estimate$sigma2_re <- sigma2_re
    }
    if (is.null(fix_alpha)) out$estimate$alpha <- unname(estim$par[ind_alpha])
    out$estimate$gamma <- unname(estim$par[ind_gamma])
    if (sst) out$estimate$psi <- unname(estim$par[ind_psi])

    out$grad_MLE <- grad.MC.log.lik(estim$par)
    hess.MLE <- hess.MC.log.lik(estim$par)
    out$covariance <- solve(-hess.MLE)
    out$log_lik <- -estim$objective
    if(return_samples) out$S_samples <- S_tot_samples
    class(out) <- c("RiskMap_dast", "RiskMap")
    return(out)
  }


##' @title Extract Parameter Estimates from a "RiskMap_dast" Model Fit
##' @description \code{coef} method for \code{dast()} fits; extends
##' \code{coef.RiskMap()} with the MDA decay parameters.
##' @param object A fitted \code{dast()} object.
##' @param ... other parameters, passed on to \code{coef.RiskMap()}.
##' @return A list as for \code{coef.RiskMap()}, plus \code{alpha} (unless
##' fixed via \code{drop}), \code{gamma}, and \code{psi} (if spatio-temporal).
##' @seealso \code{\link{dast}}, \code{\link{coef.RiskMap}}
##' @method coef RiskMap_dast
##' @export
coef.RiskMap_dast <- function(object, ...) {
  res <- NextMethod()
  estimate <- object$estimate

  if (!is.null(estimate$alpha)) res$alpha <- plogis(estimate$alpha)
  res$gamma <- exp(estimate$gamma)
  if (!is.null(estimate$psi)) res$psi <- exp(estimate$psi)

  res
}

##' @title Summarize dast() Model Fits
##' @description \code{summary} method for \code{dast()} fits; extends
##' \code{summary.RiskMap()} with the MDA decay parameters.
##' @param object A fitted \code{dast()} object.
##' @param ... other parameters, passed on to \code{summary.RiskMap()}.
##' @param conf_level Confidence level for intervals (default 0.95).
##' @return A list of class \code{c("summary.RiskMap_dast", "summary.RiskMap")}
##' with an additional \code{mda} matrix of MDA decay parameter estimates/CIs.
##' @seealso \code{\link{dast}}, \code{\link{summary.RiskMap}}
##' @method summary RiskMap_dast
##' @export
summary.RiskMap_dast <- function(object, ..., conf_level = 0.95) {
  res <- NextMethod()

  z_crit   <- qnorm(1 - (1 - conf_level) / 2)
  estimate <- unlist(object$estimate)
  se_par   <- setNames(sqrt(diag(object$covariance)), names(estimate))

  mda_rows <- list()
  if ("alpha" %in% names(estimate)) {
    est <- estimate[["alpha"]]
    se <- se_par[["alpha"]]
    mda_rows[["Drop (alpha)"]] <- c("Estimate" = plogis(est),
                                    "Lower limit" = plogis(est - z_crit * se),
                                    "Upper limit" = plogis(est + z_crit * se))
  }
  est <- estimate[["gamma"]]
  se <- se_par[["gamma"]]
  mda_rows[["Scale of the decay (gamma)"]] <- c("Estimate" = exp(est),
                                                "Lower limit" = exp(est - z_crit * se),
                                                "Upper limit" = exp(est + z_crit * se))
  res$mda       <- do.call(rbind, mda_rows)
  res$fix_alpha <- object$fix_alpha
  res$power_val <- object$power_val

  if ("psi" %in% names(estimate)) {
    est <- estimate[["psi"]]
    se <- se_par[["psi"]]
    psi_row <- c("Estimate" = exp(est),
                 "Lower limit" = exp(est - z_crit * se),
                 "Upper limit" = exp(est + z_crit * se))
    res$sp <- rbind(res$sp, "Temporal corr. scale (psi)" = psi_row)
  }

  class(res) <- c("summary.RiskMap_dast", "summary.RiskMap")
  res
}

##' @title Print Summary of dast() Model Fits
##' @description Print method for \code{"summary.RiskMap_dast"} objects;
##' extends \code{print.summary.RiskMap()} with the MDA decay parameters.
##' @param x An object of class \code{"summary.RiskMap_dast"}.
##' @param ... other parameters.
##' @return Invisibly returns \code{x}.
##' @method print summary.RiskMap_dast
##' @export
print.summary.RiskMap_dast <- function(x, ...) {
  NextMethod()

  cat("\nMDA impact function\n")
  cat("f(v) = alpha * exp(-(v/gamma)^delta),  delta fixed at ", x$power_val, "\n", sep = "")
  printCoefmat(x$mda, P.values = FALSE, has.Pvalue = FALSE)
  if (!is.null(x$fix_alpha))
    cat("Drop (alpha) fixed at ", x$fix_alpha, "\n", sep = "")

  invisible(x)
}
