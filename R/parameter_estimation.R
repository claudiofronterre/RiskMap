


##' Simulation from Generalized Linear Gaussian Process Models
##'
##' Simulates data from a fitted Generalized Linear Gaussian Process Model (GLGPM) or a specified model formula and data.
##'
##' @param n_sim Number of simulations to perform.
##' @param model_fit Fitted GLGPM model object of class `RiskMap`. If provided, overrides `formula`, `data`, `family`, `convert_to_crs` and `scale_to_km` arguments.
##' @param formula Model formula indicating the variables of the model to be simulated.
##' @param data `sf` object containing the variables in the model formula.
##' @param family Distribution family for the response variable. Must be one of `gaussian`, `binomial`, or `poisson.`
##' @param den Required for `binomial` to denote the denominator (i.e. number of trials) of the Binomial distribution.
##' For the `poisson` family, the argument is optional and is used a multiplicative term to express the mean counts.
##' @param cov_offset Offset for the covariate part of the GLGPM.
##' @param convert_to_crs CRS code to convert data to.
##' @param scale_to_km Logical; if `TRUE`, distances between locations are computed in kilometers; if `FALSE`, in meters.
##' @param sim_pars List of simulation parameters including `beta`, `sigma2`, `tau2`, `phi`, `sigma2_me`, and optionally `sigma2_re`.
##' If multiple covariates or random effects are included, the lengths of `beta` and `sigma2_re` must match the number of covariates and random effects respectively.
##' @param messages Logical; if `TRUE`, display progress and informative messages.
##'
##' @details
##' Generalized Linear Gaussian Process Models (GLGPMs) extend generalized linear models (GLMs) by incorporating spatial Gaussian processes to model spatial correlation. This function simulates data from GLGPMs using Markov Chain Monte Carlo (MCMC) methods. It supports Gaussian, binomial, and Poisson response families, utilizing a Matern correlation function to model spatial dependence.
##'
##' The simulation process involves generating spatially correlated random effects and simulating responses based on the fitted or specified model parameters. For `gaussian` family, the function simulates response values by adding measurement error.
##'
##' Additionally, GLGPMs can incorporate unstructured random effects specified through the [`re()`] term in the model formula, allowing for capturing additional variability beyond fixed and spatial covariate effects.
##'
##' @return A list containing simulated data, simulated spatial random effects (if applicable), and other simulation parameters.
##' @export
simulate_glgpm <- function(n_sim,
                      model_fit = NULL,
                      formula = NULL,
                      data = NULL,
                      family = NULL,
                      den = NULL,
                      cov_offset = NULL,
                      convert_to_crs = NULL,
                      scale_to_km = TRUE,
                      sim_pars = list(beta = NULL,
                                      sigma2 = NULL,
                                      tau2 = NULL,
                                      phi = NULL,
                                      sigma2_me = NULL,
                                      sigma2_re = NULL),
                      messages = TRUE) {

  check_positive_integer(n_sim, "n_sim")

  if(!is.null(model_fit)) {
    if(!inherits(model_fit, "RiskMap")){
      stop("'model_fit' must be of class 'RiskMap'")
    }
    if (!is.null(data) | !is.null(formula)){
      stop("if you provide 'model_fit' you should not provide 'data' or 'formula'")
    }
    formula <- as.formula(model_fit$formula)
    data <- model_fit$data_sf
    family <- model_fit$family
    convert_to_crs <- model_fit$convert_to_crs
    scale_to_km <- model_fit$scale_to_km
  }

  check_data(data)
  check_formula(formula, data)
  inter_f <- interpret.formula(formula)

  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")

  if(family != "gaussian" & family != "binomial" &
     family != "poisson") stop("'family' must be either 'gaussian', 'binomial'
                               or 'poisson'")

  mf <- model.frame(inter_f$pf,data = data, na.action = na.fail)

  # Extract covariates matrix
  D <- as.matrix(model.matrix(attr(mf,"terms"), data = data))
  n <- nrow(D)

  hr_re <- if (length(inter_f$re.spec) > 0L) {
    inter_f$re.spec$term
  } else {
    NULL
  }
  random_effects <- prepare_random_effects(data, hr_re)
  n_re <- random_effects$n_re
  ID_re <- random_effects$ID_re
  re_unique <- random_effects$re_unique

  # Number of covariates
  p <- ncol(D)

  if(!is.null(model_fit)) {
    par_hat <- coef(model_fit)

    beta <- par_hat$beta
    sigma2 <- par_hat$sigma2
    phi <- par_hat$phi

    if(isTRUE(model_fit$fix_tau2)) {
      tau2 <- par_hat$tau2
      if(is.null(model_fit$fix_var_me)) {
        sigma2_me <- par_hat$sigma2_me
      } else {
        sigma2_me <- model_fit$fix_var_me
      }
      if(n_re>0) {
        sigma2_re <- par_hat$sigma2_me
      }
    } else {
      tau2 <- model_fit$fix_tau2
      if(is.null(model_fit$fix_var_me)) {
        sigma2_me <- par_hat$sigma2_me
      } else {
        sigma2_me <- model_fit$fix_var_me
      }
      if(n_re>0) {
        sigma2_re <- par_hat$sigma2_re
      }
    }
  } else {
    # extract non-NULL names
    par_names <- names(sim_pars)[!vapply(sim_pars, is.null, logical(1))]

    # [[]] syntax avoids partial matching
    if (!"beta" %in% par_names) stop("'beta' is missing")
    beta <- sim_pars[["beta"]]
    if (length(beta)!=p) stop("the number of values provided for 'beta' must be one plus
    the number of covariates specified in the formula")
    if (!"sigma2" %in% par_names) stop("'sigma2' is missing")
    sigma2 <- sim_pars[["sigma2"]]
    if (!"phi" %in% par_names) stop("'phi' is missing")
    phi <- sim_pars[["phi"]]
    if (!"tau2" %in% par_names) stop("'tau2' is missing")
    tau2 <- sim_pars[["tau2"]]
    if (!"sigma2_me" %in% par_names) stop("'sigma2_me' is missing")
    sigma2_me <- sim_pars[["sigma2_me"]]
    if (n_re > 0) {
      if(!"sigma2_re" %in% par_names) stop("'sigma2_re' is missing")
      if(length(sim_pars[["sigma2_re"]]) != n_re) stop("the values passed to 'sigma2_re' in 'sim_pars'
      does not match the number of random effects specfied in re() in the formula")
      sigma2_re <- sim_pars[["sigma2_re"]]
    }
    if (n_re == 0 & "sigma2_re" %in% par_names){
      warning("'sigma2_re' will be ignored as no random effects are included")
    }
  }

  # Extract coordinates
  if(!is.null(convert_to_crs)) {
    if(!is.numeric(convert_to_crs)) stop("'convert_to_crs' must be a numeric object")
    data <- st_transform(data, crs = convert_to_crs)
    crs <- convert_to_crs
  }
  if(messages) message("The CRS used is ", as.list(st_crs(data))$input, "\n")

  coords_o <- st_coordinates(data)
  coords <- unique(coords_o)

  m <- nrow(coords_o)
  ID_coords <- sapply(1:m, function(i)
    which(coords_o[i,1]==coords[,1] &
            coords_o[i,2]==coords[,2]))
  s_unique <- unique(ID_coords)

  if(all(table(ID_coords)==1) & !is.null(tau2) &
     !is.null(sigma2_me) && (tau2!=0 & sigma2_me!=0)) {
    warning("When there is only one observation per location, both the nugget and measurement error cannot
         be estimated. Consider removing either one of them. ")
  }

  if(scale_to_km) {
    coords_o <- coords_o/1000
    coords <- coords/1000
    if(messages) message("Distances between locations are computed in kilometers \n")
  } else {
    if(messages) message("Distances between locations are computed in meters \n")
  }

  # Simulate S
  Sigma <- sigma2*matern_correlation(dist(coords), phi = phi, kappa = kappa,
                             return_sym_matrix = TRUE)
  diag(Sigma) <- diag(Sigma) + tau2
  Sigma_sroot <- t(chol(Sigma))
  S_sim <- t(sapply(1:n_sim, function(i) Sigma_sroot%*%rnorm(nrow(coords))))

  # Simulate random effects
  if(n_re>0) {
    re_sim <- list()
    if(!is.null(model_fit)) {
      re_names <- names(model_fit$re)
    } else {
      re_names <- inter_f$re.spec$term
    }

    dim_re <- sapply(1:n_re, function(j) length(re_unique[[j]]))
    for(i in 1:n_sim) {
      re_sim[[i]] <- list()
      for(j in 1:n_re) {
        re_sim[[i]][[paste(re_names[j])]] <- rnorm(dim_re[j])*sqrt(sigma2_re[j])
      }
    }
  }

  # Linear predictor
  # try adding cov_offset here
  eta_sim <- t(sapply(1:n_sim, function(i) D%*%beta + S_sim[i,][ID_coords]))

  if(n_re > 0) {
    for(i in 1:n_sim) {
      for(j in 1:n_re) {
        eta_sim[i,] <- eta_sim[i,] + re_sim[[i]][[paste(re_names[j])]][ID_re[,j]]
      }
    }
  }

  if(family!="gaussian") {
    if(!is.null(den))  {
      do_name <- deparse(substitute(den))
      y <- as.numeric(model.response(mf))
      units_m <- data[[do_name]]

      if (family == "binomial") check_binomial(y, units_m)
      if(is.integer(units_m)) units_m <- as.numeric(units_m)
      if(!is.numeric(units_m)) stop("the variable passed to `den` must be numeric")

    } else {
      units_m <- model_fit$units_m
    }

  }

  y_sim <- matrix(NA, nrow=n_sim, ncol=n)
  if(family=="gaussian") {

    for(i in 1:n_sim) {
      y_sim[i,] <- eta_sim[i,] + sqrt(sigma2_me)*rnorm(n)
    }
  } else {

    if(family=="binomial") {
      for(i in 1:n_sim) {
        prob_i <- exp(eta_sim[i,])/(1+exp(eta_sim[i,]))
        y_sim[i,] <- rbinom(n, size = units_m, prob = prob_i)
      }
    } else if(family=="poisson") {
      for(i in 1:n_sim) {
        mean_i <- exp(eta_sim[i,])/(1+exp(eta_sim[i,]))
        y_sim[i,] <- rpois(n,lambda = units_m*mean_i)
      }
    }
  }

  if(!is.null(model_fit)) {
    data_sim <- model_fit$data_sf
  } else {
    data_sim <- data
  }

  for(i in 1:n_sim) {
    data_sim[[paste(inter_f$response,"_sim",i,sep="")]] <- y_sim[i,]
  }
  out <- list(data_sim = data_sim,
              S_sim = S_sim,
              lin_pred_sim = eta_sim,
              beta = beta,
              sigma2 = sigma2,
              tau2 = tau2,
              phi = phi)
  if(family=="gaussian") {
    out$sigma2_me <- sigma2_me
  }
  if(n_re>0) {
    out$sigma2_re <- sigma2_re
    out$re_sim <- re_sim
  }
  return(out)
}

##' Maximization of the Integrand for Generalized Linear Gaussian Process Models
##'
##' Maximizes the integrand function for Generalized Linear Gaussian Process Models (GLGPMs), which involves the evaluation of likelihood functions with spatially correlated random effects.
##'
##' @param y Response variable vector.
##' @param units_m Units of measurement for the response variable.
##' @param mu Mean vector of the response variable.
##' @param Sigma Covariance matrix of the spatial process.
##' @param ID_coords Indices mapping response to locations.
##' @param ID_re Indices mapping response to unstructured random effects.
##' @param family Distribution family for the response variable. Must be one of 'gaussian', 'binomial', or 'poisson'.
##' @param sigma2_re Variance of the unstructured random effects.
##' @param hessian Logical; if TRUE, compute the Hessian matrix.
##' @param gradient Logical; if TRUE, compute the gradient vector.
##' @param invlink A function that defines the inverse of the link function for the distribution of the data given the random effects.
##'
##' @details
##' This function maximizes the integrand for GLGPMs using the Nelder-Mead optimization algorithm. It computes the likelihood function incorporating spatial covariance and unstructured random effects, if provided.
##'
##' The integrand includes terms for the spatial process (Sigma), unstructured random effects (sigma2_re), and the likelihood function (llik) based on the specified distribution family ('gaussian', 'binomial', or 'poisson').
##'
##' @return A list containing the mode estimate, and optionally, the Hessian matrix and gradient vector.
##' @export
maxim_integrand <- function(
    y, units_m, mu, Sigma, ID_coords, ID_re = NULL, family,
    sigma2_re = NULL, hessian = FALSE, gradient = FALSE, invlink = NULL
) {
  stopifnot(family %in% c("poisson", "binomial"))

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

  cross_sum <- function(v, grp1, n1, grp2, n2) {
    f1 <- factor(grp1, levels = seq_len(n1))
    f2 <- factor(grp2, levels = seq_len(n2))
    as.matrix(xtabs(v ~ f1 + f2))  # n1 x n2, zeros where empty
  }

  # ---------- dimensions ----------
  Sigma.inv <- solve(Sigma)
  n_loc <- nrow(Sigma)
  n <- length(y)

  if (( !is.null(ID_re) && is.null(sigma2_re)) ||
      (  is.null(ID_re) && !is.null(sigma2_re))) {
    stop("To add unstructured random effects, provide both `ID_re` and `sigma2_re`.")
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

  make_link_funs <- function(family, invlink, ncheck) {
    have_Deriv <- requireNamespace("Deriv", quietly = TRUE)
    have_numDeriv <- requireNamespace("numDeriv", quietly = TRUE)

    # Canonical defaults
    if (is.null(invlink)) {
      if (family == "poisson") {
        inv <- function(x) exp(x)
        d1  <- function(x) exp(x)
        d2  <- function(x) exp(x)
      } else {
        inv <- function(x) plogis(x)
        d1  <- function(x) { p <- inv(x); p * (1 - p) }
        d2  <- function(x) { p <- inv(x); d <- p * (1 - p); d * (1 - 2 * p) }
      }
      check_vec_fun(inv, ncheck, "canonical invlink")
      check_vec_fun(d1,  ncheck, "canonical invlink_prime")
      check_vec_fun(d2,  ncheck, "canonical invlink_second")
      return(list(inv = inv, d1 = d1, d2 = d2, name = "canonical"))
    }

    # If a bare function is given, treat it as the inverse link
    if (is.function(invlink)) {
      inv_user <- invlink
      d1_user <- NULL
      d2_user <- NULL
    } else if (is.list(invlink)) {
      inv_user <- invlink$inv %||% invlink$inv_link %||% invlink$invlink
      d1_user  <- invlink$d1 %||% invlink$inv_link_prime %||% invlink$mu_eta
      d2_user  <- invlink$d2 %||% invlink$inv_link_second
    } else {
      stop("'invlink' must be NULL, a function, or a list with components inv, d1, d2")
    }

    # Validate the inverse link
    check_vec_fun(inv_user, ncheck, "invlink")

    # Obtain missing derivatives once
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

    if (is.null(d2_user)) {
      if (have_Deriv) {
        d2_user <- Deriv::Deriv(d1_user, "eta")
      } else if (have_numDeriv) {
        d2_user <- function(eta) vapply(eta, function(z)
          numDeriv::grad(function(x) d1_user(x), z), numeric(1))
      } else {
        stop("Cannot auto-derive second derivative. Install `Deriv` or `numDeriv`, or provide `d2`.")
      }
    }
    check_vec_fun(d2_user, ncheck, "invlink_second")

    list(inv = inv_user, d1 = d1_user, d2 = d2_user, name = "custom")
  }

  linkf <- make_link_funs(family, invlink, n)
  inv_link <- linkf$inv
  inv1     <- linkf$d1
  inv2     <- linkf$d2

  # Per-observation contributions for arbitrary inverse link
  contribs <- function(eta) {
    if (family == "poisson") {
      mu  <- inv_link(eta)
      mu1 <- inv1(eta)
      mu2 <- inv2(eta)

      if (any(!is.finite(mu))) stop("invlink produced non-finite values for Poisson.")
      if (any(mu <= 0)) stop("invlink must return strictly positive means for Poisson.")

      g  <- (y / mu - units_m) * mu1
      l2 <- - y * (mu1^2) / (mu^2) + (y / mu - units_m) * mu2
      w  <- -l2

      llik <- sum(y * log(pmax(units_m * mu, .Machine$double.eps)) - units_m * mu)
      list(g = g, w = w, llik = llik)
    } else { # binomial
      p  <- inv_link(eta)
      p1 <- inv1(eta)
      p2 <- inv2(eta)

      if (any(!is.finite(p))) stop("invlink produced non-finite values for Binomial.")
      if (any(p <= 0 | p >= 1)) stop("invlink must return values in (0,1) for Binomial.")

      den <- p * (1 - p)
      g  <- (y - units_m * p) * p1 / den
      l2 <- - units_m * (p1^2) / den + (y - units_m * p) * ( p2 / den - (p1^2) * (1 - 2 * p) / (den^2) )
      w  <- -l2

      llik <- sum(y * log(pmax(p, .Machine$double.eps)) +
                    (units_m - y) * log(pmax(1 - p, .Machine$double.eps)))
      list(g = g, w = w, llik = llik)
    }
  }

  # ---------- core functions ----------
  integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]
    S_re_list <- if (n_re > 0) lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]]) else NULL

    eta <- mu + S[ID_coords]
    if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]

    cw <- contribs(eta)

    qS  <- as.numeric(crossprod(S, Sigma.inv %*% S))
    qre <- if (n_re > 0) sum(vapply(seq_len(n_re), function(i) sum(S_re_list[[i]]^2) / sigma2_re[i], numeric(1))) else 0

    -0.5 * qS - 0.5 * qre + cw$llik
  }

  grad.integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]
    S_re_list <- if (n_re > 0) lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]]) else NULL

    eta <- mu + S[ID_coords]
    if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]

    cw <- contribs(eta)
    g  <- cw$g

    out <- numeric(n_tot)
    out[1:n_loc] <- as.numeric(-Sigma.inv %*% S) + sum_by_group(g, ID_coords, n_loc)
    if (n_re > 0) {
      for (j in seq_len(n_re)) {
        out[ind_re[[j]]] <- - S_re_list[[j]] / sigma2_re[j] +
          sum_by_group(g, ID_re[, j], n_dim_re[j])
      }
    }
    out
  }

  hessian.integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]
    S_re_list <- if (n_re > 0) lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]]) else NULL

    eta <- mu + S[ID_coords]
    if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]

    cw <- contribs(eta)
    w  <- cw$w  # positive if l'' is negative

    H <- matrix(0, nrow = n_tot, ncol = n_tot)

    # S block
    H[1:n_loc, 1:n_loc] <- -Sigma.inv
    diag(H)[1:n_loc] <- diag(H)[1:n_loc] - sum_by_group(w, ID_coords, n_loc)

    # RE blocks
    if (n_re > 0) {
      for (j in seq_len(n_re)) {
        idx <- ind_re[[j]]
        # diagonal block j
        diag(H)[idx] <- diag(H)[idx] - 1 / sigma2_re[j] - sum_by_group(w, ID_re[, j], n_dim_re[j])

        # cross S vs RE j
        M <- cross_sum(w, ID_coords, n_loc, ID_re[, j], n_dim_re[j]) # n_loc x n_dim_re[j]
        H[1:n_loc, idx] <- H[1:n_loc, idx] - M
        H[idx, 1:n_loc] <- t(H[1:n_loc, idx])

        # cross RE j vs RE k
        if (j < n_re) {
          for (k in (j + 1):n_re) {
            Mk <- cross_sum(w, ID_re[, j], n_dim_re[j], ID_re[, k], n_dim_re[k])
            idxk <- ind_re[[k]]
            H[idx, idxk] <- H[idx, idxk] - Mk
            H[idxk, idx] <- t(H[idx, idxk])
          }
        }
      }
    }
    H
  }

  # ---------- optimize ----------
  estim <- nlminb(
    rep(0, n_tot),
    function(x) -integrand(x),
    function(x) -grad.integrand(x),
    function(x) -hessian.integrand(x)
  )

  out <- list(mode = estim$par, invlink_used = linkf$name)
  if (hessian) {
    out$hessian <- hessian.integrand(out$mode)
  } else {
    out$Sigma.tilde <- solve(-hessian.integrand(out$mode))
  }
  if (gradient) out$gradient <- grad.integrand(out$mode)

  out
}

