##' @title Prediction of the random effects components and covariates effects over a spatial grid
##' @description Computes predictions over a spatial grid using a fitted model from
##'   \code{\link{glgpm}}.
##' @param object A RiskMap object.
##' @param grid_pred An \code{sfc} or \code{sf} of POINT geometries, or a list thereof for joint predictions.
##' If not provided, the predictions will be generated at the geometries provided when fitting the model.
##' @param predictors Optional dataframe or list of dataframes containing predictor variables at prediction locations.
##' Must be provided if you specify `grid_pred`.
##' @param re_predictors Optional dataframe containing random effect predictors.
##' Not supported if `grid_pred` is a list.
##' @param pred_cov_offset Optional numeric vector containing covariate offsets at prediction locations.
##' Must be provided if there is an offset included in the model and not supported if `grid_pred` is a list.
##' @param control_sim Control parameters from \code{\link{set_control_mcmc}}.
##' @param type Whether the predictions are `marginal` or `joint`. `marginal` predictions are less
##' computationally expensive than `joint` predictions but cannot be used to predict areal targets.
##' If `grid_pred` is a list or random effects are included, must be set to `joint`. Defaults to `marginal`.
##' @param messages Logical; display progress messages. Defaults to `TRUE`.
##' @return An object of class \code{"RiskMap_pred"} containing:
##'   \describe{
##'     \item{mu_pred}{Fixed-effects component of the linear predictor at the
##'       prediction locations, i.e. \eqn{D_{pred}\hat{\beta}}. A numeric vector
##'       of length \code{n_pred}, or a list of such vectors when \code{grid_pred}
##'       was supplied as a list of grids. Equal to \code{0} when the model
##'       contains no covariates.}
##'     \item{grid_pred}{The locations of the predictions}
##'     \item{par_hat}{Named list of the maximum likelihood estimates returned by
##'       \code{coef()} on the fitted \code{\link{glgpm}} object, containing
##'       \code{beta} (regression coefficients), \code{sigma2} (spatial variance),
##'       \code{phi} (scale of the spatial correlation) and, where applicable,
##'       \code{tau2} (nugget), \code{sigma2_re} (variances of the unstructured
##'       random effects) and \code{sigma2_me} (measurement error variance).}
##'     \item{S_samples}{Samples from the predictive distribution of the spatial
##'       Gaussian process at the prediction locations. A matrix with
##'       \code{n_pred} rows and \code{n_samples} columns, or a list of such
##'       matrices (one per grid) when \code{grid_pred} was supplied as a list.}
##'     \item{re}{List with two elements describing the unstructured random
##'       effects: \code{D_pred}, a list of design matrices mapping the
##'       prediction locations onto the levels of each random effect, and
##'       \code{samples}, a nested list giving, for each random effect and each
##'       of its levels, a vector of \code{n_samples} draws. Both elements are
##'       \code{NULL} when the model contains no unstructured random effects.}
##'     \item{obs_loc}{Logical; \code{TRUE} when predictions were made at the
##'       observed data locations, i.e. when \code{grid_pred} was left as
##'       \code{NULL} in \code{\link{setup_prediction}}.}
##'     \item{inter_f}{The model formula after interpretation by
##'       \code{interpret.formula}, separating the fixed-effects terms, the
##'       spatial term and the unstructured random effect terms. Used internally
##'       to build the linear predictor at the prediction locations.}
##'     \item{family}{The model family}
##'     \item{cov_offset}{Covariate offsets}
##'     \item{type}{The type of predictions - `marginal` or `joint`}
##'   }
##' @importFrom Matrix solve
##' @examples
##'
##' data(italy_sim)
##' italy_subset <- italy_sim[1:100,]
##'
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_subset,
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' # using locations in dataset
##' prediction_setup <- setup_prediction(fit)
##'
##' # using new locations
##' hull <- create_convex_hull(italy_subset)
##' grid_pred <- create_grid(hull, 20)
##' prediction_setup <- setup_prediction(
##'   fit,
##'   grid_pred = grid_pred,
##'   predictors = data.frame(y = rnorm(length(grid_pred)))
##' )
##'
##' @export
setup_prediction <- function(object,
                             grid_pred = NULL,
                             predictors = NULL,
                             re_predictors = NULL,
                             pred_cov_offset = NULL,
                             control_sim = set_control_mcmc(),
                             type = "marginal",
                             messages = TRUE) {

  # ---------------------------------------------------------------------------
  # validate inputs
  # ---------------------------------------------------------------------------
  stopifnot("'object' must be of class RiskMap" = inherits(object, "RiskMap"))

  list_mode <- inherits(grid_pred, "list")
  if (list_mode) {
    if (type != "joint")
      stop("When 'grid_pred' is a list, 'type' must be 'joint'")
    if (length(grid_pred) == 0L)
      stop("'grid_pred' is a list but has length 0")
    tryCatch(
      lapply(grid_pred, check_data, type = "sfc"),
      error = function(e){
        stop("Each element of 'grid_pred' must be an 'sf' or 'sfc' object with POINT geometries")
      }
    )
  } else {
    if (!is.null(grid_pred))
      check_data(grid_pred, type = "sfc")
  }

  if (!inherits(control_sim, "RiskMap_control_mcmc"))
    stop("'control_sim' must be an output from 'set_control_mcmc()'")

  if (!type %in% c("marginal", "joint"))
    stop("'type' must be either 'marginal' or 'joint'")

  if (!is.null(grid_pred) && is.null(predictors))
    stop("'predictors' must be supplied if 'grid_pred' is supplied")

  obs_loc <- is.null(grid_pred)
  if (obs_loc) {
    if (!is.null(predictors))
      warning("You have set 'predictors' but not 'grid_pred' so 'predictors' will be ignored")
    predictors <- as.data.frame(st_drop_geometry(object$data_sf))
    grid_pred  <- st_as_sfc(object$data_sf)
  } else {
    if (list_mode) {
      grid_pred <- lapply(grid_pred, st_transform, crs = object$crs)
    } else {
      grid_pred <- st_transform(grid_pred, crs = object$crs)
    }
  }

  if (list_mode) {
    grp    <- lapply(grid_pred, st_coordinates)
    n_pred <- vapply(grp, nrow, integer(1))
  } else {
    grp    <- st_coordinates(grid_pred)
    n_pred <- nrow(grp)
  }

  object$D <- as.matrix(object$D)
  p <- ncol(object$D)

  if (!all(object$cov_offset == 0)) {
    if (obs_loc){
      pred_cov_offset <- object$cov_offset
    } else {
      if (list_mode)
        stop("Predictions including covariate offsets are not yet supported when 'pred_grid' is a list")
      if (is.null(pred_cov_offset))
        stop("'pred_cov_offset' must be specified at each prediction location")
      if (!inherits(pred_cov_offset, "numeric"))
        stop("'pred_cov_offset' must be a numeric vector")
      if (length(pred_cov_offset) != n_pred)
        stop("The length of 'pred_cov_offset' does not match the number of prediction locations")
    }
  } else {
    if (!is.null(pred_cov_offset))
      warning("You have set 'pred_cov_offset' but 'object' does not contain a cov_offset so this will be ignored")
    pred_cov_offset <- 0
  }

  # ---------------------------------------------------------------------------
  # Extract terms from object
  # ---------------------------------------------------------------------------

  par_hat <- coef(object)
  inter_f      <- interpret.formula(object$formula)
  inter_lt_f   <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~.)

  intercept_only <- if (p == 1) prod(object$D[, 1] == 1) == 1 else FALSE

  n_re <- length(object$re)
  if (n_re > 0 && type == "marginal" && !is.null(re_predictors))
    stop("Random effect predictions require 'type' to be set to 'joint'")

  if (!is.null(re_predictors) && list_mode)
    stop("Prediction of random effects with a list of prediction grids ('grid_pred') is ",
         "not yet supported - supply a single 'grid_pred' or omit 're_predictors'")

  # ---------------------------------------------------------------------------
  # Build mu_pred (fixed-effects linear predictor at prediction locations)
  # ---------------------------------------------------------------------------

  if (!is.null(predictors)) {

    .build_D_pred <- function(predictors, n_predictors, index = NULL) {
      response <- as.character(object$formula[[2]])
      re_names <- names(object$re)
      offset_names <- names(object$cov_offset)
      model_predictors <- setdiff(get_formula_terms(object$formula), c(response, re_names, offset_names))
      if (!is.data.frame(predictors))
        stop(if (is.null(index)) "'predictors' must be a data.frame"
             else sprintf("'predictors[[%d]]' must be a data.frame", index))
      if (!all(model_predictors %in% names(predictors)))
        stop("The column names in 'predictors' do not match the variables in the model formula")
      if (nrow(predictors) != n_predictors)
        stop(if (is.null(index)) "The number of rows in 'predictors' does not match the number of locations in 'grid_pred'"
             else sprintf("The number of rows in of 'predictors[[%d]]' does not match the number of locations in 'grid_pred[[%d]]'", index, index))
      mf <- model.frame(inter_lt_f$pf, data = predictors, na.action = na.fail)
      as.matrix(model.matrix(attr(mf, "terms"), data = predictors))
    }

    if (list_mode) {
      D_pred  <- mapply(.build_D_pred, predictors, n_pred, seq_along(predictors), SIMPLIFY = FALSE)
      mu_pred <- lapply(D_pred, function(D) as.numeric(D %*% par_hat$beta))
    } else {
      D_pred  <- .build_D_pred(predictors, n_pred)
      mu_pred <- as.numeric(D_pred %*% par_hat$beta)
    }

  }

  # ---------------------------------------------------------------------------
  # Random effects (unchanged from original)
  # ---------------------------------------------------------------------------
  if (n_re > 0) {
    ID_g         <- as.matrix(cbind(object$ID_coords, object$ID_re))
    re_unique    <- object$re
    n_dim_re_tot <- sapply(seq_len(n_re + 1), function(i) length(unique(ID_g[, i])))

    if (!is.null(re_predictors)) {
      if (any(is.na(re_predictors))) {
        warning("Missing values in 're_predictors'; removing affected locations")
        ind_c       <- complete.cases(re_predictors)
        re_predictors <- re_predictors[ind_c, , drop = FALSE]
        grid_pred   <- if (list_mode) lapply(grid_pred, `[`, ind_c) else grid_pred[ind_c]
        grp         <- if (list_mode) lapply(grid_pred, st_coordinates) else st_coordinates(grid_pred)
        n_pred      <- if (list_mode) vapply(grp, nrow, integer(1)) else nrow(grp)
      }
      if (!is.data.frame(re_predictors)) stop("'re_predictors' must be a data.frame")
      if (nrow(re_predictors) != n_pred)  stop("'re_predictors' rows do not match 'grid_pred'")
      if (ncol(re_predictors) != n_re)    stop("'re_predictors' columns do not match number of random effects")

      n_dim_re <- sapply(seq_len(n_re), function(i) length(unique(object$ID_re[, i])))
      D_re_pred <- lapply(seq_len(n_re), function(i) {
        mat <- matrix(0, nrow = n_pred, ncol = n_dim_re[i])
        for (j in seq_along(object$re[[i]]))
          mat[re_predictors[, i] == object$re[[i]][j], j] <- 1
        mat
      })
    } else {
      D_re_pred <- NULL
    }
  } else {
    D_re_pred <- NULL
  }

  # ---------------------------------------------------------------------------
  # Spatial quantities
  # ---------------------------------------------------------------------------
  out <- list(mu_pred = mu_pred, grid_pred = grid_pred, par_hat = par_hat)

  if (object$scale_to_km) {
    grp <- if (list_mode) lapply(grp, function(g) g / 1000) else grp / 1000
  }

  if (object$family != "gaussian" && !obs_loc) {
    if (list_mode) {
      U_pred <- lapply(grp, function(g) {
        t(sapply(seq_len(nrow(g)), function(i)
          sqrt((object$coords[, 1] - g[i, 1])^2 + (object$coords[, 2] - g[i, 2])^2)))
      })
    } else {
      U_pred <- t(sapply(seq_len(n_pred), function(i)
        sqrt((object$coords[, 1] - grp[i, 1])^2 + (object$coords[, 2] - grp[i, 2])^2)))
    }
  } else if (object$family == "gaussian" && !obs_loc) {
    if (list_mode) {
      U_pred <- lapply(grp, function(g) {
        t(sapply(seq_len(nrow(g)), function(i)
          sqrt((object$coords[object$ID_coords, 1] - g[i, 1])^2 +
                 (object$coords[object$ID_coords, 2] - g[i, 2])^2)))
      })
    } else {
      U_pred <- t(sapply(seq_len(n_pred), function(i)
        sqrt((object$coords[object$ID_coords, 1] - grp[i, 1])^2 +
               (object$coords[object$ID_coords, 2] - grp[i, 2])^2)))
    }
  }

  U <- dist(object$coords)
  R <- matern_correlation(U, phi = par_hat$phi, kappa = object$kappa, return_sym_matrix = TRUE)

  if (!obs_loc) {
    C <- if (list_mode)
      lapply(U_pred, function(u) par_hat$sigma2 * matern_correlation(u, phi = par_hat$phi, kappa = object$kappa))
    else
      par_hat$sigma2 * matern_correlation(U_pred, phi = par_hat$phi, kappa = object$kappa)
  } else {
    C <- par_hat$sigma2 * R[, object$ID_coords]
    grp <- object$coords
  }

  n_pred_spatial <- if (obs_loc) nrow(object$coords) else n_pred

  mu <- as.numeric(object$D %*% par_hat$beta)

  n_samples <- if (control_sim$linear_model) control_sim$n_sim else (control_sim$n_sim - control_sim$burnin) / control_sim$thin

  # ---------------------------------------------------------------------------
  # FIX 2: nu2 / nugget
  # STH: tau2 fixed at 0 -> near-zero for numerical stability
  # LF:  tau2 may be estimated -> use it
  # glgpm: existing behaviour
  # ---------------------------------------------------------------------------
  nu2 <- if (!is.null(object$fix_tau2)) object$fix_tau2 / par_hat$sigma2 else par_hat$tau2 / par_hat$sigma2
  if (nu2 == 0) nu2 <- 1e-10
  diag(R) <- diag(R) + nu2

  diff.y <- if (object$family == "gaussian") object$y - mu else NULL

  # ===========================================================================
  # NON-GAUSSIAN MODELS
  # ===========================================================================
  if (object$family != "gaussian") {

    Sigma     <- par_hat$sigma2 * R
    Sigma_inv <- solve(Sigma)

    if (!obs_loc) {
      A <- if (list_mode) lapply(C, function(Ci) Ci %*% Sigma_inv)
      else C %*% Sigma_inv
    }

    simulation <- laplace_sampling_mcmc(
      y = object$y, units_m = object$units_m, mu = mu, Sigma = Sigma,
      sigma2_re = par_hat$sigma2_re, invlink = object$linkf,
      ID_coords = object$ID_coords, ID_re = object$ID_re,
      family = object$family, control_mcmc = control_sim, messages = messages)

    if (obs_loc) {
      out$S_samples <- t(simulation$samples$S)
    } else {
      mu_cond_S <- if (list_mode)
        lapply(A, function(Ai) Ai %*% t(simulation$samples$S))
      else
        A %*% t(simulation$samples$S)

      if (type == "marginal") {
        sd_cond_S <- sqrt(par_hat$sigma2 - diag(A %*% t(C)))
        out$S_samples <- sapply(seq_len(n_samples), function(i)
          mu_cond_S[, i] + sd_cond_S * rnorm(n_pred))

      } else {
        if (list_mode) {
          out$S_samples <- lapply(seq_along(mu_cond_S), function(i) {
            Sp    <- par_hat$sigma2 * matern_correlation(dist(grp[[i]]), phi = par_hat$phi,
                                                         kappa = object$kappa, return_sym_matrix = TRUE)
            Sc    <- Sp - A[[i]] %*% t(C[[i]])
            Scr   <- t(chol(Sc))
            sapply(seq_len(n_samples), function(j)
              mu_cond_S[[i]][, j] + Scr %*% rnorm(nrow(mu_cond_S[[i]])))
          })
        } else {
          Sp  <- par_hat$sigma2 * matern_correlation(dist(grp), phi = par_hat$phi,
                                                     kappa = object$kappa, return_sym_matrix = TRUE)
          Sc  <- Sp - A %*% t(C)
          Scr <- t(chol(Sc))
          out$S_samples <- sapply(seq_len(n_samples), function(i)
            mu_cond_S[, i] + Scr %*% rnorm(nrow(mu_cond_S)))
        }
      }
    }

  } else {
    # =========================================================================
    # GAUSSIAN MODEL (unchanged)
    # =========================================================================
    if (!is.null(object$fix_var_me) && object$fix_var_me > 0 || is.null(object$fix_var_me)) {
      m        <- length(object$y)
      s_unique <- unique(object$ID_coords)
      ID_g     <- as.matrix(cbind(object$ID_coords, object$ID_re))
      n_dim_re_tot <- sapply(seq_len(n_re + 1), function(i) length(unique(ID_g[, i])))
      C_g      <- matrix(0, nrow = m, ncol = sum(n_dim_re_tot))
      for (i in seq_len(m)) {
        ind_s_i <- which(s_unique == ID_g[i, 1])
        C_g[i, seq_len(n_dim_re_tot[1])][ind_s_i] <- 1
      }
      if (n_re > 0) {
        for (j in seq_len(n_re)) {
          select_col <- sum(n_dim_re_tot[seq_len(j)])
          for (i in seq_len(m)) {
            ind_re_j_i <- which(re_unique[[j]] == ID_g[i, j + 1])
            C_g[i, select_col + seq_len(n_dim_re_tot[j + 1])][ind_re_j_i] <- 1
          }
        }
      }
      C_g      <- Matrix(C_g, sparse = TRUE, doDiag = FALSE)
      C_g_m    <- forceSymmetric(Matrix::t(C_g) %*% C_g)
      Sigma_g  <- matrix(0, nrow = sum(n_dim_re_tot), ncol = sum(n_dim_re_tot))
      Sigma_g_inv <- matrix(0, nrow = sum(n_dim_re_tot), ncol = sum(n_dim_re_tot))
      Sigma_g[seq_len(n_dim_re_tot[1]), seq_len(n_dim_re_tot[1])] <- par_hat$sigma2 * R
      Sigma_g_inv[seq_len(n_dim_re_tot[1]), seq_len(n_dim_re_tot[1])] <- solve(R) / par_hat$sigma2
      if (n_re > 0) {
        for (j in seq_len(n_re)) {
          sc <- sum(n_dim_re_tot[seq_len(j)])
          diag(Sigma_g[sc + seq_len(n_dim_re_tot[j+1]), sc + seq_len(n_dim_re_tot[j+1])]) <- par_hat$sigma2_re[j]
          diag(Sigma_g_inv[sc + seq_len(n_dim_re_tot[j+1]), sc + seq_len(n_dim_re_tot[j+1])]) <- 1 / par_hat$sigma2_re[j]
        }
      }
      Sigma_star     <- Sigma_g_inv + C_g_m / par_hat$sigma2_me
      Sigma_star_inv <- forceSymmetric(Matrix::solve(Sigma_star))
      B    <- -C_g %*% Sigma_star_inv %*% Matrix::t(C_g) / (par_hat$sigma2_me^2)
      diag(B) <- Matrix::diag(B) + 1 / par_hat$sigma2_me

      A <- if (list_mode) lapply(C, function(single_grid_C) single_grid_C %*% B) else C %*% B
    } else {
      Sigma     <- par_hat$sigma2 * R
      Sigma_inv <- solve(Sigma)
      A <- if (list_mode) lapply(C, function(single_grid_C) single_grid_C %*% Sigma_inv) else C %*% Sigma_inv
    }

    mu_cond_S <- if (list_mode) {
      lapply(A, function(single_grid_A) as.numeric(single_grid_A %*% diff.y))
    } else {
      as.numeric(A %*% diff.y)
    }

    if (type == "marginal") {
      if (list_mode) {
        out$S_samples <- lapply(seq_along(A), function(i) {
          sd_cond_S_i <- sqrt(par_hat$sigma2 - Matrix::diag(A[[i]] %*% t(C[[i]])))
          sapply(seq_len(n_samples), function(j)
            mu_cond_S[[i]] + sd_cond_S_i * rnorm(n_pred[i]))
        })
      } else {
        sd_cond_S <- sqrt(par_hat$sigma2 - Matrix::diag(A %*% t(C)))
        out$S_samples <- sapply(seq_len(n_samples), function(i)
          mu_cond_S + sd_cond_S * rnorm(n_pred_spatial))
      }
    } else {
      if (list_mode) {
        out$S_samples <- lapply(seq_along(A), function(i) {
          spatial_covariance_i <- par_hat$sigma2 * matern_correlation(dist(grp[[i]]), phi = par_hat$phi,
                                                                      kappa = object$kappa, return_sym_matrix = TRUE)
          conditional_covariance_i <- spatial_covariance_i - A[[i]] %*% t(C[[i]])
          cholesky_root_i <- t(chol(conditional_covariance_i))
          sapply(seq_len(n_samples), function(j)
            mu_cond_S[[i]] + cholesky_root_i %*% rnorm(n_pred[i]))
        })
      } else {
        Sp  <- par_hat$sigma2 * matern_correlation(dist(grp), phi = par_hat$phi,
                                                   kappa = object$kappa, return_sym_matrix = TRUE)
        Sc  <- Sp - A %*% t(C)
        Scr <- t(chol(Sc))
        out$S_samples <- sapply(seq_len(n_samples), function(i)
          mu_cond_S + Scr %*% rnorm(n_pred_spatial))
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Random effect posterior samples (unchanged from original)
  # ---------------------------------------------------------------------------

  if (n_re > 0 && !is.null(re_predictors)) {
    out$re          <- list()
    out$re$D_pred   <- D_re_pred
    out$re$samples  <- list()
    re_names        <- colnames(object$ID_re)
    if (object$family == "gaussian") {
      Sigma_cond_inv <- solve(Sc)
      C_Z  <- C_g[, -(seq_len(n_dim_re_tot[1]))]
      add  <- 0
      for (i in seq_along(n_dim_re_tot[-1])) {
        C_Z[, add + seq_len(n_dim_re_tot[i+1])] <- par_hat$sigma2_re[i] *
          C_Z[, add + seq_len(n_dim_re_tot[i+1])]
        add <- n_dim_re_tot[i+1]
      }
      A_Z          <- Matrix::t(C_Z) %*% B %*% t(C) %*% Sigma_cond_inv
      Sigma_Z_cond <- diag(rep(par_hat$sigma2_re, n_dim_re_tot[-1])) -
        Matrix::t(C_Z) %*% B %*% C_Z -
        A_Z %*% C %*% Matrix::t(B) %*% C_Z
      Scr_Z        <- t(chol(Sigma_Z_cond))
      mu_Z_cond    <- sapply(seq_len(n_samples), function(i)
        as.matrix(A_Z %*% (out$S_samples[, i] - mu_cond_S)))
      add <- 0
      for (i in seq_len(n_re)) {
        for (j in seq_len(n_dim_re_tot[1 + i])) {
          add <- add + 1
          ind_ij   <- which(object$ID_re[[i]] == re_unique[[i]][j])
          C_re_ij  <- matrix(0, ncol = m)
          C_re_ij[, ind_ij] <- par_hat$sigma2_re[i]
          mu_Z_cond[add, ] <- as.numeric(C_re_ij %*% B %*% diff.y) + mu_Z_cond[add, ]
        }
      }
      re_samples <- sapply(seq_len(n_samples), function(i)
        as.numeric(mu_Z_cond[, i] + Scr_Z %*% rnorm(sum(n_dim_re_tot[-1]))))
    } else {
      re_samples <- matrix(0, nrow = sum(n_dim_re_tot[-1]), ncol = n_samples)
      add <- 0
      for (i in seq_len(n_re)) {
        re_samples[seq_len(n_dim_re_tot[i+1]), ] <- t(simulation$samples[[i+1]])
        add <- add + n_dim_re_tot[i+1]
      }
    }
    add <- 0
    for (i in seq_len(n_re)) {
      for (j in seq_len(n_dim_re_tot[1 + i])) {
        add <- add + 1
        out$re$samples[[re_names[[i]]]][[as.character(re_unique[[i]][[j]])]] <- re_samples[add, ]
      }
    }
  } else {
    out$re <- list(D_pred = NULL, samples = NULL)
  }

  out$obs_loc  <- obs_loc
  if (obs_loc){
    out$ID_coords <- object$ID_coords
  } else {
    out["ID_coords"] <- list(NULL)
  }
  out$inter_f  <- inter_f
  out$family   <- object$family
  out$par_hat  <- par_hat
  out$cov_offset <- pred_cov_offset
  out$type     <- type
  class(out)   <- "RiskMap_pred"
  return(out)
}


##' @title Update Predictors for a RiskMap Prediction Object
##'
##' @description
##' This function updates the predictors of a given RiskMap prediction object. It ensures that the new predictors match the original prediction grid and updates the relevant components of the object accordingly.
##'
##' @param object A `RiskMap_pred` object, which is the output of the \code{\link{setup_prediction}} function.
##' @param predictors A data frame containing the new predictor values. The number of rows must match the prediction grid in the `object`.
##'
##' @details
##' The function performs several checks and updates:
##' \itemize{
##'   \item Ensures that `object` is of class `RiskMap_pred`.
##'   \item Ensures that the number of rows in `predictors` matches the prediction grid in `object`.
##'   \item Removes any rows with missing values in `predictors` and updates the corresponding components of the `object`.
##'   \item Updates the prediction locations, the predictive samples for the random effects, and the linear predictor.
##' }
##'
##' @return The updated `RiskMap_pred` object.
##' @export
update_predictors <- function(object, predictors) {
  if (!inherits(object, what = "RiskMap_pred", which = FALSE)) {
    stop("The object passed to 'object' must be an output of the function 'setup_prediction'")
  }

  list_mode <- is.list(object$grid_pred) && !is.null(object$grid_pred) &
    !any(class(object$grid_pred)=="sf" | class(object$grid_pred)=="sfc")
  par_hat <- object$par_hat
  p <- length(par_hat$beta)

  if (p == 1) {
    stop("No update of the predictors can be done for an intercept-only model")
  }

  inter_f <- object$inter_f
  inter_lt_f <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~ .)

  if (list_mode) {
    if (!is.list(predictors)) stop("If object$grid_pred is a list, then 'predictors' must also be a list")
    if (length(predictors) != length(object$grid_pred)) stop("Length of 'predictors' list must match length of 'grid_pred' list")

    n_pred_list <- lapply(object$grid_pred, function(x) nrow(st_coordinates(x)))
    for (i in seq_along(predictors)) {
      if (!is.data.frame(predictors[[i]])) stop(sprintf("predictors[[%d]] must be a data.frame", i))
      if (nrow(predictors[[i]]) != n_pred_list[i]) stop(sprintf("predictors[[%d]] does not match the number of prediction locations", i))
    }

    if (!is.null(object$re_predictors)) {
      for (i in seq_along(predictors)) {
        comb_pred <- data.frame(predictors[[i]], object$re_predictors[[i]])
        ind_c <- complete.cases(comb_pred)

        predictors_aux <- predictors[[i]][ind_c, , drop = FALSE]
        predictors[[i]] <- predictors_aux

        re_predictors_aux <- object$re_predictors[[i]][ind_c, , drop = FALSE]
        object$re_predictors[[i]] <- re_predictors_aux

        n_re <- length(object$re$samples)
        for (r in 1:n_re) {
          object$re$D_pred[[r]][[i]] <- object$re$D_pred[[r]][[i]][ind_c, , drop = FALSE]
        }

        object$grid_pred[[i]] <- object$grid_pred[[i]][ind_c, ]
        if (is.matrix(object$S_samples[[i]])) {
          object$S_samples[[i]] <- object$S_samples[[i]][ind_c, , drop = FALSE]
        } else {
          object$S_samples[[i]] <- object$S_samples[[i]][ind_c]
        }
      }
    } else {
      for (i in seq_along(predictors)) {
        comb_pred <- predictors[[i]]
        ind_c <- complete.cases(comb_pred)
        predictors[[i]] <- predictors[[i]][ind_c, , drop = FALSE]

        object$grid_pred[[i]] <- object$grid_pred[[i]][ind_c, ]
        if (is.matrix(object$S_samples[[i]])) {
          object$S_samples[[i]] <- object$S_samples[[i]][ind_c, , drop = FALSE]
        } else {
          object$S_samples[[i]] <- object$S_samples[[i]][ind_c]
        }
      }
    }

    object$mu_pred <- vector("list", length(predictors))
    for (i in seq_along(predictors)) {
      mf_pred <- model.frame(inter_lt_f$pf, data = predictors[[i]], na.action = na.fail)
      D_pred <- as.matrix(model.matrix(attr(mf_pred, "terms"), data = predictors[[i]]))
      if (ncol(D_pred) != p) stop("Mismatch in number of predictors")
      object$mu_pred[[i]] <- as.numeric(D_pred %*% par_hat$beta)
    }
  } else {
    grp <- st_coordinates(object$grid_pred)
    n_pred <- nrow(grp)

    if (!is.data.frame(predictors)) stop("'predictors' must be an object of class 'data.frame'")
    if (nrow(predictors) != n_pred) stop("the values provided for 'predictors' do not match the prediction grid")

    if (any(is.na(predictors))) {
      warning("There are missing values in 'predictors'; these values have been removed alongside the corresponding prediction locations")
    }

    if (!is.null(object$re_predictors)) {
      comb_pred <- data.frame(predictors, object$re_predictors)
      ind_c <- complete.cases(comb_pred)
      predictors <- predictors[ind_c, , drop = FALSE]
      object$re_predictors <- object$re_predictors[ind_c, , drop = FALSE]

      n_re <- length(object$re$samples)
      for (i in 1:n_re) {
        object$re$D_pred[[i]] <- object$re$D_pred[[i]][ind_c, , drop = FALSE]
      }
    } else {
      comb_pred <- predictors
      ind_c <- complete.cases(comb_pred)
      predictors <- predictors[ind_c, , drop = FALSE]
    }

    object$grid_pred <- object$grid_pred[ind_c, ]
    grp <- st_coordinates(object$grid_pred)
    n_pred <- nrow(grp)
    object$S_samples <- object$S_samples[ind_c, , drop = FALSE]

    mf_pred <- model.frame(inter_lt_f$pf, data = predictors, na.action = na.fail)
    D_pred <- as.matrix(model.matrix(attr(mf_pred, "terms"), data = predictors))
    if (ncol(D_pred) != p) stop("Mismatch in number of predictors")
    object$mu_pred <- as.numeric(D_pred %*% par_hat$beta)
  }

  class(object) <- "RiskMap_pred"
  return(object)
}
