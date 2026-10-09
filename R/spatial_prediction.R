##' @title Prediction of the random effects components and covariates effects over a spatial grid
##' @description Computes predictions over a spatial grid using a fitted model from
##'   \code{\link{glgpm}}.
##' @param object A RiskMap object.
##' @param grid_pred An \code{sfc} or \code{sf} of POINT geometries, or a list
##' thereof for joint predictions. Its coordinate reference system (CRS) must
##' match the CRS used to fit \code{object}; transform it explicitly with
##' \code{sf::st_transform()} if required. If not provided, predictions are
##' generated at the observed locations.
##' @param predictors Optional dataframe or list of dataframes containing predictor variables at prediction locations.
##' Must be provided if you specify `grid_pred`.
##' @param re_predictors Optional dataframe containing random effect predictors.
##' Not supported if `grid_pred` is a list.
##' @param pred_cov_offset Optional numeric vector containing covariate offsets at prediction locations.
##' Must be provided if there is an offset included in the model and not supported if `grid_pred` is a list.
##' @param control_mcmc Control parameters from \code{\link{set_control_mcmc}}.
##' @param type Whether the predictions are `"marginal"` or `"joint"`. `"marginal"` predictions are less
##' computationally expensive than `"joint"` predictions but cannot be used to predict areal targets.
##' If `grid_pred` is a list or random effects are included, must be set to `"joint"`. Defaults to `"marginal"`.
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
##'       \code{interpret_formula}, separating the fixed-effects terms, the
##'       spatial term and the unstructured random effect terms. Used internally
##'       to build the linear predictor at the prediction locations.}
##'     \item{family}{The model family}
##'     \item{cov_offset}{Covariate offsets}
##'     \item{type}{The type of predictions - `"marginal"` or `"joint"`}
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
##'   predictors = data.frame(y = rnorm(nrow(grid_pred)))
##' )
##'
##' @export
setup_prediction <- function(object,
                           grid_pred = NULL,
                           predictors = NULL,
                           re_predictors = NULL,
                           pred_cov_offset = NULL,
                           control_mcmc = set_control_mcmc(),
                           type = "marginal",
                           messages = TRUE) {

  # ---------------------------------------------------------------------------
  # validate inputs
  # ---------------------------------------------------------------------------
  stopifnot("'object' must be of class RiskMap" = inherits(object, "RiskMap"))

  list_mode <- inherits(grid_pred, "list")
  model_crs <- st_crs(object$data)

  if (is.na(model_crs)) {
    stop(
      "The fitted model does not contain a valid coordinate reference system (CRS)",
      call. = FALSE
    )
  }

  validate_grid_element <- function(x, index) {
    tryCatch(
      check_data(x, type = "sfc"),
      error = function(e) {
        reason <- sub("^'[^']+' ", "", conditionMessage(e))
        stop(sprintf("'grid_pred[[%d]]' %s", index, reason), call. = FALSE)
      }
    )
  }

  crs_label <- function(x) {
    x_crs <- st_crs(x)
    if (!is.na(x_crs$epsg)) paste0("EPSG:", x_crs$epsg) else x_crs$input
  }

  if (list_mode) {
    if (type != "joint") {
      stop("When 'grid_pred' is a list, 'type' must be 'joint'")
    }
    if (length(grid_pred) == 0L) {
      stop("'grid_pred' is a list but has length 0")
    }
    invisible(Map(validate_grid_element, grid_pred, seq_along(grid_pred)))
    crs_mismatch <- vapply(
      grid_pred,
      function(x) !isTRUE(st_crs(x) == model_crs),
      logical(1)
    )
    if (any(crs_mismatch)) {
      mismatch_indices <- which(crs_mismatch)
      mismatch_labels <- vapply(
        grid_pred[mismatch_indices],
        crs_label,
        character(1)
      )
      mismatch_details <- paste0(
        mismatch_indices,
        " (",
        mismatch_labels,
        ")",
        collapse = ", "
      )
      stop(
        "The CRS of every element of 'grid_pred' must match the fitted model CRS (",
        crs_label(model_crs), "). Mismatches: ", mismatch_details,
        ". Transform these elements explicitly with sf::st_transform().",
        call. = FALSE
      )
    }
  } else {
    if (!is.null(grid_pred)) {
      check_data(grid_pred, type = "sfc")
      if (!isTRUE(st_crs(grid_pred) == model_crs)) {
        stop(
          "The CRS of 'grid_pred' (", crs_label(grid_pred),
          ") must match the fitted model CRS (", crs_label(model_crs),
          "). Transform 'grid_pred' explicitly with sf::st_transform().",
          call. = FALSE
        )
      }
    }
  }

  if (!inherits(control_mcmc, "RiskMap_control_mcmc"))
    stop("'control_mcmc' must be an output from 'set_control_mcmc()'")

  if (!is.null(control_mcmc$seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(control_mcmc$seed)
  }

  if (!type %in% c("marginal", "joint"))
    stop("'type' must be either 'marginal' or 'joint'")

  if (!is.null(grid_pred) && is.null(predictors) && ncol(object$D) > 1)
    stop("'predictors' must be supplied if 'grid_pred' is supplied")

  obs_loc <- is.null(grid_pred)
  if (obs_loc) {
    if (!is.null(predictors))
      warning("You have set 'predictors' but not 'grid_pred' so 'predictors' will be ignored")
    predictors <- as.data.frame(st_drop_geometry(object$data))
    grid_pred  <- st_as_sfc(object$data)
  }

  if (list_mode) {
    grp    <- lapply(grid_pred, coordinates_in_units,
                     distance_units = object$distance_units)
    n_pred <- vapply(grp, nrow, integer(1))
  } else {
    grp    <- coordinates_in_units(grid_pred, object$distance_units)
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
  inter_f      <- interpret_formula(object$formula)
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
      check_complete_data(predictors, model_predictors)
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

  } else if (intercept_only) {
    mu_pred <- if (list_mode) lapply(n_pred, function(n) rep(par_hat$beta, n)) else par_hat$beta
  } else {
    mu_pred <- 0
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
        grp         <- if (list_mode) {
          lapply(grid_pred, coordinates_in_units,
                 distance_units = object$distance_units)
        } else {
          coordinates_in_units(grid_pred, object$distance_units)
        }
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
  distance_scale <- attr(object, "distance_scale") %||% 1
  fitting_coords <- object$coords / distance_scale
  fitting_phi <- par_hat$phi / distance_scale
  fitting_grp <- if (list_mode) {
    lapply(grp, `/`, distance_scale)
  } else {
    grp / distance_scale
  }

  conditioning_coords <- if (object$family == "gaussian") {
    fitting_coords[object$ID_coords, , drop = FALSE]
  } else {
    fitting_coords
  }
  marginal_batch_size <- if (!list_mode && !obs_loc && type == "marginal") {
    marginal_prediction_batch_size(n_pred, nrow(conditioning_coords))
  } else {
    n_pred
  }
  batch_marginal <- !list_mode && !obs_loc && type == "marginal" &&
    marginal_batch_size < n_pred

  if (object$family != "gaussian" && !obs_loc && !batch_marginal) {
    if (list_mode) {
      U_pred <- lapply(fitting_grp, cross_distances, second = fitting_coords)
    } else {
      U_pred <- cross_distances(fitting_grp, fitting_coords)
    }
  } else if (object$family == "gaussian" && !obs_loc && !batch_marginal) {
    if (list_mode) {
      U_pred <- lapply(fitting_grp, cross_distances,
                       second = conditioning_coords)
    } else {
      U_pred <- cross_distances(fitting_grp, conditioning_coords)
    }
  }

  U <- pairwise_distances(fitting_coords)
  R <- matern_correlation(U, phi = fitting_phi, kappa = object$kappa,
                          return_sym_matrix = TRUE)

  if (obs_loc) {
    C <- par_hat$sigma2 * R[, object$ID_coords]
    grp <- object$coords
    fitting_grp <- fitting_coords
  } else if (!batch_marginal) {
    C <- if (list_mode)
      lapply(U_pred, function(u) par_hat$sigma2 * matern_correlation(u, phi = fitting_phi, kappa = object$kappa))
    else
      par_hat$sigma2 * matern_correlation(U_pred, phi = fitting_phi, kappa = object$kappa)
  }

  n_pred_spatial <- if (obs_loc) nrow(object$coords) else n_pred

  mu <- as.numeric(object$D %*% par_hat$beta)

  n_samples <- if (control_mcmc$linear_model) control_mcmc$n_sim else (control_mcmc$n_sim - control_mcmc$burnin) / control_mcmc$thin

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
    Sigma_root <- factor_covariance(Sigma, "observation covariance")

    prediction_weights <- function(cross_covariance) {
      cholesky_prediction_weights(cross_covariance, Sigma_root)
    }

    if (!obs_loc && !batch_marginal) {
      A <- if (list_mode) {
        lapply(C, prediction_weights)
      } else {
        prediction_weights(C)
      }
    }

    simulation <- laplace_sampling_mcmc(
      y = object$y, units_m = object$units_m, mu = mu, Sigma = Sigma,
      sigma2_re = par_hat$sigma2_re, invlink = object$linkf,
      ID_coords = object$ID_coords, ID_re = object$ID_re,
      family = object$family, control_mcmc = control_mcmc, messages = messages)

    if (obs_loc) {
      out$S_samples <- t(simulation$samples$S)
    } else if (batch_marginal) {
      out$S_samples <- batched_marginal_prediction(
        fitting_grp, conditioning_coords, prediction_weights,
        t(simulation$samples$S), par_hat$sigma2, fitting_phi,
        object$kappa, n_samples, marginal_batch_size
      )
    } else {
      mu_cond_S <- if (list_mode)
        lapply(A, function(Ai) Ai %*% t(simulation$samples$S))
      else
        A %*% t(simulation$samples$S)

      if (type == "marginal") {
        sd_cond_S <- sqrt(conditional_variances(par_hat$sigma2, A, C))
        out$S_samples <- sample_independent_gaussian(
          mu_cond_S, sd_cond_S, n_samples
        )

      } else {
        if (list_mode) {
          out$S_samples <- lapply(seq_along(mu_cond_S), function(i) {
            Sp    <- par_hat$sigma2 * matern_correlation(pairwise_distances(fitting_grp[[i]]), phi = fitting_phi,
                                                 kappa = object$kappa, return_sym_matrix = TRUE)
            Sc    <- Sp - A[[i]] %*% t(C[[i]])
            Scr   <- t(factor_covariance(
              Sc, "conditional prediction covariance"
            ))
            sample_correlated_gaussian(mu_cond_S[[i]], Scr, n_samples)
          })
        } else {
          Sp  <- par_hat$sigma2 * matern_correlation(pairwise_distances(fitting_grp), phi = fitting_phi,
                                             kappa = object$kappa, return_sym_matrix = TRUE)
          Sc  <- Sp - A %*% t(C)
          Sc_root <- factor_covariance(
            Sc, "conditional prediction covariance"
          )
          Scr <- t(Sc_root)
          out$S_samples <- sample_correlated_gaussian(
            mu_cond_S, Scr, n_samples
          )
        }
      }
    }

  } else {
    # =========================================================================
    # GAUSSIAN MODEL
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
      Sigma_g  <- matrix(0, nrow = sum(n_dim_re_tot), ncol = sum(n_dim_re_tot))
      Sigma_g[seq_len(n_dim_re_tot[1]), seq_len(n_dim_re_tot[1])] <- par_hat$sigma2 * R
      if (n_re > 0) {
        for (j in seq_len(n_re)) {
          sc <- sum(n_dim_re_tot[seq_len(j)])
          diag(Sigma_g[sc + seq_len(n_dim_re_tot[j+1]), sc + seq_len(n_dim_re_tot[j+1])]) <- par_hat$sigma2_re[j]
        }
      }
      prediction_weights <- gaussian_prediction_weights(C_g, Sigma_g, ID_g,
                                                         par_hat$sigma2_me)
      if (!batch_marginal) {
        A <- if (list_mode) {
          lapply(C, prediction_weights)
        } else {
          prediction_weights(C)
        }
      }
    } else {
      Sigma     <- par_hat$sigma2 * R
      Sigma_root <- factor_covariance(Sigma, "observation covariance")
      prediction_weights <- function(cross_covariance) {
        cholesky_prediction_weights(cross_covariance, Sigma_root)
      }
      if (!batch_marginal) {
        A <- if (list_mode) {
          lapply(C, prediction_weights)
        } else {
          prediction_weights(C)
        }
      }
    }

    if (batch_marginal) {
      out$S_samples <- batched_marginal_prediction(
        fitting_grp, conditioning_coords, prediction_weights, diff.y,
        par_hat$sigma2, fitting_phi, object$kappa, n_samples,
        marginal_batch_size
      )
    } else {
      mu_cond_S <- if (list_mode) {
        lapply(A, function(single_grid_A) {
          as.numeric(single_grid_A %*% diff.y)
        })
      } else {
        as.numeric(A %*% diff.y)
      }

      if (type == "marginal") {
        if (list_mode) {
          out$S_samples <- lapply(seq_along(A), function(i) {
            sd_cond_S_i <- sqrt(conditional_variances(
              par_hat$sigma2, A[[i]], C[[i]]
            ))
            sample_independent_gaussian(
              mu_cond_S[[i]], sd_cond_S_i, n_samples
            )
          })
        } else {
          sd_cond_S <- sqrt(conditional_variances(par_hat$sigma2, A, C))
          out$S_samples <- sample_independent_gaussian(
            mu_cond_S, sd_cond_S, n_samples
          )
        }
      } else {
        if (list_mode) {
          out$S_samples <- lapply(seq_along(A), function(i) {
            spatial_covariance_i <- par_hat$sigma2 * matern_correlation(
              pairwise_distances(fitting_grp[[i]]), phi = fitting_phi,
              kappa = object$kappa, return_sym_matrix = TRUE
            )
            conditional_covariance_i <- spatial_covariance_i -
              A[[i]] %*% t(C[[i]])
            cholesky_root_i <- t(factor_covariance(
              conditional_covariance_i,
              "conditional prediction covariance"
            ))
            sample_correlated_gaussian(
              mu_cond_S[[i]], cholesky_root_i, n_samples
            )
          })
        } else {
          Sp <- par_hat$sigma2 * matern_correlation(
            pairwise_distances(fitting_grp), phi = fitting_phi,
            kappa = object$kappa, return_sym_matrix = TRUE
          )
          Sc <- Sp - A %*% t(C)
          Sc_root <- factor_covariance(
            Sc, "conditional prediction covariance"
          )
          Scr <- t(Sc_root)
          out$S_samples <- sample_correlated_gaussian(
            mu_cond_S, Scr, n_samples
          )
        }
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
      C_Z  <- C_g[, -(seq_len(n_dim_re_tot[1]))]
      add  <- 0
      for (i in seq_along(n_dim_re_tot[-1])) {
        C_Z[, add + seq_len(n_dim_re_tot[i+1])] <- par_hat$sigma2_re[i] *
          C_Z[, add + seq_len(n_dim_re_tot[i+1])]
        add <- n_dim_re_tot[i+1]
      }
      W_Z          <- prediction_weights(Matrix::t(C_Z))
      A_Z          <- cholesky_prediction_weights(W_Z %*% t(C), Sc_root)
      Sigma_Z_cond <- diag(rep(par_hat$sigma2_re, n_dim_re_tot[-1])) -
        W_Z %*% C_Z -
        A_Z %*% C %*% t(W_Z)
      Scr_Z        <- t(factor_covariance(
        Sigma_Z_cond, "conditional random-effect covariance"
      ))
      mu_Z_cond    <- A_Z %*% (out$S_samples - mu_cond_S)
      mu_Z_cond <- mu_Z_cond + as.numeric(W_Z %*% diff.y)
      re_samples <- sample_correlated_gaussian(
        mu_Z_cond, Scr_Z, n_samples
      )
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

##' @title Predictive Target Over a Regular Spatial Grid
##'
##' @description Computes predictions over a regular spatial grid using outputs from
##' \code{\link{setup_prediction}}. Custom targets can be supplied via \code{f_target}.
##'
##' @param object Output from \code{\link{setup_prediction}}, a
##'   \code{RiskMap_pred} object.
##' @param include_covariates Logical. Include covariate effects in the linear
##'   predictor. Default \code{TRUE}.
##' @param include_nugget Logical. Add a nugget draw to each spatial sample.
##'   Default \code{FALSE}.
##' @param include_cov_offset Logical. Include the covariate offset.
##'   Default \code{FALSE}.
##' @param include_re Logical. Include unstructured random effects.
##'   Default \code{FALSE}.
##' @param f_target Optional named list of functions to apply on the linear
##'   predictor samples. Each function receives a matrix
##'   (\code{n_pred x n_samples}) and returns a matrix of the same dimensions.
##'   Overrides the model-specific defaults described above.
##' @param pd_summary Optional named list of summary functions applied
##'   row-wise to each target matrix (default: mean, median, sd, 2.5% and
##'   97.5% quantiles).
##'
##' @return An object of class `RiskMap_predict_grid_target` containing:
##' \describe{
##'   \item{target}{List of the predictions for each target}
##'   \item{grid_pred}{`sfc` object containing the coordinates for the predictions}
##'   \item{f_target}{Character vector giving the names of the target functions
##'     that were applied to the linear predictor samples. These names index the
##'     first level of the \code{target} list. If \code{f_target} was supplied
##'     unnamed, the default names \code{"f_target_1"}, \code{"f_target_2"},
##'     ... are used.}
##'   \item{pd_summary}{Character vector giving the names of the summary
##'     functions applied to each target matrix. These names index the second
##'     level of the \code{target} list. If \code{pd_summary} was supplied
##'     unnamed, the default names \code{"pd_summary_1"}, \code{"pd_summary_2"},
##'     ... are used.}
##'   \item{family}{The model family}
##'   \item{lp_samples}{Samples of the linear predictor at the prediction
##'     locations, before the target functions in \code{f_target} are applied.
##'     A matrix with \code{n_pred} rows and \code{n_samples} columns, or a list
##'     of such matrices when \code{object} holds a list of prediction grids.
##'     The contributions of the covariates, the covariate offset, the nugget
##'     and the unstructured random effects are included only when the
##'     corresponding \code{include_*} arguments are set to \code{TRUE}.}
##' }
##' @seealso \code{\link{setup_prediction}}
##' @importFrom Matrix solve
##'
##' @examples
##' data(italy_sim)
##'
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' # using locations in dataset
##' prediction_setup <- setup_prediction(fit)
##'
##' predictions <- predict_grid_target(prediction_setup)
##'
##' @export
predict_grid_target <- function(object,
                             include_covariates  = TRUE,
                             include_nugget      = FALSE,
                             include_cov_offset  = FALSE,
                             include_re          = FALSE,
                             f_target            = NULL,
                             pd_summary          = NULL) {

  if (!inherits(object, "RiskMap_pred"))
    stop("'object' must be an output of setup_prediction()")

  # ---------------------------------------------------------------------------
  # list-mode detection
  # ---------------------------------------------------------------------------
  list_mode <- inherits(object$grid_pred, "list")

  if (list_mode) {
    n_pred <- vapply(object$grid_pred,
                     function(g) nrow(st_coordinates(g)), integer(1))
  } else {
    n_pred <- nrow(object$S_samples)
  }

  if (!is.null(object$par_hat$tau2) == FALSE && include_nugget)
    stop("No nugget was estimated; cannot include it in the predictive target")

  # ---------------------------------------------------------------------------
  # Default f_target — model-specific
  # ---------------------------------------------------------------------------
  if (is.null(f_target)) {

    # glgpm: identity on linear predictor
    f_target <- list(linear_target = function(x) x)
  }

  # ---------------------------------------------------------------------------
  # Default pd_summary
  # ---------------------------------------------------------------------------
  if (is.null(pd_summary)) {
    pd_summary <- list(
      mean   = mean,
      median = median,
      sd     = sd,
      lower  = function(x) quantile(x, 0.025),
      upper  = function(x) quantile(x, 0.975)
    )
  }

  n_f         <- length(f_target)
  n_summaries <- length(pd_summary)
  names_f     <- names(f_target)
  if (is.null(names_f)) names_f <- paste0("f_target_", seq_len(n_f))
  names_s     <- names(pd_summary)
  if (is.null(names_s)) names_s <- paste0("pd_summary_", seq_len(n_summaries))

  if (list_mode) {
    n_samples <- ncol(object$S_samples[[1]])
  } else {
    n_samples <- ncol(object$S_samples)
  }

  n_re <- length(object$re$samples)

  # ---------------------------------------------------------------------------
  # Covariate / offset checks
  # ---------------------------------------------------------------------------
  if (length(object$mu_pred) == 1 && object$mu_pred == 0 && include_covariates)
    stop("Covariates were not provided in setup_prediction(); rerun with 'predictors'")

  if (n_re == 0 && include_re)
    stop("Random effect categories not provided; rerun setup_prediction() with 're_predictors'")

  if (list_mode) {
    mu_target  <- if (include_covariates) object$mu_pred  else lapply(n_pred, function(n) rep(0, n))
    cov_offset <- if (include_cov_offset) object$cov_offset else lapply(n_pred, function(n) rep(0, n))
  } else {
    mu_target  <- if (include_covariates) object$mu_pred  else 0
    cov_offset <- if (include_cov_offset) object$cov_offset else 0
  }

  if (include_cov_offset && length(object$cov_offset) == 1)
    stop("No covariate offset was included in the model")

  # ---------------------------------------------------------------------------
  # Optional nugget draw
  # ---------------------------------------------------------------------------
  if (include_nugget) {
    tau2 <- object$par_hat$tau2
    if (list_mode) {
      object$S_samples <- lapply(seq_along(object$S_samples), function(i) {
        object$S_samples[[i]] +
          matrix(rnorm(n_samples * n_pred[i], sd = sqrt(tau2)), ncol = n_samples)
      })
    } else {
      object$S_samples <- object$S_samples +
        matrix(rnorm(n_samples * n_pred, sd = sqrt(tau2)), ncol = n_samples)
    }
  }

  # ---------------------------------------------------------------------------
  # Build linear predictor samples:  lp = mu + offset + S(x)
  # ---------------------------------------------------------------------------
  if (list_mode) {
    object$S_samples <- lapply(object$S_samples, function(x) {
      if (is.numeric(x) && is.vector(x)) matrix(x, nrow = 1) else x
    })

    lp_samples <- vector("list", length(object$grid_pred))
    for (i in seq_along(object$grid_pred)) {
      lp_i <- vapply(seq_len(n_samples), function(j) {
        mu_i <- if (is.matrix(mu_target[[i]])) mu_target[[i]][, j] else mu_target[[i]]
        mu_i + cov_offset[[i]] + object$S_samples[[i]][, j]
      }, numeric(n_pred[i]))
      lp_samples[[i]] <- .as_sample_matrix(
        lp_i,
        nrow_expected = n_pred[i],
        ncol_expected = n_samples,
        context = sprintf("list-mode linear predictor samples for group %d", i)
      )
    }

  } else {
    ID_coords <- if (object$obs_loc) object$ID_coords else seq_len(n_pred)

    lp_samples <- if (is.matrix(mu_target)) {
      sapply(seq_len(n_samples), function(i)
        mu_target[, i] + cov_offset + object$S_samples[ID_coords, i])
    } else {
      sapply(seq_len(n_samples), function(i)
        mu_target + cov_offset + object$S_samples[ID_coords, i])
    }
    # sapply drops dimensions when n_pred == 1; restore matrix shape
    if (!is.matrix(lp_samples))
      lp_samples <- matrix(lp_samples, nrow = length(ID_coords))
  }

  # ---------------------------------------------------------------------------
  # Unstructured random effects
  # ---------------------------------------------------------------------------
  if (include_re) {
    n_dim_re <- sapply(seq_len(n_re), function(i) length(object$re$samples[[i]]))
    for (i in seq_len(n_re)) {
      for (j in seq_len(n_dim_re[i])) {
        re_samp <- object$re$samples[[i]][[j]]   # length n_samples
        if (list_mode) {
          for (g in seq_along(lp_samples))
            lp_samples[[g]] <- lp_samples[[g]] +
              outer(object$re$D_pred[[i]][, j], re_samp)
        } else {
          lp_samples <- lp_samples +
            outer(object$re$D_pred[[i]][, j], re_samp)
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Apply f_target transformations + summaries
  # ---------------------------------------------------------------------------
  out <- list()
  out$target  <- list()

  if (list_mode) {
    group_names <- names(object$grid_pred) %||%
      paste0("group_", seq_along(object$grid_pred))

    for (i in seq_along(object$grid_pred)) {
      out$target[[group_names[i]]] <- list()

      for (fi in seq_len(n_f)) {
        target_mat <- f_target[[fi]](lp_samples[[i]])
        target_mat <- .as_sample_matrix(
          target_mat,
          nrow_expected = n_pred[i],
          ncol_expected = n_samples,
          context = sprintf("list-mode target matrix for group %d", i)
        )

        out$target[[group_names[i]]][[names_f[fi]]] <- list()
        for (si in seq_len(n_summaries))
          out$target[[group_names[i]]][[names_f[fi]]][[names_s[si]]] <-
          apply(target_mat, 1, pd_summary[[si]])
        }
    }

  } else {
    for (fi in seq_len(n_f)) {
      target_mat <- f_target[[fi]](lp_samples)
      if (!is.matrix(target_mat))
        target_mat <- matrix(target_mat, nrow = nrow(lp_samples))

      out$target[[names_f[fi]]] <- list()
      for (si in seq_len(n_summaries))
        out$target[[names_f[fi]]][[names_s[si]]] <-
        apply(target_mat, 1, pd_summary[[si]])
    }
  }

  # ---------------------------------------------------------------------------
  # Metadata
  # ---------------------------------------------------------------------------
  out$grid_pred  <- object$grid_pred
  out$f_target   <- names_f
  out$pd_summary <- names_s
  out$family     <- object$family
  out$lp_samples <- lp_samples

  class(out) <- "RiskMap_predict_grid_target"
  return(out)
}

##' Plot Method for RiskMap_predict_grid_target Objects
##'
##' Generates a map of the predicted values or summaries over the regular spatial grid
##' from an object of class 'RiskMap_predict_grid_target'.
##'
##' @param x An object of class 'RiskMap_predict_grid_target'.
##' @param target Character string specifying which target prediction to plot,
##' one of \code{x$f_target}. If \code{NULL} (the default), the first target is used.
##' @param summary Character string specifying which summary statistic to plot
##' (e.g., "mean", "sd"), one of \code{x$pd_summary}. Defaults to \code{"mean"}.
##' @param palette Either the name of a palette from [grDevices::hcl.pals()]
##' (e.g. `"viridis"`, `"Blues"`, `"Spectral"`) or a character vector of two or
##' more colours to interpolate between. Defaults to `"viridis"`.
##' @param reverse_palette Logical; if `TRUE`, reverse the order of the palette
##' colours. Defaults to `FALSE`.
##' @param north_arrow Logical; if `TRUE` (the default), add a north arrow to the map.
##' @param scale_bar Logical; if `TRUE` (the default), add a scale bar to the map.
##' @param limits Range of the colour scale: `"shared"` (the default),
##' `"independent"` or a numeric vector of length two giving custom limits.
##' For grid maps, `"shared"` and `"independent"` both use the range of the
##' plotted values; this is the range that
##' [plot.RiskMap_predict_areal_target()] shares by default.
##' @param ... Additional arguments passed to [ggplot2::scale_fill_gradientn()],
##' e.g. `name` or `breaks`.
##' @return A \code{ggplot} object representing the specified prediction target or summary statistic over the spatial grid.
##' @details
##' The grid cells are drawn with [ggplot2::geom_tile()] in the coordinate
##' reference system of the prediction grid, with axes labelled in longitude
##' and latitude. Tiles are used rather than a raster so that grids with
##' missing rows or columns, or that are not exactly regular (e.g. after
##' reprojection), are drawn at their true locations. When the north arrow or
##' scale bar is shown, space is added above or below the data so that they do
##' not cover it. The styling matches [plot.RiskMap_predict_areal_target()].
##' The returned object can be modified further with standard \pkg{ggplot2}
##' functions.
##'
##' @seealso \code{\link{predict_grid_target}}, \code{\link{plot.RiskMap_predict_areal_target}}
##'
##' @method plot RiskMap_predict_grid_target
##' @export
##'
##'
plot.RiskMap_predict_grid_target <- function(x, target = NULL, summary = "mean",
                                             palette = "viridis",
                                             reverse_palette = FALSE,
                                             north_arrow = TRUE,
                                             scale_bar = TRUE,
                                             limits = "shared", ...) {
  stopifnot("'x' must be of class RiskMap_predict_grid_target" =
              inherits(x, "RiskMap_predict_grid_target"))
  if (inherits(x$grid_pred, "list")) {
    stop("Plotting is not supported when 'grid_pred' is a list of grids")
  }
  target <- check_plot_target(x, target, summary)

  grid_coordinates <- st_coordinates(x$grid_pred)
  plot_data <- data.frame(x = grid_coordinates[, 1],
                          y = grid_coordinates[, 2],
                          value = x$target[[target]][[summary]])

  # extent of the tiles, not just their centres
  half_width <- resolution(plot_data$x, zero = FALSE) / 2
  half_height <- resolution(plot_data$y, zero = FALSE) / 2
  extent <- c(xmin = min(plot_data$x) - half_width, xmax = max(plot_data$x) + half_width,
              ymin = min(plot_data$y) - half_height, ymax = max(plot_data$y) + half_height)

  out <- ggplot(plot_data) +
    geom_tile(aes(x = .data$x, y = .data$y, fill = .data$value))

  add_map_layers(out,
                 extent = extent,
                 crs = st_crs(x$grid_pred),
                 legend_title = paste0(target, "_", summary),
                 palette = palette,
                 reverse_palette = reverse_palette,
                 limits = resolve_colour_limits(limits),
                 north_arrow = north_arrow,
                 scale_bar = scale_bar,
                 ...)
}

##' Validate the target and summary requested from a prediction object
##'
##' @param x An object of class 'RiskMap_predict_grid_target' or
##' 'RiskMap_predict_areal_target'.
##' @param target Character string naming the target, or `NULL` for the first.
##' @param summary Character string naming the summary statistic.
##' @return The name of the target to plot.
##' @noRd
check_plot_target <- function(x, target, summary) {
  if (is.null(target)) {
    target <- x$f_target[1]
  } else if (!target %in% x$f_target) {
    stop("'target' must be one of: ", paste(shQuote(x$f_target), collapse = ", "))
  }
  if (!summary %in% x$pd_summary) {
    stop("'summary' must be one of: ", paste(shQuote(x$pd_summary), collapse = ", "))
  }
  target
}

##' Convert a palette specification into a vector of colours
##'
##' @param palette Either the name of a palette in [grDevices::hcl.pals()] or a
##' character vector of two or more valid colours.
##' @param reverse_palette Logical; reverse the order of the colours.
##' @return A character vector of colours.
##' @noRd
palette_colours <- function(palette, reverse_palette = FALSE) {
  stopifnot("'palette' must be a character vector" = is.character(palette))
  check_logical(reverse_palette)

  if (length(palette) == 1) {
    # hcl.colors() ignores case, spaces and punctuation in palette names
    normalise_name <- function(name) tolower(gsub("[^[:alnum:]]", "", name))
    if (!normalise_name(palette) %in% normalise_name(grDevices::hcl.pals())) {
      stop("'palette' must be one of grDevices::hcl.pals() or a vector of two or more colours")
    }
    colours <- grDevices::hcl.colors(256, palette = palette)
  } else {
    valid_colour <- vapply(palette, function(colour) {
      !inherits(try(grDevices::col2rgb(colour), silent = TRUE), "try-error")
    }, logical(1))
    if (!all(valid_colour)) {
      stop("'palette' contains invalid colours: ",
           paste(shQuote(palette[!valid_colour]), collapse = ", "))
    }
    colours <- palette
  }
  if (reverse_palette) colours <- rev(colours)
  colours
}

##' Add the fill scale, map extent, north arrow and scale bar shared by the
##' prediction maps
##'
##' @param map A `ggplot` object.
##' @param extent Named numeric vector (`xmin`, `xmax`, `ymin`, `ymax`) giving
##' the extent of the data in `crs`.
##' @param crs The coordinate reference system of the data.
##' @param legend_title Character string used as the fill legend title.
##' @param palette,reverse_palette Passed to `palette_colours()`.
##' @param limits `NULL` or a numeric vector of length two for the colour scale.
##' @param north_arrow,scale_bar Logical; whether to add each annotation.
##' @param ... Additional arguments passed to [ggplot2::scale_fill_gradientn()].
##' @return The `ggplot` object with the layers added.
##' @noRd
add_map_layers <- function(map, extent, crs, legend_title, palette,
                           reverse_palette, limits, north_arrow, scale_bar, ...) {
  check_logical(north_arrow)
  check_logical(scale_bar)
  map <- map +
    scale_fill_gradientn(colours = palette_colours(palette, reverse_palette),
                         limits = limits, ...) +
    coord_sf(crs = crs,
             xlim = as.numeric(extent[c("xmin", "xmax")]),
             ylim = pad_map_extent(extent, north_arrow, scale_bar)) +
    labs(x = NULL, y = NULL, fill = legend_title)

  if (north_arrow) {
    map <- map +
      ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                        height = unit(1, "cm"),
                                        width = unit(1, "cm"))
  }
  if (scale_bar) {
    map <- map + ggspatial::annotation_scale(location = "bl")
  }
  map
}

##' Resolve the colour scale limits requested for a prediction map
##'
##' @param limits `"shared"`, `"independent"` or a numeric vector of length two.
##' @param shared_range Numeric values whose range defines the `"shared"`
##' limits, or `NULL` to use the range of the plotted values.
##' @return `NULL` (use the range of the plotted values) or a numeric vector of
##' length two.
##' @noRd
resolve_colour_limits <- function(limits, shared_range = NULL) {
  if (is.numeric(limits)) {
    if (length(limits) != 2 || anyNA(limits)) {
      stop("Custom 'limits' must be a numeric vector of length two")
    }
    return(limits)
  }
  if (!(is.character(limits) && length(limits) == 1 &&
        limits %in% c("shared", "independent"))) {
    stop("'limits' must be \"shared\", \"independent\" or a numeric vector of length two")
  }
  if (limits == "shared" && !is.null(shared_range)) {
    return(range(shared_range, na.rm = TRUE))
  }
  NULL
}

##' Extend the vertical extent of a map to make room for annotations
##'
##' @param extent Named numeric vector (`xmin`, `xmax`, `ymin`, `ymax`).
##' @param north_arrow Logical; add space above the data for a north arrow.
##' @param scale_bar Logical; add space below the data for a scale bar.
##' @return Numeric vector of length two giving the padded y limits.
##' @noRd
pad_map_extent <- function(extent, north_arrow, scale_bar) {
  height <- unname(extent["ymax"] - extent["ymin"])
  # roughly the height of each annotation on a typical map
  top_padding <- if (north_arrow) 0.15 * height else 0
  bottom_padding <- if (scale_bar) 0.08 * height else 0
  unname(c(extent["ymin"] - bottom_padding, extent["ymax"] + top_padding))
}

##' @title Predictive Targets over Boundaries (grid-aggregated)
##'
##' @description
##' Computes predictive targets over polygon features using joint prediction
##' samples from \code{\link{setup_prediction}}. Targets can incorporate
##' covariates, offsets, optional unstructured random effects.
##'
##' @param object Output from \code{\link{setup_prediction}} (class \code{RiskMap_pred}),
##'   typically fitted with \code{type = "joint"} so that linear predictor samples are available.
##' @param boundaries An \pkg{sf} polygon object representing regions over which predictions are aggregated.
##' @param areal_target A function that aggregates grid-cell values within each polygon to a
##'   single regional value (default \code{mean}). Examples: \code{mean}, \code{sum},
##'   a custom weighted mean, etc.
##' @param weights Optional numeric vector of weights used inside \code{areal_target}.
##'   If supplied with \code{standardize_weights = TRUE}, weights are normalized within each region.
##' @param standardize_weights Logical; standardize \code{weights} within each region (\code{FALSE} by default).
##' @param col_names Name or column index in \code{boundaries} containing region identifiers to use in outputs.
##' @param include_covariates Logical; include fitted covariate effects in the linear predictor (default \code{TRUE}).
##' @param include_nugget Logical; include the nugget (unstructured measurement error) in the linear predictor (default \code{FALSE}).
##' @param include_cov_offset Logical; include any covariate offset term (default \code{FALSE}).
##' @param return_boundaries Logical; if \code{TRUE}, return \code{boundaries} with appended summary columns
##'   defined by \code{pd_summary} (default \code{TRUE}).
##' @param include_re Logical; include unstructured random effects (RE) in the linear predictor (default \code{FALSE}).
##' @param f_target List of target functions applied to linear predictor samples (e.g.,
##'   \code{list(prev = plogis)} for prevalence on the probability scale). If \code{NULL},
##'   the identity is used.
##' @param pd_summary Named list of summary functions applied to each region's target samples
##'   (e.g., \code{list(mean = mean, sd = sd, q025 = function(x) quantile(x, 0.025), q975 = function(x) quantile(x, 0.975))}).
##'   Names are used as column suffixes in the outputs.
##' @param messages Logical; if \code{TRUE}, print progress messages while computing regional targets.
##' @param return_target_samples Logical; if \code{TRUE}, also return the raw target samples per region
##'   (default \code{FALSE}).
##'
##' @details
##' For each polygon in \code{boundaries}, grid-cell samples of the linear predictor are transformed with
##' \code{f_target}, optionally adjusted for covariates, offset, nugget and/or REs, and
##' then aggregated via \code{areal_target} (optionally weighted). The list \code{pd_summary} is applied
##' to each region's target samples to produce summary statistics.
##'
##' @return An object of class \code{RiskMap_predict_areal_target} with components:
##' \itemize{
##'   \item \code{target}: \code{data.frame} of region-level summaries (one row per region).
##'   \item \code{target_samples}: (optional) \code{list} with one element per region; each contains
##'         a \code{data.frame}/matrix of raw samples for each named target in \code{f_target},
##'         if \code{return_target_samples = TRUE}.
##'   \item \code{boundaries}: (optional) the input \code{sf} object with appended summary columns,
##'         included if \code{return_boundaries = TRUE}.
##'   \item \code{grid_summary_range}: nested \code{list} giving, for each target and
##'         summary, the range of the same summary computed for each grid cell. Used by
##'         [plot.RiskMap_predict_areal_target()] so that its colour scale matches
##'         [plot.RiskMap_predict_grid_target()].
##'   \item \code{f_target}, \code{pd_summary}, \code{grid_pred}: inputs echoed for reproducibility.
##' }
##'
##' @seealso \code{\link{setup_prediction}}, \code{\link{predict_grid_target}}
##'
##' @importFrom terra rast as.data.frame
##'
##' @examples
##' library(sf)
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
##' # joint predictions are required
##' prediction_setup <- setup_prediction(fit, type = "joint")
##'
##' # split the whole area into two
##' areal <- st_sf(geometry = st_make_grid(italy_subset, n = c(1, 2)))
##'
##' predictions <- predict_areal_target(prediction_setup, areal)
##'
##' @export
predict_areal_target <- function(object,
                            boundaries,
                            areal_target = mean,
                            weights = NULL,
                            standardize_weights = FALSE,
                            col_names = NULL,
                            include_covariates = TRUE,
                            include_nugget = FALSE,
                            include_cov_offset = FALSE,
                            return_boundaries = TRUE,
                            include_re = FALSE,
                            f_target = NULL,
                            pd_summary = NULL,
                            messages = TRUE,
                            return_target_samples = FALSE) {

  if(!inherits(object, what = "RiskMap_pred", which = FALSE)) {
    stop("The object passed to 'object' must be an output of
         the function 'setup_prediction'")
  }

  if(object$type != "joint") {
    stop("To run predictions for areas, joint predictions must be used;
         rerun 'setup_prediction' and set 'type' = \"joint\"")
  }

  check_data(boundaries, "polygon")

  list_mode <- inherits(object$grid_pred, "list")

  if (list_mode) {

    n_pred <- vapply(object$grid_pred, function(g) nrow(st_coordinates(g)), integer(1))

    if (!is.null(weights)) {
      if (!is.list(weights)) {
        stop("When 'grid_pred' passed to 'setup_prediction' is a list, 'weights' must also be a list,
             with one numeric vector per element of 'object$grid_pred'")
      }
      if (length(weights) != length(object$grid_pred)) {
        stop("Length of 'weights' must match length of 'grid_pred' passed to 'setup_prediction'")
      }
      for (i in seq_along(weights)) {
        if (!is.numeric(weights[[i]])) {
          stop(sprintf("'weights[[%d]]' must be numeric", i))
        }
        if (length(weights[[i]]) != n_pred[i]) {
          stop(sprintf("Length of 'weights[[%d]]' (%d) must equal number of locations in 'grid_pred[[%d]]' passed to 'setup_prediction' (%d).",
                       i, length(weights[[i]]), i, n_pred[i]))
        }
        if (anyNA(weights[[i]])) {
          warning(sprintf("NA values found in 'weights[[%d]]'; they will be treated as 0.", i))
          weights[[i]][is.na(weights[[i]])] <- 0
        }
      }
      no_weights <- FALSE
    } else {
      weights <- lapply(n_pred, function(n) rep(1, n))
      no_weights <- TRUE
    }

  } else {
    n_pred <- nrow(st_coordinates(object$grid_pred))
  }

  n_re <- length(object$re$samples)
  re_names <- names(object$re$samples)

  if(n_re == 0 && include_re) {
    stop("The categories of the random effects variables have not been provided;
         re-run 'setup_prediction' and provide the covariates through the argument 're_predictors'")
  }

  if(!is.null(weights)) {
    no_weights <- FALSE
  } else {
    if(list_mode) {
      weights <- sapply(unlist(n_pred), function(i) rep(1, i))
    } else {
      weights <- rep(1, n_pred)
    }
    no_weights <- TRUE
  }

  if(is.null(object$par_hat$tau2) & include_nugget) {
    stop("The nugget cannot be included in the predictive target
             because it was not included when fitting the model")
  }

  if(is.null(f_target)) {
    f_target <- list(linear_target = function(x) x)
  }

  if(is.null(pd_summary)) {
    pd_summary <- list(mean = mean, sd = sd)
  }

  n_summaries <- length(pd_summary)
  n_f <- length(f_target)
  if(list_mode) {
    n_samples <- ncol(object$S_samples[[1]])
  } else {
    n_samples <- ncol(object$S_samples)
  }

  if(!list_mode && (length(object$mu_pred) == 1 && object$mu_pred == 0 && include_covariates)) {
    stop("Covariates have not been provided; re-run setup_prediction
         and provide the covariates through the argument 'predictors'")
  }

  if(!include_covariates) {
    mu_target <- 0
  } else {
    if(is.null(object$mu_pred)) stop("the output obtained from 'setup_prediction' does not
                                     contain any covariates; if including covariates
                                     in the predictive target these should be included
                                     when running 'setup_prediction'")
    mu_target <- object$mu_pred
  }

  if(!include_cov_offset) {
    if(list_mode) {
      cov_offset <- lapply(n_pred, function(i) rep(0, i))
    } else {
      cov_offset <- 0
    }
  } else {
    if(length(object$cov_offset) == 1) {
      stop("No covariate offset was included in the model;
           set include_cov_offset = FALSE, or refit the model and include
           the covariate offset")
    }
    cov_offset <- object$cov_offset
  }

  if(include_nugget) {
    if (list_mode) {
      object$S_samples <- lapply(seq_along(object$S_samples), function(i) {
        object$S_samples[[i]] +
          matrix(rnorm(n_samples * n_pred[i], sd = sqrt(object$par_hat$tau2)),
                 nrow = n_pred[i], ncol = n_samples)
      })
    } else {
      Z_sim <- matrix(rnorm(n_samples * n_pred, sd = sqrt(object$par_hat$tau2)),
                      ncol = n_samples)
      object$S_samples <- object$S_samples + Z_sim
    }
  }

  out <- list()
  if (return_target_samples) out$target_samples <- list()

  if(list_mode) {
    object$S_samples <- lapply(object$S_samples, function(x) {
      if (is.numeric(x) && is.vector(x)) {
        matrix(x, nrow = 1)
      } else {
        x
      }
    })

    if(is.matrix(mu_target[[1]])) {
      out$lp_samples <- lapply(seq_along(object$grid_pred), function(j) {
        lp_j <- vapply(seq_len(n_samples),
                       function(i)
                         mu_target[[j]][, i] + cov_offset[[j]] +
                         object$S_samples[[j]][, i],
                       numeric(n_pred[j]))
        .as_sample_matrix(
          lp_j,
          nrow_expected = n_pred[j],
          ncol_expected = n_samples,
          context = sprintf("list-mode linear predictor samples for group %d", j)
        )
      })
    } else {
      out$lp_samples <- lapply(seq_along(object$grid_pred), function(j) {
        lp_j <- vapply(seq_len(n_samples),
                       function(i)
                         mu_target[[j]] + cov_offset[[j]] +
                         object$S_samples[[j]][, i],
                       numeric(n_pred[j]))
        .as_sample_matrix(
          lp_j,
          nrow_expected = n_pred[j],
          ncol_expected = n_samples,
          context = sprintf("list-mode linear predictor samples for group %d", j)
        )
      })
    }
  } else {
    if(is.matrix(mu_target)) {
      out$lp_samples <- sapply(1:n_samples, function(i)
        mu_target[, i] + cov_offset + object$S_samples[, i])
    } else {
      out$lp_samples <- sapply(1:n_samples, function(i)
        mu_target + cov_offset + object$S_samples[, i])
    }
  }

  if(include_re) {
    n_dim_re <- sapply(1:n_re, function(i) length(object$re$samples[[i]]))
    for(i in 1:n_re) {
      for(j in 1:n_dim_re[i]) {
        if (list_mode) {
          for(g in seq_along(out$lp_samples)) {
            D_re_g <- object$re$D_pred[[i]]
            if (is.list(D_re_g)) D_re_g <- D_re_g[[g]]
            out$lp_samples[[g]] <- out$lp_samples[[g]] +
              outer(D_re_g[, j], object$re$samples[[i]][[j]])
          }
        } else {
          for(h in 1:n_samples) {
            out$lp_samples[, h] <- out$lp_samples[, h] +
              object$re$D_pred[[i]][, j] * object$re$samples[[i]][[j]][h]
          }
        }
      }
    }
  }

  names_f <- names(f_target)
  names_s <- names(pd_summary)
  out$target <- list()

  n_reg <- nrow(boundaries)
  if(is.null(col_names)) {
    boundaries$region <- paste("reg", 1:n_reg, sep = "")
    col_names <- "region"
    names_reg <- boundaries$region
  } else {
    names_reg <- boundaries[[col_names]]
    if(n_reg != length(names_reg)) {
      stop("The names in the column identified by 'col_names' do not
         provide a unique set of names, but there are duplicates")
    }
  }

  if(list_mode) {
    boundaries <- st_transform(boundaries, crs = st_crs(object$grid_pred[[1]]))
  } else {
    boundaries <- st_transform(boundaries, crs = st_crs(object$grid_pred))
  }

  if(!list_mode) {
    inter <- st_intersects(boundaries, object$grid_pred)
    if(any(is.na(weights))) {
      warning("Missing values found in 'weights' are set to 0 \n")
      weights[is.na(weights)] <- 0
    }
  } else {
    for(i in 1:length(object$grid_pred)) {
      if(any(is.na(weights[[i]]))) {
        warning("Missing values found in 'weights' are set to 0 \n")
        weights[[i]][is.na(weights[[i]])] <- 0
      }
    }
  }

  no_comp <- NULL
  for(h in 1:n_reg) {

    if(list_mode) {
      if(messages) message("Computing predictive target for: ", boundaries[[col_names]][h])
      if(standardize_weights & !no_weights) {
        weights_h <- weights[[h]] / sum(weights[[h]])
      } else {
        weights_h <- weights[[h]]
      }
      for(i in 1:n_f) {
        target_grid_samples_i <- .as_sample_matrix(
          f_target[[i]](out$lp_samples[[h]]),
          nrow_expected = n_pred[h],
          ncol_expected = n_samples,
          context = sprintf(
            "list-mode target aggregation for group %d",
            h
          )
        )
        if (nrow(target_grid_samples_i) != length(weights_h)) {
          stop(sprintf(
            "Dimension mismatch in list-mode target aggregation for group %d: nrow(target_grid_samples_i) = %d, length(weights_h) = %d",
            h, nrow(target_grid_samples_i), length(weights_h)
          ))
        }

        target_samples_i <- apply(target_grid_samples_i, 2, function(x) areal_target(weights_h * x))

        if (return_target_samples) {
          if (is.null(out$target_samples[[ names_reg[h] ]])) out$target_samples[[ names_reg[h] ]] <- list()
          out$target_samples[[ names_reg[h] ]][[ names_f[i] ]] <- target_samples_i
        }

        out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]] <- list()
        for(j in 1:n_summaries) {
          out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]][[ paste(names_s[j]) ]] <-
            pd_summary[[j]](target_samples_i)
        }
      }
      if(messages) message(" \n")

    } else {
      if(messages) message("Computing predictive target for:", boundaries[[col_names]][h])
      if(length(inter[[h]]) == 0) {
        warning(paste("No points on the grid fall within", boundaries[[col_names]][h],
                      "and no predictions are carried out for this area"))
        no_comp <- c(no_comp, h)
      } else {
        ind_grid_h <- inter[[h]]
        if(standardize_weights & !no_weights) {
          weights_h <- weights[ind_grid_h] / sum(weights[ind_grid_h])
        } else {
          weights_h <- weights[ind_grid_h]
        }
        for(i in 1:n_f) {
          target_grid_samples_i <- as.matrix(f_target[[i]](out$lp_samples[ind_grid_h, ]))

          target_samples_i <- apply(target_grid_samples_i, 2, function(x) areal_target(weights_h * x))

          if (return_target_samples) {
            if (is.null(out$target_samples[[ names_reg[h] ]])) out$target_samples[[ names_reg[h] ]] <- list()
            out$target_samples[[ names_reg[h] ]][[ names_f[i] ]] <- target_samples_i
          }

          out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]] <- list()
          for(j in 1:n_summaries) {
            out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]][[ paste(names_s[j]) ]] <-
              pd_summary[[j]](target_samples_i)
          }
        }
      }
      if(messages) message(" \n")
    }
  }

  if(return_boundaries) {
    if(length(no_comp) > 0) {
      ind_reg <- (1:n_reg)[-no_comp]
    } else {
      ind_reg <- 1:n_reg
    }
    for(i in 1:n_f) {
      for(j in 1:n_summaries) {
        name_ij <- paste(names_f[i], "_", paste(names_s[j]), sep = "")
        boundaries[[name_ij]] <- rep(NA, n_reg)
        for(h in ind_reg) {
          which_reg <- which(boundaries[[col_names]] == names_reg[h])
          boundaries[which_reg, ][[name_ij]] <-
            out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]][[ paste(names_s[j]) ]]
        }
      }
    }
    out$boundaries <- boundaries
  }

  out$f_target <- names(f_target)
  out$pd_summary <- names(pd_summary)
  out$grid_pred <- object$grid_pred
  out$grid_summary_range <- summarise_grid_range(out$lp_samples, f_target, pd_summary)
  class(out) <- "RiskMap_predict_areal_target"
  return(out)
}


##' Range of grid-cell summaries for each predictive target
##'
##' Computes each summary in `pd_summary` for every grid cell and returns its
##' range, so that areal maps can share a colour scale with grid maps.
##'
##' @param lp_samples Matrix of linear predictor samples (cells by samples), or
##' a list of such matrices.
##' @param f_target Named list of target functions.
##' @param pd_summary Named list of summary functions.
##' @return A nested list indexed by target then summary, each a numeric
##' vector of length two.
##' @noRd
summarise_grid_range <- function(lp_samples, f_target, pd_summary) {
  if (is.list(lp_samples)) lp_samples <- do.call(rbind, lp_samples)
  lapply(f_target, function(target_function) {
    target_samples <- as.matrix(target_function(lp_samples))
    lapply(pd_summary, function(summary_function) {
      cell_summaries <- apply(target_samples, 1, summary_function)
      range(as.numeric(cell_summaries), na.rm = TRUE)
    })
  })
}

##' Plot Method for RiskMap_predict_areal_target Objects
##'
##' Generates a map of predictive target values or summaries over boundaries.
##'
##' @param x An object of class 'RiskMap_predict_areal_target' containing computed targets,
##' summaries, and associated spatial data.
##' @param target Character indicating the target type to plot (e.g., "linear_target"),
##' one of \code{x$f_target}. If \code{NULL} (the default), the first target is used.
##' @param summary Character indicating the summary type to plot (e.g., "mean", "sd"),
##' one of \code{x$pd_summary}. Defaults to \code{"mean"}.
##' @inheritParams plot.RiskMap_predict_grid_target
##' @return A \code{ggplot} object showing the plot of the specified predictive target or summary.
##' @param limits Range of the colour scale. One of:
##' - `"shared"` (the default): covers both the areal values and the same
##'   summary computed for each grid cell, so that colours match those of
##'   [plot.RiskMap_predict_grid_target()] for the same target and summary.
##' - `"independent"`: the range of the areal values only.
##' - a numeric vector of length two giving custom limits.
##' @details
##' The styling matches [plot.RiskMap_predict_grid_target()]. The returned
##' object can be modified further with standard \pkg{ggplot2} functions.
##' @seealso
##' \code{\link{predict_areal_target}}, \code{\link{plot.RiskMap_predict_grid_target}},
##' \code{\link[ggplot2]{geom_sf}}, \code{\link[ggplot2]{scale_fill_gradientn}}
##'
##' @method plot RiskMap_predict_areal_target
##' @export
plot.RiskMap_predict_areal_target <- function(x, target = NULL, summary = "mean",
                                              palette = "viridis",
                                              reverse_palette = FALSE,
                                              north_arrow = TRUE,
                                              scale_bar = TRUE,
                                              limits = "shared", ...) {
  stopifnot("'x' must be of class RiskMap_predict_areal_target" =
              inherits(x, "RiskMap_predict_areal_target"))
  if (is.null(x$boundaries)) {
    stop("'x' does not contain boundaries; rerun predict_areal_target() with 'return_boundaries = TRUE'")
  }
  target <- check_plot_target(x, target, summary)

  col_boundaries_name <- paste(target, "_", summary, sep = "")
  # "shared" covers the grid-level map of the same quantity as well
  limits <- resolve_colour_limits(
    limits,
    shared_range = c(x$boundaries[[col_boundaries_name]],
                     x$grid_summary_range[[target]][[summary]])
  )

  out <- ggplot(x$boundaries) +
    geom_sf(aes(fill = .data[[col_boundaries_name]]))

  add_map_layers(out,
                 extent = st_bbox(x$boundaries),
                 crs = st_crs(x$boundaries),
                 legend_title = col_boundaries_name,
                 palette = palette,
                 reverse_palette = reverse_palette,
                 limits = limits,
                 north_arrow = north_arrow,
                 scale_bar = scale_bar,
                 ...)
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



.anpit_area <- function(curve, u = seq(0, 1, length.out = length(curve))) {
  if (length(curve) != length(u)) stop("'curve' and 'u' must have the same length")
  if (length(curve) < 2) return(NA_real_)

  dx <- diff(u)
  y <- abs(curve - u)
  sum(dx * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}

.as_sample_matrix <- function(x, nrow_expected, ncol_expected, context) {
  if (!is.matrix(x)) {
    if (length(x) != nrow_expected * ncol_expected) {
      stop(sprintf(
        "%s: expected %d values, got %d",
        context, nrow_expected * ncol_expected, length(x)
      ))
    }
    x <- matrix(x, nrow = nrow_expected, ncol = ncol_expected)
  }

  if (nrow(x) != nrow_expected || ncol(x) != ncol_expected) {
    stop(sprintf(
      "%s: expected a %d x %d matrix, got %d x %d",
      context, nrow_expected, ncol_expected, nrow(x), ncol(x)
    ))
  }

  x
}

##' @title Plot Training/Test Splits
##'
##' @description
##' Plots every observation used by \code{\link{assess_prediction}}, coloured
##' by whether it fell in the training or test set for a given iteration,
##' faceted by iteration when there is more than one. Mimics the visual style
##' of \code{spatialsample::autoplot()} (used directly for
##' \code{method = "cluster"}), for the \code{"user"} and \code{"regularized"}
##' methods, whose splits aren't necessarily a complete partition of the data
##' and so aren't well suited to colouring by fold membership alone.
##'
##' @param data_split A list with a \code{splits} element, each entry itself a
##' list with an \code{sf} \code{data} (the training set) and an \code{sf}
##' \code{data_test} (the test set) - the structure built internally by
##' \code{assess_prediction()} for \code{method = "user"} and
##' \code{method = "regularized"}.
##' @param alpha Point transparency, passed to \code{ggplot2::geom_sf()}.
##' Defaults to \code{0.6}, matching \code{spatialsample::autoplot()}.
##' @return A \code{ggplot} object.
##' @noRd
plot_folds <- function(data_split, alpha = 0.6) {
  n_iter <- length(data_split$splits)

  combined <- do.call(rbind, lapply(seq_len(n_iter), function(i) {
    split_i <- data_split$splits[[i]]
    rbind(
      cbind(split_i$data,      set = "Training", iteration = i),
      cbind(split_i$data_test, set = "Testing",  iteration = i)
    )
  }))

  p <- ggplot(data = combined, aes(color = .data$set, fill = .data$set)) +
    geom_sf(alpha = alpha) +
    guides(colour = guide_legend("Set"), fill = guide_legend("Set")) +
    theme_minimal()

  if (n_iter > 1) {
    p <- p + facet_wrap(vars(.data$iteration))
  }

  p + coord_sf()
}

##' @title Assess Predictive Performance via Spatial Cross-Validation
##'
##' @description
##' This function evaluates the predictive performance of spatial models fitted
##' to `RiskMap` objects using cross-validation. It supports two classes of diagnostic tools:
##'
##' - **Scoring rules**, including the Continuous Ranked Probability Score (CRPS)
##'  and its scaled version (SCRPS), which quantify the sharpness and calibration
##'  of probabilistic forecasts (Bolin & Wallin, 2023);
##' - **Calibration diagnostics**, based on the Probability Integral Transform (PIT)
##' for Gaussian outcomes, Average nonrandomized PIT (AnPIT) curves for discrete
##' outcomes (e.g., Poisson or Binomial), and the area between the PIT/AnPIT curve
##' and the reference line (Giorgi *et al.* 2026).
##'
##' Cross-validation can be performed using either spatial clustering, regularized
##' subsampling with a minimum inter-point distance or a user-defined test set.
##' For each fold or subset, models can be refitted or evaluated with fixed parameters,
##' offering flexibility in model validation. The function also provides visualizations
##' of the spatial distribution of test folds.
##'
##' @param object A list of `RiskMap` objects, each representing a model fitted with `glgpm()`.
##' @param method Character; either `"cluster"`, `"regularized"` or `"user"` for the
##' cross-validation method:
##' \describe{
##'  - The `"cluster"` method uses spatial clustering as implemented
##' by the \code{spatial_clustering_cv} function from the `spatialsample` package.
##'  - The `"regularized"` method selects a subsample of the dataset by imposing a minimum distance,
##'  set by the `min_dist` argument, for a randomly selected subset of locations using
##'  the `subsample.distance` function from the `spatialEco` package.
##'  - The `"user"` method takes a user-defined test set defined by `user_split`
##'  }
##' @param fold Integer; required when `method = "cluster"` - number of folds for cross-validation.
##' @param min_dist Numeric; required when `method = "regularized"` - minimum distance in kilometers for regularized subsampling.
##' @param size Integer; the size of the test set, required when `method = "regularized"`.
##' @param user_split Required when `method = "user"`. A user-defined cross-validation split. Either:
##'   * a matrix with \code{nrow = n} (number of observations) and
##'     \code{ncol = iter} (number of iterations), where entries of \code{1}
##'     indicate membership in the test set for that iteration and \code{0}
##'     indicate training set; or
##'   * a list of length \code{iter}, where each element is either a vector of
##'     indices of the dataset to use as the test set, or a list with components
##'     \code{in_id} (training indices) and \code{out_id} (test indices).
##' @param iter Integer; number of times to repeat the cross-validation. Defaults to `1`.
##' @param metrics Character vector; one or more of `"CRPS"`, `"SCRPS"` and `"AnPIT"`,
##' to specify the predictive performance metrics to compute. When `"AnPIT"` is requested,
##' the scalar score `"AnPIT_area"` is also computed as the integrated absolute deviation
##' between the PIT/AnPIT curve and the reference line. Defaults to all metrics.
##' @param keep_par_fixed Logical; whether to keep parameters fixed across folds,
##' or re-estimate for each fold. Defaults to `TRUE`.
##' @param control_mcmc Control settings for simulation, an output from `set_control_mcmc()`.
##' If it carries a `seed`, the random splits generated for `method = "cluster"`
##' or `method = "regularized"` are reproducible; ignored for `method = "user"`,
##' which is already deterministic. The caller's random number generator state
##' is restored on exit.
##' @param plot_fold Logical; whether to plot each iteration's test and training sets. Defaults to `TRUE`.
##' @param messages Logical; whether to display progress messages. Defaults to `TRUE`.
##' @param ... Additional arguments passed to clustering or subsampling functions.
##'
##' @return A list of class `RiskMap_cross_validation`, containing:
##' \describe{
##'   \item{test_set}{A list of test sets used for validation, each of class `'sf'`.}
##'   \item{model}{A named list, one per model, each containing:
##'     \describe{
##'       \item{metric}{A list with CRPS, SCRPS, and/or AnPIT_area metrics for each fold if requested.}
##'       \item{PIT}{(if `family = "gaussian"` and `metrics` includes `"AnPIT"`) A list of PIT values for test data.}
##'       \item{AnPIT}{(if `family` is discrete and `metrics` includes `"AnPIT"`) A list of AnPIT curves for test data.}
##'     }
##'   }
##' }
##'
##' @seealso \code{\link{plot.RiskMap_cross_validation}}
##'
##' @references
##' Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of
##' proper scoring rules. *Statistical Science*, 38(1), 140–159. \doi{10.1214/22-STS864}.
##'
##' Giorgi, E., Fronterre, C. & Diggle, P. J. (2026). A decay-adjusted spatio-temporal
##' model to account for the impact of mass drug administration on neglected
##' tropical disease prevalence. *Journal of the Royal Statistical Society Series
##' A: Statistics in Society*.\doi{10.1093/jrsssa/qnag100}.
##'
##' @importFrom terra match
##' @importFrom spatialEco subsample.distance
##' @importFrom spatialsample spatial_clustering_cv autoplot
##'
##' @examples
##'
##' data(italy_sim)
##'
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' # cluster method
##' cross_validation <-
##'   assess_prediction(
##'     list(fit),
##'     method = "cluster",
##'     fold = 2
##'   )
##'
##' summary(cross_validation)
##'
##' # regularized method
##' cross_validation <-
##'   assess_prediction(
##'     list(fit),
##'     method = "regularized",
##'     size = 5,
##'     min_dist = 1
##'   )
##'
##' summary(cross_validation)
##'
##' # user method with a matrix
##'  cross_validation <-
##'   assess_prediction(
##'   list(fit),
##'   method = "user",
##'   user_split = matrix(
##'     sample(c(rep(1, 50), rep(0, 50))),
##'     ncol = 1)
##'   )
##'
##' summary(cross_validation)
##'
##' # user method with a list
##'  cross_validation <-
##'   assess_prediction(
##'   list(fit),
##'   method = "user",
##'   user_split = list(
##'     sample(100, 50))
##'   )
##'
##' summary(cross_validation)
##'
##' @export
assess_prediction <- function(object,
                              method,
                              fold = NULL,
                              min_dist = NULL,
                              size = NULL,
                              user_split = NULL,
                              iter = 1,
                              metrics = c("AnPIT", "CRPS", "SCRPS"),
                              keep_par_fixed = TRUE,
                              control_mcmc = set_control_mcmc(),
                              plot_fold = TRUE,
                              messages = TRUE,
                              ...) {

  ## ─────────────────────────── helpers ─────────────────────────── ##
  is_list_of_riskmap <- function(x) {
    is.list(x) && all(vapply(x, inherits, logical(1), what = "RiskMap"))
  }

  crps_gaussian <- function(y, mu, sigma) {
    if (sigma == 0) return(0)
    z <- (y - mu) / sigma
    2 * dnorm(z) + z * (2 * pnorm(z) - 1) - 1 / sqrt(pi)
  }
  crps_discrete <- function(y, Fk) {
    k <- seq_along(Fk) - 1
    sum((Fk - as.numeric(k >= y))^2)
  }
  exp_crps_discrete <- function(Fk, pk) {
    k <- seq_along(pk) - 1
    sum(pk * vapply(k, crps_discrete, numeric(1), Fk = Fk))
  }
  u_val <- seq(0, 1, length.out = 1000)

  ## ────────────────────── sanity checks (unchanged) ─────────────────────── ##
  if (!is_list_of_riskmap(object))
    stop("'object' must be a list of fitted models of class 'RiskMap'.")

  if (!all(metrics %in% c("CRPS", "SCRPS", "AnPIT")))
    stop("'metrics' must only contain 'CRPS', 'SCRPS' or 'AnPIT'")

  if (!method %in% c("cluster", "regularized", "user"))
    stop("'method' must be either 'cluster', 'regularized' or 'user'")

  if (method == "cluster"){
    if (is.null(fold)) stop("when 'method' is 'cluster' you must supply 'fold'")
    check_positive_integer(fold, "fold")
  }

  if (method == "regularized") {
    if (is.null(min_dist)) stop("when 'method' is 'regularized' you must supply 'min_dist'")
    if (is.null(size))     stop("when 'method' is 'regularized' you must supply 'size'")
    check_positive_number(min_dist, "")
    check_positive_integer(size, "size")
  }

  if (method == "user"){
    if (is.null(user_split)) stop("when 'method' is 'user' you must supply 'user_split'")
  }

  check_logical(keep_par_fixed)

  check_positive_integer(iter, "iter")

  if (!inherits(control_mcmc, "RiskMap_control_mcmc"))
    stop("'control_mcmc' must come from 'set_control_mcmc()'")

  check_logical(plot_fold)

  check_logical(messages)

  get_CRPS  <- "CRPS"  %in% metrics
  get_SCRPS <- "SCRPS" %in% metrics
  get_AnPIT <- "AnPIT" %in% metrics

  ## ───────────────────────────── data & splits ───────────────────────────── ##

  # add names to any unnamed models
  if (is.null(names(object))) {
    names(object) <- paste("Model", seq_along(object))
  } else {
    unnamed <- which(names(object) == "")
    if (length(unnamed) > 0) {
      names(object)[unnamed] <- paste("Model", unnamed)
    }
  }

  object1 <- object[[1]]
  data_sf <- object1$data
  n_obs   <- nrow(data_sf)
  data_geom <- st_geometry(data_sf)

  for (h in seq_along(object)) {
    fit_data <- object[[h]]$data
    if (nrow(fit_data) != n_obs) {
      stop("All models in 'object' supplied must have the same number of observations")
    }
    fit_geom <- st_geometry(fit_data)
    if (!identical(fit_geom, data_geom)) {
      stop("All models in 'object' must have data in the same row order and geometry")
    }
  }

  # Validates one vector of row indices for a 'user_split' entry: non-missing
  # numeric, non-empty, whole numbers, no duplicates, in range. `what` names
  # whichever vector the caller actually supplied, so it reads sensibly for
  # both call sites below (a bare vector of test indices, or an explicit
  # 'in_id'/'out_id').
  validate_split_vector <- function(v, what) {
    if (!is.numeric(v) || anyNA(v))
      stop(what, " must be a non-missing numeric vector of row indices")
    if (length(v) == 0)
      stop(what, " must be non-empty")
    if (any(v != round(v)))
      stop(what, " must contain whole numbers")
    if (anyDuplicated(v))
      stop(what, " must not contain duplicate indices")
    if (!all(v %in% seq_len(n_obs)))
      stop(what, " must be row indices between 1 and the number of observations")
    invisible(TRUE)
  }

  # Validates one user-supplied 'user_split' entry: either a bare vector of
  # test indices (the shorthand form - training indices are its computed
  # complement, which can never independently be invalid) or an explicit
  # list(in_id = ..., out_id = ...) pair, both directly supplied by the
  # caller. For the latter, also checks train and test don't overlap - an
  # overlap would silently leak a held-out observation into training.
  validate_split_indices <- function(out_id, index, in_id = NULL) {
    tag <- sprintf("user_split[[%d]]", index)

    if (is.null(in_id)) {
      validate_split_vector(out_id, sprintf("'%s' test indices", tag))
      in_id <- setdiff(seq_len(n_obs), out_id)
      if (length(in_id) == 0)
        stop("'", tag, "' must leave at least one observation for training")
    } else {
      validate_split_vector(in_id, sprintf("'%s's 'in_id'", tag))
      validate_split_vector(out_id, sprintf("'%s's 'out_id'", tag))
      if (length(intersect(in_id, out_id)) > 0)
        stop("'", tag, "'s 'in_id' and 'out_id' must not overlap")
    }

    list(in_id = as.integer(in_id), out_id = as.integer(out_id))
  }

  make_splits_from_user <- function(usr, n_iter_expected) {
    spl <- vector("list", n_iter_expected)
    if (is.matrix(usr)) {
      if (nrow(usr) != n_obs)
        stop("'user_split' matrix must have the same number of rows as the data in the model")
      if (ncol(usr) != n_iter_expected)
        stop("'user_split' matrix must have a number of columns equal to 'iter'")
      if (anyNA(usr))
        stop("'user_split' matrix must not contain missing values")
      if (!all(usr %in% c(0, 1)))
        stop("'user_split' matrix must only contain 0s (training) and 1s (test)")
      for (i in seq_len(n_iter_expected)) {
        out_id <- which(usr[, i] == 1)
        in_id  <- which(usr[, i] == 0)
        if (length(in_id) == 0 || length(out_id) == 0)
          stop("Column ", i, " of 'user_split' must contain at least one training (0) ",
               "and one test (1) observation")
        spl[[i]] <- list(in_id = in_id, out_id = out_id,
                         data = data_sf[in_id, ],
                         data_test = data_sf[out_id, ])
      }
    } else if (inherits(usr, "list")) {
      if (length(usr) != n_iter_expected)
        stop("'user_split' list must have the same length as 'iter'")
      for (i in seq_len(n_iter_expected)) {
        ui <- usr[[i]]
        if (is.list(ui) && !is.null(ui$in_id) && !is.null(ui$out_id)) {
          validated <- validate_split_indices(ui$out_id, i, ui$in_id)
          in_id  <- validated$in_id
          out_id <- validated$out_id
        } else if (is.numeric(ui)) {
          validated <- validate_split_indices(ui, i)
          in_id  <- validated$in_id
          out_id <- validated$out_id
        } else {
          stop("Each element of 'user_split' must be a vector of test indices or a list(in_id=..., out_id=...).")
        }
        spl[[i]] <- list(in_id = in_id, out_id = out_id,
                         data = data_sf[in_id, ],
                         data_test = data_sf[out_id, ])
      }
    } else {
      stop("'user_split' must be a matrix or a list")
    }
    list(splits = spl)
  }

  if (!is.null(control_mcmc$seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(control_mcmc$seed)
  }

  if (method == "user") {
    data_split <- make_splits_from_user(user_split, iter)
    n_iter <- iter

    if (isTRUE(plot_fold)) print(plot_folds(data_split))
  } else if (method == "cluster") {
    data_split <- spatial_clustering_cv(data = data_sf, v = fold, repeats = iter, ...)
    # ensure out_id present
    all_ids <- seq_len(nrow(data_sf))
    for (i in seq_along(data_split$splits)) {
      if (is.null(data_split$splits[[i]]$out_id) || all(is.na(data_split$splits[[i]]$out_id))) {
        in_id <- data_split$splits[[i]]$in_id
        data_split$splits[[i]]$out_id <- setdiff(all_ids, in_id)
      }
    }
    n_iter <- iter * fold
    if (plot_fold) print(autoplot(data_split))
  } else { # regularized
    data_split <- list(splits = vector("list", iter))
    for (i in seq_len(iter)) {
      locations_sf <- data_sf[!duplicated(st_as_text(data_sf$geometry)), ]
      data_split$splits[[i]] <- list()
      data_split$splits[[i]]$data_test <- subsample.distance(locations_sf, size = size, d = min_dist * 1000, ...)
      test_geom <- st_as_text(data_split$splits[[i]]$data_test$geometry)
      in_test   <- st_as_text(data_sf$geometry) %in% test_geom
      data_split$splits[[i]]$out_id <- which(in_test)
      data_split$splits[[i]]$in_id  <- which(!in_test)
      data_split$splits[[i]]$data   <- data_sf[!in_test, ]
    }
    n_iter <- iter
    if (isTRUE(plot_fold)) print(plot_folds(data_split))
  }

  ## ───────────────────────── initialise output ───────────────────────── ##
  n_models    <- length(object)
  model_names <- names(object)
  out <- list(test_set = vector("list", n_iter), model = list())

  ## ─────────────────────────── iterate over models ───────────────────────── ##
  for (h in seq_len(n_models)) {

    fit0      <- object[[h]]
    fit_data_sf <- fit0$data
    par_hat   <- coef(fit0)
    den_name  <- as.character(fit0$call$denominator)
    fam       <- fit0$family
    linkfun   <- switch(fam,
                        gaussian = identity,
                        binomial = plogis,
                        poisson  = exp)

    if (messages) {
      message(sprintf("\nModel '%s' (%s)", h, model_names[h]))
    }

    ## containers for this model
    if (get_CRPS)   CRPS  <- vector("list", n_iter)
    if (get_SCRPS) { y_CRPS <- vector("list", n_iter); SCRPS <- vector("list", n_iter) }
    if (get_AnPIT) {
      AnPIT_area <- vector("list", n_iter)
      if (fam == "gaussian") PIT <- vector("list", n_iter) else AnPIT <- vector("list", n_iter)
    }

    ## ─────────────────── iterate over CV splits (refit as needed) ─────────────────── ##
    for (i in seq_len(n_iter)) {

      in_id  <- data_split$splits[[i]]$in_id
      out_id <- data_split$splits[[i]]$out_id

      ## ----- refit or slice -----
      if (!keep_par_fixed) {

        message("\nRe-estimating model for subset ", i)
        model_crs <- st_crs(fit0$data)

        refit_args <- list(
          formula        = fit0$formula,
          data           = fit_data_sf[in_id, ],
          family         = fam,
          model_crs      = model_crs,
          distance_units = fit0$distance_units,
          control_mcmc   = control_mcmc,
          fix_var_me     = fit0$fix_var_me,
          messages       = FALSE,
          start_pars     = par_hat
        )
        ## 'denominator' is passed as an unquoted column name (NSE); include it
        ## when the original model was fitted with one
        if (length(den_name) == 1 && nzchar(den_name)) {
          refit_args$denominator <- as.name(den_name)
        }
        refit_i <- do.call(glgpm, refit_args)
      } else {
        ## quick slice without re-fitting
        refit_i <- fit0
        keep <- in_id
        refit_i$data  <- refit_i$data [keep, ]
        refit_i$units_m  <- refit_i$units_m[keep]
        keep_coord <- unique(refit_i$ID_coords[keep])
        refit_i$coords   <- refit_i$coords[keep_coord, , drop = FALSE]
        refit_i$y        <- refit_i$y      [keep]
        refit_i$D        <- refit_i$D      [keep, , drop = FALSE]
        if (!is.null(refit_i$cov_offset))
          refit_i$cov_offset <- refit_i$cov_offset[keep]
        if (!is.null(refit_i$ID_re)) {
          re_terms <- names(refit_i$ID_re)
          random_effects_i <- prepare_random_effects(refit_i$data, re_terms)
          refit_i$ID_re <- as.data.frame(random_effects_i$ID_re)
          colnames(refit_i$ID_re) <- random_effects_i$names_re
          refit_i$re <- random_effects_i$re_unique_f
        }
        ## recompute ID_coords mapping
        refit_i$ID_coords <- create_ids(refit_i$data)$ID_coords
      }

      ## ----- held-out set and offsets -----
      data_test_i <- fit_data_sf[out_id, ]
      keep_test <- complete.cases(st_drop_geometry(data_test_i))
      data_test_i <- data_test_i[keep_test, ]
      out_id_i <- out_id[keep_test]
      if (h == 1) out$test_set[[i]] <- data_test_i

      pred_coff_i <- if (is.null(fit0$cov_offset) || all(fit0$cov_offset == 0)) {
        NULL
      } else {
        fit0$cov_offset[out_id_i]
      }

      if (messages) message("\nModel: ", model_names[h], "\nSpatial prediction for subset ", i)

      ## ----- prediction over test set -----
      pred_S <- setup_prediction(
        object          = refit_i,
        grid_pred       = st_as_sfc(data_test_i),
        control_mcmc     = control_mcmc,
        predictors      = data_test_i,
        pred_cov_offset = pred_coff_i,
        type            = "marginal",
        messages        = FALSE
      )

      pred_lp  <- predict_grid_target(
        pred_S,
        include_nugget     = is.null(refit_i$fix_tau2) || refit_i$fix_tau2 != 0,
        include_cov_offset = !all(refit_i$cov_offset == 0)
      )
      eta_samp <- pred_lp$lp_samples
      if (fam == "gaussian") {
        sigma2_me <- if (is.null(refit_i$fix_var_me)) coef(refit_i)$sigma2_me else refit_i$fix_var_me
        eta_samp <- eta_samp + sqrt(sigma2_me) * rnorm(length(eta_samp))
      }
      mu_samp <- linkfun(eta_samp)

      n_pred <- nrow(eta_samp)
      n_draw <- ncol(eta_samp)
      if (get_CRPS)   CRPS [[i]] <- numeric(n_pred)
      if (get_SCRPS) { y_CRPS[[i]] <- numeric(n_pred); SCRPS[[i]] <- numeric(n_pred) }

      if (get_AnPIT) {
        if (fam == "gaussian") {
          PIT_i <- numeric(n_pred)
        } else {
          AnPIT_i <- matrix(NA_real_, nrow = length(u_val), ncol = n_pred)
          npit_fun <- function(y, u, Fk) {
            F_y1 <- if (y == 0) 0 else Fk[y]
            F_y  <- Fk[y + 1]
            ifelse(u <= F_y1, 0,
                   ifelse(u <= F_y, (u - F_y1) / (F_y - F_y1), 1))
          }
        }
      }

      units_m_i <- fit0$units_m[out_id_i]
      y_i       <- fit0$y      [out_id_i]

      for (j in seq_len(n_pred)) {
        if (fam == "gaussian") {
          mu_j <- mean(mu_samp[j, ])
          sd_j <- sd(mu_samp[j, ])
          if (get_CRPS)  CRPS[[i]][j] <- sd_j * crps_gaussian(y_i[j], mu_j, sd_j)
          if (get_SCRPS) {
            y_CRPS[[i]][j] <- sd_j / sqrt(pi)
            SCRPS [[i]][j] <- -0.5 * (1 + CRPS[[i]][j] / y_CRPS[[i]][j] + log(2 * abs(y_CRPS[[i]][j])))
          }
          if (get_AnPIT) PIT_i[j] <- pnorm(y_i[j], mean = mu_j, sd = sd_j)
        } else {
          if (fam == "binomial") {
            y_samp  <- rbinom(n_draw, size = units_m_i[j], prob = mu_samp[j, ])
            support <- 0:units_m_i[j]
          } else { # Poisson
            lambda  <- units_m_i[j] * mu_samp[j, ]
            y_samp  <- rpois(n_draw, lambda)
            support <- 0:max(max(y_samp), y_i[j], qpois(0.999, mean(lambda)))
          }
          pk <- tabulate(y_samp + 1, nbins = length(support)) / n_draw
          Fk <- cumsum(pk)
          if (get_CRPS)  CRPS[[i]][j] <- crps_discrete(y_i[j], Fk)
          if (get_SCRPS) {
            y_CRPS[[i]][j] <- exp_crps_discrete(Fk, pk)
            SCRPS [[i]][j] <- -0.5 * (1 + CRPS[[i]][j] / y_CRPS[[i]][j] + log(2 * abs(y_CRPS[[i]][j])))
          }
          if (get_AnPIT) {
            AnPIT_i[, j] <- vapply(u_val, npit_fun, numeric(1), y = y_i[j], Fk = Fk)
          }
        }
      }

      if (get_AnPIT) {
        if (fam == "gaussian") {
          PIT[[i]] <- PIT_i
          AnPIT_area[[i]] <- .anpit_area(ecdf(PIT_i)(u_val), u_val)
        } else {
          AnPIT[[i]] <- rowMeans(AnPIT_i)
          AnPIT_area[[i]] <- .anpit_area(AnPIT[[i]], u_val)
        }
      }

    } # end i loop

    ## ─────────────── finalise output for this model ─────────────── ##
    out$model[[model_names[h]]] <- list(metric = list())
    if (get_CRPS)  out$model[[model_names[h]]]$metric$CRPS  <- CRPS
    if (get_SCRPS) out$model[[model_names[h]]]$metric$SCRPS <- SCRPS
    if (get_AnPIT) {
      out$model[[model_names[h]]]$metric$AnPIT_area <- AnPIT_area
      if (fam == "gaussian") out$model[[model_names[h]]]$PIT <- PIT else out$model[[model_names[h]]]$AnPIT <- AnPIT
    }

  } # end h loop

  class(out) <- "RiskMap_cross_validation"
  return(out)
}

##' @title Plot Method for RiskMap_simulation Objects
##' @description
##' Plots the linear predictor from [simulate_glgpm()] on a regular grid,
##' with sampling locations overlaid when replicated data are available.
##'
##' @param x Output from [simulate_glgpm()] including a surface.
##' @param sim The simulation index to plot.
##' @param ... Additional graphical parameters to be passed to the plotting function of the `terra` package.
##'
##' @return A plot of the simulation results.
##'
##' @importFrom terra rast rasterize vect
##'
##' @method plot RiskMap_simulation
##' @export
plot.RiskMap_simulation <- function(x, sim, ...) {

  sf_object <- simulated_surface(x, sim)

  # Points are assumed to fall on a regular lattice (e.g. from create_grid());
  # infer its cell size from the smallest gap between distinct coordinates,
  # and pad the extent by half a cell so points land at cell centres.
  coords <- st_coordinates(sf_object)
  cellsize <- c(min(diff(sort(unique(coords[, "X"])))),
               min(diff(sort(unique(coords[, "Y"])))))
  template <- rast(
    xmin = min(coords[, "X"]) - cellsize[1] / 2,
    xmax = max(coords[, "X"]) + cellsize[1] / 2,
    ymin = min(coords[, "Y"]) - cellsize[2] / 2,
    ymax = max(coords[, "Y"]) + cellsize[2] / 2,
    resolution = cellsize,
    crs = st_crs(sf_object)$wkt
  )
  r <- rasterize(vect(sf_object), template, field = "linear_predictor")

  plot(r, main = paste("Simulation no.", sim), ...)
  if (!is.null(x$locations$data)) {
    points(st_coordinates(x$locations$data), pch = 20)
  }

}

##' @title Assess Simulations
##'
##' @description Evaluates models using joint data and surface simulations.
##'
##' @param obj_sim Output from [simulate_glgpm()] with both data and surface.
##'   The current assessment interface supports generating models without
##'   grouped `re()` effects or custom inverse links.
##' @param models A named list of fitted models returned by [glgpm()] or model
##'   specifications returned by [specify_glgpm()]. Each object defines a
##'   candidate model specification and is refitted to every simulated dataset;
##'   its estimated parameter values are not reused.
##' @param control_mcmc A control object for MCMC sampling, created with `set_control_mcmc()`. Default is `set_control_mcmc()`.
##' @param spatial_scale The scale(s) at which predictions are assessed: `"grid"`, `"area"`, or `c("grid", "area")`
##'   to compute both from a single fit-and-predict pass over the simulations.
##' @param messages Logical, if `TRUE` messages will be displayed during processing. Default is `TRUE`.
##' @param target_transform A function that converts a numeric matrix of linear
##'   predictors to a numeric matrix of scientific targets with the same
##'   dimensions. For example, use [plogis()] for prevalence from a binomial
##'   logit model or [exp()] for the mean of a Poisson log-link model.
##' @param area_summary A function that combines a numeric vector of grid-cell
##'   targets within one area and returns one finite numeric value, such as
##'   [base::mean()] or [base::sum()]. Required only when `spatial_scale`
##'   includes `"area"`.
##' @param boundaries An `sf` object containing only POLYGON or MULTIPOLYGON geometries, required if `spatial_scale` includes `"area"`.
##' @param col_names Column name in `boundaries` containing unique region names. If `NULL`, defaults to `"region"`.
##' @param pred_objective A character vector specifying objectives, either `"mse"`, `"classify"`, or both.
##' @param categories A numeric vector of thresholds defining categories for classification. Required if `pred_objective = "classify"`.
##'
##' @return A list of class `RiskMap_assess_simulation`. `pred_objective` holds one element per
##'   requested `spatial_scale` (`"grid"` and/or `"area"`), each in turn holding `mse` and/or
##'   `classify` per the requested `pred_objective`.
##'
##' @examples
##' library(sf)
##' data(italy_sim)
##' italy_subset <- italy_sim[1:30, ]
##' italy_grid <- italy_subset[!duplicated(st_coordinates(italy_subset)), ]
##' model <- specify_glgpm(
##'   y ~ pop_dens + gp(), italy_subset, family = "gaussian",
##'   parameters = list(beta = c(1, 0.001), sigma2 = 1, phi = 20,
##'                     sigma2_me = 0.1)
##' )
##' simulations <- simulate_glgpm(
##'   model, nsim = 1, what = c("data", "surface"),
##'   prediction_grid = italy_grid, seed = 1
##' )
##' boundary <- st_sf(
##'   region = "study_area",
##'   geometry = st_convex_hull(st_union(italy_grid))
##' )
##' assessment <- assess_simulation(
##'   simulations,
##'   models = list(candidate = model),
##'   spatial_scale = c("grid", "area"),
##'   target_transform = exp,
##'   area_summary = mean,
##'   boundaries = boundary,
##'   pred_objective = "mse",
##'   messages = FALSE
##' )
##'
##' @export
assess_simulation <- function(obj_sim,
                       models,
                       control_mcmc = set_control_mcmc(),
                       spatial_scale,
                       messages = TRUE,
                       target_transform = NULL,
                       area_summary = NULL,
                       boundaries = NULL,
                       col_names = NULL,
                       pred_objective = c("mse", "classify"),
                       categories = NULL) {

  if (!inherits(obj_sim, "RiskMap_simulation")) {
    stop("'obj_sim' must be output from simulate_glgpm()")
  }
  if (!is.character(pred_objective) || length(pred_objective) == 0 ||
      anyDuplicated(pred_objective) ||
      !all(pred_objective %in% c("mse", "classify"))) {
    stop("'pred_objective' must be either 'mse', 'classify' or c('mse', 'classify')")
  }
  if (!is.character(spatial_scale) || length(spatial_scale) == 0 ||
      anyDuplicated(spatial_scale) ||
      !all(spatial_scale %in% c("grid", "area"))) {
    stop("'spatial_scale' must be set to 'grid', 'area', or c('grid', 'area')")
  }

  want_grid <- "grid" %in% spatial_scale
  want_area <- "area" %in% spatial_scale
  want_mse <- "mse" %in% pred_objective
  want_classify <- "classify" %in% pred_objective

  if (!is.list(models) || !length(models) ||
      is.null(names(models)) || any(!nzchar(names(models))) ||
      anyDuplicated(names(models)) ||
      !all(vapply(models, function(model) {
        inherits(model, "RiskMap") ||
          inherits(model, "RiskMap_simulation_model")
      }, logical(1)))) {
    stop("'models' must be a non-empty, uniquely named list of objects returned by glgpm() or specify_glgpm().")
  }

  if (want_area) {
    if (is.null(boundaries)) {
      stop("if spatial_scale includes 'area' then an sf object of the area(s) must be passed to
           'boundaries'")
    }
    check_data(boundaries, "polygon")
  }

  obj_sim <- simulation_assessment_data(obj_sim)
  if (!is.function(target_transform)) {
    stop("'target_transform' must be a function that transforms the linear predictor.")
  }

  apply_target_transform <- function(x, context) {
    result <- tryCatch(
      target_transform(x),
      error = function(e) {
        stop("'target_transform' failed for ", context, ": ",
             conditionMessage(e), call. = FALSE)
      }
    )
    if (!is.numeric(result) || !is.matrix(result) ||
        !identical(dim(result), dim(x)) || anyNA(result) ||
        any(!is.finite(result))) {
      stop("'target_transform' must return a finite numeric matrix with the same dimensions as its input (failed for ",
           context, ").", call. = FALSE)
    }
    result
  }

  if (want_classify) {
    if (is.null(categories)) stop("if 'pred_objective' is 'classify', a value for 'categories' must be specified")
    if (!is.numeric(categories) || length(categories) < 3 ||
        anyNA(categories) || any(!is.finite(categories)) ||
        any(diff(categories) <= 0)) {
      stop("'categories' must contain at least three unique, strictly increasing values.")
    }
  }
  n_sim <- length(obj_sim$data_sim)
  n_models <- length(models)

  if(want_area && !is.function(area_summary)) {
    stop("'area_summary' must be a function when 'spatial_scale' includes 'area'.")
  }

  apply_area_summary <- function(x, context) {
    result <- tryCatch(
      area_summary(x),
      error = function(e) {
        stop("'area_summary' failed for ", context, ": ",
             conditionMessage(e), call. = FALSE)
      }
    )
    if (!is.numeric(result) || length(result) != 1L || is.na(result) ||
        !is.finite(result)) {
      stop("'area_summary' must return one finite numeric value (failed for ",
           context, ").", call. = FALSE)
    }
    as.numeric(result)
  }
  model_names <- names(models)

  include_covariates <- obj_sim$include_covariates
  include_cov_offset <- obj_sim$include_cov_offset
  include_nugget <- obj_sim$nugget_over_grid

  no_comp <- NULL

  # A joint prediction is required for area-level aggregation, since it needs
  # spatially correlated samples across the grid; it also carries everything
  # a marginal grid-level assessment needs, so requesting both scales still
  # only takes one fit-and-predict pass over the simulations (#109).
  type <- if (want_area) "joint" else "marginal"

  if (want_area) {
    n_reg <- nrow(boundaries)

    if(is.null(col_names)) {
      boundaries$region <- paste("reg",1:n_reg, sep="")
      col_names <- "region"
      names_reg <- boundaries$region
    } else {
      names_reg <- boundaries[[col_names]]
      if(n_reg != length(names_reg)) {
        stop("The names in the column identified by 'col_names' do not
         provide a unique set of names, but there are duplicates")
      }
    }
    boundaries <- st_transform(boundaries, st_crs(obj_sim$lp_grid_sim))
    inter <- st_intersects(boundaries, obj_sim$lp_grid_sim)
  }

  n_samples <- (control_mcmc$n_sim-control_mcmc$burnin)/control_mcmc$thin
  n_pred <- nrow(obj_sim$lp_grid_sim)

  if (want_classify) {
    category_labels <- paste0("(", utils::head(categories, -1), ",",
                              categories[-1], "]")
    categories_class <- factor(category_labels, levels = category_labels)
  }

  # One `mse`/`classify` store per requested spatial scale, all sharing the
  # same layout - only what each store is computed from differs below.
  init_objective_store <- function() {
    store <- list()
    if(want_mse) {
      store$mse <- array(NA, c(n_models, n_sim))
      rownames(store$mse) <- model_names
      colnames(store$mse) <- paste0("sim_", 1:n_sim)
    }
    if(want_classify) {
      store$classify <- setNames(vector("list", length(model_names)), model_names)
      for(i in 1:n_models) {
        store$classify[[model_names[i]]] <- list(by_cat = vector("list", n_sim),
                                                  across_cat = list())
        for(j in 1:n_sim) {
          store$classify[[model_names[i]]]$by_cat[[j]] <-
            data.frame(
              Class = categories_class,
              Sensitivity = NA,
              Specificity = NA,
              PPV = NA,
              NPV = NA,
              CC = NA
            )
        }
        store$classify[[model_names[i]]]$CC <- rep(NA, n_sim)
      }
    }
    store
  }

  out <- list(pred_objective = list())
  if(want_grid) out$pred_objective$grid <- init_objective_store()
  if(want_area) out$pred_objective$area <- init_objective_store()

  # Update one model/simulation result. The helper below keeps every category
  # in the confusion matrix, including categories absent from this simulation.
  update_classify_store <- function(store, model_name, sim_index, true_vals, samples) {
    metrics <- simulation_classification_metrics(
      true_vals, samples, categories, levels(categories_class)
    )
    store$classify[[model_name]]$by_cat[[sim_index]] <- metrics$by_cat
    store$classify[[model_name]]$CC[sim_index] <- metrics$overall_cc
    store
  }

  lp_true_sim <- as.matrix(st_drop_geometry(obj_sim$lp_grid_sim[, grepl("^lp_sim_[0-9]+$",
                                                                       names(obj_sim$lp_grid_sim))]))

  true_target_grid_sim <- apply_target_transform(
    lp_true_sim, "the simulated true surface"
  )

  if(want_area) {
    true_target_area_sim <- matrix(NA, nrow = n_reg, ncol = n_sim)
    for(i in 1:n_reg) {
      for(j in 1:n_sim) {
        if(length(inter[[i]]) == 0) {
          warning(paste("No points on the grid fall within", boundaries[[col_names]][i],
                        "and no predictions are carried out for this area"))
          no_comp <- c(no_comp, i)
        } else {
          true_target_area_sim[i,j] <- apply_area_summary(
            true_target_grid_sim[inter[[i]], j],
            paste0("area '", boundaries[[col_names]][i],
                   "' in true simulation ", j)
          )
        }
      }
    }
  }

  for(i in 1:n_models) {
    if (messages) message("Model: ", model_names[i], "\n")

    template_i <- models[[i]]
    if (template_i$family != obj_sim$family) {
      stop("Model '", model_names[i], "' uses family '", template_i$family,
           "', but the simulations use family '", obj_sim$family, "'.")
    }
    if_i <- interpret_formula(template_i$formula)
    rhs_terms <- attr(terms(if_i$pf), "term.labels")
    predictors_i <- if (length(rhs_terms) == 0) NULL else obj_sim$lp_grid_sim
    f_i <- update(template_i$formula, y ~ .)

    for(j in 1:n_sim) {
      if (messages) message("Processing simulation no. ", j)

      # Fit, predict and score one model/simulation pair at a time. Keeping
      # these objects local avoids retaining every fit and posterior sample
      # matrix until the complete assessment has finished.
      if (messages) message("Estimation")
      refit_args <- assessment_refit_args(
        template_i, f_i, obj_sim$data_sim[[j]], control_mcmc
      )
      fit_ij <- do.call(glgpm, refit_args)

      if (messages) message("Prediction over the grid")
      obj_pred_ij <- setup_prediction(
        fit_ij,
        grid_pred = st_as_sfc(obj_sim$lp_grid_sim),
        predictors = predictors_i,
        pred_cov_offset = if (is.null(if_i$offset)) NULL else
          obj_sim$lp_grid_sim[[if_i$offset]],
        control_mcmc = control_mcmc,
        type = type,
        messages = FALSE
      )

      if(length(obj_pred_ij$mu_pred) == 1 && obj_pred_ij$mu_pred == 0 &&
         include_covariates) {
        stop("Covariates have not been provided; re-run setup_prediction
         and provide the covariates through the argument 'predictors'")
      }


      if(!include_covariates) {
        mu_target <- 0
      } else {

        if(is.null(obj_pred_ij$mu_pred)) stop("the output obtained from 'setup_prediction' does not
                                     contain any covariates; if including covariates
                                     in the predictive target these should be included
                                     when running 'setup_prediction'")
        mu_target <- obj_pred_ij$mu_pred
      }

      if(!include_cov_offset) {
        cov_offset <- 0
      } else {
        if(length(obj_pred_ij$cov_offset) == 1) {
          stop("No covariate offset was included in the model;
           set include_cov_offset = FALSE, or refit the model and include
           the covariate offset")
        }
        cov_offset <- obj_pred_ij$cov_offset
      }

      if(include_nugget) {
        if(is.null(obj_pred_ij$par_hat$tau2)) stop("'include_nugget' cannot be
                                                   set to TRUE if this has not been included
                                                   in the fit of the model")
        n_samples <- ncol(obj_pred_ij$S_samples)
        Z_sim <- matrix(rnorm(n_samples*n_pred,
                              sd = sqrt(obj_pred_ij$par_hat$tau2)),
                        ncol = n_samples)
        obj_pred_ij$S_samples <- obj_pred_ij$S_samples+Z_sim
      }

      n_samples <- ncol(obj_pred_ij$S_samples)

      if(is.matrix(mu_target)) {
        lp_samples_ij <- sapply(1:n_samples,
                                function(h)
                                  mu_target[,h] + cov_offset +
                                  obj_pred_ij$S_samples[,h])
      } else {
        lp_samples_ij <- sapply(1:n_samples,
                                function(h)
                                  mu_target + cov_offset +
                                  obj_pred_ij$S_samples[,h])
      }

      # Grid-cell-level target samples are shared by both scales: grid
      # objectives use them directly, area objectives aggregate them by
      # region below.
      target_samples_ij <- apply_target_transform(
        lp_samples_ij,
        paste0("model '", model_names[i], "', simulation ", j)
      )

      if(want_grid) {
        mean_target_grid_ij <- apply(target_samples_ij, 1, mean)

        if(want_mse) {
          out$pred_objective$grid$mse[i,j] <-
            mean((mean_target_grid_ij - true_target_grid_sim[,j])^2)
        }
        if(want_classify) {
          out$pred_objective$grid <- update_classify_store(
            out$pred_objective$grid, model_names[i], j,
            true_target_grid_sim[,j], target_samples_ij)
        }
      }

      if(want_area) {
        target_area_samples_ij <- matrix(NA, nrow = n_reg, ncol = n_samples)
        mean_target_area_ij <- rep(NA,n_reg)
        for(h in 1:n_reg) {
          if(length(inter[[h]]) > 0) {
            ind_grid_h <- inter[[h]]
            target_area_samples_ij[h,] <- vapply(
              seq_len(n_samples),
              function(sample_index) {
                apply_area_summary(
                  target_samples_ij[ind_grid_h, sample_index],
                  paste0("area '", boundaries[[col_names]][h], "', model '",
                         model_names[i], "', simulation ", j,
                         ", predictive sample ", sample_index)
                )
              },
              numeric(1)
            )
            mean_target_area_ij[h] <- mean(target_area_samples_ij[h,])
          }
        }

        if(want_mse) {
          out$pred_objective$area$mse[i,j] <-
            mean((mean_target_area_ij - true_target_area_sim[,j])^2)
        }
        if(want_classify) {
          out$pred_objective$area <- update_classify_store(
            out$pred_objective$area, model_names[i], j,
            true_target_area_sim[,j], target_area_samples_ij)
        }
      }
    }
  }
  if(want_classify) {
    if(want_grid) out$pred_objective$grid$classify$Class <- categories_class
    if(want_area) out$pred_objective$area$classify$Class <- categories_class
  }
  out$n_sim <- n_sim
  out$spatial_scale <- spatial_scale
  class(out) <- "RiskMap_assess_simulation"
  return(out)
}

##' Build a refit call from a fitted assessment-model template
##'
##' The fitted coefficients are intentionally excluded. Data-dependent inputs
##' are resolved against each simulated dataset rather than copied as vectors
##' from the original fit.
##' @noRd
assessment_refit_args <- function(template, formula, data, control_mcmc) {
  fitted_template <- inherits(template, "RiskMap")
  args <- list(
    formula = formula,
    family = template$family,
    data = data,
    distance_units = template$distance_units,
    control_mcmc = control_mcmc,
    control_mcml = if (fitted_template) {
      attr(template, "control_mcml") %||% set_control_mcml()
    } else {
      set_control_mcml()
    },
    fix_var_me = if (fitted_template) template$fix_var_me else NULL,
    messages = FALSE
  )
  if (template$family != "gaussian") {
    args$denominator <- quote(units_m)
    if (fitted_template && identical(template$link_function$name, "custom")) {
      args$invlink <- template$link_function[c("inv", "d1", "d2")]
    } else if (!fitted_template && isTRUE(template$custom_link)) {
      args$invlink <- template$invlink
    }
  }
  args
}

##' Classification metrics for one simulated dataset
##'
##' @param true_values True target values.
##' @param samples Matrix of predictive target samples, one row per target.
##' @param breaks Category boundaries.
##' @param labels Fixed category labels.
##' @return A list containing category-specific metrics and overall accuracy.
##' @noRd
simulation_classification_metrics <- function(true_values, samples, breaks, labels) {
  true_class <- cut(true_values, breaks = breaks, labels = labels)
  true_class <- factor(true_class, levels = labels)

  probabilities <- do.call(cbind, lapply(seq_len(length(breaks) - 1L), function(i) {
    rowMeans(breaks[i] < samples & samples <= breaks[i + 1L])
  }))
  predicted_class <- factor(labels[max.col(probabilities, ties.method = "first")],
                            levels = labels)
  confusion <- table(true_class, predicted_class)
  total <- sum(confusion)

  divide_or_na <- function(numerator, denominator) {
    if (denominator == 0) NA_real_ else numerator / denominator
  }
  by_category <- lapply(seq_along(labels), function(i) {
    true_positive <- confusion[i, i]
    false_positive <- sum(confusion[, i]) - true_positive
    false_negative <- sum(confusion[i, ]) - true_positive
    true_negative <- total - true_positive - false_positive - false_negative
    data.frame(
      Class = labels[i],
      Sensitivity = divide_or_na(true_positive, true_positive + false_negative),
      Specificity = divide_or_na(true_negative, true_negative + false_positive),
      PPV = divide_or_na(true_positive, true_positive + false_positive),
      NPV = divide_or_na(true_negative, true_negative + false_negative),
      CC = divide_or_na(true_positive + true_negative, total)
    )
  })

  valid <- !is.na(true_class)
  overall_cc <- if (any(valid)) {
    mean(true_class[valid] == predicted_class[valid])
  } else {
    NA_real_
  }
  list(by_cat = do.call(rbind, by_category), overall_cc = overall_cc)
}

##' @title Summarize Simulation Results
##'
##' @description Summarizes the results of model evaluations from a `RiskMap_assess_simulation` object. Provides average metrics for classification by category and overall correct classification (CC) summary.
##'
##' @param object An object of class `RiskMap_assess_simulation`, as returned by `assess_simulation`.
##' @param ... Additional arguments (not used).
##'
##' @return A list containing summary data for each model:
##' - `by_cat_summary`: A data frame with average sensitivity, specificity, PPV, NPV, and CC by category.
##' - `CC_summary`: A numeric vector with mean, 2.5th percentile, and 97.5th percentile for CC across simulations.
##'
##' @method summary RiskMap_assess_simulation
##' @export
summary.RiskMap_assess_simulation <- function(object, ...) {
  stopifnot(inherits(object, "RiskMap_assess_simulation"))

  if (identical(object$n_sim, 1L)) {
    warning("The assessment contains one simulation. Point estimates are shown, ",
            "but across-simulation uncertainty cannot be estimated.", call. = FALSE)
  }

  # `object$pred_objective` holds one element per requested spatial scale
  # ("grid" and/or "area"); each is summarized the same way.
  results <- list()
  for (scale in intersect(c("grid", "area"), names(object$pred_objective))) {
    results[[scale]] <- summarize_pred_objective(object$pred_objective[[scale]])
  }

  # Assign class for S3 print method
  class(results) <- "summary.RiskMap_assess_simulation"
  return(results)
}

##' Summarize one spatial scale's `mse`/`classify` results
##'
##' @param pred_objective The `mse`/`classify` element of a
##'   `RiskMap_assess_simulation` object for a single spatial scale.
##' @return A list with `mse` and/or `classify` summaries.
##' @noRd
summarize_pred_objective <- function(pred_objective) {
  results <- list()

  # Check for "mse" in pred_objective
  if ("mse" %in% names(pred_objective)) {
    mse_data <- pred_objective$mse

    # Check if mse_data is a matrix
    if (is.matrix(mse_data)) {
      n_valid <- rowSums(is.finite(mse_data))
      mse_summary <- data.frame(
        Model = rownames(mse_data),
        n_sim = ncol(mse_data),
        n_valid = n_valid,
        MSE_mean = apply(mse_data, 1, finite_mean),
        MSE_sd = apply(mse_data, 1, finite_sd)
      )

      results$mse <- mse_summary
    } else {
      stop("mse_data must be a matrix.")
    }
  }

  # Check for "classify" in pred_objective
  if ("classify" %in% names(pred_objective)) {
    classify_data <- pred_objective$classify

    name_models <- setdiff(names(classify_data), "Class")
    results$classify <- list()
    for(model_name in name_models) {
      model_data <- classify_data[[model_name]]
      n_sim <- length(model_data$by_cat)
      metric_names <- setdiff(names(model_data$by_cat[[1]]), "Class")
      classes <- model_data$by_cat[[1]]$Class
      metric_array <- array(
        NA_real_,
        dim = c(length(classes), length(metric_names), n_sim),
        dimnames = list(classes, metric_names, paste0("sim_", seq_len(n_sim)))
      )
      for(j in seq_len(n_sim)) {
        current <- model_data$by_cat[[j]]
        row_index <- match(classes, current$Class)
        metric_array[, , j] <- as.matrix(current[row_index, metric_names, drop = FALSE])
      }
      metric_mean <- apply(metric_array, c(1, 2), finite_mean)
      metric_n_valid <- apply(is.finite(metric_array), c(1, 2), sum)
      classify_res <- data.frame(Class = classes, metric_mean,
                                 check.names = FALSE, row.names = NULL)
      n_valid <- data.frame(Class = classes, metric_n_valid,
                            check.names = FALSE, row.names = NULL)

      results$classify[[model_name]] <- list(
        classify_res = classify_res,
        n_valid = n_valid,
        n_sim = n_sim,
        cc_summary = summarize_simulation_metric(model_data$CC, n_sim)
      )
    }
  }

  results
}

##' @noRd
finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

##' @noRd
finite_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) >= 2L) stats::sd(x) else NA_real_
}

##' @noRd
summarize_simulation_metric <- function(x, n_sim = length(x)) {
  x <- x[is.finite(x)]
  n_valid <- length(x)
  list(
    mean = if (n_valid) mean(x) else NA_real_,
    sd = if (n_valid >= 2L) stats::sd(x) else NA_real_,
    lower = if (n_valid >= 2L) unname(stats::quantile(x, 0.025)) else NA_real_,
    upper = if (n_valid >= 2L) unname(stats::quantile(x, 0.975)) else NA_real_,
    n_valid = n_valid,
    n_sim = n_sim
  )
}



##' @title Print Simulation Results
##'
##' @description Prints a concise summary of simulation results from a `RiskMap_assess_simulation` object, including average metrics by category and a summary of overall correct classification (CC).
##'
##' @param x An object of class `summary.RiskMap_assess_simulation`, as returned by `summary.RiskMap_assess_simulation`.
##' @param ... Additional arguments (not used).
##'
##' @return Invisibly returns `x`.
##'
##'
##' Print Simulation Results
##'
##' Prints a concise summary of simulation results from a `summary.RiskMap_assess_simulation` object,
##' including average metrics by category and a summary of overall correct classification (CC).
##'
##' @param x An object of class `summary.RiskMap_assess_simulation`, as returned by `summary.RiskMap_assess_simulation`.
##' @param ... Additional arguments (not used).
##'
##' @return Invisibly returns `x`.
##'
##' @method print summary.RiskMap_assess_simulation
##' @export
print.summary.RiskMap_assess_simulation <- function(x, ...) {
  cat("Summary of Simulation Results\n\n")

  scale_label <- c(grid = "Grid", area = "Area")
  for (scale in intersect(c("grid", "area"), names(x))) {
    if (length(x) > 1) cat(sprintf("== %s-level results ==\n\n", scale_label[[scale]]))
    print_pred_objective(x[[scale]])
  }

  invisible(x)
}

##' Print a data frame as plain text
##'
##' Mirrors `print.data.frame()` but prints the formatted character matrix
##' directly, so RStudio notebooks keep the table inline with the surrounding
##' `cat()` output instead of rendering it as a separate data frame widget.
##'
##' @param df A data frame.
##' @return Invisibly returns `df`.
##' @noRd
print_text_table <- function(df) {
  print(as.matrix(format(df)), quote = FALSE, right = TRUE)
  invisible(df)
}

##' Print one spatial scale's `mse`/`classify` summary
##'
##' @param x An element of a `summary.RiskMap_assess_simulation` object, as
##'   returned by `summarize_pred_objective()`.
##' @return Invisibly returns `x`.
##' @noRd
print_pred_objective <- function(x) {
  if (!is.null(x$mse)) {
    cat("Mean Squared Error (MSE):\n")
    print_text_table(x$mse)
    cat("\n")
  }

  if (!is.null(x$classify)) {
    cat("Classification Results:\n")

    # Iterate over each model in classify results
    for (model_name in names(x$classify)) {
      model_data <- x$classify[[model_name]]

      cat(sprintf("\nModel: %s\n", model_name))

      cat("\nAverages across simulations by Category:\n")
      print_text_table(model_data$classify_res)

      cat("\nNumber of valid simulations by Category and metric ",
          sprintf("(out of %d):\n", model_data$n_sim), sep = "")
      print_text_table(model_data$n_valid)

      cat("\nProportion of Correct Classification (CC) across categories:\n")
      cc_summary <- model_data$cc_summary
      if (cc_summary$n_valid >= 2L) {
        cat(sprintf("Mean: %.3f, SD: %.3f, 95%% simulation interval: [%.3f, %.3f] ",
                    cc_summary$mean, cc_summary$sd,
                    cc_summary$lower, cc_summary$upper))
      } else {
        cat(sprintf("Mean: %.3f, uncertainty unavailable ", cc_summary$mean))
      }
      cat(sprintf("(n_valid = %d of %d)\n",
                  cc_summary$n_valid, cc_summary$n_sim))
    }
    cat("\n")
  }

  invisible(x)
}
