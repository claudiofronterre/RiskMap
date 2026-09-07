##' @title Assess Predictive Performance via Spatial Cross-Validation
##'
##' @description
##' This function evaluates the predictive performance of spatial models fitted to `RiskMap` objects using cross-validation. It supports two classes of diagnostic tools:
##'
##' - **Scoring rules**, including the Continuous Ranked Probability Score (CRPS) and its scaled version (SCRPS), which quantify the sharpness and calibration of probabilistic forecasts;
##' - **Calibration diagnostics**, based on the Probability Integral Transform (PIT) for Gaussian outcomes, Aggregated nonparametric PIT (AnPIT) curves for discrete outcomes (e.g., Poisson or Binomial), and the area between the PIT/AnPIT curve and the reference line.
##'
##' Cross-validation can be performed using either spatial clustering or regularized subsampling with a minimum inter-point distance. For each fold or subset, models can be refitted or evaluated with fixed parameters, offering flexibility in model validation. The function also provides visualizations of the spatial distribution of test folds.
##'
##' @param object A list of `RiskMap` objects, each representing a model fitted with `glgpm`.
##' @param method Character; either `"cluster"` or `"regularized"` for the cross-validation method. The `"cluster"` method uses
##' spatial clustering as implemented by the \code{spatial_clustering_cv} function from the `spatialEco` package, while the `"regularized"` method
##' selects a subsample of the dataset by imposing a minimum distance, set by the `min_dist` argument, for a randomly selected
##' subset of locations.
##' @param keep_par_fixed Logical; if `TRUE`, parameters are kept fixed across folds, otherwise the model is re-estimated for each fold.
##' @param iter Integer; number of times to repeat the cross-validation.
##' @param fold Integer; number of folds for cross-validation (required if `method = "cluster"`).
##' @param n_size Optional; the size of the test set, required if `method = "regularized"`.
##' @param control_sim Control settings for simulation, an output from `set_control_mcmc`.
##' @param min_dist Optional; minimum distance for regularized subsampling (required if `method = "regularized"`).
##' @param plot_fold Logical; if `TRUE`, plots each fold's test set.
##' @param messages Logical; if `TRUE`, displays progress messages.
##' @param which_metric Character vector; one or more of `"CRPS"`, `"SCRPS"`, or `"AnPIT"`, to specify the predictive performance metrics to compute. When `"AnPIT"` is requested, the scalar score `"AnPIT_area"` is also computed as the integrated absolute deviation between the PIT/AnPIT curve and the reference line.
##' @param user_split A user-defined cross-validation split. Either:
##'   * a matrix with \code{nrow = n} (number of observations) and
##'     \code{ncol = iter} (number of iterations), where entries of \code{1}
##'     indicate membership in the test set for that iteration and \code{0}
##'     indicate training set; or
##'   * a list of length \code{iter}, where each element is either a vector of
##'     test indices, or a list with components \code{in_id} (training indices)
##'     and \code{out_id} (test indices).
##'   When supplied, \code{user_split} overrides the automatic clustering or
##'   regularized distance splitting defined by \code{method}.
##' @param ... Additional arguments passed to clustering or subsampling functions.
##'
##' @return A list of class `RiskMap_cross_validation`, containing:
##' \describe{
##'   \item{test_set}{A list of test sets used for validation, each of class `'sf'`.}
##'   \item{model}{A named list, one per model, each containing:
##'     \describe{
##'       \item{score}{A list with CRPS, SCRPS, and/or AnPIT area scores for each fold if requested.}
##'       \item{PIT}{(if `family = "gaussian"` and `which_metric` includes `"AnPIT"`) A list of PIT values for test data.}
##'       \item{AnPIT}{(if `family` is discrete and `which_metric` includes `"AnPIT"`) A list of AnPIT curves for test data.}
##'     }
##'   }
##' }
##'
##' @seealso \code{\link{plot_AnPIT}}
##'
##' @references
##' Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of proper scoring rules. *Statistical Science*, 38(1), 140–159. \doi{10.1214/22-STS864}.
##'
##' @importFrom terra match
##' @importFrom gridExtra grid.arrange
##' @importFrom spatialEco subsample.distance
##' @importFrom spatialsample spatial_clustering_cv autoplot
##' @export
assess_prediction <- function(object,
                              method,
                              keep_par_fixed = TRUE,
                              iter = 1,
                              fold = NULL,
                              n_size = NULL,
                              control_sim = set_control_mcmc(),
                              min_dist = NULL,
                              plot_fold = TRUE,
                              messages = TRUE,
                              which_metric = c("AnPIT", "CRPS", "SCRPS"),
                              user_split = NULL,
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

  if(!is.null(user_split)) {
    iter <- ncol(user_split)
  }
  ## ────────────────────── sanity checks (unchanged) ─────────────────────── ##
  if (!is_list_of_riskmap(object))
    stop("'object' must be a list of fitted models of class 'RiskMap'.")

  if (!all(which_metric %in% c("CRPS", "SCRPS", "AnPIT")))
    stop("'which_metric' must only contain 'CRPS', 'SCRPS' or 'AnPIT'")

  if (is.null(user_split)) {
    if (!method %in% c("cluster", "regularized"))
      stop("'method' must be either 'cluster' or 'regularized' (unless 'user_split' is supplied).")

    if (method == "regularized") {
      if (is.null(min_dist)) stop("for 'regularized', supply 'min_dist'")
      if (is.null(n_size))   stop("for 'regularized', supply 'n_size'")
    }
    if (method == "cluster" && is.null(fold))
      stop("when 'method' is 'cluster', you must supply 'fold'")
  }

  if (!inherits(control_sim, "RiskMap_control_mcmc"))
    stop("'control_sim' must come from 'set_control_mcmc()'")

  get_CRPS  <- "CRPS"  %in% which_metric
  get_SCRPS <- "SCRPS" %in% which_metric
  get_AnPIT <- "AnPIT" %in% which_metric

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
  data_sf <- object1$data_sf
  n_obs   <- nrow(data_sf)
  data_geom <- st_as_text(st_geometry(data_sf))

  for (h in seq_along(object)) {
    fit_data <- object[[h]]$data_sf
    if (nrow(fit_data) != n_obs) {
      stop("All models supplied to 'assess_prediction()' must have the same number of observations.")
    }
    fit_geom <- st_as_text(st_geometry(fit_data))
    if (!identical(fit_geom, data_geom)) {
      stop("All models supplied to 'assess_prediction()' must have data in the same row order and geometry.")
    }
  }

  make_splits_from_user <- function(usr, n_iter_expected) {
    spl <- vector("list", n_iter_expected)
    if (is.matrix(usr)) {
      if (nrow(usr) != n_obs)
        stop("'user_split' matrix must have nrow == nrow(data).")
      if (ncol(usr) != n_iter_expected)
        stop("'user_split' matrix must have ncol == 'iter'.")
      for (i in seq_len(n_iter_expected)) {
        out_id <- which(usr[, i] != 0 & !is.na(usr[, i]))
        in_id  <- setdiff(seq_len(n_obs), out_id)
        spl[[i]] <- list(in_id = in_id, out_id = out_id,
                         data = data_sf[in_id, ],
                         data_test = data_sf[out_id, ])
      }
    } else if (is.list(usr)) {
      if (length(usr) != n_iter_expected)
        stop("'user_split' list must have length == 'iter'.")
      for (i in seq_len(n_iter_expected)) {
        ui <- usr[[i]]
        if (is.list(ui) && !is.null(ui$in_id) && !is.null(ui$out_id)) {
          in_id  <- ui$in_id
          out_id <- ui$out_id
        } else if (is.integer(ui) || is.double(ui)) {
          out_id <- as.integer(ui)
          in_id  <- setdiff(seq_len(n_obs), out_id)
        } else {
          stop("Each element of 'user_split' must be a vector of test indices or a list(in_id=..., out_id=...).")
        }
        spl[[i]] <- list(in_id = in_id, out_id = out_id,
                         data = data_sf[in_id, ],
                         data_test = data_sf[out_id, ])
      }
    } else {
      stop("'user_split' must be a matrix (nrow=n, ncol=iter) or a list.")
    }
    list(splits = spl)
  }

  if (!is.null(user_split)) {
    data_split <- make_splits_from_user(user_split, iter)
    n_iter <- iter

    if (isTRUE(plot_fold)) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("plot_fold = TRUE requires the 'ggplot2' package; skipping plots.", call. = FALSE)
      } else {
        if (n_iter == 1) {
          p <- ggplot(data_split$splits[[1]]$data_test) +
            geom_sf() +
            theme_minimal() +
            ggtitle("Test set")
          print(p)
        } else {
          plots <- lapply(seq_len(n_iter), function(i) {
            ggplot(data_split$splits[[i]]$data_test) +
              geom_sf() +
              theme_minimal() +
              ggtitle(paste("Test", i))
          })

          if (requireNamespace("gridExtra", quietly = TRUE)) {
            do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
          } else {
            warning("Optional package 'gridExtra' not installed; printing plots sequentially.", call. = FALSE)
            for (p in plots) print(p)
          }
        }
      }
    }
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
      data_split$splits[[i]]$data_test <- subsample.distance(locations_sf, size = n_size, d = min_dist * 1000, ...)
      test_geom <- st_as_text(data_split$splits[[i]]$data_test$geometry)
      in_test   <- st_as_text(data_sf$geometry) %in% test_geom
      data_split$splits[[i]]$out_id <- which(in_test)
      data_split$splits[[i]]$in_id  <- which(!in_test)
      data_split$splits[[i]]$data   <- data_sf[!in_test, ]
    }
    n_iter <- iter
    if (isTRUE(plot_fold)) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("plot_fold = TRUE requires the 'ggplot2' package; skipping plots.", call. = FALSE)
      } else if (!requireNamespace("sf", quietly = TRUE)) {
        warning("plot_fold = TRUE with geom_sf() requires the 'sf' package; skipping plots.", call. = FALSE)
      } else {
        plots <- lapply(seq_len(n_iter), function(i) {
          ggplot(data_split$splits[[i]]$data_test) +
            geom_sf() +
            theme_minimal() +
            ggtitle(paste("Subset", i))
        })

        if (n_iter > 1 && requireNamespace("gridExtra", quietly = TRUE)) {
          do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
        } else {
          # Either only one plot or gridExtra not available: print sequentially
          for (p in plots) print(p)
        }
      }
    }
  }

  ## ───────────────────────── initialise output ───────────────────────── ##
  n_models    <- length(object)
  model_names <- names(object)
  out <- list(test_set = vector("list", n_iter), model = list())

  ## ─────────────────────────── iterate over models ───────────────────────── ##
  for (h in seq_len(n_models)) {

    fit0      <- object[[h]]
    fit_data_sf <- fit0$data_sf
    par_hat   <- coef(fit0)
    den_name  <- as.character(fit0$call$den)
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
        crs_num <- if (is.numeric(fit0$crs)) {
          as.integer(fit0$crs)
        } else {
          as.integer(sub(".*:(\\d+)$", "\\1", fit0$crs))
        }

        fit0$crs <- crs_num
        ## Original path (unchanged)
        refit_i <- eval(bquote(
          glgpm(.(
            formula      = fit0$formula,
            data         = fit_data_sf[in_id, ],
            cov_offset   = .(fit0$cov_offset),
            family       = .(fam),
            crs          = .(fit0$crs),
            scale_to_km  = .(fit0$scale_to_km),
            control_mcmc = control_sim,
            fix_var_me   = .(fit0$fix_var_me),
            den          = .(as.name(den_name)),
            messages     = FALSE,
            start_pars   = par_hat
          ))
        ))
      } else {
        ## quick slice without re-fitting
        refit_i <- fit0
        keep <- in_id
        refit_i$data_sf  <- refit_i$data_sf [keep, ]
        refit_i$units_m  <- refit_i$units_m[keep]
        keep_coord <- unique(refit_i$ID_coords[keep])
        refit_i$coords   <- refit_i$coords[keep_coord, , drop = FALSE]
        refit_i$y        <- refit_i$y      [keep]
        refit_i$D        <- refit_i$D      [keep, , drop = FALSE]
        if (!is.null(refit_i$cov_offset))
          refit_i$cov_offset <- refit_i$cov_offset[keep]
        if (!is.null(refit_i$ID_re)) {
          re_terms <- names(refit_i$ID_re)
          random_effects_i <- prepare_random_effects(refit_i$data_sf, re_terms)
          refit_i$ID_re <- as.data.frame(random_effects_i$ID_re)
          colnames(refit_i$ID_re) <- random_effects_i$names_re
          refit_i$re <- random_effects_i$re_unique_f
        }
        ## recompute ID_coords mapping
        refit_i$ID_coords <- create_ids(refit_i$data_sf)$ID_coords
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
        control_sim     = control_sim,
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
    out$model[[model_names[h]]] <- list(score = list())
    if (get_CRPS)  out$model[[model_names[h]]]$score$CRPS  <- CRPS
    if (get_SCRPS) out$model[[model_names[h]]]$score$SCRPS <- SCRPS
    if (get_AnPIT) {
      out$model[[model_names[h]]]$score$AnPIT_area <- AnPIT_area
      if (fam == "gaussian") out$model[[model_names[h]]]$PIT <- PIT else out$model[[model_names[h]]]$AnPIT <- AnPIT
    }

  } # end h loop

  class(out) <- "RiskMap_cross_validation"
  return(out)
}


.anpit_area <- function(curve, u = seq(0, 1, length.out = length(curve))) {
  if (length(curve) != length(u)) stop("'curve' and 'u' must have the same length")
  if (length(curve) < 2) return(NA_real_)

  dx <- diff(u)
  y <- abs(curve - u)
  sum(dx * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}
