##' @title Predictive Targets over a Shapefile (grid-aggregated)
##'
##' @description
##' Computes predictive targets over polygon features using joint prediction
##' samples from \code{\link{setup_prediction}}. Targets can incorporate
##' covariates, offsets, optional unstructured random effects.
##'
##' @param object Output from \code{\link{setup_prediction}} (class \code{RiskMap_pred}),
##'   typically fitted with \code{type = "joint"} so that linear predictor samples are available.
##' @param shp An \pkg{sf} polygon object representing regions over which predictions are aggregated.
##' @param shp_target A function that aggregates grid-cell values within each polygon to a
##'   single regional value (default \code{mean}). Examples: \code{mean}, \code{sum},
##'   a custom weighted mean, etc.
##' @param weights Optional numeric vector of weights used inside \code{shp_target}.
##'   If supplied with \code{standardize_weights = TRUE}, weights are normalized within each region.
##' @param standardize_weights Logical; standardize \code{weights} within each region (\code{FALSE} by default).
##' @param col_names Name or column index in \code{shp} containing region identifiers to use in outputs.
##' @param include_covariates Logical; include fitted covariate effects in the linear predictor (default \code{TRUE}).
##' @param include_nugget Logical; include the nugget (unstructured measurement error) in the linear predictor (default \code{FALSE}).
##' @param include_cov_offset Logical; include any covariate offset term (default \code{FALSE}).
##' @param return_shp Logical; if \code{TRUE}, return the shapefile with appended summary columns
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
##' For each polygon in \code{shp}, grid-cell samples of the linear predictor are transformed with
##' \code{f_target}, optionally adjusted for covariates, offset, nugget and/or REs, and
##' then aggregated via \code{shp_target} (optionally weighted). The list \code{pd_summary} is applied
##' to each region's target samples to produce summary statistics.
##'
##' @return An object of class \code{RiskMap_predict_areal_target} with components:
##' \itemize{
##'   \item \code{target}: \code{data.frame} of region-level summaries (one row per region).
##'   \item \code{target_samples}: (optional) \code{list} with one element per region; each contains
##'         a \code{data.frame}/matrix of raw samples for each named target in \code{f_target},
##'         if \code{return_target_samples = TRUE}.
##'   \item \code{shp}: (optional) the input \code{sf} object with appended summary columns,
##'         included if \code{return_shp = TRUE}.
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
                                 shp,
                                 shp_target = mean,
                                 weights = NULL,
                                 standardize_weights = FALSE,
                                 col_names = NULL,
                                 include_covariates = TRUE,
                                 include_nugget = FALSE,
                                 include_cov_offset = FALSE,
                                 return_shp = TRUE,
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
    stop("To run predictions with a shape file, joint predictions must be used;
         rerun 'setup_prediction' and set 'type' = \"joint\"")
  }

  check_data(shp, "polygon")

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

  n_reg <- nrow(shp)
  if(is.null(col_names)) {
    shp$region <- paste("reg", 1:n_reg, sep = "")
    col_names <- "region"
    names_reg <- shp$region
  } else {
    names_reg <- shp[[col_names]]
    if(n_reg != length(names_reg)) {
      stop("The names in the column identified by 'col_names' do not
         provide a unique set of names, but there are duplicates")
    }
  }

  if(list_mode) {
    shp <- st_transform(shp, crs = st_crs(object$grid_pred[[1]])$input)
  } else {
    shp <- st_transform(shp, crs = st_crs(object$grid_pred)$input)
  }

  if(!list_mode) {
    inter <- st_intersects(shp, object$grid_pred)
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
      if(messages) message("Computing predictive target for: ", shp[[col_names]][h])
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

        target_samples_i <- apply(target_grid_samples_i, 2, function(x) shp_target(weights_h * x))

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
      if(messages) message("Computing predictive target for:", shp[[col_names]][h])
      if(length(inter[[h]]) == 0) {
        warning(paste("No points on the grid fall within", shp[[col_names]][h],
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

          target_samples_i <- apply(target_grid_samples_i, 2, function(x) shp_target(weights_h * x))

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

  if(return_shp) {
    if(length(no_comp) > 0) {
      ind_reg <- (1:n_reg)[-no_comp]
    } else {
      ind_reg <- 1:n_reg
    }
    for(i in 1:n_f) {
      for(j in 1:n_summaries) {
        name_ij <- paste(names_f[i], "_", paste(names_s[j]), sep = "")
        shp[[name_ij]] <- rep(NA, n_reg)
        for(h in ind_reg) {
          which_reg <- which(shp[[col_names]] == names_reg[h])
          shp[which_reg, ][[name_ij]] <-
            out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]][[ paste(names_s[j]) ]]
        }
      }
    }
  }

  out$shp <- shp
  out$f_target <- names(f_target)
  out$pd_summary <- names(pd_summary)
  out$grid_pred <- object$grid_pred
  class(out) <- "RiskMap_predict_areal_target"
  return(out)
}


##' Plot Method for RiskMap_predict_areal_target Objects
##'
##' Generates a plot of predictive target values or summaries over a shapefile.
##'
##' @param x An object of class 'RiskMap_predict_areal_target' containing computed targets,
##' summaries, and associated spatial data.
##' @param which_target Character indicating the target type to plot (e.g., "linear_target").
##' @param which_summary Character indicating the summary type to plot (e.g., "mean", "sd").
##' @param ... Additional arguments passed to 'scale_fill_distiller' in 'ggplot2'.
##' @return A \code{ggplot} object showing the plot of the specified predictive target or summary.
##' @details
##' This function plots the predictive target values or summaries over a shapefile.
##' It requires the 'ggplot2' package for plotting and 'sf' objects for spatial data.
##'
##' @seealso
##' \code{\link{predict_areal_target}}, \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_sf}},
##' \code{\link[ggplot2]{aes}}, \code{\link[ggplot2]{scale_fill_distiller}}
##'
##' @method plot RiskMap_predict_areal_target
##' @export
plot.RiskMap_predict_areal_target <- function(x, which_target = "linear_target",
                                              which_summary = "mean", ...) {
  col_shp_name <- paste(which_target,"_",which_summary,sep="")

  out <- ggplot(x$shp) +
    geom_sf(aes(fill = x$shp[[col_shp_name]])) +
    scale_fill_distiller(...)
  return(out)
}
