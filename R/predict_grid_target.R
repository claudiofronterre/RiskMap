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
##' Generates a plot of the predicted values or summaries over the regular spatial grid
##' from an object of class 'RiskMap_predict_grid_target'.
##'
##' @param x An object of class 'RiskMap_predict_grid_target'.
##' @param which_target Character string specifying which target prediction to plot.
##' @param which_summary Character string specifying which summary statistic to plot (e.g., "mean", "sd").
##' @param ... Additional arguments passed to the \code{\link[terra]{plot}} function of the \code{terra} package.
##' @return A \code{ggplot} object representing the specified prediction target or summary statistic over the spatial grid.
##' @details
##' This function requires the 'terra' package for spatial data manipulation and plotting.
##' It plots the values or summaries over a regular spatial grid, allowing for visual examination of spatial patterns.
##'
##' @seealso \code{\link{predict_grid_target}}
##'
##' @importFrom terra as.data.frame rast plot
##' @method plot RiskMap_predict_grid_target
##' @export
##'
##'
plot.RiskMap_predict_grid_target <- function(x, which_target = "linear_target", which_summary = "mean", ...) {
  t_data.frame <-
    terra::as.data.frame(cbind(st_coordinates(x$grid_pred),
                               x$target[[which_target]][[which_summary]]),
                         xy = TRUE)
  raster_out <- terra::rast(t_data.frame, crs = st_crs(x$grid_pred)$input)

  terra::plot(raster_out, ...)
}
