##' @title Plot Calibration Curves (AnPIT / PIT) from Spatial Cross-Validation
##'
##' @description
##' Produce calibration plots from a \code{RiskMap.spatial.cv} object returned by
##' \code{\link{assess_pp}}.
##' * For Binomial or Poisson models the function visualises the
##'   \emph{Aggregated normalised Probability Integral Transform} (AnPIT)
##'   curves stored in \code{$AnPIT}.
##' * For Gaussian models it detects the list \code{$PIT} and instead plots
##'   the empirical \emph{Probability Integral Transform} curve
##'   (ECDF of PIT values) on the same \eqn{u}-grid.
##'
##' A 45° dashed red line indicates perfect calibration.
##'
##' @param object       A \code{RiskMap.spatial.cv} object.
##' @param mode         One of \code{"average"} (average curve across test sets),
##'                     \code{"single"} (a specific test set),
##'                     or \code{"all"} (every test set separately).
##' @param test_set     Integer; required when \code{mode = "single"}.
##' @param model_name   Optional character string; if supplied,
##'                     only that model is plotted.
##' @param combine_panels Logical; when \code{mode = "average"}, draw
##'                       all models in a single panel (\code{TRUE})
##'                       or one panel per model (\code{FALSE}, default).
##'
##' @return A \pkg{ggplot2} object (single plot) or a \pkg{grid} object
##'   from \pkg{gridExtra} (multiple panels).
##'
##' @importFrom dplyr   group_by summarize %>%
##' @export
plot_AnPIT <- function(object,
                       mode = "average",
                       test_set = NULL,
                       model_name = NULL,
                       combine_panels = FALSE) {

  if (!inherits(object, "RiskMap.spatial.cv"))
    stop("`object` must be a 'RiskMap.spatial.cv' produced by assess_pp().")

  all_models <- names(object$model)

  if (!is.null(model_name)) {
    if (!model_name %in% all_models)
      stop("Model name '", model_name, "' not found in `object$model`.")
    all_models <- model_name
  }

  make_df <- function(mname) {
    m <- object$model[[mname]]
    if (!is.null(m$AnPIT)) {
      lapply(seq_along(m$AnPIT), function(j) {
        curve_vals <- m$AnPIT[[j]]
        if (length(curve_vals) == 0) return(NULL)
        data.frame(
          u_val   = seq(0, 1, length.out = length(curve_vals)),
          value   = curve_vals,
          test_set = j,
          model    = mname,
          type     = "AnPIT"
        )
      })
    } else if (!is.null(m$PIT)) {
      u_grid <- seq(0, 1, length.out = 1000)
      lapply(seq_along(m$PIT), function(j) {
        pit_vec <- m$PIT[[j]]
        if (length(pit_vec) == 0) return(NULL)
        data.frame(
          u_val   = u_grid,
          value   = ecdf(pit_vec)(u_grid),
          test_set = j,
          model    = mname,
          type     = "PIT"
        )
      })
    } else {
      NULL
    }
  }

  plot_data <- do.call(rbind, unlist(lapply(all_models, make_df), recursive = FALSE))

  if (is.null(plot_data) || nrow(plot_data) == 0)
    stop("No AnPIT or PIT data available for plotting.")

  y_label <- unique(plot_data$type)
  if (length(y_label) > 1) y_label <- "Calibration curve"

  id_line <- geom_abline(intercept = 0, slope = 1,
                         linetype = "dashed", colour = "red")

  if (mode == "average" && combine_panels) {
    avg <- plot_data %>%
      dplyr::group_by(.data$model, .data$u_val) %>%
      dplyr::summarize(value = mean(.data$value), .groups = "drop")

    return(
      ggplot(avg, aes(.data$u_val, .data$value, colour = .data$model)) +
        geom_line() + id_line +
        labs(title = "Average calibration curves",
             x = "", y = y_label) +
        theme_minimal() +
        guides(colour = guide_legend(title = "Model"))
    )
  }

  build_plot <- function(df, title_suffix = "") {
    ggplot(df, aes(.data$u_val, .data$value,
                   colour = if (mode == "all") as.factor(test_set) else NULL)) +
      geom_line() + id_line +
      labs(title = title_suffix, x = "", y = unique(df$type)) +
      theme_minimal() +
      guides(colour = guide_legend(title = "Test set"))
  }

  plots <- list()
  for (mname in all_models) {
    df_model <- dplyr::filter(plot_data, .data$model == mname)

    p <- switch(mode,
                average = {
                  avg <- df_model %>%
                    dplyr::group_by(.data$u_val) %>%
                    dplyr::summarize(value = mean(.data$value), .groups = "drop")
                  avg$type <- unique(df_model$type)
                  build_plot(avg, paste("Model", mname, ": average"))
                },
                single  = {
                  if (is.null(test_set))
                    stop("Provide `test_set` when mode = 'single'.")
                  df_ts <- dplyr::filter(df_model, test_set == test_set)
                  if (nrow(df_ts) == 0)
                    stop("No data for test_set ", test_set, " in model ", mname)
                  build_plot(df_ts,
                             paste("Model", mname, "- test set", test_set))
                },
                all     = build_plot(df_model,
                                     paste("Model", mname, "- all test sets")),
                stop("Invalid `mode`. Use 'average', 'single' or 'all'.")
    )

    plots[[mname]] <- p
  }

  if (length(plots) == 1) {
    plots[[1]]
  } else {
    ncol <- ifelse(length(plots) == 2, 2, 2)
    nrow <- ceiling(length(plots) / ncol)
    do.call(gridExtra::grid.arrange, c(plots, ncol = ncol, nrow = nrow))
  }
}


##' @title Plot Spatial Scores for a Specific Model and Metric
##'
##' @description This function visualizes spatial scores for a specified model and metric.
##' It combines test set data, handles duplicate locations by averaging scores,
##' and creates a customizable map using ggplot2.
##'
##' @param object A list containing test sets and model scores. The structure should include
##'   `object$test_set` (list of sf objects) and `object$model[[which_model]]$score[[which_score]]`.
##' @param which_score A string specifying the score to visualize. Must match a score computed in the model.
##' @param which_model A string specifying the model whose scores to visualize.
##' @param ... Additional arguments to customize ggplot, such as `scale_color_gradient` or `scale_color_manual`.
##' @return A ggplot object visualizing the spatial distribution of the specified score.
##' @export
plot_score <- function(object, which_score, which_model, ...) {

  # Check if "which_score" exists
  if (!which_score %in% names(object$model[[which_model]]$score)) {
    stop(paste("Error: The score", shQuote(which_score), "was not computed for model", shQuote(which_model)))
  }

  # Extract the test sets and number of test sets
  test_sets <- object$test_set
  n_test <- length(test_sets)

  # Combine the data and add the score variable
  data_full <- st_as_sf(test_sets[[1]])
  data_full$score <- object$model[[which_model]]$score[[which_score]][[1]]

  if (n_test > 1) {
    for (i in 1:n_test) {
      test_sets[[i]]$score <- object$model[[which_model]]$score[[which_score]][[i]]
      data_full <- rbind(data_full, test_sets[[i]])
    }
  }

  # Check for duplicate locations and average the score
  data_full <- data_full %>%
    mutate(geom_id = st_as_text(.data$geometry)) %>%
    group_by(.data$geom_id) %>%
    summarize(score = mean(.data$score, na.rm = TRUE),
              geometry = first(.data$geometry), .groups = "drop") %>%
    st_as_sf()


  # Create the base plot
  out <- ggplot(data = data_full) +
    geom_sf(aes(color = .data$score), size = 2) +
    ggtitle(paste("Visualizing", which_score, "for model", which_model)) +
    theme_minimal()


  return(out)
}

##' @title Plot the estimated MDA impact function
##'
##' @description
##' Generate a plot of the estimated impact of mass drug administration (MDA)
##' on infection prevalence, based on a fitted decay-adjusted spatio-temporal (DAST) model.
##' The function simulates draws from the posterior distribution of model parameters,
##' propagates them through the MDA effect function, and produces uncertainty bands
##' around the estimated impact curve.
##'
##' @param object A fitted DAST model object, returned by \code{\link{dast}}.
##' @param mda_history Specification of the MDA schedule. This can be either:
##'   \itemize{
##'     \item A numeric vector of event times (integers starting at 0, e.g. \code{c(0,1,2,6)}),
##'     \item OR a 0/1 indicator vector on the yearly grid (e.g. \code{c(1,1,1,0,0,0,1)}),
##'     where position \code{i} corresponds to year \code{i-1}.
##'   }
##'   If omitted, the default is a single MDA at time 0.
##' @param n_sim Number of posterior draws used for uncertainty quantification (default: 1000).
##' @param x_min Minimum value for the x-axis (default: \code{1e-6}).
##' @param x_max Maximum value for the x-axis (default: \code{10}).
##' @param conf_level Confidence level for the pointwise uncertainty interval (default: 0.95).
##' @param lower_f Optional lower bound for the y-axis. If not provided, computed from the data.
##' @param upper_f Optional upper bound for the y-axis. If not provided, computed from the data.
##' @param mc_cores Number of CPU cores to use for parallel simulation. Default is 1 (serial).
##' @param parallel_backend Parallelisation backend to use. Options are \code{"none"} (default),
##'   \code{"fork"} (Unix-like systems), or \code{"psock"} (cross-platform).
##' @param ... Additional arguments (currently unused).
##'
##' @details
##' The time axis is assumed to start at 0 and increase in integer steps of 1 year.
##' The argument \code{mda_history} allows the user to specify when MDAs occurred either
##' by listing the years directly or by giving a binary indicator on the yearly grid.
##' The function then evaluates the cumulative relative reduction
##' \eqn{1 - \mathrm{effect}(t)} at a dense grid of time points between \code{x_min}
##' and \code{x_max}, using the fitted parameters from the supplied DAST model.
##'
##' @return
##' A \code{ggplot2} object showing the median estimated MDA impact function
##' and the pointwise uncertainty band at the chosen confidence level.
##'
##' @export
plot_mda <- function(object,
                     mda_history  = NULL,   # numeric event times (integers, starting at 0) OR 0/1 vector on yearly grid
                     n_sim        = 1000,
                     x_min        = 1e-6,
                     x_max        = 10,
                     conf_level   = 0.95,
                     lower_f      = NULL,
                     upper_f      = NULL,
                     mc_cores     = 1,
                     parallel_backend = c("none","fork","psock"),
                     ...) {

  parallel_backend <- match.arg(parallel_backend)

  # --- Time axis for evaluation ---
  stopifnot(is.numeric(x_min), is.numeric(x_max), x_max > x_min)
  survey_times <- seq(x_min, x_max, length.out = 200)
  n_t <- length(survey_times)

  # --- MDA schedule ---
  if (is.null(mda_history)) {
    mda_times <- 0
  } else if (is.numeric(mda_history) && all(mda_history %in% c(0,1))) {
    if (length(mda_history) == 0L) {
      mda_times <- numeric(0)
    } else {
      mda_times <- which(mda_history == 1) - 1
    }
  } else if (is.numeric(mda_history)) {
    mda_times <- sort(unique(as.numeric(mda_history)))
  } else {
    stop("`mda_history` must be numeric: either integer event times (0,1,2,...) or a 0/1 vector on that yearly grid.")
  }

  # --- Extract params ---
  par_hat   <- coef(object)
  n_par     <- length(object$estimate)
  power_val <- object$power_val

  if (is.null(par_hat$alpha)) {
    ind_dast    <- n_par
    par_dast    <- log(par_hat$gamma)
    alpha_fixed <- object$fix_alpha
    has_alpha   <- FALSE
  } else {
    ind_dast    <- (n_par - 1):n_par
    par_dast    <- c(log(par_hat$alpha / (1 - par_hat$alpha)),
                     log(par_hat$gamma))
    alpha_fixed <- NA_real_
    has_alpha   <- TRUE
  }

  Sigma_par       <- as.matrix(object$covariance[ind_dast, ind_dast])
  Sigma_par_sroot <- t(chol(Sigma_par))
  par_hat_sim <- t(vapply(
    X   = seq_len(n_sim),
    FUN = function(i) par_dast + Sigma_par_sroot %*% rnorm(length(ind_dast)),
    FUN.VALUE = numeric(length(ind_dast))
  ))

  alphas <- if (has_alpha) plogis(par_hat_sim[, 1]) else rep(alpha_fixed, n_sim)
  gammas <- if (has_alpha) exp(par_hat_sim[, 2]) else exp(par_hat_sim[, 1])

  # --- Simulate effects ---
  if (length(mda_times) == 0L) {
    effects_mat <- matrix(0, nrow = n_t, ncol = n_sim)
  } else {
    intervention_mat <- matrix(1, nrow = n_t, ncol = length(mda_times))

    one_sim <- function(j) {
      eff <- compute_mda_effect(
        survey_times_data = survey_times,
        mda_times         = mda_times,
        intervention      = intervention_mat,
        alpha             = alphas[j],
        gamma             = gammas[j],
        kappa             = power_val
      )
      1 - eff
    }

    if (parallel_backend == "none" || mc_cores <= 1L) {
      eff_list <- lapply(seq_len(n_sim), one_sim)
    } else if (parallel_backend == "fork" && .Platform$OS.type == "unix") {
      mc_cores <- as.integer(max(1L, min(mc_cores, parallel::detectCores(logical = TRUE) - 1L, n_sim)))
      eff_list <- parallel::mclapply(seq_len(n_sim), one_sim,
                                     mc.cores = mc_cores, mc.preschedule = TRUE)
    } else if (parallel_backend == "psock") {
      mc_cores <- as.integer(max(1L, min(mc_cores, parallel::detectCores(logical = TRUE), n_sim)))
      cl <- parallel::makeCluster(mc_cores, type = "PSOCK")
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      parallel::clusterExport(cl,
                              varlist = c("survey_times","mda_times","intervention_mat","alphas","gammas","power_val","one_sim"),
                              envir = environment())
      eff_list <- parallel::parLapply(cl, seq_len(n_sim), one_sim)
    } else {
      eff_list <- lapply(seq_len(n_sim), one_sim)
    }

    effects_mat <- do.call(cbind, eff_list)
  }

  # --- Summaries ---
  alpha_q <- (1 - conf_level) / 2
  med   <- apply(effects_mat, 1, median,   na.rm = TRUE)
  lower <- apply(effects_mat, 1, quantile, probs = alpha_q, na.rm = TRUE)
  upper <- apply(effects_mat, 1, quantile, probs = 1 - alpha_q, na.rm = TRUE)

  in_view <- survey_times >= x_min & survey_times <= x_max
  if (!any(in_view)) in_view <- rep(TRUE, length(survey_times))

  if (is.null(lower_f)) lower_f <- min(lower[in_view], na.rm = TRUE)
  if (is.null(upper_f)) upper_f <- max(upper[in_view], na.rm = TRUE)

  plot_data <- data.frame(
    time   = survey_times,
    median = med,
    lower  = lower,
    upper  = upper
  )

  p <- ggplot(plot_data, aes(x = time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "grey70", alpha = 0.3) +
    geom_line(aes(y = median),
              color = "black", linewidth = 1) +
    labs(
      x = "Years since baseline",
      y = "Relative reduction from baseline prevalence",
      title = "MDA Impact Over Time"
    ) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(lower_f, upper_f)) +
    theme_minimal()

  # --- Add vertical dashed lines for MDA times ---
  if (length(mda_times) > 0) {
    p <- p + geom_vline(xintercept = mda_times,
                        linetype = "dashed", color = "red", alpha = 0.7)
  }

  return(p)
}

##' @title Plotting the empirical variogram
##' @description Plots the empirical variogram generated by \code{\link{variogram}}
##' @param variogram_output The output generated by the function \code{\link{variogram}}.
##' @param plot_envelope A logical value indicating if the envelope of spatial independence
##' generated using the permutation test must be displayed (\code{plot_envelope = TRUE}) or not
##' (\code{plot_envelope = FALSE}). By default \code{plot_envelope = FALSE}. Note: if \code{n_permutations = 0} when
##' running the function \code{\link{variogram}}, the function will display an error message because no envelope can be generated.
##' @param color If \code{plot_envelope = TRUE}, it sets the colour of the envelope; run \code{vignette("ggplot2-specs")} for more details on this argument.
##' @return A \code{ggplot} object representing the empirical variogram plot, optionally including the envelope of spatial independence.
##' @details This function plots the empirical variogram, which shows the spatial dependence structure of the data. If \code{plot_envelope} is set to \code{TRUE}, the plot will also include an envelope indicating the range of values under spatial independence, based on a permutation test.
##' @seealso \code{\link{variogram}}
##' @export
plot_variogram <- function(variogram_output,
                           plot_envelope = FALSE,
                           color = "royalblue1") {

  if (!inherits(variogram_output, "RiskMap_variogram")){
    stop("'variogram' must be an object of class 'RiskMap_variogram'")
  }

  basic_plot <- ggplot(data = variogram_output$variogram,
                       aes_string(x = "mid_points", y = "obs_vari")) +
    geom_point() +
    geom_line()

  if (plot_envelope) {
    if (variogram_output$n_permutations == 0){
      stop("To plot the envelope for spatial independence 'n_permutations' in `variogram()` must be greater than 0")
    }
    basic_plot <- basic_plot +
      geom_ribbon(aes(ymin = variogram_output$variogram$lower_bound,
                      ymax = variogram_output$variogram$upper_bound),
                  fill = color, alpha = 0.3)
  }

  if (variogram_output$scale_to_km) {
    x_label <- "Distance (km)"
  } else {
    x_label <- "Distance (m)"
  }

  basic_plot <- basic_plot + labs(x = x_label,
                                  y = "Variogram")

  return(basic_plot)
}


##' Plot simulated surface data for a given simulation
##'
##' This function plots the simulated surface data for a specific simulation from the result of `surf_sim`. It visualizes the linear predictor values on a raster grid along with the actual data points.
##'
##' @param surf_obj The output object from `surf_sim`, containing both simulated data (`data_sim`) and predicted grid simulations (`lp_grid_sim`).
##' @param sim The simulation index to plot.
##' @param ... Additional graphical parameters to be passed to the plotting function of the `terra` package.
##'
##' @return A plot of the simulation results.
##'
##' @importFrom stars st_rasterize
##'
##' @export
plot_sim_surf <-  function(surf_obj, sim, ...) {

  sf_object <- surf_obj$lp_grid_sim
  value_column <- paste0("lp_sim_",sim)
  r <- rast(st_rasterize(sf_object[,c("x",value_column)]))
  r[r == 0] <- NA

  plot(r, main = paste("Simulation no.", sim), ...)
  points(st_coordinates(surf_obj$data_sim[[sim]]), pch = 20)

}
