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
##'   `object$test_set` (list of sf objects) and `object$model[[model]]$score[[score]]`.
##' @param score A string specifying the score to visualize. Must match a score computed in the model.
##' @param model A string specifying the model whose scores to visualize.
##' @param ... Additional arguments to customize ggplot, such as `scale_color_gradient` or `scale_color_manual`.
##' @return A ggplot object visualizing the spatial distribution of the specified score.
##' @export
plot_score <- function(object, score, model, ...) {

  # Check if "score" exists
  if (!score %in% names(object$model[[model]]$score)) {
    stop(paste("Error: The score", shQuote(score), "was not computed for model", shQuote(model)))
  }

  # Extract the test sets and number of test sets
  test_sets <- object$test_set
  n_test <- length(test_sets)

  # Combine the data and add the score variable
  data_full <- st_as_sf(test_sets[[1]])
  data_full$score <- object$model[[model]]$score[[score]][[1]]

  if (n_test > 1) {
    for (i in 1:n_test) {
      test_sets[[i]]$score <- object$model[[model]]$score[[score]][[i]]
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
    ggtitle(paste("Visualizing", score, "for model", model)) +
    theme_minimal()


  return(out)
}


##' @title Plotting the empirical variogram
##' @description Plots the empirical variogram generated by \code{\link{variogram}}
##' @param x The output generated by the function \code{\link{variogram}}.
##' @param plot_envelope A logical value indicating if the envelope of spatial independence
##' generated using the permutation test must be displayed (\code{plot_envelope = TRUE}) or not
##' (\code{plot_envelope = FALSE}). By default \code{plot_envelope = FALSE}. Note: if \code{n_permutations = 0} when
##' running the function \code{\link{variogram}}, the function will display an error message because no envelope can be generated.
##' @param color If \code{plot_envelope = TRUE}, it sets the colour of the envelope; run \code{vignette("ggplot2-specs")} for more details on this argument.
##' @param ... Further arguments passed to or from other methods (currently unused).
##' @return A \code{ggplot} object representing the empirical variogram plot, optionally including the envelope of spatial independence.
##' @details This function plots the empirical variogram, which shows the spatial dependence structure of the data. If \code{plot_envelope} is set to \code{TRUE}, the plot will also include an envelope indicating the range of values under spatial independence, based on a permutation test.
##' @seealso \code{\link{variogram}}
##' @method plot RiskMap_variogram
##' @export
plot.RiskMap_variogram <- function(x,
                                   plot_envelope = FALSE,
                                   color = "royalblue1",
                                   ...) {
  variogram_output <- x

  basic_plot <- ggplot(data = variogram_output$variogram,
                       aes(x = .data$mid_points, y = .data$obs_vari)) +
    geom_point() +
    geom_line()

  if (plot_envelope) {
    if (variogram_output$n_permutations == 0) {
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
##' @param x The output object from `surf_sim`, containing both simulated data (`data_sim`) and predicted grid simulations (`lp_grid_sim`).
##' @param sim The simulation index to plot.
##' @param ... Additional graphical parameters to be passed to the plotting function of the `terra` package.
##'
##' @return A plot of the simulation results.
##'
##' @importFrom stars st_rasterize
##' @method plot RiskMap_sim
##' @export
plot.RiskMap_sim <- function(x, sim, ...) {

  sf_object <- x$lp_grid_sim
  value_column <- paste0("lp_sim_",sim)
  r <- rast(st_rasterize(sf_object[,c("x",value_column)]))
  r[r == 0] <- NA

  plot(r, main = paste("Simulation no.", sim), ...)
  points(st_coordinates(x$data_sim[[sim]]), pch = 20)

}

##' Plot Method for RiskMap_pred_target_grid Objects
##'
##' Generates a plot of the predicted values or summaries over the regular spatial grid
##' from an object of class 'RiskMap_pred_target_grid'.
##'
##' @param x An object of class 'RiskMap_pred_target_grid'.
##' @param target Character string specifying which target prediction to plot.
##' @param summary Character string specifying which summary statistic to plot (e.g., "mean", "sd").
##' @param ... Additional arguments passed to the \code{\link[terra]{plot}} function of the \code{terra} package.
##' @return A \code{ggplot} object representing the specified prediction target or summary statistic over the spatial grid.
##' @details
##' This function requires the 'terra' package for spatial data manipulation and plotting.
##' It plots the values or summaries over a regular spatial grid, allowing for visual examination of spatial patterns.
##'
##' @seealso \code{\link{pred_target_grid}}
##'
##' @importFrom terra as.data.frame rast plot
##' @method plot RiskMap_pred_target_grid
##' @export
##'
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterre@@lancaster.ac.uk}
plot.RiskMap_pred_target_grid <- function(x, target = "linear_target", summary = "mean", ...) {
  t_data.frame <-
    terra::as.data.frame(cbind(st_coordinates(x$grid_pred),
                               x$target[[target]][[summary]]),
                         xy = TRUE)
  raster_out <- terra::rast(t_data.frame, crs = st_crs(x$grid_pred)$input)

  terra::plot(raster_out, ...)
}
