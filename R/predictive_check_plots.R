# Plots of the assessment of models against simulated data.

##' @title Plot the assessment of models against simulated data
##'
##' @description
##' Draws the observed data against datasets simulated from each model, from the
##' result of [assess_fit()], as `plot()` of the result does.
##'
##' @param x The result of [assess_fit()].
##' @param type What to draw, each `"statistic"` (the default), the
##'   `"distribution"` of the prevalence, the rate or the response, or the
##'   `"variogram"` of the residuals from the covariates.
##' @param model Names of the models to include, all of them by default.
##' @param colour The colour of what is observed. The histograms and the bands
##'   have paler versions of it.
##' @param ... For `plot()`, the arguments of `plot_predictive_check()`. For
##'   `plot_predictive_check()`, not used, and an error if given, so that a
##'   misspelt argument is not ignored. Add \pkg{ggplot2} components to the
##'   plot with `+`.
##'
##' @details
##' # Statistics
##'
##' The histograms are those of the simulated datasets, with the models in rows,
##' and the vertical line is the observed value. The darker bars are the simulated
##' values at or above the observed one, whose proportion
##' [summary()][summary.RiskMap_predictive_check()] gives. The observed value is
##' unusual for the model when it is in a tail of the histogram, beyond nearly all
##' of the bars. It doesn't have to be in the middle. A bar at the observed value
##' is darker, so the whole histogram is when the observed value is the smallest
##' possible one, such as a proportion of zeros of 0, which need not be unusual.
##'
##' # Distribution and variogram
##'
##' The band is the 95% joint envelope of the simulated curves and the dashed
##' line their median at each value. The observed distribution function is the
##' thick line, and the observed semivariances are the points. Where the
##' observed curve is outside the band, it is orange. The band is joint, so that
##' the whole curve of a dataset from the model is inside it with a probability
##' of 95%, however many values or distance classes there are. Its edge depends
##' on the simulations, so a single point near the edge says little, and the
##' first distance classes often have few pairs. The points of the variogram are
##' not joined by lines, as in [plot_variogram()]. Each model has its own
##' points, because the residuals are those from its covariates.
##'
##' @return A \pkg{ggplot2} object, to which components can be added.
##'
##' @seealso [assess_fit()]
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:100, ]
##' fit <- glgpm(y ~ log(pop_dens) + gp(), data = example_data,
##'              family = "gaussian", messages = FALSE)
##' result <- assess_fit(fit, nsim = 199, seed = 1)
##' plot_predictive_check(result)
##' plot_predictive_check(result, type = "distribution")
##' plot_predictive_check(result, type = "variogram")
##' @export
plot_predictive_check <- function(x,
                                  type = c("statistic", "distribution",
                                           "variogram"),
                                  model = NULL, colour = "#00857F", ...) {
  check_supplied(x)
  check_plot_object(x, "RiskMap_predictive_check", "plot_predictive_check",
                    "assess_fit")
  check_dots_unused(...)
  type <- match_choice(type)
  models <- names(x$family)
  check_assessed_model(model, models, "model", several = TRUE)
  check_colour(colour)
  colours <- validation_palette(colour)
  shown <- if (is.null(model)) models else models[models %in% model]
  switch(type,
         statistic = plot_check_statistics(x, shown, colours),
         distribution = plot_check_distribution(x, shown, colours),
         variogram = plot_check_variogram(x, shown, colours))
}

##' @rdname plot_predictive_check
##' @method plot RiskMap_predictive_check
##' @export
plot.RiskMap_predictive_check <- function(x, ...) plot_predictive_check(x, ...)

##' The theme of the panels of a predictive check
##' @noRd
predictive_check_theme <- function(colours) {
  validation_theme() +
    theme(panel.border = element_rect(colour = colours$frame, fill = NA,
                                      linewidth = 0.4))
}

##' The envelope of the simulated curves and their median
##' @noRd
envelope_layers <- function(colours) {
  list(geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper),
                   fill = colours$band),
       geom_path(aes(y = .data$median), colour = colours$reference,
                 linewidth = 0.4, linetype = "22"))
}

##' The colours of observed values inside and outside the envelope
##' @noRd
scale_colour_outside <- function(colours) {
  scale_colour_manual(values = c(`FALSE` = colours$line,
                                 `TRUE` = colours$contrast), guide = "none")
}

##' Rows of some models, in their order
##' @noRd
rows_of_models <- function(data, models) {
  data <- data[data$model %in% models, ]
  data$model <- factor(data$model, levels = models)
  data
}

##' Histograms of the simulated statistics with the observed ones
##' @noRd
plot_check_statistics <- function(x, models, colours) {
  summarised <- as.data.frame(unclass(summary(x)))
  summarised <- rows_of_models(summarised, models)
  statistics <- unique(summarised$statistic)
  labels <- check_statistic_labels(statistics, x)
  label <- function(data) {
    data$statistic <- factor(data$statistic, levels = statistics, labels = labels)
    data
  }
  simulated <- rows_of_models(x$simulated, models)
  simulated <- simulated[is.finite(simulated$value), ]
  if (!nrow(simulated)) {
    abort_input("{.arg x} has no finite statistics to draw.", arg = "x")
  }
  observed <- summarised[is.finite(summarised$observed), ]
  bars <- label(statistic_bars(simulated, observed))
  out <- ggplot(bars) +
    geom_rect(aes(xmin = .data$from, xmax = .data$to, ymin = 0,
                  ymax = .data$count, fill = .data$above), colour = "white",
              linewidth = 0.15) +
    scale_fill_manual(values = c(`FALSE` = paler_colour(colours$line, 88, 0.3),
                                 `TRUE` = paler_colour(colours$line, 74, 0.47)),
                      guide = "none") +
    geom_vline(data = label(observed), aes(xintercept = .data$observed),
               colour = colours$line, linewidth = 0.8) +
    scale_y_continuous(expand = expansion(c(0, 0.06))) +
    labs(x = NULL, y = "Simulated datasets") +
    predictive_check_theme(colours) +
    theme(panel.grid = element_blank())
  if (length(models) > 1L) {
    out + facet_grid(rows = vars(.data$model), cols = vars(.data$statistic),
                     scales = "free", labeller = label_wrap_gen(22L)) +
      theme(strip.text.y = element_text(angle = -90, hjust = 0.5,
                                        margin = margin(0, 2, 0, 4)))
  } else {
    out + facet_wrap(vars(.data$statistic), scales = "free",
                     labeller = label_wrap_gen(22L))
  }
}

##' The bars of the histograms of simulated statistics
##'
##' Each panel has its own bars, about 25 of them, and each bar holds values that
##' are all below the observed one or all at or above it, as in
##' `summarise_statistic()`. A statistic that takes equally spaced values, such as
##' a proportion of locations, has bars that each hold the same number of these
##' values.
##' @param simulated,observed Data frames with the `model` and the `statistic`,
##'   and with the `value` or the `observed` value.
##' @noRd
statistic_bars <- function(simulated, observed) {
  panel_of <- function(data) paste(data$model, data$statistic, sep = "\r")
  pieces <- lapply(split(simulated, panel_of(simulated)), function(panel) {
    value <- observed$observed[match(panel_of(panel[1L, ]), panel_of(observed))]
    above <- at_or_above(panel$value, value)
    bins <- histogram_breaks(panel$value, above, value)
    n <- length(bins$breaks)
    bin <- findInterval(panel$value, bins$breaks, rightmost.closed = TRUE)
    data.frame(model = panel$model[1L], statistic = panel$statistic[1L],
               from = bins$breaks[-n], to = bins$breaks[-1L],
               count = tabulate(bin, n - 1L), above = bins$above)
  })
  out <- do.call(rbind, unname(pieces))
  out[out$count > 0, ]
}

##' The limits of the bars of a histogram, and whether each bar holds the values
##' that are `above`
##'
##' @param values The values of the histogram.
##' @param above Whether each value is at or above the observed one.
##' @param observed The observed value.
##' @noRd
histogram_breaks <- function(values, above, observed, n_bars = 25L) {
  range <- range(values)
  distinct <- sort(unique(signif(values, 10L)))
  if (length(distinct) == 1L) {
    width <- if (range[1L] == 0) 1 else 0.1 * abs(range[1L])
    return(list(breaks = range[1L] + c(-0.5, 0.5) * width, above = any(above)))
  }
  step <- min(diff(distinct))
  multiples <- diff(distinct) / step
  spaced <- diff(range) / step < 1e3 &&
    all(abs(multiples - round(multiples)) < 1e-4)
  if (spaced) {
    # The bars go from half a step below a value.
    width <- step * ceiling((round(diff(range) / step) + 1) / n_bars)
    range <- range + c(-0.5, 0.5) * step
  } else {
    width <- diff(range) / n_bars
  }
  # A bar ends between the values below the observed one and the others, and at
  # the observed one if it isn't at a value.
  divide <- range[1L]
  if (any(above) && !all(above)) {
    between <- c(max(values[!above]), min(values[above]))
    divide <- if (spaced) {
      between[2L] - 0.5 * step
    } else if (all(abs(observed - between) > 1e-6 * width)) {
      observed
    } else {
      mean(between)
    }
  }
  start <- divide - width * ceiling((divide - range[1L]) / width - 1e-6)
  breaks <- start + width * (0:ceiling((range[2L] - start) / width - 1e-6))
  # The first and the last bar don't go beyond the values.
  breaks <- c(range[1L], breaks[-c(1L, length(breaks))], range[2L])
  middle <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  list(breaks = breaks, above = any(above) & (all(above) | middle > divide))
}

##' The distribution of the data over the joint envelope of the simulated ones
##' @noRd
plot_check_distribution <- function(x, models, colours) {
  steps <- distribution_steps(rows_of_models(x$distribution, models))
  label <- switch(x$quantity, prevalence = "Prevalence", rate = "Rate",
                  x$response)
  out <- ggplot(steps, aes(x = .data$value)) +
    envelope_layers(colours) +
    geom_path(aes(y = .data$observed, group = .data$model,
                  colour = .data$outside), linewidth = 0.7) +
    scale_colour_outside(colours) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = c("0", "0.25", "0.5", "0.75", "1"),
                       expand = expansion(0.02)) +
    labs(x = label, y = paste0("Proportion of ", x$unit, "s at or below")) +
    predictive_check_theme(colours)
  out <- if (x$quantity == "prevalence") {
    out + scale_x_continuous(labels = scales::label_percent(),
                             expand = expansion(0))
  } else {
    out + scale_x_continuous(expand = expansion(0))
  }
  if (length(models) > 1L) out <- out + facet_wrap(vars(.data$model))
  out
}

##' The steps of distribution functions and of their envelopes
##'
##' Each value holds up to the next one, so that the lines and the band are
##' drawn as steps.
##' @noRd
distribution_steps <- function(distribution) {
  pieces <- lapply(split(distribution, distribution$model, drop = TRUE),
                   function(rows) {
    n <- nrow(rows)
    if (n < 2L) return(rows)
    out <- rows[rep(seq_len(n - 1L), each = 2L), ]
    out$value <- rows$value[c(rbind(seq_len(n - 1L), seq_len(n)[-1L]))]
    # The rise to the next step is outside when both of its ends are.
    out$outside[2L * seq_len(n - 1L)] <- rows$outside[-n] & rows$outside[-1L]
    # The last step goes a little further.
    end <- rows[c(n, n), ]
    end$value[2L] <- rows$value[n] + 0.02 * (rows$value[n] - rows$value[1L])
    rbind(out, end)
  })
  do.call(rbind, unname(pieces))
}

##' The variogram of the data over the joint envelope of the simulated ones
##' @noRd
plot_check_variogram <- function(x, models, colours) {
  if (is.null(x$variogram)) {
    abort_input(c("{.arg x} has no variogram, because no pair of locations is in a distance class.", "i" = "Give the distance classes in {.arg breaks} of {.fn assess_fit}."), arg = "x")
  }
  variogram <- rows_of_models(x$variogram, models)
  if (!nrow(variogram)) {
    abort_input("The covariates of {.val {models}} leave no residuals, so there is no variogram to draw.", arg = "model")
  }
  # One distance class has an interval instead of a band.
  envelope <- if (length(unique(variogram$distance)) > 1L) {
    envelope_layers(colours)
  } else {
    geom_linerange(aes(ymin = .data$lower, ymax = .data$upper),
                   colour = colours$band, linewidth = 4)
  }
  out <- ggplot(variogram, aes(x = .data$distance)) +
    envelope +
    geom_point(aes(y = .data$observed, colour = .data$outside), size = 1.7) +
    scale_colour_outside(colours) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(expand = expansion(c(0, 0.03))) +
    scale_y_continuous(expand = expansion(c(0, 0.05))) +
    labs(x = paste0("Distance (", x$distance_units, ")"),
         y = "Semivariance of the residuals") +
    predictive_check_theme(colours)
  if (length(models) > 1L) out <- out + facet_wrap(vars(.data$model))
  out
}
