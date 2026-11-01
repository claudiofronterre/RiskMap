##' Plot profile likelihoods of a fitted geostatistical model
##'
##' Draws the profiles computed by [profile.RiskMap()]. Each marginal profile
##' shows the relative profile log-likelihood, the profile log-likelihood
##' minus its maximum, against the parameter. It is 0 at the maximum, which
##' is the estimate unless a point on a profile is higher, and negative
##' elsewhere. The panel also shows the quadratic approximation
##' implied by the standard error (dashed), the cut-off for the requested
##' level and the profile interval shaded, with the estimate, the profile
##' interval and the Wald interval in the subtitle, one per line. Each joint
##' profile shows the relative log-likelihood over the grid as filled regions
##' at the 50%, 80% and 95% levels of a chi-squared distribution on two
##' degrees of freedom, the estimate, and the ridges of the two marginal
##' profiles (the value of one parameter that maximises the likelihood for
##' each value of the other). With one joint profile and both its marginal
##' profiles, the panels are arranged as a corner plot. Variances and the
##' scale are shown on a logarithmic axis.
##'
##' With `scale = "deviance"`, the panels show the profile deviance instead,
##' twice the drop in the log-likelihood from its maximum, which is minus
##' twice the relative log-likelihood. The cut-off of a profile interval at
##' level \eqn{1 - \alpha} is then the quantile \eqn{\chi^2_{1, 1 - \alpha}},
##' and on the relative log-likelihood it is \eqn{-\chi^2_{1, 1 - \alpha}/2}
##' (about -1.92 at the 95% level).
##'
##' @param x A `RiskMap_profile` object from [profile.RiskMap()].
##' @param parameters Names of the marginal profiles to draw. The default is
##'   all of them.
##' @param pairs Names of the joint profiles to draw, as `"first:second"`. The
##'   default is all of them.
##' @param scale The vertical scale of the profiles: `"loglik"` (the default)
##'   for the relative log-likelihood, or `"deviance"` for the profile
##'   deviance.
##' @param ... Not used.
##' @return Invisibly, a list of `ggplot` objects, one per panel. The combined
##'   figure is drawn.
##' @seealso [profile.RiskMap()]
##' @export
plot.RiskMap_profile <- function(x, parameters = NULL, pairs = NULL,
                                 scale = c("loglik", "deviance"), ...) {
  stopifnot(inherits(x, "RiskMap_profile"))
  scale <- match.arg(scale)
  parameters <- parameters %||% names(x$marginal)
  pairs <- pairs %||% names(x$joint)
  stopifnot("unknown marginal profile" = all(parameters %in% names(x$marginal)),
            "unknown joint profile" = all(pairs %in% names(x$joint)))
  axis <- profile_axis(scale, x$level)
  marginal <- lapply(parameters, function(name) {
    plot_marginal_profile(x, name, axis)
  })
  names(marginal) <- parameters
  joint <- lapply(pairs, function(pair) plot_joint_profile(x, pair, axis))
  names(joint) <- pairs
  arrange_profile_panels(x, marginal, joint, axis)
  invisible(c(marginal, joint))
}

##' The vertical scale of the profiles
##'
##' `value()` maps a deviance to the scale shown, and `limits` are the limits
##' of the marginal panels on that scale. The relative log-likelihood is
##' minus half the deviance, so its limits are those of the deviance halved
##' and negated.
##' @noRd
profile_axis <- function(scale, level) {
  top <- max(stats::qchisq(level, 1) * 1.8, 6)
  if (identical(scale, "deviance")) {
    list(scale = scale, value = function(deviance) deviance,
         limits = c(0, top), title = "Profile deviance",
         joint = "Profile deviance regions and profile ridges")
  } else {
    list(scale = scale, value = function(deviance) -deviance / 2,
         limits = c(-top / 2, 0), title = "Relative log-likelihood",
         joint = "Relative log-likelihood regions and profile ridges")
  }
}

##' Labels of a logarithmic axis as plain numbers, such as 0.00001 and 1,000
##' rather than 1e-05 and 1e+03
##' @noRd
plain_number_labels <- function() {
  scales::label_number(drop0trailing = TRUE, big.mark = ",")
}

##' Colours and theme shared by the profile panels
##' @noRd
profile_palette <- function() {
  list(profile = "#0072B2", wald = "#7a7f87", band = "#cfe3f2",
       cutoff = "#D55E00", regions = c("#08519c", "#3f87c2", "#9ecae1", "#f4f8fb"),
       ridge = "#1d232b")
}

##' @noRd
profile_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", size = 11),
                   plot.subtitle = ggplot2::element_text(colour = "#586271", size = 9))
}

##' Values of a parameter on the scale shown to the user
##' @noRd
profile_natural <- function(name, working) {
  if (grepl("^beta\\.", name)) working else exp(working)
}

##' One marginal profile, on the vertical scale of `axis` from
##' `profile_axis()`
##' @noRd
plot_marginal_profile <- function(x, name, axis = profile_axis("loglik", x$level)) {
  colours <- profile_palette()
  profile <- x$marginal[[name]]
  profile <- profile[is.finite(profile$deviance), ]
  estimate <- x$estimate[[name]]
  se <- x$standard_errors[[name]]
  grid <- seq(min(profile$value), max(profile$value), length.out = 200)
  # A fit without standard errors (a scale at zero, a singular information)
  # has no Wald approximation to draw.
  wald <- if (is.finite(se) && se > 0) {
    data.frame(value = profile_natural(name, grid),
               height = axis$value(((grid - estimate) / se)^2))
  } else {
    data.frame(value = numeric(), height = numeric())
  }
  shown <- data.frame(value = profile_natural(name, profile$value),
                      height = axis$value(profile$deviance))
  interval <- x$intervals[x$intervals$parameter == name, ]
  cutoff <- axis$value(stats::qchisq(x$level, 1))
  # An interval left open by the grid is shaded to the edge of the grid.
  plot <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = interval$profile_lower %|NA|% min(shown$value),
                      xmax = interval$profile_upper %|NA|% max(shown$value),
                      ymin = -Inf, ymax = Inf, fill = colours$band, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = cutoff, colour = colours$cutoff, linetype = "dotted") +
    ggplot2::geom_line(data = wald, ggplot2::aes(.data$value, .data$height),
                       colour = colours$wald, linetype = "dashed") +
    ggplot2::geom_line(data = shown, ggplot2::aes(.data$value, .data$height),
                       colour = colours$profile, linewidth = 0.9) +
    ggplot2::geom_point(data = shown, ggplot2::aes(.data$value, .data$height),
                        colour = colours$profile, size = 1.1) +
    ggplot2::geom_vline(xintercept = interval$estimate, colour = colours$ridge, linewidth = 0.3) +
    ggplot2::coord_cartesian(ylim = axis$limits) +
    ggplot2::labs(title = x$labels[[name]], x = NULL, y = axis$title,
                  subtitle = profile_interval_text(interval, x$level)) +
    profile_theme()
  if (!grepl("^beta\\.", name)) {
    plot <- plot + ggplot2::scale_x_log10(labels = plain_number_labels())
  }
  plot
}

##' Interval text for a panel subtitle
##'
##' The estimate and each interval go on a line of their own, so that the
##' subtitle fits a panel a third of a page wide.
##' @noRd
profile_interval_text <- function(interval, level) {
  format_value <- function(value) if (is.na(value)) "open" else format(signif(value, 3))
  sprintf("Estimate %s\n%d%% profile interval %s to %s\n%d%% Wald interval %s to %s",
          format_value(interval$estimate), round(100 * level),
          format_value(interval$profile_lower), format_value(interval$profile_upper),
          round(100 * level), format_value(interval$wald_lower),
          format_value(interval$wald_upper))
}

##' Replace a missing value
##' @noRd
`%|NA|%` <- function(value, replacement) if (is.na(value)) replacement else value

##' One joint profile, as filled confidence regions with the marginal ridges,
##' on the vertical scale of `axis` from `profile_axis()`
##'
##' The regions are those within the chi-squared cut-offs on two degrees of
##' freedom, which are the same on either scale. The bins of the filled
##' contours run from low to high values, so the innermost region is the
##' first bin of the deviance and the last of the log-likelihood.
##' @noRd
plot_joint_profile <- function(x, pair_name, axis = profile_axis("loglik", x$level)) {
  colours <- profile_palette()
  joint <- x$joint[[pair_name]]
  pair <- attr(joint, "pair")
  joint <- joint[is.finite(joint$deviance), ]
  levels <- c(0.5, 0.8, 0.95)
  cutoffs <- sort(axis$value(stats::qchisq(levels, 2)))
  breaks <- sort(axis$value(c(0, stats::qchisq(levels, 2),
                              max(joint$deviance, stats::qchisq(0.99, 2)) + 1)))
  regions <- c("50%", "80%", "95%", "outside")
  bins <- if (identical(axis$scale, "deviance")) regions else rev(regions)
  surface <- data.frame(first = profile_natural(pair[1], joint$first),
                        second = profile_natural(pair[2], joint$second),
                        height = axis$value(joint$deviance))
  ridges <- profile_ridges(x, pair)
  estimate <- data.frame(first = profile_natural(pair[1], x$estimate[[pair[1]]]),
                         second = profile_natural(pair[2], x$estimate[[pair[2]]]))
  plot <- ggplot2::ggplot(surface, ggplot2::aes(.data$first, .data$second)) +
    ggplot2::geom_contour_filled(
      ggplot2::aes(z = .data$height,
                   fill = ggplot2::after_stat(factor(bins[as.integer(.data$level)],
                                                     levels = regions))),
      breaks = breaks) +
    ggplot2::geom_contour(ggplot2::aes(z = .data$height), breaks = cutoffs,
                          colour = "white", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = stats::setNames(colours$regions, regions),
                               name = "Joint region", drop = FALSE,
                               guide = ggplot2::guide_legend(order = 2)) +
    ggplot2::geom_point(data = estimate, colour = "white", fill = colours$cutoff,
                        shape = 21, size = 2.6, stroke = 0.8) +
    ggplot2::labs(x = x$labels[[pair[1]]], y = x$labels[[pair[2]]],
                  title = "Joint profile", subtitle = axis$joint) +
    profile_theme() +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
  if (!is.null(ridges)) {
    plot <- plot + ggplot2::geom_path(data = ridges, ggplot2::aes(linetype = .data$ridge),
                                      colour = colours$ridge, linewidth = 0.5) +
      # Above the legend of the regions. The axes name the parameters, so
      # the entries refer to them as x and y.
      ggplot2::scale_linetype_manual(values = c("solid", "longdash"),
                                     name = "Profile ridges",
                                     guide = ggplot2::guide_legend(order = 1))
  }
  if (!grepl("^beta\\.", pair[1])) {
    plot <- plot + ggplot2::scale_x_log10(labels = plain_number_labels())
  }
  if (!grepl("^beta\\.", pair[2])) {
    plot <- plot + ggplot2::scale_y_log10(labels = plain_number_labels())
  }
  plot + ggplot2::coord_cartesian(xlim = range(surface$first), ylim = range(surface$second), expand = FALSE)
}

##' Ridges of the marginal profiles of a pair, on the scale shown to the user
##' @noRd
profile_ridges <- function(x, pair) {
  if (!all(pair %in% names(x$marginal))) return(NULL)
  path <- function(fixed, other) {
    profile <- x$marginal[[fixed]]
    profile <- profile[is.finite(profile$loglik), ]
    values <- if (fixed == pair[1]) {
      data.frame(first = profile$value, second = profile[[other]])
    } else {
      data.frame(first = profile[[other]], second = profile$value)
    }
    values$first <- profile_natural(pair[1], values$first)
    values$second <- profile_natural(pair[2], values$second)
    values$ridge <- if (fixed == pair[1]) "best y for each x" else "best x for each y"
    values
  }
  ridges <- rbind(path(pair[1], pair[2]), path(pair[2], pair[1]))
  ridges$ridge <- factor(ridges$ridge,
                         levels = c("best y for each x", "best x for each y"))
  ridges
}

##' Arrange the panels: a corner plot for one joint profile with both
##' marginals, otherwise a grid of panels
##' @noRd
arrange_profile_panels <- function(x, marginal, joint,
                                   axis = profile_axis("loglik", x$level)) {
  panels <- c(marginal, joint)
  if (length(joint) == 1L) {
    pair <- attr(x$joint[[names(joint)]], "pair")
    if (all(pair %in% names(marginal)) && requireNamespace("patchwork", quietly = TRUE)) {
      others <- marginal[setdiff(names(marginal), pair)]
      # The marginal panels share the joint panel's ranges so the axes line up.
      grid <- x$joint[[names(joint)]]
      first_range <- range(profile_natural(pair[1], grid$first))
      second_range <- range(profile_natural(pair[2], grid$second))
      first_panel <- suppressMessages(marginal[[pair[1]]] +
        ggplot2::coord_cartesian(xlim = first_range, ylim = axis$limits, expand = FALSE))
      second_panel <- suppressMessages(marginal[[pair[2]]] +
        ggplot2::coord_flip(xlim = second_range, ylim = axis$limits, expand = FALSE))
      corner <- patchwork::wrap_plots(
        first_panel, patchwork::plot_spacer(), joint[[1]], second_panel,
        ncol = 2, widths = c(2, 1), heights = c(1, 2))
      if (length(others)) {
        print(patchwork::wrap_plots(others, nrow = 1) / corner + patchwork::plot_layout(heights = c(1, 3)))
      } else {
        print(corner)
      }
      return(invisible(NULL))
    }
  }
  columns <- min(3L, length(panels))
  gridExtra::grid.arrange(grobs = panels, ncol = columns)
  invisible(NULL)
}
