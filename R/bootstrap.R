##' Parametric bootstrap of a fitted geostatistical model
##'
##' Simulates datasets from a model fitted with [glgpm()] at its estimated
##' parameters, refits the model to each, and summarises the sampling
##' distribution of the estimates with standard errors, bias and two kinds of
##' bootstrap confidence interval.
##'
##' @param object A fitted `RiskMap` object returned by [glgpm()], from a
##'   Gaussian, binomial or Poisson model, with or without a nugget and
##'   [re()] terms.
##' @param n_sim Number of bootstrap replicates. Percentile intervals at the
##'   95% level need about 1000 for stable limits.
##' @param seed Optional non-negative integer that makes the bootstrap
##'   reproducible. When `NULL` (the default), one integer is drawn from the
##'   session's random number generator, so [set.seed()] before the call also
##'   reproduces it. The seed used is returned in `settings`. The caller's
##'   random number state is otherwise left unchanged.
##' @param level Probability covered by the intervals.
##' @param cores Number of processes. With more than one, the replicates run
##'   in parallel by forking, which is not available on Windows, where they
##'   run one after another.
##' @param messages Logical, whether to report progress.
##'
##' @details
##' Each replicate simulates a response at the original locations, with the
##' original covariates, denominators, offsets and random-effect groups, from
##' the fitted model (see [simulate_glgpm()]), and refits the model to it by
##' [glgpm()], starting from the estimate, without the final conditional
##' samples and quietly. The refits reuse the design matrix, the locations
##' and the settings of `control_mcml` stored in the fit, so the model
##' formula and the data are not evaluated again.
##'
##' Replicate \eqn{i} uses one seed, \eqn{s_i = s + 10000 (i - 1)}, where
##' \eqn{s} is `seed`. It simulates its data with `simulate_glgpm(seed =` \eqn{s_i}`)`
##' and refits with `glgpm(seed =` \eqn{s_i}`)`. The Monte Carlo streams of a
##' binomial or Poisson fit are derived from its seed by offsets between 1 and
##' about 6000, so consecutive integers would make replicates share streams,
##' and seeds 10000 apart do not. A replicate depends only on its seed, so the
##' results are identical for any number of `cores` and the first replicates
##' of a larger bootstrap are those of a smaller one. Each replicate runs with
##' one BLAS and OpenMP thread when the RhpcBLASctl package is installed, and
##' the user's thread settings are restored afterwards.
##'
##' A refit that fails is recorded in `failures` with its error message, and
##' its estimates are `NA`. Warnings from the refits, for example about
##' variances estimated at zero, are recorded in `warnings`.
##'
##' `summary()` reports, for each parameter, the estimate, the median of the
##' bootstrap estimates, the bias (their mean minus the estimate), their
##' standard deviation as the standard error, the proportion below the
##' estimate, and two intervals. The percentile interval takes the quantiles
##' of the bootstrap estimates at \eqn{\alpha/2} and \eqn{1 - \alpha/2}, where
##' \eqn{\alpha = 1 -} `level`. The bias-corrected (BC) interval takes them at
##' \eqn{\Phi(2 z_0 + z_{\alpha/2})} and \eqn{\Phi(2 z_0 + z_{1 - \alpha/2})},
##' where \eqn{z_0 = \Phi^{-1}(p)} and \eqn{p} is the proportion of bootstrap
##' estimates below the estimate (Efron and Tibshirani, 1993, Chapter 14).
##' The BC interval is not available when \eqn{p} is 0 or 1.
##'
##' The difference matters for the spatial variance \eqn{\sigma^2} and the
##' scale \eqn{\phi}. Their maximum likelihood estimators are biased
##' downward, so most bootstrap estimates fall below the estimate, just as
##' the estimate tends to fall below the true value. The percentile interval
##' then shifts downward, the wrong way, because it reproduces the bias of
##' the estimator instead of correcting it. With \eqn{p > 0.5},
##' \eqn{z_0 > 0} and the BC interval moves both limits upward. Both intervals
##' are unchanged by monotone transformations, so they are the same whether
##' computed for \eqn{\phi} or \eqn{\log \phi}.
##'
##' @return An object of class `RiskMap_bootstrap`, a list with
##' \item{estimates}{A matrix of the bootstrap estimates on their natural
##'   scale, with one row per replicate (`NA` for failed refits) and one
##'   column per parameter, named as `unlist(coef(object))`.}
##' \item{estimate}{The estimate of `object`, on the same scale.}
##' \item{labels}{Readable names of the parameters.}
##' \item{failures}{A data frame with the `replicate`, `seed` and error
##'   `message` of each failed refit.}
##' \item{warnings}{A data frame with the `replicate`, `seed` and `message`
##'   of each warning from the refits.}
##' \item{settings}{A list with `n_sim`, `seed`, `level`, `cores`, the
##'   `seed_spacing` between replicates, the `family` and the
##'   `distance_units` of the model, and the `elapsed_seconds`.}
##'
##' @references
##' Efron, B. and Tibshirani, R. J. (1993). *An Introduction to the
##' Bootstrap*. Chapman & Hall, New York.
##'
##' @seealso [simulate_glgpm()], [profile.RiskMap()],
##'   [plot.RiskMap_bootstrap()]
##'
##' @examples
##' data(italy_sim)
##' fit <- glgpm(y ~ gp(), data = italy_sim[1:50, ], family = "gaussian",
##'              messages = FALSE)
##' bootstrap <- parametric_bootstrap(fit, n_sim = 20, seed = 1,
##'                                   messages = FALSE)
##' bootstrap
##'
##' \donttest{
##' data(liberia)
##' liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
##' fit <- glgpm(npos ~ log(elevation) + gp(), denominator = ntest,
##'              data = liberia, model_crs = 32629, family = "binomial",
##'              seed = 2026)
##' bootstrap <- parametric_bootstrap(fit, n_sim = 100, seed = 2026)
##' summary(bootstrap)
##' plot(bootstrap)
##' }
##' @export
parametric_bootstrap <- function(object, n_sim = 1000, seed = NULL,
                                 level = 0.95, cores = 1, messages = TRUE) {
  if (!inherits(object, "RiskMap")) {
    stop("'object' must be a 'RiskMap' object returned by 'glgpm'",
         call. = FALSE)
  }
  check_positive_integer(n_sim, "n_sim")
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)
  check_range(level, min = 0, max = 1, allow_equal = FALSE, name = "level")
  check_positive_integer(cores, "cores")
  check_logical(messages)
  control_mcml <- attr(object, "control_mcml") %||% set_control_mcml()
  seed_spacing <- bootstrap_seed_spacing(control_mcml$max_stages)
  if (n_sim > .Machine$integer.max %/% seed_spacing) {
    stop("'n_sim' must be at most ", .Machine$integer.max %/% seed_spacing,
         call. = FALSE)
  }
  cores <- as.integer(cores)
  if (cores > 1L && .Platform$OS.type == "windows") {
    if (messages) message("Parallel replicates need forking, which Windows ",
                          "does not support, so they run on one core.")
    cores <- 1L
  }

  # Without a seed, draw one from the session's stream, which advances it
  # like any random call. Nothing else changes the caller's state.
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  seed <- as.integer(seed)
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)

  # The workers make no call to sf, because GDAL and PROJ can deadlock in a
  # forked process, so everything that needs it is prepared here.
  replicate_seeds <- bootstrap_seeds(seed, n_sim, seed_spacing)
  prepared <- prepare_simulation(simulation_model(object), "data")
  refit <- make_bootstrap_refit(object)
  run_replicate <- function(i) {
    run_bootstrap_replicate(prepared, refit, replicate_seeds[i])
  }

  started <- proc.time()[["elapsed"]]
  if (messages) {
    message("Simulating and refitting ", n_sim, " replicates",
            if (cores > 1L) paste0(" on ", cores, " cores"), "...")
  }
  results <- run_bootstrap_replicates(seq_len(n_sim), run_replicate, cores,
                                      messages)
  elapsed <- proc.time()[["elapsed"]] - started

  out <- assemble_bootstrap(object, results, replicate_seeds)
  out$settings <- list(n_sim = as.integer(n_sim), seed = seed, level = level,
                       cores = cores, seed_spacing = seed_spacing,
                       family = object$family,
                       distance_units = object$distance_units,
                       elapsed_seconds = elapsed)
  n_failed <- nrow(out$failures)
  if (n_failed == n_sim) {
    stop("All ", n_sim, " refits failed. The first failure: ",
         out$failures$message[1L], call. = FALSE)
  }
  if (n_failed > 0L) {
    warning(n_failed, " of ", n_sim, " refits failed and are left out. ",
            "See the 'failures' element of the result.", call. = FALSE)
  }
  if (messages) {
    message("Done in ", format_bootstrap_time(elapsed), ". ",
            n_sim - n_failed, " of ", n_sim, " refits succeeded.")
  }
  out
}

##' Distance between the seeds of consecutive bootstrap replicates
##'
##' Larger than the largest offset `fit_adaptive_mcml()` and the final
##' sampling add to the seed of a fit, so no two replicates share a Monte
##' Carlo stream. The offsets reach 2000 plus three per stage, 6000 plus one
##' per stage, and three per stage plus one, all below 6000 plus three per
##' stage.
##' @noRd
bootstrap_seed_spacing <- function(max_stages) {
  spacing <- 10000L
  if (6000 + 3 * max_stages >= spacing) {
    stop("parametric_bootstrap() supports fits with at most ",
         (spacing - 6001) %/% 3, " stages of certified MCML ('max_stages' ",
         "in set_control_mcml())", call. = FALSE)
  }
  spacing
}

##' Seeds of the bootstrap replicates, wrapped like `fit_adaptive_mcml()`
##' wraps its streams
##' @noRd
bootstrap_seeds <- function(seed, n_sim, spacing) {
  offsets <- (seq_len(n_sim) - 1) * spacing
  as.integer(((seed - 1 + offsets) %% .Machine$integer.max) + 1)
}

##' A function refitting a model to a new response
##'
##' The refits use what the fit stores: its design matrix, unique coordinates
##' in the model's distance units, location and random-effect indices,
##' denominators, offset, correlation and nugget settings, inverse link and
##' `control_mcml`. They go straight to the estimation step of [glgpm()]
##' (`fit_glgpm_design()`), from the estimate, without the final conditional
##' samples and quietly, so the user's formula and data are not evaluated
##' again and sf is not called. The returned function takes a response and a
##' seed and returns the estimates on their natural scale, named as
##' `unlist(coef(object))`.
##' @noRd
make_bootstrap_refit <- function(object) {
  family <- object$family
  invlink <- if (family != "gaussian" &&
                 identical(object$link_function$name, "custom")) {
    object$link_function[c("inv", "d1", "d2")]
  }
  design <- as.matrix(object$D)
  re_unique <- if (!is.null(object$ID_re)) lapply(object$re, seq_along)
  start <- coef(object)
  control_mcml <- attr(object, "control_mcml") %||% set_control_mcml()

  function(response, seed) {
    fit <- fit_glgpm_design(
      y = response, D = design, coords = object$coords,
      ID_coords = object$ID_coords, ID_re = object$ID_re,
      s_unique = unique(object$ID_coords), re_unique = re_unique,
      units_m = object$units_m, cov_offset = object$cov_offset,
      kappa = object$kappa, fix_tau2 = object$fix_tau2,
      fix_var_me = object$fix_var_me, family = family, invlink = invlink,
      start_pars = start, control_mcmc = set_control_mcmc(),
      control_mcml = control_mcml, seed = seed, return_samples = FALSE,
      messages = FALSE
    )
    fit$D <- design
    fit$family <- family
    unlist(coef.RiskMap(fit))
  }
}

##' Simulate and refit one bootstrap replicate
##'
##' Draws a response from the prepared simulation (`prepare_simulation()`)
##' and refits the model to it, both with `seed`. Returns the estimates with
##' the status of the refit, its error message and its warnings, instead of
##' stopping when the refit fails.
##' @noRd
run_bootstrap_replicate <- function(prepared, refit, seed) {
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    samples <- draw_simulation(prepared, nsim = 1, seed = seed)
    response <- as.numeric(samples$data[, 1L, "response"])
    list(status = "ok", message = NA_character_,
         estimate = refit(response, seed))
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    list(status = "failed", message = conditionMessage(e), estimate = NULL)
  })
  result$warnings <- unique(warnings)
  result
}

##' Run the bootstrap replicates, in parallel by forking when `cores` is
##' above 1
##'
##' Every replicate runs with one BLAS and OpenMP thread when RhpcBLASctl is
##' installed, so that its result does not depend on `cores`. The calling
##' process switches to one thread before any replicate runs and restores the
##' user's settings afterwards. It does so before forking too, because a
##' worker forked from a process whose BLAS uses several threads can
##' deadlock. One after another, progress is reported about ten times.
##' @noRd
run_bootstrap_replicates <- function(replicates, work, cores, messages) {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    blas_threads <- RhpcBLASctl::blas_get_num_procs()
    omp_threads <- RhpcBLASctl::omp_get_max_threads()
    on.exit({
      RhpcBLASctl::blas_set_num_threads(blas_threads)
      RhpcBLASctl::omp_set_num_threads(omp_threads)
    }, add = TRUE)
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  if (cores > 1L && length(replicates) > 1L) {
    results <- parallel::mclapply(replicates, work,
                                  mc.cores = min(cores, length(replicates)),
                                  mc.preschedule = TRUE)
    # A worker that dies returns an error object in place of its results.
    lost <- !vapply(results, is.list, logical(1))
    results[lost] <- lapply(results[lost], function(result) {
      list(status = "failed", estimate = NULL, warnings = character(),
           message = paste("The worker process failed:",
                           as.character(result)))
    })
    return(results)
  }

  started <- proc.time()[["elapsed"]]
  report_every <- max(1L, length(replicates) %/% 10L)
  results <- vector("list", length(replicates))
  for (k in seq_along(replicates)) {
    results[[k]] <- work(replicates[k])
    if (messages && k %% report_every == 0L && k < length(replicates)) {
      elapsed <- proc.time()[["elapsed"]] - started
      message(k, " of ", length(replicates), " replicates done, about ",
              format_bootstrap_time(elapsed / k * (length(replicates) - k)),
              " left")
    }
  }
  results
}

##' @noRd
format_bootstrap_time <- function(seconds) {
  if (seconds < 1.5) return("1 second")
  if (seconds < 90) return(paste(round(seconds), "seconds"))
  if (seconds < 5400) return(paste(round(seconds / 60), "minutes"))
  paste(round(seconds / 3600, 1), "hours")
}

##' Collect the replicates into a `RiskMap_bootstrap` object
##' @noRd
assemble_bootstrap <- function(object, results, replicate_seeds) {
  estimate <- unlist(coef(object))
  estimates <- matrix(NA_real_, length(results), length(estimate),
                      dimnames = list(NULL, names(estimate)))
  failed <- integer()
  failure_messages <- character()
  warning_rows <- list()
  for (i in seq_along(results)) {
    result <- results[[i]]
    if (identical(result$status, "ok")) {
      estimates[i, ] <- result$estimate[names(estimate)]
    } else {
      failed <- c(failed, i)
      failure_messages <- c(failure_messages, result$message)
    }
    if (length(result$warnings)) {
      warning_rows[[length(warning_rows) + 1L]] <- data.frame(
        replicate = i, seed = replicate_seeds[i], message = result$warnings
      )
    }
  }
  warnings <- if (length(warning_rows)) {
    do.call(rbind, warning_rows)
  } else {
    data.frame(replicate = integer(), seed = integer(), message = character())
  }
  out <- list(
    estimates = estimates,
    estimate = estimate,
    labels = parameter_labels(names(estimate), object$distance_units),
    failures = data.frame(replicate = failed, seed = replicate_seeds[failed],
                          message = failure_messages),
    warnings = warnings
  )
  class(out) <- "RiskMap_bootstrap"
  out
}

##' Variances and the scale, which are positive and can be shown on a log
##' scale
##' @noRd
bootstrap_positive <- function(names) !grepl("^beta\\.", names)

##' @title Summarise a parametric bootstrap
##' @description Computes, for each parameter, the bootstrap median, bias,
##'   standard error and the percentile and bias-corrected intervals described
##'   in [parametric_bootstrap()].
##' @param object A `RiskMap_bootstrap` object from [parametric_bootstrap()].
##' @param level Probability covered by the intervals. Defaults to the level
##'   of the bootstrap.
##' @param ... Not used.
##' @return An object of class `summary.RiskMap_bootstrap`, a list with the
##'   table of results (`table`), the `level`, the numbers of replicates and
##'   of failed refits (`n_sim`, `n_failed`) and the `seed`.
##' @method summary RiskMap_bootstrap
##' @export
summary.RiskMap_bootstrap <- function(object, level = object$settings$level,
                                      ...) {
  check_range(level, min = 0, max = 1, allow_equal = FALSE, name = "level")
  alpha <- 1 - level
  limits <- c(alpha / 2, 1 - alpha / 2)
  rows <- lapply(names(object$estimate), function(name) {
    values <- object$estimates[, name]
    values <- values[!is.na(values)]
    estimate <- object$estimate[[name]]
    below <- mean(values < estimate)
    percentile <- stats::quantile(values, limits, names = FALSE)
    bias_corrected <- if (below > 0 && below < 1) {
      z0 <- stats::qnorm(below)
      stats::quantile(values, stats::pnorm(2 * z0 + stats::qnorm(limits)),
                      names = FALSE)
    } else {
      c(NA_real_, NA_real_)
    }
    data.frame(
      parameter = name, label = object$labels[[name]], estimate = estimate,
      median = stats::median(values), bias = mean(values) - estimate,
      se = stats::sd(values), below = below,
      percentile_lower = percentile[1L], percentile_upper = percentile[2L],
      bc_lower = bias_corrected[1L], bc_upper = bias_corrected[2L]
    )
  })
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  out <- list(table = table, level = level,
              n_sim = object$settings$n_sim,
              n_failed = nrow(object$failures),
              seed = object$settings$seed,
              family = object$settings$family)
  class(out) <- "summary.RiskMap_bootstrap"
  out
}

##' @title Print a parametric bootstrap
##' @description Prints the summary of a parametric bootstrap from
##'   [parametric_bootstrap()].
##' @param x A `RiskMap_bootstrap` or `summary.RiskMap_bootstrap` object.
##' @param ... Not used.
##' @return `x`, invisibly.
##' @method print RiskMap_bootstrap
##' @export
print.RiskMap_bootstrap <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

##' @rdname print.RiskMap_bootstrap
##' @method print summary.RiskMap_bootstrap
##' @export
print.summary.RiskMap_bootstrap <- function(x, ...) {
  family <- c(gaussian = "linear", binomial = "binomial",
              poisson = "Poisson")[[x$family]]
  cat("Parametric bootstrap of a ", family, " geostatistical model: ",
      x$n_sim, " replicates (seed ", x$seed, ")\n", sep = "")
  if (x$n_failed > 0L) {
    cat(x$n_failed, " refits failed and are left out\n", sep = "")
  }
  cat("Variances and the scale on their natural scale\n")
  significant <- function(values) {
    ifelse(is.na(values), "NA", vapply(signif(values, 4), format, ""))
  }
  distribution <- data.frame(
    Parameter = x$table$label,
    Estimate = significant(x$table$estimate),
    Median = significant(x$table$median),
    Bias = significant(x$table$bias),
    "Std. error" = significant(x$table$se),
    Below = sprintf("%.2f", x$table$below),
    check.names = FALSE
  )
  cat("\nBootstrap distribution (bias: mean minus estimate; below: proportion",
      "below the\nestimate)\n")
  print(distribution, row.names = FALSE, right = TRUE)
  intervals <- data.frame(
    Parameter = x$table$label,
    "Percentile lower" = significant(x$table$percentile_lower),
    "Percentile upper" = significant(x$table$percentile_upper),
    "BC lower" = significant(x$table$bc_lower),
    "BC upper" = significant(x$table$bc_upper),
    check.names = FALSE
  )
  cat("\n", round(100 * x$level), "% intervals: percentile and ",
      "bias-corrected (BC)\n", sep = "")
  print(intervals, row.names = FALSE, right = TRUE)
  if (anyNA(x$table[, c("bc_lower", "bc_upper")])) {
    cat("\nNA: every bootstrap estimate is on one side of the estimate, so",
        "the BC interval\nis not available.\n")
  }
  invisible(x)
}

##' Plot a parametric bootstrap
##'
##' Draws a histogram of the bootstrap estimates of each parameter, with the
##' estimate, the percentile interval and the bias-corrected interval (see
##' [parametric_bootstrap()]). Every parameter is shown on its natural axis,
##' unless `log_scale = TRUE`.
##'
##' @param x A `RiskMap_bootstrap` object from [parametric_bootstrap()].
##' @param parameters Names of the parameters to draw, as in the columns of
##'   `x$estimates`. The default is all of them.
##' @param bins Number of bins of each histogram.
##' @param log_scale Logical, whether to show the variances and the scale on
##'   a logarithmic axis, which spreads out bootstrap estimates near zero.
##'   Estimates at zero are then left out. The regression coefficients stay
##'   on their natural axis.
##' @param ... Not used.
##' @return A `ggplot` object for one parameter, or a `patchwork` object
##'   combining one `ggplot` per parameter.
##' @seealso [parametric_bootstrap()]
##' @method plot RiskMap_bootstrap
##' @export
plot.RiskMap_bootstrap <- function(x, parameters = NULL, bins = 30,
                                   log_scale = FALSE, ...) {
  stopifnot(inherits(x, "RiskMap_bootstrap"))
  parameters <- parameters %||% names(x$estimate)
  unknown <- setdiff(parameters, names(x$estimate))
  if (length(unknown)) {
    stop("Unknown parameter(s): ", paste(unknown, collapse = ", "),
         ". Choose from: ", paste(names(x$estimate), collapse = ", "),
         call. = FALSE)
  }
  check_positive_integer(bins, "bins")
  check_logical(log_scale)
  table <- summary(x)$table
  panels <- lapply(parameters, function(name) {
    plot_bootstrap_parameter(x$estimates[, name],
                             table[table$parameter == name, ],
                             log_scale && bootstrap_positive(name), bins,
                             round(100 * x$settings$level))
  })
  if (length(panels) == 1L) return(panels[[1L]])
  patchwork::wrap_plots(panels, ncol = min(2L, length(panels))) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

##' Histogram of the bootstrap estimates of one parameter
##' @noRd
plot_bootstrap_parameter <- function(values, row, log_scale, bins, level) {
  values <- values[!is.na(values)]
  if (log_scale) values <- values[values > 0]
  lines <- data.frame(
    value = c(row$estimate, row$percentile_lower, row$percentile_upper,
              row$bc_lower, row$bc_upper),
    marker = factor(c("Estimate", rep(paste(level, "% percentile interval",
                                            sep = ""), 2L),
                      rep(paste(level, "% BC interval", sep = ""), 2L)),
                    levels = c("Estimate",
                               paste(level, "% percentile interval", sep = ""),
                               paste(level, "% BC interval", sep = "")))
  )
  lines <- lines[is.finite(lines$value), ]
  colours <- stats::setNames(c("#1d232b", "#0072B2", "#D55E00"),
                             levels(lines$marker))
  linetypes <- stats::setNames(c("solid", "dashed", "dotdash"),
                               levels(lines$marker))
  plot <- ggplot2::ggplot(data.frame(value = values),
                          ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = "#cfd6de", colour = "white") +
    ggplot2::geom_vline(
      data = lines,
      ggplot2::aes(xintercept = .data$value, colour = .data$marker,
                   linetype = .data$marker),
      linewidth = 0.7
    ) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE, name = NULL) +
    ggplot2::scale_linetype_manual(values = linetypes, drop = FALSE,
                                   name = NULL) +
    ggplot2::labs(title = row$label,
                  x = if (log_scale) "Bootstrap estimate (log scale)" else
                    "Bootstrap estimate",
                  y = "Replicates") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", size = 11),
                   legend.position = "bottom")
  if (log_scale) {
    plot <- plot + ggplot2::scale_x_log10(labels = plain_number_labels())
  }
  plot
}
