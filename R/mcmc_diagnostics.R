##' Effective sample sizes of the MCMC samples of the random effects
##'
##' Computes the effective sample size of the Markov chain Monte Carlo (MCMC)
##' samples of the random effects stored in a fitted binomial or Poisson
##' model, for the mean of the spatial random effects over the locations, for
##' the spatial random effect at each location and for each level of every
##' unstructured random effect.
##'
##' @param object A fitted `RiskMap` object returned by [glgpm()] for a
##'   binomial or Poisson model with `return_samples = TRUE`.
##'
##' @details
##' The samples are those drawn by [glgpm()] at the fitted parameter values,
##' which are used for the fitted values and for spatial prediction. When the
##' model has a nugget, the spatial random effect at a location includes it.
##'
##' The effective sample size of a chain is the number of independent draws
##' that would estimate the mean of its target distribution as precisely as
##' the correlated MCMC samples do. It is the number of samples divided by
##' the integrated autocorrelation time, \eqn{1 + 2 \sum_{k \ge 1} \rho_k},
##' where \eqn{\rho_k} is the autocorrelation at lag \eqn{k}. The sum is
##' estimated with the initial monotone sequence estimator of Geyer (1992),
##' which adds pairs of consecutive sample autocorrelations while their sum is
##' positive and forces these sums to decrease. As in Stan, the
##' autocorrelation time is not allowed to fall below \eqn{1 / \log_{10} N}
##' for \eqn{N} samples, so the effective sample size of an antithetic chain
##' is at most \eqn{N \log_{10} N}. A chain that never moves has an effective
##' sample size of zero.
##'
##' @return A tibble with one row per chain and columns
##' \item{effect}{`"mean of S"` for the mean of the spatial random effects
##'   over the locations, `"S"` for the spatial random effect at a location,
##'   or the name of the variable of an unstructured random effect.}
##' \item{component}{The location, numbered as the rows of `object$coords`,
##'   or the level of the random effect, numbered as in `object$re`. It is
##'   `NA` for the mean.}
##' \item{ess}{The effective sample size.}
##' The number of samples is stored in the attribute `"n_samples"`.
##' @references Geyer, C. J. (1992). Practical Markov chain Monte Carlo.
##'   *Statistical Science*, 7(4), 473-483.
##' @seealso [plot_mcmc()]
##' @examples
##' \donttest{
##' data(liberia)
##' liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
##' fit <- glgpm(npos ~ log(elevation) + gp(), denominator = ntest,
##'              data = liberia, model_crs = 32629, family = "binomial",
##'              messages = FALSE)
##' ess <- summarise_mcmc(fit)
##' # The mean and the five locations with the smallest effective sample size
##' ess[1, ]
##' head(ess[order(ess$ess), ], 5)
##' }
##' @export
summarise_mcmc <- function(object) {
  blocks <- mcmc_sample_blocks(object)
  spatial <- blocks$S
  out <- data.frame(
    effect = "mean of S",
    component = NA_integer_,
    ess = mcmc_effective_size(rowMeans(spatial))
  )
  for (name in names(blocks)) {
    out <- rbind(out, data.frame(
      effect = name,
      component = seq_len(ncol(blocks[[name]])),
      ess = mcmc_effective_size(blocks[[name]])
    ))
  }
  out <- as_tibble(out)
  attr(out, "n_samples") <- nrow(spatial)
  out
}

##' Check the mixing of the MCMC samples of the random effects
##'
##' Draws diagnostics of the Markov chain Monte Carlo (MCMC) samples of the
##' random effects stored in a fitted binomial or Poisson model. For the mean
##' of the spatial random effects over the locations, and for the location
##' whose samples have the smallest effective sample size, it shows a trace
##' plot and the autocorrelation function. A histogram of the effective sample
##' sizes at all the locations shows whether some locations mix much worse
##' than the rest.
##'
##' @inheritParams summarise_mcmc
##' @param check_mean Logical. If `TRUE` (the default), the chain of the mean
##'   of the spatial random effects over the locations is shown.
##' @param component The locations whose chains are shown, numbered as the
##'   rows of `object$coords`. The default, `NULL`, shows the location with
##'   the smallest effective sample size.
##' @param ... Further arguments passed to [stats::acf()], such as `lag.max`.
##'
##' @details
##' Inspecting the mean \eqn{\bar{S} = \sum_{i=1}^n S(x_i) / n} avoids
##' checking every component of the spatial random effect, but poor mixing at
##' a few locations can be hidden by the average. The location with the
##' smallest effective sample size is the one most worth inspecting, and the
##' histogram shows how much the effective sample sizes vary between
##' locations. The effective sample sizes, of unstructured random effects too,
##' are returned by [summarise_mcmc()], whose help page describes how they are
##' estimated.
##'
##' The trace plots show the iterations of the chain at which the samples were
##' kept, after the burn-in and thinning set by [set_control_mcmc()]. The
##' autocorrelation is shown against the lag between kept samples, with
##' dashed lines at \eqn{\pm 1.96 / \sqrt{N}} for \eqn{N} samples. For models
##' fitted with versions of RiskMap that did not record the burn-in and
##' thinning, the samples are numbered in order.
##'
##' @return A `patchwork` object, which is drawn when printed and can be
##'   modified like a `ggplot` object.
##' @seealso [summarise_mcmc()], [glgpm()], [set_control_mcmc()]
##' @examples
##' \donttest{
##' data(liberia)
##' liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
##' fit <- glgpm(npos ~ log(elevation) + gp(), denominator = ntest,
##'              data = liberia, model_crs = 32629, family = "binomial",
##'              messages = FALSE)
##' plot_mcmc(fit)
##' plot_mcmc(fit, check_mean = FALSE, component = c(10, 20))
##' }
##' @export
plot_mcmc <- function(object, check_mean = TRUE, component = NULL, ...) {
  check_logical(check_mean)
  spatial <- mcmc_sample_blocks(object)$S
  spatial_ess <- mcmc_effective_size(spatial)
  mean_chain <- rowMeans(spatial)
  mean_ess <- mcmc_effective_size(mean_chain)
  n_locations <- ncol(spatial)
  smallest <- is.null(component)
  if (smallest) {
    component <- which.min(spatial_ess)
  } else if (!is.numeric(component) || length(component) == 0L ||
             !all(is.finite(component)) || any(component %% 1 != 0) ||
             any(component < 1) || any(component > n_locations) ||
             anyDuplicated(component) > 0L) {
    stop("'component' must contain distinct whole numbers between 1 and ",
         "the number of locations, ", n_locations, call. = FALSE)
  }
  component <- as.integer(component)

  chains <- lapply(component, function(location) spatial[, location])
  labels <- if (smallest) {
    sprintf("S at location %d, the smallest ESS", component)
  } else {
    sprintf("S at location %d", component)
  }
  chain_ess <- spatial_ess[component]
  if (check_mean) {
    chains <- c(list(mean_chain), chains)
    labels <- c(sprintf("Mean of S over %d locations", n_locations), labels)
    chain_ess <- c(mean_ess, chain_ess)
  }

  iterations <- retained_iterations(object, nrow(spatial))
  traces <- plot_mcmc_traces(chains, labels, iterations)
  correlations <- plot_mcmc_autocorrelations(chains, chain_ess, ...)
  histogram <- plot_mcmc_ess_histogram(spatial_ess, mean_ess, nrow(spatial))
  patchwork::wrap_plots(traces, correlations, nrow = 1L,
                        widths = c(2.3, 1)) /
    histogram +
    patchwork::plot_layout(heights = c(length(chains), 1.3))
}

##' Validate a fit for MCMC diagnostics and split its samples into blocks
##' @noRd
mcmc_sample_blocks <- function(object) {
  if (!inherits(object, "RiskMap")) {
    stop("'object' must be a 'RiskMap' object returned by 'glgpm'",
         call. = FALSE)
  }
  if (identical(object$family, "gaussian")) {
    stop("'object' is a Gaussian model, which has no MCMC samples",
         call. = FALSE)
  }
  if (is.null(object$S_samples)) {
    stop("'object' does not contain any MCMC chains. Refit the model with ",
         "glgpm(..., return_samples = TRUE)", call. = FALSE)
  }
  stored_conditioning_sample(object)$samples
}

##' Iterations of the chain at which the stored samples were kept
##'
##' Returns `NULL` for fits that do not record the burn-in and thinning.
##' @noRd
retained_iterations <- function(object, n_samples) {
  settings <- attr(object, "mcmc")
  if (is.null(settings$burnin) || is.null(settings$thin)) return(NULL)
  settings$burnin + settings$thin * seq_len(n_samples)
}

##' Effective sample sizes of the columns of a matrix of MCMC samples
##'
##' Autocovariances are computed for all lags at once with a fast Fourier
##' transform of each centred column, padded with zeros to avoid wrapping
##' around. Columns are processed in blocks to bound the memory used.
##' @noRd
mcmc_effective_size <- function(samples, block_size = 500L) {
  samples <- as.matrix(samples)
  n_samples <- nrow(samples)
  padded_length <- nextn(2L * n_samples)
  out <- numeric(ncol(samples))
  for (first in seq(1L, ncol(samples), by = block_size)) {
    columns <- seq.int(first, min(first + block_size - 1L, ncol(samples)))
    block <- samples[, columns, drop = FALSE]
    centred <- sweep(block, 2L, colMeans(block))
    centred <- rbind(centred,
                     matrix(0, padded_length - n_samples, length(columns)))
    power <- Mod(mvfft(centred))^2
    autocovariance <- Re(mvfft(power, inverse = TRUE))[seq_len(n_samples), ,
                                                      drop = FALSE]
    constant <- colSums(block != block[rep(1L, n_samples), , drop = FALSE]) == 0L
    out[columns] <- ifelse(constant, 0, apply(autocovariance, 2L,
                                              initial_monotone_ess))
  }
  out
}

##' Geyer's initial monotone sequence estimate of the effective sample size
##' @noRd
initial_monotone_ess <- function(autocovariance) {
  n_samples <- length(autocovariance)
  autocorrelation <- autocovariance / autocovariance[1L]
  n_pairs <- n_samples %/% 2L
  pair_sums <- autocorrelation[2L * seq_len(n_pairs) - 1L] +
    autocorrelation[2L * seq_len(n_pairs)]
  first_non_positive <- match(TRUE, pair_sums <= 0, nomatch = n_pairs + 1L)
  pair_sums <- cummin(pair_sums[seq_len(first_non_positive - 1L)])
  autocorrelation_time <- max(2 * sum(pair_sums) - 1, 1 / log10(n_samples))
  n_samples / autocorrelation_time
}

##' Colours and theme of the MCMC diagnostic panels
##' @noRd
mcmc_palette <- function() {
  list(chain = "#0072B2", bound = "#7a7f87", fill = "#cfe3f2",
       mean = "#D55E00")
}

##' @noRd
mcmc_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold", size = 9.5))
}

##' Format effective sample sizes for labels
##' @noRd
format_ess <- function(ess) {
  format(round(ess), big.mark = ",", trim = TRUE)
}

##' Trace plots of the chains, one facet per chain
##' @noRd
plot_mcmc_traces <- function(chains, labels, iterations) {
  n_samples <- length(chains[[1L]])
  traces <- data.frame(
    iteration = rep(iterations %||% seq_len(n_samples), length(chains)),
    value = unlist(chains, use.names = FALSE),
    chain = factor(rep(labels, each = n_samples), levels = labels)
  )
  ggplot(traces, aes(x = .data$iteration, y = .data$value)) +
    geom_line(colour = mcmc_palette()$chain, linewidth = 0.25) +
    facet_wrap(~ chain, ncol = 1L, scales = "free_y") +
    scale_x_continuous(labels = scales::label_comma()) +
    labs(x = if (is.null(iterations)) "Sample" else "Iteration", y = NULL) +
    mcmc_theme()
}

##' Autocorrelation functions of the chains, labelled with their effective
##' sample sizes
##' @noRd
plot_mcmc_autocorrelations <- function(chains, chain_ess, ...) {
  colours <- mcmc_palette()
  labels <- paste("ESS =", format_ess(chain_ess))
  # Identical sizes would merge facets, so each label is made unique.
  labels <- make.unique(labels, sep = " ")
  correlations <- do.call(rbind, Map(function(chain, label) {
    correlation <- acf(chain, plot = FALSE, ...)
    data.frame(lag = as.numeric(correlation$lag),
               autocorrelation = as.numeric(correlation$acf),
               chain = label)
  }, chains, labels))
  correlations$chain <- factor(correlations$chain, levels = labels)
  bound <- qnorm(0.975) / sqrt(length(chains[[1L]]))
  ggplot(correlations, aes(x = .data$lag, y = .data$autocorrelation)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = c(-bound, bound), colour = colours$bound,
               linetype = "dashed", linewidth = 0.3) +
    geom_segment(aes(xend = .data$lag, yend = 0), colour = colours$chain,
                 linewidth = 0.5) +
    facet_wrap(~ chain, ncol = 1L) +
    labs(x = "Lag", y = "Autocorrelation") +
    mcmc_theme()
}

##' Histogram of the effective sample sizes at the locations, with that of
##' their mean marked
##' @noRd
plot_mcmc_ess_histogram <- function(spatial_ess, mean_ess, n_samples) {
  colours <- mcmc_palette()
  sizes <- data.frame(
    ess = spatial_ess,
    panel = sprintf("Effective sample sizes at the %d locations, from %s samples",
                    length(spatial_ess), format(n_samples, big.mark = ","))
  )
  ggplot(sizes, aes(x = .data$ess)) +
    geom_histogram(bins = max(10L, ceiling(sqrt(length(spatial_ess)))),
                   fill = colours$fill, colour = colours$chain,
                   linewidth = 0.3) +
    geom_vline(xintercept = mean_ess, colour = colours$mean,
               linewidth = 0.6) +
    # The label goes on the side of the line with more room.
    annotate("text", x = mean_ess, y = Inf, vjust = 1.5,
             hjust = if (mean_ess > mean(range(spatial_ess))) 1.08 else -0.08,
             label = paste("Mean of S:", format_ess(mean_ess)),
             colour = colours$mean, size = 3.2) +
    facet_wrap(~ panel) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "Effective sample size", y = "Locations") +
    mcmc_theme()
}
