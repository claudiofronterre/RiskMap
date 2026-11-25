# Calibration curves of cross-validated predictions, for assess_prediction().

##' The probability levels of the stored calibration curves
##'
##' A hundred steps, so that the bins of a histogram of 5, 10 or 20 bins have
##' their edges among them.
##' @noRd
calibration_levels <- function() seq(0, 1, length.out = 101L)

##' The number of calibration curves simulated from a model
##'
##' With the curve of the data they are a thousand, and the 95% limits are order
##' statistics.
##' @noRd
n_simulated_curves <- function() 999L

##' The probability of the bands and of the limit of the area of a correct model
##' @noRd
calibration_level <- function() 0.95

##' The probability levels at which the areas are computed, ten times as many as
##' are stored
##' @noRd
area_levels <- function() seq(0, 1, length.out = 1001L)

##' The calibration curve of a model and its reference for a correct model
##'
##' @param pit_lower,pit_upper The ends of the PIT of the test observations.
##' @param reference The `curves` and `areas` from `calibration_reference()`.
##' @return A list with the `area` between the curve and the identity line, its
##'   `expected` value and its 95% `limit` for a correct model, the `curve` with
##'   its 95% joint band at `calibration_levels()` and the `simulated` curves, a
##'   row for each simulation. `NULL` values when nothing was scored.
##' @noRd
assess_calibration <- function(pit_lower, pit_upper, reference) {
  if (!length(pit_lower)) {
    return(list(area = NA_real_, expected = NA_real_, limit = NA_real_,
                curve = NULL, simulated = NULL))
  }
  curve <- pit_curves(t(pit_lower), t(pit_upper), area_levels())
  stored <- curve[, seq(1L, ncol(curve), by = 10L)]
  band <- calibration_band(stored, reference$curves)
  list(area = pit_area(curve, area_levels()),
       expected = mean(reference$areas), limit = simulated_limit(reference$areas),
       curve = data.frame(u = calibration_levels(), curve = stored,
                          lower = band$lower, upper = band$upper),
       simulated = reference$curves)
}

##' Calibration curves of a correct model
##'
##' The curves of observations simulated from the model, which the expected area,
##' its limit and the bands come from. They are drawn jointly within each test
##' set, from its joint predictive distribution, and independently between test
##' sets. An observation in several test sets has the same PIT value in all of
##' them, the one drawn in its first test set.
##'
##' Nothing else in the assessment depends on how they are simulated.
##' @param fit The fit that is assessed.
##' @param test_sets A `RiskMap_test_sets` object.
##' @param predictive The predictive distributions of each test set, from
##'   `score_test_set()`, `NULL` for one that has no scores.
##' @param n_curves The number of curves. There are fewer when the test sets have
##'   fewer predictive draws.
##' @param seeds A seed for each test set.
##' @return A list with the `curves` at `calibration_levels()`, a matrix with a
##'   curve in each row, and their `areas`.
##' @noRd
calibration_reference <- function(fit, test_sets, predictive, n_curves, seeds) {
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  traits <- family_traits(fit$family)
  simulate <- if (traits$discrete) simulate_count_pit else simulate_gaussian_pit
  # The simulated PIT values of each observation, from its first test set.
  uniforms <- vector("list", length(test_sets$location))
  ends <- vector("list", length(predictive))
  for (i in seq_along(predictive)) {
    if (is.null(predictive[[i]])) next
    set.seed(seeds[i])
    rows <- predictive[[i]]$rows
    simulated <- simulate(predictive[[i]], traits, uniforms[rows])
    uniforms[rows] <- simulated$uniforms
    ends[[i]] <- simulated[c("lower", "upper")]
  }
  bind <- function(name) {
    values <- do.call(cbind, lapply(ends, `[[`, name))
    values[seq_len(min(n_curves, nrow(values))), , drop = FALSE]
  }
  curves <- pit_curves(bind("lower"), bind("upper"), area_levels())
  list(curves = curves[, seq(1L, ncol(curves), by = 10L), drop = FALSE],
       areas = pit_area(curves, area_levels()))
}

##' Simulated PIT values of a Gaussian test set
##'
##' @param predictive The predictive distributions of a test set.
##' @param traits The traits of the family.
##' @param earlier A list with the PIT values that each observation has in an
##'   earlier test set, `NULL` for one that isn't in any.
##' @return A list with the `lower` and `upper` ends of the PIT, matrices with a
##'   simulation in each row and an observation in each column, and the PIT
##'   values, `uniforms`, a list with a vector for each observation.
##' @noRd
simulate_gaussian_pit <- function(predictive, traits, earlier) {
  pit <- pnorm((predictive$draws - predictive$mean) / predictive$sd)
  # A prediction without spread has all the probability at its mean.
  pit[!predictive$sd > 0, ] <- 1
  pit <- t(pit)
  for (j in which(lengths(earlier) > 0L)) pit[, j] <- earlier[[j]]
  list(lower = pit, upper = pit,
       uniforms = lapply(seq_len(ncol(pit)), function(j) pit[, j]))
}

##' Simulated PIT values of a test set of counts
##'
##' A count is drawn from the family given the mean, and has a randomised PIT
##' value. An observation in an earlier test set has instead the count with the
##' PIT value that it has there.
##' @inheritParams simulate_gaussian_pit
##' @noRd
simulate_count_pit <- function(predictive, traits, earlier) {
  draws <- predictive$draws
  counts <- matrix(traits$simulate(draws, predictive$units_m,
                                   predictive$parameters),
                   nrow = ncol(draws), ncol = nrow(draws), byrow = TRUE)
  mixing <- matrix(stats::runif(length(counts)), nrow(counts))
  lower <- upper <- counts
  uniforms <- earlier
  for (j in seq_along(earlier)) {
    distribution <- predictive$distributions[[j]]
    known <- length(earlier[[j]]) > 0L
    if (known) counts[, j] <- count_quantile(distribution, earlier[[j]])
    lower[, j] <- count_distribution_at(distribution, counts[, j] - 1)
    upper[, j] <- count_distribution_at(distribution, counts[, j])
    if (!known) {
      uniforms[[j]] <- lower[, j] + mixing[, j] * (upper[, j] - lower[, j])
    }
  }
  list(lower = lower, upper = upper, uniforms = uniforms)
}

##' The joint band of a curve, from curves simulated from the model
##'
##' The whole curve of a correct model is inside it with a probability of 95%.
##' The curve of the data is one of those that the band is computed from, as in
##' `variogram()`.
##' @param curve The curve of the data.
##' @param simulated A matrix with a simulated curve in each row.
##' @return A list with the `lower` and `upper` limits.
##' @noRd
calibration_band <- function(curve, simulated) {
  global_rank_envelope(cbind(curve, t(simulated)), calibration_level())
}

##' The value that a correct model is below with a probability of 95%
##'
##' The order statistic of the simulated values that the value of the data,
##' which is one more draw when the model is correct, exceeds with a probability
##' of 5%. It is missing with fewer than 19 simulations.
##' @noRd
simulated_limit <- function(values) {
  rank <- ceiling(calibration_level() * (length(values) + 1L))
  if (rank > length(values)) return(NA_real_)
  sort(values)[rank]
}

##' Average nonrandomised PIT curves
##'
##' At each level u, the average over the observations of a function that is
##' zero up to `lower`, one from `upper` and linear in between. Where the two
##' are equal, as for Gaussian outcomes, it steps from zero to one there.
##'
##' The average is a + b u on each interval of `levels`, so the sums of the
##' intercepts and slopes are accumulated over the levels at which each
##' observation's ramp starts and ends.
##' @param lower,upper Matrices with a curve in each row and an observation in
##'   each column.
##' @param levels Increasing probability levels, from 0 to 1.
##' @return A matrix with a curve in each row and a level in each column.
##' @noRd
pit_curves <- function(lower, upper, levels, block_size = 1e6) {
  n_curves <- nrow(lower)
  n_levels <- length(levels)
  n_observations <- ncol(lower)
  # Blocks of observations, so that many of them don't need large vectors.
  per_block <- max(1L, floor(block_size / n_curves))
  blocks <- split(seq_len(n_observations),
                  (seq_len(n_observations) - 1L) %/% per_block)
  intercepts <- slopes <- matrix(0, n_curves, n_levels + 1L)
  for (columns in blocks) {
    from <- lower[, columns, drop = FALSE]
    to <- upper[, columns, drop = FALSE]
    width <- to - from
    step <- !(width > 1e-9)
    # The first level above `lower`, and the first at or above `upper`.
    position <- c(findInterval(from, levels),
                  findInterval(to, levels, left.open = TRUE))
    cell <- rep(seq_len(n_curves), length.out = length(position)) +
      n_curves * position
    slope <- ifelse(step, 0, 1 / width)
    intercept <- ifelse(step, 0, -from / width)
    change <- rowsum(cbind(c(intercept, 1 - intercept), c(slope, -slope)), cell,
                     reorder = FALSE)
    cells <- as.integer(rownames(change))
    intercepts[cells] <- intercepts[cells] + change[, 1L]
    slopes[cells] <- slopes[cells] + change[, 2L]
  }
  # Sums over the levels up to each.
  accumulate <- function(changes) {
    total <- changes[, seq_len(n_levels), drop = FALSE]
    for (level in seq_len(n_levels)[-1L]) {
      total[, level] <- total[, level] + total[, level - 1L]
    }
    total
  }
  curves <- (accumulate(intercepts) +
               accumulate(slopes) * rep(levels, each = n_curves)) /
    n_observations
  pmin(pmax(curves, 0), 1)
}

##' The area between calibration curves and the identity line
##'
##' @param curves A matrix with a curve in each row.
##' @return A number for each curve, by the trapezoidal rule.
##' @noRd
pit_area <- function(curves, levels) {
  distance <- abs(curves - rep(levels, each = nrow(curves)))
  weights <- (c(diff(levels), 0) + c(0, diff(levels))) / 2
  as.numeric(distance %*% weights)
}

##' The calibrations of the models, as they are in the result
##'
##' @param calibrations A named list of results of `assess_calibration()`.
##' @noRd
combine_calibrations <- function(calibrations) {
  models <- names(calibrations)
  value <- function(name) {
    vapply(calibrations, function(x) x[[name]], numeric(1), USE.NAMES = FALSE)
  }
  curves <- do.call(rbind, lapply(models, function(name) {
    curve <- calibrations[[name]]$curve
    if (is.null(curve)) NULL else cbind(data.frame(model = name), curve)
  }))
  list(area = data.frame(model = models, PIT_area = value("area"),
                         PIT_area_expected = value("expected"),
                         PIT_area_limit = value("limit")),
       curves = curves, simulated = lapply(calibrations, `[[`, "simulated"))
}

##' The heights of the histogram of the PIT, with the band of a correct model
##'
##' The height of a bin is the mean of the curve's increase over the bin,
##' divided by the bin's width, as in the nonrandomised PIT histogram of Czado
##' et al. (2009). A PIT value of exactly 0 belongs to the first bin.
##' @param curves Curves at `calibration_levels()`, a matrix with a row for
##'   each.
##' @param n_bins The number of bins, a divisor of 100.
##' @return A matrix with a row for each curve and a column for each bin.
##' @noRd
pit_histogram <- function(curves, n_bins) {
  edges <- seq(1L, 101L, length.out = n_bins + 1L)
  at_edges <- curves[, edges, drop = FALSE]
  at_edges[, 1L] <- 0
  n_bins * (at_edges[, -1L, drop = FALSE] - at_edges[, -(n_bins + 1L), drop = FALSE])
}
