##' The predictive distribution of a count, averaged over the draws of its mean
##'
##' F(k), the average over the draws of the mean in `means` of the family's
##' distribution function of a count with `units_m` trials, or exposure
##' `units_m`, and the family's `parameters`. No count is simulated, so the
##' draws of the mean are the only random step in the scores of a count.
##'
##' F is evaluated at every count up to 100, at y - 1 and y, for the observed
##' count `y`, and at the largest count that is scored (`score_limit()`).
##' Between two of these counts a and b it is taken as linear. The true F and
##' the line both lie between F(a) and F(b), so the interpolation changes the
##' sum of F over the counts of the cell by at most (b - a) (F(b) - F(a)). A
##' cell where that exceeds `tolerance` times the smaller of the CRPS and its
##' expected value is halved, until none does. The grid is then fine wherever
##' the draws put probability, whether their means are close together or
##' differ by orders of magnitude, and the relative error of the CRPS and of
##' its expected value is of the order of `tolerance` or less. Cells of one
##' count are exact, as is every F when the largest count is at most 100.
##'
##' @param y The observed count.
##' @param units_m Its number of trials or exposure.
##' @param means The draws of the mean of the inverse link, checked by
##'   `check_count_means()`.
##' @param parameters The family's parameters, as `log_density()` takes them.
##' @param traits The family's traits (`family_traits()`).
##' @param tolerance The relative error allowed in a cell.
##' @return A list with `counts`, which start at -1, where F is zero, and
##'   `values`, F at each.
##' @noRd
count_distribution <- function(y, units_m, means, parameters, traits,
                               tolerance = 1e-4) {
  average <- function(counts) {
    values <- suppressWarnings(vapply(counts, function(count) {
      mean(traits$distribution_function(count, means, units_m, parameters))
    }, numeric(1)))
    # As for a number of trials that isn't whole.
    if (anyNA(values)) {
      riskmap_abort("The predictive distribution of a count with a denominator of {units_m} is not a number.")
    }
    values
  }
  limit <- traits$score_limit(y, units_m, means, parameters)
  counts <- sort(unique(c(-1, 0:min(limit, 100), y - 1, y, limit)))
  values <- cummax(c(0, average(counts[-1])))
  repeat {
    distribution <- list(counts = counts, values = values)
    lower <- counts[-length(counts)]
    upper <- counts[-1]
    middle <- floor((lower + upper) / 2)
    # Far beyond 10^15 neighbouring numbers are more than one count apart, and
    # a cell between two of them has no middle.
    wide <- middle > lower & middle < upper
    largest <- tolerance * min(count_crps(y, distribution),
                               count_expected_crps(distribution))
    added <- middle[wide & (upper - lower) * diff(values) > largest]
    if (!length(added)) return(distribution)
    position <- order(c(counts, added))
    values <- cummax(c(values, average(added))[position])
    counts <- c(counts, added)[position]
  }
}

##' The CRPS of a count, the sum over k of (F(k) - 1(k >= y))^2
##'
##' The sum over the counts of each cell of `distribution`
##' (`count_distribution()`), on which F is linear, in closed form. The
##' distribution must have y - 1 and y among its counts. Counts beyond the
##' largest are left out.
##' @noRd
count_crps <- function(y, distribution) {
  width <- diff(distribution$counts)
  rise <- diff(distribution$values)
  last <- distribution$values[-1] - (distribution$counts[-1] >= y)
  sum(width * last^2 - last * rise * (width - 1) +
        rise^2 * (width - 1) * (2 * width - 1) / (6 * width))
}

##' The expected CRPS of a count from its own predictive distribution
##'
##' The sum over k of F(k) (1 - F(k)), which equals half the expected distance
##' between two draws from F, over the cells of `distribution` as in
##' `count_crps()`.
##' @noRd
count_expected_crps <- function(distribution) {
  width <- diff(distribution$counts)
  rise <- diff(distribution$values)
  last <- distribution$values[-1]
  sum(width * last * (1 - last) + (2 * last - 1) * rise * (width - 1) / 2 -
        rise^2 * (width - 1) * (2 * width - 1) / (6 * width))
}

##' The predictive distribution function of a count at any counts
##'
##' Zero below zero, and its last value beyond the largest count scored.
##' @noRd
count_distribution_at <- function(distribution, counts) {
  stats::approx(distribution$counts, distribution$values, xout = counts,
                rule = 2)$y
}

##' The smallest count at which a predictive distribution reaches `probability`
##'
##' The inverse of `count_distribution_at()` over whole counts, and the largest
##' count scored for a probability that the distribution does not reach.
##' @noRd
count_quantile <- function(distribution, probability) {
  counts <- distribution$counts
  values <- distribution$values
  cell <- findInterval(probability, values, left.open = TRUE)
  cell <- pmin(cell, length(counts) - 1L)
  share <- (probability - values[cell]) / (values[cell + 1L] - values[cell])
  steps <- ceiling(share * (counts[cell + 1L] - counts[cell]))
  pmin(counts[cell] + steps, counts[cell + 1L])
}

##' Check that the draws of the mean of counts can be scored
##'
##' Stops if a draw of the mean is not a number, is outside the range of the
##' family, or gives an infinite expected count, so that the test set is
##' recorded as a failure instead of getting a score without that draw.
##' @param means A matrix of draws, with a row for each test observation.
##' @param units_m The denominators of the test observations.
##' @param traits The family's traits (`family_traits()`).
##' @noRd
check_count_means <- function(means, units_m, traits) {
  range <- traits$mean_range
  usable <- is.finite(means * units_m) & means >= range[1L] &
    means <= range[2L]
  n_draws <- sum(!usable)
  if (n_draws > 0L) {
    n_observations <- sum(rowSums(!usable) > 0)
    riskmap_abort("{n_draws} draw{?s} of the mean at {n_observations} test observation{?s} {cli::qty(n_draws)}{?isn't/aren't} finite and in the range of the family, so the counts can't be scored.")
  }
  invisible(TRUE)
}
