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
##' The curves of datasets simulated from a Gaussian working model of the fit
##' (`working_model()`), each put through the whole cross-validation, which the
##' expected area, its limit and the bands come from. In each simulated dataset
##' the covariance parameters of the working model are estimated again, as the
##' best of `n_candidates` points around the estimates (`candidate_points()`),
##' from all its observations or, with `refit`, from each training set, and the
##' coefficients by generalised least squares. Each test observation then has
##' the error of its prediction from the training set, in closed form, and the
##' normal distribution function at the standardised error is its PIT value
##' (`candidate_errors()`). For counts it is mapped to a count by the predictive
##' distribution of the data, which gives the ends of the PIT of that count.
##'
##' The data and the simulated datasets thus go through the same cross-validation
##' with estimated parameters, so the areas of the data and of the simulations
##' have about the same distribution for a correct model. An observation in
##' several test sets has a PIT value in each of them, from the same simulated
##' datasets, as it has for the data.
##'
##' Nothing else in the assessment depends on how they are simulated.
##' @param fit The fit that is assessed.
##' @param test_sets A `RiskMap_test_sets` object.
##' @param predictive The predictive distributions of each test set, from
##'   `score_test_set()`, `NULL` for one that has no scores.
##' @param n_curves The number of simulated datasets.
##' @param seeds The seeds of the test sets, the first of which sets the
##'   simulated datasets.
##' @param refit Whether the parameters are estimated again in each training
##'   set.
##' @param n_candidates The number of candidate points of the covariance
##'   parameters, the estimates among them.
##' @return A list with the `curves` at `calibration_levels()`, a matrix with a
##'   curve in each row, and their `areas`.
##' @noRd
calibration_reference <- function(fit, test_sets, predictive, n_curves, seeds,
                                  refit = FALSE, n_candidates = 200L) {
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  set.seed(seeds[1L])
  working <- working_model(fit)
  # The rows that each scored test set tests, and those its training set
  # leaves out, which include rows that it doesn't test.
  scored <- which(!vapply(predictive, is.null, logical(1)))
  sets <- lapply(scored, function(i) {
    list(test = predictive[[i]]$rows,
         left_out = setdiff(seq_along(fit$y), test_sets$sets[[i]]$train))
  })
  # A test observation in each row and a simulated dataset in each column.
  uniforms <- candidate_errors(working, sets, n_curves, n_candidates, refit)
  lower <- upper <- uniforms
  if (family_traits(fit$family)$discrete) {
    distributions <- unlist(lapply(predictive[scored], `[[`, "distributions"),
                            recursive = FALSE)
    for (j in seq_along(distributions)) {
      counts <- count_quantile(distributions[[j]], uniforms[j, ])
      lower[j, ] <- count_distribution_at(distributions[[j]], counts - 1)
      upper[j, ] <- count_distribution_at(distributions[[j]], counts)
    }
  }
  curves <- pit_curves(t(lower), t(upper), area_levels())
  list(curves = curves[, seq(1L, ncol(curves), by = 10L), drop = FALSE],
       areas = pit_area(curves, area_levels()))
}

##' The number of candidate points of the reference of a calibrated model
##'
##' Its time grows about with the cube of the number of observations. On one
##' thread it took 8 seconds at 500 observations and 29 at 1,000 with 200
##' candidates, and 52 at 2,001 with 100, which about halves it.
##' @noRd
reference_candidates <- function(n_observations) {
  if (n_observations > 2000L) 100L else 200L
}

##' The Gaussian working model of a fit
##'
##' The linear predictor of the fit, with its spatial term, nugget and random
##' effects, observed with the measurement error of a Gaussian fit, or else
##' with the variance 1 / w of each record, where w is the curvature of its log
##' likelihood at the conditional mode of the linear predictor given the data
##' (`conditional_weights()`), so that a count is about as informative as it
##' is in the fit.
##' @return A list with the `design` matrix, the covariance parameters `theta`,
##'   the names of those that the fit estimates, `free`, a function that gives
##'   the `covariance` of the records at any `theta`, and its `derivatives` at
##'   `theta` with respect to the logarithms of the free parameters.
##' @noRd
working_model <- function(fit) {
  parameters <- coef(fit)
  design <- as.matrix(fit$D)
  location <- fit$ID_coords
  distances <- pairwise_distances(fit$coords)
  effects <- lapply(fit$ID_re %||% list(), function(level) {
    outer(level, level, "==") * 1
  })
  effect_names <- sprintf("sigma2_re.%d", seq_along(effects))
  traits <- family_traits(fit$family)
  measurement <- if (traits$measurement_error) {
    parameters$sigma2_me %||% fit_fix_var_me(fit) %||% 0
  } else {
    0
  }
  nugget <- fitted_nugget_variance(fit, parameters)
  same_location <- if (nugget > 0) outer(location, location, "==") * 1
  variances_re <- as.numeric(parameters$sigma2_re)[seq_along(effects)]
  names(variances_re) <- effect_names
  theta <- c(sigma2 = parameters$sigma2[[1L]], phi = parameters$phi[[1L]],
             if (nugget > 0) c(tau2 = nugget),
             if (measurement > 0) c(sigma2_me = measurement), variances_re)
  held <- fit_fixed_pars(fit)
  held_names <- setdiff(names(held), "sigma2_re")
  if (!is.null(held$sigma2_re)) {
    terms <- names(held$sigma2_re)
    held_effects <- if (is.null(terms)) {
      seq_along(effects)
    } else {
      match(terms, names(parameters$sigma2_re))
    }
    held_names <- c(held_names, effect_names[held_effects])
  }
  free <- setdiff(names(theta), held_names)
  noise <- if (traits$exact) {
    numeric(length(location))
  } else {
    1 / pmax(conditional_weights(fit, parameters), 1e-4)
  }
  covariance <- function(theta) {
    out <- theta[["sigma2"]] * matern_correlation(
      distances, phi = theta[["phi"]], kappa = fit$kappa,
      return_sym_matrix = TRUE
    )[location, location]
    if ("tau2" %in% names(theta)) out <- out + theta[["tau2"]] * same_location
    for (j in seq_along(effects)) {
      out <- out + theta[[effect_names[j]]] * effects[[j]]
    }
    measurement <- if ("sigma2_me" %in% names(theta)) theta[["sigma2_me"]]
    diag(out) <- diag(out) + noise + (measurement %||% 0)
    out
  }
  # Central differences on the logarithmic scale.
  derivative <- function(name, step = 1e-4) {
    up <- down <- theta
    up[[name]] <- theta[[name]] * exp(step)
    down[[name]] <- theta[[name]] * exp(-step)
    (covariance(up) - covariance(down)) / (2 * step)
  }
  list(design = design, theta = theta, free = free, covariance = covariance,
       derivatives = lapply(free, derivative))
}

##' Candidate points of the covariance parameters of a working model
##'
##' The estimates, and `n - 1` points of a Halton sequence mapped to a normal
##' distribution on the logarithmic scale of the free parameters, whose
##' covariance is `stretch^2` times the inverse of their expected information,
##' with 0.1 added to its diagonal so that a parameter the data say little
##' about doesn't move too far. No parameter moves by more than a factor of
##' exp(4).
##' @param working A working model, from `working_model()`.
##' @param n The number of points.
##' @return A matrix with a point in each row and a parameter in each column.
##' @noRd
candidate_points <- function(working, n, stretch = 1.4) {
  out <- matrix(working$theta, n, length(working$theta), byrow = TRUE,
                dimnames = list(NULL, names(working$theta)))
  n_free <- length(working$free)
  if (!n_free || n < 2L) return(out[1L, , drop = FALSE])
  precision <- chol2inv(chol(working$covariance(working$theta)))
  scaled <- lapply(working$derivatives, function(derivative) {
    precision %*% derivative
  })
  information <- diag(0.1, n_free)
  for (i in seq_len(n_free)) {
    for (j in seq_len(n_free)) {
      information[i, j] <- information[i, j] +
        sum(scaled[[i]] * t(scaled[[j]])) / 2
    }
  }
  halton <- vapply(first_primes(n_free), function(base) {
    radical_inverse(seq_len(n - 1L), base)
  }, numeric(n - 1L))
  halton <- matrix(halton, n - 1L, n_free)
  steps <- stretch * t(chol(chol2inv(chol(information)))) %*% t(qnorm(halton))
  out[-1L, working$free] <- out[-1L, working$free, drop = FALSE] *
    exp(pmin(pmax(t(steps), -4), 4))
  out
}

##' The first `n` prime numbers
##' @noRd
first_primes <- function(n) {
  primes <- integer(0)
  candidate <- 2L
  while (length(primes) < n) {
    if (all(candidate %% primes[primes <= sqrt(candidate)] != 0L)) {
      primes <- c(primes, candidate)
    }
    candidate <- candidate + 1L
  }
  primes
}

##' The radical inverse of whole numbers in a base, the Halton sequence
##' @noRd
radical_inverse <- function(index, base) {
  vapply(index, function(i) {
    factor <- 1
    value <- 0
    while (i > 0) {
      factor <- factor / base
      value <- value + factor * (i %% base)
      i <- i %/% base
    }
    value
  }, numeric(1))
}

##' PIT values of the test observations of datasets simulated from a working
##' model
##'
##' Each dataset is simulated with the estimates, and takes the candidate point
##' with the largest likelihood profiled over the coefficients, of all its
##' observations or, with `refit`, of each training set. The prediction of the
##' rows that a training set leaves out, given those it has, is in closed form
##' with the precision matrix Q of the candidate point. Their covariance is
##' the inverse of the block of Q at these rows, and the errors of their
##' prediction are that inverse times the block of Q (y - X beta), with the
##' coefficients beta estimated from the observations the parameters are.
##' @param working A working model, from `working_model()`.
##' @param sets A list with the `test` rows of each test set and the rows its
##'   training set leaves out, `left_out`.
##' @param n_simulations The number of simulated datasets.
##' @param n_candidates The number of candidate points.
##' @param refit Whether the parameters are estimated in each training set.
##' @return A matrix with a test observation in each row, in the order of
##'   `sets`, and a simulated dataset in each column.
##' @noRd
candidate_errors <- function(working, sets, n_simulations, n_candidates,
                             refit) {
  design <- working$design
  n <- nrow(design)
  simulated <- crossprod(chol(working$covariance(working$theta)),
                         matrix(rnorm(n * n_simulations), n))
  points <- candidate_points(working, n_candidates)
  size <- lengths(lapply(sets, `[[`, "test"))
  first <- c(0L, cumsum(size))
  uniforms <- matrix(0, sum(size), n_simulations)
  best <- matrix(-Inf, if (refit) length(sets) else 1L, n_simulations)
  for (k in seq_len(nrow(points))) {
    progress_spin()
    root <- tryCatch(chol(working$covariance(points[k, ])),
                     error = function(error) NULL)
    if (is.null(root)) next
    precision <- chol2inv(root)
    log_determinant <- 2 * sum(log(diag(root)))
    QZ <- precision %*% simulated
    QX <- precision %*% design
    XQX <- crossprod(design, QX)
    XQZ <- crossprod(QX, simulated)
    ZQZ <- colSums(simulated * QZ)
    if (!refit) {
      # Q (y - X beta), with beta from all the observations.
      residual <- QZ - QX %*% solve(XQX, XQZ)
      log_likelihood <- -(log_determinant + colSums(simulated * residual)) / 2
      better <- which(log_likelihood > best[1L, ])
      if (!length(better)) next
      best[1L, better] <- log_likelihood[better]
    }
    for (j in seq_along(sets)) {
      progress_spin()
      left_out <- sets[[j]]$left_out
      test <- match(sets[[j]]$test, left_out)
      root_left_out <- chol(precision[left_out, left_out, drop = FALSE])
      conditional <- chol2inv(root_left_out)
      if (refit) {
        # The likelihood of the training set, from the blocks of Q.
        CX <- conditional %*% QX[left_out, , drop = FALSE]
        CZ <- conditional %*% QZ[left_out, , drop = FALSE]
        information <- XQX - crossprod(QX[left_out, , drop = FALSE], CX)
        score <- XQZ - crossprod(QX[left_out, , drop = FALSE], CZ)
        beta <- solve(information, score)
        log_likelihood <- -(log_determinant +
                              2 * sum(log(diag(root_left_out))) + ZQZ -
                              colSums(QZ[left_out, , drop = FALSE] * CZ) -
                              colSums(score * beta)) / 2
        better <- which(log_likelihood > best[j, ])
        if (!length(better)) next
        best[j, better] <- log_likelihood[better]
        errors <- CZ[test, better, drop = FALSE] -
          CX[test, , drop = FALSE] %*% beta[, better, drop = FALSE]
      } else {
        errors <- conditional[test, , drop = FALSE] %*%
          residual[left_out, better, drop = FALSE]
      }
      uniforms[first[j] + seq_along(test), better] <-
        pnorm(errors / sqrt(diag(conditional)[test]))
    }
  }
  uniforms
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
