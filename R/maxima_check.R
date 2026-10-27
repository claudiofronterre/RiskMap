##' Gauss-Hermite rule for the weight exp(-x^2)
##'
##' Nodes and weights by the Golub-Welsch algorithm, from the eigenvalues of
##' the Jacobi matrix and the first components of its eigenvectors.
##' @noRd
gauss_hermite_rule <- function(nodes) {
  off_diagonal <- sqrt(seq_len(nodes - 1L) / 2)
  jacobi <- matrix(0, nodes, nodes)
  jacobi[cbind(seq_len(nodes - 1L), 2:nodes)] <- off_diagonal
  jacobi[cbind(2:nodes, seq_len(nodes - 1L))] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(x = decomposition$values,
       w = sqrt(pi) * decomposition$vectors[1L, ]^2)
}

##' Log-likelihood of the model with independent location effects
##'
##' As the spatial scale goes to zero the spatial correlation between distinct
##' locations vanishes, and the likelihood tends to that of independent
##' location effects with the total variance at a location (the spatial
##' variance plus any nugget), a product of one-dimensional integrals. Each is
##' evaluated by adaptive Gauss-Hermite quadrature at the posterior mode of
##' the location's effect. The observation terms are those of the MCML engine
##' for the canonical links, so the result is on the engine's scale, with the
##' same constants. `D` is the engine's design matrix, `location` maps each
##' observation to its location, and `offset` is the covariate offset.
##'
##' The returned `log_likelihood(beta, log_variance)` gives the value, its
##' gradient in (beta, log_variance), and the posterior mean and variance of
##' each location's effect, which the Jensen bound
##' (`compute_jensen_bound()`) needs. `maximise()` maximises the value from
##' several starting variances and records the change in the maximum from 30
##' to 60 nodes as its numerical error, keeping the 30-node value and moments
##' so that the error of the Jensen bound can be judged the same way.
##' @noRd
make_independent_site_likelihood <- function(D, y, units_m, location, offset,
                                             family, nodes = 60L) {
  poisson <- identical(family, "poisson")
  location <- as.integer(location)
  n_locations <- max(location)
  single <- identical(location, seq_len(length(y)))
  by_location <- function(values) {
    if (single) return(values)
    rowsum(values, location, reorder = TRUE)
  }
  terms <- function(eta) {
    if (poisson) {
      mean <- exp(eta)
      list(value = y * eta - units_m * mean, score = y - units_m * mean,
           curvature = units_m * mean)
    } else {
      probability <- stats::plogis(eta)
      softplus <- pmax(eta, 0) + log1p(exp(-abs(eta)))
      list(value = y * eta - units_m * softplus,
           score = y - units_m * probability,
           curvature = units_m * probability * (1 - probability))
    }
  }
  rules <- list()
  rule <- function(count) {
    key <- as.character(count)
    if (is.null(rules[[key]])) rules[[key]] <<- gauss_hermite_rule(count)
    rules[[key]]
  }

  # Posterior mode of every location's effect by Newton steps, each halved
  # where it does not increase that location's log-density.
  find_modes <- function(linear, variance) {
    density_at <- function(effects) {
      drop(by_location(terms(linear + effects[location])$value)) -
        effects^2 / (2 * variance)
    }
    effects <- numeric(n_locations)
    current <- density_at(effects)
    for (iteration in seq_len(200L)) {
      observed <- terms(linear + effects[location])
      gradient <- drop(by_location(observed$score)) - effects / variance
      curvature <- drop(by_location(observed$curvature)) + 1 / variance
      if (max(abs(gradient)) < 1e-9) break
      step <- gradient / curvature
      for (halving in seq_len(50L)) {
        proposed <- density_at(effects + step)
        worse <- !(proposed >= current - 1e-12)
        if (!any(worse)) break
        step[worse] <- step[worse] / 2
      }
      effects <- effects + step
      current <- density_at(effects)
    }
    curvature <- drop(by_location(terms(linear + effects[location])$curvature)) +
      1 / variance
    list(mode = effects, scale = 1 / sqrt(curvature))
  }

  log_likelihood <- function(beta, log_variance, count = nodes) {
    variance <- exp(log_variance)
    linear <- drop(D %*% beta) + offset
    modes <- find_modes(linear, variance)
    gh <- rule(count)
    # Location effects at the nodes, one row per location.
    effects <- modes$mode + sqrt(2) * outer(modes$scale, gh$x)
    observed <- terms(linear + effects[location, , drop = FALSE])
    log_integrand <- by_location(observed$value) -
      effects^2 / (2 * variance) +
      rep(log(gh$w) + gh$x^2, each = n_locations)
    largest <- apply(log_integrand, 1L, max)
    weights <- exp(log_integrand - largest)
    totals <- rowSums(weights)
    posterior <- weights / totals
    value <- sum(largest + log(totals) + log(sqrt(2) * modes$scale)) -
      n_locations / 2 * log(2 * pi * variance)
    mean <- rowSums(posterior * effects)
    second <- rowSums(posterior * effects^2)
    score <- rowSums(observed$score * posterior[location, , drop = FALSE])
    list(value = value,
         gradient = c(drop(crossprod(D, score)),
                      sum(second) / (2 * variance) - n_locations / 2),
         mean = mean, variance = pmax(second - mean^2, 0))
  }

  maximise <- function(beta_start, log_variance_starts = c(-2, 0, 2),
                       lower = -20, upper = 20) {
    p <- length(beta_start)
    # optim() asks for the value and the gradient at the same point in turn,
    # so the last evaluation is kept.
    last <- NULL
    evaluate <- function(parameters) {
      if (is.null(last) || !identical(last$parameters, parameters)) {
        last <<- list(parameters = parameters,
                      result = log_likelihood(parameters[seq_len(p)],
                                              parameters[p + 1L]))
      }
      last$result
    }
    fits <- lapply(log_variance_starts, function(log_variance) {
      negative <- function(parameters) {
        value <- evaluate(parameters)$value
        if (is.finite(value)) -value else 1e100
      }
      gradient <- function(parameters) -evaluate(parameters)$gradient
      tryCatch(
        stats::optim(c(beta_start, log_variance), negative, gradient,
                     method = "L-BFGS-B",
                     lower = c(rep(-30, p), lower),
                     upper = c(rep(30, p), upper),
                     control = list(maxit = 500L, factr = 1e5)),
        error = function(error) NULL
      )
    })
    fits <- Filter(function(fit) !is.null(fit) && fit$value < 1e99, fits)
    if (!length(fits)) return(NULL)
    best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
    beta <- best$par[seq_len(p)]
    log_variance <- best$par[p + 1L]
    at_maximum <- log_likelihood(beta, log_variance)
    coarse <- log_likelihood(beta, log_variance, count = 30L)
    list(beta = beta, log_variance = log_variance,
         log_likelihood = at_maximum$value,
         quadrature_error = abs(at_maximum$value - coarse$value),
         mean = at_maximum$mean, variance = at_maximum$variance,
         coarse = list(log_likelihood = coarse$value, mean = coarse$mean,
                       variance = coarse$variance),
         convergence = best$convergence)
  }

  list(log_likelihood = log_likelihood, maximise = maximise)
}

##' Jensen lower bound on the log-likelihood at a positive spatial scale
##'
##' For the location effects with prior covariance `covariance`, the
##' log-likelihood is at least
##' l0 - log|C| / 2 - tr{(C^-1 - I)(V + m m')} / (2 v), with C the covariance
##' divided by v, where l0 is the log-likelihood with independent location
##' effects of variance v at the same coefficients and m and V are the
##' posterior means and variances of those effects. It is the evidence lower
##' bound with that posterior as the approximating distribution. Returns NA
##' if the covariance is not positive definite.
##' @noRd
compute_jensen_bound <- function(independent, covariance) {
  variance <- exp(independent$log_variance)
  root <- tryCatch(chol(covariance / variance), error = function(error) NULL)
  if (is.null(root)) return(NA_real_)
  inverse <- chol2inv(root)
  mean <- independent$mean
  trace <- sum(diag(inverse) * independent$variance) +
    drop(crossprod(mean, inverse %*% mean)) -
    sum(independent$variance) - sum(mean^2)
  independent$log_likelihood - sum(log(diag(root))) - trace / (2 * variance)
}

##' Settings of the screen of short scales and the exact alarms
##'
##' The screen and the alarms rest on the likelihood with independent location
##' effects (`make_independent_site_likelihood()`), which is computed once per
##' survey, when first needed, and kept. They are built only for the
##' canonical links and for models without unstructured random effects, the
##' cases in which that likelihood is a product of one-dimensional integrals
##' computed exactly. The screened scales are 1/8, 1/4, 1/2 and 1 times the
##' median distance between neighbouring locations (`spacing`), all on the
##' internal distance scale. The independent-site likelihood is maximised from
##' the coefficients of the fit that the check follows.
##' @noRd
make_short_scale_screen <- function(engine, independent_sites, distances,
                                    kappa = 0.5, messages = FALSE,
                                    distance_scale = 1) {
  layout <- engine$layout
  spacing <- median_nearest_neighbour_distance(distances)
  closest <- min(distances[distances > 0])
  cache <- new.env(parent = emptyenv())
  independent <- function(beta_start) {
    if (!exists("fit", envir = cache, inherits = FALSE)) {
      started <- proc.time()[["elapsed"]]
      fit <- tryCatch(independent_sites$maximise(beta_start),
                      error = function(error) NULL)
      if (!is.null(fit)) fit$seconds <- proc.time()[["elapsed"]] - started
      assign("fit", fit, envir = cache)
    }
    get("fit", envir = cache)
  }
  # The point at a screened scale: the coefficients and total variance of the
  # independent-site maximum, the other parameters as fitted. With a nugget
  # the total variance at a location is the spatial variance plus the nugget,
  # and the spatial variance is set so that the total matches.
  plug_in <- function(estimate, fit, log_phi) {
    point <- estimate
    point[layout$index$beta] <- fit$beta
    point[layout$index$phi] <- log_phi
    values <- engine$unpack(estimate)
    total <- exp(fit$log_variance)
    sigma2 <- if (!is.null(layout$index$nu2)) {
      total / (1 + values$nugget_ratio)
    } else {
      total - values$sigma2 * values$nugget_ratio
    }
    point[layout$index$sigma2] <- if (sigma2 > 0) log(sigma2) else -20
    point
  }
  announce <- function(event, log_phi = NULL) {
    if (!messages) return(invisible(NULL))
    message(switch(
      event,
      screen = paste0(
        "The likelihood at a spatial scale shorter than the fitted one (phi = ",
        signif(exp(log_phi) * distance_scale, 3), ") came close to the ",
        "fitted maximum, so the fit is repeated from that point to check for ",
        "a higher maximum at short scales."
      ),
      alarm = paste0(
        "An exact lower bound on the likelihood at a short spatial scale is ",
        "above the fitted maximum, so the fit is repeated from that point."
      ),
      failed = paste0(
        "The likelihood with independent location effects could not be ",
        "computed, so the fit is repeated from short-scale starting values ",
        "to check for a higher maximum."
      )
    ))
  }
  list(
    beta_index = layout$index$beta, sigma2_index = layout$index$sigma2,
    phi_index = layout$index$phi,
    log_phi = log(spacing) + log(c(1 / 8, 1 / 4, 1 / 2, 1)),
    margin = 0.25, window = 2, step = 0.5, step_limit = 2,
    # A numerical comparison tolerance: differences from an exact value
    # smaller than this are not acted on. It is separate from the Monte Carlo
    # and quadrature errors.
    alarm_multiplier = 2.58, comparison_tolerance = 1e-4,
    independent = independent, plug_in = plug_in,
    # The Jensen bound with its numerical error, the change from 30 to 60
    # quadrature nodes, which enters through the posterior moments as well as
    # the limit.
    bound = function(point, fit) {
      covariance <- engine$spatial_covariance(point)
      fine <- compute_jensen_bound(fit, covariance)
      coarse <- if (is.null(fit$coarse)) fine else
        compute_jensen_bound(utils::modifyList(fit, fit$coarse), covariance)
      c(value = fine, error = abs(fine - coarse))
    },
    # A scale at which the correlation between the two closest distinct
    # locations is below 1e-6 makes the spatial correlation numerically
    # negligible. It shows the fit is numerically independent, not that the
    # maximum lies at zero scale.
    at_limit = function(log_phi) {
      scaled <- closest / exp(log_phi)
      correlation <- if (scaled > 700) 0 else
        2^(1 - kappa) / gamma(kappa) * scaled^kappa * besselK(scaled, kappa)
      !is.finite(correlation) || correlation < 1e-6
    },
    announce = announce
  )
}

##' Evaluate the short scales after a fit at a long scale
##'
##' At each screened scale the estimated log-likelihood on the fitting bank is
##' evaluated at the point from the independent-site maximum
##' (`make_short_scale_screen()`). The gate opens, and the fit is repeated
##' from the best screened point, when a value comes within `margin` of the
##' fitted value. Where the gate is still closed, each value within `window`
##' of it is improved by one three-point step in the log spatial variance
##' (`step_log_variance()`), keeping the better of the plug-in and the step,
##' and the gate is judged again. A failed evaluation, of a plug-in point or
##' within a variance step, opens the gate. The margins are search settings,
##' not significance levels.
##' @return The table of screened scales, whether the gate `opened`, whether
##' an evaluation `failed`, the `start` for the refit, and the number of
##' evaluations.
##' @noRd
screen_short_scales <- function(objective, fit_bank, estimate, fitted_value,
                                fit, screen, lower, upper) {
  evaluations <- 0L
  value_at <- function(point) {
    evaluations <<- evaluations + 1L
    value <- tryCatch(objective(point, fit_bank), error = function(error) NA)
    if (length(value) == 1L && is.finite(value)) value else NA_real_
  }
  gate <- fitted_value - screen$margin
  points <- lapply(screen$log_phi, function(log_phi) {
    pmin(pmax(screen$plug_in(estimate, fit, log_phi), lower), upper)
  })
  plug_in <- vapply(points, value_at, numeric(1))
  failed <- anyNA(plug_in)
  best_values <- plug_in
  best_points <- points
  stepped <- rep(NA_real_, length(points))
  opened <- failed || any(plug_in >= gate, na.rm = TRUE)
  if (!opened) {
    for (near in which(plug_in >= gate - screen$window)) {
      step <- step_log_variance(value_at, points[[near]], plug_in[near],
                                screen$sigma2_index, screen$step,
                                screen$step_limit, lower, upper)
      stepped[near] <- best_values[near] <- step$value
      best_points[[near]] <- step$point
      failed <- failed || step$failed
    }
    opened <- failed || any(best_values >= gate, na.rm = TRUE)
  }
  best <- which.max(best_values)
  # With no finite value the refit starts at a quarter of the spacing.
  start <- if (length(best)) best_points[[best]] else points[[2L]]
  list(table = data.frame(log_phi = screen$log_phi, plug_in = plug_in,
                          stepped = stepped,
                          gap = best_values - fitted_value),
       opened = opened, failed = failed, start = start,
       evaluations = evaluations)
}

##' One three-point step in the log spatial variance
##'
##' Evaluates the point with the log variance moved by plus and minus
##' `half_width`, fits a parabola through the three values and evaluates its
##' vertex, moved at most `limit` (or, if the parabola is convex, a step of
##' twice `half_width` uphill). Returns the best value and point, the
##' starting point included, so the step never makes the screen worse, and
##' whether an evaluation failed.
##' @noRd
step_log_variance <- function(value_at, point, centre_value, index,
                              half_width, limit, lower, upper) {
  moved <- function(change) {
    point[index] <- min(max(point[index] + change, lower[index]),
                        upper[index])
    point
  }
  points <- list(point, moved(-half_width), moved(half_width))
  values <- c(centre_value, value_at(points[[2L]]), value_at(points[[3L]]))
  if (all(is.finite(values))) {
    curvature <- values[3L] + values[2L] - 2 * values[1L]
    change <- if (curvature < 0) {
      half_width * (values[2L] - values[3L]) / (2 * curvature)
    } else {
      sign(values[3L] - values[2L]) * 2 * half_width
    }
    change <- max(-limit, min(limit, change))
    if (abs(change) > 1e-8 && abs(abs(change) - half_width) > 1e-8) {
      points <- c(points, list(moved(change)))
      values <- c(values, value_at(points[[4L]]))
    }
  }
  best <- which.max(values)
  list(value = values[best], point = points[[best]], failed = anyNA(values))
}

##' Estimated log-likelihood on the comparison bank, with its standard error
##'
##' The bank of `n_draws` draws from `seed` is the one on which separate
##' maxima are compared (`adjudicate_maxima()`), drawn before the correction
##' and certification banks. The standard error is the delta-method one,
##' sd(w / mean(w)) / sqrt(n).
##' @noRd
estimate_on_comparison_bank <- function(engine, penalty, parameters, n_draws,
                                        seed) {
  set.seed(seed)
  bank <- matrix(stats::rnorm(n_draws * engine$latent_dimension),
                 nrow = n_draws)
  log_weights <- engine$evaluate(parameters, bank, log_weights_only = TRUE)
  if (!any(is.finite(log_weights))) return(list(value = -Inf, se = Inf))
  largest <- max(log_weights[is.finite(log_weights)])
  relative <- exp(log_weights - largest)
  relative[!is.finite(relative)] <- 0
  list(value = largest + log(mean(relative)) +
         if (is.null(penalty)) 0 else penalty(parameters),
       se = stats::sd(relative / mean(relative)) / sqrt(n_draws))
}

##' Exact alarm that a fit is not the global maximum
##'
##' The likelihood with independent location effects is the limit of the
##' likelihood as the scale goes to zero, so if it exceeds the fitted
##' log-likelihood the fit is not the global maximum. Jensen bounds at the
##' screened scales (`compute_jensen_bound()`) are exact lower bounds on the
##' likelihood there. The largest of these exact values is compared with the
##' fitted log-likelihood estimated on the comparison bank, and the alarm
##' fires when it exceeds that estimate by more than `alarm_multiplier`
##' standard errors plus its own quadrature error (for a bound, the change in
##' the bound from 30 to 60 nodes) and `comparison_tolerance`, a numerical
##' comparison tolerance. With a penalty the limit of
##' zero scale is left out, because the penalty diverges there, and the
##' penalty is added to the bounds. The `start` is the point that gave the
##' largest value, at the shortest screened scale for the limit itself.
##' @noRd
judge_exact_alarm <- function(engine, penalty, screen, fit, estimate,
                              n_draws, seed) {
  points <- lapply(screen$log_phi, function(log_phi) {
    screen$plug_in(estimate, fit, log_phi)
  })
  bounds <- vapply(points, function(point) {
    bound <- screen$bound(point, fit)
    bound[["value"]] <- bound[["value"]] +
      if (is.null(penalty)) 0 else penalty(point)
    bound
  }, numeric(2))
  values <- c(if (is.null(penalty)) fit$log_likelihood, bounds["value", ])
  errors <- c(if (is.null(penalty)) fit$quadrature_error, bounds["error", ])
  sources <- c(if (is.null(penalty)) -Inf, screen$log_phi)
  starts <- c(if (is.null(penalty)) points[1L], points)
  if (!any(is.finite(values))) return(NULL)
  top <- which.max(values)
  fitted <- estimate_on_comparison_bank(engine, penalty, estimate, n_draws,
                                        seed)
  list(statistic = values[top], log_phi = sources[top],
       quadrature_error = errors[top],
       fitted_comparison = fitted$value, fitted_comparison_se = fitted$se,
       fired = values[top] - errors[top] - screen$comparison_tolerance >
         fitted$value + screen$alarm_multiplier * fitted$se,
       start = starts[[top]])
}

##' Check for other maxima on the fitting bank
##'
##' A fitted spatial scale below the trigger of `limits$dispersed` (see
##' `make_dispersed_starts()`) is repeated from the dispersed starts, with the
##' other parameters at their fitted values. A fit at a longer scale is
##' followed by the screen of short scales (`screen_short_scales()`), which
##' repeats the fit from the best screened point when its gate opens. In both
##' cases an exact alarm (`judge_exact_alarm()`) repeats the fit from the
##' point that raised it. A fit with the spatial variance at zero is checked
##' only by the alarm. Without the screen settings (`limits$dispersed$screen`
##' is NULL) only the dispersed starts are used, and if the likelihood with
##' independent location effects cannot be computed a fit at a long scale is
##' repeated from the two dispersed starts at the shortest scale.
##' `choose_maximum()` decides which maximum to keep, on the independent bank
##' of `n_draws` draws from `seed`, or on a bank four times larger from
##' `enlarged_seed` when the first is ambiguous. That bank serves every later
##' comparison. The kept fit is numerically independent
##' (`numerically_independent` is TRUE) when its spatial correlation is
##' numerically negligible, which leaves the spatial scale unresolved.
##' A kept maximum at a positive scale below the shortest screened one whose
##' log-likelihood is not clearly above the limit at this Monte Carlo
##' precision stays a candidate like any other, and the record notes it as
##' `near_limit` for a warning. An alarm persists if its exact value still
##' exceeds the kept maximum, whatever the comparison with the limit says.
##' @return The kept `optimised` result, a `check` table of the candidates
##' (NULL if no refit was run), whether the kept maximum is `separate` from
##' the fitted one, whether the choice is `ambiguous`, the record of the
##' `screen` and alarm (NULL without screen settings) and
##' `numerically_independent`.
##' @noRd
check_other_maxima <- function(engine, objective, penalty, optimised,
                               fit_bank, lower, upper, limits, n_draws,
                               seed, enlarged_seed = NULL) {
  dispersed <- limits$dispersed
  screen <- dispersed[["screen"]]
  estimate <- optimised$fit$par
  at_zero_variance <- dispersed$sigma2_index %in% optimised$boundary
  short_fit <- !at_zero_variance &&
    estimate[dispersed$phi_index] < dispersed$trigger
  starts <- list()
  dispersed_start <- function(which) {
    start <- estimate
    start[dispersed$sigma2_index] <- dispersed$log_sigma2[which]
    start[dispersed$phi_index] <- dispersed$log_phi[which]
    start
  }
  if (short_fit) {
    if (!is.null(dispersed$announce)) {
      dispersed$announce(estimate[[dispersed$phi_index]])
    }
    starts <- lapply(seq_along(dispersed$labels), dispersed_start)
    names(starts) <- dispersed$labels
  }
  record <- NULL
  if (is.null(screen) && !is.null(dispersed$screen_unavailable)) {
    record <- list(gate = "unavailable", reason = dispersed$screen_unavailable)
  }
  if (!is.null(screen)) {
    started <- proc.time()[["elapsed"]]
    fit <- screen$independent(estimate[screen$beta_index])
    record <- list(independent = fit[c("beta", "log_variance",
                                        "log_likelihood", "quadrature_error",
                                        "convergence", "seconds")],
                   gate = "not run")
    if (is.null(fit)) {
      record$gate <- "failed"
      if (!short_fit && !at_zero_variance) {
        screen$announce("failed")
        shortest <- which(dispersed$log_phi == min(dispersed$log_phi))
        starts <- lapply(shortest, dispersed_start)
        names(starts) <- dispersed$labels[shortest]
      }
    } else {
      if (!short_fit && !at_zero_variance) {
        screened <- screen_short_scales(objective, fit_bank, estimate,
                                        -optimised$fit$objective, fit,
                                        screen, lower, upper)
        record$screen <- screened$table
        record$screen_evaluations <- screened$evaluations
        record$gate <- if (screened$opened) "opened" else "closed"
        record$failed_evaluation <- screened$failed
        if (screened$opened) {
          screen$announce("screen", screened$start[[dispersed$phi_index]])
          starts[["short-scale screen"]] <- screened$start
        }
      }
      alarm <- judge_exact_alarm(engine, penalty, screen, fit, estimate,
                                 n_draws, seed)
      if (!is.null(alarm)) {
        record$alarm <- alarm[setdiff(names(alarm), "start")]
        if (alarm$fired) {
          screen$announce("alarm")
          starts[["exact alarm"]] <- alarm$start
        }
      }
    }
    record$seconds <- proc.time()[["elapsed"]] - started
  }
  if (!length(starts)) {
    return(list(optimised = optimised, check = NULL, separate = FALSE,
                ambiguous = FALSE, screen = record,
                numerically_independent = FALSE))
  }
  refits <- lapply(starts, function(start) {
    optimise_on_bank(objective, pmin(pmax(start, lower), upper), fit_bank,
                     lower, upper, TRUE, limits)
  })
  # A refit that ends with non-finite values or the penalty value is not a
  # maximum and is left out of the comparison.
  usable <- vapply(refits, function(refit) is_valid_fit(refit$fit),
                   logical(1))
  candidates <- c(list(optimised), unname(refits[usable]))
  information <- tryCatch(
    compute_central_hessian(function(parameters) {
      -objective(parameters, fit_bank)
    }, estimate),
    error = function(error) NULL
  )
  choice <- choose_maximum(
    engine, penalty, lapply(candidates, function(candidate) {
      candidate$fit$par
    }),
    vapply(candidates, function(candidate) -candidate$fit$objective,
           numeric(1)),
    information, n_draws, seed, limits$tolerance, enlarged_seed
  )
  choice$table$start <- c("fitted", names(starts)[usable])
  choice$table$log_phi <- vapply(candidates, function(candidate) {
    candidate$fit$par[[dispersed$phi_index]]
  }, numeric(1))
  kept <- candidates[[choice$chosen]]
  # Every comparison after the choice uses the bank that decided it, the
  # enlarged one when it was drawn.
  enlarged <- any(choice$table$independent_draws > n_draws, na.rm = TRUE)
  judged <- NULL
  judge_kept <- function() {
    if (is.null(judged)) {
      judged <<- if (!enlarged && choice$chosen == 1L &&
                     is.finite(record$alarm$fitted_comparison %||% NA)) {
        list(value = record$alarm$fitted_comparison,
             se = record$alarm$fitted_comparison_se, draws = n_draws)
      } else {
        draws <- if (enlarged) 4L * n_draws else n_draws
        c(estimate_on_comparison_bank(engine, penalty, kept$fit$par, draws,
                                      if (enlarged) enlarged_seed else seed),
          list(draws = draws))
      }
    }
    judged
  }
  scale <- kept$fit$par[[dispersed$phi_index]]
  at_zero_variance <- dispersed$sigma2_index %in% kept$boundary
  # A fit whose spatial correlation is numerically negligible is reported as
  # numerically independent, with the scale unresolved.
  independent <- !is.null(screen$at_limit) && !at_zero_variance &&
    screen$at_limit(scale)
  limit <- record$independent
  if (!independent && !is.null(limit$log_likelihood) && is.null(penalty) &&
      !at_zero_variance && scale < min(screen$log_phi)) {
    # A positive short scale whose likelihood is not clearly above the limit
    # at this Monte Carlo precision. The candidate is kept and assessed as
    # usual, with a warning.
    gap <- judge_kept()$value - limit$log_likelihood
    if (gap <= limit$quadrature_error + screen$comparison_tolerance +
        screen$alarm_multiplier * judge_kept()$se) {
      record$near_limit <- list(gap = gap, se = judge_kept()$se,
                                draws = judge_kept()$draws, log_phi = scale)
    }
  }
  if (isTRUE(record$alarm$fired)) {
    # The alarm persists if the exact value still exceeds the kept maximum.
    record$alarm$persists <- record$alarm$statistic -
      record$alarm$quadrature_error - screen$comparison_tolerance >
      judge_kept()$value + screen$alarm_multiplier * judge_kept()$se
  }
  list(optimised = kept, check = choice$table, separate = choice$separate,
       ambiguous = choice$ambiguous, screen = record,
       numerically_independent = independent)
}
