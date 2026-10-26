##' Profile likelihood of a fitted geostatistical model
##'
##' Computes profile likelihoods of a binomial or Poisson model fitted by
##' [glgpm()] with certified Monte Carlo maximum likelihood: for each chosen
##' parameter, the log-likelihood maximised over all the other parameters on a
##' grid of values, and for each chosen pair, the same over a grid of pairs.
##' The Monte Carlo likelihood is evaluated on one fixed bank of `n_draws`
##' standard normal draws with the importance distribution anchored at the
##' estimate, so the profiles are smooth functions of the parameters. Use
##' [plot.RiskMap_profile()] to draw them.
##'
##' Grids are centred at the estimate on the scale on which the model is fitted
##' (the logarithm for variances and the scale) and span `span` standard errors
##' either way. If the profile deviance has not passed the 99% point of a
##' chi-squared distribution on one degree of freedom by the end of a marginal
##' grid, the grid is extended with steps that double each time, up to
##' `n_points` points on each side. Each point starts the optimisation from
##' the solution at its neighbour. Profiles can have more than one constrained
##' maximum in sparse data, so a point whose deviance is beyond the interval is
##' also solved from the estimate, and the higher solution is kept.
##'
##' Profile intervals are where the profile deviance, twice the drop in the
##' log-likelihood from its maximum, crosses the chi-squared quantile for the
##' requested level. The crossing is interpolated linearly in the signed root
##' of the deviance, which is close to linear in the parameter, after one more
##' point is solved at the interpolated crossing when the two points around it
##' are far apart. An end that the profile does not reach within the grid or
##' the bounds of the parameter is reported as open (`NA`). If a point on a
##' profile is higher than the estimate, the estimate is not the highest
##' maximum found, and a message says so.
##'
##' Each point is a maximisation of the Monte Carlo likelihood over the other
##' parameters, so the time grows with the number of points and of parameters.
##' A model that takes a minute to fit can take tens of seconds per point.
##' Joint profiles cost `n_grid^2` points each and are computed only when
##' `pairs` asks for them. With `messages = TRUE` the time left is estimated
##' after the first parameter.
##'
##' @param fitted A `RiskMap` object from [glgpm()] for a binomial or Poisson
##'   model.
##' @param parameters Names of the parameters to profile one at a time, as in
##'   `names(unlist(fitted$estimate))` (the prefix `beta.` can be left out).
##'   The default is every parameter.
##' @param pairs A list of pairs of parameter names for joint profiles, for
##'   example `list(c("sigma2", "phi"))`. The default is none.
##' @param n_points Number of points of each marginal grid before any
##'   extension, and the most points on each side of the estimate.
##' @param n_grid Number of points along each side of a joint grid.
##' @param span Half-width of the grids, in standard errors.
##' @param n_draws Number of Monte Carlo draws in the bank.
##' @param seed Seed of the bank. The random number state of the session is
##'   restored afterwards.
##' @param level Confidence level of the reported profile intervals.
##' @param messages Logical, whether to report progress.
##' @param ... Not used.
##' @return An object of class `RiskMap_profile`: a list with the marginal
##'   profiles (`marginal`), the joint profiles (`joint`), the intervals
##'   (`intervals`), the estimate and the maximum of the log-likelihood on the
##'   bank (`estimate`, `max_loglik`) and the settings.
##' @examples
##' \donttest{
##' data(liberia)
##' liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
##' fit <- glgpm(npos ~ log(elevation) + gp(), denominator = ntest,
##'              data = liberia, model_crs = 32629, family = "binomial",
##'              control_mcmc = set_control_mcmc(seed = 2026))
##' profiles <- profile(fit, parameters = c("sigma2", "phi"),
##'                     pairs = list(c("sigma2", "phi")), n_grid = 9,
##'                     n_draws = 2000)
##' profiles
##' plot(profiles)
##' }
##' @seealso [plot.RiskMap_profile()]
##' @importFrom stats profile
##' @export
profile.RiskMap <- function(fitted, parameters = NULL, pairs = list(),
                            n_points = 11L, n_grid = 11L, span = 3,
                            n_draws = 5000L, seed = 1L, level = 0.95,
                            messages = TRUE, ...) {
  stopifnot(inherits(fitted, "RiskMap"))
  if (!fitted$family %in% c("binomial", "poisson") ||
      is.null(fitted$mcml_history)) {
    stop("Profiles are available for binomial and Poisson models fitted by ",
         "certified Monte Carlo maximum likelihood.", call. = FALSE)
  }
  check_profile_settings(n_points, n_grid, span, n_draws, level)
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  if (messages) {
    message("Re-maximising the likelihood on a bank of ", n_draws, " draws...")
  }
  problem <- make_profile_problem(fitted, n_draws, seed)
  parameters <- match_profile_parameters(problem, parameters %||% problem$names)
  pairs <- lapply(pairs %||% list(), function(pair) {
    match_profile_parameters(problem, pair)
  })

  marginal <- vector("list", length(parameters))
  names(marginal) <- parameters
  for (k in seq_along(parameters)) {
    name <- parameters[k]
    if (messages) {
      message("Profiling ", problem$labels[[name]], " (", k, " of ",
              length(parameters), ")...")
    }
    started <- proc.time()[["elapsed"]]
    marginal[[name]] <- profile_one_parameter(problem, name, n_points, span,
                                              level)
    if (messages && k == 1L &&
        (length(parameters) > 1L || length(pairs) > 0L)) {
      seconds <- proc.time()[["elapsed"]] - started
      per_point <- seconds / max(1L, nrow(marginal[[name]]) - 1L)
      planned <- (length(parameters) - 1L) * (n_points + 1L) +
        length(pairs) * n_grid^2
      message("About ", format_profile_time(per_point * planned),
              " left, at about ", signif(per_point, 2), " s per point.")
    }
  }
  joint <- lapply(pairs, function(pair) {
    if (messages) {
      message("Profiling ", problem$labels[[pair[1]]], " with ",
              problem$labels[[pair[2]]], "...")
    }
    profile_two_parameters(problem, pair, n_grid, span, marginal)
  })
  names(joint) <- vapply(pairs, paste, "", collapse = ":")

  values <- c(unlist(lapply(marginal, `[[`, "loglik")),
              unlist(lapply(joint, `[[`, "loglik")))
  higher <- max(values, na.rm = TRUE) - problem$max_loglik
  if (messages && higher > 0.01) {
    message("A point on the profile is higher than the estimate by ",
            signif(higher, 3), " log-likelihood units: the estimate is not ",
            "the highest maximum found.")
  }
  max_loglik <- max(problem$max_loglik, values, na.rm = TRUE)
  marginal <- lapply(marginal, add_profile_deviance, max_loglik)
  joint <- lapply(joint, add_profile_deviance, max_loglik)

  structure(list(
    marginal = marginal, joint = joint,
    intervals = profile_intervals(problem, marginal, level),
    estimate = problem$estimate, standard_errors = problem$standard_errors,
    labels = problem$labels, max_loglik = max_loglik,
    higher_than_estimate = higher, level = level,
    settings = list(n_draws = n_draws, seed = seed, span = span)
  ), class = "RiskMap_profile")
}

##' A duration in seconds as readable text
##' @noRd
format_profile_time <- function(seconds) {
  if (seconds < 90) return(paste(max(1, round(seconds)), "seconds"))
  if (seconds < 5400) return(paste(round(seconds / 60), "minutes"))
  paste(signif(seconds / 3600, 2), "hours")
}

##' Check the settings of a profile
##' @noRd
check_profile_settings <- function(n_points, n_grid, span, n_draws, level) {
  stopifnot(
    "'n_points' must be at least 5" = is.numeric(n_points) && n_points >= 5,
    "'n_grid' must be at least 3" = is.numeric(n_grid) && n_grid >= 3,
    "'span' must be positive" = is.numeric(span) && span > 0,
    "'n_draws' must be at least 100" = is.numeric(n_draws) && n_draws >= 100,
    "'level' must be between 0 and 1" = is.numeric(level) && level > 0 &&
      level < 1
  )
  invisible(TRUE)
}

##' The Monte Carlo likelihood of a fitted model, rebuilt for profiling
##'
##' Rebuilds the engine of the fit from the data and settings kept in the fitted
##' object (the standardised design, the coordinates rescaled by the largest
##' distance, the latent structure and the Monte Carlo controls), draws one
##' bank, anchors the global proposal at the estimate and re-maximises the
##' log-likelihood on that bank. Parameters are handled on the working scale of
##' the user (coefficients, log variance, log scale in distance units), and
##' converted to the internal scale of the engine for each evaluation.
##' @noRd
make_profile_problem <- function(fitted, n_draws, seed) {
  design_scaling <- standardize_design_matrix(fitted$D)
  spatial_scaling <- scale_spatial_coordinates(fitted$coords)
  distances <- as.matrix(pairwise_distances(spatial_scaling$coordinates))
  control_mcml <- attr(fitted, "control_mcml")
  engine <- make_adaptive_mcml_engine(
    D = design_scaling$design, y = fitted$y, units_m = fitted$units_m,
    distances = distances, family = fitted$family, kappa = fitted$kappa,
    cov_offset = fitted$cov_offset, ID_coords = fitted$ID_coords,
    fix_tau2 = fitted$fix_tau2,
    ID_re = if (is.null(fitted$ID_re)) NULL else as.matrix(fitted$ID_re),
    defensive_weight = control_mcml$defensive_weight %||% 0.05,
    proposal = control_mcml$proposal %||% "global"
  )
  objective <- make_mcml_objective(engine, control_mcml, distances,
                                   spatial_scaling$distance_scale)$objective
  p <- ncol(fitted$D)
  transform <- design_scaling$coefficient_transform
  log_scale <- log(spatial_scaling$distance_scale)
  to_internal <- function(working) {
    c(solve(transform, working[seq_len(p)]), working[p + 1L],
      working[p + 2L] - log_scale, working[-seq_len(p + 2L)])
  }
  estimate <- unlist(fitted$estimate)
  internal <- to_internal(estimate)
  engine$set_anchor(internal)
  set.seed(seed)
  bank <- matrix(stats::rnorm(n_draws * engine$latent_dimension), nrow = n_draws)
  value <- function(working) {
    result <- tryCatch(objective(to_internal(working), bank),
                       error = function(error) NA_real_)
    if (is.finite(result) && result > -1e99) result else NA_real_
  }
  problem <- list(
    value = value, estimate = estimate, names = names(estimate),
    labels = profile_labels(names(estimate), fitted$distance_units),
    standard_errors = sqrt(pmax(diag(fitted$covariance), 0)),
    lower = c(rep(-30, p), -20, -20 + log_scale, rep(-20, length(estimate) - p - 2L)),
    upper = c(rep(30, p), 20, 10 + log_scale, rep(20, length(estimate) - p - 2L))
  )
  names(problem$standard_errors) <- names(estimate)
  # The tolerance of the certified fits: about 1e-4 log-likelihood units, far
  # below what a deviance needs.
  problem$control <- make_optimiser_control(abs(value(estimate)))
  best <- maximise_profile(problem, integer(), numeric(), estimate)
  problem$bank_estimate <- best$par
  problem$max_loglik <- best$value
  problem
}

##' Readable labels for the parameters of a fitted model
##' @noRd
profile_labels <- function(names, distance_units) {
  labels <- sub("^beta\\.", "", names)
  labels[names == "sigma2"] <- "Spatial variance"
  labels[names == "phi"] <- paste0("Scale (", distance_units %||% "units", ")")
  labels[names %in% c("tau2", "nu2")] <- "Nugget variance ratio"
  as.list(stats::setNames(labels, names))
}

##' Match user-supplied parameter names to those of the fit
##' @noRd
match_profile_parameters <- function(problem, parameters) {
  matched <- vapply(parameters, function(name) {
    hit <- c(name, paste0("beta.", name))
    hit <- hit[hit %in% problem$names]
    if (!length(hit)) {
      stop("Unknown parameter '", name, "'. Choose from: ",
           paste(problem$names, collapse = ", "), ".", call. = FALSE)
    }
    hit[1]
  }, "")
  unname(matched)
}

##' Maximise the log-likelihood on the profile bank with some parameters fixed
##' @noRd
maximise_profile <- function(problem, fixed, fixed_values, start) {
  free <- setdiff(seq_along(start), fixed)
  full <- function(free_values) {
    working <- start
    working[fixed] <- fixed_values
    working[free] <- free_values
    working
  }
  negative <- function(free_values) {
    result <- problem$value(full(free_values))
    if (is.finite(result)) -result else 1e100
  }
  fit <- tryCatch(stats::nlminb(
    pmin(pmax(start[free], problem$lower[free]), problem$upper[free]),
    negative, lower = problem$lower[free], upper = problem$upper[free],
    control = problem$control %||% make_optimiser_control()
  ), error = function(error) NULL)
  if (is.null(fit) || fit$objective >= 1e99) {
    return(list(par = full(start[free]), value = NA_real_, convergence = NA))
  }
  list(par = full(fit$par), value = -fit$objective,
       convergence = fit$convergence)
}

##' Marginal profile of one parameter, stepping outwards from the estimate
##'
##' Regular steps of 2 `span` standard errors over `n_points - 1` cover `span`
##' standard errors on each side. Beyond them the steps double until the
##' deviance passes the 99% point, with at most `n_points` points on a side.
##' One more point is then solved at each crossing of the `level` cut-off
##' whose neighbouring points are more than a regular step apart.
##' @noRd
profile_one_parameter <- function(problem, name, n_points, span,
                                  level = 0.95) {
  index <- match(name, problem$names)
  centre <- problem$bank_estimate[index]
  step <- 2 * span * problem$standard_errors[[name]] / (n_points - 1)
  if (!is.finite(step) || step <= 0) step <- 0.25
  regular <- (n_points - 1) %/% 2
  beyond <- stats::qchisq(0.99, 1)
  cutoff <- stats::qchisq(level, 1)
  deviance_of <- function(result) 2 * (problem$max_loglik - result$value)
  solve_at <- function(value, start) {
    result <- maximise_profile(problem, index, value, start)
    result$second_start <- FALSE
    # Constrained maxima can be multimodal, and continuation from a
    # neighbour can follow a lower branch, so a point that may decide an end
    # of the interval is also solved from the estimate.
    if (!is.finite(result$value) || deviance_of(result) > cutoff) {
      other <- maximise_profile(problem, index, value, problem$bank_estimate)
      if (is.finite(other$value) &&
          (!is.finite(result$value) || other$value > result$value + 1e-6)) {
        result <- other
        result$second_start <- TRUE
      }
    }
    result
  }
  walk <- function(direction) {
    points <- list()
    start <- problem$bank_estimate
    offset <- 0
    for (k in seq_len(n_points)) {
      offset <- offset + if (k <= regular) step else step * 2^(k - regular)
      value <- centre + direction * offset
      if (value < problem$lower[index] || value > problem$upper[index]) break
      result <- solve_at(value, start)
      points[[k]] <- profile_row(problem, index, value, result)
      if (is.finite(result$value)) start <- result$par
      if (k >= regular && is.finite(result$value) &&
          deviance_of(result) > beyond) break
    }
    points
  }
  centre_row <- profile_row(problem, index, centre,
                            list(par = problem$bank_estimate,
                                 value = problem$max_loglik, convergence = 0L,
                                 second_start = FALSE))
  profile <- do.call(rbind, c(rev(walk(-1)), list(centre_row), walk(1)))
  refine_profile_crossings(problem, profile, index, step, cutoff, solve_at)
}

##' Solve one more point at each crossing of the cut-off that lies between
##' points more than a step apart, at the crossing interpolated in the signed
##' root of the deviance
##' @noRd
refine_profile_crossings <- function(problem, profile, index, step, cutoff,
                                     solve_at) {
  profile <- profile[order(profile$value), ]
  root <- sqrt(pmax(2 * (problem$max_loglik - profile$loglik), 0))
  lowest <- which.min(ifelse(is.finite(root), root, Inf))
  added <- list()
  for (side in c(-1, 1)) {
    k <- crossing_index(root, lowest, side, sqrt(cutoff))
    if (is.na(k)) next
    inner <- k - side
    if (!is.finite(root[inner]) ||
        abs(profile$value[k] - profile$value[inner]) <= 1.01 * step) next
    value <- stats::approx(root[c(inner, k)], profile$value[c(inner, k)],
                           xout = sqrt(cutoff))$y
    start <- unlist(profile[inner, problem$names])
    added[[length(added) + 1L]] <- profile_row(problem, index, value,
                                               solve_at(value, start))
  }
  if (length(added)) profile <- rbind(profile, do.call(rbind, added))
  profile <- profile[order(profile$value), ]
  rownames(profile) <- NULL
  profile
}

##' Index of the first point beyond `threshold` on one side of `lowest`, or NA
##' @noRd
crossing_index <- function(root, lowest, side, threshold) {
  along <- if (side < 0) rev(seq_len(lowest)) else lowest:length(root)
  above <- along[is.finite(root[along]) & root[along] > threshold]
  if (!length(above)) NA_integer_ else above[1L]
}

##' One row of a profile: the fixed value, the maximised log-likelihood and the
##' other parameters at the maximum
##' @noRd
profile_row <- function(problem, index, value, result) {
  others <- stats::setNames(as.list(result$par), problem$names)
  data.frame(parameter = problem$names[index], value = value,
             loglik = result$value, convergence = result$convergence,
             second_start = isTRUE(result$second_start), others,
             check.names = FALSE)
}

##' Joint profile of two parameters on a grid, in a snake order so that each
##' point starts from its neighbour's solution. Where the marginal profile of a
##' parameter is available, its axis covers the values whose marginal deviance
##' is below the 99% point of a chi-squared distribution on two degrees of
##' freedom, which contains the joint 95% region; otherwise it spans `span`
##' standard errors either way.
##' @noRd
profile_two_parameters <- function(problem, pair, n_grid, span, marginal = list()) {
  index <- match(pair, problem$names)
  axes <- lapply(seq_along(pair), function(j) {
    profile <- marginal[[pair[j]]]
    if (!is.null(profile)) {
      deviance <- 2 * (problem$max_loglik - profile$loglik)
      inside <- profile$value[is.finite(deviance) & deviance <= stats::qchisq(0.99, 2)]
      if (length(inside) >= 2) {
        margin <- diff(range(inside)) / (2 * (n_grid - 1))
        return(seq(min(inside) - margin, max(inside) + margin, length.out = n_grid))
      }
    }
    se <- problem$standard_errors[[pair[j]]]
    if (!is.finite(se) || se <= 0) se <- 0.5
    centre <- problem$bank_estimate[index[j]]
    seq(centre - span * se, centre + span * se, length.out = n_grid)
  })
  rows <- vector("list", n_grid^2); start <- problem$bank_estimate; k <- 0
  for (i in seq_len(n_grid)) {
    order_j <- if (i %% 2 == 1) seq_len(n_grid) else rev(seq_len(n_grid))
    for (j in order_j) {
      values <- c(axes[[1]][i], axes[[2]][j])
      result <- maximise_profile(problem, index, values, start)
      if (is.finite(result$value)) start <- result$par
      k <- k + 1
      rows[[k]] <- data.frame(first = values[1], second = values[2],
                              loglik = result$value,
                              convergence = result$convergence)
    }
  }
  joint <- do.call(rbind, rows)
  attr(joint, "pair") <- pair
  joint
}

##' Profile deviance relative to the highest log-likelihood found
##' @noRd
add_profile_deviance <- function(profile, max_loglik) {
  profile$deviance <- 2 * (max_loglik - profile$loglik)
  profile
}

##' Profile and Wald intervals from the marginal profiles
##' @noRd
profile_intervals <- function(problem, marginal, level) {
  cutoff <- stats::qchisq(level, 1)
  rows <- lapply(names(marginal), function(name) {
    profile <- marginal[[name]][is.finite(marginal[[name]]$deviance), ]
    profile <- profile[order(profile$value), ]
    root <- sqrt(pmax(profile$deviance, 0))
    lowest <- which.min(root)
    # Interpolated in the signed root of the deviance, which is close to
    # linear in the parameter.
    crossing <- function(side) {
      k <- crossing_index(root, lowest, side, sqrt(cutoff))
      if (is.na(k) || !is.finite(root[k - side])) return(NA_real_)
      inner <- k - side
      stats::approx(root[c(inner, k)], profile$value[c(inner, k)],
                    xout = sqrt(cutoff))$y
    }
    estimate <- problem$estimate[[name]]
    se <- problem$standard_errors[[name]]
    z <- stats::qnorm(1 - (1 - level) / 2)
    data.frame(parameter = name, label = problem$labels[[name]],
               estimate = estimate, profile_lower = crossing(-1),
               profile_upper = crossing(1), wald_lower = estimate - z * se,
               wald_upper = estimate + z * se)
  })
  intervals <- do.call(rbind, rows)
  positive <- !grepl("^beta\\.", intervals$parameter)
  columns <- c("estimate", "profile_lower", "profile_upper", "wald_lower", "wald_upper")
  intervals[positive, columns] <- exp(intervals[positive, columns])
  rownames(intervals) <- NULL
  intervals
}

##' @export
print.RiskMap_profile <- function(x, ...) {
  cat("Profile likelihood of a RiskMap fit (", x$settings$n_draws,
      " Monte Carlo draws)\n", sep = "")
  cat(round(100 * x$level), "% intervals; variances and the scale on their ",
      "natural scale\n\n", sep = "")
  table <- x$intervals[, c("label", "estimate", "profile_lower", "profile_upper",
                           "wald_lower", "wald_upper")]
  names(table) <- c("Parameter", "Estimate", "Profile lower", "Profile upper",
                    "Wald lower", "Wald upper")
  print(table, row.names = FALSE, digits = 4)
  if (anyNA(x$intervals[, c("profile_lower", "profile_upper")])) {
    cat("\nNA: the profile does not reach the cut-off on that side within the",
        "grid or the\nbounds of the parameter, so the interval is open there.\n")
  }
  if (x$higher_than_estimate > 0.01) {
    cat("\nA point on the profile is higher than the estimate by",
        signif(x$higher_than_estimate, 3), "log-likelihood units.\n")
  }
  invisible(x)
}
