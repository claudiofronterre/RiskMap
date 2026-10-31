##' Profile likelihood of a fitted geostatistical model
##'
##' Computes profile likelihoods of a model fitted by [glgpm()]: for each
##' chosen parameter, the log-likelihood maximised over all the other
##' parameters on a grid of values, and for each chosen pair, the same over a
##' grid of pairs. Profiles are available for linear geostatistical
##' (Gaussian) models and for binomial and Poisson models fitted with
##' certified Monte Carlo maximum likelihood. For a binomial or Poisson model
##' the Monte Carlo likelihood is evaluated on one fixed bank of `n_draws`
##' standard normal draws with the importance distribution anchored at the
##' estimate, so the profiles are smooth functions of the parameters. Use
##' [plot.RiskMap_profile()] to draw them.
##'
##' The likelihood of a Gaussian model is exact, so there is no bank of draws
##' and `n_draws` and `seed` are not used. For each value of the covariance
##' parameters the regression coefficients that maximise the likelihood are
##' found by generalised least squares, so each point is a maximisation over
##' the covariance parameters only. The parameters are those of
##' `summary()`, with its labels: the regression coefficients, the spatial
##' variance and scale, the variance of the nugget (`tau2`) when it is
##' estimated, the measurement error variance (`sigma2_me`) unless it is
##' fixed by `fix_var_me`, and the variances of the random effects. Since
##' points of an exact likelihood are cheap, each end of a profile interval
##' is refined until the deviance of a solved point is within 0.01 of the
##' cut-off. Grid steps of the variances and the scale are at most 1 on the
##' logarithmic scale, so that a variance estimated near zero, whose standard
##' error is huge, is still profiled, and each point of a joint profile is
##' solved from its neighbour's solution and from the estimate. A profile
##' that stays below the cut-off down to the lower bound of a variance (about
##' \eqn{e^{-20}}) shows that the data cannot tell that variance from zero,
##' and its interval is open there.
##'
##' The smoothness \eqn{\kappa} of the Matern correlation stays fixed at the
##' value used in the fit. It can be profiled by refitting the model over a
##' set of values of \eqn{\kappa} and comparing the maximised
##' log-likelihoods (`log_lik`) of the fits.
##'
##' Grids are centred at the estimate on the scale on which the model is fitted
##' (the logarithm for variances and the scale) and span `span` standard errors
##' either way. If the profile deviance has not passed the 99% point of a
##' chi-squared distribution on one degree of freedom by the end of a marginal
##' grid, the grid is extended with steps that double each time, up to
##' `n_points` points on each side. Each point starts the optimisation from
##' the extrapolation of the two previous solutions, which follows a ridge in
##' the other parameters. The optimisation works on the scale of the fit, with
##' standardised covariates, while a profiled coefficient stays fixed in its
##' own units. Profiles can have more than one constrained
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
##' Each point is a maximisation of the likelihood over the other
##' parameters, so the time grows with the number of points and of parameters.
##' A binomial or Poisson model that takes a minute to fit can take tens of
##' seconds per point.
##' Joint profiles cost `n_grid^2` points each and are computed only when
##' `pairs` asks for them. With `messages = TRUE` and one core the time left
##' is estimated after the first parameter. `cores` runs independent parts in
##' parallel.
##'
##' @param fitted A `RiskMap` object from [glgpm()] for a Gaussian model, or
##'   for a binomial or Poisson model fitted with certified Monte Carlo
##'   maximum likelihood.
##' @param parameters Names of the parameters to profile one at a time, as in
##'   `names(unlist(fitted$estimate))` (the prefix `beta.` can be left out).
##'   For a Gaussian model the nugget is `tau2`, its variance, in place of
##'   `nu2`, its ratio to the spatial variance. The default is every
##'   parameter.
##' @param pairs A list of pairs of parameter names for joint profiles, for
##'   example `list(c("sigma2", "phi"))`. The default is none.
##' @param n_points Number of points of each marginal grid before any
##'   extension, and the most points on each side of the estimate.
##' @param n_grid Number of points along each side of a joint grid.
##' @param span Half-width of the grids, in standard errors.
##' @param n_draws Number of Monte Carlo draws in the bank. Not used for a
##'   Gaussian model.
##' @param seed Seed of the bank. The random number state of the session is
##'   restored afterwards. Not used for a Gaussian model.
##' @param level Confidence level of the reported profile intervals.
##' @param cores Number of processes. With more than one, the two sides of
##'   each marginal profile, the ends of the intervals and the joint profiles
##'   run in parallel by forking, which is not available on Windows. When
##'   the RhpcBLASctl package is installed, each process uses one BLAS and
##'   OpenMP thread, and the results are the same as on one core with one
##'   BLAS thread. With several BLAS threads on one core they agree up to
##'   rounding.
##' @param messages Logical, whether to report progress.
##' @param ... Not used.
##' @return An object of class `RiskMap_profile`: a list with the marginal
##'   profiles (`marginal`), the joint profiles (`joint`), the intervals
##'   (`intervals`), the estimate and the maximum of the log-likelihood on the
##'   bank, or of the exact likelihood (`estimate`, `max_loglik`), and the
##'   settings.
##' @examples
##' \donttest{
##' data(galicia)
##' fit <- glgpm(log(lead) ~ gp(kappa = 1.5), data = galicia,
##'              family = "gaussian", distance_units = "km", messages = FALSE)
##' profiles <- profile(fit, pairs = list(c("sigma2", "phi")))
##' profiles
##' plot(profiles)
##'
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
                            cores = 1L, messages = TRUE, ...) {
  stopifnot(inherits(fitted, "RiskMap"))
  gaussian <- identical(fitted$family, "gaussian")
  if (!gaussian && (!fitted$family %in% c("binomial", "poisson") ||
                    is.null(fitted$mcml_history))) {
    stop("Profiles are available for Gaussian models, and for binomial and ",
         "Poisson models fitted by certified Monte Carlo maximum likelihood.",
         call. = FALSE)
  }
  check_profile_settings(n_points, n_grid, span, if (!gaussian) n_draws,
                         level)
  stopifnot("'cores' must be a positive whole number" =
              is.numeric(cores) && length(cores) == 1L && cores >= 1 &&
              cores == round(cores))
  cores <- as.integer(cores)
  if (cores > 1L && .Platform$OS.type == "windows") {
    if (messages) message("Parallel profiles need forking, which Windows ",
                          "does not support, so they run on one core.")
    cores <- 1L
  }
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  if (gaussian) {
    problem <- make_gaussian_profile_problem(fitted)
  } else {
    if (messages) {
      message("Re-maximising the likelihood on a bank of ", n_draws,
              " draws...")
    }
    problem <- make_profile_problem(fitted, n_draws, seed)
  }
  parameters <- match_profile_parameters(problem, parameters %||% problem$names)
  pairs <- lapply(pairs %||% list(), function(pair) {
    match_profile_parameters(problem, pair)
  })

  if (cores > 1L) {
    marginal <- profile_in_parallel(problem, parameters, n_points, span, level,
                                    cores, messages)
  } else {
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
  }
  if (messages && length(pairs)) {
    message("Profiling ", length(pairs), " pair", if (length(pairs) > 1L) "s",
            " of parameters on a grid of ", n_grid, " by ", n_grid, "...")
  }
  joint <- run_profile_tasks(pairs, function(pair) {
    profile_two_parameters(problem, pair, n_grid, span, marginal)
  }, cores)
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
    settings = if (gaussian) list(span = span) else
      list(n_draws = n_draws, seed = seed, span = span)
  ), class = "RiskMap_profile")
}

##' Marginal profiles with the walks and then the ends run in parallel
##' @noRd
profile_in_parallel <- function(problem, parameters, n_points, span, level,
                                cores, messages) {
  plans <- lapply(parameters, function(name) {
    plan_one_parameter(problem, name, n_points, span, level)
  })
  tasks <- expand.grid(parameter = seq_along(parameters), direction = c(-1, 1))
  if (messages) {
    message("Profiling ", length(parameters), " parameter",
            if (length(parameters) > 1L) "s", " on ",
            min(cores, nrow(tasks)), " cores...")
  }
  walks <- run_profile_tasks(seq_len(nrow(tasks)), function(row) {
    plans[[tasks$parameter[row]]]$walk(tasks$direction[row])
  }, cores)
  marginal <- run_profile_tasks(seq_along(parameters), function(k) {
    side <- function(direction) {
      walks[[which(tasks$parameter == k & tasks$direction == direction)]]
    }
    plans[[k]]$finish(side(-1), side(1))
  }, cores)
  names(marginal) <- parameters
  marginal
}

##' A duration in seconds as readable text
##' @noRd
format_profile_time <- function(seconds) {
  if (seconds < 90) return(paste(max(1, round(seconds)), "seconds"))
  if (seconds < 5400) return(paste(round(seconds / 60), "minutes"))
  paste(signif(seconds / 3600, 2), "hours")
}

##' Check the settings of a profile, without `n_draws` when it is `NULL`
##' @noRd
check_profile_settings <- function(n_points, n_grid, span, n_draws, level) {
  stopifnot(
    "'n_points' must be at least 5" = is.numeric(n_points) && n_points >= 5,
    "'n_grid' must be at least 3" = is.numeric(n_grid) && n_grid >= 3,
    "'span' must be positive" = is.numeric(span) && span > 0,
    "'n_draws' must be at least 100" = is.null(n_draws) ||
      (is.numeric(n_draws) && n_draws >= 100),
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
  maps <- profile_scale_maps(transform, log_scale, p)
  to_internal <- maps$to
  from_internal <- maps$from
  estimate <- unlist(fitted$estimate)
  internal <- to_internal(estimate)
  engine$set_anchor(internal)
  set.seed(seed)
  bank <- matrix(stats::rnorm(n_draws * engine$latent_dimension), nrow = n_draws)
  value_internal <- function(internal) {
    result <- tryCatch(objective(internal, bank),
                       error = function(error) NA_real_)
    if (is.finite(result) && result > -1e99) result else NA_real_
  }
  value <- function(working) value_internal(to_internal(working))
  others <- length(estimate) - p - 2L
  problem <- list(
    value = value, estimate = estimate, names = names(estimate),
    labels = profile_labels(names(estimate), fitted$distance_units),
    standard_errors = sqrt(pmax(diag(fitted$covariance), 0)),
    lower = c(rep(-30, p), -20, -20 + log_scale, rep(-20, others)),
    upper = c(rep(30, p), 20, 10 + log_scale, rep(20, others)),
    # The optimiser works on the engine's scale (standardised covariates,
    # rescaled distances), as the fit does, with the bounds of the fit.
    internal = list(value = value_internal, to = to_internal,
                    from = from_internal, transform = transform, p = p,
                    lower = c(rep(-30, p), -20, -20, rep(-20, others)),
                    upper = c(rep(30, p), 20, 10, rep(20, others))),
    predict_starts = TRUE
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

##' Maps between the user's working scale and the scale of the fit
##'
##' The coefficients are standardised with `transform`, and the logarithm of
##' the scale is shifted by `log_scale`, the logarithm of the largest
##' distance. The other parameters are the same on both scales.
##' @noRd
profile_scale_maps <- function(transform, log_scale, p) {
  list(
    to = function(working) {
      c(solve(transform, working[seq_len(p)]), working[p + 1L],
        working[p + 2L] - log_scale, working[-seq_len(p + 2L)])
    },
    from = function(internal) {
      c(drop(transform %*% internal[seq_len(p)]), internal[p + 1L],
        internal[p + 2L] + log_scale, internal[-seq_len(p + 2L)])
    }
  )
}

##' Readable labels for the parameters of a fitted model
##'
##' Those of a Gaussian fit are the labels of `summary()`.
##' @noRd
profile_labels <- function(names, distance_units, family = NULL) {
  labels <- sub("^beta\\.", "", names)
  if (identical(family, "gaussian")) {
    labels[names == "sigma2"] <- "Spatial process var."
    labels[names == "phi"] <- paste0("Spatial corr. scale (",
                                     distance_units %||% "units", ")")
    labels[names == "tau2"] <- "Variance of the nugget"
    labels[names == "sigma2_me"] <- "Measurement error var."
    random <- grepl("^sigma2_re\\.", names)
    labels[random] <- paste0(sub("^sigma2_re\\.", "", names[random]),
                             " (random eff. var.)")
    return(as.list(stats::setNames(labels, names)))
  }
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
##'
##' The fixed parameters and the start are on the user's scale, and so is the
##' result, but the optimiser works on the engine's scale, with standardised
##' covariates and rescaled distances, as the fit does. In the user's units a
##' coefficient of a covariate with a large range is badly scaled for the
##' optimiser, which then stops well below the maximum. A fixed variance,
##' scale or nugget is fixed at its value on that scale. A fixed coefficient
##' is a linear constraint on the standardised coefficients, met exactly by
##' solving for one of them (chosen by a pivoted QR decomposition) from the
##' others. This matters for the intercept, which centring mixes with the
##' other coefficients. The fixed values are matched to the fixed parameters
##' by position, and a solved coefficient outside its bounds makes the point
##' infeasible. Since the optimiser gets no direction from an infeasible
##' point, the coefficients of the start are first moved to the nearest point
##' that meets both the constraint and the bounds (`feasible_coefficients()`).
##' A problem without `internal` (the tests' stand-ins) is optimised on the
##' user's scale. When `internal` has a `concentrate` function, as for
##' Gaussian fits, the coefficients are not searched: for each value of the
##' other parameters they are set to the solution of `concentrate()`, which
##' meets the constraint of the fixed coefficients exactly.
##'
##' A stop with a non-zero code is solved again from the same start with a
##' central-difference gradient with a step of 1e-4, keeping the better
##' solution. `budget` caps the number of evaluations, and the best point
##' evaluated so far is returned when the cap is reached. The number of
##' evaluations is returned, and `status` is `"ok"` with a finite value,
##' `"infeasible"` when no coefficients within the bounds meet the constraint,
##' and `"failed"` when a feasible point exists but no finite value was found.
##' @noRd
maximise_profile <- function(problem, fixed, fixed_values, start,
                             budget = Inf) {
  scale <- problem$internal %||% list(
    value = problem$value, to = identity, from = identity, transform = NULL,
    p = 0L, lower = problem$lower, upper = problem$upper
  )
  p <- scale$p
  concentrate <- scale$concentrate
  target <- start
  target[fixed] <- fixed_values
  internal_start <- scale$to(target)
  fixed_beta <- fixed[fixed <= p]
  fixed_other <- fixed[fixed > p]
  pivot <- integer()
  constraint <- beta_targets <- NULL
  if (length(fixed_beta)) {
    # The targets of the fixed coefficients, by position, since the fixed
    # parameters may come in any order.
    beta_targets <- fixed_values[match(fixed_beta, fixed)]
    constraint <- scale$transform[fixed_beta, , drop = FALSE]
    if (is.null(concentrate)) {
      pivot <- qr(constraint)$pivot[seq_along(fixed_beta)]
      rest <- setdiff(seq_len(p), pivot)
    }
  }
  # Coefficients that `concentrate` solves for the other parameters are not
  # searched.
  solved <- if (is.null(concentrate)) pivot else seq_len(p)
  free <- setdiff(seq_along(start), c(solved, fixed_other))
  lower <- scale$lower[free]
  upper <- scale$upper[free]
  feasibility <- "ok"
  if (length(pivot)) {
    coefficients <- seq_len(p)
    feasible <- feasible_coefficients(constraint, beta_targets,
                                      internal_start[coefficients],
                                      scale$lower[coefficients],
                                      scale$upper[coefficients])
    feasibility <- feasible$status
    if (identical(feasibility, "ok")) {
      internal_start[coefficients] <- feasible$par
    }
  }
  full_internal <- function(free_values) {
    internal <- internal_start
    internal[free] <- free_values
    if (!is.null(concentrate)) {
      internal[seq_len(p)] <- concentrate(internal, constraint, beta_targets)
    } else if (length(pivot)) {
      internal[pivot] <- solve(
        constraint[, pivot, drop = FALSE],
        beta_targets - constraint[, rest, drop = FALSE] %*% internal[rest]
      )
    }
    internal
  }
  full <- function(free_values) {
    working <- scale$from(full_internal(free_values))
    working[fixed] <- fixed_values
    working
  }
  evaluations <- 0L
  best <- list(objective = Inf, par = NULL)
  negative <- function(free_values) {
    if (evaluations >= budget) {
      stop(structure(class = c("profile_budget", "error", "condition"),
                     list(message = "budget reached", call = NULL)))
    }
    evaluations <<- evaluations + 1L
    internal <- full_internal(free_values)
    # A coefficient solved from the constraint must respect the bounds too,
    # otherwise the point is outside the parameter space.
    outside <- length(pivot) &&
      any(internal[pivot] < scale$lower[pivot] |
            internal[pivot] > scale$upper[pivot])
    result <- if (outside) NA_real_ else scale$value(internal)
    objective <- if (is.finite(result)) -result else 1e100
    if (objective < best$objective) {
      best <<- list(objective = objective, par = free_values)
    }
    objective
  }
  gradient <- function(free_values) {
    vapply(seq_along(free_values), function(j) {
      up <- down <- free_values
      up[j] <- min(free_values[j] + 1e-4, upper[j])
      down[j] <- max(free_values[j] - 1e-4, lower[j])
      (negative(up) - negative(down)) / (up[j] - down[j])
    }, numeric(1))
  }
  run <- function(from, with_gradient) {
    tryCatch(
      stats::nlminb(pmin(pmax(from, lower), upper), negative,
                    if (with_gradient) gradient,
                    lower = lower, upper = upper,
                    control = problem$control %||% make_optimiser_control()),
      profile_budget = function(condition) NULL,
      error = function(error) NULL
    )
  }
  better <- function(first, second) {
    if (is.null(second)) return(first)
    if (is.null(first) || second$objective < first$objective) second else first
  }
  from <- internal_start[free]
  if (!identical(feasibility, "ok")) {
    return(list(par = full(from), value = NA_real_, convergence = NA,
                evaluations = 0L,
                status = if (identical(feasibility, "infeasible"))
                  "infeasible" else "failed"))
  }
  # Nothing left to search, as when a joint profile fixes every covariance
  # parameter of a Gaussian fit.
  if (!length(free)) {
    objective <- negative(from)
    return(list(par = full(from),
                value = if (objective < 1e99) -objective else NA_real_,
                convergence = 0L, evaluations = evaluations,
                status = if (objective < 1e99) "ok" else "failed"))
  }
  fit <- run(from, FALSE)
  if (!is.null(fit) && fit$convergence != 0L && fit$objective < 1e99) {
    fit <- better(fit, run(from, TRUE))
  }
  if (is.null(fit)) {
    if (is.null(best$par) || best$objective >= 1e99) {
      return(list(par = full(from), value = NA_real_,
                  convergence = NA, evaluations = evaluations,
                  status = "failed"))
    }
    return(list(par = full(best$par), value = -best$objective,
                convergence = 1L, evaluations = evaluations, status = "ok"))
  }
  if (fit$objective >= 1e99) {
    return(list(par = full(from), value = NA_real_, convergence = NA,
                evaluations = evaluations, status = "failed"))
  }
  list(par = full(fit$par), value = -fit$objective,
       convergence = fit$convergence, evaluations = evaluations,
       status = "ok")
}

##' Coefficients that meet a linear constraint within bounds, near a start
##'
##' Finds the point of `{a : constraint %*% a = targets}` within slightly
##' narrowed bounds that is nearest to `start`, by Dykstra's alternating
##' projections between the affine set and the box, which converge to that
##' point whenever the two sets meet. The bounds are narrowed by a millionth
##' of their width, so that a coefficient solved again from the others stays
##' within the original bounds despite rounding. When the projections do not
##' reach a feasible point, the smallest violation of the constraint within
##' the bounds, a convex problem, tells an infeasible constraint
##' (`"infeasible"`) from one where no point was found (`"not found"`).
##' @noRd
feasible_coefficients <- function(constraint, targets, start, lower, upper) {
  inside <- function(a, low, high) all(a >= low & a <= high)
  if (inside(start, lower, upper)) return(list(par = start, status = "ok"))
  margin <- ifelse(is.finite(upper - lower), 1e-6 * (upper - lower), 0)
  low <- lower + margin
  high <- upper - margin
  normal <- tcrossprod(constraint)
  to_affine <- function(a) {
    drop(a - crossprod(constraint,
                       solve(normal, constraint %*% a - targets)))
  }
  x <- to_affine(start)
  box_correction <- affine_correction <- 0
  for (iteration in seq_len(5000L)) {
    y <- pmin(pmax(x + box_correction, low), high)
    box_correction <- x + box_correction - y
    updated <- to_affine(y + affine_correction)
    affine_correction <- y + affine_correction - updated
    change <- max(abs(updated - x))
    x <- updated
    if (change < 1e-12) break
  }
  if (inside(x, lower, upper)) return(list(par = x, status = "ok"))
  violation <- function(a) sum((constraint %*% a - targets)^2)
  slope <- function(a) {
    drop(2 * crossprod(constraint, constraint %*% a - targets))
  }
  closest <- stats::nlminb(pmin(pmax(start, lower), upper), violation, slope,
                           lower = lower, upper = upper)
  tolerance <- (1e-8 * (1 + sum(abs(targets))))^2
  list(par = start,
       status = if (closest$objective > tolerance) "infeasible" else
         "not found")
}

##' Marginal profile of one parameter, stepping outwards from the estimate
##'
##' Regular steps of 2 `span` standard errors over `n_points - 1` cover `span`
##' standard errors on each side. Beyond them the steps double until the
##' deviance passes the 99% point, with at most `n_points` points on a side.
##' Each point starts from its neighbour's solution or from the extrapolation
##' of the last two solutions (`problem$predict_starts`). Constrained maxima
##' can be multimodal, and continuation
##' can follow a lower branch, so the two points around each crossing of the
##' `level` cut-off are also solved from the estimate, with a capped number
##' of evaluations, keeping the higher solution. One more point is then solved
##' at each crossing whose neighbouring points are more than a regular step
##' apart.
##' @noRd
profile_one_parameter <- function(problem, name, n_points, span,
                                  level = 0.95) {
  plan <- plan_one_parameter(problem, name, n_points, span, level)
  plan$finish(plan$walk(-1), plan$walk(1))
}

##' The steps of a marginal profile, as functions that can run apart
##'
##' `walk(direction)` steps out on one side of the estimate, and
##' `finish(below, above)` joins the two walks, re-solves the points around
##' each crossing from the estimate and refines the crossings
##' (`profile_one_parameter()`). The two walks are independent, so they can run
##' in parallel.
##' @noRd
plan_one_parameter <- function(problem, name, n_points, span, level = 0.95) {
  index <- match(name, problem$names)
  centre <- problem$bank_estimate[index]
  step <- 2 * span * problem$standard_errors[[name]] / (n_points - 1)
  if (!is.finite(step) || step <= 0) step <- 0.25
  step <- min(step, profile_step_cap(problem, name))
  regular <- (n_points - 1) %/% 2
  beyond <- stats::qchisq(0.99, 1)
  cutoff <- stats::qchisq(level, 1)
  deviance_of <- function(result) 2 * (problem$max_loglik - result$value)
  # The start of each point: its neighbour's solution, or, with
  # `problem$predict_starts`, the linear extrapolation of the last two
  # solutions in the profiled parameter, which follows a ridge.
  start_for <- function(solutions, value) {
    last <- solutions[[length(solutions)]]
    if (!isTRUE(problem$predict_starts) || length(solutions) < 2L) {
      return(last$par)
    }
    before <- solutions[[length(solutions) - 1L]]
    ratio <- (value - last$value) / (last$value - before$value)
    predicted <- last$par + ratio * (last$par - before$par)
    pmin(pmax(predicted, problem$lower), problem$upper)
  }
  walk <- function(direction) {
    points <- list()
    solutions <- list(list(value = centre, par = problem$bank_estimate))
    offset <- 0
    for (k in seq_len(n_points)) {
      offset <- offset + if (k <= regular) step else step * 2^(k - regular)
      value <- centre + direction * offset
      if (value < problem$lower[index] || value > problem$upper[index]) break
      result <- maximise_profile(problem, index, value,
                                 start_for(solutions, value))
      points[[k]] <- profile_row(problem, index, value, result)
      # Beyond a value that no coefficients within the bounds can reach, the
      # profile has no more points on this side.
      if (identical(result$status, "infeasible")) break
      if (!is.finite(result$value)) next
      solutions <- c(solutions, list(list(value = value, par = result$par)))
      if (k >= regular && deviance_of(result) > beyond) break
    }
    points
  }
  second_start <- function(value, result) {
    other <- maximise_profile(problem, index, value, problem$bank_estimate,
                              budget = 400L)
    if (is.finite(other$value) &&
        (!is.finite(result$value) || other$value > result$value + 1e-6)) {
      other$second_start <- TRUE
      return(other)
    }
    result
  }
  finish <- function(below, above) {
    centre_row <- profile_row(problem, index, centre,
                              list(par = problem$bank_estimate,
                                   value = problem$max_loglik,
                                   convergence = 0L))
    profile <- do.call(rbind, c(rev(below), list(centre_row), above))
    profile <- recheck_profile_crossings(problem, profile, index, cutoff,
                                         second_start)
    refine_profile_crossings(problem, profile, index, step, cutoff,
                             function(value, start) {
                               second_start(value, maximise_profile(
                                 problem, index, value, start))
                             })
  }
  list(walk = walk, finish = finish)
}

##' The largest grid step of a parameter
##'
##' A variance estimated near zero, or a likelihood that is nearly flat, has
##' a huge standard error on the log scale, and steps in proportion to it
##' would leave the bounds at once. A problem can cap the steps of the
##' parameters on the log scale with `max_step` (Gaussian fits use 1). The
##' coefficients, in their own units, are not capped.
##' @noRd
profile_step_cap <- function(problem, name) {
  if (grepl("^beta\\.", name)) Inf else problem$max_step %||% Inf
}

##' Apply a function to tasks, in parallel by forking when `cores` is above 1
##'
##' Each worker limits BLAS and OpenMP to one thread when RhpcBLASctl is
##' installed, because several forked processes that each use every core
##' oversubscribe the machine. The setting dies with the worker. On Windows,
##' where R cannot fork, the tasks run one after another. An error in a
##' worker stops the profile with its message.
##' @noRd
run_profile_tasks <- function(tasks, work, cores) {
  if (cores <= 1L || length(tasks) <= 1L) return(lapply(tasks, work))
  results <- parallel::mclapply(tasks, function(task) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    }
    work(task)
  }, mc.cores = min(cores, length(tasks)), mc.preschedule = FALSE)
  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("A profile task failed: ",
         conditionMessage(attr(results[[which(failed)[1L]]], "condition")),
         call. = FALSE)
  }
  results
}

##' Solve the two points around each crossing of the cut-off again from the
##' estimate, keeping the higher solution, until the points around each
##' crossing have both been checked
##' @noRd
recheck_profile_crossings <- function(problem, profile, index, cutoff,
                                      second_start) {
  checked <- numeric()
  repeat {
    profile <- profile[order(profile$value), ]
    root <- sqrt(pmax(2 * (problem$max_loglik - profile$loglik), 0))
    lowest <- which.min(ifelse(is.finite(root), root, Inf))
    pending <- integer()
    for (side in c(-1, 1)) {
      k <- crossing_index(root, lowest, side, sqrt(cutoff))
      if (is.na(k)) next
      pending <- c(pending, k, k - side)
    }
    pending <- pending[!profile$value[pending] %in% checked &
                         pending != lowest]
    if (!length(pending)) break
    for (row in unique(pending)) {
      current <- list(value = profile$loglik[row],
                      par = unlist(profile[row, problem$names]),
                      convergence = profile$convergence[row],
                      status = profile$status[row],
                      second_start = profile$second_start[row])
      result <- second_start(profile$value[row], current)
      profile[row, ] <- profile_row(problem, index, profile$value[row], result)
      checked <- c(checked, profile$value[row])
    }
  }
  rownames(profile) <- NULL
  profile
}

##' Solve one more point at each crossing of the cut-off that lies between
##' points more than a step apart, at the crossing interpolated in the signed
##' root of the deviance
##'
##' An exact likelihood (`problem$exact`) is cheap to solve again, so its
##' crossings are refined in up to eight rounds, wherever the points are,
##' until the deviance of a point next to each crossing is within 0.01 of the
##' cut-off. The signed root of the deviance can be far from linear over a
##' step when the quadratic approximation is poor.
##' @noRd
refine_profile_crossings <- function(problem, profile, index, step, cutoff,
                                     solve_at) {
  exact <- isTRUE(problem$exact)
  for (round in seq_len(if (exact) 8L else 1L)) {
    profile <- profile[order(profile$value), ]
    root <- sqrt(pmax(2 * (problem$max_loglik - profile$loglik), 0))
    lowest <- which.min(ifelse(is.finite(root), root, Inf))
    added <- list()
    for (side in c(-1, 1)) {
      k <- crossing_index(root, lowest, side, sqrt(cutoff))
      if (is.na(k)) next
      inner <- k - side
      if (!is.finite(root[inner])) next
      if (exact) {
        if (min(abs(root[c(inner, k)]^2 - cutoff)) < 0.01) next
      } else if (abs(profile$value[k] - profile$value[inner]) <= 1.01 * step) {
        next
      }
      value <- stats::approx(root[c(inner, k)], profile$value[c(inner, k)],
                             xout = sqrt(cutoff))$y
      start <- unlist(profile[inner, problem$names])
      added[[length(added) + 1L]] <- profile_row(problem, index, value,
                                                 solve_at(value, start))
    }
    if (!length(added)) break
    profile <- rbind(profile, do.call(rbind, added))
  }
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
             status = result$status %||% "ok",
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
    half <- min(span * se, profile_step_cap(problem, pair[j]) * (n_grid - 1) / 2)
    centre <- problem$bank_estimate[index[j]]
    seq(centre - half, centre + half, length.out = n_grid)
  })
  rows <- vector("list", n_grid^2); start <- problem$bank_estimate; k <- 0
  for (i in seq_len(n_grid)) {
    order_j <- if (i %% 2 == 1) seq_len(n_grid) else rev(seq_len(n_grid))
    for (j in order_j) {
      values <- c(axes[[1]][i], axes[[2]][j])
      result <- maximise_profile(problem, index, values, start)
      if (isTRUE(problem$exact)) {
        # A neighbour's solution with a variance near zero lies on a plateau
        # where the optimiser stops, so points of a cheap exact likelihood
        # are also solved from the estimate.
        other <- maximise_profile(problem, index, values,
                                  problem$bank_estimate)
        if (is.finite(other$value) &&
            (!is.finite(result$value) || other$value > result$value)) {
          result <- other
        }
      }
      if (is.finite(result$value)) start <- result$par
      k <- k + 1
      rows[[k]] <- data.frame(first = values[1], second = values[2],
                              loglik = result$value,
                              convergence = result$convergence,
                              status = result$status)
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
  # A profile without draws is of a Gaussian fit, whose likelihood is exact.
  cat("Profile likelihood of a RiskMap fit (",
      if (is.null(x$settings$n_draws)) "exact likelihood" else
        paste(x$settings$n_draws, "Monte Carlo draws"), ")\n", sep = "")
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
