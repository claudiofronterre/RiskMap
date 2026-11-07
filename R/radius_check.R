##' Check of the radius of the certificate, run on request by `check_fit()`
##'
##' Predicts, at a fitted binomial or Poisson model, how often the radius of
##' its Monte Carlo certificate is too small. Short Markov chains on the
##' conditional distribution of the latent effects reach the rare draws of
##' large importance weight that a validation bank usually misses. They are
##' pooled with one bank of the size the fit used, and banks resampled from
##' the pool show how the certificate behaves.
##'
##' The certificate is the length of a Newton correction plus a radius that
##' allows for the Monte Carlo error of the certification bank, in statistical
##' standard errors, and the radius is estimated from the bank's own draws. For
##' the radius the fit uses (`certify_with_replicate()`), the check predicts
##' the level as run, the percentage of samples whose error, in the form the
##' radius allows for, exceeds it, which the 95% of the certificate refers to.
##' It gives its interval from resampling the banks, which leaves out the
##' error of the chains, the percentage of certification banks above their own
##' radius, the percentage of fits whose Monte Carlo error exceeds the
##' certificate, and the factor by which the radius would have to grow for
##' the level to be 5%. Outside binomial with the logit link the factor is
##' raised by a `margin` of 15%, since the predictions were low on simulated
##' Poisson surveys. The row is `"SUSPICIOUS"` when the raised factor exceeds
##' one, which is when the level, with that margin, is above 5%.
##'
##' The chains use elliptical slice sampling (Murray, Adams and MacKay, 2010)
##' with the Gaussian part of the fit's importance distribution as their
##' prior, and are pooled with the bank by multiple importance sampling with
##' the balance heuristic, with the normalising constant from the bridge
##' equation of Meng and Wong (1996). On simulated binomial and Poisson
##' surveys, one bank with 16 chains of 4,000 rounds predicted the level within
##' a factor of 1.65 and the factor within 14%. The chains can miss a separate
##' mode of the latent effects, so a low prediction is not a proof. The check
##' takes about as long as the fit. Without a `seed` its streams are 9500 to
##' 9502 after the seed of the fit, so it is reproducible, and the caller's
##' random number state is left unchanged.
##'
##' Returns the `row` of the check (`fit_check_row()`), the `panel` data of
##' each resampled bank's error and radius, and the `details`, with the
##' `level` and its `level_interval`, `own`, `exceedance`, `factor`, `margin`,
##' `widening` (the raised factor, at least one), the `certificate` and the
##' certificate with its radius so widened (`widened`), the `tolerance`, the
##' median Pareto shape of 100 resampled banks (`pareto_k`) and the bank's own
##' (`bank_pareto_k`), the split R-hat of the chains' log weight (`rhat`), the
##' `settings` and the elapsed `seconds`.
##' @noRd
fit_check_radius <- function(object, chains = 16L, rounds = 4000L,
                             banks = 1000L, seed = NULL) {
  if (!inherits(object, "RiskMap")) {
    stop("'object' must be a 'RiskMap' object returned by 'glgpm'",
         call. = FALSE)
  }
  check_positive_integer(chains, "chains")
  check_positive_integer(rounds, "rounds")
  check_positive_integer(banks, "banks")
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)
  if (rounds < 10) stop("'rounds' must be at least 10", call. = FALSE)
  if (banks < 3) stop("'banks' must be at least 3", call. = FALSE)
  history <- object$mcml_history
  if (!object$family %in% c("binomial", "poisson") || !length(history)) {
    stop("the radius can be checked only for binomial and Poisson models ",
         "fitted by certified Monte Carlo maximum likelihood", call. = FALSE)
  }
  last <- history[[length(history)]]
  if (!is.finite(last$certificate_radius %||% NA_real_)) {
    stop("the fit has no certificate whose radius could be checked ",
         "(outcome '", last$outcome %||% "unknown", "')", call. = FALSE)
  }
  started <- proc.time()[["elapsed"]]
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  streams <- make_radius_check_streams(seed, object$seed)
  problem <- rebuild_certification_problem(object, last)
  bank <- draw_radius_bank(problem, streams[[1L]])
  chain <- run_radius_chains(problem, chains, rounds, streams[[2L]])
  pool <- pool_radius_draws(bank, chain)
  resampled <- resample_radius_banks(pool, problem, banks, streams[[3L]])
  rates <- measure_radius_rule(resampled$errors, resampled$radii)
  margin <- if (is_binomial_logit(object)) 0 else 0.15
  widening <- max(1, rates$factor * (1 + margin))
  details <- list(
    level = 100 * rates$level, level_interval = 100 * rates$level_interval,
    own = 100 * rates$own, exceedance = 100 * rates$exceedance,
    factor = rates$factor, margin = margin, widening = widening,
    certificate = last$certificate,
    widened = last$certificate_distance + widening * last$certificate_radius,
    tolerance = attr(object, "control_mcml")$tolerance %||% NA_real_,
    pareto_k = stats::median(resampled$pareto_k), bank_pareto_k = bank$pareto_k,
    rhat = chain$rhat,
    settings = list(chains = as.integer(chains), rounds = as.integer(rounds),
                    banks = as.integer(banks), draws = problem$draws,
                    seed = streams[[1L]]),
    seconds = proc.time()[["elapsed"]] - started
  )
  list(row = describe_radius_check(details),
       panel = data.frame(error = sqrt(rowSums(resampled$errors^2)),
                          radius = resampled$radii),
       details = details)
}

##' Whether a fit is binomial with the logit link
##' @noRd
is_binomial_logit <- function(object) {
  identical(object$family, "binomial") &&
    (identical(object$link_function$name, "canonical") ||
       identical(object$link_function$link, "logit"))
}

##' The row of the check of the radius, from its details
##' @noRd
describe_radius_check <- function(details) {
  percent <- function(value) paste0(signif(value, 2), "%")
  found <- paste0(
    "The radius of the bound is predicted to be exceeded in ",
    percent(details$level), " of samples (", percent(details$level_interval[1]),
    " to ", percent(details$level_interval[2]), "), against 5%.",
    if (details$margin > 0) {
      paste0(" Outside binomial models with the logit link the prediction ",
             "is read with a margin of ", 100 * details$margin, "%, since it ",
             "ran low on simulated Poisson surveys.")
    }
  )
  if (details$widening <= 1) {
    return(fit_check_row("Radius of the bound", "OK",
                         paste0(found, " The radius holds.")))
  }
  fit_check_row(
    "Radius of the bound", "SUSPICIOUS",
    paste0(found, " The radius would have to widen by a factor of ",
           signif(details$widening, 3), ", which makes the bound ",
           signif(details$widened, 2), " standard errors (limit ",
           details$tolerance, ")."),
    paste0("Read the bound with its radius so widened. More draws, through ",
           "'max_samples' in set_control_mcml(), bring more of the tail of ",
           "the weights into each sample.")
  )
}

##' Seeds of the three random steps of `fit_check_radius()`
##'
##' Without a `seed` they are streams 9500 to 9502 after the seed of the fit,
##' which `fit_adaptive_mcml()` does not use, wrapped as it wraps its streams.
##' With one they are that seed and the two after it.
##' @noRd
make_radius_check_streams <- function(seed, fit_seed) {
  base <- as.numeric(if (is.null(seed)) fit_seed else seed)
  offsets <- if (is.null(seed)) 9500 + 0:2 else 0:2
  as.integer(((base - 1 + offsets) %% .Machine$integer.max) + 1)
}

##' The certification problem of a fit, rebuilt at its estimate
##'
##' Rebuilds the engine as `make_profile_problem()` does, with the link and
##' the proposal of the fit's last stage, anchors the global proposal at the
##' estimate, and recovers the information over the free parameters of that
##' stage from the covariance of the fit, on the engine's scale. Gives the
##' estimate on that scale, the upper Cholesky factor `root` of the
##' information, the `draws` of the last validation sample with its
##' `defensive` rows, and the derivatives of the log weights
##' (`make_radius_derivatives()`).
##' @noRd
rebuild_certification_problem <- function(object, last) {
  design_scaling <- standardize_design_matrix(object$D)
  spatial_scaling <- scale_spatial_coordinates(object$coords)
  distances <- as.matrix(pairwise_distances(spatial_scaling$coordinates))
  control_mcml <- attr(object, "control_mcml") %||% set_control_mcml()
  engine <- make_adaptive_mcml_engine(
    D = design_scaling$design, y = object$y, units_m = object$units_m,
    distances = distances, family = object$family, kappa = object$kappa,
    cov_offset = object$cov_offset, ID_coords = object$ID_coords,
    fix_tau2 = object$fix_tau2,
    ID_re = if (is.null(object$ID_re)) NULL else as.matrix(object$ID_re),
    invlink = refit_invlink(object),
    defensive_weight = control_mcml$defensive_weight %||% 0.05,
    proposal = control_mcml$proposal %||% "global"
  )
  if (identical(last$proposal, "laplace")) engine$use_proposal("laplace")
  coefficients <- seq_len(ncol(object$D))
  maps <- profile_scale_maps(design_scaling$coefficient_transform,
                             log(spatial_scaling$distance_scale),
                             length(coefficients))
  estimate <- maps$to(unlist(object$estimate))
  engine$set_anchor(estimate)
  # The coefficients are the only parameters whose scale is not a shift.
  to_engine <- diag(length(estimate))
  to_engine[coefficients, coefficients] <-
    solve(design_scaling$coefficient_transform)
  covariance <- object$covariance
  covariance[!is.finite(covariance)] <- 0
  covariance <- to_engine %*% covariance %*% t(to_engine)
  free <- last$free %||% seq_along(estimate)
  root <- chol(solve(covariance[free, free, drop = FALSE]))
  draws <- as.integer(last$validation_budget %||%
                        attr(object, "optimizer")$mcml_validation_budget)
  problem <- list(engine = engine, estimate = estimate, free = free,
                  root = root, draws = draws,
                  defensive = defensive_rows(draws, engine$defensive_weight))
  problem$derivatives <- make_radius_derivatives(problem)
  problem
}

##' Log weights and their derivatives at draws, as a validation bank has them
##'
##' `at()` gives the log weights of draws under the mixture with the
##' defensive component and their derivatives over the free parameters, by
##' central differences with the steps of `score_validation_bank()`, in the
##' metric of the information (rows of R^-T g). `core()` gives the log
##' weights against the Gaussian component alone, log pi~ - log q_core, the
##' log-likelihood of the chains.
##' @noRd
make_radius_derivatives <- function(problem) {
  engine <- problem$engine
  estimate <- problem$estimate
  free <- problem$free
  step <- .Machine$double.eps^(1 / 3) * pmax(1, abs(estimate[free]))
  shifted <- list(estimate)
  for (parameter in seq_along(free)) {
    upper <- lower <- estimate
    upper[free[parameter]] <- upper[free[parameter]] + step[parameter]
    lower[free[parameter]] <- lower[free[parameter]] - step[parameter]
    shifted <- c(shifted, list(upper, lower))
  }
  proposals <- lapply(shifted, function(values) {
    tryCatch(engine$proposal_at(values), error = function(error) NULL)
  })
  if (any(vapply(proposals, is.null, logical(1)))) {
    stop("the importance distribution of the fit could not be rebuilt at ",
         "its estimate", call. = FALSE)
  }
  share <- mean(problem$defensive)
  weights_at <- function(variant, draws, defensive) {
    engine$evaluate(shifted[[variant]], draws, log_weights_only = TRUE,
                    built_proposal = proposals[[variant]],
                    defensive = defensive, defensive_share = share)
  }
  list(
    at = function(draws, defensive) {
      log_weights <- weights_at(1L, draws, defensive)
      gradients <- matrix(0, nrow(draws), length(free))
      for (parameter in seq_along(free)) {
        gradients[, parameter] <- (
          weights_at(2L * parameter, draws, defensive) -
            weights_at(2L * parameter + 1L, draws, defensive)
        ) / (2 * step[parameter])
      }
      valid <- is.finite(log_weights)
      gradients[!valid, ] <- 0
      if (any(!is.finite(gradients))) {
        stop("the derivatives of the importance weights are not finite at ",
             "some draws", call. = FALSE)
      }
      list(log_weights = log_weights,
           gradients = t(forwardsolve(t(problem$root), t(gradients))))
    },
    core = function(draws) {
      engine$evaluate(estimate, draws, log_weights_only = TRUE,
                      built_proposal = proposals[[1L]],
                      defensive = rep(FALSE, nrow(draws)),
                      defensive_share = 0)
    }
  )
}

##' Apply `evaluate` to the rows of `draws` in blocks
##'
##' Blocks have the rows `score_validation_bank()` uses, so that a large
##' set of draws is never evaluated at once. `evaluate` returns a list of
##' row-wise results, vectors or matrices, which are stacked.
##' @noRd
evaluate_in_blocks <- function(draws, evaluate) {
  block_rows <- max(250L, floor(2.5e6 / ncol(draws)))
  starts <- seq.int(1L, nrow(draws), by = block_rows)
  parts <- lapply(starts, function(first) {
    evaluate(first:min(nrow(draws), first + block_rows - 1L))
  })
  stats::setNames(lapply(names(parts[[1L]]), function(name) {
    pieces <- lapply(parts, `[[`, name)
    if (is.matrix(pieces[[1L]])) do.call(rbind, pieces) else unlist(pieces)
  }), names(parts[[1L]]))
}

##' One validation bank at the estimate, with each draw kept
##'
##' Drawn as `score_validation_bank()` draws it, in the same blocks from
##' `seed`, keeping each draw's log weight and derivatives and, for its
##' Gaussian rows, the log weight against the Gaussian component.
##' @noRd
draw_radius_bank <- function(problem, seed) {
  dimension <- problem$engine$latent_dimension
  block_rows <- max(250L, floor(2.5e6 / dimension))
  set.seed(seed)
  bank <- lapply(seq.int(1L, problem$draws, by = block_rows), function(first) {
    rows <- first:min(problem$draws, first + block_rows - 1L)
    block <- matrix(stats::rnorm(length(rows) * dimension),
                    nrow = length(rows))
    values <- problem$derivatives$at(block, problem$defensive[rows])
    gaussian <- !problem$defensive[rows]
    values$log_weights_core <- rep(NA_real_, length(rows))
    values$log_weights_core[gaussian] <- problem$derivatives$core(
      block[gaussian, , drop = FALSE]
    )
    values
  })
  log_weights <- unlist(lapply(bank, `[[`, "log_weights"))
  list(log_weights = log_weights,
       log_weights_core = unlist(lapply(bank, `[[`, "log_weights_core")),
       gradients = do.call(rbind, lapply(bank, `[[`, "gradients")),
       defensive = problem$defensive,
       pareto_k = summarise_importance_weights(log_weights)$pareto_k)
}

##' Elliptical slice chains on the conditional distribution of the latent
##' effects
##'
##' In the coordinates of the Gaussian component of the proposal, where it is
##' N(0, I), the target is that prior times exp(core(z)) (Murray, Adams and
##' MacKay, 2010). All chains move at once. The first fifth of the rounds are
##' discarded and every eighth of the rest is kept. Returns the kept draws'
##' log weights and derivatives as a bank has them, their log weights against
##' the Gaussian component, and the split R-hat of the log-likelihood.
##' @noRd
run_radius_chains <- function(problem, chains, rounds, seed) {
  core <- problem$derivatives$core
  dimension <- problem$engine$latent_dimension
  discarded <- rounds %/% 5L
  thinning <- 8L
  kept <- (rounds - discarded) %/% thinning
  set.seed(seed)
  state <- matrix(stats::rnorm(chains * dimension), chains)
  log_likelihood <- core(state)
  for (attempt in seq_len(100L)) {
    restart <- which(!is.finite(log_likelihood))
    if (!length(restart)) break
    state[restart, ] <- stats::rnorm(length(restart) * dimension)
    log_likelihood[restart] <- core(state[restart, , drop = FALSE])
  }
  if (any(!is.finite(log_likelihood))) {
    stop("the chains could not be started at draws of positive weight",
         call. = FALSE)
  }
  stored <- array(0, c(kept, chains, dimension))
  traced <- matrix(NA_real_, rounds, chains)
  for (round in seq_len(rounds)) {
    moved <- move_elliptical_slice(state, log_likelihood, core)
    state <- moved$state
    log_likelihood <- moved$log_likelihood
    traced[round, ] <- log_likelihood
    after <- round - discarded
    if (after > 0L && after %% thinning == 0L) {
      stored[after %/% thinning, , ] <- state
    }
  }
  # Chain by chain within each kept round.
  draws <- matrix(aperm(stored, c(2L, 1L, 3L)), ncol = dimension)
  values <- evaluate_in_blocks(draws, function(rows) {
    block <- draws[rows, , drop = FALSE]
    c(problem$derivatives$at(block, rep(FALSE, length(rows))),
      list(log_weights_core = core(block)))
  })
  values$rhat <- compute_split_rhat(traced[-seq_len(discarded), ,
                                           drop = FALSE])
  values
}

##' One round of elliptical slice sampling for every chain
##'
##' Each row of `state` is a chain with prior N(0, I) and log-likelihood
##' `log_likelihood`. A chain that finds no point above its slice in 200
##' shrinkages stays where it is.
##' @noRd
move_elliptical_slice <- function(state, log_likelihood, core) {
  chains <- nrow(state)
  direction <- matrix(stats::rnorm(length(state)), chains)
  level <- log_likelihood + log(stats::runif(chains))
  angle <- stats::runif(chains, 0, 2 * pi)
  lower <- angle - 2 * pi
  upper <- angle
  active <- seq_len(chains)
  for (shrinkage in seq_len(200L)) {
    candidate <- state[active, , drop = FALSE] * cos(angle[active]) +
      direction[active, , drop = FALSE] * sin(angle[active])
    value <- core(candidate)
    accepted <- is.finite(value) & value > level[active]
    if (any(accepted)) {
      taken <- active[accepted]
      state[taken, ] <- candidate[accepted, , drop = FALSE]
      log_likelihood[taken] <- value[accepted]
    }
    active <- active[!accepted]
    if (!length(active)) break
    negative <- angle[active] < 0
    lower[active[negative]] <- angle[active[negative]]
    upper[active[!negative]] <- angle[active[!negative]]
    angle[active] <- stats::runif(length(active), lower[active],
                                  upper[active])
  }
  list(state = state, log_likelihood = log_likelihood)
}

##' Split R-hat of a quantity traced over rounds (rows) and chains (columns)
##' @noRd
compute_split_rhat <- function(traced) {
  half <- nrow(traced) %/% 2L
  split <- cbind(traced[seq_len(half), , drop = FALSE],
                 traced[half + seq_len(half), , drop = FALSE])
  within <- mean(apply(split, 2L, stats::var))
  between <- half * stats::var(colMeans(split))
  sqrt(((half - 1) / half * within + between / half) / within)
}

##' Log normalising constant from the bridge equation
##'
##' With `log_ratios` the log of pi~ / q at `bank_draws` draws from q and
##' `chain_draws` draws from pi, the optimal bridge of Meng and Wong (1996)
##' is the L for which the sum over all draws of 1 / (n_q + n_pi pi~ / (L q))
##' is one.
##' @noRd
solve_bridge_equation <- function(log_ratios, bank_draws, chain_draws) {
  finite <- log_ratios[is.finite(log_ratios)]
  equation <- function(log_constant) {
    log(sum(1 / (bank_draws + chain_draws * exp(log_ratios - log_constant))))
  }
  stats::uniroot(equation, c(stats::median(finite) - 200, max(finite) + 60),
                 tol = 1e-12, extendInt = "upX")$root
}

##' The law of one draw of a bank, estimated from the bank and the chains
##'
##' The Gaussian rows of the bank and the chains' draws are pooled by
##' multiple importance sampling with the balance heuristic, each with
##' probability proportional to 1 / (n_q + n_pi r), r = pi~ / (L q_core), as
##' an estimate of the law of a draw from the Gaussian component. The
##' defensive rows of the bank are kept apart, since the chains add little to
##' that component.
##' @noRd
pool_radius_draws <- function(bank, chain) {
  gaussian <- which(!bank$defensive)
  log_ratios <- c(bank$log_weights_core[gaussian], chain$log_weights_core)
  bank_draws <- length(gaussian)
  chain_draws <- length(chain$log_weights_core)
  log_constant <- solve_bridge_equation(log_ratios, bank_draws, chain_draws)
  balance <- 1 / (bank_draws + chain_draws * exp(log_ratios - log_constant))
  list(log_weights = c(bank$log_weights[gaussian], chain$log_weights),
       gradients = rbind(bank$gradients[gaussian, , drop = FALSE],
                         chain$gradients),
       probability = balance / sum(balance),
       defensive_log_weights = bank$log_weights[bank$defensive],
       defensive_gradients = bank$gradients[bank$defensive, , drop = FALSE])
}

##' Score and Monte Carlo covariance of a bank from its draws
##'
##' As `score_validation_bank()` computes them, from each draw's log weight
##' and derivatives.
##' @noRd
compute_weighted_score <- function(log_weights, gradients) {
  valid <- is.finite(log_weights)
  relative <- exp(log_weights - max(log_weights[valid]))
  relative[!valid] <- 0
  weights <- relative / sum(relative)
  score <- colSums(gradients * weights)
  centred <- sweep(gradients, 2L, score) * weights
  list(score = score, covariance = crossprod(centred))
}

##' Banks of the fit's size drawn from the pool
##'
##' Each has the bank's number of defensive rows, drawn from its defensive
##' rows, and the rest from the pooled Gaussian draws. Returns their errors,
##' centred on their mean, and their radii, both in the metric of the
##' information, and the Pareto shapes of the first `shapes` of them.
##' @noRd
resample_radius_banks <- function(pool, problem, banks, seed, shapes = 100L) {
  set.seed(seed)
  defensive_draws <- sum(problem$defensive)
  gaussian_draws <- problem$draws - defensive_draws
  available <- length(pool$defensive_log_weights)
  pareto_k <- numeric()
  resampled <- lapply(seq_len(banks), function(bank) {
    rows <- sample.int(length(pool$log_weights), gaussian_draws, TRUE,
                       prob = pool$probability)
    extra <- if (defensive_draws > 0L && available > 0L) {
      sample.int(available, defensive_draws, TRUE)
    } else {
      integer()
    }
    log_weights <- c(pool$log_weights[rows], pool$defensive_log_weights[extra])
    if (bank <= shapes) {
      pareto_k[bank] <<- summarise_importance_weights(log_weights)$pareto_k
    }
    compute_weighted_score(
      log_weights,
      rbind(pool$gradients[rows, , drop = FALSE],
            pool$defensive_gradients[extra, , drop = FALSE])
    )
  })
  scores <- do.call(rbind, lapply(resampled, `[[`, "score"))
  list(errors = sweep(scores, 2L, colMeans(scores)),
       radii = compute_validation_radii(lapply(resampled, `[[`,
                                               "covariance")),
       pareto_k = pareto_k)
}

##' Radii of several covariances already in the metric of the information
##'
##' The same as `compute_validation_radius()` with the identity as the root
##' of the information, with its fixed draws made once for all covariances.
##' @noRd
compute_validation_radii <- function(covariances, probability = 0.95,
                                     n_draws = 100000L) {
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  set.seed(424242L, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  squares <- matrix(stats::rnorm(n_draws * nrow(covariances[[1L]])),
                    n_draws)^2
  vapply(covariances, function(covariance) {
    weights <- pmax(eigen((covariance + t(covariance)) / 2, symmetric = TRUE,
                          only.values = TRUE)$values, 0)
    sqrt(stats::quantile(as.numeric(squares %*% weights), probability,
                         names = FALSE))
  }, numeric(1))
}

##' Rates of the radius of the certificate on resampled banks
##'
##' On `tuples` triples (i, j, l) of distinct banks, with i the correction
##' bank, j the certification bank and l its replicate, the radius is
##' r = max(q_j, q_l) (`certify_with_replicate()`). The estimate is off by the
##' correction bank's error G_i, the measured distance is ||G_j - G_i||, and
##' the error that the radius allows for is ||(G_i + G_j) / 2||, which bounds
##' ||G_i|| by the certificate when it is below r. Gives the rates at which
##' that error exceeds r (`level`), ||G_j|| exceeds r (`own`) and ||G_i||
##' exceeds the certificate (`exceedance`), and the 95% point of the ratio of
##' that error to r, the factor that brings `level` to 5%. The 95% interval of
##' `level` comes from `replicates` resamplings of the banks with replacement,
##' on the same triples of positions.
##' @noRd
measure_radius_rule <- function(errors, radii, tuples = 40000L,
                                replicates = 200L) {
  triples <- draw_distinct_triples(nrow(errors), tuples)
  rates_of <- function(banks) {
    correction <- errors[banks[triples[, 1L]], , drop = FALSE]
    certification <- errors[banks[triples[, 2L]], , drop = FALSE]
    radius <- pmax(radii[banks[triples[, 2L]]], radii[banks[triples[, 3L]]])
    allowed <- sqrt(rowSums(((correction + certification) / 2)^2))
    list(allowed = allowed, radius = radius, correction = correction,
         certification = certification)
  }
  all <- rates_of(seq_len(nrow(errors)))
  distance <- sqrt(rowSums((all$certification - all$correction)^2))
  levels <- vapply(seq_len(replicates), function(replicate) {
    resampled <- rates_of(sample.int(nrow(errors), replace = TRUE))
    mean(resampled$allowed > resampled$radius)
  }, numeric(1))
  list(level = mean(all$allowed > all$radius),
       level_interval = stats::quantile(levels, c(0.025, 0.975),
                                        names = FALSE),
       own = mean(sqrt(rowSums(all$certification^2)) > all$radius),
       exceedance = mean(sqrt(rowSums(all$correction^2)) >
                           distance + all$radius),
       factor = stats::quantile(all$allowed / all$radius, 0.95,
                                names = FALSE))
}

##' Ordered triples of distinct integers from 1 to `size`, uniformly
##'
##' The second is drawn from the other `size - 1` values and the third from
##' the remaining `size - 2`, by shifting past the values already taken.
##' @noRd
draw_distinct_triples <- function(size, tuples) {
  first <- sample.int(size, tuples, replace = TRUE)
  second <- sample.int(size - 1L, tuples, replace = TRUE)
  second <- second + (second >= first)
  third <- sample.int(size - 2L, tuples, replace = TRUE)
  low <- pmin(first, second)
  high <- pmax(first, second)
  third <- third + (third >= low)
  third <- third + (third >= high)
  cbind(first, second, third)
}
