test_that("profiles of a binomial fit have their grids, deviances and intervals", {
  profiles <- profile(binomial_model, parameters = c("cov", "sigma2"),
                      pairs = list(c("sigma2", "phi")), n_points = 5,
                      n_grid = 3, n_draws = 200, messages = FALSE)
  expect_s3_class(profiles, "RiskMap_profile")
  expect_named(profiles$marginal, c("beta.cov", "sigma2"))
  expect_named(profiles$joint, "sigma2:phi")
  expect_equal(nrow(profiles$joint[[1]]), 9L)
  for (marginal in profiles$marginal) {
    expect_gte(nrow(marginal), 5L)
    expect_gte(min(marginal$deviance, na.rm = TRUE), -1e-8)
  }
  intervals <- profiles$intervals
  expect_equal(intervals$parameter, c("beta.cov", "sigma2"))
  inside <- with(intervals, (is.na(profile_lower) | profile_lower <= estimate) &
                   (is.na(profile_upper) | estimate <= profile_upper))
  expect_true(all(inside))
  # Variances are reported on their natural scale.
  expect_equal(intervals$estimate[2], exp(binomial_model$estimate$sigma2))
})

test_that("profiling restores the random number state", {
  set.seed(42)
  expected <- stats::runif(1)
  set.seed(42)
  invisible(profile(binomial_model, parameters = "cov", pairs = list(),
                    n_points = 5, n_draws = 200, messages = FALSE))
  expect_equal(stats::runif(1), expected)
})

test_that("profile rejects fits and settings it cannot handle", {
  uncertified <- binomial_model
  uncertified$mcml_history <- NULL
  expect_error(profile(uncertified), "binomial and Poisson")
  expect_error(profile(binomial_model, parameters = "kappa", messages = FALSE),
               "Unknown parameter")
  expect_error(profile(binomial_model, n_points = 2), "n_points")
  expect_error(profile(binomial_model, level = 1.5), "level")
})

test_that("the profile plot draws its panels and returns them", {
  profiles <- profile(binomial_model, parameters = c("sigma2", "phi"),
                      pairs = list(c("sigma2", "phi")), n_points = 5,
                      n_grid = 3, n_draws = 200, messages = FALSE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  panels <- plot(profiles)
  expect_named(panels, c("sigma2", "phi", "sigma2:phi"))
  expect_true(all(vapply(panels, ggplot2::is_ggplot, logical(1))))
  expect_output(print(profiles), "Profile likelihood")
})

# A profile problem with an exact log-likelihood in two parameters.
make_toy_profile <- function(loglik, se = c(1, 1)) {
  names <- c("beta.a", "b")
  problem <- list(
    value = loglik, estimate = stats::setNames(c(0, 0), names), names = names,
    labels = list(beta.a = "a", b = "b"),
    standard_errors = stats::setNames(se, names),
    lower = c(-30, -30), upper = c(30, 30)
  )
  problem$control <- make_optimiser_control(1)
  best <- maximise_profile(problem, integer(), numeric(), c(0, 0))
  problem$bank_estimate <- best$par
  problem$max_loglik <- best$value
  problem
}

test_that("profile intervals are exact for a quadratic log-likelihood", {
  problem <- make_toy_profile(function(w) -0.5 * w[1]^2 - 0.5 * (w[2] - w[1])^2)
  profile <- profile_one_parameter(problem, "beta.a", 5L, 3)
  intervals <- profile_intervals(
    problem, list(beta.a = add_profile_deviance(profile, problem$max_loglik)),
    0.95)
  expect_equal(c(intervals$profile_lower, intervals$profile_upper),
               c(-1, 1) * stats::qnorm(0.975), tolerance = 1e-5)
  expect_type(profile$second_start, "logical")
})

test_that("a flat profile is extended with doubling steps and left open", {
  problem <- make_toy_profile(function(w) -0.001 * w[1]^2 - 0.5 * w[2]^2)
  profile <- profile_one_parameter(problem, "beta.a", 5L, 3)
  # Two regular steps of 1.5, then steps of 3, 6 and 12, five points a side.
  expect_equal(sort(profile$value),
               c(-24, -12, -6, -3, -1.5, 0, 1.5, 3, 6, 12, 24), tolerance = 1e-8)
  intervals <- profile_intervals(
    problem, list(beta.a = add_profile_deviance(profile, problem$max_loglik)),
    0.95)
  expect_true(is.na(intervals$profile_lower) && is.na(intervals$profile_upper))
})

test_that("a crossing between distant points is refined", {
  # The deviance is quadratic near the estimate and much steeper beyond 2,
  # so the interval ends are at 1.96 but the grid jumps over them.
  loglik <- function(w) {
    a <- abs(w[1])
    -0.5 * a^2 - 5 * max(a - 2, 0)^3 - 0.5 * w[2]^2
  }
  problem <- make_toy_profile(loglik, se = c(0.5, 1))
  profile <- profile_one_parameter(problem, "beta.a", 5L, 1)
  intervals <- profile_intervals(
    problem, list(beta.a = add_profile_deviance(profile, problem$max_loglik)),
    0.95)
  expect_equal(intervals$profile_upper, stats::qnorm(0.975), tolerance = 1e-4)
})

test_that("profiles on two cores equal those on one", {
  skip_on_os("windows")
  skip_if_not_installed("RhpcBLASctl")
  # The workers use one BLAS thread. The serial profile does too here, since
  # with several threads its sums are ordered differently and the profiles
  # agree only up to rounding.
  blas_threads <- RhpcBLASctl::blas_get_num_procs()
  omp_threads <- RhpcBLASctl::omp_get_max_threads()
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  on.exit({
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    RhpcBLASctl::omp_set_num_threads(omp_threads)
  }, add = TRUE)
  serial <- profile(binomial_model, parameters = c("cov", "sigma2"),
                    pairs = list(c("sigma2", "phi")), n_points = 5,
                    n_grid = 3, n_draws = 200, messages = FALSE)
  parallel <- profile(binomial_model, parameters = c("cov", "sigma2"),
                      pairs = list(c("sigma2", "phi")), n_points = 5,
                      n_grid = 3, n_draws = 200, cores = 2,
                      messages = FALSE)
  expect_identical(parallel$marginal, serial$marginal)
  expect_identical(parallel$joint, serial$joint)
  expect_identical(parallel$intervals, serial$intervals)
  expect_error(profile(binomial_model, cores = 0), "cores")
})

test_that("a coefficient is profiled in the user's units on the standardised scale", {
  problem <- make_profile_problem(binomial_model, 200L, 1L)
  index <- match("beta.(Intercept)", problem$names)
  value <- unname(problem$bank_estimate[index] + 0.3)
  standardised <- maximise_profile(problem, index, value,
                                   problem$bank_estimate)
  # The intercept stays fixed in the user's units, although centring mixes it
  # with the other coefficients on the standardised scale.
  expect_identical(unname(standardised$par[[index]]), value)
  expect_equal(problem$value(standardised$par), standardised$value,
               tolerance = 1e-8)
  # Started from its own solution, the optimiser stays there.
  again <- maximise_profile(problem, index, value, standardised$par)
  expect_equal(again$value, standardised$value, tolerance = 1e-3)
})

# A problem on a standardised scale with an exact log-likelihood: internal
# coefficients a1, a2 and a covariance parameter c, user coefficients
# b1 = a1 + a2 and b2 = a2, so that fixing b1 solves for a1.
make_scaled_problem <- function(bound = 30) {
  transform <- matrix(c(1, 0, 1, 1), 2)
  loglik <- function(internal) {
    -0.5 * ((internal[1] - 1)^2 + (internal[2] - 2)^2 + (internal[3] - 3)^2)
  }
  to <- function(w) c(solve(transform, w[1:2]), w[3])
  from <- function(x) c(drop(transform %*% x[1:2]), x[3])
  list(value = function(w) loglik(to(w)), names = c("beta.b1", "beta.b2", "c"),
       lower = c(-30, -30, -20), upper = c(30, 30, 20),
       internal = list(value = loglik, to = to, from = from,
                       transform = transform, p = 2L,
                       lower = c(-bound, -bound, -20),
                       upper = c(bound, bound, 20)),
       control = make_optimiser_control(1))
}

test_that("fixed coefficients and covariance parameters can come in any order", {
  problem <- make_scaled_problem()
  start <- c(3, 2, 3)
  first <- maximise_profile(problem, c(1, 3), c(4, 5), start)
  second <- maximise_profile(problem, c(3, 1), c(5, 4), start)
  expect_equal(second$value, first$value, tolerance = 1e-8)
  expect_equal(second$par, first$par, tolerance = 1e-6)
  # The reported log-likelihood is the one at the returned parameters.
  for (result in list(first, second)) {
    expect_equal(problem$value(result$par), result$value, tolerance = 1e-8)
  }
})

test_that("a coefficient solved from the constraint respects its bounds", {
  # With internal bounds of 1, fixing b1 = a1 + a2 at 2.4 needs a1 = 1.4 or
  # more, so the constraint is infeasible.
  problem <- make_scaled_problem(bound = 1)
  result <- maximise_profile(problem, 1, 2.4, c(0, 0, 3))
  expect_true(is.na(result$value))
  expect_identical(result$status, "infeasible")
})

test_that("a feasible constraint is solved from an infeasible start", {
  # Fixing b1 = a1 + a2 at 29 with b2 = a2 = -2 at the start needs a1 = 31,
  # outside the bounds of 30, although a1 = 14, a2 = 15 is feasible and is
  # the constrained maximum.
  problem <- make_scaled_problem()
  result <- maximise_profile(problem, 1, 29, c(0, -2, 3))
  expect_identical(result$status, "ok")
  expect_equal(result$value, -169, tolerance = 1e-6)
  expect_equal(unname(result$par), c(29, 15, 3), tolerance = 1e-4)
})

test_that("the nearest feasible coefficients meet the constraint and bounds", {
  constraint <- matrix(c(1, 1), 1)
  feasible <- feasible_coefficients(constraint, 29, c(31, -2), c(-30, -30),
                                    c(30, 30))
  expect_identical(feasible$status, "ok")
  expect_equal(sum(feasible$par), 29, tolerance = 1e-10)
  expect_true(all(abs(feasible$par) <= 30))
  # The projection of (31, -2) onto a1 + a2 = 29 within the bounds.
  expect_equal(feasible$par, c(30, -1), tolerance = 1e-4)
  impossible <- feasible_coefficients(constraint, 61, c(31, 30), c(-30, -30),
                                      c(30, 30))
  expect_identical(impossible$status, "infeasible")
})

test_that("a failure to find a value is not reported as infeasible", {
  problem <- make_scaled_problem()
  problem$internal$value <- function(internal) NA_real_
  result <- maximise_profile(problem, 1, 4, c(3, 2, 3))
  expect_true(is.na(result$value))
  expect_identical(result$status, "failed")
})

# A well-identified linear geostatistical model with a nugget: 40 locations
# with four observations each, so that the likelihood uses the Woodbury
# identity. The shared Gaussian fixtures have one observation at each of ten
# locations, and some of their variances are estimated near zero, so their
# profiles reach the bounds, where brute-force maximisation is no check.
replicated_data <- local({
  set.seed(2026)
  n_locations <- 40L
  locations <- cbind(stats::runif(n_locations, 0, 20000),
                     stats::runif(n_locations, 0, 20000))
  correlation <- matern_correlation(dist(locations), phi = 4000, kappa = 0.5,
                                    return_sym_matrix = TRUE)
  field <- drop(t(chol(correlation)) %*% stats::rnorm(n_locations)) +
    stats::rnorm(n_locations, sd = sqrt(0.3))
  rows <- rep(seq_len(n_locations), each = 4L)
  data <- data.frame(x = locations[rows, 1], z = locations[rows, 2],
                     cov = stats::rnorm(length(rows)))
  data$y <- 1 + 0.5 * data$cov + field[rows] +
    stats::rnorm(length(rows), sd = sqrt(0.2))
  sf::st_as_sf(data, coords = c("x", "z"), crs = 32637)
})
replicated_fit <- glgpm(y ~ cov + gp(nugget = TRUE), data = replicated_data,
                        family = "gaussian", messages = FALSE)

test_that("the profile likelihood of a Gaussian fit is the likelihood of the fit", {
  for (fitted in list(gaussian_model, gaussian_offset_model, replicated_fit)) {
    problem <- make_gaussian_profile_problem(fitted)
    expect_equal(problem$value(problem$estimate), fitted$log_lik,
                 tolerance = 1e-10)
    expect_equal(problem$max_loglik, fitted$log_lik, tolerance = 1e-6)
  }
})

test_that("the Woodbury identity and direct factorisation give one likelihood", {
  expect_false(use_direct_gaussian_covariance(160L, 40L, NULL, TRUE))
  # The second fit has a random effect, so the latent effects have more than
  # one block.
  for (fitted in list(replicated_fit, gaussian_model)) {
    problem <- make_gaussian_profile_problem(fitted)
    scaling <- scale_spatial_coordinates(fitted$coords)
    likelihood <- function(woodbury) {
      make_gaussian_likelihood(
        y = fitted$y, D = standardize_design_matrix(fitted$D)$design,
        distances = pairwise_distances(scaling$coordinates),
        kappa = fitted$kappa, ID_coords = fitted$ID_coords,
        ID_re = if (is.null(fitted$ID_re)) NULL else as.matrix(fitted$ID_re),
        fix_tau2 = fitted$fix_tau2, fix_var_me = fitted$fix_var_me,
        woodbury = woodbury
      )
    }
    latent <- likelihood(TRUE)
    direct <- likelihood(FALSE)
    at_estimate <- problem$internal$to(pmax(problem$estimate, -10))
    p <- ncol(fitted$D)
    for (shift in c(0, 0.5, -2)) {
      internal <- at_estimate
      internal[-seq_len(p)] <- internal[-seq_len(p)] + shift
      expect_equal(latent$value(internal), direct$value(internal),
                   tolerance = 1e-10)
      expect_equal(latent$coefficients(internal),
                   direct$coefficients(internal), tolerance = 1e-8)
    }
  }
})

test_that("Gaussian profiles have every parameter, with the labels of summary()", {
  profiles <- profile(gaussian_model, pairs = list(c("sigma2", "phi")),
                      n_points = 5, n_grid = 3, n_draws = 1,
                      messages = FALSE)
  expect_named(profiles$marginal, names(unlist(gaussian_model$estimate)))
  expect_identical(nrow(profiles$joint[[1]]), 9L)
  summaries <- summary(gaussian_model)
  expect_identical(profiles$intervals$label,
                   c(rownames(summaries$reg_coef), rownames(summaries$sp),
                     rownames(summaries$me), rownames(summaries$ranef)))
  inside <- with(profiles$intervals,
                 (is.na(profile_lower) | profile_lower <= estimate) &
                   (is.na(profile_upper) | estimate <= profile_upper))
  expect_true(all(inside))
  expect_null(profiles$settings$n_draws)
  printed <- capture.output(print(profiles))
  expect_match(printed[1], "exact likelihood")
  expect_false(any(grepl("Monte Carlo", printed)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_named(plot(profiles), c(names(profiles$marginal), "sigma2:phi"))
})

test_that("a Gaussian profile point maximises the likelihood over the rest", {
  problem <- make_gaussian_profile_problem(replicated_fit)
  # The likelihood maximised over all the free parameters, the coefficients
  # included, without generalised least squares.
  brute_force <- function(index, value) {
    negative <- function(free) {
      working <- problem$estimate
      working[-index] <- free
      working[index] <- value
      -problem$value(working)
    }
    best <- stats::optim(problem$estimate[-index], negative, method = "BFGS",
                         control = list(reltol = 1e-14, maxit = 1000))
    -stats::optim(best$par, negative,
                  control = list(reltol = 1e-14, maxit = 5000))$value
  }
  for (index in seq_along(problem$names)) {
    value <- problem$estimate[[index]] + problem$standard_errors[[index]]
    point <- maximise_profile(problem, index, value, problem$bank_estimate)
    expect_identical(point$status, "ok")
    expect_identical(unname(point$par[index]), value)
    expect_equal(problem$value(point$par), point$value, tolerance = 1e-10)
    # The optimiser stops within about 1e-4 log-likelihood units.
    expect_lt(abs(point$value - brute_force(index, value)), 1e-4)
  }
})

test_that("Gaussian profiles of fixed variances match refits", {
  problem <- make_gaussian_profile_problem(replicated_fit)
  measurement <- match("sigma2_me", problem$names)
  value <- problem$estimate[[measurement]] + 0.3
  point <- maximise_profile(problem, measurement, value,
                            problem$bank_estimate)
  refit <- glgpm(y ~ cov + gp(nugget = TRUE), data = replicated_data,
                 family = "gaussian", fix_var_me = exp(value),
                 messages = FALSE)
  expect_lt(abs(point$value - refit$log_lik), 1e-4)
  # The nugget is profiled as its variance, so a fixed point is a model with
  # the nugget fixed at that variance.
  nugget <- match("tau2", problem$names)
  expect_false("nu2" %in% problem$names)
  expect_equal(problem$estimate[["tau2"]],
               replicated_fit$estimate$nu2 + replicated_fit$estimate$sigma2,
               tolerance = 1e-12)
  value <- problem$estimate[[nugget]] + 1
  point <- maximise_profile(problem, nugget, value, problem$bank_estimate)
  refit <- glgpm(y ~ cov + gp(nugget = exp(value)), data = replicated_data,
                 family = "gaussian", messages = FALSE)
  expect_lt(abs(point$value - refit$log_lik), 1e-4)
})

test_that("the Wald interval of the nugget variance is that of summary()", {
  profiles <- profile(replicated_fit, parameters = "tau2", n_points = 5,
                      messages = FALSE)
  expect_identical(profiles$labels$tau2, "Variance of the nugget")
  nugget <- summary(replicated_fit)$sp["Variance of the nugget", ]
  expect_equal(unlist(profiles$intervals[c("estimate", "wald_lower",
                                           "wald_upper")]),
               nugget, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("a joint Gaussian profile can fix every covariance parameter", {
  fitted <- glgpm(y ~ gp(), data = replicated_data, family = "gaussian",
                  fix_var_me = 0.2, messages = FALSE)
  profiles <- profile(fitted, parameters = "phi",
                      pairs = list(c("sigma2", "phi")), n_points = 5,
                      n_grid = 3, messages = FALSE)
  joint <- profiles$joint[[1]]
  expect_true(all(is.finite(joint$loglik)))
  expect_true(all(joint$deviance >= -1e-8))
})

test_that("Gaussian profiles on two cores equal those on one", {
  skip_on_os("windows")
  serial <- profile(replicated_fit, parameters = c("cov", "sigma2"),
                    pairs = list(c("sigma2", "phi")), n_points = 5,
                    n_grid = 3, messages = FALSE)
  parallel <- profile(replicated_fit, parameters = c("cov", "sigma2"),
                      pairs = list(c("sigma2", "phi")), n_points = 5,
                      n_grid = 3, cores = 2, messages = FALSE)
  expect_identical(parallel$marginal, serial$marginal)
  expect_identical(parallel$joint, serial$joint)
  expect_identical(parallel$intervals, serial$intervals)
})

test_that("a poorly identified variance has an open profile interval", {
  skip_on_cran()
  # The Galicia lead data cannot tell the measurement error variance from
  # zero, so the profile stays below the cut-off down to the lower bound,
  # while the upper end is refined to the crossing of the cut-off.
  data(galicia)
  fitted <- glgpm(log(lead) ~ gp(kappa = 1.5), data = galicia,
                  family = "gaussian", distance_units = "km",
                  messages = FALSE)
  profiles <- profile(fitted, parameters = "sigma2_me", messages = FALSE)
  interval <- profiles$intervals
  expect_true(is.na(interval$profile_lower))
  expect_equal(interval$profile_upper, 0.0763, tolerance = 2e-3)
  marginal <- profiles$marginal$sigma2_me
  below <- marginal$value < log(interval$estimate)
  expect_lt(max(marginal$deviance[below]), stats::qchisq(0.95, 1))
  expect_lt(min(marginal$value), -15)
  cutoff <- stats::qchisq(0.95, 1)
  expect_lt(min(abs(marginal$deviance - cutoff)), 0.01)
})
