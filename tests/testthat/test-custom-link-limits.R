## Custom inverse links at the limits of their range and where the likelihood
## is not log-concave. Each of these failed a certified fit before the likelihood
## of a custom link was computed by bound_log_probability() and
## compute_custom_link_terms().

make_link <- function(inv, d1, d2) {
  make_inverse_link_functions("binomial", list(inv = inv, d1 = d1, d2 = d2), 3L)
}
cloglog_link <- make_link(
  function(eta) -expm1(-exp(eta)),
  function(eta) exp(eta - exp(eta)),
  function(eta) exp(eta - exp(eta)) * (1 - exp(eta))
)
probit_link <- make_link(pnorm, dnorm, function(eta) -eta * dnorm(eta))
imperfect_link <- make_link(
  function(eta) 0.01 + 0.89 * plogis(eta),
  function(eta) 0.89 * dlogis(eta),
  function(eta) 0.89 * dlogis(eta) * (1 - 2 * plogis(eta))
)
softplus_link <- make_inverse_link_functions(
  "poisson",
  list(inv = function(eta) log1p(exp(eta)), d1 = plogis, d2 = dlogis), 3L
)

test_that("a draw at a limit of the link leaves the other draws unchanged", {
  y <- c(3, 20, 0)
  units_m <- rep(20, 3)
  moderate <- c(0, 1, -1)
  rounded <- c(5, 1, -1)
  expect_identical(cloglog_link$inv(5), 1)
  log_likelihood <- compute_observation_log_likelihood(
    rbind(moderate, rounded, moderate), y, units_m, "binomial", cloglog_link
  )
  exact <- sum(y * log(cloglog_link$inv(moderate)) -
                 (units_m - y) * exp(moderate))
  expect_equal(log_likelihood[c(1L, 3L)], rep(exact, 2L), ignore_attr = TRUE,
               tolerance = 1e-12)
  # The failures of the first site are charged at the bound 2^-53.
  others <- 20 * log(cloglog_link$inv(1)) - 20 * exp(-1)
  expect_equal(log_likelihood[[2L]], 17 * (-53 * log(2)) + others,
               tolerance = 1e-12)
})

test_that("a mean at a limit costs nothing if no observation contradicts it", {
  expect_identical(
    compute_observation_log_likelihood(matrix(50), 20, 20, "binomial",
                                       cloglog_link),
    0
  )
  expect_identical(pnorm(-50), 0)
  expect_identical(
    compute_observation_log_likelihood(matrix(-50), 0, 20, "binomial",
                                       probit_link),
    0
  )
  expect_identical(softplus_link$inv(-800), 0)
  expect_identical(
    compute_observation_log_likelihood(matrix(-800), 0, 1, "poisson",
                                       softplus_link),
    0
  )
})

test_that("only the draws outside the support have zero likelihood", {
  # Built by hand, because a link like this is refused when it is supplied.
  link <- list(
    inv = function(eta) {
      ifelse(eta > 10, NaN, ifelse(eta < -10, -0.1, plogis(eta)))
    },
    d1 = dlogis, d2 = function(eta) dlogis(eta) * (1 - 2 * plogis(eta)),
    name = "custom"
  )
  log_likelihood <- compute_observation_log_likelihood(
    rbind(0, 11, -11), 3, 20, "binomial", link
  )
  expect_true(is.finite(log_likelihood[[1L]]))
  expect_identical(log_likelihood[2:3], c(-Inf, -Inf))
  expect_null(laplace_observation_terms(11, 3, 20, "binomial", link))
  expect_null(laplace_observation_terms(-11, 3, 20, "binomial", link))
})

test_that("a link that goes wrong in its tails is refused when it is supplied", {
  expect_error(
    make_inverse_link_functions(
      "binomial", function(eta) ifelse(abs(eta) > 6, NaN, pnorm(eta)), 3L
    ),
    "not a number or is outside the range of the family"
  )
  expect_error(
    make_inverse_link_functions(
      "binomial", function(eta) pnorm(eta) + (eta > 20), 3L
    ),
    "linear predictor of 20.5"
  )
  expect_error(
    make_inverse_link_functions(
      "poisson", function(eta) ifelse(eta < -10, -1, exp(eta)), 3L
    ),
    "not a number or is outside the range of the family"
  )
  # A valid link at its limits is accepted, and so is overflow for a mean.
  expect_no_error(make_inverse_link_functions("binomial", pnorm, 3L))
  expect_no_error(
    make_inverse_link_functions("poisson", function(eta) exp(exp(eta)), 3L)
  )
})

test_that("a custom link's terms are the derivatives of its log-likelihood", {
  eta <- c(-2.5, -0.7, 0.4, 1.1)
  y <- c(0, 3, 11, 20)
  units_m <- rep(20, 4)
  for (linkf in list(cloglog_link, probit_link, imperfect_link)) {
    log_likelihood <- function(point) {
      sum(y * log(linkf$inv(point)) + (units_m - y) * log1p(-linkf$inv(point)))
    }
    terms <- laplace_observation_terms(eta, y, units_m, "binomial", linkf)
    expect_equal(terms$log_likelihood, log_likelihood(eta), tolerance = 1e-12)
    expect_equal(terms$score, numDeriv::grad(log_likelihood, eta),
                 tolerance = 1e-8)
    expect_equal(terms$curvature,
                 -diag(numDeriv::hessian(log_likelihood, eta)),
                 tolerance = 1e-6)
  }
  counts <- c(0, 2, 5, 9)
  exposure <- rep(2, 4)
  log_likelihood <- function(point) {
    mean <- softplus_link$inv(point)
    sum(counts * log(exposure * mean) - exposure * mean)
  }
  terms <- laplace_observation_terms(eta, counts, exposure, "poisson",
                                     softplus_link)
  expect_equal(terms$log_likelihood, log_likelihood(eta), tolerance = 1e-12)
  expect_equal(terms$score, numDeriv::grad(log_likelihood, eta),
               tolerance = 1e-8)
  expect_equal(terms$curvature, -diag(numDeriv::hessian(log_likelihood, eta)),
               tolerance = 1e-6)
})

test_that("a custom link's terms agree with the canonical formula for logit", {
  eta <- c(-30, -3, 0, 3, 30)
  y <- c(0, 2, 9, 18, 20)
  units_m <- rep(20, 5)
  custom <- laplace_observation_terms(
    eta, y, units_m, "binomial",
    make_link(plogis, dlogis, function(point) {
      dlogis(point) * (1 - 2 * plogis(point))
    })
  )
  canonical <- laplace_observation_terms(
    eta, y, units_m, "binomial", make_inverse_link_functions("binomial", NULL, 5L)
  )
  expect_equal(custom$score, canonical$score, tolerance = 1e-10)
  expect_equal(custom$curvature, canonical$curvature, tolerance = 1e-10)
  expect_equal(custom$log_likelihood, canonical$log_likelihood,
               tolerance = 1e-10)
})

test_that("the curvature survives far into the tails of the probability", {
  # (p (1 - p))^2 underflows below a probit linear predictor of -27.
  eta <- c(-30, -37)
  expect_identical(pnorm(eta)^2, c(0, 0))
  terms <- laplace_observation_terms(eta, c(3, 3), c(20, 20), "binomial",
                                     probit_link)
  mills <- exp(dnorm(eta, log = TRUE) - pnorm(eta, log.p = TRUE))
  expect_equal(terms$score, 3 * mills, tolerance = 1e-10)
  expect_equal(terms$curvature, 3 * mills * (eta + mills), tolerance = 1e-8)
  # Without a success the likelihood is flat there, and nothing overflows.
  flat <- laplace_observation_terms(c(-40, 40), c(0, 20), c(20, 20),
                                    "binomial", probit_link)
  expect_identical(flat$score, c(0, 0))
  expect_identical(flat$curvature, c(0, 0))
})

test_that("a bank that reaches the limit estimates the likelihood", {
  y <- 18
  trials <- 20
  beta <- 2
  sigma2 <- 2
  parameters <- c(beta, log(sigma2), log(1))
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, nrow = 1), y = y, units_m = trials, distances = matrix(0),
    family = "binomial",
    invlink = cloglog_link[c("inv", "d1", "d2")]
  )
  exact <- log(integrate(
    function(effect) {
      rate <- exp(beta + effect)
      exp(y * log(-expm1(-rate)) - (trials - y) * rate) *
        dnorm(effect, sd = sqrt(sigma2))
    },
    lower = -Inf, upper = Inf, rel.tol = 1e-12
  )$value)
  set.seed(8128)
  bank <- matrix(rnorm(20000), ncol = 1)
  expect_true(any(cloglog_link$inv(beta + sqrt(sigma2) * bank) == 1))
  expect_equal(engine$evaluate(parameters, bank), exact, tolerance = 0.005)
})

test_that("a saturated bank still has a validation score", {
  set.seed(2)
  sites <- 30L
  distances <- as.matrix(dist(cbind(runif(sites), runif(sites))))
  y <- rbinom(sites, 20, 0.8)
  engine <- make_adaptive_mcml_engine(
    D = matrix(1, sites, 1L), y = y, units_m = rep(20, sites),
    distances = distances, family = "binomial",
    invlink = cloglog_link[c("inv", "d1", "d2")]
  )
  parameters <- c(0.5, log(3), log(0.2))
  set.seed(3)
  bank <- matrix(rnorm(200 * sites), 200L)
  log_weights <- engine$evaluate(parameters, bank, log_weights_only = TRUE)
  expect_true(all(is.finite(log_weights)))
  expect_true(is.finite(engine$evaluate(parameters, bank)))
  validation <- score_validation_bank(engine, parameters, 1:3, 500L, 4L)
  expect_false(is.null(validation))
  expect_true(all(is.finite(validation$score)))
})

test_that("the mode is found where the likelihood is convex", {
  sites <- 4L
  y <- c(8, 0, 0, 0)
  units_m <- rep(20, sites)
  mu <- rep(-6, sites)
  precision <- diag(0.25, sites)
  latent_map <- make_latent_map(cbind(seq_len(sites)), sites)
  start_terms <- laplace_observation_terms(mu, y, units_m, "binomial",
                                           imperfect_link)
  # The Newton step of the search does not exist at the start.
  expect_lt(start_terms$curvature[[1L]], 0)
  expect_error(chol(precision + diag(start_terms$curvature)))
  found <- adaptive_laplace_mode(mu, y, units_m, latent_map, precision,
                                 "binomial", imperfect_link, numeric(sites))
  terms <- laplace_observation_terms(mu + found$mode, y, units_m, "binomial",
                                     imperfect_link)
  expect_false(found$reached_limit)
  expect_equal(terms$score, as.numeric(precision %*% found$mode),
               tolerance = 1e-7)
  expect_gt(found$mode[[1L]], 2)
  expect_no_error(chol(found$hessian))
})

test_that("the global proposal starts from a positive definite precision", {
  sites <- 3L
  y <- c(8, 8, 0)
  units_m <- rep(20, sites)
  mu <- rep(-6, sites)
  precision <- diag(0.25, sites)
  latent_map <- make_latent_map(cbind(seq_len(sites)), sites)
  # A point where the curvature is too negative for the factors at the mode to
  # give a positive definite precision.
  curvature_at <- laplace_observation_terms(mu, y, units_m, "binomial",
                                            imperfect_link)$curvature
  expect_lt(min(curvature_at + diag(precision)), 0)
  proposal <- fit_global_proposal(
    numeric(sites), mu, y, units_m, latent_map, precision, "binomial",
    imperfect_link, canonical = FALSE
  )
  expect_false(is.null(proposal))
  expect_no_error(chol(proposal$precision))
})

test_that("a quadratic fit without enough weighted nodes is a failed refit", {
  sites <- 3L
  precision <- diag(0.25, sites)
  latent_map <- make_latent_map(cbind(seq_len(sites)), sites)
  local_mocked_bindings(update_observation_factors = function(...) {
    list(linear = c(0, NA, 0), curvature = c(1, NA, 1))
  })
  expect_null(refit_global_factors(
    c(rep(0, sites), rep(1, sites)), rep(-1, sites), c(1, 2, 3),
    rep(20, sites), latent_map, precision, "binomial", probit_link, FALSE,
    build_normal_quadrature(), function(curvature) {
      if (all(curvature >= 0)) curvature else stop("unexpected missing value")
    }
  ))
})

test_that("a search started past a limit of the link finds the mode", {
  sites <- 3L
  units_m <- rep(20, sites)
  latent_map <- make_latent_map(cbind(seq_len(sites)), sites)
  precision <- diag(0.5, sites)
  loglog_link <- make_link(
    function(eta) exp(-exp(-eta)),
    function(eta) exp(-eta - exp(-eta)),
    function(eta) exp(-eta - exp(-eta)) * (exp(-eta) - 1)
  )
  cases <- list(
    list(link = cloglog_link, y = c(15, 20, 20), mu = rep(5, sites)),
    list(link = loglog_link, y = c(3, 0, 0), mu = rep(-40, sites))
  )
  for (case in cases) {
    # The start has no slope to follow, because the mean is at its limit.
    expect_true(any(case$link$inv(case$mu) %in% c(0, 1)))
    found <- adaptive_laplace_mode(case$mu, case$y, units_m, latent_map,
                                   precision, "binomial", case$link,
                                   numeric(sites))
    terms <- laplace_observation_terms(case$mu + found$mode, case$y, units_m,
                                       "binomial", case$link)
    expect_false(found$reached_limit)
    expect_equal(terms$score, as.numeric(precision %*% found$mode),
                 tolerance = 1e-7)
    expect_gt(abs(found$mode[[1L]]), 3)
    # maxim_integrand() has its own search, from zero.
    direct <- maxim_integrand(
      y = case$y, denominator = units_m, mu = case$mu,
      Sigma = solve(precision), ID_coords = seq_len(sites),
      family = "binomial", invlink = case$link[c("inv", "d1", "d2")]
    )
    expect_equal(direct$mode, found$mode, tolerance = 1e-4)
  }
})

test_that("moving a start off a limit leaves a start away from the limits alone", {
  start <- c(0.5, -0.5, 1)
  expect_identical(
    move_start_off_limits(
      start, function(point) point, c(3, 0, 20), rep(20, 3L), "binomial",
      cloglog_link, 1:3
    ),
    start
  )
  # A site at the limit without a failure, or without a success, stays there.
  expect_identical(
    move_start_off_limits(
      c(50, -50), function(point) point, c(20, 0), c(20, 20), "binomial",
      probit_link, 1:2
    ),
    c(50, -50)
  )
})

test_that("maxim_integrand() is not stopped by a probability that rounds to one", {
  sites <- 5L
  distances <- as.matrix(dist(seq_len(sites)))
  Sigma <- 3 * exp(-distances / 2)
  y <- rep(20, sites)
  mu <- rep(4, sites)
  expect_identical(cloglog_link$inv(mu), rep(1, sites))
  found <- maxim_integrand(
    y = y, denominator = y, mu = mu, Sigma = Sigma,
    ID_coords = seq_len(sites), family = "binomial",
    invlink = cloglog_link[c("inv", "d1", "d2")]
  )
  # Every trial is positive, so the likelihood is flat and the prior decides.
  expect_equal(found$mode, rep(0, sites), tolerance = 1e-6)
})

test_that("the sampler is not stopped by a probability that rounds to one", {
  sites <- 5L
  Sigma <- diag(2, sites)
  y <- rep(20, sites)
  mu <- rep(3, sites)
  samples <- laplace_sampling_mcmc(
    y = y, units_m = y, mu = mu, Sigma = Sigma, ID_coords = seq_len(sites),
    family = "binomial", invlink = cloglog_link[c("inv", "d1", "d2")],
    control_mcmc = set_control_mcmc(n_sim = 3000, burnin = 500, thin = 1),
    seed = 11L
  )
  # All the trials are positive, so the conditional distribution above the
  # rounding point at 3.6 is the prior's, with about 37% of its mass.
  above <- mean(sweep(samples$samples$S, 2L, mu, "+") > 3.7)
  expect_gt(above, 0.2)
})

test_that("cloglog at high prevalence is certified", {
  skip_on_cran()
  set.seed(4)
  sites <- 60L
  coords <- cbind(runif(sites), runif(sites)) * 100
  field <- as.numeric(
    rnorm(sites) %*% chol(2 * exp(-as.matrix(dist(coords)) / 15))
  )
  inverse_link <- function(eta) -expm1(-exp(eta))
  survey <- sf::st_as_sf(
    data.frame(y = rbinom(sites, 20, inverse_link(-0.3 + field)), m = 20,
               east = coords[, 1] * 1000, north = coords[, 2] * 1000),
    coords = c("east", "north"), crs = 32632
  )
  fit <- glgpm(
    y ~ gp(), denominator = m, data = survey, family = "binomial",
    invlink = inverse_link, distance_units = "km", return_samples = FALSE,
    messages = FALSE, seed = 1L
  )
  expect_identical(attr(fit, "optimizer")$mcml_outcome, "certified")
})
