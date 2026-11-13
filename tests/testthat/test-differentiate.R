## differentiate() (#90): auto-derives a missing d1/d2 for a user-supplied
## invlink, preferring symbolic differentiation (Deriv::Deriv()) and falling
## back to numerical differentiation (numDeriv::grad()) when Deriv can't
## process the function body.

test_that("uses the symbolic derivative when Deriv can differentiate the function", {
  f <- function(eta) 1 / (1 + exp(-eta))
  d <- differentiate(f, 5)

  eta <- c(-2, -1, 0, 1, 2)
  p <- f(eta)
  analytic <- p * (1 - p)

  expect_equal(d(eta), analytic, tolerance = 1e-8)
  ## symbolic result should be near machine precision, not a finite-difference
  ## approximation
  expect_lt(max(abs(d(eta) - analytic)), 1e-10)
})

test_that("falls back to numerical differentiation when Deriv errors on the function", {
  ## Deriv::Deriv() has no derivative rule for approxfun()'s returned closure
  ## and raises a hard error - the case this fallback exists for
  f <- function(eta) approxfun(c(-10, 10), c(0, 1))(eta)
  expect_error(Deriv::Deriv(function(eta) f(eta), "eta"))

  d <- differentiate(f, 5)
  expect_equal(d(2), 0.05, tolerance = 1e-6)
  ## stays vectorised over eta
  expect_equal(d(c(-2, 0, 2)), c(0.05, 0.05, 0.05), tolerance = 1e-6)
})

test_that("differentiates whichever d1 was actually produced, for d2", {
  f <- function(eta) 1 / (1 + exp(-eta))
  d1 <- differentiate(f, 5)
  d2 <- differentiate(d1, 5)

  eta <- c(-2, -1, 0, 1, 2)
  p <- f(eta)
  analytic_d2 <- p * (1 - p) * (1 - 2 * p)

  expect_equal(d2(eta), analytic_d2, tolerance = 1e-6)
})

test_that("a custom invlink that Deriv cannot differentiate still fits via the numeric fallback", {
  ## approxfun()'s returned closure has no Deriv derivative rule, so this
  ## exercises the same fallback as the unit tests above, end-to-end through
  ## glgpm(). Bounded to [0, 1], so paired with "binomial", not "poisson".
  ## Built once and closed over, rather than rebuilt on every call: the
  ## numeric fallback's Richardson extrapolation calls `inv()` many times per
  ## optimizer iteration, and re-running approxfun()'s table setup that often
  ## made this test take ~30s.
  af <- approxfun(c(-50, 50), c(0.01, 0.99), rule = 2)
  invlink <- list(inv = function(eta) af(eta))

  fit <- suppressWarnings(glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial",
              denominator = denominator, invlink = invlink, control_mcmc = control_mcmc,
              seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
              messages = FALSE))

  expect_s3_class(fit, "RiskMap")
  expect_true(all(is.finite(unlist(fit$estimate))))
})

test_that("the numeric fallback matches numDeriv on each value, in one pass", {
  ## A bare pnorm takes the numeric fallback, which is why it fitted slowly.
  ## Equality with numDeriv to the last bit shows that the fallback ran.
  d <- differentiate(pnorm, 5)
  eta <- c(-30, -8.3, -1, 0, 1e-12, 0.4, 3.6, 8.3)
  expect_identical(d(eta), vapply(eta, function(z) {
    numDeriv::grad(pnorm, z, method.args = list(zero.tol = 1))
  }, numeric(1)))

  ## Every element shares each evaluation of the link: 2 x 4 calls in all.
  af <- approxfun(c(-50, 50), c(0.01, 0.99), rule = 2)
  calls <- 0L
  counted <- function(eta) {
    calls <<- calls + 1L
    af(eta)
  }
  d_counted <- differentiate(counted, 5)
  calls <- 0L
  values <- d_counted(seq(-3, 3, length.out = 200))
  expect_identical(calls, 8L)
  expect_equal(values, rep(0.98 / 100, 200), tolerance = 1e-8)
})

test_that("a function that is not elementwise is differentiated value by value", {
  af <- approxfun(c(-50, 50), c(0.01, 0.99), rule = 2)
  f <- function(eta) af(eta - mean(eta))
  d <- differentiate(f, 5)
  eta <- c(-1, 0, 2)
  expect_identical(d(eta), vapply(eta, function(z) {
    numDeriv::grad(f, z, method.args = list(zero.tol = 1))
  }, numeric(1)))
})

test_that("numerical derivatives of a link are accurate near a linear predictor of 0", {
  ## None of these is differentiated by Deriv, so d1 is numerical and d2 is
  ## the numerical derivative of d1. With numDeriv's default zero tolerance
  ## the steps fell to about 2e-9 just above 1.8e-5, and d2 was wrong by up to
  ## several hundred there (-40 at 1e-4 for the imperfect test, where it is
  ## -1.1e-5).
  links <- list(
    imperfect = list(inv = function(eta) 0.01 + 0.89 * plogis(eta),
                     d1 = function(eta) 0.89 * dlogis(eta),
                     d2 = function(eta) {
                       p <- plogis(eta)
                       0.89 * p * (1 - p) * (1 - 2 * p)
                     }),
    logit = list(inv = plogis, d1 = dlogis, d2 = function(eta) {
      p <- plogis(eta)
      p * (1 - p) * (1 - 2 * p)
    }),
    probit = list(inv = pnorm, d1 = dnorm,
                  d2 = function(eta) -eta * dnorm(eta)),
    t4 = list(inv = function(eta) pt(eta, 4), d1 = function(eta) dt(eta, 4),
              d2 = function(eta) -5 * eta * dt(eta, 4) / (4 + eta^2)),
    cauchit = list(inv = pcauchy, d1 = dcauchy,
                   d2 = function(eta) -2 * eta / (pi * (1 + eta^2)^2))
  )
  near_zero <- c(0, 1e-8, 1e-6, 2e-5, 1e-4, 1e-3, 1e-2)
  eta <- sort(unique(c(seq(-30, 30, by = 0.25), near_zero, -near_zero)))
  for (name in names(links)) {
    link <- links[[name]]
    derived <- make_inverse_link_functions("binomial", link$inv, 5)
    probability <- link$inv(eta)
    error_d1 <- abs(derived$d1(eta) - link$d1(eta))
    error_d2 <- abs(derived$d2(eta) - link$d2(eta))
    expect_lt(max(error_d1), 1e-9, label = paste(name, "d1"))
    expect_lt(max(error_d2), 1e-5, label = paste(name, "d2"))
    ## The curvature of the log-likelihood divides d2 by the probability, so
    ## below one half the error is also small relative to it.
    lower <- probability < 0.5
    expect_lt(max(error_d2[lower] / probability[lower]), 1e-4,
              label = paste(name, "d2 relative to the probability"))
  }
})

test_that("the numeric fallback stops on a link that returns NA", {
  expect_error(
    compute_richardson_derivative(function(eta) ifelse(eta > 1, NA, eta), c(0, 2)),
    "returns `NA`", class = "riskmap_error_input"
  )
})
