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

  fit <- glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial",
              denominator = denominator, invlink = invlink, control_mcmc = control_mcmc,
              messages = FALSE)

  expect_s3_class(fit, "RiskMap")
  expect_true(all(is.finite(unlist(fit$estimate))))
})
