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
  expect_identical(d(eta), vapply(eta, function(z) numDeriv::grad(pnorm, z),
                                  numeric(1)))

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
  expect_identical(d(eta), vapply(eta, function(z) numDeriv::grad(f, z),
                                  numeric(1)))
})

test_that("the numeric fallback stops on a link that returns NA", {
  expect_error(
    compute_richardson_derivative(function(eta) ifelse(eta > 1, NA, eta), c(0, 2)),
    "returns `NA`", class = "riskmap_error_input"
  )
})
