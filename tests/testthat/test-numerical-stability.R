test_that("log_mean_exp is stable and translation equivariant", {
  x <- c(-2, 0, 3)
  expect_equal(log_mean_exp(x), log(mean(exp(x))))

  large <- c(1000, 1001, 1002)
  expect_true(is.finite(log_mean_exp(large)))
  expect_equal(log_mean_exp(large) - 1000, log_mean_exp(large - 1000))
  expect_equal(log_mean_exp(rep(-Inf, 3)), -Inf)
})

test_that("normalise_log_weights remains finite on extreme scales", {
  log_weights <- c(1000, 1001, 1002)
  weights <- normalise_log_weights(log_weights)

  expect_true(all(is.finite(weights)))
  expect_equal(sum(weights), 1)
  expect_equal(weights, exp(log_weights - 1002) / sum(exp(log_weights - 1002)))
  expect_error(
    normalise_log_weights(rep(-Inf, 3)),
    "All Monte Carlo importance weights are non-finite"
  )
  expect_equal(importance_effective_sample_size(rep(0.25, 4)), 4)
  expect_equal(importance_effective_sample_size(c(1, 0, 0, 0)), 1)
})

test_that("softplus is stable in both tails", {
  ordinary <- c(-2, 0, 3)
  expect_equal(softplus(ordinary), log1p(exp(ordinary)))

  expect_equal(softplus(1000), 1000)
  expect_equal(softplus(-1000), 0)
  expect_true(all(is.finite(softplus(c(-1000, 1000)))))
})

test_that("linear_start_values uses a rank-aware least-squares solve", {
  design <- cbind(1, c(-1, 0, 1))
  response <- c(0, 1, 4)

  expect_equal(linear_start_values(design, response),
               unname(coef(lm(response ~ design[, 2]))))
  expect_error(
    linear_start_values(cbind(1, 1:3, 2 * (1:3)), response),
    "model matrix is rank deficient"
  )
})

test_that("safe_optimizer_objective penalises and counts invalid trials", {
  objective <- safe_optimizer_objective(function(x) {
    if (x < 0) stop("invalid trial")
    if (x == 0) return(Inf)
    x^2
  })

  expect_equal(objective(2), 4)
  expect_true(is.finite(objective(-1)))
  expect_true(is.finite(objective(0)))
  expect_equal(attr(objective, "diagnostics")$invalid_evaluations, 2L)
})
