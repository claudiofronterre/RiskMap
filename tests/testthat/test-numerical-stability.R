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

test_that("design standardisation preserves linear predictors", {
  design <- cbind(
    "(Intercept)" = 1,
    elevation = c(1000, 1500, 2500, 4000),
    rainfall = c(0.01, 0.04, 0.02, 0.08)
  )
  beta <- c(0.5, -0.002, 8)
  scaled <- standardize_design_matrix(design)
  scaled_beta <- solve(scaled$coefficient_transform, beta)

  expect_equal(
    as.numeric(scaled$design %*% scaled_beta),
    as.numeric(design %*% beta)
  )
  expect_equal(unname(colMeans(scaled$design[, -1, drop = FALSE])), c(0, 0),
               tolerance = 1e-14)
  expect_equal(unname(apply(scaled$design[, -1, drop = FALSE], 2, sd)), c(1, 1))

  no_intercept <- design[, -1, drop = FALSE]
  scaled_no_intercept <- standardize_design_matrix(no_intercept)
  expect_equal(unname(scaled_no_intercept$center), c(0, 0))
  expect_equal(
    as.numeric(scaled_no_intercept$design %*%
                 solve(scaled_no_intercept$coefficient_transform, beta[-1])),
    as.numeric(no_intercept %*% beta[-1])
  )

  intercept_only <- standardize_design_matrix(design[, 1, drop = FALSE])
  expect_equal(intercept_only$design, design[, 1, drop = FALSE])
  expect_equal(unname(intercept_only$coefficient_transform), diag(1))

  expect_error(
    standardize_design_matrix(cbind("(Intercept)" = 1, constant = 2)),
    "constant.*constant"
  )
  expect_equal(
    sd(standardize_design_matrix(
      cbind("(Intercept)" = 1, tiny = (1:4) * 1e-12)
    )$design[, "tiny"]),
    1
  )
})

test_that("fixed-effect results are restored with the full Jacobian", {
  transform <- matrix(c(1, 0, -2, 0.5), nrow = 2)
  working_covariance <- matrix(
    c(2, 0.2, 0.4,
      0.2, 1, 0.3,
      0.4, 0.3, 3),
    nrow = 3
  )
  result <- list(
    estimate = list(beta = c("(Intercept)" = 3, x = 4), sigma2 = 0.5),
    covariance = working_covariance,
    grad_MLE = c(0.1, -0.2, 0.3)
  )
  jacobian <- diag(3)
  jacobian[1:2, 1:2] <- transform

  restored <- restore_fixed_effect_scale(result, transform)

  expect_equal(unname(restored$estimate$beta),
               as.numeric(transform %*% c(3, 4)))
  expect_equal(unname(restored$covariance),
               jacobian %*% working_covariance %*% t(jacobian))
  expect_equal(unname(restored$grad_MLE),
               as.numeric(solve(t(jacobian), result$grad_MLE)))
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

test_that("Gaussian covariance backend selection preserves correctness overrides", {
  expect_true(use_direct_gaussian_covariance(100, 50, 0.2, FALSE))
  expect_false(use_direct_gaussian_covariance(500, 50, 0.2, FALSE))
  expect_true(use_direct_gaussian_covariance(500, 50, 0, FALSE))
  expect_true(use_direct_gaussian_covariance(500, 50, 0.2, 0.4))
})

test_that("MCML convergence requires small standardised change and overlap", {
  control <- set_control_mcml(
    max_iterations = 3,
    tolerance = 0.05,
    min_relative_ess = 0.1
  )

  expect_true(mcml_update_converged(0.04, 0.2, control))
  expect_false(mcml_update_converged(0.06, 0.2, control))
  expect_false(mcml_update_converged(0.04, 0.05, control))
  expect_false(mcml_update_converged(NA_real_, 0.2, control))
})

test_that("MCML steps are shortened when the proposal lacks overlap", {
  result <- supported_importance_step(
    reference = 0,
    proposal = 1,
    log_weight_function = function(parameter) c(0, -20 * parameter),
    min_relative_ess = 0.75
  )

  expect_lt(result$step_fraction, 1)
  expect_gte(result$relative_importance_ess, 0.75)
  expect_lt(result$proposal_relative_importance_ess, 0.75)
})
