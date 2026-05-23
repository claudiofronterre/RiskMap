test_that("summary labels spatio-temporal nugget and temporal scale correctly", {
  object <- list(
    D = matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "(Intercept)")),
    estimate = c(
      beta = 0,
      sigma2 = log(2),
      phi = log(3),
      tau2 = log(0.5 / 2),
      alpha = qlogis(0.2),
      gamma = log(4),
      psi = log(6)
    ),
    covariance = diag(rep(0.01, 7)),
    re = list(),
    family = "binomial",
    power_val = 1,
    fix_alpha = NULL,
    fix_tau2 = NULL,
    sst = TRUE,
    kappa = 0.5,
    log.lik = -1,
    cov_offset = NULL,
    call = quote(dast(y ~ gp(x, y, t), data = data))
  )
  class(object) <- "RiskMap"

  out <- summary(object)

  expect_equal(rownames(out$sp),
               c("Spatial process var.", "Spatial corr. scale",
                 "Variance of the nugget", "Temporal corr. scale"))
  expect_equal(out$sp["Spatial process var.", "Estimate"], 2)
  expect_equal(out$sp["Spatial corr. scale", "Estimate"], 3)
  expect_equal(out$sp["Variance of the nugget", "Estimate"], 0.5)
  expect_equal(out$sp["Temporal corr. scale", "Estimate"], 6)
  expect_equal(
    unname(out$sp["Temporal corr. scale", c("Lower limit", "Upper limit")]),
    exp(log(6) + c(-1, 1) * stats::qnorm(0.975) * 0.1),
    tolerance = 1e-8
  )
})
