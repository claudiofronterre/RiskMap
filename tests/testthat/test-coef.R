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

  object$re <- list("test" = 1)

  coef(object)

  object$re <- list("test" = 1,
                    "test2" = 1,
                    "test3" = 1)

  coef(object)

})
