test_that("dast fits without error", {
  fit <- dast(y ~ gp(), data = dast_data, den = denominator, time = survey_time,
             mda_times = dast_mda_times, int_mat = dast_int_mat, power_val = 1,
             control_mcmc = control_mcmc, messages = FALSE)

  expect_s3_class(fit, "RiskMap")
})
