# Tests for the summary() methods and their print() methods.
#
# The print methods only format summary objects, so their tests build the
# summaries by hand with fixed values. Snapshots of summaries from real fits
# would change with platform, BLAS or optimiser details rather than with the
# formatting being tested.

# summary.RiskMap ------------------------------------------------------------

test_that("estimates are consistent between coef and summary", {

  sum <- summary(gaussian_model)
  cof <- coef(gaussian_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_me, sum$me[1,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)

  sum <- summary(binomial_model)
  cof <- coef(binomial_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)

  sum <- summary(poisson_model)
  cof <- coef(poisson_model)

  expect_equal(cof$beta[1], sum$reg_coef[1,1], ignore_attr = TRUE)
  expect_equal(cof$beta[2], sum$reg_coef[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2, sum$sp[1,1], ignore_attr = TRUE)
  expect_equal(cof$phi, sum$sp[2,1], ignore_attr = TRUE)
  expect_equal(cof$sigma2_re, sum$ranef[1,1], ignore_attr = TRUE)
})

test_that("summary.RiskMap returns the expected components for each family", {
  gaussian_summary <- summary(gaussian_model)
  expect_s3_class(gaussian_summary, "summary.RiskMap")
  expect_setequal(names(gaussian_summary),
                  c("reg_coef", "me", "sp", "tau2", "ranef", "level",
                    "family", "kappa", "log_lik", "cov_offset_used", "aic",
                    "call", "link_name", "invlink_expression"))
  expect_identical(gaussian_summary$family, "gaussian")
  expect_identical(gaussian_summary$link_name, "identity")
  expect_identical(gaussian_summary$call, gaussian_model$call)
  expect_identical(gaussian_summary$log_lik, gaussian_model$log_lik)
  expect_identical(gaussian_summary$kappa, gaussian_model$kappa)
  expect_false(gaussian_summary$cov_offset_used)
  expect_equal(gaussian_summary$aic,
               2 * length(unlist(gaussian_model$estimate)) - 2 * gaussian_model$log_lik)
  expect_identical(rownames(gaussian_summary$ranef), "i (random eff. var.)")
  expect_identical(rownames(gaussian_summary$sp),
                   c("Spatial process var.", "Spatial corr. scale (km)"))

  binomial_summary <- summary(binomial_model)
  expect_identical(binomial_summary$link_name, "canonical (logit)")
  expect_null(binomial_summary$me)
  expect_null(binomial_summary$aic)

  poisson_summary <- summary(poisson_model)
  expect_identical(poisson_summary$link_name, "canonical (log)")
  expect_null(poisson_summary$me)
  expect_null(poisson_summary$aic)
})

test_that("summary.RiskMap computes standard errors and confidence intervals", {
  fit_summary <- summary(binomial_model)
  estimate <- unlist(binomial_model$estimate)
  standard_error <- sqrt(diag(binomial_model$covariance))
  is_beta <- startsWith(names(estimate), "beta.")
  z_crit <- qnorm(0.975)

  reg_coef <- fit_summary$reg_coef
  expect_identical(colnames(reg_coef),
                   c("Estimate", "Lower limit", "Upper limit", "StdErr",
                     "z.value", "p.value"))
  expect_identical(rownames(reg_coef), colnames(binomial_model$D))
  expect_equal(reg_coef[, "Estimate"], estimate[is_beta], ignore_attr = TRUE)
  expect_equal(reg_coef[, "StdErr"], standard_error[is_beta], ignore_attr = TRUE)
  expect_equal(reg_coef[, "Lower limit"],
               estimate[is_beta] - z_crit * standard_error[is_beta],
               ignore_attr = TRUE)
  expect_equal(reg_coef[, "z.value"],
               estimate[is_beta] / standard_error[is_beta], ignore_attr = TRUE)
  expect_equal(reg_coef[, "p.value"],
               2 * pnorm(-abs(reg_coef[, "z.value"])))

  # covariance parameters are estimated on the log scale, so their intervals
  # are log-normal
  sp_index <- match(c("sigma2", "phi"), names(estimate))
  expect_equal(fit_summary$sp[, "Estimate"], exp(estimate[sp_index]),
               ignore_attr = TRUE)
  expect_equal(fit_summary$sp[, "Lower limit"],
               exp(estimate[sp_index] - z_crit * standard_error[sp_index]),
               ignore_attr = TRUE)
  expect_equal(fit_summary$sp[, "Upper limit"],
               exp(estimate[sp_index] + z_crit * standard_error[sp_index]),
               ignore_attr = TRUE)

  re_index <- match("sigma2_re.i", names(estimate))
  expect_equal(fit_summary$ranef[, "Upper limit"],
               exp(estimate[re_index] + z_crit * standard_error[re_index]),
               ignore_attr = TRUE)
})

test_that("summary.RiskMap respects level", {
  default_summary <- summary(gaussian_model)
  narrow_summary <- summary(gaussian_model, level = 0.9)

  expect_identical(narrow_summary$level, 0.9)
  expect_equal(narrow_summary$reg_coef[, "Estimate"],
               default_summary$reg_coef[, "Estimate"])
  default_width <- default_summary$reg_coef[, "Upper limit"] -
    default_summary$reg_coef[, "Lower limit"]
  narrow_width <- narrow_summary$reg_coef[, "Upper limit"] -
    narrow_summary$reg_coef[, "Lower limit"]
  expect_equal(narrow_width / default_width,
               rep(qnorm(0.95) / qnorm(0.975), 2), ignore_attr = TRUE)
})

test_that("summary.RiskMap reports fixed and estimated nugget and measurement error", {
  # gaussian_offset_model: gp(nugget = TRUE), fix_var_me = 0, offset
  offset_summary <- summary(gaussian_offset_model)
  expect_true(offset_summary$tau2)
  expect_true(offset_summary$cov_offset_used)
  expect_identical(offset_summary$me, 0)
  expect_identical(rownames(offset_summary$sp)[3], "Variance of the nugget")
  # the nugget is estimated as a ratio to sigma2 but reported as a variance
  expect_equal(offset_summary$sp["Variance of the nugget", "Estimate"],
               coef(gaussian_offset_model)$tau2)

  # gaussian_model: no nugget, estimated measurement error
  gaussian_summary <- summary(gaussian_model)
  expect_identical(gaussian_summary$tau2, 0)
  expect_identical(rownames(gaussian_summary$me), "Measurement error var.")
  expect_equal(gaussian_summary$me[, "Estimate"],
               coef(gaussian_model)$sigma2_me, ignore_attr = TRUE)
})

# A summary.RiskMap for a gaussian model; 'changes' replaces fields, and
# setting a field to NULL removes it.
build_summary_riskmap <- function(changes = list()) {
  reg_coef <- cbind(
    Estimate      = c(1.2, -0.4),
    "Lower limit" = c(0.8, -0.9),
    "Upper limit" = c(1.6, 0.1),
    StdErr        = c(0.2, 0.25),
    z.value       = c(6, -1.6),
    p.value       = c(1.97e-09, 0.1096)
  )
  rownames(reg_coef) <- c("(Intercept)", "cov")

  me <- cbind(Estimate = 0.3, "Lower limit" = 0.1, "Upper limit" = 0.9)
  rownames(me) <- "Measurement error var."

  sp <- cbind(
    Estimate      = c(0.5, 0.2, 0.1),
    "Lower limit" = c(0.25, 0.05, 0.02),
    "Upper limit" = c(1, 0.8, 0.5)
  )
  rownames(sp) <- c("Spatial process var.", "Spatial corr. scale (km)",
                    "Variance of the nugget")

  ranef <- cbind(Estimate = 0.4, "Lower limit" = 0.2, "Upper limit" = 0.8)
  rownames(ranef) <- "i (random eff. var.)"

  summary_object <- list(
    reg_coef = reg_coef,
    me = me,
    sp = sp,
    tau2 = TRUE,
    ranef = ranef,
    level = 0.95,
    family = "gaussian",
    kappa = 0.5,
    log_lik = -12.3456,
    cov_offset_used = TRUE,
    aic = 38.6912,
    call = quote(glgpm(formula = y ~ cov + gp(nugget = TRUE) + re(i),
                       data = gaussian_data, family = "gaussian")),
    link_name = "identity",
    invlink_expression = "Inverse link function = x"
  )
  summary_object <- modifyList(summary_object, changes)
  class(summary_object) <- "summary.RiskMap"
  summary_object
}

test_that("print.summary.RiskMap formats a gaussian model with every optional section", {
  summary_object <- build_summary_riskmap()

  expect_snapshot(print(summary_object))
  expect_output(expect_invisible(print(summary_object)))
})

test_that("print.summary.RiskMap reports a fixed measurement error variance", {
  summary_object <- build_summary_riskmap(list(me = 0.1, call = NULL))

  expect_snapshot(print(summary_object))
})

test_that("print.summary.RiskMap formats binomial and poisson models", {
  sp_without_nugget <- build_summary_riskmap()$sp[1:2, ]

  binomial_summary <- build_summary_riskmap(list(
    family = "binomial", me = NULL, sp = sp_without_nugget, tau2 = 0,
    ranef = NULL, cov_offset_used = FALSE, aic = NULL, call = NULL,
    link_name = "canonical (logit)",
    invlink_expression = "Inverse link function = 1 / (1 + exp(-x))"
  ))
  expect_snapshot(print(binomial_summary))

  poisson_summary <- build_summary_riskmap(list(
    family = "poisson", me = NULL, sp = sp_without_nugget, tau2 = 0.25,
    ranef = NULL, cov_offset_used = FALSE, aic = NULL, call = NULL,
    link_name = "canonical (log)",
    invlink_expression = "Inverse link function = exp(x)"
  ))
  expect_snapshot(print(poisson_summary))
})

test_that("print.summary.RiskMap labels a summary of a real fit", {
  expect_output(print(summary(gaussian_model)), "Linear geostatistical model")
  expect_output(print(summary(gaussian_model)), "Unstructured random effects")
  expect_output(print(summary(binomial_model)), "Binomial geostatistical model")
  expect_output(print(summary(poisson_model)), "Poisson geostatistical model")
})

# summary.RiskMap_cross_validation -------------------------------------------

# A RiskMap_cross_validation with two models, two metrics and two test sets
# of sizes 1 and 3, so the overall scores are weighted by test set size.
build_cross_validation <- function() {
  list(model = list(
    model_a = list(metric = list(CRPS = list(0.2, c(0.3, 0.4, 0.5)),
                                 SCRPS = list(-0.5, c(-0.6, -0.7, -0.8)))),
    model_b = list(metric = list(CRPS = list(0.6, c(0.1, 0.2, 0.3)),
                                 SCRPS = list(-1, c(-0.2, -0.3, -0.4))))
  ))
}

test_that("summary.RiskMap_cross_validation weights test sets by their size", {
  cross_validation_summary <- summary.RiskMap_cross_validation(build_cross_validation())

  expect_s3_class(cross_validation_summary, "summary.RiskMap_cross_validation")
  expect_identical(dimnames(cross_validation_summary),
                   list(c("model_a", "model_b"), c("CRPS", "SCRPS")))
  expect_equal(cross_validation_summary["model_a", "CRPS"], (0.2 + 3 * 0.4) / 4)
  expect_equal(cross_validation_summary["model_b", "SCRPS"], (-1 + 3 * -0.3) / 4)

  test_set_means <- attr(cross_validation_summary, "test_set_means")
  expect_named(test_set_means, c("model_a", "model_b"))
  # scores per test set are stored by metric position, in column order
  expect_equal(test_set_means$model_a[[1]], c(0.2, 0.4))
  expect_equal(attr(cross_validation_summary, "overall_averages"),
               colMeans(cross_validation_summary))
  expect_true(attr(cross_validation_summary, "view_all"))

  expect_false(attr(summary.RiskMap_cross_validation(build_cross_validation(),
                                                     view_all = FALSE),
                    "view_all"))
})

test_that("summary.RiskMap_cross_validation errors when no metrics were computed", {
  no_metrics <- list(model = list(model_a = list(metric = list())))
  expect_error(summary.RiskMap_cross_validation(no_metrics),
               "No metrics of predictive performance were computed")
})

test_that("summary.RiskMap_cross_validation summarises assess_prediction output", {
  cross_validation_summary <- summary(cross_validation)
  expect_identical(rownames(cross_validation_summary), c("model_a", "model_b"))
  expect_true(all(is.finite(cross_validation_summary)))
})

test_that("print.summary.RiskMap_cross_validation formats scores by model", {
  scores <- matrix(c(0.3, 0.45, -0.8, -0.95), nrow = 2,
                   dimnames = list(c("model_a", "model_b"), c("CRPS", "SCRPS")))
  attr(scores, "test_set_means") <- list(
    model_a = list(CRPS = c(0.25, 0.35), SCRPS = c(-0.75, -0.85)),
    model_b = list(CRPS = c(0.4, 0.5), SCRPS = c(-0.9, -1))
  )
  attr(scores, "overall_averages") <- colMeans(scores)
  class(scores) <- "summary.RiskMap_cross_validation"

  attr(scores, "view_all") <- TRUE
  expect_snapshot(print(scores))
  expect_output(expect_invisible(print(scores)))

  attr(scores, "view_all") <- FALSE
  expect_snapshot(print(scores))
})

# summary.RiskMap_assess_simulation -------------------------------------------

# A RiskMap_assess_simulation with two simulations of one model; summary()
# of it is deterministic, so it is used to build the object being printed.
build_assess_simulation <- function(spatial_scale) {
  classify_result <- function(sensitivity, cc) {
    data.frame(Class = c("low", "high"), Sensitivity = sensitivity,
               Specificity = c(0.8, 0.7), PPV = c(0.6, 0.5),
               NPV = c(0.9, 0.85), CC = cc)
  }
  scale_results <- list(
    mse = matrix(c(0.12, 0.18), nrow = 1,
                 dimnames = list("model", c("sim_1", "sim_2"))),
    classify = list(
      model = list(by_cat = list(classify_result(c(0.5, 0.6), c(0.7, 0.7)),
                                 classify_result(c(0.7, 0.4), c(0.8, 0.8))),
                   CC = c(0.7, 0.8)),
      Class = factor(c("low", "high"), levels = c("low", "high"))
    )
  )
  pred_objective <- rep(list(scale_results), length(spatial_scale))
  names(pred_objective) <- spatial_scale

  structure(
    list(pred_objective = pred_objective, n_sim = 2L,
         spatial_scale = spatial_scale),
    class = "RiskMap_assess_simulation"
  )
}

test_that("classification summaries handle one simulation and metric-wise missingness #108", {
  one_result <- data.frame(
    Class = c("low", "high"),
    Sensitivity = c(0.5, NA),
    Specificity = c(NA, 0.75),
    PPV = c(0.4, NA),
    NPV = c(NA, 0.8),
    CC = c(0.6, 0.7)
  )
  one_simulation <- structure(
    list(
      pred_objective = list(
        grid = list(
          classify = list(
            model = list(by_cat = list(one_result), CC = 0.65),
            Class = factor(c("low", "high"), levels = c("low", "high"))
          )
        )
      ),
      n_sim = 1L,
      spatial_scale = "grid"
    ),
    class = "RiskMap_assess_simulation"
  )

  expect_warning(summary_one <- summary(one_simulation), "one simulation")
  model_one <- summary_one$grid$classify$model
  expect_identical(model_one$classify_res$Sensitivity, c(0.5, NA))
  expect_identical(model_one$n_valid$Sensitivity, c(1L, 0L))
  expect_identical(model_one$cc_summary$mean, 0.65)
  expect_true(is.na(model_one$cc_summary$sd))
  expect_true(is.na(model_one$cc_summary$lower))
  expect_output(print(summary_one), "uncertainty unavailable")

  second_result <- one_result
  second_result$Sensitivity <- c(NA, 0.25)
  second_result$Specificity <- c(0.5, NA)
  two_simulations <- one_simulation
  two_simulations$n_sim <- 2L
  two_simulations$pred_objective$grid$classify$model$by_cat <-
    list(one_result, second_result)
  two_simulations$pred_objective$grid$classify$model$CC <- c(0.65, 0.75)

  summary_two <- summary(two_simulations)$grid$classify$model
  expect_identical(summary_two$classify_res$Sensitivity, c(0.5, 0.25))
  expect_identical(summary_two$n_valid$Sensitivity, c(1L, 1L))
  expect_identical(summary_two$n_valid$Specificity, c(1L, 1L))
  expect_identical(summary_two$cc_summary$n_valid, 2L)
  expect_identical(summary_two$cc_summary$n_sim, 2L)
})

test_that("summary.RiskMap_assess_simulation summarises MSE and classification", {
  assessment <- build_assess_simulation(c("grid", "area"))
  assessment$pred_objective$grid$mse[1, 2] <- NA

  expect_error(summary.RiskMap_assess_simulation(list()), "RiskMap_assess_simulation")

  simulation_summary <- summary(assessment)
  expect_s3_class(simulation_summary, "summary.RiskMap_assess_simulation")
  expect_named(simulation_summary, c("grid", "area"))

  # non-finite MSE values are dropped from the mean and counted as invalid
  expect_identical(simulation_summary$grid$mse$n_sim, 2L)
  expect_identical(simulation_summary$grid$mse$n_valid, 1)
  expect_identical(simulation_summary$grid$mse$MSE_mean, 0.12)
  expect_true(is.na(simulation_summary$grid$mse$MSE_sd))
  expect_equal(simulation_summary$area$mse$MSE_mean, 0.15)
  expect_equal(simulation_summary$area$mse$MSE_sd, sd(c(0.12, 0.18)))

  model_summary <- simulation_summary$grid$classify$model
  expect_identical(model_summary$classify_res$Class, c("low", "high"))
  expect_equal(model_summary$classify_res$Sensitivity, c(0.6, 0.5))
  expect_equal(model_summary$cc_summary$mean, 0.75)
  expect_equal(model_summary$cc_summary$lower,
               unname(quantile(c(0.7, 0.8), 0.025)))

  assessment$pred_objective$grid$mse <- c(0.12, 0.18)
  expect_error(summary(assessment), "mse_data must be a matrix")
})

test_that("print.summary.RiskMap_assess_simulation formats one spatial scale", {
  summary_object <- summary(build_assess_simulation("grid"))

  expect_snapshot(print(summary_object))
  expect_output(expect_invisible(print(summary_object)))
})

test_that("print.summary.RiskMap_assess_simulation labels each spatial scale", {
  summary_object <- summary(build_assess_simulation(c("grid", "area")))

  expect_snapshot(print(summary_object))
})
