# What a fit promises (the interface review, F1, F2, F4, N1, N6, N9, O2, O4).

# The covariance on the scale of coef() by the delta method, with the
# Jacobian of the natural values in the working ones taken numerically.
numerical_vcov <- function(fit) {
  parameters <- fit_parameter_table(fit)
  rows <- estimated_rows(parameters)
  natural <- function(working) {
    values <- natural_values(parameters, working_estimate(parameters, working))
    unlist(values)[rows$reported]
  }
  jacobian <- numDeriv::jacobian(natural, unlist(fit$estimate))
  jacobian <- jacobian[, rows$index, drop = FALSE]
  covariance <- fit$covariance[rows$index, rows$index, drop = FALSE]
  out <- jacobian %*% covariance %*% t(jacobian)
  dimnames(out) <- list(rows$reported, rows$reported)
  out
}

test_that("glgpm() takes distance_crs, with seed and messages last", {
  expect_identical(
    names(formals(glgpm)),
    c("formula", "data", "family", "invlink", "denominator", "distance_crs",
      "distance_units", "start_pars", "fixed_pars", "control_mcml",
      "control_mcmc", "return_samples", "seed", "messages")
  )
  expect_identical(names(formals(summary.RiskMap)),
                   c("object", "level", "..."))
})

test_that("cores and messages come last, in that order, cores is an integer and profile() has no default seed (N6, N7)", {
  expect_identical(names(formals(check_fit)),
                   c("object", "extra", "profile", "cores", "messages"))
  expect_identical(utils::tail(names(formals(parametric_bootstrap)), 4L),
                   c("seed", "level", "cores", "messages"))
  expect_identical(utils::tail(names(formals(profile.RiskMap)), 5L),
                   c("seed", "level", "cores", "messages", "..."))
  expect_identical(formals(parametric_bootstrap)$cores, 1L)
  expect_identical(formals(check_fit)$cores, 1L)
  expect_null(formals(profile.RiskMap)$seed)
})

test_that("logLik(), nobs(), AIC() and BIC() read a fit", {
  log_likelihood <- logLik(gaussian_model)
  expect_s3_class(log_likelihood, "logLik")
  expect_identical(as.numeric(log_likelihood), gaussian_model$log_lik)
  expect_identical(attr(log_likelihood, "df"),
                   length(unlist(coef(gaussian_model))))
  expect_identical(nobs(gaussian_model), 10L)
  expect_identical(attr(log_likelihood, "nobs"), 10L)
  expect_equal(AIC(gaussian_model), summary(gaussian_model)$aic)
  expect_equal(BIC(gaussian_model),
               -2 * gaussian_model$log_lik +
                 log(10) * attr(log_likelihood, "df"))
  # A held parameter is not counted.
  held <- logLik(gaussian_offset_model)
  expect_identical(attr(held, "df"),
                   length(unlist(coef(gaussian_offset_model))) - 1L)
  expect_equal(AIC(gaussian_offset_model),
               summary(gaussian_offset_model)$aic)
  expect_identical(nobs(binomial_model), 10L)
  expect_identical(as.numeric(logLik(binomial_model)), binomial_model$log_lik)
})

test_that("vcov() is the covariance on the scale of coef()", {
  for (fit in list(gaussian_model, gaussian_offset_model, binomial_model)) {
    covariance <- vcov(fit)
    rows <- estimated_rows(fit_parameter_table(fit))
    expect_identical(rownames(covariance), rows$reported)
    expect_identical(colnames(covariance), rows$reported)
    expect_true(all(rownames(covariance) %in% names(unlist(coef(fit)))))
    expect_equal(covariance, numerical_vcov(fit), tolerance = 1e-6)
    # The coefficients are on the same scale, so their standard errors are
    # those of summary().
    beta <- grep("^beta\\.", rownames(covariance), value = TRUE)
    expect_equal(sqrt(diag(covariance)[beta]),
                 summary(fit)$reg_coef[, "StdErr"], ignore_attr = TRUE)
  }
  # The nugget's variance is the product of its ratio and the spatial
  # variance, so its covariance with the spatial variance is not zero.
  expect_true("tau2" %in% rownames(vcov(gaussian_offset_model)))
  expect_false("sigma2_me" %in% rownames(vcov(gaussian_offset_model)))
  expect_true(vcov(gaussian_offset_model)["tau2", "sigma2"] != 0)
})

test_that("vcov() maps the size of a negative binomial fit and leaves NA rows", {
  # A negative binomial fit made from a Poisson one, with the size on the
  # scale of log(1 / theta) and a covariance that covers it.
  fit <- poisson_model
  fit$family <- "negative_binomial"
  fit$estimate <- c(fit$estimate[c("beta", "sigma2", "phi")],
                    list(theta = -log(2)), fit$estimate["sigma2_re"])
  names_flat <- names(unlist(fit$estimate))
  set.seed(3)
  root <- matrix(rnorm(length(names_flat)^2), length(names_flat)) / 4
  fit$covariance <- crossprod(root) + diag(0.1, length(names_flat))
  dimnames(fit$covariance) <- list(names_flat, names_flat)
  covariance <- vcov(fit)
  expect_true("theta" %in% rownames(covariance))
  expect_equal(covariance, numerical_vcov(fit), tolerance = 1e-6)
  # The size falls as its working value rises.
  theta_row <- which(names_flat == "theta")
  expect_equal(covariance["theta", "phi"],
               -2 * coef(fit)$phi * fit$covariance[theta_row, "phi"])
  # A parameter without a standard error has an NA row and column.
  fit$covariance["phi", ] <- NA
  fit$covariance[, "phi"] <- NA
  covariance <- vcov(fit)
  expect_true(all(is.na(covariance["phi", ])))
  expect_true(all(is.na(covariance[, "phi"])))
  expect_false(anyNA(covariance[rownames(covariance) != "phi",
                                colnames(covariance) != "phi"]))
})

test_that("vcov() of a fit that holds every parameter is empty", {
  fit <- glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
               fixed_pars = list(beta = c(1, 0.3), sigma2 = 0.5, phi = 1.2,
                                 sigma2_me = 0.1),
               messages = FALSE)
  expect_identical(dim(vcov(fit)), c(0L, 0L))
  expect_identical(attr(logLik(fit), "df"), 0L)
})

test_that("the accessors and the print and summary of a fit stop on unused arguments", {
  expect_input_error(print(gaussian_model, color = TRUE), "color")
  expect_input_error(summary(gaussian_model, conf_level = 0.9), "conf_level")
  expect_error(summary(gaussian_model, conf_level = 0.9),
               "confidence level is given in")
  expect_input_error(print(summary(gaussian_model), digits = 3), "digits")
  expect_input_error(vcov(gaussian_model, complete = TRUE), "complete")
  expect_input_error(logLik(gaussian_model, REML = TRUE), "REML")
  expect_input_error(nobs(gaussian_model, use.fallback = TRUE), "use.fallback")
  expect_input_error(print(negative_binomial(), digits = 3), "digits")
  # level is matched by position now that it comes before the dots.
  expect_identical(summary(gaussian_model, 0.9)$level, 0.9)
})

test_that("every result records the version of RiskMap that made it", {
  stamp <- riskmap_version_stamp()
  expect_named(stamp, c("version", "sha"))
  expect_s3_class(stamp$version, "package_version")
  expect_identical(stamp$version, utils::packageVersion("RiskMap"))
  expect_type(stamp$sha, "character")
  expect_length(stamp$sha, 1L)
  stamped <- function(object) {
    expect_identical(attr(object, "riskmap_version"), stamp)
  }
  stamped(gaussian_model)
  stamped(binomial_model)
  stamped(predict(gaussian_intercept_model, grid, messages = FALSE))
  stamped(predict(gaussian_intercept_model, grid, boundaries = areal,
                  aggregate = "mean", messages = FALSE))
  stamped(cross_validation)
  stamped(cross_validation$test_sets)
  stamped(create_test_sets(gaussian_data))
  stamped(simulation)
  stamped(simulation$model)
  stamped(variogram(gaussian_data, "y", n_permutations = 19, seed = 1,
                    distance_units = "m"))
  stamped(check_fit(gaussian_intercept_model, messages = FALSE))
  stamped(profile(gaussian_intercept_model, parameters = "phi", n_points = 5,
                  seed = 1, messages = FALSE))
  stamped(parametric_bootstrap(gaussian_intercept_model, n_sim = 2, seed = 1,
                               messages = FALSE))
  stamped(suppressWarnings(assess_extrapolation(gaussian_intercept_model, grid,
                                                spatial_threshold = 0.8)))
  without_re <- simulate_glgpm(gaussian_intercept_model,
                               what = c("data", "surface"),
                               prediction_grid = grid,
                               sample_locations = gaussian_data, seed = 1)
  stamped(without_re)
  stamped(assess_simulation(without_re,
                            models = list(model = gaussian_intercept_model),
                            spatial_scale = "location", objective = "mse",
                            messages = FALSE))
  # The stamp is not printed.
  expect_false(any(grepl("riskmap_version", capture.output(print(gaussian_model)))))
})

test_that("to_table() has a method for each object it tabulates", {
  expect_identical(to_table(gaussian_model, digits = 2),
                   to_table(summary(gaussian_model), digits = 2))
  expect_identical(to_table(cross_validation, digits = 2),
                   to_table(summary(cross_validation), digits = 2))
  expect_identical(
    to_table(cross_validation, reference = "model_a"),
    to_table(summary(cross_validation, reference = "model_a"))
  )
  expect_input_error(to_table(list(a = 1)), "object")
  expect_error(to_table(1), "RiskMap_cross_validation")
})

test_that("the engine's functions are internal and the family class is RiskMap_family", {
  exports <- getNamespaceExports("RiskMap")
  expect_false(any(c("maxim_integrand", "matern_gradient_phi",
                     "matern_hessian_phi", "create_ids") %in% exports))
  expect_false(exists("create_ids", envir = asNamespace("RiskMap"),
                      inherits = FALSE))
  expect_s3_class(negative_binomial(), "RiskMap_family", exact = TRUE)
})
