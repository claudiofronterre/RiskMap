# The parameter table (step 2a of the parameter table of the architecture
# review) moved the arithmetic of the sites that convert, bound, start and
# label the parameters without changing it, so each function is checked with
# identical() against the arithmetic it replaced, written out here as it was.

table_with_everything <- function() {
  parameter_table(c("(Intercept)", "x"), fix_tau2 = TRUE,
                  estimated_family = "theta", re_names = c("g", "h"),
                  distance_units = "km")
}

test_that("the table has a row for each parameter, in the order of the working vector", {
  parameters <- table_with_everything()
  rows <- parameters$rows
  expect_identical(rows$name, c("beta", "beta", "sigma2", "phi", "tau2",
                                "theta", "sigma2_re", "sigma2_re"))
  expect_identical(rows$stored, c("beta.(Intercept)", "beta.x", "sigma2",
                                  "phi", "nu2", "theta", "sigma2_re.g",
                                  "sigma2_re.h"))
  expect_identical(rows$reported[5], "tau2")
  expect_identical(rows$index, unname(unlist(parameters$layout$index)))
  expect_identical(rows$block, c("coefficients", "coefficients", "covariance",
                                 "covariance", "covariance", "family",
                                 "random_effects", "random_effects"))
  expect_identical(rows$label, c("(Intercept)", "x", "Spatial process var.",
                                 "Spatial corr. scale (km)",
                                 "Variance of the nugget",
                                 "Negative binomial size",
                                 "g (random eff. var.)",
                                 "h (random eff. var.)"))
  expect_identical(rows$working_label[5], "Nugget variance ratio")
  expect_identical(rows$start, c("regression", "regression", "one",
                                 "distance_quantile", "one",
                                 "family_regression", "one", "one"))
  expect_true(all(is.na(rows$held)))
})

test_that("held parameters have rows without a position", {
  parameters <- model_parameter_table(
    "(Intercept)", fix_tau2 = 0.5, exact = TRUE, fix_var_me = 0.1,
    estimated_family = character(), family_parameters = NULL,
    re_names = NULL, distance_units = "km"
  )
  held <- parameters$rows[is.na(parameters$rows$index), ]
  expect_identical(held$name, c("tau2", "sigma2_me"))
  expect_identical(held$held, c(0.5, 0.1))
  expect_identical(held$label, c("Variance of the nugget",
                                 "Measurement error var."))
  expect_identical(held_value(parameters, "tau2"), 0.5)
  expect_null(held_value(parameters, "sigma2"))
  # A model without a nugget holds it at zero.
  none <- parameter_table("(Intercept)", fix_tau2 = 0)
  expect_identical(held_value(none, "tau2"), 0)
  expect_false(nugget_estimated(none))
  # An estimated measurement error variance is in the family block.
  estimated <- model_parameter_table(
    "(Intercept)", fix_tau2 = TRUE, exact = TRUE, fix_var_me = NULL,
    estimated_family = character(), family_parameters = NULL,
    re_names = "g", distance_units = "km"
  )
  expect_identical(estimated$layout$family, "sigma2_me")
  expect_identical(estimate_names(estimated),
                   c("beta", "sigma2", "phi", "nu2", "sigma2_me",
                     "sigma2_re"))
})

test_that("working_start() repeats the arithmetic of the MCML start", {
  start <- list(beta = c(-0.3, 0.2), sigma2 = 0.7, phi = 0.13, tau2 = 0.31,
                theta = 2.4, sigma2_re = c(0.4, 1.7))
  old <- c(start$beta, log(start$sigma2), log(start$phi),
           log(start$tau2 / start$sigma2),
           family_parameter("theta")$working(start$theta),
           log(start$sigma2_re))
  expect_identical(working_start(table_with_everything(), start), old)
  # Without a nugget or random effects, and with a held nugget.
  plain <- parameter_table(c("(Intercept)", "x"), fix_tau2 = 0.2)
  expect_identical(working_start(plain, start),
                   c(start$beta, log(start$sigma2), log(start$phi)))
})

test_that("working_start() repeats the arithmetic of the direct Gaussian start", {
  parameters <- model_parameter_table(
    c("(Intercept)", "x"), fix_tau2 = TRUE, exact = TRUE, fix_var_me = NULL,
    estimated_family = character(), family_parameters = NULL,
    re_names = c("g", "h"), distance_units = "km"
  )
  start <- list(beta = c(1, 0.5), sigma2 = 1.9, phi = 0.21, tau2 = 0.37,
                sigma2_me = 0.11, sigma2_re = c(0.6, 2.2))
  covariance_start <- c(start$sigma2, start$phi, start$tau2, start$sigma2_me,
                        start$sigma2_re)
  old <- c(start$beta, log(covariance_start[1:2]),
           log(covariance_start[3] / covariance_start[1L]),
           log(covariance_start[4]), log(covariance_start[4 + 1:2]))
  expect_identical(working_start(parameters, start), old)
})

test_that("woodbury_start() repeats the start of the Woodbury Gaussian path", {
  parameters <- model_parameter_table(
    c("(Intercept)", "x"), fix_tau2 = TRUE, exact = TRUE, fix_var_me = NULL,
    estimated_family = character(), family_parameters = NULL,
    re_names = c("g", "h"), distance_units = "km"
  )
  start <- list(beta = c(1, 0.5), sigma2 = 1.9, phi = 0.21, tau2 = 0.37,
                sigma2_me = 0.11, sigma2_re = c(0.6, 2.2))
  # pack_cov_pars() followed by the division of every slot after the scale.
  start_cov_pars <- c(start$sigma2, start$phi, start$tau2, start$sigma2_me,
                      start$sigma2_re)
  start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)] / start_cov_pars[1]
  expect_identical(woodbury_start(parameters, start),
                   c(start$beta, log(start_cov_pars)))
  # It differs from working_start() only in the measurement error and the
  # random-effect variances, which step 2b of the table settles.
  differs <- woodbury_start(parameters, start) !=
    working_start(parameters, start)
  expect_identical(which(differs),
                   unlist(parameters$layout$index[c("sigma2_me",
                                                    "sigma2_re")],
                          use.names = FALSE))
})

test_that("working_estimate() splits the working vector as the engines did", {
  parameters <- table_with_everything()
  par <- c(0.1, -0.2, log(0.8), log(0.3), -1.1, -0.7, log(0.2), log(1.3))
  # structure_estimate().
  old <- list(beta = c("(Intercept)" = 0.1, x = -0.2),
              sigma2 = log(0.8), phi = log(0.3), nu2 = -1.1, theta = -0.7,
              sigma2_re = c(g = log(0.2), h = log(1.3)))
  expect_identical(working_estimate(parameters, par), old)
  expect_identical(names(unlist(old)), parameters$rows$stored)
})

test_that("natural_values() gives what coef() gave, in coef()'s order", {
  old_coef <- function(object) {
    estimate <- object$estimate
    res <- list()
    res$beta <- estimate$beta
    names(res$beta) <- colnames(as.matrix(object$D))
    res$sigma2 <- exp(estimate$sigma2)
    res$phi <- exp(estimate$phi)
    for (name in estimated_family_parameters(object)) {
      res[[name]] <- fit_traits(object)$parameters[[name]]$natural(
        estimate[[name]]
      )
    }
    if (!is.null(estimate$nu2)) {
      res$tau2 <- exp(estimate$nu2 + estimate$sigma2)
    }
    if (!is.null(estimate$sigma2_re)) res$sigma2_re <- exp(estimate$sigma2_re)
    res
  }
  for (fit in list(gaussian_model, gaussian_offset_model,
                   gaussian_intercept_model, binomial_model, poisson_model)) {
    expect_identical(coef(fit), old_coef(fit))
  }
  # A fit with a nugget, the family's parameter and two random effects.
  fit <- list(D = cbind("(Intercept)" = 1, x = 1:3), family = "negative_binomial",
              fix_tau2 = TRUE,
              estimate = list(beta = c("(Intercept)" = 0.1, x = -0.2),
                              sigma2 = log(0.8), phi = log(0.3), nu2 = -1.1,
                              theta = -0.7,
                              sigma2_re = c(g = log(0.2), h = log(1.3))))
  class(fit) <- "RiskMap"
  expect_identical(coef(fit), old_coef(fit))
  expect_named(coef(fit), c("beta", "sigma2", "phi", "theta", "tau2",
                            "sigma2_re"))
})

test_that("natural_values() with held values repeats the MCML reference point", {
  # estimate_to_mcml_reference().
  old_reference <- function(estimate, fix_tau2) {
    out <- list(beta = unname(estimate$beta), sigma2 = exp(estimate$sigma2),
                phi = exp(estimate$phi))
    out$tau2 <- if (isTRUE(fix_tau2)) {
      exp(estimate$sigma2 + estimate$nu2)
    } else {
      fix_tau2
    }
    if (!is.null(estimate$sigma2_re)) out$sigma2_re <- exp(estimate$sigma2_re)
    out
  }
  estimate <- list(beta = c(a = 0.1), sigma2 = log(0.8), phi = log(0.3),
                   nu2 = -1.1, sigma2_re = c(g = log(0.2)))
  for (fix_tau2 in list(TRUE, 0, 0.4)) {
    if (!isTRUE(fix_tau2)) estimate$nu2 <- NULL
    parameters <- parameter_table("a", fix_tau2 = fix_tau2, re_names = "g")
    new <- natural_values(parameters, estimate, held = TRUE)
    old <- old_reference(estimate, fix_tau2)
    for (name in names(old)) {
      expect_identical(unname(new[[name]]), unname(old[[name]]))
    }
  }
})

test_that("the nugget's maps repeat the conversions they replaced", {
  sigma2 <- 0.83
  tau2 <- 0.29
  expect_identical(nugget_working(tau2, sigma2), log(tau2 / sigma2))
  expect_identical(nugget_natural(-1.3, log(sigma2)),
                   exp(-1.3 + log(sigma2)))
  # Addition is commutative in floating point, so the order in which
  # estimate_to_mcml_reference() added the two did not matter.
  expect_identical(nugget_natural(-1.3, log(sigma2)),
                   exp(log(sigma2) + -1.3))
  expect_identical(nugget_ratio_natural(tau2, sigma2), tau2 / sigma2)
  fit <- list(fix_tau2 = TRUE,
              estimate = list(sigma2 = log(sigma2), nu2 = log(tau2 / sigma2)))
  expect_identical(fit_nugget_variance(fit),
                   exp(fit$estimate$nu2 + fit$estimate$sigma2))
  expect_identical(fit_nugget_variance(list(fix_tau2 = 0.2)), 0.2)
  expect_null(fit_nugget_variance(list(fix_tau2 = FALSE, estimate = list())))
})

test_that("summary()'s Jacobian is the one it built", {
  parameters <- table_with_everything()
  nm <- parameters$rows$stored
  J <- diag(length(nm))
  dimnames(J) <- list(nm, nm)
  J["nu2", "sigma2"] <- 1
  expect_identical(parameter_jacobian(parameters), J)
  plain <- parameter_table("(Intercept)", fix_tau2 = 0)
  expect_identical(unname(parameter_jacobian(plain)), diag(3))
})

test_that("parameter_bounds() repeats the literals of the fit and the profiles", {
  parameters <- table_with_everything()
  p <- 2L
  others <- nrow(estimated_rows(parameters)) - p - 2L
  log_scale <- log(37.2)
  # glgpm_nong_adaptive() and the engine's scale of make_profile_problem().
  expect_identical(parameter_bounds(parameters), list(
    lower = c(rep(-30, p), -20, -20, rep(-20, others)),
    upper = c(rep(30, p), 20, 10, rep(20, others))
  ))
  # make_profile_problem() on the scale of the model's units.
  expect_identical(parameter_bounds(parameters, log_scale), list(
    lower = c(rep(-30, p), -20, -20 + log_scale, rep(-20, others)),
    upper = c(rep(30, p), 20, 10 + log_scale, rep(20, others))
  ))
  # make_gaussian_profile_problem(), which solves the coefficients.
  expect_identical(
    parameter_bounds(parameters, log_scale, coefficients = c(-Inf, Inf)),
    list(lower = c(rep(-Inf, p), -20, -20 + log_scale, rep(-20, others)),
         upper = c(rep(Inf, p), 20, 10 + log_scale, rep(20, others)))
  )
})

test_that("variance_index() and limit_index() give the positions the engine read", {
  parameters <- table_with_everything()
  layout <- parameters$layout
  expect_identical(variance_index(parameters),
                   unlist(layout$index[c("sigma2", "nu2", "sigma2_re")],
                          use.names = FALSE))
  expect_identical(limit_index(parameters),
                   unlist(layout$index["theta"], use.names = FALSE))
  plain <- parameter_table("(Intercept)", fix_tau2 = 0)
  expect_identical(variance_index(plain),
                   unlist(plain$layout$index[c("sigma2", "nu2", "sigma2_re")],
                          use.names = FALSE))
  expect_identical(limit_index(plain),
                   unlist(plain$layout$index[character()], use.names = FALSE))
})

test_that("the labels are those of parameter_labels() and profile_labels()", {
  parameters <- table_with_everything()
  reported <- names(unlist(natural_values(parameters, working_estimate(
    parameters, numeric(8)
  ))))
  expect_identical(table_labels(parameters, reported),
                   parameter_labels(reported, "km"))
  old_profile <- parameter_labels(parameters$rows$stored, "km")
  old_profile[parameters$rows$stored == "nu2"] <- "Nugget variance ratio"
  expect_identical(profile_labels(parameters, parameters$rows$stored),
                   as.list(old_profile))
  expect_identical(table_labels(fit_parameter_table(binomial_model),
                                names(unlist(coef(binomial_model)))),
                   parameter_labels(names(unlist(coef(binomial_model))),
                                    binomial_model$distance_units))
})

test_that("the blocks of flat names repeat the tests on names", {
  names <- c("beta.(Intercept)", "beta.sigma2", "beta.theta", "sigma2", "phi",
             "tau2", "nu2", "sigma2_me", "theta", "sigma2_re.g")
  expect_identical(bootstrap_positive(names), !grepl("^beta\\.", names))
  family <- names %in% family_parameter_names()
  block <- parameter_block(names)
  expect_identical(block == "family", family)
  expect_identical(block %in% c("covariance", "random_effects"),
                   !grepl("^beta\\.", names) & !family)
})

test_that("the default starts are those glgpm() used", {
  parameters <- table_with_everything()
  coords <- cbind(c(0, 1, 3, 7), c(0, 2, 1, 5))
  expect_identical(default_start(parameters, "sigma2"), 1)
  expect_identical(default_start(parameters, "tau2"), 1)
  expect_identical(default_start(parameters, "sigma2_re"), rep(1, 2))
  expect_identical(default_start(parameters, "phi", coords),
                   quantile(pairwise_distances(coords), 0.1))
  gaussian <- model_parameter_table(
    "(Intercept)", fix_tau2 = 0, exact = TRUE, fix_var_me = NULL,
    estimated_family = character(), family_parameters = NULL,
    re_names = NULL, distance_units = "km"
  )
  expect_identical(default_start(gaussian, "sigma2_me"), 1)
  expect_error(default_start(parameters, "beta"),
               class = "riskmap_error_internal")
})

test_that("a fit's table is the one glgpm() built, with the engines' layout", {
  # The table glgpm() builds from the model's settings, before any estimate
  # exists, and the one rebuilt from the fit must agree, since the bootstrap
  # refits with the rebuilt table and an engine that builds its own layout
  # from the same settings.
  for (fit in list(gaussian_model, gaussian_offset_model,
                   gaussian_intercept_model, binomial_model, poisson_model)) {
    exact <- isTRUE(fit_traits(fit)$exact)
    built <- model_parameter_table(
      colnames(fit$D), fix_tau2 = fit$fix_tau2, exact = exact,
      fix_var_me = fit$fix_var_me, estimated_family = character(),
      family_parameters = fit$family_parameters, re_names = names(fit$re),
      distance_units = fit$distance_units
    )
    rebuilt <- fit_parameter_table(fit)
    expect_identical(rebuilt, built)
    engine <- cov_par_layout(
      ncol(fit$D), fix_tau2 = fit$fix_tau2,
      family_parameters = if (exact) gaussian_family_block(fit$fix_var_me) else
        estimated_family_parameters(fit),
      n_re = length(fit$re)
    )
    expect_identical(rebuilt$layout[c("index", "length")],
                     engine[c("index", "length")])
    expect_identical(estimated_rows(rebuilt)$stored,
                     names(unlist(fit$estimate)))
    expect_identical(working_estimate(rebuilt, unname(unlist(fit$estimate))),
                     fit$estimate)
  }
  # A model without a nugget has none in its table, whatever it estimates.
  expect_identical(gaussian_model$fix_tau2, 0)
  expect_false(nugget_estimated(fit_parameter_table(gaussian_model)))
  expect_null(gaussian_model$estimate$nu2)
})

test_that("a start that is missing or of the wrong length stops", {
  parameters <- model_parameter_table(
    c("(Intercept)", "x"), fix_tau2 = TRUE, exact = TRUE, fix_var_me = NULL,
    estimated_family = character(), family_parameters = NULL,
    re_names = c("g", "h"), distance_units = "km"
  )
  start <- list(beta = c(1, 0.5), sigma2 = 1.9, phi = 0.21, tau2 = 0.37,
                sigma2_me = 0.11, sigma2_re = c(0.6, 2.2))
  for (make in list(working_start, woodbury_start)) {
    expect_error(make(parameters, start[names(start) != "tau2"]),
                 class = "riskmap_error_internal")
    expect_error(make(parameters, utils::modifyList(start, list(
      sigma2_re = 0.6
    ))), class = "riskmap_error_internal")
    expect_error(make(parameters, utils::modifyList(start, list(beta = 1))),
                 class = "riskmap_error_internal")
  }
})

test_that("the Jacobian places the nugget by name, whatever the order", {
  parameters <- table_with_everything()
  shuffled <- parameters$rows$stored[c(1, 2, 3, 4, 6, 5, 7, 8)]
  J <- parameter_jacobian(parameters, shuffled)
  expect_identical(dimnames(J), list(shuffled, shuffled))
  expect_identical(J["nu2", "sigma2"], 1)
  expect_identical(sum(J), 9)
})
