test_that("glgpm produces errors", {
  skip_on_cran()
  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(beta = c(0, 1), unknown = 1)),
    "start_pars"
  )

  expect_input_error(
    glgpm("not formula", data = gaussian_data, family = "gaussian"),
    "formula"
  )

  expect_error(
    glgpm(~ cov + gp(), data = gaussian_data, family = "gaussian", messages = FALSE),
    "must have a response"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data, family = "gaussian", messages = FALSE),
    "formula"
  )

  expect_error(
    glgpm(y ~ cov + gp() + re(), data = gaussian_data, family = "gaussian", messages = FALSE),
    "needs at least one variable"
  )

  saturated_data <- gaussian_data
  saturated_data$row_id <- seq_len(nrow(saturated_data))

  expect_input_error(
    glgpm(y ~ cov + gp() + re(row_id), data = saturated_data, family = "gaussian", messages = FALSE),
    "formula"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = data, family = "gaussian"),
    "data"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "invalid"),
    "family"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", invlink = function(x) x),
    "invlink"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = poisson_data, family = "poisson", invlink = "not func", messages = FALSE),
    "invlink"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", denominator = 1),
    "denominator"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial", denominator = 1),
    "denominator"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial", denominator = not_present),
    "denominator"
  )

  # Records that aren't of the family.
  wrong <- list(
    binomial = list(y = c(1.5, -1), denominator = c(2.5, 0, NA, Inf, -1)),
    poisson = list(y = c(1.5, -1), denominator = c(0, NA, Inf, -1))
  )
  for (family in names(wrong)) {
    for (column in c("y", "denominator")) {
      for (value in wrong[[family]][[column]]) {
        wrong_data <- if (family == "binomial") binomial_data else poisson_data
        wrong_data$y[1] <- 0
        wrong_data[[column]][1] <- value
        expect_input_error(
          glgpm(y ~ cov + gp(), data = wrong_data, family = family, denominator = denominator, messages = FALSE),
          if (column == "y") "data" else column
        )
      }
    }
  }

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", distance_crs = 12345),
    "distance_crs", class = "riskmap_error_input"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", distance_units = "not valid"),
    "distance_units"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = latlon_data, family = "gaussian", distance_crs = 4326, messages = FALSE),
    "distance_crs"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", return_samples = "not logical"),
    "return_samples"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", messages = "not logical"),
    "messages"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial",
          denominator = denominator, fixed_pars = list(sigma2_me = 1)),
    "fixed_pars"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          fixed_pars = list(sigma2_me = c(1, 2))),
    "fixed_pars"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          fixed_pars = list(sigma2_me = "not number")),
    "fixed_pars"
  )

  expect_input_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          fixed_pars = list(sigma2_me = -1)),
    "fixed_pars"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(invalid = 1), messages = FALSE),
    "start_pars", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(invalid = 1, silly = 2), messages = FALSE),
    "\"invalid\" and \"silly\" are not among them", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(beta = 1), messages = FALSE),
    "must have a value for each column of the design matrix, 2, not 1",
    fixed = TRUE
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          control_mcml = list(max_stages = 2), messages = FALSE),
    "set_control_mcml", class = "riskmap_error_input"
  )

  expect_error(set_control_mcml(initial_draws = 0), "positive integer")
  expect_error(set_control_mcml(max_stages = 0), "positive integer")
  expect_error(set_control_mcml(tolerance = -1), "positive number")
  expect_error(
    set_control_mcml(initial_draws = 1000, max_draws = 500),
    "at least `initial_draws`", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(beta = c("a", "b")), messages = FALSE),
    "must be finite numbers", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2 = -1), messages = FALSE),
    "sigma2", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2 = "a"), messages = FALSE),
    "sigma2", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(phi = -1), messages = FALSE),
    "phi", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data, family = "gaussian",
          fixed_pars = list(sigma2_me = 0), start_pars = list(tau2 = -1), messages = FALSE),
    "tau2", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(tau2 = 1), messages = FALSE),
    "tau2", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2_me = -1), messages = FALSE),
    "sigma2_me", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2_re = c(1, 2)), messages = FALSE),
    "\"sigma2_re\" is not among them", class = "riskmap_error_input"
  )

  expect_error(
    glgpm(y ~ cov + gp() + re(i), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2_re = c(1, 2)), messages = FALSE),
    "a value for each `re()` term", fixed = TRUE, class = "riskmap_error_input"
  )

})

expected_output <- c("estimate", "grad_MLE", "covariance", "log_lik",
                     "y", "D", "coords", "ID_coords", "re", "ID_re", "nugget",
                     "fixed_pars", "formula", "family", "distance_units",
                     "data", "input_crs", "kappa", "units_m", "cov_offset", "call",
                     "S_samples", "link_function", "mcml_history",
                     "mcml_certified", "mcml_selected_iteration",
                     "linear_predictors", "fitted_values", "seed")

test_that("glgpm produces expected output for gaussian models", {

  fit_no_re <- glgpm(y ~ cov + gp(),
                     data = gaussian_data,
                     family = "gaussian",
                     distance_units = "m",
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_null(fit_no_re$mcml_history)
  expect_length(fitted(fit_no_re), nrow(gaussian_data))
  expect_identical(fitted(fit_no_re), fit_no_re$fitted_values)
  expect_identical(fitted(fit_no_re, type = "link"),
               fit_no_re$linear_predictors)
  expect_null(fit_no_re$mcml_certified)
  expect_null(fit_no_re$mcml_selected_iteration)
  expect_identical(fit_no_re$family, "gaussian")
  expect_identical(fit_no_re$coords[,1], data$x)
  expect_identical(fit_no_re$coords[,2], data$z)
  expect_identical(fit_no_re$y, gaussian_data$y)
  expect_identical(unname(fit_no_re$D[,2]), gaussian_data$cov)
  optimizer <- attr(fit_no_re, "optimizer")
  expect_named(
    optimizer,
    c("convergence", "message", "evaluations", "invalid_evaluations",
      "max_abs_gradient", "stationary", "information_rcond",
      "information_jitter")
  )
  expect_type(optimizer$stationary, "logical")
  expect_gte(optimizer$information_rcond, 0)

  fit_re_fixed_me <- glgpm(y ~ cov + gp() + re(i),
                           data = gaussian_data,
                           family = "gaussian",
                           fixed_pars = list(sigma2_me = 0.1),
                           messages = FALSE)
  expect_true(is.finite(coef(fit_re_fixed_me)$sigma2_re[["i"]]))

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = gaussian_data,
                  family = "gaussian",
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_identical(fit_re$family, "gaussian")
  expect_length(fit_re$re, 1)

  fit_re2 <- glgpm(y ~ cov + gp() + re(i, j),
                  data = gaussian_data,
                  family = "gaussian",
                  messages = FALSE)

  expect_s3_class(fit_re2, "RiskMap")
  expect_length(fit_re2$re, 2)

  fit_kappa <- glgpm(y ~ cov + gp(kappa = 0.4),
                     data = gaussian_data,
                     family = "gaussian",
                     distance_units = "m",
                     messages = FALSE)

  expect_s3_class(fit_kappa, "RiskMap")
  expect_false(identical(coef(fit_no_re)[["sigma2"]], coef(fit_kappa)[["sigma2"]]))
  expect_false(identical(coef(fit_no_re)[["phi"]], coef(fit_kappa)[["phi"]]))
})

test_that("glgpm fits replicated gaussian models with the Woodbury path #165", {
  skip_on_cran()
  # replicating locations gives more than three observations per latent
  # effect, so glgpm_lm should not fall back to glgpm_lm_direct; six copies
  # are needed for re(i, j), which has 20 + 10 + 2 = 32 latent effects
  set.seed(2)
  replicated_data <- do.call(rbind, replicate(6, gaussian_data, simplify = FALSE))
  replicated_data$y <- replicated_data$y + rnorm(nrow(replicated_data), sd = 0.3)

  # re(i, j) exercises the cross-terms between random effects in the Hessian
  formulas <- list(y ~ cov + gp(), y ~ cov + gp() + re(i),
                   y ~ cov + gp() + re(i, j))
  for (formula in formulas) {
    woodbury_fit <- with_mocked_bindings(
      glgpm(formula, data = replicated_data, family = "gaussian",
            messages = FALSE),
      glgpm_lm_direct = function(...) stop("direct path used")
    )
    direct_fit <- with_mocked_bindings(
      glgpm(formula, data = replicated_data, family = "gaussian",
            messages = FALSE),
      use_direct_gaussian_covariance = function(...) TRUE
    )

    expect_s3_class(woodbury_fit, "RiskMap")
    expect_setequal(names(woodbury_fit), expected_output)
    # near-zero sigma2 for re(i, j) makes tiny absolute differences large relative ones
    expect_equal(woodbury_fit$estimate, direct_fit$estimate, tolerance = 1e-4)
    expect_equal(woodbury_fit$log_lik, direct_fit$log_lik, tolerance = 1e-4)
    expect_equal(woodbury_fit$covariance, direct_fit$covariance,
                 tolerance = 1e-4)
  }
})

test_that("glgpm estimates a nugget in replicated gaussian models with the Woodbury path", {
  # replicated locations separate the nugget from the measurement error, and
  # give enough observations per latent effect to avoid glgpm_lm_direct. The
  # spatial correlation reaches several locations, so that it separates the
  # nugget from the spatial variance and their covariance is well determined.
  set.seed(4)
  sites <- cbind(runif(40, 0, 10000), runif(40, 0, 10000))
  spatial <- crossprod(chol(2 * exp(-as.matrix(dist(sites)) / 3000)),
                       rnorm(40))
  site <- rep(seq_len(40), 4)
  replicated_data <- data.frame(x = sites[site, 1], z = sites[site, 2],
                                cov = rnorm(160))
  replicated_data$y <- 1 + 0.5 * replicated_data$cov + spatial[site] +
    rnorm(40, sd = 0.7)[site] + rnorm(160, sd = 0.3)
  replicated_data <- st_as_sf(replicated_data, coords = c("x", "z"),
                              crs = 32637)

  woodbury_fit <- with_mocked_bindings(
    glgpm(y ~ cov + gp(nugget = TRUE), data = replicated_data,
          family = "gaussian", messages = FALSE),
    glgpm_lm_direct = function(...) stop("direct path used")
  )
  direct_fit <- with_mocked_bindings(
    glgpm(y ~ cov + gp(nugget = TRUE), data = replicated_data,
          family = "gaussian", messages = FALSE),
    use_direct_gaussian_covariance = function(...) TRUE
  )

  expect_s3_class(woodbury_fit, "RiskMap")
  expect_setequal(names(woodbury_fit$estimate),
                  c("beta", "sigma2", "phi", "nu2", "sigma2_me"))
  expect_true(all(is.finite(unlist(woodbury_fit$estimate))))
  expect_equal(woodbury_fit$estimate, direct_fit$estimate, tolerance = 1e-6)
  expect_equal(woodbury_fit$log_lik, direct_fit$log_lik, tolerance = 1e-6)
  # the analytical Hessian (including its nu2 rows/columns) gives the covariance
  expect_equal(woodbury_fit$covariance, direct_fit$covariance,
               tolerance = 1e-4)
})

test_that("glgpm fits gaussian models combining re() with a fixed measurement error variance #117", {

  fit <- glgpm(y ~ cov + gp() + re(i),
              data = gaussian_data,
              family = "gaussian",
              fixed_pars = list(sigma2_me = 0.1),
              messages = FALSE)

  expect_s3_class(fit, "RiskMap")
  expect_setequal(names(fit), expected_output)
  expect_true(all(is.finite(unlist(fit$estimate))))

  fit_nugget <- glgpm(y ~ cov + gp(nugget = TRUE) + re(i),
                      data = gaussian_data,
                      family = "gaussian",
                      fixed_pars = list(sigma2_me = 0.1),
                      messages = FALSE)

  expect_s3_class(fit_nugget, "RiskMap")
  expect_true(all(is.finite(unlist(fit_nugget$estimate))))
})


test_that("glgpm fits gaussian models with saturated RE when fix_var_me is provided", {
  saturated_data <- gaussian_data
  saturated_data$row_id <- seq_len(nrow(saturated_data))

  fit <- glgpm(formula = y ~ cov + gp() + re(row_id), data = saturated_data,
               family = "gaussian", messages = FALSE, fixed_pars = list(sigma2_me = 0.1))

  expect_s3_class(fit, "RiskMap")
  expect_setequal(names(fit$estimate), c("beta", "sigma2", "phi", "sigma2_re"))

  expect_error(
    glgpm(formula = y ~ cov + gp(nugget = TRUE) + re(row_id),
          data = saturated_data,
          family = "gaussian",
          messages = FALSE,
          fixed_pars = list(sigma2_me = 0.1)),
    "told apart from the nugget", class = "riskmap_error_input"
  )
})

test_that("check_computation() flags a spatial scale shorter than the spacing of the locations (#206)", {
  # true range (phi) is far shorter than the spacing between locations
  set.seed(42)
  n_pts <- 20L
  coords <- cbind(runif(n_pts, 0, 10000), runif(n_pts, 0, 10000))
  cov <- rnorm(n_pts)
  sigma2 <- 1; phi <- 50; kappa <- 1
  Sigma <- sigma2 * matern_correlation(dist(coords), phi = phi, kappa = kappa,
                                       return_sym_matrix = TRUE)
  S <- as.numeric(t(chol(Sigma)) %*% rnorm(n_pts))
  denom <- sample(10:20, n_pts, replace = TRUE)
  y <- rbinom(n_pts, denom, plogis(0.2 + 0.3 * cov + S))
  no_spatial_data <- st_as_sf(
    data.frame(x = coords[, 1], z = coords[, 2], cov = cov, y = y,
               denominator = denom),
    coords = c("x", "z"), crs = 32637
  )

  fit <- glgpm(y ~ cov + gp(), data = no_spatial_data, family = "binomial",
               denominator = denominator, control_mcmc = control_mcmc,
               seed = fixture_seed, control_mcml = fixture_control_mcml,
               messages = FALSE)
  expect_true(attr(fit, "optimizer")$mcml_small_scale)
  checks <- check_computation(fit)$checks
  expect_identical(checks$status[checks$check == "Spatial scale"], "SUSPICIOUS")
})

test_that("a random-effect variance is held at zero when the data have no such effect (#206)", {
  # real, well-resolved spatial effect but no random-effect contribution at all
  set.seed(42)
  n_pts <- 20L
  coords <- cbind(runif(n_pts, 0, 10000), runif(n_pts, 0, 10000))
  i <- rep(seq_len(n_pts / 2L), each = 2)
  cov <- rnorm(n_pts)
  sigma2 <- 1; phi <- 5500; kappa <- 1
  Sigma <- sigma2 * matern_correlation(dist(coords), phi = phi, kappa = kappa,
                                       return_sym_matrix = TRUE)
  S <- as.numeric(t(chol(Sigma)) %*% rnorm(n_pts))
  denom <- sample(10:20, n_pts, replace = TRUE)
  y <- rbinom(n_pts, denom, plogis(0.2 + 0.3 * cov + S))
  no_re_data <- st_as_sf(
    data.frame(x = coords[, 1], z = coords[, 2], cov = cov, y = y,
               denominator = denom, i = i),
    coords = c("x", "z"), crs = 32637
  )

  fit <- glgpm(y ~ cov + gp() + re(i), data = no_re_data, family = "binomial",
               denominator = denominator, control_mcmc = control_mcmc,
               seed = fixture_seed, control_mcml = fixture_control_mcml,
               messages = FALSE)
  expect_identical(attr(fit, "optimizer")$mcml_boundary,
                   "random-effect variance for 'i'")
  variances <- check_computation(fit)$results$Variances$row
  expect_identical(variances$status, "OK")
  expect_match(variances$brief, "Zero is the maximum", fixed = TRUE)
})

test_that("glgpm warns when the observed information matrix is singular (#206)", {
  # three locations can't identify both a spatial variance and a nugget, as in
  # test-nugget.R, whose locations give a singular information on every platform
  three_locations <- st_as_sf(
    data.frame(x = c(1, 2, 3), z = c(0, 1, 2), y = 0:2, denominator = 4),
    coords = c("x", "z"), crs = 4326
  )

  expect_warnings(
    fit <- glgpm(y ~ gp(nugget = TRUE), data = three_locations,
                 family = "binomial", denominator = denominator,
                 distance_crs = 32631, control_mcmc = control_mcmc,
                 seed = fixture_seed, control_mcml = fixture_control_mcml,
                 messages = FALSE),
    c("`sigma2` was estimated as effectively zero",
      "observed information matrix is singular")
  )
  expect_true(attr(fit, "optimizer")$mcml_identifiability_failure)
})

test_that("glgpm is invariant to fixed-effect covariate scale", {
  scaled_data <- gaussian_data
  scaled_data$cov_large <- scaled_data$cov * 1e6

  ordinary_fit <- glgpm(
    y ~ cov + gp(),
    data = scaled_data,
    family = "gaussian",
    messages = FALSE
  )
  scaled_fit <- glgpm(
    y ~ cov_large + gp(),
    data = scaled_data,
    family = "gaussian",
    messages = FALSE
  )

  expect_equal(ordinary_fit$log_lik, scaled_fit$log_lik,
               tolerance = 1e-7)
  expect_equal(
    as.numeric(ordinary_fit$D %*% ordinary_fit$estimate$beta),
    as.numeric(scaled_fit$D %*% scaled_fit$estimate$beta),
    tolerance = 1e-7
  )
  expect_equal(
    unname(ordinary_fit$estimate$beta[["cov"]]),
    unname(scaled_fit$estimate$beta[["cov_large"]]) * 1e6,
    tolerance = 1e-7
  )
  expect_equal(
    unname(diag(ordinary_fit$covariance)[1:2]),
    unname(diag(scaled_fit$covariance)[1:2]) * c(1, 1e12),
    tolerance = 1e-6
  )
})

test_that("glgpm produces expected output for binomial models", {

  fit_no_re <- glgpm(y ~ cov + gp(),
                     data = binomial_data,
                     family = "binomial",
                     denominator = denominator,
                     control_mcmc = control_mcmc,
                     seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_identical(fit_no_re$family, "binomial")
  expect_length(fitted(fit_no_re), nrow(binomial_data))
  expect_true(all(fitted(fit_no_re) >= 0 & fitted(fit_no_re) <= 1))
  expect_length(fit_no_re$linear_predictors, nrow(binomial_data))
  expect_false(is.null(fit_no_re$S_samples))
  starting_values <- attr(fit_no_re, "starting_values")
  expect_true(starting_values$selected %in%
                c("current", "transformed_gaussian"))
  expect_true(is.finite(starting_values$current_laplace))
  expect_true(is.finite(starting_values$transformed_laplace))

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = binomial_data,
                  family = "binomial",
                  denominator = denominator,
                  control_mcmc = control_mcmc, seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_identical(fit_re$family, "binomial")
  eta_samples <- sweep(
    fit_re$S_samples[, fit_re$ID_coords, drop = FALSE], 2,
    as.numeric(fit_re$D %*% fit_re$estimate$beta + fit_re$cov_offset), "+"
  )
  n_locations <- nrow(fit_re$coords)
  re_samples <- fit_re$S_samples[, n_locations + seq_along(fit_re$re[[1]]),
                                 drop = FALSE]
  eta_samples <- eta_samples +
    re_samples[, fit_re$ID_re[[1]], drop = FALSE]
  expect_equal(fitted(fit_re), colMeans(plogis(eta_samples)))
  expect_equal(fit_re$linear_predictors, colMeans(eta_samples))
  # A model with re() terms keeps the default start without a comparison.
  re_starting_values <- attr(fit_re, "starting_values")
  expect_identical(re_starting_values$selected, "current")
  expect_identical(re_starting_values$current_laplace, NA_real_)
  expect_identical(re_starting_values$transformed_laplace, NA_real_)
  expect_match(re_starting_values$fallback_reason, "re() terms",
               fixed = TRUE)
})

test_that("responses and denominators of counts are whole to floating point", {
  skip_on_cran()
  fit <- function(data, family) {
    glgpm(y ~ cov + gp(), data = data, family = family,
          denominator = denominator, seed = fixture_seed,
          fixed_pars = list(beta = c(0.2, 0.3), sigma2 = 0.5, phi = 2),
          control_mcmc = control_mcmc, messages = FALSE)
  }
  # As a number found by a division may be.
  nearly <- function(x) x * (1 + 1e-12)
  data <- binomial_data
  whole <- fit(data, "binomial")
  data$y <- nearly(data$y)
  data$denominator <- nearly(data$denominator)
  rounded <- fit(data, "binomial")
  expect_identical(rounded$y, whole$y)
  expect_identical(rounded$units_m, whole$units_m)
  expect_identical(rounded$S_samples, whole$S_samples)
  # The simulation takes the denominators of the fit, as whole numbers.
  expect_identical(simulated_values(simulate_glgpm(rounded, nsim = 2, seed = 1)),
                   simulated_values(simulate_glgpm(whole, nsim = 2, seed = 1)))

  # The start of a Poisson model regresses the rounded responses, which can be
  # just below zero, and the exposures aren't rounded.
  data <- poisson_data
  data$denominator <- data$denominator + 0.5
  data$y[1] <- 0
  whole <- fit(data, "poisson")
  data$y <- nearly(data$y)
  data$y[1] <- 0.3 - 0.1 * 3
  data$y[2] <- data$y[2] + 1.2e-8
  expect_lt(data$y[1], 0)
  rounded <- fit(data, "poisson")
  expect_identical(rounded$y, whole$y)
  expect_identical(rounded$units_m, data$denominator)
  expect_identical(rounded$S_samples, whole$S_samples)
})

test_that("a model with re() terms keeps the default start without fitting a Gaussian model", {
  skip_on_cran()
  # The transformed-Gaussian candidate of such a model was fitted with dense
  # matrices over all the records or all the levels of the terms, and the
  # comparison could never choose it. The calls are counted, since an error
  # in either would only be caught as a failed comparison.
  calls <- 0L
  local_mocked_bindings(
    glgpm_lm = function(...) {
      calls <<- calls + 1L
      stop("the Gaussian model was fitted")
    },
    starting_laplace_log_marginal = function(...) {
      calls <<- calls + 1L
      stop("the starting values were compared")
    }
  )
  fit_with <- function(family, start_pars = list()) {
    suppressWarnings(glgpm(
      y ~ cov + gp() + re(i),
      data = if (family == "binomial") binomial_data else poisson_data,
      family = family, denominator = denominator,
      control_mcmc = control_mcmc, seed = fixture_seed,
      control_mcml = fixture_control_mcml, start_pars = start_pars,
      return_samples = FALSE, messages = FALSE
    ))
  }
  for (family in c("binomial", "poisson")) {
    fit <- fit_with(family)
    expect_identical(calls, 0L)
    expect_identical(attr(fit, "starting_values")$selected, "current")
    expect_identical(attr(fit, "starting_values")$reason, "re_terms")
    # The default start, as when it is given.
    given <- fit_with(family, list(sigma2_re = 1))
    expect_identical(fit$estimate, given$estimate)
  }
})

test_that("the transformed-Gaussian start is fitted to the totals of each location", {
  # One record per location: the records as they are.
  totals <- location_totals(1:3, c(5, 6, 7), cbind(1, c(0.1, 0.2, 0.3)),
                            c(0, 0.5, 1), c(2L, 1L, 3L), c(2L, 1L, 3L))
  expect_identical(totals$y, 1:3)
  expect_identical(totals$ID_coords, c(2L, 1L, 3L))
  expect_identical(totals$s_unique, c(2L, 1L, 3L))
  # Several records at a location: sums of the responses and denominators,
  # and the design rows and offsets averaged with the denominators as
  # weights, in the order of the locations.
  design <- cbind("(Intercept)" = 1, x = c(1, 3, 2, 4))
  totals <- location_totals(c(1, 2, 0, 1), c(2, 6, 4, 4), design,
                            c(0, 1, 0, 0), c(2L, 2L, 1L, 1L), c(2L, 1L))
  expect_identical(totals$y, c(1, 3))
  expect_identical(totals$units_m, c(8, 8))
  expect_equal(totals$D, cbind("(Intercept)" = c(1, 1),
                               x = c((2 * 4 + 4 * 4) / 8, (1 * 2 + 3 * 6) / 8)))
  expect_equal(totals$offset, c(0, 6 / 8))
  expect_identical(totals$ID_coords, 1:2)
  expect_identical(totals$s_unique, 1:2)
  # A location whose denominators are all zero averages with equal weights.
  totals <- location_totals(c(1, 0, 0, 0), c(2, 6, 0, 0), design,
                            c(0, 1, 0, 1), c(2L, 2L, 1L, 1L), c(2L, 1L))
  expect_identical(totals$units_m, c(0, 8))
  expect_equal(totals$D[1, ], c("(Intercept)" = 1, x = 3))
  expect_equal(totals$offset, c(0.5, 6 / 8))

  # In a fit, the Gaussian model has one observation per location.
  data <- binomial_data[rep(seq_len(nrow(binomial_data)), each = 2), ]
  observations <- NULL
  local_mocked_bindings(glgpm_lm = function(y, ...) {
    observations <<- length(y)
    stop("recorded")
  })
  fit <- glgpm(y ~ cov + gp(), data = data, family = "binomial",
               denominator = denominator, control_mcmc = control_mcmc,
               seed = fixture_seed, control_mcml = fixture_control_mcml,
               return_samples = FALSE, messages = FALSE)
  expect_identical(observations, nrow(binomial_data))
  expect_match(attr(fit, "starting_values")$fallback_reason, "recorded")
})

test_that("the comparison of starts reads named Laplace values", {
  # With re() terms the transformed candidate's value carried the name of the
  # term, so scores["transformed_gaussian"] was missing and the defaults were
  # always kept.
  value <- 0
  local_mocked_bindings(starting_laplace_log_marginal = function(...) {
    value <<- value + 1
    c(hh = -10 + value)
  })
  fit <- glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial",
               denominator = denominator, control_mcmc = control_mcmc,
               seed = fixture_seed, control_mcml = fixture_control_mcml,
               return_samples = FALSE, messages = FALSE)
  diagnostics <- attr(fit, "starting_values")
  expect_identical(diagnostics$selected, "transformed_gaussian")
  expect_identical(diagnostics$current_laplace, -9)
  expect_identical(diagnostics$transformed_laplace, -8)
  expect_null(diagnostics$reason)
})

test_that("transformed starts preserve an estimated nugget", {
  skip_on_cran()
  fit <- glgpm(
    y ~ cov + gp(nugget = TRUE), data = binomial_data,
    family = "binomial", denominator = denominator,
    control_mcmc = control_mcmc, seed = fixture_seed,
    control_mcml = fixture_control_mcml, messages = FALSE
  )

  diagnostics <- attr(fit, "starting_values")
  expect_true(diagnostics$selected %in%
                c("current", "transformed_gaussian"))
  expect_true(is.finite(diagnostics$current_laplace))
  expect_true(is.finite(diagnostics$transformed_laplace))
  expect_true(is.finite(fit$estimate$nu2))
})

test_that("explicit non-Gaussian starts remain authoritative", {
  skip_on_cran()
  data_frame <- sf::st_drop_geometry(binomial_data)
  glm_start <- glm(
    cbind(y, denominator - y) ~ cov,
    data = data_frame,
    family = binomial
  )
  start_pars <- list(
    beta = coef(glm_start),
    sigma2 = 1,
    phi = quantile(
      pairwise_distances(coordinates_in_units(binomial_data, "km")),
      0.1
    )
  )

  explicit_fit_2 <- glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
    start_pars = start_pars, messages = FALSE
  )
  explicit_fit <- glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
    start_pars = start_pars, messages = FALSE
  )

  expect_identical(attr(explicit_fit, "starting_values")$selected, "user")
  expect_identical(explicit_fit$estimate, explicit_fit_2$estimate)
  expect_identical(explicit_fit$log_lik, explicit_fit_2$log_lik)
})

test_that("final fitted samples support a nugget model", {
  fit <- glgpm(
    y ~ cov + gp(nugget = TRUE),
    fixed_pars = list(tau2 = 0.1),
    data = binomial_data,
    family = "binomial",
    denominator = denominator,
    control_mcmc = control_mcmc, seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
    messages = FALSE
  )

  expect_identical(ncol(fit$S_samples), nrow(fit$coords))
  expect_length(fitted(fit), nrow(binomial_data))
  expect_true(all(is.finite(fitted(fit))))
})

test_that("custom inverse links retain the current automatic start", {
  skip_on_cran()
  custom_link <- list(
    inv = plogis,
    d1 = function(x) plogis(x) * (1 - plogis(x)),
    d2 = function(x) {
      probability <- plogis(x)
      probability * (1 - probability) * (1 - 2 * probability)
    }
  )
  fit <- glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
    invlink = custom_link,
    messages = FALSE
  )

  diagnostics <- attr(fit, "starting_values")
  expect_identical(diagnostics$selected, "current")
  expect_match(diagnostics$fallback_reason, "custom inverse links")
})

test_that("a binomial model without an intercept starts from a regression without one", {
  # The start regression dropped the first column of the design matrix and
  # added an intercept, so a model without one started from the coefficients
  # of another model, shifted by one column. With a named link that start
  # went straight to MCML.
  regression_start <- function(formula, family) {
    start <- NULL
    local_mocked_bindings(select_nongaussian_start = function(current, ...) {
      start <<- current$beta
      stop("start recorded")
    })
    expect_error(glgpm(formula, data = binomial_data, family = family,
                       denominator = denominator, messages = FALSE),
                 "start recorded")
    start
  }
  glm_start <- function(formula, link) {
    coef(glm(formula, family = binomial(link),
             data = st_drop_geometry(binomial_data)))
  }
  for (link in c("logit", "probit")) {
    expect_identical(
      regression_start(y ~ 0 + factor(j) + cov + gp(), binomial(link)),
      glm_start(cbind(y, denominator - y) ~ 0 + factor(j) + cov, link)
    )
    expect_identical(
      regression_start(y ~ 0 + cov + gp(), binomial(link)),
      glm_start(cbind(y, denominator - y) ~ 0 + cov, link)
    )
    expect_identical(
      regression_start(y ~ factor(j) + cov + gp(), binomial(link)),
      glm_start(cbind(y, denominator - y) ~ factor(j) + cov, link)
    )
  }
})

test_that("glgpm produces expected output for poisson models", {

  fit_no_re <- glgpm(y ~ cov + gp(),
                     data = poisson_data,
                     family = "poisson",
                     control_mcmc = control_mcmc, seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_identical(fit_no_re$family, "poisson")
  expect_length(fitted(fit_no_re), nrow(poisson_data))
  expect_true(all(fitted(fit_no_re) > 0))

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = poisson_data,
                  family = "poisson",
                  control_mcmc = control_mcmc, seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_identical(fit_re$family, "poisson")

  fit_re_den <- glgpm(y ~ cov + gp() + re(i),
                      data = poisson_data,
                      denominator = denominator,
                      family = "poisson",
                      control_mcmc = control_mcmc, seed = fixture_seed,
                     control_mcml = fixture_control_mcml,
                      messages = FALSE)

  expect_s3_class(fit_re_den, "RiskMap")
  expect_setequal(names(fit_re_den), expected_output)
  expect_identical(fit_re_den$family, "poisson")
  expect_true(all(fitted(fit_re_den) > 0))
})

test_that("certified MCML records its certificate and bank diagnostics", {
  set.seed(9917)
  caller_seed <- .Random.seed
  fit <- glgpm(
    y ~ cov + gp(), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
    control_mcml = fixture_control_mcml, messages = FALSE
  )
  stage <- fit$mcml_history[[1L]]
  diagnostics <- attr(fit, "optimizer")

  expect_identical(fit$mcml_selected_iteration, 1L)
  expect_length(fit$mcml_history, 1L)
  expect_identical(stage$budget, 500L)
  expect_true(all(c(
    "fit_diagnostics", "first_validation", "second_validation",
    "replicate_validation", "certificate", "certificate_distance",
    "certificate_radius", "certificate_radii", "radius_from",
    "correction_size", "certified", "outcome", "timing"
  ) %in% names(stage)))
  expect_identical(stage$certificate_radius, max(stage$certificate_radii))
  expect_equal(stage$certificate,
               stage$certificate_distance + stage$certificate_radius)
  expect_identical(names(stage$replicate_validation),
                   names(stage$second_validation))
  expect_null(stage$information)
  expect_identical(diagnostics$mcml_budget, 500L)
  expect_identical(diagnostics$mcml_validation_budget, 500L)
  expect_identical(diagnostics$mcml_stages, 1L)
  expect_identical(diagnostics$mcml_certified, fit$mcml_certified)
  expect_identical(diagnostics$mcml_certificate, stage$certificate)
  expect_identical(.Random.seed, caller_seed)
})

test_that("a fit is reproduced from the seed it stores", {
  skip_on_cran()
  fit_binomial <- function(seed = NULL) {
    suppressWarnings(glgpm(
      y ~ cov + gp(), data = binomial_data, family = "binomial",
      denominator = denominator, control_mcmc = control_mcmc,
      control_mcml = fixture_control_mcml, seed = seed, messages = FALSE
    ))
  }
  set.seed(31)
  first <- fit_binomial()
  second <- fit_binomial()
  expect_type(first$seed, "integer")
  expect_false(identical(first$seed, second$seed))

  refit <- fit_binomial(first$seed)
  expect_identical(refit$seed, first$seed)
  expect_identical(refit$estimate, first$estimate)
  expect_identical(refit$S_samples, first$S_samples)
  expect_identical(refit$fitted_values, first$fitted_values)

  # Without a seed, one integer is drawn from the session's stream, which
  # advances it like any random call, and nothing else is drawn from it.
  set.seed(32)
  expected_seed <- sample.int(.Machine$integer.max, 1L)
  expected_state <- .Random.seed
  set.seed(32)
  drawn <- fit_binomial()
  expect_identical(drawn$seed, expected_seed)
  expect_identical(.Random.seed, expected_state)
  # A supplied seed leaves the session's stream unchanged, and is stored.
  fit_binomial(5L)
  expect_identical(.Random.seed, expected_state)
  expect_identical(binomial_model$seed, as.integer(fixture_seed))
})

test_that("Gaussian fits store no seed and seeds are validated", {
  expect_true("seed" %in% names(gaussian_model))
  expect_null(gaussian_model$seed)
  seeded <- glgpm(y ~ gp(), data = gaussian_data, family = "gaussian",
                  seed = 3, messages = FALSE)
  expect_null(seeded$seed)
  expect_identical(seeded$estimate, gaussian_intercept_model$estimate)
  expect_error(
    glgpm(y ~ gp(), data = gaussian_data, family = "gaussian", seed = -1,
          messages = FALSE),
    "seed", class = "riskmap_error_input"
  )
  expect_error(set_control_mcmc(seed = 1), "unused argument")
})


test_that("glgpm correctly reprojects to new CRS", {

  latlon <- st_transform(gaussian_data, 4326)

  suggested_crs <- propose_crs(latlon, messages = FALSE)
  sf_reproj <- st_transform(latlon, suggested_crs)
  scaled_coords <- coordinates_in_units(sf_reproj, "km")

  expect_message(
    fit <- glgpm(y ~ cov + gp(),
                 data = latlon,
                 family = "gaussian",
                 distance_units = "km",
                 messages = FALSE),
    "projected to EPSG"
  )

  expect_identical(fit$input_crs, st_crs(latlon))
  expect_identical(fit$coords, scaled_coords)
  expect_identical(st_crs(fit$data), st_crs(suggested_crs))
})

test_that("glgpm honours an explicitly supplied projected distance_crs", {

  latlon <- st_transform(gaussian_data, 4326)

  fit <- glgpm(y ~ cov + gp(),
               data = latlon,
               family = "gaussian",
               distance_crs = 32637,
               distance_units = "m",
               messages = FALSE)

  expect_identical(st_crs(fit$data), st_crs(32637))
  expect_identical(fit$coords, coordinates_in_units(fit$data, "m"))
})

test_that("coordinates_in_units converts projected CRS units", {

  metre_data <- gaussian_data
  foot_data <- st_transform(gaussian_data, 2263)

  expect_identical(coordinates_in_units(metre_data, "km"),
               st_coordinates(metre_data) / 1000)
  expect_identical(coordinates_in_units(foot_data, "m"),
               st_coordinates(foot_data) * 0.3048006096012192)
})


test_that("glgpm fits several re() calls as one call with all their variables", {
  skip_on_cran()
  gaussian_start <- list(sigma2_re = c(5, 10), sigma2_me = 0.01)
  gaussian_one_call <- glgpm(y ~ cov + gp() + re(i, j), data = gaussian_data,
                             family = "gaussian", start_pars = gaussian_start,
                             messages = FALSE)
  gaussian_two_calls <- glgpm(y ~ cov + gp() + re(i) + re(j), data = gaussian_data,
                              family = "gaussian", start_pars = gaussian_start,
                              messages = FALSE)
  expect_identical(names(gaussian_two_calls$re), c("i", "j"))
  expect_identical(coef(gaussian_two_calls), coef(gaussian_one_call))
  expect_identical(gaussian_two_calls$covariance, gaussian_one_call$covariance)
  expect_identical(
    predict(gaussian_two_calls, seed = 1, messages = FALSE)$predictions,
    predict(gaussian_one_call, seed = 1, messages = FALSE)$predictions
  )

  fit_binomial <- function(formula) {
    suppressWarnings(glgpm(
      formula, data = binomial_data, family = "binomial",
      denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
      control_mcml = fixture_control_mcml, messages = FALSE
    ))
  }
  binomial_one_call <- fit_binomial(y ~ cov + gp() + re(i, j))
  binomial_two_calls <- fit_binomial(y ~ cov + gp() + re(i) + re(j))
  expect_identical(names(binomial_two_calls$re), c("i", "j"))
  expect_identical(coef(binomial_two_calls), coef(binomial_one_call))
})

test_that("glgpm reads an offset alongside several re() calls and reports unusable offsets", {
  fit <- glgpm(y ~ cov + gp() + re(i) + re(j) + offset(offset), data = gaussian_data,
               family = "gaussian",
               start_pars = list(sigma2_re = c(5, 10), sigma2_me = 0.01),
               messages = FALSE)
  expect_identical(fit$cov_offset, gaussian_data$offset)

  expect_error(
    glgpm(y ~ cov + gp() + offset(log(denominator)), data = gaussian_data,
          family = "gaussian", messages = FALSE),
    "must name a column of `data`"
  )
  expect_error(
    glgpm(y ~ gp() + offset(offset) + offset(cov), data = gaussian_data,
          family = "gaussian", messages = FALSE),
    "only one `offset\\(\\)` term"
  )
})
