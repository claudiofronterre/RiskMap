test_that("glgpm produces errors", {

  # par0 is not checked

  expect_error(
    glgpm("not formula", data = gaussian_data, family = "gaussian"),
    "'formula' must be a 'formula'"
  )

  expect_error(
    glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data, family = "gaussian", messages = FALSE),
    "When there is only one observation per location"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = data, family = "gaussian"),
    "'data' must be of class 'sf'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "invalid"),
    "'family' must be either 'gaussian', 'binomial' or 'poisson'"
  )

  # need to document that invlink can be a list
  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", invlink = function(x) x),
    "'invlink' cannot be provided when 'family' is 'gaussian'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = poisson_data, family = "poisson", invlink = "not func", messages = FALSE),
    "'invlink' must be NULL, a function, or a list"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", denominator = 1),
    "'denominator' cannot be provided when 'family' is 'gaussian'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial", denominator = 1),
    "'denominator' must be provided as an unquoted column name for a column in 'data'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial", denominator = not_present),
    "the variable provided to 'denominator' is not present in 'data'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", model_crs = 12345),
    "The 'model_crs' provided is not a valid CRS"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", distance_units = "not valid"),
    "'distance_units' must be either 'km' or 'm'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = latlon_data, family = "gaussian", model_crs = 4326, messages = FALSE),
    "'model_crs' must be a projected CRS, not longitude/latitude"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", return_samples = "not logical"),
    "'return_samples' must be either TRUE or FALSE"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", messages = "not logical"),
    "'messages' must be either TRUE or FALSE"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial", fix_var_me = 1),
    "'fix_var_me' cannot be provided when 'family' is 'binomial'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", fix_var_me = c(1, 2)),
    "'fix_var_me' must be NULL or a single positive value or zero"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", fix_var_me = "not number"),
    "'fix_var_me' must be NULL or a single positive value or zero"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", fix_var_me = -1),
    "'fix_var_me' must be NULL or a single positive value or zero"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", par0 = 1),
    "'par0' cannot be provided when 'family' is 'gaussian'"
  )


  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(invalid = 1), messages = FALSE),
    "'invalid' is not a valid starting parameter"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(invalid = 1, silly = 2), messages = FALSE),
    "'invalid', 'silly' is not a valid starting parameter"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(beta = 1), messages = FALSE),
    "number of starting values provided for 'beta' do not match"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          control_mcml = list(max_iterations = 2), messages = FALSE),
    "set_control_mcml"
  )

  expect_error(set_control_mcml(max_iterations = 0), "positive integer")
  expect_error(set_control_mcml(tolerance = 0), "positive finite number")
  expect_error(set_control_mcml(min_relative_ess = 2), "between zero and one")

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(beta = c("a", "b")), messages = FALSE),
    "The starting values for 'beta' must be numeric"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2 = -1), messages = FALSE),
    "The starting value for 'sigma2' must be a single positive number"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2 = "a"), messages = FALSE),
    "The starting value for 'sigma2' must be a single positive number"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(phi = -1), messages = FALSE),
    "The starting value for 'phi' must be a single positive number"
  )

  expect_error(
    glgpm(y ~ cov + gp(nugget = TRUE), data = gaussian_data, family = "gaussian",
          fix_var_me = 0, start_pars = list(tau2 = -1), messages = FALSE),
    "The starting value for 'tau2' must be a single positive number"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(tau2 = 1), messages = FALSE),
    "The starting value for 'tau2' cannot be provided when 'nugget'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2_me = -1), messages = FALSE),
    "The starting value for 'sigma2_me' must be a single positive number"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2_re = c(1, 2)), messages = FALSE),
    "Starting values for 'sigma2_re' cannot be provided when no random effects are included in the model"
  )

  expect_error(
    glgpm(y ~ cov + gp() + re(i), data = gaussian_data, family = "gaussian",
          start_pars = list(sigma2_re = c(1, 2)), messages = FALSE),
    "The starting values for 'sigma2_re' do not match the number"
  )

})

expected_output <- c("estimate", "grad_MLE", "covariance", "log_lik",
                     "y", "D", "coords", "ID_coords", "re", "ID_re", "fix_tau2",
                     "fix_var_me", "formula", "family", "distance_units",
                     "data", "input_crs", "kappa", "units_m", "cov_offset", "call",
                     "S_samples", "link_function")
expected_nongaussian_output <- c(
  expected_output, "mcml_history", "mcml_converged"
)

test_that("glgpm produces expected output for gaussian models", {

  fit_no_re <- glgpm(y ~ cov + gp(),
                     data = gaussian_data,
                     family = "gaussian",
                     distance_units = "m",
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_equal(fit_no_re$family, "gaussian")
  expect_equal(fit_no_re$coords[,1], data$x)
  expect_equal(fit_no_re$coords[,2], data$z)
  expect_equal(fit_no_re$y, gaussian_data$y)
  expect_equal(unname(fit_no_re$D[,2]), gaussian_data$cov)
  optimizer <- attr(fit_no_re, "optimizer")
  expect_named(
    optimizer,
    c("convergence", "message", "evaluations", "invalid_evaluations",
      "max_abs_gradient", "stationary", "information_rcond",
      "information_jitter")
  )
  expect_true(is.logical(optimizer$stationary))
  expect_true(optimizer$information_rcond >= 0)

  fit_re_fixed_me <- glgpm(y ~ cov + gp() + re(i),
                           data = gaussian_data,
                           family = "gaussian",
                           fix_var_me = 0.1,
                           messages = FALSE)
  expect_true(is.finite(coef(fit_re_fixed_me)$sigma2_re[["i"]]))

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = gaussian_data,
                  family = "gaussian",
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_equal(fit_re$family, "gaussian")
  expect_length(fit_re$re, 1)
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
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_nongaussian_output)
  expect_equal(fit_no_re$family, "binomial")

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = binomial_data,
                  family = "binomial",
                  denominator = denominator,
                  control_mcmc = control_mcmc,
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_nongaussian_output)
  expect_equal(fit_re$family, "binomial")
})

test_that("glgpm produces expected output for poisson models", {

  fit_no_re <- glgpm(y ~ cov + gp(),
                     data = poisson_data,
                     family = "poisson",
                     control_mcmc = control_mcmc,
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_nongaussian_output)
  expect_equal(fit_no_re$family, "poisson")

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = poisson_data,
                  family = "poisson",
                  control_mcmc = control_mcmc,
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_nongaussian_output)
  expect_equal(fit_re$family, "poisson")

  fit_re_den <- glgpm(y ~ cov + gp() + re(i),
                      data = poisson_data,
                      denominator = denominator,
                      family = "poisson",
                      control_mcmc = control_mcmc,
                      messages = FALSE)

  expect_s3_class(fit_re_den, "RiskMap")
  expect_setequal(names(fit_re_den), expected_nongaussian_output)
  expect_equal(fit_re_den$family, "poisson")
})

test_that("iterative MCML updates its reference and records reproducible history", {
  iterative_control <- set_control_mcml(
    max_iterations = 2,
    tolerance = 1e6,
    min_relative_ess = 0
  )
  fit <- glgpm(
    y ~ cov + gp(),
    data = binomial_data,
    family = "binomial",
    denominator = denominator,
    control_mcmc = control_mcmc,
    control_mcml = iterative_control,
    messages = FALSE
  )

  expect_true(fit$mcml_converged)
  expect_length(fit$mcml_history, 2)
  expect_equal(
    vapply(fit$mcml_history, `[[`, numeric(1), "seed"),
    c(control_mcmc$seed, control_mcmc$seed + 1L)
  )
  expect_lte(fit$mcml_history[[2]]$max_parameter_change,
             iterative_control$tolerance)
  expect_true(is.finite(
    fit$mcml_history[[2]]$log_likelihood_ratio_gain
  ))
  expect_true(is.finite(
    fit$mcml_history[[2]]$max_standardized_change
  ))
  expect_equal(fit$mcml_history[[2]]$estimate, fit$estimate)
})


test_that("glgpm correctly reprojects to new CRS", {

  latlon <- st_transform(gaussian_data, 4326)

  suggested_crs <- propose_utm(latlon)
  sf_reproj <- st_transform(latlon, suggested_crs)
  scaled_coords <- coordinates_in_units(sf_reproj, "km")

  expect_message(
    fit <- glgpm(y ~ cov + gp(),
                 data = latlon,
                 family = "gaussian",
                 distance_units = "km",
                 messages = FALSE),
    "automatically reprojecting to EPSG"
  )

  expect_equal(fit$input_crs, st_crs(latlon))
  expect_equal(fit$coords, scaled_coords)
  expect_equal(st_crs(fit$data), st_crs(suggested_crs))
})

test_that("glgpm honours an explicitly supplied projected model_crs", {

  latlon <- st_transform(gaussian_data, 4326)

  fit <- glgpm(y ~ cov + gp(),
               data = latlon,
               family = "gaussian",
               model_crs = 32637,
               distance_units = "m",
               messages = FALSE)

  expect_equal(st_crs(fit$data), st_crs(32637))
  expect_equal(fit$coords, coordinates_in_units(fit$data, "m"))
})

test_that("coordinates_in_units converts projected CRS units", {

  metre_data <- gaussian_data
  foot_data <- st_transform(gaussian_data, 2263)

  expect_equal(coordinates_in_units(metre_data, "km"),
               st_coordinates(metre_data) / 1000)
  expect_equal(coordinates_in_units(foot_data, "m"),
               st_coordinates(foot_data) * 0.3048006096012192)
})

test_that("plot_mcmc produces errors as expected", {

  expect_error(
    plot_mcmc("not risk"),
    "'object' must be either")

  expect_error(
    plot_mcmc(gaussian_model),
    "'object' is a gaussian model")

  expect_error(
    plot_mcmc(binomial_model),
    "'object' does not contain any MCMC chains")

  expect_warning(
    plot_mcmc(poisson_model, component = 1),
    "if check_mean = TRUE"
  )

  expect_error(
    plot_mcmc(poisson_model, check_mean = FALSE, component = NULL),
    "When 'check_mean' = FALSE a component of the"
  )

  #ncol(poisson_model$S_samples) is 15
  # expect_error(
  #   plot_mcmc(poisson_model, check_mean = FALSE, component = 11),
  #   "'component' must be a positive integer"
  # )
  expect_error(
    plot_mcmc(poisson_model, check_mean = FALSE, component = 20),
    "'component' must be a single positive integer"
  )

  expect_error(
    plot_mcmc(poisson_model, check_mean = FALSE, component = 10.5),
    "'component' must be a single positive integer"
  )

  expect_error(
    plot_mcmc(poisson_model, check_mean = FALSE, component = 0),
    "'component' must be a single positive integer"
  )

  expect_error(
    plot_mcmc(poisson_model, check_mean = FALSE, component = "no"),
    "'component' must be a single positive integer"
  )

  expect_no_error(
    plot_mcmc(poisson_model, check_mean = FALSE, component = 1)
  )

  expect_no_error(
    plot_mcmc(poisson_model, check_mean = FALSE, component = 10)
  )

  poisson_grid <- setup_prediction(poisson_model, control_mcmc = control_mcmc)
  expect_no_error(
    plot_mcmc(poisson_grid)
  )

  expect_no_error(
    plot_mcmc(poisson_grid, check_mean = FALSE, component = 10)
  )

})
