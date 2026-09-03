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
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", den = 1),
    "'den' cannot be provided when 'family' is 'gaussian'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial", den = 1),
    "'den' must be provided as an unquoted column name for a column in 'data'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = binomial_data, family = "binomial", den = not_present),
    "the variable provided to 'den' is not present in 'data'"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", convert_to_crs = 12345),
    "The 'convert_to_crs' provided is not a valid CRS"
  )

  expect_error(
    glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian", scale_to_km = "not logical"),
    "'scale_to_km' must be either TRUE or FALSE"
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
                     "fix_var_me", "formula", "family", "crs", "scale_to_km",
                     "data_sf", "kappa", "units_m", "cov_offset", "call",
                     "S_samples", "link_function")

test_that("glgpm produces expected output for gaussian models", {

  fit_no_re <- glgpm(y ~ cov + gp(),
                     data = gaussian_data,
                     family = "gaussian",
                     scale_to_km = FALSE,
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_equal(fit_no_re$family, "gaussian")
  expect_equal(fit_no_re$coords[,1], data$x)
  expect_equal(fit_no_re$coords[,2], data$z)
  expect_equal(fit_no_re$y, gaussian_data$y)
  expect_equal(unname(fit_no_re$D[,2]), gaussian_data$cov)

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = gaussian_data,
                  family = "gaussian",
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_equal(fit_re$family, "gaussian")
})

test_that("glgpm produces expected output for binomial models", {

  fit_no_re <- glgpm(y ~ cov + gp(),
                     data = binomial_data,
                     family = "binomial",
                     den = den,
                     control_mcmc = control_mcmc,
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_equal(fit_no_re$family, "binomial")

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = binomial_data,
                  family = "binomial",
                  den = den,
                  control_mcmc = control_mcmc,
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_equal(fit_re$family, "binomial")
})

test_that("glgpm produces expected output for poisson models", {

  fit_no_re <- glgpm(y ~ cov + gp(),
                     data = poisson_data,
                     family = "poisson",
                     control_mcmc = control_mcmc,
                     messages = FALSE)

  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_equal(fit_no_re$family, "poisson")

  fit_re <- glgpm(y ~ cov + gp() + re(i),
                  data = poisson_data,
                  family = "poisson",
                  control_mcmc = control_mcmc,
                  messages = FALSE)

  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_equal(fit_re$family, "poisson")

  fit_re_den <- glgpm(y ~ cov + gp() + re(i),
                      data = poisson_data,
                      den = den,
                      family = "poisson",
                      control_mcmc = control_mcmc,
                      messages = FALSE)

  expect_s3_class(fit_re_den, "RiskMap")
  expect_setequal(names(fit_re_den), expected_output)
  expect_equal(fit_re_den$family, "poisson")
})


test_that("glgpm correctly reprojects to new CRS", {

  latlon <- st_transform(gaussian_data, 4326)

  suggested_crs <- propose_utm(latlon)
  sf_reproj <- st_transform(latlon, suggested_crs)
  scaled_coords <- st_coordinates(sf_reproj) / 1000
  fit <- glgpm(y ~ cov + gp(),
               data = latlon,
               family = "gaussian",
               scale_to_km = TRUE,
               messages = FALSE,
               convert_to_crs = suggested_crs)

  expect_equal(fit$coords, scaled_coords)
  expect_equal(fit$crs, suggested_crs)
})

test_that("check_mcmc produces errors as expected", {

  expect_error(
    check_mcmc("not risk"),
    "'object' must be either")

  expect_error(
    check_mcmc(gaussian_model),
    "'object' is a gaussian model")

  expect_error(
    check_mcmc(binomial_model),
    "'object' does not contain any MCMC chains")

  expect_warning(
    check_mcmc(poisson_model, component = 1),
    "if check_mean = TRUE"
  )

  expect_error(
    check_mcmc(poisson_model, check_mean = FALSE, component = NULL),
    "When 'check_mean' = FALSE a component of the"
  )

  #ncol(poisson_model$S_samples) is 15
  # expect_error(
  #   check_mcmc(poisson_model, check_mean = FALSE, component = 11),
  #   "'component' must be a positive integer"
  # )
  expect_error(
    check_mcmc(poisson_model, check_mean = FALSE, component = 20),
    "'component' must be a single positive integer"
  )

  expect_error(
    check_mcmc(poisson_model, check_mean = FALSE, component = 10.5),
    "'component' must be a single positive integer"
  )

  expect_error(
    check_mcmc(poisson_model, check_mean = FALSE, component = 0),
    "'component' must be a single positive integer"
  )

  expect_error(
    check_mcmc(poisson_model, check_mean = FALSE, component = "no"),
    "'component' must be a single positive integer"
  )

  expect_no_error(
    check_mcmc(poisson_model, check_mean = FALSE, component = 1)
  )

  expect_no_error(
    check_mcmc(poisson_model, check_mean = FALSE,component = 10)
  )
})

