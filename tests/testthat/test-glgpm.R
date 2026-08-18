expected_output <- c("y", "D", "coords", "ID_coords", "re", "ID_re", "fix_tau2",
                     "fix_var_me", "formula", "family", "crs", "scale_to_km",
                     "data_sf", "kappa", "units_m", "cov_offset", "call")

# currently missing from docs
# "estimate", "grad.MLE", "covariance", "log.lik", "S_samples", "linkf", "sst"

test_that("glgpm produces expected output for gaussian models", {

  fit_no_re <- glgpm(y ~ cov + gp(), data = sf_data, family = "gaussian", scale_to_km = FALSE, messages = FALSE)
  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_equal(fit_no_re$family, "gaussian")
  expect_equal(fit_no_re$coords[,1], sf_data$lon)
  expect_equal(fit_no_re$coords[,2], sf_data$lat)
  expect_equal(fit_no_re$y, sf$y)
  expect_equal(unname(fit_no_re$D[,2]), sf_data$cov)

  fit_re <- glgpm(y ~ cov + gp() + re(i), data = sf_data, family = "gaussian", messages = FALSE)
  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_equal(fit_re$family, "gaussian")
})

test_that("glgpm produces expected output for binomial models", {

  eta <- 0.2 + 0.3 * data$cov + S
  p <- plogis(eta)
  data$y <- rbinom(n, size = data$den, prob = p)
  sf <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326)

  fit_no_re <- glgpm(y ~ cov + gp(), data = sf, family = "binomial",
                 den = den, messages = FALSE)
  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_equal(fit_no_re$family, "binomial")

  fit_re <- glgpm(y ~ cov + gp() + re(i), data = sf, family = "binomial",
                     den = den, messages = FALSE)
  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_equal(fit_re$family, "binomial")
})

test_that("glgpm produces expected output for poisson models", {

  lambda <- exp(0.1 + 0.2 * data$cov + S)
  data$y <- rpois(n, lambda)
  sf <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326)

  fit_no_re <- glgpm(y ~ cov + gp(), data = sf, family = "poisson", messages = FALSE)
  expect_s3_class(fit_no_re, "RiskMap")
  expect_setequal(names(fit_no_re), expected_output)
  expect_equal(fit_no_re$family, "poisson")

  fit_re <- glgpm(y ~ cov + gp() + re(i), data = sf, family = "poisson", messages = FALSE)
  expect_s3_class(fit_re, "RiskMap")
  expect_setequal(names(fit_re), expected_output)
  expect_equal(fit_re$family, "poisson")

  fit_re_den <- glgpm(y ~ cov + gp() + re(i), den = den, data = sf, family = "poisson", messages = FALSE)
  expect_s3_class(fit_re_den, "RiskMap")
  expect_setequal(names(fit_re_den), expected_output)
  expect_equal(fit_re_den$family, "poisson")
})


test_that("glgpm correctly reprojects to new CRS", {

  suggested_crs <- propose_utm(sf_data)
  sf_reproj <- sf::st_transform(sf_data, suggested_crs)
  scaled_coords <- sf::st_coordinates(sf_reproj) / 1000
  fit <- glgpm(y ~ cov + gp(), data = sf_data, family = "gaussian", scale_to_km = TRUE, messages = FALSE, convert_to_crs = suggested_crs)

  expect_equal(fit$coords, scaled_coords)
  expect_equal(fit$crs, suggested_crs)
})
