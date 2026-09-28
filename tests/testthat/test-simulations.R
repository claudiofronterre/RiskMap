test_that("joint simulations share effects at overlapping locations", {
  data <- gaussian_data[c(1, 1, 2), ]
  model <- specify_glgpm(y ~ cov + gp(nugget = TRUE) + re(i), data,
                         "gaussian", list(beta = c(1, 0.5), sigma2 = 1,
                                          phi = 2, tau2 = 0.2, sigma2_re = 0.3,
                                          sigma2_me = 0.1))
  sim <- simulate_glgpm(model, nsim = 2, what = c("data", "surface"),
                        prediction_grid = data[3:1, ], seed = 10)
  for (component in c("spatial_effect", "nugget", "group_effect", "linear_predictor")) {
    expect_equal(sim$samples$data[1, , component], sim$samples$data[2, , component])
    expect_equal(sim$samples$data[, , component], sim$samples$surface[3:1, , component])
  }
  expect_false(identical(sim$samples$data[1, , "response"], sim$samples$data[2, , "response"]))
  expect_false("response" %in% dimnames(sim$samples$surface)[[3]])
  expect_equal(simulated_data(sim, 2)$y, as.numeric(sim$samples$data[, 2, "response"]))
  expect_equal(simulated_data(sim, 2)$cov, data$cov)
  tidy <- simulated_values(sim, "spatial_effect", simulation = 2)
  expect_equal(nrow(tidy), 6L)
  expect_equal(tidy$value, c(sim$samples$data[, 2, "spatial_effect"],
                             sim$samples$surface[, 2, "spatial_effect"]))
  expect_equal(simulated_surface(sim, 1)$mean, as.numeric(sim$samples$surface[, 1, "mean"]))
  expect_equal(nrow(simulated_values(sim, "response")), 6L)
})

test_that("simulation model validation remains informative", {
  expect_error(
    specify_glgpm("not a formula", gaussian_data, "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "'formula' must be a 'formula'"
  )
  expect_error(
    specify_glgpm(y ~ gp(), "not spatial data", "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "class 'sf'"
  )
  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "not-a-family",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "should be one of"
  )
  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian",
                  list(sigma2 = 1, phi = 1, sigma2_me = 0)),
    "Missing simulation parameters: beta"
  )
  expect_error(simulate_glgpm(gaussian_model, nsim = -1),
               "positive integer")
})

test_that("joint draws follow the specified covariance and distance units", {
  data <- gaussian_data[1:2, ]
  model <- specify_glgpm(y ~ gp(), data, "gaussian",
                         list(beta = 0, sigma2 = 1, phi = 2, sigma2_me = 0))
  sim <- simulate_glgpm(model, nsim = 5000, seed = 10)
  expected <- matern_correlation(dist(st_coordinates(data) / 1000),
                                  phi = 2, kappa = 0.5, return_sym_matrix = TRUE)
  expect_equal(cov(t(sim$samples$data[, , "spatial_effect"])), expected, tolerance = 0.05)
  model_m <- specify_glgpm(y ~ gp(), data, "gaussian",
                           list(beta = 0, sigma2 = 1, phi = 2000, sigma2_me = 0),
                           distance_units = "m")
  expect_equal(simulate_glgpm(model_m, nsim = 2, seed = 1)$samples,
               simulate_glgpm(model, nsim = 2, seed = 1)$samples)
})

test_that("fitted models retain offsets, parameters, denominators and links", {
  for (model in list(gaussian_model, gaussian_offset_model, binomial_model, poisson_model)) {
    sim <- simulate_glgpm(model, nsim = 2, seed = 4)
    x <- sim$samples$data
    fixed <- as.numeric(model$D %*% coef(model)$beta) + model$cov_offset
    expect_equal(x[, , "linear_predictor"] - x[, , "spatial_effect"] -
                   x[, , "nugget"] - x[, , "group_effect"],
                 matrix(fixed, nrow(model$data), 2))
    inverse <- if (model$family == "gaussian") identity else model$link_function$inv
    expect_equal(x[, , "mean"], inverse(x[, , "linear_predictor"]))
    if (model$family == "binomial") {
      expect_true(all(x[, , "response"] >= 0 & x[, , "response"] <= model$units_m))
    }
  }
  custom <- binomial_model
  custom$link_function$inv <- function(x) rep(0.25, length(x))
  expect_true(all(simulate_glgpm(custom, seed = 1)$samples$data[, , "mean"] == 0.25))
})

test_that("Poisson responses use exposure times the exponential mean", {
  model <- specify_glgpm(y ~ gp() + offset(offset), poisson_data, "poisson",
                         list(beta = log(3), sigma2 = 0, phi = 1),
                         denominator = denominator)
  sim <- simulate_glgpm(model, nsim = 2000, seed = 1)
  expected <- 3 * exp(poisson_data$offset)
  expect_equal(sim$samples$data[, 1, "mean"], expected)
  expect_equal(rowMeans(sim$samples$data[, , "response"]) / poisson_data$denominator,
               expected, tolerance = 0.03)
})

test_that("validation is informative and local seeds preserve RNG", {
  set.seed(12)
  rng <- .Random.seed
  first <- simulate_glgpm(gaussian_model, seed = 3)
  expect_identical(.Random.seed, rng)
  expect_identical(first$samples, simulate_glgpm(gaussian_model, seed = 3)$samples)
  expect_error(simulate_glgpm(gaussian_model, nsim = Inf), "positive integer")
  expect_error(simulate_glgpm(gaussian_model, what = "surface"), "prediction_grid")
  expect_message(
    transformed <- simulate_glgpm(gaussian_model,
                                  sample_locations = latlon_data),
    "transformed to the model CRS"
  )
  expect_identical(st_crs(transformed$locations$data), st_crs(gaussian_model$data))
  expect_error(simulate_glgpm(binomial_model, sample_locations = binomial_data[, "cov"]),
               "Random-effect|denominator")
  expect_error(simulated_data(first, 2), "valid simulation")
  expect_error(simulated_surface(first), "No surface")
  expect_error(simulate_glgpm(list(gaussian_model)), "fitted RiskMap")
})

test_that("new locations preserve factor coding without observed responses", {
  data <- gaussian_data
  data$group <- factor(rep(c("a", "b"), 5))
  model <- specify_glgpm(y ~ group + gp(), data, "gaussian",
                         list(beta = c(1, 2), sigma2 = 0, phi = 1, sigma2_me = 0))
  new <- data[2, "group", drop = FALSE]
  sim <- simulate_glgpm(model, sample_locations = new, seed = 1)
  expect_equal(simulated_data(sim)$y, 3)
  surface <- simulate_glgpm(model, what = "surface", prediction_grid = new, seed = 1)
  expect_equal(simulated_surface(surface)$mean, 3)
})

test_that("assess_simulation validates area boundaries", {
  obj_sim <- structure(list(), class = "RiskMap_simulation")
  expect_error(assess_simulation(obj_sim, models = list(model = y ~ 1),
                                 spatial_scale = "area", boundaries = gaussian_data,
                                 f_area_target = mean),
               "'boundaries' can only contain 'POLYGON' or 'MULTIPOLYGON' geometry")
})

test_that("joint output feeds the existing grid assessment", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)
  result <- assess_simulation(sim, models = list(intercept = y ~ gp()),
                              control_mcmc = control_mcmc, spatial_scale = "grid",
                              f_grid_target = identity, pred_objective = "mse",
                              messages = FALSE)
  expect_equal(dim(result$pred_objective$mse), c(1L, 2L))
  expect_true(all(is.finite(result$pred_objective$mse)))
})
