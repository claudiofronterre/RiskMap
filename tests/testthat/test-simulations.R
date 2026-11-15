test_that("joint simulations share effects at overlapping locations", {
  data <- gaussian_data[c(1, 1, 2), ]
  model <- specify_glgpm(y ~ cov + gp(nugget = TRUE) + re(i), data,
                         "gaussian", list(beta = c(1, 0.5), sigma2 = 1,
                                          phi = 2, tau2 = 0.2, sigma2_re = 0.3,
                                          sigma2_me = 0.1))
  sim <- simulate_glgpm(model, nsim = 2, what = c("data", "surface"),
                        prediction_grid = data[3:1, ], seed = 10)
  for (component in c("spatial_effect", "nugget", "group_effect", "linear_predictor")) {
    expect_identical(sim$samples$data[1, , component], sim$samples$data[2, , component])
    expect_identical(sim$samples$data[, , component], sim$samples$surface[3:1, , component])
  }
  expect_false(identical(sim$samples$data[1, , "response"], sim$samples$data[2, , "response"]))
  expect_false("response" %in% dimnames(sim$samples$surface)[[3]])
  expect_identical(simulated_data(sim, 2)$y, as.numeric(sim$samples$data[, 2, "response"]))
  expect_identical(simulated_data(sim, 2)$cov, data$cov)
  tidy <- simulated_values(sim, "spatial_effect", simulation = 2)
  expect_identical(nrow(tidy), 6L)
  expect_identical(tidy$value, c(sim$samples$data[, 2, "spatial_effect"],
                             sim$samples$surface[, 2, "spatial_effect"]))
  expect_identical(simulated_surface(sim, 1)$mean, as.numeric(sim$samples$surface[, 1, "mean"]))
  expect_identical(nrow(simulated_values(sim, "response")), 6L)
  expect_s3_class(tidy, "tbl_df")
  expect_named(tidy, c("simulation", "location_set", "location_id",
                       "component", "value"))
})

test_that("specify_glgpm takes invlink after parameters and the CRS in distance_crs", {
  expect_identical(names(formals(specify_glgpm)),
                   c("formula", "data", "family", "parameters", "invlink",
                     "denominator", "distance_crs", "distance_units"))
  parameters <- list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0.1)
  # A CRS given in distance_crs is used without a message.
  expect_no_message(
    model <- specify_glgpm(y ~ gp(), latlon_data, "gaussian", parameters,
                           distance_crs = 32637)
  )
  expect_identical(st_crs(model$data), st_crs(32637))
  expect_equal(st_coordinates(model$data), st_coordinates(gaussian_data),
               tolerance = 1e-6)
})

test_that("specify_glgpm returns errors as expected", {
  expect_input_error(
    specify_glgpm("not a formula", gaussian_data, "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "formula"
  )

  expect_error(
    specify_glgpm(y ~ cov, gaussian_data, "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "must contain a Gaussian process term"
  )

  expect_error(
    specify_glgpm(y ~ gp(), "not spatial data", "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "must be an <sf> object"
  )

  expect_input_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian", distance_crs = 4326,
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "distance_crs"
  )

  expect_message(
    specify_glgpm(y ~ gp(), latlon_data, "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "are in longitude and latitude"
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "not-a-family",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "family", class = "riskmap_error_input"
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian",
                  denominator = cov,
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "can't be given for a Gaussian model"
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian",
                  invlink = plogis,
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "can't be given for a Gaussian model"
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian",
                  list(sigma2 = 1, phi = 1, sigma2_me = 0)),
    "must have the element \"beta\""
  )

  expect_input_error(
    specify_glgpm(y ~ gp(), binomial_data, "binomial",
                  denominator = "d",
                  list(beta = 1, sigma2 = 1, phi = 1)),
                  "denominator"
  )

  expect_input_error(
    specify_glgpm(y ~ gp(), binomial_data, "binomial",
                  denominator = d,
                  list(beta = 1, sigma2 = 1, phi = 1)),
    "denominator"
  )

  expect_input_error(
    specify_glgpm(y ~ gp(), binomial_data, "binomial",
                  denominator = denominator,
                  invlink = "a",
                  list(beta = 1, sigma2 = 1, phi = 1)),
    "invlink"
  )

})

test_that("specify_glgpm returns errors as expected", {

  expect_error(
    simulate_glgpm(gaussian_model, nsim = -1),
    "positive integer"
  )

  expect_input_error(
    simulate_glgpm(gaussian_model, nsim = 1, what = "n"),
    "what"
  )

  expect_input_error(
    simulate_glgpm(gaussian_model, nsim = 1, what = "surface",
                   sample_locations = gaussian_data),
    "sample_locations"
  )

  expect_input_error(
    simulate_glgpm(gaussian_model, nsim = 1, what = "data",
                   prediction_grid = grid),
    "prediction_grid"
  )

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
    expect_identical(x[, , "mean"], inverse(x[, , "linear_predictor"]))
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
    "projected to the coordinate reference system of the model"
  )
  expect_identical(st_crs(transformed$locations$data), st_crs(gaussian_model$data))
  expect_input_error(simulate_glgpm(binomial_model, sample_locations = binomial_data[, "cov"]),
               "sample_locations")
  expect_error(simulated_data(first, 2), "must be one simulation number")
  expect_error(simulated_surface(first), "has no simulated surface")
  expect_error(simulate_glgpm(list(gaussian_model)), "<RiskMap> fit")
})

test_that("missing formula variables in prediction_grid or sample_locations are named clearly (#204)", {
  surface_missing_cov <- grid
  expect_input_error(
    simulate_glgpm(gaussian_model, what = c("data", "surface"),
                   prediction_grid = surface_missing_cov,
                   sample_locations = gaussian_data),
    "prediction_grid"
  )

  surface_missing_re <- grid
  surface_missing_re$cov <- rnorm(nrow(grid))
  expect_input_error(
    simulate_glgpm(gaussian_model, what = c("data", "surface"),
                   prediction_grid = surface_missing_re,
                   sample_locations = gaussian_data),
    "prediction_grid"
  )

  expect_input_error(
    simulate_glgpm(gaussian_model, sample_locations = gaussian_data[, "cov"]),
    "sample_locations"
  )
})

test_that("new locations preserve factor coding without observed responses", {
  data <- gaussian_data
  data$group <- factor(rep(c("a", "b"), 5))
  model <- specify_glgpm(y ~ group + gp(), data, "gaussian",
                         list(beta = c(1, 2), sigma2 = 0, phi = 1, sigma2_me = 0))
  new <- data[2, "group", drop = FALSE]
  sim <- simulate_glgpm(model, sample_locations = new, seed = 1)
  expect_identical(simulated_data(sim)$y, 3)
  surface <- simulate_glgpm(model, what = "surface", prediction_grid = new, seed = 1)
  expect_identical(simulated_surface(surface)$mean, 3)
})

test_that("assess_simulation validates area boundaries", {
  simulation <- structure(list(), class = "RiskMap_simulation")
  expect_error(
    assess_simulation(simulation, models = list(model = y ~ gp()),
                      spatial_scale = "location"),
    "`glgpm\\(\\)` or models from `specify_glgpm\\(\\)`"
  )
  expect_input_error(
    assess_simulation(simulation, models = list(model = gaussian_intercept_model),
                      spatial_scale = character()),
    "spatial_scale"
  )
  expect_input_error(
    assess_simulation(simulation, models = list(model = gaussian_intercept_model),
                      spatial_scale = c("location", "location")),
    "spatial_scale"
  )
  expect_input_error(
    assess_simulation(simulation, models = list(model = gaussian_intercept_model),
                      spatial_scale = "location", objective = character()),
    "objective"
  )
  expect_input_error(assess_simulation(simulation, models = list(model = gaussian_intercept_model),
                                 spatial_scale = "area", boundaries = gaussian_data,
                                 aggregate = "mean"),
               "boundaries")
})

test_that("assessment refits use fitted templates without reusing estimates #108", {
  simulated <- binomial_data
  simulated$y <- rev(binomial_data$y)
  simulated$units_m <- rev(binomial_data$denominator)
  formula <- update(binomial_model$formula, y ~ .)
  args <- assessment_refit_args(binomial_model, formula, simulated, control_mcmc)

  expect_identical(args$formula, formula)
  expect_identical(args$family, "binomial")
  expect_identical(args$data, simulated)
  expect_identical(args$denominator, quote(units_m))
  expect_identical(args$control_mcmc, control_mcmc)
  expect_false("start_pars" %in% names(args))
  expect_false("estimate" %in% names(args))
})

test_that("assessment refits accept specify_glgpm templates #108", {
  template <- specify_glgpm(
    y ~ cov + gp(), gaussian_data, "gaussian",
    parameters = list(beta = c(1, 0.5), sigma2 = 1, phi = 2,
                      sigma2_me = 0.1)
  )
  formula <- update(template$formula, y ~ .)
  args <- assessment_refit_args(
    template, formula, gaussian_data, control_mcmc
  )

  expect_identical(args$formula, formula)
  expect_identical(args$family, "gaussian")
  expect_identical(args$data, gaussian_data)
  expect_null(args$fixed_pars)
  expect_false("start_pars" %in% names(args))

  sim <- simulate_glgpm(
    template, nsim = 1, what = c("data", "surface"),
    prediction_grid = gaussian_data, seed = 2
  )
  result <- assess_simulation(
    sim, models = list(candidate = template), spatial_scale = "location",
    objective = "mse", messages = FALSE
  )
  expect_identical(dim(result$objective$location$mse), c(1L, 1L))
})

test_that("joint output feeds the existing grid assessment", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)
  result <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                              control_mcmc = control_mcmc,
                              seed = fixture_seed, spatial_scale = "location",
                              objective = "mse",
                              messages = FALSE)
  expect_identical(dim(result$objective$location$mse), c(1L, 2L))
  expect_true(all(is.finite(result$objective$location$mse)))
})

test_that("assess_simulation validates target and area function outputs #171", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)

  expect_input_error(
    assess_simulation(
      sim,
      models = list(intercept = gaussian_intercept_model),
      spatial_scale = "location",
      target = "plogis",
      objective = "mse",
      messages = FALSE
    ),
    "target"
  )

  expect_error(
    assess_simulation(
      sim,
      models = list(intercept = gaussian_intercept_model),
      spatial_scale = "location",
      target = mean,
      objective = "mse",
      messages = FALSE
    ),
    "same dimensions as its input"
  )

  expect_error(
    assess_simulation(
      sim,
      models = list(intercept = gaussian_intercept_model),
      spatial_scale = "location",
      target = function(x) x * NA_real_,
      objective = "mse",
      messages = FALSE
    ),
    "finite numeric matrix"
  )

  boundaries <- create_convex_hull(gaussian_data)
  expect_error(
    assess_simulation(
      sim,
      models = list(intercept = gaussian_intercept_model),
      spatial_scale = "area",
      aggregate = function(values, weights) values,
      boundaries = boundaries,
      objective = "mse",
      messages = FALSE
    ),
    "one number for each predictive draw"
  )
  expect_input_error(
    assess_simulation(
      sim,
      models = list(intercept = gaussian_intercept_model),
      spatial_scale = "area",
      aggregate = "average",
      boundaries = boundaries,
      objective = "mse",
      messages = FALSE
    ),
    "aggregate"
  )
})

test_that("assess_simulation computes grid and area objectives in one combined run #109", {
  boundaries <- create_convex_hull(gaussian_data)
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)

  combined <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                control_mcmc = control_mcmc,
                                seed = fixture_seed,
                                spatial_scale = c("location", "area"),
                                aggregate = "mean",
                                boundaries = boundaries, objective = "mse",
                                messages = FALSE)

  expect_setequal(names(combined$objective), c("location", "area"))
  expect_identical(dim(combined$objective$location$mse), c(1L, 2L))
  expect_identical(dim(combined$objective$area$mse), c(1L, 2L))
  expect_true(all(is.finite(combined$objective$location$mse)))
  expect_true(all(is.finite(combined$objective$area$mse)))

  grid_only <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                 control_mcmc = control_mcmc,
                                 seed = fixture_seed, spatial_scale = "location",
                                 objective = "mse",
                                 messages = FALSE)
  area_only <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                 control_mcmc = control_mcmc,
                                 seed = fixture_seed, spatial_scale = "area",
                                 aggregate = "mean",
                                 boundaries = boundaries, objective = "mse",
                                 messages = FALSE)

  # A function given as `aggregate` receives the locations by draws and
  # their weights, as in predict().
  custom_mean <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                   control_mcmc = control_mcmc,
                                   seed = fixture_seed, spatial_scale = "area",
                                   aggregate = function(values, weights) {
                                     colSums(values * weights) / sum(weights)
                                   },
                                   boundaries = boundaries, objective = "mse",
                                   messages = FALSE)
  expect_equal(custom_mean$objective$area$mse, area_only$objective$area$mse)

  # Area-level results are computed identically whether or not 'grid' is
  # also requested: both paths use the same (joint-type) prediction, so
  # under the shared control_mcmc seed the two runs must agree exactly.
  expect_identical(combined$objective$area$mse, area_only$objective$area$mse)

  # A standalone grid-only run instead uses cheaper marginal-type sampling,
  # since it doesn't need spatially correlated draws across the grid; a
  # combined run always uses joint-type sampling so one prediction can serve
  # both scales. The two are therefore consistent estimates of the same
  # quantity, not bitwise-identical draws.
  expect_equal(combined$objective$location$mse, grid_only$objective$location$mse,
               tolerance = 0.05)

  s <- summary(combined)
  expect_setequal(names(s), c("location", "area", "fits"))
  expect_s3_class(s, "summary.RiskMap_assess_simulation")
  # Each spatial scale has its own section, the locations first.
  expect_match(squish(capture.output(print(s))),
               "Locations -+ Mean squared error .* Areas -+ Mean squared error")

  combined_classify <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                         control_mcmc = control_mcmc,
                                         seed = fixture_seed,
                                         spatial_scale = c("location", "area"),
                                         aggregate = "mean",
                                         boundaries = boundaries,
                                         objective = c("mse", "classify"),
                                         breaks = c(-3, -1, 0, 1, 3),
                                         messages = FALSE)

  expect_length(combined_classify$objective$location$classify$intercept$by_cat, 2L)
  expect_length(combined_classify$objective$area$classify$intercept$by_cat, 2L)
  expected_classes <- c("(-3,-1]", "(-1,0]", "(0,1]", "(1,3]")
  expect_identical(
    levels(combined_classify$objective$location$classify$Class),
    expected_classes
  )

  expect_error(
    assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                      control_mcmc = control_mcmc, seed = fixture_seed,
                      spatial_scale = "location",
                      objective = "classify", breaks = c(-1, 0, 0),
                      messages = FALSE),
    "finite, strictly increasing"
  )
})

test_that("classification metrics retain fixed categories and compute accuracy #108", {
  samples <- rbind(
    rep(0.25, 4),
    rep(0.25, 4),
    rep(0.75, 4),
    rep(0.75, 4)
  )
  metrics <- simulation_classification_metrics(
    true_values = c(0.25, 0.25, 0.25, 0.75),
    samples = samples,
    breaks = c(0, 0.5, 1),
    labels = c("low", "high")
  )

  expect_identical(metrics$by_cat$Class, c("low", "high"))
  expect_identical(metrics$by_cat$Sensitivity, c(2 / 3, 1))
  expect_identical(metrics$by_cat$CC, c(3 / 4, 3 / 4))
  expect_identical(metrics$overall_cc, 3 / 4)

  no_high_predictions <- simulation_classification_metrics(
    true_values = c(0.25, 0.75),
    samples = matrix(0.25, nrow = 2, ncol = 4),
    breaks = c(0, 0.5, 1),
    labels = c("low", "high")
  )
  expect_identical(no_high_predictions$by_cat$Class, c("low", "high"))
  expect_true(is.na(no_high_predictions$by_cat$PPV[2]))
})

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
      objective = list(
        location = list(
          classify = list(
            model = list(by_cat = list(one_result), CC = 0.65),
            Class = factor(c("low", "high"), levels = c("low", "high"))
          )
        )
      ),
      n_sim = 1L,
      spatial_scale = "location"
    ),
    class = "RiskMap_assess_simulation"
  )

  expect_warning(summary_one <- summary(one_simulation), "one simulation")
  model_one <- summary_one$location$classify$model
  expect_equal(model_one$classify_res$Sensitivity, c(0.5, NA))
  expect_equal(model_one$n_valid$Sensitivity, c(1L, 0L))
  expect_equal(model_one$cc_summary$mean, 0.65)
  expect_true(is.na(model_one$cc_summary$sd))
  expect_true(is.na(model_one$cc_summary$lower))
  expect_output(print(summary_one), "uncertainty unavailable")

  second_result <- one_result
  second_result$Sensitivity <- c(NA, 0.25)
  second_result$Specificity <- c(0.5, NA)
  two_simulations <- one_simulation
  two_simulations$n_sim <- 2L
  two_simulations$objective$location$classify$model$by_cat <-
    list(one_result, second_result)
  two_simulations$objective$location$classify$model$CC <- c(0.65, 0.75)

  summary_two <- summary(two_simulations)$location$classify$model
  expect_equal(summary_two$classify_res$Sensitivity, c(0.5, 0.25))
  expect_equal(summary_two$n_valid$Sensitivity, c(1L, 1L))
  expect_equal(summary_two$n_valid$Specificity, c(1L, 1L))
  expect_equal(summary_two$cc_summary$n_valid, 2L)
  expect_equal(summary_two$cc_summary$n_sim, 2L)
})

test_that("plot.RiskMap_simulation maps the simulated mean by default (#87)", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = grid, seed = 1)

  expect_s3_class(sim, "RiskMap_simulation")
  expect_error(plot(sim, simulation = 3), "must be one simulation number")
  expect_no_error(plot(sim, simulation = 2))

  p <- plot(sim)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(p$labels$title, "Simulation 1")
  expect_equal(p$labels$fill, "Mean")
  expect_equal(p$data$value, as.numeric(sim$samples$surface[, 1, "mean"]))
  # sample locations are overlaid when data were simulated
  expect_length(p$layers, 2)
  expect_equal(nrow(p$layers[[2]]$data), nrow(sim$locations$data))

  p2 <- plot(sim, simulation = 2, quantity = "spatial_effect")
  expect_equal(p2$labels$title, "Simulation 2")
  expect_equal(p2$labels$fill, "Spatial effect")
  expect_equal(p2$data$value,
               as.numeric(sim$samples$surface[, 2, "spatial_effect"]))
  # A signed effect has a scale centred at zero.
  limits <- ggplot2::ggplot_build(p2)$plot$scales$get_scales("fill")$limits
  expect_equal(limits[[1L]], -limits[[2L]])

  linear <- plot(sim, quantity = "linear_predictor", limits = c(-5, 5))
  expect_equal(linear$labels$fill, "Linear predictor")
  expect_equal(linear$scales$get_scales("fill")$limits, c(-5, 5))

  # The model has no nugget, so the simulated nugget is zero.
  nugget <- plot(sim, quantity = "nugget")
  expect_equal(nugget$data$value, rep(0, nrow(grid)))
  expect_no_error(ggplot2::ggplot_build(nugget))
  expect_input_error(plot(sim, quantity = "response"), "quantity")
  expect_input_error(plot(sim, quantity = c("mean", "spatial_effect")),
                     "quantity")
  expect_input_error(plot(sim, limits = c(1, 0)), "limits")
})

test_that("plot.RiskMap_simulation labels the prevalence of a binomial model (#87)", {
  model <- specify_glgpm(y ~ gp(), binomial_data, "binomial",
                         list(beta = -1, sigma2 = 1, phi = 1),
                         denominator = denominator)
  sim <- simulate_glgpm(model, nsim = 1, what = "surface",
                        prediction_grid = grid, seed = 1)
  p <- plot(sim)
  expect_equal(p$labels$fill, "Prevalence")
  values <- as.numeric(sim$samples$surface[, 1, "mean"])
  expect_equal(p$data$value, values)
  built <- ggplot2::ggplot_build(p)
  labels <- built$plot$scales$get_scales("fill")$get_labels()
  expect_true(all(grepl("%$", labels[!is.na(labels)])))
})

test_that("plot.RiskMap_simulation omits points when only a surface was simulated (#183)", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 1, what = "surface",
                        prediction_grid = grid, seed = 1)

  p <- plot(sim)
  expect_length(p$layers, 1)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("simulated_surface() output has a plot method matching plot.RiskMap_simulation (#183)", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = grid, seed = 1)
  surface <- simulated_surface(sim, 2)

  expect_s3_class(surface, c("RiskMap_simulated_surface", "sf"))
  p <- plot(surface, palette = "-Blues")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(p$data$value, surface$mean)
  expect_equal(p$labels$fill, "Mean")
  expect_identical(ggplot2::ggplot_build(p)$data[[1]]$fill,
                   ggplot2::ggplot_build(plot(sim, simulation = 2,
                                              palette = "-Blues"))$data[[1]]$fill)
  # A reversed palette is the named palette in reverse.
  reversed <- ggplot2::ggplot_build(plot(surface, palette = "-Blues"))
  forward <- ggplot2::ggplot_build(plot(surface, palette = rev(
    grDevices::hcl.colors(256L, "Blues"))))
  expect_identical(reversed$data[[1]]$fill, forward$data[[1]]$fill)
  expect_input_error(plot(surface, reverse_palette = TRUE), "reverse_palette")
  # The map sets the breaks and labels of its fill scale itself.
  expect_input_error(plot(sim, breaks = c(0, 0.5, 1)), "breaks")
  expect_input_error(plot(surface, labels = NULL), "labels")
  expect_no_error(ggplot2::ggplot_build(plot(surface, transform = "reverse")))
  expect_input_error(plot(surface, palette = "not a palette"), "palette")
  expect_error(plot.RiskMap_simulated_surface(sim), "<RiskMap_simulated_surface>")
})


test_that("assess_simulation records a failed fit instead of stopping", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 3,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)
  real_glgpm <- glgpm
  calls <- new.env()
  calls$n <- 0
  local_mocked_bindings(
    glgpm = function(...) {
      calls$n <- calls$n + 1
      if (calls$n == 2) stop("fit failed")
      real_glgpm(...)
    },
    .package = "RiskMap"
  )
  expect_warning(
    result <- assess_simulation(
      sim, models = list(intercept = gaussian_intercept_model),
      control_mcmc = control_mcmc, seed = fixture_seed, spatial_scale = "location",
      objective = "mse", messages = FALSE
    ),
    "1 of 3 model fits failed"
  )

  expect_identical(result$fits$status, c("ok", "failed", "ok"))
  expect_identical(result$fits$message[2], "fit failed")
  expect_true(is.na(result$objective$location$mse[1, 2]))
  expect_true(all(is.finite(result$objective$location$mse[1, c(1, 3)])))
  expect_true(all(is.na(result$estimates$intercept[2, ])))
  expect_output(print(summary(result)),
                "1 of 3 model fits failed, so their results are missing.",
                fixed = TRUE)
  expect_output(print(result), "1 of 3 model fits failed.", fixed = TRUE)
})

test_that("assess_simulation leaves out areas without grid locations", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)
  hull <- st_geometry(create_convex_hull(gaussian_data))
  boundaries <- st_sf(
    region = c("near", "far"),
    geometry = c(hull, st_set_crs(hull + c(1e6, 1e6), st_crs(hull)))
  )
  expect_warning(
    result <- assess_simulation(
      sim, models = list(intercept = gaussian_intercept_model),
      control_mcmc = control_mcmc, seed = fixture_seed, spatial_scale = "area",
      aggregate = "mean",
      boundaries = boundaries, label = "region",
      objective = c("mse", "classify"), breaks = c(-10, 0, 10),
      messages = FALSE
    ),
    "\"far\", so it is left out of the assessment of areas"
  )
  expect_true(all(is.finite(result$objective$area$mse)))
})

test_that("assess_simulation requires unique area names", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 1,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)
  hull <- st_geometry(create_convex_hull(gaussian_data))
  boundaries <- st_sf(region = c("a", "a"), geometry = c(hull, hull))
  expect_error(
    assess_simulation(
      sim, models = list(intercept = gaussian_intercept_model),
      spatial_scale = "area",
      aggregate = "mean", boundaries = boundaries, label = "region",
      objective = "mse", messages = FALSE
    ),
    "distinct names of areas"
  )
})

test_that("assess_simulation names missing candidate covariates on the grid", {
  grid_with_gap <- gaussian_data
  grid_with_gap$cov[1] <- NA
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 1,
                        what = c("data", "surface"),
                        prediction_grid = grid_with_gap, seed = 2)
  candidate <- specify_glgpm(
    y ~ cov + gp(), gaussian_data, family = "gaussian",
    parameters = list(beta = c(1, 0.5), sigma2 = 1, phi = 2, sigma2_me = 0.1)
  )
  expect_error(
    assess_simulation(
      sim, models = list(candidate = candidate), spatial_scale = "location",
      objective = "mse", messages = FALSE
    ),
    "needs complete values of \"cov\" on the prediction grid"
  )
})

test_that("assess_simulation scores candidates without the generating model's nugget", {
  generating <- specify_glgpm(
    y ~ gp(nugget = TRUE), data = binomial_data, family = "binomial",
    denominator = denominator,
    parameters = list(beta = 0, sigma2 = 0.5, phi = 2, tau2 = 0.2)
  )
  sim <- simulate_glgpm(generating, nsim = 2, what = c("data", "surface"),
                        prediction_grid = binomial_data, seed = 3)
  candidate <- specify_glgpm(
    y ~ gp(), data = binomial_data, family = "binomial",
    denominator = denominator, parameters = list(beta = 0, sigma2 = 0.5, phi = 2)
  )
  result <- suppressWarnings(assess_simulation(
    sim, models = list(no_nugget = candidate), control_mcmc = control_mcmc,
    seed = fixture_seed,
    spatial_scale = "location", objective = "mse",
    messages = FALSE
  ))
  expect_identical(result$fits$status, c("ok", "ok"))
  expect_true(all(is.finite(result$objective$location$mse)))
})

# Runs assess_simulation() and records, for every fit in the order of the
# fits, the name of its link and its predictive draws of the linear predictor.
assess_recording_links <- function(...) {
  recorded <- new.env()
  recorded$links <- character()
  recorded$draws <- list()
  real_components <- prediction_sample_components
  local_mocked_bindings(
    prediction_sample_components = function(object, ...) {
      components <- real_components(object, ...)
      recorded$links <- c(recorded$links,
                          object$link_function$link %||% "custom")
      recorded$draws <- c(recorded$draws, list(components$link))
      components
    },
    .package = "RiskMap"
  )
  result <- suppressWarnings(assess_simulation(...))
  list(result = result, links = recorded$links, draws = recorded$draws)
}

test_that("a canonical logit assessment scores plogis of the linear predictors", {
  # The scores that target_transform = plogis gave before every model was
  # mapped through its own inverse link, computed here from the simulated
  # linear predictor and the predictive draws of the linear predictor.
  generating <- specify_glgpm(
    y ~ gp(), data = binomial_data, family = "binomial",
    denominator = denominator,
    parameters = list(beta = -0.5, sigma2 = 0.5, phi = 2)
  )
  sim <- simulate_glgpm(generating, nsim = 2, what = c("data", "surface"),
                        prediction_grid = binomial_data, seed = 3)
  run <- assess_recording_links(
    sim, models = list(logit = generating), control_mcmc = control_mcmc,
    seed = fixture_seed, spatial_scale = "location", objective = "mse",
    messages = FALSE
  )
  expect_identical(run$result$fits$status, c("ok", "ok"))
  expect_identical(run$links, c("logit", "logit"))
  truth <- plogis(sim$samples$surface[, , "linear_predictor"])
  expected <- vapply(1:2, function(j) {
    mean((rowMeans(plogis(run$draws[[j]])) - truth[, j])^2)
  }, numeric(1))
  expect_identical(unname(run$result$objective$location$mse[1, ]), expected)
})

test_that("assess_simulation scores every model through its own inverse link", {
  generating <- specify_glgpm(
    y ~ gp(), data = binomial_data, family = binomial("probit"),
    denominator = denominator,
    parameters = list(beta = -0.5, sigma2 = 0.3, phi = 2)
  )
  sim <- simulate_glgpm(generating, nsim = 2, what = c("data", "surface"),
                        prediction_grid = binomial_data, seed = 3)
  linear_predictor <- sim$samples$surface[, , "linear_predictor"]
  truth <- simulation_assessment_data(sim)$mean_sim
  expect_identical(truth, matrix(pnorm(linear_predictor),
                                 nrow(linear_predictor)))

  logit <- specify_glgpm(
    y ~ gp(), data = binomial_data, family = "binomial",
    denominator = denominator,
    parameters = list(beta = -0.8, sigma2 = 0.8, phi = 2)
  )
  targets <- list()
  record_target <- function(mean) {
    targets[[length(targets) + 1L]] <<- mean
    mean
  }
  run <- assess_recording_links(
    sim, models = list(probit = generating, logit = logit),
    control_mcmc = control_mcmc, seed = fixture_seed, spatial_scale = "location",
    target = record_target, objective = "mse", messages = FALSE
  )
  expect_identical(run$result$fits$status, rep("ok", 4))
  expect_identical(run$links, rep(c("probit", "logit"), 2))

  # The target receives the simulated prevalence first, then the predictive
  # prevalence of each fit through its own inverse link, and the scores
  # compare the two.
  expect_identical(targets[[1]], truth)
  inverse <- list(probit = function(eta) pnorm(eta), logit = plogis)
  for (fit in seq_along(run$draws)) {
    draws <- run$draws[[fit]]
    prevalence <- matrix(inverse[[run$links[fit]]](as.numeric(draws)),
                         nrow(draws))
    expect_identical(targets[[fit + 1L]], prevalence)
    model <- 2L - fit %% 2L
    simulation <- (fit + 1L) %/% 2L
    expect_identical(
      run$result$objective$location$mse[model, simulation],
      mean((rowMeans(prevalence) - truth[, simulation])^2)
    )
  }
})

test_that("assess_simulation accepts a generating model with a custom link", {
  imperfect_test <- function(eta) 0.01 + 0.89 * plogis(eta)
  generating <- specify_glgpm(
    y ~ gp(), data = binomial_data, family = "binomial",
    denominator = denominator, invlink = imperfect_test,
    parameters = list(beta = 0, sigma2 = 0.5, phi = 2)
  )
  sim <- simulate_glgpm(generating, nsim = 1, what = c("data", "surface"),
                        prediction_grid = binomial_data, seed = 3)
  linear_predictor <- sim$samples$surface[, , "linear_predictor"]
  truth <- simulation_assessment_data(sim)$mean_sim
  expect_identical(truth, matrix(imperfect_test(linear_predictor), ncol = 1))

  targets <- list()
  record_target <- function(mean) {
    targets[[length(targets) + 1L]] <<- mean
    mean
  }
  run <- assess_recording_links(
    sim, models = list(imperfect_test = generating),
    control_mcmc = control_mcmc, seed = fixture_seed, spatial_scale = "location",
    target = record_target, objective = "mse", messages = FALSE
  )
  expect_identical(run$result$fits$status, "ok")
  expect_identical(run$links, "custom")
  draws <- run$draws[[1]]
  expect_identical(targets[[2]],
                   matrix(imperfect_test(as.numeric(draws)), nrow(draws)))
  expect_true(is.finite(run$result$objective$location$mse[1, 1]))
})

test_that("assess_simulation no longer takes target_transform", {
  expect_error(
    assess_simulation(structure(list(), class = "RiskMap_simulation"),
                      models = list(model = gaussian_intercept_model),
                      spatial_scale = "location", target_transform = plogis),
    "unused argument \\(target_transform = plogis\\)"
  )
})

test_that("simulation objects print briefly", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)
  expect_output(print(sim), "Linear geostatistical model, y ~ gp().",
                fixed = TRUE)
  expect_output(print(sim), "2 simulations (seed 2)", fixed = TRUE)
  model <- specify_glgpm(
    y ~ cov + gp(), gaussian_data, family = "gaussian",
    parameters = list(beta = c(1, 0.5), sigma2 = 1, phi = 2, sigma2_me = 0.1)
  )
  expect_output(print(model), "beta +\\(Intercept\\) = 1, cov = 0.5")
  expect_output(print(model), "phi +2 km")
})

test_that("simulated surfaces can be plotted with fixed colour limits", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 1, what = "surface",
                        prediction_grid = grid, seed = 1)
  p <- plot(simulated_surface(sim), limits = c(-5, 5))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(p$scales$get_scales("fill")$limits, c(-5, 5))
  expect_no_error(plot(sim, limits = c(-5, 5)))
})

test_that("specify_glgpm combines several re() calls", {
  parameters <- list(beta = c(1, 0.5), sigma2 = 1, phi = 2,
                     sigma2_re = c(i = 0.3, j = 0.6), sigma2_me = 0.1)
  one_call <- specify_glgpm(y ~ cov + gp() + re(i, j), gaussian_data, "gaussian",
                            parameters)
  two_calls <- specify_glgpm(y ~ cov + gp() + re(i) + re(j), gaussian_data,
                             "gaussian", parameters)
  expect_identical(two_calls$re_terms, c("i", "j"))
  expect_identical(simulate_glgpm(two_calls, nsim = 2, seed = 1)$samples,
                   simulate_glgpm(one_call, nsim = 2, seed = 1)$samples)
})

test_that("simulations and assessments print counts in the singular and plural", {
  model <- list(family = "gaussian", formula = y ~ gp(),
                parameters = list(beta = c("(Intercept)" = 1), sigma2 = 1,
                                  phi = 2),
                distance_units = "km")
  sim <- structure(list(model = model, what = "surface",
                        locations = list(surface = data.frame(x = 1)),
                        nsim = 1, seed = 3),
                   class = "RiskMap_simulation")
  expect_match(squish(capture.output(print(sim))),
               "1 simulation (seed 3), of the surface at 1 location.",
               fixed = TRUE)
  assessment <- structure(list(fits = NULL, n_sim = 1, spatial_scale = "location",
                               objective = list(location = list(mse = NULL))),
                          class = "RiskMap_assess_simulation")
  expect_match(squish(capture.output(print(assessment))),
               paste("Models refitted to 1 simulation, assessed at the",
                     "locations by mean squared error."), fixed = TRUE)
  # A long row of coefficients wraps under its column.
  many <- model
  many$parameters$beta <- stats::setNames(seq(0.1, 1.2, by = 0.1),
                                          paste0("covariate_", 1:12))
  many$data <- data.frame(x = 1:10)
  class(many) <- "RiskMap_simulation_model"
  withr::local_options(cli.width = 40L)
  printed <- capture.output(print(many))
  expect_lte(max(nchar(printed)), 40L)
  expect_match(squish(printed), "covariate_12 = 1.2", fixed = TRUE)
})
