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
})

test_that("specify_glgpm returns errors as expected", {
  expect_error(
    specify_glgpm("not a formula", gaussian_data, "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "'formula' must be a 'formula'"
  )

  expect_error(
    specify_glgpm(y ~ cov, gaussian_data, "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "The 'formula' must contain a Gaussian Process term"
  )

  expect_error(
    specify_glgpm(y ~ gp(), "not spatial data", "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "class 'sf'"
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian", model_crs = 4326,
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "'model_crs' must be a projected CRS"
  )

  expect_message(
    specify_glgpm(y ~ gp(), latlon_data, "gaussian",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "'data' are in longitude/latitude "
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "not-a-family",
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "should be one of"
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian",
                  denominator = cov,
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "When 'family' is 'gaussian', neither 'denominator'"
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian",
                  invlink = plogis,
                  list(beta = 1, sigma2 = 1, phi = 1, sigma2_me = 0)),
    "When 'family' is 'gaussian', neither 'denominator'"
  )

  expect_error(
    specify_glgpm(y ~ gp(), gaussian_data, "gaussian",
                  list(sigma2 = 1, phi = 1, sigma2_me = 0)),
    "Missing simulation parameters: beta"
  )

  expect_error(
    specify_glgpm(y ~ gp(), binomial_data, "binomial",
                  denominator = "d",
                  list(beta = 1, sigma2 = 1, phi = 1)),
                  "'denominator' must be provided as an unquoted column"
  )

  expect_error(
    specify_glgpm(y ~ gp(), binomial_data, "binomial",
                  denominator = d,
                  list(beta = 1, sigma2 = 1, phi = 1)),
    "'denominator' is not present in 'data'"
  )

  expect_error(
    specify_glgpm(y ~ gp(), binomial_data, "binomial",
                  denominator = denominator,
                  invlink = "a",
                  list(beta = 1, sigma2 = 1, phi = 1)),
    "'invlink' must be a function"
  )

})

test_that("specify_glgpm returns errors as expected", {

  expect_error(
    simulate_glgpm(gaussian_model, nsim = -1),
    "positive integer"
  )

  expect_error(
    simulate_glgpm(gaussian_model, nsim = 1, what = "n"),
    "'what' must be 'data', 'surface'"
  )

  expect_error(
    simulate_glgpm(gaussian_model, nsim = 1, what = "surface",
                   sample_locations = gaussian_data),
    "'sample_locations' requires 'data' in 'what'"
  )

  expect_error(
    simulate_glgpm(gaussian_model, nsim = 1, what = "data",
                   prediction_grid = grid),
    "'prediction_grid' requires 'surface' in 'what'"
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
    "transformed to the model CRS"
  )
  expect_identical(st_crs(transformed$locations$data), st_crs(gaussian_model$data))
  expect_error(simulate_glgpm(binomial_model, sample_locations = binomial_data[, "cov"]),
               "'i' is present in the model formula but not a column in 'sample_locations'")
  expect_error(simulated_data(first, 2), "valid simulation")
  expect_error(simulated_surface(first), "No surface")
  expect_error(simulate_glgpm(list(gaussian_model)), "fitted RiskMap")
})

test_that("missing formula variables in prediction_grid or sample_locations are named clearly (#204)", {
  surface_missing_cov <- grid
  expect_error(
    simulate_glgpm(gaussian_model, what = c("data", "surface"),
                   prediction_grid = surface_missing_cov,
                   sample_locations = gaussian_data),
    "'cov', 'i' are present in the model formula but not a column in 'prediction_grid'"
  )

  surface_missing_re <- grid
  surface_missing_re$cov <- rnorm(nrow(grid))
  expect_error(
    simulate_glgpm(gaussian_model, what = c("data", "surface"),
                   prediction_grid = surface_missing_re,
                   sample_locations = gaussian_data),
    "'i' is present in the model formula but not a column in 'prediction_grid'"
  )

  expect_error(
    simulate_glgpm(gaussian_model, sample_locations = gaussian_data[, "cov"]),
    "'i' is present in the model formula but not a column in 'sample_locations'"
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
  obj_sim <- structure(list(), class = "RiskMap_simulation")
  expect_error(
    assess_simulation(obj_sim, models = list(model = y ~ gp()),
                      spatial_scale = "grid", target_transform = identity),
    "glgpm\\(\\) or specify_glgpm\\(\\)"
  )
  expect_error(
    assess_simulation(obj_sim, models = list(model = gaussian_intercept_model),
                      spatial_scale = character(), target_transform = identity),
    "'spatial_scale' must be set"
  )
  expect_error(
    assess_simulation(obj_sim, models = list(model = gaussian_intercept_model),
                      spatial_scale = c("grid", "grid"), target_transform = identity),
    "'spatial_scale' must be set"
  )
  expect_error(
    assess_simulation(obj_sim, models = list(model = gaussian_intercept_model),
                      spatial_scale = "grid", pred_objective = character(),
                      target_transform = identity),
    "'pred_objective' must be either"
  )
  expect_error(assess_simulation(obj_sim, models = list(model = gaussian_intercept_model),
                                 spatial_scale = "area", boundaries = gaussian_data,
                                 area_summary = mean),
               "'boundaries' can only contain 'POLYGON' or 'MULTIPOLYGON' geometry")
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
  expect_null(args$fix_var_me)
  expect_false("start_pars" %in% names(args))

  sim <- simulate_glgpm(
    template, nsim = 1, what = c("data", "surface"),
    prediction_grid = gaussian_data, seed = 2
  )
  result <- assess_simulation(
    sim, models = list(candidate = template), spatial_scale = "grid",
    target_transform = identity, pred_objective = "mse", messages = FALSE
  )
  expect_identical(dim(result$pred_objective$grid$mse), c(1L, 1L))
})

test_that("joint output feeds the existing grid assessment", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)
  result <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                              control_mcmc = control_mcmc, spatial_scale = "grid",
                              target_transform = identity, pred_objective = "mse",
                              messages = FALSE)
  expect_identical(dim(result$pred_objective$grid$mse), c(1L, 2L))
  expect_true(all(is.finite(result$pred_objective$grid$mse)))
})

test_that("assess_simulation validates target and area function outputs #171", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)

  expect_error(
    assess_simulation(
      sim,
      models = list(intercept = gaussian_intercept_model),
      spatial_scale = "grid",
      target_transform = mean,
      pred_objective = "mse",
      messages = FALSE
    ),
    "same dimensions as its input"
  )

  expect_error(
    assess_simulation(
      sim,
      models = list(intercept = gaussian_intercept_model),
      spatial_scale = "grid",
      target_transform = function(x) x * NA_real_,
      pred_objective = "mse",
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
      target_transform = identity,
      area_summary = identity,
      boundaries = boundaries,
      pred_objective = "mse",
      messages = FALSE
    ),
    "one finite numeric value"
  )
})

test_that("assess_simulation computes grid and area objectives in one combined run #109", {
  boundaries <- create_convex_hull(gaussian_data)
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)

  combined <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                control_mcmc = control_mcmc,
                                spatial_scale = c("grid", "area"),
                                target_transform = identity, area_summary = mean,
                                boundaries = boundaries, pred_objective = "mse",
                                messages = FALSE)

  expect_setequal(names(combined$pred_objective), c("grid", "area"))
  expect_identical(dim(combined$pred_objective$grid$mse), c(1L, 2L))
  expect_identical(dim(combined$pred_objective$area$mse), c(1L, 2L))
  expect_true(all(is.finite(combined$pred_objective$grid$mse)))
  expect_true(all(is.finite(combined$pred_objective$area$mse)))

  grid_only <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                 control_mcmc = control_mcmc, spatial_scale = "grid",
                                 target_transform = identity, pred_objective = "mse",
                                 messages = FALSE)
  area_only <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                 control_mcmc = control_mcmc, spatial_scale = "area",
                                 target_transform = identity, area_summary = mean,
                                 boundaries = boundaries, pred_objective = "mse",
                                 messages = FALSE)

  # Area-level results are computed identically whether or not 'grid' is
  # also requested: both paths use the same (joint-type) prediction, so
  # under the shared control_mcmc seed the two runs must agree exactly.
  expect_identical(combined$pred_objective$area$mse, area_only$pred_objective$area$mse)

  # A standalone grid-only run instead uses cheaper marginal-type sampling,
  # since it doesn't need spatially correlated draws across the grid; a
  # combined run always uses joint-type sampling so one prediction can serve
  # both scales. The two are therefore consistent estimates of the same
  # quantity, not bitwise-identical draws.
  expect_equal(combined$pred_objective$grid$mse, grid_only$pred_objective$grid$mse,
               tolerance = 0.05)

  s <- summary(combined)
  expect_setequal(names(s), c("grid", "area", "fits"))
  expect_s3_class(s, "summary.RiskMap_assess_simulation")
  expect_output(print(s), "Grid-level results")
  expect_output(print(s), "Area-level results")

  combined_classify <- assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                                         control_mcmc = control_mcmc,
                                         spatial_scale = c("grid", "area"),
                                         target_transform = identity, area_summary = mean,
                                         boundaries = boundaries,
                                         pred_objective = c("mse", "classify"),
                                         categories = c(-3, -1, 0, 1, 3),
                                         messages = FALSE)

  expect_length(combined_classify$pred_objective$grid$classify$intercept$by_cat, 2L)
  expect_length(combined_classify$pred_objective$area$classify$intercept$by_cat, 2L)
  expected_classes <- c("(-3,-1]", "(-1,0]", "(0,1]", "(1,3]")
  expect_identical(
    levels(combined_classify$pred_objective$grid$classify$Class),
    expected_classes
  )

  expect_error(
    assess_simulation(sim, models = list(intercept = gaussian_intercept_model),
                      control_mcmc = control_mcmc,
                      spatial_scale = "grid", target_transform = identity,
                      pred_objective = "classify", categories = c(-1, 0, 0),
                      messages = FALSE),
    "unique, strictly increasing"
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
  two_simulations$pred_objective$grid$classify$model$by_cat <-
    list(one_result, second_result)
  two_simulations$pred_objective$grid$classify$model$CC <- c(0.65, 0.75)

  summary_two <- summary(two_simulations)$grid$classify$model
  expect_equal(summary_two$classify_res$Sensitivity, c(0.5, 0.25))
  expect_equal(summary_two$n_valid$Sensitivity, c(1L, 1L))
  expect_equal(summary_two$n_valid$Specificity, c(1L, 1L))
  expect_equal(summary_two$cc_summary$n_valid, 2L)
  expect_equal(summary_two$cc_summary$n_sim, 2L)
})

test_that("plot.RiskMap_simulation plots the simulated surface (#87)", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = grid, seed = 1)

  expect_s3_class(sim, "RiskMap_simulation")
  expect_error(plot(sim, simulation = 3), "must select valid simulation numbers")
  expect_no_error(plot(sim, simulation = 2))

  p <- plot(sim)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(p$labels$title, "Simulation 1")
  expect_equal(p$labels$fill, "linear_predictor")
  expect_equal(p$data$value,
               as.numeric(sim$samples$surface[, 1, "linear_predictor"]))
  # sample locations are overlaid when data were simulated
  expect_length(p$layers, 2)
  expect_equal(nrow(p$layers[[2]]$data), nrow(sim$locations$data))

  p2 <- plot(sim, simulation = 2)
  expect_equal(p2$labels$title, "Simulation 2")
  expect_equal(p2$data$value,
               as.numeric(sim$samples$surface[, 2, "linear_predictor"]))
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
  p <- plot(surface, palette = "Blues", reverse_palette = TRUE)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(p$data$value, surface$linear_predictor)
  expect_equal(p$labels$fill, "linear_predictor")
  expect_identical(ggplot2::ggplot_build(p)$data[[1]]$fill,
                   ggplot2::ggplot_build(plot(sim, simulation = 2, palette = "Blues",
                                              reverse_palette = TRUE))$data[[1]]$fill)
  expect_error(plot.RiskMap_simulated_surface(sim), "must be of class RiskMap_simulated_surface")
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
      control_mcmc = control_mcmc, spatial_scale = "grid",
      target_transform = identity, pred_objective = "mse", messages = FALSE
    ),
    "1 of 3 model fits failed"
  )

  expect_identical(result$fits$status, c("ok", "failed", "ok"))
  expect_identical(result$fits$message[2], "fit failed")
  expect_true(is.na(result$pred_objective$grid$mse[1, 2]))
  expect_true(all(is.finite(result$pred_objective$grid$mse[1, c(1, 3)])))
  expect_true(all(is.na(result$estimates$intercept[2, ])))
  expect_output(print(summary(result)),
                "Failed model fits \\(results missing\\): 1 of 3")
  expect_output(print(result), "Failed model fits: 1 of 3")
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
      control_mcmc = control_mcmc, spatial_scale = "area",
      target_transform = identity, area_summary = mean,
      boundaries = boundaries, col_names = "region",
      pred_objective = c("mse", "classify"), categories = c(-10, 0, 10),
      messages = FALSE
    ),
    "'far', which is left out of the area-level assessment"
  )
  expect_true(all(is.finite(result$pred_objective$area$mse)))
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
      spatial_scale = "area", target_transform = identity,
      area_summary = mean, boundaries = boundaries, col_names = "region",
      pred_objective = "mse", messages = FALSE
    ),
    "unique region names"
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
      sim, models = list(candidate = candidate), spatial_scale = "grid",
      target_transform = identity, pred_objective = "mse", messages = FALSE
    ),
    "needs complete values of 'cov' on the prediction grid"
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
    spatial_scale = "grid", target_transform = plogis, pred_objective = "mse",
    messages = FALSE
  ))
  expect_identical(result$fits$status, c("ok", "ok"))
  expect_true(all(is.finite(result$pred_objective$grid$mse)))
})

test_that("simulation objects print briefly", {
  sim <- simulate_glgpm(gaussian_intercept_model, nsim = 2,
                        what = c("data", "surface"),
                        prediction_grid = gaussian_data, seed = 2)
  expect_output(print(sim), "Simulations from a gaussian geostatistical model")
  expect_output(print(sim), "Number of simulations: 2 \\(seed 2\\)")
  model <- specify_glgpm(
    y ~ cov + gp(), gaussian_data, family = "gaussian",
    parameters = list(beta = c(1, 0.5), sigma2 = 1, phi = 2, sigma2_me = 0.1)
  )
  expect_output(print(model), "beta: \\(Intercept\\) = 1, cov = 0.5")
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
