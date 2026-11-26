# Each plot_*() function is the plot() method of its class under its own name.
# plot() and the function give equal objects, the function refuses an object
# of another class, and an argument the function doesn't have stops.

plot_cases <- function() {
  variogram_result <- variogram(gaussian_data, variable = "y", n_bins = 3L,
                                n_permutations = 19L, seed = 1)
  prediction <- predict(gaussian_intercept_model, newdata = grid,
                        seed = 1, messages = FALSE)
  extrapolation <- assess_extrapolation(gaussian_intercept_model, grid,
                                        spatial_threshold = 0.8)
  checked <- check_computation(binomial_model, messages = FALSE)
  profiles <- profile(gaussian_intercept_model, parameters = "phi",
                      n_points = 5L, seed = 1, messages = FALSE)
  bootstrap <- parametric_bootstrap(gaussian_intercept_model, n_sim = 3,
                                    seed = 1, messages = FALSE)
  list(
    variogram = list(fun = plot_variogram, x = variogram_result,
                     args = list(show_pairs = TRUE), made_by = "variogram"),
    prediction = list(fun = plot_prediction, x = prediction,
                      args = list(summaries = "upper"), made_by = "predict"),
    extrapolation = list(fun = plot_extrapolation, x = extrapolation,
                         args = list(indicator = "spatial_unexplained"),
                         made_by = "assess_extrapolation"),
    check = list(fun = plot_computation, x = checked,
                 args = list(which = "certificate"), made_by = "check_computation"),
    profile = list(fun = plot_profile, x = profiles,
                   args = list(scale = "deviance"), made_by = "profile"),
    bootstrap = list(fun = plot_bootstrap, x = bootstrap,
                     args = list(parameters = "phi", log_scale = TRUE),
                     made_by = "parametric_bootstrap"),
    predictive_check = list(
      fun = plot_predictive_check, made_by = "assess_fit",
      x = assess_fit(gaussian_intercept_model, nsim = 19, seed = 1),
      args = list(type = "variogram", colour = "darkblue")
    ),
    cross_validation = list(fun = plot_cross_validation, x = cross_validation,
                            args = list(type = "score", score = "SCRPS"),
                            made_by = "assess_prediction"),
    test_sets = list(fun = plot_test_sets, x = cross_validation$test_sets,
                     args = list(type = "distance"),
                     made_by = "create_test_sets"),
    simulation = list(fun = plot_simulation, x = simulation,
                      args = list(quantity = "spatial_effect"),
                      made_by = "simulate_glgpm"),
    simulated_surface = list(fun = plot_simulated_surface,
                             x = simulated_surface(simulation),
                             args = list(quantity = "linear_predictor"),
                             made_by = "simulated_surface")
  )
}

cases <- plot_cases()

test_that("plot() and the plot_*() function of a class give the same plot", {
  skip_on_cran()
  for (case in cases) {
    expect_equal(case$fun(case$x), plot(case$x))
    expect_equal(do.call(case$fun, c(list(case$x), case$args)),
                 do.call(plot, c(list(case$x), case$args)))
    # The arguments change the plot, so the comparison above is not empty.
    # A plot keeps the frames it was built in, which hold the arguments, so
    # these are left out of the comparison.
    expect_failure(expect_equal(case$fun(case$x),
                                do.call(case$fun, c(list(case$x), case$args)),
                                ignore_function_env = TRUE,
                                ignore_formula_env = TRUE))
  }
})

test_that("each plot_*() function refuses the result of another function", {
  for (name in names(cases)) {
    case <- cases[[name]]
    wrong <- if (name == "variogram") cases$prediction$x else cases$variogram$x
    error <- expect_input_error(case$fun(wrong), "x")
    expect_match(conditionMessage(error),
                 paste0("needs the result of `", case$made_by, "()`"),
                 fixed = TRUE)
    expect_input_error(case$fun(), "x")
  }
})

test_that("plot_*() and plot() stop on an argument they don't use", {
  # The maps pass their other arguments to the fill scale, which refuses them.
  checked <- setdiff(names(cases),
                     c("prediction", "simulation", "simulated_surface"))
  for (name in checked) {
    case <- cases[[name]]
    expect_input_error(case$fun(case$x, not_an_argument = 1),
                       "not_an_argument")
    expect_input_error(plot(case$x, not_an_argument = 1), "not_an_argument")
  }
  for (name in c("prediction", "simulation", "simulated_surface")) {
    case <- cases[[name]]
    expect_error(case$fun(case$x, not_an_argument = 1), "not_an_argument")
    expect_error(plot(case$x, not_an_argument = 1), "not_an_argument")
  }
})

test_that("plot() of a variogram or a prediction names the replacements of former arguments", {
  variogram_result <- cases$variogram$x
  expect_match(conditionMessage(expect_input_error(
    plot(variogram_result, plot_envelope = FALSE), "plot_envelope"
  )), "envelope")
  expect_match(conditionMessage(expect_input_error(
    plot(variogram_result, color = "red"), "color"
  )), "colour")
  expect_input_error(plot(cases$prediction$x, summary = "mean"), "summary")
})
