test_that("setup_prediction produces errors as expected", {

  expect_error(
    setup_prediction("not model"),
    "'object' must be of class RiskMap")

  expect_error(
    setup_prediction(gaussian_model,
                   grid_pred = "not sf"),
    "'grid_pred' must be of class 'sf' or 'sfc'")

  expect_error(
    setup_prediction(gaussian_model,
                   grid_pred = list("not sf"),
                   type = "joint"),
    "Each element of 'grid_pred' must be an 'sf' or 'sfc'")

  expect_error(
    setup_prediction(gaussian_model,
                   grid_pred = list(gaussian_data),
                   type = "marginal"),
    "When 'grid_pred' is a list, 'type' must be 'joint'")

  expect_error(
    setup_prediction(gaussian_model,
                   control_sim = "not mcmc"),
    "'control_sim' must be an output from 'set_control_sim")

  expect_error(
    setup_prediction(gaussian_model,
                   type = "not type"),
    "'type' must be either 'marginal' or 'joint'")

  expect_error(
    setup_prediction(gaussian_offset_model,
                   grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid)))),
    "'pred_cov_offset' must be specified at each prediction location"
  )

  expect_error(
    setup_prediction(gaussian_offset_model,
                   grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid))),
                   pred_cov_offset = "not n"),
    "'pred_cov_offset' must be a numeric vector"
  )

  expect_error(
    setup_prediction(gaussian_offset_model,
                   grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid))),
                   pred_cov_offset = rnorm(length(grid) - 1)),
    "The length of 'pred_cov_offset' does"
  )

  expect_error(
    setup_prediction(gaussian_offset_model,
                   grid_pred = list(grid, grid),
                   predictors = list(data.frame(cov = rnorm(length(grid))),
                                     data.frame(cov = rnorm(length(grid)))),
                   pred_cov_offset = rnorm(length(grid) - 1),
                   type = "joint"),
    "Predictions including covariate offsets"
  )

  expect_error(
    setup_prediction(gaussian_model,
                   grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid))),
                   re_predictors = 1),
    "Random effect predictions require 'type' to be set to 'joint'"
  )

  expect_error(
    setup_prediction(gaussian_model,
                   grid_pred = grid),
    "'predictors' must be supplied if 'grid_pred' is supplied"
  )

  expect_error(
    setup_prediction(gaussian_model,
                   grid_pred = grid,
                   predictors = "not df"),
    "'predictors' must be a data.frame"
  )

  expect_error(
    setup_prediction(gaussian_model,
                   grid,
                   predictors = data.frame(cov = rnorm(length(grid) - 1))),
    "The number of rows in 'predictors' does not match the number of locations in 'grid_pred'"
  )

  expect_error(
    setup_prediction(gaussian_model,
                   grid,
                   predictors = data.frame(not_cov = rnorm(length(grid)))),
    "The column names in 'predictors' do not match the variables in the model formula"
  )

  expect_error(
    setup_prediction(gaussian_model,
                   grid_pred = list(grid, grid),
                   predictors = list(data.frame(cov = rnorm(length(grid))),
                                     data.frame(cov = rnorm(length(grid)))),
                   re_predictors = list(data.frame(i = 1:5),
                                        data.frame(i = 1:5)),
                   type = "joint"),
    "Prediction of random effects with a list of prediction"
  )

  expect_error(
    setup_prediction(gaussian_model,
                   list(grid, grid),
                   type = "joint",
                   predictors = list(data.frame(not_cov = rnorm(length(grid))),
                                     data.frame(not_cov = rnorm(length(grid))))
                   ),
    "The column names in 'predictors' do not match the variables in the model formula"
  )

})

test_that("setup_prediction generates warnings", {
  old_opt <- options(warn = 2)  # Convert warnings to errors to force stop
  on.exit(options(old_opt))     # Restore original setting

  expect_error(
    setup_prediction(gaussian_model,
                   grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid))),
                   pred_cov_offset = 1:n),
    "You have set 'pred_cov_offset' but 'object'"
  )

  expect_error(
    setup_prediction(gaussian_model, predictors = "not df"),
    "You have set 'predictors' but not 'grid_pred' so 'predictors' will be ignored"
  )

})

test_that("setup_prediction produces expected output", {

  expected_output <- c("mu_pred", "grid_pred", "par_hat", "S_samples", "re", "obs_loc",
                       "ID_coords", "inter_f", "family", "cov_offset", "type")

  result <- setup_prediction(gaussian_model)
  expect_setequal(names(result), expected_output)

  result <- setup_prediction(binomial_model, control_sim = control_mcmc)
  expect_setequal(names(result), expected_output)

  result <- setup_prediction(poisson_model, control_sim = control_mcmc)
  expect_setequal(names(result), expected_output)

  result <- setup_prediction(gaussian_offset_model)
  expect_setequal(names(result), expected_output)

  result <- setup_prediction(gaussian_model,
                           grid_pred = grid,
                           predictors = data.frame(cov = rnorm(length(grid))))
  expect_setequal(names(result), expected_output)

  result <- setup_prediction(gaussian_model,
                           grid_pred = grid,
                           predictors = data.frame(cov = rnorm(length(grid))),
                           re_predictors = data.frame(i = 1:5),
                           type = "joint")
  expect_setequal(names(result), expected_output)

  result <- setup_prediction(gaussian_model,
                           grid_pred = list(grid, grid),
                           predictors = list(data.frame(cov = rnorm(length(grid))),
                                             data.frame(cov = rnorm(length(grid)))),
                           type = "joint")
  expect_setequal(names(result), expected_output)

  result <- setup_prediction(binomial_model,
                           grid_pred = list(grid, grid),
                           predictors = list(data.frame(cov = rnorm(length(grid))),
                                             data.frame(cov = rnorm(length(grid)))),
                           control_sim = control_mcmc,
                           type = "joint")
  expect_setequal(names(result), expected_output)

  result <- setup_prediction(poisson_model,
                           grid_pred = list(grid, grid),
                           predictors = list(data.frame(cov = rnorm(length(grid))),
                                             data.frame(cov = rnorm(length(grid)))),
                           control_sim = control_mcmc,
                           type = "joint")
  expect_setequal(names(result), expected_output)



})

