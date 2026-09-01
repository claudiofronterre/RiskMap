hull <- convex_hull_sf(gaussian_data)
grid <- create_grid(hull, 3, propose_utm(hull))

test_that("pred_over_grid produces errors as expected", {

  expect_error(
    pred_over_grid("not model"),
    "'object' must be of class RiskMap")

  expect_error(
    pred_over_grid(gaussian_model, grid_pred = "not sf"),
    "'grid_pred' must be of class 'sf' or 'sfc'")

  expect_error(
    pred_over_grid(gaussian_model, grid_pred = list("not sf"), type = "joint"),
    "Each element of 'grid_pred' must be an 'sf' or 'sfc'")

  expect_error(
    pred_over_grid(gaussian_model, grid_pred = list(gaussian_data), type = "marginal"),
    "When 'grid_pred' is a list, 'type' must be 'joint'")

  expect_error(
    pred_over_grid(gaussian_model, control_sim = "not mcmc"),
    "'control_sim' must be an output from 'set_control_sim")

  expect_error(
    pred_over_grid(gaussian_model, type = "not type"),
    "'type' must be either 'marginal' or 'joint'")

  expect_error(
    pred_over_grid(gaussian_offset_model, grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid)))),
    "'pred_cov_offset' must be specified at each prediction location"
  )

  expect_error(
    pred_over_grid(gaussian_offset_model, grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid))),
                   pred_cov_offset = "not n"),
    "'pred_cov_offset' must be a numeric vector"
  )

  expect_error(
    pred_over_grid(gaussian_offset_model, grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid))),
                   pred_cov_offset = rnorm(length(grid) - 1)),
    "The length of 'pred_cov_offset' does"
  )

  expect_error(
    pred_over_grid(gaussian_model, grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid))),
                   re_predictors = 1),
    "Random effect predictions require 'type' to be set to 'joint'"
  )

  expect_error(
    pred_over_grid(gaussian_model, grid_pred = grid),
    "'predictors' must be supplied if 'grid_pred' is supplied"
  )

  expect_error(
    pred_over_grid(gaussian_model, grid_pred = grid, predictors = "not df"),
    "'predictors' must be a data.frame"
  )

  expect_error(
    pred_over_grid(gaussian_model, grid, predictors = data.frame(cov = rnorm(length(grid) - 1))),
    "The number of rows in 'predictors' does not match the number of locations in 'grid_pred'"
  )

  expect_error(
    pred_over_grid(gaussian_model, grid, predictors = data.frame(not_cov = rnorm(length(grid)))),
    "The column names in 'predictors' do not match the variables in the model formula"
  )

  # expect_error(
  #   pred_over_grid(gaussian_model,
  #                  list(grid, grid),
  #                  type = "joint",
  #                  predictors = list(data.frame(not_cov = rnorm(length(grid))),
  #                                    data.frame(not_cov = rnorm(length(grid))))
  #                  ),
  #   "The column names in 'predictors' do not match the variables in the model formula"
  # )
  #
  # expect_no_error(
  #   pred_over_grid(gaussian_model,
  #                  list(grid, grid),
  #                  type = "joint",
  #                  predictors = list(data.frame(cov = rnorm(length(grid))),
  #                                    data.frame(cov = rnorm(length(grid))))
  #   )
  # )

})

test_that("pred_over_grid generates warnings", {
  old_opt <- options(warn = 2)  # Convert warnings to errors to force stop
  on.exit(options(old_opt))     # Restore original setting

  expect_error(
    pred_over_grid(gaussian_model,
                   grid_pred = grid,
                   predictors = data.frame(not_cov = rnorm(length(grid))),
                   pred_cov_offset = 1:n),
    "You have set 'pred_cov_offset' but 'object'"
  )

  expect_error(
    pred_over_grid(gaussian_model, predictors = "not df"),
    "You have set 'predictors' but not 'grid_pred' so 'predictors' will be ignored"
  )

})

test_that("pred_over_grid produces expected output", {

  expected_output <- c("mu_pred", "grid_pred", "par_hat", "S_samples", "re", "obs_loc",
                       "ID_coords", "inter_f", "family", "cov_offset", "type")

  result <- pred_over_grid(gaussian_model)
  expect_setequal(names(result), expected_output)

  result <- pred_over_grid(binomial_model, control_sim = control_mcmc)
  expect_setequal(names(result), expected_output)

  result <- pred_over_grid(gaussian_offset_model)
  expect_setequal(names(result), expected_output)

  result <- pred_over_grid(gaussian_model,
                           grid_pred = grid,
                           predictors = data.frame(cov = rnorm(length(grid))))
  expect_setequal(names(result), expected_output)

  result <- pred_over_grid(gaussian_model,
                           grid_pred = grid,
                           predictors = data.frame(cov = rnorm(length(grid))),
                           re_predictors = data.frame(i = 1:5),
                           type = "joint")
  expect_setequal(names(result), expected_output)

})

test_that("pred_over_grid functions for poisson models and in list mode", {
  skip()

  result <- pred_over_grid(poisson_model, control_sim = control_mcmc)
  expect_setequal(names(result), expected_output)

  # Error in C %*% B : non-conformable arguments
  # C is a list of matrices
  result <- pred_over_grid(gaussian_model,
                           grid_pred = list(grid, grid),
                           predictors = list(data.frame(cov = rnorm(length(grid))),
                                             data.frame(cov = rnorm(length(grid)))),
                           # cannot handle these as lists
                           # re_predictors = list(data.frame(i = 1:5),
                           #                      data.frame(i = 1:5)),
                           type = "joint")
  expect_setequal(names(result), expected_output)

  pred_target_grid(result)

  # fails as n_pred is a vector
  result <- pred_over_grid(gaussian_offset_model, list(grid, grid), pred_cov_offset = rnorm(length(grid)), type = "joint")

  # mu_pred is zero?
  result <- pred_over_grid(no_offset, grid)


})
