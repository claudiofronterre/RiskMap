test_that("predict_grid_target produces expected output with default arguments", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- setup_prediction(gaussian_model)
  gaussian_offset_grid <- setup_prediction(gaussian_offset_model)
  binomial_grid <- setup_prediction(binomial_model, control_sim = control_mcmc)
  poisson_grid <- setup_prediction(poisson_model, control_sim = control_mcmc)

  result <- predict_grid_target(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(gaussian_offset_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(poisson_grid)
  expect_setequal(names(result), expected_output)

})

test_that("predict_grid_target produces expected output when grid is provided", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- setup_prediction(gaussian_model,
                                  grid_pred = grid,
                                  predictors = data.frame(cov = rnorm(length(grid))),
                                  re_predictors = data.frame(i = 1:5),
                                  type = "joint")

  binomial_grid <- setup_prediction(binomial_model,
                                  grid_pred = grid,
                                  predictors = data.frame(cov = rnorm(length(grid))),
                                  re_predictors = data.frame(i = 1:5),
                                  control_sim = control_mcmc,
                                  type = "joint")

  poisson_grid <- setup_prediction(poisson_model,
                                 grid_pred = grid,
                                 predictors = data.frame(cov = rnorm(length(grid))),
                                 re_predictors = data.frame(i = 1:5),
                                 control_sim = control_mcmc,
                                 type = "joint")

  result <- predict_grid_target(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(poisson_grid)
  expect_setequal(names(result), expected_output)

})


test_that("predict_grid_target produces expected output when in list mode", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- setup_prediction(gaussian_model,
                                  grid_pred = list(grid, grid),
                                  predictors = list(data.frame(cov = rnorm(length(grid))),
                                                    data.frame(cov = rnorm(length(grid)))),
                                  type = "joint")

  binomial_grid <- setup_prediction(binomial_model,
                                  grid_pred = list(grid, grid),
                                  predictors = list(data.frame(cov = rnorm(length(grid))),
                                                    data.frame(cov = rnorm(length(grid)))),
                                  control_sim = control_mcmc,
                                  type = "joint")

  poisson_grid <- setup_prediction(poisson_model,
                                 grid_pred = list(grid, grid),
                                 predictors = list(data.frame(cov = rnorm(length(grid))),
                                                   data.frame(cov = rnorm(length(grid)))),
                                 control_sim = control_mcmc,
                                 type = "joint")

  result <- predict_grid_target(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- predict_grid_target(poisson_grid)
  expect_setequal(names(result), expected_output)
})

test_that("predict_grid_target handles one-pixel groups in list mode", {
  grid_pred <- list(
    group_one = sf::st_as_sf(data.frame(x = 0, y = 0), coords = c("x", "y"), crs = 4326),
    group_two = sf::st_as_sf(data.frame(x = c(1, 2), y = c(1, 2)), coords = c("x", "y"), crs = 4326)
  )

  object <- list(
    grid_pred = grid_pred,
    S_samples = list(
      matrix(c(1, 2, 3), nrow = 1),
      matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
    ),
    mu_pred = list(c(0), c(0, 0)),
    cov_offset = list(c(0), c(0, 0)),
    re = list(samples = list()),
    par_hat = list(),
    family = "gaussian"
  )
  class(object) <- "RiskMap_pred"

  out <- predict_grid_target(
    object,
    f_target = list(identity_target = function(x) x),
    pd_summary = list(mean = mean, sd = sd)
  )

  expect_equal(out$target$group_one$identity_target$mean, 2)
  expect_equal(out$target$group_one$identity_target$sd, sd(c(1, 2, 3)))
  expect_equal(out$target$group_two$identity_target$mean, c(2, 5))
  expect_equal(out$target$group_two$identity_target$sd, c(sd(c(1, 2, 3)),
                                                          sd(c(4, 5, 6))))
  expect_equal(dim(out$lp_samples[[1]]), c(1, 3))
  expect_equal(dim(out$lp_samples[[2]]), c(2, 3))
})

