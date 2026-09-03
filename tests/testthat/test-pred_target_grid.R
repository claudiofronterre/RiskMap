test_that("pred_target_grid produces expected output with default arguments", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- pred_over_grid(gaussian_model)
  gaussian_offset_grid <- pred_over_grid(gaussian_offset_model)
  binomial_grid <- pred_over_grid(binomial_model, control_sim = control_mcmc)
  poisson_grid <- pred_over_grid(poisson_model, control_sim = control_mcmc)

  result <- pred_target_grid(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- pred_target_grid(gaussian_offset_grid)
  expect_setequal(names(result), expected_output)

  result <- pred_target_grid(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- pred_target_grid(poisson_grid)
  expect_setequal(names(result), expected_output)

})

test_that("pred_target_grid produces expected output when grid is provided", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- pred_over_grid(gaussian_model,
                                  grid_pred = grid,
                                  predictors = data.frame(cov = rnorm(length(grid))),
                                  re_predictors = data.frame(i = 1:5),
                                  type = "joint")

  binomial_grid <- pred_over_grid(binomial_model,
                                  grid_pred = grid,
                                  predictors = data.frame(cov = rnorm(length(grid))),
                                  re_predictors = data.frame(i = 1:5),
                                  control_sim = control_mcmc,
                                  type = "joint")

  poisson_grid <- pred_over_grid(poisson_model,
                                 grid_pred = grid,
                                 predictors = data.frame(cov = rnorm(length(grid))),
                                 re_predictors = data.frame(i = 1:5),
                                 control_sim = control_mcmc,
                                 type = "joint")

  result <- pred_target_grid(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- pred_target_grid(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- pred_target_grid(poisson_grid)
  expect_setequal(names(result), expected_output)

})


test_that("pred_target_grid produces expected output when in list mode", {

  expected_output <- c("target", "grid_pred", "f_target", "pd_summary", "family", "lp_samples")

  gaussian_grid <- pred_over_grid(gaussian_model,
                                  grid_pred = list(grid, grid),
                                  predictors = list(data.frame(cov = rnorm(length(grid))),
                                                    data.frame(cov = rnorm(length(grid)))),
                                  type = "joint")

  binomial_grid <- pred_over_grid(binomial_model,
                                  grid_pred = list(grid, grid),
                                  predictors = list(data.frame(cov = rnorm(length(grid))),
                                                    data.frame(cov = rnorm(length(grid)))),
                                  control_sim = control_mcmc,
                                  type = "joint")

  poisson_grid <- pred_over_grid(poisson_model,
                                 grid_pred = list(grid, grid),
                                 predictors = list(data.frame(cov = rnorm(length(grid))),
                                                   data.frame(cov = rnorm(length(grid)))),
                                 control_sim = control_mcmc,
                                 type = "joint")

  result <- pred_target_grid(gaussian_grid)
  expect_setequal(names(result), expected_output)

  result <- pred_target_grid(binomial_grid)
  expect_setequal(names(result), expected_output)

  result <- pred_target_grid(poisson_grid)
  expect_setequal(names(result), expected_output)
})

test_that("pred_target_grid handles one-pixel groups in list mode", {
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

  out <- pred_target_grid(
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

