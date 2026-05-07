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
  class(object) <- "RiskMap.pred.re"

  out <- pred_target_grid(
    object,
    f_target = list(identity_target = function(x) x),
    pd_summary = list(mean = mean, sd = sd)
  )

  expect_equal(out$target$group_one$identity_target$mean, 2)
  expect_equal(out$target$group_one$identity_target$sd, stats::sd(c(1, 2, 3)))
  expect_equal(out$target$group_two$identity_target$mean, c(2, 5))
  expect_equal(out$target$group_two$identity_target$sd, c(stats::sd(c(1, 2, 3)),
                                                          stats::sd(c(4, 5, 6))))
})
