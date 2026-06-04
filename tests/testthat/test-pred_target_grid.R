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
  expect_equal(dim(out$lp_samples[[1]]), c(1, 3))
  expect_equal(dim(out$lp_samples[[2]]), c(2, 3))
})

test_that("pred_target_shp preserves posterior samples for one-pixel list-mode regions", {
  grid_pred <- list(
    group_one = sf::st_as_sf(data.frame(x = 0, y = 0), coords = c("x", "y"), crs = 4326),
    group_two = sf::st_as_sf(data.frame(x = c(1, 2), y = c(1, 2)), coords = c("x", "y"), crs = 4326)
  )
  shp <- sf::st_sf(
    region = c("group_one", "group_two"),
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 4326)
  )
  object <- list(
    type = "joint",
    grid_pred = grid_pred,
    S_samples = list(
      matrix(c(0.01, 0.02, 0.03), nrow = 1),
      matrix(c(0.10, 0.10, 0.10, 0.40, 0.40, 0.40), nrow = 2, byrow = TRUE)
    ),
    mu_pred = list(c(0), c(0, 0)),
    cov_offset = list(c(0), c(0, 0)),
    re = list(samples = list()),
    par_hat = list()
  )
  class(object) <- "RiskMap.pred.re"

  out <- pred_target_shp(
    object,
    shp = shp,
    shp_target = sum,
    weights = list(1, c(0.25, 0.75)),
    standardize_weights = FALSE,
    col_names = "region",
    f_target = list(identity_target = identity),
    pd_summary = list(mean = mean),
    return_shp = FALSE,
    return_target_samples = TRUE,
    messages = FALSE
  )

  expect_equal(out$target_samples$group_one$identity_target, c(0.01, 0.02, 0.03))
  expect_equal(out$target$group_one$identity_target$mean, 0.02)
  expect_equal(out$target_samples$group_two$identity_target, rep(0.325, 3))
})

test_that("pred_target_shp errors on wrong list-mode target orientation", {
  grid_pred <- list(
    group_one = sf::st_as_sf(data.frame(x = c(0, 1), y = c(0, 1)), coords = c("x", "y"), crs = 4326)
  )
  shp <- sf::st_sf(
    region = "group_one",
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326)
  )
  object <- list(
    type = "joint",
    grid_pred = grid_pred,
    S_samples = list(matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06), nrow = 2)),
    mu_pred = list(c(0, 0)),
    cov_offset = list(c(0, 0)),
    re = list(samples = list()),
    par_hat = list()
  )
  class(object) <- "RiskMap.pred.re"

  expect_error(
    pred_target_shp(
      object,
      shp = shp,
      weights = list(c(0.5, 0.5)),
      col_names = "region",
      f_target = list(bad_target = function(x) t(x)),
      pd_summary = list(mean = mean),
      return_shp = FALSE,
      messages = FALSE
    ),
    "expected a 2 x 3 matrix"
  )
})
