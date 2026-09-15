test_that("simulate_surface accepts a valid custom CRS without an EPSG code", {
  custom_crs <- st_crs(
    "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
  )
  pred_grid <- st_as_sf(
    data.frame(
      x = c(500000, 501000, 500000, 501000),
      y = c(1000000, 1000000, 1001000, 1001000),
      cov = 0
    ),
    coords = c("x", "y"),
    crs = custom_crs
  )
  sample_data <- pred_grid[1:2, ]
  sample_data$units_m <- 1
  sampling_f <- function() sample_data

  expect_true(is.na(st_crs(pred_grid)$epsg))
  result <- simulate_surface(
    n_sim = 1,
    pred_grid = pred_grid,
    formula = y ~ cov + gp(),
    sampling_f = sampling_f,
    family = "binomial",
    scale_to_km = FALSE,
    par0 = list(beta = c(0, 0), sigma2 = 1, phi = 1000),
    messages = FALSE
  )

  expect_s3_class(result, "RiskMap_simulation")
  expect_true(isTRUE(result$crs == custom_crs))
  expect_true(isTRUE(st_crs(result$lp_grid_sim) == custom_crs))
})
