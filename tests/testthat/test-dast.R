test_that("dast validates intervention matrix dimensions", {
  df <- data.frame(
    y = c(1, 2, 3),
    m = c(10, 10, 10),
    time = c(1, 2, 3),
    x = c(0, 1, 2),
    z = c(0, 1, 0)
  )

  data <- st_as_sf(df, coords = c("x", "z"), crs = 4326)

  expect_error(
    dast(
      y ~ gp(),
      data = data,
      den = m,
      time = time,
      mda_times = c(0.5, 1.5),
      int_mat = matrix(0, nrow = 2, ncol = 2),
      power_val = 1
    ),
    "'int_mat' must have 3 rows",
    fixed = TRUE
  )

  expect_error(
    dast(
      y ~ gp(),
      data = data,
      den = m,
      time = time,
      mda_times = c(0.5, 1.5),
      int_mat = matrix(0, nrow = 4, ncol = 2),
      power_val = 1
    ),
    "'int_mat' must have 3 rows",
    fixed = TRUE
  )

  expect_error(
    dast(
      y ~ gp(),
      data = data,
      den = m,
      time = time,
      mda_times = c(0.5, 1.5),
      int_mat = matrix(0, nrow = 3, ncol = 1),
      power_val = 1
    ),
    "'int_mat' must have 2 columns",
    fixed = TRUE
  )
})
