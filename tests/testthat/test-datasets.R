test_that("tz_covariates has the covariates of tz_malaria on the same scales", {
  data("tz_malaria", package = "RiskMap", envir = environment())
  data("tz_covariates", package = "RiskMap", envir = environment())
  for (name in names(sf::st_drop_geometry(tz_covariates))) {
    expect_false(any(tz_covariates[[name]] == -9999, na.rm = TRUE))
  }
  expect_true(all(tz_covariates$ITN >= 0 & tz_covariates$ITN <= 100,
                  na.rm = TRUE))
  # Each cluster against the nearest point of the grid, which for most
  # clusters is a few kilometres away. Values on one scale differ little at
  # the median.
  nearest <- sf::st_nearest_feature(tz_malaria, tz_covariates)
  for (name in c("EVI", "Temperature", "Precipitation", "Population", "ITN",
                 "NTL")) {
    cluster <- tz_malaria[[name]]
    cell <- tz_covariates[[name]][nearest]
    kept <- is.finite(cluster) & is.finite(cell)
    cluster <- cluster[kept]
    cell <- cell[kept]
    expect_lt(abs(stats::median(cell - cluster)),
              0.25 * stats::IQR(cluster), label = name)
    expect_gt(stats::median(cell) / stats::median(cluster), 0.5,
              label = name)
    expect_lt(stats::median(cell) / stats::median(cluster), 2, label = name)
  }
})
