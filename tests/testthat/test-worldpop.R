fake_raster <- function(path, value) {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
                   ymax = 2, crs = "EPSG:4326")
  terra::values(r) <- value
  terra::writeRaster(r, path, overwrite = TRUE)
  path
}

fake_download <- function(urls, dest_dir, messages) {
  paths <- file.path(dest_dir, basename(urls))
  for (i in seq_along(paths)) fake_raster(paths[i], i)
  paths
}

age_catalogue <- function(dataset, iso3) {
  groups <- c(0, 1, seq(5, 80, by = 5))
  files <- sprintf("https://example.org/lbr_%s_%d_2011.tif",
                   rep(c("f", "m"), each = length(groups)), groups)
  data.frame(popyear = c("2010", "2011"),
             files = I(list("https://example.org/old.tif", files)))
}

test_that("get_worldpop selects the age groups and sums them", {
  local_mocked_bindings(worldpop_catalogue = age_catalogue,
                        worldpop_download = fake_download)
  adults <- get_worldpop("lbr", 2011, ages = c(20, Inf),
                         dest_dir = withr::local_tempdir(), messages = FALSE)
  expect_s4_class(adults, "SpatRaster")
  expect_identical(names(adults), "population")
  # 13 groups (20 to 80+) for each of two sexes, with values 1 to 26
  expect_equal(unique(terra::values(adults)[, 1]), sum(1:26))

  women <- get_worldpop("LBR", 2011, ages = c(15, 49), sex = "female",
                        dest_dir = withr::local_tempdir(), messages = FALSE)
  expect_equal(unique(terra::values(women)[, 1]), sum(1:7))

  expect_error(get_worldpop("LBR", 2011, ages = c(18, Inf),
                            messages = FALSE),
               "must start and end on WorldPop age groups")
  expect_error(get_worldpop("LBR", 2012, ages = c(20, Inf),
                            messages = FALSE),
               "no unconstrained")
})

test_that("get_worldpop chooses products and checks arguments", {
  product <- RiskMap:::worldpop_product
  expect_identical(product(2014, NULL, "100m", FALSE, FALSE)$dataset,
                   "pop/wpgp")
  expect_identical(product(2014, NULL, "1km", FALSE, TRUE)$dataset,
                   "pop/wpicuadj1km")
  expect_identical(product(2011, c(20, Inf), "100m", FALSE, FALSE)$dataset,
                   "age_structures/aswpgp")
  expect_true(is.na(product(2011, c(20, Inf), "1km", FALSE, FALSE)$dataset))
  expect_identical(product(2020, c(0, 4), "1km", TRUE, FALSE)$dataset,
                   "age_structures/G2_CN_Age_R25A_1km")
  expect_error(product(2011, NULL, "100m", TRUE, FALSE), "2015 to 2030")
  expect_error(product(2025, NULL, "100m", FALSE, FALSE), "2000 to 2020")
  expect_error(product(2011, c(20, Inf), "100m", FALSE, TRUE),
               "whole population")

  urls <- RiskMap:::worldpop_files(product(2011, c(20, Inf), "1km", FALSE,
                                           FALSE),
                                   "LBR", 2011, c(20, Inf), "male")
  expect_length(urls, 13L)
  expect_match(urls[1], "Global_2000_2020_1km/unconstrained/2011/LBR/lbr_m_20_2011_1km.tif",
               fixed = TRUE)

  expect_error(get_worldpop("Liberia", 2011), "ISO3")
  expect_error(get_worldpop("LBR", 2011.5), "whole number")
  expect_error(get_worldpop("LBR", 2011, ages = 20), "c\\(from, to\\)")
})

test_that("get_worldpop uses the both-sexes layers when available", {
  catalogue <- function(dataset, iso3) {
    files <- sprintf("https://example.org/lbr_%s_%02d_2015_CN_1km.tif",
                     rep(c("f", "m", "t"), each = 3), rep(c(0, 1, 5), 3))
    data.frame(popyear = "2015", files = I(list(files)))
  }
  local_mocked_bindings(worldpop_catalogue = catalogue)
  urls <- RiskMap:::worldpop_files(
    list(dataset = "age_structures/G2_CN_Age_R25A_1km"), "LBR", 2015,
    c(0, 4), "both"
  )
  expect_identical(basename(urls),
                   c("lbr_t_00_2015_CN_1km.tif", "lbr_t_01_2015_CN_1km.tif"))
})

test_that("get_worldpop downloads from WorldPop", {
  skip_on_cran()
  skip_if_offline("hub.worldpop.org")
  total <- get_worldpop("LBR", 2011, resolution = "1km",
                        dest_dir = withr::local_tempdir(), messages = FALSE)
  expect_gt(terra::global(total, "sum", na.rm = TRUE)[[1]], 3e6)
})
