lonlat_points <- function(lon, lat) {
  st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"),
           crs = 4326)
}

lonlat_box <- function(west, east, south, north) {
  ring <- rbind(c(west, south), c(east, south), c(east, north),
                c(west, north), c(west, south))
  st_sf(geometry = st_sfc(st_polygon(list(ring)), crs = 4326))
}

# Corners and centre of the Democratic Republic of the Congo's extent, too
# wide for a UTM zone to keep distances within 0.5%
wide_points <- lonlat_points(lon = c(12, 31, 12, 31, 21.5),
                             lat = c(-13, -13, 5, 5, -4))

test_that("propose_crs checks its arguments", {
  expect_input_error(propose_crs("a"), "data")
  expect_input_error(propose_crs(st_set_crs(gaussian_data, NA)), "data")
  expect_input_error(propose_crs(gaussian_data, epsg_tolerance = -0.1),
                     "epsg_tolerance")
  expect_input_error(propose_crs(gaussian_data, epsg_tolerance = 1),
                     "epsg_tolerance")
  expect_input_error(propose_crs(gaussian_data, messages = NA), "messages")
  expect_input_error(
    propose_crs(st_sf(geometry = st_sfc(st_point(), crs = 4326))),
    "data"
  )
})

test_that("propose_crs gives the UTM zone of a single location", {
  expect_message(epsg <- propose_crs(lonlat_points(33, 5)),
                 "EPSG:32636", class = "riskmap_message")
  expect_identical(epsg, 32636L)

  # identical locations are a single location
  epsg <- propose_crs(lonlat_points(c(33, 33, 33), c(5, 5, 5)),
                      messages = FALSE)
  expect_identical(epsg, 32636L)
  distortion <- crs_distortion(lonlat_points(33, 5), epsg)
  expect_lt(distortion$max_error_pct, 0.05)
})

test_that("propose_crs takes the zone at the centre, not the modal zone", {
  # most locations are in zone 36, but the centre of the region is in zone 37
  points <- lonlat_points(lon = c(35.5, 35.6, 35.7, 35.8, 40), lat = 5)
  expect_identical(propose_crs(points, messages = FALSE), 32637L)
})

test_that("propose_crs counts the Equator as north", {
  expect_identical(propose_crs(lonlat_points(c(33, 34), 0), messages = FALSE),
                   32636L)
  south <- propose_crs(lonlat_points(c(33, 34), c(0.1, -0.3)),
                       messages = FALSE)
  expect_identical(south, 32736L)
  # the hemisphere does not change distances
  points <- lonlat_points(c(33, 34), c(0.1, -0.3))
  expect_equal(crs_distortion(points, 32636), crs_distortion(points, 32736),
               tolerance = 1e-4)
})

test_that("propose_crs handles data across the antimeridian", {
  near <- lonlat_points(c(179.6, -179.6), c(-17, -17.5))
  epsg <- propose_crs(near, messages = FALSE)
  expect_true(epsg %in% c(32701L, 32760L))
  expect_lt(crs_distortion(near, epsg)$max_error_pct, 0.2)

  wide <- lonlat_points(c(166, 174, -175, -166, 170, -170),
                        c(-12, -24, -14, -22, -18, -18))
  crs <- propose_crs(wide, messages = FALSE)
  expect_s3_class(crs, "crs")
  expect_lt(crs_distortion(wide, crs)$max_error_pct, 0.5)
})

test_that("propose_crs gives a centred projection for a wide region", {
  expect_message(crs <- propose_crs(wide_points), "custom",
                 class = "riskmap_message")
  expect_s3_class(crs, "crs")
  expect_match(crs$Name, "^RiskMap ")
  expect_true(is.na(crs$epsg))
  expect_false(st_is_longlat(crs))

  # the extent of the Democratic Republic of the Congo is about 2,900 km
  # across its diagonal, so the best is about 0.65%
  distortion <- crs_distortion(wide_points, crs)
  expect_lt(distortion$max_error_pct, 0.7)
  expect_lt(distortion$max_anisotropy_pct, 1e-4)
  # balanced about 1, up to the rounding of the scale factor
  expect_lt(abs(distortion$min_scale + distortion$max_scale - 2), 1e-4)

  # the UTM zone at the centre distorts distances by more than 1%
  expect_gt(crs_distortion(wide_points, 32734)$max_error_pct, 1)

  # the result can be used to project the data
  projected <- st_transform(wide_points, crs)
  expect_false(st_is_longlat(projected))

  # a larger tolerance accepts the best UTM zone, with a warning that a
  # centred projection is better
  expect_warning(
    epsg <- propose_crs(wide_points, epsg_tolerance = 0.5, messages = FALSE),
    "smaller", class = "riskmap_warning_crs"
  )
  expect_true(is.integer(epsg))
})

test_that("the name of a centred projection identifies it", {
  crs <- propose_crs(wide_points, messages = FALSE)
  narrower <- lonlat_points(lon = c(13, 30, 13, 30, 21.5),
                            lat = c(-12, -12, 4, 4, -4))
  other <- propose_crs(narrower, messages = FALSE)
  expect_false(crs == other)
  expect_false(identical(crs$Name, other$Name))
  expect_match(crs$Name, " k0 0\\.9[0-9]{3}$")
})

test_that("propose_crs does not depend on the axis order or on digits", {
  crs <- propose_crs(wide_points, messages = FALSE)
  axis_order <- st_axis_order(TRUE)
  withr::defer(st_axis_order(axis_order))
  expect_identical(propose_crs(wide_points, messages = FALSE)$wkt, crs$wkt)
  expect_identical(propose_crs(lonlat_points(c(33, 34), c(5, 6)),
                               messages = FALSE), 32636L)
  st_axis_order(axis_order)

  points <- lonlat_points(c(118, 129, 118, 129, 123.6), c(-13, -13, 5, 5, -4))
  expected <- propose_crs(points, messages = FALSE)$wkt
  withr::local_options(digits = 3)
  expect_identical(propose_crs(points, messages = FALSE)$wkt, expected)
})

test_that("propose_crs is reproducible and shared by similar extents", {
  crs <- propose_crs(wide_points, messages = FALSE)
  expect_identical(propose_crs(wide_points, messages = FALSE)$wkt, crs$wkt)
  shifted <- lonlat_points(lon = c(12.05, 31, 12, 31.05, 21.5),
                           lat = c(-13, -13.05, 5, 5, -4))
  expect_identical(propose_crs(shifted, messages = FALSE)$wkt, crs$wkt)
  expect_no_warning(propose_crs(wide_points, messages = FALSE))
})

test_that("propose_crs handles polygons and projected data", {
  box <- lonlat_box(12, 31, -13, 5)
  crs <- propose_crs(box, messages = FALSE)
  expect_s3_class(crs, "crs")
  expect_lt(crs_distortion(box, crs)$max_error_pct, 0.7)

  small <- lonlat_box(-10, -8, 6, 8)
  expect_identical(propose_crs(small, messages = FALSE), 32629L)

  # an sfc in a projected CRS gives the same as its longitude and latitude
  projected <- st_transform(st_geometry(wide_points), 3857)
  expect_identical(propose_crs(projected, messages = FALSE)$wkt,
                   propose_crs(wide_points, messages = FALSE)$wkt)
  expect_identical(propose_crs(gaussian_data, messages = FALSE),
                   propose_crs(latlon_data, messages = FALSE))
})

test_that("propose_crs reports no infinite distortion for global data", {
  global <- lonlat_points(rep(seq(-180, 120, 60), 3), rep(c(-60, 0, 60), each = 6))
  texts <- character()
  result <- tryCatch(
    withCallingHandlers(
      propose_crs(global),
      riskmap_message = function(m) {
        texts <<- c(texts, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      riskmap_warning_crs = function(w) {
        texts <<- c(texts, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    riskmap_error_input = function(e) e
  )
  expect_false(any(grepl("Inf", texts, fixed = TRUE)))
})

test_that("propose_crs warns for a region too wide for a planar CRS", {
  africa <- lonlat_points(lon = c(-17, 51, 15, 20, 35, -10),
                          lat = c(15, 11, -30, 30, -25, 5))
  expect_warning(crs <- propose_crs(africa, messages = FALSE),
                 "too wide", class = "riskmap_warning_crs")
  expect_gt(crs_distortion(africa, crs)$max_error_pct, 1)
})

test_that("crs_distortion measures known projections", {
  # on the central meridian of a UTM zone the scale is 0.9996
  meridian <- lonlat_points(c(39, 39), c(1, 2))
  distortion <- crs_distortion(meridian, 32637)
  expect_equal(distortion$min_scale, 0.9996, tolerance = 1e-5)
  expect_equal(distortion$max_error_pct, 0.04, tolerance = 1e-2)
  expect_lt(distortion$max_anisotropy_pct, 1e-4)
  expect_identical(names(distortion),
                   c("max_error_pct", "max_anisotropy_pct", "min_scale",
                     "max_scale", "n_locations"))

  # an equal-area projection is not conformal
  laea <- "+proj=laea +lat_0=-4 +lon_0=21.5 +datum=WGS84 +units=m"
  expect_gt(crs_distortion(wide_points, laea)$max_anisotropy_pct, 0.1)

  # coordinates in feet are converted
  feet <- "+proj=tmerc +lat_0=0 +lon_0=39 +k_0=0.9996 +x_0=0 +y_0=0 +datum=WGS84 +units=us-ft"
  expect_equal(crs_distortion(meridian, feet), distortion, tolerance = 1e-6)

  # the CRS of projected data is used by default
  expect_identical(crs_distortion(st_transform(meridian, 32637)),
                   distortion)
})

test_that("crs_distortion checks its arguments", {
  expect_input_error(crs_distortion(latlon_data), "crs")
  expect_input_error(crs_distortion(gaussian_data, crs = 4326), "crs")
  expect_input_error(crs_distortion(gaussian_data, crs = "nonsense"), "crs")
})

test_that("a CRS that distorts distances by more than 1% is reported", {
  # Pseudo-Mercator at the latitudes of Europe
  europe <- lonlat_points(c(0, 10, 0, 10), c(45, 45, 55, 55))
  expect_warning(
    summarise_distance(europe, distance_crs = 3857),
    "distance_crs", class = "riskmap_warning_crs"
  )
  expect_warning(
    summarise_distance(st_transform(europe, 3857)),
    "propose_crs", class = "riskmap_warning_crs"
  )
  expect_no_warning(summarise_distance(gaussian_data))
  # refits inside RiskMap don't repeat the warning
  expect_no_warning(
    without_crs_check(summarise_distance(st_transform(europe, 3857)))
  )
})

test_that("create_grid checks the CRS of a raster", {
  skip_if_not_installed("terra")
  # Pseudo-Mercator over Europe
  counts <- terra::rast(terra::ext(0, 1e6, 5.6e6, 7.4e6), resolution = 1e5,
                        crs = "EPSG:3857")
  terra::values(counts) <- seq_len(terra::ncell(counts))
  names(counts) <- "count"
  expect_warning(create_grid(raster = counts), "raster",
                 class = "riskmap_warning_crs")
})

test_that("longitude and latitude are projected with the CRS reported", {
  # The messages, unwrapped
  messages_of <- function(code) {
    texts <- character()
    value <- withCallingHandlers(code, riskmap_message = function(m) {
      text <- gsub("\\s+", " ", cli::ansi_strip(conditionMessage(m)))
      texts <<- c(texts, text)
      invokeRestart("muffleMessage")
    })
    list(value = value, messages = texts)
  }

  result <- messages_of(summarise_distance(wide_points))
  expect_length(result$messages, 1)
  expect_match(result$messages, "projected to a custom", fixed = TRUE)
  expect_match(result$messages, "for computing distances in km.", fixed = TRUE)
  expect_match(result$messages, "No UTM zone keeps distances within 0.5% over them.",
               fixed = TRUE)
  expect_match(result$messages, "Set `distance_crs` to choose the projection.",
               fixed = TRUE)
  # The distortion is left to the warning above 1% and to crs_distortion().
  expect_no_match(result$messages, "distorted", fixed = TRUE)

  result <- messages_of(create_grid(lonlat_box(-10, -8, 6, 8), spacing = 50))
  expect_length(result$messages, 1)
  expect_match(result$messages,
               "projected to EPSG:32629 (WGS 84 / UTM zone 29N) to create the grid, with distances in km.",
               fixed = TRUE)
  expect_match(result$messages, "Project `boundaries` yourself", fixed = TRUE)
  expect_identical(st_crs(result$value), st_crs(32629))
})

test_that("a fit in longitude and latitude says where it is fitted once", {
  messages_of <- function(code) {
    texts <- character()
    value <- withCallingHandlers(code, riskmap_message = function(m) {
      texts <<- c(texts, gsub("\\s+", " ", cli::ansi_strip(conditionMessage(m))))
      invokeRestart("muffleMessage")
    })
    list(value = value, messages = texts)
  }
  fit_to <- function(data) {
    messages_of(suppressWarnings(glgpm(
      y ~ cov + gp(), data = data, family = "binomial",
      denominator = denominator, control_mcml = fixture_control_mcml,
      seed = fixture_seed, messages = TRUE
    )))
  }
  lonlat <- fit_to(sf::st_transform(binomial_data, 4326))
  crs_messages <- grep("EPSG|fitted in", lonlat$messages, value = TRUE)
  expect_length(crs_messages, 1L)
  expect_match(crs_messages, "for fitting, with distances in km.", fixed = TRUE)
  expect_match(crs_messages, "Set `model_crs` to choose the projection.",
               fixed = TRUE)
  # The fit records where its starting values came from, and summary() says
  # it, so the fit doesn't announce it.
  # This fit takes the transformed-Gaussian start, the case that had a
  # message, so its absence here shows the message is gone.
  expect_identical(attr(lonlat$value, "starting_values")$selected,
                   "transformed_gaussian")
  expect_false(any(grepl("starting values", lonlat$messages, fixed = TRUE)))
  expect_match(squish(capture.output(print(summary(lonlat$value)))),
               describe_starting_values("transformed_gaussian"),
               fixed = TRUE)
  # Projected data get the line on where the model is fitted instead.
  projected <- fit_to(binomial_data)
  expect_identical(grep("EPSG|fitted in", projected$messages, value = TRUE),
                   "i The model is fitted in EPSG:32637, with distances in km.")
})
