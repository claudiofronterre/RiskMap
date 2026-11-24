test_that("extrapolation indicators separate spatial and covariate novelty", {
  locations <- st_coordinates(gaussian_data)
  cov_range <- range(gaussian_data$cov)
  centre <- colMeans(locations)
  newdata <- st_as_sf(
    data.frame(
      x = c(locations[1, 1], centre[[1]], locations[1, 1], 1e7, 5e6),
      z = c(locations[1, 2], centre[[2]], locations[1, 2], 1e7, 5e6),
      cov = c(gaussian_data$cov[1], mean(cov_range),
              cov_range[2] + diff(cov_range), cov_range[1], NA)
    ),
    coords = c("x", "z"), crs = 32637
  )
  # The fixture's short scale leaves a sampled location with no spatial
  # information from the others, so the classes are checked with a fixed
  # threshold.
  expect_warning(
    extrapolation <- assess_extrapolation(gaussian_model, newdata,
                                          spatial_threshold = 0.8),
    "missing values, so"
  )

  expect_s3_class(extrapolation, c("RiskMap_extrapolation", "sf"))
  expect_identical(nrow(extrapolation), 5L)
  expect_equal(extrapolation$covariates_outside[1:4], c(0L, 0L, 1L, 0L))
  expect_equal(extrapolation$covariate_novelty[3], 1)
  expect_identical(extrapolation$most_novel_covariate[3], "cov")
  expect_true(is.na(extrapolation$most_novel_covariate[1]))
  expect_lt(extrapolation$spatial_unexplained[1],
            extrapolation$spatial_unexplained[2])
  expect_equal(extrapolation$spatial_unexplained[4], 1)
  expect_identical(
    as.character(extrapolation$extrapolation[c(1, 3, 4)]),
    c("Supported", "Covariate", "Spatial")
  )
  expect_true(all(is.na(st_drop_geometry(extrapolation)[5, c(
    "covariates_outside", "covariate_novelty", "covariate_leverage",
    "spatial_unexplained", "extrapolation"
  )])))
  expect_lte(max(extrapolation$covariate_leverage[1:2]), 1 + 1e-8)
  expect_gt(extrapolation$covariate_leverage[3], 1)
})

test_that("unexplained spatial variance matches the kriging formula", {
  parameters <- coef(gaussian_intercept_model)
  observed <- gaussian_intercept_model$coords
  ratio <- parameters$sigma2_me / parameters$sigma2
  correlation <- matern_correlation(as.matrix(dist(observed)),
                                    phi = parameters$phi,
                                    kappa = gaussian_intercept_model$kappa)
  diag(correlation) <- 1 + ratio
  point <- observed[1, , drop = FALSE] + c(0.5, -0.3)
  cross <- matern_correlation(
    sqrt(colSums((t(observed) - as.numeric(point))^2)),
    phi = parameters$phi, kappa = gaussian_intercept_model$kappa
  )
  expected <- 1 - drop(cross %*% solve(correlation, cross))

  factor <- crs_to_distance_factor(gaussian_data,
                                   gaussian_intercept_model$distance_units)
  newdata <- st_as_sf(data.frame(x = point[1] / factor, z = point[2] / factor),
                      coords = c("x", "z"), crs = 32637)
  extrapolation <- assess_extrapolation(gaussian_intercept_model, newdata,
                                        spatial_threshold = 0.8)
  expect_equal(extrapolation$spatial_unexplained, expected, tolerance = 1e-8)
  expect_identical(extrapolation$covariates_outside, 0L)
})

test_that("extrapolation maps and argument checks", {
  boundary <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
    c(0, 0), c(10000, 0), c(10000, 10000), c(0, 10000), c(0, 0)
  ))), crs = 32637))
  grid <- create_grid(boundary, spacing = 2, distance_units = "km")
  grid$cov <- seq(-3, 3, length.out = nrow(grid))
  # The fixture's locations share almost no correlation, which leaves the
  # default threshold at 1, so the maps use a fixed one.
  extrapolation <- assess_extrapolation(gaussian_model, grid,
                                        spatial_threshold = 0.8)
  expect_false(is.null(attr(extrapolation, "grid_metadata")))
  for (indicator in c("extrapolation", "spatial_unexplained",
                      "covariate_novelty", "covariate_leverage")) {
    map <- plot(extrapolation, indicator = indicator)
    expect_s3_class(map, "ggplot")
    expect_s3_class(map$layers[[1L]]$geom, "GeomSf")
  }
  expect_identical(ggplot2::get_labs(plot(extrapolation))$fill,
                   "Extrapolation")
  # Every class gets a coloured key, including classes absent from the map
  expect_true(plot(extrapolation)$layers[[1L]]$show.legend)
  expect_true(is.na(
    plot(extrapolation,
         indicator = "spatial_unexplained")$layers[[1L]]$show.legend
  ))
  with_missing <- grid
  with_missing$cov[1] <- NA
  expect_warning(missing_map <- plot(assess_extrapolation(
    gaussian_model, with_missing, spatial_threshold = 0.8
  )),
                 "missing values, so")
  expect_identical(ggplot2::get_guide_data(missing_map, "fill")$.label,
                   c("Supported", "Spatial", "Covariate",
                     "Spatial and covariate"))
  expect_true("grey90" %in% ggplot2::ggplot_build(missing_map)$data[[1]]$fill)
  skip_if_not_installed("patchwork")
  both <- plot(extrapolation,
               indicator = c("spatial_unexplained", "extrapolation"))
  expect_length(c(both$patches$plots, list(both)), 2L)
  expect_identical(both$patches$layout$ncol, 2L)
  expect_input_error(plot(extrapolation, indicator = "distance"), "indicator")
  # An argument that plot() doesn't use stops instead of being ignored.
  expect_input_error(plot(extrapolation, indicators = "extrapolation"),
                     "indicators", "Check the spelling")
  expect_input_error(assess_extrapolation(gaussian_model, grid,
                                          spatial_threshold = 2),
                     "spatial_threshold")
  expect_input_error(assess_extrapolation(list(), grid), "object")
})

test_that("an extrapolation map at points not made by create_grid() draws every point", {
  extrapolation <- assess_extrapolation(gaussian_model, gaussian_data,
                                        spatial_threshold = 0.8)
  expect_null(attr(extrapolation, "grid_metadata"))
  map <- plot(extrapolation, indicator = "spatial_unexplained")
  expect_no_warning(drawn <- ggplot2::layer_grob(map, 1L)[[1L]])
  expect_length(drawn$x, nrow(gaussian_data))
})

test_that("unexplained spatial variance uses fixed variance components", {
  fixed_model <- glgpm(y ~ gp(), data = gaussian_data, family = "gaussian",
                       fixed_pars = list(sigma2_me = 0.2), messages = FALSE)
  estimated_model <- fixed_model
  estimated_model$fixed_pars <- NULL
  estimated_model$estimate$sigma2_me <- log(0.2)
  share <- function(model) {
    assess_extrapolation(model, gaussian_data,
                         spatial_threshold = 0.8)$spatial_unexplained
  }
  expect_equal(share(fixed_model), share(estimated_model), tolerance = 1e-12)

  nugget_model <- gaussian_intercept_model
  nugget_model$nugget <- "fixed"
  nugget_model$fixed_pars <- list(tau2 = 0.25)
  expect_true(all(share(nugget_model) > share(gaussian_intercept_model)))
})

test_that("the unexplained share includes the binomial sampling noise", {
  # The fixture's scale leaves its locations almost uncorrelated, so a longer
  # one gives a leave-one-out threshold well below 1.
  binomial_model$estimate$phi <- log(4)
  parameters <- coef(binomial_model)
  sigma2 <- parameters$sigma2
  observed <- binomial_model$coords
  correlation <- matern_correlation(as.matrix(dist(observed)),
                                    phi = parameters$phi,
                                    kappa = binomial_model$kappa)
  # The conditional mode of the effects at the estimates, found by
  # maxim_integrand(), and the binomial information m p (1 - p) there.
  ratio <- nugget_ratio_natural(fitted_nugget_variance(binomial_model),
                                sigma2)
  mode <- maxim_integrand(
    y = binomial_model$y, denominator = binomial_model$units_m,
    mu = as.numeric(binomial_model$D %*% parameters$beta),
    Sigma = sigma2 * (correlation + diag(max(ratio, 1e-10), nrow(observed))),
    ID_coords = binomial_model$ID_coords, ID_re = binomial_model$ID_re,
    sigma2_re = parameters$sigma2_re, family = "binomial"
  )$mode
  eta <- as.numeric(binomial_model$D %*% parameters$beta) +
    mode[binomial_model$ID_coords] +
    mode[nrow(observed) + binomial_model$ID_re[, 1L]]
  p <- stats::plogis(eta)
  information <- rowsum(binomial_model$units_m * p * (1 - p),
                        binomial_model$ID_coords)[, 1L]
  noisy <- correlation
  diag(noisy) <- 1 + ratio + 1 / (information * sigma2)
  point <- observed[1, , drop = FALSE] + c(0.4, 0.2)
  cross <- matern_correlation(
    sqrt(colSums((t(observed) - as.numeric(point))^2)),
    phi = parameters$phi, kappa = binomial_model$kappa
  )
  expected <- 1 - drop(cross %*% solve(noisy, cross))
  # Each location left out in turn, from the other locations' data
  left_out <- vapply(seq_len(nrow(observed)), function(j) {
    1 - drop(correlation[j, -j] %*%
               solve(noisy[-j, -j], correlation[-j, j]))
  }, numeric(1))

  factor <- crs_to_distance_factor(binomial_data,
                                   binomial_model$distance_units)
  newdata <- st_as_sf(data.frame(x = point[1] / factor, z = point[2] / factor,
                                 cov = binomial_data$cov[1]),
                      coords = c("x", "z"), crs = 32637)
  extrapolation <- assess_extrapolation(binomial_model, newdata)
  expect_equal(extrapolation$spatial_unexplained, expected, tolerance = 1e-5)
  threshold <- stats::quantile(left_out, 0.95, names = FALSE)
  expect_lt(threshold, 0.95)
  expect_lt(threshold, max(left_out) - 0.01)
  expect_equal(attr(extrapolation, "spatial_threshold"), threshold,
               tolerance = 1e-6)
  at_sites <- assess_extrapolation(binomial_model, binomial_data)
  expect_identical(as.character(at_sites$extrapolation),
                   rep("Supported", nrow(binomial_data)))
  expect_identical(attr(extrapolation, "spatial_rule"), "leave_one_out")
  expect_true(attr(extrapolation, "sampling_noise"))
  expect_identical(
    as.character(extrapolation$extrapolation),
    if (expected > threshold) "Spatial" else "Supported"
  )

  # Without the noise, the share is that of the locations and the nugget,
  # which is lower everywhere.
  without <- assess_extrapolation(binomial_model, binomial_data,
                                  sampling_noise = FALSE)
  with_noise <- assess_extrapolation(binomial_model, binomial_data)
  expect_false(attr(without, "sampling_noise"))
  expect_true(all(without$spatial_unexplained <
                    with_noise$spatial_unexplained))
  expect_equal(without$spatial_unexplained, rep(0, nrow(binomial_data)),
               tolerance = 1e-6)
})

test_that("a Gaussian share is the kriging variance with the records averaged", {
  repeated <- rbind(gaussian_data, gaussian_data[1:3, ])
  withr::with_seed(3, repeated$y <- repeated$y + stats::rnorm(nrow(repeated),
                                                              sd = 0.3))
  fit <- glgpm(y ~ gp(), data = repeated, family = "gaussian",
               messages = FALSE)
  expect_identical(nrow(fit$coords), nrow(gaussian_data))
  parameters <- coef(fit)
  observed <- fit$coords
  records <- observed[fit$ID_coords, ]
  point <- observed[2, , drop = FALSE] + c(0.3, 0.1)
  correlation <- function(a, b) {
    distances <- sqrt(pmax(outer(rowSums(a^2), rowSums(b^2), "+") -
                             2 * a %*% t(b), 0))
    matern_correlation(distances, phi = parameters$phi, kappa = fit$kappa)
  }
  same_location <- outer(fit$ID_coords, fit$ID_coords, "==")
  covariance <- parameters$sigma2 * correlation(records, records) +
    fitted_nugget_variance(fit, parameters) * same_location +
    diag(parameters$sigma2_me, nrow(records))
  cross <- parameters$sigma2 * correlation(point, records)
  expected <- 1 - drop(cross %*% solve(covariance, t(cross))) /
    parameters$sigma2

  factor <- crs_to_distance_factor(repeated, fit$distance_units)
  newdata <- st_as_sf(data.frame(x = point[1] / factor, z = point[2] / factor),
                      coords = c("x", "z"), crs = 32637)
  expect_equal(assess_extrapolation(fit, newdata)$spatial_unexplained,
               expected, tolerance = 1e-8)
})

test_that("documents that print data frames as tables show the print", {
  # Stands for the methods that R Markdown and Quarto use for data frames.
  knit_print.data.frame <- function(x, ...) cat("a table\n")
  expect_output(knitr::knit_print(data.frame(x = 1)), "a table")
  extrapolation <- assess_extrapolation(binomial_model, binomial_data,
                                        spatial_threshold = 0.5)
  expect_identical(capture.output(knitr::knit_print(extrapolation)),
                   capture.output(print(extrapolation)))
  # A subset without the classes is a table like any other.
  expect_output(knitr::knit_print(extrapolation[, "covariates_outside"]),
                "a table")
  listed <- glgpm_parameters(y ~ cov + gp(), data = binomial_data,
                             family = "binomial")
  expect_identical(capture.output(knitr::knit_print(listed)),
                   capture.output(print(listed)))
})

test_that("the print of indicators states the classes and the thresholds", {
  # A longer scale than the fixture's, so that the threshold is below 1
  binomial_model$estimate$phi <- log(4)
  extrapolation <- assess_extrapolation(binomial_model, binomial_data)
  printed <- capture.output(print(extrapolation))
  expect_match(printed[1L], "Extrapolation")
  expect_true("Indicators at 10 locations." %in% printed)
  share <- paste0(format(100 * attr(extrapolation, "spatial_threshold"),
                         digits = 3), "%")
  expect_match(squish(printed), paste0(
    "Spatial: more than ", share, " of the spatial variance left ",
    "unexplained, the 95th percentile of the shares at the sampled ",
    "locations when their own data are left out."
  ), fixed = TRUE)
  printed_again <- capture.output(shown <- withVisible(print(extrapolation)))
  expect_false(shown$visible)
  expect_identical(shown$value, extrapolation)
  given <- assess_extrapolation(binomial_model, binomial_data,
                                spatial_threshold = 0.5,
                                sampling_noise = FALSE)
  expect_identical(attr(given, "spatial_threshold"), 0.5)
  expect_identical(attr(given, "spatial_rule"), "given")
  printed <- squish(capture.output(print(given)))
  expect_match(printed, "more than 50% of the spatial variance left unexplained",
               fixed = TRUE)
  expect_match(printed, "leaves out the sampling noise", fixed = TRUE)
  # Locations that share almost no correlation set the threshold at 1, so
  # nothing can be flagged, which assess_extrapolation() warns of and the
  # print says.
  far <- st_as_sf(data.frame(x = 1e7, z = 1e7, cov = 0), coords = c("x", "z"),
                  crs = 32637)
  expect_warning(isolated <- assess_extrapolation(gaussian_model, far),
                 "no location can be classed as spatial extrapolation",
                 class = "riskmap_warning")
  expect_equal(attr(isolated, "spatial_threshold"), 1, tolerance = 1e-12)
  expect_identical(as.character(isolated$extrapolation), "Supported")
  printed <- squish(capture.output(print(isolated)))
  expect_match(printed, paste0(
    "Spatial: none, since the threshold is 100%, the 95th percentile of the ",
    "shares at the sampled locations when their own data are left out."
  ), fixed = TRUE)
  expect_match(printed, "unexplained at more than 5% of them", fixed = TRUE)
  expect_no_match(squish(capture.output(print(extrapolation))),
                  "Spatial: none", fixed = TRUE)
  # A threshold that is given, or a percentile below 1, gives no warning.
  expect_no_warning(assess_extrapolation(gaussian_model, far,
                                         spatial_threshold = 0.8))
  expect_no_warning(assess_extrapolation(binomial_model, binomial_data))
  # A threshold below 1 never prints as 100%.
  expect_identical(
    vapply(c(0.5, 0.8181409, 0.99978, 0.9999999, 1),
           format_extrapolation_share, character(1)),
    c("50%", "81.8%", "99.98%", "99.99999%", "100%")
  )
  # A subset without the classes prints as an sf object, with the arguments
  # of the sf method.
  expect_output(print(extrapolation[1:2, "covariates_outside"]),
                "Simple feature collection")
  expect_output(print(extrapolation[, "covariates_outside"], n = 2),
                "First 2 features")
  # Otherwise an argument that print() doesn't use stops.
  expect_input_error(print(extrapolation, n = 2), "n")
  expect_input_error(assess_extrapolation(binomial_model, binomial_data,
                                          spatial_threshold = 1),
                     "spatial_threshold")
  expect_input_error(assess_extrapolation(binomial_model, binomial_data,
                                          sampling_noise = NA),
                     "sampling_noise")
})
