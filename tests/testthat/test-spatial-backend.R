test_that("compiled pairwise distances match stats::dist", {
  coordinates <- rbind(c(1e7, -1e7),
                       c(1e7 + 3, -1e7 + 4),
                       c(1e7 - 8, -1e7 + 15),
                       c(1e7 + 20, -1e7 - 7))

  observed <- pairwise_distances(coordinates)
  expected <- dist(coordinates)

  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-14)
  expect_equal(attr(observed, "Size"), attr(expected, "Size"))
})

test_that("compiled cross-distances match direct Euclidean distances", {
  observed_locations <- rbind(c(0, 0), c(3, 4), c(-2, 5))
  prediction_locations <- rbind(c(1, 1), c(8, -3))
  expected <- outer(seq_len(nrow(prediction_locations)),
                    seq_len(nrow(observed_locations)),
                    Vectorize(function(i, j) {
                      sqrt(sum((prediction_locations[i, ] -
                        observed_locations[j, ])^2))
                    }))

  expect_equal(cross_distances(prediction_locations, observed_locations),
               expected, tolerance = 1e-14)
})

test_that("compiled distances reject invalid coordinate matrices", {
  expect_error(pairwise_distances(matrix(c(1, NA), ncol = 1)),
               "finite values")
  expect_error(cross_distances(matrix(1:4, ncol = 2),
                               matrix(1:6, ncol = 3)),
               "same number of columns")
})

test_that("half-integer Matérn kernels agree with the general definition", {
  distances <- c(0, 1e-8, 0.1, 1, 10, 1000)
  phi <- 2.3

  general_matern <- function(kappa) {
    scaled <- distances / phi
    output <- ifelse(distances > 0,
                     2^(1 - kappa) / gamma(kappa) * scaled^kappa *
                       besselK(scaled, kappa),
                     1)
    output[distances > 600 * phi] <- 0
    output
  }

  for (kappa in c(0.5, 1.5, 2.5)) {
    expect_equal(matern_correlation(distances, phi, kappa),
                 general_matern(kappa), tolerance = 1e-12)
  }
})

test_that("compiled Matérn kernels preserve cross-covariance dimensions", {
  distances <- matrix(seq(0, 2, length.out = 12), nrow = 3)

  observed <- matern_correlation(distances, phi = 0.7, kappa = 1.5)

  expect_equal(dim(observed), dim(distances))
})

test_that("distance scaling leaves Matérn correlations unchanged", {
  coordinates <- rbind(c(500000, 5700000),
                       c(503000, 5704000),
                       c(515000, 5698000))
  phi <- 7500
  scaled <- scale_spatial_coordinates(coordinates, phi)

  original_correlation <- matern_correlation(pairwise_distances(coordinates),
                                              phi, 1.5)
  scaled_correlation <- matern_correlation(
    pairwise_distances(scaled$coordinates), scaled$phi, 1.5
  )

  expect_equal(scaled_correlation, original_correlation, tolerance = 5e-14)
  expect_equal(restore_spatial_range(scaled$phi, scaled$distance_scale), phi)
})

test_that("distance scaling rejects coincident observed locations", {
  expect_error(scale_spatial_coordinates(matrix(1, nrow = 3, ncol = 2)),
               "at least two distinct locations")
})
