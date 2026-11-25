# Calibration curves of cross-validated predictions.

# The curve of each row by its definition.
direct_curves <- function(lower, upper, levels) {
  t(vapply(seq_len(nrow(lower)), function(i) {
    vapply(levels, function(u) {
      width <- upper[i, ] - lower[i, ]
      mean(ifelse(width > 0, pmin(pmax((u - lower[i, ]) / width, 0), 1),
                  as.numeric(lower[i, ] <= u)))
    }, numeric(1))
  }, numeric(length(levels))))
}

test_that("calibration curves average the PIT functions of the observations", {
  set.seed(40)
  levels <- calibration_levels()
  # Counts, with ends at 0 and 1 and on the levels, and equal ends.
  lower <- matrix(runif(60), 4)
  upper <- pmin(lower + matrix(runif(60, 0, 0.3), 4), 1)
  lower[1, 1:3] <- c(0, 0.25, 0.9)
  upper[1, 1:3] <- c(0.1, 0.5, 1)
  lower[2, 1:2] <- upper[2, 1:2] <- c(0.3, 1)
  expect_equal(pit_curves(lower, upper, levels),
               direct_curves(lower, upper, levels), tolerance = 1e-12)
  # In blocks of observations, as for many of them.
  expect_equal(pit_curves(lower, upper, levels, block_size = 9),
               direct_curves(lower, upper, levels), tolerance = 1e-12)
  # Gaussian, where the curve is the empirical distribution function.
  pit <- matrix(c(runif(48), 0, 1), 1)
  expect_equal(as.numeric(pit_curves(pit, pit, levels)),
               stats::ecdf(pit)(levels), tolerance = 1e-12)
  # One observation.
  expect_equal(pit_curves(matrix(0.2), matrix(0.4), c(0, 0.2, 0.3, 0.4, 1)),
               matrix(c(0, 0, 0.5, 1, 1), 1))
})

test_that("pit_area() is the area between the curve and the identity line", {
  u <- seq(0, 1, length.out = 1001)
  curves <- rbind(u, 0, u^2, 1 - u)
  expect_equal(pit_area(curves, u), c(0, 1 / 2, 1 / 6, 1 / 2),
               tolerance = 1e-6)
})

test_that("the PIT histogram has the increases of the curve", {
  set.seed(41)
  levels <- calibration_levels()
  pit <- c(runif(50), 0, 1)
  curve <- pit_curves(t(pit), t(pit), levels)
  for (bins in c(5, 10, 20)) {
    heights <- pit_histogram(curve, bins)
    edges <- seq(0, 1, length.out = bins + 1)
    counts <- graphics::hist(pit, breaks = edges, plot = FALSE)$counts
    expect_equal(as.numeric(heights), bins * counts / length(pit))
    expect_equal(sum(heights) / bins, 1)
  }
  # The heights of counts are the averages of those of the observations.
  expect_equal(
    as.numeric(pit_histogram(pit_curves(t(0.1), t(0.5), levels), 5)),
    c(1.25, 2.5, 1.25, 0, 0)
  )
})

# The reference from PIT values without ends.
reference_of <- function(pit) {
  curves <- pit_curves(pit, pit, area_levels())
  list(curves = curves[, seq(1, 1001, by = 10)],
       areas = pit_area(curves, area_levels()))
}

test_that("assess_calibration() gives the area, its reference and the band", {
  set.seed(42)
  n <- 80
  pit <- runif(n)
  simulated <- matrix(runif(999 * n), 999)
  out <- assess_calibration(pit, pit, reference_of(simulated))
  u <- seq(0, 1, length.out = 1001)
  expect_equal(out$area, pit_area(t(stats::ecdf(pit)(u)), u))
  areas <- apply(simulated, 1, function(x) pit_area(t(stats::ecdf(x)(u)), u))
  expect_equal(out$expected, mean(areas))
  expect_equal(out$limit, sort(areas)[950])
  # About 0.31 / sqrt(n) for independent observations.
  expect_equal(out$expected, 0.31 / sqrt(n), tolerance = 0.05)
  expect_identical(dim(out$simulated), c(999L, 101L))
  expect_identical(out$curve$u, calibration_levels())
  expect_equal(out$curve$curve, stats::ecdf(pit)(out$curve$u))
  # The band is the joint envelope of the curve of the data and the simulated.
  envelope <- global_rank_envelope(cbind(out$curve$curve, t(out$simulated)), 0.95)
  expect_equal(out$curve$lower, envelope$lower)
  expect_equal(out$curve$upper, envelope$upper)
  inside <- apply(out$simulated, 1, function(curve) {
    all(curve >= out$curve$lower & curve <= out$curve$upper)
  })
  expect_equal(mean(inside), 0.95, tolerance = 0.02)
  empty <- assess_calibration(numeric(0), numeric(0), NULL)
  expect_true(is.na(empty$area))
  expect_true(is.na(empty$expected))
  expect_true(is.na(empty$limit))
  expect_null(empty$curve)
})

test_that("a correct model exceeds the limit of the area with a probability of 5%", {
  # The value of the data is one of a hundred exchangeable values.
  values <- c(3, 1, 2, rep(0, 16))
  expect_identical(simulated_limit(values), 3)
  expect_identical(simulated_limit(as.numeric(1:99)), 95)
  # Too few simulations for a limit.
  expect_identical(simulated_limit(as.numeric(1:18)), NA_real_)
  expect_identical(simulated_limit(as.numeric(1:19)), 19)
})

test_that("the bands and the limit have their coverage for a correct model", {
  skip_on_cran()
  # Dependent PIT values, as of one test set.
  set.seed(43)
  n <- 30
  root <- chol(0.5 * exp(-as.matrix(dist(seq_len(n))) / 8) + diag(0.5, n))
  draw <- function(k) pnorm(matrix(rnorm(k * n), k) %*% root)
  outside <- replicate(400, {
    pit <- draw(1)
    simulated <- draw(199)
    out <- assess_calibration(as.numeric(pit), as.numeric(pit),
                              reference_of(simulated))
    height <- as.numeric(pit_histogram(t(out$curve$curve), 10))
    bars <- calibration_band(height, pit_histogram(out$simulated, 10))
    c(area = out$area > out$limit,
      curve = any(out$curve$curve < out$curve$lower |
                    out$curve$curve > out$curve$upper),
      bars = any(height < bars$lower | height > bars$upper))
  })
  # 5%, with a standard error of 1.1%.
  expect_true(all(rowMeans(outside) > 0.02 & rowMeans(outside) < 0.08))
})

test_that("a Gaussian observation has the PIT value of its first test set", {
  predictive <- function(rows, draws) {
    list(rows = rows, draws = draws, mean = c(0, 0), sd = c(1, 1))
  }
  traits <- family_traits("gaussian")
  draws <- matrix(c(0.1, 0.2, 0.3, 0.4), 2)
  first <- simulate_gaussian_pit(predictive(c(2L, 4L), draws), traits,
                                 vector("list", 2))
  expect_equal(first$lower, pnorm(t(draws)))
  expect_identical(first$upper, first$lower)
  expect_equal(first$uniforms[[2]], pnorm(draws[2, ]))
  second <- simulate_gaussian_pit(predictive(c(4L, 1L), draws + 1), traits,
                                  list(first$uniforms[[2]], NULL))
  expect_equal(second$lower, cbind(pnorm(draws[2, ]), pnorm(draws[2, ] + 1)))
})

test_that("calibration_reference() has the curves and the areas", {
  set.seed(44)
  predictive <- lapply(list(1:3, 3:5), function(rows) {
    list(rows = rows, draws = matrix(rnorm(3 * 50), 3), mean = rep(0, 3),
         sd = rep(1, 3))
  })
  fit <- list(family = "gaussian")
  test_sets <- list(location = 1:6)
  set.seed(1)
  state <- .Random.seed
  out <- calibration_reference(fit, test_sets, predictive, 40L, c(1L, 2L))
  expect_identical(.Random.seed, state)
  expect_identical(dim(out$curves), c(40L, 101L))
  # The third observation has the same PIT value in both test sets.
  pit <- pnorm(rbind(predictive[[1]]$draws,
                     predictive[[1]]$draws[3, ], predictive[[2]]$draws[2:3, ]))
  direct <- t(apply(pit[, 1:40], 2, function(x) stats::ecdf(x)(area_levels())))
  expect_equal(out$curves, direct[, seq(1, 1001, by = 10)])
  expect_equal(out$areas, pit_area(direct, area_levels()))
  # A test set without scores isn't in it, and few draws give fewer curves.
  one <- calibration_reference(fit, test_sets, c(predictive[1], list(NULL)),
                               999L, c(1L, 2L))
  expect_identical(dim(one$curves), c(50L, 101L))
  expect_equal(one$curves[1:40, ],
               t(apply(pit[1:3, 1:40], 2, function(x) {
                 stats::ecdf(x)(calibration_levels())
               })))
})
