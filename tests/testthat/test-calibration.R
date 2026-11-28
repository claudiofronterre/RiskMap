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

# The test sets of a fit, with the predictive distributions that the reference
# takes, from counts drawn from the means of the fit.
reference_inputs <- function(fit, sets) {
  test_sets <- read_test_sets(sets, fit)
  traits <- family_traits(fit$family)
  units_m <- assessed_denominators(fit)
  predictive <- lapply(test_sets$sets, function(set) {
    rows <- set$test
    if (!traits$discrete) return(list(rows = rows))
    link <- as.numeric(fit$D[rows, , drop = FALSE] %*% coef(fit)$beta)
    means <- traits$inverse_link(link + matrix(rnorm(length(rows) * 50, 0, 0.5),
                                               length(rows)))
    c(list(rows = rows),
      score_counts(fit$y[rows], units_m[rows], means, fit_family_values(fit),
                   traits)$predictive)
  })
  list(test_sets = test_sets, predictive = predictive)
}

test_that("the working model has the covariance of the fit", {
  fit <- gaussian_model
  p <- coef(fit)
  working <- working_model(fit)
  expect_identical(working$free, c("sigma2", "phi", "sigma2_me", "sigma2_re.1"))
  spatial <- p$sigma2 * matern_correlation(
    as.matrix(dist(fit$coords))[fit$ID_coords, fit$ID_coords], phi = p$phi,
    kappa = fit$kappa
  )
  group <- outer(fit$data$i, fit$data$i, "==")
  covariance <- spatial + p$sigma2_re[[1]] * group + diag(p$sigma2_me, fixture_n)
  expect_equal(working$covariance(working$theta), covariance,
               ignore_attr = TRUE)
  # The derivatives are on the logarithmic scale.
  expect_equal(working$derivative("sigma2"), spatial, ignore_attr = TRUE,
               tolerance = 1e-6)
  expect_equal(working$derivative("sigma2_me"), diag(p$sigma2_me, fixture_n),
               ignore_attr = TRUE, tolerance = 1e-6)
  expect_identical(colnames(working$design), c("(Intercept)", "cov"))
  # A nugget, a measurement error held at zero, and a held scale.
  nugget <- working_model(gaussian_offset_model)
  expect_identical(nugget$free, c("sigma2", "phi", "tau2"))
  held <- glgpm(y ~ gp(nugget = TRUE), data = gaussian_data,
                family = "gaussian", fixed_pars = list(phi = 2, tau2 = 0.1),
                messages = FALSE)
  expect_identical(working_model(held)$free, c("sigma2", "sigma2_me"))
  # Held coefficients are known, so they aren't estimated.
  slope <- glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
                 fixed_pars = list(beta = c(cov = 0.5), phi = 2),
                 messages = FALSE)
  expect_identical(colnames(working_model(slope)$design), "(Intercept)")
  # Counts are observed with the inverse of their curvature.
  counts <- working_model(binomial_model)
  p <- coef(binomial_model)
  expect_equal(diag(counts$covariance(counts$theta)),
               p$sigma2 + p$sigma2_re[[1]] + 1 / conditional_weights(binomial_model))
})

test_that("candidate points move the free parameters around the estimates", {
  working <- working_model(gaussian_offset_model)
  points <- candidate_points(working, 50)
  expect_identical(dim(points), c(50L, 3L))
  expect_identical(points[1, ], working$theta)
  ratio <- log(points[-1, ] / rep(working$theta, each = 49))
  expect_true(all(abs(ratio) <= 4 + 1e-12))
  expect_true(all(apply(ratio, 2, stats::sd) > 0.1))
  # The same points every time, a Halton sequence.
  expect_identical(candidate_points(working, 50), points)
  expect_identical(candidate_points(working, 1), points[1, , drop = FALSE])
  # Without free parameters, only the estimates.
  held <- glgpm(y ~ gp(), data = gaussian_intercept_data, family = "gaussian",
                fixed_pars = list(sigma2 = 1, phi = 2, sigma2_me = 0.1),
                messages = FALSE)
  expect_identical(nrow(candidate_points(working_model(held), 50)), 1L)
  expect_identical(first_primes(5), c(2L, 3L, 5L, 7L, 11L))
  expect_equal(radical_inverse(1:4, 2), c(0.5, 0.25, 0.75, 0.125))
})

test_that("the leave-out errors of simulated datasets are those of direct refits", {
  working <- working_model(gaussian_model)
  n <- nrow(working$design)
  # A test set that leaves out a row it doesn't test.
  sets <- list(list(test = c(2, 5), left_out = c(2, 5, 7)),
               list(test = 1:3, left_out = 1:3), list(test = 9, left_out = 9))
  n_simulations <- 4
  n_candidates <- 6
  direct <- function(refit) {
    set.seed(3)
    simulated <- t(chol(working$covariance(working$theta))) %*%
      matrix(rnorm(n * n_simulations), n)
    points <- candidate_points(working, n_candidates)
    covariances <- lapply(seq_len(n_candidates), function(k) {
      working$covariance(points[k, ])
    })
    # The likelihood profiled over the coefficients, by generalised least
    # squares.
    profile <- function(V, z, rows) {
      X <- working$design[rows, , drop = FALSE]
      precision <- solve(V[rows, rows])
      beta <- solve(t(X) %*% precision %*% X, t(X) %*% precision %*% z[rows])
      residual <- z[rows] - X %*% beta
      list(beta = beta,
           value = -(determinant(V[rows, rows])$modulus +
                       t(residual) %*% precision %*% residual) / 2)
    }
    do.call(rbind, lapply(sets, function(set) {
      train <- setdiff(seq_len(n), set$left_out)
      vapply(seq_len(n_simulations), function(b) {
        z <- simulated[, b]
        rows <- if (refit) train else seq_len(n)
        values <- vapply(covariances, function(V) {
          as.numeric(profile(V, z, rows)$value)
        }, numeric(1))
        V <- covariances[[which.max(values)]]
        residual <- z - working$design %*% profile(V, z, rows)$beta
        weights <- V[set$test, train] %*% solve(V[train, train])
        error <- residual[set$test] - weights %*% residual[train]
        variance <- V[set$test, set$test] - weights %*% V[train, set$test]
        pnorm(as.numeric(error) / sqrt(diag(as.matrix(variance))))
      }, numeric(length(set$test)))
    }))
  }
  for (refit in c(FALSE, TRUE)) {
    set.seed(3)
    out <- candidate_errors(working, sets, n_simulations, n_candidates, refit)
    expect_equal(out, matrix(direct(refit), ncol = n_simulations),
                 tolerance = 1e-8)
  }
})

test_that("held parameters are known in the simulated datasets", {
  # Every parameter held, so the errors are those of the prediction with them.
  fit <- glgpm(y ~ cov + gp(), data = gaussian_data, family = "gaussian",
               fixed_pars = list(beta = c(1, 0.5), sigma2 = 1, phi = 2,
                                 sigma2_me = 0.1),
               messages = FALSE)
  working <- working_model(fit)
  expect_identical(ncol(working$design), 0L)
  expect_length(working$free, 0L)
  sets <- list(list(test = c(2, 5), left_out = c(2, 5, 7)))
  for (refit in c(FALSE, TRUE)) {
    set.seed(4)
    out <- candidate_errors(working, sets, 3L, 20L, refit)
    set.seed(4)
    V <- working$covariance(working$theta)
    z <- t(chol(V)) %*% matrix(rnorm(3 * fixture_n), fixture_n)
    train <- setdiff(seq_len(fixture_n), c(2, 5, 7))
    weights <- V[c(2, 5), train] %*% solve(V[train, train])
    sd <- sqrt(diag(V[c(2, 5), c(2, 5)] - weights %*% V[train, c(2, 5)]))
    expect_equal(out, pnorm((z[c(2, 5), ] - weights %*% z[train, ]) / sd),
                 tolerance = 1e-8)
  }
})

test_that("a coefficient that a training set can't estimate is left out", {
  # A column that is zero on the training rows of the first test set.
  working <- working_model(gaussian_model)
  sets <- list(list(test = c(2, 5), left_out = c(2, 5, 7)),
               list(test = 1:3, left_out = 1:3))
  wider <- working
  wider$design <- cbind(working$design, as.numeric(seq_len(fixture_n) %in% c(2, 5, 7)))
  set.seed(5)
  out <- candidate_errors(wider, sets, 4L, 6L, TRUE)
  set.seed(5)
  narrower <- candidate_errors(working, sets, 4L, 6L, TRUE)
  expect_equal(out[1:2, ], narrower[1:2, ], tolerance = 1e-8)
  expect_false(isTRUE(all.equal(out[3:5, ], narrower[3:5, ])))
})

test_that("calibration_reference() is reproducible and keeps the session's numbers", {
  inputs <- reference_inputs(gaussian_model, list(1:4, 5:7))
  reference <- function(seeds, predictive = inputs$predictive,
                        test_sets = inputs$test_sets) {
    calibration_reference(gaussian_model, test_sets, predictive, 50L, seeds,
                          n_candidates = 10L)
  }
  set.seed(1)
  state <- .Random.seed
  out <- reference(c(3L, 4L))
  expect_identical(.Random.seed, state)
  expect_identical(dim(out$curves), c(50L, 101L))
  expect_identical(reference(c(3L, 4L)), out)
  expect_false(identical(reference(c(5L, 4L)), out))
  # The PIT values of candidate_errors(), as a curve for each dataset.
  set.seed(3)
  uniforms <- candidate_errors(
    working_model(gaussian_model),
    list(list(test = 1:4, left_out = 1:4), list(test = 5:7, left_out = 5:7)),
    50L, 10L, FALSE
  )
  curves <- pit_curves(t(uniforms), t(uniforms), area_levels())
  expect_equal(out$curves, curves[, seq(1, 1001, by = 10)])
  expect_equal(out$areas, pit_area(curves, area_levels()))
  # A test set without scores isn't in it.
  one <- reference(c(3L, 4L), list(inputs$predictive[[1]], NULL))
  curves <- pit_curves(t(uniforms[1:4, ]), t(uniforms[1:4, ]), area_levels())
  expect_equal(one$curves, curves[, seq(1, 1001, by = 10)])
  # An observation in two test sets has the same PIT values in both, so twice
  # the same test set gives the curves of one.
  twice <- reference_inputs(gaussian_model, list(1:4, 1:4))
  expect_equal(reference(c(3L, 4L), twice$predictive, twice$test_sets)$curves,
               one$curves)
})

test_that("counts have the PIT of counts from their predictive distributions", {
  skip_on_cran()
  fit <- function(family, data, ...) {
    glgpm(
      y ~ gp(), data = data, family = family, control_mcmc = control_mcmc,
      seed = fixture_seed, control_mcml = fixture_control_mcml,
      messages = FALSE, ...
    )
  }
  # Binary outcomes simulated from y ~ gp().
  bernoulli_data <- simulate_fixture(30, function(columns) {
    rbinom(fixture_n, 1, plogis(0.2 + fixture_field))
  }, denominator = rep(1L, fixture_n))
  fits <- list(
    binomial = binomial_model, poisson = poisson_model,
    negative_binomial = fit(negative_binomial(), negative_binomial_data,
                            denominator = denominator,
                            fixed_pars = list(theta = 2)),
    bernoulli = fit("binomial", bernoulli_data, denominator = denominator)
  )
  for (name in names(fits)) {
    set.seed(27)
    inputs <- reference_inputs(fits[[name]], list(c(1, 4, 6), 7:10))
    out <- calibration_reference(fits[[name]], inputs$test_sets,
                                 inputs$predictive, 40L, c(2L, 1L),
                                 n_candidates = 8L)
    set.seed(2)
    uniforms <- candidate_errors(
      working_model(fits[[name]]),
      lapply(inputs$test_sets$sets, function(set) {
        list(test = set$test, left_out = set$test)
      }), 40L, 8L, FALSE
    )
    distributions <- c(inputs$predictive[[1]]$distributions,
                       inputs$predictive[[2]]$distributions)
    lower <- upper <- uniforms
    for (j in seq_along(distributions)) {
      counts <- count_quantile(distributions[[j]], uniforms[j, ])
      # The count whose PIT contains the uniform.
      lower[j, ] <- count_distribution_at(distributions[[j]], counts - 1)
      upper[j, ] <- count_distribution_at(distributions[[j]], counts)
      expect_true(all(lower[j, ] < uniforms[j, ] + 1e-12 &
                        uniforms[j, ] <= upper[j, ] + 1e-12), label = name)
    }
    if (name == "bernoulli") expect_true(all(upper[lower > 0] == 1))
    curves <- pit_curves(t(lower), t(upper), area_levels())
    expect_equal(out$curves, curves[, seq(1, 1001, by = 10)], label = name)
    expect_true(all(is.finite(out$areas) & out$areas > 0), label = name)
  }
})

test_that("counts drawn at known parameters have a calibrated reference", {
  skip_on_cran()
  # With every parameter held the simulated PIT values are uniform, so the
  # counts have their predictive distributions, whose average curve is the
  # identity line.
  held <- function(family, data) {
    glgpm(
      y ~ gp(), data = data, family = family, denominator = denominator,
      fixed_pars = list(beta = 0.2, sigma2 = 0.5, phi = 2),
      control_mcmc = control_mcmc, seed = fixture_seed, messages = FALSE
    )
  }
  for (family in c("binomial", "poisson")) {
    fit <- held(family, if (family == "binomial") binomial_data else poisson_data)
    set.seed(29)
    inputs <- reference_inputs(fit, list(1:5, 6:10))
    out <- calibration_reference(fit, inputs$test_sets, inputs$predictive,
                                 999L, c(3L, 1L), n_candidates = 8L)
    expect_lt(max(abs(colMeans(out$curves) - calibration_levels())), 0.03,
              label = family)
  }
})

test_that("the reference has fewer candidates for many observations", {
  expect_identical(reference_candidates(2000), 200L)
  expect_identical(reference_candidates(2001), 100L)
  expect_no_message(inform_reference_time(1500))
  expect_message(inform_reference_time(1501),
                 "With 1501 observations, simulating the calibration of a correct model can take a minute or more")
})
