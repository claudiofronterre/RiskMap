# The scores of binomial, Poisson and negative binomial predictions, from their
# distribution function averaged over the draws of the mean.

# F at every count, and the scores by their definitions as sums over the counts.
every_count <- function(family, y, units_m, means, parameters, limit) {
  distribution <- switch(
    family,
    binomial = function(k) pbinom(k, units_m, means),
    poisson = function(k) ppois(k, units_m * means),
    negative_binomial = function(k) {
      pnbinom(k, size = parameters$theta, mu = units_m * means)
    }
  )
  values <- vapply(0:limit, function(k) mean(distribution(k)), numeric(1))
  list(values = values, crps = sum((values - (0:limit >= y))^2),
       expected = sum(values * (1 - values)))
}

scaled_crps <- function(crps, expected) {
  0.5 * (1 + crps / expected + log(2 * expected))
}

# A list with a family, its parameters, denominator, observed count and draws of
# the mean.
count_case <- function(family, y, units_m, means, parameters = NULL) {
  list(family = family, y = y, units_m = units_m, means = means,
       parameters = parameters)
}

score_case <- function(case, ...) {
  count_distribution(case$y, case$units_m, case$means, case$parameters,
                     family_traits(case$family), ...)
}

set.seed(20)
small_cases <- list(
  binomial = count_case("binomial", 4, 20, plogis(rnorm(300, -1, 0.8))),
  no_events = count_case("binomial", 0, 20, plogis(rnorm(300, -4, 0.8))),
  one_trial = count_case("binomial", 1, 1, plogis(rnorm(300, -1, 0.8))),
  one_trial_none = count_case("binomial", 0, 1, plogis(rnorm(300, -1, 0.8))),
  poisson = count_case("poisson", 4, 2.5, exp(rnorm(300, 0, 0.4))),
  zero_count = count_case("poisson", 0, 0.3, exp(rnorm(300, 0, 0.4))),
  negative_binomial = count_case("negative_binomial", 3, 0.8,
                                 exp(rnorm(300, 0, 0.4)), list(theta = 4))
)

test_that("scores of small counts are the sums over every count", {
  for (case in small_cases) {
    distribution <- score_case(case)
    limit <- max(distribution$counts)
    expect_lte(limit, 100)
    expect_identical(distribution$counts, as.numeric(-1:limit))
    direct <- every_count(case$family, case$y, case$units_m, case$means,
                          case$parameters, limit)
    expect_equal(distribution$values, c(0, direct$values), tolerance = 1e-14)
    expect_equal(count_crps(case$y, distribution), direct$crps,
                 tolerance = 1e-14)
    expect_equal(count_expected_crps(distribution), direct$expected,
                 tolerance = 1e-14)
    # The counts beyond the largest one have nearly no probability.
    expect_lt(1 - direct$values[limit + 1], 1e-9)
  }
})

test_that("scores of counts agree with sets of simulated counts", {
  skip_on_cran()
  set.seed(21)
  cases <- small_cases[c("binomial", "poisson", "negative_binomial")]
  for (case in cases) {
    traits <- family_traits(case$family)
    distribution <- score_case(case)
    # One count for each draw of the mean, many times, as the scores were
    # computed before. The mean over many sets of the distribution function is
    # F, and the mean of the CRPS exceeds that of F by the variance of the
    # simulated distribution function.
    simulated <- replicate(2000, {
      traits$simulate(case$means, rep(case$units_m, length(case$means)),
                      case$parameters)
    })
    expect_equal(mean(simulated <= case$y),
                 count_distribution_at(distribution, case$y),
                 tolerance = 0.01)
    # The CRPS as E|X - y| - E|X - X'| / 2, and the expected one as the latter.
    other <- replicate(2000, {
      traits$simulate(case$means, rep(case$units_m, length(case$means)),
                      case$parameters)
    })
    # Two independent draws from F have draws of the mean chosen at random.
    spread <- mean(abs(sample(simulated) - sample(other))) / 2
    expect_equal(count_expected_crps(distribution), spread, tolerance = 0.01)
    expect_equal(count_crps(case$y, distribution),
                 mean(abs(simulated - case$y)) - spread, tolerance = 0.02)
  }
})

test_that("large counts are scored as at every count", {
  skip_on_cran()
  set.seed(22)
  draws <- 200
  cases <- list(
    trials = count_case("binomial", 300, 1000, plogis(rnorm(draws, -1, 0.8))),
    many_trials = count_case("binomial", 900, 20000,
                             plogis(rnorm(draws, -3, 0.3))),
    equal_means = count_case("poisson", 9900, 1, rep(1e4, draws)),
    close_means = count_case("poisson", 9900, 2,
                             5000 * exp(rnorm(draws, 0, 0.002))),
    spread_means = count_case("poisson", 700, 1, 1000 * exp(rnorm(draws))),
    # Half of the draws have a mean of 0.1 and the others one of 10^5, and one
    # draw has a mean far from the rest.
    two_means = count_case("poisson", 60, 1,
                           rep(c(0.1, 1e5), each = draws / 2)),
    one_large_mean = count_case("poisson", 1, 1,
                                c(rep(0.1, draws - 1), 1e5)),
    wide_means = count_case("poisson", 5, 1, exp(rnorm(draws, 2, 3))),
    heavy_tail = count_case("negative_binomial", 60, 1,
                            500 * exp(rnorm(draws, 0, 0.5)), list(theta = 0.5))
  )
  for (name in names(cases)) {
    case <- cases[[name]]
    distribution <- score_case(case)
    limit <- max(distribution$counts)
    expect_gt(limit, 900)
    expect_lt(length(distribution$counts), 2000)
    expect_true(all(c(-1:100, case$y - 1, case$y) %in% distribution$counts))
    expect_false(is.unsorted(distribution$values))
    direct <- every_count(case$family, case$y, case$units_m, case$means,
                          case$parameters, limit)
    crps <- count_crps(case$y, distribution)
    expected <- count_expected_crps(distribution)
    expect_equal(crps, direct$crps, tolerance = 5e-5, label = name)
    expect_equal(expected, direct$expected, tolerance = 2e-4, label = name)
    expect_equal(scaled_crps(crps, expected),
                 scaled_crps(direct$crps, direct$expected), tolerance = 5e-5,
                 label = name)
    # F is exact at y - 1 and y, from which the PIT comes.
    expect_equal(count_distribution_at(distribution, case$y - 1:0),
                 direct$values[case$y + 0:1], tolerance = 1e-12)
  }
})

test_that("the sum over a cell is that of the line through its ends", {
  counts <- c(-1, 0, 1, 9, 10, 40)
  values <- c(0, 0.1, 0.15, 0.6, 0.7, 1)
  distribution <- list(counts = counts, values = values)
  line <- approx(counts, values, xout = 0:40)$y
  expect_equal(count_crps(10, distribution),
               sum((line - (0:40 >= 10))^2), tolerance = 1e-14)
  expect_equal(count_expected_crps(distribution), sum(line * (1 - line)),
               tolerance = 1e-14)
  expect_equal(count_distribution_at(distribution, c(-3, -1, 0, 5, 40, 90)),
               c(0, 0, 0.1, 0.375, 1, 1))
  # The smallest whole count that reaches each probability.
  probability <- c(1e-9, 0.1, 0.1 + 1e-9, 0.15, 0.2, 0.6, 0.65, 0.7, 0.713, 1)
  expect_identical(count_quantile(distribution, probability),
                   vapply(probability, function(p) {
                     as.numeric((0:40)[line >= p - 1e-15][1])
                   }, numeric(1)))
  # A probability that the distribution does not reach.
  distribution$values[6] <- 1 - 1e-11
  expect_identical(count_quantile(distribution, 1 - 1e-13), 40)
})

test_that("an observed count beyond the predicted ones is scored", {
  for (case in list(
    count_case("poisson", 60, 1, c(0.4, 0.6)),
    count_case("poisson", 5000, 1, c(0.4, 0.6)),
    count_case("negative_binomial", 500, 1, c(0.4, 0.6), list(theta = 3))
  )) {
    distribution <- score_case(case)
    expect_identical(max(distribution$counts), case$y)
    direct <- every_count(case$family, case$y, case$units_m, case$means,
                          case$parameters, case$y)
    expect_equal(count_crps(case$y, distribution), direct$crps,
                 tolerance = 1e-10)
    expect_gt(direct$crps, case$y - 10)
    expect_equal(count_distribution_at(distribution, case$y - 1), 1)
  }
})

test_that("a mean at the limit of its range is scored", {
  traits <- family_traits("binomial")
  # A prevalence of one, or of zero, in every draw is a point mass, whose CRPS
  # is the absolute error.
  for (trials in c(1, 10, 1e5)) {
    for (prevalence in c(0, 1)) {
      distribution <- count_distribution(1, trials, rep(prevalence, 5), NULL,
                                         traits)
      expect_identical(count_crps(1, distribution),
                       abs(1 - trials * prevalence))
      expect_identical(count_expected_crps(distribution), 0)
    }
  }
  distribution <- count_distribution(2, 4, c(0, 1), NULL, traits)
  expect_identical(distribution$values, c(0, 0.5, 0.5, 0.5, 0.5, 1))
  expect_identical(count_crps(2, distribution), 1)
  zero <- count_distribution(3, 1, c(0, 0), NULL, family_traits("poisson"))
  expect_identical(count_crps(3, zero), 3)
})

test_that("counts too large to be told apart stop, and those below are scored", {
  traits <- family_traits("poisson")
  large <- count_distribution(3, 1, c(1, 1e15), NULL, traits)
  expect_equal(count_crps(3, large), 1e15 / 4, tolerance = 1e-3)
  for (mean in c(1e16, 1e30, 1e180, 1e308)) {
    expect_no_warning(expect_error(
      count_distribution(5, 1, c(1, 1, mean), NULL, traits),
      "too large to be scored", class = "riskmap_error"
    ))
    expect_no_warning(expect_error(
      count_distribution(5, 1, c(1, 1, mean), list(theta = 0.5),
                         family_traits("negative_binomial")),
      "too large to be scored", class = "riskmap_error"
    ))
  }
})

test_that("draws of the mean that can't be scored stop", {
  means <- rbind(c(0.2, 0.3), c(0.4, 0.5), c(0.6, 0.7))
  binomial <- family_traits("binomial")
  poisson <- family_traits("poisson")
  expect_true(check_count_means(means, c(5, 1, 2), binomial))
  expect_true(check_count_means(1000 * means, c(5, 1, 2), poisson))
  for (wrong in c(NaN, NA, Inf, -0.1, 1.1)) {
    means[2, 1] <- wrong
    expect_error(check_count_means(means, c(5, 1, 2), binomial),
                 "1 draw of the mean at 1 test observation isn't finite",
                 class = "riskmap_error")
  }
  means[3, ] <- NaN
  expect_error(check_count_means(means, c(5, 1, 2), binomial),
               "3 draws of the mean at 2 test observations aren't finite")
  # An expected count that overflows.
  expect_error(check_count_means(matrix(1e308, 1, 2), 10, poisson),
               "2 draws of the mean at 1 test observation aren't finite")
  # A denominator that isn't a number of trials.
  for (trials in c(2.5, 250.5)) {
    expect_no_warning(expect_error(
      count_distribution(1, trials, c(0.2, 0.3), NULL, binomial),
      "must be whole numbers", class = "riskmap_error"
    ))
  }
})

# A fit and a mocked prediction whose draws of the linear predictor are given.
# `location` gives the location of each row, one for each by default.
count_fit <- function(family, y, units_m, parameters = NULL,
                      inverse = family_traits(family)$inverse_link,
                      location = seq_along(y)) {
  n <- length(y)
  data <- st_sf(y = y, x1 = seq_len(n),
                geometry = st_sfc(lapply(location, function(i) {
                  st_point(c(1000 * i, 0))
                }), crs = 32637))
  fit <- list(
    formula = y ~ x1, data = data, input_crs = st_crs(data), family = family,
    estimate = list(beta = c(0, 0), sigma2 = 0, phi = 0),
    D = matrix(1, nrow = n, ncol = 2), re = list(), cov_offset = NULL,
    nugget = "none", fixed_pars = parameters,
    sst = FALSE, units_m = units_m, y = y,
    ID_coords = match(location, unique(location)),
    coords = unique(st_coordinates(data)), distance_units = "m",
    call = list(denominator = quote(units_m)),
    link_function = list(inv = inverse, name = "custom")
  )
  class(fit) <- "RiskMap"
  fit
}

mock_draws <- function(draws) {
  function(fit, newdata, control_mcmc, seed) draws[newdata$x1, , drop = FALSE]
}

test_that("assess_prediction() scores counts by the average distribution function", {
  set.seed(23)
  draws <- matrix(rnorm(4 * 500, -0.5, 0.7), 4)
  local_mocked_bindings(predict_test_set = mock_draws(draws),
                        calibration_reference = independent_reference)
  cases <- list(
    binomial = list(units_m = c(10, 1, 1, 30), y = c(3, 1, 0, 30),
                    parameters = NULL, mean = plogis),
    poisson = list(units_m = c(1, 2, 0.5, 4), y = c(0, 9, 1, 3),
                   parameters = NULL, mean = exp),
    negative_binomial = list(units_m = c(1, 2, 0.5, 4), y = c(0, 9, 1, 3),
                             parameters = list(theta = 1.5), mean = exp)
  )
  for (family in names(cases)) {
    case <- cases[[family]]
    fit <- count_fit(family, case$y, case$units_m, case$parameters)
    out <- assess_prediction(list(m = fit), list(2:4), messages = FALSE)
    direct <- lapply(2:4, function(i) {
      every_count(family, case$y[i], case$units_m[i], case$mean(draws[i, ]),
                  case$parameters,
                  if (family == "binomial") case$units_m[i] else 2000)
    })
    crps <- vapply(direct, `[[`, numeric(1), "crps")
    expected <- vapply(direct, `[[`, numeric(1), "expected")
    expect_equal(out$scores$CRPS, crps, tolerance = 1e-7, label = family)
    expect_equal(out$scores$SCRPS, scaled_crps(crps, expected),
                 tolerance = 1e-5, label = family)
    ends <- vapply(1:3, function(i) {
      c(0, direct[[i]]$values)[case$y[i + 1] + 1:2]
    }, numeric(2))
    expect_equal(out$scores$pit_lower, ends[1, ], tolerance = 1e-6,
                 label = family)
    expect_equal(out$scores$pit_upper, ends[2, ], tolerance = 1e-6,
                 label = family)
    curve <- function(u) {
      rowMeans(apply(ends, 2, function(end) {
        pmin(pmax((u - end[1]) / diff(end), 0), 1)
      }))
    }
    curves <- out$calibration$curves
    expect_equal(curves$curve, curve(curves$u), tolerance = 1e-6,
                 label = family)
    u <- seq(0, 1, length.out = 1001)
    expect_equal(out$calibration$area$PIT_area, pit_area(t(curve(u)), u),
                 tolerance = 1e-6, label = family)
  }
})

test_that("the CRPS of a binary outcome is the Brier score", {
  set.seed(28)
  means <- plogis(matrix(rnorm(5 * 300), 5))
  y <- c(0, 1, 1, 0, 1)
  scored <- score_counts(y, rep(1, 5), means, NULL, family_traits("binomial"))
  expect_equal(scored$scores$CRPS, (rowMeans(means) - y)^2)
})

test_that("a rare event has a SCRPS", {
  # The prevalence is 10^-6 in every draw, so a count simulated for each draw
  # is 0 in every one of them, a distribution without spread.
  local_mocked_bindings(
    predict_test_set = mock_draws(matrix(qlogis(1e-6), 2, 100)),
    calibration_reference = independent_reference
  )
  out <- assess_prediction(count_fit("binomial", c(0, 0), c(10, 10)), list(2),
                           messages = FALSE)
  values <- pbinom(0:10, 10, 1e-6)
  crps <- sum((values - 1)^2)
  expected <- sum(values * (1 - values))
  expect_equal(out$scores$CRPS, crps)
  expect_equal(out$scores$SCRPS, scaled_crps(crps, expected))
})

test_that("a test set with a draw of the mean that isn't finite is a failure", {
  draws <- matrix(0, 3, 50)
  draws[2, 7] <- 800
  # exp(800) is infinite, and the custom link returns no number above 100,
  # draws that a simulated count would have left out.
  broken <- function(eta) ifelse(eta > 100, NaN, plogis(eta))
  fits <- list(poisson = count_fit("poisson", c(1, 2, 3), rep(1, 3)),
               negative_binomial = count_fit("negative_binomial", c(1, 2, 3),
                                             rep(1, 3), list(theta = 2)),
               binomial = count_fit("binomial", c(1, 2, 3), rep(5, 3),
                                    inverse = broken))
  local_mocked_bindings(predict_test_set = mock_draws(draws),
                        calibration_reference = independent_reference)
  for (name in names(fits)) {
    expect_no_warning(expect_warning(
      out <- assess_prediction(fits[name], list(1:2, 3), messages = FALSE),
      "1 of 2 test-set assessments failed"
    ))
    expect_identical(out$failures$test_set, 1L)
    expect_match(out$failures$message,
                 "1 draw of the mean at 1 test observation isn't finite")
    expect_identical(is.finite(out$scores$CRPS), c(FALSE, FALSE, TRUE))
    expect_true(is.finite(out$calibration$area$PIT_area))
  }
})

test_that("a fit has the same count scores alone and in any position of a list", {
  skip_on_cran()
  fit_family <- function(family, data, ...) {
    suppressWarnings(glgpm(
      y ~ gp(), data = data, family = family, denominator = denominator,
      control_mcmc = control_mcmc, seed = fixture_seed,
      control_mcml = fixture_control_mcml, messages = FALSE, ...
    ))
  }
  fits <- list(
    binomial = fit_family("binomial", binomial_data),
    poisson = fit_family("poisson", poisson_data),
    negative_binomial = fit_family(negative_binomial(), poisson_data,
                                   fixed_pars = list(theta = 2))
  )
  assess <- function(models, seed = fixture_seed) {
    assess_prediction(models, list(1:3, c(2, 5:7)),
                      control_mcmc = control_mcmc, seed = seed,
                      messages = FALSE)
  }
  for (name in names(fits)) {
    fit <- fits[[name]]
    set.seed(24)
    state <- .Random.seed
    pair <- assess(list(a = fit, b = fit))
    expect_identical(.Random.seed, state)
    scores <- split(pair$scores[-1], pair$scores$model)
    expect_equal(scores$a, scores$b, ignore_attr = TRUE, tolerance = 0,
                 label = name)
    simulated <- pair$calibration$simulated
    expect_identical(simulated$a, simulated$b, label = name)
    alone <- assess(list(a = fit))
    expect_identical(alone$scores[-1], scores$a, label = name)
    expect_identical(alone$calibration$simulated$a, simulated$a, label = name)
    expect_identical(nrow(alone$failures), 0L)
    # Without a seed only the seeds of the test sets come from the session, so
    # its state does not depend on the number of models.
    set.seed(24)
    assess(list(a = fit), seed = NULL)
    state <- .Random.seed
    set.seed(24)
    assess(list(a = fit, b = fit), seed = NULL)
    expect_identical(.Random.seed, state)
  }
})

test_that("large counts, scored on a grid, have the PIT of every count", {
  skip_on_cran()
  set.seed(26)
  means <- 400 * exp(matrix(rnorm(40 * 100, 0, 0.5), 40))
  observed <- rpois(40, means[, 1])
  fit <- count_fit("poisson", c(observed, 1), rep(1, 41))
  local_mocked_bindings(predict_test_set = mock_draws(log(means)),
                        calibration_reference = independent_reference)
  out <- assess_prediction(fit, list(1:40), seed = 1, messages = FALSE)
  exact <- function(shift) vapply(1:40, function(i) {
    mean(ppois(observed[i] + shift, means[i, ]))
  }, numeric(1))
  expect_equal(out$scores$pit_lower, exact(-1), tolerance = 1e-3)
  expect_equal(out$scores$pit_upper, exact(0), tolerance = 1e-3)
})

# Binary outcomes ------------------------------------------------------------

binary_y <- c(0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1)
binary_sets <- list(1:4, 5:8, 9:12)
binary_draws <- function() {
  withr::with_seed(29, matrix(rnorm(12 * 300, -0.3, 0.8), 12))
}

test_that("binary models are told apart by their scores", {
  local_mocked_bindings(predict_test_set = mock_draws(binary_draws()),
                        calibration_reference = independent_reference)
  binary <- function(fit, test_sets = binary_sets) {
    binary_models(assess_prediction(list(m = fit), test_sets, messages = FALSE))
  }
  expect_identical(binary(count_fit("binomial", binary_y, rep(1, 12))),
                   c(m = TRUE))
  # Individual-level data, with three rows at each location.
  rows <- count_fit("binomial", binary_y, rep(1, 12),
                    location = rep(1:4, each = 3))
  expect_identical(binary(rows, create_test_sets(rows)), c(m = TRUE))
  # Several trials, even when every count is 0 or the number of trials, or in
  # one observation only.
  expect_identical(binary(count_fit("binomial", 3 * binary_y, rep(3, 12))),
                   c(m = FALSE))
  expect_identical(binary(count_fit("binomial", rep(0, 12), rep(2, 12))),
                   c(m = FALSE))
  expect_identical(binary(count_fit("binomial", binary_y, c(2, rep(1, 11)))),
                   c(m = FALSE))
  expect_identical(binary(count_fit("poisson", binary_y, rep(1, 12))),
                   c(m = FALSE))
})

test_that("assess_prediction() says once that the area of binary outcomes is not a check", {
  local_mocked_bindings(predict_test_set = mock_draws(binary_draws()),
                        calibration_reference = independent_reference)
  bernoulli <- count_fit("binomial", binary_y, rep(1, 12))
  said <- function(models) {
    capture_messages(assess_prediction(models, binary_sets))
  }
  messages <- said(list(m = bernoulli))
  expect_length(grep("binary outcomes, the PIT area", messages), 1L)
  expect_length(grep("Brier score", messages), 1L)
  expect_no_message(assess_prediction(bernoulli, binary_sets, messages = FALSE))
  # Each model of a mixed assessment on its own.
  mixed <- said(list(b = bernoulli, p = count_fit("poisson", binary_y,
                                                  rep(1, 12))))
  expect_length(grep("binary outcomes \\(.b.\\), the PIT area", mixed), 1L)
  expect_no_match(said(list(m = count_fit("binomial", 3 * binary_y,
                                          rep(3, 12)))),
                  "binary|Brier")
})

test_that("the methods of a cross-validation give no limit or band for binary outcomes", {
  local_mocked_bindings(predict_test_set = mock_draws(binary_draws()),
                        calibration_reference = independent_reference)
  models <- list(b = count_fit("binomial", binary_y, rep(1, 12)),
                 p = count_fit("poisson", binary_y, rep(1, 12)))
  out <- assess_prediction(models, binary_sets, messages = FALSE)
  # The result keeps the limit and the band.
  expect_true(all(is.finite(out$calibration$area$PIT_area_limit)))
  expect_true(all(is.finite(out$calibration$curves$upper)))
  scores <- summary(out)
  expect_identical(attr(scores, "binary"), c(b = TRUE, p = FALSE))
  expect_identical(scores$PIT_area_limit, out$calibration$area$PIT_area_limit)
  limits <- cv_summary_table(scores, 3L)$Limit
  expect_identical(limits[1], "none (binary)")
  expect_match(limits[2], "^0\\.[0-9]{3}$")
  printed <- squish(capture.output(print(out)))
  expect_match(printed, "none (binary)", fixed = TRUE)
  expect_match(printed, "binary outcomes (b)", fixed = TRUE)
  expect_match(squish(to_table(out)), "none (binary)", fixed = TRUE)
  geoms <- function(plot) {
    vapply(plot$layers, function(layer) class(layer$geom)[1], character(1))
  }
  p <- plot_cross_validation(out)
  band <- ggplot2::layer_data(p, match("GeomRibbon", geoms(p)))
  expect_identical(unique(as.integer(band$PANEL)), 2L)
  label <- ggplot2::layer_data(p, match("GeomText", geoms(p)))$label
  expect_no_match(label[1], "Below")
  expect_match(label[2], "Below")
  alone <- plot_cross_validation(out, model = "b")
  expect_false("GeomRibbon" %in% geoms(alone))
})

test_that("a binomial model with several trials keeps its limit and band", {
  local_mocked_bindings(predict_test_set = mock_draws(binary_draws()),
                        calibration_reference = independent_reference)
  fit <- count_fit("binomial", 3 * binary_y, rep(3, 12))
  out <- assess_prediction(list(m = fit), binary_sets, messages = FALSE)
  scores <- summary(out)
  expect_identical(attr(scores, "binary"), c(m = FALSE))
  expect_match(cv_summary_table(scores, 3L)$Limit, "^0\\.[0-9]{3}$")
  expect_match(squish(capture.output(print(scores))),
               "below the limit with a probability of 95%", fixed = TRUE)
  for (style in c("curve", "histogram")) {
    p <- plot_cross_validation(out, style = style)
    layers <- vapply(p$layers, function(layer) class(layer$geom)[1],
                     character(1))
    expect_true("GeomRibbon" %in% layers, label = style)
    expect_match(ggplot2::layer_data(p, match("GeomText", layers))$label,
                 "Below", label = style)
  }
})
