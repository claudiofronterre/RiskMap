# assess_fit() and the methods of its result.

# Data with two records at each location, and fits with every parameter held at
# given values, which need no estimation.
check_data_twice <- gaussian_data[rep(1:10, each = 2), ]
check_data_twice$cov <- check_data_twice$cov + rep(c(0, 1), 10)
set.seed(31)
check_data_twice$y <- 1 + rnorm(20)
check_data_twice$count <- rpois(20, 2)
check_data_twice$binary <- rbinom(20, 1, 0.5)
check_data_twice$z <- rnorm(20)
check_data_twice$exposure <- rep(c(0.5, 2), 10)

held_fit <- function(formula, family, parameters, denominator = NULL) {
  arguments <- list(
    formula, data = check_data_twice, family = family,
    fixed_pars = c(list(sigma2 = 0.5, phi = 2), parameters),
    control_mcmc = control_mcmc, seed = 1, messages = FALSE
  )
  if (!is.null(denominator)) arguments$denominator <- as.name(denominator)
  suppressWarnings(do.call(glgpm, arguments))
}
gaussian_twice <- held_fit(y ~ cov + gp(), "gaussian",
                           parameters = list(beta = c(1, 0.2),
                                             sigma2_me = 0.5))
binary_twice <- held_fit(binary ~ gp(), "binomial",
                         parameters = list(beta = 0))
poisson_twice <- held_fit(count ~ cov + gp(), "poisson", denominator = "exposure",
                          parameters = list(beta = c(0.5, 0)))
negative_binomial_twice <- held_fit(count ~ gp(), negative_binomial(),
                                    denominator = "exposure",
                                    parameters = list(beta = 0.5, theta = 2))

checked <- assess_fit(list(a = gaussian_model, b = gaussian_intercept_model),
                      nsim = 99, seed = 7, messages = FALSE)

test_that("assess_fit() describes the observed and the simulated data", {
  expect_s3_class(checked, "RiskMap_predictive_check")
  expect_identical(names(checked$family), c("a", "b"))
  expect_identical(checked$n_locations, 10L)
  expect_identical(checked$observed$statistic,
                   rep(c("mean", "sd", "min", "max"), 2))
  y <- gaussian_data$y
  expect_equal(checked$observed$value, rep(c(mean(y), sd(y), range(y)), 2))
  expect_identical(nrow(checked$simulated), 2L * 4L * 99L)

  # The simulated statistics are those of the datasets of simulate_glgpm().
  simulations <- simulate_glgpm(gaussian_model, nsim = 99, seed = 7)
  fifth <- simulated_data(simulations, simulation = 5)$y
  expect_equal(
    subset(checked$simulated, model == "a" & simulation == 5)$value,
    c(mean(fifth), sd(fifth), range(fifth))
  )
})

test_that("assess_fit() takes the variogram of the residuals from the covariates", {
  shown <- subset(checked$variogram, model == "a")
  residuals <- resid(lm(y ~ cov, gaussian_data))
  data <- gaussian_data
  data$residuals <- residuals
  expected <- variogram(data, "residuals", n_permutations = 0)
  kept <- expected$variogram$n_pairs > 0
  expect_equal(checked$breaks, expected$breaks)
  expect_equal(shown$distance, expected$variogram$distance[kept])
  expect_equal(shown$n_pairs, expected$variogram$n_pairs[kept])
  expect_equal(shown$observed, expected$variogram$semivariance[kept])
  # Without covariates it is that of the response.
  expect_equal(
    subset(checked$variogram, model == "b")$observed,
    variogram(gaussian_data, "y", n_permutations = 0)$variogram$semivariance[kept]
  )

  # The regression is done again in each simulated dataset, and the envelope
  # is the one of variogram(), with the observed curve first.
  simulations <- simulate_glgpm(gaussian_model, nsim = 99, seed = 7)
  curves <- vapply(0:99, function(k) {
    y <- if (k == 0) gaussian_data$y else simulated_data(simulations, k)$y
    data$residuals <- resid(lm(y ~ gaussian_data$cov))
    variogram(data, "residuals", n_permutations = 0)$variogram$semivariance[kept]
  }, numeric(sum(kept)))
  envelope <- global_rank_envelope(curves, 0.95)
  expect_equal(shown$lower, envelope$lower)
  expect_equal(shown$upper, envelope$upper)
  expect_equal(shown$median, apply(curves[, -1], 1, median))
  expect_identical(shown$outside, shown$observed < shown$lower |
                     shown$observed > shown$upper)
})

test_that("the offset is taken out of the residuals", {
  with_offset <- assess_fit(gaussian_offset_model, nsim = 19, seed = 1)
  data <- gaussian_data
  data$residuals <- resid(lm(y - offset ~ cov, gaussian_data))
  expected <- variogram(data, "residuals", n_permutations = 0)$variogram
  expect_equal(with_offset$variogram$observed,
               expected$semivariance[expected$n_pairs > 0])
})

test_that("observations at a location are added up, or averaged", {
  index <- rep(1:10, each = 2)
  locations <- check_data_twice[2 * (1:10), ]
  semivariances <- function(values) {
    locations$values <- values
    variogram(locations, "values", n_permutations = 0)$variogram |>
      subset(n_pairs > 0, select = "semivariance", drop = TRUE)
  }
  statistics <- function(check) {
    stats::setNames(check$observed$value, check$observed$statistic)
  }
  covariate <- as.numeric(tapply(check_data_twice$cov, index, mean))

  # Gaussian, the statistics and the distribution are of the observations.
  check <- assess_fit(gaussian_twice, nsim = 19, seed = 1)
  y <- check_data_twice$y
  expect_equal(statistics(check), c(mean = mean(y), sd = sd(y), min = min(y),
                                    max = max(y)))
  expect_equal(check$variogram$observed,
               semivariances(resid(lm(as.numeric(tapply(y, index, mean)) ~
                                        covariate))))
  expect_identical(check$n_observations, 20L)
  expect_identical(check$n_locations, 10L)

  # Denominators of 1. The positives at a location are added up.
  check <- assess_fit(binary_twice, nsim = 19, seed = 1)
  positive <- as.numeric(tapply(check_data_twice$binary, index, sum))
  expect_equal(statistics(check),
               c(mean_prevalence = mean(positive / 2),
                 sd_prevalence = sd(positive / 2), zeros = mean(positive == 0)))
  expect_equal(check$variogram$observed, semivariances(elogit(positive, 2)))
  expect_identical(check$quantity, "prevalence")
  expect_true(all(check$distribution$value %in% c(0, 0.5, 1)))

  # Counts with exposures. The covariates are weighted by them.
  exposure <- as.numeric(tapply(check_data_twice$exposure, index, sum))
  count <- as.numeric(tapply(check_data_twice$count, index, sum))
  weighted <- as.numeric(tapply(check_data_twice$cov * check_data_twice$exposure,
                                index, sum)) / exposure
  log_rate <- log((count + 0.5) / exposure)
  check <- assess_fit(list(poisson = poisson_twice,
                           negbin = negative_binomial_twice),
                      nsim = 19, seed = 1)
  expect_identical(check$quantity, "rate")
  # The statistics and the distribution are of the observations.
  rate <- check_data_twice$count / check_data_twice$exposure
  expect_equal(subset(check$observed, model == "negbin")$value,
               c(mean(rate), sd(rate), mean(rate == 0)))
  shown <- subset(check$distribution, model == "negbin")
  expect_equal(shown$observed, ecdf(rate)(shown$value))
  expect_equal(subset(check$variogram, model == "negbin")$observed,
               semivariances(log_rate))
  expect_equal(subset(check$variogram, model == "poisson")$observed,
               semivariances(resid(lm(log_rate ~ weighted))))
})

test_that("the distribution is the empirical one, with its joint envelope", {
  shown <- subset(checked$distribution, model == "b")
  expect_equal(shown$observed, ecdf(gaussian_data$y)(shown$value))
  expect_identical(min(shown$observed), as.numeric(min(shown$value) >=
                                                     min(gaussian_data$y)) * 0.1)
  expect_identical(max(shown$observed), 1)
  expect_false(is.unsorted(shown$value, strictly = TRUE))
  simulations <- simulate_glgpm(gaussian_intercept_model, nsim = 99, seed = 7)
  curves <- vapply(0:99, function(k) {
    y <- if (k == 0) gaussian_data$y else simulated_data(simulations, k)$y
    ecdf(y)(shown$value)
  }, numeric(nrow(shown)))
  envelope <- global_rank_envelope(curves, 0.95)
  expect_equal(shown$lower, envelope$lower)
  expect_equal(shown$upper, envelope$upper)
  # Models share their values.
  expect_identical(subset(checked$distribution, model == "a")$value, shown$value)
})

test_that("many values are summarised by quantiles, few are kept", {
  set.seed(4)
  values <- list(a = matrix(rnorm(5000), 50), b = matrix(rnorm(5000, 0, 3), 50))
  at <- distribution_values(values)
  expect_lte(length(at), 202L)
  expect_true(all(range(values$a[, 1]) >= min(at) & range(values$a[, 1]) <= max(at)))
  expect_identical(distribution_values(list(matrix(c(0, 1, 0, 1, 1, 0), 2))),
                   c(0, 1))
})

test_that("statistics can be a function or a named list of functions", {
  describe <- function(data) {
    c(mean = mean(data$y), correlation = cor(data$y, data$cov))
  }
  check <- assess_fit(gaussian_model, describe, nsim = 19, seed = 7)
  expect_equal(check$observed$value,
               c(mean(gaussian_data$y), cor(gaussian_data$y, gaussian_data$cov)))
  simulated <- simulated_data(simulate_glgpm(gaussian_model, nsim = 19, seed = 7),
                              simulation = 19)
  expect_equal(subset(check$simulated, simulation == 19)$value,
               c(mean(simulated$y), cor(simulated$y, simulated$cov)))
  expect_identical(check_statistic_labels("mean", check), "mean")

  # A list, a data frame and a list of functions.
  same <- list(function(data) as.list(describe(data)),
               function(data) as.data.frame(as.list(describe(data))),
               list(mean = function(data) mean(data$y),
                    correlation = function(data) c(rho = cor(data$y, data$cov))))
  for (statistics in same) {
    other <- assess_fit(gaussian_model, statistics, nsim = 19, seed = 7)
    expect_equal(other$simulated, check$simulated)
  }
  # One number without a name, and a logical one.
  one <- assess_fit(gaussian_model, function(data) all(data$y > 0), nsim = 19,
                    seed = 7)
  expect_identical(one$observed$statistic, "statistic")
  expect_identical(one$observed$value, 1)
})

test_that("statistics that can't be used are errors", {
  run <- function(statistics) {
    assess_fit(gaussian_model, statistics, nsim = 19, seed = 7)
  }
  expect_input_error(run("mean"), "statistics", "a function, or a list")
  expect_input_error(run(list(mean, sd)), "statistics", "distinct names")
  expect_input_error(run(list()), "statistics")
  expect_input_error(run(function(data) c(1, 2)), "statistics", "named numbers")
  expect_input_error(run(function(data) c(a = "x")), "statistics", "named numbers")
  expect_input_error(run(function(data) data.frame(a = 1:2)), "statistics")
  expect_input_error(run(function(data) c(a = 1, a = 2)), "statistics", "repeated")
  expect_input_error(run(function(data) NULL), "statistics")
  failed <- expect_input_error(run(function(data) stop("no luck")), "statistics",
                               "failed on the observed data")
  expect_match(conditionMessage(failed$parent), "no luck")
  # Other names for a simulated dataset.
  first <- TRUE
  changing <- function(data) {
    names <- if (first) "a" else "b"
    first <<- FALSE
    stats::setNames(1, names)
  }
  expect_input_error(run(list(changing = function(data) c(1, changing(data)))),
                     "statistics")
  first <- TRUE
  expect_input_error(run(function(data) c(changing(data), other = 1)),
                     "statistics", "same named numbers")
})

test_that("a transformed response needs the default statistics", {
  data <- gaussian_data
  data$positive <- exp(data$y)
  fit <- glgpm(log(positive) ~ gp(), data = data, family = "gaussian",
               messages = FALSE)
  check <- assess_fit(fit, nsim = 19, seed = 1)
  expect_identical(check$response, "log(positive)")
  expect_equal(check$observed$value[1], mean(data$y))
  expect_input_error(
    assess_fit(fit, function(data) c(a = mean(data$positive)), nsim = 19),
    "models", "plain column name"
  )
})

test_that("statistics that aren't finite are left out and reported", {
  set.seed(1)
  uncertain <- function(data) {
    c(sd = sd(data$y), log_min = suppressWarnings(log(min(data$y))))
  }
  expect_warning(
    check <- assess_fit(gaussian_model, uncertain, nsim = 99, seed = 7),
    "left out"
  )
  summarised <- summary(check)
  simulated <- subset(check$simulated, statistic == "log_min")$value
  finite <- simulated[is.finite(simulated)]
  expect_gt(length(finite), 0)
  expect_lt(length(finite), 99)
  row <- subset(as.data.frame(unclass(summarised)), statistic == "log_min")
  expect_identical(row$n_not_finite, 99L - length(finite))
  expect_equal(row$median, median(finite))
  expect_equal(row$at_or_above, mean(finite >= row$observed))
  expect_identical(summarised$n_not_finite[1], 0L)
  expect_output(print(check), "aren't finite are left out")

  # A statistic that is never finite.
  expect_warning(
    expect_warning(never <- assess_fit(gaussian_model, function(data) c(x = NaN),
                                       nsim = 19, seed = 7), "observed data"),
    "simulated datasets"
  )
  expect_true(all(is.na(summary(never)[c("observed", "median", "at_or_above")])))
  expect_output(print(never), "NA")
  expect_input_error(plot_predictive_check(never), "x", "no finite statistics")
  # A statistic that can be drawn is, without the one that can't.
  some <- suppressWarnings(
    assess_fit(gaussian_model, function(data) c(x = NaN, mean = mean(data$y)),
               nsim = 19, seed = 7)
  )
  expect_identical(unique(as.character(plot_predictive_check(some)$data$statistic)),
                   "mean")
})

test_that("the warnings of statistics are raised once for each model", {
  noisy <- function(data) {
    warning("careful")
    c(mean = mean(data$y))
  }
  warnings <- capture_warnings(
    assess_fit(list(a = gaussian_model, b = gaussian_intercept_model), noisy,
               nsim = 19, seed = 7)
  )
  expect_length(warnings, 2)
  expect_match(warnings, "gave 20 warnings")
})

test_that("summary() and to_table() say where the observed values lie", {
  summarised <- summary(checked)
  expect_s3_class(summarised, "summary.RiskMap_predictive_check")
  expect_named(summarised, c("model", "statistic", "observed", "median", "lower",
                             "upper", "at_or_above", "n_not_finite"))
  expect_identical(summarised$model, rep(c("a", "b"), 4))
  expect_identical(summarised$statistic,
                   rep(c("mean", "sd", "min", "max"), each = 2))
  simulated <- subset(checked$simulated, model == "b" & statistic == "sd")$value
  row <- as.data.frame(unclass(summarised))[4, ]
  expect_equal(unlist(row[c("median", "lower", "upper")], use.names = FALSE),
               quantile(simulated, c(0.5, 0.025, 0.975), names = FALSE))
  expect_equal(row$at_or_above, mean(simulated >= sd(gaussian_data$y)))

  expect_s3_class(to_table(checked), "knitr_kable")
  expect_identical(to_table(checked), to_table(summarised))
  expect_match(to_table(checked), "95% interval", all = FALSE)
  expect_input_error(to_table(summarised[1:3]), "object")
  expect_input_error(to_table(checked, digits = -1), "digits")
  expect_error(summary(checked, level = 0.9), class = "riskmap_error_input")
})

test_that("a value equal to the observed one counts as at or above it", {
  row <- summarise_statistic(0.3, c(0.1 + 0.2, 0.3, 0.2, 0.5))
  expect_identical(row$at_or_above, 0.75)
})

test_that("an assessment of fits prints", {
  # With held parameters the numbers don't depend on an optimiser.
  constant <- held_fit(y ~ gp(), "gaussian",
                       parameters = list(beta = 1, sigma2_me = 0.5))
  two <- assess_fit(list(constant = constant, line = gaussian_twice), nsim = 99,
                    seed = 1)
  expect_snapshot(print(two))
  expect_snapshot(print(summary(assess_fit(constant, nsim = 19, seed = 1))))
  # Some columns of the summary print as a data frame.
  expect_output(print(summary(checked)[1:3]), "statistic +observed")
  expect_error(print(checked, digits = 2), class = "riskmap_error_input")
})

test_that("the seed sets the simulations, and leaves the caller's state", {
  set.seed(11)
  before <- .Random.seed
  again <- assess_fit(list(a = gaussian_model, b = gaussian_intercept_model),
                      nsim = 99, seed = 7)
  expect_identical(.Random.seed, before)
  expect_identical(again, checked)

  # The draws of a model don't depend on the other models.
  alone <- assess_fit(list(b = gaussian_intercept_model), nsim = 99, seed = 7)
  expect_identical(alone$simulated$value,
                   subset(checked$simulated, model == "b")$value)

  # Statistics with random numbers have the same ones for every model.
  random <- function(data) c(u = runif(1))
  check <- assess_fit(list(a = gaussian_model, b = gaussian_intercept_model),
                      random, nsim = 19, seed = 7)
  expect_identical(subset(check$simulated, model == "a")$value,
                   subset(check$simulated, model == "b")$value)
  expect_identical(.Random.seed, before)

  # Without a seed the simulations come from the caller's state.
  first <- assess_fit(gaussian_model, nsim = 19)
  second <- assess_fit(gaussian_model, nsim = 19)
  expect_false(identical(first$simulated, second$simulated))
  set.seed(11)
  expect_identical(assess_fit(gaussian_model, nsim = 19)$simulated,
                   first$simulated)
})

test_that("the simulations are returned only if asked for", {
  expect_null(checked$simulations)
  check <- assess_fit(list(a = gaussian_model), nsim = 19, seed = 7,
                      return_simulations = TRUE)
  expect_named(check$simulations, "a")
  expect_equal(check$simulations$a$samples,
               simulate_glgpm(gaussian_model, nsim = 19, seed = 7)$samples)
})

test_that("the distance classes can be given", {
  check <- assess_fit(gaussian_model, breaks = c(0, 3, 6, 9), nsim = 19, seed = 1)
  expect_identical(check$breaks, c(0, 3, 6, 9))
  expect_identical(nrow(check$variogram), 3L)
  run <- function(breaks) {
    assess_fit(gaussian_model, breaks = breaks, nsim = 19, seed = 1)
  }
  expect_input_error(run(c(1, 0)), "breaks", "increasing")
  expect_input_error(run("a"), "breaks")
  expect_input_error(run(c(-1, 1)), "breaks", "non-negative")
  expect_input_error(run(c(500, 600)), "breaks", "No pair")
  expect_warning(run(c(0, 5, 500, 600)), "no pairs")
})

test_that("few locations have no variogram, or one class of it", {
  few <- gaussian_data[c(1, 2), ]
  fit <- gaussian_intercept_model
  fit$data <- fit$data[1:2, ]
  fit$y <- fit$y[1:2]
  fit$D <- fit$D[1:2, , drop = FALSE]
  fit$cov_offset <- fit$cov_offset[1:2]
  expect_warning(check <- assess_fit(fit, nsim = 19, seed = 1), "no variogram")
  expect_null(check$variogram)
  expect_output(print(check), "2\\s+locations")
  expect_input_error(plot_predictive_check(check, type = "variogram"), "x",
                     "no variogram")
  one_class <- assess_fit(fit, breaks = c(0, 1e3), nsim = 19, seed = 1)
  expect_identical(nrow(one_class$variogram), 1L)
  expect_no_message(ggplot2::ggplotGrob(
    plot_predictive_check(one_class, type = "variogram")
  ))
  # As many coefficients as locations leave no residuals.
  fit$D <- cbind(1, c(0, 1))
  none <- assess_fit(fit, breaks = c(0, 1e3), nsim = 19, seed = 1)
  expect_identical(none$variogram$observed, 0)
  expect_false(none$variogram$outside)
})

test_that("a location without anyone examined is left out of rates", {
  fit <- binary_twice
  fit$units_m[1:2] <- 0
  fit$y[1:2] <- 0
  check <- assess_fit(fit, nsim = 19, seed = 1)
  positive <- as.numeric(tapply(fit$y[-(1:2)], rep(1:9, each = 2), sum))
  expect_equal(check$observed$value[1], mean(positive / 2))
  expect_true(all(is.finite(check$simulated$value)))
  expect_true(all(is.finite(check$variogram$observed)))
})

test_that("assess_fit() checks its arguments", {
  expect_input_error(assess_fit(), "models")
  expect_input_error(assess_fit(1), "models")
  expect_input_error(assess_fit(list(gaussian_model, binomial_model)), "models")
  expect_input_error(assess_fit(gaussian_model, nsim = 18), "nsim", "at least 19")
  expect_input_error(assess_fit(gaussian_model, nsim = 1.5), "nsim")
  expect_input_error(assess_fit(gaussian_model, seed = -1), "seed")
  expect_input_error(assess_fit(gaussian_model, return_simulations = NA),
                     "return_simulations")
  expect_input_error(assess_fit(gaussian_model, messages = 1), "messages")
  # Families of different quantities, for the same counts.
  counts <- check_data_twice
  fit <- poisson_twice
  fit$family <- "gaussian"
  fit$units_m <- poisson_twice$units_m
  expect_input_error(assess_fit(list(a = poisson_twice, b = fit)), "models",
                     "share one distribution")
})

test_that("plot_predictive_check() draws the statistics", {
  p <- plot_predictive_check(checked)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetGrid")
  expect_identical(levels(p$data$statistic), c("Mean", "SD", "Minimum", "Maximum"))
  # Every simulated dataset is in a bar of its panel, and the darker bars are
  # those at or above the observed value.
  panel <- paste(p$data$model, p$data$statistic)
  expect_true(all(tapply(p$data$count, panel, sum) == 99))
  summarised <- summary(checked)
  expect_setequal(ggplot2::layer_data(p, 2)$xintercept, summarised$observed)
  above <- tapply(p$data$count * p$data$above, panel, sum)
  expect_equal(
    as.numeric(above[paste(summarised$model,
                           check_statistic_labels(summarised$statistic, checked))]),
    99 * summarised$at_or_above
  )
  expect_length(p$layers, 2L)
  expect_identical(ggplot2::layer_data(plot(checked), 1), ggplot2::layer_data(p, 1))
  one <- plot_predictive_check(checked, model = "b")
  expect_s3_class(one$facet, "FacetWrap")
  expect_identical(levels(one$data$model), "b")
})

test_that("plot_predictive_check() draws the distribution and the variogram", {
  p <- plot_predictive_check(checked, type = "variogram")
  points <- ggplot2::layer_data(p, 3)
  expect_equal(points$y, checked$variogram$observed)
  band <- ggplot2::layer_data(p, 1)
  expect_equal(band$ymax, checked$variogram$upper)
  expect_identical(p$labels$x, "Distance (km)")
  expect_identical(p$labels$y, "Semivariance of the residuals")

  # What is outside has another colour.
  altered <- checked
  altered$variogram$outside[1] <- TRUE
  colours <- ggplot2::layer_data(plot(altered, type = "variogram"), 3)$colour
  expect_identical(colours[1], validation_palette()$contrast)
  expect_true(all(colours[-1] == validation_palette()$line))
  green <- plot(checked, type = "variogram", colour = "darkblue")
  expect_true(all(ggplot2::layer_data(green, 3)$colour == "darkblue"))

  p <- plot_predictive_check(checked, type = "distribution", model = "a")
  expect_s3_class(p$facet, "FacetNull")
  expect_identical(p$labels$x, "y")
  expect_match(p$labels$y, "observations")
  line <- ggplot2::layer_data(p, 3)
  expect_equal(range(line$y), range(subset(checked$distribution,
                                           model == "a")$observed))
  prevalence <- plot(assess_fit(binary_twice, nsim = 19, seed = 1),
                     type = "distribution")
  expect_identical(prevalence$labels$x, "Prevalence")
  expect_match(prevalence$labels$y, "locations")
  rate <- plot(assess_fit(poisson_twice, nsim = 19, seed = 1),
               type = "distribution")
  expect_identical(rate$labels$x, "Rate")
})

test_that("plot_predictive_check() checks its arguments", {
  expect_input_error(plot_predictive_check(), "x")
  expect_input_error(plot_predictive_check(gaussian_model), "x", "assess_fit")
  expect_input_error(plot_predictive_check(checked, type = "pit"), "type")
  expect_input_error(plot_predictive_check(checked, model = "c"), "model")
  expect_input_error(plot_predictive_check(checked, colour = "nocolour"), "colour")
  expect_error(plot_predictive_check(checked, bins = 3),
               class = "riskmap_error_input")
})

test_that("the bars of a histogram hold equally spaced values alike", {
  # Proportions of 90 locations, up to 12 of them.
  values <- rep(0:12, 1:13) / 90
  breaks <- histogram_breaks(values)
  expect_equal(diff(breaks), rep(1 / 90, 13))
  bars <- statistic_bars(
    data.frame(model = "a", statistic = "s", value = values),
    data.frame(model = "a", statistic = "s", observed = 3 / 90)
  )
  expect_identical(bars$count, 1:13)
  expect_identical(bars$above, 0:12 >= 3)
  # Up to 80 of them, four in a bar, one of which starts at the observed one.
  values <- (0:80) / 90
  breaks <- histogram_breaks(values)
  expect_equal(diff(breaks), c(rep(4, 20), 1) / 90)
  breaks <- histogram_breaks(values, 6 / 90)
  expect_true(any(abs(breaks - 5.5 / 90) < 1e-6))
  # The bars start and end with the values.
  expect_equal(range(breaks), c(-0.5, 80.5) / 90)
  bars <- statistic_bars(
    data.frame(model = "a", statistic = "s", value = values),
    data.frame(model = "a", statistic = "s", observed = 6 / 90)
  )
  expect_identical(sum(bars$count[bars$above]), 75L)

  # Values that aren't equally spaced.
  set.seed(1)
  values <- rnorm(100)
  expect_length(histogram_breaks(values), 26L)
  expect_identical(range(histogram_breaks(values, 0.3)), range(values))
  for (observed in c(-10, -1, 0.3, 1, 10)) {
    bars <- statistic_bars(
      data.frame(model = "a", statistic = "s", value = values),
      data.frame(model = "a", statistic = "s", observed = observed)
    )
    expect_identical(sum(bars$count), 100L)
    expect_identical(sum(bars$count[bars$above]), sum(values >= observed))
  }
  # Only one value, and an observed one that is missing.
  expect_equal(histogram_breaks(rep(2, 5)), c(1.9, 2.1))
  expect_equal(histogram_breaks(0), c(-0.5, 0.5))
  bars <- statistic_bars(
    data.frame(model = "a", statistic = "s", value = c(2, 2)),
    data.frame(model = "a", statistic = "s", observed = 2)
  )
  expect_true(bars$above)
  bars <- statistic_bars(
    data.frame(model = "a", statistic = "s", value = c(2, 3)),
    data.frame(model = character(), statistic = character(),
               observed = numeric())
  )
  expect_false(any(bars$above))
})
