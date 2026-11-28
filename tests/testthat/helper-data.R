# The shared fixtures of the tests (#206). Every dataset is simulated from the
# model that the tests fit to it, at the same 20 locations in a 10 km square
# and with the same spatial field, so that each fit estimates effects its data
# contain.
#
# - `gaussian_data`, `binomial_data` and `poisson_data` come from
#   `y ~ cov + gp()`, with a measurement error for the Gaussian data and the
#   exposure `denominator` for the counts.
# - `gaussian_re_data`, `binomial_re_data` and `poisson_re_data` come from
#   `y ~ cov + gp() + re(i)`, and `gaussian_model`, `binomial_model` and
#   `poisson_model` are fitted to them.
# - `gaussian_two_re_data` and `binomial_two_re_data` come from
#   `y ~ cov + gp() + re(i, j)`, with four levels of `j`.
# - `gaussian_offset_data` comes from the model of `gaussian_offset_model` and
#   `gaussian_intercept_data` from `y ~ gp()`, that of
#   `gaussian_intercept_model`.
# - `binomial_nugget_data` comes from `y ~ cov + gp(nugget = TRUE)`, with five
#   times the denominators, without which the nugget is not identified.
# - `negative_binomial_data` comes from `y ~ cov + gp()` with a negative
#   binomial family of size 2 and the exposure `denominator`.
#
# Each dataset draws its response with its own seed, chosen so that these fits
# are clean, with a scale longer than the median distance between neighbouring
# locations (1.1 km), variances away from zero and every check of
# `check_computation()` OK. `simulate_fixture()` simulates data for another
# model at the same locations.
set.seed(1)
fixture_n <- 20L
fixture_groups <- fixture_n / 2L
fixture_columns <- withr::with_seed(1, {
  data.frame(x = runif(fixture_n, 0, 10000), z = runif(fixture_n, 0, 10000),
             cov = rnorm(fixture_n),
             denominator = sample(5:20, fixture_n, replace = TRUE),
             offset = rnorm(fixture_n),
             i = rep(seq_len(fixture_groups), each = 2),
             j = rep(1:2, each = fixture_groups))
})
# A field of variance 1 and scale 5.5 km at the locations.
fixture_field <- withr::with_seed(2, {
  correlation <- matern_correlation(dist(fixture_columns[, c("x", "z")]),
                                    phi = 5500, kappa = 1,
                                    return_sym_matrix = TRUE)
  as.numeric(t(chol(correlation)) %*% rnorm(fixture_n))
})

# The locations with the response `y` that `response(columns)` draws with
# `seed`. Columns given in `...` replace those of `fixture_columns`.
simulate_fixture <- function(seed, response, ...) {
  columns <- utils::modifyList(fixture_columns, list(...))
  columns$y <- withr::with_seed(seed, response(columns))
  st_as_sf(columns, coords = c("x", "z"), crs = 32637)
}
# Independent effects of the levels, one for each row.
fixture_effect <- function(levels, variance = 1) {
  rnorm(max(levels), sd = sqrt(variance))[levels]
}

gaussian_data <- simulate_fixture(20, function(columns) {
  3 + 0.5 * columns$cov + fixture_field + rnorm(fixture_n, sd = sqrt(0.3))
})
gaussian_re_data <- simulate_fixture(10, function(columns) {
  3 + 0.5 * columns$cov + fixture_field + fixture_effect(columns$i) +
    rnorm(fixture_n, sd = sqrt(0.3))
})
gaussian_two_re_data <- simulate_fixture(60, function(columns) {
  3 + 0.5 * columns$cov + fixture_field + fixture_effect(columns$i) +
    fixture_effect(columns$j) + rnorm(fixture_n, sd = sqrt(0.3))
}, j = rep(1:4, length.out = fixture_n))
gaussian_offset_data <- simulate_fixture(14, function(columns) {
  3 + 0.5 * columns$cov + columns$offset + fixture_field +
    rnorm(fixture_n, sd = sqrt(0.3))
})
gaussian_intercept_data <- simulate_fixture(38, function(columns) {
  3 + fixture_field + rnorm(fixture_n, sd = sqrt(0.3))
})
latlon_data <- st_transform(gaussian_data, crs = 4326)

binomial_data <- simulate_fixture(33, function(columns) {
  rbinom(fixture_n, columns$denominator,
         plogis(0.2 + 0.3 * columns$cov + fixture_field))
})
binomial_re_data <- simulate_fixture(18, function(columns) {
  rbinom(fixture_n, columns$denominator,
         plogis(0.2 + 0.3 * columns$cov + fixture_field +
                  fixture_effect(columns$i)))
})
binomial_two_re_data <- simulate_fixture(31, function(columns) {
  rbinom(fixture_n, columns$denominator,
         plogis(0.2 + 0.3 * columns$cov + fixture_field +
                  fixture_effect(columns$i) + fixture_effect(columns$j)))
}, j = rep(1:4, length.out = fixture_n))
binomial_nugget_data <- simulate_fixture(20, function(columns) {
  rbinom(fixture_n, columns$denominator,
         plogis(0.2 + 0.3 * columns$cov + fixture_field +
                  rnorm(fixture_n, sd = sqrt(0.5))))
}, denominator = 5L * fixture_columns$denominator)

poisson_data <- simulate_fixture(35, function(columns) {
  rpois(fixture_n, columns$denominator *
          exp(-1 + 0.2 * columns$cov + fixture_field))
})
poisson_re_data <- simulate_fixture(3, function(columns) {
  rpois(fixture_n, columns$denominator *
          exp(-1 + 0.2 * columns$cov + fixture_field +
                fixture_effect(columns$i)))
})
negative_binomial_data <- simulate_fixture(27, function(columns) {
  rnbinom(fixture_n, size = 2, mu = columns$denominator *
            exp(-1 + 0.2 * columns$cov + fixture_field))
})

# reduce iterations to speed up fits
control_mcmc <- set_control_mcmc(n_draws = 80,
                                 burnin = 200)
fixture_seed <- 12345
fixture_control_mcml <- set_control_mcml(
  initial_draws = 500, validation_draws = 500, max_draws = 500,
  max_stages = 1
)

gaussian_model <- glgpm(y ~ cov + gp() + re(i),
                        data = gaussian_re_data,
                        family = "gaussian",
                        messages = FALSE)

# gaussian_model without the covariate, a second model of the same data.
gaussian_no_cov_model <- glgpm(y ~ gp() + re(i),
                               data = gaussian_re_data,
                               family = "gaussian",
                               messages = FALSE)

gaussian_offset_model <- glgpm(y ~ cov + gp(nugget = TRUE) + offset(offset),
                               data = gaussian_offset_data,
                               family = "gaussian",
                               fixed_pars = list(sigma2_me = 0),
                               messages = FALSE)

gaussian_intercept_model <- glgpm(y ~ gp(),
                                  data = gaussian_intercept_data,
                                  family = "gaussian",
                                  messages = FALSE)

binomial_model <- glgpm(
  y ~ cov + gp() + re(i), data = binomial_re_data, family = "binomial",
  denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
  control_mcml = fixture_control_mcml, messages = FALSE
)

poisson_model <- glgpm(
  y ~ cov + gp() + re(i), data = poisson_re_data, family = "poisson",
  denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
  control_mcml = fixture_control_mcml, return_samples = TRUE,
  messages = FALSE
)

hull <- create_convex_hull(gaussian_data)
grid <- create_grid(hull, 3)
squares <- st_make_grid(hull, n = c(2, 2))
areal <- st_sf(geometry = squares)

prediction_data <- grid
prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
prediction_data$offset <- seq_len(nrow(prediction_data)) / 100
prediction_data$i <- rep(gaussian_model$re[[1]],
                         length.out = nrow(prediction_data))
prediction_data$population <- seq_len(nrow(prediction_data))

# A reference of independent uniform PIT values, for tests of the scores of
# fits that are mocks without a working model (`working_model()`).
independent_reference <- function(fit, test_sets, predictive, n_curves, ...) {
  rows <- unlist(lapply(predictive, `[[`, "rows"))
  pit <- withr::with_seed(1, matrix(runif(n_curves * length(rows)), n_curves))
  curves <- pit_curves(pit, pit, area_levels())
  list(curves = curves[, seq(1L, ncol(curves), by = 10L), drop = FALSE],
       areas = pit_area(curves, area_levels()))
}

cross_validation <- assess_prediction(
  list(model_a = gaussian_model, model_b = gaussian_no_cov_model),
  list(c(1, 3, 5, 7), c(2, 4, 6, 8, 9, 10)),
  messages = FALSE
)

# A cross-validation in which `models` are binomial models of binary outcomes.
# Their scores become those of one trial, a PIT interval from 0 for y = 0 or to
# 1 for y = 1 with the Brier score as the CRPS, with y the parity of the row.
# The CRPS is computed from F(0) = 1 - p as `count_crps()` computes it. The
# calibration is left as it is.
binary_scores <- function(x, models) {
  for (model in models) {
    rows <- x$scores$model == model
    y <- x$scores$row[rows] %% 2
    below <- 1 - seq(0.2, 0.6, length.out = sum(rows))
    x$scores$pit_lower[rows] <- ifelse(y == 0, 0, below)
    x$scores$pit_upper[rows] <- ifelse(y == 0, below, 1)
    x$scores$CRPS[rows] <- (below - (y == 0))^2
    x$family[model] <- "binomial"
  }
  x
}

surface <- grid
surface$cov <- rnorm(nrow(grid))
surface$i <- rep(1, nrow(grid))

simulation <- simulate_glgpm(gaussian_model,
                             what = c("data", "surface"),
                             prediction_grid = surface,
                             sample_locations = gaussian_re_data
                             )
