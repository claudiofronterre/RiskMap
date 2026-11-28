# A Gaussian fit must not depend on the units of the response. Multiplied by a
# constant, the response has coefficients multiplied by it, variances by its
# square, the same scale of correlation, and a log-likelihood that differs by
# the number of observations times the constant's logarithm.

# Locations in kilometres, as data with a projected CRS in metres.
units_data <- function(locations, columns) {
  columns$easting <- 1000 * locations[, 1] + 500000
  columns$northing <- 1000 * locations[, 2]
  st_as_sf(columns, coords = c("easting", "northing"), crs = 32633)
}

# Records at 30 locations with a covariate, an offset, a cluster effect, a
# nugget and measurement error.
simulate_units <- function(records, seed = 11) {
  set.seed(seed)
  locations <- cbind(runif(30, 0, 100), runif(30, 0, 100))
  location <- rep(seq_len(30), each = records)
  n <- length(location)
  spatial <- crossprod(chol(exp(-as.matrix(dist(locations)) / 20)), rnorm(30))
  columns <- data.frame(x = rnorm(n), known = rnorm(n, 10, 3),
                        cluster = sample(rep_len(1:6, n)))
  columns$y <- 2 + 0.7 * columns$x + spatial[location] +
    rnorm(30, 0, 0.5)[location] + rnorm(6, 0, 0.5)[columns$cluster] +
    rnorm(n, 0, 0.3)
  list(locations = locations[location, ], location = location,
       columns = columns)
}

# The records with the offset `known` added to the response, for a model with
# that offset.
with_offset <- function(simulated) {
  simulated$columns$y <- simulated$columns$y + simulated$columns$known
  simulated
}

in_units <- function(simulated, scale) {
  columns <- simulated$columns
  columns[c("y", "known")] <- scale * columns[c("y", "known")]
  units_data(simulated$locations, columns)
}

# The numbers of a fit in units in which the response was multiplied by scale.
from_units <- function(fit, scale) {
  factor_of <- function(names) {
    ifelse(grepl("^beta", names), scale,
           ifelse(grepl("^sigma2|^tau2", names), scale^2, 1))
  }
  estimates <- unlist(coef(fit))
  errors <- sqrt(diag(vcov(fit)))
  list(
    estimates = estimates / factor_of(names(estimates)),
    errors = errors / factor_of(names(errors)),
    log_likelihood = as.numeric(logLik(fit)) + nobs(fit) * log(scale),
    fitted = fitted(fit) / scale
  )
}

# The log-likelihood, up to the constant and with the coefficients at their
# best, written out without RiskMap.
dense_log_likelihood <- function(simulated, terms, scale) {
  y <- scale * simulated$columns$y
  design <- cbind(1, simulated$columns$x)
  distances <- as.matrix(dist(simulated$locations))
  same_location <- outer(simulated$location, simulated$location, "==")
  same_cluster <- outer(simulated$columns$cluster, simulated$columns$cluster,
                        "==")
  function(par) {
    value <- exp(stats::setNames(par, terms))
    covariance <- value[["sigma2"]] * exp(-distances / value[["phi"]]) +
      diag(value[["sigma2_me"]], length(y))
    if ("tau2" %in% terms) {
      covariance <- covariance + value[["tau2"]] * same_location
    }
    if ("sigma2_re" %in% terms) {
      covariance <- covariance + value[["sigma2_re"]] * same_cluster
    }
    root <- tryCatch(chol(covariance), error = function(error) NULL)
    if (is.null(root)) return(-1e10)
    whitened <- lm.fit(backsolve(root, design, transpose = TRUE),
                       backsolve(root, y, transpose = TRUE))$residuals
    -sum(log(diag(root))) - sum(whitened^2) / 2
  }
}

# Its maximum, from three scales of correlation and with the variance of the
# response divided equally between the terms.
dense_maximum <- function(simulated, terms, scale) {
  objective <- dense_log_likelihood(simulated, terms, scale)
  start <- stats::setNames(rep(scale^2 * stats::var(simulated$columns$y) /
                                 (length(terms) - 1), length(terms)), terms)
  max(vapply(c(5, 20, 60), function(phi) {
    start[["phi"]] <- phi
    -optim(log(start), function(par) -objective(par), method = "BFGS",
           control = list(maxit = 500, reltol = 1e-12))$value
  }, numeric(1)))
}

test_that("a Gaussian fit does not depend on the units of the response", {
  skip_on_cran()
  single <- simulate_units(1)
  repeated <- simulate_units(2)
  # A seed at which the fits of six records estimate the nugget away from zero.
  many <- simulate_units(6, seed = 8)
  models <- list(
    plain = list(y ~ gp(), single),
    covariate = list(y ~ x + gp(), single),
    no_intercept = list(y ~ 0 + x + gp(kappa = 1.5), single),
    offset = list(y ~ x + offset(known) + gp(), with_offset(single)),
    random_effect = list(y ~ x + gp() + re(cluster), single),
    nugget = list(y ~ x + gp(nugget = TRUE), repeated),
    nugget_re = list(y ~ x + gp(nugget = TRUE) + re(cluster), repeated),
    # Six records at a location take the Woodbury path.
    woodbury = list(y ~ x + gp(nugget = TRUE) + re(cluster), many),
    no_measurement_error = list(y ~ x + gp(nugget = TRUE), single,
                                fixed = function(k) list(sigma2_me = 0)),
    held_nugget = list(y ~ x + gp(nugget = TRUE), repeated,
                       fixed = function(k) list(tau2 = 0.25 * k^2)),
    held_error = list(y ~ x + gp(nugget = TRUE), many,
                      fixed = function(k) list(sigma2_me = 0.3 * k^2)),
    held_others = list(y ~ x + gp() + re(cluster), single,
                       fixed = function(k) {
                         list(sigma2 = k^2, sigma2_re = 0.25 * k^2)
                       }),
    held_coefficients = list(y ~ x + gp(), single, fixed = function(k) {
      list(beta = c("(Intercept)" = 1.5 * k))
    }),
    held_slope = list(y ~ x + gp(), single, fixed = function(k) {
      list(beta = c(x = 0.5 * k), phi = 15)
    }),
    all_held = list(y ~ x + gp(), single, fixed = function(k) {
      list(beta = k * c(2, 0.7), sigma2 = k^2, phi = 20, sigma2_me = 0.3 * k^2)
    }),
    some_started = list(y ~ x + gp(nugget = TRUE), repeated,
                        start = function(k) {
                          list(beta = c(x = 2 * k), sigma2_me = 0.1 * k^2)
                        }),
    all_started = list(y ~ x + gp(), single, start = function(k) {
      list(beta = c(k, k), sigma2 = 3 * k^2, phi = 30, sigma2_me = k^2)
    })
  )
  fit_in <- function(model, scale) {
    glgpm(model[[1]], data = in_units(model[[2]], scale), family = "gaussian",
          fixed_pars = if (!is.null(model$fixed)) model$fixed(scale),
          start_pars = if (!is.null(model$start)) model$start(scale),
          messages = FALSE)
  }
  for (name in names(models)) {
    reference <- from_units(fit_in(models[[name]], 1), 1)
    for (scale in c(1e-6, 1e3, 1e6)) {
      changed <- from_units(fit_in(models[[name]], scale), scale)
      label <- paste(name, "times", scale)
      expect_equal(changed$log_likelihood, reference$log_likelihood,
                   tolerance = 1e-8, label = label)
      expect_equal(changed$estimates, reference$estimates, tolerance = 1e-6,
                   label = label)
      expect_equal(changed$errors, reference$errors, tolerance = 1e-5,
                   label = label)
      expect_equal(changed$fitted, reference$fitted, tolerance = 1e-6,
                   label = label)
    }
  }
})

test_that("a Gaussian fit does not depend on the mean of the response", {
  simulated <- simulate_units(1)
  fit <- glgpm(y ~ x + gp(), data = in_units(simulated, 1),
               family = "gaussian", messages = FALSE)
  simulated$columns$y <- simulated$columns$y + 1e6
  moved <- glgpm(y ~ x + gp(), data = in_units(simulated, 1),
                 family = "gaussian", messages = FALSE)
  # An offset that takes the mean away again.
  simulated$columns$known <- 1e6
  offset <- glgpm(y ~ x + offset(known) + gp(), data = in_units(simulated, 1),
                  family = "gaussian", messages = FALSE)
  for (other in list(moved, offset)) {
    expect_equal(as.numeric(logLik(other)), as.numeric(logLik(fit)),
                 tolerance = 1e-8)
    expect_equal(unlist(coef(other))[-1], unlist(coef(fit))[-1],
                 tolerance = 1e-6)
  }
  expect_equal(coef(moved)$beta[[1]] - 1e6, coef(fit)$beta[[1]],
               tolerance = 1e-5)
  expect_equal(coef(offset)$beta[[1]], coef(fit)$beta[[1]], tolerance = 1e-5)
})

test_that("a Gaussian fit in any units reaches the maximum found without RiskMap", {
  skip_on_cran()
  single <- simulate_units(1)
  repeated <- simulate_units(2)
  many <- simulate_units(4)
  cases <- list(
    list(y ~ x + gp(), single, c("sigma2", "phi", "sigma2_me")),
    list(y ~ x + gp(nugget = TRUE) + re(cluster), repeated,
         c("sigma2", "phi", "tau2", "sigma2_me", "sigma2_re")),
    list(y ~ x + gp(nugget = TRUE), many,
         c("sigma2", "phi", "tau2", "sigma2_me"))
  )
  for (case in cases) {
    highest <- dense_maximum(case[[2]], case[[3]], 1)
    for (scale in c(1e-3, 1e3)) {
      fit <- glgpm(case[[1]], data = in_units(case[[2]], scale),
                   family = "gaussian", messages = FALSE)
      expect_equal(as.numeric(logLik(fit)) + nobs(fit) * log(scale), highest,
                   tolerance = 1e-6)
      # RiskMap's log-likelihood is the one written out above.
      objective <- dense_log_likelihood(case[[2]], case[[3]], scale)
      expect_equal(objective(log(unlist(coef(fit)[case[[3]]]))),
                   as.numeric(logLik(fit)), tolerance = 1e-8)
    }
  }
})

test_that("a response with a large variance does not stop far from the maximum", {
  skip_on_cran()
  # With every variance starting at 1, 55% of these data sets multiplied by
  # 100, and 90% multiplied by 1000, stopped more than 0.1 below the maximum,
  # by up to 17.8, most of them with no spatial variance or at a scale far
  # from the distances. Some of the data sets have their maximum with no
  # measurement error or at a scale below the spacing of the locations, which
  # the search must reach too.
  for (seed in seq_len(20)) {
    set.seed(seed)
    locations <- cbind(runif(50, 0, 100), runif(50, 0, 100))
    spatial <- crossprod(chol(exp(-as.matrix(dist(locations)) / 20)),
                         rnorm(50))
    columns <- data.frame(x = 0, y = 1 + as.numeric(spatial) +
                            rnorm(50, 0, 0.2), known = 0)
    simulated <- list(locations = locations, location = seq_len(50),
                      columns = columns)
    fits <- lapply(c(1, 100, 1000), function(scale) {
      from_units(glgpm(y ~ gp(), data = in_units(simulated, scale),
                       family = "gaussian", messages = FALSE), scale)
    })
    reached <- vapply(fits, `[[`, numeric(1), "log_likelihood")
    expect_equal(reached, rep(reached[1], 3), tolerance = 1e-8)
    # The search from the simulated values, without RiskMap. The design there
    # has a column of zeros, which lm.fit() drops.
    terms <- c("sigma2", "phi", "sigma2_me")
    objective <- dense_log_likelihood(simulated, terms, 1)
    best <- optim(log(c(1, 20, 0.04)), function(par) -objective(par),
                  method = "BFGS", control = list(maxit = 500))
    expect_gt(reached[1], -best$value - 1e-4)
  }
})

test_that("the variances of a Gaussian model start at its residual variance", {
  simulated <- simulate_units(2)
  data <- in_units(simulated, 50)
  start_of <- function(formula, ...) {
    captured <- NULL
    try(with_mocked_bindings(
      glgpm(formula, data = data, family = "gaussian", messages = FALSE, ...),
      fit_glgpm_design = function(start_pars, ...) {
        captured <<- start_pars
        stop("captured the start")
      }
    ), silent = TRUE)
    captured
  }
  regression <- lm(y ~ x, data = data)
  variance <- mean(residuals(regression)^2)
  formula <- y ~ x + gp(nugget = TRUE) + re(cluster)
  start <- start_of(formula)
  expect_equal(start$beta, unname(coef(regression)), tolerance = 1e-10)
  expect_equal(unname(unlist(start[c("sigma2", "tau2", "sigma2_me",
                                     "sigma2_re")])),
               rep(variance, 4), tolerance = 1e-10)
  # A value given or held replaces its own and no other.
  start <- start_of(formula, start_pars = list(sigma2_me = 3))
  expect_identical(start$sigma2_me, 3)
  expect_equal(start$sigma2, variance, tolerance = 1e-10)
  start <- start_of(formula, start_pars = list(sigma2_re = c(cluster = 2)))
  expect_identical(start$sigma2_re, 2)
  expect_equal(start$sigma2, variance, tolerance = 1e-10)
  # The regression takes the held coefficients as known, and the offset.
  start <- start_of(y ~ x + offset(known) + gp(),
                    fixed_pars = list(beta = c(x = 10)))
  given_slope <- lm(y ~ 1 + offset(10 * x + known), data = data)
  expect_equal(start$beta, c(coef(given_slope)[[1]], 10), tolerance = 1e-10)
  expect_equal(start$sigma2, mean(residuals(given_slope)^2),
               tolerance = 1e-10)

  parameters <- glgpm_parameters(formula, data = data, family = "gaussian")
  expect_identical(unique(parameters$start[!parameters$name %in%
                                             c("beta", "phi")]),
                   "residual_variance")
  held <- glgpm_parameters(formula, data = data, family = "gaussian",
                           fixed_pars = list(tau2 = 1, sigma2_re = 1))
  expect_identical(held$start, parameters$start)
  # The variances of the other families start at 1.
  counts <- glgpm_parameters(y ~ x + gp(nugget = TRUE) + re(i),
                             data = transform(poisson_data, x = cov),
                             family = "poisson")
  expect_identical(unique(counts$start[!counts$name %in% c("beta", "phi")]),
                   "one")
})

test_that("the starting regression is fitted on the standardised design", {
  # The columns of a cubic in the year are collinear to lm.fit() as they are,
  # and independent once they are centred and scaled.
  simulated <- simulate_units(1)
  simulated$columns$year <- 2000 + 20 * seq(0, 1, length.out = 30)
  data <- in_units(simulated, 1)
  formula <- y ~ x + year + I(year^2) + I(year^3) + gp()
  expect_lt(qr(model.matrix(~ year + I(year^2) + I(year^3), data))$rank, 4)
  fit <- glgpm(formula, data = data, family = "gaussian", messages = FALSE)
  started <- glgpm(formula, data = data, family = "gaussian",
                   messages = FALSE,
                   start_pars = list(beta = c(mean(data$y), 0, 0, 0, 0)))
  expect_equal(as.numeric(logLik(started)), as.numeric(logLik(fit)),
               tolerance = 1e-6)
  # The same model in the year since 2010.
  data$year <- data$year - 2010
  centred <- glgpm(formula, data = data, family = "gaussian",
                   messages = FALSE)
  expect_equal(as.numeric(logLik(centred)), as.numeric(logLik(fit)),
               tolerance = 1e-6)
})

test_that("a response without residuals is fitted as it is", {
  parameters <- model_parameter_table(
    "(Intercept)", fix_tau2 = 0, exact = TRUE, fix_var_me = NULL,
    estimated_family = character(), family_parameters = NULL,
    re_names = NULL, distance_units = "km"
  )
  design <- matrix(1, 7, 1, dimnames = list(NULL, "(Intercept)"))
  regression <- gaussian_regression(rep(1 / 3, 7), 0, design, parameters)
  expect_identical(regression$variance, 0)
  # Nor does what is left by subtracting an offset that differs from the
  # response by a constant.
  offset <- 1e6 * c(1.1, 0.7, 1.3, 0.2, 0.9, 1.7, 0.4) / 3
  expect_gt(stats::sd(offset + 0.1 - offset), 0)
  expect_identical(
    gaussian_regression(offset + 0.1, offset, design, parameters)$variance, 0
  )
  expect_identical(default_start(parameters, "sigma2", NULL, 0), 1)
  scaling <- scale_gaussian_response(regression)
  expect_equal(scaling$response, rep(0, 7))
  expect_identical(scaling$variance(0.4), 0.4)
  expect_null(scaling$variance(NULL))
  expect_true(scaling$variance(TRUE))
})

test_that("the profiles of a Gaussian fit do not depend on the units", {
  skip_on_cran()
  simulated <- simulate_units(2)
  intervals <- lapply(c(1, 1e-6, 1e6), function(scale) {
    fit <- glgpm(y ~ x + gp(nugget = TRUE) + re(cluster),
                 data = in_units(simulated, scale),
                 family = "gaussian", messages = FALSE)
    profiled <- profile(fit, parameters = c("x", "sigma2", "sigma2_me"),
                        n_points = 11, messages = FALSE)$intervals
    as.matrix(profiled[c("profile_lower", "profile_upper")]) /
      c(scale, scale^2, scale^2)
  })
  expect_true(all(is.finite(intervals[[1]])))
  for (changed in intervals[-1]) {
    expect_equal(changed, intervals[[1]], tolerance = 1e-3)
  }
})
