# A minimal Gaussian fit, whose test sets are predicted by a mock.
make_assess_prediction_fit <- function(data, covariate) {
  fit <- list(
    formula = stats::as.formula(paste("y ~", covariate)),
    data = data,
    input_crs = sf::st_crs(data),
    family = "gaussian",
    estimate = list(beta = c(0, 0), sigma2 = 0, phi = 0),
    D = matrix(1, nrow = nrow(data), ncol = 2),
    re = list(),
    cov_offset = NULL,
    nugget = "none",
    fixed_pars = list(sigma2_me = 0),
    sst = FALSE,
    units_m = rep(1, nrow(data)),
    y = data$y,
    ID_coords = seq_len(nrow(data)),
    coords = sf::st_coordinates(data),
    distance_units = "m",
    call = list(denominator = quote(units_m)),
    model_id = covariate
  )
  class(fit) <- "RiskMap"
  fit
}

# Data at a few points of a line, with any other columns.
line_data <- function(...) {
  columns <- data.frame(...)
  points <- lapply(seq_len(nrow(columns)), function(i) st_point(c(1000 * i, 0)))
  st_sf(columns, geometry = st_sfc(points, crs = 32637))
}

# Draws of the linear predictor for any test set, those of `draws` at the rows
# with these values of `x1`.
mock_draws <- function(draws) {
  function(fit, newdata, control_mcmc, seed) draws[newdata$x1, , drop = FALSE]
}

mock_constant <- function(value) {
  function(fit, newdata, control_mcmc, seed) {
    matrix(value, nrow(newdata), 4000)
  }
}

assess <- function(...) assess_prediction(..., messages = FALSE)

test_that("assess_prediction() checks its arguments", {
  fit <- gaussian_intercept_model
  expect_input_error(assess_prediction(), "models")
  expect_input_error(assess(fit), "test_sets")
  expect_input_error(assess("fit", list(1:2)), "models", "a list of fits")
  expect_input_error(assess(list(), list(1:2)), "models")
  expect_input_error(assess(list(fit, 1), list(1:2)), "models")
  expect_input_error(assess(list(a = fit, a = fit), list(1:2)), "models",
                     "\"a\" is repeated")
  expect_input_error(assess(fit, "cluster"), "test_sets")
  expect_input_error(assess(fit, list(1:2), refit = "no"), "refit")
  expect_input_error(assess(fit, list(1:2), control_mcmc = "no"),
                     "control_mcmc")
  expect_input_error(assess(fit, list(1:2), seed = 1.5), "seed")
  expect_input_error(assess_prediction(fit, list(1:2), messages = "no"),
                     "messages")
})

test_that("the former arguments of assess_prediction() say what replaces them", {
  fit <- gaussian_intercept_model
  former <- list(
    method = list("cluster", "create_test_sets()"),
    fold = list(3, "create_test_sets()"),
    min_dist = list(3, "create_test_sets()"),
    size = list(3, "create_test_sets()"),
    iter = list(3, "create_test_sets()"),
    plot_fold = list(TRUE, "plot_test_sets()"),
    user_split = list(list(1:2), "in `test_sets`"),
    metrics = list("CRPS", "All metrics are always computed"),
    keep_par_fixed = list(FALSE, "`refit`"),
    which_metric = list("CRPS", "All metrics are always computed"),
    n_size = list(3, "create_test_sets()"),
    control_sim = list(1, "`control_mcmc`")
  )
  for (name in names(former)) {
    call <- c(list(fit), stats::setNames(former[[name]][1], name))
    expect_input_error(do.call(assess_prediction, call), name,
                       former[[name]][[2]])
  }
  # The former name of `models`, without any of the present arguments.
  expect_input_error(assess_prediction(object = list(fit), method = "cluster"),
                     "object", "It is now `models`")
  # Arguments with the same advice give it once.
  message <- tryCatch(assess(fit, method = "cluster", fold = 3),
                      error = conditionMessage)
  expect_length(gregexpr("Make the test sets", message)[[1]], 1L)
  # Also when test_sets is given, and a misspelt argument too.
  expect_input_error(assess(fit, list(1:2), keep_par_fixed = TRUE),
                     "keep_par_fixed", "`refit`")
  expect_input_error(assess(fit, list(1:2), refits = TRUE), "refits")
})

test_that("a fit or an unnamed list of fits has names in the results", {
  local_mocked_bindings(predict_test_set = mock_constant(0))
  fit <- gaussian_intercept_model
  expect_identical(unique(assess(fit, list(1:2))$scores$model), "Model 1")
  expect_identical(unique(assess(list(fit, fit), list(1:2))$scores$model),
                   c("Model 1", "Model 2"))
  expect_identical(unique(assess(list(a = fit, fit), list(1:2))$scores$model),
                   c("a", "Model 2"))
})

test_that("models of different observations can't be compared", {
  fit <- function(data, ...) glgpm(y ~ gp(), data = data, family = "gaussian",
                                   messages = FALSE, ...)
  first <- gaussian_intercept_model
  refused <- function(other, ...) {
    expect_input_error(assess(list(a = first, b = other), list(1:2)), "models",
                       ...)
  }
  refused(fit(gaussian_data[1:9, ]), "is fitted to 9 of them, and \"a\" to 10")
  refused(fit(gaussian_data[10:1, ]), "same rows in the same order")
  moved <- gaussian_data
  st_geometry(moved)[1] <- st_geometry(moved)[1] + c(1, 0)
  refused(fit(st_set_crs(moved, st_crs(gaussian_data))),
          "same rows in the same order")
  # Defect 17, a response on another scale.
  other <- gaussian_data
  other$y <- exp(other$y)
  refused(fit(other), "The responses of \"b\" and \"a\" differ")
  # Rounding in the last digits isn't a difference.
  rounded <- gaussian_data
  rounded$y <- (rounded$y / 3) * 3
  local_mocked_bindings(predict_test_set = mock_constant(0))
  expect_s3_class(assess(list(first, fit(rounded)), list(1:2)),
                  "RiskMap_cross_validation")
})

test_that("models of different denominators can't be compared", {
  data <- line_data(y = c(1, 2, 3), n = c(5, 5, 5))
  binomial <- function(units_m) {
    fit <- make_assess_prediction_fit(data, "1")
    fit$family <- "binomial"
    fit$units_m <- units_m
    fit$link_function <- list(inv = plogis, name = "logit")
    fit
  }
  expect_input_error(
    assess(list(a = binomial(c(5, 5, 5)), b = binomial(c(5, 5, 6))), list(1:2)),
    "models", "same denominators"
  )
  # A model without denominators has ones.
  poisson <- binomial(c(1, 1, 1))
  poisson$family <- "poisson"
  poisson$link_function <- list(inv = exp, name = "log")
  without <- poisson
  without$units_m <- NULL
  local_mocked_bindings(predict_test_set = mock_constant(0))
  out <- assess(list(a = poisson, b = without), list(1:2))
  expect_identical(out$scores$CRPS[1:2], out$scores$CRPS[3:4])
})

test_that("each model predicts from its own data", {
  seen <- list()
  local_mocked_bindings(
    predict_test_set = function(fit, newdata, ...) {
      seen[[fit$model_id]] <<- names(newdata)
      matrix(1, nrow(newdata), 2)
    }
  )
  out <- assess(
    list(model_x1 = make_assess_prediction_fit(
           line_data(y = 1:3, x1 = c(10, 20, 30)), "x1"),
         model_x2 = make_assess_prediction_fit(
           line_data(y = 1:3, x2 = c(100, 200, 300)), "x2")),
    list(2:3)
  )
  expect_setequal(seen$x1, c("y", "x1", "geometry"))
  expect_setequal(seen$x2, c("y", "x2", "geometry"))
})

test_that("the result has a row for each model and test observation", {
  local_mocked_bindings(predict_test_set = mock_constant(0))
  sets <- list(c(3, 1), list(train = 5:10, test = c(2, 3)))
  out <- assess(list(a = gaussian_intercept_model, b = gaussian_offset_model),
                sets)
  expect_s3_class(out, "RiskMap_cross_validation")
  expect_named(out, c("scores", "calibration", "test_sets", "refit", "family",
                      "failures"), ignore.order = TRUE)
  expect_named(out$scores, c("model", "test_set", "row", "CRPS", "SCRPS",
                             "pit_lower", "pit_upper"))
  expect_identical(out$scores$model, rep(c("a", "b"), each = 4))
  expect_identical(out$scores$test_set, rep(c(1L, 1L, 2L, 2L), 2))
  expect_identical(out$scores$row, rep(c(3L, 1L, 2L, 3L), 2))
  expect_identical(row.names(out$scores), as.character(1:8))
  expect_s3_class(out$test_sets, "RiskMap_test_sets")
  expect_identical(out$test_sets$sets[[2]]$train, 5:10)
  expect_false(out$refit)
  expect_identical(out$family, c(a = "gaussian", b = "gaussian"))
  expect_named(out$failures, c("model", "test_set", "message"))
  expect_identical(nrow(out$failures), 0L)
  expect_named(out$calibration, c("area", "curves", "simulated"))
  expect_identical(out$calibration$area$model, c("a", "b"))
  expect_identical(unique(out$calibration$curves$model), c("a", "b"))
  expect_identical(nrow(out$calibration$curves), 202L)
  expect_named(out$calibration$simulated, c("a", "b"))
})

test_that("the scores of Gaussian models are those of their normal predictions", {
  # Every predictive distribution has mean 0 and standard deviation 1.
  data <- line_data(y = c(0, 0.5, -1, 2), x1 = 1:4)
  local_mocked_bindings(
    predict_test_set = mock_draws(matrix(c(-1, 0, 1), 4, 3, byrow = TRUE))
  )
  out <- assess(make_assess_prediction_fit(data, "x1"), list(1, 2:4), seed = 3)
  y <- data$y
  crps <- 2 * dnorm(y) + y * (2 * pnorm(y) - 1) - 1 / sqrt(pi)
  expect_equal(out$scores$CRPS, crps)
  # The negative of Bolin and Wallin's scaled CRPS, a loss (#230).
  expect_equal(out$scores$SCRPS,
               0.5 * (1 + crps * sqrt(pi) + log(2 / sqrt(pi))))
  expect_equal(out$scores$pit_lower, pnorm(y))
  expect_identical(out$scores$pit_upper, out$scores$pit_lower)
  u <- seq(0, 1, length.out = 1001)
  expect_equal(out$calibration$area$PIT_area,
               pit_area(t(stats::ecdf(pnorm(y))(u)), u))
})

test_that("a Gaussian prediction without spread has the absolute error", {
  data <- line_data(y = c(0, 0.5, -1, 0), x1 = 1:4)
  local_mocked_bindings(predict_test_set = mock_constant(0.25))
  out <- assess(make_assess_prediction_fit(data, "x1"), list(1:3), seed = 1)
  expect_equal(out$scores$CRPS, abs(data$y[1:3] - 0.25))
  expect_identical(out$scores$SCRPS, rep(NA_real_, 3))
  expect_identical(out$scores$pit_upper, c(0, 1, 0))
})

test_that("a Gaussian prediction that isn't finite is a failure", {
  data <- line_data(y = c(0, 0.5, -1, 0), x1 = 1:4)
  fit <- make_assess_prediction_fit(data, "x1")
  draws <- matrix(rnorm(40), 4)
  draws[2, 3] <- Inf
  local_mocked_bindings(predict_test_set = mock_draws(draws))
  expect_warning(out <- assess(fit, list(1, 2:3)), "1 of 2 test-set")
  expect_identical(out$failures$test_set, 2L)
  expect_match(out$failures$message, "at 1 test observation isn't finite")
  expect_identical(is.na(out$scores$CRPS), c(FALSE, TRUE, TRUE))
  # One draw has no standard deviation.
  local_mocked_bindings(predict_test_set = mock_draws(matrix(0, 4, 1)))
  expect_warning(out <- assess(fit, list(1, 2:3)), "2 of 2 test-set")
  expect_match(out$failures$message[2], "at 2 test observations isn't finite")
})

test_that("a sharper calibrated prediction has the lower scores (#230)", {
  fit <- make_assess_prediction_fit(line_data(y = c(0, 0), x1 = 1:2), "x1")
  scores_with <- function(spread) {
    local_mocked_bindings(
      predict_test_set = mock_draws(matrix(spread * c(-1, 0, 1), 2, 3,
                                           byrow = TRUE))
    )
    assess(fit, list(2))$scores
  }
  expect_lt(scores_with(0.5)$CRPS, scores_with(2)$CRPS)
  expect_lt(scores_with(0.5)$SCRPS, scores_with(2)$SCRPS)
})

# A minimal binomial fit with a probit link.
make_probit_fit <- function() {
  fit <- make_assess_prediction_fit(
    line_data(y = c(3, 9, 6), n = c(10, 10, 10)), "1"
  )
  fit$family <- "binomial"
  fit$units_m <- fit$data$n
  fit$estimate <- list(beta = 0, sigma2 = 0, phi = 0)
  fit$link_function <- list(inv = function(x) pnorm(x), name = "custom")
  fit
}

test_that("counts are scored with the inverse link of the model", {
  local_mocked_bindings(predict_test_set = mock_constant(1))
  out <- assess(list(probit = make_probit_fit()), list(2:3), seed = 1)
  exact_crps <- function(y, p) {
    F_k <- pbinom(0:10, 10, p)
    sum((F_k - as.numeric(0:10 >= y))^2)
  }
  expect_equal(out$scores$CRPS,
               vapply(c(9, 6), exact_crps, numeric(1), p = pnorm(1)))
  expect_equal(out$scores$pit_lower, pbinom(c(8, 5), 10, pnorm(1)))
  expect_equal(out$scores$pit_upper, pbinom(c(9, 6), 10, pnorm(1)))
})

test_that("a failed refit leaves the scores of its test set missing", {
  captured <- new.env()
  fit <- make_probit_fit()
  local_mocked_bindings(
    glgpm = function(..., data, invlink = NULL) {
      captured$invlink <- invlink
      if (nrow(data) == 1L) stop("refit failed")
      fit
    },
    predict_test_set = mock_constant(1)
  )
  expect_warning(
    out <- assess(list(probit = fit), list(2:3, 1), refit = TRUE),
    "1 of 2 test-set assessments failed, and their scores are missing"
  )
  expect_equal(captured$invlink$inv(1), pnorm(1))
  expect_identical(out$failures$model, "probit")
  expect_identical(out$failures$message, "refit failed")
  expect_identical(out$failures$test_set, 1L)
  expect_identical(out$scores$row, c(2L, 3L, 1L))
  expect_identical(is.na(out$scores$CRPS), c(TRUE, TRUE, FALSE))
  expect_true(is.finite(out$calibration$area$PIT_area))
  expect_true(out$refit)
  # Nothing was scored.
  expect_warning(out <- assess(list(probit = fit), list(2:3), refit = TRUE),
                 "1 of 1")
  expect_true(is.na(out$calibration$area$PIT_area))
  expect_true(is.na(out$calibration$area$PIT_area_expected))
  expect_null(out$calibration$curves)
  expect_true(is.na(summary(out)["probit", "CRPS"]))
  expect_input_error(plot(out), "x", "no scores")
  expect_input_error(plot(out, type = "score"), "x", "no CRPS")
})

test_that("only the observations that every model can use are scored", {
  data <- line_data(y = c(1, 2, 3), x1 = c(1, 2, NA), unused = c(1, NA, NA))
  local_mocked_bindings(predict_test_set = mock_constant(0))
  expect_warning(
    out <- assess(list(a = make_assess_prediction_fit(data, "x1"),
                       b = make_assess_prediction_fit(data, "1")), list(2:3)),
    "1 observation with missing values"
  )
  expect_identical(out$scores$row, c(2L, 2L))
  expect_input_error(assess(make_assess_prediction_fit(data, "x1"), list(3)),
                     "test_sets", "no test observation without missing values")
  # A test set can be left without observations.
  expect_warning(
    out <- assess(make_assess_prediction_fit(data, "x1"), list(1:2, 3)),
    "1 observation with missing values"
  )
  expect_identical(out$scores$test_set, c(1L, 1L))
  expect_identical(nrow(out$failures), 0L)
})

test_that("the test sets have the seeds they had, shared by the models", {
  seeds <- integer(0)
  local_mocked_bindings(
    predict_test_set = function(fit, newdata, control_mcmc, seed) {
      seeds <<- c(seeds, seed)
      matrix(0, nrow(newdata), 2)
    }
  )
  fit <- gaussian_intercept_model
  set.seed(1)
  state <- .Random.seed
  assess(list(fit, fit), list(1, 2, 3), seed = 9)
  expect_identical(.Random.seed, state)
  set.seed(9)
  expect_identical(seeds, rep(sample.int(.Machine$integer.max, 3), 2))
})

test_that("the seed reproduces an assessment and keeps the session's numbers", {
  skip_on_cran()
  set.seed(1)
  state <- .Random.seed
  run <- function(...) {
    assess(list(a = gaussian_model, b = binomial_like), list(1:4, c(2, 7)), ...)
  }
  binomial_like <- gaussian_offset_model
  first <- run(seed = 7)
  expect_identical(.Random.seed, state)
  expect_identical(run(seed = 7), first)
  expect_false(identical(run(seed = 8)$scores, first$scores))
  # Without a seed the numbers are the session's.
  set.seed(3)
  one <- run()
  expect_false(identical(.Random.seed, state))
  set.seed(3)
  expect_identical(run(), one)
})

test_that("a model has the same scores in any position", {
  skip_on_cran()
  fit <- gaussian_offset_model
  pair <- assess(list(a = fit, b = fit), list(1:3, 4:6), seed = 5)
  columns <- c("CRPS", "SCRPS", "pit_lower", "pit_upper")
  of <- function(model) {
    x <- pair$scores[pair$scores$model == model, columns]
    row.names(x) <- NULL
    x
  }
  expect_identical(of("a"), of("b"))
  expect_identical(pair$calibration$simulated$a, pair$calibration$simulated$b)
  alone <- assess(list(b = fit), list(1:3, 4:6), seed = 5)
  expect_identical(alone$scores[columns], of("b"))
  expect_identical(alone$calibration$area$PIT_area_expected,
                   pair$calibration$area$PIT_area_expected[2])
})

test_that("assess_prediction() says how the parameters are set", {
  local_mocked_bindings(predict_test_set = mock_constant(0))
  fits <- list(m0 = gaussian_intercept_model, m1 = gaussian_offset_model)
  said <- function(models) {
    paste(capture_messages(assess_prediction(models, list(1:2))), collapse = " ")
  }
  expect_match(said(fits), "fixed at the estimates from all the data")
  expect_match(said(fits), "different numbers of parameters")
  expect_no_match(said(list(a = fits$m0, b = fits$m0)), "different numbers")
  expect_no_message(assess(fits, list(1:2)))
  local_mocked_bindings(glgpm = function(...) gaussian_intercept_model)
  expect_message(assess_prediction(fits, list(1:2, 3), refit = TRUE),
                 "re-estimated in each training set \\(4 model fits\\)")
})

test_that("a location with training and test observations is warned about", {
  data <- gaussian_data[c(1:10, 1:10), ]
  data$y <- data$y + rep(c(-0.3, 0.3, 0.1, -0.2), 5)
  fit <- glgpm(y ~ gp(), data = data, family = "gaussian", messages = FALSE)
  with_nugget <- glgpm(y ~ gp(nugget = TRUE), data = data, family = "gaussian",
                       fixed_pars = list(tau2 = 0.1), messages = FALSE)
  local_mocked_bindings(predict_test_set = mock_constant(0))
  expect_warning(assess(list(fit, with_nugget), list(c(2, 12), 11, c(4, 14))),
                 "In test set 2, some locations also have training")
  expect_no_warning(assess(fit, list(c(2, 12), 11)))
  expect_no_warning(assess(with_nugget, list(c(2, 12), c(1, 11))))
})

# The joint predictive mean and covariance of rows of a Gaussian fit, from the
# rest of its data.
exact_prediction <- function(fit, test, train = NULL) {
  p <- coef(fit)
  n <- length(fit$y)
  train <- train %||% setdiff(seq_len(n), test)
  spatial <- p$sigma2 * matern_correlation(
    dist(fit$coords[fit$ID_coords, ]), phi = p$phi, kappa = fit$kappa,
    return_sym_matrix = TRUE
  )
  same_site <- outer(fit$ID_coords, fit$ID_coords, "==")
  latent <- spatial + fitted_nugget_variance(fit) * same_site
  for (term in names(fit$re)) {
    latent <- latent + p$sigma2_re[[term]] *
      outer(fit$data[[term]], fit$data[[term]], "==")
  }
  observed <- latent + diag(p$sigma2_me %||% 0, n)
  mean <- as.numeric(fit$D %*% p$beta) + (fit$cov_offset %||% 0)
  weights <- observed[test, train] %*% solve(observed[train, train])
  list(mean = mean[test] +
         as.numeric(weights %*% (fit$y[train] - mean[train])),
       latent = latent[test, test] - weights %*% latent[train, test],
       outcome = observed[test, test] - weights %*% observed[train, test])
}

normal_crps <- function(y, mean, sd) {
  z <- (y - mean) / sd
  sd * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
}

# Twelve sites with three records each, in four groups of three sites, with a
# nugget, a random effect and a measurement error, at known parameters.
grouped_fit <- function() {
  set.seed(50)
  site <- rep(1:12, each = 3)
  coordinates <- cbind(runif(12, 0, 5e4), runif(12, 0, 5e4))[site, ]
  data <- st_as_sf(data.frame(
    x = coordinates[, 1], z = coordinates[, 2], covariate = rnorm(36),
    group = factor(rep(letters[1:4], each = 9)), y = rnorm(36, sd = 1.7)
  ), coords = c("x", "z"), crs = 32637)
  glgpm(y ~ covariate + gp(nugget = TRUE) + re(group), data = data,
        family = "gaussian", messages = FALSE,
        fixed_pars = list(beta = c(0.2, 0.5), sigma2 = 1, phi = 12, tau2 = 0.3,
                          sigma2_me = 0.25, sigma2_re = c(group = 0.8)))
}

test_that("test sets are predicted jointly, with the nugget and random effects", {
  skip_on_cran()
  fit <- grouped_fit()
  control <- set_control_mcmc(n_draws = 20000, burnin = 0, thin = 1)
  # The sites of a group, and one site of a group that the training data have.
  test <- 1:12
  exact <- exact_prediction(fit, test)
  link <- predict_test_set(condition_on_rows(fit, setdiff(1:36, test)),
                           fit$data[test, ], control, seed = 1)
  expect_identical(dim(link), c(12L, 20000L))
  expect_equal(rowMeans(link), exact$mean, tolerance = 0.03)
  expect_equal(stats::cov(t(link)), exact$latent, tolerance = 0.03,
               ignore_attr = TRUE)
})

test_that("the scores of Gaussian models with re() are those of their closed form", {
  skip_on_cran()
  fit <- grouped_fit()
  control <- set_control_mcmc(n_draws = 5000, burnin = 0, thin = 1)
  # A site of each group, whose other sites are in training, and whole groups,
  # whose levels are new.
  by_site <- lapply(0:2, function(k) which((rep(1:12, each = 3) - 1) %% 3 == k))
  by_group <- lapply(0:3, function(k) 9 * k + 1:9)
  for (sets in list(by_site, by_group)) {
    out <- suppressWarnings(assess(fit, sets, control_mcmc = control, seed = 2))
    expect_identical(nrow(out$failures), 0L)
    exact <- lapply(sets, function(test) {
      prediction <- exact_prediction(fit, test)
      sd <- sqrt(diag(prediction$outcome))
      data.frame(crps = normal_crps(fit$y[test], prediction$mean, sd),
                 pit = pnorm(fit$y[test], prediction$mean, sd))
    })
    exact <- do.call(rbind, exact)
    expect_equal(out$scores$CRPS, exact$crps, tolerance = 0.03)
    expect_equal(out$scores$pit_lower, exact$pit, tolerance = 0.03)
  }
})

test_that("the scores of counts with re() are those of their closed form", {
  skip_on_cran()
  # Almost no spatial variance, so that the effect of a group given its counts
  # has a density in one dimension.
  set.seed(51)
  data <- st_as_sf(data.frame(x = runif(12, 0, 5e4), z = runif(12, 0, 5e4),
                              group = rep(letters[1:3], each = 4),
                              y = c(9, 12, 7, 10, 0, 1, 0, 2, 3, 3, 5, 2),
                              n = rep(c(2, 1), 6)),
                   coords = c("x", "z"), crs = 32637)
  fit <- suppressWarnings(glgpm(
    y ~ gp() + re(group), data = data, family = "poisson", denominator = n,
    messages = FALSE, return_samples = FALSE, seed = 1,
    fixed_pars = list(beta = 0.8, sigma2 = 1e-6, phi = 1,
                      sigma2_re = c(group = 0.7))
  ))
  out <- assess(fit, list(c(1, 2, 5, 6), 9:12), seed = 3,
                control_mcmc = set_control_mcmc(n_draws = 4000, burnin = 2000,
                                                thin = 10))
  expect_identical(nrow(out$failures), 0L)
  distribution <- function(row, train) {
    same <- train[data$group[train] == data$group[row]]
    density <- function(effect) vapply(effect, function(e) {
      prod(dpois(data$y[same], data$n[same] * exp(0.8 + e))) *
        dnorm(e, sd = sqrt(0.7))
    }, numeric(1))
    grid <- seq(-5, 5, length.out = 2001)
    weight <- density(grid) / sum(density(grid))
    vapply(0:200, function(k) {
      sum(weight * ppois(k, data$n[row] * exp(0.8 + grid)))
    }, numeric(1))
  }
  exact <- vapply(seq_len(nrow(out$scores)), function(i) {
    row <- out$scores$row[i]
    values <- distribution(row, out$test_sets$sets[[out$scores$test_set[i]]]$train)
    c(sum((values - (0:200 >= data$y[row]))^2), values[data$y[row] + 1])
  }, numeric(2))
  expect_equal(out$scores$CRPS, exact[1, ], tolerance = 0.04)
  expect_equal(out$scores$pit_upper, exact[2, ], tolerance = 0.04)
})

test_that("the reference of the calibration has joint draws", {
  # Draws that are the same for every observation give curves that step from 0
  # to 1 at a uniform, whose expected area is 1/3, and independent draws about
  # 0.31 / sqrt(n).
  set.seed(52)
  n <- 50
  # The last row is for training.
  fit <- make_assess_prediction_fit(
    line_data(y = rep(0, n + 1), x1 = seq_len(n + 1)), "x1"
  )
  expected_with <- function(draws, sets = list(1:n)) {
    local_mocked_bindings(predict_test_set = mock_draws(draws))
    assess(fit, sets, seed = 1)$calibration$area$PIT_area_expected
  }
  shared <- matrix(rnorm(1000), n, 1000, byrow = TRUE)
  expect_equal(expected_with(shared), 1 / 3, tolerance = 0.05)
  independent <- matrix(rnorm(n * 1000), n)
  expect_equal(expected_with(independent), 0.31 / sqrt(n), tolerance = 0.1)
  # The test sets are independent of each other. Here they are half the data
  # each, which move together within each half.
  halves <- rbind(shared[1:25, ], matrix(rnorm(1000), 25, 1000, byrow = TRUE))
  two <- expected_with(halves, list(1:25, 26:50))
  expect_gt(two, 0.2)
  expect_lt(two, 0.3)
})

test_that("an observation in two test sets counts as one in the reference", {
  set.seed(53)
  n <- 12
  fit <- make_assess_prediction_fit(
    line_data(y = rnorm(n + 1), x1 = seq_len(n + 1)), "x1"
  )
  draws <- matrix(rnorm(n * 500, 0.5, 0.8), n)
  local_mocked_bindings(predict_test_set = mock_draws(draws))
  once <- assess(fit, list(1:n), seed = 1)$calibration
  for (sets in list(list(1:n, 1:n), list(1:n, 1:6, 7:n))) {
    again <- assess(fit, sets, seed = 1)$calibration
    expect_equal(again$area, once$area)
    expect_lt(max(abs(again$simulated[[1]] - once$simulated[[1]])), 1e-12)
  }

  # Counts with the same means in both test sets are the same counts.
  traits <- family_traits("poisson")
  simulate <- function(rows, earlier, shift = 0) {
    means <- exp(draws[rows, , drop = FALSE] + shift)
    scored <- score_counts(rep(2, length(rows)), rep(1, length(rows)), means,
                           NULL, traits)
    simulate_count_pit(c(list(draws = means[, 1:200, drop = FALSE]),
                         scored$predictive), traits, earlier)
  }
  first <- simulate(1:n, vector("list", n))
  expect_true(all(first$upper > first$lower))
  second <- simulate(c(7, 3), first$uniforms[c(7, 3)])
  expect_equal(second$lower, first$lower[, c(7, 3)])
  expect_equal(second$upper, first$upper[, c(7, 3)])
  expect_identical(second$uniforms, first$uniforms[c(7, 3)])
  # With other means they have the same PIT value.
  moved <- simulate(4L, first$uniforms[4], shift = 1)
  expect_true(all(moved$lower <= first$uniforms[[4]] &
                    first$uniforms[[4]] <= moved$upper))
  expect_false(isTRUE(all.equal(moved$lower, first$lower[, 4, drop = FALSE])))
})

test_that("assess_prediction() refits the models of every family", {
  skip_on_cran()
  for (fit in list(gaussian_model, binomial_model, poisson_model)) {
    out <- suppressWarnings(assess(
      list(model = fit), list(2), refit = TRUE, control_mcmc = control_mcmc,
      seed = fixture_seed
    ))
    expect_identical(nrow(out$failures), 0L)
    expect_true(is.finite(out$scores$CRPS))
    expect_true(out$refit)
  }
})

test_that("refits keep a custom CRS", {
  custom_data <- st_transform(
    gaussian_data,
    st_crs("+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs")
  )
  fit <- glgpm(y ~ cov + gp(), data = custom_data, family = "gaussian",
               messages = FALSE)
  expect_true(is.na(st_crs(fit$data)$epsg))
  out <- assess(fit, list(1:2), refit = TRUE, control_mcmc = control_mcmc,
                seed = fixture_seed)
  expect_identical(nrow(out$failures), 0L)
})
