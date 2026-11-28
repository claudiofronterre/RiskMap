# Predictions compared with closed forms and independent calculations, at given
# parameters so that nothing is estimated.

held_gaussian <- list(beta = c(1, 0.5), sigma2 = 1, phi = 20, sigma2_me = 0.25)
held_re <- c(g = 1.5, h = 0.8, k = 0.4)

closed_form_sf <- function(coordinates, values) {
  values$east <- coordinates[, 1] * 1000 + 5e5
  values$north <- coordinates[, 2] * 1000
  sf::st_as_sf(values, coords = c("east", "north"), crs = 32633)
}

# Records at 40 locations, with several at some of them, and three groupings
# that cross each other and the locations.
make_records <- function(seed = 11, n = 60, n_locations = 40) {
  set.seed(seed)
  coordinates <- cbind(runif(n_locations, 0, 100), runif(n_locations, 0, 100))
  location <- c(seq_len(n_locations), sample(n_locations, n - n_locations))
  groups <- function(n_groups) {
    match(x <- c(seq_len(n_groups), sample(n_groups, n - n_groups, TRUE)),
          unique(x))
  }
  g <- groups(6)
  h <- groups(5)
  k <- groups(4)
  z <- rnorm(n)
  spatial <- as.numeric(t(chol(exp(-as.matrix(dist(coordinates)) / 20))) %*%
                          rnorm(n_locations))
  y <- 1 + 0.5 * z + spatial[location] + rnorm(6, 0, sqrt(1.5))[g] +
    rnorm(5, 0, sqrt(0.8))[h] + rnorm(4, 0, sqrt(0.4))[k] + rnorm(n, 0, 0.5)
  list(coordinates = coordinates[location, ],
       data = closed_form_sf(coordinates[location, ],
                             data.frame(y, z, g, h, k)))
}

fit_held <- function(data, terms = character(), nugget = FALSE, ...) {
  formula <- reformulate(
    c("z", if (nugget) "gp(nugget = TRUE)" else "gp()",
      if (length(terms)) paste0("re(", paste(terms, collapse = ", "), ")")),
    "y"
  )
  held <- utils::modifyList(held_gaussian, list(...))
  if (length(terms)) held$sigma2_re <- held_re[terms]
  if (nugget && is.null(held$tau2)) held$tau2 <- 0.3
  glgpm(formula, data = data, family = "gaussian",
        fixed_pars = held, messages = FALSE)
}

# The distribution, given the data, of the linear predictor and of its spatial
# and random-effect parts at the rows of `newdata`.
exact_gaussian <- function(fit, coordinates, new_coordinates, newdata,
                           terms = character(), new_nugget = FALSE) {
  parameters <- coef(fit)
  data <- sf::st_drop_geometry(fit$data)
  new <- sf::st_drop_geometry(newdata)
  distances <- function(a, b) {
    sqrt(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2)
  }
  spatial <- function(a, b) {
    parameters$sigma2 * exp(-distances(a, b) / parameters$phi)
  }
  random <- function(first, second) {
    total <- 0
    for (term in terms) {
      total <- total + parameters$sigma2_re[[term]] *
        outer(as.character(first[[term]]), as.character(second[[term]]), "==")
    }
    total
  }
  tau2 <- parameters$tau2
  if (is.null(tau2)) tau2 <- 0
  nugget <- function(a) tau2 * (distances(a, a) == 0)
  V <- spatial(coordinates, coordinates) + nugget(coordinates) +
    random(data, data) + diag(parameters$sigma2_me, nrow(data))
  residual <- fit$y - as.numeric(fit$D %*% parameters$beta)
  cross <- list(spatial = spatial(new_coordinates, coordinates),
                random = random(new, data) + matrix(0, nrow(new), nrow(data)))
  prior <- list(spatial = spatial(new_coordinates, new_coordinates),
                random = random(new, new) + matrix(0, nrow(new), nrow(new)))
  conditional <- function(first, second) {
    - first %*% solve(V, t(second))
  }
  fixed <- as.numeric(parameters$beta[[1]] + parameters$beta[[2]] * new$z)
  out <- lapply(stats::setNames(nm = names(cross)), function(part) {
    list(mean = as.numeric(cross[[part]] %*% solve(V, residual)),
         covariance = prior[[part]] +
           conditional(cross[[part]], cross[[part]]))
  })
  both <- cross$spatial + cross$random
  out$link <- list(
    mean = fixed + as.numeric(both %*% solve(V, residual)),
    covariance = prior$spatial + prior$random + conditional(both, both) +
      if (new_nugget) nugget(new_coordinates) else 0
  )
  out
}

# Draws within Monte Carlo error of a Gaussian distribution, with 5 standard
# errors for the means and for the entries of the covariance matrix.
expect_gaussian_draws <- function(draws, exact) {
  n <- ncol(draws)
  variances <- diag(exact$covariance)
  expect_lt(max(abs(rowMeans(draws) - exact$mean) /
                  sqrt(pmax(variances, 1e-8) / n)), 5)
  error <- sqrt((outer(variances, variances) + exact$covariance^2) / n) + 1e-8
  expect_lt(max(abs(stats::cov(t(draws)) - exact$covariance) / error), 5)
}

draws_of <- function(prediction) {
  samples <- prediction$samples
  link <- prediction$predictions$fixed + prediction$predictions$offset +
    samples$spatial
  for (part in c("random_effects", "nugget")) {
    if (!is.null(samples[[part]])) link <- link + samples[[part]]
  }
  list(link = link, spatial = samples$spatial,
       random = samples$random_effects, nugget = samples$nugget)
}

many_draws <- set_control_mcmc(n_draws = 10000, burnin = 0, thin = 1)

new_records <- function(seed = 5) {
  set.seed(seed)
  coordinates <- cbind(runif(4, 0, 100), runif(4, 0, 100))
  list(coordinates = coordinates,
       data = closed_form_sf(coordinates, data.frame(
         z = rnorm(4), g = c(1, 2, 2, 5), h = c(1, 1, 3, 4), k = c(2, 1, 4, 4)
       )))
}

test_that("the labels of Gaussian random effects don't change predictions", {
  skip_on_cran()
  records <- make_records()
  new <- new_records()
  data <- records$data
  reference <- fit_held(data, "g")
  exact <- exact_gaussian(reference, records$coordinates, new$coordinates,
                          new$data, "g")
  predict_with <- function(fit, newdata, ...) {
    predict(fit, newdata, seed = 1, control_mcmc = many_draws,
            keep_samples = TRUE, messages = FALSE, ...)
  }
  plain <- predict_with(reference, new$data, joint = TRUE,
                        include = "random_effects")
  expect_gaussian_draws(draws_of(plain)$spatial, exact$spatial)

  other <- c(3, 6, 1, 5, 2, 4)
  labels <- list(
    character = function(code) letters[code],
    factor = function(code) factor(letters[other[code]], levels = letters[1:6]),
    from_101 = function(code) code + 100,
    permuted = function(code) other[code],
    numeric = function(code) other[code] / 10
  )
  for (kind in names(labels)) {
    labelled <- data
    labelled$g <- labels[[kind]](data$g)
    new_labelled <- new$data
    new_labelled$g <- labels[[kind]](new$data$g)
    fit <- fit_held(labelled, "g")
    # The same noise, so the spatial draws are the same whatever the labels.
    marginal <- predict_with(fit, new_labelled)
    expect_equal(marginal$samples$spatial,
                 predict_with(reference, new$data)$samples$spatial,
                 tolerance = 1e-8, info = kind)
    joint <- predict_with(fit, new_labelled, joint = TRUE,
                          include = "random_effects")
    expect_equal(joint$samples$spatial, plain$samples$spatial,
                 tolerance = 1e-8, info = kind)
    # The levels of a factor are in another order, so the draws of the
    # effects are others of the same distribution.
    expect_gaussian_draws(draws_of(joint)$link, exact$link)
  }
})

test_that("the labels of random effects don't change count predictions", {
  skip_on_cran()
  records <- make_records()
  data <- records$data
  set.seed(4)
  eta <- as.numeric(data$y - mean(data$y))
  data$trials <- 20
  data$binomial <- rbinom(nrow(data), 20, plogis(eta))
  data$counts <- rpois(nrow(data), exp(1 + eta / 2))
  held <- list(beta = c(0, 0.5), sigma2 = 1, phi = 20, sigma2_re = c(g = 1.5))
  new <- new_records()$data
  families <- list(binomial = "binomial", poisson = "poisson",
                   negative_binomial = negative_binomial())
  for (name in names(families)) {
    predict_with <- function(labels) {
      labelled <- data
      labelled$g <- labels(data$g)
      newdata <- new
      newdata$g <- labels(new$g)
      labelled$response <- if (name == "binomial") data$binomial else
        data$counts
      values <- held
      if (name == "negative_binomial") values$theta <- 2
      fit <- glgpm(
        response ~ z + gp() + re(g), data = labelled, family = families[[name]],
        denominator = trials, fixed_pars = values, seed = 3,
        control_mcmc = set_control_mcmc(n_draws = 100, burnin = 100, thin = 1),
        messages = FALSE
      )
      predict(fit, newdata, joint = TRUE, include = "random_effects", seed = 1,
              keep_samples = TRUE, messages = FALSE)$samples
    }
    reference <- predict_with(identity)
    other <- c(3, 6, 1, 5, 2, 4)
    expect_equal(predict_with(function(code) letters[code]), reference)
    expect_equal(predict_with(function(code) code + 100), reference)
    expect_equal(predict_with(function(code) other[code]), reference)
  }
})

test_that("joint Gaussian predictions with random effects are exact", {
  skip_on_cran()
  records <- make_records()
  new <- new_records()
  area <- sf::st_sf(geometry = sf::st_as_sfc(
    sf::st_bbox(sf::st_buffer(new$data, 1000))
  ))
  for (terms in list("g", c("g", "h"), c("g", "h", "k"))) {
    for (nugget in c(FALSE, TRUE)) {
      fit <- fit_held(records$data, terms, nugget = nugget)
      include <- c("random_effects", if (nugget) "nugget")
      exact <- exact_gaussian(fit, records$coordinates, new$coordinates,
                              new$data, terms, new_nugget = nugget)
      prediction <- predict(
        fit, new$data, joint = TRUE, include = include, seed = 2,
        control_mcmc = many_draws, keep_samples = TRUE, messages = FALSE
      )
      draws <- draws_of(prediction)
      expect_gaussian_draws(draws$spatial, exact$spatial)
      expect_gaussian_draws(draws$random, exact$random)
      expect_gaussian_draws(draws$link, exact$link)

      areal <- predict(
        fit, new$data, boundaries = area, aggregate = "mean",
        include = include, seed = 2, control_mcmc = many_draws,
        keep_samples = TRUE, messages = FALSE
      )
      mean_of <- function(part) {
        list(mean = mean(part$mean), covariance = matrix(mean(part$covariance)))
      }
      expect_gaussian_draws(areal$samples$mean$link, mean_of(exact$link))
      expect_gaussian_draws(areal$samples$mean$random_effects,
                            mean_of(exact$random))
    }
  }
})

test_that("random effects are kept when the measurement error is held at zero", {
  records <- make_records(n = 40)
  new <- new_records()
  area <- sf::st_sf(geometry = sf::st_as_sfc(
    sf::st_bbox(sf::st_buffer(new$data, 1000))
  ))
  for (terms in list("g", c("g", "h", "k"))) {
    fit <- fit_held(records$data, terms, sigma2_me = 0)
    exact <- exact_gaussian(fit, records$coordinates, new$coordinates,
                            new$data, terms)
    # Marginal draws have the same means and variances.
    marginal <- predict(fit, new$data, seed = 2, control_mcmc = many_draws,
                        keep_samples = TRUE, messages = FALSE)$samples$spatial
    expect_gaussian_draws(marginal, list(
      mean = exact$spatial$mean,
      covariance = diag(diag(exact$spatial$covariance))
    ))
    prediction <- predict(
      fit, new$data, joint = TRUE, include = "random_effects", seed = 2,
      control_mcmc = many_draws, keep_samples = TRUE, messages = FALSE
    )
    expect_gaussian_draws(draws_of(prediction)$link, exact$link)
    # The data determine the sums of the effects of their records, so the
    # covariance of the effects needs a jitter, which warns.
    expect_warnings(
      areal <- predict(
        fit, new$data, boundaries = area, aggregate = "mean",
        include = "random_effects", seed = 2, control_mcmc = many_draws,
        keep_samples = TRUE, messages = FALSE
      ),
      "needed a diagonal jitter of 1e-12 to be factorised"
    )
    expect_gaussian_draws(areal$samples$mean$link, list(
      mean = mean(exact$link$mean),
      covariance = matrix(mean(exact$link$covariance))
    ))
  }
})

test_that("joint predictions share the effects of a location and of a level", {
  records <- make_records()
  new <- new_records()
  # Two locations with two and three records, and one with one, in mixed order
  rows <- c(1, 2, 1, 3, 2, 2)
  coordinates <- new$coordinates[rows, ]
  newdata <- closed_form_sf(coordinates, data.frame(
    z = c(0.2, -0.4, 1, 0, 0.3, 0.1), g = c(1, 2, 5, 5, 2, 6),
    h = c(1, 1, 2, 4, 5, 5), k = 1
  ))
  fit <- fit_held(records$data, c("g", "h"), nugget = TRUE)
  prediction <- expect_no_warning(predict(
    fit, newdata, joint = TRUE, include = c("random_effects", "nugget"),
    seed = 2, control_mcmc = many_draws, keep_samples = TRUE, messages = FALSE
  ))
  draws <- draws_of(prediction)
  exact <- exact_gaussian(fit, records$coordinates, coordinates, newdata,
                          c("g", "h"), new_nugget = TRUE)
  expect_gaussian_draws(draws$link, exact$link)
  expect_identical(draws$nugget[c(1, 2, 2), ], draws$nugget[c(3, 5, 6), ])
  expect_identical(draws$spatial[c(1, 2, 2), ], draws$spatial[c(3, 5, 6), ])
  expect_lt(abs(cor(draws$nugget[1, ], draws$nugget[2, ])), 0.05)

  area <- sf::st_sf(geometry = sf::st_as_sfc(
    sf::st_bbox(sf::st_buffer(newdata, 1000))
  ))
  areal <- expect_no_warning(predict(
    fit, newdata, boundaries = area, aggregate = "mean",
    include = c("random_effects", "nugget"), seed = 2,
    control_mcmc = many_draws, keep_samples = TRUE, messages = FALSE
  ))
  expect_gaussian_draws(areal$samples$mean$link, list(
    mean = mean(exact$link$mean),
    covariance = matrix(mean(exact$link$covariance))
  ))
})

test_that("a level that isn't in the data stops a prediction with its effect", {
  records <- make_records()
  fit <- fit_held(records$data, c("g", "h"))
  newdata <- new_records()$data
  newdata$g[2] <- 99
  expect_input_error(
    predict(fit, newdata, joint = TRUE, include = "random_effects",
            messages = FALSE),
    "newdata"
  )
  # It doesn't matter when the effect isn't included.
  expect_s3_class(predict(fit, newdata, joint = TRUE, messages = FALSE),
                  "RiskMap_prediction")
})

test_that("effects of new levels are drawn from their distribution inside", {
  skip_on_cran()
  # What the assessment of predictions needs, where a test set can have levels
  # of its own and several records at a location.
  records <- make_records()
  new <- new_records()
  rows <- c(1, 2, 1, 3, 2)
  coordinates <- new$coordinates[rows, ]
  newdata <- closed_form_sf(coordinates, data.frame(
    z = c(0.2, -0.4, 1, 0, 0.3), g = c(1, 98, 99, 2, 98), h = c(98, 1, 1, 2, 99)
  ))
  include <- c("random_effects", "nugget")
  components <- function(fit, control) {
    input <- RiskMap:::validate_prediction_data(
      fit, newdata, include, allow_new_levels = TRUE
    )
    set.seed(3)
    RiskMap:::prediction_sample_components(
      fit, input$model_data, "joint", include, control, messages = FALSE
    )
  }
  fit <- fit_held(records$data, c("g", "h"), nugget = TRUE)
  exact <- exact_gaussian(fit, records$coordinates, coordinates, newdata,
                          c("g", "h"), new_nugget = TRUE)
  expect_gaussian_draws(components(fit, many_draws)$link, exact$link)

  data <- records$data
  data$g <- letters[data$g]
  newdata$g <- c("a", "new", "other", "b", "new")
  data$trials <- 20
  data$binomial <- rbinom(nrow(data), 20, plogis(data$y - mean(data$y)))
  data$counts <- rpois(nrow(data), exp(1 + (data$y - mean(data$y)) / 2))
  families <- list(binomial = "binomial", poisson = "poisson",
                   negative_binomial = negative_binomial())
  for (name in names(families)) {
    data$response <- if (name == "binomial") data$binomial else data$counts
    held <- list(beta = c(0, 0.5), sigma2 = 1, phi = 20, tau2 = 0.3,
                 sigma2_re = c(g = 1.5, h = 0.8))
    if (name == "negative_binomial") held$theta <- 2
    fit <- glgpm(
      response ~ z + gp(nugget = TRUE) + re(g, h), data = data,
      family = families[[name]], denominator = trials, fixed_pars = held,
      seed = 3, messages = FALSE,
      control_mcml = set_control_mcml(integrate_noise = FALSE),
      control_mcmc = set_control_mcmc(n_draws = 1500, burnin = 500, thin = 5)
    )
    draws <- components(fit, NULL)
    n <- ncol(draws$link)
    expect_identical(draws$nugget[c(1, 2), ], draws$nugget[c(3, 5), ])
    expect_identical(draws$spatial[c(1, 2), ], draws$spatial[c(3, 5), ])
    expect_equal(draws$link, draws$fixed + draws$spatial +
                   draws$random_effects + draws$nugget)
    # The effects given the data of the known levels, and for the new ones
    # draws from their distribution, shared by the rows of one level.
    samples <- fit$S_samples
    n_locations <- nrow(fit$coords)
    g <- function(level) samples[, n_locations + match(level, fit$re$g)]
    h <- function(level) {
      samples[, n_locations + length(fit$re$g) + match(level, fit$re$h)]
    }
    known <- rbind(draws$random_effects[4, ] - g("b") - h(2),
                   draws$random_effects[1, ] - g("a"),
                   draws$random_effects[2, ] - h(1),
                   draws$random_effects[3, ] - h(1))
    expect_equal(known[1, ], numeric(n))
    expect_lt(abs(var(known[2, ]) / 0.8 - 1), 5 * sqrt(2 / n))
    expect_lt(abs(var(known[3, ]) / 1.5 - 1), 5 * sqrt(2 / n))
    expect_lt(abs(var(known[4, ]) / 1.5 - 1), 5 * sqrt(2 / n))
    expect_lt(max(abs(cor(t(known[2:4, ]))[upper.tri(diag(3))])), 5 / sqrt(n))
    fifth <- draws$random_effects[5, ] - known[3, ]
    expect_lt(abs(var(fifth) / 0.8 - 1), 5 * sqrt(2 / n))
    expect_lt(max(abs(cor(fifth, t(known[2:4, ])))), 5 / sqrt(n))
  }
})

test_that("predict() without newdata gives each record its location's effect", {
  skip_on_cran()
  records <- make_records()
  data <- records$data
  order <- sample(nrow(data))
  data <- data[order, ]
  data$trials <- 20
  set.seed(8)
  data$binomial <- rbinom(nrow(data), 20, plogis(data$y - mean(data$y)))
  data$counts <- rpois(nrow(data), exp(1 + (data$y - mean(data$y)) / 2))
  control <- set_control_mcmc(n_draws = 60, burnin = 100, thin = 1)
  families <- list(binomial = "binomial", poisson = "poisson",
                   negative_binomial = negative_binomial())
  for (name in names(families)) {
    data$response <- if (name == "binomial") data$binomial else data$counts
    held <- list(beta = c(0, 0.5), sigma2 = 1, phi = 20)
    if (name == "negative_binomial") held$theta <- 2
    for (with_re in c(FALSE, TRUE)) {
      formula <- if (with_re) response ~ z + gp() + re(g) else response ~ z + gp()
      if (with_re) held$sigma2_re <- c(g = 1.5)
      fit <- glgpm(
        formula, data = data, family = families[[name]], denominator = trials,
        fixed_pars = held, seed = 3, control_mcmc = control, messages = FALSE
      )
      prediction <- expect_no_warning(predict(
        fit, seed = 1, joint = with_re,
        include = if (with_re) "random_effects", keep_samples = TRUE,
        messages = FALSE
      ))
      expect_equal(prediction$predictions$link_mean, fitted(fit, "link"))
      expect_identical(nrow(prediction$samples$spatial), nrow(data))
      expect_identical(
        prediction$samples$spatial,
        t(fit$S_samples[, fit$ID_coords])
      )
      exceedance <- predict(fit, seed = 1, summaries = c("mean", "exceedance"),
                            thresholds = 0.3, messages = FALSE)
      expect_identical(nrow(exceedance$predictions), nrow(data))
    }
  }

  # With `standardise`
  standardised <- predict(fit, standardise = data.frame(z = c(0, 1)), seed = 1,
                          messages = FALSE)
  profiles <- rowMeans(vapply(c(0, 1), function(z) {
    with_z <- fit$S_samples[, fit$ID_coords] + 0.5 * z
    colMeans(exp(with_z))
  }, numeric(nrow(data))))
  expect_equal(standardised$predictions$response_mean, profiles)
})

test_that("predict() without newdata follows the records at every draw", {
  skip_on_cran()
  data <- make_records()$data
  data <- data[sample(nrow(data)), ]
  data$counts <- rpois(nrow(data), exp(1 + (data$y - mean(data$y)) / 2))
  covariance <- list(sigma2 = 1, phi = 20)
  fits <- list(
    glgpm(y ~ z + gp(), data = data, family = "gaussian",
          fixed_pars = c(covariance, sigma2_me = 0.25), messages = FALSE),
    glgpm(
      counts ~ z + gp(), data = data, family = "poisson", denominator = NULL,
      fixed_pars = covariance, seed = 3, messages = FALSE,
      control_mcmc = set_control_mcmc(n_draws = 60, burnin = 100, thin = 1)
    )
  )
  for (fit in fits) {
    prediction <- expect_no_warning(predict(
      fit, parameter_uncertainty = TRUE, parameter_draws = 3, seed = 1,
      keep_samples = TRUE, messages = FALSE
    ))
    samples <- prediction$samples
    expect_identical(nrow(samples$spatial), nrow(data))
    first <- match(fit$ID_coords, fit$ID_coords)
    expect_identical(samples$spatial, samples$spatial[first, ])
    # The coefficients differ between the draws, and the records of a location
    # in their covariate.
    expect_identical(nrow(samples$fixed), nrow(data))
    each <- ncol(samples$fixed) / 3
    slopes <- apply(samples$fixed[, c(1, 1 + each, 1 + 2 * each)], 2,
                    function(fixed) coef(lm(fixed ~ data$z))[[2]])
    expect_length(unique(round(slopes, 10)), 3L)
    expect_equal(prediction$predictions$link_mean,
                 rowMeans(samples$fixed + samples$spatial))
  }
})

test_that("Gaussian predictions without newdata follow the records", {
  skip_on_cran()
  records <- make_records()
  order <- sample(nrow(records$data))
  data <- records$data[order, ]
  coordinates <- records$coordinates[order, ]
  for (terms in list(character(), c("g", "h"))) {
    fit <- fit_held(data, terms)
    exact <- exact_gaussian(fit, coordinates, coordinates, data, terms)
    for (joint in c(FALSE, TRUE)) {
      prediction <- expect_no_warning(predict(
        fit, joint = joint, seed = 1, control_mcmc = many_draws,
        include = if (joint && length(terms)) "random_effects",
        keep_samples = TRUE, messages = FALSE
      ))
      draws <- draws_of(prediction)
      target <- exact[[if (joint) "link" else "spatial"]]
      values <- if (joint) draws$link else draws$spatial
      covariance <- target$covariance
      # Records at one location share their draws.
      if (!joint) covariance <- covariance * outer(fit$ID_coords, fit$ID_coords, "==")
      expect_gaussian_draws(values, list(mean = target$mean,
                                         covariance = covariance))
    }
    # The rows are those of the data, and the means the fitted values.
    expect_equal(rowMeans(draws$link), fitted(fit), tolerance = 0.02)
  }
})

test_that("the stored samples are drawn at the held values", {
  skip_on_cran()
  set.seed(21)
  n <- 50
  coordinates <- cbind(runif(n, 0, 100), runif(n, 0, 100))
  spatial <- as.numeric(t(chol(exp(-as.matrix(dist(coordinates)) / 20))) %*%
                          rnorm(n))
  data <- closed_form_sf(coordinates, data.frame(
    z = 1000 + rnorm(n, 0, 100), trials = 20, g = rep(1:10, 5)
  ))
  eta <- -1 + (data$z - 1000) / 200 + spatial + rnorm(10, 0, 0.7)[data$g]
  data$binomial <- rbinom(n, 20, plogis(eta))
  data$counts <- rnbinom(n, size = 3, mu = exp(1 + eta / 2))
  all_held <- list(beta = c("(Intercept)" = -6, z = 0.005), sigma2 = 0.9,
                   phi = 18, tau2 = 0.2, sigma2_re = c(g = 0.4))
  # Each kind of hold alone, the two coefficients together and one at a time,
  # and every parameter.
  holds <- list(
    phi = all_held["phi"], beta = all_held["beta"],
    slope = list(beta = c(z = 0.005)),
    intercept = list(beta = c("(Intercept)" = -6)),
    sigma2 = all_held["sigma2"], tau2 = all_held["tau2"],
    sigma2_re = all_held["sigma2_re"],
    several = list(phi = 18, beta = c(z = 0.005), sigma2_re = c(g = 0.4)),
    all = all_held
  )
  control_mcml <- set_control_mcml(
    initial_draws = 300, validation_draws = 300, max_draws = 300,
    max_stages = 1, integrate_noise = FALSE
  )
  control <- set_control_mcmc(n_draws = 50, burnin = 100, thin = 1)
  check_stored <- function(fit) {
    parameters <- coef(fit)
    Sigma <- parameters$sigma2 *
      exp(-as.matrix(dist(coordinates)) / parameters$phi)
    diag(Sigma) <- diag(Sigma) + parameters$tau2
    direct <- RiskMap:::laplace_sampling_mcmc(
      y = fit$y, units_m = fit$units_m,
      mu = as.numeric(fit$D %*% parameters$beta), Sigma = Sigma,
      ID_coords = fit$ID_coords, ID_re = fit$ID_re,
      sigma2_re = parameters$sigma2_re, family = fit$family,
      control_mcmc = control,
      seed = RiskMap:::final_sampling_seed(fit$seed, length(fit$mcml_history)),
      family_parameters = parameters["theta"]
    )
    expect_equal(unname(fit$S_samples), unname(do.call(cbind, direct$samples)),
                 tolerance = 1e-5)
  }
  for (hold in names(holds)) {
    fit <- glgpm(
      binomial ~ z + gp(nugget = TRUE) + re(g), data = data,
      family = "binomial", denominator = trials, fixed_pars = holds[[hold]],
      control_mcml = control_mcml, control_mcmc = control, seed = 5,
      messages = FALSE
    )
    check_stored(fit)
    # The fitted values and the default prediction are from them.
    expect_equal(fitted(fit, "link"), colMeans(
      fit$S_samples[, fit$ID_coords] +
        fit$S_samples[, n + fit$ID_re$g]
    ) + as.numeric(fit$D %*% coef(fit)$beta))
  }
  # Every parameter of the families with counts, with the size held
  for (family in list("poisson", negative_binomial())) {
    held <- all_held
    held$beta <- c("(Intercept)" = -1, z = 0.002)
    if (!is.character(family)) held$theta <- 3
    check_stored(glgpm(
      counts ~ z + gp(nugget = TRUE) + re(g), data = data, family = family,
      fixed_pars = held, control_mcml = control_mcml, control_mcmc = control,
      seed = 5, messages = FALSE
    ))
  }
})

test_that("a Gaussian variance estimated at zero keeps its estimate", {
  skip_on_cran()
  set.seed(3)
  n <- 45
  coordinates <- cbind(runif(n, 0, 100), runif(n, 0, 100))
  root <- t(chol(exp(-as.matrix(dist(coordinates)) / 20)))
  location <- rep(seq_len(n), 2)
  draw <- function(seed) {
    set.seed(seed)
    closed_form_sf(coordinates[location, ], data.frame(
      y = 1 + as.numeric(root %*% rnorm(n))[location] + rnorm(2 * n, 0, 0.3)
    ))
  }
  # The first of these data sets whose nugget is estimated at zero
  for (seed in 1:20) {
    data <- draw(seed)
    fit <- glgpm(y ~ gp(nugget = TRUE), data = data,
                 family = "gaussian", messages = FALSE)
    if (coef(fit)$tau2 < 1e-4) break
  }
  expect_lt(coef(fit)$tau2, 1e-4)
  expect_gt(sqrt(diag(fit$covariance))[["nu2"]], 100)
  without <- glgpm(y ~ gp(), data = data, family = "gaussian",
                   messages = FALSE)
  set.seed(1)
  draws <- RiskMap:::draw_parameter_values(fit, 500)
  # The nugget keeps its estimate, as its variance.
  expect_equal(draws[, "nu2"] + draws[, "sigma2"],
               rep(log(coef(fit)$tau2), 500))
  # The other parameters are drawn as in the model without a nugget.
  set.seed(1)
  reference <- RiskMap:::draw_parameter_values(without, 500)
  expect_equal(draws[, colnames(reference)], reference, tolerance = 0.01)

  newdata <- new_records()$data
  for (seed in 1:5) {
    prediction <- expect_no_warning(predict(
      fit, newdata, parameter_uncertainty = TRUE, parameter_draws = 20,
      seed = seed, summaries = c("mean", "sd"), messages = FALSE
    ))
    plain <- predict(without, newdata, parameter_uncertainty = TRUE,
                     parameter_draws = 20, seed = seed,
                     summaries = c("mean", "sd"), messages = FALSE)
    expect_equal(prediction$predictions$link_sd, plain$predictions$link_sd,
                 tolerance = 0.02)
  }
})

test_that("the other parameters are drawn given those that are kept", {
  fit <- glgpm(
    y ~ z + gp(nugget = TRUE), data = make_records()$data, family = "gaussian",
    messages = FALSE
  )
  # A nugget near zero with a wide likelihood, correlated with the spatial
  # variance and scale, on the scale with the variance of the nugget
  fit$estimate$sigma2 <- log(0.8)
  fit$estimate$nu2 <- log(0.004 / 0.8)
  names_flat <- names(unlist(fit$estimate))
  set.seed(1)
  root <- diag(c(0.2, 0.1, 0.4, 0.6, 12, 0.3))
  root[5, 3] <- 1.5
  root[5, 4] <- 1.2
  root[5, 1] <- 0.4
  reported <- crossprod(root)
  dimnames(reported) <- list(names_flat, names_flat)
  to_working <- diag(6)
  dimnames(to_working) <- dimnames(reported)
  to_working["nu2", "sigma2"] <- -1
  fit$covariance <- to_working %*% reported %*% t(to_working)
  conditional <- solve(solve(reported)[-5, -5])
  expect_gt(reported["sigma2", "sigma2"], 5 * conditional[3, 3])
  set.seed(2)
  draws <- expect_no_warning(RiskMap:::draw_parameter_values(fit, 20000))
  expect_equal(draws[, "nu2"] + draws[, "sigma2"], rep(log(0.004), 20000))
  expect_equal(unname(cov(draws[, -5])), unname(conditional), tolerance = 0.05)

  # A spatial variance and a nugget that only their sum determines
  reported <- diag(c(0.04, 0.01, 1e6, 1e6, 1e6, 0.09))
  reported[3, 5] <- reported[5, 3] <- -1e6
  dimnames(reported) <- list(names_flat, names_flat)
  conditional <- RiskMap:::conditional_covariance_of(
    reported, !names_flat %in% c("sigma2", "phi", "nu2")
  )
  expect_equal(unname(conditional), diag(c(0.04, 0.01, 0.09)))
  fit$estimate$nu2 <- 0
  fit$covariance <- to_working %*% reported %*% t(to_working)
  expect_warning(RiskMap:::draw_parameter_values(fit, 10),
                 "don't determine `sigma2`, `phi`, and `tau2`",
                 class = "riskmap_warning_fit")
})

test_that("a spatial variance at zero leaves the nugget its uncertainty", {
  fit <- glgpm(
    y ~ z + gp(nugget = TRUE), data = make_records()$data, family = "gaussian",
    messages = FALSE
  )
  # A spatial variance at zero, with the likelihood flat in it and in the
  # scale. The ratio of the nugget to it moves with it, so that only the
  # variance of the nugget is determined.
  fit$estimate$sigma2 <- log(1e-8)
  fit$estimate$nu2 <- log(0.4 / 1e-8)
  names_flat <- names(unlist(fit$estimate))
  reported <- diag(c(0.04, 0.01, 1e8, 1e8, 0.25, 0.09))
  dimnames(reported) <- list(names_flat, names_flat)
  to_working <- diag(6)
  dimnames(to_working) <- dimnames(reported)
  to_working["nu2", "sigma2"] <- -1
  fit$covariance <- to_working %*% reported %*% t(to_working)
  set.seed(2)
  # The spatial variance is zero for the prediction, so there is no warning.
  draws <- expect_no_warning(RiskMap:::draw_parameter_values(fit, 4000))
  expect_equal(draws[, "sigma2"], rep(log(1e-8), 4000))
  expect_equal(draws[, "phi"], rep(fit$estimate$phi, 4000))
  log_tau2 <- draws[, "nu2"] + draws[, "sigma2"]
  expect_equal(mean(log_tau2), log(0.4), tolerance = 0.05)
  expect_equal(sd(log_tau2), 0.5, tolerance = 0.05)
  expect_equal(unname(apply(draws, 2, sd)[c(1, 2, 6)]), c(0.2, 0.1, 0.3),
               tolerance = 0.05)
})
