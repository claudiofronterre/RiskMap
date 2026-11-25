# Columns of the user's data ---------------------------------------------------

# The data and a grid of new locations with one covariate under several names.
named_data <- gaussian_data[c("cov", "denominator")]
named_data$y <- named_data$units_m <- named_data$cov
named_grid <- st_as_sf(
  data.frame(x = c(1000, 3000, 5000, 7000, 9000), z = c(2000, 8000, 5000, 1000, 6000),
             cov = c(-1, -0.5, 0, 0.5, 1)),
  coords = c("x", "z"), crs = 32637
)
named_grid$y <- named_grid$units_m <- named_grid$cov
named_areas <- st_sf(name = c("west", "east"),
                     geometry = st_make_grid(named_grid, n = c(2, 1)))

assess_named <- function(formula, family, parameters, ...) {
  model <- if (family == "gaussian") {
    specify_glgpm(formula, named_data, family, parameters)
  } else {
    specify_glgpm(formula, named_data, family, parameters,
                  denominator = denominator)
  }
  sim <- simulate_glgpm(model, nsim = 2, what = c("data", "surface"),
                        prediction_grid = named_grid, seed = 3)
  suppressWarnings(assess_simulation(
    sim, models = list(model = model), spatial_scale = c("location", "area"),
    aggregate = "mean", boundaries = named_areas, label = "name",
    objective = "mse", control_mcmc = control_mcmc, seed = fixture_seed,
    messages = FALSE, ...
  ))
}

test_that("assess_simulation keeps covariates called y and units_m", {
  expect_same_assessment <- function(x, y) {
    expect_identical(x$fits$status, c("ok", "ok"))
    expect_identical(x$objective, y$objective)
    expect_identical(unname(x$estimates$model), unname(y$estimates$model))
  }
  gaussian <- list(beta = c(1, 0.5), sigma2 = 1, phi = 2, sigma2_me = 0.1)
  expect_same_assessment(assess_named(cases ~ y + gp(), "gaussian", gaussian),
                         assess_named(cases ~ cov + gp(), "gaussian", gaussian))
  expect_same_assessment(
    assess_named(cases ~ I(y^2) + gp(), "gaussian", gaussian),
    assess_named(cases ~ I(cov^2) + gp(), "gaussian", gaussian)
  )
  binomial <- list(beta = c(-0.5, 0.5), sigma2 = 0.5, phi = 2)
  reference <- assess_named(cases ~ cov + gp(), "binomial", binomial)
  expect_same_assessment(assess_named(cases ~ units_m + gp(), "binomial", binomial),
                         reference)
  expect_same_assessment(assess_named(cases ~ y + gp(), "binomial", binomial),
                         reference)
})

test_that("assess_simulation asks for a covariate called y on the grid", {
  gaussian <- list(beta = c(1, 0.5), sigma2 = 1, phi = 2, sigma2_me = 0.1)
  generating <- specify_glgpm(cases ~ cov + gp(), named_data, "gaussian", gaussian)
  candidate <- specify_glgpm(cases ~ y + gp(), named_data, "gaussian", gaussian)
  sim <- simulate_glgpm(generating, nsim = 1, what = c("data", "surface"),
                        prediction_grid = named_grid["cov"], seed = 3)
  expect_input_error(
    assess_simulation(sim, models = list(candidate = candidate),
                      spatial_scale = "location", objective = "mse",
                      messages = FALSE),
    "models", "has no column \"y\""
  )
})

test_that("a covariate called as the response of the generating model is kept", {
  # The response of a specification need not exist in the data, so the data
  # can have a covariate of that name, which another model uses.
  gaussian <- list(beta = c(1, 0.5), sigma2 = 1, phi = 2, sigma2_me = 0.1)
  generating <- specify_glgpm(y ~ cov + gp(), named_data, "gaussian", gaussian)
  sim <- simulate_glgpm(generating, nsim = 2, what = c("data", "surface"),
                        prediction_grid = named_grid, seed = 3)
  result <- suppressWarnings(assess_simulation(
    sim,
    models = list(
      named_y = specify_glgpm(cases ~ y + gp(), named_data, "gaussian", gaussian),
      named_cov = specify_glgpm(cases ~ cov + gp(), named_data, "gaussian", gaussian)
    ),
    spatial_scale = "location", objective = "mse", control_mcmc = control_mcmc,
    seed = fixture_seed, messages = FALSE
  ))
  mse <- result$objective$location$mse
  expect_identical(mse["named_y", ], mse["named_cov", ])
})

test_that("assess_simulation refuses a candidate with a transformed response", {
  parameters <- list(beta = 1, sigma2 = 1, phi = 2, sigma2_me = 0.1)
  generating <- specify_glgpm(y ~ gp(), gaussian_data, "gaussian", parameters)
  candidate <- specify_glgpm(log(y) ~ gp(), gaussian_data, "gaussian", parameters)
  sim <- simulate_glgpm(generating, nsim = 1, what = c("data", "surface"),
                        prediction_grid = named_grid, seed = 3)
  expect_input_error(
    assess_simulation(sim, models = list(plain = generating, logged = candidate),
                      spatial_scale = "location", objective = "mse",
                      messages = FALSE),
    "models", "log(y)"
  )
  # A candidate may call its response differently.
  other_name <- specify_glgpm(cases ~ gp(), gaussian_data, "gaussian", parameters)
  result <- suppressWarnings(assess_simulation(
    sim, models = list(y = generating, cases = other_name),
    spatial_scale = "location", objective = "mse", seed = fixture_seed,
    messages = FALSE
  ))
  mse <- result$objective$location$mse
  expect_identical(mse["y", ], mse["cases", ])
})

# Classes ----------------------------------------------------------------------

test_that("every value belongs to a class, as in plot_prediction()", {
  breaks <- c(0, 0.2, 0.4)
  values <- c(-1, 0, 0.1, 0.2, 0.21, 0.4, 5)
  classes <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L)
  expect_identical(simulation_class(values, breaks), classes)

  scale <- prediction_fill_scale(colours = grDevices::hcl.colors(256L),
                                 limits = NULL, values = values,
                                 breaks = breaks, role = "level",
                                 probability = FALSE)
  colours <- scale$map(values)
  expect_identical(match(colours, unique(colours)), classes)

  # A matrix of draws keeps its shape.
  draws <- matrix(values[1:6], 2)
  expect_identical(simulation_class(draws, breaks), matrix(classes[1:6], 2))
})

test_that("values outside the breaks and at the first break count in the classification", {
  # Two values in the low class, one at the first break, and three in the high
  # class, two above the last break, with draws that equal the truth.
  truth <- c(0, 0.1, 0.3, 0.9, 5)
  metrics <- simulation_classification_metrics(
    truth, matrix(truth, 5, 50), c(0, 0.2, 0.4), c("low", "high")
  )
  expect_identical(metrics$overall_cc, 1)
  expect_identical(metrics$by_cat$Sensitivity, c(1, 1))

  always_low <- simulation_classification_metrics(
    truth, matrix(0.1, 5, 50), c(0, 0.2, 0.4), c("low", "high")
  )
  expect_identical(always_low$overall_cc, 2 / 5)
  expect_identical(always_low$by_cat$PPV, c(2 / 5, NA))
  expect_identical(always_low$by_cat$CC, c(2 / 5, 2 / 5))
})

test_that("draws beyond the breaks are classed in the end classes", {
  breaks <- c(-10, 0, 10)
  above <- simulation_classification_metrics(c(-1, 1, 5), matrix(30, 3, 50),
                                             breaks, c("low", "high"))
  expect_identical(above$overall_cc, 2 / 3)
  expect_identical(above$by_cat$Sensitivity, c(0, 1))
  # 6 draws far above the last break and 4 in the low class.
  mixed <- simulation_classification_metrics(
    1, matrix(rep(c(30, -5), c(6, 4)), 1), breaks, c("low", "high")
  )
  expect_identical(mixed$overall_cc, 1)
  # A tie goes to the lowest class.
  tied <- simulation_classification_metrics(
    -1, matrix(rep(c(30, -5), c(5, 5)), 1), breaks, c("low", "high")
  )
  expect_identical(tied$overall_cc, 1)
})

# The nugget at sampled sites --------------------------------------------------

# Sites with repeated observations in groups, and a grid that mixes sampled
# sites with new locations.
nugget_sites <- local({
  set.seed(4)
  coordinates <- cbind(runif(8, 0, 10000), runif(8, 0, 10000))
  site <- c(1:8, 1:8, 1:4)
  st_as_sf(data.frame(x = coordinates[site, 1], z = coordinates[site, 2],
                      site = site, group = rep(1:4, 5),
                      cov = coordinates[site, 1] / 10000,
                      m = rep(c(20, 30), 10), y = rnorm(20)),
           coords = c("x", "z"), crs = 32637)
})
nugget_grid <- local({
  new <- st_as_sf(data.frame(x = c(1500, 4500, 6500, 8500, 2500),
                             z = c(1500, 8500, 2500, 5500, 5500),
                             site = NA, cov = c(0.15, 0.45, 0.65, 0.85, 0.25)),
                  coords = c("x", "z"), crs = 32637)
  rbind(nugget_sites[c(2, 7), c("site", "cov")], new[1:2, ],
        nugget_sites[5, c("site", "cov")], new[3:5, ])
})
# Two blocks with sampled sites and one without.
nugget_groups <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L)

# The mean and covariance of the linear predictor at `grid` given the data of
# a Gaussian fit to `nugget_sites`, with the nugget of a sampled site the one
# of its observations.
exact_gaussian_prediction <- function(fit, grid) {
  pars <- coef(fit)
  correlation <- function(a, b) {
    matern_correlation(cross_distances(a, b), phi = pars$phi, kappa = fit$kappa)
  }
  sites <- st_coordinates(nugget_sites) / 1000
  points <- st_coordinates(grid) / 1000
  same_site <- outer(grid$site, nugget_sites$site, `==`)
  same_site[is.na(same_site)] <- FALSE
  observed <- pars$sigma2 * correlation(sites, sites) +
    pars$tau2 * outer(nugget_sites$site, nugget_sites$site, `==`) +
    diag(pars$sigma2_me, nrow(nugget_sites))
  if (length(pars$sigma2_re)) {
    observed <- observed + pars$sigma2_re[[1]] *
      outer(nugget_sites$group, nugget_sites$group, `==`)
  }
  cross <- pars$sigma2 * correlation(points, sites) + pars$tau2 * same_site
  weights <- cross %*% solve(observed)
  prior <- pars$sigma2 * correlation(points, points) +
    diag(pars$tau2, nrow(grid))
  list(
    mean = as.numeric(cbind(1, grid$cov) %*% pars$beta +
                        weights %*% (fit$y - fit$D %*% pars$beta)),
    covariance = prior - weights %*% t(cross)
  )
}

# Compares the draws with their exact distribution, in units of the Monte
# Carlo standard errors of their means and covariances. Every location is in
# one of `groups`, and draws are joint within each.
expect_exact_draws <- function(draws, exact, groups) {
  n_draws <- ncol(draws)
  variances <- diag(exact$covariance)
  mean_error <- (rowMeans(draws) - exact$mean) / sqrt(variances / n_draws)
  expect_lt(max(abs(mean_error)), 4)
  standard_error <- sqrt((outer(variances, variances) + exact$covariance^2) /
                           n_draws)
  covariance_error <- (cov(t(draws)) - exact$covariance) / standard_error
  expect_lt(max(abs(covariance_error[outer(groups, groups, `==`)])), 4.5)
}

nugget_gaussian_fit <- function(formula, ...) {
  suppressWarnings(glgpm(
    formula, data = nugget_sites, family = "gaussian", messages = FALSE,
    fixed_pars = list(beta = c(0.5, 1), sigma2 = 1, phi = 2, tau2 = 1.5,
                      sigma2_me = 0.2, ...)
  ))
}

test_that("draws of a Gaussian model at sampled sites condition on their nugget", {
  fit <- nugget_gaussian_fit(y ~ cov + gp(nugget = TRUE))
  exact <- exact_gaussian_prediction(fit, nugget_grid)
  control <- set_control_mcmc(n_draws = 5000, burnin = 0, thin = 1)
  set.seed(11)
  joint <- assessment_response_draws(fit, nugget_grid, "nugget", control,
                                     nugget_groups)
  expect_exact_draws(joint, exact, nugget_groups)
  set.seed(12)
  marginal <- assessment_response_draws(fit, nugget_grid, "nugget", control)
  expect_exact_draws(marginal, exact, seq_len(nrow(nugget_grid)))

  # The groups of a candidate with a random effect are part of the data.
  grouped <- nugget_gaussian_fit(y ~ cov + gp(nugget = TRUE) + re(group),
                                 sigma2_re = 0.7)
  set.seed(13)
  draws <- assessment_response_draws(grouped, nugget_grid, "nugget", control,
                                     nugget_groups)
  expect_exact_draws(draws, exact_gaussian_prediction(grouped, nugget_grid),
                     nugget_groups)
})

test_that("draws of a Gaussian model without measurement error are its data", {
  sites <- nugget_sites[1:8, ]
  fit <- suppressWarnings(glgpm(
    y ~ gp(nugget = TRUE), data = sites, family = "gaussian", messages = FALSE,
    fixed_pars = list(beta = 0.5, sigma2 = 1, phi = 2, tau2 = 1.5,
                      sigma2_me = 0)
  ))
  control <- set_control_mcmc(n_draws = 20, burnin = 0, thin = 1)
  for (groups in list(NULL, nugget_groups)) {
    expect_no_warning(
      draws <- assessment_response_draws(fit, nugget_grid, "nugget", control,
                                         groups)
    )
    expect_equal(draws[c(1, 2, 5), ], matrix(sites$y[c(2, 7, 5)], 3, 20),
                 tolerance = 1e-5)
    expect_gt(min(apply(draws[-c(1, 2, 5), ], 1, sd)), 0.5)
  }

  # With groups, the data leave only their variance.
  grouped <- suppressWarnings(glgpm(
    y ~ gp(nugget = TRUE) + re(group), data = sites, family = "gaussian",
    messages = FALSE,
    fixed_pars = list(beta = 0.5, sigma2 = 1, phi = 2, tau2 = 1.5,
                      sigma2_me = 0, sigma2_re = 0.7)
  ))
  expect_no_warning(
    draws <- assessment_response_draws(grouped, nugget_grid, "nugget", control)
  )
  # Sites 2 and 5 are in different groups, and 1 and 5 in the same.
  group_effect <- sites$y[c(2, 7, 5)] - draws[c(1, 2, 5), ]
  expect_gt(min(apply(group_effect, 1, sd)), 0.1)
  same_group <- rbind(nugget_sites[c(1, 5), c("site", "cov")])
  expect_no_warning(
    draws <- assessment_response_draws(grouped, same_group, "nugget", control)
  )
  expect_equal(sites$y[1] - draws[1, ], sites$y[5] - draws[2, ],
               tolerance = 1e-4)
})

test_that("joint draws share the nugget of a repeated location", {
  # Its nugget is that of one point in the simulated truth too.
  control <- set_control_mcmc(n_draws = 50, burnin = 100, thin = 2)
  nugget_sites$cases <- rep(c(3, 12, 7, 9), 5)
  fits <- list(
    nugget_gaussian_fit(y ~ cov + gp(nugget = TRUE)),
    suppressWarnings(glgpm(
      cases ~ cov + gp(nugget = TRUE), data = nugget_sites, family = "binomial",
      denominator = m, messages = FALSE,
      fixed_pars = list(beta = c(0.2, 0.1), sigma2 = 0.3, phi = 2, tau2 = 1.5)
    ))
  )
  for (fit in fits) {
    for (rows in list(c(3, 3, 4), c(1, 1, 3, 3, 4))) {
      grid <- nugget_grid[rows, ]
      # Only the fixed part differs between the rows of a repeated location.
      grid$cov[which(rows == 3)] <- c(0.15, 0.95)
      set.seed(17)
      # The prediction engine warns of the jitter that equal rows need.
      draws <- suppressWarnings(
        assessment_response_draws(fit, grid, "nugget", control,
                                  rep(1L, length(rows))),
        classes = "riskmap_warning_fit"
      )
      link <- if (fit$family == "gaussian") draws else qlogis(draws)
      difference <- link[which(rows == 3), ]
      slope <- coef(fit)$beta[[2]]
      expect_equal(difference[2, ] - difference[1, ], rep(0.8 * slope, 50),
                   tolerance = 1e-4)
      expect_gt(sd(link[1, ] - link[length(rows), ]), 1)
    }
  }
})

test_that("draws of a binomial model at sampled sites are those of the fit's latent effects", {
  # Sites so far apart that their effects are independent, which makes the
  # exact distribution at each an integral over one effect.
  far <- st_as_sf(data.frame(x = c(0, 1e6, 2e6, 3e6), z = 0,
                             cases = c(1, 19, 12, 4), m = 20),
                  coords = c("x", "z"), crs = 32637)
  fit <- suppressWarnings(glgpm(
    cases ~ gp(nugget = TRUE), data = far, family = "binomial",
    denominator = m, messages = FALSE,
    fixed_pars = list(beta = 0.2, sigma2 = 0.3, phi = 2, tau2 = 1.5)
  ))
  new <- st_as_sf(data.frame(x = 5e5, z = 0), coords = c("x", "z"), crs = 32637)
  grid <- rbind(far[c(2, 1), "geometry"], new)
  exact <- vapply(c(2, 1), function(i) {
    weight <- function(w) dbinom(far$cases[i], 20, plogis(0.2 + w)) *
      dnorm(w, 0, sqrt(1.8))
    integrate(function(w) plogis(0.2 + w) * weight(w), -10, 10)$value /
      integrate(weight, -10, 10)$value
  }, numeric(1))
  exact <- c(exact, integrate(function(w) plogis(0.2 + w) * dnorm(w, 0, sqrt(1.8)),
                              -10, 10)$value)
  control <- set_control_mcmc(n_draws = 4000, burnin = 1000, thin = 5)
  for (groups in list(NULL, c(1L, 1L, 2L))) {
    set.seed(14)
    draws <- assessment_response_draws(fit, grid, "nugget", control, groups)
    expect_equal(rowMeans(draws), exact, tolerance = 0.02)
  }
})

test_that("draws away from the sampled sites are those of the prediction engine", {
  control <- set_control_mcmc(n_draws = 50, burnin = 100, thin = 2)
  gaussian <- nugget_gaussian_fit(y ~ cov + gp(nugget = TRUE))
  nugget_sites$cases <- rep(c(3, 12, 7, 9), 5)
  binomial <- suppressWarnings(glgpm(
    cases ~ cov + gp(nugget = TRUE) + re(group), data = nugget_sites,
    family = "binomial", denominator = m, messages = FALSE,
    control_mcml = set_control_mcml(integrate_noise = FALSE),
    fixed_pars = list(beta = c(0.2, 0.1), sigma2 = 0.3, phi = 2, tau2 = 1.5,
                      sigma2_re = 0.2)
  ))
  sampled <- !is.na(nugget_grid$site)
  for (fit in list(gaussian, binomial)) {
    for (groups in list(NULL, nugget_groups)) {
      set.seed(15)
      plain <- prediction_sample_components(
        fit, nugget_grid, if (is.null(groups)) "marginal" else "joint",
        "nugget", control, messages = FALSE, groups = groups
      )$response
      set.seed(15)
      draws <- assessment_response_draws(fit, nugget_grid, "nugget", control,
                                         groups)
      same <- if (is.null(groups)) !sampled else groups == 3L
      expect_identical(draws[same, ], plain[same, ])
      expect_false(any(draws[sampled, ] == plain[sampled, ]))
      # Nothing changes without a sampled site or without a nugget to predict.
      for (case in list(list(include = "nugget", rows = !sampled),
                        list(include = character(), rows = rep(TRUE, 8)))) {
        set.seed(16)
        plain <- prediction_sample_components(
          fit, nugget_grid[case$rows, ],
          if (is.null(groups)) "marginal" else "joint", case$include, control,
          messages = FALSE, groups = groups[case$rows]
        )$response
        set.seed(16)
        draws <- assessment_response_draws(fit, nugget_grid[case$rows, ],
                                           case$include, control,
                                           groups[case$rows])
        expect_identical(draws, plain)
      }
    }
  }
})

test_that("assess_simulation predicts the realised nugget at sampled sites", {
  set.seed(1)
  n <- 60
  sites <- st_as_sf(data.frame(x = runif(n, 0, 1e5), z = runif(n, 0, 1e5)),
                    coords = c("x", "z"), crs = 32633)
  parameters <- list(beta = 0, sigma2 = 1, phi = 20, tau2 = 1, sigma2_me = 0.05)
  model <- specify_glgpm(y ~ gp(nugget = TRUE), sites, "gaussian", parameters)
  candidate <- suppressWarnings(glgpm(
    y ~ gp(nugget = TRUE), data = transform(sites, y = rnorm(n)),
    family = "gaussian", messages = FALSE, fixed_pars = parameters
  ))
  sim <- simulate_glgpm(model, nsim = 10, what = c("data", "surface"),
                        prediction_grid = sites, seed = 2)
  result <- suppressWarnings(assess_simulation(
    sim, models = list(model = candidate), spatial_scale = "location",
    objective = "mse", seed = 1, messages = FALSE
  ))
  # The mean squared errors of the exact conditional means of the surface with
  # its nugget, 0.048 on average.
  correlation <- exp(-as.matrix(dist(st_coordinates(sites) / 1000)) / 20)
  covariance <- correlation + diag(n)
  weights <- covariance %*% solve(covariance + diag(0.05, n))
  exact <- colMeans((weights %*% sim$samples$data[, , "response"] -
                       sim$samples$surface[, , "mean"])^2)
  expect_equal(as.numeric(result$objective$location$mse), exact,
               tolerance = 0.03)
})

# An offset that differs between the data and the grid -------------------------

test_that("assess_simulation adds an offset that is zero at the samples", {
  assess_offset <- function(value, formula, family, parameters) {
    nugget_sites$exposure <- 0
    nugget_sites$count <- rep(c(3, 12, 7, 9), 5)
    nugget_grid$exposure <- value
    model <- specify_glgpm(formula, nugget_sites, family, parameters)
    candidate <- suppressWarnings(glgpm(
      formula, data = nugget_sites, family = family, messages = FALSE,
      fixed_pars = parameters
    ))
    sim <- simulate_glgpm(model, nsim = 2, what = c("data", "surface"),
                          prediction_grid = nugget_grid, seed = 3)
    result <- suppressWarnings(assess_simulation(
      sim, models = list(model = candidate),
      spatial_scale = c("location", "area"), aggregate = "mean",
      boundaries = named_areas, label = "name", objective = "mse",
      control_mcmc = control_mcmc, seed = fixture_seed, messages = FALSE
    ))
    list(data = sim$samples$data[, , "response"],
         surface = sim$samples$surface[, , "mean"],
         mse = lapply(result$objective, function(scale) unname(scale$mse)))
  }
  parameters <- list(beta = c(0.5, 1), sigma2 = 1, phi = 2)
  gaussian <- c(parameters, sigma2_me = 0.2)
  cases <- list(
    list(count ~ cov + gp() + offset(exposure), "gaussian", gaussian),
    # Three of the grid locations are sampled sites, whose draws use their nugget.
    list(count ~ cov + gp(nugget = TRUE) + offset(exposure), "gaussian",
         c(gaussian, tau2 = 1.5)),
    list(count ~ cov + gp(nugget = TRUE) + offset(exposure), "poisson",
         c(parameters, tau2 = 1.5))
  )
  for (case in cases) {
    zero <- do.call(assess_offset, c(0, case))
    one <- do.call(assess_offset, c(1, case))
    expect_identical(names(zero$mse), c("location", "area"))
    expect_identical(one$data, zero$data)
    # An offset of 1 adds 1 to the surface and its predictions with an identity
    # link, and multiplies them by exp(1) with a log link.
    if (case[[2]] == "gaussian") {
      expect_equal(one$surface, zero$surface + 1)
      expect_equal(one$mse, zero$mse)
    } else {
      expect_equal(one$surface, exp(1) * zero$surface)
      expect_equal(one$mse, lapply(zero$mse, `*`, exp(2)))
    }
  }
})

# Edge cases -------------------------------------------------------------------

edge_model <- specify_glgpm(y ~ gp(), gaussian_data, "gaussian",
                            list(beta = 1, sigma2 = 1, phi = 2, sigma2_me = 0.1))
edge_simulation <- simulate_glgpm(edge_model, nsim = 2,
                                  what = c("data", "surface"),
                                  prediction_grid = named_grid, seed = 2)

test_that("a model can't be called Class when models are classified", {
  models <- list(Class = edge_model)
  expect_input_error(
    assess_simulation(edge_simulation, models, spatial_scale = "location",
                      objective = "classify", breaks = c(-10, 0, 10),
                      messages = FALSE),
    "models", "Class"
  )
  result <- suppressWarnings(assess_simulation(
    edge_simulation, models, spatial_scale = "location", objective = "mse",
    seed = 1, messages = FALSE
  ))
  expect_true(all(is.finite(result$objective$location$mse)))
})

test_that("assess_simulation stops when no area has grid locations", {
  far <- st_sf(name = "far",
               geometry = st_set_crs(st_geometry(named_areas)[1] + c(1e6, 1e6),
                                     st_crs(named_areas)))
  expect_input_error(
    assess_simulation(edge_simulation, list(model = edge_model),
                      spatial_scale = c("location", "area"), aggregate = "mean",
                      boundaries = far, label = "name", objective = "mse",
                      messages = FALSE),
    "boundaries"
  )
})

test_that("summary() and print() of an assessment report fits with warnings", {
  real_glgpm <- glgpm
  local_mocked_bindings(
    glgpm = function(..., seed) {
      fit <- real_glgpm(..., seed = seed)
      warning("no maximum")
      fit
    },
    .package = "RiskMap"
  )
  expect_warning(
    result <- assess_simulation(
      edge_simulation, list(model = edge_model), spatial_scale = "location",
      objective = "mse", seed = 1, messages = FALSE
    ),
    "2 of 2 model fits gave warnings"
  )
  expect_identical(result$fits$status, c("ok", "ok"))
  expect_identical(summary(result)$fits, c(n = 2L, failed = 0L, warned = 2L))
  expect_match(squish(capture.output(print(summary(result)))),
               "2 of 2 model fits gave warnings", fixed = TRUE)
  expect_match(squish(capture.output(print(result))),
               "2 of 2 model fits gave warnings", fixed = TRUE)
})

test_that("a surface on points that are not a regular grid is mapped as points", {
  # A regular grid in longitude and latitude isn't one in the model's CRS.
  corners <- st_bbox(latlon_data)
  latlon_grid <- st_as_sf(
    expand.grid(x = seq(corners[["xmin"]], corners[["xmax"]], length.out = 6),
                z = seq(corners[["ymin"]], corners[["ymax"]], length.out = 5)),
    coords = c("x", "z"), crs = 4326
  )
  sim <- suppressMessages(simulate_glgpm(
    edge_model, nsim = 1, what = c("data", "surface"),
    prediction_grid = latlon_grid, seed = 2
  ))
  map <- plot_simulation(sim)
  expect_s3_class(map$layers[[1]]$geom, "GeomPoint")
  drawn <- ggplot2::layer_data(map, 1)
  expect_identical(nrow(drawn), 30L)
  expect_false(anyNA(drawn$fill))
  expect_gt(length(unique(drawn$fill)), 1L)
  expect_length(map$layers, 2)

  # A regular grid is drawn with tiles, whatever its gaps and spacings.
  regular <- st_as_sf(expand.grid(x = seq(0, 9000, 1000), z = seq(0, 9000, 3000)),
                      coords = c("x", "z"), crs = 32637)
  for (rows in list(1:40, seq(1, 40, by = 2), 3)) {
    surface <- simulate_glgpm(edge_model, what = "surface",
                              prediction_grid = regular[rows, ], seed = 2)
    expect_s3_class(plot_simulation(surface)$layers[[1]]$geom, "GeomTile")
  }
})
