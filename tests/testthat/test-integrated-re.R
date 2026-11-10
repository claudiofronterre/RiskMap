# Random effects within locations integrated out (integrate_noise), with the
# rules of plan_integrated_records() and the nested patterns of site_table.h.

# The binomial log-likelihood without its constant.
binomial_kernel <- function(y, m, eta) {
  y * eta - m * (pmax(eta, 0) + log1p(exp(-abs(eta))))
}

# The Gauss-Hermite rule with `nodes` nodes for the standard normal
# distribution, from the eigenvalues of the Jacobi matrix of the Hermite
# polynomials (Golub and Welsch, 1969), whose weights sum to one.
normal_rule <- function(nodes) {
  jacobi <- matrix(0, nodes, nodes)
  off <- sqrt(seq_len(nodes - 1L))
  jacobi[cbind(seq_len(nodes - 1L), 2:nodes)] <- off
  jacobi[cbind(2:nodes, seq_len(nodes - 1L))] <- off
  decomposition <- eigen(jacobi, symmetric = TRUE)
  order <- rev(seq_len(nodes))
  list(nodes = decomposition$values[order],
       weights = decomposition$vectors[1L, order]^2)
}

# log integral N(e; 0, vo) prod_u integral N(d; 0, vi) G_u(x + e + d) dd de by
# Gauss-Hermite rules with `nodes` nodes on each Gaussian, independent of the
# package's tables. Accurate to rounding for households of a few tested.
nested_by_quadrature <- function(y, m, x, vo, vi, nodes = 60L) {
  rule <- normal_rule(nodes)
  vapply(x, function(point) {
    outer <- point + sqrt(vo) * rule$nodes
    log_inner <- vapply(seq_along(y), function(u) {
      values <- outer(outer, sqrt(vi) * rule$nodes, `+`)
      kernel <- binomial_kernel(y[u], m[u], values)
      top <- apply(kernel, 1L, max)
      top + log(drop(exp(kernel - top) %*% rule$weights))
    }, numeric(nodes))
    total <- rowSums(log_inner)
    top <- max(total)
    top + log(sum(rule$weights * exp(total - top)))
  }, numeric(1))
}

test_that("each level and location follows the rules of the plan", {
  # Location 1 has one household. Location 2 has three households alike in
  # covariates. Location 3 has two households that differ in a covariate.
  # Household 8 is recorded at locations 4 and 5.
  location <- c(1L, 2L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)
  household <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 8L, 9L)
  x <- c(0, 1, 1, 1, 0, 1, 0, 0, 0, 0)
  plan <- plan_integrated_records(
    y = c(1, 0, 2, 1, 0, 1, 3, 1, 0, 2), units_m = rep(5, 10),
    D = cbind(1, x), cov_offset = numeric(10), ID_coords = location,
    ID_re = cbind(household = household), n_loc = 5L,
    allowed = c(nugget = TRUE, "re(household)" = TRUE)
  )
  expect_identical(plan$nugget_integrated, c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(plan$ID_coords, c(1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L))
  expect_identical(plan$kind, c(1L, 2L, 1L, 1L, 1L, 0L, 0L, 1L))
  expect_identical(plan$classes, matrix(c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
                                        3L, dimnames = list(NULL, c(
                                          "nugget", "re(household)"))))
  expect_identical(plan$a, c(1L, 3L, 2L, 2L, 2L, 1L, 1L, 2L))
  expect_identical(plan$b[2L], 2L)
  # The nested record sums its households and keeps their pairs.
  expect_identical(plan$y[2L], 3)
  expect_identical(plan$patterns[[2L]], list(y = c(0, 2, 1), m = c(5, 5, 5)))
  expect_identical(unname(plan$ID_re[2L, 1L]), NA_integer_)
  # Household 8 alone stays sampled, as the first random-effect column.
  expect_identical(plan$sampled_levels, 1L)
  expect_identical(plan$latent_columns[, 2L], c(0L, 0L, 0L, 0L, 0L, 6L, 6L, 0L))
  expect_identical(plan$inner, list(2L))
  expect_identical(plan$counts$integrated, c(2L, 8L))
  expect_identical(plan$counts$levels, c(5L, 9L))
  expect_identical(plan$integrated_terms, c("nugget", "re(household)"))
  # Without the households the nugget is integrated only where a location has
  # one unit.
  nugget_only <- plan_integrated_records(
    y = c(1, 0, 2, 1, 0, 1, 3, 1, 0, 2), units_m = rep(5, 10),
    D = cbind(1, x), cov_offset = numeric(10), ID_coords = location,
    ID_re = cbind(household = household), n_loc = 5L,
    allowed = c(nugget = TRUE, "re(household)" = FALSE)
  )
  expect_identical(nugget_only$nugget_integrated,
                   c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(nugget_only$sampled_levels, 9L)
  expect_identical(nugget_only$integrated_terms, "nugget")
})

test_that("an effect with one level per location is an outer term", {
  # Villages with one location each, two households in each, and no nugget.
  plan <- plan_integrated_records(
    y = c(1, 0, 2, 1), units_m = rep(5, 4), D = matrix(1, 4L),
    cov_offset = numeric(4), ID_coords = c(1L, 1L, 2L, 2L),
    ID_re = cbind(village = c(1L, 1L, 2L, 2L), household = 1:4), n_loc = 2L,
    allowed = c(nugget = FALSE, "re(village)" = TRUE, "re(household)" = TRUE)
  )
  expect_identical(plan$kind, c(2L, 2L))
  expect_identical(plan$classes[plan$a[1L], ],
                   c(nugget = FALSE, "re(village)" = TRUE,
                     "re(household)" = FALSE))
  expect_identical(plan$classes[plan$b[1L], ],
                   c(nugget = FALSE, "re(village)" = FALSE,
                     "re(household)" = TRUE))
  expect_identical(plan$sampled_levels, c(0L, 0L))
  expect_identical(plan$counts$integrated, c(0L, 2L, 4L))
  expect_identical(plan$ID_re, cbind(village = 1:2, household = NA_integer_))
})

test_that("a nested pattern matches numerical integration", {
  y <- c(0, 2, 1)
  m <- c(5, 4, 6)
  sites <- make_site_table(sum(y), sum(m), "binomial", kind = 2L, a = 1L,
                           b = 2L, patterns = list(list(y = y, m = m)))
  x <- c(-3.1, -0.6, 1.3)
  for (vo in c(0.1, 0.8)) for (vi in c(0.15, 0.5)) {
    values <- cpp_site_evaluate(sites$pointer, c(vo, vi), sites$units,
                                matrix(x, ncol = 1L), FALSE, 0L)$h
    expect_equal(as.numeric(values), nested_by_quadrature(y, m, x, vo, vi),
                 tolerance = 1e-7)
  }
})

test_that("a nested pattern with large variances matches numerical integration", {
  # The outer integral over the shared effect by numerical integration, with
  # the units' integrals at vi from their record table, which the tests of the
  # record table check against numerical integration of their own.
  outer_reference <- function(x, records, vo, vi, n) {
    log_integrand <- function(eta) {
      inner <- cpp_site_evaluate(records$pointer, vi, records$units,
                                 matrix(rep(eta, each = n), n), TRUE, 0L)$h
      colSums(matrix(inner, n)) + stats::dnorm(eta, x, sqrt(vo), log = TRUE)
    }
    breaks <- sort(unique(c(-100, -50, 0, 30, 80,
                            x + c(-12, -6, -3, 0, 3, 6, 12) * sqrt(vo))))
    top <- max(log_integrand(seq(min(breaks), max(breaks), length.out = 2001)))
    pieces <- mapply(function(lower, upper) {
      stats::integrate(function(eta) exp(log_integrand(eta) - top), lower,
                       upper, rel.tol = 1e-11, subdivisions = 2000L,
                       stop.on.error = FALSE)$value
    }, c(-Inf, breaks), c(breaks, Inf))
    log(sum(pieces)) + top
  }
  y <- c(0, 2, 1)
  m <- c(5, 4, 6)
  probit <- make_inverse_link_functions("binomial", stats::pnorm, 3L)
  for (linkf in list(NULL, probit)) {
    sites <- make_site_table(sum(y), sum(m), "binomial", linkf = linkf,
                             kind = 2L, a = 1L, b = 2L,
                             patterns = list(list(y = y, m = m)))
    records <- make_site_table(y, m, "binomial", linkf = linkf)
    # A large outer or inner variance, or both, inside the grid and beyond.
    for (pair in list(c(1e4, 0.5), c(0.8, 1e4), c(50, 100), c(1e4, 1e4))) {
      x <- c(-15.789, 120)
      values <- cpp_site_evaluate(sites$pointer, pair, sites$units,
                                  matrix(x, nrow = 1L), TRUE, 0L)$h
      expected <- vapply(x, outer_reference, numeric(1), records = records,
                         vo = pair[1], vi = pair[2], n = length(y))
      expect_equal(as.numeric(values), expected, tolerance = 1e-8)
    }
  }
  # Both variances at the upper bound of the search for a log variance, 20,
  # where the walk took more than a minute for these values with the logit
  # link, and 1e4 with the probit, where it took more than 60 s.
  skip_on_cran()
  x <- matrix(seq(-300, 300, length.out = 20), nrow = 1L)
  sites <- make_site_table(sum(y), sum(m), "binomial", kind = 2L, a = 1L,
                           b = 2L, patterns = list(list(y = y, m = m)))
  expect_lt(system.time(cpp_site_evaluate(sites$pointer, rep(exp(20), 2L),
                                          sites$units, x, TRUE, 0L))[["elapsed"]], 3)
  sites <- make_site_table(sum(y), sum(m), "binomial", linkf = probit,
                           kind = 2L, a = 1L, b = 2L,
                           patterns = list(list(y = y, m = m)))
  expect_lt(system.time(cpp_site_evaluate(sites$pointer, c(1e4, 1e4),
                                          sites$units, x, TRUE, 0L))[["elapsed"]], 15)
})

test_that("a nested pattern's derivatives are those of its interpolant", {
  sites <- make_site_table(3, 15, "binomial", kind = 2L, a = 1L, b = 2L,
                           patterns = list(list(y = c(0, 2, 1),
                                                m = c(5, 5, 5))))
  x <- seq(-3, 2, by = 0.137)
  at <- function(x) {
    cpp_site_evaluate(sites$pointer, c(0.4, 0.6), sites$units,
                      matrix(x, ncol = 1L), FALSE, 2L)
  }
  step <- 1e-5
  centre <- at(x)
  above <- at(x + step)
  below <- at(x - step)
  expect_equal(centre$h1, (above$h - below$h) / (2 * step), tolerance = 1e-7)
  expect_equal(centre$h2, (above$h1 - below$h1) / (2 * step),
               tolerance = 1e-7)
  info <- cpp_site_table_info(sites$pointer)
  expect_identical(info$patterns, 1L)
  expect_identical(info$nested_tables_built, 1L)
})

test_that("a nested pattern with a zero variance reduces to records", {
  y <- c(0, 2, 1)
  m <- c(5, 5, 5)
  nested <- make_site_table(3, 15, "binomial", kind = 2L, a = 1L, b = 2L,
                            patterns = list(list(y = y, m = m)))
  records <- make_site_table(c(3, y), c(15, m), "binomial")
  x <- c(-2.2, 0.3)
  # One column per record and one row per point.
  evaluate <- function(sites, classes, eta) {
    cpp_site_evaluate(sites$pointer, classes, sites$units, eta, FALSE, 2L)
  }
  # Households held at zero: the summed record with the whole nugget.
  held <- evaluate(nested, c(0.7, 0), matrix(x, ncol = 1L))
  summed <- evaluate(make_site_table(3, 15, "binomial"), 0.7,
                     matrix(x, ncol = 1L))
  expect_equal(held, summed)
  # The nugget held at zero: the households' own integrals at x.
  alone <- evaluate(nested, c(0, 0.6), matrix(x, ncol = 1L))
  each <- evaluate(make_site_table(y, m, "binomial"), 0.6,
                   matrix(x, length(x), 3L))
  expect_equal(as.numeric(alone$h), rowSums(each$h))
  expect_equal(as.numeric(alone$h1), rowSums(each$h1))
  # Both held: the plain kernel of the summed record.
  plain <- evaluate(nested, c(0, 0), matrix(x, ncol = 1L))
  expect_equal(as.numeric(plain$h), binomial_kernel(3, 15, x))
})

test_that("a nested pattern's slope in its inner variance is its derivative", {
  y <- c(0, 2, 1)
  m <- c(5, 5, 5)
  sites <- make_site_table(3, 15, "binomial", kind = 2L, a = 1L, b = 2L,
                           patterns = list(list(y = y, m = m)))
  x <- matrix(c(-2.4, -0.5, 0.8), ncol = 1L)
  h_at <- function(vi) {
    as.numeric(nested_by_quadrature(y, m, as.numeric(x), 0.5, vi, 80L))
  }
  for (vi in c(0, 0.3)) {
    slopes <- cpp_site_inner_slopes(sites$pointer, c(0.5, vi), sites$units,
                                    x)
    step <- 1e-4
    numerical <- if (vi > 0) (h_at(vi + step) - h_at(vi - step)) / (2 * step) else
      (h_at(step) - h_at(0 + 1e-12)) / step
    expect_equal(as.numeric(slopes), numerical, tolerance = 2e-3)
  }
})

# A survey of households nested in locations, households alike in covariates
# within a location unless `household_covariate`.
household_survey <- function(n = 30L, households = 3L, seed = 4L,
                             household_covariate = FALSE) {
  set.seed(seed)
  coords <- cbind(runif(n), runif(n)) * 50
  distances <- as.matrix(dist(coords))
  location <- rep(seq_len(n), each = households)
  x <- if (household_covariate) rnorm(n * households) else rnorm(n)[location]
  S <- drop(t(chol(1.2 * exp(-distances / 8) + diag(1e-9, n))) %*% rnorm(n))
  eta <- -1.3 + 0.3 * x + S[location] + rnorm(n, 0, sqrt(0.5))[location] +
    rnorm(n * households, 0, sqrt(0.6))
  m <- rep(5, n * households)
  list(D = cbind(1, x), y = rbinom(length(m), m, plogis(eta)), m = m,
       location = location, distances = distances, coords = coords,
       ID_re = cbind(household = seq_along(m)),
       parameters = c(-1.3, 0.3, log(1.2), log(8), log(0.5 / 1.2), log(0.6)))
}

household_engine <- function(survey, integrate, nugget = TRUE) {
  make_adaptive_mcml_engine(survey$D, survey$y, survey$m, survey$distances,
                            family = "binomial", ID_coords = survey$location,
                            fix_tau2 = nugget, ID_re = survey$ID_re,
                            integrate_noise = integrate)
}

estimate_with <- function(engine, parameters, n_draws, seed) {
  set.seed(seed)
  draws <- matrix(rnorm(n_draws * engine$latent_dimension),
                  ncol = engine$latent_dimension)
  log_weights <- engine$evaluate(parameters, draws, log_weights_only = TRUE)
  top <- max(log_weights)
  weights <- exp(log_weights - top)
  c(estimate = top + log(mean(weights)),
    se = stats::sd(weights) / mean(weights) / sqrt(length(weights)))
}

test_that("the sampled and integrated forms estimate the same likelihood", {
  cases <- list(
    nested = list(survey = household_survey(), nugget = TRUE),
    households = list(survey = household_survey(), nugget = FALSE),
    covariate = list(survey = household_survey(household_covariate = TRUE),
                     nugget = TRUE)
  )
  for (name in names(cases)) {
    survey <- cases[[name]]$survey
    nugget <- cases[[name]]$nugget
    parameters <- if (nugget) survey$parameters else survey$parameters[-5L]
    sampled <- estimate_with(household_engine(survey, FALSE, nugget),
                             parameters, 20000L, 21L)
    engine <- household_engine(survey, TRUE, nugget)
    expect_true(engine$integrated)
    integrated <- estimate_with(engine, parameters, 20000L, 22L)
    expect_lt(abs(sampled[["estimate"]] - integrated[["estimate"]]),
              4 * sqrt(sampled[["se"]]^2 + integrated[["se"]]^2),
              label = name)
  }
})

test_that("the integrated form matches an exact likelihood of three locations", {
  # Three locations with three households each, the nugget and the households
  # integrated out, against a tensor Gauss-Hermite rule over the field and a
  # quadrature of each location's nested integral.
  survey <- household_survey(n = 3L, seed = 9L)
  engine <- household_engine(survey, TRUE)
  expect_identical(engine$latent_dimension, 3L)
  parameters <- survey$parameters
  sigma2 <- exp(parameters[3L])
  covariance <- sigma2 * exp(-survey$distances / exp(parameters[4L])) +
    diag(1e-9, 3L)
  vo <- sigma2 * exp(parameters[5L])
  vi <- exp(parameters[6L])
  mu <- drop(survey$D %*% parameters[1:2])
  # Each location's log-likelihood on a fine grid of its location effect.
  grid <- seq(-9, 9, by = 0.01)
  splines <- lapply(seq_len(3L), function(l) {
    rows <- which(survey$location == l)
    stats::splinefun(grid, nested_by_quadrature(
      survey$y[rows], survey$m[rows], mu[rows[1L]] + grid, vo, vi, 50L
    ))
  })
  rule <- normal_rule(40L)
  z <- as.matrix(expand.grid(rule$nodes, rule$nodes, rule$nodes))
  w <- Reduce(`*`, expand.grid(rule$weights, rule$weights, rule$weights))
  S <- z %*% t(t(chol(covariance)))
  total <- splines[[1L]](S[, 1L]) + splines[[2L]](S[, 2L]) +
    splines[[3L]](S[, 3L])
  top <- max(total)
  exact <- top + log(sum(w * exp(total - top)))
  estimate <- estimate_with(engine, parameters, 40000L, 23L)
  expect_lt(abs(estimate[["estimate"]] - exact), 4 * estimate[["se"]] + 1e-4)
})

test_that("a household variance held at zero gives the model without it", {
  survey <- household_survey()
  integrated <- household_engine(survey, TRUE)
  without <- make_adaptive_mcml_engine(
    survey$D, survey$y, survey$m, survey$distances, family = "binomial",
    ID_coords = survey$location, fix_tau2 = TRUE, integrate_noise = TRUE
  )
  expect_identical(integrated$hold_indices, c(5L, 6L))
  held <- replace(survey$parameters, 6L, -20)
  set.seed(24)
  draws <- matrix(rnorm(500 * 30), ncol = 30)
  expect_equal(
    integrated$evaluate(held, draws, log_weights_only = TRUE),
    without$evaluate(survey$parameters[1:5], draws, log_weights_only = TRUE)
  )
  expect_identical(
    hold_noise_at_zero(integrated, survey$parameters, boundary = 6L), held
  )
  expect_identical(
    hold_noise_at_zero(integrated, survey$parameters, boundary = 3L),
    survey$parameters
  )
})

test_that("both forms estimate the slope of a household variance at zero", {
  survey <- household_survey()
  held <- replace(survey$parameters, 6L, -20)
  slope <- function(engine, seed) {
    engine$set_anchor(held)
    proposal <- engine$proposal_at(held)
    set.seed(seed)
    draws <- matrix(rnorm(20000 * engine$latent_dimension),
                    ncol = engine$latent_dimension)
    defensive <- defensive_rows(nrow(draws), engine$defensive_weight)
    log_weights <- engine$evaluate(held, draws, log_weights_only = TRUE)
    terms <- engine$variance_slope_terms(draws, proposal, defensive, 6L)
    weights <- exp(log_weights - max(log_weights))
    weights <- weights / sum(weights)
    estimate <- sum(weights * terms)
    c(estimate = estimate,
      se = sqrt(sum(weights^2 * (terms - estimate)^2)))
  }
  sampled <- slope(household_engine(survey, FALSE), 25L)
  integrated <- slope(household_engine(survey, TRUE), 26L)
  expect_true(is.finite(integrated[["estimate"]]))
  expect_lt(abs(sampled[["estimate"]] - integrated[["estimate"]]),
            4 * sqrt(sampled[["se"]]^2 + integrated[["se"]]^2))
})

test_that("the levels integrated out are described term by term", {
  counts <- data.frame(term = c("nugget", "re(household)", "re(district)"),
                       integrated = c(28L, 1195L, 0L),
                       levels = c(30L, 1200L, 10L))
  expect_identical(
    describe_integration_counts(counts),
    paste0(" The nugget was integrated out at 28 of 30 locations and sampled",
           " at the other 2. The random effect re(household) was integrated",
           " out at 1,195 of 1,200 levels and sampled at the other 5. The",
           " random effect re(district) stayed sampled.")
  )
  expect_null(describe_integration_counts(NULL))
  counts$integrated <- 0L
  expect_null(describe_integration_counts(counts))
})

test_that("sampled households with a large variance suggest integration", {
  fit <- list(fix_tau2 = FALSE, ID_coords = c(1L, 1L, 2L, 2L),
              ID_re = cbind(household = 1:4),
              estimate = list(sigma2 = 0, sigma2_re = list(household = log(0.5))))
  sampled <- list(mcml_integrate_noise = FALSE)
  expect_true(suggests_integration(fit, sampled))
  fit$estimate$sigma2_re$household <- log(0.05)
  expect_false(suggests_integration(fit, sampled))
  # Districts spanning both locations do not.
  fit$ID_re <- cbind(district = c(1L, 2L, 1L, 2L))
  fit$estimate$sigma2_re <- list(district = log(0.5))
  expect_false(suggests_integration(fit, sampled))
})

test_that("a fit integrates households out and reports them", {
  skip_on_cran()
  survey <- household_survey()
  data <- sf::st_as_sf(data.frame(
    y = survey$y, m = survey$m, x = survey$D[, 2L],
    household = survey$ID_re[, 1L],
    east = 1000 * survey$coords[survey$location, 1L],
    north = 1000 * survey$coords[survey$location, 2L]
  ), coords = c("east", "north"), crs = 32636)
  fit <- suppressWarnings(glgpm(
    y ~ x + gp(nugget = TRUE) + re(household), data = data,
    family = "binomial", denominator = m, distance_units = "km",
    control_mcmc = control_mcmc,
    control_mcml = set_control_mcml(integrate_noise = TRUE, max_stages = 1),
    seed = fixture_seed, messages = FALSE
  ))
  optimizer <- attr(fit, "optimizer")
  expect_true(optimizer$mcml_integrate_noise)
  expect_identical(optimizer$mcml_integrated_terms,
                   c("nugget", "re(household)"))
  expect_match(check_fit(fit)$draws,
               "The random effect re\\(household\\) was integrated out at 90 of 90 levels.")
  expect_identical(summarise_mcml_fit(fit)$integration$integrated,
                   c(30L, 90L))
})

test_that("an automatic fit with households is the fit with TRUE", {
  skip_on_cran()
  survey <- household_survey()
  data <- sf::st_as_sf(data.frame(
    y = survey$y, m = survey$m, x = survey$D[, 2L],
    household = survey$ID_re[, 1L],
    east = 1000 * survey$coords[survey$location, 1L],
    north = 1000 * survey$coords[survey$location, 2L]
  ), coords = c("east", "north"), crs = 32636)
  fit_with <- function(integrate_noise) {
    control <- set_control_mcml(integrate_noise = integrate_noise,
                                max_stages = 1)
    # A threshold no variance reaches would sample a nugget that switches.
    control$integrate_threshold <- 1e6
    suppressWarnings(glgpm(
      y ~ x + gp(nugget = TRUE) + re(household), data = data,
      family = "binomial", denominator = m, distance_units = "km",
      control_mcmc = control_mcmc, control_mcml = control,
      seed = fixture_seed, messages = FALSE
    ))
  }
  automatic <- fit_with("auto")
  integrated <- fit_with(TRUE)
  expect_identical(coef(automatic), coef(integrated))
  expect_identical(automatic$log_lik, integrated$log_lik)
  expect_identical(attr(automatic, "optimizer")$mcml_certificate,
                   attr(integrated, "optimizer")$mcml_certificate)
  expect_identical(summarise_mcml_fit(automatic)$integration$integrated,
                   c(30L, 90L))
})
