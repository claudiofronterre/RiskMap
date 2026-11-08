# The likelihood of a record with its nugget integrated out, by adaptive
# integration about the mode of the integrand, as a reference for the table.
integrated_by_integrate <- function(y, m, x, v, poisson) {
  log_integrand <- function(e) {
    eta <- x + e
    kernel <- if (poisson) y * eta - m * exp(eta) else
      y * eta - m * (pmax(eta, 0) + log1p(exp(-abs(eta))))
    kernel + stats::dnorm(e, 0, sqrt(v), log = TRUE)
  }
  reach <- 40 * sqrt(v) + 80
  mode <- stats::optimize(log_integrand, c(-reach, reach), maximum = TRUE,
                          tol = 1e-12)$maximum
  top <- log_integrand(mode)
  width <- min(sqrt(v), 1)
  log(stats::integrate(function(e) exp(log_integrand(e) - top),
                       mode - 60 * width, mode + 60 * width, rel.tol = 1e-12,
                       subdivisions = 2000L)$value) + top
}

test_that("the table of integrated sites matches numerical integration", {
  for (family in c("binomial", "poisson")) {
    poisson <- family == "poisson"
    y <- if (poisson) c(0, 3, 40) else c(0, 3, 20)
    m <- if (poisson) c(1, 1, 2) else c(20, 20, 20)
    sites <- make_site_table(y, m, family)
    # Both rules, adaptive Gauss-Hermite up to v = 0.2 and the trapezium above,
    # inside the grid and beyond it.
    for (v in c(0.01, 0.15, 0.5, 12)) {
      for (x in c(-55, -6, -0.4, 1.7)) {
        values <- cpp_site_evaluate(sites$pointer, v, sites$pair,
                                    matrix(x, 3L, 1L), TRUE, 0L)$h
        reference <- vapply(seq_along(y), function(i) {
          integrated_by_integrate(y[i], m[i], x, v, poisson)
        }, numeric(1))
        expect_equal(as.numeric(values), reference, tolerance = 1e-7)
      }
    }
  }
})

test_that("the table's derivatives are those of its interpolant", {
  sites <- make_site_table(c(2, 7), c(20, 20), "binomial")
  x <- seq(-3, 2, by = 0.137)
  at <- function(x) {
    cpp_site_evaluate(sites$pointer, 0.8, sites$pair, cbind(x, x), FALSE, 2L)
  }
  step <- 1e-5
  centre <- at(x)
  above <- at(x + step)
  below <- at(x - step)
  expect_equal(centre$h1, (above$h - below$h) / (2 * step), tolerance = 1e-7)
  expect_equal(centre$h2, (above$h1 - below$h1) / (2 * step),
               tolerance = 1e-7)
})

test_that("the table fills only the cells it needs and keeps 16 tables", {
  sites <- make_site_table(c(1, 5), c(10, 10), "binomial")
  cpp_site_evaluate(sites$pointer, 0.3, sites$pair, matrix(0.1, 2L, 1L), TRUE,
                    0L)
  info <- cpp_site_table_info(sites$pointer)
  expect_identical(info$tables_built, 1L)
  expect_identical(info$cells_filled, 4)
  for (v in seq(0.1, 2, length.out = 20)) {
    cpp_site_evaluate(sites$pointer, v, sites$pair, matrix(0.1, 2L, 1L), TRUE,
                      0L)
  }
  # The first table was dropped, so v = 0.3 is built again.
  before <- cpp_site_table_info(sites$pointer)$tables_built
  cpp_site_evaluate(sites$pointer, 0.3, sites$pair, matrix(0.1, 2L, 1L), TRUE,
                    0L)
  expect_identical(cpp_site_table_info(sites$pointer)$tables_built,
                   before + 1L)
})

test_that("a negative pair or a zero variance takes the plain kernel", {
  sites <- make_site_table(c(3, 3), c(20, 20), "binomial",
                           plain = c(FALSE, TRUE))
  expect_identical(sites$pair, c(1L, -1L))
  eta <- matrix(c(-1.2, 0.4), ncol = 1L)
  plain <- 3 * eta - 20 * log1p(exp(eta))
  values <- cpp_site_evaluate(sites$pointer, 0.5, sites$pair, eta, TRUE, 0L)$h
  expect_equal(values[2], plain[2])
  expect_false(isTRUE(all.equal(values[1], plain[1])))
  held <- cpp_site_evaluate(sites$pointer, 0, sites$pair, eta, TRUE, 0L)$h
  expect_equal(as.numeric(held), as.numeric(plain))
})

test_that("identical records at a location are summed and others flagged", {
  D <- cbind(1, c(0.5, 0.5, 1, 2))
  records <- prepare_integrated_records(
    y = c(1, 2, 3, 4), units_m = c(10, 10, 10, 10), D = D,
    cov_offset = numeric(4), ID_coords = c(1L, 1L, 2L, 2L), ID_re = NULL,
    n_loc = 3L
  )
  expect_identical(records$y, c(3, 3, 4))
  expect_identical(records$units_m, c(20, 10, 10))
  expect_identical(records$ID_coords, c(1L, 2L, 2L))
  expect_identical(records$several, c(FALSE, TRUE, FALSE))
})

# A small survey with a nugget, for both forms.
integrated_survey <- function(family, n = 40L, seed = 3L) {
  set.seed(seed)
  coords <- cbind(runif(n), runif(n)) * 50
  distances <- as.matrix(dist(coords))
  x <- rnorm(n)
  S <- drop(t(chol(1.2 * exp(-distances / 8) + diag(1e-9, n))) %*% rnorm(n))
  eta <- -1.5 + 0.3 * x + S + rnorm(n, 0, sqrt(0.7))
  units <- if (family == "binomial") rep(20, n) else rep(1, n)
  y <- if (family == "binomial") rbinom(n, 20, plogis(eta)) else
    rpois(n, exp(eta + 1))
  list(D = cbind(1, x), y = y, units = units, distances = distances,
       parameters = c(-1.5, 0.3, log(1.2), log(8), log(0.7 / 1.2)))
}

estimate_log_likelihood <- function(engine, parameters, draws) {
  log_weights <- engine$evaluate(parameters, draws, log_weights_only = TRUE)
  top <- max(log_weights)
  weights <- exp(log_weights - top)
  c(estimate = top + log(mean(weights)),
    se = stats::sd(weights) / mean(weights) / sqrt(length(weights)))
}

test_that("both forms estimate the same likelihood", {
  for (family in c("binomial", "poisson")) {
    survey <- integrated_survey(family)
    make <- function(integrate) {
      make_adaptive_mcml_engine(survey$D, survey$y, survey$units,
                                survey$distances, family = family,
                                fix_tau2 = TRUE, integrate_noise = integrate)
    }
    sampled <- make(FALSE)
    integrated <- make(TRUE)
    expect_false(sampled$integrated)
    expect_true(integrated$integrated)
    expect_identical(integrated$integrated_locations, 40L)
    set.seed(11)
    draws <- matrix(rnorm(20000 * 40), ncol = 40)
    a <- estimate_log_likelihood(sampled, survey$parameters, draws)
    b <- estimate_log_likelihood(integrated, survey$parameters, draws)
    expect_lt(abs(a[["estimate"]] - b[["estimate"]]),
              4 * sqrt(a[["se"]]^2 + b[["se"]]^2))
    # The covariance of the location effects with the whole nugget is the same.
    expect_equal(integrated$spatial_covariance(survey$parameters),
                 sampled$spatial_covariance(survey$parameters))
  }
})

test_that("the nugget stays sampled where a location keeps several records", {
  survey <- integrated_survey("binomial")
  # A second record at location 1, with another covariate value.
  D <- rbind(survey$D, c(1, 2))
  y <- c(survey$y, 5)
  units <- c(survey$units, 20)
  ID_coords <- c(seq_len(40), 1L)
  make <- function(integrate) {
    make_adaptive_mcml_engine(D, y, units, survey$distances,
                              family = "binomial", ID_coords = ID_coords,
                              fix_tau2 = TRUE, integrate_noise = integrate)
  }
  sampled <- make(FALSE)
  integrated <- make(TRUE)
  expect_identical(integrated$integrated_locations, 39L)
  set.seed(12)
  draws <- matrix(rnorm(20000 * 40), ncol = 40)
  a <- estimate_log_likelihood(sampled, survey$parameters, draws)
  b <- estimate_log_likelihood(integrated, survey$parameters, draws)
  expect_lt(abs(a[["estimate"]] - b[["estimate"]]),
            4 * sqrt(a[["se"]]^2 + b[["se"]]^2))
})

test_that("a nugget held at zero gives the model without a nugget", {
  survey <- integrated_survey("binomial")
  integrated <- make_adaptive_mcml_engine(
    survey$D, survey$y, survey$units, survey$distances, family = "binomial",
    fix_tau2 = TRUE, integrate_noise = TRUE
  )
  without <- make_adaptive_mcml_engine(
    survey$D, survey$y, survey$units, survey$distances, family = "binomial",
    fix_tau2 = FALSE
  )
  expect_identical(integrated$hold_log_ratio, -20)
  held <- replace(survey$parameters, 5L, -20)
  set.seed(13)
  draws <- matrix(rnorm(500 * 40), ncol = 40)
  expect_equal(
    integrated$evaluate(held, draws, log_weights_only = TRUE),
    without$evaluate(survey$parameters[1:4], draws, log_weights_only = TRUE)
  )
  expect_identical(
    hold_noise_at_zero(integrated, survey$parameters, boundary = 5L), held
  )
  expect_identical(
    hold_noise_at_zero(integrated, survey$parameters, boundary = 3L),
    survey$parameters
  )
  expect_identical(
    hold_noise_at_zero(without, survey$parameters[1:4], boundary = 3L),
    survey$parameters[1:4]
  )
})

test_that("integrate_noise is ignored without a nugget", {
  survey <- integrated_survey("binomial")
  without <- make_adaptive_mcml_engine(
    survey$D, survey$y, survey$units, survey$distances, family = "binomial",
    fix_tau2 = FALSE, integrate_noise = TRUE
  )
  expect_false(without$integrated)
})

test_that("the table of a custom link matches numerical integration", {
  links <- list(
    binomial = list(probit = stats::pnorm,
                    cloglog = function(eta) -expm1(-exp(eta)),
                    cauchit = stats::pcauchy),
    poisson = list(softplus = function(eta) log1p(exp(eta)))
  )
  for (family in names(links)) for (name in names(links[[family]])) {
    inverse <- links[[family]][[name]]
    linkf <- make_inverse_link_functions(family, inverse, 3L)
    y <- if (family == "poisson") c(0, 4, 30) else c(0, 3, 20)
    m <- if (family == "poisson") c(1, 1, 2) else c(20, 20, 20)
    sites <- make_site_table(y, m, family, linkf = linkf)
    for (v in c(0.004, 0.3, 4)) for (x in c(-4, -0.7, 1.2)) {
      values <- cpp_site_evaluate(sites$pointer, v, sites$pair,
                                  matrix(x, 3L, 1L), TRUE, 0L)$h
      reference <- vapply(seq_along(y), function(i) {
        log_integrand <- function(e) {
          mean <- inverse(x + e)
          kernel <- if (family == "poisson") {
            y[i] * log(m[i] * mean) - m[i] * mean
          } else {
            y[i] * log(pmax(mean, 2^-53)) +
              (m[i] - y[i]) * log(pmax(1 - mean, 2^-53))
          }
          kernel + stats::dnorm(e, 0, sqrt(v), log = TRUE)
        }
        mode <- stats::optimize(log_integrand, c(-30, 30), maximum = TRUE,
                                tol = 1e-12)$maximum
        top <- log_integrand(mode)
        width <- min(sqrt(v), 1)
        log(stats::integrate(function(e) exp(log_integrand(e) - top),
                             mode - 60 * width, mode + 60 * width,
                             rel.tol = 1e-12, subdivisions = 2000L)$value) +
          top
      }, numeric(1))
      expect_equal(as.numeric(values), reference, tolerance = 1e-6,
                   label = paste(family, name, v, x))
    }
  }
})

test_that("both forms estimate the same likelihood with a custom link", {
  survey <- integrated_survey("binomial")
  make <- function(integrate) {
    make_adaptive_mcml_engine(survey$D, survey$y, survey$units,
                              survey$distances, family = "binomial",
                              fix_tau2 = TRUE, integrate_noise = integrate,
                              invlink = stats::pnorm)
  }
  sampled <- make(FALSE)
  integrated <- make(TRUE)
  parameters <- survey$parameters
  parameters[1L] <- -0.9
  set.seed(14)
  draws <- matrix(rnorm(20000 * 40), ncol = 40)
  a <- estimate_log_likelihood(sampled, parameters, draws)
  b <- estimate_log_likelihood(integrated, parameters, draws)
  expect_lt(abs(a[["estimate"]] - b[["estimate"]]),
            4 * sqrt(a[["se"]]^2 + b[["se"]]^2))
})

test_that("a fit with integrate_noise reports the integrated locations", {
  skip_on_cran()
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(nugget = TRUE), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    control_mcml = set_control_mcml(integrate_noise = TRUE, max_stages = 1),
    seed = fixture_seed, messages = FALSE
  ))
  optimizer <- attr(fit, "optimizer")
  expect_true(optimizer$mcml_integrate_noise)
  expect_identical(optimizer$mcml_integrated_locations, 10L)
  expect_null(optimizer$laplace_start)
  expect_match(check_fit(fit)$draws,
               "The nugget was integrated out at 10 of 10 locations.")
  expect_true(is.finite(optimizer$mcml_history[[1]]$refit_change))
})

test_that("check_fit suggests integrating a sampled nugget out", {
  bank <- list(effective_draws = 40, pareto_k = 0.9)
  expect_match(fit_check_weights(bank, suggest_integration = TRUE)$fix,
               "integrate_noise = TRUE")
  expect_no_match(fit_check_weights(bank)$fix, "integrate_noise")
})

test_that("integration is suggested only for a nugget variance of 0.1 or more", {
  fit <- function(tau2, fix_tau2 = TRUE, sigma2 = 1.5) {
    list(fix_tau2 = fix_tau2, estimate = list(
      sigma2 = log(sigma2),
      nu2 = if (isTRUE(fix_tau2)) log(tau2 / sigma2)
    ))
  }
  sampled <- list(mcml_integrate_noise = FALSE)
  expect_true(suggests_integration(fit(0.5), sampled))
  expect_true(suggests_integration(fit(0.1 + 1e-9), sampled))
  expect_false(suggests_integration(fit(0.05), sampled))
  # A nugget held at zero.
  expect_false(suggests_integration(fit(1.5 * exp(-20)), sampled))
  expect_false(suggests_integration(fit(0.5),
                                    list(mcml_integrate_noise = TRUE)))
  expect_true(suggests_integration(fit(NA, fix_tau2 = 0.2), sampled))
  expect_false(suggests_integration(fit(NA, fix_tau2 = 0.02), sampled))
  expect_false(suggests_integration(fit(NA, fix_tau2 = FALSE), sampled))
})

test_that("a global iteration stopped near its fixed point is accepted only when asked", {
  stopped <- function(change) {
    function(mode, mu, ...) {
      list(centre = mode, precision = diag(length(mode)), root = NULL,
           factors = list(mean = mode), iterations = 501L, change = change,
           converged = FALSE)
    }
  }
  find_mode <- function() list(mode = c(0, 0), hessian = diag(2))
  build <- function(change, tolerance) {
    local_mocked_bindings(fit_global_proposal = stopped(change))
    build_gaussian_proposal("global", 0, 0, 1, NULL, diag(2), "binomial",
                            NULL, TRUE, NULL, find_mode, tolerance)
  }
  accepted <- build(5e-3, 1e-2)
  expect_false(accepted$failed)
  expect_identical(accepted$stalled_change, 5e-3)
  expect_false(is.null(accepted$factors))
  far <- build(5e-2, 1e-2)
  expect_true(far$failed)
  expect_identical(far$failure$reason, "no convergence")
  expect_null(far$stalled_change)
  # Without a tolerance, as for the sampled form, the rule of the release.
  release <- build(5e-3, NULL)
  expect_true(release$failed)
  expect_null(release$stalled_change)
})

test_that("only the integrated engine accepts a stalled global iteration", {
  survey <- integrated_survey("binomial")
  make <- function(integrate) {
    make_adaptive_mcml_engine(survey$D, survey$y, survey$units,
                              survey$distances, family = "binomial",
                              fix_tau2 = TRUE, integrate_noise = integrate)
  }
  sampled <- make(FALSE)
  integrated <- make(TRUE)
  real <- fit_global_proposal
  local_mocked_bindings(fit_global_proposal = function(...) {
    global <- real(...)
    if (!is.null(global)) {
      global$converged <- FALSE
      global$change <- 5e-3
    }
    global
  })
  sampled$set_anchor(survey$parameters)
  expect_identical(sampled$proposal_failures(), 1L)
  expect_identical(sampled$stalled_anchors()$count, 0L)
  integrated$set_anchor(survey$parameters)
  expect_identical(integrated$proposal_failures(), 0L)
  expect_identical(integrated$stalled_anchors(), list(count = 1L, change = 5e-3))
  mark <- integrated$mark_proposal_failures()
  integrated$reset_proposal_failures()
  expect_identical(integrated$stalled_anchors()$count, 0L)
  integrated$restore_proposal_failures(mark)
  expect_identical(integrated$stalled_anchors()$count, 1L)
})

test_that("a sampled fit keeps the release's rule when the global iteration stalls", {
  skip_on_cran()
  real <- fit_global_proposal
  local_mocked_bindings(fit_global_proposal = function(...) {
    global <- real(...)
    if (!is.null(global)) {
      global$converged <- FALSE
      global$change <- 5e-3
    }
    global
  })
  fit_with <- function(integrate) {
    suppressWarnings(glgpm(
      y ~ cov + gp(nugget = TRUE), data = binomial_data, family = "binomial",
      denominator = denominator, control_mcmc = control_mcmc,
      control_mcml = set_control_mcml(integrate_noise = integrate,
                                      max_stages = 1),
      seed = fixture_seed, messages = FALSE
    ))
  }
  sampled <- attr(fit_with(FALSE), "optimizer")
  expect_gt(sampled$mcml_history[[1]]$proposal_failures, 0L)
  expect_true(isTRUE(sampled$mcml_history[[1]]$proposal_fallback))
  expect_null(sampled$mcml_stalled_anchors)
  expect_null(sampled$mcml_history[[1]]$stalled_anchors)
  integrated <- attr(fit_with(TRUE), "optimizer")
  expect_identical(integrated$mcml_history[[1]]$proposal_failures, 0L)
  expect_null(integrated$mcml_history[[1]]$proposal_fallback)
  expect_gt(integrated$mcml_history[[1]]$stalled_anchors, 0L)
  expect_identical(integrated$mcml_history[[1]]$stalled_change, 5e-3)
  expect_identical(integrated$mcml_stalled_anchors$stages, 1L)
})

test_that("stalled anchors are summarised over the stages", {
  expect_null(summarise_stalled_anchors(list(list(), list())))
  history <- list(list(), list(stalled_anchors = 2L, stalled_change = 4e-3),
                  list(stalled_anchors = 1L, stalled_change = 7e-3))
  expect_identical(summarise_stalled_anchors(history),
                   list(stages = 2:3, count = 3, change = 7e-3))
})

test_that("set_control_mcml accepts an automatic choice of the form", {
  expect_identical(set_control_mcml(integrate_noise = "auto")$integrate_noise,
                   "auto")
  expect_error(set_control_mcml(integrate_noise = "yes"),
               "TRUE, FALSE or \"auto\"")
  expect_error(set_control_mcml(integrate_noise = NA), "TRUE, FALSE or")
})

test_that("a switching engine takes its form from the nugget variance", {
  survey <- integrated_survey("binomial")
  make <- function(integrate) {
    make_adaptive_mcml_engine(
      survey$D, survey$y, survey$units, survey$distances, family = "binomial",
      fix_tau2 = TRUE, integrate_noise = integrate
    )
  }
  sampled <- make(FALSE)
  integrated <- make(TRUE)
  engine <- make_switching_mcml_engine(sampled, integrated)
  expect_identical(engine$form(), "sampled")
  expect_false(engine$integrated)
  expect_null(engine$hold_log_ratio)
  expect_identical(engine$integrated_locations, 0L)
  expect_identical(engine$latent_dimension, sampled$latent_dimension)
  # The survey's nugget variance is 0.7.
  large <- survey$parameters
  small <- replace(large, 5L, log(0.05 / 1.2))
  set.seed(15)
  draws <- matrix(rnorm(300 * 40), ncol = 40)
  engine$choose_form(large)
  expect_identical(engine$form(), "integrated")
  expect_true(engine$integrated)
  expect_identical(engine$hold_log_ratio, -20)
  expect_identical(engine$integrated_locations, 40L)
  expect_identical(engine$evaluate(large, draws, log_weights_only = TRUE),
                   integrated$evaluate(large, draws, log_weights_only = TRUE))
  engine$choose_form(small)
  expect_identical(engine$form(), "sampled")
  expect_identical(engine$evaluate(small, draws, log_weights_only = TRUE),
                   sampled$evaluate(small, draws, log_weights_only = TRUE))
  # A nugget held at zero takes the sampled form.
  engine$choose_form(replace(large, 5L, -20))
  expect_identical(engine$form(), "sampled")
  expect_identical(make_switching_mcml_engine(sampled, integrated,
                                              threshold = 0.04)$
                     choose_form(small), "integrated")
  engine$use_proposal("laplace")
  expect_identical(c(sampled$proposal_type(), integrated$proposal_type()),
                   c("laplace", "laplace"))
  engine$use_proposal("global")
  expect_identical(
    engine$mode_search_counts()$searches,
    sampled$mode_search_counts()$searches +
      integrated$mode_search_counts()$searches
  )
})

test_that("an automatic fit records the form of each stage", {
  skip_on_cran()
  fit <- suppressWarnings(glgpm(
    y ~ cov + gp(nugget = TRUE), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcmc = control_mcmc,
    control_mcml = set_control_mcml(integrate_noise = "auto", max_stages = 2),
    seed = fixture_seed, messages = FALSE
  ))
  optimizer <- attr(fit, "optimizer")
  forms <- optimizer$mcml_noise_forms
  expect_length(forms, optimizer$mcml_stages)
  expect_true(all(forms %in% c("sampled", "integrated")))
  expect_identical(optimizer$mcml_integrate_noise,
                   identical(utils::tail(forms, 1L), "integrated"))
  # The first stage takes its form from the Laplace estimate.
  laplace <- optimizer$laplace_start$estimate
  expect_false(is.null(laplace))
  expect_identical(forms[1L], if (exp(laplace[3L] + laplace[5L]) >= 0.1) {
    "integrated"
  } else {
    "sampled"
  })
  expect_false(summarise_mcml_fit(fit)$suggest_integration)
  expect_match(check_fit(fit)$draws, "The automatic choice")
})

test_that("the forms of an automatic fit are described by stage", {
  expect_null(describe_noise_forms(NULL))
  expect_identical(describe_noise_forms("integrated"),
                   " The automatic choice integrated the nugget out.")
  expect_identical(describe_noise_forms(c("sampled", "sampled")),
                   " The automatic choice sampled the nugget in every stage.")
  expect_identical(
    describe_noise_forms(c("sampled", "integrated", "integrated")),
    paste0(" The automatic choice sampled the nugget in stage 1 and ",
           "integrated it out in stages 2 and 3.")
  )
  expect_identical(
    describe_noise_forms(c("integrated", "sampled", "integrated", "sampled")),
    paste0(" The automatic choice sampled the nugget in stages 2 and 4 and ",
           "integrated it out in stages 1 and 3.")
  )
})
