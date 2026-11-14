# Groups of units whose linear predictors differ, such as the children of a
# household who differ in sex, integrated out (integrate_noise), with rule (c)
# of plan_integrated_records() and the group and location tables of
# site_table.h.

# A location with two households of children who differ in sex and age and a
# third household of one child.
group_survey <- function() {
  list(location = rep(1L, 5L), household = c(1L, 1L, 1L, 2L, 3L),
       D = cbind(1, sex = c(0, 1, 1, 0, 1), age = c(0.2, -1, 0.5, 0, 1.1)),
       y = c(1, 0, 1, 0, 1), m = c(1, 1, 2, 1, 1))
}

# The plan, table and members of a survey's records at the coefficients beta.
group_table <- function(survey, allowed, family = make_family("binomial"),
                        beta = c(-0.6, 0.5, -0.3)) {
  n_loc <- max(survey$location)
  plan <- plan_integrated_records(
    survey$y, survey$m, survey$D, numeric(length(survey$y)), survey$location,
    cbind(household = survey$household), n_loc, allowed
  )
  sites <- make_site_table(plan$y, plan$units_m, family, kind = plan$kind,
                           a = plan$a, b = plan$b, patterns = plan$patterns,
                           groups = plan$groups)
  cpp_site_set_relative(sites$pointer, relative_predictors(sites, beta))
  members <- lapply(seq_along(plan$kind), function(r) {
    groups <- plan$groups[[r]] %||% list(list(
      y = plan$y[r], m = plan$units_m[r], D = matrix(0, 1L, length(beta)),
      offset = 0
    ))
    lapply(groups, function(g) {
      list(y = g$y, m = g$m, r = as.numeric(g$D %*% beta + g$offset))
    })
  })
  list(plan = plan, sites = sites, members = members,
       mu = as.numeric(plan$D %*% beta))
}

# The Gauss-Hermite rule with `nodes` nodes for the standard normal
# distribution (Golub and Welsch, 1969).
standard_rule <- function(nodes) {
  jacobi <- matrix(0, nodes, nodes)
  off <- sqrt(seq_len(nodes - 1L))
  jacobi[cbind(seq_len(nodes - 1L), 2:nodes)] <- off
  jacobi[cbind(2:nodes, seq_len(nodes - 1L))] <- off
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = decomposition$values, weights = decomposition$vectors[1L, ]^2)
}

log_sum_exp_columns <- function(a) {
  top <- apply(a, 2L, max)
  top + log(colSums(exp(a - rep(top, each = nrow(a)))))
}

# log of the integral over N(z; 0, vo) of the product over groups of the
# integral over N(u; 0, vi) of the product of the members' likelihoods at
# x + r + z + u, by Gauss-Hermite rules on both Gaussians, independent of the
# package's tables, for the log-likelihood `kernel(y, m, eta)` of a member. A
# variance of zero takes the value at zero.
nested_groups_by_quadrature <- function(groups, x, vo, vi, kernel,
                                        nodes = 80L) {
  rule <- standard_rule(nodes)
  log_w <- log(rule$weights)
  outer_nodes <- sqrt(vo) * rule$nodes
  inner_nodes <- sqrt(vi) * rule$nodes
  shifts <- outer(inner_nodes, outer_nodes, `+`)
  vapply(x, function(point) {
    # One row per inner node, one column per outer node.
    total <- Reduce(`+`, lapply(groups, function(g) {
      values <- Reduce(`+`, lapply(seq_along(g$y), function(j) {
        kernel(g$y[j], g$m[j], point + g$r[j] + shifts)
      }))
      log_sum_exp_columns(values + log_w)
    }))
    log_sum_exp_columns(matrix(total + log_w, ncol = 1L))
  }, numeric(1))
}

logit_kernel <- function(y, m, eta) {
  y * eta - m * (pmax(eta, 0) + log1p(exp(-abs(eta))))
}

link_kernel <- function(inverse) {
  function(y, m, eta) {
    p <- inverse(eta)
    y * log(pmax(p, 2^-53)) + (m - y) * log(pmax(1 - p, 2^-53))
  }
}

test_that("children who differ in covariates form groups by the rules", {
  # Location 1 has households 1 (two children who differ in sex) and 2 (one
  # child). Location 2 has one household of two children who differ. Location
  # 3 has one child. At location 4 every child has the same covariates, so
  # the two children of household 5 are one unit and the location follows
  # rule (b). Household 6 is recorded at locations 5 and 6.
  location <- c(1L, 1L, 1L, 2L, 2L, 3L, 4L, 4L, 4L, 5L, 6L, 6L)
  household <- c(1L, 1L, 2L, 3L, 3L, 4L, 5L, 5L, 7L, 6L, 6L, 8L)
  sex <- c(0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0)
  y <- c(1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1)
  plan <- plan_integrated_records(
    y, rep(1, 12), cbind(1, sex = sex), numeric(12), location,
    cbind(household = household), 6L,
    c(nugget = TRUE, "re(household)" = TRUE)
  )
  # Location 6 keeps its nugget sampled, since household 6 spans two
  # locations.
  expect_identical(plan$kind, c(4L, 3L, 1L, 2L, 1L, 0L, 1L))
  expect_identical(plan$ID_coords, c(1L, 2L, 3L, 4L, 5L, 6L, 6L))
  expect_identical(plan$nugget_integrated,
                   c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(lengths(plan$groups), c(2L, 1L, 0L, 0L, 0L, 0L, 0L))
  expect_identical(plan$groups[[1L]][[1L]]$D, cbind(0, sex = c(0, 1)))
  # The location of one household is one group under the nugget and the
  # household together, the class of a record with one unit.
  expect_identical(plan$a[2L], plan$a[3L])
  expect_identical(unname(plan$classes[plan$a[2L], ]), c(TRUE, TRUE))
  expect_identical(unname(plan$classes[plan$a[1L], ]), c(TRUE, FALSE))
  expect_identical(unname(plan$classes[plan$b[1L], ]), c(FALSE, TRUE))
  expect_identical(unname(plan$ID_re[, 1L]), c(NA, 3L, 4L, NA, 6L, 6L, 8L))
  expect_identical(plan$inner, list(c(1L, 4L)))
  expect_identical(plan$sampled_levels, 1L)
  expect_identical(plan$counts$integrated, c(5L, 7L))
  # Without the nugget, each household of several children is a group of its
  # own and a household of one child a record of rule (a).
  without <- plan_integrated_records(
    y, rep(1, 12), cbind(1, sex = sex), numeric(12), location,
    cbind(household = household), 6L,
    c(nugget = FALSE, "re(household)" = TRUE)
  )
  expect_identical(without$kind, c(3L, 1L, 3L, 1L, 1L, 1L, 0L, 0L, 1L))
  expect_identical(unname(without$classes[without$a[1L], ]), c(FALSE, TRUE))
  expect_identical(without$inner, list(integer()))
  expect_identical(without$counts$integrated, c(0L, 7L))
})

test_that("rule (c) leaves a location sampled where two terms vary", {
  # Two households at one location, with an effect of each child, so that
  # both the households and the children vary there. Each child's effect is
  # integrated out of its own record by rule (a), and household 2's too.
  plan <- plan_integrated_records(
    c(1, 0, 1), rep(1, 3), cbind(1, c(0, 1, 0)), numeric(3), rep(1L, 3),
    cbind(household = c(1L, 1L, 2L), child = 1:3), 1L,
    c(nugget = TRUE, "re(household)" = TRUE, "re(child)" = TRUE)
  )
  expect_identical(plan$kind, c(1L, 1L, 1L))
  expect_identical(plan$nugget_integrated, FALSE)
  expect_identical(plan$sampled_levels, c(1L, 0L))
  # With the household shared, the children's effects are the groups.
  shared <- plan_integrated_records(
    c(1, 0), rep(1, 2), cbind(1, c(0, 1)), numeric(2), rep(1L, 2),
    cbind(household = c(1L, 1L), child = 1:2), 1L,
    c(nugget = TRUE, "re(household)" = TRUE, "re(child)" = TRUE)
  )
  expect_identical(shared$kind, 4L)
  expect_identical(unname(shared$classes[shared$a, ]), c(TRUE, TRUE, FALSE))
  expect_identical(unname(shared$classes[shared$b, ]), c(FALSE, FALSE, TRUE))
})

test_that("the plan of units that share their predictors is unchanged", {
  # With the same covariate for every child of a location, rule (c) does not
  # apply and rules (a) and (b) give the plan.
  location <- c(1L, 1L, 1L, 2L, 2L, 3L)
  household <- c(1L, 1L, 2L, 3L, 4L, 5L)
  x <- c(0.3, 0.3, 0.3, -1, -1, 2)
  plan <- plan_integrated_records(
    c(1, 0, 1, 0, 0, 1), rep(1, 6), cbind(1, x), numeric(6), location,
    cbind(household = household), 3L,
    c(nugget = TRUE, "re(household)" = TRUE)
  )
  expect_identical(plan$kind, c(2L, 2L, 1L))
  expect_true(all(vapply(plan$groups, is.null, logical(1))))
})

test_that("groups and locations of groups match numerical integration", {
  survey <- group_survey()
  x <- c(-2.1, -0.3, 1.2)
  for (outer in c(TRUE, FALSE)) {
    table <- group_table(survey, c(nugget = outer, "re(household)" = TRUE))
    plan <- table$plan
    expect_identical(plan$kind, if (outer) 4L else c(3L, 1L, 1L))
    for (variances in list(c(0.1, 0.15), c(0.8, 0.5), c(0.3, 2.5))) {
      vo <- if (outer) variances[1L] else 0
      vi <- variances[2L]
      classes <- as.numeric((plan$classes + 0) %*% c(vo, vi))
      eta <- matrix(rep(table$mu, each = length(x)) + x, length(x))
      values <- cpp_site_evaluate(table$sites$pointer, classes,
                                  table$sites$units, eta, FALSE, 2L)
      for (r in seq_along(plan$kind)) {
        reference <- nested_groups_by_quadrature(table$members[[r]], eta[, r],
                                                 vo, vi, logit_kernel)
        expect_equal(values$h[, r], reference, tolerance = 1e-8)
      }
      # Derivatives in the linear predictor, by central differences.
      step <- 1e-4
      up <- cpp_site_evaluate(table$sites$pointer, classes, table$sites$units,
                              eta + step, FALSE, 2L)
      down <- cpp_site_evaluate(table$sites$pointer, classes,
                                table$sites$units, eta - step, FALSE, 2L)
      expect_equal(values$h1, (up$h - down$h) / (2 * step), tolerance = 1e-6)
      expect_equal(values$h2, (up$h1 - down$h1) / (2 * step),
                   tolerance = 1e-6)
    }
  }
})

test_that("a location of groups has the slope of its inner variance", {
  survey <- group_survey()
  table <- group_table(survey, c(nugget = TRUE, "re(household)" = TRUE))
  eta <- matrix(table$mu + c(-1.4, 0.6), ncol = 1L)
  members <- table$members[[1L]]
  step <- 1e-5
  for (variances in list(c(0.1, 0.3), c(0.9, 1.6), c(0.7, 0), c(0, 0.6),
                         c(0, 0))) {
    classes <- as.numeric((table$plan$classes + 0) %*% variances)
    slopes <- cpp_site_inner_slopes(table$sites$pointer, classes,
                                    table$sites$units, eta)
    value_at <- function(vi) {
      nested_groups_by_quadrature(members, eta, variances[1L], vi,
                                  logit_kernel)
    }
    # At an inner variance of zero the slope is one-sided.
    numeric_slope <- if (variances[2L] > 0) {
      (value_at(variances[2L] + step) - value_at(variances[2L] - step)) /
        (2 * step)
    } else {
      (value_at(step) - value_at(0)) / step
    }
    expect_equal(as.numeric(slopes), numeric_slope, tolerance = 1e-4)
    # Held variances give the integrals of what remains.
    values <- cpp_site_evaluate(table$sites$pointer, classes,
                                table$sites$units, eta, FALSE, 0L)$h
    expect_equal(as.numeric(values), value_at(variances[2L]),
                 tolerance = 1e-8)
  }
})

test_that("groups under a custom link match numerical integration", {
  survey <- group_survey()
  links <- list(probit = stats::pnorm,
                imperfect = function(eta) 0.01 + 0.89 * stats::plogis(eta))
  x <- c(-1.1, 0.4)
  for (name in names(links)) {
    linkf <- if (name == "probit") {
      make_inverse_link_functions("binomial", named_inverse_link("probit"),
                                  5L)
    } else {
      make_inverse_link_functions("binomial", links[[name]], 5L)
    }
    for (outer in c(TRUE, FALSE)) {
      table <- group_table(survey, c(nugget = outer, "re(household)" = TRUE),
                           family = make_family("binomial", linkf))
      for (variances in list(c(0.1, 0.3), c(0.9, 1.6))) {
        vo <- if (outer) variances[1L] else 0
        classes <- as.numeric((table$plan$classes + 0) %*%
                                c(vo, variances[2L]))
        eta <- matrix(rep(table$mu, each = length(x)) + x, length(x))
        values <- cpp_site_evaluate(table$sites$pointer, classes,
                                    table$sites$units, eta, FALSE, 0L)$h
        for (r in seq_along(table$plan$kind)) {
          reference <- nested_groups_by_quadrature(
            table$members[[r]], eta[, r], vo, variances[2L],
            link_kernel(links[[name]])
          )
          expect_equal(values[, r], reference, tolerance = 1e-7)
        }
      }
    }
  }
})

test_that("group tables are a function of the coefficients alone", {
  survey <- group_survey()
  table <- group_table(survey, c(nugget = TRUE, "re(household)" = TRUE))
  sites <- table$sites
  eta <- matrix(table$mu + seq(-2, 2, by = 0.37), ncol = 1L)
  first <- relative_predictors(sites, c(-0.6, 0.5, -0.3))
  second <- relative_predictors(sites, c(-0.6, 0.9, 0.2))
  evaluate <- function() {
    cpp_site_evaluate(sites$pointer, c(0.4, 0.7), sites$units, eta, FALSE,
                      2L)
  }
  at_first <- evaluate()
  cpp_site_set_relative(sites$pointer, second)
  at_second <- evaluate()
  expect_false(isTRUE(all.equal(at_first$h, at_second$h)))
  cpp_site_set_relative(sites$pointer, first)
  expect_identical(evaluate(), at_first)
  # A fresh table, filled in another order, gives the same values.
  fresh <- group_table(survey, c(nugget = TRUE, "re(household)" = TRUE))
  backwards <- cpp_site_evaluate(fresh$sites$pointer, c(0.4, 0.7),
                                 fresh$sites$units,
                                 eta[rev(seq_len(nrow(eta))), , drop = FALSE],
                                 FALSE, 2L)
  expect_identical(backwards$h[rev(seq_len(nrow(eta))), , drop = FALSE],
                   at_first$h)
  info <- cpp_site_table_info(sites$pointer)
  # The location's three households and the group of all its children.
  expect_identical(info$groups, 4L)
  expect_identical(info$group_locations, 1L)
  expect_gt(info$location_cells_filled, 0)
})

test_that("a household keeps its integral when shared coefficients change", {
  # Two households at one location with a household covariate, whose
  # children differ in sex.
  survey <- list(location = rep(1L, 4L), household = c(1L, 1L, 2L, 2L),
                 D = cbind(1, sex = c(0, 1, 0, 1), wealth = c(0, 0, 1, 1)),
                 y = c(1, 0, 0, 1), m = rep(1, 4))
  table <- group_table(survey, c(nugget = TRUE, "re(household)" = TRUE),
                       beta = c(-0.6, 0.5, 0.2))
  sites <- table$sites
  eta <- matrix(table$mu + c(-0.5, 0.5), ncol = 1L)
  evaluate <- function() {
    cpp_site_evaluate(sites$pointer, c(0.4, 0.7), sites$units, eta, FALSE,
                      0L)$h
  }
  evaluate()
  before <- cpp_site_table_info(sites$pointer)
  # A new coefficient of the household covariate moves the second household
  # relative to the first but changes neither household's own integral.
  cpp_site_set_relative(sites$pointer,
                        relative_predictors(sites, c(-0.6, 0.5, 0.9)))
  evaluate()
  after <- cpp_site_table_info(sites$pointer)
  expect_identical(after$group_tables_built, before$group_tables_built)
  expect_identical(after$location_tables_built,
                   before$location_tables_built + 1L)
  # A new coefficient of sex changes both.
  cpp_site_set_relative(sites$pointer,
                        relative_predictors(sites, c(-0.6, 0.1, 0.9)))
  evaluate()
  expect_identical(cpp_site_table_info(sites$pointer)$group_tables_built,
                   after$group_tables_built + 2L)
})

test_that("a table with groups needs its relative predictors", {
  survey <- group_survey()
  plan <- plan_integrated_records(
    survey$y, survey$m, survey$D, numeric(5), survey$location,
    cbind(household = survey$household), 1L,
    c(nugget = TRUE, "re(household)" = TRUE)
  )
  sites <- make_site_table(plan$y, plan$units_m, make_family("binomial"),
                           kind = plan$kind, a = plan$a, b = plan$b,
                           groups = plan$groups)
  expect_error(cpp_site_evaluate(sites$pointer, c(0.4, 0.7), sites$units,
                                 matrix(0), TRUE, 0L), "not set")
  expect_error(cpp_site_set_relative(sites$pointer, 1), "each member")
  expect_error(make_integrated_kernel(sites, c(0.4, 0.7)),
               class = "riskmap_error_internal")
  kernel <- make_integrated_kernel(sites, c(0.4, 0.7), c(-0.6, 0.5, -0.3))
  expect_true(is.finite(kernel$terms(-0.4)$log_likelihood))
})

# A small survey of households of children who differ in sex, for both forms
# of the engine.
child_survey <- function(n = 25L, seed = 5L) {
  set.seed(seed)
  coords <- cbind(runif(n), runif(n)) * 50
  distances <- as.matrix(dist(coords))
  S <- drop(t(chol(exp(-distances / 8) + diag(1e-9, n))) %*% rnorm(n))
  households <- 1L + rpois(n, 2)
  hh_location <- rep(seq_len(n), households)
  children <- 1L + rpois(length(hh_location), 1.5)
  household <- rep(seq_along(hh_location), children)
  location <- hh_location[household]
  sex <- rbinom(length(household), 1, 0.5)
  eta <- -1 + 0.4 * sex + S[location] + rnorm(n, 0, 0.6)[location] +
    rnorm(length(hh_location), 0, 0.8)[household]
  list(D = cbind(1, sex), y = rbinom(length(eta), 1, plogis(eta)),
       units = rep(1, length(eta)), distances = distances,
       ID_coords = location, ID_re = cbind(household = household),
       parameters = c(-1, 0.4, 0, log(8), log(0.36), log(0.64)))
}

test_that("both forms estimate the same likelihood with groups", {
  survey <- child_survey()
  for (nugget in c(TRUE, FALSE)) {
    parameters <- if (nugget) survey$parameters else survey$parameters[-5L]
    make <- function(integrate) {
      make_adaptive_mcml_engine(survey$D, survey$y, survey$units,
                                survey$distances, family = "binomial",
                                ID_coords = survey$ID_coords,
                                ID_re = survey$ID_re, fix_tau2 = nugget,
                                integrate_noise = integrate)
    }
    sampled <- make(FALSE)
    integrated <- make(TRUE)
    expect_true(integrated$integrated)
    # Every household is integrated out, so the latent effects are the
    # locations.
    expect_identical(integrated$latent_dimension, 25L)
    expect_identical(integrated$integration_counts$integrated[
      integrated$integration_counts$term == "re(household)"],
      max(survey$ID_re))
    set.seed(13)
    draws_sampled <- matrix(rnorm(20000 * sampled$latent_dimension), 20000)
    draws_integrated <- matrix(rnorm(20000 * 25), 20000)
    estimate <- function(engine, draws) {
      log_weights <- engine$evaluate(parameters, draws,
                                     log_weights_only = TRUE)
      top <- max(log_weights)
      weights <- exp(log_weights - top)
      c(estimate = top + log(mean(weights)),
        se = stats::sd(weights) / mean(weights) / sqrt(length(weights)))
    }
    a <- estimate(sampled, draws_sampled)
    b <- estimate(integrated, draws_integrated)
    expect_lt(abs(a[["estimate"]] - b[["estimate"]]),
              4 * sqrt(a[["se"]]^2 + b[["se"]]^2))
    # The compiled weights and those in R agree.
    proposal <- integrated$proposal_at(parameters)
    expect_true(is.finite(sum(proposal$centre)))
  }
})
