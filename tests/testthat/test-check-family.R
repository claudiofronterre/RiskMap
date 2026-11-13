# The checks of a family's likelihood, mean and simulator (#253, step 5). Every
# family and link of the table is checked here, so a family or link added to
# riskmap_families is checked as soon as it is there.

family_link <- function(family, link) {
  if (identical(link, family_traits(family)$default_link)) return(NULL)
  make_inverse_link_functions(family, c(named_inverse_link(link), link = link),
                              1L)
}

# The values of the family's parameters at which it is checked.
check_parameters <- function(family) {
  name <- names(family_traits(family)$parameters)
  if (identical(name, "sigma2_me")) return(list(list(sigma2_me = 0.5)))
  if (identical(name, "theta")) {
    return(list(list(theta = 0.5), list(theta = 50)))
  }
  list(NULL)
}

test_that("every family and link of the table passes its checks", {
  for (family in names(riskmap_families)) {
    traits <- family_traits(family)
    expect_true(all(traits$concave_links %in% traits$links))
    expect_true(all(traits$radius_links %in% traits$links))
    for (link in traits$links) {
      for (parameters in check_parameters(family)) {
        checked <- check_family(family, family_link(family, link), parameters)
        label <- paste(family, link, format(unlist(parameters)))
        expect(length(checked$problems) == 0L,
               paste0(label, ": ", paste(checked$problems, collapse = " ")))
        # A link on the list of concave links has no convex record on the
        # grid, and the cauchit link, which is not on it, has one.
        expect(is.null(checked$convex_at) ==
                 (link %in% traits$concave_links || traits$exact),
               paste(label, "is convex on the grid but on the list, or the",
                     "reverse"))
      }
    }
  }
})

test_that("the checks find a wrong score, curvature, density or simulator", {
  terms <- canonical_binomial_terms
  traits <- family_traits
  local_mocked_bindings(canonical_binomial_terms = function(...) {
    result <- terms(...)
    result$score <- result$score * 1.01
    result
  })
  expect_match(check_family("binomial")$problems,
               "The score differs from a central difference", all = FALSE)
  local_mocked_bindings(canonical_binomial_terms = function(...) {
    result <- terms(...)
    result$curvature <- result$curvature + 0.01
    result
  })
  expect_match(check_family("binomial")$problems,
               "The curvature differs from a central difference", all = FALSE)
  local_mocked_bindings(canonical_binomial_terms = terms)

  local_mocked_bindings(family_traits = function(family) {
    result <- traits(family)
    result$log_density <- function(y, mean, units_m, parameters) {
      dbinom(y, units_m, pmin(1.05 * mean, 1), log = TRUE)
    }
    result
  })
  expect_match(check_family("binomial")$problems,
               "is not the log density plus a term free", all = FALSE)
  local_mocked_bindings(family_traits = function(family) {
    result <- traits(family)
    result$simulate <- function(means, units_m, parameters) {
      rpois(length(means), lambda = 1.05 * as.numeric(means) * units_m)
    }
    result
  })
  expect_match(check_family("poisson")$problems,
               "The draws of the simulator do not have the mean", all = FALSE)
  # The Gaussian family is checked at its identity link, with means of -1, 0
  # and 1, so a simulator that keeps means within [0, 1] is found.
  local_mocked_bindings(family_traits = function(family) {
    result <- traits(family)
    result$simulate <- function(means, units_m, parameters) {
      pmin(pmax(as.numeric(means), 0), 1) +
        rnorm(length(means), sd = sqrt(parameters$sigma2_me))
    }
    result
  })
  expect_match(check_family("gaussian",
                            parameters = list(sigma2_me = 0.5))$problems,
               "The draws of the simulator do not have the mean", all = FALSE)
  local_mocked_bindings(family_traits = traits)

  # A link whose mean leaves the range of the family.
  above_one <- list(inv = function(eta) 1.01 * plogis(eta),
                    d1 = function(eta) 1.01 * dlogis(eta),
                    d2 = function(eta) 1.01 * dlogis(eta) *
                      (1 - 2 * plogis(eta)),
                    name = "custom", link = NULL)
  expect_match(check_family("binomial", above_one)$problems,
               "The mean is outside the range of the family", all = FALSE)
})

test_that("the checks find where a link given as a function is convex", {
  imperfect <- make_inverse_link_functions(
    "binomial", function(eta) 0.01 + 0.89 * plogis(eta), 5L
  )
  # log p is convex where the floor of 0.01 takes over from 0.89 F.
  convex_at <- check_family("binomial", imperfect)$convex_at
  expect_gt(convex_at$eta, -6)
  expect_lt(convex_at$eta, -3)
  expect_identical(c(convex_at$y, convex_at$units_m), c(1, 1))
  # pnorm is log-concave, and the t link with 4 degrees of freedom is not.
  # Deriv knows neither pnorm given bare nor plogis, so their derivatives
  # are numerical, and the curvature of a record with failures, a difference
  # of two terms near one, loses its sign far into the upper tail, which the
  # second differences of the log-likelihood do not depend on. A cloglog
  # written with a cancellation is rounded in its lower tail.
  concave <- list(pnorm, plogis, function(eta) plogis(3 * eta),
                  function(eta) 1 - exp(-exp(eta)))
  for (inverse in concave) {
    link <- make_inverse_link_functions("binomial", inverse, 5L)
    expect_null(check_family("binomial", link)$convex_at)
  }
  t_link <- make_inverse_link_functions("binomial",
                                        function(eta) pt(eta, 4), 5L)
  expect_false(is.null(check_family("binomial", t_link)$convex_at))
})

test_that("the checks leave the random number state as it was", {
  set.seed(4)
  state <- .Random.seed
  check_family("negative_binomial", parameters = list(theta = 2))
  expect_identical(.Random.seed, state)
})
