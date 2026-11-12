## Links named by a family object, as in family = binomial("probit"), read by
## read_family() on entry to glgpm(), specify_glgpm() and maxim_integrand().

fit_binomial <- function(family, invlink = NULL) {
  glgpm(y ~ cov + gp(), data = binomial_data, family = family,
        invlink = invlink, denominator = denominator,
        control_mcmc = control_mcmc, seed = fixture_seed,
        control_mcml = fixture_control_mcml, messages = FALSE)
}

# What a fit holds apart from its call and formula, whose environments differ,
# its link functions, which are closures made afresh by each fit, and the
# times recorded anywhere in its history.
fit_without_timing <- function(fit) {
  drop_timing <- function(x) {
    if (!is.list(x)) return(x)
    if (!is.null(names(x))) {
      x <- x[!names(x) %in% c("timing", "elapsed_seconds", "seconds")]
    }
    x[] <- lapply(x, drop_timing)
    x
  }
  drop_timing(unclass(fit)[setdiff(names(fit), c("call", "formula",
                                                 "link_function"))])
}

# The same links as users would write them in 'invlink'. Bare pnorm and
# pcauchy are compiled, so their derivatives are found numerically.
link_functions <- list(
  probit = pnorm,
  cloglog = function(eta) 1 - exp(-exp(eta)),
  cauchit = pcauchy
)
named_fits <- lapply(stats::setNames(nm = names(link_functions)),
                     function(link) fit_binomial(binomial(link)))

test_that("strings, family functions and canonical family objects read alike", {
  for (family in list("binomial", binomial, binomial(), binomial("logit"))) {
    expect_identical(read_family(family),
                     list(family = "binomial", invlink = NULL))
  }
  expect_identical(read_family(poisson()),
                   list(family = "poisson", invlink = NULL))
  expect_identical(read_family(gaussian),
                   list(family = "gaussian", invlink = NULL))
  expect_identical(read_family(binomial(), plogis),
                   list(family = "binomial", invlink = plogis))
})

test_that("named links carry their name and exact derivatives", {
  eta <- c(-30, -8, -2.5, -0.7, 0.3, 1.9, 3)
  for (link in names(link_functions)) {
    read <- read_family(binomial(link))
    expect_identical(read$family, "binomial")
    expect_identical(read$invlink$link, link)
    expect_equal(read$invlink$inv(eta), binomial(link)$linkinv(eta),
                 tolerance = 1e-15)
    expect_equal(read$invlink$d1(eta),
                 compute_richardson_derivative(read$invlink$inv, eta),
                 tolerance = 1e-6)
    expect_equal(read$invlink$d2(eta),
                 compute_richardson_derivative(read$invlink$d1, eta),
                 tolerance = 1e-6)
  }
  # The second derivative of cloglog stays finite where exp(eta) overflows.
  expect_identical(named_inverse_link("cloglog")$d2(c(-1000, 710, 1000)),
                   c(0, 0, 0))
})

test_that("canonical family objects fit exactly as the strings do", {
  binomial_object <- suppressWarnings(glgpm(
    y ~ cov + gp() + re(i), data = binomial_data, family = binomial(),
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    messages = FALSE
  ))
  expect_identical(fit_without_timing(binomial_object),
                   fit_without_timing(binomial_model))
  expect_identical(binomial_object$link_function$name, "canonical")
  expect_identical(binomial_object$link_function$link, "logit")

  poisson_object <- suppressWarnings(glgpm(
    y ~ cov + gp() + re(i), data = poisson_data, family = poisson,
    denominator = denominator, control_mcmc = control_mcmc,
    seed = fixture_seed, control_mcml = fixture_control_mcml,
    return_samples = TRUE, messages = FALSE
  ))
  expect_identical(fit_without_timing(poisson_object),
                   fit_without_timing(poisson_model))
  expect_identical(poisson_object$link_function$link, "log")
})

test_that("a named link gives the estimates of its inverse given in invlink", {
  for (link in names(link_functions)) {
    named <- named_fits[[link]]
    supplied <- fit_binomial("binomial", link_functions[[link]])
    expect_identical(named$link_function$link, link)
    expect_null(supplied$link_function$link)
    expect_true(named$mcml_certified)
    expect_true(supplied$mcml_certified)
    # The derivatives only steer the importance sampler, so the two fits
    # differ by Monte Carlo error alone. Each estimate is within its
    # certificate, in standard errors, of the maximum likelihood estimate, so
    # the two are within the sum of the certificates of each other.
    standard_errors <- sqrt(diag(named$covariance))
    moved <- abs(unlist(named$estimate) - unlist(supplied$estimate)) /
      standard_errors[names(unlist(named$estimate))]
    expect_lt(max(moved), attr(named, "optimizer")$mcml_certificate +
                attr(supplied, "optimizer")$mcml_certificate)
    # With 500 draws and a relative effective sample size above 0.8, each
    # Monte Carlo log-likelihood has a standard error of about 0.01 to 0.02,
    # so 0.05 is two to four standard errors of their difference.
    expect_lt(abs(named$log_lik - supplied$log_lik), 0.05)
  }
})

test_that("a named link is never differentiated", {
  local_mocked_bindings(differentiate = function(...) {
    stop("differentiate() was called")
  }, .package = "RiskMap")
  # The fit includes its final conditional samples, and predict() draws
  # afresh, so every sampler is given the exact derivatives.
  expect_no_error(fit <- fit_binomial(binomial("probit")))
  expect_identical(fit_without_timing(fit),
                   fit_without_timing(named_fits$probit))
  expect_no_error(predict(fit, control_mcmc = control_mcmc, messages = FALSE))
  expect_error(fit_binomial("binomial", pnorm), "differentiate\\(\\) was called")
})

# The message of an input error, which must name `invlink` as the way out.
refused_message <- function(object) {
  error <- expect_input_error({{ object }}, "family")
  expect_match(conditionMessage(error), "`invlink`", fixed = TRUE)
  conditionMessage(error)
}

test_that("links that can leave the range of the family are refused", {
  expect_match(refused_message(fit_binomial(binomial("log"))),
               "binomial(\"log\")` can't be used, because the log link can give probabilities above 1",
               fixed = TRUE)
  expect_match(
    refused_message(glgpm(y ~ cov + gp(), data = poisson_data,
                          family = poisson("identity"), messages = FALSE)),
    "because the identity link can give negative means", fixed = TRUE
  )
  expect_match(
    refused_message(glgpm(y ~ cov + gp(), data = poisson_data,
                          family = poisson("sqrt"), messages = FALSE)),
    "a positive and a negative linear predictor", fixed = TRUE
  )
})

test_that("quasi and other families are refused", {
  for (family in list(quasibinomial(), quasipoisson(), quasi())) {
    expect_match(refused_message(fit_binomial(family)),
                 "quasi-likelihood family, which has no likelihood")
  }
  expect_match(refused_message(fit_binomial(Gamma())),
               "Gaussian, binomial or Poisson family, not \"Gamma\"")
  error <- expect_input_error(fit_binomial(gaussian("log")), "family")
  expect_match(conditionMessage(error),
               "The log link isn't available for the gaussian family")
  other_link <- make.link("inverse")
  expect_match(refused_message(fit_binomial(binomial(link = other_link))),
               "The inverse link isn't available for the binomial family")
  expect_input_error(fit_binomial(1), "family")
})

test_that("a link is given either in family or in invlink", {
  error <- expect_input_error(fit_binomial(binomial("probit"), pnorm),
                              "invlink")
  expect_match(conditionMessage(error), "not both", fixed = TRUE)
  expect_input_error(fit_binomial("binomial", "probit"), "invlink")
  error <- expect_input_error(fit_binomial("binomial", make.link("cloglog")),
                              "invlink")
  expect_match(conditionMessage(error), "binomial(\"cloglog\")", fixed = TRUE)
  # The internal list form of invlink takes one set of names.
  expect_error(fit_binomial("binomial", list(inv_link = plogis)),
               "takes only the elements `inv`, `d1`, `d2` and `link`")
})

test_that("print() and summary() name the link", {
  fit <- named_fits$cloglog
  expect_output(print(fit),
                "Binomial geostatistical model with cloglog link, 10 obs")
  fit_summary <- summary(fit)
  expect_identical(fit_summary$link_name, "cloglog")
  expect_identical(fit_summary$invlink_expression,
                   "Inverse link function = 1 - exp(-exp(x))")
  expect_output(print(fit_summary), "Link: cloglog")
  expect_output(print(named_fits$probit), "with probit link")
  expect_output(print(fit_binomial("binomial", link_functions$cloglog)),
                "with a custom link")
  expect_identical(summary(named_fits$cauchit)$link_name, "cauchit")
  expect_output(print(binomial_model),
                "Binomial geostatistical model, 10 observations")
})

test_that("check_fit() notes a link not known to be log-concave", {
  link_rows <- function(fit) {
    checks <- check_fit(fit)$checks
    checks[checks$check == "Link", ]
  }
  for (fit in list(binomial_model, poisson_model, named_fits$probit,
                   named_fits$cloglog)) {
    expect_identical(nrow(link_rows(fit)), 0L)
  }
  cauchit <- link_rows(named_fits$cauchit)
  expect_identical(cauchit$status, "OK")
  expect_match(cauchit$finding, paste0(
    "^The cauchit link is not among those known to make the likelihood of ",
    "each location log-concave in its latent effect, so the latent effects ",
    "can have more than one conditional mode"
  ))
  expect_identical(cauchit$fix, "")
  # The row is OK, so summary() does not count it among the checks to review.
  expect_identical(summarise_mcml_fit(named_fits$cauchit)$checks_to_review,
                   sum(check_fit(named_fits$cauchit)$checks$status != "OK"))

  supplied <- link_rows(fit_binomial("binomial", link_functions$cloglog))
  expect_match(supplied$finding, paste0(
    "^The link was given as a function, so it is not known whether.*as the ",
    "logit, probit and cloglog links do\\. If it does not, the latent effects"
  ))
  expect_match(fit_check_link(list(name = "custom"), "poisson")$finding,
               "as the log link does\\.")
  # A refit receives the name of its link in the list form, so it is known.
  expect_null(fit_check_link(list(name = "custom", link = "probit"),
                             "binomial"))
})

test_that("a named link starts from a regression with its own link", {
  aux_data <- data.frame(y = binomial_data$y,
                         units_m = binomial_data$denominator,
                         cov = binomial_data$cov)
  # The regression the start was taken from before it regressed on the design
  # matrix, which it reproduces exactly for a model with an intercept.
  regress <- function(glm_family) {
    coef(glm(cbind(y, units_m - y) ~ ., offset = rep(0, nrow(aux_data)),
             data = aux_data, family = glm_family))
  }
  D <- model.matrix(~ cov, data = aux_data)
  start <- function(link) {
    coef(binomial_start_regression(aux_data$y, aux_data$units_m, D,
                                   rep(0, nrow(aux_data)), link))
  }
  for (link in names(link_functions)) {
    expect_identical(start(link), regress(binomial(link)))
  }
  expect_identical(start(NULL), regress(binomial()))
  expect_identical(start("logit"), regress(binomial()))
  # A regression that cannot be fitted with the link falls back on logit.
  expect_identical(start("no such link"), regress(binomial()))
})

test_that("refits, the bootstrap and predict() keep the named link", {
  fit <- named_fits$probit
  carried <- refit_invlink(fit)
  expect_identical(carried$link, "probit")
  refit <- fit_binomial("binomial", carried)
  expect_identical(refit$link_function$link, "probit")
  expect_identical(fit_without_timing(refit), fit_without_timing(fit))
  expect_identical(environment(make_bootstrap_refit(fit))$invlink, carried)
  expect_null(refit_invlink(binomial_model))
  expect_null(refit_invlink(gaussian_model))

  prediction <- predict(fit, keep_samples = TRUE, messages = FALSE)
  expect_equal(prediction$predictions$response_mean,
               rowMeans(pnorm(prediction$samples$spatial +
                                prediction$predictions$fixed)),
               tolerance = 1e-12)
})

test_that("specify_glgpm() reads family objects", {
  parameters <- list(beta = c(0.2, 0.3), sigma2 = 0.5, phi = 0.3)
  specify <- function(family) {
    specify_glgpm(y ~ cov + gp(), binomial_data, family, parameters,
                  denominator = denominator)
  }
  probit <- specify(binomial("probit"))
  expect_identical(probit$invlink$link, "probit")
  expect_true(probit$custom_link)
  expect_output(print(probit), "Family: binomial \\(probit link\\)")
  draws <- simulate_glgpm(probit, nsim = 2, seed = 1)$samples$data
  expect_identical(draws[, , "mean"], pnorm(draws[, , "linear_predictor"]))

  expect_identical(
    simulate_glgpm(specify(binomial()), nsim = 2, seed = 1)$samples,
    simulate_glgpm(specify("binomial"), nsim = 2, seed = 1)$samples
  )
  expect_output(print(specify(binomial())), "Family: binomial\n")

  from_fit <- simulation_model(named_fits$probit)
  expect_identical(from_fit$invlink$link, "probit")
  expect_true(from_fit$custom_link)
  expect_input_error(specify(poisson("identity")), "family")
})

test_that("maxim_integrand() reads family objects", {
  coords <- cbind(c(0, 1, 2), c(0, 0, 1))
  sigma <- exp(-as.matrix(dist(coords)))
  mode <- function(family, invlink = NULL) {
    maxim_integrand(y = c(3, 5, 1), denominator = c(10, 10, 10),
                    mu = rep(-0.5, 3), Sigma = sigma, ID_coords = 1:3,
                    family = family, invlink = invlink)
  }
  expect_identical(mode(binomial()), mode("binomial"))
  probit <- mode(binomial("probit"))
  expect_identical(probit$invlink_used, "custom")
  expect_equal(probit$mode, mode("binomial", function(eta) pnorm(eta))$mode,
               tolerance = 1e-8)
  expect_input_error(mode(gaussian()), "family")
})
