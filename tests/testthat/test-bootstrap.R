test_that("a refit to the observed response is a glgpm() fit from the estimate", {
  refit <- make_bootstrap_refit(binomial_model)
  refitted <- suppressWarnings(refit(binomial_model$y, 77L))
  direct <- suppressWarnings(glgpm(
    y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
    denominator = denominator, control_mcml = fixture_control_mcml,
    start_pars = coef(binomial_model), seed = 77L, return_samples = FALSE,
    messages = FALSE
  ))
  expect_identical(refitted, unlist(coef(direct)))

  # The offset, the nugget and the measurement error variance are carried
  # over from the fit.
  gaussian_refit <- make_bootstrap_refit(gaussian_offset_model)
  gaussian_direct <- glgpm(
    y ~ cov + gp(nugget = TRUE) + offset(offset), data = gaussian_data,
    family = "gaussian", fix_var_me = 0,
    start_pars = coef(gaussian_offset_model), messages = FALSE
  )
  expect_identical(gaussian_refit(gaussian_offset_model$y, 1L),
                   unlist(coef(gaussian_direct)))
})

test_that("replicates have their own seeds and do not depend on cores", {
  set.seed(4)
  caller_state <- .Random.seed
  first <- parametric_bootstrap(gaussian_model, n_sim = 4, seed = 11,
                                messages = FALSE)
  expect_identical(.Random.seed, caller_state)
  expect_s3_class(first, "RiskMap_bootstrap")
  expect_identical(dim(first$estimates), c(4L, 6L))
  expect_identical(colnames(first$estimates),
                   names(unlist(coef(gaussian_model))))
  expect_identical(first$estimate, unlist(coef(gaussian_model)))
  expect_identical(first$settings$seed, 11L)
  expect_identical(bootstrap_seeds(11L, 3L, 10000L),
                   c(11L, 10011L, 20011L))
  expect_identical(bootstrap_seed_spacing(1333L), 10000L)
  expect_error(bootstrap_seed_spacing(1334L), "at most 1333 stages")

  smaller <- parametric_bootstrap(gaussian_model, n_sim = 2, seed = 11,
                                  messages = FALSE)
  expect_identical(smaller$estimates, first$estimates[1:2, ])

  # The first replicate simulates and refits with the seed itself.
  simulated <- simulate_glgpm(gaussian_model, seed = 11)
  refitted <- suppressWarnings(make_bootstrap_refit(gaussian_model)(
    simulated$samples$data[, 1L, "response"], 11L
  ))
  expect_identical(first$estimates[1L, ], refitted)

  skip_on_os("windows")
  parallel <- parametric_bootstrap(gaussian_model, n_sim = 4, seed = 11,
                                   cores = 2, messages = FALSE)
  expect_identical(parallel$estimates, first$estimates)
})

test_that("without a seed, one is drawn from the session's stream", {
  set.seed(5)
  expected_seed <- sample.int(.Machine$integer.max, 1L)
  expected_state <- .Random.seed
  set.seed(5)
  drawn <- parametric_bootstrap(gaussian_intercept_model, n_sim = 2,
                                messages = FALSE)
  expect_identical(drawn$settings$seed, expected_seed)
  expect_identical(.Random.seed, expected_state)
})

test_that("binomial bootstraps with random effects use separate streams", {
  boot <- suppressWarnings(parametric_bootstrap(binomial_model, n_sim = 2,
                                                seed = 3, messages = FALSE))
  expect_identical(colnames(boot$estimates),
                   names(unlist(coef(binomial_model))))
  expect_false(anyNA(boot$estimates))
  expect_identical(boot$settings$seed_spacing, 10000L)
  expect_contains(names(boot$warnings), c("replicate", "seed", "message"))
})

test_that("failed refits are recorded instead of stopping", {
  original_glgpm <- glgpm
  local_mocked_bindings(glgpm = function(..., seed) {
    if (identical(seed, 10001L)) stop("simulated failure")
    original_glgpm(..., seed = seed)
  })
  expect_warning(
    boot <- parametric_bootstrap(gaussian_intercept_model, n_sim = 3,
                                 seed = 1, messages = FALSE),
    "1 of 3 refits failed"
  )
  expect_identical(boot$failures,
                   data.frame(replicate = 2L, seed = 10001L,
                              message = "simulated failure"))
  expect_true(all(is.na(boot$estimates[2L, ])))
  expect_false(anyNA(boot$estimates[-2L, ]))
  expect_output(print(boot), "1 refits failed")

  local_mocked_bindings(glgpm = function(...) stop("always fails"))
  expect_error(
    parametric_bootstrap(gaussian_intercept_model, n_sim = 2, seed = 1,
                         messages = FALSE),
    "All 2 refits failed. The first failure: always fails"
  )
})

test_that("summary computes percentile and bias-corrected intervals", {
  values <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) / 10
  boot <- structure(list(
    estimates = cbind(phi = values, sigma2 = values),
    estimate = c(phi = 0.75, sigma2 = 0.05),
    labels = c(phi = "Scale", sigma2 = "Variance"),
    failures = data.frame(replicate = integer(), seed = integer(),
                          message = character()),
    warnings = data.frame(),
    settings = list(n_sim = 10L, seed = 1L, level = 0.9,
                    family = "binomial")
  ), class = "RiskMap_bootstrap")
  table <- summary(boot)$table
  phi <- table[table$parameter == "phi", ]
  expect_equal(phi$median, median(values))
  expect_equal(phi$bias, mean(values) - 0.75)
  expect_equal(phi$se, sd(values))
  expect_equal(phi$below, 0.7)
  expect_equal(c(phi$percentile_lower, phi$percentile_upper),
               quantile(values, c(0.05, 0.95), names = FALSE))
  z0 <- qnorm(0.7)
  expect_equal(c(phi$bc_lower, phi$bc_upper),
               quantile(values, pnorm(2 * z0 + qnorm(c(0.05, 0.95))),
                        names = FALSE))
  # Most estimates below the estimate move the BC interval upward.
  expect_gt(phi$bc_lower, phi$percentile_lower)
  expect_gt(phi$bc_upper, phi$percentile_upper)
  # Every estimate above the estimate leaves the BC interval undefined.
  sigma2 <- table[table$parameter == "sigma2", ]
  expect_true(is.na(sigma2$bc_lower) && is.na(sigma2$bc_upper))
  expect_output(print(boot), "BC interval\\s+is not available")
  expect_equal(summary(boot, level = 0.5)$level, 0.5)
})

test_that("plot draws one histogram per parameter", {
  boot <- parametric_bootstrap(gaussian_intercept_model, n_sim = 5, seed = 2,
                               messages = FALSE)
  single <- plot(boot, parameters = "phi")
  expect_s3_class(single, "ggplot")
  expect_identical(single$labels$title,
                   paste0("Spatial corr. scale (",
                          gaussian_intercept_model$distance_units, ")"))
  skip_if_not_installed("patchwork")
  expect_s3_class(plot(boot), "patchwork")
  expect_error(plot(boot, parameters = "kappa"), "Unknown parameter")
})

test_that("parametric_bootstrap validates its arguments", {
  expect_error(parametric_bootstrap(list()), "'RiskMap' object")
  expect_error(parametric_bootstrap(gaussian_model, n_sim = 0),
               "'n_sim' must be a single positive integer")
  expect_error(parametric_bootstrap(gaussian_model, seed = -1),
               "'seed' must be a single non-negative integer")
  expect_error(parametric_bootstrap(gaussian_model, level = 1), "level")
  expect_error(parametric_bootstrap(gaussian_model, cores = 0),
               "'cores' must be a single positive integer")
  expect_error(parametric_bootstrap(gaussian_model, messages = NA),
               "'messages' must be either TRUE or FALSE")
  expect_error(parametric_bootstrap(gaussian_model, n_sim = 300000),
               "'n_sim' must be at most")
})
