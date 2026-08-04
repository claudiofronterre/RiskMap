test_that("random effects are encoded consistently across data containers", {
  base_data <- data.frame(
    region = factor(
      c("north", "south", "north", "south"),
      levels = c("north", "south", "unused")
    ),
    province = c(10, 30, 10, 20),
    x = seq_len(4),
    y = seq_len(4)
  )
  tibble_data <- tibble::as_tibble(base_data)
  sf_data <- sf::st_as_sf(
    base_data,
    coords = c("x", "y"),
    crs = 32634
  )
  sf_tibble <- sf::st_as_sf(
    tibble_data,
    coords = c("x", "y"),
    crs = 32634
  )

  expected <- RiskMap:::prepare_random_effects(
    base_data,
    c("region", "province")
  )

  expect_equal(
    RiskMap:::prepare_random_effects(tibble_data, c("region", "province")),
    expected
  )
  expect_equal(
    RiskMap:::prepare_random_effects(sf_data, c("region", "province")),
    expected
  )
  expect_equal(
    RiskMap:::prepare_random_effects(sf_tibble, c("region", "province")),
    expected
  )
  expect_equal(
    expected$ID_re,
    matrix(
      c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 3L),
      ncol = 2L,
      dimnames = list(NULL, c("region", "province"))
    )
  )
  expect_equal(expected$re_unique$region, seq_len(3))
  expect_equal(expected$re_unique_f$region, c("north", "south", "unused"))
  expect_equal(expected$re_unique$province, seq_len(3))
  expect_equal(expected$re_unique_f$province, c(10, 30, 20))
})

test_that("supported random-effect types retain their labels and ordering", {
  data <- data.frame(
    ordered_group = ordered(
      c("medium", "low", "high"),
      levels = c("low", "medium", "high")
    ),
    character_group = c("beta", "alpha", "beta"),
    integer_group = c(20L, 10L, 20L)
  )

  result <- RiskMap:::prepare_random_effects(
    data,
    c("ordered_group", "character_group", "integer_group")
  )

  expect_equal(result$ID_re[, "ordered_group"], c(2L, 1L, 3L))
  expect_equal(
    result$re_unique_f$ordered_group,
    c("low", "medium", "high")
  )
  expect_equal(result$ID_re[, "character_group"], c(1L, 2L, 1L))
  expect_equal(result$re_unique_f$character_group, c("beta", "alpha"))
  expect_equal(result$ID_re[, "integer_group"], c(1L, 2L, 1L))
  expect_equal(result$re_unique_f$integer_group, c(20L, 10L))
})

test_that("a single random effect retains a matrix index", {
  result <- RiskMap:::prepare_random_effects(
    data.frame(group = factor(c("a", "b", "a"))),
    "group"
  )

  expect_equal(dim(result$ID_re), c(3L, 1L))
  expect_equal(colnames(result$ID_re), "group")
  expect_equal(result$ID_re[, 1], c(1L, 2L, 1L))
})

test_that("one-row data retain all random-effect columns", {
  result <- RiskMap:::prepare_random_effects(
    data.frame(
      region = factor("north"),
      province = 10
    ),
    c("region", "province")
  )

  expect_equal(dim(result$ID_re), c(1L, 2L))
  expect_equal(colnames(result$ID_re), c("region", "province"))
  expect_equal(result$ID_re[1, ], c(region = 1L, province = 1L))
})

test_that("no random effects return the established null representation", {
  result <- RiskMap:::prepare_random_effects(data.frame(value = 1:3), NULL)

  expect_equal(result$n_re, 0L)
  expect_null(result$names_re)
  expect_null(result$ID_re)
  expect_null(result$re_unique)
  expect_null(result$re_unique_f)
})

test_that("invalid random-effect variables fail informatively", {
  expect_error(
    RiskMap:::prepare_random_effects(
      data.frame(group = c("a", NA)),
      "group"
    ),
    "contains missing values",
    fixed = TRUE
  )
  expect_error(
    RiskMap:::prepare_random_effects(
      data.frame(group = c(1, Inf)),
      "group"
    ),
    "contains non-finite values",
    fixed = TRUE
  )
  expect_error(
    RiskMap:::prepare_random_effects(
      data.frame(group = c(TRUE, FALSE)),
      "group"
    ),
    "must be a factor, character, integer, or numeric vector",
    fixed = TRUE
  )
  expect_error(
    RiskMap:::prepare_random_effects(data.frame(value = 1:2), "group"),
    "not found in `data`",
    fixed = TRUE
  )
})

test_that("glgpm simulations support random effects in tibble-backed sf data", {
  data <- tibble::tibble(
    y = 0,
    group = factor(c("a", "a", "b", "b")),
    x = c(0, 1, 0, 1),
    z = c(0, 0, 1, 1)
  )
  data <- sf::st_as_sf(data, coords = c("x", "z"), crs = 32634)

  set.seed(1)
  result <- suppressWarnings(
    glgpm_sim(
      n_sim = 2,
      formula = y ~ gp(kappa = 0.5, nugget = 0) + re(group),
      data = data,
      family = "gaussian",
      sim_pars = list(
        beta = 0,
        sigma2 = 1,
        tau2 = 0,
        phi = 1,
        sigma2_me = 0.1,
        sigma2_re = 0.2
      ),
      messages = FALSE
    )
  )

  expect_length(result$re_sim, 2L)
  expect_length(result$re_sim[[1]]$group, 2L)
  expect_equal(result$sigma2_re, 0.2)
})

test_that("dast simulations support random effects in tibble-backed sf data", {
  data <- tibble::tibble(
    y = c(1, 2, 1, 2),
    m = rep(10, 4),
    survey_time = seq_len(4),
    group = factor(c("a", "a", "b", "b")),
    x = c(0, 1, 0, 1),
    z = c(0, 0, 1, 1)
  )
  data <- sf::st_as_sf(data, coords = c("x", "z"), crs = 32634)
  m <- data$m

  set.seed(1)
  result <- suppressWarnings(
    dast_sim(
      n_sim = 2,
      formula = y ~ gp(kappa = 0.5, nugget = 0) + re(group),
      data = data,
      den = m,
      time = "survey_time",
      mda_times = 2,
      int_mat = matrix(0, nrow = 4, ncol = 1),
      power_val = 1,
      sim_pars = list(
        beta = 0,
        sigma2 = 1,
        tau2 = 0,
        phi = 1,
        psi = NULL,
        sigma2_re = 0.2,
        alpha = 0.1,
        gamma = 0.2
      ),
      messages = FALSE
    )
  )

  expect_length(result$re_sim, 2L)
  expect_length(result$re_sim[[1]]$group, 2L)
  expect_equal(result$sigma2_re, 0.2)
})
