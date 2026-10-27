fake_areal_prediction <- function(weights = "people") {
  boundaries <- st_sf(
    name = c("North", "South"),
    geometry = st_sfc(
      st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))),
      st_polygon(list(rbind(c(0, 1), c(1, 1), c(1, 2), c(0, 2), c(0, 1)))),
      crs = 32631
    )
  )
  rows <- function(aggregate, summary, values, threshold = NA) {
    data.frame(boundary_id = 1:2, aggregate = aggregate,
               quantity = "response", summary = summary,
               threshold = threshold, value = values)
  }
  summary <- rbind(
    rows("mean", "mean", c(0.10, 0.30)),
    rows("mean", "lower", c(0.05, 0.20)),
    rows("mean", "upper", c(0.15, 0.40)),
    rows("mean", "exceedance", c(0.02, 0.90), threshold = 0.2),
    rows("sum", "mean", c(1234.5, 98765)),
    rows("sum", "lower", c(1000, 80000)),
    rows("sum", "upper", c(1500, 120000)),
    rows("sum", "exceedance", c(1, 1), threshold = 0.2)
  )
  structure(
    list(predictions = NULL, boundaries = boundaries,
         boundary_summary = summary, family = "binomial", level = 0.95,
         weights = weights),
    class = "RiskMap_prediction"
  )
}

table_lines <- function(table) as.character(table)

test_that("to_table formats polygon predictions", {
  table <- table_lines(to_table(fake_areal_prediction()))
  expect_match(table[1], "Area")
  expect_match(table[1], "Weighted mean prevalence (%, 95% interval)",
               fixed = TRUE)
  expect_match(table[1], "Probability weighted mean prevalence > 20% (%)",
               fixed = TRUE)
  expect_match(table[1], "Weighted sum of prevalence (95% interval)",
               fixed = TRUE)
  # Exceedance thresholds on prevalence do not apply to sums
  expect_false(any(grepl("sum of prevalence >", table, fixed = TRUE)))
  expect_match(table[3], "North")
  expect_match(table[3], "10.0 (5.0 to 15.0)", fixed = TRUE)
  expect_match(table[3], "1,230 (1,000 to 1,500)", fixed = TRUE)
  expect_match(table[4], "98,800 (80,000 to 120,000)", fixed = TRUE)
  expect_match(table[4], "90.0", fixed = TRUE)
})

test_that("to_table selects, sorts and renames columns", {
  prediction <- fake_areal_prediction(weights = NULL)
  sorted <- table_lines(to_table(prediction, sort_by = "mean"))
  expect_match(sorted[3], "South")
  expect_match(sorted[1], "Mean prevalence (%, 95% interval)", fixed = TRUE)

  means <- table_lines(to_table(prediction, aggregate = "mean",
                                summaries = "mean", digits = 0,
                                col_names = c("Area", "Prevalence (%)")))
  expect_match(means[1], "Prevalence (%)", fixed = TRUE)
  expect_match(means[3], "\\| *10\\|")

  expect_error(to_table(prediction, col_names = "Area"), "must have 4")
  expect_error(to_table(prediction, aggregate = "max"), "not calculated")
  expect_error(to_table(prediction, sort_by = "median"), "sort_by")
  expect_error(to_table(prediction, label = "missing"), "label")
  prediction$boundaries <- NULL
  expect_error(to_table(prediction), "boundaries")
})

test_that("to_table keeps working for fitted models, with col_names", {
  table <- table_lines(to_table(gaussian_intercept_model, digits = 2))
  expect_true(length(table) > 2)
  n_columns <- ncol(summary(gaussian_intercept_model)$reg_coef[, 1:3])
  renamed <- table_lines(to_table(gaussian_intercept_model,
                                  col_names = c("Est", "Low", "High")))
  expect_match(renamed[1], "Est")
})
