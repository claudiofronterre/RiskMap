table_test_model <- function() {
  structure(
    list(
      reg_coef = matrix(
        c(0.49, 0.4, 0.6, 0.01,
          -1.25, -1.5, -1.0, 0.02),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("(Intercept)", "x"),
                        c("Estimate", "Lower limit", "Upper limit", "p-value"))
      ),
      sp = matrix(c(2, 1.5, 2.5),
                  nrow = 1,
                  dimnames = list("Spatial process var.",
                                  c("Estimate", "Lower limit", "Upper limit"))),
      ranef = NULL,
      me = NULL
    ),
    class = "summary.RiskMap"
  )
}

# The summary of a cross-validation, with a difference from a reference.
table_test_cv <- function(crps = c(0.49, 0.25), reference = NULL) {
  out <- data.frame(CRPS = crps, SCRPS = c(0.25, 0.3), PIT_area = c(0.05, 0.04),
                    PIT_area_expected = c(0.03, 0.035),
                    PIT_area_limit = c(0.06, 0.07),
                    row.names = c("model_a", "model_b"))
  if (!is.null(reference)) {
    out$CRPS_difference <- c(0, -0.24)
    out$CRPS_se <- c(0, 0.05)
    out$SCRPS_difference <- c(0, 0.05)
    out$SCRPS_se <- c(0, 0.01)
  }
  structure(out, class = c("summary.RiskMap_cross_validation", "data.frame"),
            reference = reference,
            n_observations = c(model_a = 20, model_b = 20))
}

table_test_rows <- function(...) {
  rows <- as.character(to_table(..., format = "pipe"))
  lapply(strsplit(rows[grepl("^\\|model", rows)], "|", fixed = TRUE),
         function(cells) trimws(cells[-1]))
}

test_that("to_table formats model summaries directly", {
  object <- table_test_model()
  table <- to_table(object, digits = 3, format = "pipe")
  rendered <- paste(as.character(table), collapse = "\n")

  expect_s3_class(table, "knitr_kable")
  expect_match(rendered, "0.490", fixed = TRUE)
  expect_match(rendered, "Spatial process var.", fixed = TRUE)
  expect_match(rendered, "Lower limit", fixed = TRUE)
})

test_that("to_table names the intercept and adds no row for a held parameter", {
  # A fit with only an intercept, and one with a held measurement error, whose summary
  # holds the fixed measurement error apart from the table (#253, step 2).
  rendered <- paste(as.character(to_table(gaussian_intercept_model,
                                          format = "pipe")), collapse = "\n")
  expect_match(rendered, "(Intercept)", fixed = TRUE)
  table <- to_table(gaussian_offset_model, format = "pipe")
  rows <- as.character(table)[-(1:2)]
  expect_length(rows, nrow(summary(gaussian_offset_model)$reg_coef) +
                  nrow(summary(gaussian_offset_model)$sp))
  expect_false(any(grepl("Measurement error", rows)))
})

test_that("to_table formats cross-validation summaries directly", {
  table <- to_table(table_test_cv(), format = "pipe")
  expect_s3_class(table, "knitr_kable")
  expect_match(as.character(table)[1], "Model +\\| +CRPS\\| +SCRPS\\| +PIT area\\| +Limit")
  expect_identical(table_test_rows(table_test_cv()),
                   list(c("model_a", "0.490", "0.250", "0.050", "0.060"),
                        c("model_b", "0.250", "0.300", "0.040", "0.070")))
  compared <- table_test_cv(reference = "model_a")
  expect_identical(table_test_rows(compared, digits = 2),
                   list(c("model_a", "0.49", "", "0.25", "", "0.05", "0.06"),
                        c("model_b", "0.25", "-0.24 (0.05)", "0.30",
                          "+0.05 (0.01)", "0.04", "0.07")))
  expect_input_error(to_table(table_test_cv(), reference = "model_a"),
                     "reference", "summary()")
})

# The cells of the row with this name, after the name itself
table_test_cells <- function(object, name, digits = 3) {
  rows <- as.character(to_table(object, digits = digits, format = "pipe"))
  cells <- strsplit(rows[startsWith(rows, paste0("|", name))], "|",
                    fixed = TRUE)[[1]]
  trimws(cells[-(1:2)])
}

# Those of the row of `x` with this estimate and limits
table_test_row <- function(values, digits = 3) {
  object <- table_test_model()
  object$reg_coef[2, 1:3] <- values
  table_test_cells(object, "x ", digits)
}

test_that("to_table uses scientific notation only when fixed notation is unsuitable", {
  expect_identical(table_test_row(c(1e100, 1.5, 2)),
                   c("1.000e+100", "1.500e+00", "2.000e+00"))
  expect_identical(table_test_row(c(1e6, 1.5, 2)),
                   c("1.000e+06", "1.500e+00", "2.000e+00"))
  expect_identical(table_test_row(c(999999, 1.5, 2)),
                   c("999999.000", "1.500", "2.000"))
  # 1.993e-08 takes nine decimals and 2.5e-06 seven, but 2.5e-05 six
  expect_identical(table_test_row(c(1.993e-8, 0, 0.5)),
                   c("1.993e-08", "0.000e+00", "5.000e-01"))
  expect_identical(table_test_row(c(2.5e-6, 0.1, 0.5)),
                   c("2.500e-06", "1.000e-01", "5.000e-01"))
  expect_identical(table_test_row(c(2.5e-5, 0.1, 0.5)),
                   c("0.000025", "0.100000", "0.500000"))
  # The limit doesn't apply to the decimals that are asked for
  expect_identical(table_test_row(c(0.25, 1e-8, 0.5), digits = 8),
                   c("0.25000000", "0.00000001", "0.50000000"))
  expect_identical(table_test_row(c(0.25, 1e-9, 0.5), digits = 8),
                   c("2.50000000e-01", "1.00000000e-09", "5.00000000e-01"))
  # A value of a million or more decides, whatever the smallest is
  expect_identical(table_test_row(c(2e6, 0.0004, 0.5)),
                   c("2.000e+06", "4.000e-04", "5.000e-01"))
})

test_that("to_table gives small values two significant digits", {
  expect_identical(table_test_row(c(-0.002243, -0.004232, -0.0002506)),
                   c("-0.00224", "-0.00423", "-0.00025"))
  # All the values keep the decimals of the smallest, which isn't a zero
  expect_identical(table_test_row(c(0, -0.00042, 12.3)),
                   c("0.00000", "-0.00042", "12.30000"))
  expect_identical(table_test_row(c(0.00042, NA, Inf)),
                   c("0.00042", "NA", "Inf"))
  # A value that rounds up across a power of ten keeps two digits
  expect_identical(table_test_row(c(0.000996, 0.5, 0.75)),
                   c("0.0010", "0.5000", "0.7500"))
  expect_identical(table_test_row(c(0.0000996, 0.5, 0.75)),
                   c("0.00010", "0.50000", "0.75000"))
  expect_identical(table_test_row(c(9.96e-6, 0.5, 0.75)),
                   c("0.000010", "0.500000", "0.750000"))
  # The smallest value that `digits` decimals show with a digit that isn't 0
  expect_identical(table_test_row(c(0.001, 0.5, 0.75)),
                   c("0.001", "0.500", "0.750"))
  expect_identical(table_test_row(c(0.49, 2, 3), digits = 0),
                   c("0.49", "2.00", "3.00"))
  # Zeros and missing values change nothing by themselves
  expect_identical(table_test_row(c(0, NA, 0.5)), c("0.000", "NA", "0.500"))
  expect_identical(table_test_row(c(NA_real_, NA, NA)), c("NA", "NA", "NA"))
})

test_that("to_table changes only the rows that need it", {
  object <- table_test_model()
  for (small in c(-0.0002506, 1.993e-8, 1e100)) {
    object$reg_coef[2, 1] <- small
    expect_identical(table_test_cells(object, "(Intercept)"),
                     c("0.490", "0.400", "0.600"))
    expect_identical(table_test_cells(object, "Spatial process var."),
                     c("2.000", "1.500", "2.500"))
  }
})

test_that("to_table shows one metric of a cross-validation in one notation", {
  expect_identical(
    vapply(table_test_rows(table_test_cv(crps = c(0.49, 0.0004))), `[[`,
           character(1), 2L),
    c("0.49000", "0.00040")
  )
  expect_identical(
    vapply(table_test_rows(table_test_cv(crps = c(0.5, 1e-8))), `[[`,
           character(1), 2L),
    c("5.000e-01", "1.000e-08")
  )
})

test_that("to_table validates its inputs", {
  object <- table_test_model()

  expect_error(to_table(object, digits = -1), "non-negative integer")
  expect_error(to_table(object, digits = 1.5), "non-negative integer")
  expect_error(to_table(1), "<RiskMap> fit")
})
