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

table_test_cv <- function() {
  structure(
    matrix(c(0.49, 0.25),
           nrow = 1,
           dimnames = list("model_a", c("crps", "scrps"))),
    class = c("summary.RiskMap_cross_validation", "matrix", "array")
  )
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

test_that("to_table formats cross-validation summaries directly", {
  object <- table_test_cv()
  table <- to_table(object, digits = 3, format = "pipe")
  rendered <- paste(as.character(table), collapse = "\n")

  expect_s3_class(table, "knitr_kable")
  expect_match(rendered, "model_a", fixed = TRUE)
  expect_match(rendered, "0.490", fixed = TRUE)
  expect_match(rendered, "CRPS", fixed = TRUE)
})

test_that("to_table validates its inputs", {
  object <- table_test_model()

  expect_error(to_table(object, digits = -1), "non-negative integer")
  expect_error(to_table(object, digits = 1.5), "non-negative integer")
  expect_error(to_table(1), "RiskMap model")
})
