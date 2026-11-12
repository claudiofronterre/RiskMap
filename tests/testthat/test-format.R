checks_to_print <- function() {
  structure(list(
    checks = data.frame(
      check = c("Variances", "Monte Carlo accuracy", "Importance weights"),
      status = c("SUSPICIOUS", "OK", "OK"),
      finding = c(paste("The nugget variance (tau2) was estimated as",
                        "effectively zero, so its standard error is not",
                        "reported."),
                  "The error of the estimate is within 0.059 standard errors.",
                  "The certification draws had 1,081 effective draws."),
      fix = c("Consider refitting without it, `gp(nugget = FALSE)`.", "", ""),
      brief = c("", "Within 0.059 standard errors.",
                "1,081 effective draws, tail shape 0.46."),
      stringsAsFactors = FALSE
    ),
    draws = paste("Fitted in 2 stages, with 2,000 fitting and 5,000",
                  "validation draws at the last stage."),
    not_computed = NULL
  ), class = "RiskMap_check")
}

test_that("printed output is not a message, so hiding messages keeps it", {
  expect_no_message(printed <- capture.output(print(binomial_model)))
  hidden <- suppressMessages(capture.output(print(binomial_model)))
  expect_identical(hidden, printed)
  expect_gt(length(printed), 10L)
})

test_that("wrapping never breaks code between backticks", {
  lines <- out_wrap("Consider refitting without it, `gp(nugget = FALSE)`.",
                    width = 40L)
  expect_identical(as.character(lines),
                   c("Consider refitting without it,", "`gp(nugget = FALSE)`."))
  expect_identical(as.character(out_wrap("plain text", width = 40L)),
                   "plain text")
})

test_that("bullets carry cli's symbols and wrap under their text", {
  withr::local_options(cli.width = 40L)
  lines <- out_bullets(c(
    v = "Done.",
    "!" = "A sentence long enough to need a second line under the first."
  ))
  expect_identical(lines[1L], paste(cli::symbol$tick, "Done."))
  expect_true(startsWith(lines[2L], "! A sentence"))
  expect_true(startsWith(lines[3L], "  "))
  expect_lte(max(nchar(lines)), 40L)
  expect_identical(out_bullets(character()), character())
  expect_identical(out_as_bullets(NULL, "i"), character())
})

test_that("a function is linked to its help page only where links work", {
  expect_identical(out_fun("check_fit"), "`check_fit()`")
  withr::local_options(cli.hyperlink = TRUE, cli.hyperlink_help = TRUE)
  linked <- out_fun("summary", "summary.RiskMap")
  expect_match(linked, "x-r-help:RiskMap::summary.RiskMap", fixed = TRUE)
  expect_match(cli::ansi_strip(linked), "summary()", fixed = TRUE)
})

test_that("the checks print with their symbols, in Unicode and in ASCII", {
  checks <- checks_to_print()
  expect_snapshot(print(checks))
  testthat::local_reproducible_output(unicode = TRUE)
  expect_snapshot(print(checks))
})
