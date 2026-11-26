# Expect an input error about one argument, without matching the whole
# wording. `arg` is compared with the field `arg` of the condition. When
# `message` is given, the message must contain it, which tells apart several
# checks of the same argument.
expect_input_error <- function(object, arg, message = NULL) {
  condition <- expect_error({{ object }}, class = "riskmap_error_input")
  if (inherits(condition, "riskmap_error_input")) {
    expect_identical(condition$arg, arg)
    if (!is.null(message)) {
      expect_match(conditionMessage(condition), message, fixed = TRUE)
    }
  }
  invisible(condition)
}

# Expect exactly the warnings a call gives, as `expect_warning()` catches only
# the first that matches and lets the others through. Each element of
# `messages` must be contained in one warning, and every warning must contain
# one of them. cli wraps messages at the console width, so runs of whitespace
# are read as one space before matching.
expect_warnings <- function(object, messages) {
  warnings <- capture_warnings({{ object }})
  warnings <- gsub("\\s+", " ", warnings)
  contains <- function(warning, message) grepl(message, warning, fixed = TRUE)
  given <- if (length(warnings) == 0) "none" else paste("*", warnings, collapse = "\n")
  for (message in messages) {
    expect(
      any(contains(warnings, message)),
      sprintf("No warning contained \"%s\". Warnings given:\n%s", message, given)
    )
  }
  unexpected <- warnings[!vapply(warnings, function(warning) {
    any(vapply(messages, contains, logical(1), warning = warning))
  }, logical(1))]
  expect(
    length(unexpected) == 0,
    sprintf("Unexpected warnings:\n%s", paste("*", unexpected, collapse = "\n"))
  )
  invisible(warnings)
}
