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
