# Expect an input error about one argument, without matching the wording.
# `arg` is compared with the field `arg` of the condition.
expect_input_error <- function(object, arg) {
  condition <- expect_error({{ object }}, class = "riskmap_error_input")
  if (inherits(condition, "riskmap_error_input")) {
    expect_identical(condition$arg, arg)
  }
  invisible(condition)
}
