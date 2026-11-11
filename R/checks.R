# Checks of the arguments of exported functions.
#
# Each check raises an error of class `riskmap_error_input` whose field `arg`
# names the argument, and returns its input invisibly (`match_choice()`
# returns the matched choice). They follow the conventions of rlang's own
# type checks. The argument name defaults to the expression passed as `x`,
# and `call` to the function that called the check. The error names the
# user-facing function the user called (`riskmap_abort()`), so a check made
# inside an internal function still names the exported one.

#' Raise an input error for an argument
#' @param message The message, as in `riskmap_abort()`.
#' @param arg The name of the argument, stored in the field `arg`.
#' @param call Where the error was raised, as in `riskmap_abort()`. A check
#'   helper passes the `call` it was given.
#' @param ... Other fields, such as `parent`.
#' @param .envir The environment in which the message is interpolated.
#' @noRd
abort_input <- function(message, arg, call = rlang::caller_env(), ...,
                        .envir = parent.frame()) {
  riskmap_abort(message, ..., arg = arg, class = "riskmap_error_input",
                call = call, .envir = .envir)
}

#' Describe a supplied value for an error message
#'
#' A single atomic value is shown as itself, for example `0.5`, `NA` or
#' `"yes"`. Anything else is described by its type, for example "a double
#' vector" or "a data frame". Used after "not" in input errors.
#' @param x The value supplied by the user.
#' @return A character string, already formatted by cli.
#' @noRd
describe_value <- function(x) {
  if (is.factor(x) && length(x) == 1L) x <- as.character(x)
  if (is.atomic(x) && length(x) == 1L) {
    cli::format_inline("{.val {x}}")
  } else {
    cli::format_inline("{.obj_type_friendly {x}}")
  }
}

#' Check that a value is a single positive integer
#' @param x The value to check.
#' @param arg The name of the argument, used in the message and stored in the
#'   field `arg`.
#' @param allow_null Whether `NULL` is accepted.
#' @param allow_zero Whether zero is accepted.
#' @param call The call named in the error.
#' @noRd
check_positive_integer <- function(x, arg = rlang::caller_arg(x),
                                   allow_null = FALSE, allow_zero = FALSE,
                                   call = rlang::caller_env()) {
  if (is.null(x) && allow_null) return(invisible(x))
  description <- if (allow_zero) "non-negative" else "positive"
  invalid <- !is.numeric(x) || length(x) != 1L || is.na(x) ||
    !is.finite(x) || x %% 1 != 0 || x < as.integer(!allow_zero) ||
    x > .Machine$integer.max
  if (invalid) {
    supplied <- describe_value(x)
    or_null <- if (allow_null) " or {.code NULL}" else ""
    abort_input(
      paste0("{.arg {arg}} must be a single {description} integer", or_null, ", not {supplied}."),
      arg = arg, call = call
    )
  }
  invisible(x)
}

#' Check that a value is a single positive number
#'
#' @param x The value to check.
#' @param type The kind of value. The default, `"starting "`, is for an
#'   element of a list of starting values such as `start_pars[["sigma2"]]`,
#'   and the message names the element. Use `""` for an argument.
#' @param arg The expression supplied as `x`. The field `arg` stores the
#'   argument it comes from, `start_pars` for `start_pars[["sigma2"]]`.
#' @param call The call named in the error.
#' @noRd
check_positive_number <- function(x, type = "starting ",
                                  arg = rlang::caller_arg(x),
                                  call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    # The element of a list, as in start_pars[["sigma2"]], or the argument
    name <- gsub('.*\\[\\["([^"]+)"\\]\\].*', "\\1", arg)
    # The argument that holds the value, such as `start_pars`.
    argument <- sub("\\[.*", "", arg)
    supplied <- describe_value(x)
    problem <- if (nzchar(type)) {
      "The {trimws(type)} value for {.code {name}} must be a single positive number, not {supplied}."
    } else {
      "{.arg {name}} must be a single positive number, not {supplied}."
    }
    abort_input(problem, arg = argument, call = call)
  }
  invisible(x)
}

#' Check that a value is a single number in a range
#'
#' @param x The value to check.
#' @param min,max The endpoints. An infinite endpoint is left out of the
#'   message.
#' @param allow_equal Whether the endpoints are accepted.
#' @param arg The name of the argument.
#' @param call The call named in the error.
#' @noRd
check_range <- function(x, min = -Inf, max = Inf, allow_equal = TRUE,
                        arg = rlang::caller_arg(x),
                        call = rlang::caller_env()) {
  valid_number <- is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x)
  outside <- if (!valid_number) {
    TRUE
  } else if (allow_equal) {
    x < min || x > max
  } else {
    x <= min || x >= max
  }
  if (!valid_number || outside) {
    description <- describe_range(min, max, allow_equal)
    supplied <- describe_value(x)
    abort_input(
      "{.arg {arg}} must be a single {description}, not {supplied}.",
      arg = arg, call = call
    )
  }
  invisible(x)
}

#' Describe the numbers that `check_range()` accepts
#'
#' An infinite bound is left out, so that `min = 0, max = Inf` gives
#' "positive finite number" rather than "finite number strictly between 0 and
#' Inf".
#' @param min,max,allow_equal As in `check_range()`.
#' @return A character string to follow "must be a single".
#' @noRd
describe_range <- function(min, max, allow_equal) {
  has_lower <- is.finite(min)
  has_upper <- is.finite(max)
  if (has_lower && has_upper) {
    interval <- if (allow_equal) "between" else "strictly between"
    return(paste("finite number", interval, min, "and", max))
  }
  if (has_lower && min == 0) {
    return(if (allow_equal) "non-negative finite number" else "positive finite number")
  }
  if (has_upper && max == 0) {
    return(if (allow_equal) "non-positive finite number" else "negative finite number")
  }
  if (has_lower) {
    comparison <- if (allow_equal) "greater than or equal to" else "greater than"
    return(paste("finite number", comparison, min))
  }
  if (has_upper) {
    comparison <- if (allow_equal) "less than or equal to" else "less than"
    return(paste("finite number", comparison, max))
  }
  "finite number"
}

#' Check that a value is a single number between 0 and 1, inclusive
#' @inheritParams check_range
#' @noRd
check_zero_one <- function(x, arg = rlang::caller_arg(x),
                           call = rlang::caller_env()) {
  check_range(x, min = 0, max = 1, arg = arg, call = call)
}

#' Check that a value is `TRUE` or `FALSE`
#' @inheritParams check_range
#' @noRd
check_logical <- function(x, arg = rlang::caller_arg(x),
                          call = rlang::caller_env()) {
  if (!isTRUE(x) && !isFALSE(x)) {
    supplied <- describe_value(x)
    abort_input(
      "{.arg {arg}} must be {.code TRUE} or {.code FALSE}, not {supplied}.",
      arg = arg, call = call
    )
  }
  invisible(x)
}

#' Check that a value inherits from one of some classes
#'
#' @param x The value to check.
#' @param class The accepted classes. A value that inherits from any of them
#'   passes.
#' @param allow_null Whether `NULL` is accepted.
#' @inheritParams check_range
#' @noRd
check_class <- function(x, class, arg = rlang::caller_arg(x),
                        allow_null = FALSE, call = rlang::caller_env()) {
  if (is.null(x) && allow_null) return(invisible(x))
  if (!inherits(x, class)) {
    # One {.cls} per class, because {.cls} joins a vector into one hierarchy
    classes <- vapply(class, function(name) cli::format_inline("{.cls {name}}"),
                      character(1), USE.NAMES = FALSE)
    abort_input(
      "{.arg {arg}} must be an object of class {.or {classes}}, not {.obj_type_friendly {x}}.",
      arg = arg, call = call
    )
  }
  invisible(x)
}

#' Check that a control argument was made by its function
#'
#' For `control_mcmc` and `control_mcml`, made by [set_control_mcmc()] and
#' [set_control_mcml()].
#' @param maker The name of the function that makes it.
#' @noRd
check_control <- function(x, maker = c("set_control_mcmc", "set_control_mcml"),
                          arg = rlang::caller_arg(x),
                          call = rlang::caller_env()) {
  maker <- match.arg(maker)
  class <- switch(maker, set_control_mcmc = "RiskMap_control_mcmc",
                  set_control_mcml = "RiskMap_control_mcml")
  if (!inherits(x, class)) {
    abort_input(
      "{.arg {arg}} must be made by {.fn {maker}}, not {.obj_type_friendly {x}}.",
      arg = arg, call = call
    )
  }
  invisible(x)
}

#' Match a value to one of some choices
#'
#' A replacement for [match.arg()] that raises a RiskMap input error. As with
#' `match.arg()`, the choices default to the default value of the argument in
#' the calling function, the first choice is returned when the argument was
#' not supplied, and a unique abbreviation is accepted.
#'
#' @param x The value to match.
#' @param choices The accepted values. By default, the default value of the
#'   argument `arg` in the calling function.
#' @inheritParams check_range
#' @return The matched choice.
#' @noRd
match_choice <- function(x, choices = NULL, arg = rlang::caller_arg(x),
                         call = rlang::caller_env()) {
  if (is.null(choices)) {
    caller_formals <- formals(sys.function(sys.parent()))
    # An argument without a default has the empty symbol as its formal, which
    # can't be stored in a variable, so it is tested where it is extracted.
    if (!arg %in% names(caller_formals) || !is.call(caller_formals[[arg]])) {
      riskmap_bug("{.fn match_choice} needs {.arg choices} for {.code {arg}}.")
    }
    choices <- eval(caller_formals[[arg]], envir = call)
  }
  # NULL or the whole default give the first choice, as with match.arg().
  if (is.null(x) || identical(x, choices)) return(choices[[1L]])
  index <- if (is.character(x) && length(x) == 1L && !is.na(x)) {
    pmatch(x, choices)
  } else {
    NA_integer_
  }
  if (is.na(index)) {
    supplied <- describe_value(x)
    abort_input(
      "{.arg {arg}} must be one of {.or {.val {choices}}}, not {supplied}.",
      arg = arg, call = call
    )
  }
  choices[[index]]
}

#' Check that a value is a valid coordinate reference system
#'
#' @param crs The value to check, as accepted by [sf::st_crs()].
#' @inheritParams check_range
#' @return `crs` invisibly. The error has the warning or error from
#'   [sf::st_crs()] as its parent.
#' @noRd
check_crs <- function(crs, arg = rlang::caller_arg(crs),
                      call = rlang::caller_env()) {
  invalid_crs <- function(condition) {
    abort_input(
      "{.arg {arg}} must be a valid coordinate reference system.",
      arg = arg, call = call, parent = condition
    )
  }
  tryCatch(st_crs(crs), warning = invalid_crs, error = invalid_crs)
  invisible(crs)
}
