#' Conditions raised by RiskMap
#'
#' @description
#' Errors, warnings and messages from RiskMap carry classes, so that code can
#' handle them with [tryCatch()] or [withCallingHandlers()] without matching
#' their text, which may be reworded between versions.
#'
#' - `riskmap_error` is the class of every error that RiskMap raises on
#'   purpose. A more specific class says what went wrong.
#'   - `riskmap_error_input`, an argument has a value that the function can't
#'     use. The field `arg` of the condition names the argument.
#'   - `riskmap_error_worldpop`, WorldPop can't be reached or has no data for
#'     the request.
#'   - `riskmap_error_internal`, a problem in RiskMap itself. Please report it.
#' - `riskmap_warning` is the class of every warning.
#' - `riskmap_message` is the class of every message. Functions that report
#'   on their progress take an argument `messages`, whose default is the
#'   option `RiskMap.messages`, so that `options(RiskMap.messages = FALSE)`
#'   silences them all.
#'
#' RiskMap is moving to these classes one source file at a time, so some of
#' its conditions don't carry them yet.
#'
#' @examples
#' # Which argument was rejected
#' tryCatch(
#'   get_worldpop("Liberia", 2020),
#'   riskmap_error_input = function(condition) condition$arg
#' )
#' @name riskmap_conditions
NULL

#' Raise an error, warning or message
#'
#' Wrappers of [cli::cli_abort()], [cli::cli_warn()] and [cli::cli_inform()]
#' that add the classes `riskmap_error`, `riskmap_warning` and
#' `riskmap_message` after any class given. Every condition in RiskMap is
#' raised through these, so that a change to all of them is made here.
#' @param message The message, with cli inline markup and bullets. It is
#'   interpolated in `.envir`.
#' @param ... Fields stored in the condition, such as `arg` or `parent`.
#' @param class More specific classes, placed before the RiskMap class.
#' @param call The call named in the error. The default is the function that
#'   called `riskmap_abort()`. A check helper passes the `call` it was given.
#' @param .envir The environment in which the message is interpolated.
#' @noRd
riskmap_abort <- function(message, ..., class = NULL,
                          call = rlang::caller_env(),
                          .envir = parent.frame()) {
  cli::cli_abort(message, ..., class = c(class, "riskmap_error"),
                 call = call, .envir = .envir)
}

#' @rdname riskmap_abort
#' @noRd
riskmap_warn <- function(message, ..., class = NULL,
                         .envir = parent.frame()) {
  cli::cli_warn(message, ..., class = c(class, "riskmap_warning"),
                .envir = .envir)
}

#' @rdname riskmap_abort
#' @noRd
riskmap_inform <- function(message, ..., class = NULL,
                           .envir = parent.frame()) {
  cli::cli_inform(message, ..., class = c(class, "riskmap_message"),
                  .envir = .envir)
}

#' Raise an error for a state that RiskMap should never reach
#'
#' For checks of RiskMap's own logic, not of what the user supplied. The error
#' names the internal function and asks for a report.
#' @inheritParams riskmap_abort
#' @noRd
riskmap_bug <- function(message, ..., call = rlang::caller_env(),
                        .envir = parent.frame()) {
  riskmap_abort(
    c(message,
      "i" = "This is a problem in RiskMap. Please report it at {.url https://github.com/claudiofronterre/RiskMap/issues}."),
    ...,
    class = "riskmap_error_internal",
    call = call,
    .envir = .envir
  )
}
