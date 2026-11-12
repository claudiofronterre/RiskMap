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
#'   - `riskmap_error_fit`, a model can't be fitted or a computation on a
#'     fit fails, for example because the search for the mode of the latent
#'     effects fails or every Monte Carlo weight is zero.
#'   - `riskmap_error_worldpop`, WorldPop can't be reached or has no data for
#'     the request.
#'   - `riskmap_error_internal`, a problem in RiskMap itself. Please report it.
#' - `riskmap_warning` is the class of every warning.
#'   - `riskmap_warning_fit`, a fit may not be reliable, for example because
#'     the optimiser did not converge or the Monte Carlo error of the estimate
#'     is not certified.
#' - `riskmap_message` is the class of every message. Functions that report
#'   on their progress take an argument `messages`, whose default is the
#'   option `RiskMap.messages`, so that `options(RiskMap.messages = FALSE)`
#'   silences them all.
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
#' @param call Where the error was raised, by default the function that
#'   called `riskmap_abort()`. The error names the outermost function at or
#'   above it on the call stack that users call, an exported function or a
#'   registered S3 method (`user_facing_frame()`), so that internal functions
#'   need not pass `call` down. A call or `NULL` is used as it is.
#' @param user_facing Whether to look for that function. `riskmap_bug()`
#'   names the internal function instead, for the report.
#' @param .envir The environment in which the message is interpolated.
#' @noRd
riskmap_abort <- function(message, ..., class = NULL,
                          call = rlang::caller_env(), user_facing = TRUE,
                          .envir = parent.frame()) {
  if (user_facing) call <- user_facing_frame(call)
  cli::cli_abort(message, ..., class = c(class, "riskmap_error"),
                 call = call, .envir = .envir)
}

#' The frame of the user-facing RiskMap function the user called
#'
#' Looks at the call stack from the outermost frame down to the frame `env`
#' and returns the first frame whose function is exported by RiskMap or
#' registered as one of its S3 methods. An error raised deep inside `glgpm()`
#' then names `glgpm()`, and one raised in a refit by `glgpm()` inside
#' `assess_prediction()` names `assess_prediction()`, the function the user
#' called. A frame entered through `do.call()` with a function, whose call has
#' no name to show, is passed over, and if it is the only one the error names
#' no function (`NULL`). Returns `env` itself when it is not an environment,
#' is not on the stack, or has no user-facing function above it, as when an
#' internal function is called directly in a test.
#' @noRd
user_facing_frame <- function(env) {
  if (!is.environment(env)) return(env)
  frames <- sys.frames()
  position <- which(vapply(frames, identical, logical(1), env))
  if (!length(position)) return(env)
  namespace <- topenv(environment(user_facing_frame))
  user_facing <- user_facing_functions(namespace)
  unnamed <- FALSE
  for (i in seq_len(max(position))) {
    fn <- sys.function(i)
    if (!identical(environment(fn), namespace)) next
    if (!any(vapply(user_facing, identical, logical(1), fn))) next
    if (has_function_name(sys.call(i))) return(frames[[i]])
    unnamed <- TRUE
  }
  if (unnamed) NULL else env
}

#' Whether a call names its function, as `glgpm(...)` or `RiskMap::glgpm(...)`
#' @noRd
has_function_name <- function(call) {
  head <- call[[1L]]
  is.symbol(head) ||
    (is.call(head) && (identical(head[[1L]], quote(`::`)) ||
                         identical(head[[1L]], quote(`:::`))))
}

#' The exported functions and registered S3 methods of a namespace
#'
#' The S3 methods are read from the registrations of the namespace, because
#' a method for a generic of another package, such as `predict.RiskMap()`, is
#' stored in the method table of that package.
#' @noRd
user_facing_functions <- function(namespace) {
  methods <- getNamespaceInfo(namespace, "S3methods")[, 3L]
  names <- unique(c(getNamespaceExports(namespace), methods))
  found <- mget(names, envir = namespace, ifnotfound = list(NULL))
  Filter(is.function, found)
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

#' Raise a condition again, unchanged
#'
#' For a handler that catches errors in general but must let some, such as
#' input errors, reach the user as they are, and for warnings caught in a
#' forked process and raised again in the main one. A warning is raised as a
#' warning and anything else as an error.
#' @noRd
riskmap_rethrow <- function(condition) {
  if (inherits(condition, "warning")) {
    warning(condition)
    return(invisible(NULL))
  }
  stop(condition)
}

#' The message of a condition as one plain line, for storing in a result
#'
#' Without the colours and line breaks that cli adds for the console, which
#' depend on the session and would be kept in the stored text.
#' @noRd
condition_text <- function(condition) {
  old <- options(cli.condition_width = Inf, cli.num_colors = 1)
  on.exit(options(old), add = TRUE)
  text <- cli::ansi_strip(conditionMessage(condition))
  gsub("[[:space:]]*\n[[:space:]]*", " ", trimws(text))
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
    user_facing = FALSE,
    .envir = .envir
  )
}
