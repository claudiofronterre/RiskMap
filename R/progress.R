#' Progress of long computations
#'
#' @description
#' RiskMap shows how far a long computation has gone and about how long it
#' has left.
#'
#' - [parametric_bootstrap()], [profile.RiskMap()], [assess_prediction()] and
#'   [assess_simulation()] show a progress bar with the work done and an
#'   estimate of the time left.
#' - [predict.RiskMap()] shows a bar while it predicts at many locations in
#'   blocks, and while it propagates parameter uncertainty over its parameter
#'   draws.
#' - [glgpm()] shows a line for each stage of Monte Carlo maximum likelihood,
#'   with the time the stage took. While a stage runs, the line has a spinner
#'   and the time so far, which move each time the likelihood is evaluated.
#'   A step computed in one call to compiled code, such as the Langevin
#'   sampler, shows its line without them.
#'
#' Progress is shown when the function's argument `messages` is `TRUE` and
#' the option `RiskMap.progress` is `TRUE`. The option's default is
#' [rlang::is_interactive()], so progress appears in the console but not in
#' scripts run with `Rscript`, in documents rendered by knitr or in tests,
#' even when these are started from the console. Set
#' `options(RiskMap.progress = TRUE)` to show it in the log of a script too.
#' A log can't redraw a line, so there a bar writes a line each time another
#' tenth of the work is done, and each stage of a fit writes a line when it
#' starts and another, with its time, when it ends.
#' `options(RiskMap.progress = FALSE)` hides progress everywhere and leaves
#' the other messages.
#'
#' The bars are drawn by the cli package, so its options apply as well. For
#' example, a bar appears only once a computation has run for
#' `getOption("cli.progress_show_after", 2)` seconds.
#'
#' @examples
#' # Show progress in the log of a script
#' options(RiskMap.progress = TRUE)
#'
#' # Hide it everywhere
#' options(RiskMap.progress = FALSE)
#' @name riskmap_progress
NULL

#' Whether progress is shown
#'
#' When `messages` is `TRUE` and the option `RiskMap.progress` (by default
#' [rlang::is_interactive()], which is `FALSE` while knitr renders and while
#' testthat runs) allows it.
#' @noRd
progress_enabled <- function(messages = TRUE) {
  isTRUE(messages) &&
    isTRUE(getOption("RiskMap.progress", rlang::is_interactive()))
}

#' A progress bar
#'
#' Every bar in RiskMap is made here, so that their look and their behaviour
#' in a log are decided in one place. The bar belongs to the frame `.envir`,
#' by default the caller's, and closes when that frame exits if
#' `progress_done()` has not closed it before. Returns `NULL`, which
#' `progress_tick()` and `progress_done()` accept and ignore, when progress
#' is not shown.
#'
#' In a terminal that can redraw a line, cli redraws the bar a few times a
#' second. In a log, the bar writes a line each time another tenth of
#' `total` is done, without the bar itself, and only once the computation
#' has run for `cli.progress_show_after` seconds.
#' @param name What is being done, such as "Refitting replicates".
#' @param total The amount of work, in the units of `progress_tick()`.
#' @param messages The `messages` argument of the calling function.
#' @param unit `"count"` shows the work done as `done/total`, for tasks
#'   counted one by one. `"percent"` shows a percentage, for work measured in
#'   approximate units, such as the expected number of profile points.
#' @noRd
progress_bar <- function(name, total, messages = TRUE,
                         unit = c("count", "percent"),
                         .envir = parent.frame()) {
  unit <- match.arg(unit)
  if (!progress_enabled(messages) || !isTRUE(total > 0)) return(NULL)
  dynamic <- cli::is_dynamic_tty()
  amount <- switch(unit,
                   count = "{cli::pb_current}/{cli::pb_total}",
                   percent = "{cli::pb_percent}")
  format <- paste0("{cli::pb_name}", if (dynamic) "{cli::pb_bar} ", amount,
                   " | ETA {cli::pb_eta}")
  bar <- new.env(parent = emptyenv())
  bar$id <- cli::cli_progress_bar(name = name, total = total, format = format,
                                  auto_terminate = FALSE, .envir = .envir)
  bar$envir <- .envir
  bar$total <- total
  bar$current <- 0
  bar$dynamic <- dynamic
  bar$next_report <- total / 10
  bar$started <- proc.time()[["elapsed"]]
  bar
}

#' Advance a progress bar by `amount`, never past its total
#' @noRd
progress_tick <- function(bar, amount = 1) {
  if (is.null(bar) || !isTRUE(amount > 0)) return(invisible(bar))
  bar$current <- min(bar$current + amount, bar$total)
  if (bar$dynamic) {
    cli::cli_progress_update(id = bar$id, set = bar$current,
                             .envir = bar$envir)
    return(invisible(bar))
  }
  tenth <- bar$total / 10
  if (bar$current < bar$next_report - 1e-8 * bar$total) {
    return(invisible(bar))
  }
  bar$next_report <- (floor(bar$current / tenth + 1e-8) + 1) * tenth
  elapsed <- proc.time()[["elapsed"]] - bar$started
  if (elapsed >= getOption("cli.progress_show_after", 2)) {
    cli::cli_progress_update(id = bar$id, set = bar$current, force = TRUE,
                             .envir = bar$envir)
  }
  invisible(bar)
}

#' Redraw a bar in a terminal, so that its time left moves between ticks
#' @noRd
progress_refresh <- function(bar) {
  if (!is.null(bar) && bar$dynamic) {
    cli::cli_progress_update(id = bar$id, set = bar$current,
                             .envir = bar$envir)
  }
  invisible(bar)
}

#' Close a progress bar
#' @noRd
progress_done <- function(bar) {
  if (is.null(bar)) return(invisible(NULL))
  cli::cli_progress_done(id = bar$id, .envir = bar$envir)
  invisible(NULL)
}

#' A line saying which step of a computation is running
#'
#' In a terminal, cli shows the step while it runs and then a tick with the
#' time it took, when the next step of `.envir` starts, when
#' `progress_step_done()` is called, or when the frame `.envir` exits. The
#' line starts with an information sign, and once the computation calls
#' `progress_spin()` it shows a spinner and the time so far instead, so that a
#' long step shows it is still running. In a log, where cli would add a blank
#' line after each step, the step writes a line when it starts and another,
#' with its time, when the next step of `.envir` starts or
#' `progress_step_done()` is called. `message` is formatted with cli inline
#' markup in `.envir` when the step starts, so the line that ends it reads the
#' same.
#' @param inform When progress is hidden but `messages` is `TRUE`, write
#'   `message` as an ordinary message instead, with a full stop, for a step
#'   that used to be a message.
#' @noRd
progress_step <- function(message, messages = TRUE, inform = FALSE,
                          .envir = parent.frame()) {
  if (!progress_enabled(messages)) {
    if (inform && isTRUE(messages)) {
      text <- cli::format_inline(message, .envir = .envir)
      riskmap_inform("{text}.")
    }
    return(invisible(NULL))
  }
  text <- cli::format_inline(message, .envir = .envir)
  if (cli::is_dynamic_tty()) {
    open_spinning_step(text, .envir)
  } else {
    end_log_step(.envir)
    riskmap_inform(c("i" = "{text}"))
    log_steps$open[[length(log_steps$open) + 1L]] <-
      list(envir = .envir, text = text, started = proc.time()[["elapsed"]])
  }
  invisible(NULL)
}

#' The step that `progress_spin()` moves, and the process that opened it
#' @noRd
spinning_step <- new.env(parent = emptyenv())

#' Open a step in a terminal, as `cli::cli_progress_step()` does
#'
#' The running line starts with cli's information sign until the first
#' `progress_spin()`, and then with a spinner, followed at the end of the line
#' by the time so far. The line that ends the step is a tick, or a cross after
#' an error, with the time the step took, as cli's own steps end. cli builds
#' each line when it draws it, through the functions in `extra`, so that the
#' text is shortened to the console against the time shown at that moment.
#' The text is passed as data, so braces in it are printed as they are.
#' @noRd
open_spinning_step <- function(text, .envir) {
  running <- function() {
    if (isTRUE(cli::pb_extra$spinning)) {
      step_line(cli::pb_spin, text, cli::pb_elapsed)
    } else {
      step_line(cli::col_cyan(cli::symbol$info), text)
    }
  }
  ended <- function(sign) {
    force(sign)
    function() step_line(sign, text, cli::pb_elapsed)
  }
  restore <- options(cli.progress_show_after = 0)
  on.exit(options(restore), add = TRUE)
  id <- cli::cli_progress_bar(
    type = "custom", format = "{cli::pb_extra$running()}",
    format_done = "{cli::pb_extra$done()}",
    format_failed = "{cli::pb_extra$failed()}",
    clear = FALSE,
    extra = list(spinning = FALSE, running = running,
                 done = ended(cli::col_green(cli::symbol$tick)),
                 failed = ended(cli::col_red(cli::symbol$cross))),
    .envir = .envir
  )
  cli::cli_progress_update(id = id, force = TRUE, .envir = .envir)
  spinning_step$id <- id
  spinning_step$process <- Sys.getpid()
  spinning_step$spinning <- FALSE
  invisible(id)
}

#' One line of a step, shortened to fit the console with its time
#'
#' cli cuts a line wider than the console at its end, which would take the
#' time with it, so the text is shortened instead, against the width of the
#' time shown, which grows from "[5ms]" to "[1h 21m 58.7s]". `time` is `NULL`
#' for a line without one.
#' @noRd
step_line <- function(sign, text, time = NULL, width = cli::console_width()) {
  clock <- if (is.null(time)) "" else cli::col_grey(paste0(" [", time, "]"))
  room <- width - cli::ansi_nchar(sign) - 1L - cli::ansi_nchar(clock)
  paste0(sign, " ", cli::ansi_strtrim(text, max(1L, room)), clock)
}

#' Move the spinner of the running step
#'
#' Called where a long step does its work again and again, such as each
#' evaluation of the Monte Carlo likelihood. It does nothing when no step is
#' running in a terminal, or in a forked process, which must not draw over
#' the line of its parent. cli redraws the line only a few times a second, so
#' a call costs a few microseconds.
#' @noRd
progress_spin <- function() {
  id <- spinning_step$id
  if (is.null(id) || !identical(spinning_step$process, Sys.getpid())) {
    return(invisible(NULL))
  }
  first <- !isTRUE(spinning_step$spinning)
  spinning_step$spinning <- TRUE
  # The step may have closed with its frame, and cli then no longer knows it.
  tryCatch(
    if (first) {
      cli::cli_progress_update(id = id, extra = list(spinning = TRUE),
                               force = TRUE)
    } else {
      cli::cli_progress_update(id = id)
    },
    error = function(error) spinning_step$id <- NULL
  )
  invisible(NULL)
}

#' Close the current step of `.envir`
#' @noRd
progress_step_done <- function(messages = TRUE, .envir = parent.frame()) {
  if (!progress_enabled(messages)) return(invisible(NULL))
  if (cli::is_dynamic_tty()) {
    cli::cli_progress_done(.envir = .envir)
    spinning_step$id <- NULL
  } else {
    end_log_step(.envir)
  }
  invisible(NULL)
}

#' Steps open in a log, each with its frame, text and start time
#' @noRd
log_steps <- new.env(parent = emptyenv())
log_steps$open <- list()

#' Write the line that ends the open step of `.envir` in a log
#'
#' Steps whose frame has exited, after an error, are dropped without a line.
#' @noRd
end_log_step <- function(.envir) {
  frames <- sys.frames()
  on_stack <- function(step) {
    any(vapply(frames, identical, logical(1), step$envir))
  }
  open <- Filter(on_stack, log_steps$open)
  mine <- vapply(open, function(step) identical(step$envir, .envir),
                 logical(1))
  for (step in open[mine]) {
    took <- format_step_time(proc.time()[["elapsed"]] - step$started)
    text <- step$text
    riskmap_inform(c("v" = "{text} [{took}]"))
  }
  log_steps$open <- open[!mine]
  invisible(NULL)
}

#' A step's time as cli shows it, as in 338ms, 2.2s or 3m 5s
#' @noRd
format_step_time <- function(seconds) {
  if (seconds < 1) return(paste0(round(seconds * 1000), "ms"))
  if (seconds < 60) return(paste0(format(round(seconds, 1), nsmall = 1), "s"))
  paste0(seconds %/% 60, "m ", round(seconds %% 60), "s")
}
