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
#'   blocks.
#' - [glgpm()] shows a line for each stage of Monte Carlo maximum likelihood,
#'   with the time the stage took.
#'
#' Progress is shown when the function's argument `messages` is `TRUE` and
#' the option `RiskMap.progress` is `TRUE`. The option's default is
#' [interactive()], so progress appears in the console but not in scripts run
#' with `Rscript`, in rendered documents or in tests. Set
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
#' [interactive()]) allows it.
#' @noRd
progress_enabled <- function(messages = TRUE) {
  isTRUE(messages) && isTRUE(getOption("RiskMap.progress", interactive()))
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

#' Close a progress bar
#' @noRd
progress_done <- function(bar) {
  if (is.null(bar)) return(invisible(NULL))
  cli::cli_progress_done(id = bar$id, .envir = bar$envir)
  invisible(NULL)
}

#' A line saying which step of a computation is running
#'
#' cli shows the step while it runs and then a tick with the time it took,
#' when the next step of `.envir` starts, when `progress_step_done()` is
#' called, or when the frame `.envir` exits. In a log the step writes a line
#' when it starts and another, with its time, when it ends. `message` is
#' formatted with cli inline markup in `.envir` when the step starts, so
#' the line that ends it reads the same.
#' @noRd
progress_step <- function(message, messages = TRUE, .envir = parent.frame()) {
  if (!progress_enabled(messages)) return(invisible(NULL))
  text <- cli::format_inline(message, .envir = .envir)
  # Already formatted, so cli must not interpolate it again.
  literal <- gsub("}", "}}", gsub("{", "{{", text, fixed = TRUE), fixed = TRUE)
  cli::cli_progress_step(literal, .envir = .envir)
  invisible(NULL)
}

#' Close the current step of `.envir`
#' @noRd
progress_step_done <- function(messages = TRUE, .envir = parent.frame()) {
  if (progress_enabled(messages)) cli::cli_progress_done(.envir = .envir)
  invisible(NULL)
}
