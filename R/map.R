#' Apply a function to tasks, in forked processes when `cores` is above 1
#'
#' The one loop of RiskMap for repeated, independent work, such as bootstrap
#' replicates, the walks of a profile or simulation studies. On one core the
#' tasks run in order. With more, each task runs in a process forked by
#' [parallel::mcparallel()], at most `cores` at a time, and the results are
#' collected as the processes finish, so that a progress bar moves while
#' they run. Windows can't fork, so there the tasks run on one core. A task
#' must not depend on the tasks before it, and one that draws random numbers
#' must set its own seed, so that the results don't depend on `cores` or on
#' the order in which the tasks finish.
#'
#' Forked processes must not call sf, because GDAL and PROJ can deadlock in a
#' fork, so the caller prepares anything that needs it.
#' @param tasks A vector or list, one element per task.
#' @param work A function of one task. With `reports_progress`, a function
#'   of the task and of `report(amount)`, which it calls as it goes to move
#'   the bar by part of the task's weight.
#' @param cores Number of processes.
#' @param progress A bar from `progress_bar()`, or `NULL`.
#' @param weights The share of the bar of each task, 1 each by default. When
#'   a task finishes, the bar moves by what remains of its weight.
#' @param reports_progress Whether `work` takes `report`. In a forked process
#'   `report` does nothing, and the bar moves when the task finishes.
#' @param errors `"stop"` stops at the first task that fails, and the other
#'   processes, with an error whose cause is the task's error. `"return"`
#'   puts the error in place of the task's result. A process that ends
#'   without returning a result, for example because the system ran out of
#'   memory, fails in the same way.
#' @param one_thread Whether the tasks run with one BLAS and OpenMP thread,
#'   when the RhpcBLASctl package is installed. The settings are restored
#'   afterwards. Always the case with more than one core, because a process
#'   forked while BLAS runs several threads can deadlock, and processes that
#'   each use every core slow one another down.
#' @param call The call named in an error, by default the caller.
#' @return A list of the results, in the order of `tasks`.
#' @noRd
riskmap_map <- function(tasks, work, cores = 1L, progress = NULL,
                        weights = NULL, reports_progress = FALSE,
                        errors = c("stop", "return"), one_thread = FALSE,
                        call = rlang::caller_env()) {
  errors <- match.arg(errors)
  n <- length(tasks)
  if (n == 0L) return(list())
  weights <- weights %||% rep(1, n)
  if (length(weights) != n) {
    riskmap_bug("{.arg weights} must have one value per task.")
  }
  if (.Platform$OS.type == "windows") cores <- 1L
  cores <- min(as.integer(cores), n)
  if (one_thread || cores > 1L) {
    restore_threads <- use_one_thread()
    on.exit(restore_threads(), add = TRUE)
  }

  # How much of each task's weight the bar has shown, so that partial
  # reports never move it by more than the weight.
  shown <- numeric(n)
  advance <- function(k, amount) {
    amount <- min(amount, weights[k] - shown[k])
    if (amount > 0) {
      shown[k] <<- shown[k] + amount
      progress_tick(progress, amount)
    }
  }
  results <- vector("list", n)
  finish <- function(k, outcome) {
    force(outcome)
    advance(k, weights[k])
    if (outcome$ok) {
      results[k] <<- list(outcome$value)
    } else if (errors == "return") {
      results[k] <<- list(outcome$error)
    } else {
      abort_task(k, n, outcome$error, call)
    }
  }

  if (cores <= 1L) {
    for (k in seq_len(n)) {
      task_work <- if (reports_progress) {
        function(task) work(task, function(amount = 1) advance(k, amount))
      } else {
        work
      }
      outcome <- run_task(task_work, tasks[[k]], errors, k, n, call)
      finish(k, outcome)
    }
  } else {
    forked_work <- if (reports_progress) {
      function(task) work(task, function(amount = 1) invisible(NULL))
    } else {
      work
    }
    run_forked(tasks, forked_work, cores, finish)
  }
  results
}

#' Run one task, as `list(ok = TRUE, value =)` or `list(ok = FALSE, error =)`
#'
#' With `errors = "stop"` an error is raised from a calling handler, so that
#' its backtrace still reaches into the task.
#' @noRd
run_task <- function(work, task, errors, k, n, call) {
  if (errors == "stop") {
    return(withCallingHandlers(
      list(ok = TRUE, value = work(task)),
      error = function(error) abort_task(k, n, error, call)
    ))
  }
  tryCatch(list(ok = TRUE, value = work(task)),
           error = function(error) list(ok = FALSE, error = error))
}

#' @noRd
abort_task <- function(k, n, error, call) {
  riskmap_abort("Task {k} of {n} failed.", parent = error, call = call)
}

#' Run tasks in forked processes, at most `cores` at a time
#'
#' `finish(k, outcome)` is called in the parent for each task as its process
#' finishes. If it raises an error, or the user interrupts, the processes
#' still running are stopped.
#' @noRd
run_forked <- function(tasks, work, cores, finish) {
  jobs <- list()
  on.exit(stop_jobs(jobs), add = TRUE)
  next_task <- 1L
  while (next_task <= length(tasks) || length(jobs) > 0L) {
    while (length(jobs) < cores && next_task <= length(tasks)) {
      job <- start_job(work, tasks[[next_task]])
      job$task <- next_task
      jobs[[as.character(job$pid)]] <- job
      next_task <- next_task + 1L
    }
    # A process that ended without a result is collected as NULL, and
    # mccollect() warns about it, which the outcome below reports instead.
    collected <- suppressWarnings(
      parallel::mccollect(unname(jobs), wait = FALSE, timeout = 1)
    )
    for (pid in names(collected)) {
      k <- jobs[[pid]]$task
      outcome <- collected[[pid]]
      jobs[[pid]] <- NULL
      if (!is.list(outcome) || is.null(outcome$ok)) {
        outcome <- list(ok = FALSE, error = rlang::error_cnd(
          "riskmap_error",
          message = paste("The process running the task ended without",
                          "returning a result, for example because the",
                          "system ran out of memory.")
        ))
      }
      finish(k, outcome)
    }
  }
  invisible(NULL)
}

#' Fork a process running one task, which sends back its outcome
#' @noRd
start_job <- function(work, task) {
  parallel::mcparallel(
    tryCatch(list(ok = TRUE, value = work(task)),
             error = function(error) list(ok = FALSE, error = error)),
    mc.set.seed = FALSE
  )
}

#' Stop forked processes that are still running and collect them
#' @noRd
stop_jobs <- function(jobs) {
  if (!length(jobs)) return(invisible(NULL))
  pids <- vapply(jobs, function(job) as.integer(job$pid), integer(1))
  tools::pskill(pids, tools::SIGTERM)
  suppressWarnings(parallel::mccollect(unname(jobs), wait = TRUE))
  invisible(NULL)
}

#' Use one BLAS and OpenMP thread, returning a function that restores the
#' settings. Nothing changes without RhpcBLASctl.
#' @noRd
use_one_thread <- function() {
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    return(function() invisible(NULL))
  }
  blas_threads <- RhpcBLASctl::blas_get_num_procs()
  omp_threads <- RhpcBLASctl::omp_get_max_threads()
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  function() {
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    RhpcBLASctl::omp_set_num_threads(omp_threads)
    invisible(NULL)
  }
}
