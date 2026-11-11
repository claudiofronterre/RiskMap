fail_third <- function(k) if (k == 3) stop("The third task fails.") else k

test_that("riskmap_map returns the results in the order of the tasks", {
  expect_identical(riskmap_map(1:5, function(k) k^2), as.list((1:5)^2))
  expect_identical(riskmap_map(list(), identity), list())
  expect_identical(riskmap_map(1:2, function(k) NULL), list(NULL, NULL))
  expect_identical(riskmap_map(list(c("a", "b"), "c"), length), list(2L, 1L))
})

test_that("riskmap_map gives the same results on several cores", {
  skip_on_os("windows")
  draw <- function(seed) {
    set.seed(seed)
    runif(3)
  }
  expect_identical(riskmap_map(1:6, draw, cores = 2), riskmap_map(1:6, draw))
  # The first task finishes last, and its result still comes first.
  slow_first <- function(k) {
    Sys.sleep(if (k == 1) 0.5 else 0)
    k
  }
  expect_identical(riskmap_map(1:4, slow_first, cores = 2), as.list(1:4))
})

test_that("a failing task stops riskmap_map, with its error as the cause", {
  error <- expect_error(riskmap_map(1:4, fail_third), "Task 3 of 4 failed",
                        class = "riskmap_error")
  expect_match(conditionMessage(error$parent), "third task fails")

  skip_on_os("windows")
  error <- expect_error(riskmap_map(1:4, fail_third, cores = 2),
                        "Task 3 of 4 failed", class = "riskmap_error")
  expect_match(conditionMessage(error$parent), "third task fails")
})

test_that("riskmap_map can return errors in place of results", {
  results <- riskmap_map(1:4, fail_third, errors = "return")
  expect_s3_class(results[[3]], "error")
  expect_identical(results[-3], list(1L, 2L, 4L))

  skip_on_os("windows")
  expect_identical(riskmap_map(1:4, fail_third, cores = 2, errors = "return"),
                   results)
  # A process killed while running its task returns no result.
  killed_second <- function(k) {
    if (k == 2) tools::pskill(Sys.getpid(), tools::SIGKILL)
    k
  }
  results <- riskmap_map(1:3, killed_second, cores = 2, errors = "return")
  expect_s3_class(results[[2]], "riskmap_error")
  expect_match(conditionMessage(results[[2]]), "ended without returning")
  expect_identical(results[-2], list(1L, 3L))
})

test_that("progress is shown when messages and the option allow it", {
  withr::local_options(RiskMap.progress = FALSE)
  expect_false(progress_enabled(TRUE))
  expect_null(progress_bar("Working", 10))

  withr::local_options(RiskMap.progress = TRUE)
  expect_true(progress_enabled(TRUE))
  expect_false(progress_enabled(FALSE))
  expect_null(progress_bar("Working", 10, messages = FALSE))
  expect_null(progress_bar("Working", 0))
  # Functions that take no bar accept NULL.
  expect_null(progress_tick(NULL))
  expect_null(progress_done(NULL))
})

test_that("a task moves the bar by at most its weight", {
  withr::local_options(RiskMap.progress = TRUE, cli.dynamic = FALSE,
                       cli.progress_show_after = Inf)
  report_six_times <- function(k, report) {
    for (i in 1:6) report()
    k
  }
  bar <- progress_bar("Working", 8, unit = "percent")
  riskmap_map(1:2, function(k, report) {
    report()
    expect_equal(bar$current, 4 * (k - 1) + 1)
    report_six_times(k, report)
  }, progress = bar, weights = c(4, 4), reports_progress = TRUE)
  expect_equal(bar$current, 8)

  skip_on_os("windows")
  bar <- progress_bar("Working", 8, unit = "percent")
  riskmap_map(1:2, report_six_times, cores = 2, progress = bar,
              weights = c(4, 4), reports_progress = TRUE)
  expect_equal(bar$current, 8)
})

test_that("in a log a bar writes a line for each tenth of the work", {
  withr::local_options(RiskMap.progress = TRUE, cli.dynamic = FALSE,
                       cli.progress_show_after = 0)
  run <- function() {
    bar <- progress_bar("Refitting replicates", 20)
    riskmap_map(1:20, identity, progress = bar)
    progress_done(bar)
  }
  lines <- trimws(capture_messages(run()))
  lines <- lines[nzchar(lines)]
  expect_length(lines, 10L)
  expect_match(lines[1], "^Refitting replicates 2/20 \\| ETA")
  expect_match(lines[10], "^Refitting replicates 20/20 \\| ETA")
})

test_that("in a log a step writes a line when it starts and when it ends", {
  withr::local_options(RiskMap.progress = TRUE, cli.dynamic = FALSE)
  run <- function() {
    for (budget in c(1000, 2000)) {
      progress_step("Stage {budget / 1000}, with {budget} draws {{1}}")
    }
    progress_step_done()
  }
  lines <- trimws(capture_messages(run()))
  lines <- lines[nzchar(lines)]
  expect_length(lines, 4L)
  expect_match(lines[1], "Stage 1, with 1000 draws {1}", fixed = TRUE)
  expect_match(lines[2], "Stage 1, with 1000 draws {1} [", fixed = TRUE)
  expect_match(lines[4], "Stage 2, with 2000 draws {1} [", fixed = TRUE)
  expect_no_message(progress_step("Stage 1", messages = FALSE))
  withr::local_options(RiskMap.progress = FALSE)
  expect_no_message(progress_step("Stage 1"))
})
