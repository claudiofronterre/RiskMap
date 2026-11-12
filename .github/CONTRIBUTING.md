# Contributing to RiskMap

This guide outlines how to propose a change to RiskMap. It is currently a draft document and aspirational as the current code base may not comply with all conventions.

## Branches

- The main branch should be consistent with the current version on CRAN.
- The dev branch contains updates that have been merged, but not yet released.
- Feature branches should be created by branching off from dev. If the branch relates to an open issue, then include this in the branch name e.g. `40-fix_something`.
- Once a feature is complete, open a pull request and ask another contributor to review the changes.
- If small changes are required during review, make them yourself, but if you are uncertain add a comment in the review.

## Code style

- New code should follow the [tidyverse style guide](https://style.tidyverse.org). 
- Use spaces around operators e.g. `like == this` not `like==this`.
- Use verbose variable names to make it clear what they represent. Do not use acronyms or single letters.
- Make function names verbs e.g. `fit_model()` not `model()`.
- Use `snake_case` for function and variable names. Do not use `dot.notation`
- Where it may not be comprehensible from the code, add comments to explain what is happening, but try to make code self-explanatory through sensible function and variable names.
- Simple data manipulation should use base R, but for more complex manipulations use dplyr.
- If using variables present in the data in dplyr, append `.data$` in front of them to avoid generating warnings about undefined variables. See [Programming with dplyr](https://dplyr.tidyverse.org/articles/programming.html) for more details.
- Outputs should be assigned classes with `class()` and if necessary, assign multiple classes e.g. to model families.
- The classes should be used in other functions to determine what code to run rather than inspecting objects to see what the contain and inferring from that.

## Defensive programming

- Code should be written defensively so that users are provided with informative errors.
- Check the arguments of an exported function when it starts, with the helpers in `R/checks.R` (see below). Internal functions can then trust what they are given.
- If writing C++ functions, arguments should be checked on the R side, so do not need checking in the C++.

## Errors, warnings and messages

RiskMap reports conditions with [cli](https://cli.r-lib.org), following the [error messages chapter](https://style.tidyverse.org/errors.html) of the tidyverse style guide. The package is being converted, so older code still uses base R. `R/worldpop.R` and `R/checks.R` show the conventions below.

### Raising a condition

Raise every error, warning and message with `riskmap_abort()`, `riskmap_warn()` or `riskmap_inform()` from `R/conditions.R`, never with `stop()`, `warning()`, `message()`, `stopifnot()` or cli directly. They wrap `cli::cli_abort()`, `cli::cli_warn()` and `cli::cli_inform()`, so messages take cli markup, and they add the classes `riskmap_error`, `riskmap_warning` and `riskmap_message`. A change that should apply to every condition is then made in one place.

`test-conditions.R` enforces this. It reads every function in the installed package and fails if one raises a condition directly without being listed in `tests/testthat/conditions-to-convert.txt`. When you convert a function, remove it from that list. The test also fails if a listed function no longer needs to be there, so the list can only shrink.

For a state that RiskMap should never reach, which would be a bug in RiskMap rather than a mistake by the user, use `riskmap_bug()`. Its error asks the user to report it.

### Checking arguments

Use the helpers in `R/checks.R` for the common checks rather than writing the test and the message again.

- `check_logical()` for `TRUE` or `FALSE`, as for `messages`.
- `check_positive_integer()`, with `allow_zero` and `allow_null`.
- `check_positive_number()`, for a starting value or, with `type = ""`, an argument.
- `check_range()` and `check_zero_one()` for a number in an interval.
- `check_class()` for an object of a given class, instead of `stopifnot(inherits())`.
- `match_choice()` instead of `match.arg()`. It takes the choices from the default of the argument, accepts a unique abbreviation as `match.arg()` does, and raises a RiskMap error.
- `check_crs()` for a coordinate reference system.

Each helper takes the argument name from `rlang::caller_arg()`, so that `check_logical(messages)` reports `messages`, and stores it in the field `arg` of the error. For an input error that no helper covers, call `abort_input()` with the message and `arg`.

### Naming the function that went wrong

An error names the function the user called, not the internal function that found the problem. `riskmap_abort()` finds it on the call stack, as the outermost function at or above the one raising the error that RiskMap exports or registers as an S3 method. An error raised deep inside `glgpm()` therefore reads "Error in `glgpm()`", and internal functions don't need to pass `call` down. A check helper takes `call = rlang::caller_env()` and passes it on, so that the search starts at the function that called the check, as below. `riskmap_bug()` names the internal function instead, because a bug report needs it.

```r
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
```

When an error from another package is the cause, pass it as `parent`, so that both are shown. `check_crs()` does this with the error from `sf::st_crs()`.

### Writing the message

- Start with a problem statement. Use "must" when the requirement is clear, as in "`n_sim` must be a single positive integer, not 0.5.", and "can't" when the request is impossible, as in "Constrained estimates can't be requested for 2035."
- Add bullets after the problem statement when they help. A bullet named `"x"` says what went wrong and a bullet named `"i"` gives context or a hint. Two or three bullets are usually enough.

  ```r
  abort_input(
    c("{.arg ages} must start and end on WorldPop age groups.",
      "x" = "You supplied {.code c({ages[1]}, {ages[2]})}.",
      "i" = "Use, for example, {.code c(20, Inf)} or {.code c(15, 49)}."),
    arg = "ages", call = call
  )
  ```

- Write full sentences that start with a capital letter and end with a full stop.
- Mark up names and values with inline styles instead of quoting them by hand.
  - `{.arg n_sim}` for an argument, shown as `` `n_sim` ``.
  - `{.val {x}}` for a value. Strings are quoted and numbers coloured.
  - `{.cls RiskMap}` for a class.
  - `{.fn glgpm}` for a function, shown as `` `glgpm()` ``.
  - `{.file {path}}` for a file and `{.url {url}}` for a web address.
  - `{.code TRUE}` for code, and also for an element of a list such as `start_pars`, as in `{.code sigma2}`. The `{.field}` style is meant for such names but shows them without quotes in plain text.
  - `{.obj_type_friendly {x}}` for the type of an object, such as "a data frame". The internal `describe_value()` shows a single value as itself and anything else by its type, for the "not ..." part of an input error.
- Pluralise with `{?}`. The number is taken from the value interpolated just before, so `"{n} file{?s}"` gives "1 file" and "3 files", and `{?is/are}` works the same way. A vector is collapsed with commas and "and", or with "or" inside `{.or }`.
- Interpolated values are not parsed again for markup, so braces in user data are safe. Write a literal brace in the template as `{{` or `}}`.
- Keep each line of a message on one line of source, even past 80 characters, so that it can be found with grep.

### Condition classes

The classes are documented for users in `?riskmap_conditions` (`R/conditions.R`). Every error has class `riskmap_error`, every warning `riskmap_warning` and every message `riskmap_message`. A more specific class goes before these when code or tests need to tell a condition apart. Invalid arguments have class `riskmap_error_input`, which every helper and `abort_input()` raise. Other failures have a class named after their subject, such as `riskmap_error_worldpop` when WorldPop can't be reached or has no data for a request, and warnings follow the same pattern, `riskmap_warning_<subject>`. A class is easy to add and hard to rename, because users' code may catch it, so choose names with care and add each new one to `?riskmap_conditions`.

### Testing conditions

- Test an input error with `expect_input_error(code, "n_bins")`, from `tests/testthat/helper-conditions.R`. It checks the class `riskmap_error_input` and that the field `arg` is `"n_bins"`, without matching the wording.
- Test another error by its class, as in `expect_error(code, class = "riskmap_error_worldpop")`.
- Use `expect_snapshot(error = TRUE)` when the wording itself matters, for example for an error with hints. A snapshot records the message as users see it, including the function named. Review changes with `testthat::snapshot_review()`, and commit new snapshots, because a snapshot missing on CI fails there.
- Do not match whole sentences, so that the wording can improve without breaking tests.

### Messages and the `messages` argument

Functions that report progress or choices made for the user take `messages`. When it is `TRUE` they report with `riskmap_inform()`. Use `riskmap_warn()` instead for a result the user should doubt. Do not use `cat()` or `print()` for progress, because those can't be silenced like messages. Keep them for `print()` methods.

The default of `messages` is the option `RiskMap.messages`, so that a script or a book chapter can silence RiskMap once with `options(RiskMap.messages = FALSE)`.

```r
get_worldpop <- function(country, year, ..., messages = getOption("RiskMap.messages", TRUE))
```

Internal functions keep taking `messages` from their caller without a default. Because messages are conditions, users can also silence them with `suppressMessages()`, and `options(rlib_message_verbosity = "quiet")` silences every message raised through cli or rlang, in RiskMap and in other packages.

## Long computations

### Repeated work

Run repeated, independent work, such as bootstrap replicates, profile walks or simulations, through `riskmap_map()` in `R/map.R`, not through `lapply()`, `parallel::mclapply()` or a loop. A loop that can't be written that way, such as the test sets of `assess_prediction()`, which share state between models, moves a bar with `progress_tick()` instead. On one core it runs the tasks in order. With `cores` above 1 it runs each task in a forked process, at most `cores` at a time, and collects the results as the processes finish, so that a progress bar moves while they run. It stops the other processes when a task fails, reports a process that died, raises the warnings of a forked process again in the main one, uses one BLAS thread while forking, and runs on one core on Windows. Two rules keep the results independent of `cores`.

- A task that draws random numbers sets its own seed, derived from the task.
- A task makes no call to sf, because GDAL and PROJ can deadlock in a forked process. Prepare everything that needs sf before `riskmap_map()`.

When a task fails, `riskmap_map()` stops with an error that names the task, from `labels`, such as "The profile of phi below the estimate failed.", keeps the task's error as its cause and keeps its RiskMap classes and `arg`, so that a handler for `riskmap_error_fit` still catches it. With `errors = "return"` a failed task's error takes the place of its result, for callers that record failures instead of stopping. Store the text of a condition in a result with `condition_text()`, which leaves out the colours and line breaks cli adds for the console.

### Progress

Show progress with the helpers in `R/progress.R`, never with `message()` inside a loop.

- `progress_bar(name, total, messages)` makes a bar for work whose amount is known. Pass it to `riskmap_map()` as `progress`, or move it with `progress_tick(bar)`, and close it with `progress_done(bar)`. When the units are only estimates, such as the expected number of profile points, use `unit = "percent"` and give `riskmap_map()` the `weights` of the tasks.
- `progress_step(message, messages)` shows a line for each stage of work whose number of stages isn't known in advance, such as the stages of Monte Carlo maximum likelihood. `progress_step_done(messages)` closes the last one.

```r
bar <- progress_bar("Refitting replicates", n_sim, messages)
results <- riskmap_map(seq_len(n_sim), run_replicate, cores = cores,
                       progress = bar)
progress_done(bar)
```

Both helpers show progress only when `messages` is `TRUE` and the option `RiskMap.progress` allows it (`?riskmap_progress`), so callers don't check either. Progress is for the console. Say what a long computation is about to do and what it found with `riskmap_inform()`, which stays in logs where progress is off.

## Debugging

- To view data internal to a function call `debug(<function name>)` and then call the function to enter the debugger.
- Once complete call `undebug(<function name>)` to stop debugging.
- If you know the location of code that is causing a problem, then add `browser()` above the code chunk, reload the functions with `devtools::load_all()` or `Ctrl + Shift + L` and then call the function. Remember to remove `browser()` afterwards.

## Importing functions

- Generally when functions from other packages are called, they should be scoped using `::` e.g. `dplyr::filter` in which case it is not necessary to use `@importFrom dplyr filter` in the documentation.
- If many functions from one package are used, consider adding `@import <package>` to `R/RiskMap-Package.R` and then scoping is not required. 
- Currently `ggplot2`, `sf` and `stats` are imported.

## Function length

- Functions should be as short as possible so that they are easy to comprehend and easier to test.
- Where functions are over 100 lines long, they should be split into smaller functions.

## Deprecating parameters and functions

- If functions are renamed or removed, they should initially be left either as aliases to the updated function that warn with `riskmap_warn()` or, if no updated function exists, that stop with `riskmap_abort()`.
- Similarly, if function parameters are renamed or removed, they should be left with a default `NULL` parameter and if a user uses the parameter they should be informed that it is deprecated.
- One year after a version containing the deprecated versions has been released on CRAN, the functions and parameters can be removed entirely.

## Documentation

- We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation. This means you can use `` `code` `` to describe variables, `` [function] `` to link to functions and `` - `` to create bulleted lists.
- Make sure that you document the package prior to submitting a PR to update documentation.
- Exported functions should contain at least one example which should run quickly.
- Unexported functions should still be documented but add `@NoRd` to prevent them being included in the manual.

## Testing

- We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. Each function should have unit tests to verify that they function as expected and fail gracefully when incorrect inputs are provided.
- Test data is included in `tests/testthat/helper-data.R` which is loaded when `devtools::load_all()` is used and available to use in all the tests.
- Tests affecting a change made should be run locally prior to submitting a PR, but the whole suite will be run whenever changes are pushed.
