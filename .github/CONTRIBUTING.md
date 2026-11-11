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
- Start each exported function by checking that objects are of the correct class using `stopifnot(inherits())`.
- If objects must be of a certain value i.e. one of certain strings or numbers within certain bounds, check those next.
- There are various `check_*` functions to help with this - use or adapt those rather than adding large chunks of checking code to functions.
- If writing C++ functions, arguments should be checked on the R side, so do not need checking in the C++. 

## Errors, warnings and messages

RiskMap reports conditions with [cli](https://cli.r-lib.org), following the [error messages chapter](https://style.tidyverse.org/errors.html) of the tidyverse style guide. Use `cli::cli_abort()`, `cli::cli_warn()` and `cli::cli_inform()` instead of `stop()`, `warning()` and `message()`. They format inline markup, wrap long lines and name the function that went wrong. The package is being converted, so older code still uses the base functions. `R/worldpop.R` and the `check_*()` helpers in `R/auxiliary.R` show the conventions below.

### Writing the message

- Start with a problem statement. Use "must" when the requirement is clear, as in "`n_sim` must be a single positive integer, not 0.5.", and "can't" when the request is impossible, as in "Constrained estimates can't be requested for 2035."
- Add bullets after the problem statement when they help. A bullet named `"x"` says what went wrong and a bullet named `"i"` gives context or a hint. Two or three bullets are usually enough.

  ```r
  cli::cli_abort(
    c("{.arg ages} must start and end on WorldPop age groups.",
      "x" = "You supplied {.code c({ages[1]}, {ages[2]})}.",
      "i" = "Use, for example, {.code c(20, Inf)} or {.code c(15, 49)}."),
    class = "riskmap_error_input",
    call = call
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
- Pluralise with `{?}`. The number is taken from the value interpolated just before, so `"{n} file{?s}"` gives "1 file" and "3 files", and `{?is/are}` works the same way. A vector is collapsed with commas and "and".
- Interpolated values are not parsed again for markup, so braces in user data are safe. Write a literal brace in the template as `{{` or `}}`.
- Keep each line of a message on one line of source, even past 80 characters, so that it can be found with grep.

### Naming the function that went wrong

By default `cli_abort()` names the function that called it. That is right for an exported function that checks its own arguments. A check helper, or an internal function that checks arguments for an exported one, must name the exported function instead. It takes `call = rlang::caller_env()` and passes it to `cli_abort()`, and an internal function that calls a helper on behalf of an exported one passes its own `call` on. Without this the error names a function the user never called, such as `extract_distances()`. Helpers take the argument name from `rlang::caller_arg()`, so that `check_logical(messages)` reports `messages`.

```r
check_logical <- function(x, name = rlang::caller_arg(x),
                          call = rlang::caller_env()) {
  if (!isTRUE(x) && !isFALSE(x)) {
    supplied <- describe_value(x)
    cli::cli_abort(
      "{.arg {name}} must be {.code TRUE} or {.code FALSE}, not {supplied}.",
      class = "riskmap_error_input",
      call = call
    )
  }
  invisible(TRUE)
}
```

When an error from another package is the cause, pass it as `parent`, so that both are shown. `check_crs()` does this with the error from `sf::st_crs()`.

### Condition classes

Give a condition a class when code or tests need to tell it apart from others. Invalid arguments have class `riskmap_error_input`, which every `check_*()` helper raises. Other failures have a class named after their subject, such as `riskmap_error_worldpop` when WorldPop can't be reached or has no data for a request. Warnings follow the same pattern, `riskmap_warning_<subject>`. A class is easy to add and hard to rename, because users' code may catch it, so choose names with care.

### Testing conditions

- Test an error by its class, as in `expect_error(fit(), class = "riskmap_error_input")`. Add a short pattern, such as the argument name, when the test is about which argument was rejected.
- Use `expect_snapshot(error = TRUE)` when the wording itself matters, for example for an error with hints. A snapshot records the message as users see it, including the function named. Review changes with `testthat::snapshot_review()`.
- Do not match whole sentences, so that the wording can improve without breaking tests.

### Messages and the `messages` argument

Functions that report progress or choices made for the user take `messages`. When it is `TRUE` they report with `cli::cli_inform()`. Use `cli::cli_warn()` instead for a result the user should doubt. Do not use `cat()` or `print()` for progress, because those can't be silenced like messages. Keep them for `print()` methods.

The default of `messages` is the option `RiskMap.messages`, so that a script or a book chapter can silence RiskMap once with `options(RiskMap.messages = FALSE)`.

```r
get_worldpop <- function(country, year, ..., messages = getOption("RiskMap.messages", TRUE))
```

Internal functions keep taking `messages` from their caller without a default. Because cli messages are conditions, users can also silence them with `suppressMessages()`, and `options(rlib_message_verbosity = "quiet")` silences every message raised with `cli::cli_inform()` or `rlang::inform()`, in RiskMap and in other packages.

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

- If functions are renamed or removed, they should initially be left either as aliases to the updated function with a `warning()` message or if no updated function exists, a `stop()` message.
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
