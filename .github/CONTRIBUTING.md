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

## Defensive programming

- Code should be written defensively so that users are provided with informative errors.
- Start each exported function by checking that objects are of the correct class using `stopifnot(inherits())`.
- If objects must be of a certain value i.e. one of certain strings or numbers within certain bounds, check those next. 

## Debugging

- To view data internal to a function call `debug(<function name>)` and then call the function to enter the debugger.
- Once complete call `undebug(<function name>)` to stop debugging.
- If you know the location of code that is causing a problem, then add `browser()` above the code chunk, reload the functions with `devtools::load_all()` or `Ctrl + Shift + L` and then call the function. Remember to remove `browser()` afterwards.

## Importing functions

- Generally when functions from other packages are called, they should be scoped using `::` e.g. `dplyr::filter` in which case it is not necessary to use `@importFrom dplyr filter` in the documentation.
- If many functions from one package are used, consider using `@import <package>` instead and then scoping is not required. Example packages where this approach may be sensible are `stats`, `dplyr`, `ggplot2` and `sf` 

## Function length

- Functions should be as short as possible so that they are easy to comprehend and easier to test.
- Where functions are over 100 lines long, they should be split into smaller functions.

## Deprecating parameters and functions

- If functions are renamed or removed, they should initially be left either as aliases to the updated function with a `warning()` message or if no updated function exists, a `stop()` message.
- Similarly, if function parameters are renamed or removed, they should be left with a default `NULL` parameter and if a user uses the parameter they should be informed that it is deprecated.
- One year after a version containing the deprecated versions has been released on CRAN, the functions and parameters can be removed entirely.

## Documentation

- We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation. Make sure that you document the package prior to submitting a PR to update documentation.
- Exported functions should contain at least one example which should run quickly.
- Unexported functions should still be documented but add `@NoRd` to prevent them being included in the manual.

## Testing

- We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. Each function should have unit tests to verify that they function as expected and fail gracefully when incorrect inputs are provided.
- Test data is included in `tests/testthat/helper-data.R` which is loaded when `devtools::load_all()` is used and available to use in all the tests.
- Tests affecting a change made should be run locally prior to submitting a PR, but the whole suite will be run whenever changes are pushed.
