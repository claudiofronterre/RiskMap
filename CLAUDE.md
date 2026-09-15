# CLAUDE.md

Guidance for Claude Code working in this repository.

## Read first

`.github/CONTRIBUTING.md` covers code style, defensive programming, documentation,
and testing conventions for this package — read it before making changes and
follow it. This file only covers things CONTRIBUTING.md doesn't.

## Useful commands

- `Rscript -e 'devtools::load_all(quiet = TRUE)'` — load the package for manual
  testing (equivalent to `Ctrl+Shift+L` in RStudio).
- `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_dir("tests/testthat")'`
  — run the full test suite.
- `Rscript -e 'devtools::document(quiet = TRUE)'` — regenerate `man/*.Rd` and
  `NAMESPACE` from roxygen comments. Run this after editing any `##'` block.
- Running tests can leave a stray `tests/testthat/Rplots.pdf` behind from
  plotting code — safe to delete, not a real change to include in a diff.

## Test fixtures

`tests/testthat/helper-data.R` is sourced automatically before the test suite
runs and defines the shared fixtures used across most test files:
`gaussian_data`, `binomial_data`, `poisson_data` (small, `n = 10`, synthetic
datasets) and pre-fitted models `gaussian_model`, `gaussian_offset_model`,
`gaussian_intercept_model`, `binomial_model`, `poisson_model`, plus
`control_mcmc` (a reduced-iteration `set_control_mcmc()` for fast fits). Prefer
reusing these over creating new fixtures, unless a test specifically needs
different data (e.g. a covariate name collision, or a formula shape none of
the fixtures cover). If you ever need to produce examples of bugs, use these
same objects.

## Commit messages

Reference the issue number where one exists, e.g. `fix thing #92`.
