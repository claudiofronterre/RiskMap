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
runs and defines the shared fixtures used across most test files. Each
dataset has 20 locations and is simulated from the model that the tests fit
to it. `gaussian_data`, `binomial_data` and `poisson_data` come from
`y ~ cov + gp()`, the `*_re_data` from `y ~ cov + gp() + re(i)`, the
`*_two_re_data` from `y ~ cov + gp() + re(i, j)`, `binomial_nugget_data` from
`y ~ cov + gp(nugget = TRUE)`, `negative_binomial_data` from a negative
binomial `y ~ cov + gp()`, and `gaussian_offset_data` and
`gaussian_intercept_data` from their models. Fit a model to data that contain
every effect it estimates (#206), and simulate data for another model with
`simulate_fixture()`. The pre-fitted models are `gaussian_model`,
`gaussian_no_cov_model`, `gaussian_offset_model`, `gaussian_intercept_model`,
`binomial_model` and `poisson_model`, plus `control_mcmc` and
`fixture_control_mcml` for fast fits.
Prefer reusing these over creating new fixtures, unless a test specifically
needs different data (e.g. a covariate name collision, or a formula shape none
of the fixtures cover). If you ever need to produce examples of bugs, use
these same objects.

## Golden battery

`tools/golden/` (outside the built package, see `.Rbuildignore`) runs a fixed set of
fits, predictions, profiles, bootstraps and engine evaluations and compares two
installs of RiskMap task by task. A PR that is meant to change no number must leave
every task identical bit for bit, and a PR that moves fits shows by how much with the
declared paths or the tolerance mode of `compare.R`. Run both arms with the same
battery, from installs of clean `git archive`s (never `load_all()`), with one BLAS
thread and on one machine type, on Wayland with `tools/golden/hpc/push.sh` and
`submit.sh`. Say in the PR which arms were compared and what the comparison found.
`tools/golden/README.md` has the details and the current reference.

## Commit messages

Reference the issue number where one exists, e.g. `fix thing #92`.
