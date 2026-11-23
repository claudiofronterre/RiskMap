## Submission

This is RiskMap 2.0.0, a major release that redesigns the interface of 1.0.0.
`NEWS.md` lists each function of 1.0.0 with its replacement.

## Maintainer

The maintainer changes from Emanuele Giorgi <e.giorgi@bham.ac.uk> to me,
Claudio Fronterre <c.fronterre@bham.ac.uk>. Emanuele Giorgi stays an author and
will confirm the change by email.

## Test environments

- Ubuntu 24.04 (local), R 4.6.2.
- GitHub Actions, with R release on macOS, Windows and Ubuntu, and R-devel and
  R 4.5 (oldrel-1), the version in `Depends`, on Ubuntu. These run the full test
  suite.
- R-hub. On gcc-asan and m1-san the examples and the full test suite ran
  without reports from the address or undefined behaviour sanitizers. Under
  valgrind the examples of the first 22 of the 47 help pages with examples ran
  without errors before R-hub's limit of 6 hours. Valgrind runs of the tests
  alone, one over all the test files and five over parts of them, also stopped
  at that limit. Up to there they reported no invalid reads or writes and no use
  of uninitialised values. Valgrind looks for leaks when a process exits, which
  only those forked by parallel did, and there it found nothing definitely lost.
  rchk reports only Rcpp's own protection helpers (`Armor.h` and `Shield.h`). On
  clang-asan and clang-ubsan the dependency terra could not be installed in
  R-hub's containers, so RiskMap was not checked there.

## R CMD check results

0 errors | 0 warnings | 1 note

The note is about the change of maintainer. It also flags two words in
DESCRIPTION as possibly misspelled. "Fronterre" is my surname and
"geostatistics" is the field of the package.

On my laptop the check gives two more notes, about the flag
`-mno-omit-leaf-frame-pointer`, which comes from Ubuntu's build of R, and about
HTML Tidy, which is not installed there.

The installed size is about 17 MB, of which 15 MB is the compiled code in
`libs`.

`R CMD check --as-cran` takes about 6 minutes on my laptop, with 2 minutes of
tests and 2 minutes of examples, including those in `\donttest{}`. Each
example outside `\donttest{}` runs in under 3 seconds. Tests that take more
than a second are skipped on CRAN and run on GitHub Actions.

## Reverse dependencies

RiskMap has one reverse dependency on CRAN, MBGapp 0.1.0, which imports it.
`R CMD check` of MBGapp passes with RiskMap 2.0.0, with the same NOTE as with
1.0.0 (namespaces in Imports not imported from), since its R code and tests
don't call RiskMap. Its Shiny app does call RiskMap, and those calls stop
working with 2.0.0, because they use functions and arguments of 1.0.0 that
2.0.0 removed or renamed.

- `set_control_sim()`, `convex_hull_sf()`, `pred_over_grid()` and
  `pred_target_grid()` no longer exist.
- `glgpm()` no longer takes `crs`, `convert_to_crs` and `scale_to_km`, and
  `create_grid()` no longer takes `spat_res` and `grid_crs`.

TODO whether the maintainer of MBGapp was told.
