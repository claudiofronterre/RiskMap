## Submission

This is RiskMap 2.0.0, a major release that redesigns the interface of 1.0.0.
`NEWS.md` lists each function of 1.0.0 with its replacement.

## Maintainer

The maintainer changes from Emanuele Giorgi <e.giorgi@bham.ac.uk> to me,
Claudio Fronterre <c.fronterre@bham.ac.uk>. Emanuele Giorgi stays an author and
will confirm the change by email.

## Test environments

- Ubuntu 24.04 (local), R 4.6.2.
- GitHub Actions, R release on macOS, Windows and Ubuntu, and R-devel and
  R 4.5 (oldrel-1) on Ubuntu.
- R-hub, the clang-asan, clang-ubsan, gcc-asan, valgrind, rchk and m1-san
  platforms. TODO results.
- R 4.1.3 on Ubuntu 20.04 (rocker/geospatial), for the minimum version of R in
  `Depends`. TODO result.

## R CMD check results

0 errors | 0 warnings | TODO notes

TODO the NOTEs.

`R CMD check --as-cran` takes TODO minutes on my laptop, with TODO seconds of
tests and TODO seconds of examples, including those in `\donttest{}`.

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
