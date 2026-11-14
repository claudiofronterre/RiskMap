# Progress of long computations

RiskMap shows how far a long computation has gone and about how long it
has left.

- [`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md),
  [`profile.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md),
  [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
  and
  [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md)
  show a progress bar with the work done and an estimate of the time
  left.

- [`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
  shows a bar while it predicts at many locations in blocks, and while
  it propagates parameter uncertainty over its parameter draws.

- [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  shows a line for each stage of Monte Carlo maximum likelihood, with
  the time the stage took. While a stage runs, the line has a spinner
  and the time so far, which move each time the likelihood is evaluated.
  A step computed in one call to compiled code, such as the Langevin
  sampler, shows its line without them.

Progress is shown when the function's argument `messages` is `TRUE` and
the option `RiskMap.progress` is `TRUE`. The option's default is
[`rlang::is_interactive()`](https://rlang.r-lib.org/reference/is_interactive.html),
so progress appears in the console but not in scripts run with
`Rscript`, in documents rendered by knitr or in tests, even when these
are started from the console. Set `options(RiskMap.progress = TRUE)` to
show it in the log of a script too. A log can't redraw a line, so there
a bar writes a line each time another tenth of the work is done, and
each stage of a fit writes a line when it starts and another, with its
time, when it ends. `options(RiskMap.progress = FALSE)` hides progress
everywhere and leaves the other messages.

The bars are drawn by the cli package, so its options apply as well. For
example, a bar appears only once a computation has run for
`getOption("cli.progress_show_after", 2)` seconds.

## Examples

``` r
# Show progress in the log of a script
options(RiskMap.progress = TRUE)

# Hide it everywhere
options(RiskMap.progress = FALSE)
```
