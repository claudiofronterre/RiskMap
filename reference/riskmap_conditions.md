# Conditions raised by RiskMap

Errors, warnings and messages from RiskMap carry classes, so that code
can handle them with
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) or
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html)
without matching their text, which may be reworded between versions.

- `riskmap_error` is the class of every error that RiskMap raises on
  purpose. A more specific class says what went wrong.

  - `riskmap_error_input`, an argument has a value that the function
    can't use. The field `arg` of the condition names the argument.

  - `riskmap_error_fit`, a model can't be fitted or a computation on a
    fit fails, for example because the search for the mode of the latent
    effects fails or every Monte Carlo weight is zero.

  - `riskmap_error_worldpop`, WorldPop can't be reached or has no data
    for the request.

  - `riskmap_error_internal`, a problem in RiskMap itself. Please report
    it.

- `riskmap_warning` is the class of every warning.

  - `riskmap_warning_fit`, a fit may not be reliable, for example
    because the optimiser did not converge or the Monte Carlo error of
    the estimate is not certified.

  - `riskmap_warning_crs`, the projected coordinate reference system
    used distorts distances by more than 1% over the data (see
    [`propose_crs()`](https://claudiofronterre.github.io/RiskMap/reference/propose_crs.md)
    and
    [`crs_distortion()`](https://claudiofronterre.github.io/RiskMap/reference/crs_distortion.md)).

- `riskmap_message` is the class of every message. Functions that report
  on their progress take an argument `messages`, whose default is the
  option `RiskMap.messages`, so that `options(RiskMap.messages = FALSE)`
  silences them all.

## Examples

``` r
# Which argument was rejected
tryCatch(
  get_worldpop("Liberia", 2020),
  riskmap_error_input = function(condition) condition$arg
)
#> [1] "country"
```
