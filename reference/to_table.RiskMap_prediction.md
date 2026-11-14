# Table of polygon predictions

Formats the polygon summaries returned by
[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
as a presentation-ready table, with one row per polygon.

## Usage

``` r
# S3 method for class 'RiskMap_prediction'
to_table(
  object,
  label = NULL,
  aggregate = NULL,
  summaries = NULL,
  quantity = "response",
  sort_by = NULL,
  digits = 1,
  col_names = NULL,
  ...
)
```

## Arguments

- object:

  A `RiskMap_prediction` object from
  [`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
  with `boundaries`.

- label:

  Name of the column of the boundaries used to name the rows. By
  default, the first character or factor column, or the polygon number
  when there is none.

- aggregate:

  Aggregations to include, such as `"mean"` or `"sum"`. By default, all
  those calculated.

- summaries:

  Summaries to include, such as `"mean"`, `"exceedance"` or `"sd"`. By
  default, all those calculated.

- quantity:

  Predictive quantity to tabulate. Defaults to `"response"`.

- sort_by:

  Name of a summary, such as `"mean"`, used to sort the rows in
  decreasing order, using the first aggregation that has it. By default,
  rows keep the order of the boundaries.

- digits:

  Number of decimal places for values that are not counts.

- col_names:

  Optional column names, replacing the default headers.

- ...:

  Additional arguments passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html).

## Value

An object of class `knitr_kable`.

## Details

Each aggregation gives one or more columns. When the predictive mean and
both interval limits are available, they are combined in one column as
"mean (lower to upper)". Exceedance probabilities are omitted for sums
of a probability, because the thresholds refer to the probability
itself. Probabilities, including prevalence averaged over a polygon and
exceedance probabilities, are shown as percentages. Sums are shown to
three significant digits with thousands separators. The default headers
name the quantity as in the maps drawn by
[`plot.RiskMap_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_prediction.md).

## See also

[`predict.RiskMap()`](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md),
[`plot.RiskMap_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/plot.RiskMap_prediction.md)

## Examples

``` r
data(italy_sim)
example_data <- italy_sim[1:50, ]
fit <- glgpm(y ~ gp(), data = example_data, family = "gaussian",
             messages = FALSE)
areas <- sf::st_sf(name = c("West", "East"),
                   geometry = sf::st_make_grid(example_data, n = c(2, 1)))
prediction <- predict(fit, newdata = example_data, boundaries = areas,
                      aggregate = "mean",
                      summaries = c("mean", "lower", "upper"),
                      messages = FALSE)
#> Warning: The conditional prediction covariance of 2 of 2 blocks needed a diagonal
#> jitter, of up to 1e-12, to be factorised.
to_table(prediction)
#> 
#> 
#> |Area | Mean y (95% interval)|
#> |:----|---------------------:|
#> |West |    10.0 (9.8 to 10.1)|
#> |East |      9.6 (9.4 to 9.8)|
```
