# Format RiskMap Results as a Table

Converts a fitted "RiskMap" model, cross-validation results or polygon
predictions into a table that renders directly in Quarto, R Markdown,
HTML, LaTeX and the R console.

## Usage

``` r
to_table(object, ...)

# Default S3 method
to_table(object, digits = 3, col_names = NULL, ...)
```

## Arguments

- object:

  An object of class "RiskMap" resulting from a call to
  [`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
  a "summary.RiskMap" object, a "summary.RiskMap_cross_validation"
  object, or a "RiskMap_prediction" object (see
  [`to_table.RiskMap_prediction`](https://claudiofronterre.github.io/RiskMap/reference/to_table.RiskMap_prediction.md)).

- ...:

  Additional arguments passed to
  [`kable`](https://rdrr.io/pkg/knitr/man/kable.html).

- digits:

  A non-negative integer giving the number of decimal places used to
  display numeric results.

- col_names:

  Optional column names, replacing the default headers.

## Value

An object of class "knitr_kable" that can be rendered directly.

## Details

This function creates a presentation-ready summary table from a fitted
"RiskMap" model or cross-validation results for multiple models. Use
[`coef`](https://rdrr.io/r/stats/coef.html) or
[`summary`](https://rdrr.io/r/base/summary.html) when numeric results
are required for further analysis. Numeric values use fixed notation
with `digits` decimal places, except when scientific notation is needed
to represent very large or very small values clearly. Values that are
compared share one notation: the estimate and confidence limits of a
parameter, and one metric across models.

When the input is a "RiskMap" model object, the table includes:

- Regression coefficients with their estimates, confidence intervals,
  and p-values.

- Parameters for the spatial process.

- Random effect variances.

- The parameters of the family that were estimated, such as the
  measurement error variance of a Gaussian model.

When the input is a cross-validation summary object
("summary.RiskMap_cross_validation"), the table includes:

- A row for each model being compared.

- Performance metrics such as CRPS and SCRPS for each model.

## See also

[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`summary.RiskMap_cross_validation`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md)

## Examples

``` r
data(italy_sim)

fit <- glgpm(
  formula = y ~ gp(),
  data = italy_sim[1:100,],
  family = "gaussian",
  messages = FALSE
)
to_table(fit)
#> 
#> 
#> |                         | Estimate| Lower limit| Upper limit|
#> |:------------------------|--------:|-----------:|-----------:|
#> |(Intercept)              |    9.323|       7.701|      10.946|
#> |Spatial process var.     |    4.849|       1.820|      12.924|
#> |Spatial corr. scale (km) |   39.015|       2.777|     548.084|
#> |Measurement error var.   |    0.179|       0.133|       0.239|
to_table(fit, digits = 1)
#> 
#> 
#> |                         | Estimate| Lower limit| Upper limit|
#> |:------------------------|--------:|-----------:|-----------:|
#> |(Intercept)              |      9.3|         7.7|        10.9|
#> |Spatial process var.     |      4.8|         1.8|        12.9|
#> |Spatial corr. scale (km) |     39.0|         2.8|       548.1|
#> |Measurement error var.   |      0.2|         0.1|         0.2|
```
