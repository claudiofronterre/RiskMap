# Format RiskMap Results as a Table

Converts a fitted "RiskMap" model, cross-validation results or polygon
predictions into a table that renders directly in Quarto, R Markdown,
HTML, LaTeX and the R console.

## Usage

``` r
to_table(object, ...)

# S3 method for class 'RiskMap'
to_table(object, digits = 3, col_names = NULL, ...)

# S3 method for class 'summary.RiskMap'
to_table(object, digits = 3, col_names = NULL, ...)

# Default S3 method
to_table(object, ...)

# S3 method for class 'RiskMap_cross_validation'
to_table(object, digits = 3, col_names = NULL, reference = NULL, ...)

# S3 method for class 'summary.RiskMap_cross_validation'
to_table(object, digits = 3, col_names = NULL, ...)

# S3 method for class 'RiskMap_predictive_check'
to_table(object, digits = 3, col_names = NULL, ...)

# S3 method for class 'summary.RiskMap_predictive_check'
to_table(object, digits = 3, col_names = NULL, ...)
```

## Arguments

- object:

  A fit from
  [`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
  or its
  [summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md),
  an assessment of fits from
  [`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md)
  or its
  [summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_predictive_check.md),
  a cross-validation from
  [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
  or its
  [summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md),
  or a prediction over `boundaries` from
  [predict()](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
  (see
  [`to_table.RiskMap_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/to_table.RiskMap_prediction.md)).

- ...:

  Additional arguments passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html).

- digits:

  A non-negative integer giving the number of decimal places used to
  display numeric results. Small values get more (see Details).

- col_names:

  Optional column names, replacing the default headers.

- reference:

  For a cross-validation, optional name of the model that the others are
  compared with, as in its
  [summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md).

## Value

An object of class "knitr_kable" that can be rendered directly.

## Details

This function creates a presentation-ready summary table from a fitted
"RiskMap" model, an assessment of fits or cross-validation results for
multiple models. Use
[coef()](https://claudiofronterre.github.io/RiskMap/reference/coef.RiskMap.md)
or
[summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md)
when numeric results are required for further analysis.

Numeric values have `digits` decimal places. Values that are compared,
the estimate and confidence limits of a parameter or one metric across
models, share one notation and one number of decimal places. When one of
them is smaller in size than `10^(-digits)`, they get as many decimal
places as the smallest needs for two significant digits. They are in
scientific notation when that takes more than six decimal places, or
when one of them is a million or more in size. Zeros and missing values
don't count.

The table of a fit or its summary has

- the regression coefficients with their estimates and confidence
  intervals,

- the parameters of the spatial process,

- the variances of the random effects,

- the parameters of the family that were estimated, such as the
  measurement error variance of a Gaussian model.

The table of a cross-validation or its summary has a row for each model,
with its scores, their differences from a `reference` model with their
standard errors in brackets, and its `PIT_area` with the limit of a
correct model. A model of binary outcomes has "none (binary)" as its
limit, because the area is not a check of calibration there (see
[`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)).

## See also

[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`summary.RiskMap_cross_validation()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_cross_validation.md)

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
