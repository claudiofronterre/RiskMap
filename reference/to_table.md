# Generate LaTeX Tables from RiskMap Model Fits and Validation

Converts a fitted "RiskMap" model or cross-validation results into an
`xtable` object, formatted for easy export to LaTeX or HTML.

## Usage

``` r
to_table(object, ...)
```

## Arguments

- object:

  An object of class "RiskMap" resulting from a call to
  [`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
  or a summary object of class "summary.RiskMap.spatial.cv" containing
  cross-validation results.

- ...:

  Additional arguments to be passed to
  [`xtable`](https://rdrr.io/pkg/xtable/man/xtable.html) for
  customization.

## Value

An object of class "xtable", which contains the formatted table as a
`data.frame` and several attributes specifying table formatting options.

## Details

This function creates a summary table from a fitted "RiskMap" model or
cross-validation results for multiple models, returning it as an
`xtable` object.

When the input is a "RiskMap" model object, the table includes:

- Regression coefficients with their estimates, confidence intervals,
  and p-values.

- Parameters for the spatial process.

- Random effect variances.

- Measurement error variance, if applicable.

When the input is a cross-validation summary object
("summary.RiskMap.spatial.cv"), the table includes:

- A row for each model being compared.

- Performance metrics such as CRPS and SCRPS for each model.

The resulting `xtable` object can be further customized with additional
formatting options and printed as a LaTeX or HTML table for reports or
publications.

## See also

[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[`xtable`](https://rdrr.io/pkg/xtable/man/xtable.html),
[`summary.RiskMap.spatial.cv`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.spatial.cv.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
