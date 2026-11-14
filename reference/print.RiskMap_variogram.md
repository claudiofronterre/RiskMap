# Print Method for RiskMap_variogram Objects

Prints the distance classes and the envelope settings of an empirical
variogram computed by
[`variogram`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md),
followed by its table.

## Usage

``` r
# S3 method for class 'RiskMap_variogram'
print(x, digits = 3, ...)
```

## Arguments

- x:

  an object of class `RiskMap_variogram`.

- digits:

  the number of significant digits to print. By default `digits = 3`.

- ...:

  further arguments passed to
  [`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html).

## Value

`x`, invisibly.

## See also

[`variogram`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md)
