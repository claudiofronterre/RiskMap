# Plot Method for RiskMap_pred_target_grid Objects

Generates a plot of the predicted values or summaries over the regular
spatial grid from an object of class 'RiskMap_pred_target_grid'.

## Usage

``` r
# S3 method for class 'RiskMap_pred_target_grid'
plot(x, which_target = "linear_target", which_summary = "mean", ...)
```

## Arguments

- x:

  An object of class 'RiskMap_pred_target_grid'.

- which_target:

  Character string specifying which target prediction to plot.

- which_summary:

  Character string specifying which summary statistic to plot (e.g.,
  "mean", "sd").

- ...:

  Additional arguments passed to the
  [`plot`](https://rspatial.github.io/terra/reference/plot.html)
  function of the `terra` package.

## Value

A `ggplot` object representing the specified prediction target or
summary statistic over the spatial grid.

## Details

This function requires the 'terra' package for spatial data manipulation
and plotting. It plots the values or summaries over a regular spatial
grid, allowing for visual examination of spatial patterns.

## See also

[`pred_target_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_target_grid.md)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
