# Plot Method for RiskMap_pred_target_shp Objects

Generates a plot of predictive target values or summaries over a
shapefile.

## Usage

``` r
# S3 method for class 'RiskMap_pred_target_shp'
plot(x, which_target = "linear_target", which_summary = "mean", ...)
```

## Arguments

- x:

  An object of class 'RiskMap_pred_target_shp' containing computed
  targets, summaries, and associated spatial data.

- which_target:

  Character indicating the target type to plot (e.g., "linear_target").

- which_summary:

  Character indicating the summary type to plot (e.g., "mean", "sd").

- ...:

  Additional arguments passed to 'scale_fill_distiller' in 'ggplot2'.

## Value

A `ggplot` object showing the plot of the specified predictive target or
summary.

## Details

This function plots the predictive target values or summaries over a
shapefile. It requires the 'ggplot2' package for plotting and 'sf'
objects for spatial data.

## See also

[`pred_target_shp`](https://claudiofronterre.github.io/RiskMap/reference/pred_target_shp.md),
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html),
[`geom_sf`](https://ggplot2.tidyverse.org/reference/ggsf.html),
[`aes`](https://ggplot2.tidyverse.org/reference/aes.html),
[`scale_fill_distiller`](https://ggplot2.tidyverse.org/reference/scale_brewer.html)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
