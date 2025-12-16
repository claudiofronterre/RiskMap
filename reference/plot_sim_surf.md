# Plot simulated surface data for a given simulation

This function plots the simulated surface data for a specific simulation
from the result of \`surf_sim\`. It visualizes the linear predictor
values on a raster grid along with the actual data points.

## Usage

``` r
plot_sim_surf(surf_obj, sim, ...)
```

## Arguments

- surf_obj:

  The output object from \`surf_sim\`, containing both simulated data
  (\`data_sim\`) and predicted grid simulations (\`lp_grid_sim\`).

- sim:

  The simulation index to plot.

- ...:

  Additional graphical parameters to be passed to the plotting function
  of the \`terra\` package.

## Value

A plot of the simulation results.
