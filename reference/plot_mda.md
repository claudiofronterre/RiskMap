# Plot the estimated MDA impact function

Generate a plot of the estimated impact of mass drug administration
(MDA) on infection prevalence, based on a fitted decay-adjusted
spatio-temporal (DAST) model. The function simulates draws from the
posterior distribution of model parameters, propagates them through the
MDA effect function, and produces uncertainty bands around the estimated
impact curve.

## Usage

``` r
plot_mda(
  object,
  mda_history = NULL,
  n_sim = 1000,
  x_min = 1e-06,
  x_max = 10,
  conf_level = 0.95,
  lower_f = NULL,
  upper_f = NULL,
  mc_cores = 1,
  parallel_backend = c("none", "fork", "psock"),
  ...
)
```

## Arguments

- object:

  A fitted DAST model object, returned by
  [`dast`](https://claudiofronterre.github.io/RiskMap/reference/dast.md).

- mda_history:

  Specification of the MDA schedule. This can be either:

  - A numeric vector of event times (integers starting at 0, e.g.
    `c(0,1,2,6)`),

  - OR a 0/1 indicator vector on the yearly grid (e.g.
    `c(1,1,1,0,0,0,1)`), where position `i` corresponds to year `i-1`.

  If omitted, the default is a single MDA at time 0.

- n_sim:

  Number of posterior draws used for uncertainty quantification
  (default: 1000).

- x_min:

  Minimum value for the x-axis (default: `1e-6`).

- x_max:

  Maximum value for the x-axis (default: `10`).

- conf_level:

  Confidence level for the pointwise uncertainty interval (default:
  0.95).

- lower_f:

  Optional lower bound for the y-axis. If not provided, computed from
  the data.

- upper_f:

  Optional upper bound for the y-axis. If not provided, computed from
  the data.

- mc_cores:

  Number of CPU cores to use for parallel simulation. Default is 1
  (serial).

- parallel_backend:

  Parallelisation backend to use. Options are `"none"` (default),
  `"fork"` (Unix-like systems), or `"psock"` (cross-platform).

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot2` object showing the median estimated MDA impact function and
the pointwise uncertainty band at the chosen confidence level.

## Details

The time axis is assumed to start at 0 and increase in integer steps of
1 year. The argument `mda_history` allows the user to specify when MDAs
occurred either by listing the years directly or by giving a binary
indicator on the yearly grid. The function then evaluates the cumulative
relative reduction \\1 - \mathrm{effect}(t)\\ at a dense grid of time
points between `x_min` and `x_max`, using the fitted parameters from the
supplied DAST model.
