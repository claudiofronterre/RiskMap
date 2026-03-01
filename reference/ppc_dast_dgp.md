# Internal Posterior Predictive Validation for DAST Fits

Internal utility to validate DAST model fit quality by comparing
observed outcomes to simulations generated from the fitted DAST
data-generating process. The function computes discrepancy summaries at
global, temporal, and IU aggregation levels, and optionally returns
diagnostic plots.

## Usage

``` r
ppc_dast_dgp(
  y_obs,
  m,
  y_sim,
  year,
  iu,
  m_bins = c(0, 23.5, 30.5, 49.5, Inf),
  iu_min_total_m = 50,
  iu_year_min_total_m = 25,
  prob = 0.9,
  n_overlay = 50,
  seed = 1L,
  make_plots = TRUE,
  prevalence_estimand = c("both", "pooled", "unweighted")
)
```

## Arguments

- y_obs:

  Numeric vector of observed positive counts.

- m:

  Numeric vector of denominators (examined counts).

- y_sim:

  Numeric matrix of simulated counts with dimension
  `length(y_obs) x nsim` (rows are observations, columns are simulation
  draws).

- year:

  Observation year (numeric/integer/factor/character coercible to
  numeric).

- iu:

  IU identifier vector of same length as `y_obs`.

- m_bins:

  Breaks used to create denominator groups for grouped PPC plots.

- iu_min_total_m:

  Minimum total examined count required to include an IU in IU-level
  summaries.

- iu_year_min_total_m:

  Minimum total examined count required to include an IU-year cell in
  IU-year summaries.

- prob:

  Central predictive interval probability used in coverage summaries and
  envelopes.

- n_overlay:

  Number of simulation draws to overlay in PPC plots.

- seed:

  Integer random seed used for randomized PIT and draw subsampling.

- make_plots:

  Logical; if `TRUE`, return PPC plots when plotting dependencies are
  available.

- prevalence_estimand:

  Which prevalence estimand to report and plot: `"pooled"`,
  `"unweighted"` or `"both"`.

## Value

An object of class `"ppc_dast_dgp"` with components:

- inputs:

  Input metadata and key settings.

- summaries:

  A nested list with `global`, `temporal`, and `iu` summary outputs.

- plots:

  A nested list with `global`, `temporal`, and `iu` diagnostic plots
  (possibly empty).

Year-level and IU-level prevalence trajectories use the
examined-weighted prevalence (pooled prevalence), computed as \\\sum_i
y_i / \sum_i m_i\\ within each aggregation stratum. This is different
from the unweighted mean of location-level prevalences.

For IU-level summaries and plots, an IU is included when it has at least
one observed record and total examined count \\\sum_i m_i \ge
iu\\min\\total\\m\\.

## Details

Interpretation guide:

- `summaries$global$pvals`: posterior predictive p-values (PPP) for
  global discrepancy statistics. Values close to 0 or 1 indicate
  mismatch between data and model-implied simulations for the
  corresponding summary.

- `summaries$global$pit_summary`: PIT mean/variance; departures from
  roughly Uniform(0,1) moments suggest calibration issues.

- `summaries$global$z_summary`: quantiles of standardized residuals
  \\(y_i-\mu_i)/\sqrt{v_i}\\ based on simulation moments. Heavy tails or
  skewness can indicate poor dispersion/shape fit.

- `summaries$temporal$pvals`: year-level trajectory mismatch summaries.

- `summaries$iu$pvals`: IU and IU-year support-aware mismatch summaries.

- `summaries$iu$iu_support`: IU-level support table (number of surveys,
  total examined, observed pooled/unweighted prevalence and inclusion
  flag).

- `summaries$iu$diagnostics`: counts summarising IU inclusion and zero
  vs non-zero observed records, including out-of-interval IU counts.

The two-sided PPP is computed as a median-centered tail area: \$\$PPP =
\frac{1}{S}\sum\_{s=1}^{S}
\mathbf{1}\left(\left\|T_s^{rep}-\mathrm{median}(T^{rep})\right\| \ge
\left\|T^{obs}-\mathrm{median}(T^{rep})\right\|\right)\$\$ This is
robust for asymmetric predictive distributions and does not assume
normality of discrepancy statistics.

When `make_plots = TRUE`, bayesplot is required for PPC plot primitives.
If bayesplot is unavailable, the function returns numeric summaries and
sets `out$plot_note` with a warning message. Prevalence values are
computed internally on the 0-1 scale but shown in plots as percentages
for readability.

The IU mean prevalence PPC plot is a forest-style diagnostic: observed
IU prevalence is shown as points, simulated central intervals are shown
as horizontal ranges, and observed point fill color encodes the number
of survey records contributing to each IU mean. A companion plot
`plots$iu$iu_mean_prev_misfit` shows only IUs where observed prevalence
falls outside the simulated central interval.

## Examples

``` r
if (FALSE) { # \dontrun{
out <- dast_sim(n_sim = 250, model_fit = dast_with_pen, messages = FALSE)
ppc_all <- ppc_dast_dgp(
  y_obs = dast_with_pen$data_sf$positive,
  m = dast_with_pen$data_sf$examined,
  y_sim = out$data_sim,
  year = dast_with_pen$data_sf$year,
  iu = dast_with_pen$data_sf$IUs_NAME,
  prob = 0.9
)
print(ppc_all)
ppc_all$plots$temporal$year_traj
ppc_all$plots$temporal$year_diffs
ppc_all$plots$iu$iu_mean_prev
ppc_all$plots$iu$iu_mean_prev_misfit
ppc_all$plots$iu$iu_year_cell_ecdf
ppc_all$plots$global$pit_ecdf
ppc_all$plots$global$stat_hist$binom_deviance_plugin_mean
} # }
```
