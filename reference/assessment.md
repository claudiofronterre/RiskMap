# Which function answers which question about a model

Several functions of RiskMap check or assess a model, and each answers a
different question. The first questions come before the model is fitted
and just after, and the last ones when it is used to predict.

|  |  |
|----|----|
| **Question** | **Function** |
| Are the data spatially correlated? | [`variogram()`](https://claudiofronterre.github.io/RiskMap/reference/variogram.md) |
| Did the computation of the fit work? | [`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md), and [`plot_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/plot_mcmc.md) and [`summarise_mcmc()`](https://claudiofronterre.github.io/RiskMap/reference/summarise_mcmc.md) for its samples |
| How uncertain are the estimates? | [summary()](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap.md), [profile()](https://claudiofronterre.github.io/RiskMap/reference/profile.RiskMap.md) and [`parametric_bootstrap()`](https://claudiofronterre.github.io/RiskMap/reference/parametric_bootstrap.md) |
| Could the model have produced the data? | [`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md) |
| How well does the model predict, and which of several models predicts better? | [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md) |
| Where do the predictions go beyond what the data support? | [`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md) |
| How well would the predictions of a model do if the truth were known? | [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md) |

## Details

- [`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md)
  is about the computation and not about the model. It says whether the
  maximum of the likelihood was found, and how large the Monte Carlo
  error of the estimates is. A wrong model can pass it.

- [`assess_fit()`](https://claudiofronterre.github.io/RiskMap/reference/assess_fit.md)
  compares the data with datasets simulated from the model, to see
  whether the model misses a feature of them. It says nothing about
  predictions, and a model that passes is not shown to be right.

- [`assess_prediction()`](https://claudiofronterre.github.io/RiskMap/reference/assess_prediction.md)
  leaves observations out, predicts them from the others and scores the
  predictions, by their sharpness and calibration. It compares models,
  and the model with the best scores may still miss a feature of the
  data.

- [`assess_extrapolation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_extrapolation.md)
  needs no outcome at the prediction locations. It shows where their
  covariates are outside the range of the data and where the data say
  little about the spatial process.

- [`assess_simulation()`](https://claudiofronterre.github.io/RiskMap/reference/assess_simulation.md)
  needs datasets simulated with a known surface, from
  [`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md),
  and so it is for planning a survey or comparing methods, and not for
  assessing a fit to data.

The order above is a useful one. A fit whose computation failed
shouldn't be assessed further, and the scores of a model that misses a
feature of the data are easier to understand when the feature is known.

## See also

[`glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md),
[predict()](https://claudiofronterre.github.io/RiskMap/reference/predict.RiskMap.md)
and
[`simulate_glgpm()`](https://claudiofronterre.github.io/RiskMap/reference/simulate_glgpm.md).
