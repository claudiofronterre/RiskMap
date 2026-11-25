##' Which function answers which question about a model
##'
##' @name assessment
##' @description
##' Several functions of RiskMap check or assess a model, and each answers a
##' different question. The first questions come before the model is fitted and
##' just after, and the last ones when it is used to predict.
##'
##' \tabular{ll}{
##'   **Question** \tab **Function** \cr
##'   Are the data spatially correlated? \tab [variogram()] \cr
##'   Did the computation of the fit work? \tab [check_fit()], and
##'     [plot_mcmc()] and [summarise_mcmc()] for its samples \cr
##'   How uncertain are the estimates? \tab [summary()][summary.RiskMap()],
##'     [profile()][profile.RiskMap()] and [parametric_bootstrap()] \cr
##'   Could the model have produced the data? \tab [assess_fit()] \cr
##'   How well does the model predict, and which of several models predicts
##'     better? \tab [assess_prediction()] \cr
##'   Where do the predictions go beyond what the data support? \tab
##'     [assess_extrapolation()] \cr
##'   How well would the predictions of a model do if the truth were known?
##'     \tab [assess_simulation()]
##' }
##'
##' @details
##' - [check_fit()] is about the computation and not about the model. It says
##'   whether the maximum of the likelihood was found, and how large the Monte
##'   Carlo error of the estimates is. A wrong model can pass it.
##' - [assess_fit()] compares the data with datasets simulated from the model, to
##'   see whether the model misses a feature of them. It says nothing about
##'   predictions, and a model that passes is not shown to be right.
##' - [assess_prediction()] leaves observations out, predicts them from the
##'   others and scores the predictions, by their sharpness and calibration. It
##'   compares models, and the model with the best scores may still miss a
##'   feature of the data.
##' - [assess_extrapolation()] needs no outcome at the prediction locations. It
##'   shows where their covariates are outside the range of the data and where
##'   the data say little about the spatial process.
##' - [assess_simulation()] needs datasets simulated with a known surface, from
##'   [simulate_glgpm()], and so it is for planning a survey or comparing
##'   methods, and not for assessing a fit to data.
##'
##' The order above is a useful one. A fit whose computation failed shouldn't be
##' assessed further, and the scores of a model that misses a feature of the data
##' are easier to understand when the feature is known.
##'
##' @seealso [glgpm()], [predict()][predict.RiskMap()] and [simulate_glgpm()].
NULL
