% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_prediction.R
\name{surf_sim}
\alias{surf_sim}
\title{Simulate surface data based on a spatial model}
\usage{
surf_sim(
  n_sim,
  pred_grid,
  formula,
  sampling_f,
  family,
  scale_to_km = TRUE,
  control_mcmc = set_control_sim(),
  par0,
  nugget_over_grid = FALSE,
  include_covariates = TRUE,
  fix_var_me = NULL,
  messages = TRUE
)
}
\arguments{
\item{n_sim}{The number of simulations to run.}

\item{pred_grid}{A spatial object (either `sf` or `data.frame`) representing the prediction grid where the simulation will take place.}

\item{formula}{A formula object specifying the model to be fitted. It should include both fixed effects and random effects if applicable.}

\item{sampling_f}{A function that returns a sampled dataset (of class `sf` or `data.frame`) to simulate data from.}

\item{family}{A character string specifying the family of the model. Must be one of "gaussian", "binomial", or "poisson".}

\item{scale_to_km}{A logical indicating whether the coordinates should be scaled to kilometers. Defaults to `TRUE`.}

\item{control_mcmc}{A list of control parameters for MCMC (not used in this implementation but can be expanded later).}

\item{par0}{A list containing initial parameter values for the simulation, including `beta`, `sigma2`, `phi`, `tau2`, and `sigma2_me`.}

\item{nugget_over_grid}{A logical indicating whether to include a nugget effect over the entire prediction grid.}

\item{include_covariates}{A logical indicateing if the covariates (or the intercept if no covariates are used) should be included in the linear
predictor. By default \code{include_covariates = TRUE}}

\item{fix_var_me}{A parameter to fix the variance of the random effects for the measurement error. Defaults to `NULL`.}

\item{messages}{A logical value indicating whether to print messages during the simulation. Defaults to `TRUE`.}
}
\value{
A list containing the simulated data (\code{data_sim}), the linear predictors (\code{lp_grid_sim}),
a logical value indicating if covariates have been included in the linear predictor (\code{include_covariates}),
a logical value indicating if the nugget has been included into the simulations of the linear predictor over the grid
(\code{nugget_over_grid}), a logical  indicating if a covariate offset has been included in the linear predictor (\code{include_cov_offset}),
the model parameters set for the simulation (\code{par0}) and the family used in the model (\code{family}).
}
\description{
This function simulates surface data based on a user-defined formula and other parameters. It allows for simulation of spatial data with various model families (Gaussian, Binomial, or Poisson). The simulation involves creating spatially correlated random fields and generating outcomes for data points in a given prediction grid.
}
\author{
Emanuele Giorgi \email{e.giorgi@lancaster.ac.uk}
}
