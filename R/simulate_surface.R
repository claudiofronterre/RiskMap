##' Simulate surface data based on a spatial model
##'
##' This function simulates surface data based on a user-defined formula and other parameters. It allows for simulation of spatial data with various model families (Gaussian, Binomial, or Poisson). The simulation involves creating spatially correlated random fields and generating outcomes for data points in a given prediction grid.
##'
##' @param n_sim The number of simulations to run.
##' @param pred_grid An `sf` object representing the prediction grid where the simulation will take place.
##' @param formula A formula object specifying the model to be fitted. It should include both fixed effects and random effects if applicable.
##' @param sampling_f A function that returns a sampled dataset (of class `sf`) to simulate data from.
##' @param family A character string specifying the family of the model. Must be one of "gaussian", "binomial", or "poisson".
##' @param scale_to_km A logical indicating whether the coordinates should be scaled to kilometers. Defaults to `TRUE`.
##' @param control_mcmc A list of control parameters for MCMC (not used in this implementation but can be expanded later).
##' @param par0 A list containing initial parameter values for the simulation, including `beta`, `sigma2`, `phi`, `tau2`, and `sigma2_me`.
##' @param include_covariates A logical indicating if the covariates (or the intercept if no covariates are used) should be included in the linear
##' predictor. By default \code{include_covariates = TRUE}
##' @param nugget_over_grid A logical indicating whether to include a nugget effect over the entire prediction grid.
##' @param fix_var_me A parameter to fix the variance of the random effects for the measurement error. Defaults to `NULL`.
##' @param messages A logical value indicating whether to print messages during the simulation. Defaults to `TRUE`.
##'
##' @return A list containing the simulated data (\code{data_sim}), the linear predictors (\code{lp_grid_sim}),
##' a logical value indicating if covariates have been included in the linear predictor (\code{include_covariates}),
##' a logical value indicating if the nugget has been included into the simulations of the linear predictor over the grid
##' (\code{nugget_over_grid}), a logical  indicating if a covariate offset has been included in the linear predictor (\code{include_cov_offset}),
##' the model parameters set for the simulation (\code{par0}) and the family used in the model (\code{family}).
##'
##' @export
simulate_surface <- function(n_sim,
                             pred_grid,
                             formula,
                             sampling_f,
                             family,
                             scale_to_km = TRUE,
                             control_mcmc = set_control_mcmc(),
                             par0, nugget_over_grid = FALSE,
                             include_covariates = TRUE,
                             fix_var_me = NULL,
                             messages = TRUE) {


  if(!inherits(formula,
               what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula'
                                     object indicating the variables of the
                                     model to be fitted")
  }
  inter_f <- interpret.formula(formula)
  include_cov_offset <- !is.null(inter_f$offset)
  if(!inherits(pred_grid, "sf")) {
    stop("'pred_grid' must be an 'sf'
          object indicating the variables of the
          model to be fitted")
  }

  sim_crs <- st_crs(pred_grid)$epsg
  if (is.na(sim_crs)) {
    stop("'pred_grid' must have a valid coordinate reference system (CRS) set.")
  }

  if(!inherits(sampling_f, "function")){
    stop("'sampling_f' must be an object of class 'function'")
  }

  data_test <- sampling_f()
  if(!inherits(data_test, "sf")) {
    stop("The object return by 'sampling_f' must be of an 'sf' object")
  }

  if (!"units_m" %in% colnames(data_test)) {
    stop("The object returned by 'sampling_f' must contain a column named 'units_m'.")
  }


  data_sim <- list()
  coords_sim <- list()
  for(i in 1:n_sim) {
    data_sim[[i]] <- sampling_f()
    coords_sim[[i]] <- st_coordinates(data_sim[[i]])
    if(scale_to_km) coords_sim[[i]] <- coords_sim[[i]]/1000
    if (st_crs(data_sim[[i]]) != st_crs(pred_grid)) {
      pred_grid <- st_transform(pred_grid, st_crs(data_sim[[i]]))
      if(i==1) {
        sim_crs_data <- st_crs(data_sim[[i]])
      }
    } else {
      if(i==1) {
        sim_crs_data <- sim_crs
      }
    }

    # Find nearest neighbor in 'pred_grid' for each feature in 'data'
    nearest_indices <- st_nearest_feature(data_sim[[i]], pred_grid)

    # Extract variables from the nearest features in 'pred_grid'
    pred_grid_vars <- st_drop_geometry(pred_grid[nearest_indices, ])

    # Bind the extracted variables to 'data'
    data_sim[[i]] <- cbind(data_sim[[i]], pred_grid_vars)
  }

  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")

  if(family != "gaussian" & family != "binomial" &
     family != "poisson") stop("'family' must be either 'gaussian', 'binomial'
                               or 'poisson'")

  inter_lt_f <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~.)
  D <- list()
  n <- list()
  for(i in 1:n_sim) {
    mf <- model.frame(inter_lt_f$pf,data=data_sim[[i]], na.action = na.fail)


    # Extract covariates matrix
    D[[i]] <- as.matrix(model.matrix(attr(mf,"terms"),data=data_sim[[i]]))
    n[[i]] <- nrow(D[[i]])
  }


  if(length(inter_f$re.spec) > 0) {
    stop("In the current impletementation of 'simulate_surface' the addition of random effects
        with re() is not supported")
  }
  mf_pred <- model.frame(inter_lt_f$pf,data=pred_grid, na.action = na.fail)
  D_pred <- as.matrix(model.matrix(attr(mf_pred,"terms"),data=pred_grid))
  if(ncol(D_pred)!=ncol(D[[1]])) stop("the provided variables in 'grid_pred' do not match the number of explanatory variables used to fit the model.")

  # Number of covariates
  p <- ncol(D[[1]])

  beta <- par0$beta

  if(length(beta)!=p) stop("The values passed to 'beta' do not match the variables in 'grid_pred'")

  sigma2 <- par0$sigma2
  phi <- par0$phi

  if(is.null(par0$tau2)) {
    tau2 <- 0
  } else {
    tau2 <- par0$tau2
  }

  if(is.null(fix_var_me)) {
    sigma2_me <- 0
  } else {
    sigma2_me <- par0$sigma2_me
  }

  ID_coords <- sapply(1:n_sim, function(i) 1:n[[i]])

  grid_pred <- st_coordinates(pred_grid)
  if(scale_to_km) grid_pred <- grid_pred/1000
  n_pred <- nrow(grid_pred)

  # Simulate on the grid
  # Simulate S
  Sigma <- sigma2*matern_correlation(dist(grid_pred), phi = phi, kappa = kappa,
                                     return_sym_matrix = TRUE)
  if(nugget_over_grid) {
    diag(Sigma) <- diag(Sigma) + tau2
  }
  Sigma_sroot <- t(chol(Sigma))

  S_sim <- sapply(1:n_sim, function(i) Sigma_sroot%*%rnorm(n_pred))

  sim_columns <- paste0("sim_", 1:n_sim)
  S_sim_df <- as.data.frame(S_sim)
  colnames(S_sim_df) <- sim_columns

  # Combine simulations with grid_pred
  S_sim <- cbind(grid_pred, S_sim_df)
  S_sim <- st_sf(S_sim, geometry = st_geometry(pred_grid))

  coords_tot <- list()
  out <- list()

  for(i in 1:n_sim) {
    coords_tot[[i]] <- rbind(coords_sim[[i]], grid_pred)
    ind_data <- 1:n[[i]]
    ind_pred <- (n[[i]]+1):(n[[i]]+n_pred)

    D_tot <- rbind(D[[i]], D_pred)

    # Find the nearest features in grid_pred_with_sim for each location in data_sim[[i]]
    nearest_indices <- st_nearest_feature(data_sim[[i]],
                                          S_sim)

    # Extract the i-th simulation values
    sim_column <- paste0("sim_", i)  # Column name for the i-th simulation
    S_sim_data <- S_sim[nearest_indices, sim_column, drop = TRUE]

    if(tau2>0 & !nugget_over_grid) {
      S_sim_data <- S_sim_data + sqrt(tau2)*rnorm(n[[i]])
    }
    # Linear predictor
    if(include_covariates) {
      eta_sim_tot <- D_tot%*%beta +
        c(S_sim_data, S_sim[[sim_column]])
    } else {
      eta_sim_tot <- c(S_sim_data, S_sim[[sim_column]])
    }

    data_sim[[i]]$y <- NA

    if(family=="gaussian") {

      data_sim[[i]]$y <- eta_sim_tot[1:n[[i]]] + sqrt(sigma2_me)*rnorm(n)

    } else if(family=="binomial") {
      prob_i <- exp(eta_sim_tot[1:n[[i]]])/(1+exp(eta_sim_tot[1:n[[i]]]))
      data_sim[[i]]$y <- rbinom(n[[i]], size = data_sim[[i]]$units_m,
                                prob = prob_i)
    } else if(family=="poisson") {
      mean_i <- data_sim[[i]]$units_m*exp(eta_sim_tot[1:n[[i]]])
      data_sim[[i]]$y <- rpois(n[[i]], lambda = mean_i)
    }

    data_sim[[i]]$lp_data <- S_sim_data

    pred_grid[[paste0("lp_",sim_column)]] <- eta_sim_tot[-(1:n[[i]])]

  }

  out$data_sim <- data_sim
  out$crs <- sim_crs_data
  out$lp_grid_sim <- pred_grid
  out$include_covariates <- include_covariates
  out$nugget_over_grid <- nugget_over_grid
  out$include_cov_offset <- include_cov_offset
  out$par0 <- par0
  out$family <- family
  class(out) <- "RiskMap_simulation"
  return(out)
}

##' Plot simulated surface data for a given simulation
##'
##' This function plots the simulated surface data for a specific simulation from the result of `simulate_surface`. It visualizes the linear predictor values on a raster grid along with the actual data points.
##'
##' @param surf_obj The output object from `simulate_surface`, containing both simulated data (`data_sim`) and predicted grid simulations (`lp_grid_sim`).
##' @param sim The simulation index to plot.
##' @param ... Additional graphical parameters to be passed to the plotting function of the `terra` package.
##'
##' @return A plot of the simulation results.
##'
##' @importFrom stars st_rasterize
##'
##' @export
plot_sim_surf <-  function(surf_obj, sim, ...) {

  sf_object <- surf_obj$lp_grid_sim
  value_column <- paste0("lp_sim_",sim)
  r <- rast(st_rasterize(sf_object[,c("x",value_column)]))
  r[r == 0] <- NA

  plot(r, main = paste("Simulation no.", sim), ...)
  points(st_coordinates(surf_obj$data_sim[[sim]]), pch = 20)

}
