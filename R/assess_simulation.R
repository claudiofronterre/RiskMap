##' @title Assess Simulations
##'
##' @description This function evaluates the performance of models based on simulation results from the `simulate_surface` function.
##'
##' @param obj_sim An object of class `RiskMap_simulation`, obtained as an output from the `simulate_surface` function.
##' @param models A named list of models to be evaluated.
##' @param control_mcmc A control object for MCMC sampling, created with `set_control_mcmc()`. Default is `set_control_mcmc()`.
##' @param spatial_scale The scale at which predictions are assessed, either `"grid"` or `"area"`.
##' @param messages Logical, if `TRUE` messages will be displayed during processing. Default is `TRUE`.
##' @param f_grid_target A function for processing grid-level predictions.
##' @param f_area_target A function for processing area-level predictions.
##' @param shp A shapefile of class `sf` for area-level analysis, required if `spatial_scale = "area"`.
##' @param col_names Column name in `shp` containing unique region names. If `NULL`, defaults to `"region"`.
##' @param pred_objective A character vector specifying objectives, either `"mse"`, `"classify"`, or both.
##' @param categories A numeric vector of thresholds defining categories for classification. Required if `pred_objective = "classify"`.
##'
##' @return A list of class `RiskMap_assess_simulation` containing model evaluation results.
##'
##' @export
assess_simulation <- function(obj_sim,
                       models,
                       control_mcmc = set_control_mcmc(),
                       spatial_scale,
                       messages = TRUE,
                       f_grid_target = NULL,
                       f_area_target = NULL,
                       shp = NULL, col_names = NULL,
                       pred_objective = c("mse","classify"),
                       categories= NULL) {

  if (!inherits(obj_sim, "RiskMap_simulation")) {
    stop("'obj_sim' must be an object of class 'RiskMap_simulation' obtained as an output from the 'simulate_surface' function")
  }
  if (length(setdiff(pred_objective, c("mse","classify")))>0) {
    stop(paste("Invalid value for pred_objective. Allowed values are:", paste(c("mse","classify"), collapse = ", ")))
  }
  if(spatial_scale != "grid" & spatial_scale != "area") {
    stop("'spatial_scale' must be set to 'grid' or 'area'")
  }
  if(spatial_scale=="area" & is.null(shp)) {
    stop("if spatial_scale='area' then a shape file of the area(s) must be passed to
         'shp'")
    check_data(shp, "polygon")
  }

  # Determine the binomial denominator column, if relevant to the family
  units_m <- NULL
  if (obj_sim$family == "binomial") {
    stopifnot(!is.null(obj_sim$data_sim[[1]]$units_m))
    units_m <- obj_sim$data_sim[[1]]$units_m
  }
  if(any(pred_objective=="classify")) {
    if(is.null(categories)) stop("if pred_objective='class', a value for 'categories' must be specified")
    if (length(categories) < 3) {
      stop("Categories vector must contain at least three unique, strictly increasing values.")
    }
  }
  n_sim <- length(obj_sim$data_sim)
  n_models <- length(models)

  if(spatial_scale == "area" & is.null(f_area_target)) {
    stop("If 'spatial_scale' is set to 'area', then 'f_area_target' must be provided")
  }
  model_names <- names(models)

  include_covariates <- obj_sim$include_covariates
  include_cov_offset <- obj_sim$include_cov_offset
  include_nugget <- obj_sim$nugget_over_grid

  fits <- list()
  preds <- list()

  no_comp <- NULL

  if(spatial_scale=="grid") {
    type <- "marginal"
  } else if(spatial_scale=="area") {
    type <- "joint"
    n_reg <- nrow(shp)
    if(is.null(shp)) stop("If spatial_scale='area', then 'shp' must be specified")
    if(!inherits(shp,
                 what = c("sf"), which = FALSE)) {
      stop("The object passed to 'shp' must be an object of class 'sf'")
    }

    if(is.null(col_names)) {
      shp$region <- paste("reg",1:n_reg, sep="")
      col_names <- "region"
      names_reg <- shp$region
    } else {
      names_reg <- shp[[col_names]]
      if(n_reg != length(names_reg)) {
        stop("The names in the column identified by 'col_names' do not
         provide a unique set of names, but there are duplicates")
      }
    }
    shp <- st_transform(shp, st_crs(obj_sim$lp_grid_sim))
    inter <- st_intersects(shp, obj_sim$lp_grid_sim)
  }

  for(i in 1:n_models) {
    if(messages) message("Model: ", paste(model_names[i]),"\n")

    if_i <- interpret.formula(models[[i]])
    rhs_terms <- attr(terms(if_i$pf), "term.labels")
    # Check if there are any covariates
    if (length(rhs_terms) == 0) {
      predictors_i <- NULL
    } else {
      predictors_i <- obj_sim$lp_grid_sim
    }
    for(j in 1:n_sim) {
      if(messages) message("Processing simulation no.", j)
      f_i <- update(models[[i]], y ~ .)
      if(messages) message("Estimation")
      fits[[paste(model_names[i])]][[j]] <- glgpm(formula = f_i,
                                                  den = units_m,
                                                  family = obj_sim$family,
                                                  data = obj_sim$data_sim[[j]],
                                                  control_mcmc = control_mcmc,
                                                  messages = FALSE)

      if(messages) message("Prediction over the grid")
      preds[[paste(model_names[i])]][[j]] <-
        setup_prediction(fits[[paste(model_names[i])]][[j]],
                       grid_pred = st_as_sfc(obj_sim$lp_grid_sim),
                       predictors = predictors_i,
                       type = type, messages = FALSE)
    }
  }

  n_samples <- (control_mcmc$n_sim-control_mcmc$burnin)/control_mcmc$thin
  n_pred <- nrow(obj_sim$lp_grid_sim)


  out <- list(pred_objective = list())

  if(any(pred_objective=="mse")) {
    out$pred_objective$mse <- array(NA, c(n_models, n_sim))
    rownames(out$pred_objective$mse) <- model_names
    colnames(out$pred_objective$mse) <- paste0("sim_",1:n_sim)
  }

  if(any(pred_objective=="classify")) {

    # Ensure categories are unique and strictly increasing
    categories <- unique(sort(categories))


    # Assign classification to the output object
    out$pred_objective$classify <- setNames(vector("list", length(model_names)), model_names)

    # Correctly generate breaks and labels
    breaks <- categories  # Use categories directly as breaks
    categories_class <- factor(paste0("(", utils::head(categories, -1), ",",
                                      categories[-1], "]"))  # Labels to match intervals



    for(i in 1:n_models) {
      out$pred_objective$classify[[model_names[i]]] <- list(by_cat = list(),
                                                            across_cat = list())
      out$pred_objective$classify[[model_names[i]]]$by_cat <- vector("list", n_sim)
      for(j in 1:n_sim) {
        out$pred_objective$classify[[paste(model_names[i])]]$by_cat[[j]] <-
          data.frame(
            Class = categories_class,
            Sensitivity = NA,
            Specificity = NA,
            PPV = NA,
            NPV = NA,
            CC = NA
          )
      }
      out$pred_objective$classify[[model_names[i]]]$CC <- rep(NA,n_sim)
    }
  }
  lp_true_sim <- st_drop_geometry(obj_sim$lp_grid_sim[, grepl("lp_sim",
                                                              names(obj_sim$lp_grid_sim))])

  if(spatial_scale == "grid") {
    true_target_sim <- f_grid_target(lp_true_sim)
  } else if(spatial_scale == "area") {
    true_target_sim <- matrix(NA, nrow = n_reg, ncol = n_sim)
    true_target_grid_sim <- f_grid_target(lp_true_sim)
    for(i in 1:n_reg) {
      for(j in 1:n_sim) {
        if(length(inter[[i]])==0) {
          warning(paste("No points on the grid fall within", shp[[col_names]][i],
                        "and no predictions are carried out for this area"))
          no_comp <- c(no_comp, i)
        } else {
          true_target_sim[i,j] <- f_area_target(true_target_grid_sim[inter[[i]],j])
        }
      }
    }
  }



  for(i in 1:n_models) {
    for(j in 1:n_sim) {
      obj_pred_ij <- preds[[paste(model_names[i])]][[j]]
      if(length(obj_pred_ij$mu_pred)==1 && obj_pred_ij$mu_pred==0 &&
         include_covariates) {
        stop("Covariates have not been provided; re-run setup_prediction
         and provide the covariates through the argument 'predictors'")
      }


      if(!include_covariates) {
        mu_target <- 0
      } else {

        if(is.null(obj_pred_ij$mu_pred)) stop("the output obtained from 'setup_prediction' does not
                                     contain any covariates; if including covariates
                                     in the predictive target these should be included
                                     when running 'setup_prediction'")
        mu_target <- obj_pred_ij$mu_pred
      }

      if(!include_cov_offset) {
        cov_offset <- 0
      } else {
        if(length(obj_pred_ij$cov_offset)==1) {
          stop("No covariate offset was included in the model;
           set include_cov_offset = FALSE, or refit the model and include
           the covariate offset")
        }
        cov_offset <- obj_pred_ij$cov_offset
      }

      if(include_nugget) {
        if(is.null(obj_pred_ij$par_hat$tau2)) stop("'include_nugget' cannot be
                                                   set to TRUE if this has not been included
                                                   in the fit of the model")
        Z_sim <- matrix(rnorm(n_samples*n_pred,
                              sd = sqrt(obj_pred_ij$par_hat$tau2)),
                        ncol = n_samples)
        obj_pred_ij$S_samples <- obj_pred_ij$S_samples+Z_sim
      }

      if(is.matrix(mu_target)) {
        lp_samples_ij <- sapply(1:n_samples,
                                function(h)
                                  mu_target[,h] + cov_offset +
                                  obj_pred_ij$S_samples[,h])
      } else {
        lp_samples_ij <- sapply(1:n_samples,
                                function(h)
                                  mu_target + cov_offset +
                                  obj_pred_ij$S_samples[,h])
      }

      target_samples_ij <- f_grid_target(lp_samples_ij)

      if(spatial_scale == "grid") {
        mean_target_ij <- apply(target_samples_ij, 1, mean)
      } else if(spatial_scale == "area") {
        target_area_samples_ij <- matrix(NA, nrow = n_reg, ncol = n_samples)
        mean_target_ij <- rep(NA,n_reg)
        for(h in 1:n_reg) {
          if(length(inter[[h]])==0) {
            warning(paste("No points on the grid fall within", shp[[col_names]][h],
                          "and no predictions are carried out for this area"))
            no_comp <- c(no_comp, h)
          } else {
            ind_grid_h <- inter[[h]]
            target_area_samples_ij[h,] <-  apply(target_samples_ij[ind_grid_h,], 2,
                                                 f_area_target)
            mean_target_ij[h] <- mean(target_area_samples_ij[h,])
          }
        }
      }

      if(any(pred_objective=="mse")) {
        out$pred_objective$mse[i,j] <- mean((mean_target_ij-true_target_sim[,j])^2)
      }

      if(any(pred_objective=="classify")) {
        true_class_ij <- cut(true_target_sim[,j], breaks = categories)
        n_categories <- length(categories)-1
        if(spatial_scale == "grid") {
          prob_cat_ij <- matrix(0, nrow=n_pred, ncol = n_categories)
        } else if(spatial_scale == "area") {
          prob_cat_ij <- matrix(0, nrow=n_reg, ncol = n_categories)
        }
        for(h in 1:(n_categories)) {
          if(spatial_scale == "grid") {
            prob_cat_ij[,h] <- apply(categories[h] < target_samples_ij &
                                       categories[h+1] > target_samples_ij, 1, mean)
          } else if(spatial_scale == "area") {
            prob_cat_ij[,h] <- apply(categories[h] < target_area_samples_ij &
                                       categories[h+1] > target_area_samples_ij, 1, mean)
          }
        }

        pred_class_ij <- apply(prob_cat_ij, 1, function(x) categories_class[which.max(x)])

        # Define the confusion matrix
        conf_matrix <- table(true_class_ij, pred_class_ij)

        # Calculate metrics for each class
        for (h in 1:nrow(conf_matrix)) {
          TP <- conf_matrix[h, h]  # True Positives: diagonal entry
          FP <- sum(conf_matrix[, h]) - conf_matrix[h, h]  # False Positives: column sum minus diagonal
          FN <- sum(conf_matrix[h, ]) - conf_matrix[h, h]  # False Negatives: row sum minus diagonal
          TN <- sum(conf_matrix) - sum(conf_matrix[h, ]) - sum(conf_matrix[, h]) + conf_matrix[h, h]


          # Handle cases where division by zero could occur
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$Sensitivity <-
            ifelse((TP + FN) == 0, NA, TP / (TP + FN))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$Specificity <-
            ifelse((TN + FP) == 0, NA, TN / (TN + FP))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$PPV <-
            ifelse((TP + FP) == 0, NA, TP / (TP + FP))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$NPV <-
            ifelse((TN + FN) == 0, NA, TN / (TN + FN))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$CC <-
            ifelse((TP + FN) == 0, NA, (TP ) / (TP + FN))
        }
        out$pred_objective$classify[[paste(model_names[[i]])]]$CC[j] <-
          mean(true_class_ij==pred_class_ij)
      }
    }
  }
  if(any(pred_objective=="classify")) out$pred_objective$classify$Class <- categories_class
  class(out) <- "RiskMap_assess_simulation"
  return(out)
}

##' @title Summarize Simulation Results
##'
##' @description Summarizes the results of model evaluations from a `RiskMap_assess_simulation` object. Provides average metrics for classification by category and overall correct classification (CC) summary.
##'
##' @param object An object of class `RiskMap_assess_simulation`, as returned by `assess_simulation`.
##' @param ... Additional arguments (not used).
##'
##' @return A list containing summary data for each model:
##' - `by_cat_summary`: A data frame with average sensitivity, specificity, PPV, NPV, and CC by category.
##' - `CC_summary`: A numeric vector with mean, 2.5th percentile, and 97.5th percentile for CC across simulations.
##'
##' @method summary RiskMap_assess_simulation
##' @export
summary.RiskMap_assess_simulation <- function(object, ...) {
  stopifnot(inherits(object, "RiskMap_assess_simulation"))

  # Initialize results
  results <- list()

  # Check for "mse" in pred_objective
  if ("mse" %in% names(object$pred_objective)) {
    mse_data <- object$pred_objective$mse

    # Check if mse_data is a matrix
    if (is.matrix(mse_data)) {
      # Compute mean and SD for each model (row)
      mse_summary <- data.frame(
        Model = rownames(mse_data),
        MSE_mean = rowMeans(mse_data, na.rm = TRUE),
        MSE_sd = apply(mse_data, 1, sd, na.rm = TRUE)
      )

      results$mse <- mse_summary
    } else {
      stop("mse_data must be a matrix.")
    }
  }

  # Check for "classify" in pred_objective
  if ("classify" %in% names(object$pred_objective)) {
    classify_data <- object$pred_objective$classify

    # Loop over each model (e.g., M1, M2)
    n_models <- length(classify_data)-1
    name_models <- names(classify_data)[1:n_models]
    results$classify <- list()
    for(i in 1:n_models) {

      model_data <- classify_data[[i]]

      n_sim <- length(model_data$by_cat)
      res_class <- model_data$by_cat[[1]][,-1]
      den <- 0
      for(j in 2:n_sim) {
        if(!any(is.na(model_data$by_cat[[j]][,-1]))) {
          den <- den + 1
          res_class <- res_class+model_data$by_cat[[j]][,-1]
        }
      }
      res_class <- data.frame(res_class/den)
      res_class$Class <- model_data$by_cat[[1]][,1]

      cc_summary <- list(mean = mean(model_data$CC, na.rm = TRUE),
                         lower = quantile(model_data$CC, 0.025, na.rm = TRUE),
                         upper = quantile(model_data$CC, 0.975, na.rm = TRUE))

      results$classify[[paste(name_models[i])]] <- list(classify_res = res_class,
                                            cc_summary = list(mean = mean(model_data$CC, na.rm = TRUE),
                                                                     lower = quantile(model_data$CC, 0.025, na.rm = TRUE),
                                                                     upper = quantile(model_data$CC, 0.975, na.rm = TRUE)))
    }
  }

  # Assign class for S3 print method
  class(results) <- "summary.RiskMap_assess_simulation"
  return(results)
}



##' @title Print Simulation Results
##'
##' @description Prints a concise summary of simulation results from a `RiskMap_assess_simulation` object, including average metrics by category and a summary of overall correct classification (CC).
##'
##' @param x An object of class `summary.RiskMap_assess_simulation`, as returned by `summary.RiskMap_assess_simulation`.
##' @param ... Additional arguments (not used).
##'
##' @return Invisibly returns `x`.
##'
##'
##' Print Simulation Results
##'
##' Prints a concise summary of simulation results from a `summary.RiskMap_assess_simulation` object,
##' including average metrics by category and a summary of overall correct classification (CC).
##'
##' @param x An object of class `summary.RiskMap_assess_simulation`, as returned by `summary.RiskMap_assess_simulation`.
##' @param ... Additional arguments (not used).
##'
##' @return Invisibly returns `x`.
##'
##' @method print summary.RiskMap_assess_simulation
##' @export
print.summary.RiskMap_assess_simulation <- function(x, ...) {
  cat("Summary of Simulation Results\n\n")

  if (!is.null(x$mse)) {
    cat("Mean Squared Error (MSE):\n")
    print(x$mse)
    cat("\n")
  }

  if (!is.null(x$classify)) {
    cat("Classification Results:\n")

    # Iterate over each model in classify results
    for (model_name in names(x$classify)) {
      model_data <- x$classify[[model_name]]

      cat(sprintf("\nModel: %s\n", model_name))

      cat("\nAverages across simulations by Category:\n")
      print(model_data$classify_res)

      cat("\nProportion of Correct Classification (CC) across categories:\n")
      cc_summary <- model_data$cc_summary
      cat(sprintf("Mean: %.3f, 95%% CI: [%.3f, %.3f]\n",
                  cc_summary$mean, cc_summary$lower, cc_summary$upper))
    }
    cat("\n")
  }

  invisible(x)
}
