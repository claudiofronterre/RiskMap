##' @title Estimation of Generalized Linear Gaussian Process Models
##' @description Fits generalized linear Gaussian process models to spatial data,
##' incorporating spatial Gaussian processes with a Matern correlation function.
##' Supports Gaussian, binomial, and Poisson response families.
##' @param formula A formula object specifying the model to be fitted.
##' The formula should include fixed effects and spatial effects specified using
##' [gp()] and optionally, random effects specified using [re()]. It can also
##' include one `offset()` term that names a column of `data`, e.g. `offset(log_exposure)`.
##' @param data An sf object containing the variables in the model.
##' @param family The distribution of the response given the random effects,
##' with its link. Either one of the strings `"gaussian"`, `"binomial"` and
##' `"poisson"`, which take the identity, logit and log links, or a family
##' object as in [stats::glm()], such as `binomial("probit")`. Binomial
##' responses take the logit, probit, complementary log-log (`"cloglog"`) and
##' Cauchy (`"cauchit"`) links, Poisson responses the log link and Gaussian
##' responses the identity link. `binomial("log")` and `poisson("identity")`
##' are not accepted, because they can give means outside the range of the
##' family, nor is `poisson("sqrt")`, whose inverse gives the same mean to a
##' positive and a negative linear predictor. Quasi-likelihood families, such
##' as [stats::quasibinomial()], have no likelihood to maximise and are not
##' accepted either. Any other link can be given in `invlink`. A named link
##' gives the same estimates as its inverse given in `invlink`, and is faster
##' when that inverse cannot be differentiated symbolically, as for
##' `invlink = pnorm`.
##' @param invlink A function that maps the linear predictor to the mean of
##' the response given the random effects, for a binomial or Poisson model
##' whose link is not one of those that `family` names, such as a link for an
##' imperfect diagnostic test. Its first and second derivatives are found
##' symbolically where possible and numerically otherwise. They only steer the
##' importance sampler, so they affect the speed of a fit and not its
##' estimates. Not applicable for Gaussian models, and not to be combined with
##' a family object that names a link other than logit or log.
##' @param denominator Optional denominator or exposure for binomial or Poisson distributions.
##' Passed as a bare/unquoted column name present in `data`.
##' If not provided, defaults to `1` for binomial models.
##' @param model_crs Optional CRS (e.g. an EPSG code) used internally for model fitting.
##' If `data` is already in a projected CRS and `model_crs` is not provided, the data are used as-is.
##' If `data` are in longitude/latitude and `model_crs` is not provided, the data are automatically
##' reprojected to an appropriate UTM zone (see [propose_utm()]) and a message reports the conversion used.
##' If provided, `model_crs` must be a projected (not longitude/latitude) CRS.
##' @param distance_units Character string, either `"km"` or `"m"`, indicating the units used
##' for internal coordinate distances and spatial parameters such as `phi`. The coordinates in
##' the returned `data` retain the linear units declared by its CRS.
##' Defaults to `"km"`.
##' @param control_mcmc Control parameters for MCMC sampling for binomial or Poisson models.
##' Must be an object of class `RiskMap_control_mcmc` as returned by [set_control_mcmc()].
##' @param control_mcml Control parameters for adaptive Monte Carlo maximum
##' likelihood estimation. Must be returned by [set_control_mcml()].
##' @param seed Optional non-negative integer that seeds every random step of
##' a binomial or Poisson fit: the Monte Carlo samples of the estimation and
##' the final conditional samples. When `NULL` (the default), one integer is
##' drawn from the session's random number generator, so [set.seed()] before
##' the call also reproduces the fit. The seed used is returned as `seed`.
##' The caller's random number state is otherwise left unchanged. Ignored for
##' Gaussian models, which have no random step.
##' @param return_samples Logical indicating whether to draw and retain
##' conditional samples at the fitted parameters. These samples are used to
##' compute fitted values. If `FALSE`, final sampling is skipped and fitted
##' values are not computed. Defaults to `TRUE`. Ignored for Gaussian models.
##' @param messages Logical indicating whether to print messages, and a
##'   line for each stage of Monte Carlo maximum likelihood as described in
##'   [riskmap_progress]. Defaults to the option `RiskMap.messages`, or
##'   `TRUE` if it is not set.
##' @param fix_var_me Optional fixed value for the measurement error variance when fitting a Gaussian model.
##' When not provided, the value will be estimated, but cannot be if each location only has one sample and
##' the `nugget` term in [gp()] is also set to `TRUE`.
##' @param start_pars Optional list of starting values for model parameters:
##' \describe{
##'   \item{beta}{regression coefficients}
##'   \item{sigma2}{spatial process variance}
##'   \item{tau2}{nugget effect variance}
##'   \item{phi}{spatial correlation scale}
##'   \item{sigma2_me}{measurement error variance - only for Gaussian models}
##'   \item{sigma2_re}{random effects variances - only when random effects are included}
##' }
##' @details
##' Generalized linear Gaussian process models extend generalized linear models (GLMs) by incorporating
##' spatial Gaussian processes to account for spatial correlation in the data. This function fits GLGPMs
##' using maximum likelihood methods, allowing for Gaussian, binomial, and Poisson response families.
##' In the case of the Binomial and Poisson families, a Monte Carlo maximum likelihood algorithm is used.
##'
##' The spatial Gaussian process is modeled with a Matern correlation function, which is flexible and
##' commonly used in geostatistical modeling. The function supports both spatial covariates and
##' unstructured random effects, providing a comprehensive framework to analyze spatially correlated
##' data across different response distributions.
##'
##' Additionally, the function allows for the inclusion of unstructured random effects, specified through
##' the [re()] term in the model formula. These random effects can capture unexplained variability
##' at specific locations beyond the fixed and spatial covariate effects, enhancing the model's flexibility
##' in capturing complex spatial patterns.
##'
##' The `model_crs` argument can be used to reproject the spatial coordinates to a different, projected
##' CRS for model fitting. If `data` are in longitude/latitude and `model_crs` is not supplied, the
##' coordinates are automatically reprojected to an appropriate UTM zone and a message details the
##' conversion; automatic projection is avoided for `model_crs` explicitly supplied in longitude/latitude,
##' which raises an error.
##' The `distance_units` argument controls whether distances between locations (and so parameters such
##' as `phi`) are expressed in kilometers (`"km"`, the default) or meters (`"m"`).
##'
##' The `control_mcmc` argument specifies the control parameters for MCMC sampling.
##' This argument must be an object returned by [set_control_mcmc()].
##' Non-Gaussian models are fitted by certified Monte Carlo maximum likelihood
##' (MCML), controlled by `control_mcml` (see [set_control_mcml()]). The
##' estimate is corrected with an independent validation sample and certified
##' with a second, fresh one, which bounds its Monte Carlo error in statistical
##' standard errors, with an allowance for that sample's own error taken from
##' it and from a replicate of it (see [set_control_mcml()]). The samples
##' are enlarged only when the certificate exceeds
##' the tolerance. Every stage is recorded in `mcml_history`, and a warning is
##' issued if the tolerance is not met within the configured limits. Warnings
##' are kept for what changes how a result is read, and [check_fit()] goes through
##' all the checks on the computation, with what to do about any that need a look.
##' The estimation samples and the final conditional samples are drawn from
##' separate streams derived from `seed`, so `glgpm(..., seed = fit$seed)`
##' reproduces a fit exactly.
##'
##' The `start_pars` argument allows starting values to be supplied explicitly.
##' Otherwise, binomial and Poisson models compare two automatic candidates
##' and use the one with the better Laplace approximation.
##'
##' @return An object of class `RiskMap` containing the fitted model and relevant information:
##'
##' \item{estimate}{Estimated parameters as a named list: \code{beta} (named
##' regression coefficients on the scale of the supplied covariates), \code{sigma2},
##' \code{phi}, \code{nu2} \eqn{= \tau^2/\sigma^2} when the nugget is
##' estimated, \code{sigma2_me} for Gaussian models with an estimated
##' measurement error variance, and \code{sigma2_re} (a named vector, one
##' entry per random effect) when random effects are included. Keeping these
##' in separate list elements - rather than one flat named vector - avoids a
##' covariate whose name collides with a parameter name (e.g. a covariate
##' literally called \code{"sigma2"}) corrupting name-based lookups.
##' \code{covariance} is estimated as one joint matrix over every parameter;
##' its rows/columns are named to match \code{unlist(estimate)}.}
##' \item{grad_MLE}{Gradient of the maximum likelihood function}
##' \item{covariance}{Covariance}
##' \item{log_lik}{Log likelihood}
##' \item{y}{Response variable}
##' \item{D}{Covariate matrix}
##' \item{coords}{Unique spatial coordinates}
##' \item{ID_coords}{Index of coordinates}
##' \item{re}{Random effects if `re()` is included in `formula`}
##' \item{ID_re}{Index of random effects if `re()` is included in `formula`}
##' \item{fix_tau2}{Fixed nugget effect variance}
##' \item{fix_var_me}{Fixed measurement error variance}
##' \item{formula}{Model formula}
##' \item{family}{Response family}
##' \item{link_function}{For binomial and Poisson models, the inverse link
##' (`inv`), its first and second derivatives (`d1` and `d2`) and the name of
##' the link (`link`), such as `"probit"`, which is `NULL` for an inverse given
##' in `invlink`. `NULL` for Gaussian models.}
##' \item{distance_units}{Units (`"km"` or `"m"`) used for internal coordinate distances and spatial parameters}
##' \item{data}{The `sf` data used for model fitting, in the CRS used internally}
##' \item{input_crs}{The CRS of the input data}
##' \item{kappa}{Spatial correlation parameter}
##' \item{units_m}{Distribution offset if `family` is `binomial` or `poisson`}
##' \item{cov_offset}{Covariate offset}
##' \item{call}{Matched call}
##' \item{S_samples}{For non-Gaussian models, fresh conditional MCMC samples
##' evaluated at the fitted parameters when `return_samples` is `TRUE`. These
##' are not the reference samples used for MCML estimation. The first block is
##' the location-level Gaussian effect (including the nugget when present),
##' followed by any unstructured random-effect blocks.}
##' \item{linear_predictors}{Conditional fitted linear predictors. For
##' non-Gaussian models these are Monte Carlo means over the final-parameter
##' latent samples and are `NULL` when `return_samples = FALSE`; for Gaussian
##' models they use the conditional latent mean.}
##' \item{fitted_values}{Conditional fitted response means. Binomial models
##' return probabilities; Poisson models return expected counts and therefore
##' include any supplied exposure in `denominator`. These are `NULL` for
##' non-Gaussian models when `return_samples = FALSE`.}
##' \item{mcml_history}{For non-Gaussian models, one entry per stage of
##' certified Monte Carlo maximum likelihood (see [set_control_mcml()]),
##' containing the sizes of the fitting and validation samples, the
##' importance distribution used, the optimum of the fitting sample
##' (`fitted_estimate`) and the corrected estimate (`estimate`), the optimizer
##' result, the size of the correction, the certificate and its two parts
##' (`certificate_distance` and `certificate_radius`), effective sample sizes
##' and Pareto tail diagnostics of each sample, the outcome of the stage and
##' its timing. A stage where the global importance distribution failed also
##' holds why, in `proposal_failure` or, at the start of the first stage,
##' `start_failure`: the `reason`, the number of `refits` and the last `change`
##' of its factors. Estimates use the same parameterisation as `estimate`, with
##' regression coefficients and spatial scale restored to the supplied
##' covariate and distance units. `NULL` for Gaussian models.}
##' \item{mcml_certified}{Whether the certificate of the returned estimate is
##' within the tolerance. The certificate is an approximate 95% bound, in
##' statistical standard errors, on the Monte Carlo error of the estimate, that
##' is its distance from the maximum of the exact likelihood that the fit has
##' located. It says nothing about whether that maximum is the highest.
##' `NULL` for Gaussian models.}
##' \item{mcml_selected_iteration}{The stage whose estimate is returned.
##' `NULL` for Gaussian models.}
##' \item{seed}{The seed used for the random steps of a binomial or Poisson
##' fit, whether supplied or drawn. `NULL` for Gaussian models.}
##'
##' @examples
##'
##' data(italy_sim)
##'
##' # intercept only
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' summary(fit)
##'
##' # add a covariate
##' fit <- glgpm(
##'   formula = y ~ gp() + pop_dens,
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' summary(fit)
##'
##' # add a random effect
##' fit <- glgpm(
##'   formula = y ~ gp() + re(province),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' summary(fit)
##'
##' # estimate the nugget
##' fit <- glgpm(
##'   formula = y ~ gp(nugget = TRUE),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' summary(fit)
##'
##' # binomial model
##' data(loaloa)
##'
##' # reduce iterations to speed up the fit, and set a
##' # seed for reproducibility
##' control_mcmc <- set_control_mcmc(n_sim = 100,
##'                                  burnin = 20)
##'
##' fit <- glgpm(NO_INF ~ gp(),
##'              data = loaloa[1:50,],
##'              denominator = NO_EXAM,
##'              family = "binomial",
##'              control_mcmc = control_mcmc,
##'              seed = 42,
##'              messages = FALSE)
##'
##' summary(fit)
##'
##' # family = binomial("probit") fits the same model with a probit link
##'
##' @seealso [set_control_mcmc()], [set_control_mcml()],
##'   \code{\link{summary.RiskMap}}, \code{\link{to_table}}
##' @export
glgpm <- function(formula,
                 data,
                 family,
                 invlink = NULL,
                 denominator = NULL,
                 model_crs = NULL,
                 distance_units = c("km", "m"),
                 control_mcmc = set_control_mcmc(),
                 control_mcml = set_control_mcml(),
                 seed = NULL,
                 return_samples = TRUE,
                 messages = getOption("RiskMap.messages", TRUE),
                 fix_var_me = NULL,
                 start_pars = list(beta = NULL,
                                   sigma2 = NULL,
                                   tau2 = NULL,
                                   phi = NULL,
                                   sigma2_me = NULL,
                                   sigma2_re = NULL)) {
  check_supplied(formula)
  check_supplied(data)
  check_supplied(family)

  check_data(data)
  check_formula(formula, data)
  family_link <- read_family(family, invlink)
  family <- family_link$family
  invlink <- family_link$invlink
  family_parameters <- family_link$parameters
  model <- family_traits(family)

  distance_units <- match_choice(distance_units)
  check_logical(return_samples)
  check_logical(messages)
  check_control(control_mcml, "set_control_mcml")
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)

  if (model$exact) {
    if (!is.null(invlink)) {
      abort_input("{.arg invlink} can't be given for a Gaussian model, which has the identity link.", arg = "invlink")
    }
    # `denominator` names a column, so it is read without being evaluated,
    # but a variable holding NULL passed on by another function is accepted.
    given <- !is.null(substitute(denominator)) &&
      !is.null(tryCatch(denominator, error = function(error) TRUE))
    if (given) {
      abort_input("{.arg denominator} can't be given for a Gaussian model.", arg = "denominator")
    }
    if (!is.null(fix_var_me)) check_range(fix_var_me, min = 0)
  } else if (!is.null(fix_var_me)) {
    abort_input("{.arg fix_var_me} can be given only for a Gaussian model, not a {family} one.", arg = "fix_var_me")
  }

  inter_f <- interpret_formula(formula)
  kappa <- inter_f$gp_spec$kappa
  mf <- model.frame(inter_f$pf, data = data, na.action = na.fail)

  # Extract outcome data
  y <- as.numeric(model.response(mf))
  n <- length(y)

  # Extract covariates matrix
  D <- as.matrix(model.matrix(attr(mf, "terms"), data = data))
  original_D <- D

  if(is.null(inter_f$offset)) {
    cov_offset <- rep(0, nrow(data))
  } else {
    cov_offset <- data[[inter_f$offset]]
  }

  # Define denominators for Binomial and Poisson distributions
  if (!model$exact) {
    sub_den <- substitute(denominator)
    if (is.null(sub_den)){
      units_m <- rep(1, nrow(data))
      if (model$warn_missing_denominator) {
        riskmap_warn("{.arg denominator} isn't given, so it is taken to be 1 for every observation.")
      }
    } else {
      if (!is.symbol(sub_den)) {
        abort_input("{.arg denominator} must be the unquoted name of a column of {.arg data}.", arg = "denominator")
      }
      do_name <- deparse(sub_den)
      if (!do_name %in% names(data)) {
        abort_input(c("{.arg denominator} must name a column of {.arg data}.", "x" = "There is no column {.val {do_name}}."), arg = "denominator")
      }
      units_m <- data[[do_name]]
    }
    if (!is.numeric(units_m)) {
      abort_input("{.arg denominator} must name a numeric column, not {.obj_type_friendly {units_m}}.", arg = "denominator")
    }
    if (is.integer(units_m)) units_m <- as.numeric(units_m)
    model$validate(y, units_m)
    check_control(control_mcmc, "set_control_mcmc")
  }

  if (length(inter_f$re_spec) > 0L) {
    hr_re <- inter_f$re_spec$term
  } else {
    hr_re <- NULL
  }

  random_effects <- prepare_random_effects(data, hr_re)
  n_re <- random_effects$n_re
  names_re <- random_effects$names_re
  ID_re <- random_effects$ID_re
  if (!is.null(ID_re)) {
    ID_re <- as.data.frame(ID_re)
  }
  re_unique <- random_effects$re_unique
  re_unique_f <- random_effects$re_unique_f


  # transform crs and extract coordinates
  input_crs <- st_crs(data)
  if(!is.null(model_crs)) {
    check_crs(model_crs)
    data <- st_transform(data, crs = model_crs)
    if(st_is_longlat(data)) {
      abort_input("{.arg model_crs} must be a projected coordinate reference system, not longitude and latitude.", arg = "model_crs")
    }
  } else if(st_is_longlat(data)) {
    auto_crs <- propose_utm(data)
    data <- st_transform(data, crs = auto_crs)
    # Shown even with `messages = FALSE`, because it changes the model.
    riskmap_inform(c("{.arg data} are in longitude and latitude, so they are projected to EPSG:{auto_crs} for fitting.", "i" = "Set {.arg model_crs} to choose the projection."))
  }

  coords_o <- coordinates_in_units(
    data, distance_units, if (is.null(model_crs)) "data" else "model_crs"
  )
  coordinate_index <- index_coordinate_rows(coords_o)
  coords <- coordinate_index$coordinates
  ID_coords <- coordinate_index$index

  m <- nrow(coords_o)
  s_unique <- unique(ID_coords)

  fix_tau2 <- inter_f$gp_spec$nugget

  if(all(table(ID_coords) == 1) &&
     model$measurement_error &&
     isTRUE(fix_tau2) &&
     is.null(fix_var_me)){
    abort_input(c("The nugget and the measurement error variance can't both be estimated with one observation per location.", "i" = "Set {.code nugget = FALSE} or a fixed nugget in {.fn gp}, or give {.arg fix_var_me}."), arg = "formula")
  }

  # For the Gaussian family, an observation-level random effect is perfectly
  # confounded with an estimated measurement-error variance. When there is
  # only one observation per location, it is also confounded with an estimated
  # nugget. Fixed variance components do not create this identifiability issue.
  if (model$measurement_error && n_re > 0) {
    n_levels_re <- vapply(re_unique, length, integer(1))
    is_saturated_re <- n_levels_re == n
    confounded_with_me <- is.null(fix_var_me)
    confounded_with_nugget <- isTRUE(fix_tau2) && all(table(ID_coords) == 1)

    if (any(is_saturated_re) && (confounded_with_me || confounded_with_nugget)) {
      confounded_component <- if (confounded_with_me && confounded_with_nugget) {
        "the measurement error and nugget"
      } else if (confounded_with_me) {
        "the measurement error"
      } else {
        "the nugget"
      }
      saturated <- names_re[is_saturated_re]
      n_saturated <- length(saturated)
      abort_input(c("The random effect{?s} {.val {saturated}} {?has/have} one level per observation.", "x" = "In a Gaussian model {cli::qty(n_saturated)}{?such an effect/these effects} can't be told apart from {confounded_component}.", "i" = "Drop {cli::qty(n_saturated)}{?it/them} from {.arg formula}, or fix the variance of {confounded_component}."), arg = "formula")
    }
  }

  if (messages) {
    crs <- as.list(st_crs(data))$input
    riskmap_inform("The model is fitted in {crs}, with distances in {distance_units}.")
  }

  valid_start_pars <- c("beta", "sigma2", "phi", "tau2", "sigma2_re", "sigma2_me")
  invalid <- names(start_pars)[!names(start_pars) %in% valid_start_pars]
  if (length(invalid)) {
    abort_input(c("{.arg start_pars} must name parameters among {.val {valid_start_pars}}.", "x" = "{.val {invalid}} {?is/are} not among them."), arg = "start_pars")
  }
  user_supplied_start <- any(vapply(
    start_pars,
    function(value) !is.null(value),
    logical(1)
  ))

  if(is.null(start_pars[["beta"]])) {
    # From the records the engine fits, the negative binomial family's as
    # binomial ones (engine_records()).
    records <- engine_records(family, y, if (!model$exact) units_m,
                              cov_offset, family_parameters, invlink)
    start_pars[["beta"]] <- model$start_coefficients(
      y = records$y, units_m = records$units_m, D = D,
      cov_offset = records$offset, link = if (is.list(invlink)) invlink$link,
      formula = inter_f$pf, data = data
    )
  } else {
    if (!is.numeric(start_pars[["beta"]])) {
      abort_input("{.code start_pars$beta} must be numeric, not {.obj_type_friendly {start_pars[['beta']]}}.", arg = "start_pars")
    }
    if (length(start_pars[["beta"]]) != ncol(D)) {
      abort_input("{.code start_pars$beta} must have {ncol(D)} value{?s}, one for each column of the design matrix including the intercept, not {length(start_pars[['beta']])}.", arg = "start_pars")
    }
  }

  if(is.null(start_pars[["sigma2"]])) {
    start_pars[["sigma2"]] <- 1
  } else {
    check_positive_number(start_pars[["sigma2"]])
  }

  if(is.null(start_pars[["phi"]])) {
    start_pars[["phi"]] <- quantile(pairwise_distances(coords), 0.1)
  } else {
    check_positive_number(start_pars[["phi"]])
  }

  if(isTRUE(fix_tau2)) {
    if(is.null(start_pars[["tau2"]])) {
      start_pars[["tau2"]] <- 1
    } else {
      check_positive_number(start_pars[["tau2"]])
    }
  } else {
    if(!is.null(start_pars[["tau2"]]))
      abort_input("A starting value for {.code tau2} can't be given when {.fn gp} has {.code nugget = FALSE}.", arg = "start_pars")
  }

  if(n_re > 0) {
    if(is.null(start_pars[["sigma2_re"]])) {
      start_pars[["sigma2_re"]] <- rep(1,n_re)
    } else {
      if (length(start_pars[["sigma2_re"]]) != n_re) {
        abort_input("{.code start_pars$sigma2_re} must have {n_re} value{?s}, one for each random effect, not {length(start_pars[['sigma2_re']])}.", arg = "start_pars")
      }
      if (!is.numeric(start_pars[["sigma2_re"]]) ||
          !all(is.finite(start_pars[["sigma2_re"]]) &
               start_pars[["sigma2_re"]] > 0)) {
        abort_input("The starting values in {.code start_pars$sigma2_re} must be positive finite numbers.", arg = "start_pars")
      }
    }
  } else {
    if (!is.null(start_pars[["sigma2_re"]])) {
      abort_input("Starting values for {.code sigma2_re} can't be given for a model without random effects.", arg = "start_pars")
    }
  }

  starting_diagnostics <- list(
    selected = if (user_supplied_start) "user" else "current",
    current_laplace = NA_real_,
    transformed_laplace = NA_real_,
    fallback_reason = NULL
  )
  if (!model$exact && !user_supplied_start) {
    selection <- select_nongaussian_start(
      current = start_pars,
      y = y,
      units_m = units_m,
      D = D,
      coords = coords,
      ID_coords = ID_coords,
      family = family,
      kappa = kappa,
      cov_offset = cov_offset,
      invlink = invlink,
      fix_tau2 = fix_tau2,
      n_re = n_re,
      ID_re = ID_re,
      s_unique = s_unique,
      re_unique = re_unique,
      messages = messages,
      family_parameters = family_parameters
    )
    start_pars <- selection$start
    starting_diagnostics <- selection$diagnostics
  }

  res <- fit_glgpm_design(
    y = y, D = D, coords = coords, ID_coords = ID_coords, ID_re = ID_re,
    s_unique = s_unique, re_unique = re_unique,
    units_m = if (!model$exact) units_m, cov_offset = cov_offset,
    kappa = kappa, fix_tau2 = fix_tau2,
    fix_var_me = fix_var_me, family = family, invlink = invlink,
    start_pars = start_pars, control_mcmc = control_mcmc,
    control_mcml = control_mcml, seed = seed,
    return_samples = return_samples, messages = messages,
    family_parameters = family_parameters
  )

  res$y <- y
  res$D <- original_D
  res$coords <- coords
  if (!model$exact) {
    attr(res, "starting_values") <- starting_diagnostics
  }
  res$ID_coords <- ID_coords
  if(n_re > 0) {
    res$re <- re_unique_f
    res$ID_re <- as.data.frame(ID_re)
    colnames(res$ID_re) <- names_re
  } else {
    res["re"] <- list(NULL)
    res["ID_re"] <- list(NULL)
  }
  res$fix_tau2 <- fix_tau2
  res["fix_var_me"] <- list(fix_var_me)
  res$formula <- formula
  res$family <- family
  # The family's own parameters, held at the values given with it.
  if (!is.null(family_parameters)) res$family_parameters <- family_parameters
  res$distance_units <- distance_units
  # Retain fitting controls as metadata so a fitted object can serve as a
  # complete specification when refitted to simulated or held-out data.
  attr(res, "control_mcml") <- control_mcml
  res$data <- data
  res$input_crs <- input_crs
  res$kappa <- kappa
  if (!model$exact) res$units_m <- units_m
  res$cov_offset <- cov_offset
  res$call <- match.call()
  if (model$exact) {
    gaussian_fitted <- gaussian_conditional_fitted(res)
    res$linear_predictors <- gaussian_fitted
    res$fitted_values <- gaussian_fitted
  }
  return(res)
}

##' Fit a geostatistical model to its design
##'
##' The estimation step of [glgpm()], from validated inputs and complete
##' starting values: the response `y`, the design matrix `D` in the supplied
##' covariate units, the unique coordinates `coords` in the model's distance
##' units, the location and random-effect indices, the denominators and
##' offsets, and the correlation and nugget settings. Coordinates and
##' covariates are rescaled for fitting and the estimate is returned in the
##' supplied units, with the seed used (`NULL` for Gaussian models) and the
##' scalings as attributes. It makes no call to sf, so it can run in forked
##' processes, where GDAL and PROJ can deadlock. `family_parameters` are the
##' family's own (`read_family()`), NULL for a family without them. `call` is
##' the call named in input errors, [glgpm()] when it calls this function.
##' @noRd
fit_glgpm_design <- function(y, D, coords, ID_coords, ID_re, s_unique,
                             re_unique, units_m, cov_offset, kappa, fix_tau2,
                             fix_var_me, family, invlink, start_pars,
                             control_mcmc, control_mcml, seed,
                             return_samples, messages,
                             family_parameters = NULL,
                             call = rlang::caller_env()) {
  exact <- family_traits(family)$exact
  n_re <- if (is.null(ID_re)) 0L else ncol(ID_re)
  spatial_scaling <- scale_spatial_coordinates(coords, start_pars[["phi"]])
  fitting_coords <- spatial_scaling$coordinates
  fitting_start_phi <- spatial_scaling$phi

  design_scaling <- standardize_design_matrix(D)
  D <- design_scaling$design
  start_pars[["beta"]] <- as.numeric(
    solve(design_scaling$coefficient_transform, start_pars[["beta"]])
  )

  if (exact) {
    if(is.null(fix_var_me)) {
      if(is.null(start_pars[["sigma2_me"]])) {
        start_pars[["sigma2_me"]] <- 1
      } else {
        check_positive_number(start_pars[["sigma2_me"]], call = call)
      }
    }
    gaussian_layout <- cov_par_layout(
      ncol(D), fix_tau2 = fix_tau2,
      family_parameters = gaussian_family_block(fix_var_me), n_re = n_re
    )
    res <- glgpm_lm(y = y-cov_offset, D, fitting_coords,
            kappa = kappa,
            ID_coords, ID_re, s_unique, re_unique,
            fix_var_me, fix_tau2,
            start_beta = start_pars[["beta"]],
            start_cov_pars = pack_cov_pars(gaussian_layout,
                                           sigma2    = start_pars[["sigma2"]],
                                           phi       = fitting_start_phi,
                                           nu2       = start_pars[["tau2"]],
                                           family = start_pars["sigma2_me"],
                                           sigma2_re = start_pars[["sigma2_re"]]),
            messages = messages)
    res["mcml_history"] <- list(NULL)
    res["mcml_certified"] <- list(NULL)
    res["mcml_selected_iteration"] <- list(NULL)
  } else {
    # Without a seed, draw one from the session's stream, which advances it
    # like any random call. Every random step of the fit derives from it.
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
    seed <- as.integer(seed)
    fitting_reference <- start_pars
    fitting_reference$phi <- fitting_start_phi
    res <- glgpm_nong_adaptive(
      y = y, D = D, coords = fitting_coords, units_m = units_m,
      kappa = kappa, reference_pars = fitting_reference,
      cov_offset = cov_offset, ID_coords = ID_coords, ID_re = ID_re,
      fix_tau2 = fix_tau2, family = family, invlink = invlink,
      control_mcml = control_mcml, seed = seed,
      messages = messages, distance_scale = spatial_scaling$distance_scale,
      family_parameters = family_parameters
    )
    optimizer <- attr(res, "optimizer")
    mcml_history <- optimizer$mcml_history
    res$mcml_history <- mcml_history
    res$mcml_certified <- optimizer$mcml_certified
    res$mcml_selected_iteration <- length(mcml_history)

    if (return_samples) {
      # The MCML draws target importance distributions along the optimisation
      # path, not the conditional distribution at the estimate. When
      # requested, draw once more at the final estimate and use that draw for
      # fitted values as well as user-facing diagnostics.
      progress_step("Drawing conditional samples at the estimate", messages)
      final_draw <- sample_at_nongaussian_fit(
        estimate = res$estimate,
        y = y,
        units_m = units_m,
        D = D,
        coords = fitting_coords,
        cov_offset = cov_offset,
        ID_coords = ID_coords,
        ID_re = ID_re,
        fix_tau2 = fix_tau2,
        family = family,
        link_function = res$link_function,
        kappa = kappa,
        control_mcmc = control_mcmc,
        seed = final_sampling_seed(seed, length(mcml_history)),
        messages = messages,
        family_parameters = family_parameters
      )
      progress_step_done(messages)
      res$linear_predictors <- final_draw$linear_predictors
      res$fitted_values <- final_draw$fitted_values
      res$S_samples <- final_draw$samples
      attr(res, "mcmc") <- final_draw$diagnostics
    } else {
      res["S_samples"] <- list(NULL)
      res["linear_predictors"] <- list(NULL)
      res["fitted_values"] <- list(NULL)
    }
  }

  res <- restore_fixed_effect_scale(
    res,
    design_scaling$coefficient_transform
  )
  if (!is.null(res$mcml_history)) {
    res$mcml_history <- restore_mcml_history_scale(
      res$mcml_history,
      design_scaling$coefficient_transform,
      spatial_scaling$distance_scale
    )
    # The optimiser keeps the same records of the check for other maxima.
    optimizer <- attr(res, "optimizer")
    if (!is.null(optimizer$mcml_maxima_check)) {
      optimizer$mcml_maxima_check <- res$mcml_history[[1L]]$maxima_check
    }
    if (!is.null(optimizer$mcml_maxima_screen)) {
      optimizer$mcml_maxima_screen <- res$mcml_history[[1L]]$maxima_screen
    }
    attr(res, "optimizer") <- optimizer
  }
  # The fitting engines optimise log(phi) on scaled coordinates. An additive
  # shift restores the original distance units without changing its Hessian.
  res$estimate$phi <- res$estimate$phi +
    log(spatial_scaling$distance_scale)
  attr(res, "design_scaling") <- design_scaling[c("center", "scale")]
  attr(res, "distance_scale") <- spatial_scaling$distance_scale
  res["seed"] <- list(if (!exact) seed)
  res
}

##' Conditional fitted mean for a Gaussian RiskMap model
##' @noRd
gaussian_conditional_fitted <- function(object) {
  pars <- coef(object)
  n <- length(object$y)
  n_loc <- nrow(object$coords)
  spatial_design <- matrix(0, n, n_loc)
  spatial_design[cbind(seq_len(n), object$ID_coords)] <- 1

  correlation <- matern_correlation(
    pairwise_distances(object$coords), phi = pars$phi,
    kappa = object$kappa, return_sym_matrix = TRUE
  )
  tau2 <- pars$tau2 %||%
    if (is.numeric(object$fix_tau2)) object$fix_tau2 else 0
  spatial_covariance <- pars$sigma2 * correlation
  diag(spatial_covariance) <- diag(spatial_covariance) + tau2

  latent_designs <- list(spatial_design)
  latent_covariances <- list(spatial_covariance)
  if (!is.null(object$ID_re)) {
    ID_re <- as.matrix(object$ID_re)
    for (j in seq_len(ncol(ID_re))) {
      design <- matrix(0, n, max(ID_re[, j]))
      design[cbind(seq_len(n), ID_re[, j])] <- 1
      latent_designs[[j + 1L]] <- design
      latent_covariances[[j + 1L]] <- diag(
        pars$sigma2_re[j], nrow = ncol(design)
      )
    }
  }

  Z <- do.call(cbind, latent_designs)
  G <- as.matrix(Matrix::bdiag(latent_covariances))
  measurement_variance <- pars$sigma2_me %||% object$fix_var_me %||% 0
  observed_covariance <- Z %*% G %*% t(Z)
  diag(observed_covariance) <- diag(observed_covariance) +
    measurement_variance
  root <- factor_covariance(observed_covariance, "Gaussian fitted covariance")
  fixed <- as.numeric(object$D %*% pars$beta + object$cov_offset)
  latent_mean <- G %*% t(Z) %*%
    solve_from_cholesky(root, object$y - fixed)
  as.numeric(fixed + Z %*% latent_mean)
}


##' Canonical working-scale covariance-parameter layout
##'
##' Single source of truth for where each non-beta parameter sits in the flat
##' numeric vector used by the fitting engines: `sigma2`, `phi`,
##' optionally `nu2` (only when the nugget is estimated), then the family
##' block, the family's own parameters that are estimated (the measurement
##' error variance `sigma2_me` of a Gaussian model whose variance is not
##' fixed, #253 step 2), and finally one entry per unstructured random
##' effect. Fixed-effect coefficients always occupy `1:p` and are handled
##' separately.
##'
##' @param p Number of fixed-effect coefficients.
##' @param fix_tau2 Logical; `TRUE` when the nugget is estimated, so `nu2` is
##' a free parameter.
##' @param family_parameters The names of the family's own parameters that are
##' estimated, in the order of `riskmap_families`.
##' @param n_re Number of unstructured random effects (`0` for none).
##' @return A list with `index` (a named list of integer positions: `beta`,
##' `sigma2`, `phi`, and - only when applicable - `nu2`, one for each family
##' parameter, `sigma2_re`), `family`, the names of the family block, and
##' `length` (the total vector length, including `beta`).
##' @noRd
cov_par_layout <- function(p, fix_tau2, family_parameters = character(),
                           n_re = 0L) {
  index <- list(beta = seq_len(p), sigma2 = p + 1L, phi = p + 2L)
  idx   <- p + 2L

  if (isTRUE(fix_tau2)) {
    idx <- idx + 1L
    index$nu2 <- idx
  }
  for (name in family_parameters) {
    idx <- idx + 1L
    index[[name]] <- idx
  }
  if (n_re > 0) {
    index$sigma2_re <- idx + seq_len(n_re)
    idx <- idx + n_re
  }

  list(index = index, family = family_parameters, length = idx)
}

##' The family parameters a Gaussian model estimates
##'
##' `"sigma2_me"` unless `fix_var_me` holds the measurement error variance,
##' for the family block of `cov_par_layout()`.
##' @noRd
gaussian_family_block <- function(fix_var_me) {
  if (is.null(fix_var_me)) "sigma2_me" else character()
}

##' Pack working-scale starting values using a `cov_par_layout()` layout
##'
##' @param layout A layout from `cov_par_layout()`.
##' @param sigma2,phi Starting values (always required).
##' @param nu2,sigma2_re Starting values for the parameters the layout
##' includes; only consulted when `layout$index` has a matching entry.
##' @param family Starting values of the family block, a named list, of which
##' only the parameters the layout's family block includes are consulted, so
##' passing e.g. `sigma2_me` when the layout has no `sigma2_me` slot is
##' silently ignored (matching this family/specification not using it).
##' @return A numeric vector of the non-beta working-scale parameters, in
##' `layout`'s order - the same contract `start_cov_pars` has always had.
##' @noRd
pack_cov_pars <- function(layout, sigma2, phi, nu2 = NULL, family = list(),
                          sigma2_re = NULL) {
  p   <- length(layout$index$beta)
  out <- numeric(layout$length - p)

  out[layout$index$sigma2 - p] <- sigma2
  out[layout$index$phi - p]    <- phi
  if (!is.null(layout$index$nu2))       out[layout$index$nu2 - p]       <- nu2
  for (name in layout$family) out[layout$index[[name]] - p] <- family[[name]]
  if (!is.null(layout$index$sigma2_re)) out[layout$index$sigma2_re - p] <- sigma2_re

  out
}

##' Structure the raw parameter vector returned by the optimizer
##'
##' Splits the working-scale (i.e. not yet exponentiated) vector of estimates
##' returned by the fitting engines into a named list, using the same
##' `cov_par_layout()` those fitting engines build their `ind_*` positions
##' from, so this can never drift out of step with them (see #135/#117).
##'
##' A flat named vector (as used prior to #92) can't safely be keyed by
##' parameter name: a covariate literally named e.g. `"sigma2"` collides with
##' the spatial variance parameter's own name, silently corrupting name-based
##' lookups. Structuring `estimate` as a list instead - `estimate$beta` for
##' the regression coefficients, `estimate$sigma2` for the spatial variance,
##' etc. - keeps those namespaces separate. `coef.RiskMap()`/`summary.RiskMap()`
##' read from this list directly; `summary.RiskMap()` additionally needs a
##' *flat* vector to line up against `covariance` (which is estimated as one
##' joint matrix over every parameter) - `unlist()` provides that, and
##' disambiguates the same way (e.g. `"beta.sigma2"` vs `"sigma2"`).
##'
##' @noRd
structure_estimate <- function(par, beta_names, fix_tau2,
                               family_parameters = character(),
                               re_names = NULL) {
  p      <- length(beta_names)
  layout <- cov_par_layout(p, fix_tau2 = fix_tau2,
                           family_parameters = family_parameters,
                           n_re = length(re_names))

  beta <- par[layout$index$beta]
  names(beta) <- beta_names

  out <- list(beta = beta,
              sigma2 = unname(par[layout$index$sigma2]),
              phi    = unname(par[layout$index$phi]))

  if (!is.null(layout$index$nu2))       out$nu2       <- unname(par[layout$index$nu2])
  for (name in layout$family) out[[name]] <- unname(par[layout$index[[name]]])
  if (!is.null(layout$index$sigma2_re)) {
    sigma2_re <- unname(par[layout$index$sigma2_re])
    names(sigma2_re) <- re_names
    out$sigma2_re <- sigma2_re
  }

  out
}

##' Compute the nugget-to-spatial variance ratio
##'
##' The parser stores `TRUE` when the nugget is estimated and a numeric value
##' when it is fixed. Keeping this distinction in one helper prevents the
##' likelihood and its derivatives from interpreting the same model
##' differently.
##'
##' @noRd
nugget_ratio <- function(fix_tau2, sigma2, log_nu2 = NULL) {
  if (isTRUE(fix_tau2)) {
    if (is.null(log_nu2)) {
      riskmap_bug("{.arg log_nu2} is required when the nugget is estimated.")
    }

    return(exp(log_nu2))
  }

  as.numeric(fix_tau2) / sigma2
}

##' Convert log weights to probabilities without overflow or underflow
##'
##' The common shift cancels during normalisation, while keeping every
##' exponentiated value no greater than one.
##'
##' @noRd
normalise_log_weights <- function(log_weights) {
  shift <- max(log_weights)
  if (!is.finite(shift)) {
    riskmap_abort("All Monte Carlo importance weights are non-finite.", class = "riskmap_error_fit")
  }

  weights <- exp(log_weights - shift)
  weights / sum(weights)
}

##' Obtain stable Gaussian regression starting values
##'
##' QR decomposition avoids squaring the condition number as the normal
##' equations do, and lets us reject a rank-deficient model matrix explicitly.
##'
##' @noRd
linear_start_values <- function(design, response) {
  fit <- lm.fit(x = design, y = response)
  if (fit$rank < ncol(design)) {
    abort_input(c("The design matrix is rank deficient, so starting values for the regression coefficients can't be found.", "i" = "Check {.arg formula} for covariates that are linear combinations of others."), arg = "formula")
  }

  unname(fit$coefficients)
}

##' Standardise the fixed-effect design without changing its column space
##'
##' With an intercept, non-intercept columns are centred and scaled. Without
##' an intercept they are only scaled, because centring would implicitly add
##' an intercept to the fitted model. `coefficient_transform` maps coefficients
##' on the standardised scale back to the scale supplied by the user.
##'
##' @noRd
standardize_design_matrix <- function(design) {
  design <- as.matrix(design)
  p <- ncol(design)
  coefficient_names <- colnames(design)
  intercept <- which(coefficient_names == "(Intercept)")

  if (length(intercept) > 1L) {
    abort_input("The design matrix has more than one intercept column.", arg = "formula")
  }

  columns_to_scale <- setdiff(seq_len(p), intercept)
  center <- setNames(numeric(p), coefficient_names)
  scale <- setNames(rep(1, p), coefficient_names)

  if (length(columns_to_scale) > 0L) {
    if (length(intercept) == 1L) {
      center[columns_to_scale] <-
        colMeans(design[, columns_to_scale, drop = FALSE])
    }

    centered <- sweep(
      design[, columns_to_scale, drop = FALSE],
      2L,
      center[columns_to_scale],
      FUN = "-"
    )
    scale[columns_to_scale] <- sqrt(colSums(centered^2) /
                                      max(1, nrow(design) - 1L))

    invalid <- !is.finite(scale[columns_to_scale]) |
      scale[columns_to_scale] == 0
    if (any(invalid)) {
      invalid_names <- coefficient_names[columns_to_scale][invalid]
      abort_input("The design matrix has {?a /}constant or numerically constant column{?s}, {.val {invalid_names}}.", arg = "formula")
    }
  }

  coefficient_transform <- diag(p)
  dimnames(coefficient_transform) <-
    list(coefficient_names, coefficient_names)
  if (length(columns_to_scale) > 0L) {
    coefficient_transform[cbind(columns_to_scale, columns_to_scale)] <-
      1 / scale[columns_to_scale]
    if (length(intercept) == 1L) {
      coefficient_transform[intercept, columns_to_scale] <-
        -center[columns_to_scale] / scale[columns_to_scale]
    }
  }

  list(
    design = design %*% coefficient_transform,
    coefficient_transform = coefficient_transform,
    center = center,
    scale = scale
  )
}

##' Restore fixed-effect results to the scale supplied by the user
##'
##' Applies the same linear transformation to the estimates and their full
##' covariance matrix, including cross-covariances with covariance parameters.
##' The score is transformed with the inverse transpose Jacobian.
##'
##' The Jacobian is block diagonal, so each block is transformed on its own.
##' A parameter without a standard error, such as a variance estimated at
##' zero, has missing values in its row and column of the covariance matrix,
##' and a product with the full Jacobian would spread them to every entry.
##'
##' @noRd
restore_fixed_effect_scale <- function(result, coefficient_transform) {
  p <- nrow(coefficient_transform)
  parameter_names <- names(unlist(result$estimate))
  beta_index <- seq_len(p)
  other_index <- setdiff(seq_along(parameter_names), beta_index)

  beta_names <- names(result$estimate$beta)
  result$estimate$beta <- as.numeric(
    coefficient_transform %*% result$estimate$beta
  )
  names(result$estimate$beta) <- beta_names

  covariance <- result$covariance
  covariance[beta_index, beta_index] <- coefficient_transform %*%
    covariance[beta_index, beta_index, drop = FALSE] %*%
    t(coefficient_transform)
  cross <- coefficient_transform %*%
    covariance[beta_index, other_index, drop = FALSE]
  covariance[beta_index, other_index] <- cross
  covariance[other_index, beta_index] <- t(cross)
  result$covariance <- covariance
  dimnames(result$covariance) <- list(parameter_names, parameter_names)

  result$grad_MLE[beta_index] <- solve(
    t(coefficient_transform), result$grad_MLE[beta_index]
  )
  result$grad_MLE <- as.numeric(result$grad_MLE)
  names(result$grad_MLE) <- parameter_names
  result
}

##' Restore MCML history estimates to user-facing covariate and distance scales
##' @noRd
restore_mcml_history_scale <- function(history, coefficient_transform,
                                       distance_scale) {
  lapply(history, function(entry) {
    for (field in intersect(c("estimate", "fitted_estimate"), names(entry))) {
      if (is.list(entry[[field]])) {
        entry[[field]]$beta <- as.numeric(
          coefficient_transform %*% entry[[field]]$beta
        )
        names(entry[[field]]$beta) <- colnames(coefficient_transform)
        entry[[field]]$phi <- entry[[field]]$phi + log(distance_scale)
      } else {
        p <- ncol(coefficient_transform)
        entry[[field]][seq_len(p)] <- as.numeric(
          coefficient_transform %*% entry[[field]][seq_len(p)]
        )
        entry[[field]][p + 2L] <-
          entry[[field]][p + 2L] + log(distance_scale)
      }
    }
    if (!is.null(entry$maxima_check)) {
      entry$maxima_check$log_phi <- entry$maxima_check$log_phi +
        log(distance_scale)
    }
    screen <- entry$maxima_screen
    if (!is.null(screen$screen)) {
      screen$screen$log_phi <- screen$screen$log_phi + log(distance_scale)
    }
    if (!is.null(screen$alarm)) {
      screen$alarm$log_phi <- screen$alarm$log_phi + log(distance_scale)
      for (table in intersect(c("bounds", "shares"), names(screen$alarm))) {
        screen$alarm[[table]]$log_phi <- screen$alarm[[table]]$log_phi +
          log(distance_scale)
      }
    }
    if (!is.null(screen$near_limit)) {
      screen$near_limit$log_phi <- screen$near_limit$log_phi +
        log(distance_scale)
    }
    if (!is.null(screen$independent$beta)) {
      screen$independent$beta <- as.numeric(
        coefficient_transform %*% screen$independent$beta
      )
      names(screen$independent$beta) <- colnames(coefficient_transform)
    }
    if (!is.null(screen)) entry$maxima_screen <- screen
    entry
  })
}

##' Wrap an optimiser objective with controlled invalid-trial handling
##'
##' Covariance parameters can temporarily define a non-positive-definite
##' matrix. Such a trial is outside the valid parameter space, so return a
##' finite penalty rather than allowing a linear-algebra error to abort the
##' complete fit. The counter is retained for post-fit diagnostics.
##'
##' @noRd
safe_optimizer_objective <- function(objective) {
  diagnostics <- new.env(parent = emptyenv())
  diagnostics$invalid_evaluations <- 0L
  penalty <- sqrt(.Machine$double.xmax)

  wrapped <- function(par) {
    value <- tryCatch(objective(par), error = function(error) NA_real_)
    if (length(value) != 1L || !is.finite(value)) {
      diagnostics$invalid_evaluations <-
        diagnostics$invalid_evaluations + 1L
      return(penalty)
    }

    value
  }
  attr(wrapped, "diagnostics") <- diagnostics
  wrapped
}

##' Collect diagnostics needed to assess an optimiser result
##'
##' @noRd
collect_optimizer_diagnostics <- function(estimate, objective, gradient,
                                          information, information_root) {
  diagnostics <- estimate[c("convergence", "message", "evaluations")]
  diagnostics$invalid_evaluations <-
    attr(objective, "diagnostics")$invalid_evaluations
  diagnostics$max_abs_gradient <- max(abs(gradient))
  diagnostics$stationary <- diagnostics$max_abs_gradient <= 1e-3
  diagnostics$information_rcond <- rcond(information)
  diagnostics$information_jitter <- attr(information_root, "jitter")
  diagnostics
}

##' Report an optimiser result that should not be treated as converged
##'
##' @noRd
warn_unconverged_optimizer <- function(diagnostics) {
  if (diagnostics$convergence != 0 && !diagnostics$stationary) {
    largest_score <- format(diagnostics$max_abs_gradient, scientific = TRUE)
    riskmap_warn(c("Model optimisation did not converge (code {diagnostics$convergence}).", "i" = "The optimiser reported {.val {diagnostics$message}}.", "i" = "The largest absolute score is {largest_score}."), class = "riskmap_warning_fit")
  }
}

##' Evaluate log(1 + exp(x)) without overflow
##'
##' @noRd
softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

##' Laplace approximation used to compare automatic non-Gaussian starts
##'
##' Constants that do not depend on the candidate are retained for clarity.
##' The comparison is only used to choose an initial posterior reference; it is
##' not reported as the fitted model likelihood.
##'
##' @noRd
starting_laplace_log_marginal <- function(candidate, y, units_m, D, coords,
                                          ID_coords, family, kappa,
                                          cov_offset, ID_re = NULL,
                                          fix_tau2 = FALSE) {
  distances <- pairwise_distances(coords)
  covariance <- candidate$sigma2 * matern_correlation(
    distances,
    phi = candidate$phi,
    kappa = kappa,
    return_sym_matrix = TRUE
  )
  nugget_variance <- if (isTRUE(fix_tau2)) {
    candidate$tau2
  } else if (is.numeric(fix_tau2) && fix_tau2 > 0) {
    fix_tau2
  } else {
    0
  }
  if (nugget_variance > 0) {
    diag(covariance) <- diag(covariance) + nugget_variance
  }
  covariance_root <- factor_covariance(
    covariance,
    "starting-value spatial covariance"
  )
  mu <- as.numeric(D %*% candidate$beta) + cov_offset
  laplace <- maxim_integrand(
    y = y,
    denominator = units_m,
    mu = mu,
    Sigma = covariance,
    ID_coords = ID_coords,
    ID_re = ID_re,
    sigma2_re = candidate$sigma2_re,
    family = family
  )
  mode <- laplace$mode
  eta <- mu + mode[ID_coords]
  if (!is.null(ID_re)) {
    offset <- nrow(covariance)
    for (j in seq_len(ncol(ID_re))) {
      n_levels <- max(ID_re[, j])
      indices <- offset + seq_len(n_levels)
      eta <- eta + mode[indices][ID_re[, j]]
      offset <- offset + n_levels
    }
  }
  log_likelihood <- family_traits(family)$start_log_likelihood(eta, y,
                                                                units_m)
  if (is.null(log_likelihood)) return(-Inf)
  log_prior <- -0.5 * (
    nrow(covariance) * log(2 * pi) +
      log_determinant_from_cholesky(covariance_root) +
      sum(mode[seq_len(nrow(covariance))] * solve_from_cholesky(
        covariance_root, mode[seq_len(nrow(covariance))]
      ))
  )
  if (!is.null(ID_re)) {
    offset <- nrow(covariance)
    for (j in seq_len(ncol(ID_re))) {
      n_levels <- max(ID_re[, j])
      random_mode <- mode[offset + seq_len(n_levels)]
      variance <- candidate$sigma2_re[j]
      log_prior <- log_prior - 0.5 * (
        n_levels * log(2 * pi * variance) +
          sum(random_mode^2) / variance
      )
      offset <- offset + n_levels
    }
  }
  laplace_root <- factor_covariance(
    laplace$Sigma.tilde,
    "starting-value Laplace covariance"
  )
  log_likelihood + log_prior +
    0.5 * length(mode) * log(2 * pi) +
    0.5 * log_determinant_from_cholesky(laplace_root)
}

##' Select between the existing and transformed-Gaussian starting values
##'
##' The conservative eligibility conditions match the model structures used in
##' validation. All failures retain the existing start so that automatic
##' initialisation cannot make a previously valid fit fail. Both starts are
##' judged on the records the engine fits (`engine_records()`), with the
##' family's own `family_parameters`, so a negative binomial model is started
##' as the binomial model of its records.
##'
##' @noRd
select_nongaussian_start <- function(current, y, units_m, D, coords,
                                     ID_coords, family, kappa, cov_offset,
                                     invlink, fix_tau2, n_re,
                                     ID_re = NULL, s_unique = NULL,
                                     re_unique = NULL,
                                     messages = TRUE,
                                     family_parameters = NULL) {
  diagnostics <- list(
    selected = "current",
    current_laplace = NA_real_,
    transformed_laplace = NA_real_,
    fallback_reason = NULL
  )
  ineligible_reason <- if (!is.null(invlink)) {
    "custom inverse links are not supported by the transformed initializer"
  } else {
    NULL
  }
  if (!is.null(ineligible_reason)) {
    diagnostics$fallback_reason <- ineligible_reason
    return(list(start = current, diagnostics = diagnostics))
  }

  records <- engine_records(family, y, units_m, cov_offset,
                            family_parameters)
  family <- records$family
  units_m <- records$units_m
  cov_offset <- records$offset
  traits <- family_traits(family)
  transformed_fit <- tryCatch({
    transformed_response <- traits$empirical_link(y, units_m)
    # With one observation per location, transformed sampling variation and a
    # Gaussian nugget are otherwise confounded. The Gaussian engine currently
    # accepts a scalar measurement variance, so fix it at the median
    # delta-method variance when the nugget is estimated.
    transformed_measurement_variance <- if (isTRUE(fix_tau2)) {
      traits$transform_variance(y, units_m)
    } else {
      NULL
    }
    glgpm_lm(
      y = transformed_response - cov_offset,
      D = D,
      coords = coords,
      kappa = kappa,
      ID_coords = ID_coords,
      ID_re = ID_re,
      s_unique = s_unique,
      re_unique = re_unique,
      fix_var_me = transformed_measurement_variance,
      fix_tau2 = fix_tau2,
      start_beta = current$beta,
      start_cov_pars = c(
        current$sigma2,
        current$phi,
        if (isTRUE(fix_tau2)) current$tau2,
        if (is.null(transformed_measurement_variance)) 1,
        current$sigma2_re
      ),
      messages = FALSE
    )
  }, error = identity)
  if (inherits(transformed_fit, "error")) {
    diagnostics$fallback_reason <- paste(
      "transformed Gaussian fit failed:",
      condition_text(transformed_fit)
    )
    return(list(start = current, diagnostics = diagnostics))
  }

  transformed <- current
  transformed$beta <- unname(transformed_fit$estimate$beta)
  transformed$sigma2 <- exp(transformed_fit$estimate$sigma2)
  transformed$phi <- exp(transformed_fit$estimate$phi)
  if (isTRUE(fix_tau2)) {
    transformed$tau2 <- transformed$sigma2 *
      exp(transformed_fit$estimate$nu2)
  }
  if (n_re > 0L) {
    transformed$sigma2_re <- exp(transformed_fit$estimate$sigma2_re)
  }
  candidate_values <- c(
    transformed$beta,
    transformed$sigma2,
    transformed$phi,
    transformed$tau2,
    transformed$sigma2_re
  )
  if (any(!is.finite(candidate_values)) ||
      transformed$sigma2 <= 0 || transformed$phi <= 0 ||
      (!is.null(transformed$tau2) && transformed$tau2 <= 0) ||
      (!is.null(transformed$sigma2_re) &&
       any(transformed$sigma2_re <= 0))) {
    diagnostics$fallback_reason <-
      "transformed Gaussian fit returned invalid starting values"
    return(list(start = current, diagnostics = diagnostics))
  }

  scores <- tryCatch(
    c(
      current = starting_laplace_log_marginal(
        current, y, units_m, D, coords, ID_coords, family, kappa, cov_offset,
        ID_re, fix_tau2
      ),
      transformed_gaussian = starting_laplace_log_marginal(
        transformed, y, units_m, D, coords, ID_coords, family, kappa,
        cov_offset, ID_re, fix_tau2
      )
    ),
    error = identity
  )
  if (inherits(scores, "error")) {
    diagnostics$fallback_reason <- paste(
      "starting-value Laplace comparison failed:",
      condition_text(scores)
    )
    return(list(start = current, diagnostics = diagnostics))
  }
  diagnostics$current_laplace <- unname(scores["current"])
  diagnostics$transformed_laplace <-
    unname(scores["transformed_gaussian"])
  if (is.finite(scores["transformed_gaussian"]) &&
      (!is.finite(scores["current"]) ||
       scores["transformed_gaussian"] > scores["current"])) {
    diagnostics$selected <- "transformed_gaussian"
    if (messages) {
      riskmap_inform("Using transformed-Gaussian automatic starting values.")
    }
    return(list(start = transformed, diagnostics = diagnostics))
  }

  if (!is.finite(scores["current"])) {
    diagnostics$fallback_reason <-
      "neither automatic starting candidate had a finite Laplace value"
  } else if (!is.finite(scores["transformed_gaussian"])) {
    diagnostics$fallback_reason <-
      "transformed Gaussian candidate had a non-finite Laplace value"
  }
  list(start = current, diagnostics = diagnostics)
}

##' @importFrom Matrix Matrix forceSymmetric
##' @noRd
glgpm_lm <- function(y, D, coords, kappa, ID_coords, ID_re, s_unique, re_unique,
                     fix_var_me, fix_tau2, start_beta, start_cov_pars, messages) {

  latent_dimension <- nrow(coords)
  if (!is.null(ID_re)) {
    latent_dimension <- latent_dimension +
      sum(vapply(seq_len(ncol(ID_re)), function(j) max(ID_re[, j]), numeric(1)))
  }
  if (use_direct_gaussian_covariance(
    length(y), latent_dimension, fix_var_me, fix_tau2
  )) {
    return(glgpm_lm_direct(
      y, D, coords, kappa, ID_coords, ID_re,
      fix_var_me, fix_tau2, start_beta, start_cov_pars, messages
    ))
  }

  m <- length(y)
  p <- ncol(D)
  U <- pairwise_distances(coords)
  if(is.null(ID_re)) {
    n_re <- 0
  } else {
    n_re <- ncol(ID_re)
  }

  ID_g <- as.matrix(cbind(ID_coords, ID_re))

  n_dim_re <- sapply(1:(n_re+1), function(i) length(unique(ID_g[,i])))
  C_g <- matrix(0, nrow = m, ncol = sum(n_dim_re))

  for(i in 1:m) {
    ind_s_i <- which(s_unique==ID_g[i,1])
    C_g[i,1:n_dim_re[1]][ind_s_i] <- 1
  }

  if(n_re>0) {
    for(j in 1:n_re) {
      select_col <- sum(n_dim_re[1:j])

      for(i in 1:m) {
        ind_re_j_i <- which(re_unique[[j]]==ID_g[i,j+1])
        C_g[i,select_col+1:n_dim_re[j+1]][ind_re_j_i] <- 1
      }
    }
  }
  C_g <- Matrix(C_g, sparse = TRUE, doDiag = FALSE)


  C_g_m <- Matrix::t(C_g)%*%C_g
  C_g_m <- forceSymmetric(C_g_m)


  layout <- cov_par_layout(p, fix_tau2 = fix_tau2,
                           family_parameters = gaussian_family_block(fix_var_me),
                           n_re = n_re)
  ind_beta   <- layout$index$beta
  ind_sigma2 <- layout$index$sigma2
  ind_phi    <- layout$index$phi
  if (!is.null(layout$index$nu2))       ind_nu2       <- layout$index$nu2
  if (!is.null(layout$index$sigma2_me)) ind_omega2    <- layout$index$sigma2_me
  if (!is.null(layout$index$sigma2_re)) ind_sigma2_re <- layout$index$sigma2_re


  log.lik <- function(par) {
    beta <- par[ind_beta]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[ind_phi])
    if(!isTRUE(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    if(n_re>0) {
      sigma2_re <- exp(par[ind_sigma2_re])
    }

    if(!is.null(fix_var_me)) {
      omega2 <- fix_var_me
    } else {
      omega2 <- exp(par[ind_omega2])
    }


    R <- matern_correlation(U,phi = phi, kappa=kappa,return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2

    Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g_inv <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- sigma2*R
    Sigma_g_inv[1:n_dim_re[1], 1:n_dim_re[1]] <-
      solve(R)/sigma2
    if(n_re > 0) {
      for(j in 1:n_re) {
        select_col <- sum(n_dim_re[1:j])

        diag(Sigma_g[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          sigma2_re[j]

        diag(Sigma_g_inv[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          1/sigma2_re[j]

      }
    }

    Sigma_g_inv <- Matrix(Sigma_g_inv, sparse = TRUE)
    Sigma_g_inv <- forceSymmetric(Sigma_g_inv)

    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- as.numeric(Matrix::t(C_g)%*%diff.y)
    Sigma_star <- Sigma_g_inv+C_g_m/omega2
    Sigma_star_inv <- forceSymmetric(solve(Sigma_star))

    q.f.y <- as.numeric(sum(diff.y^2)/omega2)
    q.f.y_tilde <- as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde/
                                (omega2^2))
    Sigma_g_C_g_m <- Sigma_g%*%C_g_m
    Sigma_tilde <- Sigma_g_C_g_m/omega2
    Matrix::diag(Sigma_tilde) <- Matrix::diag(Sigma_tilde) + 1
    log_det <- as.numeric(m*log(omega2)+Matrix::determinant(Sigma_tilde)$modulus)

    out <- -0.5*(log_det+q.f.y-q.f.y_tilde)
    return(out)
  }

  D.tilde <- t(D)%*%C_g
  U <- pairwise_distances(coords)

  grad.log.lik <- function(par) {
    beta <- par[ind_beta]
    sigma2 <- exp(par[ind_sigma2])
    phi <- exp(par[ind_phi])
    if(!isTRUE(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    if(n_re>0) {
      sigma2_re <- exp(par[ind_sigma2_re])
    }
    if(!is.null(fix_var_me)) {
      omega2 <- fix_var_me
    } else {
      omega2 <- exp(par[ind_omega2])
    }

    n_p <- length(par)
    g <- rep(0, n_p)

    R <- matern_correlation(U,phi = phi, kappa=kappa,return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2

    Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g_inv <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- sigma2*R
    R.inv <- solve(R)
    Sigma_g_inv[1:n_dim_re[1], 1:n_dim_re[1]] <-
      R.inv/sigma2
    if(n_re > 0) {
      for(j in 1:n_re) {
        select_col <- sum(n_dim_re[1:j])

        diag(Sigma_g[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          sigma2_re[j]

        diag(Sigma_g_inv[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          1/sigma2_re[j]

      }
    }

    Sigma_g_inv <- Matrix(Sigma_g_inv, sparse = TRUE)
    Sigma_g_inv <- forceSymmetric(Sigma_g_inv)

    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- as.numeric(Matrix::t(C_g)%*%diff.y)
    Sigma_star <- Sigma_g_inv+C_g_m/omega2
    Sigma_star_inv <- forceSymmetric(solve(Sigma_star))
    M_aux <- D.tilde%*%Sigma_star_inv


    g[ind_beta] <- t(D)%*%diff.y/omega2-M_aux%*%diff.y.tilde/(omega2^2)

    der_Sigma_g_inv_sigma2 <- matrix(0, nrow = sum(n_dim_re),
                                     ncol = sum(n_dim_re))

    der_Sigma_g_inv_sigma2[1:n_dim_re[1], 1:n_dim_re[1]] <-
      -R.inv/sigma2^2
    der_sigma2_aux <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2
    Sigma_g_C_g_m <- Sigma_g%*%C_g_m
    Sigma_tilde <- Sigma_g_C_g_m/omega2
    Matrix::diag(Sigma_tilde) <- Matrix::diag(Sigma_tilde) + 1
    Sigma_tilde_inv <- solve(Sigma_tilde)
    der_sigma2_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    der_sigma2_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- R
    der_sigma2_Sigma_g <- der_sigma2_Sigma_g%*%C_g_m/omega2
    der_sigma2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_sigma2_Sigma_g))
    g[ind_sigma2] <- (-0.5*der_sigma2_trace-0.5*t(diff.y.tilde)%*%
                        der_sigma2_aux%*%Sigma_star_inv%*%
                        diff.y.tilde/(omega2^2))*sigma2

    der_R_phi <- matrix(0, nrow = sum(n_dim_re),
                        ncol = sum(n_dim_re))
    M.der.phi <- matern_gradient_phi(U, phi, kappa)
    der_R_phi[1:n_dim_re[1], 1:n_dim_re[1]] <-
      M.der.phi*sigma2
    der_Sigma_g_inv_phi <- Sigma_g_inv%*%der_R_phi%*%Sigma_g_inv
    der_phi_aux <- -Sigma_star_inv%*%der_Sigma_g_inv_phi
    der_phi_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    der_phi_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- sigma2*M.der.phi
    der_phi_Sigma_g <- der_phi_Sigma_g%*%C_g_m/omega2
    der_phi_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_phi_Sigma_g))
    g[p+2] <- (-0.5*der_phi_trace-0.5*t(diff.y.tilde)%*%
                 der_phi_aux%*%Sigma_star_inv%*%
                 diff.y.tilde/(omega2^2))*phi
    if(isTRUE(fix_tau2)) {
      der_R_nu2 <- matrix(0, nrow = sum(n_dim_re),
                          ncol = sum(n_dim_re))
      diag(der_R_nu2[1:n_dim_re[1], 1:n_dim_re[1]]) <-
        sigma2
      der_Sigma_g_inv_nu2_aux <- Sigma_g_inv%*%der_R_nu2%*%Sigma_g_inv
      der_nu2_aux <- -Sigma_star_inv%*%der_Sigma_g_inv_nu2_aux
      der_nu2_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
      diag(der_nu2_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]]) <- sigma2
      der_nu2_Sigma_g <- der_nu2_Sigma_g%*%C_g_m/omega2
      der_nu2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_nu2_Sigma_g))
      g[p+3] <- (-0.5*der_nu2_trace-0.5*t(diff.y.tilde)%*%
                   der_nu2_aux%*%Sigma_star_inv%*%
                   diff.y.tilde/(omega2^2))*nu2
    }


    if(is.null(fix_var_me)) {
      der_omega2_q.f.y <- -as.numeric(sum(diff.y^2)/omega2^2)
      der_omega2_Sigma_star <- -C_g_m/omega2^2
      der_omega2_Sigma_tilde <- -Sigma_g_C_g_m/omega2^2
      der_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_omega2_Sigma_tilde))
      M_beta_omega2 <- -Sigma_star_inv%*%der_omega2_Sigma_star%*%Sigma_star_inv

      num1_omega2 <- as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde)
      der_num1_omega2 <- as.numeric(t(diff.y.tilde)%*%
                                      M_beta_omega2%*%diff.y.tilde)
      der_omega2_q.f.y_tilde <- -2*num1_omega2/(omega2^3)+
        der_num1_omega2/(omega2^2)

      g[ind_omega2] <- (-0.5*(m/omega2+der_omega2_trace+
                                der_omega2_q.f.y-der_omega2_q.f.y_tilde))*omega2
    }

    if(n_re>0) {
      der_sigma2_re_trace <- list()
      sigma2_re_trace_aux <- list()
      der_sigma2_re_Sigma_tilde <- list()
      der_sigma2_re_Sigma_g <- list()
      der_Sigma_g_inv_sigma2_re_aux <- list()
      der_sigma2_re_aux <- list()
      M_beta_sigma2_re <- list()
      for(i in 1:n_re) {

        select_col <- sum(n_dim_re[1:i])
        der_Sigma_g_inv_sigma2_re_aux[[i]] <- matrix(0, nrow = sum(n_dim_re),
                                                     ncol = sum(n_dim_re))
        diag(der_Sigma_g_inv_sigma2_re_aux[[i]][select_col+1:n_dim_re[i+1],
                                                select_col+1:n_dim_re[i+1]])  <-
          -1/sigma2_re[i]^2
        der_sigma2_re_aux[[i]] <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2_re_aux[[i]]

        M_beta_sigma2_re[[i]] <- der_sigma2_re_aux[[i]]%*%Sigma_star_inv

        der_sigma2_re_Sigma_g[[i]] <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
        diag(der_sigma2_re_Sigma_g[[i]][select_col+1:n_dim_re[i+1],
                                        select_col+1:n_dim_re[i+1]]) <- 1
        der_sigma2_re_Sigma_tilde[[i]] <- der_sigma2_re_Sigma_g[[i]]%*%C_g_m/omega2
        sigma2_re_trace_aux[[i]] <- Sigma_tilde_inv%*%der_sigma2_re_Sigma_tilde[[i]]
        der_sigma2_re_trace[[i]] <- sum(Matrix::diag(sigma2_re_trace_aux[[i]]))
        g[ind_sigma2_re[i]] <- (-0.5*der_sigma2_re_trace[[i]]-0.5*t(diff.y.tilde)%*%
                                  M_beta_sigma2_re[[i]]%*%
                                  diff.y.tilde/(omega2^2))*sigma2_re[i]

      }
    }
    return(g)
  }

  DtD <- t(D)%*%D



  hessian.log.lik <- function(par) {
    beta <- par[ind_beta]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[ind_phi])
    if(!isTRUE(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    if(n_re>0) {
      sigma2_re <- exp(par[ind_sigma2_re])
    }
    if(!is.null(fix_var_me)) {
      omega2 <- fix_var_me
    } else {
      omega2 <- exp(par[ind_omega2])
    }
    n_p <- length(par)
    g <- rep(0, n_p)

    R <- matern_correlation(U,phi = phi, kappa=kappa,return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2

    Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g_inv <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- sigma2*R
    R.inv <- solve(R)
    Sigma_g_inv[1:n_dim_re[1], 1:n_dim_re[1]] <-
      R.inv/sigma2
    if(n_re > 0) {
      for(j in 1:n_re) {
        select_col <- sum(n_dim_re[1:j])

        diag(Sigma_g[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          sigma2_re[j]

        diag(Sigma_g_inv[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          1/sigma2_re[j]

      }
    }

    Sigma_g_inv <- Matrix(Sigma_g_inv, sparse = TRUE)
    Sigma_g_inv <- forceSymmetric(Sigma_g_inv)

    Sigma_g_C_g_m <- Sigma_g%*%C_g_m
    Sigma_tilde <- Sigma_g_C_g_m/omega2
    Matrix::diag(Sigma_tilde) <- Matrix::diag(Sigma_tilde) + 1
    Sigma_tilde_inv <- solve(Sigma_tilde)


    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- as.numeric(Matrix::t(C_g)%*%diff.y)
    Sigma_star <- Sigma_g_inv+C_g_m/omega2
    Sigma_star_inv <- forceSymmetric(solve(Sigma_star))
    M_aux <- D.tilde%*%Sigma_star_inv

    H <- matrix(0, n_p, n_p)

    # beta - beta
    H[ind_beta, ind_beta] <- as.matrix(-DtD/omega2+
                                         +M_aux%*%Matrix::t(D.tilde)/(omega2^2))

    # beta - sigma2
    der_Sigma_g_inv_sigma2_aux <- matrix(0, nrow = sum(n_dim_re),
                                         ncol = sum(n_dim_re))
    der_Sigma_g_inv_sigma2_aux[1:n_dim_re[1], 1:n_dim_re[1]] <-
      -R.inv/sigma2^2
    der_sigma2_aux <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2_aux
    M_beta_sigma2 <- der_sigma2_aux%*%Sigma_star_inv
    H[ind_beta, ind_sigma2] <-
      H[ind_sigma2, ind_beta] <- as.numeric(D.tilde%*%M_beta_sigma2%*%
                                              diff.y.tilde/(omega2^2))*sigma2

    # beta - phi

    # Derivatives for phi
    der_R_phi <- matrix(0, nrow = sum(n_dim_re),
                        ncol = sum(n_dim_re))
    M.der.phi <- matern_gradient_phi(U, phi, kappa)
    der_R_phi[1:n_dim_re[1], 1:n_dim_re[1]] <-
      M.der.phi*sigma2
    der_Sigma_g_inv_phi_aux <- -Sigma_g_inv%*%der_R_phi%*%Sigma_g_inv
    der_phi_aux <- Sigma_star_inv%*%der_Sigma_g_inv_phi_aux
    M_beta_phi <- der_phi_aux%*%Sigma_star_inv

    H[ind_beta, ind_phi] <-
      H[ind_phi, ind_beta] <- as.numeric(D.tilde%*%M_beta_phi%*%
                                           diff.y.tilde/(omega2^2))*phi

    # beta - nu2
    if(isTRUE(fix_tau2)) {
      # Derivatives for nu2
      der_R_nu2 <- matrix(0, nrow = sum(n_dim_re),
                          ncol = sum(n_dim_re))
      diag(der_R_nu2[1:n_dim_re[1], 1:n_dim_re[1]]) <-
        sigma2
      der_Sigma_g_inv_nu2_aux <- -Sigma_g_inv%*%der_R_nu2%*%Sigma_g_inv
      der_nu2_aux <- Sigma_star_inv%*%der_Sigma_g_inv_nu2_aux
      M_beta_nu2 <- der_nu2_aux%*%Sigma_star_inv

      H[ind_beta, ind_nu2] <-
        H[ind_nu2, ind_beta] <- as.numeric(D.tilde%*%M_beta_nu2%*%
                                             diff.y.tilde/(omega2^2))*nu2
    }


    if(is.null(fix_var_me)) {
      # beta - omega2
      der_omega2_Sigma_star <- -C_g_m/omega2^2
      M_beta_omega2 <- -Sigma_star_inv%*%
        der_omega2_Sigma_star%*%
        Sigma_star_inv
      H[ind_beta, ind_omega2] <-
        H[ind_omega2, ind_beta] <-
        -(t(D)%*%diff.y/(omega2^2)+
            -2*as.numeric(D.tilde%*%Sigma_star_inv%*%diff.y.tilde/
                            (omega2^3))+
            as.numeric(D.tilde%*%M_beta_omega2%*%diff.y.tilde/
                         (omega2^2)))*omega2
    }

    # beta - sigma2_re
    if(n_re > 0) {
      M_beta_sigma2_re <- list()
      der_sigma2_re_aux <- list()
      der_Sigma_g_inv_sigma2_re_aux <- list()
      for(i in 1:n_re) {
        select_col <- sum(n_dim_re[1:i])

        der_Sigma_g_inv_sigma2_re_aux[[i]] <- matrix(0, nrow = sum(n_dim_re),
                                                     ncol = sum(n_dim_re))
        diag(der_Sigma_g_inv_sigma2_re_aux[[i]][select_col+1:n_dim_re[i+1],
                                                select_col+1:n_dim_re[i+1]])  <-
          -1/sigma2_re[i]^2
        der_sigma2_re_aux[[i]] <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2_re_aux[[i]]

        M_beta_sigma2_re[[i]] <- der_sigma2_re_aux[[i]]%*%Sigma_star_inv
        H[ind_beta, ind_sigma2_re[i]] <-
          H[ind_sigma2_re[i], ind_beta] <-
          as.numeric(D.tilde%*%M_beta_sigma2_re[[i]]%*%
                       diff.y.tilde/(omega2^2))*sigma2_re[i]
      }
    }

    # sigma2 - sigma2
    # Derivatives for sigma2
    der2_Sigma_g_inv_sigma2_aux <- matrix(0, nrow = sum(n_dim_re),
                                          ncol = sum(n_dim_re))
    der2_Sigma_g_inv_sigma2_aux[1:n_dim_re[1], 1:n_dim_re[1]] <-
      2*R.inv/sigma2^3
    der_R_sigma2 <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    der_R_sigma2[1:n_dim_re[1], 1:n_dim_re[1]] <- R
    der_sigma2_Sigma_g <- der_R_sigma2%*%C_g_m/omega2
    sigma2_trace_aux <- Sigma_tilde_inv%*%der_sigma2_Sigma_g
    der_sigma2_trace <- sum(Matrix::diag(sigma2_trace_aux))
    der2_sigma2_trace <- sum(Matrix::diag(-sigma2_trace_aux%*%sigma2_trace_aux))
    der.sigma2 <- (-0.5*der_sigma2_trace-0.5*t(diff.y.tilde)%*%
                     M_beta_sigma2%*%
                     diff.y.tilde/(omega2^2))*sigma2
    M2_sigma2 <- Sigma_star_inv%*%(2*der_Sigma_g_inv_sigma2_aux%*%
                                     der_sigma2_aux-der2_Sigma_g_inv_sigma2_aux)%*%
      Sigma_star_inv
    H[ind_sigma2, ind_sigma2] <-as.numeric(
      der.sigma2+
        (-0.5*der2_sigma2_trace+0.5*t(diff.y.tilde)%*%
           M2_sigma2%*%
           diff.y.tilde/(omega2^2))*sigma2^2)

    # sigma2 - phi
    der_R_sigma2_phi <- der_R_phi/sigma2
    der_phi_Sigma_g <- der_R_phi%*%C_g_m/omega2
    der_sigma2_phi_Sigma_g <- der_phi_Sigma_g/sigma2
    phi_trace_aux <- Sigma_tilde_inv%*%der_phi_Sigma_g
    der_phi_trace <- sum(Matrix::diag(phi_trace_aux))
    der_sigma2_phi_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_sigma2_phi_Sigma_g-
                                               sigma2_trace_aux%*%phi_trace_aux))

    der_Sigma_g_inv_sigma2_phi_aux <-
      -Sigma_g_inv%*%(
        der_R_sigma2%*%Sigma_g_inv%*%der_R_phi+
          der_R_phi%*%Sigma_g_inv%*%der_R_sigma2-
          der_R_sigma2_phi
      )%*%
      Sigma_g_inv

    M2_sigma2_phi <- -Sigma_star_inv%*%(
      -der_Sigma_g_inv_sigma2_aux%*%der_phi_aux+
        -der_Sigma_g_inv_phi_aux%*%der_sigma2_aux-
        der_Sigma_g_inv_sigma2_phi_aux
    )%*%Sigma_star_inv

    H[ind_sigma2, ind_phi] <-
      H[ind_phi, ind_sigma2] <- as.numeric(
        (-0.5*der_sigma2_phi_trace+0.5*t(diff.y.tilde)%*%
           M2_sigma2_phi%*%
           diff.y.tilde/(omega2^2))*sigma2*phi)

    # sigma2 - nu2
    if(isTRUE(fix_tau2)) {
      der_R_sigma2_nu2 <- der_R_nu2/sigma2
      der_nu2_Sigma_g <- der_R_nu2%*%C_g_m/omega2
      der_sigma2_nu2_Sigma_g <- der_nu2_Sigma_g/sigma2
      nu2_trace_aux <- Sigma_tilde_inv%*%der_nu2_Sigma_g
      der_nu2_trace <- sum(Matrix::diag(nu2_trace_aux))
      der_sigma2_nu2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_sigma2_nu2_Sigma_g-
                                                 sigma2_trace_aux%*%nu2_trace_aux))

      der_Sigma_g_inv_sigma2_nu2_aux <-
        -Sigma_g_inv%*%(
          der_R_sigma2%*%Sigma_g_inv%*%der_R_nu2+
            der_R_nu2%*%Sigma_g_inv%*%der_R_sigma2-
            der_R_sigma2_nu2
        )%*%
        Sigma_g_inv

      M2_sigma2_nu2 <- -Sigma_star_inv%*%(
        -der_Sigma_g_inv_sigma2_aux%*%der_nu2_aux+
          -der_Sigma_g_inv_nu2_aux%*%der_sigma2_aux-
          der_Sigma_g_inv_sigma2_nu2_aux
      )%*%Sigma_star_inv

      H[ind_sigma2, ind_nu2] <-
        H[ind_nu2, ind_sigma2] <- as.numeric(
          (-0.5*der_sigma2_nu2_trace+0.5*t(diff.y.tilde)%*%
             M2_sigma2_nu2%*%
             diff.y.tilde/(omega2^2))*sigma2*nu2)
    }

    if(is.null(fix_var_me)) {
      # sigma2 - omega2
      M2_sigma2_omega2 <- -Sigma_star_inv%*%
        (der_omega2_Sigma_star%*%Sigma_star_inv%*%der_Sigma_g_inv_sigma2_aux+
           der_Sigma_g_inv_sigma2_aux%*%Sigma_star_inv%*%der_omega2_Sigma_star)%*%
        Sigma_star_inv
      der_omega2_Sigma_tilde <- -Sigma_g_C_g_m/omega2^2
      omega2_trace_aux <- Sigma_tilde_inv%*%der_omega2_Sigma_tilde
      der_sigma2_omega2_Sigma_g <- -der_sigma2_Sigma_g/omega2
      der_sigma2_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_sigma2_omega2_Sigma_g-
                                                    sigma2_trace_aux%*%omega2_trace_aux))
      der_sigma2_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%
                                                       M_beta_sigma2%*%
                                                       diff.y.tilde/
                                                       (omega2^3))+
        as.numeric(t(diff.y.tilde)%*%
                     M2_sigma2_omega2%*%diff.y.tilde/
                     (omega2^2))

      H[ind_sigma2, ind_omega2] <-
        H[ind_omega2, ind_sigma2] <-
        (-0.5*(der_sigma2_omega2_trace+
                 der_sigma2_omega2_q.f.y_tilde))*omega2*sigma2
    }


    # sigma2 - sigma2_re
    if(n_re>0) {
      sigma2_re_trace_aux <- list()
      der_sigma2_re_Sigma_g <- list()
      der_sigma2_re_Sigma_tilde <- list()
      for(i in 1:n_re) {
        select_col <- sum(n_dim_re[1:i])


        der_sigma2_re_Sigma_g[[i]] <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
        diag(der_sigma2_re_Sigma_g[[i]][select_col+1:n_dim_re[i+1],
                                        select_col+1:n_dim_re[i+1]]) <- 1
        der_sigma2_re_Sigma_tilde[[i]] <- der_sigma2_re_Sigma_g[[i]]%*%C_g_m/omega2
        sigma2_re_trace_aux[[i]] <- Sigma_tilde_inv%*%der_sigma2_re_Sigma_tilde[[i]]
        der_sigma2_sigma2_re_trace <- sum(Matrix::diag(-sigma2_trace_aux%*%sigma2_re_trace_aux[[i]]))

        M2_sigma2_sigma2_re <- -Sigma_star_inv%*%(
          -der_Sigma_g_inv_sigma2_aux%*%der_sigma2_re_aux[[i]]+
            -der_Sigma_g_inv_sigma2_re_aux[[i]]%*%der_sigma2_aux
        )%*%Sigma_star_inv

        H[ind_sigma2, ind_sigma2_re[i]] <-
          H[ind_sigma2_re[i], ind_sigma2] <- as.numeric(
            (-0.5*der_sigma2_sigma2_re_trace+0.5*t(diff.y.tilde)%*%
               M2_sigma2_sigma2_re%*%
               diff.y.tilde/(omega2^2))*sigma2*sigma2_re[i])
      }
    }

    # phi - phi
    der2_R_phi <- matrix(0, nrow = sum(n_dim_re),
                         ncol = sum(n_dim_re))
    M.der2.phi <- matern_hessian_phi(U, phi, kappa)
    der2_R_phi[1:n_dim_re[1], 1:n_dim_re[1]] <-
      M.der2.phi*sigma2

    der2_Sigma_g_inv_phi_aux <- Sigma_g_inv%*%(
      2*der_R_phi%*%Sigma_g_inv%*%der_R_phi-
        der2_R_phi)%*%Sigma_g_inv
    der2_phi_Sigma_g <- der2_R_phi%*%C_g_m/omega2
    der2_phi_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der2_phi_Sigma_g
                                       -phi_trace_aux%*%phi_trace_aux))
    der.phi <- (-0.5*der_phi_trace-0.5*t(diff.y.tilde)%*%
                  M_beta_phi%*%
                  diff.y.tilde/(omega2^2))*phi
    M2_phi <- Sigma_star_inv%*%(2*der_Sigma_g_inv_phi_aux%*%
                                  der_phi_aux-
                                  der2_Sigma_g_inv_phi_aux)%*%
      Sigma_star_inv
    H[ind_phi, ind_phi] <-as.numeric(
      der.phi+
        (-0.5*der2_phi_trace+0.5*t(diff.y.tilde)%*%
           M2_phi%*%
           diff.y.tilde/(omega2^2))*phi^2)

    # phi - nu2
    if(isTRUE(fix_tau2)) {
      der_phi_nu2_trace <- sum(Matrix::diag(-phi_trace_aux%*%nu2_trace_aux))

      der_Sigma_g_inv_phi_nu2_aux <-
        -Sigma_g_inv%*%(
          der_R_phi%*%Sigma_g_inv%*%der_R_nu2+
            der_R_nu2%*%Sigma_g_inv%*%der_R_phi
        )%*%
        Sigma_g_inv

      M2_phi_nu2 <- -Sigma_star_inv%*%(
        -der_Sigma_g_inv_phi_aux%*%der_nu2_aux+
          -der_Sigma_g_inv_nu2_aux%*%der_phi_aux-
          der_Sigma_g_inv_phi_nu2_aux
      )%*%Sigma_star_inv

      H[ind_phi, ind_nu2] <-
        H[ind_nu2, ind_phi] <- as.numeric(
          (-0.5*der_phi_nu2_trace+0.5*t(diff.y.tilde)%*%
             M2_phi_nu2%*%
             diff.y.tilde/(omega2^2))*phi*nu2)
    }

    if(is.null(fix_var_me)) {
      # phi - omega2
      M2_phi_omega2 <- -Sigma_star_inv%*%
        (der_omega2_Sigma_star%*%Sigma_star_inv%*%der_Sigma_g_inv_phi_aux+
           der_Sigma_g_inv_phi_aux%*%Sigma_star_inv%*%der_omega2_Sigma_star)%*%
        Sigma_star_inv
      der_phi_omega2_Sigma_g <- -der_phi_Sigma_g/omega2
      der_phi_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_phi_omega2_Sigma_g-
                                                 phi_trace_aux%*%omega2_trace_aux))
      der_phi_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%
                                                    M_beta_phi%*%
                                                    diff.y.tilde/
                                                    (omega2^3))+
        as.numeric(t(diff.y.tilde)%*%
                     M2_phi_omega2%*%diff.y.tilde/
                     (omega2^2))

      H[ind_phi, ind_omega2] <-
        H[ind_omega2, ind_phi] <-
        (-0.5*(der_phi_omega2_trace+
                 der_phi_omega2_q.f.y_tilde))*omega2*phi
    }


    #phi - sigma2_re
    if(n_re>0) {
      for(i in 1:n_re) {
        der_phi_sigma2_re_trace <- sum(Matrix::diag(-phi_trace_aux%*%sigma2_re_trace_aux[[i]]))

        M2_phi_sigma2_re <- -Sigma_star_inv%*%(
          -der_Sigma_g_inv_phi_aux%*%der_sigma2_re_aux[[i]]+
            -der_Sigma_g_inv_sigma2_re_aux[[i]]%*%der_phi_aux
        )%*%Sigma_star_inv

        H[ind_phi, ind_sigma2_re[i]] <-
          H[ind_sigma2_re[i], ind_phi] <- as.numeric(
            (-0.5*der_phi_sigma2_re_trace+0.5*t(diff.y.tilde)%*%
               M2_phi_sigma2_re%*%
               diff.y.tilde/(omega2^2))*phi*sigma2_re[i])
      }
    }

    if(isTRUE(fix_tau2)) {
      # nu2 - nu2
      der2_Sigma_g_inv_nu2_aux <- Sigma_g_inv%*%(
        2*der_R_nu2%*%Sigma_g_inv%*%der_R_nu2)%*%Sigma_g_inv

      der2_nu2_trace <- sum(Matrix::diag(-nu2_trace_aux%*%nu2_trace_aux))
      der.nu2 <- (-0.5*der_nu2_trace-0.5*t(diff.y.tilde)%*%
                    M_beta_nu2%*%
                    diff.y.tilde/(omega2^2))*nu2
      M2_nu2 <- Sigma_star_inv%*%(2*der_Sigma_g_inv_nu2_aux%*%
                                    der_nu2_aux-
                                    der2_Sigma_g_inv_nu2_aux)%*%
        Sigma_star_inv
      H[ind_nu2, ind_nu2] <-as.numeric(
        der.nu2+
          (-0.5*der2_nu2_trace+0.5*t(diff.y.tilde)%*%
             M2_nu2%*%
             diff.y.tilde/(omega2^2))*nu2^2)

      if(is.null(fix_var_me)) {
        # nu2 - omega2
        M2_nu2_omega2 <- -Sigma_star_inv%*%
          (der_omega2_Sigma_star%*%Sigma_star_inv%*%der_Sigma_g_inv_nu2_aux+
             der_Sigma_g_inv_nu2_aux%*%Sigma_star_inv%*%der_omega2_Sigma_star)%*%
          Sigma_star_inv
        der_nu2_omega2_Sigma_g <- -der_nu2_Sigma_g/omega2
        der_nu2_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_nu2_omega2_Sigma_g-
                                                   nu2_trace_aux%*%omega2_trace_aux))
        der_nu2_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%
                                                      M_beta_nu2%*%
                                                      diff.y.tilde/
                                                      (omega2^3))+
          as.numeric(t(diff.y.tilde)%*%
                       M2_nu2_omega2%*%diff.y.tilde/
                       (omega2^2))

        H[ind_nu2, ind_omega2] <-
          H[ind_omega2, ind_nu2] <-
          (-0.5*(der_nu2_omega2_trace+
                   der_nu2_omega2_q.f.y_tilde))*omega2*nu2
      }


      #nu2 - sigma2_re
      if(n_re>0) {
        for(i in 1:n_re) {
          der_nu2_sigma2_re_trace <- sum(Matrix::diag(-nu2_trace_aux%*%sigma2_re_trace_aux[[i]]))

          M2_nu2_sigma2_re <- -Sigma_star_inv%*%(
            -der_Sigma_g_inv_nu2_aux%*%der_sigma2_re_aux[[i]]+
              -der_Sigma_g_inv_sigma2_re_aux[[i]]%*%der_nu2_aux
          )%*%Sigma_star_inv

          H[ind_nu2, ind_sigma2_re[i]] <-
            H[ind_sigma2_re[i], ind_nu2] <- as.numeric(
              (-0.5*der_nu2_sigma2_re_trace+0.5*t(diff.y.tilde)%*%
                 M2_nu2_sigma2_re%*%
                 diff.y.tilde/(omega2^2))*nu2*sigma2_re[i])
        }
      }
    }

    if(is.null(fix_var_me)) {
      #omega2 - omega2
      der_omega2_q.f.y <- -as.numeric(sum(diff.y^2)/omega2^2)
      der2_omega2_q.f.y <- 2*as.numeric(sum(diff.y^2)/omega2^3)

      omega2_trace_aux <- Sigma_tilde_inv%*%der_omega2_Sigma_tilde
      der_omega2_trace <- sum(Matrix::diag(omega2_trace_aux))
      der2_omega2_Sigma_tilde <- 2*Sigma_g_C_g_m/omega2^3
      der2_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der2_omega2_Sigma_tilde-
                                              omega2_trace_aux%*%omega2_trace_aux))

      num1_omega2 <- as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde)
      der_num1_omega2 <- as.numeric(t(diff.y.tilde)%*%
                                      M_beta_omega2%*%diff.y.tilde)
      der2_omega2_Sigma_star <- 2*C_g_m/omega2^3
      M2_beta_omega2 <- Sigma_star_inv%*%
        (2*der_omega2_Sigma_star%*%Sigma_star_inv%*%der_omega2_Sigma_star-
           der2_omega2_Sigma_star)%*%
        Sigma_star_inv
      der2_num1_omega2 <- as.numeric(t(diff.y.tilde)%*%
                                       M2_beta_omega2%*%diff.y.tilde)

      der_omega2_q.f.y_tilde <- -2*num1_omega2/(omega2^3)+
        der_num1_omega2/(omega2^2)
      der2_omega2_q.f.y_tilde <- -2*(der_num1_omega2*omega2-3*num1_omega2)/
        (omega2^4)+
        (der2_num1_omega2*omega2-
           2*der_num1_omega2)/(omega2^3)

      #     g[ind_omega2] <- (-0.5*(m/omega2+der_omega2_trace+
      #     der_omega2_q.f.y-der_omega2_q.f.y_tilde))*omega2

      der.omega2 <- (-0.5*(m/omega2+der_omega2_trace+
                             der_omega2_q.f.y-der_omega2_q.f.y_tilde))*omega2

      H[ind_omega2, ind_omega2] <- der.omega2+
        -0.5*(-m/omega2^2+der2_omega2_trace+
                der2_omega2_q.f.y-der2_omega2_q.f.y_tilde)*(omega2^2)

      #omega2 - sigma2_re
      if(n_re>0) {
        for(i in 1:n_re) {
          der_omega2_sigma2_re_Sigma_g <- -der_sigma2_re_Sigma_tilde[[i]]/omega2
          der_omega2_sigma2_re_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_omega2_sigma2_re_Sigma_g-
                                                           omega2_trace_aux%*%sigma2_re_trace_aux[[i]]))

          M2_omega2_sigma2_re <- -Sigma_star_inv%*%
            (der_omega2_Sigma_star%*%Sigma_star_inv%*% der_Sigma_g_inv_sigma2_re_aux[[i]]+
               der_Sigma_g_inv_sigma2_re_aux[[i]]%*%Sigma_star_inv%*%der_omega2_Sigma_star)%*%
            Sigma_star_inv

          der_sigma2_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%
                                                           M_beta_sigma2_re[[i]]%*%
                                                           diff.y.tilde/
                                                           (omega2^3))+
            as.numeric(t(diff.y.tilde)%*%
                         M2_omega2_sigma2_re%*%diff.y.tilde/
                         (omega2^2))

          H[ind_omega2, ind_sigma2_re[i]] <-
            H[ind_sigma2_re[i], ind_omega2] <- as.numeric(
              -0.5*(der_omega2_sigma2_re_trace+der_sigma2_omega2_q.f.y_tilde)*
                omega2*sigma2_re[i])

        }
      }
    }
    # sigma2_re - sigma2_re
    if(n_re > 0) {
      der_sigma2_re_trace <- list()
      der.sigma2_re <- list()
      M2_sigma2_re <- list()
      der2_Sigma_g_inv_sigma2_re_aux <- list()
      der2_sigma2_re_trace <- list()
      for(i in 1:n_re) {
        select_col <- sum(n_dim_re[1:i])
        der2_Sigma_g_inv_sigma2_re_aux[[i]] <- matrix(0, nrow = sum(n_dim_re),
                                                      ncol = sum(n_dim_re))
        diag(der2_Sigma_g_inv_sigma2_re_aux[[i]][select_col+1:n_dim_re[i+1],
                                                 select_col+1:n_dim_re[i+1]]) <-
          2/sigma2_re[i]^3
        der_sigma2_re_trace[[i]] <- sum(Matrix::diag(sigma2_re_trace_aux[[i]]))
        der.sigma2_re[[i]] <- (-0.5*der_sigma2_re_trace[[i]]-0.5*t(diff.y.tilde)%*%
                                 M_beta_sigma2_re[[i]]%*%
                                 diff.y.tilde/(omega2^2))*sigma2_re[i]
        M2_sigma2_re[[i]] <- Sigma_star_inv%*%(
          2*der_Sigma_g_inv_sigma2_re_aux[[i]]%*%
            der_sigma2_re_aux[[i]]-der2_Sigma_g_inv_sigma2_re_aux[[i]])%*%
          Sigma_star_inv
        der2_sigma2_re_trace[[i]] <- sum(Matrix::diag(-sigma2_re_trace_aux[[i]]%*%
                                                        sigma2_re_trace_aux[[i]]))
        H[ind_sigma2_re[i], ind_sigma2_re[i]] <-as.numeric(
          der.sigma2_re[[i]]+
            (-0.5*der2_sigma2_re_trace[[i]]+0.5*t(diff.y.tilde)%*%
               M2_sigma2_re[[i]]%*%
               diff.y.tilde/(omega2^2))*sigma2_re[i]^2)
        if(i < n_re) {
          for(j in (i+1):n_re) {

            der_sigma2_re_ij_trace <- sum(Matrix::diag(-sigma2_re_trace_aux[[i]]%*%
                                                         sigma2_re_trace_aux[[j]]))

            M2_sigma2_re_ij <- -Sigma_star_inv%*%(
              -der_Sigma_g_inv_sigma2_re_aux[[i]]%*%der_sigma2_re_aux[[j]]+
                -der_Sigma_g_inv_sigma2_re_aux[[j]]%*%der_sigma2_re_aux[[i]]
            )%*%Sigma_star_inv

            H[ind_sigma2_re[i], ind_sigma2_re[j]] <-
              H[ind_sigma2_re[j], ind_sigma2_re[i]] <- as.numeric(
                (-0.5*der_sigma2_re_ij_trace+0.5*t(diff.y.tilde)%*%
                   M2_sigma2_re_ij%*%
                   diff.y.tilde/(omega2^2))*sigma2_re[i]*sigma2_re[j])
          }
        }
      }
    }


    return(H)
  }

  start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)]/start_cov_pars[1]
  start_par <- c(start_beta, log(start_cov_pars))

  out <- list()
  objective <- safe_optimizer_objective(function(x) -log.lik(x))
  estim <- nlminb(start_par,
                  objective,
                  function(x) -grad.log.lik(x),
                  function(x) -hessian.log.lik(x),
                  control = list(trace = 1 * messages))

  out$estimate <- structure_estimate(
    estim$par,
    beta_names = colnames(D),
    fix_tau2   = fix_tau2,
    family_parameters = gaussian_family_block(fix_var_me),
    re_names   = if (n_re > 0) names(ID_re) else NULL
  )
  out$grad_MLE <- grad.log.lik(estim$par)
  hess.MLE <- hessian.log.lik(estim$par)
  information <- -hess.MLE
  information_root <- factor_covariance(
    information,
    "observed information matrix"
  )
  out$covariance <- chol2inv(information_root)
  flat_names <- names(unlist(out$estimate))
  dimnames(out$covariance) <- list(flat_names, flat_names)
  out$log_lik <- -estim$objective
  out["link_function"] <- list(NULL)
  out["units_m"] <- list(NULL)
  out["S_samples"] <- list(NULL)
  optimizer_diagnostics <- collect_optimizer_diagnostics(
    estim, objective, out$grad_MLE, information, information_root
  )
  warn_unconverged_optimizer(optimizer_diagnostics)
  attr(out, "optimizer") <- optimizer_diagnostics

  class(out) <- "RiskMap"
  return(out)
}


##' Differentiate a vector-in/vector-out function w.r.t. `eta`
##'
##' Used to auto-derive a missing `d1`/`d2` from a user-supplied inverse link
##' (or from an already-derived `d1`, for `d2`). Tries symbolic
##' differentiation via `Deriv::Deriv()` first, since the result is fast and
##' stable when evaluated repeatedly. `Deriv::Deriv()` raises a hard error
##' when it cannot process `f`'s body (e.g. a call it has no derivative rule
##' for), so its construction is wrapped in `tryCatch()`, and a derivative
##' that fails validation is dropped. If either fails, falls back to
##' numerical differentiation by the method of `numDeriv::grad()`, applied to
##' all elements of `eta` at once (`compute_richardson_derivative()`).
##'
##' @param f Function of a single argument `eta`, returning a numeric vector.
##' @param ncheck Length of the probe vector used to validate a candidate
##' derivative before accepting it.
##' @return A function of `eta` approximating the derivative of `f`.
##' @noRd
differentiate <- function(f, ncheck) {
  is_valid <- function(d) {
    out <- tryCatch(d(rep(0, ncheck)), error = function(e) NULL)
    is.numeric(out) && length(out) == ncheck && all(is.finite(out))
  }

  symbolic <- tryCatch({
    wrapped <- function(eta) f(eta)
    Deriv::Deriv(wrapped, "eta")
  }, error = function(e) NULL)
  if (!is.null(symbolic) && !is_valid(symbolic)) symbolic <- NULL

  if (!is.null(symbolic)) return(symbolic)

  # An inverse link acts element by element, so every element can share each
  # evaluation of `f`. A function that does not is caught on a probe and
  # differentiated one value at a time, as numDeriv::grad() requires.
  shared <- function(eta) compute_richardson_derivative(f, eta)
  one_at_a_time <- function(eta) {
    vapply(eta, function(z) numDeriv::grad(f, z), numeric(1))
  }
  probe <- c(-2.5, -0.7, 0, 0.3, 1.9)
  agrees <- tryCatch(
    isTRUE(all.equal(shared(probe), one_at_a_time(probe), tolerance = 1e-10)),
    error = function(e) FALSE
  )
  if (agrees) shared else one_at_a_time
}

##' Derivative of an elementwise function by Richardson extrapolation
##'
##' The default method of `numDeriv::grad()` (central differences with
##' `rounds` halvings of the step and Richardson extrapolation, with its
##' default settings) applied to every element of `eta` at once. For a
##' function that acts element by element it gives the same values as
##' `numDeriv::grad()` called on each element, with `2 * rounds` calls of `f`
##' in all in place of `2 * rounds` per element.
##' @noRd
compute_richardson_derivative <- function(f, eta, eps = 1e-4, d = 1e-4,
                                          zero_tolerance = sqrt(
                                            .Machine$double.eps / 7e-7
                                          ),
                                          rounds = 4L, ratio = 2) {
  step <- abs(d * eta) + eps * (abs(eta) < zero_tolerance)
  estimates <- matrix(NA_real_, rounds, length(eta))
  for (round in seq_len(rounds)) {
    estimates[round, ] <- (f(eta + step) - f(eta - step)) / (2 * step)
    if (anyNA(estimates[round, ])) {
      abort_input("{.arg invlink} returns {.code NA} near the values being differentiated.", arg = "invlink")
    }
    step <- step / ratio
  }
  for (m in seq_len(rounds - 1L)) {
    estimates <- (estimates[2:(rounds + 1L - m), , drop = FALSE] * (4^m) -
                    estimates[1:(rounds - m), , drop = FALSE]) / (4^m - 1)
  }
  c(estimates)
}

##' The binomial regression that gives the starting regression coefficients
##'
##' A named link starts from a regression with its own link, since the
##' coefficients of a logit regression are on another scale (an intercept of
##' -2.94 where a cauchit regression gives -6.31). Should that regression fail
##' or not converge, as one with a cauchit link can, the logit regression is
##' used, as it is for the canonical link and for a link given as a function.
##' The regression is on the columns of the model's design matrix, so it has
##' an intercept only when the model has one, and its coefficients are in the
##' order of those columns.
##' @param y,units_m The positives and the denominators.
##' @param D The design matrix of the model.
##' @param link The name of the link, NULL for a link given as a function.
##' @param quasi Whether to regress with the quasi-binomial family, for the
##'   records of a negative binomial model, whose trials need not be whole.
##' @noRd
binomial_start_regression <- function(y, units_m, D, cov_offset,
                                      link = NULL, quasi = FALSE) {
  regress <- function(glm_family) {
    glm.fit(x = D, y = cbind(y, units_m - y), offset = cov_offset,
            family = glm_family)
  }
  if (!is.null(link) && !identical(link, "logit")) {
    own <- tryCatch(suppressWarnings(regress(binomial(link))),
                    error = function(e) NULL)
    if (isTRUE(own$converged)) return(own)
  }
  # The quasi-binomial family gives the same estimates without glm()'s
  # warning about trials that are not whole, as those of the negative
  # binomial records are (engine_records()).
  regress(if (quasi) stats::quasibinomial() else binomial())
}

##' What a refit of a fitted model passes as `invlink`
##'
##' `NULL` for a Gaussian fit or a canonical link, so that the refit takes the
##' canonical path again, and otherwise the stored inverse, its derivatives
##' and the name of the link, in the list form of `invlink`. The likelihood
##' that `profile()` rebuilds from a fit takes its link from here too.
##' @noRd
refit_invlink <- function(object) {
  link_function <- object$link_function
  if (is.null(link_function) ||
      identical(link_function$name, "canonical")) {
    return(NULL)
  }
  list(inv = link_function$inv, d1 = link_function$d1,
       d2 = link_function$d2, link = link_function$link)
}

##' Maximization of the Integrand for Generalized Linear Gaussian Process Models
##'
##' Finds the conditional mode of the latent effects in a binomial or Poisson
##' generalized linear Gaussian process model (GLGPM), that is, the maximiser of
##' the joint log-density of the responses and the latent spatial and
##' unstructured random effects, viewed as a function of the latent effects.
##'
##' @param y Response variable vector: numbers of positives for
##'   `family = "binomial"`, or counts for `family = "poisson"`.
##' @param denominator Numeric vector with one positive value per response: the
##'   number of trials for `family = "binomial"`, or the exposure (the expected
##'   count is `denominator` times the inverse link of the linear predictor) for
##'   `family = "poisson"`. This is the same quantity as the `denominator`
##'   argument of [glgpm()].
##' @param mu Vector of the fixed part of the linear predictor, one value per
##'   response.
##' @param Sigma Covariance matrix of the spatial process.
##' @param ID_coords Indices mapping response to locations.
##' @param ID_re Indices mapping response to unstructured random effects.
##' @param family The distribution of the response, `"binomial"` or
##'   `"poisson"`, or a binomial or Poisson family object with one of the
##'   links that [glgpm()] accepts, such as `binomial("probit")`.
##' @param sigma2_re Variance of the unstructured random effects.
##' @param hessian Logical; if TRUE, compute the Hessian matrix.
##' @param gradient Logical; if TRUE, compute the gradient vector.
##' @param invlink A function that maps the linear predictor to the mean of
##'   the response given the random effects, for a link that `family` does
##'   not name, as in [glgpm()].
##'
##' @details
##' The mode is found by Newton-type optimisation ([stats::nlminb()]) using the
##' analytic gradient and Hessian of the integrand. The integrand includes terms
##' for the spatial process (`Sigma`), the unstructured random effects
##' (`sigma2_re`), and the log-likelihood of the responses for the specified
##' family. With the canonical links, the terms of a response whose
##' probability rounds to 0 or 1, or whose mean underflows to 0, are computed
##' in their exact canonical form, so that the optimisation is not stopped
##' there.
##'
##' @return A list with the mode (`mode`), the name of the inverse link used
##'   (`invlink_used`) and either the Hessian of the integrand at the mode
##'   (`hessian`, if `hessian = TRUE`) or the inverse of its negative
##'   (`Sigma.tilde`, otherwise). If `gradient = TRUE`, the gradient at the mode
##'   is also returned (`gradient`).
##' @examples
##' coords <- cbind(c(0, 1, 2), c(0, 0, 1))
##' Sigma <- exp(-as.matrix(dist(coords)))
##' mode <- maxim_integrand(
##'   y = c(3, 5, 1), denominator = c(10, 10, 10), mu = rep(-0.5, 3),
##'   Sigma = Sigma, ID_coords = 1:3, family = "binomial"
##' )
##' mode$mode
##' @export
maxim_integrand <- function(
    y, denominator, mu, Sigma, ID_coords, ID_re = NULL, family,
    sigma2_re = NULL, hessian = FALSE, gradient = FALSE, invlink = NULL
) {
  check_supplied(y)
  check_supplied(mu)
  check_supplied(Sigma)
  check_supplied(ID_coords)
  check_supplied(family)
  if (missing(denominator)) {
    abort_input("{.arg denominator} must be given.", arg = "denominator")
  }
  family_link <- read_family(family, invlink, families = engine_families())
  family <- family_link$family
  invlink <- family_link$invlink
  # The integrand and its derivatives below refer to the denominator as
  # `units_m`, the name used by the package's internal fitting functions.
  # The records are those the engine fits, the negative binomial family's as
  # binomial ones, whose offset goes into `mu`.
  records <- engine_records(family, y, denominator, mu,
                            family_link$parameters, invlink)
  family <- records$family
  units_m <- records$units_m
  mu <- records$offset

  # ---------- utilities ----------
  sum_by_group <- function(v, grp, nlev) {
    sums <- rowsum(v, grp, reorder = FALSE)
    out <- numeric(nlev)
    out[as.integer(rownames(sums))] <- sums[, 1L]
    out
  }

  cross_sum <- function(v, grp1, n1, grp2, n2) {
    f1 <- factor(grp1, levels = seq_len(n1))
    f2 <- factor(grp2, levels = seq_len(n2))
    as.matrix(xtabs(v ~ f1 + f2))  # n1 x n2, zeros where empty
  }

  # ---------- dimensions ----------
  Sigma.inv <- solve(Sigma)
  n_loc <- nrow(Sigma)
  n <- length(y)

  if (( !is.null(ID_re) && is.null(sigma2_re)) ||
      (  is.null(ID_re) && !is.null(sigma2_re))) {
    abort_input("{.arg ID_re} and {.arg sigma2_re} must be given together.", arg = if (is.null(ID_re)) "ID_re" else "sigma2_re")
  }

  if (is.null(ID_re)) {
    n_re <- 0L
    n_dim_re <- integer(0)
    ind_re <- list()
  } else {
    n_re <- length(sigma2_re)
    n_dim_re <- vapply(seq_len(n_re), function(i) length(unique(ID_re[, i])), integer(1))
    ind_re <- vector("list", n_re)
    add_i <- 0L
    for (i in seq_len(n_re)) {
      ind_re[[i]] <- (add_i + n_loc + 1):(add_i + n_loc + n_dim_re[i])
      if (i < n_re) add_i <- sum(n_dim_re[1:i])
    }
  }
  n_tot <- n_loc + if (n_re > 0) sum(n_dim_re) else 0L

  linkf <- make_inverse_link_functions(family, invlink, n)
  model <- make_family(family, linkf)

  # Per-observation contributions for arbitrary inverse link
  contribs <- function(eta) {
    terms <- model$terms(eta, y, units_m)
    if (is.null(terms)) {
      if (!model$canonical) {
        riskmap_abort("{.arg invlink} gives a mean outside the support of the family, or a derivative that isn't finite, at the current linear predictor.", class = "riskmap_error_fit")
      }
      if (any(!is.finite(eta))) {
        riskmap_abort("The linear predictor isn't finite.", class = "riskmap_error_fit")
      }
      # Only the canonical Poisson mean, exp(eta), overflows.
      riskmap_abort("The Poisson mean isn't finite at the current linear predictor.", class = "riskmap_error_fit")
    }
    list(g = terms$score, w = terms$curvature, llik = terms$log_likelihood)
  }

  # ---------- core functions ----------
  integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]
    S_re_list <- if (n_re > 0) lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]]) else NULL

    eta <- mu + S[ID_coords]
    if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]

    cw <- contribs(eta)

    qS  <- as.numeric(crossprod(S, Sigma.inv %*% S))
    qre <- if (n_re > 0) sum(vapply(seq_len(n_re), function(i) sum(S_re_list[[i]]^2) / sigma2_re[i], numeric(1))) else 0

    -0.5 * qS - 0.5 * qre + cw$llik
  }

  grad.integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]
    S_re_list <- if (n_re > 0) lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]]) else NULL

    eta <- mu + S[ID_coords]
    if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]

    cw <- contribs(eta)
    g  <- cw$g

    out <- numeric(n_tot)
    out[1:n_loc] <- as.numeric(-Sigma.inv %*% S) + sum_by_group(g, ID_coords, n_loc)
    if (n_re > 0) {
      for (j in seq_len(n_re)) {
        out[ind_re[[j]]] <- - S_re_list[[j]] / sigma2_re[j] +
          sum_by_group(g, ID_re[, j], n_dim_re[j])
      }
    }

    out
  }

  hessian.integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]
    S_re_list <- if (n_re > 0) lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]]) else NULL

    eta <- mu + S[ID_coords]
    if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]

    cw <- contribs(eta)
    w  <- cw$w  # positive if l'' is negative

    H <- matrix(0, nrow = n_tot, ncol = n_tot)

    # S block
    H[1:n_loc, 1:n_loc] <- -Sigma.inv
    diag(H)[1:n_loc] <- diag(H)[1:n_loc] - sum_by_group(w, ID_coords, n_loc)

    # RE blocks
    if (n_re > 0) {
      for (j in seq_len(n_re)) {
        idx <- ind_re[[j]]
        # diagonal block j
        diag(H)[idx] <- diag(H)[idx] - 1 / sigma2_re[j] - sum_by_group(w, ID_re[, j], n_dim_re[j])

        # cross S vs RE j
        M <- cross_sum(w, ID_coords, n_loc, ID_re[, j], n_dim_re[j]) # n_loc x n_dim_re[j]
        H[1:n_loc, idx] <- H[1:n_loc, idx] - M
        H[idx, 1:n_loc] <- t(H[1:n_loc, idx])

        # cross RE j vs RE k
        if (j < n_re) {
          for (k in (j + 1):n_re) {
            Mk <- cross_sum(w, ID_re[, j], n_dim_re[j], ID_re[, k], n_dim_re[k])
            idxk <- ind_re[[k]]
            H[idx, idxk] <- H[idx, idxk] - Mk
            H[idxk, idx] <- t(H[idx, idxk])
          }
        }
      }
    }
    H
  }

  # ---------- optimize ----------
  start <- rep(0, n_tot)
  if (model$moves_start) {
    start <- move_start_off_limits(
      start,
      function(S_tot) {
        eta <- mu + S_tot[ID_coords]
        if (n_re > 0) for (i in seq_len(n_re)) {
          eta <- eta + S_tot[ind_re[[i]]][ID_re[, i]]
        }
        eta
      },
      y, units_m, family, linkf, ID_coords
    )
  }
  estim <- nlminb(
    start,
    function(x) -integrand(x),
    function(x) -grad.integrand(x),
    function(x) -hessian.integrand(x)
  )

  out <- list(mode = estim$par, invlink_used = linkf$name)
  if (hessian) {
    out$hessian <- hessian.integrand(out$mode)
  } else {
    out$Sigma.tilde <- solve(-hessian.integrand(out$mode))
  }
  if (gradient) out$gradient <- grad.integrand(out$mode)

  out
}

##' @title Laplace-sampling MCMC for Generalized Linear Gaussian Process Models
##'
##' @description
##' Draws the latent effects of a binomial or Poisson model given the data
##' (the spatial effects at the locations, then the levels of each `re()`
##' term) with a Metropolis-adjusted Langevin chain preconditioned by the
##' Laplace approximation at their conditional mode, with the step size
##' tuned during burn-in.
##'
##' @param y Numeric vector of responses of length \eqn{n}: the number of
##'   positives for \code{family = "binomial"}, or counts for
##'   \code{family = "poisson"}.
##' @param units_m Numeric vector of positive values, one per response: the
##'   number of trials for \code{family = "binomial"}, or the exposure for
##'   \code{family = "poisson"} (the \code{denominator} of [glgpm()]).
##' @param mu Numeric vector of length \eqn{n} providing the fixed part of the
##'   linear predictor for each observation.
##' @param Sigma Numeric positive-definite covariance matrix for the latent spatial
##'   process \eqn{S} at the unique locations referenced by \code{ID_coords}.
##' @param ID_coords Integer vector of length \eqn{n} mapping each response in
##'   \code{y} to a row/column of \code{Sigma} (i.e., the index of the corresponding location).
##' @param ID_re Optional matrix or data.frame with one column per unstructured
##'   random effect (RE). Each column is an integer vector of length \eqn{n}
##'   mapping observations in \code{y} to RE levels (e.g., cluster, survey, etc.).
##'   Use \code{NULL} to exclude REs.
##' @param sigma2_re Optional named numeric vector of RE variances. Names must
##'   match the column names of \code{ID_re}. Ignored if \code{ID_re = NULL}.
##' @param family Character string: \code{"binomial"} or \code{"poisson"}.
##' @param control_mcmc List of control parameters:
##'   \describe{
##'     \item{n_sim}{Total number of MCMC iterations (including burn-in).}
##'     \item{burnin}{Number of initial iterations to discard.}
##'     \item{thin}{Thinning interval for saving samples.}
##'     \item{h}{Initial step size for the Gaussian proposal. Defaults to \eqn{1.65 / n_\mathrm{tot}^{1/6}} if not supplied.}
##'     \item{c1.h, c2.h}{Positive tuning constants for adaptive step-size updates.}
##'   }
##' @param seed Optional integer seed for the chain. When supplied, the
##'   caller's random number state is restored on exit.
##' @param invlink Optional inverse-link function. If \code{NULL}, defaults are used:
##'   \code{identity} (gaussian), \code{plogis} (binomial), and \code{exp} (poisson).
##' @param start Optional starting value for the search for the mode, in the
##'   order of the latent effects, such as the mode at nearby parameter values.
##' @param average_step Logical. If `TRUE`, the step kept after burn-in is the
##'   geometric mean of the tuned step over the second half of burn-in, and
##'   otherwise its last value.
##' @param start_from_draw Logical. If `TRUE`, the chain starts from a draw of
##'   the Laplace approximation, and otherwise at the mode, from which a chain
##'   in many dimensions accepts almost nothing until its step has shrunk.
##' @param adapt_step Logical. If `FALSE`, the step is not tuned during burn-in
##'   and stays at `control_mcmc$h`, so that burn-in only discards draws.
##' @param family_parameters The family's own parameters (`read_family()`),
##'   NULL for a family without them. The chain samples the records the
##'   engine fits (`engine_records()`), the negative binomial family's as
##'   binomial ones.
##'
##' @details
##' The mode is found by Newton's method (`find_sampler_mode()`). The Laplace
##' precision there, Q = P + A' W A for prior precision P, record weights W
##' and the map A from latent effects to records, is factored through the
##' Schur complement of the diagonal block of the `re()` term with most
##' levels (`sampler_latent_structure()`), so that an iteration of the chain
##' (`cpp_latent_mala()`) costs two triangular solves of the size of the
##' locations and passes over the records, whatever the number of levels.
##' Every link runs in the compiled chain, a custom link through one call of
##' its R functions per iteration.
##'
##' Dimensions must be consistent:
##' \code{length(y) = n}, \code{nrow(Sigma) = ncol(Sigma) = n_loc},
##' and \code{length(ID_coords) = n} with entries in \eqn{1,\dots,n_\mathrm{loc}}.
##' If \code{ID_re} is provided, each column must have length \eqn{n}; when
##' \code{sigma2_re} is supplied, it must be named and match \code{colnames(ID_re)}.
##'
##' @return An object of class \code{"RiskMap_mcmc"} with components:
##' \describe{
##'   \item{samples}{A list containing posterior draws. Always includes
##'                 \code{$S} (latent spatial field). If \code{ID_re} is supplied,
##'                 each unstructured RE is returned under \code{$<re_name>}.}
##'   \item{tuning_par}{Numeric vector of step sizes (\code{h}) used over iterations.}
##'   \item{acceptance_prob}{Numeric vector of cumulative Metropolis–Hastings acceptance rates.}
##'   \item{acceptance}{Logical vector recording whether each proposal was accepted.}
##'   \item{acceptance_rate}{Named vector giving acceptance rates during burn-in
##'                         and retained-sample iterations.}
##'   \item{mode}{The conditional mode of the latent effects.}
##' }
##'
##' @section Default links:
##' The default inverse links are: identity (gaussian), logistic (binomial),
##' and exponential (poisson). Supply \code{invlink} to override.
##'
##' @seealso \code{\link{maxim_integrand}}
##' @noRd
laplace_sampling_mcmc <- function(y,
                                  units_m,
                                  mu,
                                  Sigma,
                                  ID_coords,
                                  ID_re = NULL,
                                  sigma2_re = NULL,
                                  family,
                                  control_mcmc,
                                  seed = NULL,
                                  invlink = NULL,
                                  start = NULL,
                                  average_step = TRUE,
                                  start_from_draw = TRUE,
                                  adapt_step = TRUE,
                                  family_parameters = NULL) {

  traits <- check_engine_family(family)

  n <- length(y)
  n_loc <- nrow(Sigma)
  if (!is.numeric(y) || n < 1L || any(!is.finite(y)) ||
      any(y < 0) || any(y != as.integer(y))) {
    abort_input("The response must be non-negative whole numbers.", arg = "data")
  }
  if (!is.numeric(units_m) || length(units_m) != n ||
      any(!is.finite(units_m)) || any(units_m <= 0)) {
    abort_input("The denominators must be positive finite numbers, one for each response.", arg = "denominator")
  }
  if (traits$upper_limit && any(y > units_m)) {
    abort_input("Binomial responses can't exceed their denominators.", arg = "denominator")
  }
  if (!is.numeric(Sigma) || length(dim(Sigma)) != 2L ||
      n_loc != ncol(Sigma) || n_loc < 1L) {
    riskmap_bug("{.arg Sigma} must be a non-empty square numeric matrix.")
  }
  if (!is.numeric(mu) || length(mu) != n || any(!is.finite(mu))) {
    riskmap_bug("{.arg mu} must contain one finite value per response.")
  }
  if (length(ID_coords) != n || anyNA(ID_coords) ||
      any(ID_coords != as.integer(ID_coords)) ||
      any(ID_coords < 1L | ID_coords > n_loc)) {
    riskmap_bug("{.arg ID_coords} must contain one valid location index per response.")
  }

  # set seed if it exists and reset on exit
  if (!is.null(seed)){
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(seed)
  }

  # ---------- dimensions ----------
  if (( !is.null(ID_re) && is.null(sigma2_re)) ||
      (  is.null(ID_re) && !is.null(sigma2_re))) {
    riskmap_bug("{.arg ID_re} and {.arg sigma2_re} must be given together.")
  }

  if (is.null(ID_re)) {
    n_re <- 0L
    n_dim_re <- integer(0)
  } else {
    ID_re <- as.matrix(ID_re)
    if (nrow(ID_re) != n || ncol(ID_re) < 1L || anyNA(ID_re) ||
        any(ID_re != as.integer(ID_re)) || any(ID_re < 1L)) {
      riskmap_bug("{.arg ID_re} must contain positive integer indices and one row per response.")
    }
    n_re <- length(sigma2_re)
    if (n_re != ncol(ID_re) || any(!is.finite(sigma2_re)) ||
        any(sigma2_re <= 0)) {
      riskmap_bug("{.arg sigma2_re} must contain one positive finite variance per random effect.")
    }
    if (!is.null(names(sigma2_re)) && !is.null(colnames(ID_re)) &&
        !identical(names(sigma2_re), colnames(ID_re))) {
      riskmap_bug("The names of {.arg sigma2_re} must match the columns of {.arg ID_re}.")
    }
    n_dim_re <- vapply(seq_len(n_re), function(i) length(unique(ID_re[, i])), integer(1))
    for (i in seq_len(n_re)) {
      if (!identical(sort(unique(ID_re[, i])), seq_len(n_dim_re[i]))) {
        riskmap_bug("Each column of {.arg ID_re} must use consecutive indices starting at 1.")
      }
    }
  }
  if (!is.null(start) && (!is.numeric(start) ||
      length(start) != n_loc + sum(n_dim_re) || any(!is.finite(start)))) {
    riskmap_bug("{.arg start} must contain one finite value for each latent effect.")
  }

  # The records the engine samples from, the negative binomial family's as
  # binomial ones, whose offset goes into the linear predictors.
  records <- engine_records(family, y, units_m, mu, family_parameters,
                            invlink)
  family <- records$family
  units_m <- records$units_m
  mu <- records$offset
  linkf <- make_inverse_link_functions(family, invlink, n)
  model <- make_family(family, linkf)
  structure <- sampler_latent_structure(Sigma, ID_coords, ID_re, sigma2_re)
  found <- tryCatch(
    find_sampler_mode(
      structure, y, units_m, mu, family, linkf,
      start = if (!is.null(start)) start[order(structure$position)]
    ),
    error = function(error) NULL
  )
  if (is.null(found)) {
    # The search of maxim_integrand() is slower but differs from Newton's
    # method, so a survey where Newton's method fails still gets a centre.
    found <- fallback_sampler_mode(structure, y, units_m, mu, Sigma,
                                   ID_coords, ID_re, sigma2_re, family,
                                   invlink, linkf)
  }
  link_terms <- if (!model$compiled) {
    # One call of the link per iteration, on all the linear predictors. A
    # point where the link fails or is not finite is rejected by the chain.
    function(eta) {
      terms <- tryCatch(
        compute_custom_link_terms(eta, y, units_m, family, linkf,
                                  with_curvature = FALSE),
        error = function(error) NULL
      )
      if (is.null(terms) || !is.finite(terms$log_likelihood)) return(NULL)
      list(terms$log_likelihood, terms$score)
    }
  }
  # The compiled chain stops with base errors, which are raised again with
  # the fit class. RiskMap's own errors, from a custom link, pass as they are.
  chain <- tryCatch(cpp_latent_mala(
    y_r = as.numeric(y), units_r = as.numeric(units_m), mu_r = as.numeric(mu),
    columns_r = structure$columns, n_u = structure$n_u, n_d = structure$n_d,
    prior_u_r = structure$prior_u, prior_d_r = structure$prior_d,
    weights_r = found$weights, mode_r = found$mode, family = family,
    link_terms = link_terms, n_sim = control_mcmc$n_sim,
    burnin = control_mcmc$burnin, thin = control_mcmc$thin,
    step_size = control_mcmc$h %||% (1.65 / structure$dimension^(1 / 6)),
    adaptation_scale = if (adapt_step) control_mcmc$c1.h else 0,
    adaptation_exponent = control_mcmc$c2.h,
    average_step = average_step, start_from_draw = start_from_draw
  ), error = function(error) {
    if (inherits(error, "riskmap_error")) riskmap_rethrow(error)
    riskmap_abort("The Langevin sampler of the latent effects failed.", class = "riskmap_error_fit", parent = error)
  })

  samples <- chain$samples[, structure$position, drop = FALSE]
  out_sim <- list(samples = list(S = samples[, seq_len(n_loc), drop = FALSE]))
  if (n_re > 0) {
    re_names <- colnames(ID_re) %||% paste0("re", seq_len(n_re))
    first <- n_loc + c(0L, cumsum(n_dim_re))
    for (i in seq_len(n_re)) {
      out_sim$samples[[re_names[i]]] <-
        samples[, first[i] + seq_len(n_dim_re[i]), drop = FALSE]
    }
  }
  out_sim$tuning_par <- chain$tuning
  out_sim$acceptance_prob <- chain$cumulative_acceptance
  out_sim$acceptance <- as.logical(chain$acceptance)
  out_sim$acceptance_rate <- c(
    burnin = if (control_mcmc$burnin > 0L) {
      mean(out_sim$acceptance[seq_len(control_mcmc$burnin)])
    } else NA_real_,
    sampling = mean(out_sim$acceptance[seq.int(
      control_mcmc$burnin + 1L, control_mcmc$n_sim
    )])
  )
  out_sim$mode <- found$mode[structure$position]
  out_sim$invlink_used <- linkf$name
  class(out_sim) <- "RiskMap_mcmc"
  out_sim
}

##' Latent effects of the prediction sampler in the order of its factor
##'
##' The sampler orders the latent effects as the locations, then the levels of
##' every `re()` term but the one with most levels, then the levels of that
##' term, whose block of the Laplace precision is diagonal (`n_d` of them). The
##' first `n_u` have the dense prior precision `prior_u` (the inverse of
##' `Sigma` for the locations, `1 / sigma2_re` on the diagonal for the levels),
##' and the last `n_d` the prior precision `prior_d`. `columns` holds, for each
##' record, the 0-based positions of the latent effects in its linear
##' predictor, with the diagonal term last, and `position` maps the usual
##' order (locations, then each term's levels in turn) to this one.
##' @noRd
sampler_latent_structure <- function(Sigma, ID_coords, ID_re, sigma2_re) {
  n_loc <- nrow(Sigma)
  precision <- chol2inv(factor_covariance(Sigma, "spatial covariance matrix"))
  n_re <- if (is.null(ID_re)) 0L else ncol(ID_re)
  levels <- if (n_re > 0) as.integer(apply(ID_re, 2L, max)) else integer(0)
  diagonal_term <- if (n_re > 0) which.max(levels) else 0L
  dense_terms <- setdiff(seq_len(n_re), diagonal_term)
  n_u <- n_loc + sum(levels[dense_terms])
  n_d <- if (n_re > 0) levels[diagonal_term] else 0L
  # First position (0-based) of each term, in the usual and in this order.
  usual_first <- n_loc + c(0L, cumsum(levels))[seq_len(n_re)]
  sampler_first <- integer(n_re)
  sampler_first[dense_terms] <- n_loc + c(0L, cumsum(levels[dense_terms]))[
    seq_along(dense_terms)]
  if (n_re > 0) sampler_first[diagonal_term] <- n_u
  position <- seq_len(n_loc)
  for (j in seq_len(n_re)) {
    position[usual_first[j] + seq_len(levels[j])] <-
      sampler_first[j] + seq_len(levels[j])
  }
  columns <- matrix(as.integer(ID_coords) - 1L, ncol = 1L)
  for (j in c(dense_terms, if (n_re > 0) diagonal_term)) {
    columns <- cbind(columns, sampler_first[j] + as.integer(ID_re[, j]) - 1L)
  }
  storage.mode(columns) <- "integer"
  dimnames(columns) <- NULL
  prior_u <- matrix(0, n_u, n_u)
  prior_u[seq_len(n_loc), seq_len(n_loc)] <- precision
  for (j in dense_terms) {
    index <- sampler_first[j] + seq_len(levels[j])
    prior_u[cbind(index, index)] <- 1 / sigma2_re[[j]]
  }
  prior_d <- if (n_re > 0) rep(1 / sigma2_re[[diagonal_term]], n_d) else
    numeric(0)
  list(columns = columns, n_u = n_u, n_d = n_d, n_loc = n_loc,
       dimension = n_u + n_d, prior_u = prior_u, prior_d = prior_d,
       position = position)
}

##' Conditional mode of the latent effects for the prediction sampler
##'
##' Newton's method with step halving on log f(y | u) - u' P u / 2, in the
##' order and with the structured factor of `sampler_latent_structure()`, so
##' that each step costs a factorisation of the size of the locations and
##' passes over the records, whatever the number of `re()` levels. It follows
##' the rules of `find_latent_mode()`. For a link that is not
##' log-concave, curvatures that leave the Hessian without a Cholesky factor
##' are set to zero, a search from a nonzero `start` that fails is repeated
##' from zero, and a custom link's search starts off the limits of the link.
##' Returns the mode and the record weights of the Laplace precision there,
##' which are positive definite by the same rule.
##' @noRd
find_sampler_mode <- function(structure, y, units_m, mu, family, linkf,
                              start = NULL, tolerance = 1e-10,
                              maximum_iterations = 100L) {
  latent_map <- make_latent_map(structure$columns + 1L, structure$dimension)
  dense <- seq_len(structure$n_u)
  diagonal <- structure$n_u + seq_len(structure$n_d)
  model <- make_family(family, linkf)
  prior_product <- function(latent) {
    c(as.numeric(structure$prior_u %*% latent[dense]),
      structure$prior_d * latent[diagonal])
  }
  terms_at <- function(latent) {
    terms <- model$terms(
      mu + latent_linear_predictor(latent_map, latent), y, units_m
    )
    finite <- !is.null(terms) && is.finite(terms$log_likelihood) &&
      all(is.finite(terms$score)) && all(is.finite(terms$curvature))
    if (finite) terms
  }
  objective <- function(latent, terms) {
    terms$log_likelihood - 0.5 * sum(latent * prior_product(latent))
  }
  solve_at <- function(weights, right_hand_side) {
    cpp_structured_solve(structure$columns, structure$n_u, structure$n_d,
                         structure$prior_u, structure$prior_d, weights,
                         right_hand_side)
  }
  # The record weights of the Newton step: the curvatures, or for a link
  # that is not log-concave with negative ones set to zero if the Hessian
  # has no Cholesky factor.
  usable_weights <- function(curvature) {
    if (model$assumes_concave ||
        !is.null(solve_at(curvature, numeric(structure$dimension)))) {
      return(curvature)
    }
    pmax(curvature, 0)
  }
  search_from <- function(latent) {
    if (model$moves_start) {
      latent <- model$move_start(
        latent,
        function(point) mu + latent_linear_predictor(latent_map, point),
        y, units_m, latent_map$columns[, 1L]
      )
    }
    terms <- terms_at(latent)
    if (is.null(terms)) return(NULL)
    current <- objective(latent, terms)
    if (!is.finite(current)) return(NULL)
    reached_limit <- TRUE
    for (iteration in seq_len(maximum_iterations)) {
      gradient <- latent_gradient(latent_map, terms$score) -
        prior_product(latent)
      step <- solve_at(usable_weights(terms$curvature), gradient)
      if (is.null(step) || !all(is.finite(step))) return(NULL)
      accepted <- FALSE
      for (halving in 0:30) {
        candidate <- latent + step / 2^halving
        candidate_terms <- terms_at(candidate)
        if (is.null(candidate_terms)) next
        value <- objective(candidate, candidate_terms)
        if (is.finite(value) && value >= current - 1e-12 * abs(current)) {
          accepted <- TRUE
          break
        }
      }
      if (!accepted) {
        # As in find_latent_mode(), a point where the promised increase
        # is within rounding is a mode that rounding keeps from confirming.
        if (sum(gradient * step) / 2 > 1e-12 * abs(current)) return(NULL)
        reached_limit <- FALSE
        break
      }
      change <- max(abs(candidate - latent))
      latent <- candidate
      terms <- candidate_terms
      current <- value
      if (change <= tolerance * max(1, max(abs(latent)))) {
        reached_limit <- FALSE
        break
      }
    }
    list(mode = latent, weights = usable_weights(terms$curvature),
         iterations = iteration, reached_limit = reached_limit)
  }
  found <- NULL
  restarted <- FALSE
  if (!is.null(start) && !all(start == 0)) {
    found <- tryCatch(search_from(start), error = function(error) NULL)
    restarted <- is.null(found)
  }
  if (is.null(found)) found <- search_from(numeric(structure$dimension))
  if (is.null(found)) {
    riskmap_abort("The search for the mode of the latent effects failed.", class = "riskmap_error_fit")
  }
  found$restarted <- restarted
  found
}

##' Centre of the prediction sampler from maxim_integrand()
##'
##' Used when `find_sampler_mode()` fails. Returns the mode found by
##' `nlminb()` in the order of the sampler, with the record weights of the
##' Laplace precision there and negative curvatures set to zero, which keeps
##' that precision positive definite.
##' @noRd
fallback_sampler_mode <- function(structure, y, units_m, mu, Sigma, ID_coords,
                                  ID_re, sigma2_re, family, invlink, linkf) {
  maximum <- maxim_integrand(
    y = y, denominator = units_m, mu = mu, Sigma = Sigma,
    ID_coords = ID_coords, ID_re = ID_re, sigma2_re = sigma2_re,
    family = family, invlink = invlink, hessian = TRUE
  )
  mode <- maximum$mode[order(structure$position)]
  latent_map <- make_latent_map(structure$columns + 1L, structure$dimension)
  terms <- make_family(family, linkf)$terms(
    mu + latent_linear_predictor(latent_map, mode), y, units_m
  )
  if (is.null(terms) || !all(is.finite(terms$curvature))) {
    riskmap_abort("The search for the mode of the latent effects failed.", class = "riskmap_error_fit")
  }
  list(mode = mode, weights = pmax(terms$curvature, 0), iterations = NA_integer_,
       reached_limit = NA, restarted = TRUE)
}

##' Set control parameters for certified MCML estimation
##'
##' Non-Gaussian models are fitted by Monte Carlo maximum likelihood with
##' importance sampling. The estimated log-likelihood is maximised on a
##' fitting sample of draws, the estimate is corrected by one Newton step with
##' the score from an independent validation sample, followed by further steps
##' on that sample when the fitting sample predicts that one step leaves a
##' residual, and the corrected estimate is certified with a second, fresh
##' validation sample. The radius of the certificate, its allowance for the
##' Monte Carlo error of that sample's score, is the larger of the radii of
##' that sample and of a replicate of the same size drawn at the same
##' estimate. The radius of one sample is estimated from its own draws, and is
##' too small when the sample misses the rare draws of very large weight that
##' dominate its error, which the replicate misses independently. The
##' certificate is an approximate 95% bound, in statistical standard errors,
##' on the Monte Carlo error of the estimate, that is its distance from the
##' maximum of the exact likelihood that the fit has located. It is not an
##' exact finite-sample guarantee, and it says nothing about whether that
##' maximum is the highest (see `check_maxima`). If the certificate exceeds
##' `tolerance`, the samples are enlarged to the size predicted to meet it and
##' the fit is repeated.
##'
##' @param initial_samples Positive integer giving the number of draws in the
##' fitting sample at the first stage. Defaults to 1,000.
##' @param validation_samples Positive integer giving the number of draws in
##' each validation sample, and in the replicate of the second, at the first
##' stage. Later stages enlarge the
##' validation samples in the same proportion as the fitting sample. Defaults
##' to 5,000.
##' @param max_samples Positive integer giving the largest permitted fitting
##' sample. Defaults to 16,000.
##' @param max_stages Positive integer giving the maximum number of stages.
##' Defaults to 4.
##' @param tolerance Positive number giving the largest accepted certificate,
##' in statistical standard errors. Defaults to 0.1.
##' @param proposal The importance distribution for the latent effects.
##' `"global"` (the default) fits a Gaussian approximation to the whole
##' conditional distribution of the latent effects, which keeps far more of
##' the draws effective in sparse data than the Laplace approximation;
##' `"laplace"` uses the Laplace approximation at the conditional mode. If the
##' global approximation cannot be built at some parameter value, the stage is
##' repeated with the Laplace approximation.
##' @param defensive_weight Number in \[0, 1) giving the share of importance
##' draws taken from the prior distribution of the latent effects instead of
##' the Gaussian approximation. This defensive mixture keeps the variance of
##' the importance weights finite, which a Gaussian proposal alone does not in
##' most spatial datasets, so that Monte Carlo errors are valid. Defaults to
##' 0.05. Set to 0 to use the Gaussian approximation alone.
##' @param laplace_start Logical or `NULL`. If `TRUE`, the Laplace
##' approximation to the likelihood is first maximised from the starting
##' values, and Monte Carlo maximum likelihood starts from that estimate. If
##' the global importance distribution cannot be built at that estimate, as
##' can happen in sparse data when it has a very large nugget, the fit starts
##' from the starting values instead. `NULL` (the default) means `TRUE`, except
##' where `integrate_noise = TRUE` integrates a term out, where it means
##' `FALSE`, because there the Laplace approximation's spurious maximum moves
##' to the spatial effect and the start from it costs more without certifying
##' more fits. With `integrate_noise = "auto"` and a nugget as the only term
##' that can be integrated out, the Laplace estimate chooses the form of the
##' first stage, which with the nugget integrated out starts from the starting
##' values.
##' @param integrate_noise Logical, or `"auto"`. If `TRUE`, a binomial or
##' Poisson model integrates its independent Gaussian terms out of the
##' likelihood numerically wherever the data allow it, and Monte Carlo maximum
##' likelihood samples the spatial effect and the terms that stay. The model,
##' its parameters and its likelihood are the same, so the option changes only
##' how the same maximum is found. Sparse surveys with a large nugget, or with
##' household effects, then keep far more effective draws. The terms are the
##' nugget, whose levels are the locations, and each `re()` term. Records at a
##' location with the same covariates, offset and random-effect levels are
##' first combined into one unit. Then three rules apply.
##' \itemize{
##'   \item A level of an `re()` term that holds a single unit, such as a
##'   household with one record or with records that do not differ, is
##'   integrated out of that unit's likelihood, and so is the nugget of a
##'   location with a single unit.
##'   \item At a location with several units the nugget is integrated out as
##'   well, through an integral over the nugget of the units' own integrals,
##'   when the units share their covariates and offset and differ only in the
##'   levels of `re()` terms that each hold a single unit, as households of a
##'   village with covariates measured at the village do.
##'   \item At any other location the nugget stays sampled.
##' }
##' Every other level stays sampled. That includes the levels of an `re()`
##' term that span several locations, such as villages or districts with
##' several sampled points, and levels that hold units with different
##' covariates, such as a household of individuals of different ages. A term
##' can be integrated out at some of its levels and sampled at others.
##' `check_fit()` reports, for each term, at how many levels it was integrated
##' out and at how many it stayed sampled. It works with any link, including a
##' custom `invlink`. Where the nugget variance is close to zero the
##' likelihood is nearly flat along it, and there the Monte Carlo estimate of
##' its effect is noisier with the nugget integrated out than sampled, so fewer
##' such fits are certified. With `"auto"` each stage of the fit integrates out
##' the nugget if its variance at the start of the stage is at least 0.1 and
##' samples it otherwise. The first stage takes the variance from the Laplace
##' estimate, or from the starting values if `laplace_start = FALSE`. Every
##' validation bank is then drawn with the nugget sampled and integrated out,
##' and the two estimates of the score are combined with weights inverse to
##' their Monte Carlo covariances, so the certificate rests on the more precise
##' form, or on both. A bank whose importance weights have a heavy tail, with a
##' Pareto k above 0.8, can understate its covariance and is left out of the
##' combination. A model with an `re()` term that can be integrated out is
##' fitted with `"auto"` as with `TRUE`, its `re()` terms and its nugget
##' integrated out at every stage wherever the rules above allow. On simulated
##' surveys of households that certified every fit, and in a fraction of the
##' time of sampling, also where the nugget variance was zero, where switching
##' the nugget alone sometimes failed. Defaults to `"auto"`, which on simulated surveys with a nugget certified
##' more fits than either form alone. A model without a nugget or `re()` terms
##' is fitted in the same way whatever the value.
##' @param check_maxima Logical. The certificate establishes Monte Carlo
##' precision at the maximum found, not that no higher maximum exists. In
##' sparse surveys an optimiser can stop at a maximum with a spatial scale far
##' below the spacing between locations while a higher one exists at a longer
##' scale, or the reverse. If `TRUE` (the default), a fitted scale below half
##' the median distance between neighbouring locations is refitted from four
##' dispersed starting values on the same Monte Carlo draws. After a fit at a
##' longer scale the likelihood is evaluated at four shorter scales, with the
##' coefficients and variance of the model with independent location effects
##' (the limit of the model as the scale goes to zero, computed exactly by
##' quadrature), and the fit is repeated from the best of them if it comes
##' within 0.25 log-likelihood units of the fit. Exact lower bounds on the
##' likelihood at short scales force a further search if they exceed the
##' fitted log-likelihood. Separate maxima are compared on independent draws
##' and the highest is kept, with a warning if it differs. A fit whose spatial
##' correlation is numerically negligible is reported as numerically
##' independent, with the scale unresolved and without standard errors or
##' certificate. The screen of short scales and the exact bounds are
##' available for the canonical links and models without unstructured random
##' effects.
##' @param pc_prior Logical. If `TRUE`, the spatial variance and scale are
##' penalised by penalised-complexity priors (Fuglstad et al., 2019) and the
##' estimate maximises the penalised likelihood on the log scale of the
##' covariance parameters. This shrinks estimates away from degenerate
##' solutions in small or sparse datasets. It is recommended for small or
##' low-prevalence surveys, especially before prediction with
##' `parameter_uncertainty = TRUE` in [predict.RiskMap()].
##' When the spatial variance is large and well identified, the penalty shrinks
##' it somewhat towards zero, so the default is `FALSE`.
##' @param pc_sigma_threshold,pc_sigma_probability The penalty sets
##' P(sigma > `pc_sigma_threshold`) = `pc_sigma_probability`, where sigma is the
##' spatial standard deviation on the scale of the linear predictor. Default
##' 1 and 0.05.
##' @param pc_scale_threshold,pc_scale_probability The penalty sets
##' P(phi < `pc_scale_threshold`) = `pc_scale_probability`, where phi is the
##' spatial scale in the units of the coordinates. By default the threshold is
##' the median distance between neighbouring locations and the probability is
##' 0.05.
##' @return A control object of class `RiskMap_control_mcml`.
##' @examples
##' control_mcml <- set_control_mcml(max_samples = 8000)
##' @seealso [glgpm()], [set_control_mcmc()]
##' @export
set_control_mcml <- function(initial_samples = 1000L,
                             validation_samples = 5000L,
                             max_samples = 16000L,
                             max_stages = 4L,
                             tolerance = 0.1,
                             proposal = c("global", "laplace"),
                             defensive_weight = 0.05,
                             laplace_start = NULL,
                             integrate_noise = "auto",
                             check_maxima = TRUE,
                             pc_prior = FALSE,
                             pc_sigma_threshold = 1,
                             pc_sigma_probability = 0.05,
                             pc_scale_threshold = NULL,
                             pc_scale_probability = 0.05) {
  check_positive_integer(initial_samples, "initial_samples")
  check_positive_integer(validation_samples, "validation_samples")
  check_positive_integer(max_samples, "max_samples")
  check_positive_integer(max_stages, "max_stages")
  check_positive_number(tolerance, "")
  proposal <- match_choice(proposal)
  if (!is.numeric(defensive_weight) || length(defensive_weight) != 1L ||
      !is.finite(defensive_weight) || defensive_weight < 0 ||
      defensive_weight >= 1) {
    supplied <- describe_value(defensive_weight)
    abort_input("{.arg defensive_weight} must be a single number from 0 up to but not including 1, not {supplied}.", arg = "defensive_weight")
  }
  if (!is.null(laplace_start)) check_logical(laplace_start)
  if (!isTRUE(integrate_noise) && !isFALSE(integrate_noise) &&
      !identical(integrate_noise, "auto")) {
    supplied <- describe_value(integrate_noise)
    abort_input("{.arg integrate_noise} must be {.code TRUE}, {.code FALSE} or {.val auto}, not {supplied}.", arg = "integrate_noise")
  }
  check_logical(check_maxima)
  check_logical(pc_prior)
  check_positive_number(pc_sigma_threshold, "")
  check_range(pc_sigma_probability, 0, 1, allow_equal = FALSE,
              arg = "pc_sigma_probability")
  if (!is.null(pc_scale_threshold)) {
    check_positive_number(pc_scale_threshold, "")
  }
  check_range(pc_scale_probability, 0, 1, allow_equal = FALSE,
              arg = "pc_scale_probability")
  if (max_samples < initial_samples) {
    abort_input("{.arg max_samples} must be at least {.arg initial_samples} ({initial_samples}), not {max_samples}.", arg = "max_samples")
  }

  structure(
    list(
      initial_samples = as.integer(initial_samples),
      validation_samples = as.integer(validation_samples),
      max_samples = as.integer(max_samples),
      max_stages = as.integer(max_stages),
      tolerance = tolerance,
      proposal = proposal,
      defensive_weight = defensive_weight,
      laplace_start = laplace_start,
      integrate_noise = integrate_noise,
      check_maxima = check_maxima,
      pc_prior = pc_prior,
      pc_sigma_threshold = pc_sigma_threshold,
      pc_sigma_probability = pc_sigma_probability,
      pc_scale_threshold = pc_scale_threshold,
      pc_scale_probability = pc_scale_probability
    ),
    class = "RiskMap_control_mcml"
  )
}

##' Use a distinct reproducible stream for final conditional sampling
##' @noRd
final_sampling_seed <- function(seed, mcml_stages) {
  # Each certified stage consumes three independent streams: the fitting
  # bank and two validation banks. Final fitted-value sampling starts after
  # them.
  offset <- 3L * mcml_stages + 1L
  as.integer(((seed - 1 + offset) %% .Machine$integer.max) + 1)
}

##' Convert a working-scale estimate into the next MCML reference point
##' @noRd
estimate_to_mcml_reference <- function(estimate, fix_tau2) {
  out <- list(
    beta = unname(estimate$beta),
    sigma2 = exp(estimate$sigma2),
    phi = exp(estimate$phi)
  )
  if (isTRUE(fix_tau2)) {
    out$tau2 <- exp(estimate$sigma2 + estimate$nu2)
  } else {
    out$tau2 <- fix_tau2
  }
  if (!is.null(estimate$sigma2_re)) {
    out$sigma2_re <- exp(estimate$sigma2_re)
  }
  out
}

##' Draw latent effects and fitted values at a non-Gaussian fitted estimate
##'
##' MCML samples belong to its importance-sampling reference distribution.
##' This helper deliberately makes a separate draw at the fitted parameters.
##' The first block is the location-level Gaussian effect; when a nugget is
##' present it is included in that block. Subsequent blocks are unstructured
##' random effects, in the same layout historically used by `S_samples`.
##' The records are the model's, with `family_parameters` for a family that
##' has them, and the fitted values are their means under `link_function`,
##' the link of the model's mean.
##'
##' @noRd
sample_at_nongaussian_fit <- function(estimate, y, units_m, D, coords,
                                      cov_offset, ID_coords, ID_re,
                                      fix_tau2, family, link_function, kappa,
                                      control_mcmc, seed = NULL, messages,
                                      family_parameters = NULL) {
  fitted_pars <- estimate_to_mcml_reference(estimate, fix_tau2)
  mu <- as.numeric(D %*% fitted_pars$beta + cov_offset)
  distances <- pairwise_distances(coords)
  Sigma <- fitted_pars$sigma2 * matern_correlation(
    distances, phi = fitted_pars$phi, kappa = kappa,
    return_sym_matrix = TRUE
  )
  diag(Sigma) <- diag(Sigma) + fitted_pars$tau2

  simulation <- laplace_sampling_mcmc(
    y = y, units_m = units_m, mu = mu, Sigma = Sigma,
    ID_coords = ID_coords, ID_re = ID_re,
    sigma2_re = fitted_pars$sigma2_re,
    family = family,
    # Canonical links are evaluated inside the compiled chain. Custom links
    # pass the inverse link with the derivatives the fit used, so that they
    # are not derived again.
    invlink = if (identical(link_function$name, "canonical")) {
      NULL
    } else {
      link_function[c("inv", "d1", "d2")]
    },
    control_mcmc = control_mcmc, seed = seed,
    family_parameters = family_parameters
  )

  eta <- sweep(
    simulation$samples$S[, ID_coords, drop = FALSE],
    2L, mu, "+"
  )
  if (!is.null(ID_re)) {
    ID_re <- as.matrix(ID_re)
    for (j in seq_len(ncol(ID_re))) {
      eta <- eta + simulation$samples[[j + 1L]][, ID_re[, j], drop = FALSE]
    }
  }

  sample_blocks <- simulation$samples
  response_mean <- colMeans(
    matrix(link_function$inv(as.numeric(eta)), nrow = nrow(eta))
  )
  response_mean <- family_traits(family)$fitted_mean(response_mean, units_m)

  list(
    samples = do.call(cbind, sample_blocks),
    linear_predictors = colMeans(eta),
    fitted_values = response_mean,
    # The sampling settings let plot_mcmc() show the iteration of each sample.
    diagnostics = c(
      simulation[c(
        "tuning_par", "acceptance_prob", "acceptance", "acceptance_rate",
        "invlink_used"
      )],
      unclass(control_mcmc)[c("n_sim", "burnin", "thin")]
    )
  )
}

##' Set Control Parameters for Simulation
##'
##' This function sets control parameters for running simulations, supporting MCMC methods
##' for glgpm models.
##'
##' @param n_sim Integer. The total number of simulations to run. Default is
##'   22000, which with the default burn-in and thinning keeps 2,000 draws. With
##'   `linear_model = TRUE` it is the number of independent draws, 2,000 by
##'   default.
##' @param burnin Integer. The number of initial simulations to discard (burn-in/warmup period). Default is 2000.
##' @param thin Integer. The interval at which simulations are recorded (thinning interval, MCMC only). Default is 10.
##' @param h Numeric. An optional positive step size for Langevin MCMC.
##' @param c1.h Numeric. Positive scale for the burn-in Robbins-Monro step-size
##'   adaptation. Default is 1.
##' @param c2.h Numeric. Exponent controlling the diminishing adaptation gain.
##'   Must be larger than 0.5 and no larger than 1. Default is 0.6.
##' @param linear_model Logical. If TRUE, sets up parameters for a linear model. Default is FALSE.
##'
##' @details
##' If \code{linear_model = TRUE}, only \code{n_sim} is required.
##'
##' During burn-in the step size of the Langevin chain is adapted towards an
##' acceptance rate of 0.574, starting from \code{h} (by default
##' \eqn{1.65 / d^{1/6}} for \eqn{d} latent effects). The step kept after
##' burn-in is the geometric mean of the adapted step over the second half of
##' burn-in.
##'
##' Seeds are not part of the control settings. Pass `seed` to [glgpm()],
##' [predict.RiskMap()], [assess_prediction()] or [assess_simulation()].
##'
##' @return A list of control parameters with class "RiskMap_control_mcmc". Contents depend on \code{sampler}:
##' \itemize{
##'   \item For "mcmc": n_sim, burnin, thin, h, c1.h, c2.h, linear_model
##' }
##'
##' @examples
##' # Default parameters (MCMC)
##' control_mcmc <- set_control_mcmc()
##'
##' # Custom MCMC parameters
##' control_mcmc <- set_control_mcmc(n_sim = 15000, burnin = 3000, thin = 20)
##'
##' @seealso \code{\link{glgpm}}
##' @importFrom Matrix Matrix forceSymmetric
##' @export
set_control_mcmc <- function(n_sim = 22000,
                            burnin = 2000,
                            thin = 10,
                            h = NULL,
                            c1.h = 1,
                            c2.h = 0.6,
                            linear_model = FALSE){

  # =============================================================================
  # LINEAR MODEL (simple case for both samplers)
  # =============================================================================

  if (linear_model) {
    # Draws of a linear model are independent, so none is discarded.
    if (missing(n_sim)) n_sim <- 2000
    check_positive_integer(n_sim, "n_sim")
    res <- list(
      n_sim = n_sim,
      linear_model = linear_model
    )
    class(res) <- "RiskMap_control_mcmc"
    return(res)
  }

  # =============================================================================
  # MCMC SAMPLER (Langevin for glgpm)
  # =============================================================================

  # Validate MCMC parameters
  check_positive_integer(n_sim, "n_sim")
  check_positive_integer(burnin, "burnin", allow_zero = TRUE)
  check_positive_integer(thin, "thin")

  if (n_sim <= burnin) {
    abort_input("{.arg n_sim} must be larger than {.arg burnin} ({burnin}), not {n_sim}.", arg = "n_sim")
  }
  if ((n_sim - burnin) %% thin != 0) {
    abort_input("{.arg thin} must divide {.code n_sim - burnin} ({n_sim - burnin}) exactly, not {thin}.", arg = "thin")
  }
  if (!is.null(h)) check_positive_number(h, "")
  check_positive_number(c1.h, "")
  check_zero_one(c2.h, "c2.h")
  if (c2.h <= 0.5) {
    abort_input("{.arg c2.h} must be larger than 0.5 and at most 1, not {c2.h}.", arg = "c2.h")
  }

  res <- list(
    n_sim = n_sim,
    burnin = burnin,
    thin = thin,
    h = h,
    c1.h = c1.h,
    c2.h = c2.h,
    linear_model = FALSE
  )

  class(res) <- "RiskMap_control_mcmc"
  return(res)
}
