#' Specify a model for simulation
#'
#' Defines a data-generating model without fitting it. Parameters are supplied
#' on their natural scale and remain fixed across simulations.
#'
#' @param formula Model formula containing `gp()` and optionally `re()` and
#'   `offset()`, as in [glgpm()]. The response need not exist in `data`.
#' @param data An `sf` point object containing covariates and grouping variables.
#' @param family One of `"gaussian"`, `"binomial"`, or `"poisson"`.
#' @param parameters A named list containing `beta`, `sigma2` and `phi`.
#'   `beta` has one entry per design-matrix column, including the intercept
#'   when present. Supply `tau2` when `gp(nugget = TRUE)`, `sigma2_me` for
#'   Gaussian measurement error, and `sigma2_re` for each `re()` term.
#'   Variances may be zero. `phi` is expressed in `distance_units`.
#' @param denominator Optional column name, supplied as a character string,
#'   for binomial trial totals or Poisson exposure. Defaults to one.
#' @param distance_units Units of spatial distances and `phi`, `"km"` or `"m"`.
#' @param invlink Optional inverse-link function for binomial or Poisson
#'   models. Defaults to `plogis` or `exp`, respectively.
#' @return A validated `RiskMap_simulation_model` for [simulate_glgpm()].
#' @details `data` must have a projected CRS. Transform longitude/latitude
#'   coordinates explicitly before specifying the model. Its geometry retains
#'   the CRS units, independently of the requested distance units.
#' @export
specify_glgpm <- function(formula, data, family, parameters,
                         denominator = NULL, distance_units = c("km", "m"),
                         invlink = NULL) {
  distance_units <- match.arg(distance_units)
  simulation_locations(data, st_crs(data), "data")
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a model formula.", call. = FALSE)
  }
  family <- match.arg(family, c("gaussian", "binomial", "poisson"))
  custom_link <- !is.null(invlink)
  inter <- interpret_formula(formula)
  if (is.null(inter$gp_spec)) {
    stop("'formula' must contain gp().", call. = FALSE)
  }
  if (!is.null(denominator) &&
      (!is.character(denominator) || length(denominator) != 1L ||
       is.na(denominator) || !denominator %in% names(data))) {
    stop("'denominator' must name a column in 'data'.", call. = FALSE)
  }
  if (family == "gaussian" && (!is.null(denominator) || !is.null(invlink))) {
    stop("Gaussian models use the identity link and no denominator.", call. = FALSE)
  }
  if (is.null(invlink)) {
    invlink <- switch(family, gaussian = identity, binomial = plogis, poisson = exp)
  }
  if (!is.function(invlink)) {
    stop("'invlink' must be a function.", call. = FALSE)
  }
  fixed_terms <- terms(update(inter$pf, NULL ~ .))
  mf <- model.frame(fixed_terms, data = data, na.action = na.fail)
  design <- model.matrix(fixed_terms, mf)
  # Keep prediction information for transformed covariates and factor contrasts.
  fixed_terms <- attr(mf, "terms")
  factor_columns <- vapply(mf, function(x) is.factor(x) || is.character(x), logical(1))
  xlevels <- lapply(mf[factor_columns], function(x) levels(as.factor(x)))
  re_terms <- inter$re_spec$term
  prepare_random_effects(data, re_terms)
  parameters <- simulation_parameters(parameters, colnames(design),
                                      re_terms, family, inter$gp_spec$nugget)
  out <- list(formula = formula, data = data, family = family,
              parameters = parameters, denominator = denominator,
              distance_units = distance_units, invlink = invlink,
              fixed_terms = fixed_terms, contrasts = attr(design, "contrasts"),
              xlevels = xlevels, re_terms = re_terms, offset = inter$offset,
              kappa = inter$gp_spec$kappa, response = inter$response,
              custom_link = custom_link)
  class(out) <- "RiskMap_simulation_model"
  # Validate data-dependent settings before accepting the specification.
  simulation_inputs(out, data, TRUE, "data")
  out
}

#' Simulate replicated data and spatial surfaces
#'
#' Generates unconditional realizations from a fitted or specified model,
#' keeping its parameters fixed. This supports parametric bootstrap, replicated
#' data checks and simulation studies with a known spatial surface.
#'
#' @param object A fitted `RiskMap` object or a model from [specify_glgpm()].
#' @param nsim Number of independent simulations. One is useful for illustrating
#'   the mechanism. Reliable bootstrap intervals and performance estimates need
#'   many simulations and an assessment of Monte Carlo uncertainty.
#' @param what `"data"`, `"surface"`, or `c("data", "surface")`.
#' @param sample_locations Optional `sf` point data for simulated responses.
#'   Defaults to the model's original data. New data must contain the model's
#'   covariates, grouping variables, offsets and any denominator column.
#' @param prediction_grid `sf` point data containing the covariates, grouping
#'   variables and offsets for the requested surface. Required for `"surface"`.
#' @param seed Optional non-negative integer seed. When supplied, the caller's
#'   random-number state is restored on exit.
#' @return A `RiskMap_simulation` object with location data stored once and
#'   numeric simulation arrays. Use [simulated_data()] to obtain a dataset
#'   ready for refitting, [simulated_surface()] for a surface, and
#'   [simulated_values()] for a tidy table of selected simulation components.
#' @details
#' The spatial process is drawn jointly at the union of the requested locations
#' using the model's Matern covariance. There is no nearest-grid approximation
#' and no conditioning on the observed responses. This differs from conditional
#' spatial prediction with [setup_prediction()].
#'
#' The returned components distinguish the spatial effect `S`, the nugget and
#' grouped effects, the full linear predictor (including offsets), its
#' inverse-link `mean`, and the `response`. For binomial models `mean` is a
#' probability, and for Poisson models it is a rate per unit exposure.
#' Denominators or exposures enter response generation. Gaussian measurement
#' error enters only the response. A surface does not include sampled responses.
#'
#' Repeated coordinates share the spatial effect and location-level nugget.
#' Matching group labels share a newly simulated group effect, also across
#' sample and grid locations. Measurement error and response sampling are
#' independent across observations conditional on the latent effects.
#'
#' Locations must use the model's projected CRS. Transform them explicitly
#' when necessary. Exact joint simulation uses a dense covariance matrix, so
#' memory grows quadratically with the number of unique locations.
#'
#' @examples
#' library(sf)
#' locations <- st_as_sf(data.frame(x = c(0, 1000, 2000), y = 0),
#'                       coords = c("x", "y"), crs = 32629)
#' model <- specify_glgpm(response ~ gp(), locations, "gaussian",
#'                        parameters = list(beta = 1, sigma2 = 1, phi = 2,
#'                                          sigma2_me = 0.1))
#' sim <- simulate_glgpm(model, nsim = 2, what = c("data", "surface"),
#'                       prediction_grid = locations, seed = 1)
#' simulated_data(sim, simulation = 1)
#' simulated_values(sim, component = "spatial_effect")
#' @export
simulate_glgpm <- function(object, nsim = 1, what = "data",
                          sample_locations = NULL, prediction_grid = NULL,
                          seed = NULL) {
  if (is.numeric(nsim) && any(!is.finite(nsim))) {
    stop("'nsim' must be a single positive integer.", call. = FALSE)
  }
  check_positive_integer(nsim, "nsim")
  if (!is.character(what) || !length(what) || anyNA(what) ||
      any(!what %in% c("data", "surface")) || anyDuplicated(what)) {
    stop("'what' must be 'data', 'surface', or c('data', 'surface').", call. = FALSE)
  }
  model <- simulation_model(object)
  if (!"data" %in% what && !is.null(sample_locations)) {
    stop("'sample_locations' requires 'data' in 'what'.", call. = FALSE)
  }
  if (!"surface" %in% what && !is.null(prediction_grid)) {
    stop("'prediction_grid' requires 'surface' in 'what'.", call. = FALSE)
  }
  locations <- list()
  original <- is.null(sample_locations)
  if ("data" %in% what) {
    locations$data <- if (original) model$data else sample_locations
  }
  if ("surface" %in% what) {
    if (is.null(prediction_grid)) {
      stop("Provide 'prediction_grid' when requesting a surface.", call. = FALSE)
    }
    locations$surface <- prediction_grid
  }
  for (name in names(locations)) {
    simulation_locations(locations[[name]], st_crs(model$data), name)
  }
  inputs <- lapply(names(locations), function(name) {
    simulation_inputs(model, locations[[name]], name == "data", name,
                      original = name == "data" && original)
  })
  names(inputs) <- names(locations)
  sizes <- vapply(locations, nrow, integer(1))
  coords <- do.call(rbind, lapply(locations, coordinates_in_units,
                                distance_units = model$distance_units))
  # Deduplicate exact coordinates before factorising the joint covariance.
  coord_table <- data.frame(x = coords[, 1], y = coords[, 2])
  unique_coords <- unique(coord_table)
  unique_coords$coord_id <- seq_len(nrow(unique_coords))
  coord_table$row_id <- seq_len(nrow(coord_table))
  mapping <- merge(coord_table, unique_coords, by = c("x", "y"), sort = FALSE)
  ids <- mapping$coord_id[order(mapping$row_id)]
  n_unique <- nrow(unique_coords)
  pars <- model$parameters
  if (pars$sigma2 > 0) {
    covariance <- pars$sigma2 * matern_correlation(dist(unique_coords[, c("x", "y")]),
                                                  phi = pars$phi, kappa = model$kappa,
                                                  return_sym_matrix = TRUE)
    root <- tryCatch(t(chol(covariance)), error = function(e) {
      stop("The spatial covariance could not be factorised. Check spatial ",
           "parameters and nearly coincident locations. ", conditionMessage(e),
           call. = FALSE)
    })
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed) || seed < 0 || seed != floor(seed) ||
        seed > .Machine$integer.max) {
      stop("'seed' must be a non-negative integer or NULL.", call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  spatial <- matrix(0, n_unique, nsim)
  if (pars$sigma2 > 0) {
    spatial <- root %*% matrix(rnorm(n_unique * nsim), n_unique, nsim)
  }
  nugget <- matrix(rnorm(n_unique * nsim, sd = sqrt(pars$tau2)), n_unique, nsim)
  spatial <- spatial[ids, , drop = FALSE]
  nugget <- nugget[ids, , drop = FALSE]
  grouped <- matrix(0, sum(sizes), nsim)
  for (term in model$re_terms) {
    labels <- unlist(lapply(locations, function(x) as.character(x[[term]])))
    groups <- match(labels, unique(labels))
    draws <- matrix(rnorm(length(unique(labels)) * nsim,
                          sd = sqrt(pars$sigma2_re[[term]])),
                    nrow = length(unique(labels)), ncol = nsim)
    grouped <- grouped + draws[groups, , drop = FALSE]
  }
  eta <- spatial + nugget + grouped + unlist(lapply(inputs, `[[`, "fixed"),
                                             use.names = FALSE)
  mu <- model$invlink(as.numeric(eta))
  if (!is.numeric(mu) || length(mu) != length(eta) || any(!is.finite(mu)) ||
      (model$family == "binomial" && any(mu < 0 | mu > 1)) ||
      (model$family == "poisson" && any(mu < 0))) {
    stop("The inverse link returned invalid means for the model family.", call. = FALSE)
  }
  mu <- matrix(mu, nrow(eta), nsim)
  samples <- list()
  start <- 0L
  for (name in names(locations)) {
    rows <- start + seq_len(sizes[[name]])
    components <- c("spatial_effect", "nugget", "group_effect",
                    "linear_predictor", "mean")
    values <- list(spatial[rows, , drop = FALSE], nugget[rows, , drop = FALSE],
                   grouped[rows, , drop = FALSE], eta[rows, , drop = FALSE],
                   mu[rows, , drop = FALSE])
    if (name == "data") {
      means <- mu[rows, , drop = FALSE]
      exposure <- inputs[[name]]$denominator
      response <- switch(model$family,
                         gaussian = as.numeric(means) + rnorm(length(means),
                                                             sd = sqrt(pars$sigma2_me)),
                         binomial = rbinom(length(means), size = exposure,
                                           prob = as.numeric(means)),
                         poisson = rpois(length(means), lambda = as.numeric(means) * exposure))
      if (any(!is.finite(response))) {
        stop("Response simulation produced non-finite values. Check means and exposures.",
             call. = FALSE)
      }
      values[[length(values) + 1L]] <- response
      components <- c(components, "response")
    }
    samples[[name]] <- array(unlist(values, use.names = FALSE),
                             dim = c(length(rows), nsim, length(components)),
                             dimnames = list(NULL, NULL, components))
    start <- start + sizes[[name]]
  }
  out <- list(locations = locations, samples = samples, nsim = nsim,
              model = model, what = names(locations), seed = seed)
  class(out) <- "RiskMap_simulation"
  out
}

#' Extract simulated values in tidy form
#'
#' @param object Output from [simulate_glgpm()].
#' @param component Components to extract. Defaults to all available components.
#' @param simulation Simulation numbers to extract. Defaults to all simulations.
#' @param location_set `"data"`, `"surface"`, or both. Defaults to available sets.
#' @return A tibble with `simulation`, `location_set`, `location_id`, `component`
#'   and `value`. Location identifiers refer to rows in `object$locations`.
#' @export
simulated_values <- function(object, component = NULL, simulation = NULL,
                             location_set = NULL) {
  simulation_selection(object, simulation)
  if (is.null(simulation)) simulation <- seq_len(object$nsim)
  if (is.null(location_set)) location_set <- names(object$samples)
  if (!is.character(location_set) || !length(location_set) || anyNA(location_set) ||
      any(!location_set %in% names(object$samples)) || anyDuplicated(location_set)) {
    stop("'location_set' must select available data or surface locations.", call. = FALSE)
  }
  available <- unique(unlist(lapply(object$samples[location_set], function(x) dimnames(x)[[3]])))
  if (!is.null(component) && (!is.character(component) || !length(component) ||
      anyNA(component) || any(!component %in% available) || anyDuplicated(component))) {
    stop("'component' must select available simulation components.", call. = FALSE)
  }
  tables <- lapply(location_set, function(name) {
    x <- object$samples[[name]]
    selected <- if (is.null(component)) dimnames(x)[[3]] else intersect(component, dimnames(x)[[3]])
    index <- expand.grid(location_id = seq_len(dim(x)[1]), simulation = simulation,
                          component = selected, KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)
    index$location_set <- rep(name, nrow(index))
    index$value <- as.numeric(x[, simulation, selected, drop = FALSE])
    index[, c("simulation", "location_set", "location_id", "component", "value")]
  })
  as_tibble(do.call(rbind, tables))
}

#' Extract a simulated dataset or surface
#'
#' @param object Output from [simulate_glgpm()].
#' @param simulation A single simulation number, defaulting to one.
#' @return `simulated_data()` returns the sample `sf` data with its response
#'   column replaced by the selected simulation, ready for refitting.
#'   `simulated_surface()` returns the grid `sf` data with columns for the
#'   spatial effect, nugget, group effect, linear predictor and mean.
#' @export
simulated_data <- function(object, simulation = 1) {
  simulation_selection(object, simulation, single = TRUE)
  if (is.null(object$samples$data)) {
    stop("No responses were simulated. Include 'data' in 'what'.", call. = FALSE)
  }
  out <- object$locations$data
  response <- object$model$response
  if (is.null(response) || length(response) != 1L ||
      !identical(make.names(response), response)) {
    stop("Refitting requires a simple response name in the model formula. ",
         "Use simulated_values() to extract responses.", call. = FALSE)
  }
  out[[response]] <- as.numeric(object$samples$data[, simulation, "response"])
  out
}

#' @rdname simulated_data
#' @export
simulated_surface <- function(object, simulation = 1) {
  simulation_selection(object, simulation, single = TRUE)
  if (is.null(object$samples$surface)) {
    stop("No surface was simulated. Include 'surface' in 'what'.", call. = FALSE)
  }
  out <- object$locations$surface
  components <- dimnames(object$samples$surface)[[3]]
  if (any(components %in% names(out))) {
    stop("Grid columns conflict with simulation component names. ",
         "Use simulated_values() to keep them separate.", call. = FALSE)
  }
  for (name in components) {
    out[[name]] <- as.numeric(object$samples$surface[, simulation, name])
  }
  out
}

#' @noRd
simulation_selection <- function(object, simulation, single = FALSE) {
  if (!inherits(object, "RiskMap_simulation") || is.null(object$samples)) {
    stop("'object' must be output from simulate_glgpm().", call. = FALSE)
  }
  if (is.null(simulation) && !single) return(invisible(NULL))
  if (!is.numeric(simulation) || !length(simulation) || anyNA(simulation) ||
      any(!is.finite(simulation)) || any(simulation != floor(simulation)) ||
      any(simulation < 1 | simulation > object$nsim) || anyDuplicated(simulation) ||
      (single && length(simulation) != 1L)) {
    stop("'simulation' must select valid simulation numbers",
         if (single) " (one at a time)" else "", ".", call. = FALSE)
  }
}

#' @noRd
simulation_locations <- function(data, crs, label) {
  if (!inherits(data, "sf") || nrow(data) == 0L ||
      any(st_is_empty(data)) || any(st_geometry_type(data) != "POINT")) {
    stop("'", label, "' must be a non-empty sf object containing POINT geometries.", call. = FALSE)
  }
  if (is.na(st_crs(data)) || is.na(crs)) {
    stop("'", label, "' must have a known CRS.", call. = FALSE)
  }
  if (st_is_longlat(data) || !isTRUE(st_crs(data) == crs)) {
    stop("'", label, "' must use the model's projected CRS. ",
         "Use st_transform() before simulation.", call. = FALSE)
  }
  coords <- st_coordinates(data)
  if (ncol(coords) != 2L || any(!is.finite(coords))) {
    stop("'", label, "' must contain finite two-dimensional coordinates.", call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
simulation_parameters <- function(parameters, beta_names, re_terms, family, nugget) {
  if (!is.list(parameters) || is.null(names(parameters)) || anyNA(names(parameters)) ||
      any(names(parameters) == "") || anyDuplicated(names(parameters))) {
    stop("'parameters' must be a named list with unique names.", call. = FALSE)
  }
  unknown <- setdiff(names(parameters), c("beta", "sigma2", "phi", "tau2",
                                         "sigma2_me", "sigma2_re"))
  if (length(unknown)) {
    stop("Unknown simulation parameters: ", paste(unknown, collapse = ", "), ".",
         call. = FALSE)
  }
  required <- c("beta", "sigma2", "phi", if (isTRUE(nugget)) "tau2",
                 if (family == "gaussian") "sigma2_me", if (length(re_terms)) "sigma2_re")
  missing <- required[vapply(required, function(x) is.null(parameters[[x]]), logical(1))]
  if (length(missing)) {
    stop("Missing simulation parameters: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(parameters$tau2)) {
    parameters$tau2 <- if (is.numeric(nugget)) nugget else 0
  }
  if (is.null(parameters$sigma2_me)) parameters$sigma2_me <- 0
  for (name in c("sigma2", "phi", "tau2", "sigma2_me")) {
    value <- parameters[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || (name == "phi" && value == 0)) {
      stop("'", name, "' must be finite and ",
           if (name == "phi") "positive." else "non-negative.", call. = FALSE)
    }
  }
  if (!isTRUE(nugget)) {
    fixed <- if (is.numeric(nugget)) nugget else 0
    if (parameters$tau2 != fixed) {
      stop("'tau2' must agree with the nugget setting in gp().", call. = FALSE)
    }
  }
  if (family != "gaussian" && parameters$sigma2_me != 0) {
    stop("'sigma2_me' applies only to Gaussian models.", call. = FALSE)
  }
  beta <- parameters$beta
  if (!is.numeric(beta) || length(beta) != length(beta_names) || any(!is.finite(beta))) {
    stop("'beta' must have one finite value per design-matrix column, including ",
         "the intercept when present. Expected: ", paste(beta_names, collapse = ", "),
         ".", call. = FALSE)
  }
  if (!is.null(names(beta))) {
    if (anyDuplicated(names(beta)) || !setequal(names(beta), beta_names)) {
      stop("Names of 'beta' must match the design-matrix columns.", call. = FALSE)
    }
    beta <- beta[beta_names]
  }
  parameters$beta <- setNames(beta, beta_names)
  if (length(re_terms)) {
    re <- parameters$sigma2_re
    if (!is.numeric(re) || length(re) != length(re_terms) || any(!is.finite(re)) || any(re < 0)) {
      stop("'sigma2_re' must contain one non-negative variance per re() term.", call. = FALSE)
    }
    if (!is.null(names(re))) {
      if (anyDuplicated(names(re)) || !setequal(names(re), re_terms)) {
        stop("Names of 'sigma2_re' must match the re() terms.", call. = FALSE)
      }
      re <- re[re_terms]
    }
    parameters$sigma2_re <- setNames(re, re_terms)
  } else if (length(parameters$sigma2_re)) {
    stop("'sigma2_re' requires re() terms in the formula.", call. = FALSE)
  }
  parameters
}

#' @noRd
simulation_model <- function(object) {
  if (inherits(object, "RiskMap_simulation_model")) return(object)
  if (!inherits(object, "RiskMap")) {
    stop("'object' must be a fitted RiskMap model or output from specify_glgpm().", call. = FALSE)
  }
  pars <- coef(object)
  if (is.null(pars$tau2)) pars$tau2 <- object$fix_tau2
  if (object$family == "gaussian" && !is.null(object$fix_var_me)) {
    pars$sigma2_me <- object$fix_var_me
  }
  denominator <- if (!is.null(object$call$den)) as.character(object$call$den) else NULL
  invlink <- if (object$family == "gaussian") NULL else object$link_function$inv
  model <- specify_glgpm(object$formula, object$data, object$family, pars,
                         denominator = denominator, distance_units = object$distance_units,
                         invlink = invlink)
  model$original_denominator <- object$units_m
  model$original_offset <- object$cov_offset
  model$custom_link <- identical(object$link_function$name, "custom")
  model
}

#' @noRd
simulation_inputs <- function(model, data, response, label, original = FALSE) {
  mf <- tryCatch(model.frame(model$fixed_terms, data = data, na.action = na.fail,
                             xlev = model$xlevels), error = function(e) {
    stop("Invalid covariates in '", label, "': ", conditionMessage(e), call. = FALSE)
  })
  design <- model.matrix(model$fixed_terms, mf, contrasts.arg = model$contrasts)
  if (!identical(colnames(design), names(model$parameters$beta)) || any(!is.finite(design))) {
    stop("The design matrix for '", label, "' does not match the model.", call. = FALSE)
  }
  prepare_random_effects(data, model$re_terms)
  offset <- rep(0, nrow(data))
  if (!is.null(model$offset)) {
    if (!model$offset %in% names(data)) {
      stop("Missing offset column '", model$offset, "' in '", label, "'.", call. = FALSE)
    }
    offset <- data[[model$offset]]
  }
  if (original && !is.null(model$original_offset)) offset <- model$original_offset
  if (!is.numeric(offset) || length(offset) != nrow(data) || any(!is.finite(offset))) {
    stop("Offsets in '", label, "' must be finite numeric values.", call. = FALSE)
  }
  den <- rep(1, nrow(data))
  if (response && !is.null(model$denominator)) {
    if (!model$denominator %in% names(data)) {
      stop("Missing denominator column '", model$denominator, "' in '", label, "'.", call. = FALSE)
    }
    den <- data[[model$denominator]]
  }
  if (response && original && !is.null(model$original_denominator)) {
    den <- model$original_denominator
  }
  if (response && (!is.numeric(den) || length(den) != nrow(data) ||
      any(!is.finite(den)) || any(den < 0) ||
      (model$family == "binomial" && any(den != floor(den))))) {
    stop("Denominators in '", label, "' must be finite, non-negative",
         if (model$family == "binomial") " integers." else " numbers.", call. = FALSE)
  }
  list(fixed = as.numeric(design %*% model$parameters$beta) + offset,
       denominator = den)
}

#' @noRd
simulation_assessment_data <- function(object) {
  simulation_selection(object, NULL)
  if (!all(c("data", "surface") %in% object$what)) {
    stop("Assessment requires what = c('data', 'surface').", call. = FALSE)
  }
  model <- object$model
  if (length(model$re_terms) || isTRUE(model$custom_link)) {
    stop("The current assess_simulation() interface does not support generating ",
         "models with re() effects or custom inverse links. Simulated values ",
         "remain available through the simulation accessors.", call. = FALSE)
  }
  den <- simulation_inputs(model, object$locations$data, TRUE, "data",
                            original = identical(object$locations$data, model$data))$denominator
  data <- lapply(seq_len(object$nsim), function(i) {
    x <- simulated_data(object, i)
    x$y <- as.numeric(object$samples$data[, i, "response"])
    x$units_m <- den
    x
  })
  grid <- object$locations$surface
  for (i in seq_len(object$nsim)) {
    grid[[paste0("lp_sim_", i)]] <- as.numeric(object$samples$surface[, i, "linear_predictor"])
  }
  list(data_sim = data, lp_grid_sim = grid, family = model$family,
       distance_units = model$distance_units, include_covariates = TRUE,
       include_cov_offset = !is.null(model$offset),
       nugget_over_grid = model$parameters$tau2 > 0)
}
