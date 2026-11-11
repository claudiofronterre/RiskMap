#' Specify a model for simulation
#'
#' Defines a data-generating model without fitting it. Parameters are supplied
#' on their natural scale and remain fixed across simulations.
#'
#' @inheritParams glgpm
#' @param parameters A named list containing `beta`, `sigma2` and `phi`.
#'   `beta` has one entry per design-matrix column, including the intercept
#'   when present. Supply `tau2` when `gp(nugget = TRUE)`, `sigma2_me` for
#'   Gaussian measurement error, and `sigma2_re` for each `re()` term.
#'   Variances may be zero. `phi` is expressed in `distance_units`.
#' @return A validated `RiskMap_simulation_model` for [simulate_glgpm()]. It
#'   contains the formula, transformed model data, family, fixed parameters,
#'   denominator name, distance units, inverse link, fixed-effect design
#'   information, random-effect terms, offset, Matern smoothness and response.
#' @details The response in `formula` need not exist in `data`. CRS handling
#'   follows [glgpm()]: `model_crs` is used when supplied, otherwise
#'   longitude/latitude data are transformed to an automatically selected UTM
#'   CRS. A message reports automatic transformations.
#' @examples
#' library(sf)
#' locations <- st_as_sf(data.frame(x = c(0, 1000, 2000), y = 0),
#'                       coords = c("x", "y"), crs = 32629)
#'
#' grid <- st_as_sf(data.frame(x = c(500, 1500), y = 500),
#'                  coords = c("x", "y"), crs = 32629)
#'
#' model <- specify_glgpm(response ~ gp(),
#'                        data = locations,
#'                        family = "gaussian",
#'                        parameters = list(beta = 1,
#'                                          sigma2 = 1,
#'                                          phi = 2,
#'                                          sigma2_me = 0.1))
#'
#' @export
specify_glgpm <- function(formula, data, family, parameters,
                         denominator = NULL, model_crs = NULL,
                         distance_units = c("km", "m"), invlink = NULL) {
  distance_units <- match_choice(distance_units)
  check_data(data)
  check_formula(formula, data, response_required = FALSE)
  if (!is.null(model_crs)) {
    check_crs(model_crs)
    data <- st_transform(data, crs = model_crs)
    if (st_is_longlat(data)) {
      abort_input("{.arg model_crs} must be a projected coordinate reference system, not longitude and latitude.", arg = "model_crs")
    }
  } else if (st_is_longlat(data)) {
    model_crs <- propose_utm(data)
    data <- st_transform(data, crs = model_crs)
    riskmap_inform(c("{.arg data} are in longitude and latitude, so they are projected to EPSG:{model_crs} for simulation.", "i" = "Set {.arg model_crs} to choose the projection."))
  }
  family_link <- read_family(family, invlink)
  family <- family_link$family
  invlink <- family_link$invlink
  custom_link <- !is.null(invlink)
  inter <- interpret_formula(formula)

  sub_denominator <- substitute(denominator)
  denominator_name <- NULL
  if (!is.null(sub_denominator)) {
    if (!is.symbol(sub_denominator)) {
      abort_input("{.arg denominator} must be the unquoted name of a column of {.arg data}.", arg = "denominator")
    }
    denominator_name <- deparse(sub_denominator)
    if (!denominator_name %in% names(data)) {
      abort_input(c("{.arg denominator} must name a column of {.arg data}.", "x" = "There is no column {.val {denominator_name}}."), arg = "denominator")
    }
  }
  if (family == "gaussian" &&
      (!is.null(denominator_name) || !is.null(invlink))) {
    abort_input("{.arg denominator} and {.arg invlink} can't be given for a Gaussian model.", arg = if (is.null(invlink)) "denominator" else "invlink")
  }
  if (is.null(invlink)) {
    invlink <- switch(family, gaussian = identity, binomial = plogis, poisson = exp)
  }
  check_invlink(inverse_link_function(invlink), family)
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
              parameters = parameters, denominator = denominator_name,
              distance_units = distance_units, invlink = invlink,
              fixed_terms = fixed_terms, contrasts = attr(design, "contrasts"),
              xlevels = xlevels, re_terms = re_terms, offset = inter$offset,
              kappa = inter$gp_spec$kappa, response = inter$response,
              custom_link = custom_link)
  model_coords <- unique(coordinates_in_units(data, distance_units))
  attr(out, "distance_scale") <- if (parameters$sigma2 > 0) {
    scale_spatial_coordinates(model_coords, parameters$phi)$distance_scale
  } else {
    1
  }
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
#' @return A `RiskMap_simulation` object containing `locations`, the transformed
#'   `sf` data used for each requested location set; `samples`, numeric arrays
#'   indexed by location, simulation and component; `nsim`; the simulation
#'   `model`; the requested `what`; and `seed`. Use [simulated_data()] to obtain
#'   a dataset ready for refitting, [simulated_surface()] for a surface, and
#'   [simulated_values()] for a tidy table of selected components.
#' @details
#' The spatial process is drawn jointly at the union of the requested locations
#' using the model's Matern covariance. There is no nearest-grid approximation
#' and no conditioning on the observed responses. This differs from conditional
#' spatial prediction with [predict.RiskMap()].
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
#' Locations are transformed to the model CRS when necessary, with a message
#' reporting the conversion. Exact joint simulation uses a dense covariance
#' matrix, so memory grows quadratically with the number of unique locations.
#'
#' @examples
#' library(sf)
#' locations <- st_as_sf(data.frame(x = c(0, 1000, 2000), y = 0),
#'                       coords = c("x", "y"), crs = 32629)
#'
#' grid <- st_as_sf(data.frame(x = c(500, 1500), y = 500),
#'                  coords = c("x", "y"), crs = 32629)
#'
#' model <- specify_glgpm(response ~ gp(),
#'                        data = locations,
#'                        family = "gaussian",
#'                        parameters = list(beta = 1,
#'                                          sigma2 = 1,
#'                                          phi = 2,
#'                                          sigma2_me = 0.1))
#'
#' sim <- simulate_glgpm(model,
#'                       nsim = 2,
#'                       what = c("data", "surface"),
#'                       prediction_grid = grid,
#'                       seed = 1)
#'
#' simulated_data(sim, simulation = 1)
#'
#' simulated_values(sim, component = "spatial_effect")
#' @export
simulate_glgpm <- function(object, nsim = 1, what = "data",
                          sample_locations = NULL, prediction_grid = NULL,
                          seed = NULL) {
  check_positive_integer(nsim, "nsim")
  if (!is.character(what) || !length(what) || anyNA(what) ||
      any(!what %in% c("data", "surface")) || anyDuplicated(what)) {
    abort_input("{.arg what} must be {.val data}, {.val surface} or both.", arg = "what")
  }
  model <- simulation_model(object)
  if (!"data" %in% what && !is.null(sample_locations)) {
    abort_input("{.arg sample_locations} is used only when {.arg what} includes {.val data}.", arg = "sample_locations")
  }
  if (!"surface" %in% what && !is.null(prediction_grid)) {
    abort_input("{.arg prediction_grid} is used only when {.arg what} includes {.val surface}.", arg = "prediction_grid")
  }
  prepared <- prepare_simulation(model, what, sample_locations,
                                 prediction_grid)
  samples <- draw_simulation(prepared, nsim, seed)
  out <- list(locations = prepared$locations, samples = samples, nsim = nsim,
              model = model, what = names(prepared$locations), seed = seed)
  class(out) <- "RiskMap_simulation"
  out
}

#' Prepare the locations, inputs and covariance factor of a simulation
#'
#' Everything `simulate_glgpm()` needs before drawing: the locations in the
#' model CRS, the fixed part of the linear predictor and the denominators at
#' each location set, the indices of the unique coordinates, the group labels
#' of the `re()` terms and the factor of the spatial covariance. This step
#' uses sf. `draw_simulation()` does not, so the draws can run in forked
#' processes.
#'
#' @param model A `RiskMap_simulation_model`.
#' @inheritParams simulate_glgpm
#' @return A list used by `draw_simulation()`.
#' @noRd
prepare_simulation <- function(model, what, sample_locations = NULL,
                               prediction_grid = NULL) {
  locations <- list()
  original <- is.null(sample_locations)
  if ("data" %in% what) {
    locations$data <- if (original) model$data else sample_locations
  }
  if ("surface" %in% what) {
    if (is.null(prediction_grid)) {
      abort_input("{.arg prediction_grid} must be given to simulate a surface.", arg = "prediction_grid")
    }
    locations$surface <- prediction_grid
  }
  for (name in names(locations)) {
    locations[[name]] <- simulation_locations(locations[[name]],
                                               st_crs(model$data), name)
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
  indexed <- index_coordinate_rows(coords)
  unique_coords <- data.frame(x = indexed$coordinates[, 1],
                              y = indexed$coordinates[, 2])
  ids <- indexed$index
  n_unique <- nrow(unique_coords)
  pars <- model$parameters
  if (pars$sigma2 > 0) {
    distance_scale <- attr(model, "distance_scale") %||% 1
    fitting_coords <- unique_coords[, c("x", "y")] / distance_scale
    fitting_phi <- pars$phi / distance_scale
    covariance <- pars$sigma2 * matern_correlation(
      pairwise_distances(fitting_coords),
      phi = fitting_phi, kappa = model$kappa,
      return_sym_matrix = TRUE
    )
    root <- tryCatch(t(factor_covariance(
      covariance, "spatial covariance of the simulation locations"
    )), error = function(e) {
      riskmap_abort(c("The spatial covariance of the simulation locations can't be factorised.", "i" = "Check the spatial parameters and look for nearly coincident locations."), class = "riskmap_error_fit", parent = e)
    })
  }
  group_labels <- lapply(model$re_terms, function(term) {
    unlist(lapply(locations, function(x) as.character(x[[term]])))
  })
  names(group_labels) <- model$re_terms
  list(model = model, locations = locations, inputs = inputs, sizes = sizes,
       ids = ids, n_unique = n_unique,
       root = if (model$parameters$sigma2 > 0) root,
       group_labels = group_labels)
}

#' Draw simulations from a prepared simulation
#'
#' @param prepared Output from `prepare_simulation()`.
#' @inheritParams simulate_glgpm
#' @return The `samples` of a `RiskMap_simulation`: one array per location
#'   set, indexed by location, simulation and component.
#' @param call The call named in input errors.
#' @noRd
draw_simulation <- function(prepared, nsim, seed = NULL,
                            call = rlang::caller_env()) {
  model <- prepared$model
  inputs <- prepared$inputs
  sizes <- prepared$sizes
  ids <- prepared$ids
  n_unique <- prepared$n_unique
  root <- prepared$root
  pars <- model$parameters
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE,
                         call = call)
  if (!is.null(seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
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
    labels <- prepared$group_labels[[term]]
    groups <- match(labels, unique(labels))
    draws <- matrix(rnorm(length(unique(labels)) * nsim,
                          sd = sqrt(pars$sigma2_re[[term]])),
                    nrow = length(unique(labels)), ncol = nsim)
    grouped <- grouped + draws[groups, , drop = FALSE]
  }
  eta <- spatial + nugget + grouped + unlist(lapply(inputs, `[[`, "fixed"),
                                             use.names = FALSE)
  mu <- check_invlink(inverse_link_function(model$invlink), model$family, eta)
  mu <- matrix(mu, nrow(eta), nsim)
  samples <- list()
  start <- 0L
  for (name in names(sizes)) {
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
        riskmap_abort(c("The simulated responses aren't all finite.", "i" = "Check the means and the denominators."), class = "riskmap_error_fit")
      }
      values[[length(values) + 1L]] <- response
      components <- c(components, "response")
    }
    samples[[name]] <- array(unlist(values, use.names = FALSE),
                             dim = c(length(rows), nsim, length(components)),
                             dimnames = list(NULL, NULL, components))
    start <- start + sizes[[name]]
  }
  samples
}

#' Extract simulated values in tidy form
#'
#' @param object Output from [simulate_glgpm()].
#' @param component Components to extract. Choose from `"spatial_effect"`,
#'   `"nugget"`, `"group_effect"`, `"linear_predictor"`, `"mean"` and, for
#'   data simulations, `"response"`. Defaults to all available components.
#' @param simulation Simulation numbers to extract. Defaults to all simulations.
#' @param location_set `"data"`, `"surface"`, or both. Defaults to available sets.
#' @return A tibble with `simulation`, `location_set`, `location_id`, `component`
#'   and `value`. Location identifiers refer to rows in `object$locations`.
#' @examples
#' library(sf)
#' locations <- st_as_sf(data.frame(x = c(0, 1000, 2000), y = 0),
#'                       coords = c("x", "y"), crs = 32629)
#'
#' grid <- st_as_sf(data.frame(x = c(500, 1500), y = 500),
#'                  coords = c("x", "y"), crs = 32629)
#'
#' model <- specify_glgpm(response ~ gp(),
#'                        data = locations,
#'                        family = "gaussian",
#'                        parameters = list(beta = 1,
#'                                          sigma2 = 1,
#'                                          phi = 2,
#'                                          sigma2_me = 0.1))
#'
#' sim <- simulate_glgpm(model,
#'                       nsim = 2,
#'                       what = c("data", "surface"),
#'                       prediction_grid = grid,
#'                       seed = 1)
#'
#' simulated_values(sim, component = "spatial_effect")
#' @export
simulated_values <- function(object, component = NULL, simulation = NULL,
                             location_set = NULL) {
  check_simulation(object, simulation)
  if (is.null(simulation)) simulation <- seq_len(object$nsim)
  if (is.null(location_set)) location_set <- names(object$samples)
  if (!is.character(location_set) || !length(location_set) || anyNA(location_set) ||
      any(!location_set %in% names(object$samples)) || anyDuplicated(location_set)) {
    available_sets <- names(object$samples)
    abort_input("{.arg location_set} must select among the simulated location sets, {.val {available_sets}}.", arg = "location_set")
  }
  available <- unique(unlist(lapply(object$samples[location_set], function(x) dimnames(x)[[3]])))
  if (!is.null(component) && (!is.character(component) || !length(component) ||
      anyNA(component) || any(!component %in% available) || anyDuplicated(component))) {
    abort_input("{.arg component} must select among the simulated components, {.val {available}}.", arg = "component")
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
#' @param simulation A single simulation number, defaulting to one. These
#'   accessors return one `sf` object at a time; use [simulated_values()] to
#'   extract multiple simulations in one tidy table.
#' @return `simulated_data()` returns the sample `sf` data with its response
#'   column replaced by the selected simulation, ready for refitting.
#'   `simulated_surface()` returns the grid `sf` data with columns for the
#'   spatial effect, nugget, group effect, linear predictor and mean, with
#'   class `RiskMap_simulated_surface` so that it can be mapped with
#'   [plot()][plot.RiskMap_simulated_surface()].
#' @examples
#' library(sf)
#' locations <- st_as_sf(data.frame(x = c(0, 1000, 2000), y = 0),
#'                       coords = c("x", "y"), crs = 32629)
#'
#' grid <- st_as_sf(data.frame(x = c(500, 1500), y = 500),
#'                  coords = c("x", "y"), crs = 32629)
#'
#' model <- specify_glgpm(response ~ gp(),
#'                        data = locations,
#'                        family = "gaussian",
#'                        parameters = list(beta = 1,
#'                                          sigma2 = 1,
#'                                          phi = 2,
#'                                          sigma2_me = 0.1))
#'
#' sim <- simulate_glgpm(model,
#'                       nsim = 2,
#'                       what = c("data", "surface"),
#'                       prediction_grid = grid,
#'                       seed = 1)
#'
#' simulated_data(sim, simulation = 1)
#' @export
simulated_data <- function(object, simulation = 1) {
  check_simulation(object, simulation, single = TRUE)
  if (is.null(object$samples$data)) {
    abort_input(c("{.arg object} has no simulated responses.", "i" = "Simulate with {.arg what} including {.val data}."), arg = "object")
  }
  out <- object$locations$data
  response <- object$model$response
  if (is.null(response) || length(response) != 1L ||
      !identical(make.names(response), response)) {
    abort_input(c("The response of the model formula must be a plain column name to be put in a data set.", "i" = "Use {.fn simulated_values} to extract the responses."), arg = "object")
  }
  out[[response]] <- as.numeric(object$samples$data[, simulation, "response"])
  out
}

#' @rdname simulated_data
#' @export
simulated_surface <- function(object, simulation = 1) {
  check_simulation(object, simulation, single = TRUE)
  if (is.null(object$samples$surface)) {
    abort_input(c("{.arg object} has no simulated surface.", "i" = "Simulate with {.arg what} including {.val surface}."), arg = "object")
  }
  out <- object$locations$surface
  components <- dimnames(object$samples$surface)[[3]]
  if (any(components %in% names(out))) {
    clashes <- intersect(components, names(out))
    abort_input(c("The prediction grid has column{?s} {.val {clashes}} named like simulated components.", "i" = "Use {.fn simulated_values} to keep them apart."), arg = "object")
  }
  for (name in components) {
    out[[name]] <- as.numeric(object$samples$surface[, simulation, name])
  }
  class(out) <- c("RiskMap_simulated_surface", class(out))
  out
}

#' Print a simulation model or simulations
#'
#' Prints a short description of a model specification from
#' [specify_glgpm()] or of simulations from [simulate_glgpm()], instead of
#' every stored value.
#'
#' @param x Output from [specify_glgpm()] or [simulate_glgpm()].
#' @param ... Not used.
#' @return Invisibly returns `x`.
#' @method print RiskMap_simulation_model
#' @export
print.RiskMap_simulation_model <- function(x, ...) {
  cat("Geostatistical model for simulation\n")
  link <- if (is.list(x$invlink)) x$invlink$link
  cat(sprintf("Family: %s%s\n", x$family,
              if (!is.null(link)) {
                sprintf(" (%s link)", link)
              } else if (isTRUE(x$custom_link)) {
                " (custom inverse link)"
              } else {
                ""
              }))
  cat(sprintf("Formula: %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf("Locations: %d\n", nrow(x$data)))
  print_simulation_parameters(x$parameters, x$distance_units)
  invisible(x)
}

#' @rdname print.RiskMap_simulation_model
#' @method print RiskMap_simulation
#' @export
print.RiskMap_simulation <- function(x, ...) {
  cat(sprintf("Simulations from a %s geostatistical model\n", x$model$family))
  cat(sprintf("Formula: %s\n", paste(deparse(x$model$formula), collapse = " ")))
  cat(sprintf("Number of simulations: %d%s\n", x$nsim,
              if (is.null(x$seed)) "" else sprintf(" (seed %s)", x$seed)))
  for (name in x$what) {
    cat(sprintf("Simulated %s at %d locations\n",
                if (name == "data") "data" else "surface",
                nrow(x$locations[[name]])))
  }
  print_simulation_parameters(x$model$parameters, x$model$distance_units)
  cat("Extract them with simulated_data(), simulated_surface() or",
      "simulated_values().\n")
  invisible(x)
}

#' Print the parameters of a simulation model
#' @noRd
print_simulation_parameters <- function(parameters, distance_units) {
  format_values <- function(values) {
    labels <- names(values) %||% paste0("[", seq_along(values), "]")
    paste(sprintf("%s = %s", labels, signif(values, 4)), collapse = ", ")
  }
  cat("Parameters:\n")
  cat(sprintf("  beta: %s\n", format_values(parameters$beta)))
  spatial <- c(sigma2 = parameters$sigma2, phi = parameters$phi)
  cat(sprintf("  %s (phi in %s)\n", format_values(spatial), distance_units))
  others <- unlist(parameters[intersect(c("tau2", "sigma2_me"),
                                        names(parameters))])
  others <- others[others != 0]
  if (length(others)) cat(sprintf("  %s\n", format_values(others)))
  if (length(parameters$sigma2_re)) {
    cat(sprintf("  sigma2_re: %s\n", format_values(parameters$sigma2_re)))
  }
  invisible(NULL)
}

#' Validate a simulation object and requested simulation indices
#'
#' @param object A prospective `RiskMap_simulation` object.
#' @param simulation Simulation indices, or `NULL` when permitted.
#' @param single Whether exactly one simulation must be selected.
#' @return `TRUE` invisibly, or an informative error.
#' @noRd
check_simulation <- function(object, simulation, single = FALSE) {
  if (!inherits(object, "RiskMap_simulation") || is.null(object$samples)) {
    abort_input("{.arg object} must be made by {.fn simulate_glgpm}, not {.obj_type_friendly {object}}.", arg = "object")
  }
  if (is.null(simulation) && !single) return(invisible(NULL))
  if (!is.numeric(simulation) || !length(simulation) || anyNA(simulation) ||
      any(!is.finite(simulation)) || any(simulation != floor(simulation)) ||
      any(simulation < 1 | simulation > object$nsim) || anyDuplicated(simulation) ||
      (single && length(simulation) != 1L)) {
    problem <- if (single) {
      "{.arg simulation} must be one simulation number between 1 and {object$nsim}."
    } else {
      "{.arg simulation} must be distinct simulation numbers between 1 and {object$nsim}."
    }
    abort_input(problem, arg = "simulation")
  }
}

#' Validate and align simulation locations
#'
#' @param data An `sf` point object.
#' @param crs The projected model CRS.
#' @param label A user-facing name for the location set.
#' @return `data`, transformed to `crs` when necessary.
#' @noRd
simulation_locations <- function(data, crs, label) {
  arg <- simulation_location_arg(label)
  check_data(data)
  if (nrow(data) == 0L || any(st_is_empty(data))) {
    abort_input("{.arg {arg}} must have at least one point and no empty geometries.", arg = arg)
  }
  if (is.na(crs)) {
    abort_input("The simulation model must have a coordinate reference system.", arg = "object")
  }
  if (!isTRUE(st_crs(data) == crs)) {
    data <- st_transform(data, crs)
    riskmap_inform("{.arg {arg}} is projected to the coordinate reference system of the model.")
  }
  coords <- st_coordinates(data)
  if (ncol(coords) != 2L || any(!is.finite(coords))) {
    abort_input("{.arg {arg}} must have finite two-dimensional coordinates.", arg = arg)
  }
  data
}

#' The argument that gives a set of simulation locations
#' @noRd
simulation_location_arg <- function(label) {
  switch(label, data = "sample_locations", surface = "prediction_grid", label)
}

#' Validate simulation parameters
#'
#' @param parameters Named parameter list supplied by the user or a fitted model.
#' @param beta_names Expected fixed-effect coefficient names.
#' @param re_terms Random-effect grouping terms.
#' @param family Model family.
#' @param nugget Nugget setting from `gp()`.
#' @return A completed and consistently ordered parameter list.
#' @noRd
simulation_parameters <- function(parameters, beta_names, re_terms, family, nugget) {
  if (!is.list(parameters) || is.null(names(parameters)) || anyNA(names(parameters)) ||
      any(names(parameters) == "") || anyDuplicated(names(parameters))) {
    abort_input("{.arg parameters} must be a list with distinct names.", arg = "parameters")
  }
  unknown <- setdiff(names(parameters), c("beta", "sigma2", "phi", "tau2",
                                         "sigma2_me", "sigma2_re"))
  if (length(unknown)) {
    abort_input("{.arg parameters} has unknown element{?s} {.val {unknown}}.", arg = "parameters")
  }
  required <- c("beta", "sigma2", "phi", if (isTRUE(nugget)) "tau2",
                 if (family == "gaussian") "sigma2_me", if (length(re_terms)) "sigma2_re")
  missing <- required[vapply(required, function(x) is.null(parameters[[x]]), logical(1))]
  if (length(missing)) {
    abort_input("{.arg parameters} must have the element{?s} {.val {missing}}.", arg = "parameters")
  }
  if (is.null(parameters$tau2)) {
    parameters$tau2 <- if (is.numeric(nugget)) nugget else 0
  }
  if (is.null(parameters$sigma2_me)) parameters$sigma2_me <- 0
  for (name in c("sigma2", "phi", "tau2", "sigma2_me")) {
    value <- parameters[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || (name == "phi" && value == 0)) {
      kind <- if (name == "phi") "positive" else "non-negative"
      abort_input("{.code parameters${name}} must be a single finite {kind} number.", arg = "parameters")
    }
  }
  if (!isTRUE(nugget)) {
    fixed <- if (is.numeric(nugget)) nugget else 0
    if (parameters$tau2 != fixed) {
      abort_input("{.code parameters$tau2} must agree with the nugget set in {.fn gp}.", arg = "parameters")
    }
  }
  if (family != "gaussian" && parameters$sigma2_me != 0) {
    abort_input("{.code parameters$sigma2_me} can be given only for a Gaussian model.", arg = "parameters")
  }
  beta <- parameters$beta
  if (!is.numeric(beta) || length(beta) != length(beta_names) || any(!is.finite(beta))) {
    abort_input(c("{.code parameters$beta} must have one finite value for each column of the design matrix, including any intercept.", "i" = "The columns are {.val {beta_names}}."), arg = "parameters")
  }
  if (!is.null(names(beta))) {
    if (anyDuplicated(names(beta)) || !setequal(names(beta), beta_names)) {
      abort_input("The names of {.code parameters$beta} must be those of the columns of the design matrix, {.val {beta_names}}.", arg = "parameters")
    }
    beta <- beta[beta_names]
  }
  parameters$beta <- setNames(beta, beta_names)
  if (length(re_terms)) {
    re <- parameters$sigma2_re
    if (!is.numeric(re) || length(re) != length(re_terms) || any(!is.finite(re)) || any(re < 0)) {
      abort_input("{.code parameters$sigma2_re} must have one non-negative variance for each {.fn re} term.", arg = "parameters")
    }
    if (!is.null(names(re))) {
      if (anyDuplicated(names(re)) || !setequal(names(re), re_terms)) {
        abort_input("The names of {.code parameters$sigma2_re} must be those of the {.fn re} terms, {.val {re_terms}}.", arg = "parameters")
      }
      re <- re[re_terms]
    }
    parameters$sigma2_re <- setNames(re, re_terms)
  } else if (length(parameters$sigma2_re)) {
    abort_input("{.code parameters$sigma2_re} can be given only for a model with {.fn re} terms.", arg = "parameters")
  }
  parameters
}

#' Convert a fitted or specified model to a simulation specification
#'
#' @param object A fitted `RiskMap` model or `RiskMap_simulation_model`.
#' @return A `RiskMap_simulation_model`.
#' @noRd
simulation_model <- function(object) {
  if (inherits(object, "RiskMap_simulation_model")) return(object)
  if (!inherits(object, "RiskMap")) {
    abort_input("{.arg object} must be a {.cls RiskMap} fit or a model from {.fn specify_glgpm}, not {.obj_type_friendly {object}}.", arg = "object")
  }
  pars <- coef(object)
  if (is.null(pars$tau2)) pars$tau2 <- object$fix_tau2
  if (object$family == "gaussian" && !is.null(object$fix_var_me)) {
    pars$sigma2_me <- object$fix_var_me
  }
  denominator <- if (!is.null(object$call$denominator)) {
    as.character(object$call$denominator)
  } else {
    NULL
  }
  # A named or custom link carries its derivatives and name, and a canonical
  # link its stored inverse.
  invlink <- refit_invlink(object) %||% object$link_function$inv
  model <- specify_glgpm(object$formula, object$data, object$family, pars,
                         distance_units = object$distance_units,
                         invlink = invlink)
  attr(model, "distance_scale") <- attr(object, "distance_scale") %||%
    attr(model, "distance_scale")
  model$denominator <- denominator
  simulation_inputs(model, model$data, TRUE, "data", original = TRUE)
  model$original_denominator <- object$units_m
  model$original_offset <- object$cov_offset
  model$custom_link <- identical(object$link_function$name, "custom")
  model
}

#' The inverse link of a simulation model as a function
#'
#' A simulation model stores `invlink` as a function, or in the list form
#' that also carries derivatives and the name of a link (`read_family()`).
#' @noRd
inverse_link_function <- function(invlink) {
  if (is.list(invlink)) invlink$inv else invlink
}

#' Check that every variable the model formula references (covariates,
#' `gp()`'s explicit coordinate columns, and `re()` grouping variables) is
#' present in one simulation location set
#'
#' @param formula The model formula.
#' @param data Location data to check.
#' @param label `"data"` or `"surface"`, naming the location set internally;
#'   mapped to the corresponding `simulate_glgpm()` argument name for the
#'   error message.
#' @return Invisible `TRUE`, or an error naming the missing variable(s).
#' @noRd
check_simulation_variables <- function(formula, data, label) {
  required <- setdiff(unique(get_formula_terms(formula)), all.vars(formula[[2L]]))
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    argument <- simulation_location_arg(label)
    abort_input("{.arg {argument}} must have the column{?s} {.val {missing}} of the model formula.", arg = argument)
  }
  invisible(TRUE)
}

#' Construct model inputs at one simulation location set
#'
#' @param model A `RiskMap_simulation_model`.
#' @param data Location data in the model CRS.
#' @param response Whether responses will be generated at these locations.
#' @param label User-facing location-set name.
#' @param original Whether these are the model's original data.
#' @return A list containing fixed linear-predictor values and denominators.
#' @noRd
simulation_inputs <- function(model, data, response, label, original = FALSE) {
  check_simulation_variables(model$formula, data, label)
  mf <- tryCatch(model.frame(model$fixed_terms, data = data, na.action = na.fail,
                             xlev = model$xlevels), error = function(e) {
    arg <- simulation_location_arg(label)
    abort_input("{.arg {arg}} has invalid covariates.", arg = arg, parent = e)
  })
  design <- model.matrix(model$fixed_terms, mf, contrasts.arg = model$contrasts)
  if (!identical(colnames(design), names(model$parameters$beta)) || any(!is.finite(design))) {
    arg <- simulation_location_arg(label)
    abort_input("{.arg {arg}} gives a design matrix that doesn't match the model.", arg = arg)
  }
  prepare_random_effects(data, model$re_terms)
  offset <- rep(0, nrow(data))
  if (!is.null(model$offset)) {
    if (!model$offset %in% names(data)) {
      arg <- simulation_location_arg(label)
      abort_input("{.arg {arg}} must have the offset column {.val {model$offset}}.", arg = arg)
    }
    offset <- data[[model$offset]]
  }
  if (original && !is.null(model$original_offset)) offset <- model$original_offset
  if (!is.numeric(offset) || length(offset) != nrow(data) || any(!is.finite(offset))) {
    arg <- simulation_location_arg(label)
    abort_input("The offsets in {.arg {arg}} must be finite numbers.", arg = arg)
  }
  den <- rep(1, nrow(data))
  if (response && !is.null(model$denominator)) {
    if (!model$denominator %in% names(data)) {
      arg <- simulation_location_arg(label)
      abort_input("{.arg {arg}} must have the denominator column {.val {model$denominator}}.", arg = arg)
    }
    den <- data[[model$denominator]]
  }
  if (response && original && !is.null(model$original_denominator)) {
    den <- model$original_denominator
  }
  if (response && (!is.numeric(den) || length(den) != nrow(data) ||
      any(!is.finite(den)) || any(den < 0) ||
      (model$family == "binomial" && any(den != floor(den))))) {
    arg <- simulation_location_arg(label)
    kind <- if (model$family == "binomial") "whole numbers" else "numbers"
    abort_input("The denominators in {.arg {arg}} must be finite non-negative {kind}.", arg = arg)
  }
  list(fixed = as.numeric(design %*% model$parameters$beta) + offset,
       denominator = den)
}

#' Adapt unified simulations for the assessment interface
#'
#' The true surface is the simulated mean, which `draw_simulation()` computes
#' through the generating model's own inverse link, whether canonical, named
#' or custom, so the assessment needs no knowledge of that link.
#'
#' @param object A `RiskMap_simulation` containing data and surface draws.
#' @return A list with the simulated datasets (`data_sim`), the prediction
#'   grid (`grid`), the simulated mean at each grid location with one column
#'   per simulation (`mean_sim`), and the family, distance units and whether
#'   the generating model has a nugget.
#' @noRd
simulation_assessment_data <- function(object) {
  check_simulation(object, NULL)
  if (!all(c("data", "surface") %in% object$what)) {
    abort_input("{.arg obj_sim} must be simulated with {.code what = c(\"data\", \"surface\")}.", arg = "obj_sim")
  }
  model <- object$model
  if (length(model$re_terms)) {
    abort_input(c("{.fn assess_simulation} can't assess simulations from a model with {.fn re} terms.", "i" = "The simulated values are available through {.fn simulated_values}."), arg = "obj_sim")
  }
  den <- simulation_inputs(model, object$locations$data, TRUE, "data",
                            original = identical(object$locations$data, model$data))$denominator
  data <- lapply(seq_len(object$nsim), function(i) {
    x <- simulated_data(object, i)
    x$y <- as.numeric(object$samples$data[, i, "response"])
    x$units_m <- den
    x
  })
  mean_sim <- object$samples$surface[, , "mean"]
  dim(mean_sim) <- c(nrow(object$locations$surface), object$nsim)
  list(data_sim = data, grid = object$locations$surface, mean_sim = mean_sim,
       family = model$family, distance_units = model$distance_units,
       nugget_over_grid = model$parameters$tau2 > 0)
}
