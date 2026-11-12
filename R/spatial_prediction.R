##' Spatial prediction from a fitted RiskMap model
##'
##' Predicts the latent spatial process, linear predictor and response-scale
##' distribution at observed or new locations. New locations and every variable
##' required by the fitted model are supplied together in one `sf` object,
##' keeping their rows aligned.
##'
##' @param object A fitted `RiskMap` object returned by [glgpm()].
##' @param newdata An `sf` object with point geometries and the variables needed
##'   by the fitted model. When `NULL`, predictions are made at the observed
##'   locations.
##' @param boundaries Optional `sf` polygon object. When supplied, prediction
##'   locations are assigned to polygons and their joint samples are aggregated
##'   within each polygon. Grids returned by [create_grid()] use exact
##'   fractional cell coverage; other point objects use point-in-polygon
##'   assignment.
##' @param aggregate Character vector naming one or more within-polygon
##'   operations: `"mean"`, `"median"`, `"sum"`, `"min"`, or `"max"`. It must
##'   be supplied when `boundaries` is used. Several operations reuse the same
##'   joint predictive samples. Custom operations can be supplied in a named
##'   list, optionally mixed with built-ins, for example
##'   `list(mean = "mean", q90 = function(values, weights) ...)`. A custom
##'   function receives a location-by-draw matrix and the location weights and
##'   must return one numeric value per draw.
##' @param type Character string specifying `"auto"`, `"marginal"`, or
##'   `"joint"` prediction. `"auto"` uses marginal prediction for locations and
##'   joint prediction when `boundaries` are supplied.
##' @param include `NULL`, or a character vector specifying optional components
##'   to add to the complete linear predictor and response. Supported values are
##'   `"random_effects"`, for conditional predictions at known groups, and
##'   `"nugget"`, for a new realization of observation-level predictive
##'   variation at locations supplied through `newdata`. Fixed effects, the
##'   spatial process and any offset are always returned separately.
##' @param transform Optional named list of functions defining derived predictive
##'   targets from the complete linear-predictor draws. Each function receives a
##'   location-by-draw matrix and must return a numeric matrix with the same
##'   dimensions. The ordinary response-scale prediction is always calculated,
##'   so this argument is needed only for additional scientific targets, for
##'   example `list(stunting = function(link) pnorm((-2 - link) / sd_child))`.
##' @param standardize Optional data frame describing a common reference
##'   population over which response-scale predictions are averaged. Its columns
##'   must be fixed-effect predictors from the fitted model and override the
##'   corresponding columns in `newdata`. Rows receive equal weight unless an
##'   optional numeric `.weight` column is supplied. Predictor values kept
##'   together on each row are standardized jointly, preserving their empirical
##'   association. This answers what would be predicted if every location had
##'   the same reference-population distribution. `.weight` weights reference
##'   people or profiles and is distinct from the polygon `weights` argument.
##' @param summaries Character vector naming the summaries calculated from each
##'   predictive distribution. Supported values are `"mean"`, `"median"`,
##'   `"sd"`, `"cv"`, `"lower"`, `"upper"`, `"interval_width"`, and
##'   `"exceedance"`. The coefficient of variation is calculated only on the
##'   response scale and is returned as `NA` when the predictive mean is
##'   non-positive or numerically too close to zero for a stable ratio. Custom
##'   summaries can be supplied in a named list, optionally mixed with
##'   built-ins, for example `list(mean = "mean", iqr = IQR)`. Each custom
##'   function receives the predictive draws and must return one numeric value.
##' @param thresholds Numeric thresholds used when `summaries` includes
##'   `"exceedance"`. One exceedance probability is returned for each threshold
##'   on the response scale and for every target supplied through `transform`.
##' @param level Number strictly between zero and one giving the probability
##'   covered by the `"lower"` and `"upper"` interval limits and used to
##'   calculate `"interval_width"`.
##' @param weights `NULL` for equally weighted polygon summaries, or the name of
##'   a numeric column in `newdata` containing non-negative location weights.
##'   For `aggregate = "mean"` or `"median"`, weights are interpreted within
##'   each polygon. A weighted median is returned when weights are supplied.
##'   For `aggregate = "sum"`, their original scale is retained. Weights do not
##'   affect `"min"` or `"max"`.
##' @param control_mcmc Optional output from [set_control_mcmc()]. When omitted,
##'   conditional samples stored in a non-Gaussian fit are reused. Supplying a
##'   control object requests a fresh conditional run.
##' @param parameter_uncertainty Logical. If `FALSE` (the default), predictions
##'   are made at the estimated parameters (plug-in prediction). If `TRUE`,
##'   they also reflect the uncertainty in the estimated parameters; see
##'   Details.
##' @param parameter_draws Positive integer giving the number of parameter
##'   values drawn when `parameter_uncertainty = TRUE`. Default 50. Conditional
##'   samples stored in the fit are not used in that case, because they
##'   condition on the estimate only; `control_mcmc`, if supplied, sets only
##'   the total number of predictive draws.
##' @param seed Optional integer seed for random draws made during prediction,
##'   including a fresh conditional MCMC run. The caller's random-number state
##'   is restored on exit.
##' @param keep_samples Logical. If `TRUE`, retain the minimal latent samples
##'   needed for custom post-processing. Predictive samples are generated to
##'   calculate summaries regardless of this setting, but are discarded by
##'   default to limit memory use.
##' @param messages Logical. Display messages, and a progress bar while
##'   predicting at many locations, as described in [riskmap_progress].
##'   Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.
##' @param ... Additional arguments reserved for future prediction methods.
##'
##' @details
##' At new locations, the default prediction represents the underlying spatial
##' surface, with linear predictor
##' \deqn{d(x)^T \hat\beta + S(x) + \mathrm{offset}(x).}
##' The offset is included automatically when it was specified in the fitted
##' model. Nugget and unstructured random effects are excluded unless requested
##' explicitly because they generally represent observation-level or
##' group-specific variation rather than the underlying spatial surface. The
##' observed-location distinction for a fitted nugget is described below.
##'
##' For non-Gaussian fits, stored final-parameter conditional samples are reused
##' when `control_mcmc = NULL`. If the model was fitted with
##' `return_samples = FALSE`, the function draws a new conditional sample using
##' [set_control_mcmc()] defaults. Supplying `control_mcmc` always requests a
##' fresh conditional run. For Gaussian models, prediction is based on their
##' exact conditional Gaussian distribution.
##'
##' Plug-in prediction treats the estimated parameters as known, so its
##' intervals can be too narrow when the parameters, in particular the spatial
##' variance and scale, are poorly estimated, as in small or low-prevalence
##' surveys. With `parameter_uncertainty = TRUE`, `parameter_draws` parameter
##' values are drawn from the approximate sampling distribution of the
##' estimate: a multivariate normal distribution with the estimated covariance
##' matrix, on the scale on which the model is estimated (log variances and
##' log scale). For each value, conditional samples of the latent effects are
##' generated at that value and the prediction is repeated; the draws are then
##' pooled, so every summary, aggregation and derived target reflects both the
##' latent and the parameter uncertainty. For binomial and Poisson models the
##' conditional samples at each parameter value come from a short run of the
##' same Langevin (MALA) sampler used for plug-in prediction, in coordinates
##' whitened by the Laplace approximation at that value, started from a draw
##' of that approximation and using the step size tuned at the estimate (100
##' burn-in iterations, then thinning by 10). The total number of predictive
##' draws equals that of plug-in prediction, 2,000 by default. Parameters
##' without a finite positive variance, such as those held at a boundary, keep
##' their estimated value. When the covariance matrix of the estimates is not
##' available, an error is returned; refitting with
##' `set_control_mcml(pc_prior = TRUE)` usually resolves this. The location
##' column `fixed` then reports the average fixed-effect contribution.
##'
##' At the observed locations, the conditioned location-level Gaussian effect
##' already includes the fitted nugget when the model has one. Consequently,
##' `include = "nugget"` is available only for predictions supplied through
##' `newdata`, where it represents a new independent nugget realization.
##'
##' `newdata` and `boundaries` are transformed internally to the fitted model
##' CRS without changing row order; returned geometries remain in their input
##' CRS. Rows with missing required predictors are retained with `NA`
##' predictions and an informative warning.
##'
##' When `boundaries` are supplied, one prediction run produces both the
##' location-level output in `predictions` and the areal output in
##' `boundary_summary`. Joint samples are generated in the smallest independent
##' blocks that preserve dependence within every polygon. Polygons sharing a
##' fractionally covered grid cell are placed in the same block so that the
##' cell's draw is reused rather than simulated twice. This avoids constructing
##' a full joint covariance over unrelated prediction locations. Locations
##' outside all polygons remain in `predictions` but do not contribute to areal
##' aggregation; polygons with no usable locations remain in `boundaries` with
##' `NA` summaries.
##'
##' With location weights \eqn{w_i}, `"mean"` computes
##' \deqn{\sum_i w_i T_i / \sum_i w_i,}
##' so population values are converted internally to within-polygon fractions.
##' `"sum"` instead computes \eqn{\sum_i w_i T_i} and never normalizes the
##' weights. Thus, a prevalence prediction weighted by population counts gives
##' population-weighted prevalence under `"mean"` and expected cases under
##' `"sum"`. Population densities must first be multiplied by the represented
##' cell areas. When `newdata` was created by [create_grid()], each supplied
##' weight is multiplied by the exact fraction of its square cell covered by
##' each polygon. Cells crossing a boundary can therefore contribute to more
##' than one polygon. For an arbitrary `sf` point object, cell footprints cannot
##' be inferred safely and each point is assigned to the polygon containing it.
##'
##' For location-level prediction, output columns use `<quantity>_<summary>`
##' names. For example, `response_mean`, `response_lower`, and `response_upper`
##' describe the response-scale prediction, while the corresponding
##' `spatial_*` and `link_*` columns describe the spatial process and complete
##' linear predictor. `fixed` is the estimated fixed-effects contribution and
##' `offset` is reported separately. When requested through `include`,
##' `random_effects_*` and `nugget_*` report those contributions separately as
##' well as incorporating them into `link_*` and `response_*`. Binomial
##' response predictions are probabilities and Poisson response predictions
##' are rates per unit exposure.
##'
##' Derived targets supplied through `transform` are evaluated draw by draw
##' before any location-level summary or polygon aggregation. Their output uses
##' the target name as the quantity prefix. With multiple thresholds,
##' location-level exceedance columns are named, for example,
##' `response_exceedance_0.2` and `response_exceedance_0.5`; polygon output
##' stores thresholds in the tidy `threshold` column.
##'
##' When `standardize` is supplied, every conditional spatial draw is reused
##' across all rows of the reference population. The response and any derived
##' targets are evaluated for each reference profile and then averaged using
##' `.weight`. This is direct standardization (empirical or counterfactual
##' averaging), not simulation of uncertain predictors. The reported `fixed`
##' and `link_*` quantities use the weighted mean fixed-effect contribution;
##' with a nonlinear link, transforming that mean is not the same as the
##' standardized `response_*` quantity.
##'
##' @return A `RiskMap_prediction` object containing:
##' \describe{
##'   \item{predictions}{The aligned prediction points or grid cells with their
##'     summaries.}
##'   \item{boundaries}{The supplied polygons, or `NULL` when no polygon
##'     aggregation was requested.}
##'   \item{boundary_summary}{For polygon prediction, a tidy table with one row per
##'     boundary, aggregation, component, summary and threshold. `boundary_id`
##'     is the row number of the corresponding polygon in `boundaries`.}
##'   \item{samples}{`NULL` when `keep_samples = FALSE`; otherwise a list
##'     containing spatial samples and any requested random-effect or nugget
##'     samples for location-level prediction, and, with parameter
##'     uncertainty, the fixed-effect contribution of each draw (`fixed`).
##'     For polygon aggregation it is nested first by aggregation and then by
##'     component, with one matrix row per polygon.}
##'   \item{include}{Optional components included in the link and response.}
##'   \item{transform}{Names of additional derived predictive targets.}
##'   \item{standardized}{Whether response predictions were averaged over a
##'     reference population supplied through `standardize`.}
##'   \item{aggregate}{Names of the within-polygon aggregations.}
##'   \item{summaries}{Names of the calculated summaries.}
##'   \item{thresholds}{Thresholds used for exceedance probabilities.}
##'   \item{level}{Probability covered by the interval limits.}
##'   \item{weights}{Name of the weights column used for the polygon
##'     summaries, or `NULL`.}
##'   \item{type}{Whether prediction was marginal or joint.}
##'   \item{parameter_draws}{The number of parameter draws when parameter
##'     uncertainty was propagated, otherwise `NULL`.}
##'   \item{family}{The fitted response family.}
##'   \item{response_name}{Name of the response variable, used to label
##'     maps of Gaussian predictions.}
##'   \item{observed_locations}{Whether `newdata` was omitted.}
##'   \item{coverage}{`NULL` for location-only prediction. For polygon
##'     aggregation, `"exact"` when fractional grid-cell coverage was used and
##'     `"point"` for point-in-polygon assignment. This is output metadata, not
##'     a user-supplied argument.}
##' }
##'
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:50, ]
##' fit <- glgpm(y ~ pop_dens + gp(), data = example_data, family = "gaussian",
##'              messages = FALSE)
##'
##' # Calculate an additional predictive target from every linear-predictor draw
##' prediction <- predict(
##'   fit,
##'   newdata = example_data[1:5, ],
##'   transform = list(squared_link = function(link) link^2),
##'   messages = FALSE
##' )
##' prediction$predictions
##'
##' # Also reflect the uncertainty in the estimated parameters
##' uncertain_prediction <- predict(
##'   fit,
##'   newdata = example_data[1:5, ],
##'   parameter_uncertainty = TRUE,
##'   parameter_draws = 20,
##'   messages = FALSE
##' )
##' uncertain_prediction$predictions
##'
##' # Use the same reference distribution of population density at each location
##' reference_population <- data.frame(pop_dens = example_data$pop_dens)
##' standardized_prediction <- predict(
##'   fit,
##'   newdata = example_data[1:5, ],
##'   standardize = reference_population,
##'   messages = FALSE
##' )
##' standardized_prediction$predictions
##'
##' # Aggregate joint predictions within two example areas
##' areas <- sf::st_sf(
##'   area = 1:2,
##'   geometry = sf::st_make_grid(example_data, n = c(2, 1))
##' )
##' areal_prediction <- predict(
##'   fit,
##'   newdata = example_data,
##'   boundaries = areas,
##'   aggregate = "mean",
##'   messages = FALSE
##' )
##' areal_prediction$boundary_summary
##'
##' @seealso [glgpm()], [fitted()], [set_control_mcmc()]
##' @method predict RiskMap
##' @export
predict.RiskMap <- function(
    object,
    newdata = NULL,
    boundaries = NULL,
    aggregate = NULL,
    type = c("auto", "marginal", "joint"),
    include = NULL,
    transform = NULL,
    standardize = NULL,
    summaries = c("mean", "lower", "upper"),
    thresholds = NULL,
    level = 0.95,
    weights = NULL,
    control_mcmc = NULL,
    parameter_uncertainty = FALSE,
    parameter_draws = 50L,
    seed = NULL,
    keep_samples = FALSE,
    messages = getOption("RiskMap.messages", TRUE),
    ...) {

  if (is.null(include)) include <- character()
  check_logical(parameter_uncertainty)
  check_positive_integer(parameter_draws, "parameter_draws")
  parameter_draws <- if (parameter_uncertainty) {
    as.integer(parameter_draws)
  } else {
    NULL
  }
  aggregate <- normalise_prediction_operations(
    aggregate, c("mean", "median", "sum", "min", "max"),
    argument_name = "aggregate",
    allow_null = TRUE
  )
  summaries <- normalise_prediction_operations(
    summaries,
    c("mean", "median", "sd", "cv", "lower", "upper", "interval_width",
      "exceedance"),
    argument_name = "summaries"
  )
  transform <- normalise_prediction_transforms(transform)
  standardize <- validate_prediction_standardize(object, standardize)

  type <- validate_prediction_options(
    object = object,
    type = type,
    include = include,
    boundaries = boundaries,
    aggregate = aggregate,
    summaries = summaries,
    thresholds = thresholds,
    level = level,
    weights = weights,
    control_mcmc = control_mcmc,
    seed = seed,
    keep_samples = keep_samples,
    messages = messages
  )
  if (!is.null(seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(seed)
  }
  observed_locations <- is.null(newdata)
  if (observed_locations && "nugget" %in% include) {
    abort_input(c("{.arg include} can't contain {.val nugget} at the observed locations.", "i" = "There the effect of each location is drawn given the data, and it already includes the nugget."), arg = "include")
  }
  prediction_input <- validate_prediction_data(
    object, newdata, include, standardize
  )
  if (!is.null(boundaries)) {
    if (observed_locations) {
      abort_input("{.arg newdata} must be given to aggregate predictions over {.arg boundaries}.", arg = "newdata")
    }
    return(predict_over_boundaries(
      object = object,
      prediction_input = prediction_input,
      boundaries = boundaries,
      aggregate = aggregate,
      include = include,
      transform = transform,
      standardize = standardize,
      summaries = summaries,
      thresholds = thresholds,
      level = level,
      weights = weights,
      control_mcmc = control_mcmc,
      parameter_draws = parameter_draws,
      seed = seed,
      keep_samples = keep_samples,
      messages = messages
    ))
  }
  prediction_data <- prediction_input$model_data
  if (!nrow(prediction_data)) {
    output_data <- format_missing_location_predictions(
      prediction_input, include, transform, summaries, thresholds
    )
    out <- list(
      predictions = output_data,
      boundaries = NULL,
      boundary_summary = NULL,
      samples = NULL,
      include = include,
      transform = names(transform),
      standardized = !is.null(standardize),
      aggregate = NULL,
      summaries = names(summaries),
      thresholds = thresholds,
      level = level,
      type = type,
      family = object$family,
      response_name = prediction_response_name(object),
      observed_locations = observed_locations,
      parameter_draws = parameter_draws,
      coverage = NULL,
      weights = NULL
    )
    class(out) <- "RiskMap_prediction"
    return(out)
  }
  sample_components <- prediction_sample_components(
    object, prediction_data, type, include, control_mcmc, messages,
    observed_locations = observed_locations, transform = transform,
    standardize = standardize, parameter_draws = parameter_draws,
    mcmc_seed = seed
  )
  output_data <- format_location_predictions(
    prediction_input, sample_components, include, transform, summaries,
    thresholds, level
  )
  valid_rows <- prediction_input$valid_rows
  samples <- if (keep_samples) {
    retained <- sample_components[c("spatial", "random_effects", "nugget")]
    retained$fixed <- sample_components$fixed_samples
    retained <- retained[!vapply(retained, is.null, logical(1))]
    lapply(retained, restore_prediction_rows,
           valid_rows = valid_rows, n_rows = nrow(output_data))
  } else {
    NULL
  }
  out <- list(
    predictions = output_data,
    boundaries = NULL,
    boundary_summary = NULL,
    samples = samples,
    include = include,
    transform = names(transform),
    standardized = !is.null(standardize),
    aggregate = NULL,
    summaries = names(summaries),
    thresholds = thresholds,
    level = level,
    type = type,
    family = object$family,
    response_name = prediction_response_name(object),
    observed_locations = observed_locations,
    parameter_draws = parameter_draws,
    coverage = NULL,
    weights = NULL
  )
  class(out) <- "RiskMap_prediction"
  out
}

#' Validate the unified prediction options
#' @param call The call named in input errors, so that they name
#'   `predict.RiskMap()`.
#' @noRd
validate_prediction_options <- function(object, type, include,
                                        boundaries, aggregate, summaries,
                                        thresholds, level, weights,
                                        control_mcmc, seed, keep_samples,
                                        messages, call = rlang::caller_env()) {
  type <- match_choice(type, c("auto", "marginal", "joint"), arg = "type",
                       call = call)
  if (!is.null(boundaries)) {
    check_data(boundaries, geometry = "polygon", type = "sf")
    if (type == "marginal") {
      abort_input("{.arg type} must be {.val joint} when {.arg boundaries} are given.", arg = "type", call = call)
    }
    type <- "joint"
    if (is.null(aggregate)) {
      choices <- c("mean", "median", "sum", "min", "max")
      abort_input("{.arg aggregate} must be given with {.arg boundaries}, as one or more of {.or {.val {choices}}}.", arg = "aggregate", call = call)
    }
  } else {
    if (!is.null(aggregate)) {
      abort_input("{.arg aggregate} is used only with {.arg boundaries}.", arg = "aggregate", call = call)
    }
    if (type == "auto") type <- "marginal"
  }
  allowed_include <- c("random_effects", "nugget")
  if (!is.character(include) || anyNA(include) ||
      any(!include %in% allowed_include) || anyDuplicated(include)) {
    abort_input("{.arg include} must contain distinct values among {.val {allowed_include}}.", arg = "include", call = call)
  }
  if (prediction_operation_requested(summaries, "exceedance")) {
    if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds) ||
        any(!is.finite(thresholds)) || anyDuplicated(thresholds)) {
      abort_input("{.arg thresholds} must be distinct finite numbers when {.arg summaries} includes {.val exceedance}.", arg = "thresholds", call = call)
    }
  } else if (!is.null(thresholds)) {
    abort_input("{.arg thresholds} is used only when {.arg summaries} includes {.val exceedance}.", arg = "thresholds", call = call)
  }
  check_range(level, min = 0, max = 1, allow_equal = FALSE, arg = "level",
              call = call)
  if (!is.null(control_mcmc)) {
    check_control(control_mcmc, "set_control_mcmc", call = call)
  }
  check_positive_integer(seed, "seed", allow_null = TRUE,
                         allow_zero = FALSE, call = call)
  if (!is.null(weights) &&
      (!is.character(weights) || length(weights) != 1L || is.na(weights))) {
    supplied <- describe_value(weights)
    abort_input("{.arg weights} must be {.code NULL} or the name of one column of {.arg newdata}, not {supplied}.", arg = "weights", call = call)
  }
  if (!is.null(weights) && is.null(boundaries)) {
    abort_input("{.arg weights} is used only with {.arg boundaries}.", arg = "weights", call = call)
  }
  check_logical(keep_samples, call = call)
  check_logical(messages, call = call)
  type
}

#' Normalise built-in and custom prediction operations
#' @noRd
normalise_prediction_operations <- function(x, allowed, argument_name,
                                            allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    abort_input("{.arg {argument_name}} can't be {.code NULL}.", arg = argument_name)
  }
  if (is.character(x)) {
    if (!length(x) || anyNA(x) || any(!x %in% allowed) || anyDuplicated(x)) {
      abort_input("{.arg {argument_name}} must contain distinct values among {.val {allowed}}, or named functions.", arg = argument_name)
    }
    out <- as.list(x)
    names(out) <- x
    return(out)
  }
  if (!is.list(x) || !length(x)) {
    abort_input("{.arg {argument_name}} must be a character vector, or a named list of the names of built-in operations and functions.", arg = argument_name)
  }
  operation_names <- names(x)
  if (is.null(operation_names) || anyNA(operation_names) ||
      any(!nzchar(operation_names)) || anyDuplicated(operation_names)) {
    abort_input("Every entry of {.arg {argument_name}} must have a distinct, non-empty name.", arg = argument_name)
  }
  valid <- vapply(x, function(operation) {
    is.function(operation) ||
      (is.character(operation) && length(operation) == 1L &&
         !is.na(operation) && operation %in% allowed)
  }, logical(1))
  if (!all(valid)) {
    abort_input("Each entry of {.arg {argument_name}} must be a function or one of {.or {.val {allowed}}}.", arg = argument_name)
  }
  x
}

#' Validate derived predictive targets
#' @noRd
normalise_prediction_transforms <- function(transform) {
  if (is.null(transform)) return(list())
  valid_names <- !is.null(names(transform)) && !anyNA(names(transform)) &&
    all(nzchar(names(transform))) && !anyDuplicated(names(transform))
  if (!is.list(transform) || !length(transform) || !valid_names ||
      any(!vapply(transform, is.function, logical(1)))) {
    abort_input("{.arg transform} must be {.code NULL} or a list of functions with distinct, non-empty names.", arg = "transform")
  }
  reserved <- c("fixed", "offset", "spatial", "random_effects", "nugget",
                "link", "response")
  conflicts <- intersect(names(transform), reserved)
  if (length(conflicts)) {
    abort_input(c("The names in {.arg transform} can't be those of built-in prediction quantities.", "x" = "{.val {conflicts}} {cli::qty(length(conflicts))}{?is a built-in quantity/are built-in quantities}."), arg = "transform")
  }
  transform
}

#' Validate a reference population for direct standardization
#' @noRd
validate_prediction_standardize <- function(object, standardize) {
  if (is.null(standardize)) return(NULL)
  if (!is.data.frame(standardize) || inherits(standardize, "sf") ||
      !nrow(standardize)) {
    abort_input("{.arg standardize} must be a non-empty data frame without geometry, not {.obj_type_friendly {standardize}}.", arg = "standardize")
  }
  if (anyDuplicated(names(standardize)) || anyNA(names(standardize)) ||
      any(!nzchar(names(standardize)))) {
    abort_input("{.arg standardize} must have distinct, non-empty column names.", arg = "standardize")
  }
  reference_names <- setdiff(names(standardize), ".weight")
  response <- as.character(object$formula[[2L]])
  inter_f <- interpret_formula(object$formula)
  offset_names <- inter_f$offset %||% character()
  fixed_names <- setdiff(
    get_formula_terms(object$formula),
    c(response, names(object$re), offset_names)
  )
  if (!length(reference_names) || any(!reference_names %in% fixed_names)) {
    invalid <- setdiff(reference_names, fixed_names)
    problem <- "{.arg standardize} must have one or more columns, each a covariate of the model."
    if (length(invalid)) {
      problem <- c(problem, "x" = "{.val {invalid}} {cli::qty(length(invalid))}{?is not a covariate/are not covariates}.")
    }
    abort_input(problem, arg = "standardize")
  }
  if (!all(complete.cases(standardize[reference_names]))) {
    abort_input("{.arg standardize} has missing values.", arg = "standardize")
  }
  weights <- if (".weight" %in% names(standardize)) {
    standardize$.weight
  } else {
    rep(1, nrow(standardize))
  }
  if (!is.numeric(weights) || anyNA(weights) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    abort_input("The column {.code .weight} of {.arg standardize} must hold finite, non-negative values with a positive total.", arg = "standardize")
  }
  list(
    data = standardize,
    weights = weights / sum(weights),
    predictors = reference_names
  )
}

#' Remove named prediction operations
#' @noRd
drop_prediction_operations <- function(operations, built_ins) {
  keep <- !vapply(operations, function(operation) {
    is.character(operation) && operation %in% built_ins
  }, logical(1))
  operations[keep]
}

#' Check whether a built-in prediction operation was requested
#' @noRd
prediction_operation_requested <- function(operations, built_in) {
  any(vapply(operations, function(operation) {
    is.character(operation) && identical(operation, built_in)
  }, logical(1)))
}

#' Validate and align prediction locations
#' @noRd
validate_prediction_data <- function(object, newdata, include,
                                     standardize = NULL) {
  prediction_data <- if (is.null(newdata)) object$data else newdata
  data_arg <- if (is.null(newdata)) "object" else "newdata"
  check_data(prediction_data, geometry = "point", type = "sf", arg = data_arg)
  # The compiled distances need finite planar coordinates.
  coordinates <- st_coordinates(prediction_data)
  if (ncol(coordinates) != 2L) {
    abort_input("The points of {.arg {data_arg}} must have two coordinates, without Z or M values.", arg = data_arg)
  }
  if (any(!is.finite(coordinates))) {
    abort_input("The points of {.arg {data_arg}} must have finite coordinates, and none may be empty.", arg = data_arg)
  }
  output_data <- prediction_data
  if (!isTRUE(st_crs(prediction_data) == st_crs(object$data))) {
    riskmap_inform("{.arg newdata} is projected to the coordinate reference system of the fit.")
    prediction_data <- st_transform(prediction_data, st_crs(object$data))
  }
  inter_f <- interpret_formula(object$formula)
  response <- as.character(object$formula[[2]])
  required <- setdiff(get_formula_terms(object$formula), response)
  if (!"random_effects" %in% include) {
    required <- setdiff(required, names(object$re))
  }
  standardized_names <- if (is.null(standardize)) character() else
    setdiff(names(standardize$data), ".weight")
  required <- setdiff(required, standardized_names)
  missing_columns <- setdiff(required, names(prediction_data))
  if (length(missing_columns) > 0L) {
    abort_input("{.arg newdata} must have the column{?s} {.val {missing_columns}}, which the model uses.", arg = "newdata")
  }
  for (name in standardized_names) {
    prediction_data[[name]] <- standardize$data[[name]][1L]
  }
  if ("random_effects" %in% include) {
    for (name in names(object$re)) {
      supplied <- prediction_data[[name]]
      unknown <- unique(supplied[!is.na(supplied) &
                                   !supplied %in% object$re[[name]]])
      if (length(unknown) > 0L) {
        abort_input(c("{.arg newdata} has levels of the random effect {.val {name}} that the fit doesn't know.", "x" = "The unknown {cli::qty(length(unknown))}level{?s} {?is/are} {.val {unknown}}."), arg = "newdata")
      }
    }
  }
  valid_rows <- if (length(required) > 0L) {
    complete.cases(st_drop_geometry(prediction_data)[required])
  } else {
    rep(TRUE, nrow(prediction_data))
  }
  if (any(!valid_rows)) {
    riskmap_warn("{sum(!valid_rows)} prediction row{?s} {?has/have} missing values, so {?its/their} predictions are {.code NA}.")
  }
  list(
    model_data = prediction_data[valid_rows, , drop = FALSE],
    output_data = output_data,
    valid_rows = valid_rows
  )
}

#' Restore prediction matrices to the original row alignment
#' @noRd
restore_prediction_rows <- function(values, valid_rows, n_rows) {
  out <- matrix(NA_real_, nrow = n_rows, ncol = ncol(values))
  out[valid_rows, ] <- values
  out
}

#' Assign prediction support to polygon boundaries
#' @noRd
assign_prediction_locations <- function(model_data, boundaries) {
  metadata <- attr(model_data, "grid_metadata", exact = TRUE)
  if (inherits(model_data, "RiskMap_grid") && !is.null(metadata)) {
    return(assign_prediction_grid_cells(model_data, boundaries, metadata))
  }

  membership <- st_within(model_data, boundaries)
  membership_count <- lengths(membership)
  if (any(membership_count > 1L)) {
    abort_input("The polygons of {.arg boundaries} must not overlap, but {sum(membership_count > 1L)} prediction location{?s} {?lies/lie} in more than one.", arg = "boundaries")
  }
  boundary_candidates <- which(membership_count == 0L)
  if (length(boundary_candidates)) {
    inclusive <- st_intersects(model_data[boundary_candidates, , drop = FALSE],
                               boundaries)
    on_boundary <- lengths(inclusive) > 0L
    membership[boundary_candidates[on_boundary]] <- lapply(
      inclusive[on_boundary], function(index) index[[1L]]
    )
    membership_count <- lengths(membership)
  }
  assigned <- membership_count == 1L
  prediction_row <- which(assigned)
  list(
    data = model_data[prediction_row, , drop = FALSE],
    prediction_row = prediction_row,
    area_id = vapply(membership[assigned], `[[`, integer(1), 1L),
    coverage = rep(1, sum(assigned)),
    outside = which(!assigned),
    method = "point"
  )
}

#' Calculate exact polygon coverage for cells created by create_grid
#' @noRd
assign_prediction_grid_cells <- function(model_data, boundaries, metadata) {
  if (!identical(metadata$cell_shape, "square") ||
      !is.numeric(metadata$cellsize) || length(metadata$cellsize) != 1L ||
      !is.finite(metadata$cellsize) || metadata$cellsize <= 0 ||
      is.null(metadata$grid_crs)) {
    abort_input(c("{.arg newdata} has invalid grid metadata.", "i" = "Make the grid with {.fn create_grid}, or drop its class {.cls RiskMap_grid}."), arg = "newdata")
  }
  if (!nrow(model_data)) {
    return(list(
      data = model_data,
      prediction_row = integer(),
      area_id = integer(),
      coverage = numeric(),
      outside = integer(),
      method = "exact"
    ))
  }
  grid_crs <- st_crs(metadata$grid_crs)
  grid_boundaries <- st_transform(boundaries, grid_crs)
  overlaps <- st_overlaps(grid_boundaries)
  duplicates <- st_equals(grid_boundaries)
  if (any(lengths(overlaps) > 0L) || any(lengths(duplicates) > 1L)) {
    abort_input("The polygons of {.arg boundaries} must not overlap.", arg = "boundaries")
  }

  grid_points <- st_transform(model_data, grid_crs)
  coordinates <- st_coordinates(grid_points)
  half_cell <- metadata$cellsize / 2
  x_limits <- range(coordinates[, 1]) + c(-half_cell, half_cell)
  y_limits <- range(coordinates[, 2]) + c(-half_cell, half_cell)
  n_columns <- round(diff(x_limits) / metadata$cellsize)
  n_rows <- round(diff(y_limits) / metadata$cellsize)
  # Reconstruct the original lattice, leaving cells absent from newdata empty.
  cell_raster <- terra::rast(
    nrows = n_rows, ncols = n_columns,
    xmin = x_limits[1], xmax = x_limits[2],
    ymin = y_limits[1], ymax = y_limits[2],
    crs = grid_crs$wkt
  )
  raster_cells <- terra::cellFromXY(cell_raster, coordinates)
  if (anyDuplicated(raster_cells)) {
    abort_input("{.arg newdata} has the same grid cell more than once.", arg = "newdata")
  }
  raster_values <- rep(NA_real_, terra::ncell(cell_raster))
  raster_values[raster_cells] <- seq_len(nrow(model_data))
  terra::values(cell_raster) <- raster_values
  names(cell_raster) <- "prediction_row"

  overlap <- terra::extract(
    cell_raster,
    terra::vect(grid_boundaries),
    cells = TRUE,
    exact = TRUE
  )
  keep <- !is.na(overlap$prediction_row) & overlap$fraction > 0
  overlap <- overlap[keep, , drop = FALSE]
  prediction_row <- as.integer(overlap$prediction_row)
  list(
    data = model_data[prediction_row, , drop = FALSE],
    prediction_row = prediction_row,
    area_id = as.integer(overlap$ID),
    coverage = as.numeric(overlap$fraction),
    outside = setdiff(seq_len(nrow(model_data)), unique(prediction_row)),
    method = "exact"
  )
}

#' Construct independent joint-prediction blocks for polygon aggregation
#' @noRd
prediction_boundary_groups <- function(assignment, n_locations) {
  prediction_row <- assignment$prediction_row
  area_id <- assignment$area_id
  areas <- sort(unique(area_id))
  parent <- seq_len(max(c(areas, 0L)))
  find_root <- function(x) {
    while (parent[[x]] != x) x <- parent[[x]]
    x
  }
  if (length(areas)) {
    shared <- split(area_id, prediction_row)
    for (members in shared) {
      members <- unique(members)
      if (length(members) > 1L) {
        root <- find_root(members[[1L]])
        for (member in members[-1L]) parent[[find_root(member)]] <- root
      }
    }
  }
  groups <- rep(NA_integer_, n_locations)
  if (length(prediction_row)) {
    location_areas <- split(area_id, prediction_row)
    rows <- as.integer(names(location_areas))
    groups[rows] <- vapply(
      location_areas, function(x) find_root(x[[1L]]), integer(1)
    )
  }
  outside <- which(is.na(groups))
  if (length(outside)) {
    groups[outside] <- max(c(groups, 0L), na.rm = TRUE) + seq_along(outside)
  }
  match(groups, unique(groups))
}

#' Predict and aggregate jointly within polygon boundaries
#' @noRd
predict_over_boundaries <- function(object, prediction_input, boundaries,
                                    aggregate, include, transform, standardize,
                                    summaries,
                                    thresholds,
                                    level, weights,
                                    control_mcmc, seed, keep_samples,
                                    messages, parameter_draws = NULL) {
  model_data <- prediction_input$model_data
  output_boundaries <- boundaries
  if (!isTRUE(st_crs(boundaries) == st_crs(object$data))) {
    riskmap_inform("{.arg boundaries} are projected to the coordinate reference system of the fit.")
    boundaries <- st_transform(boundaries, st_crs(object$data))
  }

  assignment <- assign_prediction_locations(model_data, boundaries)
  assigned_data <- assignment$data
  area_id <- assignment$area_id
  coverage <- assignment$coverage
  if (length(assignment$outside)) {
    unit <- if (assignment$method == "exact") "cell" else "location"
    n_outside <- length(assignment$outside)
    riskmap_warn("{n_outside} prediction {unit}{cli::qty(n_outside)}{?s} {?lies/lie} outside every polygon, so {?it doesn't/they don't} count in the summaries of the areas.")
  }

  if (is.null(weights)) {
    weight_values <- coverage
  } else {
    if (!weights %in% names(assigned_data)) {
      abort_input("{.arg weights} must name a column of {.arg newdata}, and {.val {weights}} doesn't.", arg = "weights")
    }
    supplied_weights <- assigned_data[[weights]]
    if (!is.numeric(supplied_weights) || anyNA(supplied_weights) ||
        any(!is.finite(supplied_weights)) || any(supplied_weights < 0)) {
      abort_input("The column {.val {weights}} named by {.arg weights} must hold finite, non-negative values.", arg = "weights")
    }
    # Allocate each cell's supplied weight in proportion to polygon coverage.
    weight_values <- supplied_weights * coverage
  }

  area_indices <- split(seq_len(nrow(assigned_data)), area_id)
  nonempty_areas <- as.integer(names(area_indices))
  single_location <- lengths(area_indices) == 1L
  if (any(single_location)) {
    unit <- if (assignment$method == "exact") "cell" else "location"
    n_single <- sum(single_location)
    riskmap_warn("{n_single} polygon{?s} {?contains/contain} only one prediction {unit}, so {cli::qty(n_single)}{?its/their} summaries rest on a single value.")
  }
  zero_weight <- vapply(area_indices, function(index) {
    sum(weight_values[index]) <= 0
  }, logical(1))
  weighted_operations <- vapply(aggregate, function(operation) {
    is.character(operation) && operation %in% c("mean", "median", "sum")
  }, logical(1))
  if (!is.null(weights) && any(weighted_operations) && any(zero_weight)) {
    riskmap_warn("{sum(zero_weight)} polygon{?s} {?has/have} a total weight of zero, so {?its/their} weighted means, medians and sums are {.code NA}.")
  }

  if (!nrow(model_data)) {
    out <- list(
      predictions = format_missing_location_predictions(
        prediction_input, include, transform, summaries, thresholds
      ),
      boundaries = output_boundaries,
      boundary_summary = missing_boundary_summaries(
        nrow(output_boundaries), names(aggregate), include, names(transform),
        summaries, thresholds
      ),
      samples = NULL, include = include, transform = names(transform),
      standardized = !is.null(standardize), aggregate = names(aggregate),
      summaries = names(summaries), thresholds = thresholds, level = level,
      type = "joint", family = object$family,
      response_name = prediction_response_name(object),
      observed_locations = FALSE,
      parameter_draws = parameter_draws,
      coverage = assignment$method,
      weights = weights
    )
    class(out) <- "RiskMap_prediction"
    return(out)
  }

  prediction_groups <- prediction_boundary_groups(
    assignment, nrow(model_data)
  )
  sample_components <- prediction_sample_components(
    object, model_data, "joint", include, control_mcmc, messages,
    groups = prediction_groups, transform = transform,
    standardize = standardize, parameter_draws = parameter_draws,
    mcmc_seed = seed
  )
  predictions <- format_location_predictions(
    prediction_input, sample_components, include, transform, summaries,
    thresholds, level
  )
  components <- list(
    fixed = sample_components$fixed_samples %||%
      matrix(sample_components$fixed, ncol = 1L),
    offset = matrix(sample_components$offset, ncol = 1L),
    spatial = sample_components$spatial,
    random_effects = sample_components$random_effects,
    nugget = sample_components$nugget,
    link = sample_components$link,
    response = sample_components$response
  )
  components <- c(components, sample_components$transformed)
  components <- components[!vapply(components, is.null, logical(1))]

  n_areas <- nrow(output_boundaries)
  retained <- NULL
  summary_rows <- list()
  row_number <- 0L
  for (group in seq_along(area_indices)) {
    index <- area_indices[[group]]
    area <- nonempty_areas[[group]]
    group_weights <- weight_values[index]
    prediction_rows <- assignment$prediction_row[index]
    group_components <- lapply(
      components, function(values) values[prediction_rows, , drop = FALSE]
    )

    for (operation_name in names(aggregate)) {
      operation <- aggregate[[operation_name]]
      aggregated <- lapply(
        group_components, aggregate_prediction_values,
        operation = operation, weights = group_weights
      )
      for (quantity in names(aggregated)) {
        values <- aggregated[[quantity]]
        row_number <- row_number + 1L
        summary_rows[[row_number]] <- if (quantity %in% c("fixed", "offset")) {
          data.frame(
            boundary_id = area,
            aggregate = operation_name,
            quantity = quantity,
            summary = "estimate",
            threshold = NA_real_,
            # Constant across draws for plug-in prediction; the average over
            # parameter draws when parameter uncertainty is propagated.
            value = mean(values),
            stringsAsFactors = FALSE
          )
        } else {
          tidy_prediction_summaries(
            values = values,
            boundary_id = area,
            aggregate = operation_name,
            quantity = quantity,
            summaries = if (quantity %in%
                              c("response", names(transform))) summaries else
              drop_prediction_operations(summaries, c("cv", "exceedance")),
            thresholds = if (quantity %in%
                               c("response", names(transform))) thresholds else NULL,
            level = level
          )
        }
      }
      if (keep_samples) {
        if (is.null(retained[[operation_name]])) {
          retained[[operation_name]] <- list()
        }
        for (quantity in setdiff(names(aggregated), c("fixed", "offset"))) {
          if (is.null(retained[[operation_name]][[quantity]])) {
            retained[[operation_name]][[quantity]] <- matrix(
              NA_real_, nrow = n_areas, ncol = ncol(aggregated[[quantity]])
            )
          }
          retained[[operation_name]][[quantity]][area, ] <-
            aggregated[[quantity]]
        }
      }
    }
  }
  missing_areas <- setdiff(seq_len(n_areas), nonempty_areas)
  if (length(missing_areas)) {
    quantities <- c(
      "fixed", "offset", "spatial",
      if ("random_effects" %in% include) "random_effects",
      if ("nugget" %in% include) "nugget",
      "link", "response"
    )
    quantities <- c(quantities, names(transform))
    for (area in missing_areas) {
      for (operation_name in names(aggregate)) {
        for (quantity in quantities) {
          row_number <- row_number + 1L
          if (quantity %in% c("fixed", "offset")) {
            summary_rows[[row_number]] <- data.frame(
              boundary_id = area, aggregate = operation_name,
              quantity = quantity,
              summary = "estimate", threshold = NA_real_, value = NA_real_,
              stringsAsFactors = FALSE
            )
          } else {
            requested <- if (quantity %in%
                               c("response", names(transform))) summaries else
              drop_prediction_operations(summaries, c("cv", "exceedance"))
            summary_rows[[row_number]] <- missing_prediction_summaries(
              boundary_id = area,
              aggregate = operation_name,
              quantity = quantity,
              summaries = requested,
              thresholds = if (quantity %in%
                                 c("response", names(transform))) thresholds else NULL
            )
          }
        }
      }
    }
  }

  out <- list(
    predictions = predictions,
    boundaries = output_boundaries,
    boundary_summary = if (length(summary_rows))
      do.call(rbind, summary_rows) else empty_prediction_summary(),
    samples = retained,
    include = include,
    transform = names(transform),
    standardized = !is.null(standardize),
    aggregate = names(aggregate),
    summaries = names(summaries),
    thresholds = thresholds,
    level = level,
    type = "joint",
    family = object$family,
    response_name = prediction_response_name(object),
    observed_locations = FALSE,
    parameter_draws = parameter_draws,
    coverage = assignment$method,
    weights = weights
  )
  class(out) <- "RiskMap_prediction"
  out
}

#' Run the conditional-prediction engine from aligned sf data
#' @noRd
prepare_prediction <- function(object, prediction_data, prediction, include,
                               control_mcmc, messages,
                               observed_locations = FALSE, groups = NULL,
                               mcmc_seed = NULL) {
  draw_prediction_samples(
    object = object,
    prediction_data = prediction_data,
    prediction = prediction,
    include = include,
    control_mcmc = control_mcmc,
    groups = groups,
    observed_locations = observed_locations,
    mcmc_seed = mcmc_seed,
    messages = messages
  )
}

#' Construct link and response samples from prediction components
#' @noRd
build_prediction_components <- function(object, prediction_setup, include,
                                        transform = list(), standardize = NULL,
                                        prediction_data = NULL) {
  spatial <- prediction_setup$S_samples
  if (!is.matrix(spatial)) {
    spatial <- matrix(spatial, nrow = nrow(prediction_setup$grid_pred))
  }
  n_locations <- nrow(spatial)
  n_samples <- ncol(spatial)
  fixed <- rep(prediction_setup$mu_pred, length.out = n_locations)
  offset <- rep(prediction_setup$cov_offset, length.out = n_locations)
  random_effects <- prediction_random_effect_samples(
    prediction_setup, n_locations, n_samples,
    include = "random_effects" %in% include
  )
  nugget <- prediction_nugget_samples(
    object, prediction_setup, n_locations, n_samples,
    include = "nugget" %in% include
  )
  link <- sweep(spatial, 1L, fixed + offset, "+")
  if (!is.null(random_effects)) link <- link + random_effects
  if (!is.null(nugget)) link <- link + nugget
  response <- if (object$family == "gaussian") {
    link
  } else {
    matrix(object$link_function$inv(as.numeric(link)), nrow = n_locations)
  }
  transformed <- lapply(names(transform), function(name) {
    values <- transform[[name]](link)
    if (!is.numeric(values) || !is.matrix(values) ||
        !identical(dim(values), dim(link))) {
      abort_input("The function {.val {name}} in {.arg transform} must return a numeric matrix of the same dimensions as its input.", arg = "transform")
    }
    values
  })
  names(transformed) <- names(transform)
  if (!is.null(standardize)) {
    standardized <- standardize_prediction_targets(
      object = object,
      prediction_data = prediction_data,
      standardize = standardize,
      spatial = spatial,
      offset = offset,
      random_effects = random_effects,
      nugget = nugget,
      transform = transform,
      beta = prediction_setup$par_hat$beta
    )
    fixed <- standardized$fixed
    link <- standardized$link
    response <- standardized$response
    transformed <- standardized$transformed
  }
  list(fixed = fixed, offset = offset, spatial = spatial,
       random_effects = random_effects, nugget = nugget,
       link = link, response = response, transformed = transformed)
}

#' Prediction components, optionally averaged over parameter uncertainty
#'
#' With `parameter_draws = NULL` this is plug-in prediction at the estimate.
#' Otherwise the parameters are drawn from the approximate sampling
#' distribution of the estimate and, for each draw, conditional samples of the
#' latent effects are generated at that parameter value and passed through the
#' ordinary prediction machinery. The draws are pooled, so every summary,
#' aggregation and transformation reflects both sources of uncertainty.
#' `mcmc_seed` seeds a fresh conditional MCMC run of plug-in prediction, as
#' in `draw_prediction_samples()`.
#' @noRd
prediction_sample_components <- function(object, prediction_data, prediction,
                                         include, control_mcmc, messages,
                                         observed_locations = FALSE,
                                         groups = NULL, transform = list(),
                                         standardize = NULL,
                                         parameter_draws = NULL,
                                         mcmc_seed = NULL) {
  components_at <- function(fit, control) {
    setup <- prepare_prediction(
      fit, prediction_data, prediction, include, control, messages = FALSE,
      observed_locations = observed_locations, groups = groups
    )
    build_prediction_components(
      fit, setup, include, transform, standardize = standardize,
      prediction_data = prediction_data
    )
  }
  if (is.null(parameter_draws)) {
    setup <- prepare_prediction(
      object, prediction_data, prediction, include, control_mcmc, messages,
      observed_locations = observed_locations, groups = groups,
      mcmc_seed = mcmc_seed
    )
    return(build_prediction_components(
      object, setup, include, transform, standardize = standardize,
      prediction_data = prediction_data
    ))
  }

  n_total <- prediction_sample_size(object, control_mcmc)
  per_draw <- max(1L, ceiling(n_total / parameter_draws))
  parameters <- draw_parameter_values(object, parameter_draws)
  if (messages) {
    riskmap_inform("Propagating parameter uncertainty over {parameter_draws} parameter draws.")
  }
  sampler <- if (object$family == "gaussian") NULL else
    parameter_conditioning_sampler(object)
  pieces <- vector("list", parameter_draws)
  acceptance <- rep(NA_real_, parameter_draws)
  # A failure at one parameter value leaves that draw out. An input error
  # would fail at every value, so it reaches the user as it is.
  first_failure <- NULL
  leave_out <- function(error) {
    if (inherits(error, "riskmap_error_input")) riskmap_rethrow(error)
    if (is.null(first_failure)) first_failure <<- error
    NULL
  }
  for (k in seq_len(parameter_draws)) {
    fit <- object
    fit$estimate <- utils::relist(parameters[k, ], skeleton = object$estimate)
    if (is.null(sampler)) {
      control <- set_control_mcmc(n_sim = per_draw, linear_model = TRUE)
    } else {
      draws <- tryCatch(sampler(fit$estimate, per_draw), error = leave_out)
      if (is.null(draws)) next
      fit$S_samples <- draws$samples
      acceptance[k] <- draws$acceptance
      control <- NULL
    }
    pieces[[k]] <- tryCatch(components_at(fit, control), error = leave_out)
  }
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) < parameter_draws / 2) {
    riskmap_abort("Prediction failed at most parameter draws, so parameter uncertainty can't be propagated.", class = "riskmap_error_fit", parent = first_failure)
  }
  if (length(pieces) < parameter_draws) {
    failed <- parameter_draws - length(pieces)
    riskmap_warn("{failed} of {parameter_draws} parameter draws are left out, because prediction failed at those values.", class = "riskmap_warning_fit")
  }
  if (any(acceptance < 0.1, na.rm = TRUE)) {
    riskmap_warn("The Langevin sampler accepted few proposals at some parameter draws, so the predictive summaries may be noisy.", class = "riskmap_warning_fit")
  }
  pool_prediction_components(pieces)
}

#' Pool prediction components drawn at different parameter values
#'
#' The fixed-effect contribution differs between parameter draws, so it is
#' returned both as its average (`fixed`) and draw by draw (`fixed_samples`).
#' @noRd
pool_prediction_components <- function(pieces) {
  bind <- function(name) {
    values <- lapply(pieces, `[[`, name)
    if (is.null(values[[1L]])) NULL else do.call(cbind, values)
  }
  fixed_samples <- do.call(cbind, lapply(pieces, function(piece) {
    matrix(piece$fixed, nrow = length(piece$fixed), ncol = ncol(piece$spatial))
  }))
  transformed <- lapply(names(pieces[[1L]]$transformed), function(name) {
    do.call(cbind, lapply(pieces, function(piece) piece$transformed[[name]]))
  })
  names(transformed) <- names(pieces[[1L]]$transformed)
  list(fixed = rowMeans(fixed_samples), fixed_samples = fixed_samples,
       offset = pieces[[1L]]$offset, spatial = bind("spatial"),
       random_effects = bind("random_effects"), nugget = bind("nugget"),
       link = bind("link"), response = bind("response"),
       transformed = transformed)
}

#' Number of predictive draws used by plug-in prediction
#' @noRd
prediction_sample_size <- function(object, control_mcmc) {
  if (!is.null(control_mcmc)) {
    if (isTRUE(control_mcmc$linear_model)) return(control_mcmc$n_sim)
    return((control_mcmc$n_sim - control_mcmc$burnin) / control_mcmc$thin)
  }
  if (object$family != "gaussian" && !is.null(object$S_samples)) {
    return(nrow(object$S_samples))
  }
  # As many draws as plug-in prediction takes by default for every family.
  control_mcmc <- set_control_mcmc()
  (control_mcmc$n_sim - control_mcmc$burnin) / control_mcmc$thin
}

#' Draw parameters from the approximate sampling distribution of the estimate
#'
#' Draws are multivariate normal on the scale on which the model is estimated
#' (log variances and log scale), with the covariance matrix of the estimate.
#' Parameters without a finite positive variance, such as those held at a
#' boundary, keep their estimated value.
#' @noRd
draw_parameter_values <- function(object, n_draws) {
  estimate <- unlist(object$estimate)
  covariance <- object$covariance
  failure <- function() {
    riskmap_abort(c("The covariance matrix of the estimates isn't available or isn't positive definite, so parameter uncertainty can't be propagated.", "i" = "For a binomial or Poisson model, refitting with {.code control_mcml = set_control_mcml(pc_prior = TRUE)} usually solves this."), class = "riskmap_error_fit")
  }
  if (is.null(covariance)) failure()
  covariance <- as.matrix(covariance)
  if (!identical(dim(covariance), rep(length(estimate), 2L))) failure()
  free <- which(is.finite(diag(covariance)) & diag(covariance) > 0)
  block <- covariance[free, free, drop = FALSE]
  root <- if (length(free) && all(is.finite(block))) {
    tryCatch(chol((block + t(block)) / 2), error = function(error) NULL)
  }
  if (is.null(root)) failure()
  draws <- matrix(estimate, nrow = n_draws, ncol = length(estimate),
                  byrow = TRUE, dimnames = list(NULL, names(estimate)))
  draws[, free] <- draws[, free, drop = FALSE] +
    matrix(stats::rnorm(n_draws * length(free)), n_draws) %*% root
  draws
}

#' Conditional samples of the latent effects at given parameter values
#'
#' Returns a function drawing samples from the conditional distribution of the
#' latent effects (locations, then random-effect levels) given the data, at
#' any parameter value, with the same Langevin sampler and the same covariance
#' as plug-in prediction. Each short chain runs in the coordinates whitened by
#' the Laplace approximation at that parameter value and starts from a draw of
#' it, so that it is in its stationary spread from the start, and the search
#' for the mode starts from the mode at the previous parameter value. The step
#' is not tuned. It is the one tuned at the estimate, or else the default,
#' which suits the whitened chain as well, so the burn-in only discards the
#' first draws.
#' @noRd
parameter_conditioning_sampler <- function(object, burnin = 100L, thin = 10L) {
  distance_scale <- attr(object, "distance_scale") %||% 1
  fitting_coords <- object$coords / distance_scale
  distances <- pairwise_distances(fitting_coords)
  invlink <- if (identical(object$link_function$name, "canonical")) {
    NULL
  } else {
    object$link_function[c("inv", "d1", "d2")]
  }
  design <- as.matrix(object$D)
  step_size <- utils::tail(attr(object, "mcmc")$tuning_par, 1L)
  if (!length(step_size) || !is.finite(step_size) || step_size <= 0) {
    step_size <- NULL
  }
  last_mode <- NULL
  function(estimate, n_samples) {
    fit <- object
    fit$estimate <- estimate
    parameters <- coef(fit)
    simulation <- laplace_sampling_mcmc(
      y = object$y, units_m = object$units_m,
      mu = as.numeric(design %*% parameters$beta + object$cov_offset),
      Sigma = conditioning_covariance(fit, parameters, distances,
                                      distance_scale),
      ID_coords = object$ID_coords, ID_re = object$ID_re,
      sigma2_re = parameters$sigma2_re, family = object$family,
      control_mcmc = set_control_mcmc(
        n_sim = burnin + n_samples * thin, burnin = burnin, thin = thin,
        h = step_size
      ),
      invlink = invlink, start = last_mode, adapt_step = FALSE
    )
    last_mode <<- simulation$mode
    list(samples = do.call(cbind, simulation$samples),
         acceptance = simulation$acceptance_rate[["sampling"]])
  }
}

#' Covariance of the location effects that prediction conditions on
#'
#' The spatial covariance at the locations plus the nugget, which is at least
#' 1e-10 times the spatial variance, as in `prepare_prediction()`.
#' `distances` are between the locations in the fitting units, which are the
#' model's divided by `distance_scale`.
#' @noRd
conditioning_covariance <- function(object, parameters, distances,
                                    distance_scale) {
  correlation <- matern_correlation(
    distances, phi = parameters$phi / distance_scale, kappa = object$kappa,
    return_sym_matrix = TRUE
  )
  ratio <- fitted_nugget_variance(object, parameters) / parameters$sigma2
  if (ratio == 0) ratio <- 1e-10
  diag(correlation) <- diag(correlation) + ratio
  parameters$sigma2 * correlation
}

#' Average predictive targets over a common reference population
#' @noRd
standardize_prediction_targets <- function(object, prediction_data,
                                           standardize, spatial, offset,
                                           random_effects, nugget, transform,
                                           beta) {
  if (is.null(prediction_data)) {
    riskmap_bug("Prediction data are required for standardisation.")
  }
  plain_data <- as.data.frame(st_drop_geometry(prediction_data))
  n_locations <- nrow(spatial)
  n_samples <- ncol(spatial)
  base_link <- sweep(spatial, 1L, offset, "+")
  if (!is.null(random_effects)) base_link <- base_link + random_effects
  if (!is.null(nugget)) base_link <- base_link + nugget

  fixed <- numeric(n_locations)
  response <- matrix(0, nrow = n_locations, ncol = n_samples)
  transformed <- lapply(transform, function(x) {
    matrix(0, nrow = n_locations, ncol = n_samples)
  })
  # Construct model matrices in bounded batches. This retains vectorised formula
  # evaluation (including transformed terms and interactions) without creating
  # the potentially enormous location-by-reference-population cross-product.
  profiles_per_batch <- max(1L, floor(250000 / n_locations))
  batches <- split(
    seq_len(nrow(standardize$data)),
    ceiling(seq_len(nrow(standardize$data)) / profiles_per_batch)
  )
  for (batch in batches) {
    profile_data <- plain_data[rep(seq_len(n_locations), times = length(batch)),
                               , drop = FALSE]
    for (name in standardize$predictors) {
      profile_data[[name]] <- rep(standardize$data[[name]][batch],
                                  each = n_locations)
    }
    batch_fixed <- matrix(
      prediction_fixed_effect(object, profile_data, beta = beta),
      nrow = n_locations, ncol = length(batch)
    )
    for (column in seq_along(batch)) {
      profile <- batch[[column]]
      profile_fixed <- batch_fixed[, column]
      profile_link <- sweep(base_link, 1L, profile_fixed, "+")
      weight <- standardize$weights[[profile]]
      fixed <- fixed + weight * profile_fixed
      profile_response <- if (object$family == "gaussian") {
        profile_link
      } else {
        matrix(object$link_function$inv(as.numeric(profile_link)),
               nrow = n_locations)
      }
      response <- response + weight * profile_response
      for (name in names(transform)) {
        values <- transform[[name]](profile_link)
        if (!is.numeric(values) || !is.matrix(values) ||
            !identical(dim(values), dim(profile_link))) {
          abort_input("The function {.val {name}} in {.arg transform} must return a numeric matrix of the same dimensions as its input.", arg = "transform")
        }
        transformed[[name]] <- transformed[[name]] + weight * values
      }
    }
  }
  # This is the expected linear predictor over the reference population. For
  # nonlinear links, its inverse is intentionally not used as the standardized
  # response; response averaging occurs above on the scientifically relevant
  # response scale.
  link <- sweep(base_link, 1L, fixed, "+")
  list(fixed = fixed, link = link, response = response,
       transformed = transformed)
}

#' Build the fixed-effect contribution for prediction data
#' @noRd
prediction_fixed_effect <- function(object, predictors, beta = coef(object)$beta) {
  design <- prediction_design_matrix(object, predictors)
  as.numeric(design %*% beta)
}

#' Build a prediction design matrix using the fitted factor structure
#' @noRd
prediction_design_matrix <- function(object, predictors) {
  inter_f <- interpret_formula(object$formula)
  fixed_formula <- update(inter_f$pf, NULL ~ .)
  fitted_frame <- model.frame(fixed_formula, data = object$data,
                              na.action = na.fail)
  fixed_terms <- attr(fitted_frame, "terms")
  model_frame <- model.frame(
    fixed_terms, data = predictors, na.action = na.fail,
    xlev = stats::.getXlevels(fixed_terms, fitted_frame)
  )
  design <- as.matrix(model.matrix(
    fixed_terms, data = model_frame,
    contrasts.arg = attr(object$D, "contrasts")
  ))
  if (!identical(colnames(design), colnames(object$D)) ||
      ncol(design) != ncol(object$D)) {
    abort_input(c("{.arg newdata} gives a design matrix that doesn't match the one of the fit.", "i" = "Check that its covariates have the types and factor levels of the data the model was fitted to."), arg = "newdata")
  }
  design
}

#' Assemble prediction-location random-effect samples
#' @noRd
prediction_random_effect_samples <- function(prediction_setup, n_locations,
                                             n_samples, include) {
  if (!include) return(NULL)
  result <- matrix(0, nrow = n_locations, ncol = n_samples)
  for (term in seq_along(prediction_setup$re$samples)) {
    term_samples <- prediction_setup$re$samples[[term]]
    for (group in seq_along(term_samples)) {
      result <- result + outer(
        prediction_setup$re$D_pred[[term]][, group], term_samples[[group]]
      )
    }
  }
  result
}

#' Draw an optional prediction-location nugget component
#' @noRd
prediction_nugget_samples <- function(object, prediction_setup, n_locations,
                                      n_samples, include) {
  if (!include) return(NULL)
  tau2 <- fitted_nugget_variance(
    object, prediction_setup$par_hat
  )
  if (is.null(tau2) || !is.finite(tau2) || tau2 <= 0) {
    abort_input("{.arg include} can't contain {.val nugget}, because the fit has no positive nugget variance.", arg = "include")
  }
  sample_independent_gaussian(
    mean = rep(0, n_locations),
    standard_deviation = rep(sqrt(tau2), n_locations),
    n_samples = n_samples
  )
}

#' Obtain the fitted or fixed nugget variance
#' @noRd
fitted_nugget_variance <- function(object, par_hat = coef(object)) {
  if (!is.null(par_hat$tau2)) return(as.numeric(par_hat$tau2))
  if (is.numeric(object$fix_tau2) && !isTRUE(object$fix_tau2)) {
    return(as.numeric(object$fix_tau2))
  }
  0
}

#' Aggregate a location-by-sample matrix within one polygon
#' @noRd
aggregate_prediction_values <- function(values, operation, weights) {
  if (!is.matrix(values)) values <- as.matrix(values)
  if (is.function(operation)) {
    result <- operation(values, weights)
    if (!is.numeric(result) || length(result) != ncol(values)) {
      abort_input("A function in {.arg aggregate} must return one number for each predictive draw.", arg = "aggregate")
    }
    return(matrix(as.numeric(result), nrow = 1L))
  }
  if (operation == "mean") {
    total_weight <- sum(weights)
    if (total_weight <= 0) {
      return(matrix(NA_real_, nrow = 1L, ncol = ncol(values)))
    }
    return(matrix(
      as.numeric(crossprod(weights / total_weight, values)), nrow = 1L
    ))
  }
  if (operation == "sum") {
    if (sum(weights) <= 0) {
      return(matrix(NA_real_, nrow = 1L, ncol = ncol(values)))
    }
    return(matrix(as.numeric(crossprod(weights, values)), nrow = 1L))
  }
  if (operation == "median") {
    if (sum(weights) <= 0) {
      return(matrix(NA_real_, nrow = 1L, ncol = ncol(values)))
    }
    return(matrix(
      apply(values, 2L, weighted_prediction_median, weights = weights),
      nrow = 1L
    ))
  }
  reducer <- if (operation == "min") min else max
  matrix(apply(values, 2L, reducer), nrow = 1L)
}

#' Calculate a deterministic weighted median for one predictive draw
#' @noRd
weighted_prediction_median <- function(values, weights) {
  keep <- weights > 0
  values <- values[keep]
  weights <- weights[keep]
  ordering <- order(values)
  values <- values[ordering]
  cumulative <- cumsum(weights[ordering]) / sum(weights)
  index <- which(cumulative >= 0.5)[1L]
  at_midpoint <- abs(cumulative[index] - 0.5) <= sqrt(.Machine$double.eps)
  if (at_midpoint && index < length(values)) {
    return(mean(values[c(index, index + 1L)]))
  }
  values[index]
}

#' Calculate a stable coefficient of variation
#' @noRd
prediction_coefficient_of_variation <- function(values) {
  predictive_mean <- mean(values)
  numerical_tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(values)))
  if (!is.finite(predictive_mean) ||
      predictive_mean <= numerical_tolerance) {
    return(NA_real_)
  }
  sd(values) / predictive_mean
}

#' Return tidy summaries for one aggregated predictive distribution
#' @noRd
tidy_prediction_summaries <- function(values, boundary_id, aggregate,
                                      quantity, summaries, thresholds, level) {
  summary_functions <- prediction_summary_functions(level)
  draws <- as.numeric(values)
  ordinary <- drop_prediction_operations(summaries, "exceedance")
  rows <- lapply(names(ordinary), function(name) {
    operation <- ordinary[[name]]
    value <- if (all(is.na(draws))) {
      NA_real_
    } else if (is.function(operation)) {
      operation(draws)
    } else {
      summary_functions[[operation]](draws)
    }
    if (!is.numeric(value) || length(value) != 1L) {
      abort_input("The function {.val {name}} in {.arg summaries} must return a single number.", arg = "summaries")
    }
    data.frame(
      boundary_id = boundary_id,
      aggregate = aggregate,
      quantity = quantity,
      summary = name,
      threshold = NA_real_,
      value = as.numeric(value),
      stringsAsFactors = FALSE
    )
  })
  if (prediction_operation_requested(summaries, "exceedance")) {
    rows <- c(rows, lapply(thresholds, function(threshold) {
      data.frame(
        boundary_id = boundary_id,
        aggregate = aggregate,
        quantity = quantity,
        summary = "exceedance",
        threshold = threshold,
        value = mean(draws > threshold),
        stringsAsFactors = FALSE
      )
    }))
  }
  do.call(rbind, rows)
}

#' Empty tidy polygon-prediction summary
#' @noRd
empty_prediction_summary <- function() {
  data.frame(
    boundary_id = integer(),
    aggregate = character(),
    quantity = character(),
    summary = character(),
    threshold = numeric(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Add requested missing summaries when no prediction row is usable
#' @noRd
format_missing_location_predictions <- function(prediction_input, include,
                                                transform, summaries,
                                                thresholds) {
  predictions <- prediction_input$output_data
  predictions$fixed <- NA_real_
  predictions$offset <- NA_real_
  ordinary <- setdiff(names(summaries), c("cv", "exceedance"))
  requested_names <- c(
    paste("spatial", ordinary, sep = "_"),
    if ("random_effects" %in% include)
      paste("random_effects", ordinary, sep = "_"),
    if ("nugget" %in% include) paste("nugget", ordinary, sep = "_"),
    paste("link", ordinary, sep = "_"),
    prediction_summary_names("response", summaries, thresholds),
    unlist(lapply(names(transform), prediction_summary_names,
                  summaries = summaries, thresholds = thresholds),
           use.names = FALSE)
  )
  for (name in requested_names) predictions[[name]] <- NA_real_
  predictions
}

#' Add location-level summaries to the original prediction support
#' @noRd
format_location_predictions <- function(prediction_input, sample_components,
                                        include, transform, summaries,
                                        thresholds, level) {
  predictions <- prediction_input$output_data
  valid_rows <- prediction_input$valid_rows
  predictions$fixed <- NA_real_
  predictions$offset <- NA_real_
  predictions$fixed[valid_rows] <- sample_components$fixed
  predictions$offset[valid_rows] <- sample_components$offset
  component_summaries <- c(
    summarise_prediction_matrix(
      sample_components$spatial, "spatial",
      drop_prediction_operations(summaries, c("cv", "exceedance")), level
    ),
    if (!is.null(sample_components$random_effects))
      summarise_prediction_matrix(
        sample_components$random_effects, "random_effects",
        drop_prediction_operations(summaries, c("cv", "exceedance")), level
      ),
    if (!is.null(sample_components$nugget))
      summarise_prediction_matrix(
        sample_components$nugget, "nugget",
        drop_prediction_operations(summaries, c("cv", "exceedance")), level
      ),
    summarise_prediction_matrix(
      sample_components$link, "link",
      drop_prediction_operations(summaries, c("cv", "exceedance")), level
    ),
    summarise_prediction_matrix(
      sample_components$response, "response", summaries, level, thresholds
    ),
    unlist(lapply(names(transform), function(name) {
      summarise_prediction_matrix(
        sample_components$transformed[[name]], name, summaries, level,
        thresholds
      )
    }), recursive = FALSE)
  )
  for (name in names(component_summaries)) {
    predictions[[name]] <- NA_real_
    predictions[[name]][valid_rows] <- component_summaries[[name]]
  }
  predictions
}

#' Return tidy missing summaries for a polygon without usable locations
#' @noRd
missing_prediction_summaries <- function(boundary_id, aggregate, quantity,
                                         summaries, thresholds) {
  ordinary <- drop_prediction_operations(summaries, "exceedance")
  rows <- lapply(names(ordinary), function(name) {
    data.frame(
      boundary_id = boundary_id, aggregate = aggregate, quantity = quantity,
      summary = name, threshold = NA_real_, value = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  if (prediction_operation_requested(summaries, "exceedance")) {
    rows <- c(rows, lapply(thresholds, function(threshold) {
      data.frame(
        boundary_id = boundary_id, aggregate = aggregate, quantity = quantity,
        summary = "exceedance", threshold = threshold, value = NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }
  do.call(rbind, rows)
}

#' Return missing summaries for boundaries without usable prediction rows
#' @noRd
missing_boundary_summaries <- function(n_boundaries, aggregate, include,
                                       transform, summaries, thresholds) {
  rows <- list()
  row_number <- 0L
  quantities <- c(
    "fixed", "offset", "spatial",
    if ("random_effects" %in% include) "random_effects",
    if ("nugget" %in% include) "nugget",
    "link", "response", transform
  )
  for (boundary_id in seq_len(n_boundaries)) {
    for (operation in aggregate) {
      for (quantity in quantities) {
        row_number <- row_number + 1L
        rows[[row_number]] <- if (quantity %in% c("fixed", "offset")) {
          data.frame(
            boundary_id = boundary_id, aggregate = operation,
            quantity = quantity, summary = "estimate", threshold = NA_real_,
            value = NA_real_, stringsAsFactors = FALSE
          )
        } else {
          requested <- if (quantity %in% c("response", transform)) summaries else
            drop_prediction_operations(summaries, c("cv", "exceedance"))
          missing_prediction_summaries(
            boundary_id, operation, quantity, requested,
            if (quantity %in% c("response", transform)) thresholds else NULL
          )
        }
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else empty_prediction_summary()
}

#' Built-in summaries of one predictive distribution
#' @noRd
prediction_summary_functions <- function(level) {
  alpha <- (1 - level) / 2
  list(
    mean = mean,
    median = median,
    sd = sd,
    cv = prediction_coefficient_of_variation,
    lower = function(x) quantile(x, alpha, names = FALSE),
    upper = function(x) quantile(x, 1 - alpha, names = FALSE),
    interval_width = function(x) {
      diff(quantile(x, c(alpha, 1 - alpha), names = FALSE))
    }
  )
}

#' Calculate built-in summaries for every row of a prediction matrix
#'
#' Gives, for each row, the values of `prediction_summary_functions()`
#' without calling them row by row. Each complete row is partially sorted
#' once, and the median and the type 7 quantiles behind `"lower"`, `"upper"`
#' and `"interval_width"` are formed from those order statistics with the
#' arithmetic of `median()` and `quantile()`, so they are identical to the
#' row-by-row values. `"mean"`, `"sd"` and `"cv"` use `rowMeans()` and
#' `rowSums()`. Like `mean()` and `sd()`, these are within about one unit in
#' the last place of the exact mean and standard deviation of the draws, but
#' they round differently, so the two can differ by one unit in the last
#' place. Rows
#' with missing values are passed to the row-by-row functions, so they give
#' the same `NA` values and errors as those functions.
#' @noRd
prediction_row_summaries <- function(values, operations, level) {
  operations <- unique(operations)
  if (!length(operations)) return(list())
  complete <- !is.na(rowSums(values))
  if (!all(complete)) {
    by_row <- prediction_summary_functions(level)
    incomplete <- values[!complete, , drop = FALSE]
    result <- lapply(setNames(nm = operations), function(operation) {
      value <- rep(NA_real_, nrow(values))
      value[!complete] <- apply(incomplete, 1L, by_row[[operation]])
      value
    })
    complete_result <- prediction_row_summaries(
      values[complete, , drop = FALSE], operations, level
    )
    for (operation in operations) {
      result[[operation]][complete] <- complete_result[[operation]]
    }
    return(result)
  }

  result <- list()
  n_draws <- ncol(values)
  if (any(c("mean", "sd", "cv") %in% operations)) {
    result$mean <- rowMeans(values)
  }
  if (any(c("sd", "cv") %in% operations)) {
    result$sd <- if (n_draws > 1L) {
      deviations <- values - result$mean
      sqrt(rowSums(deviations * deviations) / (n_draws - 1L))
    } else {
      rep(NA_real_, nrow(values))
    }
  }
  if ("cv" %in% operations) {
    magnitude <- abs(values)
    largest <- magnitude[cbind(
      seq_len(nrow(values)), max.col(magnitude, ties.method = "first")
    )]
    tolerance <- sqrt(.Machine$double.eps) * pmax(1, largest)
    result$cv <- result$sd / result$mean
    result$cv[!is.finite(result$mean) | result$mean <= tolerance] <- NA_real_
  }

  needs_quantiles <- any(c("lower", "upper", "interval_width") %in% operations)
  needs_median <- "median" %in% operations
  if (needs_quantiles || needs_median) {
    # Positions of the order statistics used by median() and by
    # quantile(type = 7) at the two interval probabilities
    alpha <- (1 - level) / 2
    index <- 1 + max(n_draws - 1, 0) * c(alpha, 1 - alpha)
    below <- floor(index)
    above <- ceiling(index)
    half <- (n_draws + 1L) %/% 2L
    middle <- if (n_draws %% 2L == 1L) half else half + 0L:1L
    positions <- unique(c(
      if (needs_quantiles) c(below, above),
      if (needs_median) middle
    ))
    draws <- t(values)
    n_statistics <- length(positions) + needs_median
    statistics <- vapply(seq_len(ncol(draws)), function(i) {
      sorted <- sort.int(draws[, i], partial = positions)
      c(
        sorted[positions],
        if (needs_median) {
          if (length(middle) == 1L) sorted[middle] else mean(sorted[middle])
        }
      )
    }, numeric(n_statistics))
    dim(statistics) <- c(n_statistics, ncol(draws))
    if (needs_median) result$median <- statistics[n_statistics, ]
    if (needs_quantiles) {
      type_seven <- function(j) {
        low <- statistics[match(below[j], positions), ]
        high <- statistics[match(above[j], positions), ]
        weight <- index[j] - below[j]
        interpolate <- index[j] > below[j] & high != low
        low[interpolate] <- (1 - weight) * low[interpolate] +
          weight * high[interpolate]
        low
      }
      result$lower <- type_seven(1L)
      result$upper <- type_seven(2L)
      result$interval_width <- result$upper - result$lower
    }
  }
  result[operations]
}

#' Summarise a location-by-sample prediction matrix
#' @noRd
summarise_prediction_matrix <- function(values, quantity, summaries, level,
                                        thresholds = NULL) {
  ordinary <- drop_prediction_operations(summaries, "exceedance")
  built_in <- vapply(ordinary, is.character, logical(1))
  row_summaries <- prediction_row_summaries(
    values, unlist(ordinary[built_in], use.names = FALSE), level
  )
  result <- lapply(names(ordinary), function(name) {
    operation <- ordinary[[name]]
    value <- if (is.function(operation)) {
      apply(values, 1L, operation)
    } else {
      row_summaries[[operation]]
    }
    if (!is.numeric(value) || length(value) != nrow(values)) {
      abort_input("The function {.val {name}} in {.arg summaries} must return a single number.", arg = "summaries")
    }
    as.numeric(value)
  })
  names(result) <- paste(quantity, names(ordinary), sep = "_")
  if (prediction_operation_requested(summaries, "exceedance")) {
    exceedance <- lapply(thresholds, function(threshold) {
      rowMeans(values > threshold)
    })
    names(exceedance) <- paste(
      quantity, "exceedance", prediction_threshold_label(thresholds), sep = "_"
    )
    result <- c(result, exceedance)
  }
  result
}

#' Construct location-level prediction column names
#' @noRd
prediction_summary_names <- function(quantity, summaries, thresholds = NULL) {
  ordinary <- names(drop_prediction_operations(summaries, "exceedance"))
  result <- paste(quantity, ordinary, sep = "_")
  if (prediction_operation_requested(summaries, "exceedance")) {
    result <- c(result, paste(
      quantity, "exceedance", prediction_threshold_label(thresholds), sep = "_"
    ))
  }
  result
}

#' Format thresholds for stable column names
#' @noRd
prediction_threshold_label <- function(thresholds) {
  format(thresholds, scientific = FALSE, trim = TRUE)
}

#' Restore the block structure of stored fitted conditional samples
#' @noRd
stored_conditioning_sample <- function(object) {
  block_sizes <- c(nrow(object$coords), lengths(object$re))
  if (ncol(object$S_samples) != sum(block_sizes)) {
    abort_input("The conditional samples stored in {.arg object} don't match its model.", arg = "object")
  }
  block_ends <- cumsum(block_sizes)
  block_starts <- c(1L, utils::head(block_ends, -1L) + 1L)
  blocks <- Map(
    function(first, last) object$S_samples[, first:last, drop = FALSE],
    block_starts, block_ends
  )
  names(blocks) <- c("S", names(object$re))
  list(samples = blocks)
}

#' Reconstruct linear-predictor draws from a prediction object
#'
#' The public prediction object keeps deterministic contributions in `predictions`
#' and stochastic contributions in `samples`. Keeping that representation
#' compact avoids duplicating large matrices, while this helper gives internal
#' assessment routines the complete draws they need.
#' @noRd
prediction_link_samples <- function(prediction, include_fixed = TRUE,
                                    include_offset = TRUE) {
  samples <- prediction$samples
  if (is.null(samples) || is.null(samples$spatial)) {
    riskmap_bug("The prediction has no samples.")
  }

  link <- samples$spatial
  if (include_fixed) {
    # With parameter uncertainty the fixed effects differ between draws.
    link <- if (is.null(samples$fixed)) {
      sweep(link, 1L, prediction$predictions$fixed, `+`)
    } else {
      link + samples$fixed
    }
  }
  if (include_offset) {
    link <- sweep(link, 1L, prediction$predictions$offset, `+`)
  }
  if (!is.null(samples$random_effects)) {
    link <- link + samples$random_effects
  }
  if (!is.null(samples$nugget)) {
    link <- link + samples$nugget
  }
  link
}

gaussian_conditioning_distribution <- function(object, spatial_covariance,
                                               par_hat, residual) {
  n_observations <- length(object$y)
  n_spatial <- nrow(spatial_covariance)
  random_ids <- as.matrix(object$ID_re)
  n_random_terms <- if (is.null(object$ID_re)) 0L else ncol(random_ids)

  spatial_design <- matrix(0, nrow = n_observations, ncol = n_spatial)
  spatial_design[cbind(seq_len(n_observations), object$ID_coords)] <- 1
  designs <- list(spatial_design)
  covariances <- list(spatial_covariance)
  if (n_random_terms) {
    for (term in seq_len(n_random_terms)) {
      design <- matrix(0, nrow = n_observations,
                       ncol = length(object$re[[term]]))
      design[cbind(seq_len(n_observations), random_ids[, term])] <- 1
      designs[[term + 1L]] <- design
      covariances[[term + 1L]] <- diag(
        par_hat$sigma2_re[[term]], nrow = ncol(design)
      )
    }
  }

  incidence <- do.call(cbind, designs)
  latent_covariance <- as.matrix(Matrix::bdiag(covariances))
  measurement_variance <- par_hat$sigma2_me %||% object$fix_var_me %||% 0
  ids <- as.matrix(cbind(object$ID_coords, object$ID_re))
  prediction_weights <- gaussian_prediction_weights(
    Matrix::Matrix(incidence, sparse = TRUE), latent_covariance, ids,
    measurement_variance
  )
  cross_covariance <- latent_covariance %*% t(incidence)
  gain <- prediction_weights(cross_covariance)
  conditional_mean <- as.numeric(gain %*% residual)
  conditional_covariance <- latent_covariance -
    gain %*% incidence %*% latent_covariance
  random_rows <- if (n_random_terms) {
    split(seq.int(n_spatial + 1L, nrow(latent_covariance)),
          rep(seq_len(n_random_terms), lengths(object$re)))
  } else {
    list()
  }
  list(mean = conditional_mean, covariance = conditional_covariance,
       spatial_rows = seq_len(n_spatial), random_rows = random_rows)
}

#' Draw the fitted latent state for an exact Gaussian model
#'
#' The shared columns are reused by every prediction block, preserving the
#' posterior dependence between the fitted spatial state and group effects.
#' @noRd
gaussian_conditioning_samples <- function(object, spatial_covariance,
                                          par_hat, residual, n_samples) {
  distribution <- gaussian_conditioning_distribution(
    object, spatial_covariance, par_hat, residual
  )
  conditional_root <- factor_covariance(
    distribution$covariance, "Gaussian fitted latent-state covariance"
  )
  samples <- sample_correlated_gaussian(
    distribution$mean, t(conditional_root), n_samples
  )
  list(
    spatial = samples[distribution$spatial_rows, , drop = FALSE],
    random_effects = lapply(distribution$random_rows, function(rows) {
      samples[rows, , drop = FALSE]
    })
  )
}

#' Generate conditional spatial-process samples from aligned prediction data
#'
#' Prediction variables and geometry remain in one row-aligned `sf` object.
#' For areal predictions, `groups` defines the independent joint covariance
#' blocks without splitting the data into parallel lists.
#' A fresh conditional MCMC run of a non-Gaussian model is seeded with
#' `mcmc_seed`, after which the caller's random number state is restored, so
#' the predictive draws that follow start from that state.
#' @noRd
draw_prediction_samples <- function(object, prediction_data, prediction,
                                    include, control_mcmc = NULL,
                                    groups = NULL,
                                    observed_locations = FALSE,
                                    mcmc_seed = NULL,
                                    messages = TRUE) {

  control_mcmc_supplied <- !is.null(control_mcmc)
  if (is.null(control_mcmc)) control_mcmc <- set_control_mcmc()
  type <- prediction
  plain_data <- as.data.frame(st_drop_geometry(prediction_data))
  grid_pred <- if (observed_locations) NULL else st_geometry(prediction_data)
  predictors <- if (observed_locations) NULL else plain_data
  inter_f <- interpret_formula(object$formula)
  pred_cov_offset <- NULL
  if (!observed_locations && !is.null(inter_f$offset)) {
    if (!inter_f$offset %in% names(plain_data)) {
      abort_input("{.arg newdata} must have the offset column {.val {inter_f$offset}}.", arg = "newdata")
    }
    pred_cov_offset <- plain_data[[inter_f$offset]]
  }
  re_predictors <- NULL
  if ("random_effects" %in% include) {
    random_names <- colnames(object$ID_re)
    if (!length(random_names)) {
      abort_input("{.arg include} can't contain {.val random_effects}, because the model has no {.fn re} terms.", arg = "include")
    }
    if (prediction != "joint") {
      abort_input("{.arg type} must be {.val joint} to predict the random effects.", arg = "type")
    }
    prepare_random_effects(plain_data, random_names)
    re_predictors <- plain_data[random_names]
  }

  # ---------------------------------------------------------------------------
  # validate inputs
  # ---------------------------------------------------------------------------
  if (!inherits(object, "RiskMap")) riskmap_bug("{.arg object} must be a {.cls RiskMap} fit.")

  model_crs <- st_crs(object$data)

  if (is.na(model_crs)) {
    abort_input("The data of {.arg object} must have a coordinate reference system.", arg = "object")
  }

  crs_label <- function(x) {
    x_crs <- st_crs(x)
    if (!is.na(x_crs$epsg)) paste0("EPSG:", x_crs$epsg) else x_crs$input
  }

  if (!is.null(grid_pred)) {
    check_data(grid_pred, type = "sfc")
    if (!isTRUE(st_crs(grid_pred) == model_crs)) {
      riskmap_bug("The prediction locations are in {crs_label(grid_pred)}, not in {crs_label(model_crs)} as the fit.")
    }
  }

  check_control(control_mcmc, "set_control_mcmc")

  if (!type %in% c("marginal", "joint")) {
    riskmap_bug("{.arg type} must be {.val marginal} or {.val joint}.")
  }

  grouped_joint <- !is.null(groups)
  if (grouped_joint) {
    if (type != "joint") {
      riskmap_bug("Grouped prediction needs {.code type = \"joint\"}.")
    }
    if (is.null(grid_pred) || length(groups) != length(grid_pred) ||
        anyNA(groups)) {
      riskmap_bug("{.arg groups} must give a group for every prediction location.")
    }
    group_indices <- unname(split(seq_along(groups), groups))
  } else {
    group_indices <- list(seq_len(if (is.null(grid_pred)) {
      nrow(object$coords)
    } else {
      length(grid_pred)
    }))
  }

  if (!is.null(grid_pred) && is.null(predictors) && ncol(object$D) > 1) {
    riskmap_bug("Covariates must be given with the prediction locations.")
  }

  obs_loc <- is.null(grid_pred)
  if (obs_loc) {
    if (!is.null(predictors)) {
      riskmap_bug("Covariates were given without prediction locations.")
    }
    predictors <- as.data.frame(st_drop_geometry(object$data))
    grid_pred  <- st_as_sfc(object$data)
  }

  grp <- coordinates_in_units(grid_pred, object$distance_units)
  n_pred <- nrow(grp)

  object$D <- as.matrix(object$D)
  p <- ncol(object$D)

  if (!all(object$cov_offset == 0)) {
    if (obs_loc){
      pred_cov_offset <- object$cov_offset
    } else {
      if (is.null(pred_cov_offset)) {
        riskmap_bug("The offset must be given at the prediction locations.")
      }
      if (!is.numeric(pred_cov_offset)) {
        abort_input("The offset column {.val {inter_f$offset}} of {.arg newdata} must be numeric.", arg = "newdata")
      }
      if (length(pred_cov_offset) != n_pred) {
        riskmap_bug("The offset must have one value for each prediction location.")
      }
    }
  } else {
    if (!is.null(pred_cov_offset)) {
      riskmap_warn("The offset of {.arg newdata} is ignored, because the offset of the fit is zero everywhere.")
    }
    pred_cov_offset <- 0
  }

  # ---------------------------------------------------------------------------
  # Extract terms from object
  # ---------------------------------------------------------------------------

  par_hat <- coef(object)
  inter_f      <- interpret_formula(object$formula)

  intercept_only <- if (p == 1) prod(object$D[, 1] == 1) == 1 else FALSE

  n_re <- length(object$re)
  if (n_re > 0 && type == "marginal" && !is.null(re_predictors)) {
    riskmap_bug("Random effects are predicted only with {.code type = \"joint\"}.")
  }

  # ---------------------------------------------------------------------------
  # Build mu_pred (fixed-effects linear predictor at prediction locations)
  # ---------------------------------------------------------------------------

  if (!is.null(predictors)) {

    .build_D_pred <- function(predictors, n_predictors) {
      response <- as.character(object$formula[[2]])
      re_names <- names(object$re)
      offset_names <- names(object$cov_offset)
      model_predictors <- setdiff(get_formula_terms(object$formula), c(response, re_names, offset_names))

      if (!is.data.frame(predictors) ||
          !all(model_predictors %in% names(predictors)) ||
          nrow(predictors) != n_predictors) {
        riskmap_bug("The covariates must be a data frame with the columns of the model and one row per prediction location.")
      }
      check_complete_data(predictors, model_predictors)
      prediction_design_matrix(object, predictors)
    }

    D_pred <- .build_D_pred(predictors, n_pred)
    mu_pred <- as.numeric(D_pred %*% par_hat$beta)

  } else if (intercept_only) {
    mu_pred <- par_hat$beta
  } else {
    mu_pred <- 0
  }

  # ---------------------------------------------------------------------------
  # Random effects (unchanged from original)
  # ---------------------------------------------------------------------------
  if (n_re > 0) {
    ID_g         <- as.matrix(cbind(object$ID_coords, object$ID_re))
    re_unique    <- object$re
    n_dim_re_tot <- sapply(seq_len(n_re + 1), function(i) length(unique(ID_g[, i])))

    if (!is.null(re_predictors)) {
      if (any(is.na(re_predictors))) {
        riskmap_warn("Locations with missing levels of a random effect are left out.")
        ind_c       <- complete.cases(re_predictors)
        re_predictors <- re_predictors[ind_c, , drop = FALSE]
        grid_pred <- grid_pred[ind_c]
        grp <- coordinates_in_units(grid_pred, object$distance_units)
        n_pred <- nrow(grp)
      }
      if (!is.data.frame(re_predictors) || nrow(re_predictors) != n_pred ||
          ncol(re_predictors) != n_re) {
        riskmap_bug("The random-effect levels must be a data frame with one column per random effect and one row per prediction location.")
      }

      n_dim_re <- sapply(seq_len(n_re), function(i) length(unique(object$ID_re[, i])))
      build_re_design <- function(values, size) {
        lapply(seq_len(n_re), function(i) {
          mat <- matrix(0, nrow = size, ncol = n_dim_re[i])
          for (j in seq_along(object$re[[i]])) {
            mat[values[, i] == object$re[[i]][j], j] <- 1
          }
          mat
        })
      }
      D_re_pred <- build_re_design(re_predictors, n_pred)
    } else {
      D_re_pred <- NULL
    }
  } else {
    D_re_pred <- NULL
  }

  # ---------------------------------------------------------------------------
  # Spatial quantities
  # ---------------------------------------------------------------------------
  out <- list(mu_pred = mu_pred, grid_pred = grid_pred, par_hat = par_hat)
  distance_scale <- attr(object, "distance_scale") %||% 1
  fitting_coords <- object$coords / distance_scale
  fitting_phi <- par_hat$phi / distance_scale
  fitting_grp <- grp / distance_scale

  conditioning_coords <- if (object$family == "gaussian") {
    fitting_coords[object$ID_coords, , drop = FALSE]
  } else {
    fitting_coords
  }
  marginal_batch_size <- if (!obs_loc && type == "marginal") {
    marginal_prediction_batch_size(n_pred, nrow(conditioning_coords))
  } else {
    n_pred
  }
  batch_marginal <- !obs_loc && type == "marginal" &&
    marginal_batch_size < n_pred

  if (object$family != "gaussian" && !obs_loc && !batch_marginal &&
      !grouped_joint) {
    U_pred <- cross_distances(fitting_grp, fitting_coords)
  } else if (object$family == "gaussian" && !obs_loc && !batch_marginal &&
             !grouped_joint) {
    U_pred <- cross_distances(fitting_grp, conditioning_coords)
  }

  U <- pairwise_distances(fitting_coords)
  R <- matern_correlation(U, phi = fitting_phi, kappa = object$kappa,
                          return_sym_matrix = TRUE)

  if (obs_loc) {
    C <- par_hat$sigma2 * R[, object$ID_coords]
    grp <- object$coords
    fitting_grp <- fitting_coords
  } else if (!batch_marginal && !grouped_joint) {
    C <- par_hat$sigma2 * matern_correlation(
      U_pred, phi = fitting_phi, kappa = object$kappa
    )
  }

  n_pred_spatial <- if (obs_loc) nrow(object$coords) else n_pred

  mu <- as.numeric(object$D %*% par_hat$beta + object$cov_offset)

  reuse_fitted_samples <- object$family != "gaussian" &&
    !control_mcmc_supplied && !is.null(object$S_samples)
  n_samples <- if (reuse_fitted_samples) {
    nrow(object$S_samples)
  } else if (control_mcmc$linear_model) {
    control_mcmc$n_sim
  } else {
    (control_mcmc$n_sim - control_mcmc$burnin) / control_mcmc$thin
  }

  # ---------------------------------------------------------------------------
  # FIX 2: nu2 / nugget
  # STH: tau2 fixed at 0 -> near-zero for numerical stability
  # LF:  tau2 may be estimated -> use it
  # glgpm: existing behaviour
  # ---------------------------------------------------------------------------
  tau2 <- fitted_nugget_variance(object, par_hat)
  nu2 <- tau2 / par_hat$sigma2
  if (nu2 == 0) nu2 <- 1e-10
  diag(R) <- diag(R) + nu2

  diff.y <- if (object$family == "gaussian") object$y - mu else NULL
  shared_gaussian_re_samples <- NULL

  # ===========================================================================
  # NON-GAUSSIAN MODELS
  # ===========================================================================
  if (object$family != "gaussian") {

    Sigma     <- par_hat$sigma2 * R
    Sigma_root <- factor_covariance(Sigma, "observation covariance")

    prediction_weights <- function(cross_covariance) {
      cholesky_prediction_weights(cross_covariance, Sigma_root)
    }

    if (!obs_loc && !batch_marginal && !grouped_joint) {
      A <- prediction_weights(C)
    }

    if (reuse_fitted_samples) {
      simulation <- stored_conditioning_sample(object)
      if (messages) {
        riskmap_inform("Reusing the conditional samples stored in the fit.")
      }
    } else {
      simulation <- laplace_sampling_mcmc(
        y = object$y, units_m = object$units_m, mu = mu, Sigma = Sigma,
        sigma2_re = par_hat$sigma2_re,
        invlink = if (identical(object$link_function$name, "canonical")) {
          NULL
        } else {
          object$link_function[c("inv", "d1", "d2")]
        },
        ID_coords = object$ID_coords, ID_re = object$ID_re,
        family = object$family, control_mcmc = control_mcmc,
        seed = mcmc_seed
      )
    }

    if (obs_loc) {
      out$S_samples <- t(simulation$samples$S)
    } else if (batch_marginal) {
      bar <- progress_bar("Predicting at locations", nrow(fitting_grp),
                          messages)
      out$S_samples <- batched_marginal_prediction(
        fitting_grp, conditioning_coords, prediction_weights,
        t(simulation$samples$S), par_hat$sigma2, fitting_phi,
        object$kappa, n_samples, marginal_batch_size, progress = bar
      )
      progress_done(bar)
    } else if (grouped_joint) {
      out$S_samples <- grouped_joint_prediction(
        prediction_coordinates = fitting_grp,
        conditioning_coordinates = conditioning_coords,
        group_indices = group_indices,
        weight_function = prediction_weights,
        conditional_signal = t(simulation$samples$S),
        marginal_variance = par_hat$sigma2,
        phi = fitting_phi,
        kappa = object$kappa,
        n_samples = n_samples
      )
    } else {
      mu_cond_S <- A %*% t(simulation$samples$S)

      if (type == "marginal") {
        sd_cond_S <- sqrt(conditional_variances(par_hat$sigma2, A, C))
        out$S_samples <- sample_independent_gaussian(
          mu_cond_S, sd_cond_S, n_samples
        )

      } else {
        Sp  <- par_hat$sigma2 * matern_correlation(pairwise_distances(fitting_grp), phi = fitting_phi,
                                           kappa = object$kappa, return_sym_matrix = TRUE)
        Sc  <- Sp - A %*% t(C)
        Sc_root <- factor_covariance(
          Sc, "conditional prediction covariance"
        )
        Scr <- t(Sc_root)
        out$S_samples <- sample_correlated_gaussian(
          mu_cond_S, Scr, n_samples
        )
      }
    }

  } else {
    # =========================================================================
    # GAUSSIAN MODEL
    # =========================================================================
    blockwise_random_effects <- grouped_joint && !is.null(re_predictors)
    if (blockwise_random_effects) {
      spatial_covariance <- par_hat$sigma2 * R
      latent_samples <- gaussian_conditioning_samples(
        object = object,
        spatial_covariance = spatial_covariance,
        par_hat = par_hat,
        residual = diff.y,
        n_samples = n_samples
      )
      spatial_root <- factor_covariance(
        spatial_covariance, "Gaussian fitted spatial covariance"
      )
      spatial_weights <- function(cross_covariance) {
        cholesky_prediction_weights(cross_covariance, spatial_root)
      }
      out$S_samples <- grouped_joint_prediction(
        prediction_coordinates = fitting_grp,
        conditioning_coordinates = fitting_coords,
        group_indices = group_indices,
        weight_function = spatial_weights,
        conditional_signal = latent_samples$spatial,
        marginal_variance = par_hat$sigma2,
        phi = fitting_phi,
        kappa = object$kappa,
        n_samples = n_samples
      )
      shared_gaussian_re_samples <- latent_samples$random_effects
    } else if (!is.null(object$fix_var_me) && object$fix_var_me > 0 || is.null(object$fix_var_me)) {
      m        <- length(object$y)
      s_unique <- unique(object$ID_coords)
      ID_g     <- as.matrix(cbind(object$ID_coords, object$ID_re))
      n_dim_re_tot <- sapply(seq_len(n_re + 1), function(i) length(unique(ID_g[, i])))
      C_g      <- matrix(0, nrow = m, ncol = sum(n_dim_re_tot))
      for (i in seq_len(m)) {
        ind_s_i <- which(s_unique == ID_g[i, 1])
        C_g[i, seq_len(n_dim_re_tot[1])][ind_s_i] <- 1
      }
      if (n_re > 0) {
        for (j in seq_len(n_re)) {
          select_col <- sum(n_dim_re_tot[seq_len(j)])
          for (i in seq_len(m)) {
            ind_re_j_i <- which(re_unique[[j]] == ID_g[i, j + 1])
            C_g[i, select_col + seq_len(n_dim_re_tot[j + 1])][ind_re_j_i] <- 1
          }
        }
      }
      C_g      <- Matrix(C_g, sparse = TRUE, doDiag = FALSE)
      Sigma_g  <- matrix(0, nrow = sum(n_dim_re_tot), ncol = sum(n_dim_re_tot))
      Sigma_g[seq_len(n_dim_re_tot[1]), seq_len(n_dim_re_tot[1])] <- par_hat$sigma2 * R
      if (n_re > 0) {
        for (j in seq_len(n_re)) {
          sc <- sum(n_dim_re_tot[seq_len(j)])
          diag(Sigma_g[sc + seq_len(n_dim_re_tot[j+1]), sc + seq_len(n_dim_re_tot[j+1])]) <- par_hat$sigma2_re[j]
        }
      }
      # coef() omits sigma2_me when the measurement error variance is fixed
      measurement_variance <- par_hat$sigma2_me %||% object$fix_var_me
      prediction_weights <- gaussian_prediction_weights(C_g, Sigma_g, ID_g,
                                                         measurement_variance)
      if (!batch_marginal && !grouped_joint) {
        A <- prediction_weights(C)
      }
    } else {
      Sigma     <- par_hat$sigma2 * R
      Sigma_root <- factor_covariance(Sigma, "observation covariance")
      prediction_weights <- function(cross_covariance) {
        cholesky_prediction_weights(cross_covariance, Sigma_root)
      }
      if (!batch_marginal && !grouped_joint) {
        A <- prediction_weights(C)
      }
    }

    if (blockwise_random_effects) {
      # Samples were generated above from one shared fitted latent state.
    } else if (batch_marginal) {
      bar <- progress_bar("Predicting at locations", nrow(fitting_grp),
                          messages)
      out$S_samples <- batched_marginal_prediction(
        fitting_grp, conditioning_coords, prediction_weights, diff.y,
        par_hat$sigma2, fitting_phi, object$kappa, n_samples,
        marginal_batch_size, progress = bar
      )
      progress_done(bar)
    } else if (grouped_joint) {
      out$S_samples <- grouped_joint_prediction(
        prediction_coordinates = fitting_grp,
        conditioning_coordinates = conditioning_coords,
        group_indices = group_indices,
        weight_function = prediction_weights,
        conditional_signal = diff.y,
        marginal_variance = par_hat$sigma2,
        phi = fitting_phi,
        kappa = object$kappa,
        n_samples = n_samples
      )
    } else {
      mu_cond_S <- as.numeric(A %*% diff.y)

      if (type == "marginal") {
        sd_cond_S <- sqrt(conditional_variances(par_hat$sigma2, A, C))
        out$S_samples <- sample_independent_gaussian(
          mu_cond_S, sd_cond_S, n_samples
        )
      } else {
        Sp <- par_hat$sigma2 * matern_correlation(
          pairwise_distances(fitting_grp), phi = fitting_phi,
          kappa = object$kappa, return_sym_matrix = TRUE
        )
        Sc <- Sp - A %*% t(C)
        Sc_root <- factor_covariance(
          Sc, "conditional prediction covariance"
        )
        Scr <- t(Sc_root)
        out$S_samples <- sample_correlated_gaussian(
          mu_cond_S, Scr, n_samples
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Random effect posterior samples (unchanged from original)
  # ---------------------------------------------------------------------------

  if (n_re > 0 && !is.null(re_predictors)) {
    out$re          <- list()
    out$re$D_pred   <- D_re_pred
    out$re$samples  <- list()
    re_names        <- colnames(object$ID_re)
    if (object$family == "gaussian" && !is.null(shared_gaussian_re_samples)) {
      re_samples <- do.call(rbind, shared_gaussian_re_samples)
    } else if (object$family == "gaussian") {
      C_Z  <- C_g[, -(seq_len(n_dim_re_tot[1]))]
      add  <- 0
      for (i in seq_along(n_dim_re_tot[-1])) {
        C_Z[, add + seq_len(n_dim_re_tot[i+1])] <- par_hat$sigma2_re[i] *
          C_Z[, add + seq_len(n_dim_re_tot[i+1])]
        add <- n_dim_re_tot[i+1]
      }
      W_Z          <- prediction_weights(Matrix::t(C_Z))
      A_Z          <- cholesky_prediction_weights(W_Z %*% t(C), Sc_root)
      Sigma_Z_cond <- diag(rep(par_hat$sigma2_re, n_dim_re_tot[-1])) -
        W_Z %*% C_Z -
        A_Z %*% C %*% t(W_Z)
      Scr_Z        <- t(factor_covariance(
        Sigma_Z_cond, "conditional random-effect covariance"
      ))
      mu_Z_cond    <- A_Z %*% (out$S_samples - mu_cond_S)
      mu_Z_cond <- mu_Z_cond + as.numeric(W_Z %*% diff.y)
      re_samples <- sample_correlated_gaussian(
        mu_Z_cond, Scr_Z, n_samples
      )
    } else {
      re_samples <- matrix(0, nrow = sum(n_dim_re_tot[-1]), ncol = n_samples)
      add <- 0
      for (i in seq_len(n_re)) {
        rows <- add + seq_len(n_dim_re_tot[i + 1])
        re_samples[rows, ] <- t(simulation$samples[[i + 1]])
        add <- add + n_dim_re_tot[i+1]
      }
    }
    add <- 0
    for (i in seq_len(n_re)) {
      for (j in seq_len(n_dim_re_tot[1 + i])) {
        add <- add + 1
        out$re$samples[[re_names[[i]]]][[as.character(re_unique[[i]][[j]])]] <- re_samples[add, ]
      }
    }
  } else {
    out$re <- list(D_pred = NULL, samples = NULL)
  }

  out$obs_loc  <- obs_loc
  if (obs_loc){
    out$ID_coords <- object$ID_coords
  } else {
    out["ID_coords"] <- list(NULL)
  }
  out$inter_f  <- inter_f
  out$family   <- object$family
  out$par_hat  <- par_hat
  out$cov_offset <- pred_cov_offset
  out$type     <- type
  class(out)   <- "RiskMap_pred"
  return(out)
}

.anpit_area <- function(curve, u = seq(0, 1, length.out = length(curve))) {
  if (length(curve) != length(u)) {
    riskmap_bug("{.arg curve} and {.arg u} must have the same length.")
  }
  if (length(curve) < 2) return(NA_real_)

  dx <- diff(u)
  y <- abs(curve - u)
  sum(dx * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}

##' @title Plot Training/Test Splits
##'
##' @description
##' Plots every observation used by \code{\link{assess_prediction}}, coloured
##' by whether it fell in the training or test set for a given iteration,
##' faceted by iteration when there is more than one. Mimics the visual style
##' of \code{spatialsample::autoplot()} (used directly for
##' \code{method = "cluster"}), for the \code{"user"} and \code{"regularized"}
##' methods, whose splits aren't necessarily a complete partition of the data
##' and so aren't well suited to colouring by fold membership alone.
##'
##' @param data_split A list with a \code{splits} element, each entry itself a
##' list with an \code{sf} \code{data} (the training set) and an \code{sf}
##' \code{data_test} (the test set) - the structure built internally by
##' \code{assess_prediction()} for \code{method = "user"} and
##' \code{method = "regularized"}.
##' @param alpha Point transparency, passed to \code{ggplot2::geom_sf()}.
##' Defaults to \code{0.6}, matching \code{spatialsample::autoplot()}.
##' @return A \code{ggplot} object.
##' @noRd
plot_folds <- function(data_split, alpha = 0.6) {
  n_iter <- length(data_split$splits)

  combined <- do.call(rbind, lapply(seq_len(n_iter), function(i) {
    split_i <- data_split$splits[[i]]
    rbind(
      cbind(split_i$data,      set = "Training", iteration = i),
      cbind(split_i$data_test, set = "Testing",  iteration = i)
    )
  }))

  p <- ggplot(data = combined, aes(color = .data$set, fill = .data$set)) +
    geom_sf(alpha = alpha) +
    guides(colour = guide_legend("Set"), fill = guide_legend("Set")) +
    theme_minimal()

  if (n_iter > 1) {
    p <- p + facet_wrap(vars(.data$iteration))
  }

  p + coord_sf()
}

##' @title Assess Predictive Performance via Spatial Cross-Validation
##'
##' @description
##' This function evaluates the predictive performance of spatial models fitted
##' to `RiskMap` objects using cross-validation. It supports two classes of diagnostic tools:
##'
##' - **Scoring rules**, including the Continuous Ranked Probability Score (CRPS)
##'  and its scaled version (SCRPS), which quantify the sharpness and calibration
##'  of probabilistic forecasts (Bolin & Wallin, 2023);
##' - **Calibration diagnostics**, based on the Probability Integral Transform (PIT)
##' for Gaussian outcomes, Average nonrandomized PIT (AnPIT) curves for discrete
##' outcomes (e.g., Poisson or Binomial), and the area between the PIT/AnPIT curve
##' and the reference line (Giorgi *et al.* 2026).
##'
##' Cross-validation can be performed using either spatial clustering, regularized
##' subsampling with a minimum inter-point distance or a user-defined test set.
##' For each fold or subset, models can be refitted or evaluated with fixed parameters,
##' offering flexibility in model validation. The function also provides visualizations
##' of the spatial distribution of test folds.
##'
##' @param object A list of `RiskMap` objects, each representing a model fitted with `glgpm()`.
##' @param method Character; either `"cluster"`, `"regularized"` or `"user"` for the
##' cross-validation method:
##' \describe{
##'  - The `"cluster"` method uses spatial clustering as implemented
##' by the \code{spatial_clustering_cv} function from the `spatialsample` package.
##'  - The `"regularized"` method selects a subsample of the dataset by imposing a minimum distance,
##'  set by the `min_dist` argument, for a randomly selected subset of locations using
##'  the `subsample.distance` function from the `spatialEco` package.
##'  - The `"user"` method takes a user-defined test set defined by `user_split`
##'  }
##' @param fold Integer; required when `method = "cluster"` - number of folds for cross-validation.
##' @param min_dist Numeric; required when `method = "regularized"` - minimum distance in kilometers for regularized subsampling.
##' @param size Integer; the size of the test set, required when `method = "regularized"`.
##' @param user_split Required when `method = "user"`. A user-defined cross-validation split. Either:
##'   * a matrix with \code{nrow = n} (number of observations) and
##'     \code{ncol = iter} (number of iterations), where entries of \code{1}
##'     indicate membership in the test set for that iteration and \code{0}
##'     indicate training set; or
##'   * a list of length \code{iter}, where each element is either a vector of
##'     indices of the dataset to use as the test set, or a list with components
##'     \code{in_id} (training indices) and \code{out_id} (test indices).
##' @param iter Integer; number of times to repeat the cross-validation. Defaults to `1`.
##' @param metrics Character vector; one or more of `"CRPS"`, `"SCRPS"` and `"PIT"`,
##' to specify the predictive performance metrics to compute. `"PIT"` computes
##' the PIT for Gaussian models or the AnPIT curve for discrete models, and
##' also the scalar score `"PIT_area"`, the integrated absolute deviation
##' between the PIT/AnPIT curve and the reference line. Defaults to all metrics.
##' @param keep_par_fixed Logical. If `TRUE` (the default), the parameters are
##' fixed at the estimates of the fitted model, which used all the data, and
##' only the predictions are conditioned on the training set. If `FALSE`, each
##' model is re-estimated on each training set, which costs one fit per test
##' set and model. Fixed parameters are much faster, and with many locations
##' and few parameters the scores change little. The test observations still
##' contribute to the fixed estimates, however, so the scores are optimistic,
##' and more so for models with more parameters. Use `FALSE` to compare models
##' with different numbers of parameters. The choice is stored in the result
##' and printed by `summary()`.
##' @param control_mcmc Control settings for simulation, an output from `set_control_mcmc()`.
##' @param seed Optional non-negative integer that makes the assessment
##' reproducible. It sets the random splits for `method = "cluster"` and
##' `method = "regularized"`, and a seed for each test set, which is passed to
##' [glgpm()] when the models are refitted and to [predict.RiskMap()]. All
##' models share the seed of a test set, so that their predictions use common
##' random numbers. The caller's random number state is restored on exit.
##' When `NULL` (the default), the seeds are drawn from the session's random
##' number generator.
##' @param plot_fold Logical; whether to plot each iteration's test and training sets. Defaults to `TRUE`.
##' @param messages Logical, whether to report on the assessment, with a
##'   progress bar as described in [riskmap_progress].
##'   Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.
##' @param ... Additional arguments passed to clustering or subsampling functions.
##'
##' @return A list of class `RiskMap_cross_validation`, containing:
##' \describe{
##'   \item{test_set}{A list of test sets used for validation, each of class `'sf'`.}
##'   \item{model}{A named list, one per model, each containing:
##'     \describe{
##'       \item{metric}{A list with CRPS, SCRPS, and/or PIT_area metrics for each fold if requested.}
##'       \item{PIT}{(if `family = "gaussian"` and `metrics` includes `"PIT"`) A list of PIT values for test data.}
##'       \item{AnPIT}{(if `family` is discrete and `metrics` includes `"PIT"`) A list of AnPIT curves for test data.}
##'     }
##'   }
##'   \item{keep_par_fixed}{The value of `keep_par_fixed`.}
##'   \item{failures}{A data frame with one row for each model and test set
##'   whose refit or prediction failed, with columns `model`, `test_set` and
##'   `message`. The scores of these test sets are `NA`.}
##' }
##'
##' @details
##' A failed refit or prediction does not stop the assessment: the scores of
##' that test set are missing and a warning reports the failures. Test
##' observations with missing values in a variable used by any of the models
##' are left out. The predictive distributions of the test observations do not
##' yet include random effects, so the scores of models with `re()` terms are
##' not reliable, and a warning says so.
##'
##' @seealso \code{\link{plot.RiskMap_cross_validation}}
##'
##' @references
##' Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of
##' proper scoring rules. *Statistical Science*, 38(1), 140–159. \doi{10.1214/22-STS864}.
##'
##' Giorgi, E., Fronterre, C. & Diggle, P. J. (2026). A decay-adjusted spatio-temporal
##' model to account for the impact of mass drug administration on neglected
##' tropical disease prevalence. *Journal of the Royal Statistical Society Series
##' A: Statistics in Society*.\doi{10.1093/jrsssa/qnag100}.
##'
##' @importFrom spatialEco subsample.distance
##' @importFrom spatialsample spatial_clustering_cv autoplot
##'
##' @examples
##'
##' data(italy_sim)
##'
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' # cluster method
##' cross_validation <-
##'   assess_prediction(
##'     list(fit),
##'     method = "cluster",
##'     fold = 2
##'   )
##'
##' summary(cross_validation)
##'
##' # regularized method
##' cross_validation <-
##'   assess_prediction(
##'     list(fit),
##'     method = "regularized",
##'     size = 5,
##'     min_dist = 1
##'   )
##'
##' summary(cross_validation)
##'
##' # user method with a matrix
##'  cross_validation <-
##'   assess_prediction(
##'   list(fit),
##'   method = "user",
##'   user_split = matrix(
##'     sample(c(rep(1, 50), rep(0, 50))),
##'     ncol = 1)
##'   )
##'
##' summary(cross_validation)
##'
##' # user method with a list
##'  cross_validation <-
##'   assess_prediction(
##'   list(fit),
##'   method = "user",
##'   user_split = list(
##'     sample(100, 50))
##'   )
##'
##' summary(cross_validation)
##'
##' @export
assess_prediction <- function(object,
                              method,
                              fold = NULL,
                              min_dist = NULL,
                              size = NULL,
                              user_split = NULL,
                              iter = 1,
                              metrics = c("PIT", "CRPS", "SCRPS"),
                              keep_par_fixed = TRUE,
                              control_mcmc = set_control_mcmc(),
                              seed = NULL,
                              plot_fold = TRUE,
                              messages = getOption("RiskMap.messages", TRUE),
                              ...) {

  ## ─────────────────────────── helpers ─────────────────────────── ##
  is_list_of_riskmap <- function(x) {
    is.list(x) && all(vapply(x, inherits, logical(1), what = "RiskMap"))
  }

  crps_gaussian <- function(y, mu, sigma) {
    if (sigma == 0) return(0)
    z <- (y - mu) / sigma
    2 * dnorm(z) + z * (2 * pnorm(z) - 1) - 1 / sqrt(pi)
  }
  crps_discrete <- function(y, Fk) {
    k <- seq_along(Fk) - 1
    sum((Fk - as.numeric(k >= y))^2)
  }
  # E_{Y ~ F} CRPS(F, Y) = sum_k F(k) (1 - F(k)), which equals E|X - X'| / 2
  exp_crps_discrete <- function(Fk) sum(Fk * (1 - Fk))
  # SCRPS of Bolin and Wallin (2023), undefined for a degenerate predictive
  scrps <- function(crps, exp_crps) {
    if (exp_crps > 0) -0.5 * (1 + crps / exp_crps + log(2 * exp_crps)) else NA_real_
  }
  u_val <- seq(0, 1, length.out = 1000)
  # Nonrandomized PIT curve of one discrete observation (Czado et al., 2009)
  npit_curve <- function(y, Fk) {
    F_below <- if (y == 0) 0 else Fk[y]
    F_at <- Fk[y + 1]
    if (F_at > F_below) {
      pmin(pmax((u_val - F_below) / (F_at - F_below), 0), 1)
    } else {
      as.numeric(u_val > F_at)
    }
  }

  ## ────────────────────── sanity checks (unchanged) ─────────────────────── ##
  if (!is_list_of_riskmap(object)) {
    abort_input("{.arg object} must be a list of {.cls RiskMap} fits.", arg = "object")
  }
  allowed_metrics <- c("CRPS", "SCRPS", "PIT")
  if (!is.null(metrics) &&
      (!is.character(metrics) || !all(metrics %in% allowed_metrics))) {
    abort_input("{.arg metrics} must contain only {.val {allowed_metrics}}.", arg = "metrics")
  }
  method <- match_choice(method, c("cluster", "regularized", "user"))
  if (method == "cluster") {
    if (is.null(fold)) {
      abort_input("{.arg fold} must be given when {.arg method} is {.val cluster}.", arg = "fold")
    }
    check_positive_integer(fold, "fold")
  }
  if (method == "regularized") {
    if (is.null(min_dist)) {
      abort_input("{.arg min_dist} must be given when {.arg method} is {.val regularized}.", arg = "min_dist")
    }
    if (is.null(size)) {
      abort_input("{.arg size} must be given when {.arg method} is {.val regularized}.", arg = "size")
    }
    check_positive_number(min_dist, "")
    check_positive_integer(size, "size")
  }
  if (method == "user" && is.null(user_split)) {
    abort_input("{.arg user_split} must be given when {.arg method} is {.val user}.", arg = "user_split")
  }

  check_logical(keep_par_fixed)

  check_positive_integer(iter, "iter")

  check_control(control_mcmc, "set_control_mcmc")

  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)

  check_logical(plot_fold)

  check_logical(messages)

  get_CRPS  <- "CRPS"  %in% metrics
  get_SCRPS <- "SCRPS" %in% metrics
  get_PIT <- "PIT" %in% metrics

  ## ───────────────────────────── data & splits ───────────────────────────── ##

  # add names to any unnamed models
  if (is.null(names(object))) {
    names(object) <- paste("Model", seq_along(object))
  } else {
    unnamed <- which(names(object) == "")
    if (length(unnamed) > 0) {
      names(object)[unnamed] <- paste("Model", unnamed)
    }
  }

  object1 <- object[[1]]
  data_sf <- object1$data
  n_obs   <- nrow(data_sf)
  data_geom <- st_geometry(data_sf)

  for (h in seq_along(object)) {
    fit_data <- object[[h]]$data
    if (nrow(fit_data) != n_obs) {
      abort_input("The models in {.arg object} must be fitted to the same number of observations.", arg = "object")
    }
    fit_geom <- st_geometry(fit_data)
    if (!identical(fit_geom, data_geom)) {
      abort_input("The models in {.arg object} must be fitted to data with the same rows in the same order.", arg = "object")
    }
  }

  # Validates one vector of row indices for a 'user_split' entry: non-missing
  # numeric, non-empty, whole numbers, no duplicates, in range. `what` names
  # whichever vector the caller actually supplied, so it reads sensibly for
  # both call sites below (a bare vector of test indices, or an explicit
  # 'in_id'/'out_id').
  validate_split_vector <- function(v, what) {
    if (!is.numeric(v) || anyNA(v)) {
      abort_input("{what} must be numeric row numbers without missing values.", arg = "user_split")
    }
    if (length(v) == 0) abort_input("{what} must not be empty.", arg = "user_split")
    if (any(v != round(v))) {
      abort_input("{what} must be whole numbers.", arg = "user_split")
    }
    if (anyDuplicated(v)) {
      abort_input("{what} must not repeat a row.", arg = "user_split")
    }
    if (!all(v %in% seq_len(n_obs))) {
      abort_input("{what} must be row numbers between 1 and {n_obs}.", arg = "user_split")
    }
    invisible(TRUE)
  }

  # Validates one user-supplied 'user_split' entry: either a bare vector of
  # test indices (the shorthand form - training indices are its computed
  # complement, which can never independently be invalid) or an explicit
  # list(in_id = ..., out_id = ...) pair, both directly supplied by the
  # caller. For the latter, also checks train and test don't overlap - an
  # overlap would silently leak a held-out observation into training.
  validate_split_indices <- function(out_id, index, in_id = NULL) {
    tag <- sprintf("user_split[[%d]]", index)

    if (is.null(in_id)) {
      validate_split_vector(out_id, cli::format_inline("The test rows in {.code {tag}}"))
      in_id <- setdiff(seq_len(n_obs), out_id)
      if (length(in_id) == 0) {
        abort_input("{.code {tag}} must leave at least one observation for training.", arg = "user_split")
      }
    } else {
      validate_split_vector(in_id, cli::format_inline("The {.code in_id} of {.code {tag}}"))
      validate_split_vector(out_id, cli::format_inline("The {.code out_id} of {.code {tag}}"))
      if (length(intersect(in_id, out_id)) > 0) {
        abort_input("The {.code in_id} and {.code out_id} of {.code {tag}} must not overlap.", arg = "user_split")
      }
    }

    list(in_id = as.integer(in_id), out_id = as.integer(out_id))
  }

  make_splits_from_user <- function(usr, n_iter_expected) {
    spl <- vector("list", n_iter_expected)
    if (is.matrix(usr)) {
      if (nrow(usr) != n_obs) {
        abort_input("A matrix {.arg user_split} must have one row per observation ({n_obs}), not {nrow(usr)}.", arg = "user_split")
      }
      if (ncol(usr) != n_iter_expected) {
        abort_input("A matrix {.arg user_split} must have {.arg iter} ({n_iter_expected}) columns, not {ncol(usr)}.", arg = "user_split")
      }
      if (anyNA(usr)) {
        abort_input("A matrix {.arg user_split} must not have missing values.", arg = "user_split")
      }
      if (!all(usr %in% c(0, 1))) {
        abort_input("A matrix {.arg user_split} must hold only 0 for training and 1 for test.", arg = "user_split")
      }
      for (i in seq_len(n_iter_expected)) {
        out_id <- which(usr[, i] == 1)
        in_id  <- which(usr[, i] == 0)
        if (length(in_id) == 0 || length(out_id) == 0) {
          abort_input("Column {i} of {.arg user_split} must have at least one training observation (0) and one test observation (1).", arg = "user_split")
        }
        spl[[i]] <- list(in_id = in_id, out_id = out_id,
                         data = data_sf[in_id, ],
                         data_test = data_sf[out_id, ])
      }
    } else if (inherits(usr, "list")) {
      if (length(usr) != n_iter_expected) {
        abort_input("A list {.arg user_split} must have {.arg iter} ({n_iter_expected}) elements, not {length(usr)}.", arg = "user_split")
      }
      for (i in seq_len(n_iter_expected)) {
        ui <- usr[[i]]
        if (is.list(ui) && !is.null(ui$in_id) && !is.null(ui$out_id)) {
          validated <- validate_split_indices(ui$out_id, i, ui$in_id)
          in_id  <- validated$in_id
          out_id <- validated$out_id
        } else if (is.numeric(ui)) {
          validated <- validate_split_indices(ui, i)
          in_id  <- validated$in_id
          out_id <- validated$out_id
        } else {
          abort_input("Each element of {.arg user_split} must be a vector of test rows or a list with {.code in_id} and {.code out_id}.", arg = "user_split")
        }
        spl[[i]] <- list(in_id = in_id, out_id = out_id,
                         data = data_sf[in_id, ],
                         data_test = data_sf[out_id, ])
      }
    } else {
      abort_input("{.arg user_split} must be a matrix or a list, not {.obj_type_friendly {usr}}.", arg = "user_split")
    }
    list(splits = spl)
  }

  if (!is.null(seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(seed)
  }

  if (method == "user") {
    data_split <- make_splits_from_user(user_split, iter)
    n_iter <- iter

    if (isTRUE(plot_fold)) print(plot_folds(data_split))
  } else if (method == "cluster") {
    data_split <- spatial_clustering_cv(data = data_sf, v = fold, repeats = iter, ...)
    # ensure out_id present
    all_ids <- seq_len(nrow(data_sf))
    for (i in seq_along(data_split$splits)) {
      if (is.null(data_split$splits[[i]]$out_id) || all(is.na(data_split$splits[[i]]$out_id))) {
        in_id <- data_split$splits[[i]]$in_id
        data_split$splits[[i]]$out_id <- setdiff(all_ids, in_id)
      }
    }
    n_iter <- iter * fold
    if (plot_fold) print(autoplot(data_split))
  } else { # regularized
    data_split <- list(splits = vector("list", iter))
    for (i in seq_len(iter)) {
      locations_sf <- data_sf[!duplicated(st_as_text(data_sf$geometry)), ]
      data_split$splits[[i]] <- list()
      data_split$splits[[i]]$data_test <- withCallingHandlers(
        subsample.distance(locations_sf, size = size, d = min_dist * 1000, ...),
        message = function(m) if (!messages) invokeRestart("muffleMessage")
      )
      test_geom <- st_as_text(data_split$splits[[i]]$data_test$geometry)
      in_test   <- st_as_text(data_sf$geometry) %in% test_geom
      data_split$splits[[i]]$out_id <- which(in_test)
      data_split$splits[[i]]$in_id  <- which(!in_test)
      data_split$splits[[i]]$data   <- data_sf[!in_test, ]
    }
    n_iter <- iter
    if (isTRUE(plot_fold)) print(plot_folds(data_split))
  }

  ## ─────────────────────── checks on splits and models ─────────────────────── ##
  if (method == "regularized") {
    n_short <- sum(vapply(data_split$splits,
                          function(s) nrow(s$data_test) < size, logical(1)))
    if (n_short > 0) {
      riskmap_warn("{n_short} of {n_iter} test set{cli::qty(n_short)}{?s} {?has/have} fewer than {.arg size} = {size} locations, because no more locations are at least {.arg min_dist} apart.")
    }
  }

  # A held-out observation at a location that also has training observations
  # is predicted with a new nugget draw, ignoring the nugget seen in training.
  has_nugget <- vapply(object, function(fit) fitted_nugget_variance(fit) > 0,
                       logical(1))
  if (method == "user" && any(has_nugget)) {
    ids <- object1$ID_coords
    shared <- which(vapply(data_split$splits,
                           function(s) any(ids[s$out_id] %in% ids[s$in_id]),
                           logical(1)))
    if (length(shared) > 0) {
      riskmap_warn(c("In {cli::qty(length(shared))}test set{?s} {shared}, some locations also have training observations.", "i" = "Their predictions ignore the nugget seen in training, so the scores of models with a nugget aren't reliable there."))
    }
  }

  with_re <- names(object)[vapply(object, function(fit) !is.null(fit$ID_re),
                                  logical(1))]
  if (length(with_re) > 0) {
    riskmap_warn("The predictive distributions of the test observations don't include the random effects, so the scores of {.val {with_re}} aren't reliable.")
  }

  n_par <- vapply(object, function(fit) length(unlist(coef(fit))), integer(1))
  if (messages) {
    if (keep_par_fixed) {
      riskmap_inform(c("Parameters are fixed at the estimates from all the data.", "i" = "Set {.code keep_par_fixed = FALSE} to re-estimate them in each training set."))
      if (length(unique(n_par)) > 1) {
        riskmap_inform("The models have different numbers of parameters. With fixed parameters the comparison tends to favour the model with more of them, so {.code keep_par_fixed = FALSE} is recommended.")
      }
    } else {
      n_fits <- n_iter * length(object)
      riskmap_inform("Parameters are re-estimated in each training set ({n_fits} model fits).")
    }
  }

  # Test observations are scored only if every model can use them
  complete_rows <- Reduce(`&`, lapply(object, function(fit) {
    model_vars <- intersect(c(all.vars(fit$formula),
                              as.character(fit$call$denominator)),
                            names(fit$data))
    complete.cases(st_drop_geometry(fit$data)[, model_vars, drop = FALSE])
  }))
  if (!all(complete_rows)) {
    riskmap_warn("{sum(!complete_rows)} observation{?s} with missing values in the model variables {?is/are} not used as test observations.")
  }

  # One seed per test set, shared by the models so that their predictions use
  # common random numbers. The measurement error of Gaussian models is drawn
  # from a separate stream, independent of the spatial draws.
  fold_seeds  <- sample.int(.Machine$integer.max, n_iter)
  noise_seeds <- sample.int(.Machine$integer.max, n_iter)

  ## ───────────────────────── initialise output ───────────────────────── ##
  n_models    <- length(object)
  model_names <- names(object)
  out <- list(test_set = vector("list", n_iter), model = list(),
              keep_par_fixed = keep_par_fixed)
  failures <- data.frame(model = character(), test_set = integer(),
                         message = character())

  for (i in seq_len(n_iter)) {
    out_id <- data_split$splits[[i]]$out_id
    out$test_set[[i]] <- data_sf[out_id[complete_rows[out_id]], ]
  }

  bar <- progress_bar(if (keep_par_fixed) "Scoring test sets" else
                        "Refitting and scoring test sets",
                      n_models * n_iter, messages)

  ## ─────────────────────────── iterate over models ───────────────────────── ##
  for (h in seq_len(n_models)) {

    fit0      <- object[[h]]
    fit_data_sf <- fit0$data
    den_name  <- as.character(fit0$call$denominator)
    traits    <- family_traits(fit0$family)
    # A Gaussian fit stores no link, so its identity comes from the family.
    linkfun   <- fit0$link_function$inv %||% traits$inverse_link

    ## containers for this model
    if (get_CRPS)   CRPS  <- vector("list", n_iter)
    if (get_SCRPS) { y_CRPS <- vector("list", n_iter); SCRPS <- vector("list", n_iter) }
    if (get_PIT) {
      PIT_area <- vector("list", n_iter)
      if (!traits$discrete) PIT <- vector("list", n_iter) else AnPIT <- vector("list", n_iter)
    }

    ## Scores of one test set: refit or slice, predict, score
    score_test_set <- function(i, in_id, out_id_i, data_test_i) {
      ## ----- refit or slice -----
      if (!keep_par_fixed) {

        refit_args <- list(
          formula        = fit0$formula,
          data           = fit_data_sf[in_id, ],
          family         = fit0$family,
          model_crs      = st_crs(fit0$data),
          distance_units = fit0$distance_units,
          control_mcmc   = control_mcmc,
          control_mcml   = attr(fit0, "control_mcml") %||%
            set_control_mcml(),
          seed           = fold_seeds[i],
          fix_var_me     = fit0$fix_var_me,
          return_samples = FALSE,
          messages       = FALSE
        )
        ## 'denominator' is passed as an unquoted column name (NSE); include it
        ## when the original model was fitted with one
        if (length(den_name) == 1 && nzchar(den_name)) {
          refit_args$denominator <- as.name(den_name)
        }
        refit_args$invlink <- refit_invlink(fit0)
        refit_i <- do.call(glgpm, refit_args)
      } else {
        ## quick slice without re-fitting
        refit_i <- fit0
        keep <- in_id
        refit_i$data  <- refit_i$data [keep, ]
        refit_i$units_m  <- refit_i$units_m[keep]
        keep_coord <- unique(refit_i$ID_coords[keep])
        refit_i$coords   <- refit_i$coords[keep_coord, , drop = FALSE]
        refit_i$y        <- refit_i$y      [keep]
        refit_i$D        <- refit_i$D      [keep, , drop = FALSE]
        if (!is.null(refit_i$cov_offset))
          refit_i$cov_offset <- refit_i$cov_offset[keep]
        if (!is.null(refit_i$ID_re)) {
          re_terms <- names(refit_i$ID_re)
          random_effects_i <- prepare_random_effects(refit_i$data, re_terms)
          refit_i$ID_re <- as.data.frame(random_effects_i$ID_re)
          colnames(refit_i$ID_re) <- random_effects_i$names_re
          refit_i$re <- random_effects_i$re_unique_f
        }
        ## recompute ID_coords mapping
        refit_i$ID_coords <- create_ids(refit_i$data)$ID_coords
      }

      ## ----- prediction over test set -----
      included_components <- if (
        is.null(refit_i$fix_tau2) || refit_i$fix_tau2 != 0
      ) "nugget" else character()
      prediction_i <- predict(
        refit_i,
        newdata = data_test_i,
        type = "marginal",
        include = included_components,
        summaries = "mean",
        control_mcmc = control_mcmc,
        seed = fold_seeds[i],
        keep_samples = TRUE,
        messages = FALSE
      )
      eta_samp <- prediction_link_samples(prediction_i)
      if (traits$measurement_error) {
        sigma2_me <- if (is.null(refit_i$fix_var_me)) coef(refit_i)$sigma2_me else refit_i$fix_var_me
        restore_noise_seed <- preserve_random_seed()
        set.seed(noise_seeds[i])
        noise <- rnorm(length(eta_samp))
        restore_noise_seed()
        eta_samp <- eta_samp + sqrt(sigma2_me) * noise
      }
      mu_samp <- linkfun(eta_samp)

      n_pred <- nrow(eta_samp)
      n_draw <- ncol(eta_samp)
      crps_i   <- numeric(n_pred)
      y_crps_i <- numeric(n_pred)
      scrps_i  <- numeric(n_pred)
      if (!traits$discrete) {
        pit_i <- numeric(n_pred)
      } else {
        anpit_i <- matrix(NA_real_, nrow = length(u_val), ncol = n_pred)
      }

      units_m_i <- fit0$units_m[out_id_i]
      y_i       <- fit0$y      [out_id_i]

      for (j in seq_len(n_pred)) {
        if (!traits$discrete) {
          mu_j <- mean(mu_samp[j, ])
          sd_j <- sd(mu_samp[j, ])
          if (sd_j > 0) {
            crps_i[j]   <- sd_j * crps_gaussian(y_i[j], mu_j, sd_j)
            y_crps_i[j] <- sd_j / sqrt(pi)
            pit_i[j]    <- pnorm(y_i[j], mean = mu_j, sd = sd_j)
          } else {
            # a point mass: the CRPS is the absolute error
            crps_i[j]   <- abs(y_i[j] - mu_j)
            y_crps_i[j] <- 0
            pit_i[j]    <- as.numeric(y_i[j] >= mu_j)
          }
        } else {
          y_samp  <- traits$simulate(mu_samp[j, ], units_m_i[j], NULL)
          support <- traits$score_support(y_samp, y_i[j], units_m_i[j],
                                          mu_samp[j, ])
          pk <- tabulate(y_samp + 1, nbins = length(support)) / n_draw
          Fk <- cumsum(pk)
          crps_i[j]   <- crps_discrete(y_i[j], Fk)
          y_crps_i[j] <- exp_crps_discrete(Fk)
          anpit_i[, j] <- npit_curve(y_i[j], Fk)
        }
        scrps_i[j] <- scrps(crps_i[j], y_crps_i[j])
      }

      if (!traits$discrete) {
        list(CRPS = crps_i, y_CRPS = y_crps_i, SCRPS = scrps_i, PIT = pit_i,
             PIT_area = .anpit_area(ecdf(pit_i)(u_val), u_val))
      } else {
        anpit <- rowMeans(anpit_i)
        list(CRPS = crps_i, y_CRPS = y_crps_i, SCRPS = scrps_i, AnPIT = anpit,
             PIT_area = .anpit_area(anpit, u_val))
      }
    }

    ## ─────────────────── iterate over CV splits (refit as needed) ─────────────────── ##
    for (i in seq_len(n_iter)) {

      in_id  <- data_split$splits[[i]]$in_id
      out_id <- data_split$splits[[i]]$out_id
      out_id_i <- out_id[complete_rows[out_id]]
      data_test_i <- fit_data_sf[out_id_i, ]

      # A failed fit or prediction leaves this test set's scores missing
      # instead of stopping the whole assessment
      scores_i <- tryCatch(
        score_test_set(i, in_id, out_id_i, data_test_i),
        error = function(e) {
          failures[nrow(failures) + 1, ] <<- list(model_names[h], i,
                                                  condition_text(e))
          n_test <- length(out_id_i)
          list(CRPS = rep(NA_real_, n_test), y_CRPS = rep(NA_real_, n_test),
               SCRPS = rep(NA_real_, n_test), PIT = numeric(0),
               AnPIT = numeric(0), PIT_area = NA_real_)
        }
      )

      if (get_CRPS) CRPS[[i]] <- scores_i$CRPS
      if (get_SCRPS) {
        y_CRPS[[i]] <- scores_i$y_CRPS
        SCRPS[[i]]  <- scores_i$SCRPS
      }
      if (get_PIT) {
        PIT_area[[i]] <- scores_i$PIT_area
        if (!traits$discrete) PIT[[i]] <- scores_i$PIT else AnPIT[[i]] <- scores_i$AnPIT
      }
      progress_tick(bar)

    } # end i loop

    ## ─────────────── finalise output for this model ─────────────── ##
    out$model[[model_names[h]]] <- list(metric = list())
    if (get_CRPS)  out$model[[model_names[h]]]$metric$CRPS  <- CRPS
    if (get_SCRPS) out$model[[model_names[h]]]$metric$SCRPS <- SCRPS
    if (get_PIT) {
      out$model[[model_names[h]]]$metric$PIT_area <- PIT_area
      if (!traits$discrete) out$model[[model_names[h]]]$PIT <- PIT else out$model[[model_names[h]]]$AnPIT <- AnPIT
    }

  } # end h loop
  progress_done(bar)

  out$failures <- failures
  if (nrow(failures) > 0) {
    n_assessed <- n_iter * n_models
    first <- failures[1, ]
    riskmap_warn(c("{nrow(failures)} of {n_assessed} test-set assessments failed, and their scores are missing.", "i" = "The first failure, for model {.val {first$model}} and test set {first$test_set}, was {.val {first$message}}.", "i" = "See the {.code failures} element of the result."))
  }

  class(out) <- "RiskMap_cross_validation"
  return(out)
}

##' @title Plot Method for RiskMap_simulation Objects
##' @description
##' Maps the linear predictor of one simulation from [simulate_glgpm()] over
##' its prediction grid, with the sampling locations overlaid when data were
##' also simulated.
##'
##' @param x Output from [simulate_glgpm()] including a surface.
##' @param simulation A single simulation number to plot, defaulting to one.
##' @inheritParams plot.RiskMap_simulated_surface
##'
##' @return A \code{ggplot} object.
##' @details
##' The returned object can be modified further
##' with standard \pkg{ggplot2} functions.
##'
##' @seealso [simulated_surface()], [plot.RiskMap_simulated_surface()]
##' @method plot RiskMap_simulation
##' @export
plot.RiskMap_simulation <- function(x, simulation = 1, palette = "viridis",
                                    reverse_palette = FALSE, ...) {
  out <- plot(simulated_surface(x, simulation),
              palette = palette,
              reverse_palette = reverse_palette,
              ...) +
    ggtitle(paste("Simulation", simulation))

  if (!is.null(x$locations$data)) {
    # geom_point() rather than geom_sf(), which would replace the map's coord_sf()
    sample_coordinates <- st_coordinates(st_transform(x$locations$data,
                                                      st_crs(x$locations$surface)))
    out <- out + geom_point(data = data.frame(x = sample_coordinates[, 1],
                                              y = sample_coordinates[, 2]),
                            aes(x = .data$x, y = .data$y),
                            inherit.aes = FALSE, size = 1)
  }
  out
}

##' Map values over a regular grid of points
##'
##' Draws each grid point as a tile, styled with `add_map_layers()`. Shared by
##' the plot methods for grid predictions and simulated surfaces.
##'
##' @param grid An `sf` or `sfc` object of points on a regular grid.
##' @param values Numeric vector of values to map, one per point in `grid`.
##' @param legend_title Character string used as the fill legend title.
##' @param palette,reverse_palette,limits,north_arrow,scale_bar,... Passed to
##' `add_map_layers()`.
##' @param call The call named in input errors, so that they name the plot
##' method.
##' @return A `ggplot` object.
##' @noRd
plot_grid <- function(grid, values, legend_title, palette, reverse_palette,
                      limits, north_arrow, scale_bar, ...,
                      call = rlang::caller_env()) {
  grid_coordinates <- st_coordinates(grid)
  plot_data <- data.frame(x = grid_coordinates[, 1],
                          y = grid_coordinates[, 2],
                          value = values)

  tile_width <- grid_spacing(plot_data$x)
  tile_height <- grid_spacing(plot_data$y)
  # a single row or column has no spacing of its own, so borrow the other
  # axis's to keep the tiles square
  if (is.na(tile_width) && is.na(tile_height)) {
    tile_width <- tile_height <- 1
  } else if (is.na(tile_width)) {
    tile_width <- tile_height
  } else if (is.na(tile_height)) {
    tile_height <- tile_width
  }

  # extent of the tiles, not just their centres
  extent <- c(xmin = min(plot_data$x) - tile_width / 2,
              xmax = max(plot_data$x) + tile_width / 2,
              ymin = min(plot_data$y) - tile_height / 2,
              ymax = max(plot_data$y) + tile_height / 2)

  out <- ggplot(plot_data) +
    geom_tile(aes(x = .data$x, y = .data$y, fill = .data$value),
              width = tile_width, height = tile_height)

  add_map_layers(out,
                 extent = extent,
                 crs = st_crs(grid),
                 legend_title = legend_title,
                 palette = palette,
                 reverse_palette = reverse_palette,
                 limits = limits,
                 north_arrow = north_arrow,
                 scale_bar = scale_bar,
                 ...,
                 call = call)
}

##' Spacing between the points of a regular grid along one axis
##'
##' @param coordinates Numeric vector of coordinates along one axis.
##' @return The smallest gap between distinct coordinates, or `NA` when they
##' all share one value.
##' @noRd
grid_spacing <- function(coordinates) {
  if (length(unique(coordinates)) < 2) return(NA_real_)
  resolution(coordinates, zero = FALSE)
}

##' Convert a palette specification into a vector of colours
##'
##' @param palette Either the name of a palette in [grDevices::hcl.pals()] or a
##' character vector of two or more valid colours.
##' @param reverse_palette Logical; reverse the order of the colours.
##' @param call The call named in input errors.
##' @return A character vector of colours.
##' @noRd
palette_colours <- function(palette, reverse_palette = FALSE,
                            call = rlang::caller_env()) {
  if (!is.character(palette)) {
    abort_input("{.arg palette} must be a character vector, not {.obj_type_friendly {palette}}.", arg = "palette", call = call)
  }
  check_logical(reverse_palette, call = call)

  if (length(palette) == 1) {
    # hcl.colors() ignores case, spaces and punctuation in palette names
    normalise_name <- function(name) tolower(gsub("[^[:alnum:]]", "", name))
    if (!normalise_name(palette) %in% normalise_name(grDevices::hcl.pals())) {
      abort_input("{.arg palette} must be the name of a palette in {.fn grDevices::hcl.pals} or two or more colours, not {.val {palette}}.", arg = "palette", call = call)
    }
    colours <- grDevices::hcl.colors(256, palette = palette)
  } else {
    valid_colour <- vapply(palette, function(colour) {
      !inherits(try(grDevices::col2rgb(colour), silent = TRUE), "try-error")
    }, logical(1))
    if (!all(valid_colour)) {
      invalid <- palette[!valid_colour]
      abort_input("{.arg palette} has {?an /}invalid colour{?s}, {.val {invalid}}.", arg = "palette", call = call)
    }
    colours <- palette
  }
  if (reverse_palette) colours <- rev(colours)
  colours
}

##' Add the fill scale, map extent, north arrow and scale bar shared by the
##' prediction maps
##'
##' @param map A `ggplot` object.
##' @param extent Named numeric vector (`xmin`, `xmax`, `ymin`, `ymax`) giving
##' the extent of the data in `crs`.
##' @param crs The coordinate reference system of the data.
##' @param legend_title Character string used as the fill legend title.
##' @param palette,reverse_palette Passed to `palette_colours()`.
##' @param limits `NULL` or a numeric vector of length two for the colour scale.
##' @param north_arrow,scale_bar Logical; whether to add each annotation.
##' @param ... Additional arguments passed to [ggplot2::scale_fill_gradientn()].
##' @param call The call named in input errors.
##' @return The `ggplot` object with the layers added.
##' @noRd
add_map_layers <- function(map, extent, crs, legend_title, palette,
                           reverse_palette, limits, north_arrow, scale_bar, ...,
                           call = rlang::caller_env()) {
  check_logical(north_arrow, call = call)
  check_logical(scale_bar, call = call)
  map <- map +
    scale_fill_gradientn(colours = palette_colours(palette, reverse_palette,
                                                   call = call),
                         limits = limits, ...) +
    coord_sf(crs = crs,
             xlim = as.numeric(extent[c("xmin", "xmax")]),
             ylim = pad_map_extent(extent, north_arrow, scale_bar)) +
    labs(x = NULL, y = NULL, fill = legend_title)

  if (north_arrow) {
    map <- map +
      ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                        height = unit(1, "cm"),
                                        width = unit(1, "cm"))
  }
  if (scale_bar) {
    map <- map + ggspatial::annotation_scale(location = "bl")
  }
  map
}

##' Extend the vertical extent of a map to make room for annotations
##'
##' @param extent Named numeric vector (`xmin`, `xmax`, `ymin`, `ymax`).
##' @param north_arrow Logical; add space above the data for a north arrow.
##' @param scale_bar Logical; add space below the data for a scale bar.
##' @return Numeric vector of length two giving the padded y limits.
##' @noRd
pad_map_extent <- function(extent, north_arrow, scale_bar) {
  height <- unname(extent["ymax"] - extent["ymin"])
  # roughly the height of each annotation on a typical map
  top_padding <- if (north_arrow) 0.15 * height else 0
  bottom_padding <- if (scale_bar) 0.08 * height else 0
  unname(c(extent["ymin"] - bottom_padding, extent["ymax"] + top_padding))
}

##' @title Plot Method for RiskMap_simulated_surface Objects
##' @description
##' Maps the linear predictor of a surface extracted with
##' [simulated_surface()] over its prediction grid.
##'
##' @param x Output from [simulated_surface()].
##' @param palette Either the name of a palette from [grDevices::hcl.pals()]
##' (e.g. `"viridis"`, `"Blues"`, `"Spectral"`) or a character vector of two or
##' more colours to interpolate between. Defaults to `"viridis"`.
##' @param reverse_palette Logical; if `TRUE`, reverse the order of the palette
##' colours. Defaults to `FALSE`.
##' @param limits `NULL` (the default) or a numeric vector of length two giving
##' the range of the colour scale, for example to map several simulations on
##' a common scale.
##' @param ... Additional arguments passed to [ggplot2::scale_fill_gradientn()],
##' e.g. `name` or `breaks`.
##'
##' @return A \code{ggplot} object.
##' @details
##' The returned object can be modified further
##' with standard \pkg{ggplot2} functions.
##'
##' @seealso [simulated_surface()], [plot.RiskMap_simulation()]
##' @method plot RiskMap_simulated_surface
##' @export
plot.RiskMap_simulated_surface <- function(x, palette = "viridis",
                                           reverse_palette = FALSE,
                                           limits = NULL, ...) {
  check_class(x, "RiskMap_simulated_surface")

  plot_grid(x,
            values = x$linear_predictor,
            legend_title = "linear_predictor",
            palette = palette,
            reverse_palette = reverse_palette,
            limits = limits,
            north_arrow = FALSE,
            scale_bar = FALSE,
            ...)
}

##' @title Assess Simulations
##'
##' @description Evaluates models using joint data and surface simulations.
##'
##' @param obj_sim Output from [simulate_glgpm()] with both data and surface,
##'   from a generating model without grouped `re()` effects. The generating
##'   model can have any link.
##' @param models A named list of fitted models returned by [glgpm()] or model
##'   specifications returned by [specify_glgpm()]. Each object defines a
##'   candidate model specification and is refitted to every simulated dataset;
##'   its estimated parameter values are not reused.
##' @param control_mcmc A control object for MCMC sampling, created with `set_control_mcmc()`. Default is `set_control_mcmc()`.
##' @param seed Optional non-negative integer that makes the assessment
##'   reproducible; see Details. The caller's random number state is restored
##'   on exit. When `NULL` (the default), the seeds of the simulations are
##'   drawn from the session's random number generator.
##' @param spatial_scale The scale(s) at which predictions are assessed: `"grid"`, `"area"`, or `c("grid", "area")`
##'   to compute both from a single fit-and-predict pass over the simulations.
##' @param messages Logical, whether to report on the assessment, with a
##'   progress bar as described in [riskmap_progress].
##'   Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.
##' @param target A function that converts a numeric matrix of means to a
##'   numeric matrix of the scientific target, with the same dimensions. The
##'   mean is the probability for a binomial model, the mean per unit exposure
##'   for a Poisson model and the value itself for a Gaussian model. The
##'   default, [identity()], scores the mean itself, such as prevalence. Other
##'   targets are functions of the mean, such as [qlogis()] to score
##'   prevalence on the logit scale, or `function(mu) pnorm((-2 - mu) / sd)`
##'   for the proportion of children stunted from a Gaussian model of
##'   height-for-age z-scores with residual standard deviation `sd`. See
##'   Details.
##' @param area_summary A function that combines a numeric vector of grid-cell
##'   targets within one area and returns one finite numeric value, such as
##'   [base::mean()] or [base::sum()]. Required only when `spatial_scale`
##'   includes `"area"`.
##' @param boundaries An `sf` object containing only POLYGON or MULTIPOLYGON geometries, required if `spatial_scale` includes `"area"`.
##' @param col_names Column name in `boundaries` containing unique region names. If `NULL`, defaults to `"region"`.
##' @param pred_objective A character vector specifying objectives, either `"mse"`, `"classify"`, or both.
##' @param categories A numeric vector of thresholds defining categories for classification. Required if `pred_objective = "classify"`.
##'
##' @return A list of class `RiskMap_assess_simulation`:
##' \describe{
##'   \item{pred_objective}{One element per requested `spatial_scale`
##'   (`"grid"` and/or `"area"`), each holding `mse`, a matrix with one row per
##'   model and one column per simulation, and/or `classify`, the
##'   classification metrics of each model and simulation.}
##'   \item{fits}{A data frame with one row per model and simulation, giving
##'   the `status` of the fit (`"ok"` or `"failed"`), the error `message` of a
##'   failed fit and any `warnings`.}
##'   \item{estimates}{A list with, for each model, a matrix of the parameter
##'   estimates, one row per simulation.}
##'   \item{n_sim, spatial_scale}{The number of simulations and the spatial
##'   scales assessed.}
##' }
##'
##' @details
##' Every model is scored on the scale of the mean, each through its own
##' inverse link, so the links of the generating model and of the candidates
##' need not agree. The true target at each grid location is `target()`
##' applied to the simulated mean, which is the generating model's inverse
##' link of the simulated linear predictor. That linear predictor includes the
##' covariates, any offset and, when the generating model has one, the
##' realised nugget. Each candidate is refitted to every simulated dataset,
##' its predictive draws of the linear predictor are mapped to the mean
##' through the inverse link of the refitted model, and `target()` is then
##' applied to those draws. A probit candidate on simulations from a logit
##' model is thus scored on the same prevalence as the truth. Candidates
##' predict with their own covariates and offset, and those with a nugget
##' include it in their predictions when the generating model has one.
##'
##' The target receives the mean rather than the linear predictor, which is
##' what the `transform` argument of [predict.RiskMap()] receives, because
##' models with different links share no linear-predictor scale. Every link is
##' monotone, so a target defined on the linear predictor of a given link can
##' be written as a function of the mean, as [qlogis()] gives the logit scale.
##'
##' The target of an area is `area_summary()` applied to the targets at the
##' grid locations that fall within it, with equal weight, unlike the exact
##' fractional coverage used by [predict.RiskMap()] with `boundaries`. Areas
##' without grid locations are left out of the area-level assessment, with a
##' warning. Predictive draws are joint within each area, which is all that
##' area-level targets need.
##'
##' Each simulation has its own seed, derived from `seed` or, without one,
##' from the current random number generator state, and shared by all models
##' so that their fits and predictions use common random numbers. A fit or prediction that fails does not stop the
##' assessment: its results are missing, and the failure is recorded in
##' `fits` and reported in a warning, as are warnings from the fits.
##'
##' @examples
##' library(sf)
##' data(italy_sim)
##' italy_subset <- italy_sim[1:30, ]
##' italy_grid <- italy_subset[!duplicated(st_coordinates(italy_subset)), ]
##' model <- specify_glgpm(
##'   y ~ pop_dens + gp(), italy_subset, family = "gaussian",
##'   parameters = list(beta = c(1, 0.001), sigma2 = 1, phi = 20,
##'                     sigma2_me = 0.1)
##' )
##' simulations <- simulate_glgpm(
##'   model, nsim = 1, what = c("data", "surface"),
##'   prediction_grid = italy_grid, seed = 1
##' )
##' boundary <- st_sf(
##'   region = "study_area",
##'   geometry = st_convex_hull(st_union(italy_grid))
##' )
##' assessment <- assess_simulation(
##'   simulations,
##'   models = list(candidate = model),
##'   spatial_scale = c("grid", "area"),
##'   area_summary = mean,
##'   boundaries = boundary,
##'   pred_objective = "mse",
##'   messages = FALSE
##' )
##'
##' @export
assess_simulation <- function(obj_sim,
                       models,
                       control_mcmc = set_control_mcmc(),
                       seed = NULL,
                       spatial_scale,
                       messages = getOption("RiskMap.messages", TRUE),
                       target = identity,
                       area_summary = NULL,
                       boundaries = NULL,
                       col_names = NULL,
                       pred_objective = c("mse", "classify"),
                       categories = NULL) {

  if (!inherits(obj_sim, "RiskMap_simulation")) {
    abort_input("{.arg obj_sim} must be made by {.fn simulate_glgpm}, not {.obj_type_friendly {obj_sim}}.", arg = "obj_sim")
  }
  if (!is.character(pred_objective) || length(pred_objective) == 0 ||
      anyDuplicated(pred_objective) ||
      !all(pred_objective %in% c("mse", "classify"))) {
    abort_input("{.arg pred_objective} must be {.val mse}, {.val classify} or both.", arg = "pred_objective")
  }
  if (!is.character(spatial_scale) || length(spatial_scale) == 0 ||
      anyDuplicated(spatial_scale) ||
      !all(spatial_scale %in% c("grid", "area"))) {
    abort_input("{.arg spatial_scale} must be {.val grid}, {.val area} or both.", arg = "spatial_scale")
  }
  check_control(control_mcmc, "set_control_mcmc")
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)

  want_grid <- "grid" %in% spatial_scale
  want_area <- "area" %in% spatial_scale
  want_mse <- "mse" %in% pred_objective
  want_classify <- "classify" %in% pred_objective

  if (!is.list(models) || !length(models) ||
      is.null(names(models)) || any(!nzchar(names(models))) ||
      anyDuplicated(names(models)) ||
      !all(vapply(models, function(model) {
        inherits(model, "RiskMap") ||
          inherits(model, "RiskMap_simulation_model")
      }, logical(1)))) {
    abort_input("{.arg models} must be a non-empty list, with distinct names, of fits from {.fn glgpm} or models from {.fn specify_glgpm}.", arg = "models")
  }

  if (want_area) {
    if (is.null(boundaries)) {
      abort_input("{.arg boundaries} must be given when {.arg spatial_scale} includes {.val area}.", arg = "boundaries")
    }
    check_data(boundaries, "polygon")
  }

  obj_sim <- simulation_assessment_data(obj_sim)
  if (!is.function(target)) {
    abort_input("{.arg target} must be a function of the mean, such as {.fn identity} for prevalence.", arg = "target")
  }

  apply_target <- function(x, context) {
    result <- tryCatch(
      target(x),
      error = function(e) {
        abort_input("{.arg target} failed for {context}.", arg = "target", parent = e)
      }
    )
    if (!is.numeric(result) || !is.matrix(result) ||
        !identical(dim(result), dim(x)) || anyNA(result) ||
        any(!is.finite(result))) {
      abort_input("{.arg target} must return a finite numeric matrix of the same dimensions as its input, and didn't for {context}.", arg = "target")
    }
    result
  }

  if (want_classify) {
    if (is.null(categories)) {
      abort_input("{.arg categories} must be given when {.arg pred_objective} includes {.val classify}.", arg = "categories")
    }
    if (!is.numeric(categories) || length(categories) < 3 ||
        anyNA(categories) || any(!is.finite(categories)) ||
        any(diff(categories) <= 0)) {
      abort_input("{.arg categories} must be at least three finite, strictly increasing numbers.", arg = "categories")
    }
  }
  n_sim <- length(obj_sim$data_sim)
  n_models <- length(models)

  if(want_area && !is.function(area_summary)) {
    abort_input("{.arg area_summary} must be a function when {.arg spatial_scale} includes {.val area}.", arg = "area_summary")
  }

  apply_area_summary <- function(x, context) {
    result <- tryCatch(
      area_summary(x),
      error = function(e) {
        abort_input("{.arg area_summary} failed for {context}.", arg = "area_summary", parent = e)
      }
    )
    if (!is.numeric(result) || length(result) != 1L || is.na(result) ||
        !is.finite(result)) {
      abort_input("{.arg area_summary} must return one finite number, and didn't for {context}.", arg = "area_summary")
    }
    as.numeric(result)
  }
  model_names <- names(models)
  formulas <- lapply(models, function(template) update(template$formula, y ~ .))
  grid <- obj_sim$grid
  include_nugget <- obj_sim$nugget_over_grid

  # Every candidate needs complete covariates at every grid location, so that
  # its predictions line up with the true surface
  for (i in seq_len(n_models)) {
    template_i <- models[[i]]
    if (template_i$family != obj_sim$family) {
      abort_input("The model {.val {model_names[i]}} has the {template_i$family} family, but the simulations have the {obj_sim$family} family.", arg = "models")
    }
    re_names <- if (inherits(template_i, "RiskMap")) names(template_i$re) else
      template_i$re_terms
    required <- setdiff(get_formula_terms(formulas[[i]]), c("y", re_names))
    missing_columns <- setdiff(required, names(grid))
    if (length(missing_columns) > 0) {
      abort_input("The prediction grid of {.arg obj_sim} has no column{?s} {.val {missing_columns}}, which the model {.val {model_names[i]}} uses.", arg = "models")
    }
    n_incomplete <- sum(!complete.cases(st_drop_geometry(grid)[required]))
    if (n_incomplete > 0) {
      abort_input(c("The model {.val {model_names[i]}} needs complete values of {.val {required}} on the prediction grid, but {n_incomplete} grid location{?s} {?has/have} missing values.", "i" = "Remove them from the prediction grid before simulating."), arg = "models")
    }
  }

  if (want_area) {
    n_reg <- nrow(boundaries)

    if(is.null(col_names)) {
      boundaries$region <- paste0("reg", 1:n_reg)
      col_names <- "region"
    } else if (!col_names %in% names(boundaries)) {
      abort_input("{.arg col_names} must name a column of {.arg boundaries}, and {.val {col_names}} doesn't.", arg = "col_names")
    }
    names_reg <- boundaries[[col_names]]
    if (anyDuplicated(names_reg)) {
      abort_input("The column {.val {col_names}} of {.arg boundaries} must hold distinct names of areas.", arg = "col_names")
    }
    boundaries <- st_transform(boundaries, st_crs(grid))
    inter <- st_intersects(boundaries, grid)
    empty <- lengths(inter) == 0
    if (any(empty)) {
      empty_names <- names_reg[empty]
      riskmap_warn("No grid locations fall within {.val {empty_names}}, so {cli::qty(length(empty_names))}{?it is/they are} left out of the assessment of areas.")
    }
    areas <- which(!empty)
    # Joint predictive draws are needed only within each area; areas that
    # share grid locations form one block
    groups <- prediction_boundary_groups(
      list(prediction_row = unlist(inter[areas]),
           area_id = rep(areas, lengths(inter[areas]))),
      nrow(grid)
    )
  }

  if (want_classify) {
    category_labels <- paste0("(", utils::head(categories, -1), ",",
                              categories[-1], "]")
    categories_class <- factor(category_labels, levels = category_labels)
  }

  true_target_grid_sim <- apply_target(
    obj_sim$mean_sim, "the simulated true surface"
  )

  # Applies 'area_summary' to the targets of each area, one column per draw
  summarise_areas <- function(target, context) {
    t(vapply(areas, function(h) {
      values <- tryCatch(
        apply(target[inter[[h]], , drop = FALSE], 2, area_summary),
        error = function(e) {
          abort_input("{.arg area_summary} failed for the area {.val {names_reg[h]}} in {context}.", arg = "area_summary", parent = e)
        }
      )
      if (!is.numeric(values) || length(values) != ncol(target) ||
          any(!is.finite(values))) {
        abort_input("{.arg area_summary} must return one finite number, and didn't for the area {.val {names_reg[h]}} in {context}.", arg = "area_summary")
      }
      values
    }, numeric(ncol(target))))
  }

  if (want_area) {
    true_target_area_sim <- matrix(NA_real_, nrow = n_reg, ncol = n_sim)
    true_target_area_sim[areas, ] <- summarise_areas(
      true_target_grid_sim, "the true simulations"
    )
  }

  # One `mse`/`classify` store per requested spatial scale, all sharing the
  # same layout - only what each store is computed from differs below.
  init_objective_store <- function() {
    store <- list()
    if(want_mse) {
      store$mse <- array(NA, c(n_models, n_sim))
      rownames(store$mse) <- model_names
      colnames(store$mse) <- paste0("sim_", 1:n_sim)
    }
    if(want_classify) {
      store$classify <- setNames(vector("list", length(model_names)), model_names)
      for(i in 1:n_models) {
        store$classify[[model_names[i]]] <- list(by_cat = vector("list", n_sim),
                                                  across_cat = list())
        for(j in 1:n_sim) {
          store$classify[[model_names[i]]]$by_cat[[j]] <-
            data.frame(
              Class = categories_class,
              Sensitivity = NA,
              Specificity = NA,
              PPV = NA,
              NPV = NA,
              CC = NA
            )
        }
        store$classify[[model_names[i]]]$CC <- rep(NA, n_sim)
      }
    }
    store
  }

  out <- list(pred_objective = list())
  if(want_grid) out$pred_objective$grid <- init_objective_store()
  if(want_area) out$pred_objective$area <- init_objective_store()

  # Scores of one model on one simulation at one spatial scale. The
  # classification keeps every category in the confusion matrix, including
  # categories absent from this simulation.
  score_targets <- function(true_vals, samples) {
    scores <- list()
    if (want_mse) scores$mse <- mean((rowMeans(samples) - true_vals)^2)
    if (want_classify) {
      scores$classify <- simulation_classification_metrics(
        true_vals, samples, categories, levels(categories_class)
      )
    }
    scores
  }

  store_scores <- function(store, i, j, scores) {
    if (want_mse) store$mse[i, j] <- scores$mse
    if (want_classify) {
      store$classify[[model_names[i]]]$by_cat[[j]] <- scores$classify$by_cat
      store$classify[[model_names[i]]]$CC[j] <- scores$classify$overall_cc
    }
    store
  }

  # One seed per simulation, shared by the models so that their fits and
  # predictions use common random numbers. Each simulation's results then do
  # not depend on the others, which also allows running them separately.
  if (!is.null(seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(seed)
  }
  sim_seeds <- sample.int(.Machine$integer.max, n_sim)
  restore_after_tasks <- preserve_random_seed()
  on.exit(restore_after_tasks(), add = TRUE, after = FALSE)
  # The fits of a simulation take their seed from the first draw of its
  # stream, and its predictions use that stream from the start.
  fit_seeds <- vapply(sim_seeds, function(sim_seed) {
    set.seed(sim_seed)
    sample.int(.Machine$integer.max, 1L)
  }, integer(1))

  # Fits, predicts and scores one model on one simulation. A failure leaves
  # that pair's results missing instead of stopping the assessment, and
  # warnings are recorded with the fit.
  assess_pair <- function(i, j) {
    set.seed(sim_seeds[j])
    fit_warnings <- character()
    result <- tryCatch(withCallingHandlers({
      refit_args <- assessment_refit_args(
        models[[i]], formulas[[i]], obj_sim$data_sim[[j]], control_mcmc
      )
      refit_args$seed <- fit_seeds[j]
      refit_args$return_samples <- FALSE
      fit_ij <- do.call(glgpm, refit_args)

      include_ij <- if (include_nugget && fitted_nugget_variance(fit_ij) > 0) {
        "nugget"
      } else {
        character()
      }
      prediction_data <- suppressMessages(
        validate_prediction_data(fit_ij, grid, include_ij)
      )$model_data
      components <- prediction_sample_components(
        fit_ij, prediction_data, if (want_area) "joint" else "marginal",
        include_ij, control_mcmc, messages = FALSE,
        groups = if (want_area) groups else NULL
      )
      context <- paste0("model '", model_names[i], "', simulation ", j)
      # The response draws are the linear-predictor draws mapped through the
      # refitted model's own inverse link (the identity for a Gaussian fit),
      # so each model reaches the scale of the simulated mean through its own
      # link.
      target_ij <- apply_target(components$response, context)

      scores <- list(status = "ok", message = NA_character_,
                     estimate = unlist(coef(fit_ij)))
      if (want_grid) {
        scores$grid <- score_targets(true_target_grid_sim[, j], target_ij)
      }
      if (want_area) {
        scores$area <- score_targets(true_target_area_sim[areas, j],
                                     summarise_areas(target_ij, context))
      }
      scores
    }, warning = function(w) {
      fit_warnings <<- c(fit_warnings, condition_text(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      list(status = "failed", message = condition_text(e))
    })
    result$warnings <- unique(fit_warnings)
    result
  }

  assess_one_simulation <- function(j, report) {
    lapply(seq_len(n_models), function(i) {
      result <- assess_pair(i, j)
      report()
      result
    })
  }

  if (messages) {
    riskmap_inform("Refitting {n_models} model{?s} to {n_sim} simulation{?s}.")
  }
  bar <- progress_bar("Refitting models to simulations", n_sim * n_models,
                      messages)
  results <- riskmap_map(seq_len(n_sim), assess_one_simulation,
                         progress = bar, weights = rep(n_models, n_sim),
                         reports_progress = TRUE,
                         labels = paste("Simulation", seq_len(n_sim)))
  progress_done(bar)

  fits <- data.frame(
    model = rep(model_names, times = n_sim),
    simulation = rep(seq_len(n_sim), each = n_models),
    status = NA_character_, message = NA_character_, warnings = NA_character_
  )
  estimates <- setNames(vector("list", n_models), model_names)
  for (j in seq_len(n_sim)) {
    results_j <- results[[j]]
    for (i in seq_len(n_models)) {
      result_ij <- results_j[[i]]
      row <- (j - 1L) * n_models + i
      fits$status[row] <- result_ij$status
      fits$message[row] <- result_ij$message %||% NA_character_
      if (length(result_ij$warnings) > 0) {
        fits$warnings[row] <- paste(result_ij$warnings, collapse = " | ")
      }
      if (!identical(result_ij$status, "ok")) next

      if (is.null(estimates[[i]])) {
        estimates[[i]] <- matrix(
          NA_real_, n_sim, length(result_ij$estimate),
          dimnames = list(paste0("sim_", seq_len(n_sim)),
                          names(result_ij$estimate))
        )
      }
      estimates[[i]][j, ] <- result_ij$estimate
      if (want_grid) {
        out$pred_objective$grid <- store_scores(out$pred_objective$grid, i, j,
                                                result_ij$grid)
      }
      if (want_area) {
        out$pred_objective$area <- store_scores(out$pred_objective$area, i, j,
                                                result_ij$area)
      }
    }
  }

  n_failed <- sum(fits$status != "ok")
  if (n_failed > 0) {
    first <- which(fits$status != "ok")[1]
    riskmap_warn(c("{n_failed} of {nrow(fits)} model fits failed, and their results are missing.", "i" = "The first failure, for model {.val {fits$model[first]}} and simulation {fits$simulation[first]}, was {.val {fits$message[first]}}.", "i" = "See the {.code fits} element of the result."))
  }
  n_warned <- sum(!is.na(fits$warnings))
  if (n_warned > 0) {
    first <- which(!is.na(fits$warnings))[1]
    riskmap_warn(c("{n_warned} of {nrow(fits)} model fits gave warnings.", "i" = "The first, for model {.val {fits$model[first]}} and simulation {fits$simulation[first]}, was {.val {fits$warnings[first]}}.", "i" = "See the {.code fits} element of the result."))
  }

  if(want_classify) {
    if(want_grid) out$pred_objective$grid$classify$Class <- categories_class
    if(want_area) out$pred_objective$area$classify$Class <- categories_class
  }
  out$fits <- fits
  out$estimates <- estimates
  out$n_sim <- n_sim
  out$spatial_scale <- spatial_scale
  class(out) <- "RiskMap_assess_simulation"
  return(out)
}

##' Build a refit call from a fitted assessment-model template
##'
##' The fitted coefficients are intentionally excluded. Data-dependent inputs
##' are resolved against each simulated dataset rather than copied as vectors
##' from the original fit.
##' @noRd
assessment_refit_args <- function(template, formula, data, control_mcmc) {
  fitted_template <- inherits(template, "RiskMap")
  args <- list(
    formula = formula,
    family = template$family,
    data = data,
    distance_units = template$distance_units,
    control_mcmc = control_mcmc,
    control_mcml = if (fitted_template) {
      attr(template, "control_mcml") %||% set_control_mcml()
    } else {
      set_control_mcml()
    },
    fix_var_me = if (fitted_template) template$fix_var_me else NULL,
    messages = FALSE
  )
  if (!isTRUE(fit_traits(template)$exact)) {
    args$denominator <- quote(units_m)
    if (fitted_template) {
      args$invlink <- refit_invlink(template)
    } else if (isTRUE(template$custom_link)) {
      args$invlink <- template$invlink
    }
  }
  args
}

##' Classification metrics for one simulated dataset
##'
##' @param true_values True target values.
##' @param samples Matrix of predictive target samples, one row per target.
##' @param breaks Category boundaries.
##' @param labels Fixed category labels.
##' @return A list containing category-specific metrics and overall accuracy.
##' @noRd
simulation_classification_metrics <- function(true_values, samples, breaks, labels) {
  true_class <- cut(true_values, breaks = breaks, labels = labels)
  true_class <- factor(true_class, levels = labels)

  probabilities <- do.call(cbind, lapply(seq_len(length(breaks) - 1L), function(i) {
    rowMeans(breaks[i] < samples & samples <= breaks[i + 1L])
  }))
  predicted_class <- factor(labels[max.col(probabilities, ties.method = "first")],
                            levels = labels)
  confusion <- table(true_class, predicted_class)
  total <- sum(confusion)

  divide_or_na <- function(numerator, denominator) {
    if (denominator == 0) NA_real_ else numerator / denominator
  }
  by_category <- lapply(seq_along(labels), function(i) {
    true_positive <- confusion[i, i]
    false_positive <- sum(confusion[, i]) - true_positive
    false_negative <- sum(confusion[i, ]) - true_positive
    true_negative <- total - true_positive - false_positive - false_negative
    data.frame(
      Class = labels[i],
      Sensitivity = divide_or_na(true_positive, true_positive + false_negative),
      Specificity = divide_or_na(true_negative, true_negative + false_positive),
      PPV = divide_or_na(true_positive, true_positive + false_positive),
      NPV = divide_or_na(true_negative, true_negative + false_negative),
      CC = divide_or_na(true_positive + true_negative, total)
    )
  })

  valid <- !is.na(true_class)
  overall_cc <- if (any(valid)) {
    mean(true_class[valid] == predicted_class[valid])
  } else {
    NA_real_
  }
  list(by_cat = do.call(rbind, by_category), overall_cc = overall_cc)
}

##' @title Print an Assessment of Simulations
##'
##' @description Prints a short description of an assessment from
##' [assess_simulation()]. Use [summary.RiskMap_assess_simulation()] for its
##' results.
##'
##' @param x An object of class `RiskMap_assess_simulation`.
##' @param ... Not used.
##'
##' @return Invisibly returns `x`.
##'
##' @method print RiskMap_assess_simulation
##' @export
print.RiskMap_assess_simulation <- function(x, ...) {
  model_names <- unique(x$fits$model %||% character())
  cat(sprintf("Assessment of %d model%s over %d simulations\n",
              length(model_names), if (length(model_names) == 1) "" else "s",
              x$n_sim))
  if (length(model_names)) {
    cat(sprintf("Models: %s\n", paste(model_names, collapse = ", ")))
  }
  objectives <- unique(unlist(lapply(x$pred_objective, names)))
  cat(sprintf("Spatial scale: %s; objectives: %s\n",
              paste(x$spatial_scale, collapse = ", "),
              paste(objectives, collapse = ", ")))
  if (!is.null(x$fits)) {
    cat(sprintf("Failed model fits: %d of %d\n",
                sum(x$fits$status != "ok"), nrow(x$fits)))
  }
  cat("Use summary() to see the results.\n")
  invisible(x)
}

##' @title Summarize Simulation Results
##'
##' @description Summarizes the results of model evaluations from a `RiskMap_assess_simulation` object. Provides average metrics for classification by category and overall correct classification (CC) summary.
##'
##' @param object An object of class `RiskMap_assess_simulation`, as returned by `assess_simulation`.
##' @param ... Additional arguments (not used).
##'
##' @return A list of class `summary.RiskMap_assess_simulation` with one
##' element per assessed spatial scale (`grid` and/or `area`), each holding:
##' - `mse`: a data frame with, for each model, the number of simulations
##'   (`n_sim`), the number with a result (`n_valid`), and the mean and
##'   standard deviation of the MSE across simulations (`MSE_mean`, `MSE_sd`);
##' - `classify`: for each model, `classify_res`, the average sensitivity,
##'   specificity, PPV, NPV and CC by category across simulations, `n_valid`,
##'   the number of simulations with a defined value of each metric, `n_sim`,
##'   and `cc_summary`, the mean, standard deviation and 2.5% and 97.5%
##'   quantiles of the overall proportion of correct classification.
##'
##' It also holds `fits`, the number of model fits and of failed fits.
##'
##' @method summary RiskMap_assess_simulation
##' @export
summary.RiskMap_assess_simulation <- function(object, ...) {
  check_class(object, "RiskMap_assess_simulation")

  if (identical(object$n_sim, 1L)) {
    riskmap_warn("The assessment has one simulation, so point estimates are shown without their uncertainty across simulations.")
  }

  # `object$pred_objective` holds one element per requested spatial scale
  # ("grid" and/or "area"); each is summarized the same way.
  results <- list()
  for (scale in intersect(c("grid", "area"), names(object$pred_objective))) {
    results[[scale]] <- summarize_pred_objective(object$pred_objective[[scale]])
  }
  if (!is.null(object$fits)) {
    results$fits <- c(n = nrow(object$fits),
                      failed = sum(object$fits$status != "ok"))
  }

  # Assign class for S3 print method
  class(results) <- "summary.RiskMap_assess_simulation"
  return(results)
}

##' Summarize one spatial scale's `mse`/`classify` results
##'
##' @param pred_objective The `mse`/`classify` element of a
##'   `RiskMap_assess_simulation` object for a single spatial scale.
##' @return A list with `mse` and/or `classify` summaries.
##' @noRd
summarize_pred_objective <- function(pred_objective) {
  results <- list()

  # Check for "mse" in pred_objective
  if ("mse" %in% names(pred_objective)) {
    mse_data <- pred_objective$mse

    # Check if mse_data is a matrix
    if (is.matrix(mse_data)) {
      n_valid <- rowSums(is.finite(mse_data))
      mse_summary <- data.frame(
        Model = rownames(mse_data),
        n_sim = ncol(mse_data),
        n_valid = n_valid,
        MSE_mean = apply(mse_data, 1, finite_mean),
        MSE_sd = apply(mse_data, 1, finite_sd),
        row.names = NULL
      )

      results$mse <- mse_summary
    } else {
      abort_input("The mean squared errors in {.arg object} must be a matrix.", arg = "object")
    }
  }

  # Check for "classify" in pred_objective
  if ("classify" %in% names(pred_objective)) {
    classify_data <- pred_objective$classify

    name_models <- setdiff(names(classify_data), "Class")
    results$classify <- list()
    for(model_name in name_models) {
      model_data <- classify_data[[model_name]]
      n_sim <- length(model_data$by_cat)
      metric_names <- setdiff(names(model_data$by_cat[[1]]), "Class")
      classes <- model_data$by_cat[[1]]$Class
      metric_array <- array(
        NA_real_,
        dim = c(length(classes), length(metric_names), n_sim),
        dimnames = list(classes, metric_names, paste0("sim_", seq_len(n_sim)))
      )
      for(j in seq_len(n_sim)) {
        current <- model_data$by_cat[[j]]
        row_index <- match(classes, current$Class)
        metric_array[, , j] <- as.matrix(current[row_index, metric_names, drop = FALSE])
      }
      metric_mean <- apply(metric_array, c(1, 2), finite_mean)
      metric_n_valid <- apply(is.finite(metric_array), c(1, 2), sum)
      classify_res <- data.frame(Class = classes, metric_mean,
                                 check.names = FALSE, row.names = NULL)
      n_valid <- data.frame(Class = classes, metric_n_valid,
                            check.names = FALSE, row.names = NULL)

      results$classify[[model_name]] <- list(
        classify_res = classify_res,
        n_valid = n_valid,
        n_sim = n_sim,
        cc_summary = summarize_simulation_metric(model_data$CC, n_sim)
      )
    }
  }

  results
}

##' @noRd
finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

##' @noRd
finite_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) >= 2L) sd(x) else NA_real_
}

##' @noRd
summarize_simulation_metric <- function(x, n_sim = length(x)) {
  x <- x[is.finite(x)]
  n_valid <- length(x)
  list(
    mean = if (n_valid) mean(x) else NA_real_,
    sd = if (n_valid >= 2L) sd(x) else NA_real_,
    lower = if (n_valid >= 2L) unname(quantile(x, 0.025)) else NA_real_,
    upper = if (n_valid >= 2L) unname(quantile(x, 0.975)) else NA_real_,
    n_valid = n_valid,
    n_sim = n_sim
  )
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
##' @method print summary.RiskMap_assess_simulation
##' @export
print.summary.RiskMap_assess_simulation <- function(x, ...) {
  cat("Summary of Simulation Results\n\n")
  if (!is.null(x$fits) && x$fits[["failed"]] > 0) {
    cat(sprintf("Failed model fits (results missing): %d of %d\n\n",
                x$fits[["failed"]], x$fits[["n"]]))
  }

  scale_label <- c(grid = "Grid", area = "Area")
  scales <- intersect(c("grid", "area"), names(x))
  for (scale in scales) {
    if (length(scales) > 1) cat(sprintf("== %s-level results ==\n\n", scale_label[[scale]]))
    print_pred_objective(x[[scale]])
  }

  invisible(x)
}

##' Print a data frame as plain text
##'
##' Mirrors `print.data.frame()` but prints the formatted character matrix
##' directly, so RStudio notebooks keep the table inline with the surrounding
##' `cat()` output instead of rendering it as a separate data frame widget.
##'
##' @param df A data frame.
##' @return Invisibly returns `df`.
##' @noRd
print_text_table <- function(df) {
  print(as.matrix(format(df)), quote = FALSE, right = TRUE)
  invisible(df)
}

##' Print one spatial scale's `mse`/`classify` summary
##'
##' @param x An element of a `summary.RiskMap_assess_simulation` object, as
##'   returned by `summarize_pred_objective()`.
##' @return Invisibly returns `x`.
##' @noRd
print_pred_objective <- function(x) {
  if (!is.null(x$mse)) {
    cat("Mean Squared Error (MSE):\n")
    print_text_table(x$mse)
    cat("\n")
  }

  if (!is.null(x$classify)) {
    cat("Classification Results:\n")

    # Iterate over each model in classify results
    for (model_name in names(x$classify)) {
      model_data <- x$classify[[model_name]]

      cat(sprintf("\nModel: %s\n", model_name))

      cat("\nAverages across simulations by Category:\n")
      print_text_table(model_data$classify_res)

      cat("\nNumber of valid simulations by Category and metric ",
          sprintf("(out of %d):\n", model_data$n_sim), sep = "")
      print_text_table(model_data$n_valid)

      cat("\nProportion of Correct Classification (CC) across categories:\n")
      cc_summary <- model_data$cc_summary
      if (cc_summary$n_valid >= 2L) {
        cat(sprintf("Mean: %.3f, SD: %.3f, 95%% simulation interval: [%.3f, %.3f] ",
                    cc_summary$mean, cc_summary$sd,
                    cc_summary$lower, cc_summary$upper))
      } else {
        cat(sprintf("Mean: %.3f, uncertainty unavailable ", cc_summary$mean))
      }
      cat(sprintf("(n_valid = %d of %d)\n",
                  cc_summary$n_valid, cc_summary$n_sim))
    }
    cat("\n")
  }

  invisible(x)
}

#' Name of the response variable, used to label prediction maps
#' @noRd
prediction_response_name <- function(object) {
  response <- all.vars(object$formula[[2L]])
  if (length(response)) response[[1L]] else NA_character_
}
