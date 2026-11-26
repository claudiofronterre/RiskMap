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
##'   locations, with one row for each record of the data.
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
##' @param joint `NULL`, `TRUE` or `FALSE`. With `TRUE` the predictive draws
##'   at the locations are joint, so they keep the dependence between
##'   locations, and with `FALSE` each location is drawn on its own (marginal
##'   prediction), which is faster for many locations. `NULL`, the default,
##'   draws jointly when `boundaries` are given and marginally otherwise.
##'   Summaries over `boundaries` need joint draws, so `FALSE` can't be used
##'   with them. At the observed locations of a binomial, Poisson or negative
##'   binomial fit the draws are the joint conditional samples whatever
##'   `joint` is.
##' @param include `NULL`, or a character vector specifying optional components
##'   to add to the complete linear predictor and response. Supported values are
##'   `"random_effects"`, for conditional predictions at known groups, and
##'   `"nugget"`, for a new realization of observation-level predictive
##'   variation at locations supplied through `newdata`. Fixed effects, the
##'   spatial process and any offset are always returned separately. A level of
##'   a random effect that isn't in the data of the fit stops a prediction that
##'   includes `"random_effects"`, and doesn't matter to one that doesn't. The
##'   nugget is one value for each location, so with joint draws rows of
##'   `newdata` at the same location share it, as they share the spatial
##'   process.
##' @param transform Optional named list of functions defining derived predictive
##'   targets from the complete linear-predictor draws. Each function receives a
##'   location-by-draw matrix and must return a numeric matrix with the same
##'   dimensions. The ordinary response-scale prediction is always calculated,
##'   so this argument is needed only for additional scientific targets, for
##'   example `list(stunting = function(link) pnorm((-2 - link) / sd_child))`.
##' @param standardise Optional data frame describing a common reference
##'   population over which response-scale predictions are averaged. Its columns
##'   must be fixed-effect predictors from the fitted model and override the
##'   corresponding columns in `newdata`. Rows receive equal weight unless an
##'   optional numeric `.weight` column is supplied. Predictor values kept
##'   together on each row are standardised jointly, preserving their empirical
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
##' @param ... Not used, and an error if given, so that a misspelt or former
##'   argument such as `type` is not ignored.
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
##' latent and the parameter uncertainty. For binomial, Poisson and negative binomial models the
##' conditional samples at each parameter value come from a short run of the
##' same Langevin (MALA) sampler used for plug-in prediction, in coordinates
##' whitened by the Laplace approximation at that value, started from a draw
##' of that approximation and using the step size tuned at the estimate (100
##' burn-in iterations, then thinning by 10). The total number of predictive
##' draws equals that of plug-in prediction, 2,000 by default. Parameters
##' without a finite positive variance, such as those held at a boundary, keep
##' their estimated value. So does a variance, a scale or a size whose
##' logarithm has a standard error above 3, as a variance of a Gaussian model
##' estimated at zero has, because draws on the log scale would span more than
##' five orders of magnitude. The other parameters are then drawn given it, and
##' a warning names it unless it is a variance close to zero. When the
##' covariance matrix of the estimates is not
##' available, an error is returned; refitting with
##' `set_control_mcml(pc_prior = TRUE)` usually resolves this. The location
##' column `fixed` then reports the average fixed-effect contribution.
##'
##' At the observed locations of a binomial, Poisson or negative binomial fit
##' with a nugget, the effect of each location is drawn given the data and
##' already includes the nugget. At those of a Gaussian fit the prediction is of
##' the spatial process \eqn{S(x)}, without the nugget, as it is at new
##' locations, and [fitted()][fitted.RiskMap()] gives the fitted values with the
##' nugget and the random effects of the data. For every family,
##' `include = "nugget"` is available only for predictions supplied through
##' `newdata`, where it represents a new independent nugget realization.
##'
##' `newdata` and `boundaries` are transformed internally to the fitted model
##' CRS without changing row order; returned geometries remain in their input
##' CRS. Rows with missing required predictors, or with predictors that give the
##' model a covariate that isn't finite, such as `log(0)`, are retained with `NA`
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
##' When `standardise` is supplied, every conditional spatial draw is reused
##' across all rows of the reference population. The response and any derived
##' targets are evaluated for each reference profile and then averaged using
##' `.weight`. This is direct standardisation (empirical or counterfactual
##' averaging), not simulation of uncertain predictors. The reported `fixed`
##' and `link_*` quantities use the weighted mean fixed-effect contribution;
##' with a nonlinear link, transforming that mean is not the same as the
##' standardised `response_*` quantity.
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
##'     samples for location-level prediction, with one row for each row of
##'     `predictions`, and, with parameter
##'     uncertainty, the fixed-effect contribution of each draw (`fixed`).
##'     For polygon aggregation it is nested first by aggregation and then by
##'     component, with one matrix row per polygon.}
##'   \item{include}{Optional components included in the link and response.}
##'   \item{transform}{Names of additional derived predictive targets.}
##'   \item{standardised}{Whether response predictions were averaged over a
##'     reference population supplied through `standardise`.}
##'   \item{aggregate}{Names of the within-polygon aggregations.}
##'   \item{summaries}{Names of the calculated summaries.}
##'   \item{thresholds}{Thresholds used for exceedance probabilities.}
##'   \item{level}{Probability covered by the interval limits.}
##'   \item{weights}{Name of the weights column used for the polygon
##'     summaries, or `NULL`.}
##'   \item{joint}{`TRUE` when the draws were joint and `FALSE` when they were
##'     marginal.}
##'   \item{parameter_draws}{The number of parameter draws when parameter
##'     uncertainty was propagated, otherwise `NULL`.}
##'   \item{family}{The fitted response family.}
##'   \item{link}{The link of the fit when it is not the canonical one, its
##'     name or `"custom"`, and `NULL` for the canonical link.}
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
##' standardised_prediction <- predict(
##'   fit,
##'   newdata = example_data[1:5, ],
##'   standardise = reference_population,
##'   messages = FALSE
##' )
##' standardised_prediction$predictions
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
##' @seealso The methods for a prediction are
##'   [print()][print.RiskMap_prediction()],
##'   [to_table()][to_table.RiskMap_prediction()] and [plot_prediction()],
##'   which `plot()` of the prediction also calls. [assess_extrapolation()]
##'   shows where predictions go beyond the data, and
##'   [fitted()][fitted.RiskMap()] gives the fitted values at the data.
##' @method predict RiskMap
##' @export
predict.RiskMap <- function(
    object,
    newdata = NULL,
    boundaries = NULL,
    aggregate = NULL,
    joint = NULL,
    include = NULL,
    transform = NULL,
    standardise = NULL,
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
  check_dots_unused(..., .hints = c(
    type = "Use {.arg joint} to choose between joint and marginal draws. The response, link and spatial effect are always predicted.",
    standardize = "It is spelt {.arg standardise}."
  ))

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
  standardise <- validate_prediction_standardise(object, standardise)

  type <- validate_prediction_options(
    object = object,
    joint = joint,
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
    reason <- if (fit_is_exact(object)) {
      "The prediction there is of the spatial process, without the nugget of the data, which {.fn fitted} includes. For a new realisation of it, give the locations as {.arg newdata}."
    } else {
      "There the effect of each location is drawn given the data, and it already includes the nugget."
    }
    abort_input(c("{.arg include} can't contain {.val nugget} at the observed locations.", "i" = reason), arg = "include")
  }
  prediction_input <- validate_prediction_data(
    object, newdata, include, standardise, messages = messages,
    check_finite = TRUE
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
      standardise = standardise,
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
      standardised = !is.null(standardise),
      aggregate = NULL,
      summaries = names(summaries),
      thresholds = thresholds,
      level = level,
      joint = type == "joint",
      family = object$family,
      link = fit_link(object),
      response_name = prediction_response_name(object),
      observed_locations = observed_locations,
      parameter_draws = parameter_draws,
      coverage = NULL,
      weights = NULL
    )
    class(out) <- "RiskMap_prediction"
    return(stamp_version(out))
  }
  # At the observed locations of a non-Gaussian fit the draws are the joint
  # conditional samples.
  joint_draws <- type == "joint" ||
    (observed_locations && !fit_is_exact(object))
  sample_components <- prediction_sample_components(
    object, prediction_data, type, include, control_mcmc, messages,
    observed_locations = observed_locations, transform = transform,
    standardise = standardise, parameter_draws = parameter_draws,
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
    standardised = !is.null(standardise),
    aggregate = NULL,
    summaries = names(summaries),
    thresholds = thresholds,
    level = level,
    joint = joint_draws,
    family = object$family,
    link = fit_link(object),
    response_name = prediction_response_name(object),
    observed_locations = observed_locations,
    parameter_draws = parameter_draws,
    coverage = NULL,
    weights = NULL
  )
  class(out) <- "RiskMap_prediction"
  stamp_version(out)
}

#' Validate the unified prediction options
#' @param call The call named in input errors, so that they name
#'   `predict.RiskMap()`.
#' @noRd
validate_prediction_options <- function(object, joint, include,
                                        boundaries, aggregate, summaries,
                                        thresholds, level, weights,
                                        control_mcmc, seed, keep_samples,
                                        messages, call = rlang::caller_env()) {
  if (!is.null(joint) && !isTRUE(joint) && !isFALSE(joint)) {
    supplied <- describe_value(joint)
    abort_input("{.arg joint} must be {.code NULL}, {.code TRUE} or {.code FALSE}, not {supplied}.", arg = "joint", call = call)
  }
  if (!is.null(boundaries)) {
    check_data(boundaries, geometry = "polygon", type = "sf")
    if (isFALSE(joint)) {
      abort_input("{.arg joint} can't be {.code FALSE} when {.arg boundaries} are given, because summaries over areas need joint draws.", arg = "joint", call = call)
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
    type <- if (isTRUE(joint)) "joint" else "marginal"
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

#' Validate a reference population for direct standardisation
#' @noRd
validate_prediction_standardise <- function(object, standardise) {
  if (is.null(standardise)) return(NULL)
  if (!is.data.frame(standardise) || inherits(standardise, "sf") ||
      !nrow(standardise)) {
    abort_input("{.arg standardise} must be a non-empty data frame without geometry, not {.obj_type_friendly {standardise}}.", arg = "standardise")
  }
  if (anyDuplicated(names(standardise)) || anyNA(names(standardise)) ||
      any(!nzchar(names(standardise)))) {
    abort_input("{.arg standardise} must have distinct, non-empty column names.", arg = "standardise")
  }
  reference_names <- setdiff(names(standardise), ".weight")
  response <- as.character(object$formula[[2L]])
  inter_f <- interpret_formula(fit_formula(object))
  offset_names <- inter_f$offset %||% character()
  fixed_names <- setdiff(
    get_formula_terms(object$formula),
    c(response, names(object$re), offset_names)
  )
  if (!length(reference_names) || any(!reference_names %in% fixed_names)) {
    invalid <- setdiff(reference_names, fixed_names)
    problem <- "{.arg standardise} must have one or more columns, each a covariate of the model."
    if (length(invalid)) {
      problem <- c(problem, "x" = "{.val {invalid}} {cli::qty(length(invalid))}{?is not a covariate/are not covariates}.")
    }
    abort_input(problem, arg = "standardise")
  }
  if (!all(complete.cases(standardise[reference_names]))) {
    abort_input("{.arg standardise} has missing values.", arg = "standardise")
  }
  weights <- if (".weight" %in% names(standardise)) {
    standardise$.weight
  } else {
    rep(1, nrow(standardise))
  }
  if (!is.numeric(weights) || anyNA(weights) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    abort_input("The column {.code .weight} of {.arg standardise} must hold finite, non-negative values with a positive total.", arg = "standardise")
  }
  # Each profile with the other covariates of a record of the data
  profiles <- st_drop_geometry(object$data)[rep(1L, nrow(standardise)), ,
                                            drop = FALSE]
  profiles <- as.data.frame(profiles)
  profiles[reference_names] <- standardise[reference_names]
  infinite <- !prediction_design_is_finite(object, profiles)
  if (any(infinite)) {
    abort_input("{sum(infinite)} row{?s} of {.arg standardise} {?gives/give} the model a covariate that isn't a finite number, such as the logarithm of zero or of a negative number.", arg = "standardise")
  }
  list(
    data = standardise,
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
#'
#' Levels of the random effects that aren't in the fit stop, unless
#' `allow_new_levels` is `TRUE`, for the assessment of predictions, whose test
#' sets can have levels of their own. Their effects are then drawn from the
#' distribution of the random effect (`draw_prediction_samples()`).
#' With `check_finite`, rows of `newdata` whose covariates give a term that
#' isn't finite are left out as those with missing values are.
#' @noRd
validate_prediction_data <- function(object, newdata, include,
                                     standardise = NULL,
                                     allow_new_levels = FALSE,
                                     messages = TRUE, check_finite = FALSE) {
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
    if (messages) {
      riskmap_inform("{.arg newdata} is projected to the coordinate reference system of the fit.")
    }
    prediction_data <- st_transform(prediction_data, st_crs(object$data))
  }
  inter_f <- interpret_formula(fit_formula(object))
  response <- as.character(object$formula[[2]])
  required <- setdiff(get_formula_terms(object$formula), response)
  if (!"random_effects" %in% include) {
    required <- setdiff(required, names(object$re))
  }
  standardised_names <- if (is.null(standardise)) character() else
    setdiff(names(standardise$data), ".weight")
  required <- setdiff(required, standardised_names)
  missing_columns <- setdiff(required, names(prediction_data))
  if (length(missing_columns) > 0L) {
    abort_input("{.arg newdata} must have the column{?s} {.val {missing_columns}}, which the model uses.", arg = "newdata")
  }
  for (name in standardised_names) {
    prediction_data[[name]] <- standardise$data[[name]][1L]
  }
  if ("random_effects" %in% include && !allow_new_levels) {
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
  if (check_finite && !is.null(newdata) && any(valid_rows)) {
    finite <- prediction_design_is_finite(
      object, st_drop_geometry(prediction_data)[valid_rows, , drop = FALSE]
    )
    if (any(!finite)) {
      riskmap_warn("{sum(!finite)} prediction row{?s} {?gives/give} the model a covariate that isn't a finite number, such as the logarithm of zero or of a negative number, so {?its/their} predictions are {.code NA}.")
      valid_rows[valid_rows] <- finite
    }
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
#'
#' `boundaries` are in the CRS of the fit and `supplied` as the user gave them.
#' @noRd
assign_prediction_locations <- function(model_data, boundaries,
                                        supplied = boundaries) {
  metadata <- attr(model_data, "grid_metadata", exact = TRUE)
  if (inherits(model_data, "RiskMap_grid") && !is.null(metadata)) {
    return(assign_prediction_grid_cells(model_data, boundaries, metadata,
                                        supplied))
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

#' Whether any polygons overlap by more than rounding
#'
#' An overlap of up to a millionth of the smaller polygon is rounding.
#' @noRd
polygons_overlap <- function(polygons) {
  pairs <- st_overlaps(polygons)
  shapes <- st_geometry(polygons)
  any(vapply(seq_along(pairs), function(first) {
    any(vapply(pairs[[first]][pairs[[first]] > first], function(second) {
      tryCatch({
        shared <- sum(as.numeric(
          st_area(st_intersection(shapes[first], shapes[second]))
        ))
        shared > 1e-6 * min(as.numeric(st_area(shapes[c(first, second)])))
      }, error = function(error) TRUE)
    }, logical(1)))
  }, logical(1)))
}

#' Project polygons along their edges
#'
#' Projecting polygons moves their vertices and keeps the edges between them
#' straight, so a polygon with a vertex on the edge of its neighbour, as at a
#' T-junction, overlaps it once projected, or leaves a gap. Edges with vertices
#' added along them, in the plane of the supplied coordinates, stay together.
#' @noRd
project_polygon_edges <- function(polygons, crs) {
  shapes <- st_set_crs(st_geometry(polygons), NA)
  box <- st_bbox(shapes)
  longest <- max(box[["xmax"]] - box[["xmin"]], box[["ymax"]] - box[["ymin"]])
  shapes <- st_segmentize(shapes, dfMaxLength = longest / 1000)
  st_geometry(polygons) <- st_set_crs(shapes, st_crs(polygons))
  st_transform(polygons, crs)
}

#' Calculate exact polygon coverage for cells created by create_grid
#' @noRd
assign_prediction_grid_cells <- function(model_data, boundaries, metadata,
                                         supplied = boundaries) {
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
  # Polygons that only their projection makes overlap are projected again
  # along their edges.
  overlapping <- any(lengths(st_overlaps(grid_boundaries)) > 0L)
  if (overlapping) {
    grid_boundaries <- project_polygon_edges(supplied, grid_crs)
    overlapping <- polygons_overlap(grid_boundaries)
  }
  duplicates <- st_equals(grid_boundaries)
  if (overlapping || any(lengths(duplicates) > 1L)) {
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
                                    aggregate, include, transform, standardise,
                                    summaries,
                                    thresholds,
                                    level, weights,
                                    control_mcmc, seed, keep_samples,
                                    messages, parameter_draws = NULL) {
  model_data <- prediction_input$model_data
  output_boundaries <- boundaries
  if (!isTRUE(st_crs(boundaries) == st_crs(object$data))) {
    if (messages) {
      riskmap_inform("{.arg boundaries} are projected to the coordinate reference system of the fit.")
    }
    boundaries <- st_transform(boundaries, st_crs(object$data))
  }
  # Points on a grid that has lost its cells, in the CRS of newdata, are said
  # to be counted whole (#214)
  inform_point_assignment(prediction_input$output_data)

  assignment <- assign_prediction_locations(model_data, boundaries,
                                            output_boundaries)
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
      standardised = !is.null(standardise), aggregate = names(aggregate),
      summaries = names(summaries), thresholds = thresholds, level = level,
      joint = TRUE, family = object$family, link = fit_link(object),
      response_name = prediction_response_name(object),
      observed_locations = FALSE,
      parameter_draws = parameter_draws,
      coverage = assignment$method,
      weights = weights
    )
    class(out) <- "RiskMap_prediction"
    return(stamp_version(out))
  }

  prediction_groups <- prediction_boundary_groups(
    assignment, nrow(model_data)
  )
  sample_components <- prediction_sample_components(
    object, model_data, "joint", include, control_mcmc, messages,
    groups = prediction_groups, transform = transform,
    standardise = standardise, parameter_draws = parameter_draws,
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
    standardised = !is.null(standardise),
    aggregate = names(aggregate),
    summaries = names(summaries),
    thresholds = thresholds,
    level = level,
    joint = TRUE,
    family = object$family,
    link = fit_link(object),
    response_name = prediction_response_name(object),
    observed_locations = FALSE,
    parameter_draws = parameter_draws,
    coverage = assignment$method,
    weights = weights
  )
  class(out) <- "RiskMap_prediction"
  stamp_version(out)
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
                                        transform = list(), standardise = NULL,
                                        prediction_data = NULL) {
  spatial <- prediction_setup$S_samples
  if (!is.matrix(spatial)) {
    spatial <- matrix(spatial, nrow = nrow(prediction_setup$grid_pred))
  }
  # Samples of the distinct locations give each row those of its location.
  location <- prediction_setup$location_index
  if (!is.null(location)) spatial <- spatial[location, , drop = FALSE]
  n_locations <- nrow(spatial)
  n_samples <- ncol(spatial)
  fixed <- prediction_row_values(prediction_setup$mu_pred, n_locations)
  offset <- prediction_row_values(prediction_setup$cov_offset, n_locations)
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
  response <- if (fit_is_exact(object)) {
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
  if (!is.null(standardise)) {
    standardised <- standardise_prediction_targets(
      object = object,
      prediction_data = prediction_data,
      standardise = standardise,
      spatial = spatial,
      offset = offset,
      random_effects = random_effects,
      nugget = nugget,
      transform = transform,
      beta = prediction_setup$par_hat$beta
    )
    fixed <- standardised$fixed
    link <- standardised$link
    response <- standardised$response
    transformed <- standardised$transformed
  }
  list(fixed = fixed, offset = offset, spatial = spatial,
       random_effects = random_effects, nugget = nugget,
       link = link, response = response, transformed = transformed)
}

#' A value for each prediction row, from one for each or a single one
#' @noRd
prediction_row_values <- function(values, n_rows) {
  if (!length(values) %in% c(1L, n_rows)) {
    riskmap_bug("The prediction has {length(values)} values for {n_rows} rows.")
  }
  rep(values, length.out = n_rows)
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
                                         standardise = NULL,
                                         parameter_draws = NULL,
                                         mcmc_seed = NULL) {
  components_at <- function(fit, control) {
    setup <- prepare_prediction(
      fit, prediction_data, prediction, include, control, messages = FALSE,
      observed_locations = observed_locations, groups = groups
    )
    build_prediction_components(
      fit, setup, include, transform, standardise = standardise,
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
      object, setup, include, transform, standardise = standardise,
      prediction_data = prediction_data
    ))
  }

  n_total <- prediction_sample_size(object, control_mcmc)
  per_draw <- max(1L, ceiling(n_total / parameter_draws))
  parameters <- draw_parameter_values(object, parameter_draws)
  if (messages) {
    riskmap_inform("Propagating parameter uncertainty over {parameter_draws} parameter draws.")
  }
  sampler <- if (fit_is_exact(object)) NULL else
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
  bar <- progress_bar("Propagating parameter uncertainty", parameter_draws,
                      messages)
  for (k in seq_len(parameter_draws)) {
    fit <- object
    fit$estimate <- utils::relist(parameters[k, ], skeleton = object$estimate)
    if (is.null(sampler)) {
      # The latent effects of a Gaussian model are drawn independently.
      control <- set_control_mcmc(n_draws = per_draw, burnin = 0, thin = 1)
    } else {
      draws <- tryCatch(sampler(fit$estimate, per_draw), error = leave_out)
      if (is.null(draws)) {
        progress_tick(bar)
        next
      }
      fit$S_samples <- draws$samples
      acceptance[k] <- draws$acceptance
      control <- NULL
    }
    pieces[[k]] <- tryCatch(components_at(fit, control), error = leave_out)
    progress_tick(bar)
  }
  progress_done(bar)
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
  if (!is.null(control_mcmc)) return(control_mcmc$n_draws)
  if (!fit_is_exact(object) && !is.null(object$S_samples)) {
    return(nrow(object$S_samples))
  }
  # As many draws as plug-in prediction takes by default for every family.
  set_control_mcmc()$n_draws
}

#' Draw parameters from the approximate sampling distribution of the estimate
#'
#' Draws are multivariate normal on the scale on which the model is estimated
#' (log variances and log scale), with the covariance matrix of the estimate.
#' Parameters without a finite positive variance, such as those held at a
#' boundary, keep their estimated value, and so do those that the data don't
#' determine (`undetermined_parameters()`), such as a variance of a Gaussian
#' model estimated at zero.
#' @noRd
draw_parameter_values <- function(object, n_draws) {
  estimate <- unlist(object$estimate)
  covariance <- object$covariance
  failure <- function() {
    riskmap_abort(c("The covariance matrix of the estimates isn't available or isn't positive definite, so parameter uncertainty can't be propagated.", "i" = "For a binomial, Poisson or negative binomial model, refitting with {.code control_mcml = set_control_mcml(pc_prior = TRUE)} usually solves this."), class = "riskmap_error_fit")
  }
  if (!nrow(estimated_rows(fit_parameter_table(object)))) {
    abort_input(c("{.arg parameter_uncertainty} can't be {.code TRUE} for a fit that holds every parameter.", "i" = "The parameters given in {.arg fixed_pars} have no uncertainty to propagate."), arg = "parameter_uncertainty")
  }
  if (is.null(covariance)) failure()
  covariance <- as.matrix(covariance)
  if (!identical(dim(covariance), rep(length(estimate), 2L))) failure()
  free <- which(is.finite(diag(covariance)) & diag(covariance) > 0)
  block <- covariance[free, free, drop = FALSE]
  if (!length(free) || any(!is.finite(block))) failure()
  draws <- matrix(estimate, nrow = n_draws, ncol = length(estimate),
                  byrow = TRUE, dimnames = list(NULL, names(estimate)))
  undetermined <- undetermined_parameters(object, block)
  kept <- !undetermined$parameters
  if (!any(kept)) {
    riskmap_abort("The data don't determine any of the estimated parameters, so parameter uncertainty can't be propagated.", class = "riskmap_error_fit")
  }
  # With some parameters kept at their estimates, the others are drawn given
  # them, with the variance of the nugget in place of its ratio to the spatial
  # variance, so that a spatial variance estimated at zero leaves the nugget
  # its own uncertainty.
  if (!all(kept)) {
    warn_undetermined(object, names(estimate)[free][!kept])
    block <- conditional_covariance_of(undetermined$covariance, kept)
  }
  root <- tryCatch(chol((block + t(block)) / 2), error = function(error) NULL)
  if (is.null(root)) failure()
  moves <- matrix(0, nrow = n_draws, ncol = length(free))
  moves[, kept] <- matrix(stats::rnorm(n_draws * sum(kept)), n_draws) %*% root
  if (!all(kept)) moves <- moves %*% t(solve(undetermined$jacobian))
  draws[, free] <- draws[, free, drop = FALSE] + moves
  draws
}

#' The parameters that the data don't determine
#'
#' Those of the covariance, the random effects and the family whose logarithm
#' has a standard error above 3, the nugget by its variance. A Gaussian fit
#' lets the logarithm of a variance estimated at zero drift along an almost
#' flat likelihood, with a standard error of thousands, where a count fit
#' holds it. The information on the logarithm is then below 1/9, about the
#' log-likelihood that a variance of zero loses, and the draws of a normal
#' distribution on the log scale would span more than five orders of
#' magnitude. So would draws of the scale where the likelihood is flat in it.
#' The limit of 3 is a provisional choice. The book's fits are all below 1.
#' @param covariance The covariance of the free parameters of the fit, on the
#'   scale on which it is estimated.
#' @return A list with `parameters`, `TRUE` for each undetermined one, the
#'   `jacobian` to the scale with the variance of the nugget and the
#'   `covariance` there.
#' @noRd
undetermined_parameters <- function(object, covariance) {
  table <- fit_parameter_table(object)
  rows <- vector_rows(table)
  free <- match(colnames(covariance), rows$stored)
  jacobian <- parameter_jacobian(table)[free, free, drop = FALSE]
  reported <- jacobian %*% covariance %*% t(jacobian)
  reported <- (reported + t(reported)) / 2
  list(parameters = rows$block[free] != "coefficients" &
         diag(reported) > 3^2,
       jacobian = jacobian, covariance = reported)
}

#' The covariance of the parameters `kept` given the others
#'
#' The others can have enormous variances and be confounded with each other,
#' so their block is inverted on the scale of their correlations, without the
#' directions in which it is singular.
#' @noRd
conditional_covariance_of <- function(covariance, kept) {
  scale <- sqrt(diag(covariance)[!kept])
  correlation <- covariance[!kept, !kept, drop = FALSE] / outer(scale, scale)
  decomposition <- eigen(correlation, symmetric = TRUE)
  used <- decomposition$values > 1e-10 * max(decomposition$values)
  cross <- sweep(covariance[kept, !kept, drop = FALSE], 2L, scale, "/") %*%
    decomposition$vectors[, used, drop = FALSE]
  covariance[kept, kept, drop = FALSE] -
    cross %*% (t(cross) / decomposition$values[used])
}

#' Warn that parameters keep their estimates when parameters are drawn
#'
#' Not for a variance below a hundredth of the total, which is zero for the
#' prediction, nor for the scale of such a spatial variance or a parameter of
#' the family at its limit, as no warning is given for those the fit holds
#' there. The others have an uncertainty that the prediction leaves out.
#' @param names The parameters, as the estimate names them.
#' @noRd
warn_undetermined <- function(object, names) {
  table <- fit_parameter_table(object)
  rows <- vector_rows(table)
  rows <- rows[match(names, rows$stored), , drop = FALSE]
  values <- unlist(coef(object))
  variances <- values[grepl("^(sigma2|tau2|sigma2_me|sigma2_re\\..*)$",
                            names(values))]
  negligible <- names(variances)[variances < 0.01 * sum(variances)]
  if ("sigma2" %in% negligible) negligible <- c(negligible, "phi")
  named <- setdiff(
    rows$reported[rows$block != "family" | rows$name == "sigma2_me"],
    negligible
  )
  if (length(named)) {
    riskmap_warn(c("The data don't determine {.code {named}}, so {?it keeps its estimate/they keep their estimates} and the prediction leaves out {?its/their} uncertainty.", "i" = "The standard {cli::qty(length(named))}error{?s} on the log scale {?is/are} above 3."), class = "riskmap_warning_fit")
  }
  invisible(NULL)
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
        n_draws = n_samples, burnin = burnin, thin = thin, h = step_size
      ),
      invlink = invlink, start = last_mode, adapt_step = FALSE,
      # At the drawn estimate, so a size the fit estimates is drawn with it.
      family_parameters = fit_family_values(fit)
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
  ratio <- nugget_ratio_natural(fitted_nugget_variance(object, parameters),
                                parameters$sigma2)
  if (ratio == 0) ratio <- 1e-10
  diag(correlation) <- diag(correlation) + ratio
  parameters$sigma2 * correlation
}

#' Average predictive targets over a common reference population
#' @noRd
standardise_prediction_targets <- function(object, prediction_data,
                                           standardise, spatial, offset,
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
    seq_len(nrow(standardise$data)),
    ceiling(seq_len(nrow(standardise$data)) / profiles_per_batch)
  )
  for (batch in batches) {
    profile_data <- plain_data[rep(seq_len(n_locations), times = length(batch)),
                               , drop = FALSE]
    for (name in standardise$predictors) {
      profile_data[[name]] <- rep(standardise$data[[name]][batch],
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
      weight <- standardise$weights[[profile]]
      fixed <- fixed + weight * profile_fixed
      profile_response <- if (fit_is_exact(object)) {
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
  # nonlinear links, its inverse is intentionally not used as the standardised
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

#' Whether each row of prediction data gives a finite design row
#'
#' A row with log(0) or log(-1) doesn't. A term such as `ns(log(x), 2)` stops
#' where `log(x)` isn't finite, so when the design can't be built the rows where
#' a part of a term isn't finite are left out. Every row is taken to be finite
#' when the design can't be built without them either, so that the prediction
#' reports why.
#' @noRd
prediction_design_is_finite <- function(object, predictors) {
  finite_rows <- function(rows) {
    design <- tryCatch(
      suppressWarnings(prediction_design_matrix(
        object, predictors[rows, , drop = FALSE], na_action = na.pass
      )),
      error = function(error) NULL
    )
    if (!is.null(design)) unname(rowSums(!is.finite(design)) == 0)
  }
  every <- rep(TRUE, nrow(predictors))
  finite <- finite_rows(every)
  if (!is.null(finite)) return(finite)
  parts <- prediction_parts_are_finite(object, predictors)
  if (!any(parts)) return(parts)
  finite <- if (!all(parts)) finite_rows(parts)
  if (is.null(finite)) return(every)
  parts[parts] <- finite
  parts
}

#' Whether the parts of the model's terms are finite at each row
#'
#' The parts are the calls and variables inside the variables of the terms,
#' such as `log(x)` in `ns(log(x), 2)`, that give a number for each row.
#' @noRd
prediction_parts_are_finite <- function(object, predictors) {
  formula_parts_are_finite(
    update(interpret_formula(fit_formula(object))$pf, NULL ~ .), predictors
  )
}

#' Build a prediction design matrix using the fitted factor structure
#' @noRd
prediction_design_matrix <- function(object, predictors,
                                     na_action = na.fail) {
  inter_f <- interpret_formula(fit_formula(object))
  fixed_formula <- update(inter_f$pf, NULL ~ .)
  fitted_frame <- model.frame(fixed_formula, data = object$data,
                              na.action = na.fail)
  fixed_terms <- attr(fitted_frame, "terms")
  model_frame <- model.frame(
    fixed_terms, data = predictors, na.action = na_action,
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
    levels <- prediction_setup$re$levels[[term]]
    result <- result +
      prediction_setup$re$samples[[term]][levels, , drop = FALSE]
  }
  result
}

#' The levels of the random effects at prediction rows
#'
#' Each row's level of each `re()` term, as its position among the levels of
#' the fit, in the order of `object$re`, which are the codes of
#' `object$ID_re`. A level that isn't in the fit comes after them.
#' @param values The columns of the `re()` terms, one row per prediction row.
#' @return A list with a vector of positions for each term, `levels`, and the
#'   number of levels beyond those of the fit, `n_new`.
#' @noRd
prediction_random_effect_levels <- function(object, values) {
  levels <- lapply(names(object$re), function(term) {
    known <- match(values[[term]], object$re[[term]])
    missing <- is.na(known)
    labels <- as.character(values[[term]])[missing]
    known[missing] <- length(object$re[[term]]) +
      match(labels, unique(labels))
    known
  })
  names(levels) <- names(object$re)
  n_new <- vapply(seq_along(levels), function(term) {
    max(levels[[term]], 0L) - length(object$re[[term]])
  }, numeric(1))
  list(levels = levels, n_new = pmax(n_new, 0))
}

#' Draw an optional prediction-location nugget component
#'
#' The nugget is one value for each location, so the rows of a joint
#' prediction at one location share their draws, as they share the spatial
#' effect (`location_index`).
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
  location <- prediction_setup$location_index
  n_distinct <- if (is.null(location)) n_locations else max(location)
  nugget <- sample_independent_gaussian(
    mean = rep(0, n_distinct),
    standard_deviation = rep(sqrt(tau2), n_distinct),
    n_samples = n_samples
  )
  if (is.null(location)) nugget else nugget[location, , drop = FALSE]
}

#' Obtain the fitted or fixed nugget variance
#' @noRd
fitted_nugget_variance <- function(object, par_hat = coef(object)) {
  if (!is.null(par_hat$tau2)) return(as.numeric(par_hat$tau2))
  held <- fit_fix_tau2(object)
  if (is.numeric(held) && !isTRUE(held)) {
    return(as.numeric(held))
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
  ordinary <- drop_prediction_operations(summaries, c("cv", "exceedance"))
  requested_names <- c(
    prediction_summary_names("spatial", ordinary),
    if ("random_effects" %in% include)
      prediction_summary_names("random_effects", ordinary),
    if ("nugget" %in% include) prediction_summary_names("nugget", ordinary),
    prediction_summary_names("link", ordinary),
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
  names(result) <- prediction_summary_names(quantity, ordinary)
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
  # No names when only the exceedance probabilities were asked for.
  result <- paste(rep(quantity, length(ordinary)), ordinary, sep = "_")
  if (prediction_operation_requested(summaries, "exceedance")) {
    result <- c(result, paste(
      quantity, "exceedance", prediction_threshold_label(thresholds), sep = "_"
    ))
  }
  result
}

#' Format thresholds for stable column names
#'
#' Each threshold is written on its own, so that its label doesn't depend on
#' the others, with more digits only where two thresholds would have the same
#' label.
#' @noRd
prediction_threshold_label <- function(thresholds) {
  for (digits in c(7L, 15L, 17L)) {
    labels <- vapply(thresholds, format, character(1), digits = digits,
                     scientific = FALSE)
    if (!anyDuplicated(labels)) break
  }
  labels
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

#' The latent effects of a Gaussian model and their observations
#'
#' The location effects, with covariance `spatial_covariance`, and then the
#' levels of each `re()` term, in the order of `object$re`, whose codes are in
#' `object$ID_re`.
#' @return A list with the observation-by-effect `incidence`, the `covariance`
#'   of the effects, the location and group indices of the observations
#'   (`ids`), the number of levels of each term (`n_levels`) and the
#'   `measurement_variance`.
#' @noRd
gaussian_latent_structure <- function(object, spatial_covariance, par_hat) {
  n_observations <- length(object$y)
  indicators <- function(index, n_columns) {
    design <- matrix(0, nrow = n_observations, ncol = n_columns)
    design[cbind(seq_len(n_observations), index)] <- 1
    design
  }
  n_levels <- unname(lengths(object$re))
  random_ids <- if (length(n_levels)) as.matrix(object$ID_re)
  designs <- list(indicators(object$ID_coords, nrow(spatial_covariance)))
  covariances <- list(spatial_covariance)
  for (term in seq_along(n_levels)) {
    designs[[term + 1L]] <- indicators(random_ids[, term], n_levels[term])
    covariances[[term + 1L]] <- diag(
      par_hat$sigma2_re[[term]], nrow = n_levels[term]
    )
  }
  list(
    incidence = do.call(cbind, designs),
    covariance = as.matrix(Matrix::bdiag(covariances)),
    ids = as.matrix(cbind(object$ID_coords, object$ID_re)),
    n_levels = n_levels,
    measurement_variance = par_hat$sigma2_me %||% fit_fix_var_me(object) %||% 0
  )
}

gaussian_conditioning_distribution <- function(object, spatial_covariance,
                                               par_hat, residual) {
  n_spatial <- nrow(spatial_covariance)
  latent <- gaussian_latent_structure(object, spatial_covariance, par_hat)
  incidence <- latent$incidence
  latent_covariance <- latent$covariance
  prediction_weights <- gaussian_prediction_weights(
    Matrix::Matrix(incidence, sparse = TRUE), latent_covariance, latent$ids,
    latent$measurement_variance
  )
  cross_covariance <- latent_covariance %*% t(incidence)
  gain <- prediction_weights(cross_covariance)
  conditional_mean <- as.numeric(gain %*% residual)
  conditional_covariance <- latent_covariance -
    gain %*% incidence %*% latent_covariance
  random_rows <- if (length(latent$n_levels)) {
    split(seq.int(n_spatial + 1L, nrow(latent_covariance)),
          rep(seq_along(latent$n_levels), latent$n_levels))
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

#' Draw the random effects of a Gaussian model given its spatial samples
#'
#' Given the data, the spatial effects at the prediction locations and the
#' random effects are jointly Gaussian, with covariance `-C V^-1 C_Z'`, where
#' `C` and `C_Z` are their covariances with the observations and `V` the
#' covariance of the observations. The random effects are drawn from their
#' distribution given the data and the spatial samples.
#' @param latent The model's latent effects (`gaussian_latent_structure()`).
#' @param prediction_weights The function giving `C V^-1`.
#' @param residual The observations less their fixed effects and offset.
#' @param spatial_samples,spatial_mean The samples of the spatial effects and
#'   their conditional mean.
#' @param spatial_cross_covariance The covariance of the spatial effects with
#'   the observations, `C`.
#' @param spatial_root The factor of the conditional covariance of the spatial
#'   effects.
#' @return A list with a level-by-draw matrix for each `re()` term.
#' @noRd
gaussian_random_effect_samples <- function(latent, prediction_weights,
                                           residual, spatial_samples,
                                           spatial_mean,
                                           spatial_cross_covariance,
                                           spatial_root) {
  n_effects <- ncol(latent$incidence)
  rows <- seq.int(n_effects - sum(latent$n_levels) + 1L, n_effects)
  variance <- diag(latent$covariance)[rows]
  # One row for each level and one column for each observation
  cross_covariance <- variance * t(latent$incidence[, rows, drop = FALSE])
  weights <- prediction_weights(cross_covariance)
  with_spatial <- weights %*% t(spatial_cross_covariance)
  spatial_weights <- cholesky_prediction_weights(with_spatial, spatial_root)
  conditional_covariance <- diag(variance, nrow = length(variance)) -
    weights %*% t(cross_covariance) - spatial_weights %*% t(with_spatial)
  root <- factor_covariance(
    conditional_covariance, "conditional random-effect covariance"
  )
  conditional_mean <- as.numeric(weights %*% residual) -
    spatial_weights %*% (spatial_samples - spatial_mean)
  samples <- sample_correlated_gaussian(
    conditional_mean, t(root), ncol(spatial_samples)
  )
  term <- rep(seq_along(latent$n_levels), latent$n_levels)
  lapply(seq_along(latent$n_levels), function(k) {
    samples[term == k, , drop = FALSE]
  })
}

#' Generate conditional spatial-process samples from aligned prediction data
#'
#' Prediction variables and geometry remain in one row-aligned `sf` object.
#' For areal predictions, `groups` defines the independent joint covariance
#' blocks without splitting the data into parallel lists.
#' A fresh conditional MCMC run of a non-Gaussian model is seeded with
#' `mcmc_seed`, after which the caller's random number state is restored, so
#' the predictive draws that follow start from that state.
#'
#' At the observed locations, and in joint predictions with several rows at some
#' location, the spatial samples have a row for each distinct location, and
#' `location_index` gives the location of each prediction row, so that the rows
#' at one location share its samples (`build_prediction_components()`).
#' Otherwise they have a row for each prediction row. The effect of a level of
#' an `re()` term that isn't in the fit is drawn from its distribution,
#' independently of the data.
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
  inter_f <- interpret_formula(fit_formula(object))
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
      abort_input("{.arg joint} must be {.code TRUE} to predict the random effects.", arg = "joint")
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

  # The distinct location of each prediction row, where the rows at one
  # location have the same draws. A marginal prediction draws each row alone.
  location_index <- if (obs_loc) {
    object$ID_coords
  } else if (type == "joint") {
    distinct <- index_coordinate_rows(grp)
    if (nrow(distinct$coordinates) < n_pred) distinct$index
  }
  if (!obs_loc && !is.null(location_index)) {
    grp <- distinct$coordinates
    groups <- groups[!duplicated(location_index)]
  }
  group_indices <- if (grouped_joint) {
    unname(split(seq_along(groups), groups))
  }

  object$D <- as.matrix(object$D)
  p <- ncol(object$D)

  if (obs_loc) {
    pred_cov_offset <- object$cov_offset %||% 0
  } else if (is.null(inter_f$offset)) {
    pred_cov_offset <- 0
  } else {
    if (!is.numeric(pred_cov_offset)) {
      abort_input("The offset column {.val {inter_f$offset}} of {.arg newdata} must be numeric.", arg = "newdata")
    }
    if (length(pred_cov_offset) != n_pred) {
      riskmap_bug("The offset must have one value for each prediction location.")
    }
  }

  # ---------------------------------------------------------------------------
  # Extract terms from object
  # ---------------------------------------------------------------------------

  par_hat <- coef(object)
  inter_f      <- interpret_formula(fit_formula(object))

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
  # The levels of the random effects at the prediction rows
  # ---------------------------------------------------------------------------
  random_levels <- if (n_re > 0 && !is.null(re_predictors)) {
    if (!is.data.frame(re_predictors) || nrow(re_predictors) != n_pred ||
        ncol(re_predictors) != n_re) {
      riskmap_bug("The random-effect levels must be a data frame with one column per random effect and one row per prediction location.")
    }
    prediction_random_effect_levels(object, re_predictors)
  }

  # ---------------------------------------------------------------------------
  # Spatial quantities
  # ---------------------------------------------------------------------------
  out <- list(mu_pred = mu_pred, grid_pred = grid_pred, par_hat = par_hat)
  distance_scale <- attr(object, "distance_scale") %||% 1
  fitting_coords <- object$coords / distance_scale
  fitting_phi <- par_hat$phi / distance_scale
  fitting_grp <- grp / distance_scale

  conditioning_coords <- if (fit_is_exact(object)) {
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

  if (!fit_is_exact(object) && !obs_loc && !batch_marginal &&
      !grouped_joint) {
    U_pred <- cross_distances(fitting_grp, fitting_coords)
  } else if (fit_is_exact(object) && !obs_loc && !batch_marginal &&
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

  reuse_fitted_samples <- !fit_is_exact(object) &&
    !control_mcmc_supplied && !is.null(object$S_samples)
  n_samples <- if (reuse_fitted_samples) {
    nrow(object$S_samples)
  } else {
    control_mcmc$n_draws
  }

  # ---------------------------------------------------------------------------
  # FIX 2: nu2 / nugget
  # STH: tau2 fixed at 0 -> near-zero for numerical stability
  # LF:  tau2 may be estimated -> use it
  # glgpm: existing behaviour
  # ---------------------------------------------------------------------------
  tau2 <- fitted_nugget_variance(object, par_hat)
  nu2 <- nugget_ratio_natural(tau2, par_hat$sigma2)
  if (nu2 == 0) nu2 <- 1e-10
  diag(R) <- diag(R) + nu2

  diff.y <- if (fit_is_exact(object)) object$y - mu else NULL
  shared_gaussian_re_samples <- NULL

  # ===========================================================================
  # NON-GAUSSIAN MODELS
  # ===========================================================================
  if (!fit_is_exact(object)) {

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
      progress_step(paste0("Running {format(mcmc_iterations(control_mcmc), ",
                           "big.mark = ',')} iterations of the Langevin ",
                           "sampler"),
                    messages)
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
        seed = mcmc_seed, family_parameters = fit_family_values(object)
      )
      progress_step_done(messages)
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
    # Without measurement error or random effects, the observations are the
    # effects of the locations themselves.
    held_variance <- fit_fix_var_me(object)
    observed_exactly <- !is.null(held_variance) && held_variance == 0 &&
      n_re == 0L
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
    } else if (!observed_exactly) {
      latent <- gaussian_latent_structure(object, par_hat$sigma2 * R, par_hat)
      prediction_weights <- gaussian_prediction_weights(
        Matrix(latent$incidence, sparse = TRUE, doDiag = FALSE),
        latent$covariance, latent$ids, latent$measurement_variance
      )
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
  # Samples of the random effects, one matrix for each term
  # ---------------------------------------------------------------------------

  out["re"] <- list(NULL)
  if (!is.null(random_levels)) {
    re_samples <- if (!is.null(shared_gaussian_re_samples)) {
      shared_gaussian_re_samples
    } else if (fit_is_exact(object)) {
      gaussian_random_effect_samples(
        latent, prediction_weights, residual = diff.y,
        spatial_samples = out$S_samples, spatial_mean = mu_cond_S,
        spatial_cross_covariance = C, spatial_root = Sc_root
      )
    } else {
      lapply(seq_len(n_re), function(term) {
        t(simulation$samples[[term + 1L]])
      })
    }
    # The levels that aren't in the fit, independent of the data.
    for (term in which(random_levels$n_new > 0)) {
      n_new <- random_levels$n_new[term]
      re_samples[[term]] <- rbind(re_samples[[term]], sample_independent_gaussian(
        mean = rep(0, n_new),
        standard_deviation = rep(sqrt(par_hat$sigma2_re[[term]]), n_new),
        n_samples = n_samples
      ))
    }
    out$re <- list(levels = random_levels$levels, samples = re_samples)
  }

  out["location_index"] <- list(location_index)
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

##' @title Map a simulation
##' @description
##' Maps one component of one simulation from [simulate_glgpm()] over its
##' prediction grid, by default the mean, with the sampling locations overlaid
##' when data were also simulated. `plot()` of the simulation does the same.
##'
##' @param x The result of [simulate_glgpm()], including a surface.
##' @param simulation A single simulation number to plot, defaulting to one.
##' @param ... For `plot()`, the arguments of `plot_simulation()`. For
##'   `plot_simulation()`, arguments passed to the ggplot2 fill scale, as in
##'   [plot_simulated_surface()].
##' @inheritParams plot_simulated_surface
##'
##' @return A `ggplot` object.
##' @details
##' The colours and the legend follow [plot_prediction()]. A prediction grid
##' that is regular in the CRS of the model is drawn with tiles, and other
##' prediction locations as points. The returned object can be modified further
##' with standard \pkg{ggplot2} functions.
##'
##' @seealso [simulated_surface()], [plot_simulated_surface()]
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:50, ]
##' fit <- glgpm(y ~ gp(), data = example_data, family = "gaussian",
##'              messages = FALSE)
##' grid <- create_grid(create_convex_hull(example_data), spacing = 20)
##' simulation <- simulate_glgpm(fit, nsim = 2, what = c("data", "surface"),
##'                              prediction_grid = grid, seed = 1)
##' plot_simulation(simulation, simulation = 2)
##' # The same map
##' plot(simulation, simulation = 2)
##' @export
plot_simulation <- function(x, simulation = 1, quantity = "mean",
                            palette = NULL, limits = NULL, ...) {
  check_supplied(x)
  check_plot_object(x, "RiskMap_simulation", "plot_simulation",
                    "simulate_glgpm")
  surface <- simulated_surface(x, simulation)
  out <- plot_simulated_surface(surface, quantity = quantity,
                                palette = palette, limits = limits, ...) +
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

##' @rdname plot_simulation
##' @method plot RiskMap_simulation
##' @export
plot.RiskMap_simulation <- function(x, ...) plot_simulation(x, ...)

##' Map values over a grid of points
##'
##' Draws each point of a regular grid as a tile, with the fill scale and legend
##' of the prediction maps of [plot_prediction()]. Other points, such as those of
##' a grid that is regular in another CRS, are drawn as points, as those maps do.
##'
##' @param grid An `sf` or `sfc` object of points.
##' @param values Numeric vector of values to map, one per point in `grid`.
##' @param legend_title Character string used as the fill legend title.
##' @param colours,limits,role,probability,... Passed to
##'   `prediction_fill_scale()`.
##' @return A `ggplot` object.
##' @noRd
plot_grid <- function(grid, values, legend_title, colours, limits, role,
                      probability, ...) {
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

  tiles <- if (tiles_cover(grid_coordinates, tile_width * tile_height)) {
    geom_tile(aes(x = .data$x, y = .data$y, fill = .data$value),
              width = tile_width, height = tile_height)
  } else {
    geom_point(aes(x = .data$x, y = .data$y, fill = .data$value),
               shape = 21, colour = "transparent", size = 1.8)
  }

  ggplot(plot_data) +
    tiles +
    prediction_fill_scale(colours = colours, limits = limits, values = values,
                          breaks = NULL, role = role,
                          probability = probability, ...) +
    # the extent of the tiles, not just their centres
    coord_sf(crs = st_crs(grid),
             xlim = range(plot_data$x) + c(-1, 1) * tile_width / 2,
             ylim = range(plot_data$y) + c(-1, 1) * tile_height / 2) +
    # Set through labs() so that users can replace it with labs(fill = ...).
    labs(x = NULL, y = NULL, fill = legend_title) +
    theme_void() +
    theme(legend.position = "top", legend.justification = "center")
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

##' Whether tiles on the points of a grid are large enough to map it
##'
##' The smallest gap between the coordinates of points that are not on a
##' regular grid is far below the distance between neighbours, so tiles of that
##' size can't be seen.
##'
##' @param coordinates Matrix of the coordinates of the points.
##' @param tile_area The area of a tile.
##' @return `TRUE` or `FALSE`.
##' @noRd
tiles_cover <- function(coordinates, tile_area) {
  coordinates <- index_coordinate_rows(coordinates)$coordinates
  n_points <- nrow(coordinates)
  if (n_points < 2L) return(TRUE)
  # The distance of at most 200 points, spread over the rows, to their nearest
  # neighbours, which is quick on large grids.
  points <- unique(round(seq(1, n_points, length.out = 200L)))
  squared_distance <- vapply(points, function(i) {
    min((coordinates[-i, 1L] - coordinates[i, 1L])^2 +
          (coordinates[-i, 2L] - coordinates[i, 2L])^2)
  }, numeric(1))
  tile_area > 0.25 * stats::median(squared_distance)
}

##' @title Map a simulated surface
##' @description
##' Maps one component of a surface extracted with [simulated_surface()] over
##' its prediction grid, by default the mean. `plot()` of the surface does the
##' same.
##'
##' @param x The result of [simulated_surface()].
##' @param quantity The simulated component to map, one of `"mean"` (the
##'   default), `"linear_predictor"`, `"spatial_effect"`, `"nugget"` and
##'   `"group_effect"`, as in [simulated_values()]. The mean is the prevalence
##'   for a binomial model, the mean per unit exposure for a Poisson or
##'   negative binomial model and the mean of the response for a Gaussian
##'   model. The nugget and the group effect are zero when the model has
##'   none.
##' @param palette Colour palette. `NULL` selects the defaults of
##'   [plot_prediction()], a calm-to-alarm sequential palette for the
##'   mean and a palette centred at zero for the linear predictor and the
##'   effects. A character vector of colours or the name of an HCL palette
##'   accepted by [grDevices::hcl.colors()] can be supplied. Prefix a named
##'   palette with `"-"` to reverse it, for example `"-Viridis"`.
##' @param limits `NULL` (the default) or a numeric vector of two increasing
##'   numbers giving the range of the colour scale, for example to map several
##'   simulations on a common scale.
##' @param ... For `plot()`, the arguments of `plot_simulated_surface()`. For
##'   `plot_simulated_surface()`, arguments passed to the ggplot2 fill scale,
##'   for example `transform = "sqrt"` for skewed values. The map sets
##'   `colours`, `limits`, `breaks`, `labels`, `oob`, `na.value` and `guide`
##'   itself. To change them, add another fill scale to the map with `+`.
##'
##' @return A `ggplot` object.
##' @details
##' The colours and the legend follow [plot_prediction()], and the
##' legend title names the mapped component. A prediction grid that is regular
##' in the CRS of the model is drawn with tiles, and other prediction locations
##' as points. The returned object can be modified further with standard
##' \pkg{ggplot2} functions.
##'
##' @seealso [simulated_surface()], [plot_simulation()]
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:50, ]
##' fit <- glgpm(y ~ gp(), data = example_data, family = "gaussian",
##'              messages = FALSE)
##' grid <- create_grid(create_convex_hull(example_data), spacing = 20)
##' simulation <- simulate_glgpm(fit, nsim = 2, what = c("data", "surface"),
##'                              prediction_grid = grid, seed = 1)
##' surface <- simulated_surface(simulation, simulation = 2)
##' plot_simulated_surface(surface, quantity = "spatial_effect")
##' # The same map
##' plot(surface, quantity = "spatial_effect")
##' @export
plot_simulated_surface <- function(x, quantity = "mean", palette = NULL,
                                   limits = NULL, ...) {
  check_supplied(x)
  check_plot_object(x, "RiskMap_simulated_surface", "plot_simulated_surface",
                    "simulated_surface")
  # The fill scale would reject these less clearly.
  dot_names <- names(match.call(expand.dots = FALSE)$...)
  if ("reverse_palette" %in% dot_names) {
    abort_input("{.arg reverse_palette} is no longer an argument. Prefix a named palette with {.code -} to reverse it, as in {.code palette = \"-Viridis\"}, or give a reversed vector of colours.", arg = "reverse_palette")
  }
  set_by_map <- intersect(dot_names, c("colours", "colors", "limits", "breaks",
                                       "labels", "oob", "na.value", "guide"))
  if (length(set_by_map)) {
    abort_input(c("{.arg {set_by_map}} can't be given, because the map sets {cli::qty(length(set_by_map))}{?it/them} itself.", "i" = "Add another fill scale to the map with {.code +} to change {cli::qty(length(set_by_map))}{?it/them}."), arg = set_by_map[[1L]])
  }
  components <- c("mean", "linear_predictor", "spatial_effect", "nugget",
                  "group_effect")
  # A surface has every component, the nugget and the group effect being zero
  # when the model has none.
  available <- intersect(components, names(x))
  if (!is.character(quantity) || length(quantity) != 1L || is.na(quantity) ||
      !quantity %in% available) {
    supplied <- describe_value(quantity)
    abort_input("{.arg quantity} must be one of {.or {.val {available}}}, not {supplied}.", arg = "quantity")
  }
  values <- x[[quantity]]
  family <- attr(x, "family", exact = TRUE)
  traits <- if (!is.null(family)) family_traits(family)
  role <- if (quantity == "mean") "level" else "signed"
  limits <- check_plot_limits(limits)
  # Default limits end on round values so the legend labels both ends.
  colour_limits <- if (role == "signed") {
    limits %||% (c(-1, 1) * max(abs(pretty(range(values, na.rm = TRUE), n = 4))))
  } else {
    limits %||% range(pretty(range(values, na.rm = TRUE), n = 4))
  }
  colours <- if (is.null(palette)) prediction_role_palette(role) else
    prediction_plot_palette(palette)
  label <- switch(quantity,
                  mean = traits$response_quantity %||% "mean",
                  linear_predictor = "linear predictor",
                  spatial_effect = "spatial effect",
                  nugget = "nugget effect",
                  group_effect = "group effect")
  plot_grid(x,
            values = values,
            legend_title = paste0(toupper(substring(label, 1L, 1L)),
                                  substring(label, 2L)),
            colours = colours,
            limits = colour_limits,
            role = role,
            probability = quantity == "mean" && isTRUE(traits$probability),
            ...)
}

##' @rdname plot_simulated_surface
##' @method plot RiskMap_simulated_surface
##' @export
plot.RiskMap_simulated_surface <- function(x, ...) {
  plot_simulated_surface(x, ...)
}

##' @title Assess Simulations
##'
##' @description Evaluates models using joint data and surface simulations.
##'
##' @param simulation Output from [simulate_glgpm()] with both data and surface,
##'   from a generating model without grouped `re()` effects. The generating
##'   model can have any link.
##' @param models A named list of fitted models returned by [glgpm()] or model
##'   specifications returned by [specify_glgpm()]. Each object defines a
##'   candidate model specification and is refitted to every simulated dataset;
##'   its estimated parameter values are not reused. Parameters that a fit holds
##'   with `fixed_pars` stay held in its refits, as does the size `theta` of a
##'   negative binomial specification. The response of a candidate must be a
##'   plain column name, and its other columns can have any name. When the
##'   models are classified, none can be called `"Class"`.
##' @param spatial_scale The scale(s) at which predictions are assessed:
##'   `"location"`, the locations of the prediction grid, `"area"`, or
##'   `c("location", "area")` to compute both from a single fit-and-predict
##'   pass over the simulations.
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
##' @param aggregate How the targets at the locations within an area are
##'   combined into the target of the area, as in [predict.RiskMap()]: one
##'   of `"mean"`, `"median"`, `"sum"`, `"min"` and `"max"`, or a function
##'   that receives a location-by-draw matrix and the location weights and
##'   returns one number for each draw. Locations have equal weights.
##'   Required only when `spatial_scale` includes `"area"`.
##' @param boundaries An `sf` object containing only POLYGON or MULTIPOLYGON geometries, required if `spatial_scale` includes `"area"`.
##' @param label The name of the column of `boundaries` that holds distinct
##'   names of the areas. If `NULL`, the areas are named `"reg1"`, `"reg2"`
##'   and so on, in a column `"region"`.
##' @param objective A character vector specifying objectives, either `"mse"`, `"classify"`, or both.
##' @param breaks A numeric vector of at least three finite, increasing class
##'   boundaries for classification, as in [plot_prediction()]. A class includes
##'   its upper boundary, the lowest class also includes the first boundary, and
##'   values beyond the boundaries belong to the end classes, so every value has
##'   a class. Required if `objective` includes `"classify"`.
##' @param control_mcmc A control object for MCMC sampling, created with `set_control_mcmc()`. Default is `set_control_mcmc()`.
##' @param seed Optional non-negative integer that makes the assessment
##'   reproducible; see Details. The caller's random number state is restored
##'   on exit. When `NULL` (the default), the seeds of the simulations are
##'   drawn from the session's random number generator.
##' @param messages Logical, whether to report on the assessment, with a
##'   progress bar as described in [riskmap_progress].
##'   Defaults to the option `RiskMap.messages`, or `TRUE` if it is not set.
##'
##' @return A list of class `RiskMap_assess_simulation`:
##' \describe{
##'   \item{objective}{One element per requested `spatial_scale`
##'   (`"location"` and/or `"area"`), each holding `mse`, a matrix with one row per
##'   model and one column per simulation, and/or `classify`, the
##'   classification metrics of each model and simulation.}
##'   \item{fits}{A data frame with one row per model and simulation, giving
##'   the `status` of the fit (`"ok"` or `"failed"`), the error `message` of a
##'   failed fit and any `warnings`. A fit with a warning, such as one that found
##'   no maximum of the likelihood, is `"ok"` and has results.}
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
##' include it in their predictions when the generating model has one. At a
##' prediction location that is also a sampled location, the true nugget is the
##' one of the data, so the prediction is of the spatial effect and that nugget
##' given the data. At any other location the nugget is a new one, as in
##' [predict.RiskMap()], and when areas are assessed the locations with the same
##' coordinates share it. Binomial, Poisson and negative binomial candidates are
##' refitted with the denominators or exposures of the generating model, whatever
##' their own.
##'
##' For classification, the true class of a location or area is the class of its
##' true target, and the predicted class is the most probable one, which holds
##' the largest share of the predictive draws, with ties going to the lowest
##' class. `breaks` divide the truth and the draws in the same way.
##'
##' The target receives the mean rather than the linear predictor, which is
##' what the `transform` argument of [predict.RiskMap()] receives, because
##' models with different links share no linear-predictor scale. Every link is
##' monotone, so a target defined on the linear predictor of a given link can
##' be written as a function of the mean, as [qlogis()] gives the logit scale.
##'
##' The target of an area is `aggregate` applied to the targets at the grid
##' locations that fall within it, with equal weight, unlike the exact
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
##'   spatial_scale = c("location", "area"),
##'   aggregate = "mean",
##'   boundaries = boundary,
##'   objective = "mse",
##'   messages = FALSE
##' )
##'
##' @seealso The methods for the result are
##'   [print()][print.RiskMap_assess_simulation()] and
##'   [summary()][summary.RiskMap_assess_simulation()]. [specify_glgpm()] and
##'   [simulate_glgpm()] make the simulations. [assessment] says which function
##'   answers which question about a model.
##' @export
assess_simulation <- function(simulation,
                              models,
                              spatial_scale,
                              target = identity,
                              aggregate = NULL,
                              boundaries = NULL,
                              label = NULL,
                              objective = c("mse", "classify"),
                              breaks = NULL,
                              control_mcmc = set_control_mcmc(),
                              seed = NULL,
                              messages = getOption("RiskMap.messages", TRUE)) {
  check_supplied(simulation)
  check_supplied(models)
  check_supplied(spatial_scale)

  if (!inherits(simulation, "RiskMap_simulation")) {
    abort_input("{.arg simulation} must be made by {.fn simulate_glgpm}, not {.obj_type_friendly {simulation}}.", arg = "simulation")
  }
  if (!is.character(objective) || length(objective) == 0 ||
      anyDuplicated(objective) ||
      !all(objective %in% c("mse", "classify"))) {
    abort_input("{.arg objective} must be {.val mse}, {.val classify} or both.", arg = "objective")
  }
  if (!is.character(spatial_scale) || length(spatial_scale) == 0 ||
      anyDuplicated(spatial_scale) ||
      !all(spatial_scale %in% c("location", "area"))) {
    abort_input("{.arg spatial_scale} must be {.val location}, {.val area} or both.", arg = "spatial_scale")
  }
  check_control(control_mcmc, "set_control_mcmc")
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)

  want_grid <- "location" %in% spatial_scale
  want_area <- "area" %in% spatial_scale
  want_mse <- "mse" %in% objective
  want_classify <- "classify" %in% objective

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

  sims <- simulation_assessment_data(simulation)
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
    if (is.null(breaks)) {
      abort_input("{.arg breaks} must be given when {.arg objective} includes {.val classify}.", arg = "breaks")
    }
    if (!is.numeric(breaks) || length(breaks) < 3 ||
        anyNA(breaks) || any(!is.finite(breaks)) ||
        any(diff(breaks) <= 0)) {
      abort_input("{.arg breaks} must be at least three finite, strictly increasing numbers.", arg = "breaks")
    }
  }
  n_sim <- length(sims$data_sim)
  n_models <- length(models)

  if (want_area) {
    aggregations <- c("mean", "median", "sum", "min", "max")
    if (is.null(aggregate)) {
      abort_input("{.arg aggregate} must be given when {.arg spatial_scale} includes {.val area}, as one of {.or {.val {aggregations}}} or a function.", arg = "aggregate")
    }
    if (!is.function(aggregate) &&
        (!is.character(aggregate) || length(aggregate) != 1L ||
           is.na(aggregate) || !aggregate %in% aggregations)) {
      supplied <- describe_value(aggregate)
      abort_input("{.arg aggregate} must be one of {.or {.val {aggregations}}} or a function, not {supplied}.", arg = "aggregate")
    }
  }
  model_names <- names(models)
  if (want_classify && "Class" %in% model_names) {
    abort_input("{.arg models} can't have a model called {.val Class}, the name under which the result holds the classes.", arg = "models")
  }
  # Each candidate is refitted to the simulated response, which is in a column
  # of its own, so that a covariate can have any name.
  formulas <- lapply(model_names, function(name) {
    formula <- fit_formula(models[[name]])
    if (!is.symbol(formula[[2L]])) {
      response <- paste(deparse(formula[[2L]]), collapse = " ")
      abort_input(c("The response of the model {.val {name}} must be a plain column name, not {.code {response}}.", "i" = "Every model is refitted to the simulated response itself."), arg = "models")
    }
    update(formula, .riskmap_response ~ .)
  })
  grid <- sims$grid
  include_nugget <- sims$nugget_over_grid

  # Every candidate needs complete covariates at every grid location, so that
  # its predictions line up with the true surface
  for (i in seq_len(n_models)) {
    template_i <- models[[i]]
    if (template_i$family != sims$family) {
      abort_input("The model {.val {model_names[i]}} has the {template_i$family} family, but the simulations have the {sims$family} family.", arg = "models")
    }
    re_names <- if (inherits(template_i, "RiskMap")) names(template_i$re) else
      template_i$re_terms
    required <- setdiff(get_formula_terms(formulas[[i]]),
                        c(".riskmap_response", re_names))
    missing_columns <- setdiff(required, names(grid))
    if (length(missing_columns) > 0) {
      abort_input("The prediction grid of {.arg simulation} has no {cli::qty(length(missing_columns))}column{?s} {.val {missing_columns}}, which the model {.val {model_names[i]}} uses.", arg = "models")
    }
    n_incomplete <- sum(!complete.cases(st_drop_geometry(grid)[required]))
    if (n_incomplete > 0) {
      abort_input(c("The model {.val {model_names[i]}} needs complete values of {.val {required}} on the prediction grid, but {n_incomplete} grid location{?s} {?has/have} missing values.", "i" = "Remove them from the prediction grid before simulating."), arg = "models")
    }
    n_infinite <- sum(!fixed_effects_are_finite(formulas[[i]], grid))
    if (n_infinite > 0) {
      abort_input(c("The model {.val {model_names[i]}} has a covariate that isn't a finite number, such as the logarithm of zero or of a negative number, at {n_infinite} location{?s} of the prediction grid.", "i" = "Remove {cli::qty(n_infinite)}{?it/them} from the prediction grid before simulating."), arg = "models")
    }
  }

  if (want_area) {
    n_reg <- nrow(boundaries)

    if(is.null(label)) {
      boundaries$region <- paste0("reg", 1:n_reg)
      label <- "region"
    } else if (!is.character(label) || length(label) != 1L || is.na(label) ||
               !label %in% names(boundaries)) {
      supplied <- describe_value(label)
      abort_input("{.arg label} must name a column of {.arg boundaries}, not {supplied}.", arg = "label")
    }
    names_reg <- boundaries[[label]]
    if (anyDuplicated(names_reg)) {
      abort_input("The column {.val {label}} of {.arg boundaries} must hold distinct names of areas.", arg = "label")
    }
    boundaries <- st_transform(boundaries, st_crs(grid))
    inter <- st_intersects(boundaries, grid)
    empty <- lengths(inter) == 0
    if (all(empty)) {
      abort_input("No location of the prediction grid of {.arg simulation} falls within the areas of {.arg boundaries}.", arg = "boundaries")
    }
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
    category_labels <- paste0("(", utils::head(breaks, -1), ",",
                              breaks[-1], "]")
    categories_class <- factor(category_labels, levels = category_labels)
  }

  true_target_grid_sim <- apply_target(
    sims$mean_sim, "the simulated true surface"
  )

  # Applies 'aggregate' to the targets of each area, one column per draw, with
  # equal weights for its locations
  summarise_areas <- function(target, context) {
    t(vapply(areas, function(h) {
      locations <- inter[[h]]
      values <- tryCatch(
        aggregate_prediction_values(target[locations, , drop = FALSE],
                                    aggregate, rep(1, length(locations))),
        error = function(e) {
          abort_input("{.arg aggregate} failed for the area {.val {names_reg[h]}} in {context}.", arg = "aggregate", parent = e)
        }
      )
      if (!is.numeric(values) || length(values) != ncol(target) ||
          any(!is.finite(values))) {
        abort_input("{.arg aggregate} must give one finite number for each draw, and didn't for the area {.val {names_reg[h]}} in {context}.", arg = "aggregate")
      }
      as.numeric(values)
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

  out <- list(objective = list())
  if(want_grid) out$objective$location <- init_objective_store()
  if(want_area) out$objective$area <- init_objective_store()

  # Scores of one model on one simulation at one spatial scale. The
  # classification keeps every category in the confusion matrix, including
  # categories absent from this simulation.
  score_targets <- function(true_vals, samples) {
    scores <- list()
    if (want_mse) scores$mse <- mean((rowMeans(samples) - true_vals)^2)
    if (want_classify) {
      scores$classify <- simulation_classification_metrics(
        true_vals, samples, breaks, levels(categories_class)
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
        models[[i]], formulas[[i]], sims$data_sim[[j]], control_mcmc
      )
      refit_args$seed <- fit_seeds[j]
      refit_args$return_samples <- FALSE
      fit_ij <- without_crs_check(do.call(glgpm, refit_args))

      include_ij <- if (include_nugget && fitted_nugget_variance(fit_ij) > 0) {
        "nugget"
      } else {
        character()
      }
      prediction_data <- suppressMessages(
        validate_prediction_data(fit_ij, grid, include_ij)
      )$model_data
      # The response draws are the linear-predictor draws mapped through the
      # refitted model's own inverse link (the identity for a Gaussian fit),
      # so each model reaches the scale of the simulated mean through its own
      # link.
      response_draws <- assessment_response_draws(
        fit_ij, prediction_data, include_ij, control_mcmc,
        groups = if (want_area) groups
      )
      context <- paste0("model '", model_names[i], "', simulation ", j)
      target_ij <- apply_target(response_draws, context)

      scores <- list(status = "ok", message = NA_character_,
                     estimate = estimated_coef(fit_ij))
      if (want_grid) {
        scores$location <- score_targets(true_target_grid_sim[, j], target_ij)
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
        out$objective$location <- store_scores(out$objective$location, i, j,
                                                result_ij$location)
      }
      if (want_area) {
        out$objective$area <- store_scores(out$objective$area, i, j,
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
    if(want_grid) out$objective$location$classify$Class <- categories_class
    if(want_area) out$objective$area$classify$Class <- categories_class
  }
  out$fits <- fits
  out$estimates <- estimates
  out$n_sim <- n_sim
  out$spatial_scale <- spatial_scale
  class(out) <- "RiskMap_assess_simulation"
  stamp_version(out)
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
    family = fit_family_argument(template),
    data = data,
    distance_units = template$distance_units,
    control_mcmc = control_mcmc,
    control_mcml = if (fitted_template) {
      attr(template, "control_mcml") %||% set_control_mcml()
    } else {
      set_control_mcml()
    },
    # A template from specify_glgpm() holds the family's parameters it
    # gives, as a size of the negative binomial family.
    fixed_pars = if (fitted_template) {
      fit_fixed_pars(template)
    } else {
      template$family_parameters
    },
    messages = FALSE
  )
  if (!family_traits(template$family)$exact) {
    args$denominator <- quote(.riskmap_denominator)
    if (fitted_template) {
      args$invlink <- refit_invlink(template)
    } else if (isTRUE(template$custom_link)) {
      args$invlink <- template$invlink
    }
  }
  args
}

##' Predictive draws of the mean of one refit in an assessment of simulations
##'
##' Draws from `prediction_sample_components()`, except at the locations of the
##' grid that are also sampled sites of the fit when the nugget is included.
##' The simulated truth there has the nugget that the observations of the site
##' have, so the spatial effect and the nugget are replaced by draws of the
##' latent effect of the site, which is their sum, given the data. The nugget
##' of any other location is a new one, as in [predict.RiskMap()], which joint
##' draws share among the locations with the same coordinates.
##'
##' @param fit A `RiskMap` fit.
##' @param grid The prediction locations, in the CRS of the fit.
##' @param include `"nugget"`, or `character()` for draws without it.
##' @param control_mcmc The control of the draws.
##' @param groups The blocks within which draws are joint, or `NULL` for
##'   marginal draws.
##' @return A matrix of draws of the mean, one row per location.
##' @noRd
assessment_response_draws <- function(fit, grid, include, control_mcmc,
                                      groups = NULL) {
  coordinates <- coordinates_in_units(grid, fit$distance_units)
  # The index of each distinct point, the sites of the fit first, and the site
  # of each location, NA where the location isn't sampled.
  point <- index_coordinate_rows(rbind(fit$coords, coordinates))$index[
    nrow(fit$coords) + seq_len(nrow(grid))
  ]
  site <- ifelse(point > nrow(fit$coords), NA, point)
  sampled <- if ("nugget" %in% include) which(!is.na(site)) else integer()

  exact <- fit_is_exact(fit)
  if (length(sampled) && !exact) {
    # The conditional samples of the latent effects, which the predictions at
    # the other locations then reuse, so that the draws are joint.
    site_effects <- prepare_prediction(
      fit, fit$data, "marginal", character(), control_mcmc, messages = FALSE,
      observed_locations = TRUE
    )$S_samples
    # The columns of the random effects, which are not predicted, aren't used.
    fit$S_samples <- cbind(t(site_effects),
                           matrix(0, ncol(site_effects), sum(lengths(fit$re))))
    control_mcmc <- NULL
  }
  components <- prediction_sample_components(
    fit, grid, if (is.null(groups)) "marginal" else "joint", include,
    control_mcmc, messages = FALSE, groups = groups
  )
  # Draws that are joint share a nugget among the locations with the same
  # coordinates, as the truth does.
  first <- match(point, point)
  repeated <- if ("nugget" %in% include && !is.null(groups)) {
    which(first != seq_along(first))
  } else {
    integer()
  }
  if (!length(sampled) && !length(repeated)) return(components$response)

  link <- components$link
  link[repeated, ] <- link[repeated, , drop = FALSE] -
    components$nugget[repeated, , drop = FALSE] +
    components$nugget[first[repeated], , drop = FALSE]
  if (exact && length(sampled)) {
    drawn <- gaussian_site_draws(fit, coordinates, groups, sampled,
                                 n_samples = ncol(link))
    site_effects <- drawn$site_effects
    rows <- drawn$rows
    link[rows, ] <- link[rows, , drop = FALSE] + drawn$spatial -
      components$spatial[rows, , drop = FALSE]
  }
  if (length(sampled)) {
    link[sampled, ] <- components$fixed[sampled] + components$offset[sampled] +
      site_effects[site[sampled], , drop = FALSE]
  }
  if (exact) return(link)
  matrix(fit$link_function$inv(as.numeric(link)), nrow = nrow(link))
}

##' Draws of the latent effects of a Gaussian fit at its sites
##'
##' The latent effect of a site is the sum of the spatial effect and the nugget
##' there. Its draws are joint over the sites given the data. With `groups`, the
##' spatial effect of each block that holds a sampled site is drawn again given
##' them, so that the draws remain joint within the block.
##'
##' @param coordinates The prediction locations, in the distance units of the
##'   fit.
##' @param sampled The rows of `coordinates` that are sampled sites.
##' @param n_samples The number of draws.
##' @inheritParams assessment_response_draws
##' @return A list of `site_effects`, one row per site, and of the draws
##'   `spatial` of the spatial effect at `rows`, the rows that share a block
##'   with a sampled site.
##' @noRd
gaussian_site_draws <- function(fit, coordinates, groups, sampled, n_samples) {
  pars <- coef(fit)
  distance_scale <- attr(fit, "distance_scale") %||% 1
  site_coordinates <- fit$coords / distance_scale
  phi <- pars$phi / distance_scale
  covariance <- pars$sigma2 * matern_correlation(
    pairwise_distances(site_coordinates), phi = phi, kappa = fit$kappa,
    return_sym_matrix = TRUE
  )
  diag(covariance) <- diag(covariance) + fitted_nugget_variance(fit, pars)
  residual <- fit$y - as.numeric(as.matrix(fit$D) %*% pars$beta +
                                   fit$cov_offset)
  if (is.null(fit$ID_re)) fit$ID_re <- matrix(0L, length(fit$y), 0L)
  distribution <- gaussian_conditioning_distribution(fit, covariance, pars,
                                                     residual)
  # Without measurement error the data leave the effects of the sites no
  # variance, or only that of the groups, so a jitter is expected.
  measurement_variance <- pars$sigma2_me %||% fit_fix_var_me(fit) %||% 0
  conditional_root <- factor_covariance(
    distribution$covariance[distribution$spatial_rows,
                            distribution$spatial_rows, drop = FALSE],
    "conditional covariance of the sites", warn = measurement_variance > 0
  )
  site_effects <- sample_correlated_gaussian(
    distribution$mean[distribution$spatial_rows], t(conditional_root),
    n_samples
  )
  if (is.null(groups)) {
    return(list(site_effects = site_effects, rows = integer(),
                spatial = matrix(0, 0, n_samples)))
  }

  rows <- which(groups %in% groups[sampled])
  root <- factor_covariance(covariance, "covariance of the sites")
  spatial <- grouped_joint_prediction(
    prediction_coordinates = coordinates[rows, , drop = FALSE] / distance_scale,
    conditioning_coordinates = site_coordinates,
    group_indices = unname(split(seq_along(rows), groups[rows])),
    weight_function = function(cross_covariance) {
      cholesky_prediction_weights(cross_covariance, root)
    },
    conditional_signal = site_effects, marginal_variance = pars$sigma2,
    phi = phi, kappa = fit$kappa, n_samples = n_samples
  )
  list(site_effects = site_effects, rows = rows, spatial = spatial)
}

##' The class of each value
##'
##' Every value has a class. A class includes its upper break, the lowest class
##' also includes the first break, and values beyond the breaks belong to the
##' end classes, as in the maps of [plot_prediction()].
##'
##' @param values A numeric vector or matrix.
##' @param breaks The class boundaries.
##' @return The numbers of the classes, in the shape of `values`.
##' @noRd
simulation_class <- function(values, breaks) {
  classes <- findInterval(values, breaks, left.open = TRUE, all.inside = TRUE)
  dim(classes) <- dim(values)
  classes
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
  true_class <- factor(labels[simulation_class(true_values, breaks)],
                       levels = labels)

  sample_class <- simulation_class(samples, breaks)
  probabilities <- do.call(cbind, lapply(seq_along(labels), function(i) {
    rowMeans(sample_class == i)
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

  list(by_cat = do.call(rbind, by_category),
       overall_cc = mean(true_class == predicted_class))
}

##' @title Print an Assessment of Simulations
##'
##' @description Prints a short description of an assessment from
##' [assess_simulation()]. Use [summary.RiskMap_assess_simulation()] for its
##' results.
##'
##' @param x An object of class `RiskMap_assess_simulation`.
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##'
##' @return Invisibly returns `x`.
##'
##' @method print RiskMap_assess_simulation
##' @export
print.RiskMap_assess_simulation <- function(x, ...) {
  check_dots_unused(...)
  model_names <- unique(x$fits$model %||% character())
  n_models <- length(model_names)
  objectives <- unique(unlist(lapply(x$objective, names)))
  objective_names <- c(mse = "mean squared error",
                       classify = "classification")
  objectives <- ifelse(objectives %in% names(objective_names),
                       objective_names[objectives], objectives)
  scales <- c(location = "at the locations", area = "in areas")[x$spatial_scale]
  models <- if (n_models) {
    paste0(n_models, if (n_models == 1L) " model" else " models", " (",
           paste(model_names, collapse = ", "), ") refitted to ")
  } else {
    "Models refitted to "
  }
  failed <- if (!is.null(x$fits)) sum(x$fits$status != "ok") else 0L
  warned <- sum(!is.na(x$fits$warnings))
  out_print(c(
    out_rule("Assessment of simulations"),
    out_text(paste0(models, out_count(x$n_sim, "simulation"),
                    ", assessed ", paste(scales[!is.na(scales)],
                                         collapse = " and "),
                    " by ", paste(objectives, collapse = " and "), ".")),
    out_bullets(c(
      if (failed > 0L) {
        c("!" = paste0(failed, " of ", nrow(x$fits), " model fits failed."))
      },
      if (warned > 0L) {
        c("!" = paste0(warned, " of ", nrow(x$fits),
                       " model fits gave warnings."))
      }
    )),
    "",
    out_bullets(c(i = paste0(out_fun("summary",
                                     "summary.RiskMap_assess_simulation"),
                             " gives the results.")))
  ))
  invisible(x)
}

##' @title Summarize Simulation Results
##'
##' @description Summarizes the results of model evaluations from a `RiskMap_assess_simulation` object. Provides average metrics for classification by category and overall correct classification (CC) summary.
##'
##' @param object An object of class `RiskMap_assess_simulation`, as returned by `assess_simulation`.
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##'
##' @return A list of class `summary.RiskMap_assess_simulation` with one
##' element per assessed spatial scale (`location` and/or `area`), each holding:
##' - `mse`: a data frame with, for each model, the number of simulations
##'   (`n_sim`), the number with a result (`n_valid`), and the mean and
##'   standard deviation of the MSE across simulations (`MSE_mean`, `MSE_sd`);
##' - `classify`: for each model, `classify_res`, the average sensitivity,
##'   specificity, PPV, NPV and CC by category across simulations, `n_valid`,
##'   the number of simulations with a defined value of each metric, `n_sim`,
##'   and `cc_summary`, the mean, standard deviation and 2.5% and 97.5%
##'   quantiles of the overall proportion of correct classification.
##'
##' It also holds `fits`, the number of model fits (`n`), of failed fits
##' (`failed`) and of fits with warnings (`warned`).
##'
##' @method summary RiskMap_assess_simulation
##' @export
summary.RiskMap_assess_simulation <- function(object, ...) {
  check_dots_unused(...)
  check_class(object, "RiskMap_assess_simulation")

  if (identical(object$n_sim, 1L)) {
    riskmap_warn("The assessment has one simulation, so point estimates are shown without their uncertainty across simulations.")
  }

  # `object$objective` holds one element per requested spatial scale
  # ("location" and/or "area"); each is summarized the same way.
  results <- list()
  for (scale in intersect(c("location", "area"), names(object$objective))) {
    results[[scale]] <- summarise_objective(object$objective[[scale]])
  }
  if (!is.null(object$fits)) {
    results$fits <- c(n = nrow(object$fits),
                      failed = sum(object$fits$status != "ok"),
                      warned = sum(!is.na(object$fits$warnings)))
  }

  # Assign class for S3 print method
  class(results) <- "summary.RiskMap_assess_simulation"
  return(results)
}

##' Summarize one spatial scale's `mse`/`classify` results
##'
##' @param objective The `mse`/`classify` element of a
##'   `RiskMap_assess_simulation` object for a single spatial scale.
##' @return A list with `mse` and/or `classify` summaries.
##' @noRd
summarise_objective <- function(objective) {
  results <- list()

  # Check for "mse" in objective
  if ("mse" %in% names(objective)) {
    mse_data <- objective$mse

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

  # Check for "classify" in objective
  if ("classify" %in% names(objective)) {
    classify_data <- objective$classify

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
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##'
##' @return Invisibly returns `x`.
##'
##' @method print summary.RiskMap_assess_simulation
##' @export
print.summary.RiskMap_assess_simulation <- function(x, ...) {
  check_dots_unused(...)
  failed <- if (!is.null(x$fits) && x$fits[["failed"]] > 0) {
    c("!" = paste0(x$fits[["failed"]], " of ", x$fits[["n"]],
                   " model fits failed, so their results are missing."))
  }
  warned <- if ("warned" %in% names(x$fits) && x$fits[["warned"]] > 0) {
    c("!" = paste0(x$fits[["warned"]], " of ", x$fits[["n"]],
                   " model fits gave warnings, which the `fits` of the ",
                   "assessment hold."))
  }
  scale_label <- c(location = "Locations", area = "Areas")
  scales <- intersect(c("location", "area"), names(x))
  results <- unlist(lapply(scales, function(scale) {
    c(if (length(scales) > 1L) c("", out_rule(scale_label[[scale]])),
      format_objective(x[[scale]]))
  }))
  out_print(c(out_rule("Results of the simulation study"),
              out_bullets(c(failed, warned)),
              results))
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

##' The lines of one spatial scale's `mse`/`classify` summary
##'
##' @param x An element of a `summary.RiskMap_assess_simulation` object, as
##'   returned by `summarise_objective()`.
##' @return A character vector of lines.
##' @noRd
format_objective <- function(x) {
  table <- function(df) out_capture(print_text_table(df))
  mse <- if (!is.null(x$mse)) {
    c("", out_heading("Mean squared error"), table(x$mse))
  }
  classify <- unlist(lapply(names(x$classify), function(model_name) {
    model_data <- x$classify[[model_name]]
    cc <- model_data$cc_summary
    correct <- if (cc$n_valid >= 2L) {
      sprintf(paste0("Mean %.3f, SD %.3f and 95%% simulation interval %.3f ",
                     "to %.3f, from %d of %d simulations."),
              cc$mean, cc$sd, cc$lower, cc$upper, cc$n_valid, cc$n_sim)
    } else {
      sprintf(paste0("Mean %.3f, with uncertainty unavailable from %d of %d ",
                     "simulations."), cc$mean, cc$n_valid, cc$n_sim)
    }
    c("", out_heading(paste("Classification by", model_name)),
      out_text("Averages over the simulations by category"),
      table(model_data$classify_res),
      out_text(paste0("Valid simulations by category and metric, of ",
                      model_data$n_sim)),
      table(model_data$n_valid),
      out_text("Proportion of correct classification across categories"),
      out_text(correct, indent = 2L))
  }))
  c(mse, classify)
}

#' Name of the response variable, used to label prediction maps
#' @noRd
prediction_response_name <- function(object) {
  response <- all.vars(object$formula[[2L]])
  if (length(response)) response[[1L]] else NA_character_
}
