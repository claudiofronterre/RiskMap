# Checks of models against data simulated from them, assess_fit(), and the
# methods of its result. The plots are in predictive_check_plots.R.

##' @title Assess models against data simulated from them
##'
##' @description
##' Checks whether a fitted model could have produced the data. Datasets are
##' simulated from the model at the observed locations, and features of the
##' simulated datasets are compared with the same features of the observed
##' one. The features are
##' - **statistics** of a dataset, such as the standard deviation of the
##'   prevalence or the proportion of zero counts,
##' - the **distribution** of the prevalence, the rate or the response,
##' - the **variogram** of the residuals from the covariates, which shows
##'   whether the model describes the spatial correlation that they leave.
##'
##' All of them are compared by default, and `statistics` says which.
##'
##' @param models A fit from [glgpm()], or a list of fits to the same data. The
##'   names of the list name the models in the results.
##' @param statistics What to compare, the names of statistics, among `"mean"`,
##'   `"sd"`, `"min"`, `"max"`, `"zeros"`, `"distribution"` and `"variogram"`, a
##'   function that takes a dataset, as [simulated_data()] returns it, and gives
##'   named numbers, or a list of names and functions. The default, `NULL`, is
##'   numbers that depend on the family (see Details), the distribution and the
##'   variogram.
##' @param breaks Optional breakpoints of the distance classes of the variogram,
##'   in the distance units of the (first) model, as in [variogram()], whose
##'   default is used otherwise, 14 classes of equal width up to a third of the
##'   largest distance between locations. An error if `statistics` doesn't include
##'   the variogram.
##' @param nsim Number of datasets simulated from each model, at least 19.
##'   Defaults to 1,000.
##' @param seed Optional non-negative integer seed. When supplied, the caller's
##'   random-number state is restored on exit. All models are simulated with
##'   the same seed, so that their surfaces share random numbers.
##' @param return_simulations Logical, whether the result holds the simulations
##'   of each model, from [simulate_glgpm()]. Defaults to `FALSE`.
##' @param messages Logical, whether to show a progress bar, as described in
##'   [riskmap_progress]. Defaults to the option `RiskMap.messages`, or `TRUE`
##'   if it is not set.
##'
##' @return A list of class `RiskMap_predictive_check`, containing
##' \describe{
##'   \item{observed}{A data frame with the observed `value` of each `statistic`
##'     that is a number for each `model`.}
##'   \item{simulated}{A data frame with the `value` of each `statistic` in each
##'     `simulation` from each `model`.}
##'   \item{distribution}{A data frame with the distribution function of the
##'     observed data for each `model`, at the values `value`, the `median` of
##'     those of the simulated datasets, the limits of the joint envelope,
##'     `lower` and `upper`, and whether the observed one is `outside` it. `NULL`
##'     if not in `statistics`.}
##'   \item{variogram}{A data frame with the same columns for the semivariance of
##'     the residuals in each distance class, the mean `distance` of its pairs of
##'     locations and their number, `n_pairs`. `NULL` if not in `statistics` or
##'     when no pair of locations is in a class, and without the models whose
##'     covariates leave no residuals.}
##'   \item{nsim, seed, breaks, distance_units, level}{The settings. `level` is
##'     the coverage of the joint envelopes, 0.95.}
##'   \item{family, response, quantity, unit}{The family of each model, the
##'     response of the formula, the quantity whose distribution is calculated,
##'     `"prevalence"`, `"rate"` or `"response"`, and whether each value of it is
##'     that of an `"observation"` or of a `"location"`.}
##'   \item{n_locations, n_observations}{The numbers of locations and of
##'     observations.}
##'   \item{simulations}{The simulations of each model, a named list of results
##'     of [simulate_glgpm()], if `return_simulations = TRUE`.}
##' }
##'
##' @details
##' # The simulated datasets
##'
##' Each dataset has a new surface, new values of the nugget and the random
##' effects, and new responses at the observed locations, with the observed
##' covariates, offsets and denominators. The parameters stay at their
##' estimates. [simulate_glgpm()] simulates them.
##'
##' The surface is simulated again, and not predicted from the data, so that the
##' model must reproduce a feature through its structure. A surface that is
##' conditioned on the data has adapted to them. Datasets simulated from it
##' resemble the observed one more than the model itself does, so a check of them
##' can miss what is wrong with the model.
##'
##' # Statistics
##'
##' The numbers with a name are the mean, the standard deviation, the minimum
##' and the maximum of the empirical prevalence, the rate or the Gaussian
##' response, and the proportion of zeros. The default ones are `"mean"`, `"sd"`
##' and `"zeros"`, and for Gaussian models `"mean"`, `"sd"`, `"min"` and `"max"`.
##' In the result the names say what the statistic is of, such as
##' `sd_prevalence`.
##'
##' These numbers are those of the observations. The exception is a binomial
##' model whose denominators are all 1, because the prevalence of a single
##' individual is only 0 or 1. The positives and the denominators are then added
##' up at each location, and the numbers are those of the locations. With
##' other denominators the observations are not added up, because that would
##' hide the variation between those at the same location.
##'
##' A function can calculate anything from a dataset, such as a correlation with
##' a covariate, and can return a vector of named numbers, a list of them or a
##' data frame with one row. The names of a list are those of its statistics, so
##' that a function in it that returns one number has the name that it has in the
##' list. Simulated datasets in which a statistic isn't finite are left out of its
##' summary, with a warning.
##'
##' # Distribution and variogram
##'
##' These are curves and not numbers, and each has a plot of its own. They are
##' calculated by default, and with a given `statistics` only if they are in it.
##'
##' The distribution is that of the same values as the numbers with a name,
##' the prevalence, the rate or the response of each observation, or the
##' prevalence of each location when the denominators of a binomial model are all
##' 1.
##'
##' The variogram is that of the residuals from the covariates. The response is
##' first transformed, to the empirical logit for binomial models and to the
##' logarithm of the rate with 0.5 added to the count for models of counts,
##' whatever their link. The response of Gaussian models is not transformed. There
##' is one value for each location, where the observations are added up, or
##' averaged for Gaussian models. The residuals are the differences between these
##' values, less any offset, and their least squares fit on the covariates of the
##' model. An offset under a link other than the canonical one is one of the
##' covariates of this fit. The regression is done again for each simulated
##' dataset. The variogram of
##' the values themselves would be dominated by a strong trend, which the
##' covariates explain, and would show little of the spatial correlation. Each
##' model has its own residuals and so its own observed variogram.
##'
##' The joint envelope is drawn as in [variogram()], from the observed and the
##' simulated curves ordered by their extreme rank lengths, at a level of 95%.
##' With the parameters at their true values, the whole observed curve is inside
##' it with a probability of at least 95%.
##'
##' # What a check can tell
##'
##' A feature far outside what the model simulates is one that the model misses.
##' The features that tell most are those that the model has no parameter for. A
##' feature that the model was fitted to, such as the mean, can hardly fail,
##' whatever is wrong with the rest of the model.
##'
##' A model that passes is not shown to be right. The parameters are estimated
##' from the same data, so the simulated datasets are centred near the observed
##' one, and the check is conservative. A correct model is outside an envelope
##' less often than 5% of the time. With few locations, or when the spatial
##' correlation extends over much of the study area, datasets simulated from the
##' same model differ a lot from each other, and the check has little power.
##'
##' The edge of an envelope depends on the simulations, so that another `seed`
##' moves it, and a single point near the edge, on either side of it, says
##' little. The first distance classes often have few pairs, so a single point
##' outside the envelope there needs care.
##'
##' Nothing here is a p-value or a mark that a model passes or fails. An observed
##' statistic is unusual for the model when it is outside the interval that holds
##' 95% of the simulated ones. The proportion of simulated values at or above the
##' observed one says where it lies among them. It says little when many
##' simulated values are equal to the observed one. When the data have no zeros,
##' for example, the proportion of zeros of every simulated dataset is at or
##' above theirs, also for a model that seldom gives any.
##'
##' @seealso [plot_predictive_check()], which `plot()` of the result also calls,
##'   [summary()][summary.RiskMap_predictive_check()], [to_table()], and
##'   [simulate_glgpm()] for the simulations. [assessment] says which function
##'   answers which question about a model.
##'
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:100, ]
##' models <- list(
##'   intercept = glgpm(y ~ gp(), data = example_data, family = "gaussian",
##'                     messages = FALSE),
##'   covariate = glgpm(y ~ log(pop_dens) + gp(), data = example_data,
##'                     family = "gaussian", messages = FALSE)
##' )
##'
##' result <- assess_fit(models, nsim = 199, seed = 1)
##' summary(result)
##' plot_predictive_check(result, type = "variogram")
##'
##' # Numbers by name, one of one's own, and the variogram
##' describe <- function(data) {
##'   c(correlation = cor(data$y, log(data$pop_dens)))
##' }
##' result <- assess_fit(models, list("sd", describe, "variogram"), nsim = 199,
##'                      seed = 1)
##' plot_predictive_check(result)
##'
##' @export
assess_fit <- function(models, statistics = NULL, breaks = NULL, nsim = 1000,
                         seed = NULL, return_simulations = FALSE,
                         messages = getOption("RiskMap.messages", TRUE)) {
  check_supplied(models)
  models <- read_assessed_models(models)
  quantity <- check_same_quantity(models)
  statistics <- read_statistics(statistics, quantity)
  with_variogram <- "variogram" %in% statistics$curves
  if (!is.null(breaks)) {
    check_distance_breaks(breaks)
    if (!with_variogram) {
      abort_input("{.arg breaks} are the distance classes of the variogram, which {.arg statistics} doesn't include.", arg = "breaks")
    }
  }
  check_positive_integer(nsim, "nsim")
  # The fewest datasets that a 95% envelope can be drawn from.
  if (nsim < 19) {
    abort_input("{.arg nsim} must be at least 19 for a 95% envelope, not {nsim}.", arg = "nsim")
  }
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)
  check_logical(return_simulations)
  check_logical(messages)

  first <- models[[1L]]
  records <- check_records(first, quantity)
  variogram <- if (with_variogram) {
    check_variogram_classes(records, breaks, first$distance_units)
  }
  if (!is.null(seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
  }
  # One seed for all models, whose surfaces then share random numbers.
  simulation_seed <- seed %||% sample.int(.Machine$integer.max, 1L)

  bar <- progress_bar("Simulating and describing datasets",
                      length(models) * nsim, messages)
  described <- Map(function(fit, name) {
    describe_simulations(fit, name, statistics, records, nsim, simulation_seed,
                         !is.null(seed), return_simulations, bar)
  }, models, names(models))
  progress_done(bar)

  values <- lapply(described, `[[`, "statistics")
  warn_not_finite(values)
  out <- list(
    observed = observed_statistics(values),
    simulated = simulated_statistics(values),
    distribution = if ("distribution" %in% statistics$curves) {
      check_distributions(lapply(described, `[[`, "values"))
    },
    variogram = check_variograms(lapply(described, `[[`, "residuals"),
                                 variogram$classes),
    nsim = nsim, seed = seed, breaks = variogram$breaks,
    distance_units = first$distance_units, level = check_level(),
    family = vapply(models, `[[`, character(1), "family"),
    response = paste(deparse(first$formula[[2L]]), collapse = ""),
    quantity = quantity, unit = records$unit,
    statistic_labels = statistic_labels(statistics$numbers, quantity),
    types = c(if (length(statistics$numbers)) "statistic", statistics$curves),
    n_locations = records$n_locations, n_observations = length(first$y)
  )
  if (return_simulations) {
    out$simulations <- lapply(described, `[[`, "simulation")
  }
  class(out) <- "RiskMap_predictive_check"
  stamp_version(out)
}

##' What there is to draw of a predictive check
##'
##' What was asked for, without a variogram that the data don't have.
##' @param x The check.
##' @noRd
drawn_types <- function(x) {
  x$types[x$types != "variogram" | !is.null(x$variogram)]
}

##' Why a predictive check doesn't have the variogram that was asked for
##' @param breaks The breakpoints of its distance classes.
##' @noRd
missing_variogram <- function(breaks) {
  if (is.null(breaks)) {
    "no pair of locations is in a distance class"
  } else {
    "the covariates leave no residuals"
  }
}

##' The coverage of the joint envelopes of a predictive check
##' @noRd
check_level <- function() 0.95

##' The quantity whose distribution a family describes
##' @noRd
check_quantity <- function(family) {
  family_traits(family)$response_quantity %||% "response"
}

##' The quantity of the models, which must be the same for all of them
##' @noRd
check_same_quantity <- function(models, call = rlang::caller_env()) {
  quantities <- vapply(models, function(fit) check_quantity(fit$family),
                       character(1))
  if (length(unique(quantities)) > 1L) {
    families <- unique(vapply(models, `[[`, character(1), "family"))
    abort_input(c("The models in {.arg models} must be all for prevalences, all for rates or all Gaussian, so that they share one distribution.", "x" = "Their families are {.val {families}}."), arg = "models", call = call)
  }
  quantities[[1L]]
}

##' The statistics that `assess_fit()` has names for, which are numbers or curves
##' @noRd
number_names <- function() c("mean", "sd", "min", "max", "zeros")

##' @rdname number_names
##' @noRd
curve_names <- function() c("distribution", "variogram")

##' @rdname number_names
##' @noRd
statistic_names <- function() c(number_names(), curve_names())

##' Read `statistics` as a list of names and functions, and the curves
##'
##' @param quantity The quantity of the models, which sets the default.
##' @return A list with the `numbers`, a list of the names of statistics and the
##'   functions, whose names are the names given by the user, or empty, and the
##'   names of the `curves`.
##' @noRd
read_statistics <- function(statistics, quantity,
                            call = rlang::caller_env()) {
  if (is.null(statistics)) {
    numbers <- if (quantity == "response") {
      c("mean", "sd", "min", "max")
    } else {
      c("mean", "sd", "zeros")
    }
    statistics <- c(numbers, curve_names())
  }
  if (is.function(statistics)) statistics <- list(statistics)
  supplied <- describe_value(statistics)
  is_name <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
  usable <- (is.list(statistics) || is.character(statistics)) &&
    length(statistics)
  statistics <- as.list(statistics)
  if (!usable || !all(vapply(statistics, function(x) {
    is.function(x) || is_name(x)
  }, logical(1)))) {
    abort_input(c("{.arg statistics} must be the names of statistics, a function, or a list of names and functions, not {supplied}.", "i" = "The names are {.val {statistic_names()}}.", "i" = "A function takes a dataset and gives named numbers."), arg = "statistics", call = call)
  }
  named <- unlist(Filter(is_name, statistics), use.names = FALSE)
  unknown <- setdiff(named, statistic_names())
  if (length(unknown)) {
    abort_input(c("{.arg statistics} must name statistics among {.val {statistic_names()}}, and {.val {unknown}} {?isn't one/aren't among them}.", "i" = "Give a function of a dataset for another statistic."), arg = "statistics", call = call)
  }
  names(statistics) <- names(statistics) %||% rep("", length(statistics))
  names(statistics)[is.na(names(statistics))] <- ""
  is_curve <- vapply(statistics, function(x) {
    is.character(x) && x %in% curve_names()
  }, logical(1))
  renamed <- unlist(statistics[is_curve & nzchar(names(statistics))],
                    use.names = FALSE)
  if (length(renamed)) {
    abort_input("{.val {renamed}} can't have another name in {.arg statistics}, because {?it is/they are} drawn in {?a plot of its/plots of their} own.", arg = "statistics", call = call)
  }
  check_distinct_statistics(character(),
                            statistic_labels(statistics[!is_curve], quantity),
                            call)
  list(numbers = statistics[!is_curve],
       curves = intersect(curve_names(), unlist(statistics[is_curve])))
}

##' The observations, their locations and denominators, for a predictive check
##'
##' The prevalence of one person is only 0 or 1, and so are its distribution and
##' statistics. With denominators that are all 1, those of the locations are
##' compared.
##' @return A list with the `quantity`, the `unit` whose values are compared,
##'   the `index` of the location of each observation, the number of locations,
##'   `n_locations`, those with a positive denominator, `used`, the `totals` of
##'   their denominators, their `coordinates`, the `denominators` of the
##'   observations and their `weights` in their location.
##' @noRd
check_records <- function(fit, quantity) {
  indexed <- index_coordinate_rows(
    coordinates_in_units(fit$data, fit$distance_units)
  )
  denominators <- assessed_denominators(fit)
  totals <- as.numeric(rowsum(denominators, indexed$index))
  # A location where nobody was examined has no prevalence or rate.
  used <- totals > 0
  individual <- quantity == "prevalence" &&
    all(denominators[denominators > 0] == 1)
  list(quantity = quantity,
       unit = if (individual) "location" else "observation",
       index = indexed$index, n_locations = length(totals),
       used = used, totals = totals[used],
       coordinates = indexed$coordinates[used, , drop = FALSE],
       denominators = denominators,
       weights = denominators / totals[indexed$index])
}

##' The distance classes of the variogram of a predictive check
##'
##' @return A list with the `breaks` and the `classes`, which are `NULL` when no
##'   pair of locations is in a class.
##' @noRd
check_variogram_classes <- function(records, breaks, distance_units,
                                    call = rlang::caller_env()) {
  n <- nrow(records$coordinates)
  none <- list(breaks = NULL, classes = NULL)
  d <- if (n > 1L) as.numeric(pairwise_distances(records$coordinates))
  if (!is.null(breaks)) {
    if (is.null(d)) {
      abort_input("{.arg breaks} can't be given for data at one location, which have no variogram.", arg = "breaks", call = call)
    }
  } else {
    breaks <- if (!is.null(d)) default_distance_breaks(d)
    if (is.null(d) || min(d) > max(breaks)) {
      riskmap_warn(c("No pair of locations is within a third of the largest distance between them, so there is no variogram.", "i" = "Give the distance classes in {.arg breaks}."))
      return(none)
    }
  }
  list(breaks = breaks,
       classes = distance_classes(d, n, breaks, distance_units,
                                  "Check {.arg breaks}.", call = call))
}

##' Simulate from one model and describe each dataset
##'
##' @param statistics The result of `read_statistics()`.
##' @param simulation_seed The seed of the simulations.
##' @param reset Whether the functions of `statistics` start from the seed, so
##'   that any random numbers they use are the same for every model.
##' @return A list with the `statistics`, the `values` whose distribution is
##'   compared and the `residuals` whose variogram is, `NULL` if there are none.
##'   Each is a matrix with a dataset in each column, the observed one first.
##'   `simulation` is the simulation, if asked for.
##' @noRd
describe_simulations <- function(fit, name, statistics, records, nsim,
                                 simulation_seed, reset, return_simulations,
                                 bar) {
  simulation <- simulate_glgpm(fit, nsim = nsim, what = "data",
                               seed = simulation_seed)
  responses <- cbind(fit$y, matrix(simulation$samples$data[, , "response"],
                                   nrow = length(fit$y)))
  if (reset) set.seed(simulation_seed)
  described <- describe_locations(responses, records)
  out <- list(values = described$values)
  if ("variogram" %in% statistics$curves) {
    out$residuals <- residuals_from_covariates(described$field, fit, records)
    if (is.null(out$residuals)) {
      riskmap_warn("Model {.val {name}} has as many coefficients as there are locations, so it has no residuals and no variogram.")
    }
  }
  out$statistics <- calculate_statistics(
    statistics$numbers, out$values, responses, simulation$locations$data,
    simulation$model$response, records$quantity, name, bar
  )
  if (return_simulations) out$simulation <- simulation
  out
}

##' The values whose distribution is compared, and the responses transformed for
##' the variogram
##'
##' @param responses A matrix with an observation in each row and a dataset in
##'   each column.
##' @return A list with the `values`, and the `field`, with a location in each
##'   row.
##' @noRd
describe_locations <- function(responses, records) {
  totals <- rowsum(responses, records$index)[records$used, , drop = FALSE]
  dimnames(totals) <- NULL
  if (records$quantity == "response") {
    # The denominators of Gaussian observations are 1, so the totals are the
    # numbers of observations.
    return(list(values = responses, field = totals / records$totals))
  }
  values <- if (records$unit == "location") {
    totals / records$totals
  } else {
    used <- records$denominators > 0
    responses[used, , drop = FALSE] / records$denominators[used]
  }
  field <- if (records$quantity == "rate") {
    log((totals + 0.5) / records$totals)
  } else {
    log((totals + 0.5) / (records$totals - totals + 0.5))
  }
  list(values = values, field = field)
}

##' Residuals from the least squares fit on the covariates of a model
##'
##' One projection for all the datasets, so that it costs a matrix product. A
##' covariate that varies between the observations at a location has its mean
##' there, weighted by the denominators. An offset is subtracted. Under a link
##' other than the canonical one, it isn't on the scale of the field, and is one
##' of the covariates instead.
##' @param field A matrix with a location in each row and a dataset in each
##'   column.
##' @return A matrix of the size of `field`, or `NULL` when the covariates leave
##'   no residuals.
##' @noRd
residuals_from_covariates <- function(field, fit, records) {
  at_locations <- function(rows) {
    rows <- rowsum(as.matrix(rows) * records$weights, records$index)
    unname(rows[records$used, , drop = FALSE])
  }
  design <- at_locations(fit$D)
  offset <- fit$cov_offset
  if (!is.null(offset) && any(offset != 0)) {
    canonical <- is.null(fit$link_function) ||
      identical(fit$link_function$name, "canonical")
    if (!canonical) {
      design <- cbind(design, at_locations(offset))
    } else if (records$quantity == "rate") {
      # The counts at a location are added up, and so are their exposures.
      largest <- stats::ave(offset, records$index, FUN = max)
      field <- field - as.numeric(at_locations(largest)) -
        log(as.numeric(at_locations(exp(offset - largest))))
    } else {
      field <- field - as.numeric(at_locations(offset))
    }
  }
  if (!ncol(design)) return(field)
  # The columns of a polynomial of a year, say, are nearly collinear.
  decomposition <- qr(design, tol = 1e-12)
  if (decomposition$rank >= nrow(design)) return(NULL)
  qr.resid(decomposition, field)
}

##' Calculate the statistics, in the order in which they were given
##'
##' @param values A matrix with a location or an observation in each row and a
##'   dataset in each column.
##' @return A matrix with a statistic in each row and a dataset in each column.
##' @noRd
calculate_statistics <- function(statistics, values, responses, data, response,
                                 quantity, name, bar,
                                 call = rlang::caller_env()) {
  if (!length(statistics)) {
    progress_tick(bar, ncol(responses) - 1L)
    return(matrix(numeric(), 0L, ncol(responses)))
  }
  own <- vapply(statistics, is.function, logical(1))
  if (any(own)) {
    given <- apply_statistics(statistics[own], responses, data, response, name,
                              statistic_labels(statistics, quantity), bar, call)
    item <- which(own)[given$item]
  } else {
    progress_tick(bar, ncol(responses) - 1L)
    given <- item <- NULL
  }
  named <- named_statistics(statistics[!own], values, quantity)
  rbind(given$values, named)[order(c(item, which(!own))), , drop = FALSE]
}

##' Check that the statistics have different names and different labels
##'
##' The tables and the plots tell the statistics by their labels.
##' @param own The names of the statistics of the functions, which are also their
##'   labels.
##' @param labels The labels of the statistics given by name, from
##'   `statistic_labels()`.
##' @noRd
check_distinct_statistics <- function(own, labels, call) {
  repeated <- function(x) unique(x[duplicated(x)])
  same <- repeated(c(names(labels), own))
  if (!length(same)) same <- repeated(c(unname(labels), own))
  if (length(same)) {
    abort_input("The statistics must have different names, and {.val {same}} {?is/are} repeated.", arg = "statistics", call = call)
  }
  invisible(NULL)
}

##' Calculate the statistics given by name
##'
##' @param statistics A list of names of statistics, whose own names are the
##'   names given by the user, or empty.
##' @return As `calculate_statistics()`.
##' @noRd
named_statistics <- function(statistics, values, quantity) {
  if (!length(statistics)) return(NULL)
  mean <- colMeans(values)
  rows <- lapply(statistics, function(statistic) {
    switch(statistic, mean = mean,
           sd = sqrt(colSums((values - rep(mean, each = nrow(values)))^2) /
                       (nrow(values) - 1)),
           min = apply(values, 2L, min), max = apply(values, 2L, max),
           zeros = colMeans(values == 0))
  })
  out <- do.call(rbind, unname(rows))
  rownames(out) <- names(named_statistic_labels(statistics, quantity))
  out
}

##' The labels of the statistics given by name
##'
##' @return The labels, whose names are those of the statistics in the result.
##' @noRd
named_statistic_labels <- function(statistics, quantity) {
  key <- unlist(statistics, use.names = FALSE)
  name <- key
  label <- c(mean = "Mean", sd = "SD", min = "Minimum", max = "Maximum",
             zeros = "Proportion of zeros")[key]
  described <- quantity != "response" & key != "zeros"
  name[described] <- paste(key, quantity, sep = "_")[described]
  label[described] <- paste0(label, ifelse(key == "sd", " of ", " "),
                             quantity)[described]
  # A name given by the user is kept as it is.
  own <- nzchar(names(statistics))
  name[own] <- label[own] <- names(statistics)[own]
  stats::setNames(label, name)
}

##' The labels of those statistics that `assess_fit()` has names for
##' @noRd
statistic_labels <- function(statistics, quantity) {
  named <- Filter(is.character, statistics)
  if (length(named)) named_statistic_labels(named, quantity) else character()
}

##' Calculate the statistics given by the user
##'
##' @param data The data of the model, whose response is replaced.
##' @param response The name of the column of the response.
##' @param labels The labels of the statistics given by name.
##' @return A list with the `values`, a matrix with a statistic in each row and a
##'   dataset in each column, and the `item` of `statistics` that gave each row.
##' @noRd
apply_statistics <- function(statistics, responses, data, response, name,
                             labels, bar, call = rlang::caller_env()) {
  if (is.null(response) || length(response) != 1L ||
      !identical(make.names(response), response)) {
    abort_input(c("With a function in {.arg statistics}, the response of the model formulas must be a plain column name.", "i" = "Put a transformed response in a column of the data before fitting the models."), arg = "models", call = call)
  }
  warned <- NULL
  item <- NULL
  evaluate <- function(k) {
    data[[response]] <- as.numeric(responses[, k])
    dataset <- if (k == 1L) "the observed data" else
      paste("simulated dataset", k - 1L)
    withCallingHandlers(
      pieces <- unname(Map(function(statistic, label) {
        value <- tryCatch(statistic(data), error = function(error) {
          abort_input("{.arg statistics} failed on {dataset} of model {.val {name}}.", arg = "statistics", call = call, parent = error)
        })
        read_statistic_values(value, label, dataset, name, call)
      }, statistics, names(statistics))),
      warning = function(warning) {
        if (inherits(warning, "riskmap_warning")) return()
        warned <<- c(warned, conditionMessage(warning))
        tryInvokeRestart("muffleWarning")
      }
    )
    if (k == 1L) item <<- rep(seq_along(pieces), lengths(pieces))
    unlist(pieces)
  }
  observed <- evaluate(1L)
  check_distinct_statistics(names(observed), labels, call)
  out <- matrix(observed, nrow = length(observed), ncol = ncol(responses),
                dimnames = list(names(observed), NULL))
  for (k in seq_len(ncol(responses))[-1L]) {
    value <- evaluate(k)
    if (!identical(names(value), names(observed))) {
      abort_input("{.arg statistics} must give the same named numbers for every dataset, and those of simulated dataset {k - 1L} of model {.val {name}} differ from those of the observed data.", arg = "statistics", call = call)
    }
    out[, k] <- value
    progress_tick(bar)
  }
  if (length(warned)) {
    riskmap_warn(c("{.arg statistics} gave {length(warned)} warning{?s} for the datasets of model {.val {name}}.", "i" = "The first was {.val {warned[1L]}}."))
  }
  list(values = out, item = item)
}

##' Read what a function of `statistics` returned as named numbers
##' @noRd
read_statistic_values <- function(value, label, dataset, name, call) {
  if (is.data.frame(value) && nrow(value) != 1L) value <- NULL
  if (is.list(value)) {
    single <- all(lengths(value) == 1L)
    value <- if (single) unlist(value)
  }
  if (is.logical(value)) value <- value + 0
  if (length(value) == 1L && nzchar(label)) names(value) <- label
  if (length(value) == 1L && is.null(names(value))) names(value) <- "statistic"
  named <- !is.null(names(value)) && !anyNA(names(value)) &&
    all(nzchar(names(value)))
  if (!is.numeric(value) || !length(value) || !named) {
    abort_input(c("{.arg statistics} must give named numbers, such as {.code c(mean = 0.2, sd = 0.1)}.", "x" = "It didn't for {dataset} of model {.val {name}}."), arg = "statistics", call = call)
  }
  stats::setNames(as.numeric(value), names(value))
}

##' Warn once about the statistics that aren't finite in some datasets
##'
##' @param values A named list of matrices of statistics, one for each model.
##' @noRd
warn_not_finite <- function(values) {
  if (!nrow(values[[1L]])) return(invisible(NULL))
  counts <- do.call(rbind, lapply(names(values), function(name) {
    data.frame(model = name, statistic = rownames(values[[name]]),
               observed = !is.finite(values[[name]][, 1L]),
               n = rowSums(!is.finite(values[[name]][, -1L, drop = FALSE])))
  }))
  if (any(counts$observed)) {
    statistics <- unique(counts$statistic[counts$observed])
    riskmap_warn("{cli::qty(length(statistics))}The statistic{?s} {.val {statistics}} {?isn't/aren't} finite in the observed data.")
  }
  some <- counts[counts$n > 0, ]
  if (nrow(some)) {
    nsim <- ncol(values[[1L]]) - 1L
    most <- some[which.max(some$n), ]
    statistics <- unique(some$statistic)
    riskmap_warn(c("{cli::qty(length(statistics))}The statistic{?s} {.val {statistics}} {?isn't/aren't} finite in some simulated datasets, which are left out of {?its summary/their summaries}.", "i" = "The most are {most$n} of {nsim}, for {.val {most$statistic}} in model {.val {most$model}}."))
  }
  invisible(NULL)
}

##' The observed statistics as a data frame
##' @noRd
observed_statistics <- function(values) {
  do.call(rbind, lapply(names(values), function(name) {
    data.frame(model = rep(name, nrow(values[[name]])),
               statistic = as.character(rownames(values[[name]])),
               value = unname(values[[name]][, 1L]))
  }))
}

##' The simulated statistics as a data frame
##' @noRd
simulated_statistics <- function(values) {
  do.call(rbind, lapply(names(values), function(name) {
    simulated <- values[[name]][, -1L, drop = FALSE]
    data.frame(model = rep(name, length(simulated)),
               simulation = rep(seq_len(ncol(simulated)), each = nrow(simulated)),
               statistic = rep(as.character(rownames(simulated)), ncol(simulated)),
               value = as.numeric(simulated))
  }))
}

##' Curves of the observed and simulated data as a table with their envelope
##'
##' @param curves A matrix with a curve in each column, the observed one first.
##' @noRd
check_envelope <- function(curves) {
  envelope <- global_rank_envelope(curves, check_level())
  observed <- curves[, 1L]
  data.frame(observed = observed,
             median = apply(curves[, -1L, drop = FALSE], 1L, stats::median),
             lower = envelope$lower, upper = envelope$upper,
             outside = observed < envelope$lower | observed > envelope$upper)
}

##' The variograms of the models with their envelopes
##'
##' @param residuals A named list of matrices with a location in each row and a
##'   dataset in each column, or `NULL` for a model without residuals.
##' @noRd
check_variograms <- function(residuals, classes) {
  if (is.null(classes)) return(NULL)
  shown <- classes$n_pairs > 0L
  do.call(rbind, lapply(names(residuals), function(name) {
    if (is.null(residuals[[name]])) return(NULL)
    curves <- binned_semivariances(residuals[[name]], classes)
    cbind(data.frame(model = name, distance = classes$distance[shown],
                     n_pairs = classes$n_pairs[shown]),
          check_envelope(curves[shown, , drop = FALSE]))
  }))
}

##' The distributions of the models with their envelopes
##'
##' @param values A named list of matrices with a location or an observation in
##'   each row and a dataset in each column.
##' @noRd
check_distributions <- function(values) {
  at <- distribution_values(values)
  do.call(rbind, lapply(names(values), function(name) {
    curves <- apply(values[[name]], 2L, function(x) {
      findInterval(at, sort(x)) / length(x)
    })
    cbind(data.frame(model = name, value = at),
          check_envelope(matrix(curves, nrow = length(at))))
  }))
}

##' The values at which distributions are compared
##'
##' Every value when there are few of them, and otherwise quantiles of all the
##' values together, which are closer to each other where there are more data.
##' They go from the 0.5% to the 99.5% quantile, or further if the observed
##' values do, because beyond them every distribution function is near zero or
##' one.
##' @noRd
distribution_values <- function(values, n_values = 200L) {
  pooled <- sort(unlist(values, use.names = FALSE))
  # The quantiles of a million values are enough.
  kept <- round(seq(1, length(pooled), length.out = min(length(pooled), 1e6)))
  quantiles <- function(probabilities) {
    stats::quantile(pooled[kept], probabilities, type = 1L, names = FALSE)
  }
  limits <- range(values[[1L]][, 1L], quantiles(c(0.005, 0.995)))
  distinct <- unique(pooled)
  if (length(distinct) > n_values) {
    distinct <- unique(c(quantiles((seq_len(n_values) - 0.5) / n_values), limits))
  }
  sort(distinct[distinct >= limits[1L] & distinct <= limits[2L]])
}

##' @title Summarise the assessment of models against simulated data
##'
##' @description Gives, for each model and statistic, the observed value, the
##'   median and the central 95% interval of the simulated values, and the
##'   proportion of simulated values at or above the observed one.
##'
##' @param object The result of [assess_fit()].
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##'
##' @return A data frame of class `summary.RiskMap_predictive_check` with one
##'   row for each statistic and model, and the columns `model`, `statistic`,
##'   `observed`, `median`, `lower`, `upper`, `at_or_above`, the proportion of
##'   simulated values at or above the observed one, and `n_not_finite`, the
##'   number of simulated datasets in which the statistic isn't finite, which
##'   are left out.
##'
##' @details
##' The observed value is unusual for the model when it is outside the interval.
##'
##' The proportion at or above the observed value is not a p-value. The
##' parameters are estimated from the same data, so that it is closer to 0.5
##' than a p-value would be, and a statistic that the model was fitted to can
##' hardly have an extreme one. It counts the simulated values equal to the
##' observed one, and where there are many of them it can be 1, or near 1, for a
##' value that isn't unusual. See [assess_fit()].
##'
##' The printed summary also says in how many distance classes the observed
##' variogram is outside its joint envelope, where the observed distribution is,
##' and, for binomial models, whether the prevalences are those of observations or
##' of locations.
##'
##' @seealso [assess_fit()], and [to_table()] for a table to render.
##' @method summary RiskMap_predictive_check
##' @export
summary.RiskMap_predictive_check <- function(object, ...) {
  check_dots_unused(...)
  models <- names(object$family)
  # Without statistics that are numbers, there are no rows.
  statistics <- unique(object$observed$statistic)
  rows <- expand.grid(model = models, statistic = statistics,
                      stringsAsFactors = FALSE)
  # A statistic may be calculated for some of the models only.
  key <- function(data) paste(data$model, data$statistic, sep = "\r")
  rows <- rows[key(rows) %in% key(object$observed), ]
  simulated <- split(object$simulated$value, key(object$simulated))
  values <- lapply(seq_len(nrow(rows)), function(i) {
    observed <- object$observed$value[match(key(rows[i, ]),
                                            key(object$observed))]
    summarise_statistic(observed, simulated[[key(rows[i, ])]])
  })
  none <- summarise_statistic(NA_real_, numeric())[0L, ]
  out <- cbind(rows, do.call(rbind, c(list(none), values)))
  row.names(out) <- NULL
  attr(out, "check") <- object[c("nsim", "seed", "level", "quantity", "response",
                                 "statistic_labels", "variogram",
                                 "distribution", "types", "breaks", "unit")]
  class(out) <- c("summary.RiskMap_predictive_check", "data.frame")
  out
}

##' One row of the summary of a predictive check
##' @noRd
summarise_statistic <- function(observed, simulated) {
  finite <- simulated[is.finite(simulated)]
  interval <- if (length(finite)) {
    stats::quantile(finite, c(0.5, 0.025, 0.975), names = FALSE)
  } else {
    rep(NA_real_, 3L)
  }
  above <- if (is.finite(observed) && length(finite)) {
    mean(at_or_above(finite, observed))
  } else {
    NA_real_
  }
  data.frame(observed = observed, median = interval[1L], lower = interval[2L],
             upper = interval[3L], at_or_above = above,
             n_not_finite = length(simulated) - length(finite))
}

##' Whether simulated values are at or above the observed one
##'
##' Values equal but for rounding are equal. None is above an observed value that
##' isn't finite.
##' @noRd
at_or_above <- function(simulated, observed) {
  if (!is.finite(observed)) return(rep(FALSE, length(simulated)))
  tolerance <- 1e-10 * max(abs(c(observed, simulated)))
  simulated >= observed - tolerance
}

##' @title Print the assessment of models against simulated data
##'
##' @description Prints the models, the number of simulated datasets and, for
##'   each model and statistic, where the observed value lies among the simulated
##'   ones, with what the joint envelopes show.
##'
##' @param x The result of [assess_fit()], or its summary.
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored. Some columns of the summary print as a data frame, and
##'   the arguments go to its print method.
##' @return Invisibly returns `x`.
##' @seealso [assess_fit()]
##' @method print RiskMap_predictive_check
##' @export
print.RiskMap_predictive_check <- function(x, ...) {
  check_dots_unused(...)
  out_print(format_predictive_check(x))
  invisible(x)
}

##' @rdname print.RiskMap_predictive_check
##' @method print summary.RiskMap_predictive_check
##' @export
print.summary.RiskMap_predictive_check <- function(x, ...) {
  # Some of the columns print as a data frame.
  if (!is_whole_check_summary(x)) return(NextMethod())
  check_dots_unused(...)
  out_print(c(out_rule("Simulated and observed data"),
              format_check_summary(x)))
  invisible(x)
}

# As for knit_print.RiskMap_extrapolation(), so that a document that prints
# data frames as tables prints the summary as the console does.
##' @exportS3Method knitr::knit_print
knit_print.summary.RiskMap_predictive_check <- function(x, ...) {
  if (!is_whole_check_summary(x)) return(NextMethod())
  knitr::normal_print(x)
}

##' Whether a summary of a predictive check has what its print uses
##' @noRd
is_whole_check_summary <- function(x) {
  all(c("model", "statistic", "observed", "median", "lower", "upper",
        "at_or_above", "n_not_finite") %in% names(x)) &&
    is.list(attr(x, "check", exact = TRUE))
}

##' The lines that print a predictive check
##' @noRd
format_predictive_check <- function(x) {
  n_models <- length(x$family)
  setting <- paste0(
    out_count(n_models, "model"), " assessed on ", format_count(x$nsim),
    " simulated datasets", if (n_models > 1L) " each",
    if (!is.null(x$seed)) paste0(" (seed ", x$seed, ")"), ", with ",
    out_count(x$n_observations, "observation"), " at ",
    out_count(x$n_locations, "location"), "."
  )
  drawn <- c(statistic = "the statistics", distribution = "the distribution",
             variogram = "the variogram")[drawn_types(x)]
  hint <- if (length(drawn)) {
    c(i = paste0(out_fun("plot_predictive_check"), " draws ",
                 join_words(drawn), "."))
  }
  c(out_rule("Simulated and observed data"),
    out_text(setting),
    "",
    format_check_summary(summary(x)),
    if (length(hint)) c("", out_bullets(hint)))
}

##' The lines of the table of a predictive check and of its envelopes
##' @noRd
format_check_summary <- function(x) {
  settings <- attr(x, "check")
  about_statistics <- if (nrow(x)) {
    paste0("At or above is the percentage of simulated values that are equal to ",
           "the observed one or above it. An observed value outside the 95% ",
           "interval is unusual for the model, and 50% is not a target.")
  }
  notes <- c(about_statistics,
             format_check_unit(settings),
             format_check_variogram(settings),
             format_check_distribution(settings$distribution, settings))
  c(if (nrow(x)) c(format_check_tables(x), ""),
    if (length(notes)) out_bullets(c(out_as_bullets(notes, "i"), format_check_missing(x))))
}

##' One table for each statistic, when there are several models, or one for all
##' the statistics of a model
##' @noRd
format_check_tables <- function(x) {
  table <- check_summary_table(x, digits = 3L)
  lines <- function(rows) {
    out_capture(print(rows, row.names = FALSE, right = TRUE))
  }
  if (is.null(table$Model)) return(lines(table))
  statistics <- unique(table$Statistic)
  unlist(lapply(seq_along(statistics), function(k) {
    c(if (k > 1L) "", out_heading(statistics[k]),
      lines(table[table$Statistic == statistics[k], names(table) != "Statistic"]))
  }))
}

##' The simulated values that aren't finite, as a bullet
##' @noRd
format_check_missing <- function(x) {
  missing <- x[x$n_not_finite > 0, ]
  if (!nrow(missing)) return(NULL)
  statistics <- check_statistic_labels(missing$statistic, attr(x, "check"))
  c("!" = paste0("Simulated values that aren't finite are left out, ",
                 paste0(missing$n_not_finite, " of ", statistics, " (",
                        missing$model, ")", collapse = ", "), "."))
}

##' Whether the prevalences are those of observations or of locations
##'
##' Functions of a dataset see the observations in any case.
##' @noRd
format_check_unit <- function(settings) {
  compared <- length(settings$statistic_labels) || !is.null(settings$distribution)
  if (settings$quantity != "prevalence" || !compared) return(NULL)
  paste0("The prevalences are those of the ", settings$unit, "s",
         if (settings$unit == "location") {
           ", because every denominator is 1"
         }, ".")
}

##' Where the variograms are outside their envelopes
##' @noRd
format_check_variogram <- function(settings) {
  variogram <- settings$variogram
  if (!"variogram" %in% settings$types) return(NULL)
  if (is.null(variogram)) {
    return(paste0("There is no variogram, because ",
                  missing_variogram(settings$breaks), "."))
  }
  models <- unique(variogram$model)
  outside <- vapply(models, function(name) {
    sum(variogram$outside[variogram$model == name])
  }, integer(1))
  n_classes <- sum(variogram$model == models[1L])
  counts <- if (length(models) == 1L) {
    as.character(outside)
  } else {
    join_words(paste0(outside, " (", models, ")"))
  }
  paste0("The observed variogram of the residuals is outside the ",
         100 * settings$level, "% joint envelope in ", counts, " of ", n_classes,
         " distance classes.")
}

##' Whether the distributions are inside their envelopes
##' @noRd
format_check_distribution <- function(distribution, settings) {
  if (is.null(distribution)) return(NULL)
  models <- unique(distribution$model)
  outside <- vapply(models, function(name) {
    any(distribution$outside[distribution$model == name])
  }, logical(1))
  several <- length(models) > 1L
  envelope <- paste0("the ", 100 * settings$level, "% joint envelope")
  values <- switch(settings$quantity, prevalence = "prevalences",
                   rate = "rates", "values")
  of_models <- function(which) {
    if (several) paste0(" for ", join_words(models[which]))
  }
  if (!any(outside)) {
    return(paste0("The observed distribution is inside ", envelope, "."))
  }
  paste0("The observed distribution is outside ", envelope, " at some ", values,
         of_models(outside),
         if (!all(outside)) paste0(", and inside it", of_models(!outside)), ".")
}

##' The labels of the statistics, those given by name spelt out
##' @noRd
check_statistic_labels <- function(statistics, settings) {
  labels <- settings$statistic_labels[statistics]
  unname(ifelse(is.na(labels), statistics, labels))
}

##' The summary of a predictive check as a table of text
##' @noRd
check_summary_table <- function(x, digits) {
  settings <- attr(x, "check")
  values <- as.matrix(as.data.frame(unclass(x)[c("observed", "median", "lower",
                                                "upper")]))
  text <- matrix("", nrow(values), ncol(values))
  for (i in seq_len(nrow(values))) {
    text[i, ] <- format_compared_values(values[i, ], digits)
  }
  text[!is.finite(values)] <- "NA"
  share <- ifelse(is.finite(x$at_or_above),
                  paste0(formatC(100 * x$at_or_above, format = "f", digits = 1L),
                         "%"), "NA")
  out <- data.frame(
    Statistic = check_statistic_labels(x$statistic, settings),
    Model = x$model, Observed = text[, 1L], Median = text[, 2L],
    Interval = paste(text[, 3L], "to", text[, 4L]), Above = share
  )
  names(out)[5:6] <- c("95% interval", "At or above")
  if (length(unique(x$model)) == 1L) out$Model <- NULL
  out
}

##' @rdname to_table
##' @export
to_table.RiskMap_predictive_check <- function(object, digits = 3,
                                              col_names = NULL, ...) {
  to_table(summary(object), digits = digits, col_names = col_names, ...)
}

##' @rdname to_table
##' @export
to_table.summary.RiskMap_predictive_check <- function(object, digits = 3,
                                                      col_names = NULL, ...) {
  check_positive_integer(digits, "digits", allow_zero = TRUE)
  if (!is_whole_check_summary(object)) {
    abort_input("{.arg object} must be the whole summary of a predictive check, not a part of it.", arg = "object")
  }
  if (!nrow(object)) {
    abort_input(c("{.arg object} has no statistics to tabulate.", "i" = "Include some in {.arg statistics} of {.fn assess_fit}."), arg = "object")
  }
  render_results_table(check_summary_table(object, digits), row_names = FALSE,
                       col_names = col_names, ...)
}
