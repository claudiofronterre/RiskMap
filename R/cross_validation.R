# Cross-validation of geostatistical models, assess_prediction(), and the
# methods of its result. The test sets are in test_sets.R, the calibration in
# calibration.R and the plots in cross_validation_plots.R.

##' @title Assess predictions by cross-validation
##'
##' @description
##' Compares the predictions of one or more models with observations left out
##' of the data. Each test set is predicted from its training set, and every
##' test observation is scored against its predictive distribution.
##' - **Scoring rules**, the continuous ranked probability score (CRPS) and the
##'   scaled CRPS (SCRPS) of Bolin and Wallin (2023), measure the sharpness and
##'   the calibration of the predictive distributions together. Both are
##'   losses, so lower values are better.
##' - **Calibration** is assessed with the probability integral transform (PIT)
##'   of Gaussian outcomes and its nonrandomised version for counts, through
##'   the area between the calibration curve and the identity line
##'   (Giorgi *et al.* 2026).
##'
##' @param models A fit from [glgpm()], or a list of fits to the same data. The
##'   names of the list name the models in the results.
##' @param test_sets The test sets, from [create_test_sets()], or a list with
##'   one element for each test set, the row numbers of its test observations
##'   or a list of the training rows, `train`, and the test rows, `test`.
##' @param refit Logical, whether the parameters are estimated again in each
##'   training set. If `FALSE` (the default), they are fixed at the estimates
##'   from all the data, and only the predictions are conditioned on the
##'   training set. See Details.
##' @param control_mcmc Control settings of the predictive draws, from
##'   [set_control_mcmc()].
##' @param seed Optional non-negative integer that makes the assessment
##'   reproducible. It sets a seed for each test set, which is passed to
##'   [glgpm()] when the models are refitted and to
##'   [predict()][predict.RiskMap()]. All models share the seed of a test set,
##'   so that their predictions use common random numbers. The caller's random
##'   number state is restored on exit. When `NULL` (the default), the seeds
##'   are drawn from the session's random number generator.
##' @param messages Logical, whether to report on the assessment, with a
##'   progress bar as described in [riskmap_progress]. Defaults to the option
##'   `RiskMap.messages`, or `TRUE` if it is not set.
##' @param ... Not used, and an error if given, so that a misspelt or former
##'   argument such as `method` is not ignored.
##'
##' @return A list of class `RiskMap_cross_validation`, containing
##' \describe{
##'   \item{scores}{A data frame with one row for each model and test
##'     observation. `model`, `test_set` and `row` identify it, with `row` the
##'     row of the data, `CRPS` and `SCRPS` are its scores, and `pit_lower`
##'     and `pit_upper` are its predictive distribution function just below
##'     the observation and at the observation, which are equal for Gaussian
##'     models.}
##'   \item{calibration}{A list with `area`, a data frame with one row for each
##'     model, its `PIT_area`, `PIT_area_expected` and `PIT_area_limit`,
##'     `curves`, a data frame with the calibration curve of each model and the
##'     95% joint band of a correct model (`lower` and `upper`) at probability
##'     levels `u`, and `simulated`, the simulated curves the band comes from, a
##'     matrix for each model with a row for each simulation and a column for
##'     each level `u`.}
##'   \item{test_sets}{The test sets, of class `RiskMap_test_sets`.}
##'   \item{refit}{The value of `refit`.}
##'   \item{family}{The family of each model.}
##'   \item{failures}{A data frame with one row for each model and test set
##'     whose refit, prediction or scores failed, with columns `model`,
##'     `test_set` and `message`. The scores of these test sets are `NA`.}
##' }
##'
##' @details
##' # Fixed and refitted parameters
##'
##' Fixed parameters are much faster, and with many locations and few
##' parameters the scores change little. The test observations still
##' contribute to the fixed estimates, however, so the scores are optimistic,
##' and more so for models with more parameters. Use `refit = TRUE` to compare
##' models with different numbers of parameters. It costs one fit for each
##' test set and model.
##'
##' # Predictive distributions
##'
##' The observations of a test set are predicted jointly, with the nugget and
##' the random effects when the model has them. A random effect has its
##' distribution given the training data, and a level that the training data
##' don't have is drawn from the distribution of the random effect. The scores
##' of an observation are those of its own predictive distribution, and the
##' joint draws are for the reference of its calibration (see Calibration). Test
##' observations with missing values in a variable used by any of the models
##' are left out.
##'
##' The memory of a test set grows with the square of its number of locations
##' and the time of its joint draws with the cube, so large datasets are better
##' divided into more test sets.
##'
##' For a Gaussian model, the predictive distribution is Gaussian, with the
##' mean and the standard deviation of the draws, which include the
##' measurement error. The predictive distribution of a binomial, Poisson or
##' negative binomial count \eqn{Y} is the distribution of the family given
##' the mean, averaged over the \eqn{J} draws of the mean from its own
##' predictive distribution,
##' \deqn{F(k) = \frac{1}{J} \sum_{j=1}^{J} P(Y \le k \mid \mu_j), \quad
##'   k = 0, 1, 2, \ldots}
##' Here \eqn{\mu_j} is the \eqn{j}th draw of the prevalence or rate, and
##' \eqn{P(Y \le k \mid \mu_j)} also depends on the denominator and, for a
##' negative binomial model, on its size.
##'
##' # Scores
##'
##' Every score is a loss, so lower values are better. The CRPS is in the
##' units of the response, which for binomial, Poisson and negative binomial
##' models are counts (for example numbers of positive individuals), not
##' prevalences or rates. The SCRPS is reported as
##' \deqn{\frac{1}{2}\left[1 + \frac{\mathrm{CRPS}(F, y)}{c} + \log(2c)\right],
##'   \quad c = \frac{1}{2} E|X - X'|,}
##' with \eqn{X} and \eqn{X'} independent draws from the predictive
##' distribution \eqn{F}. This is the negative of the scaled CRPS of Bolin and
##' Wallin (2023), whose scores are higher for better predictions, so that the
##' CRPS and the SCRPS share one orientation. The SCRPS can be negative.
##'
##' For counts,
##' \deqn{\mathrm{CRPS}(F, y) = \sum_{k \ge 0} [F(k) - 1(k \ge y)]^2, \quad
##'   c = \sum_{k \ge 0} F(k) [1 - F(k)].}
##' The scores of counts have no random step beyond the draws of the mean, as
##' no counts are simulated for them. A model therefore gets the same scores
##' whatever its position in `models`. Where counts above 100 are predicted,
##' \eqn{F} is computed on a grid of counts and is linear in between. The grid
##' is fine enough for a relative error of about \eqn{10^{-4}} or less in the
##' CRPS and in \eqn{c}.
##'
##' # Calibration
##'
##' The calibration curve of \eqn{n} test observations is
##' \deqn{\frac{1}{n} \sum_{i=1}^{n} G_i(u), \quad 0 \le u \le 1.}
##' For a Gaussian outcome, \eqn{G_i(u)} is 1 if \eqn{F_i(y_i) \le u} and 0
##' otherwise, so the curve is the empirical distribution function of the PIT
##' values. For a count it is the average nonrandomised PIT (AnPIT), with
##' \eqn{G_i(u)} equal to 0 up to \eqn{F_i(y_i - 1)}, to 1 from
##' \eqn{F_i(y_i)}, and linear in between (Czado *et al.* 2009). `PIT_area` is
##' the area between the curve and the identity line, from 0 for a perfectly
##' calibrated curve to at most 0.5.
##'
##' # How to read the area
##'
##' The area of a correct model is not zero, because its observations are finite
##' in number, and because neighbouring test observations are above or below
##' their predictions together. Its distribution is found by simulation.
##' 1. Test observations are drawn from the model, jointly for each test set from
##'    its joint predictive distribution, once for each of 999 draws of
##'    `control_mcmc`, evenly spaced, or for every draw if there are fewer. An
##'    observation in several test sets has the same PIT value in all of them,
##'    the one drawn in its first test set.
##' 2. Each draw gives a calibration curve and an area, computed as for the data.
##' 3. `PIT_area_expected` is the mean of these areas, and `PIT_area_limit` the
##'    value that 95% of them are below, which needs at least 19 draws. The band of [plot_cross_validation()] is
##'    the joint envelope of the curves, computed as in [variogram()], which the
##'    whole curve of a correct model is inside with a probability of 95%.
##'
##' An area below `PIT_area_limit` is one that a correct model would give, and an
##' area above it points to miscalibration. Below the limit, a smaller area, or
##' one closer to `PIT_area_expected`, doesn't mean a better model. It mostly
##' reflects chance, the number of observations and the dependence between them,
##' so the areas of models shouldn't be ranked. The CRPS and the SCRPS are for
##' comparing models. The area is also not sensitive to every departure. Too many
##' observations in both tails change it little, and the curve and the histogram
##' of [plot_cross_validation()] show them.
##'
##' The simulation draws from the model with its parameters at the values used for
##' the prediction, and it treats the test sets as independent, although their
##' prediction errors are not. A test observation is predicted from neighbours
##' that are predicted in other test sets, so their errors partly cancel. The
##' areas of a correct model therefore tend to be smaller than the simulated
##' ones, more so with leave-one-out and random test sets and with repeated
##' divisions, and most with leave-one-out for binary outcomes. The limit and
##' the band are cautious there, so a curve inside the band is not evidence of
##' good calibration.
##'
##' # Failures
##'
##' A failed refit or prediction does not stop the assessment. The scores of
##' that test set are missing and a warning reports the failures. So are those
##' of a test set of counts that can't be scored, because a draw of the mean is
##' not finite, is outside the range of the family or gives counts above
##' \eqn{2^{53}}, or because a binomial denominator is not a whole number.
##'
##' @seealso [create_test_sets()] makes the test sets. The methods for a
##'   cross-validation are [print()][print.RiskMap_cross_validation()],
##'   [summary()][summary.RiskMap_cross_validation()], [to_table()] and
##'   [plot_cross_validation()], which `plot()` of the result also calls.
##'
##' @references
##' Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of
##' proper scoring rules. *Statistical Science*, 38(1), 140–159. \doi{10.1214/22-STS864}.
##'
##' Czado, C., Gneiting, T., & Held, L. (2009). Predictive model assessment
##' for count data. *Biometrics*, 65(4), 1254–1261.
##' \doi{10.1111/j.1541-0420.2009.01191.x}.
##'
##' Giorgi, E., Fronterre, C. & Diggle, P. J. (2026). A decay-adjusted spatio-temporal
##' model to account for the impact of mass drug administration on neglected
##' tropical disease prevalence. *Journal of the Royal Statistical Society Series
##' A: Statistics in Society*.\doi{10.1093/jrsssa/qnag100}.
##'
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:100, ]
##' fit <- glgpm(y ~ gp(), data = example_data, family = "gaussian",
##'              messages = FALSE)
##' covariate_fit <- glgpm(y ~ log(pop_dens) + gp(), data = example_data,
##'                        family = "gaussian", messages = FALSE)
##'
##' test_sets <- create_test_sets(example_data, method = "cluster",
##'                               n_sets = 5, seed = 1)
##' result <- assess_prediction(
##'   list(intercept = fit, covariate = covariate_fit),
##'   test_sets,
##'   seed = 1
##' )
##' summary(result, reference = "intercept")
##' plot_cross_validation(result)
##'
##' # Test sets given as rows of the data
##' result <- assess_prediction(fit, list(1:50, 51:100), seed = 1)
##' summary(result)
##'
##' @export
assess_prediction <- function(models,
                              test_sets,
                              refit = FALSE,
                              control_mcmc = set_control_mcmc(),
                              seed = NULL,
                              messages = getOption("RiskMap.messages", TRUE),
                              ...) {
  check_supplied(models)
  # Before `test_sets`, which a call with the former arguments doesn't have.
  check_dots_unused(..., .hints = former_assessment_arguments())
  check_supplied(test_sets)
  models <- read_assessed_models(models)
  check_logical(refit)
  check_control(control_mcmc, "set_control_mcmc")
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)
  check_logical(messages)
  test_sets <- read_test_sets(test_sets, models[[1L]])
  n_sets <- length(test_sets$sets)

  warn_shared_locations(models, test_sets)
  if (messages) inform_assessment(models, n_sets, refit)
  complete_rows <- complete_model_rows(models)
  if (!all(complete_rows)) {
    riskmap_warn("{sum(!complete_rows)} observation{?s} with missing values in the model variables {?is/are} not used as test observations.")
  }

  if (!is.null(seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(seed)
  }
  # One seed per test set, shared by the models so that their predictions use
  # common random numbers. The measurement error of Gaussian models, and the
  # observations simulated for the calibration, are drawn from separate
  # streams, so that they don't change the spatial draws.
  seeds <- data.frame(
    prediction = sample.int(.Machine$integer.max, n_sets),
    noise = sample.int(.Machine$integer.max, n_sets),
    reference = sample.int(.Machine$integer.max, n_sets)
  )

  bar <- progress_bar(if (refit) "Refitting and scoring test sets" else
                        "Scoring test sets",
                      length(models) * n_sets, messages)
  assessed <- Map(function(fit, name) {
    assess_one_model(fit, name, test_sets, complete_rows, refit, control_mcmc,
                     seeds, bar)
  }, models, names(models))
  progress_done(bar)

  failures <- do.call(rbind, lapply(assessed, `[[`, "failures"))
  warn_failed_assessments(failures, length(models) * n_sets)
  out <- list(
    scores = do.call(rbind, lapply(assessed, `[[`, "scores")),
    calibration = combine_calibrations(lapply(assessed, `[[`, "calibration")),
    test_sets = test_sets, refit = refit,
    family = vapply(models, `[[`, character(1), "family"),
    failures = failures
  )
  row.names(out$scores) <- NULL
  row.names(out$failures) <- NULL
  class(out) <- "RiskMap_cross_validation"
  stamp_version(out)
}

##' What to use instead of the former arguments of `assess_prediction()`
##' @noRd
former_assessment_arguments <- function() {
  made <- "Make the test sets with {.fn create_test_sets} and give them in {.arg test_sets}."
  c(method = made, fold = made, min_dist = made, size = made,
    iter = made, plot_fold = "Use {.fn plot_test_sets} to see the test sets.",
    user_split = "Give the list of test sets in {.arg test_sets}, with {.code train} and {.code test} for {.code in_id} and {.code out_id}.",
    metrics = "All metrics are always computed.",
    which_metric = "All metrics are always computed.",
    keep_par_fixed = "Use {.arg refit}, which is {.code TRUE} to re-estimate the parameters in each training set.",
    object = "It is now {.arg models}.",
    n_size = made,
    control_sim = "It is now {.arg control_mcmc}.")
}

##' The models to assess, as a named list of fits to the same data
##' @noRd
read_assessed_models <- function(models, call = rlang::caller_env()) {
  if (inherits(models, "RiskMap")) models <- list(models)
  if (!is.list(models) || !length(models) ||
      !all(vapply(models, inherits, logical(1), what = "RiskMap"))) {
    abort_input("{.arg models} must be a fit from {.fn glgpm} or a list of fits, not {.obj_type_friendly {models}}.", arg = "models", call = call)
  }
  unnamed <- if (is.null(names(models))) {
    rep(TRUE, length(models))
  } else {
    is.na(names(models)) | names(models) == ""
  }
  names(models)[unnamed] <- paste("Model", which(unnamed))
  if (anyDuplicated(names(models))) {
    repeated <- unique(names(models)[duplicated(names(models))])
    abort_input("The names of {.arg models} must differ, and {.val {repeated}} {?is/are} repeated.", arg = "models", call = call)
  }
  check_same_data(models, call)
  models
}

##' Check that the models are fitted to the same observations
##'
##' The same rows, at the same locations, with the same responses and the same
##' denominators. Scores of different outcomes can't be compared.
##' @noRd
check_same_data <- function(models, call) {
  first <- models[[1L]]
  for (name in names(models)[-1L]) {
    fit <- models[[name]]
    if (nrow(fit$data) != nrow(first$data)) {
      abort_input(c("The models in {.arg models} must be fitted to the same observations.", "x" = "Model {.val {name}} is fitted to {nrow(fit$data)} of them, and {.val {names(models)[1L]}} to {nrow(first$data)}."), arg = "models", call = call)
    }
    if (!identical(st_geometry(fit$data), st_geometry(first$data))) {
      abort_input("The models in {.arg models} must be fitted to data with the same rows in the same order, in the same coordinate reference system.", arg = "models", call = call)
    }
    if (!same_values(fit$y, first$y)) {
      abort_input(c("The models in {.arg models} must have the same response.", "x" = "The responses of {.val {name}} and {.val {names(models)[1L]}} differ, so their scores can't be compared."), arg = "models", call = call)
    }
    if (!same_values(assessed_denominators(fit), assessed_denominators(first))) {
      abort_input(c("The models in {.arg models} must have the same denominators.", "x" = "The denominators of {.val {name}} and {.val {names(models)[1L]}} differ, so their scores can't be compared."), arg = "models", call = call)
    }
  }
  invisible(TRUE)
}

##' The denominators of the responses of a fit, one for each response if it
##' has none
##' @noRd
assessed_denominators <- function(fit) {
  units_m <- as.numeric(fit$units_m)
  if (!length(units_m)) rep(1, length(fit$y)) else units_m
}

##' Whether two numeric vectors are the same but for rounding
##' @noRd
same_values <- function(x, y) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), check.attributes = FALSE,
                   tolerance = 1e-10))
}

##' Warn about test locations that also have training data
##'
##' A test observation at a location that also has training observations is
##' predicted with a new nugget draw, ignoring the nugget seen in training.
##' @noRd
warn_shared_locations <- function(models, test_sets) {
  with_nugget <- vapply(models, function(fit) fitted_nugget_variance(fit) > 0,
                        logical(1))
  if (!any(with_nugget)) return(invisible(NULL))
  location <- test_sets$location
  shared <- which(vapply(test_sets$sets, function(set) {
    any(location[set$test] %in% location[set$train])
  }, logical(1)))
  if (!length(shared)) return(invisible(NULL))
  riskmap_warn(c("In {cli::qty(length(shared))}test set{?s} {shared}, some locations also have training observations.", "i" = "Their predictions ignore the nugget seen in training, so the scores of models with a nugget aren't reliable there."))
}

##' Say how the parameters are set
##' @noRd
inform_assessment <- function(models, n_sets, refit) {
  if (refit) {
    n_fits <- n_sets * length(models)
    riskmap_inform("Parameters are re-estimated in each training set ({n_fits} model fits).")
    return(invisible(NULL))
  }
  riskmap_inform(c("Parameters are fixed at the estimates from all the data.", "i" = "Set {.code refit = TRUE} to re-estimate them in each training set."))
  n_parameters <- vapply(models, function(fit) length(unlist(fit$estimate)),
                         integer(1))
  if (length(unique(n_parameters)) > 1L) {
    riskmap_inform(c("!" = "The models have different numbers of parameters. With fixed parameters the comparison tends to favour the model with more of them, so {.code refit = TRUE} is recommended."))
  }
}

##' Whether each row of the data can be used by every model
##' @noRd
complete_model_rows <- function(models) {
  Reduce(`&`, lapply(models, function(fit) {
    variables <- intersect(c(all.vars(fit$formula),
                             as.character(fit$call$denominator)),
                           names(fit$data))
    complete.cases(st_drop_geometry(fit$data)[, variables, drop = FALSE])
  }))
}

##' Warn about the test sets whose refit, prediction or scores failed
##' @noRd
warn_failed_assessments <- function(failures, n_assessed) {
  if (nrow(failures) == 0L) return(invisible(NULL))
  first <- failures[1L, ]
  riskmap_warn(c("{nrow(failures)} of {n_assessed} test-set assessments failed, and their scores are missing.", "i" = "The first failure, for model {.val {first$model}} and test set {first$test_set}, was {.val {first$message}}.", "i" = "See the {.code failures} element of the result."))
}

##' Score one model on every test set
##'
##' A failed fit or prediction leaves the scores of its test set missing
##' instead of stopping the assessment.
##' @return A list with `scores`, `failures` and `calibration`.
##' @noRd
assess_one_model <- function(fit, name, test_sets, complete_rows, refit,
                             control_mcmc, seeds, bar) {
  failures <- data.frame(model = character(), test_set = integer(),
                         message = character())
  n_curves <- n_simulated_curves()
  assessed <- lapply(seq_along(test_sets$sets), function(i) {
    rows <- test_sets$sets[[i]]$test
    rows <- rows[complete_rows[rows]]
    result <- if (length(rows)) tryCatch(
      score_test_set(fit, test_sets$sets[[i]]$train, rows, refit, control_mcmc,
                     seeds[i, ], n_curves),
      error = function(error) {
        failures[nrow(failures) + 1L, ] <<- list(name, i, condition_text(error))
        NULL
      }
    )
    progress_tick(bar)
    missing <- rep(NA_real_, length(rows))
    scores <- result$scores %||%
      data.frame(CRPS = missing, SCRPS = missing, pit_lower = missing,
                 pit_upper = missing)
    list(scores = cbind(data.frame(model = rep(name, length(rows)),
                                   test_set = rep(i, length(rows)),
                                   row = rows), scores),
         predictive = result$predictive)
  })
  scores <- do.call(rbind, lapply(assessed, `[[`, "scores"))
  scored <- !is.na(scores$pit_lower)
  reference <- if (any(scored)) {
    calibration_reference(fit, test_sets, lapply(assessed, `[[`, "predictive"),
                          n_curves, seeds$reference)
  }
  list(
    scores = scores, failures = failures,
    calibration = assess_calibration(scores$pit_lower[scored],
                                     scores$pit_upper[scored], reference)
  )
}

##' Predict and score one test set
##'
##' @param train,rows The training and test rows of the data.
##' @param seeds The seeds of the test set.
##' @param n_curves The number of draws kept for the calibration.
##' @return A list with the `scores` of the test observations and their
##'   `predictive` distributions, for `calibration_reference()`. These are a list
##'   with the `rows`, `draws`, joint draws of the outcomes of a Gaussian test set
##'   or of the means of a test set of counts with an observation in each row,
##'   and the distribution of each observation, its `mean` and `sd`, or its
##'   `distributions`, with `units_m` and the `parameters` of the family.
##' @noRd
score_test_set <- function(fit, train, rows, refit, control_mcmc, seeds,
                           n_curves) {
  conditioned <- if (refit) {
    refit_on_rows(fit, train, control_mcmc, seeds$prediction)
  } else {
    condition_on_rows(fit, train)
  }
  traits <- family_traits(fit$family)
  parameters <- fit_family_values(conditioned)
  link <- predict_test_set(conditioned, fit$data[rows, ], control_mcmc,
                           seeds$prediction)
  if (traits$measurement_error) {
    link <- link + measurement_noise(conditioned, length(link), seeds$noise)
  }
  # A Gaussian fit stores no link, so its identity comes from the family.
  means <- (fit$link_function$inv %||% traits$inverse_link)(link)
  dim(means) <- dim(link)
  y <- fit$y[rows]
  scored <- if (traits$discrete) {
    score_counts(y, assessed_denominators(fit)[rows], means, parameters, traits)
  } else {
    score_gaussian(y, means)
  }
  pick <- unique(round(seq(1, ncol(means), length.out = n_curves)))
  scored$predictive <- c(
    list(rows = rows, draws = means[, pick, drop = FALSE]), scored$predictive
  )
  scored
}

##' Refit a model to some of the rows of its data
##' @noRd
refit_on_rows <- function(fit, rows, control_mcmc, seed) {
  arguments <- list(
    formula        = fit_formula(fit),
    data           = fit$data[rows, ],
    family         = fit_family_argument(fit),
    distance_crs   = st_crs(fit$data),
    distance_units = fit$distance_units,
    control_mcmc   = control_mcmc,
    control_mcml   = attr(fit, "control_mcml") %||% set_control_mcml(),
    seed           = seed,
    fixed_pars     = fit_fixed_pars(fit),
    return_samples = FALSE,
    messages       = FALSE
  )
  # `denominator` is the unquoted name of a column, given only if the model
  # was fitted with one.
  denominator <- as.character(fit$call$denominator)
  if (length(denominator) == 1L && nzchar(denominator)) {
    arguments$denominator <- as.name(denominator)
  }
  arguments$invlink <- refit_invlink(fit)
  # The CRS of the fit was checked when it was fitted.
  without_crs_check(do.call(glgpm, arguments))
}

##' A fit with its estimates and only some of the rows of its data
##' @noRd
condition_on_rows <- function(fit, rows) {
  locations <- unique(fit$ID_coords[rows])
  fit$data <- fit$data[rows, ]
  fit$units_m <- fit$units_m[rows]
  fit$coords <- fit$coords[locations, , drop = FALSE]
  fit$y <- fit$y[rows]
  fit$D <- fit$D[rows, , drop = FALSE]
  if (!is.null(fit$cov_offset)) fit$cov_offset <- fit$cov_offset[rows]
  if (!is.null(fit$ID_re)) {
    terms <- names(fit$ID_re)
    random_effects <- prepare_random_effects(fit$data, terms)
    fit$ID_re <- as.data.frame(random_effects$ID_re)
    colnames(fit$ID_re) <- random_effects$names_re
    fit$re <- random_effects$re_unique_f
  }
  fit$ID_coords <- index_coordinate_rows(st_coordinates(fit$data))$index
  fit
}

##' Joint draws of the linear predictor of a test set
##'
##' With the nugget and the random effects of the model. A test observation
##' of a group that the training data don't have gets a draw from the
##' distribution of the random effect, shared by its group.
##' @return A matrix with a test observation in each row and a draw in each
##'   column.
##' @noRd
predict_test_set <- function(fit, newdata, control_mcmc, seed) {
  include <- c(
    character(),
    if (fitted_nugget_variance(fit) > 0) "nugget",
    if (length(fit$re)) "random_effects"
  )
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  set.seed(seed)
  # As predict() does, but for groups that aren't in the fit.
  input <- validate_prediction_data(fit, newdata, include,
                                    allow_new_levels = TRUE, messages = FALSE)
  prediction_sample_components(
    fit, input$model_data, "joint", include, control_mcmc, messages = FALSE,
    mcmc_seed = seed
  )$link
}

##' Draws of the measurement error of a Gaussian model
##' @noRd
measurement_noise <- function(fit, n_draws, seed) {
  variance <- coef(fit)$sigma2_me %||% fit_fix_var_me(fit)
  restore_seed <- preserve_random_seed()
  on.exit(restore_seed(), add = TRUE)
  set.seed(seed)
  sqrt(variance) * rnorm(n_draws)
}

##' Score a Gaussian test set
##'
##' @param outcomes Joint draws of the test observations, one in each row.
##' @noRd
score_gaussian <- function(y, outcomes) {
  mean <- rowMeans(outcomes)
  sd <- sqrt(rowSums((outcomes - mean)^2) / (ncol(outcomes) - 1))
  spread <- sd > 0
  z <- ifelse(spread, (y - mean) / sd, 0)
  # The CRPS of a point mass is the absolute error
  crps <- ifelse(spread, sd * (2 * dnorm(z) + z * (2 * pnorm(z) - 1) -
                                 1 / sqrt(pi)),
                 abs(y - mean))
  pit <- ifelse(spread, pnorm(z), as.numeric(y >= mean))
  list(scores = data.frame(CRPS = crps, SCRPS = scaled_crps(crps, sd / sqrt(pi)),
                           pit_lower = pit, pit_upper = pit),
       predictive = list(mean = mean, sd = sd))
}

##' Score a test set of counts
##'
##' @param means Joint draws of the means of the test observations, one in
##'   each row.
##' @noRd
score_counts <- function(y, units_m, means, parameters, traits) {
  check_count_means(means, units_m, traits)
  distributions <- lapply(seq_along(y), function(j) {
    count_distribution(y[j], units_m[j], means[j, ], parameters, traits)
  })
  value <- function(f) vapply(seq_along(y), f, numeric(1))
  crps <- value(function(j) count_crps(y[j], distributions[[j]]))
  expected <- value(function(j) count_expected_crps(distributions[[j]]))
  below <- value(function(j) count_distribution_at(distributions[[j]], y[j] - 1))
  at <- value(function(j) count_distribution_at(distributions[[j]], y[j]))
  list(scores = data.frame(CRPS = crps, SCRPS = scaled_crps(crps, expected),
                           pit_lower = below, pit_upper = at),
       predictive = list(distributions = distributions, units_m = units_m,
                         parameters = parameters))
}

##' The scaled CRPS, from the CRPS and its expected value
##'
##' The negative of the SCRPS of Bolin and Wallin (2023), so that lower is
##' better as for the CRPS. It is missing for a predictive distribution with
##' no spread.
##' @noRd
scaled_crps <- function(crps, expected) {
  ifelse(expected > 0,
         0.5 * (1 + crps / expected + log(2 * pmax(expected, 0))), NA_real_)
}

##' @title Summarise a cross-validation
##'
##' @description Gives the average scores of each model over all test
##'   observations, their differences from a reference model, and the areas of
##'   the calibration curves.
##'
##' @param object The result of [assess_prediction()].
##' @param reference Optional name of the model that the others are compared
##'   with.
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##'
##' @details
##' Lower is better for both scores, and a negative difference means a better
##' score than that of the reference model. A `PIT_area` below `PIT_area_limit` is
##' one that a correct model would give (see [assess_prediction()]).
##'
##' An observation in several test sets is first given the average of its
##' scores, so that it counts once. The scores of a model are then averaged
##' over the observations. The differences from the reference are those of
##' the observations that both models have scores for, and their standard error is
##' the standard deviation of these differences divided by the square root of
##' their number. It treats the observations as independent, so it is only
##' approximate, and it is too small when the differences at neighbouring
##' locations are alike.
##'
##' The scores of each test set are in the `scores` element of the result of
##' [assess_prediction()].
##'
##' @return A data frame of class `summary.RiskMap_cross_validation` with one
##'   row for each model, named after it, and the columns `CRPS` and `SCRPS`,
##'   `PIT_area`, `PIT_area_expected` and `PIT_area_limit`, and, when
##'   `reference` is given, `CRPS_difference`, `CRPS_se`, `SCRPS_difference` and
##'   `SCRPS_se`. Its
##'   attributes include `reference`, `refit`, `n_observations` and `n_failed`,
##'   the number of failed test-set assessments.
##'
##' @seealso [assess_prediction()], and [to_table()] for a table to render.
##'
##' @export
##' @method summary RiskMap_cross_validation
summary.RiskMap_cross_validation <- function(object, reference = NULL, ...) {
  check_dots_unused(..., .hints = c(
    view_all = "The scores of each test set are in the {.code scores} element of the result of {.fn assess_prediction}."
  ))
  models <- object$calibration$area$model
  check_assessed_model(reference, models, "reference")
  scores <- lapply(stats::setNames(nm = c("CRPS", "SCRPS")), function(score) {
    observation_scores(object, score)
  })
  columns <- lapply(scores, function(values) {
    out <- data.frame(mean = mean_scores(values))
    if (!is.null(reference)) {
      differences <- values - values[, reference]
      out$difference <- mean_scores(differences)
      out$se <- apply(differences, 2L, function(x) {
        x <- x[is.finite(x)]
        if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
      })
    }
    out
  })
  out <- do.call(cbind, unname(Map(function(column, score) {
    names(column) <- c(score, if (!is.null(reference)) {
      paste0(score, c("_difference", "_se"))
    })
    column
  }, columns, names(columns))))
  out <- cbind(out, object$calibration$area[c("PIT_area", "PIT_area_expected",
                                              "PIT_area_limit")])
  row.names(out) <- models
  attr(out, "reference") <- reference
  attr(out, "refit") <- object$refit
  attr(out, "n_observations") <- nrow(scores$CRPS)
  attr(out, "n_failed") <- nrow(object$failures)
  attr(out, "counts") <- all(families_are_discrete(object$family))
  class(out) <- c("summary.RiskMap_cross_validation", "data.frame")
  out
}

##' Whether the families of the models are discrete
##' @noRd
families_are_discrete <- function(families) {
  vapply(families, function(family) family_traits(family)$discrete, logical(1))
}

##' Check that a model is one of those assessed
##' @noRd
check_assessed_model <- function(model, models, arg, several = FALSE,
                                 call = rlang::caller_env()) {
  if (is.null(model)) return(invisible(NULL))
  if (!is.character(model) || anyNA(model) || !length(model) ||
      (!several && length(model) > 1L)) {
    what <- if (several) "names of models" else "the name of one model"
    abort_input(c(paste0("{.arg {arg}} must be ", what, " of the cross-validation."), "i" = "Its models are {.val {models}}."), arg = arg, call = call)
  }
  absent <- setdiff(model, models)
  if (length(absent)) {
    abort_input(c("{.arg {arg}} must name {cli::qty(length(model))}{?a model/models} of the cross-validation.", "x" = "{.val {absent}} {cli::qty(length(absent))}{?is/are} not among them.", "i" = "Its models are {.val {models}}."), arg = arg, call = call)
  }
  invisible(model)
}

##' The scores of the observations, averaged over the test sets each is in
##'
##' @return A matrix with a row for each test observation, named after its row
##'   of the data, and a column for each model.
##' @noRd
observation_scores <- function(object, score) {
  scores <- object$scores
  models <- object$calibration$area$model
  rows <- sort(unique(scores$row))
  out <- tapply(scores[[score]],
                list(factor(scores$row, levels = rows),
                     factor(scores$model, levels = models)),
                function(values) {
                  values <- values[is.finite(values)]
                  if (length(values)) mean(values) else NA_real_
                })
  matrix(out, nrow = length(rows), dimnames = list(rows, models))
}

##' Column means of the scores there are
##' @noRd
mean_scores <- function(values) {
  apply(values, 2L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  })
}

##' @title Print the summary of a cross-validation
##'
##' @description Prints the scores of each model averaged over all test
##'   observations, with what they mean, from
##'   [summary()][summary.RiskMap_cross_validation()].
##'
##' @param x A summary of class `summary.RiskMap_cross_validation`.
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored. Some columns of the summary print as a data frame, and
##'   the arguments go to its print method.
##'
##' @return Invisibly returns `x`.
##' @export
##' @method print summary.RiskMap_cross_validation
print.summary.RiskMap_cross_validation <- function(x, ...) {
  # Some of the columns print as a data frame.
  if (!is_whole_cv_summary(x)) return(NextMethod())
  check_dots_unused(...)
  out_print(format_cv_summary(x))
  invisible(x)
}

# As for knit_print.RiskMap_extrapolation(), so that a document that prints
# data frames as tables prints the summary as the console does.
##' @exportS3Method knitr::knit_print
knit_print.summary.RiskMap_cross_validation <- function(x, ...) {
  if (!is_whole_cv_summary(x)) return(NextMethod())
  knitr::normal_print(x)
}

##' Whether a summary of a cross-validation has the columns and attributes that
##' its print uses, which a subset of its columns lost
##' @noRd
is_whole_cv_summary <- function(x) {
  all(c("CRPS", "SCRPS", "PIT_area", "PIT_area_expected",
        "PIT_area_limit") %in% names(x)) &&
    is.numeric(attr(x, "n_observations", exact = TRUE))
}

##' The lines that print the summary of a cross-validation
##' @noRd
format_cv_summary <- function(x) {
  c(out_rule("Cross-validation scores"),
    format_cv_settings(attr(x, "refit"), attr(x, "n_failed")),
    "",
    format_cv_averages(x))
}

##' How the parameters were set and how many assessments failed
##' @noRd
format_cv_settings <- function(refit, n_failed, before = NULL) {
  parameters <- if (isTRUE(refit)) {
    "Parameters are re-estimated in each training set."
  } else if (isFALSE(refit)) {
    "Parameters are fixed at the estimates from all the data."
  }
  failed <- if (isTRUE(n_failed > 0)) {
    c("!" = paste0(n_failed, " test-set assessment",
                   if (n_failed == 1L) " failed, and its scores are" else
                     "s failed, and their scores are", " missing."))
  }
  sentences <- c(before, parameters)
  c(if (length(sentences)) out_text(paste(sentences, collapse = " ")),
    out_bullets(failed))
}

##' The scores of each model over all test observations, with what they mean
##' @noRd
format_cv_averages <- function(x) {
  reference <- attr(x, "reference")
  counts <- attr(x, "counts")
  notes <- c(
    paste0("Lower is better for the CRPS and the SCRPS. The CRPS is in the ",
           "units of the response", if (isTRUE(counts)) ", which are counts",
           "."),
    "The PIT area of a correct model is below the limit with a probability of 95%.",
    if (!is.null(reference)) {
      paste0("Differences are from ", reference, ", and negative ones are in ",
             "favour of the other model. Their standard errors, in brackets, ",
             "treat the observations as independent.")
    }
  )
  n_observations <- attr(x, "n_observations")
  heading <- paste0("Scores over ", format_count(n_observations),
                    " test observation", if (n_observations != 1L) "s")
  c(out_heading(heading),
    out_capture(print(cv_summary_table(x, digits = 4L), row.names = FALSE,
                      right = TRUE)),
    out_bullets(stats::setNames(notes, rep("i", length(notes)))))
}

##' The summary of a cross-validation as a table of text
##'
##' One column for each score, followed by its difference from the reference
##' with its standard error in brackets, and the area of the calibration curve
##' with its limit.
##' @noRd
cv_summary_table <- function(x, digits) {
  reference <- attr(x, "reference")
  values <- as.data.frame(unclass(x), row.names = row.names(x))
  blank_missing <- function(values) {
    ifelse(is.finite(values), format_compared_values(values, digits), "NA")
  }
  columns <- list(Model = row.names(x))
  for (score in c("CRPS", "SCRPS")) {
    columns[[score]] <- blank_missing(values[[score]])
    if (is.null(reference)) next
    difference <- values[[paste0(score, "_difference")]]
    se <- values[[paste0(score, "_se")]]
    formatted <- format_compared_values(c(difference, se), digits)
    formatted <- matrix(formatted, ncol = 2L)
    text <- paste0(ifelse(difference > 0, "+", ""), formatted[, 1L],
                   " (", formatted[, 2L], ")")
    text[!is.finite(difference)] <- "NA"
    text[row.names(x) == reference] <- ""
    columns[[paste(score, "difference")]] <- text
  }
  areas <- cbind(values$PIT_area, values$PIT_area_limit)
  formatted <- areas
  formatted[] <- format_compared_values(areas, digits)
  formatted[!is.finite(areas)] <- "NA"
  columns[["PIT area"]] <- formatted[, 1L]
  columns[["Limit"]] <- formatted[, 2L]
  as.data.frame(columns, check.names = FALSE, optional = TRUE)
}

##' @title Print a cross-validation
##' @description Gives the models, the test sets, how the parameters were set
##'   and the scores of each model over all test observations.
##' @param x The result of [assess_prediction()].
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##' @return Invisibly returns `x`.
##' @seealso [summary()][summary.RiskMap_cross_validation()] compares the
##'   models with one of them, and [plot_cross_validation()] draws the
##'   calibration and the scores.
##' @method print RiskMap_cross_validation
##' @export
print.RiskMap_cross_validation <- function(x, ...) {
  check_dots_unused(...)
  check_class(x, "RiskMap_cross_validation")
  out_print(format_cross_validation(x))
  invisible(x)
}

##' The number of test sets and their sizes, as in "2 test sets of 80 and 20
##' observations" or "10 test sets of 50 observations each"
##' @noRd
describe_test_sets <- function(sizes) {
  n <- length(sizes)
  sets <- paste0(n, if (n == 1L) " test set" else " test sets")
  if (length(unique(sizes)) == 1L) {
    return(paste0(sets, " of ", format_count(sizes[1L]),
                  if (sizes[1L] == 1L) " observation" else " observations",
                  if (n > 1L) " each"))
  }
  observations <- if (n <= 3L) {
    paste(paste(format_count(utils::head(sizes, -1L)), collapse = ", "),
          "and", format_count(utils::tail(sizes, 1L)))
  } else {
    paste(format_count(min(sizes)), "to", format_count(max(sizes)))
  }
  paste0(sets, " of ", observations, " observations")
}

##' The lines that print the results of a cross-validation
##' @noRd
format_cross_validation <- function(x) {
  scores <- summary(x)
  n_models <- nrow(scores)
  sizes <- lengths(lapply(x$test_sets$sets, `[[`, "test"))
  setting <- paste0(n_models, if (n_models == 1L) " model" else " models",
                    " assessed on ", describe_test_sets(sizes), ".")
  hint <- paste0(
    if (n_models > 1L) {
      paste0(out_fun("summary", "summary.RiskMap_cross_validation"),
             " with ", "`reference` compares the models with one of them, and ")
    },
    out_fun("plot", "plot.RiskMap_cross_validation"),
    " draws the calibration curves or maps the scores."
  )
  c(out_rule("Cross-validation"),
    format_cv_settings(x$refit, nrow(x$failures), before = setting),
    "",
    format_cv_averages(scores),
    "",
    out_bullets(c(i = hint)))
}

##' @rdname to_table
##' @export
to_table.RiskMap_cross_validation <- function(object, digits = 3,
                                              col_names = NULL,
                                              reference = NULL, ...) {
  to_table(summary(object, reference = reference), digits = digits,
           col_names = col_names, ...)
}

##' @rdname to_table
##' @export
to_table.summary.RiskMap_cross_validation <- function(object, digits = 3,
                                                      col_names = NULL,
                                                      ...) {
  check_positive_integer(digits, "digits", allow_zero = TRUE)
  if ("reference" %in% names(list(...))) {
    abort_input("{.arg reference} can't be given for a summary. Give it to {.fn summary}.", arg = "reference")
  }
  if (!is_whole_cv_summary(object)) {
    abort_input("{.arg object} must be the whole summary of a cross-validation, not a part of it.", arg = "object")
  }
  render_results_table(cv_summary_table(object, digits), row_names = FALSE,
                       col_names = col_names, ...)
}
