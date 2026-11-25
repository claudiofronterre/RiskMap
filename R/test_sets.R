# Test sets for cross-validation, for assess_prediction().

##' @title Create test sets for cross-validation
##'
##' @description
##' Splits the locations of a dataset into training and test sets for
##' [assess_prediction()]. Each test set is predicted from the rest of the
##' data, its training set. Records at the same location are always in the
##' same set.
##'
##' @param data An `sf` object with point geometries, the data the models are
##'   fitted to, or a fit from [glgpm()].
##' @param method How the test sets are made.
##'   - `"leave_one_out"` (the default) has one location in each test set.
##'   - `"random"` divides the locations at random into `n_sets` sets of about
##'     equal size.
##'   - `"cluster"` divides the locations into `n_sets` spatial clusters, by
##'     k-means on their coordinates.
##'   - `"spaced"` draws `n_sets` test sets of `size` locations that are at
##'     least `min_dist` apart.
##' @param n_sets The number of test sets, 10 by default. With
##'   `"leave_one_out"` it is the number of locations, and it can't be given.
##' @param size For `"spaced"` only, the number of locations in each test set.
##' @param min_dist For `"spaced"` only, the smallest distance between two
##'   locations of a test set, in the units of distance of the model, which
##'   are kilometres for an `sf` object and `distance_units` for a fit.
##' @param repeats The number of times that a `"random"` or `"cluster"`
##'   division is repeated, which gives `repeats` times `n_sets` test sets.
##' @param seed Optional non-negative integer that makes the test sets
##'   reproducible. The caller's random number state is restored on exit.
##'
##' @details
##' Choose the test sets so that the test locations are about as far from
##' their training data as the locations to be predicted are from all the data
##' (Milà *et al.* 2022). [plot_test_sets()] with `type = "distance"` compares
##' these distances. Leave-one-out and random test sets suit predictions
##' between the locations of the data, and clusters suit predictions in areas
##' without data. The distance `min_dist` of `"spaced"` is between the
##' locations of a test set, which keeps them spread over the region. It isn't
##' a distance from the training data.
##'
##' Distances and clusters are computed on projected coordinates, as in
##' [glgpm()]. Data in longitude and latitude are projected to the CRS that
##' [propose_crs()] gives, and a fit has the CRS and the units of distance of
##' its model.
##'
##' The clusters of k-means are the best of 20 random starts. For `"spaced"`,
##' the locations are visited in random order and one is kept if it is at
##' least `min_dist` from those kept so far, until `size` locations are kept
##' or every location has been visited. A test set can therefore have fewer
##' than `size` locations, and a warning says so.
##'
##' # Test sets made by another package
##'
##' [assess_prediction()] also takes a list with one element for each test
##' set, the row numbers of the test observations or a list of the training
##' rows, `train`, and the test rows, `test`. Observations in neither are left
##' out, for example those in a buffer around the test locations. For the
##' folds of the rsample or the spatialsample package, this list is
##'
##' ```r
##' lapply(folds$splits, function(split) {
##'   list(train = as.integer(split),
##'        test = as.integer(split, data = "assessment"))
##' })
##' ```
##'
##' and for `folds_list` of the blockCV package it is
##'
##' ```r
##' lapply(folds$folds_list, function(fold) {
##'   list(train = fold[[1]], test = fold[[2]])
##' })
##' ```
##'
##' @return An object of class `RiskMap_test_sets`, a list with
##' \describe{
##'   \item{sets}{A list with one element for each test set, each a list of
##'     the row numbers of the training data, `train`, and of the test data,
##'     `test`.}
##'   \item{repetition}{The repeat each test set belongs to.}
##'   \item{method}{The method, or `"user"` for test sets given to
##'     [assess_prediction()] as a list.}
##'   \item{geometry}{The locations of the data, in their own CRS.}
##'   \item{coordinates, distance_units}{The projected coordinates of the rows
##'     of the data, in the units of distance.}
##'   \item{crs}{The projected CRS.}
##'   \item{location}{The distinct location of each row of the data.}
##'   \item{min_dist, size}{As given.}
##' }
##'
##' @references
##' Milà, C., Mateu, J., Pebesma, E. & Meyer, H. (2022). Nearest neighbour
##' distance matching Leave-One-Out Cross-Validation for map validation.
##' *Methods in Ecology and Evolution*, 13(6), 1304-1316.
##' \doi{10.1111/2041-210X.13851}.
##'
##' @seealso [plot_test_sets()], which `plot()` of the test sets also calls,
##'   and [assess_prediction()], which scores the predictions.
##'
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:100, ]
##'
##' # One location at a time
##' create_test_sets(example_data)
##'
##' # Five spatial clusters
##' clusters <- create_test_sets(example_data, method = "cluster", n_sets = 5,
##'                              seed = 1)
##' clusters
##' plot(clusters)
##'
##' # Ten test sets of five locations at least 50 km apart
##' create_test_sets(example_data, method = "spaced", n_sets = 10, size = 5,
##'                  min_dist = 50, seed = 1)
##' @export
create_test_sets <- function(data,
                             method = c("leave_one_out", "random", "cluster",
                                        "spaced"),
                             n_sets = NULL,
                             size = NULL,
                             min_dist = NULL,
                             repeats = 1,
                             seed = NULL) {
  check_supplied(data)
  method <- match_choice(method)
  check_test_set_arguments(method, n_sets, size, min_dist, repeats)
  check_positive_integer(seed, "seed", allow_null = TRUE, allow_zero = TRUE)
  located <- locate_test_data(
    data, announce = method %in% c("cluster", "spaced")
  )
  n_locations <- nrow(located$locations)
  if (n_locations < 2L) {
    abort_input("{.arg data} must have at least two locations, to have training and test data.", arg = "data")
  }
  n_sets <- if (method == "leave_one_out") n_locations else n_sets %||% 10
  check_test_set_number(method, n_sets, size, n_locations)

  if (!is.null(seed)) {
    restore_seed <- preserve_random_seed()
    on.exit(restore_seed(), add = TRUE)
    set.seed(seed)
  }
  divide <- switch(
    method,
    leave_one_out = function() seq_len(n_locations),
    random = function() sample(rep_len(seq_len(n_sets), n_locations)),
    cluster = function() cluster_locations(located$locations, n_sets)
  )
  locations <- if (method == "spaced") {
    lapply(seq_len(n_sets), function(k) {
      spaced_locations(located$locations, size, min_dist)
    })
  } else {
    unlist(lapply(seq_len(repeats), function(k) {
      unname(split(seq_len(n_locations), divide()))
    }), recursive = FALSE)
  }
  if (method == "spaced") warn_short_test_sets(lengths(locations), size)

  sets <- lapply(locations, function(tested) {
    in_test <- located$location %in% tested
    list(train = which(!in_test), test = which(in_test))
  })
  new_test_sets(sets, method, located,
                repetition = rep(seq_len(repeats), each = n_sets),
                min_dist = min_dist, size = size)
}

##' The object returned by `create_test_sets()`
##'
##' @param sets The test sets, each a list of row numbers `train` and `test`.
##' @param located The located data (`locate_test_data()`).
##' @noRd
new_test_sets <- function(sets, method, located,
                          repetition = rep(1L, length(sets)),
                          min_dist = NULL, size = NULL) {
  out <- list(sets = sets, repetition = as.integer(repetition),
              method = method, geometry = located$geometry,
              coordinates = located$coordinates,
              distance_units = located$distance_units, crs = located$crs,
              location = located$location, min_dist = min_dist, size = size)
  class(out) <- "RiskMap_test_sets"
  stamp_version(out)
}

##' Check the arguments of `create_test_sets()` that depend on the method
##' @noRd
check_test_set_arguments <- function(method, n_sets, size, min_dist, repeats,
                                     call = rlang::caller_env()) {
  check_positive_integer(n_sets, "n_sets", allow_null = TRUE, call = call)
  check_positive_integer(size, "size", allow_null = TRUE, call = call)
  if (!is.null(min_dist)) {
    check_positive_number(min_dist, "", arg = "min_dist", call = call)
  }
  check_positive_integer(repeats, "repeats", call = call)
  if (method == "leave_one_out" && !is.null(n_sets)) {
    abort_input("{.arg n_sets} can't be given when {.arg method} is {.val leave_one_out}, which has one test set for each location.", arg = "n_sets", call = call)
  }
  spaced <- method == "spaced"
  for (arg in c("size", "min_dist")) {
    given <- !is.null(if (arg == "size") size else min_dist)
    if (spaced && !given) {
      abort_input("{.arg {arg}} must be given when {.arg method} is {.val spaced}.", arg = arg, call = call)
    }
    if (!spaced && given) {
      abort_input("{.arg {arg}} can be given only when {.arg method} is {.val spaced}.", arg = arg, call = call)
    }
  }
  if (repeats > 1 && !method %in% c("random", "cluster")) {
    hint <- if (spaced) {
      "Give the number of test sets in {.arg n_sets}."
    } else {
      "Leave-one-out test sets are the same every time."
    }
    abort_input(c("{.arg repeats} can be above 1 only when {.arg method} is {.val random} or {.val cluster}.", "i" = hint), arg = "repeats", call = call)
  }
  invisible(TRUE)
}

##' Check the number and the size of the test sets against the locations
##' @noRd
check_test_set_number <- function(method, n_sets, size, n_locations,
                                  call = rlang::caller_env()) {
  if (method %in% c("random", "cluster")) {
    if (n_sets < 2) {
      abort_input("{.arg n_sets} must be at least 2, so that every test set has training data.", arg = "n_sets", call = call)
    }
    if (n_sets > n_locations) {
      abort_input("{.arg n_sets} must be at most the number of locations of {.arg data}, {n_locations}, not {n_sets}.", arg = "n_sets", call = call)
    }
  }
  if (method == "spaced" && size >= n_locations) {
    abort_input("{.arg size} must be below the number of locations of {.arg data}, {n_locations}, so that there are training data.", arg = "size", call = call)
  }
  invisible(TRUE)
}

##' The locations of the data on the plane in which distances are computed
##'
##' As in `glgpm()`, the projected CRS of the data, or the one proposed for
##' data in longitude and latitude, with distances in kilometres. A fit has
##' the CRS and the units of its model.
##' @param data An `sf` object or a fit.
##' @param announce Whether to say how data in longitude and latitude are
##'   projected, when the test sets depend on it.
##' @return A list with the `geometry` of the rows in their own CRS, their
##'   projected `coordinates` in `distance_units`, the projected `crs`, the
##'   distinct `locations` and the `location` of each row.
##' @noRd
locate_test_data <- function(data, announce = FALSE,
                             call = rlang::caller_env()) {
  if (inherits(data, "RiskMap")) {
    distance_units <- data$distance_units
    geometry <- st_geometry(data$data)
    projected <- geometry
  } else {
    if (!inherits(data, "sf")) {
      abort_input("{.arg data} must be an {.cls sf} object or a fit from {.fn glgpm}, not {.obj_type_friendly {data}}.", arg = "data", call = call)
    }
    check_data(data, arg = "data")
    distance_units <- "km"
    geometry <- st_geometry(data)
    projected <- if (isTRUE(st_is_longlat(geometry))) {
      project_test_data(geometry, announce)
    } else {
      check_crs_distortion(geometry, NULL, "data")
      geometry
    }
  }
  coordinates <- coordinates_in_units(projected, distance_units)
  indexed <- index_coordinate_rows(coordinates)
  list(geometry = geometry, coordinates = coordinates,
       distance_units = distance_units, crs = st_crs(projected),
       locations = indexed$coordinates, location = indexed$index)
}

##' Locations in longitude and latitude, in the CRS that `propose_crs()` gives
##' @noRd
project_test_data <- function(geometry, announce) {
  choice <- choose_crs(geometry)
  if (announce) {
    target <- choice$label
    riskmap_inform(c("{.arg data} are in longitude and latitude, so distances between their locations are computed in {target}, in km.", "i" = "Give a fit as {.arg data} for the projection and the units of its model."))
  }
  if (choice$error > crs_distortion_limit) warn_region_too_wide(choice, "data")
  st_transform(geometry, choice$crs)
}

##' Spatial clusters of locations, by k-means
##'
##' The best of 20 random starts. With as many clusters as locations, where
##' `kmeans()` stops, each location is a cluster. The clusters are numbered
##' from north to south, so that the numbers mean the same at every seed that
##' finds the same clusters.
##' @param coordinates The projected coordinates of the distinct locations.
##' @return The cluster of each location.
##' @noRd
cluster_locations <- function(coordinates, n_clusters) {
  n_locations <- nrow(coordinates)
  if (n_clusters == n_locations) {
    cluster <- seq_len(n_locations)
  } else {
    # A start that stops before it converges still divides the locations.
    cluster <- suppressWarnings(
      stats::kmeans(coordinates, n_clusters, iter.max = 100L, nstart = 20L)
    )$cluster
  }
  northing <- tapply(coordinates[, 2L], cluster, mean)
  match(cluster, as.integer(names(northing))[order(-northing)])
}

##' Locations drawn at random that are at least `min_dist` apart
##'
##' Visits the locations in random order and keeps one if it is at least
##' `min_dist` from those kept so far, until `size` are kept or none is left.
##' @return The row numbers of the locations kept.
##' @noRd
spaced_locations <- function(coordinates, size, min_dist) {
  available <- sample.int(nrow(coordinates))
  kept <- integer(0)
  while (length(kept) < size && length(available)) {
    candidate <- available[1L]
    kept <- c(kept, candidate)
    distance <- sqrt(
      (coordinates[available, 1L] - coordinates[candidate, 1L])^2 +
        (coordinates[available, 2L] - coordinates[candidate, 2L])^2
    )
    available <- available[distance >= min_dist & available != candidate]
  }
  kept
}

##' Warn about spaced test sets with fewer locations than asked for
##' @noRd
warn_short_test_sets <- function(n_kept, size) {
  n_short <- sum(n_kept < size)
  if (n_short == 0L) return(invisible(NULL))
  n_sets <- length(n_kept)
  riskmap_warn(c(
    "{n_short} of {n_sets} test set{cli::qty(n_sets)}{?s} {cli::qty(n_short)}{?has/have} fewer than {.arg size} = {size} locations.",
    "i" = "No other location is at least {.arg min_dist} from the ones drawn. The smallest test set has {min(n_kept)} location{?s}."
  ))
}

##' @title Print test sets
##' @description Gives the number and the sizes of the test sets and how far
##'   the test locations are from their training data.
##' @param x Test sets from [create_test_sets()].
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##' @return Invisibly returns `x`.
##' @method print RiskMap_test_sets
##' @export
print.RiskMap_test_sets <- function(x, ...) {
  check_dots_unused(...)
  check_class(x, "RiskMap_test_sets")
  out_print(format_test_sets(x))
  invisible(x)
}

##' The lines that print test sets
##' @noRd
format_test_sets <- function(x) {
  n_rows <- length(x$location)
  n_locations <- max(x$location)
  data_text <- paste0(
    format_count(n_rows), if (n_rows == 1L) " observation" else " observations",
    if (n_locations < n_rows) {
      paste0(" at ", format_count(n_locations), " locations")
    }
  )
  sizes <- lengths(lapply(x$sets, `[[`, "test"))
  n_repeats <- max(x$repetition)
  how <- switch(
    x$method,
    leave_one_out = "Leave-one-out, with one location in each test set.",
    random = "Locations divided at random.",
    cluster = "Spatial clusters of locations, by k-means.",
    spaced = paste0("Locations at least ", format(x$min_dist), " ",
                    x$distance_units, " apart."),
    user = "Test sets given as rows of the data."
  )
  repeated <- if (n_repeats > 1L) {
    paste0(" The division is repeated ", n_repeats, " times.")
  }
  distance <- test_set_distances(x)$distance
  quantiles <- format(round(stats::quantile(distance, c(0.5, 1)), 1),
                      trim = TRUE)
  c(out_rule("Test sets for cross-validation"),
    out_text(paste0(describe_test_sets(sizes), " from ", data_text, ". ",
                    how, repeated)),
    out_text(paste0("The distance from a test location to its nearest ",
                    "training location is ", quantiles[1L], " ",
                    x$distance_units, " in the median and at most ",
                    quantiles[2L], " ", x$distance_units, ".")),
    "",
    out_bullets(c(i = paste0(
      out_fun("plot", "plot.RiskMap_test_sets"),
      " maps the test sets, and ", out_fun("assess_prediction"),
      " scores the predictions of models on them."
    ))))
}

##' The distances of the test locations from their training data
##'
##' @param test_sets A `RiskMap_test_sets` object.
##' @return A data frame with one row for each test set and distinct test
##'   location, with the `distance` to the nearest training location in the
##'   units of the test sets.
##' @noRd
test_set_distances <- function(test_sets) {
  coordinates <- test_sets$coordinates
  location <- test_sets$location
  pieces <- lapply(seq_along(test_sets$sets), function(k) {
    set <- test_sets$sets[[k]]
    tested <- set$test[!duplicated(location[set$test])]
    training <- set$train[!duplicated(location[set$train])]
    data.frame(test_set = k, row = tested,
               distance = nearest_distances(
                 coordinates[tested, , drop = FALSE],
                 coordinates[training, , drop = FALSE]
               ))
  })
  do.call(rbind, pieces)
}

##' The distance from each row of `from` to the nearest row of `to`
##'
##' In blocks, so that many prediction locations don't need one large matrix.
##' @noRd
nearest_distances <- function(from, to, block_size = 5000L) {
  blocks <- split(seq_len(nrow(from)),
                  (seq_len(nrow(from)) - 1L) %/% block_size)
  unlist(lapply(blocks, function(rows) {
    distances <- cross_distances(from[rows, , drop = FALSE], to)
    do.call(pmin, as.data.frame(as.matrix(distances)))
  }), use.names = FALSE)
}

##' Test sets given by the user, as an object of class `RiskMap_test_sets`
##'
##' @param test_sets A list with one element for each test set, the row numbers
##'   of its test observations or a list of row numbers `train` and `test`.
##' @param fit The fit whose data the rows belong to.
##' @noRd
read_test_sets <- function(test_sets, fit, call = rlang::caller_env()) {
  located <- locate_test_data(fit)
  n_rows <- length(located$location)
  if (inherits(test_sets, "RiskMap_test_sets")) {
    check_test_set_rows(test_sets, fit, located, call)
    return(test_sets)
  }
  if (is.character(test_sets)) {
    abort_input(c("{.arg test_sets} must be made by {.fn create_test_sets}, or be a list of row numbers, not {.val {test_sets}}.", "i" = "Use {.code create_test_sets(data, method = \"cluster\")} for spatial clusters, for example."), arg = "test_sets", call = call)
  }
  if (!is.list(test_sets) || is.data.frame(test_sets)) {
    hint <- if (is.numeric(test_sets) && !is.matrix(test_sets)) {
      "Put the rows of one test set in a list of its own, as {.code list(rows)}."
    } else if (is.matrix(test_sets)) {
      "For a matrix with a 1 for each test observation, use {.code apply(x == 1, 2, which, simplify = FALSE)}."
    }
    abort_input(c("{.arg test_sets} must be made by {.fn create_test_sets}, or be a list of row numbers, not {.obj_type_friendly {test_sets}}.", "i" = hint), arg = "test_sets", call = call)
  }
  if (any(c("train", "test") %in% names(test_sets))) {
    abort_input(c("{.arg test_sets} must be a list with one element for each test set.", "i" = "Put one test set given by {.code train} and {.code test} in a list of its own, as {.code list(list(train = , test = ))}."), arg = "test_sets", call = call)
  }
  if (!length(test_sets)) {
    abort_input("{.arg test_sets} must have at least one test set.", arg = "test_sets", call = call)
  }
  sets <- lapply(seq_along(test_sets), function(k) {
    read_one_test_set(test_sets[[k]], k, n_rows, call)
  })
  new_test_sets(sets, "user", located)
}

##' One test set given by the user, as a list of `train` and `test` rows
##' @noRd
read_one_test_set <- function(set, index, n_rows, call) {
  name <- sprintf("test_sets[[%d]]", index)
  if (is.list(set)) {
    if (any(c("in_id", "out_id") %in% names(set))) {
      abort_input("The training and test rows of {.code {name}} are now called {.code train} and {.code test}.", arg = "test_sets", call = call)
    }
    if (is.null(set$test)) {
      abort_input("{.code {name}} must be a vector of test rows or a list with {.code train} and {.code test}.", arg = "test_sets", call = call)
    }
    test <- check_test_rows(set$test, paste0("test rows of ", name), n_rows, call)
    train <- if (is.null(set$train)) {
      setdiff(seq_len(n_rows), test)
    } else {
      check_test_rows(set$train, paste0("training rows of ", name), n_rows, call)
    }
    if (length(intersect(train, test))) {
      abort_input("The {.code train} and {.code test} rows of {.code {name}} must not overlap.", arg = "test_sets", call = call)
    }
  } else if (is.numeric(set)) {
    test <- check_test_rows(set, paste0("test rows of ", name), n_rows, call)
    train <- setdiff(seq_len(n_rows), test)
  } else {
    abort_input("{.code {name}} must be a vector of test rows or a list with {.code train} and {.code test}, not {.obj_type_friendly {set}}.", arg = "test_sets", call = call)
  }
  if (!length(train)) {
    abort_input("{.code {name}} must leave at least one observation for training.", arg = "test_sets", call = call)
  }
  list(train = as.integer(train), test = as.integer(test))
}

##' Check a vector of row numbers of the data
##' @noRd
check_test_rows <- function(rows, what, n_rows, call) {
  problem <- if (!is.numeric(rows) || anyNA(rows)) {
    "numeric row numbers without missing values"
  } else if (!length(rows)) {
    "at least one row number"
  } else if (any(rows != round(rows)) || any(rows < 1) || any(rows > n_rows)) {
    "row numbers between 1 and {n_rows}"
  } else if (anyDuplicated(rows)) {
    "row numbers without repeats"
  }
  if (!is.null(problem)) {
    abort_input(paste0("The {what} must be ", problem, "."), arg = "test_sets", call = call)
  }
  rows
}

##' Check that test sets belong to the data of a fit
##'
##' The same number of rows at the same locations, to within a metre, which
##' allows for a transformation from another CRS.
##' @noRd
check_test_set_rows <- function(test_sets, fit, located, call) {
  if (length(test_sets$location) != length(located$location)) {
    abort_input("{.arg test_sets} must be made for the data of the models, which have {length(located$location)} rows, not {length(test_sets$location)}.", arg = "test_sets", call = call)
  }
  geometry <- test_sets$geometry
  if (!identical(geometry, located$geometry)) {
    geometry <- st_transform(geometry, st_crs(fit$data))
    coordinates <- coordinates_in_units(geometry, "m")
    expected <- coordinates_in_units(located$geometry, "m")
    if (!all(abs(coordinates - expected) < 1)) {
      abort_input("{.arg test_sets} must be made for the data of the models, with the same rows in the same order.", arg = "test_sets", call = call)
    }
  }
  invisible(TRUE)
}
