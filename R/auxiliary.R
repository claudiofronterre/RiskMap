`%||%` <- function(a, b) if (!is.null(a)) a else b

##' @title Convex Hull of an sf Object
##'
##' @description Computes the convex hull of an `sf` object, returning the boundaries of the smallest polygon that can enclose all geometries in the input.
##'
##' @param sf_object An `sf` data frame object containing geometries.
##'
##' @return An `sf` object representing the convex hull of the input geometries.
##'
##' @details The convex hull is the smallest convex polygon that encloses all points in the input `sf` object. This function computes the convex hull by first uniting all geometries in the input using `st_union()`, and then applying `st_convex_hull()` to obtain the polygonal boundary. The result is returned as an `sf` object containing the convex hull geometry.
##'
##' @seealso \code{\link[sf]{st_convex_hull}}, \code{\link[sf]{st_union}}
##'
##' @examples
##' library(sf)
##'
##' # Create example sf object
##' points <- st_sfc(st_point(c(0,0)), st_point(c(1,1)), st_point(c(2,2)), st_point(c(0,2)))
##' sf_points <- st_sf(geometry = points)
##'
##' # Calculate the convex hull
##' convex_hull_result <- create_convex_hull(sf_points)
##'
##' # Plot the result
##' plot(sf_points, col = 'blue', pch = 19)
##' plot(convex_hull_result, add = TRUE, border = 'red')
##' @export
create_convex_hull <- function(sf_object) {
  check_supplied(sf_object)
  # Check if the input is an sf object
  if (!inherits(sf_object, "sf")) {
    abort_input("{.arg sf_object} must be an {.cls sf} object, not {.obj_type_friendly {sf_object}}.", arg = "sf_object")
  }

  # Get the geometry from the sf object
  geometry <- st_geometry(sf_object)

  # Compute the convex hull
  convex_hull <- st_convex_hull(st_union(geometry))

  # Return the convex hull as an sf object
  return(st_sf(geometry = convex_hull))
}

##' @title Empirical Logit Transformation for Binomial Data
##' @description Computes the empirical logit transformation for binomial counts.
##' @param y A numeric vector of observed positives.
##' @param m A numeric vector of denominators, i.e. the numbers tested.
##' @return A numeric vector containing the empirical logit transformation
##' \eqn{\log((y + 0.5) / (m - y + 0.5))}.
##' @details The empirical logit is often used as a finite transformation for
##' binomial data, including cases where \eqn{y = 0} or \eqn{y = m}.
##' @examples
##' y <- c(0, 3, 7, 10)
##' m <- c(10, 10, 10, 10)
##' elogit(y, m)
##' @export
elogit <- function(y, m) {
  check_supplied(y)
  check_supplied(m)
  if (!is.numeric(y)) abort_input("{.arg y} must be numeric, not {.obj_type_friendly {y}}.", arg = "y")
  if (!is.numeric(m)) abort_input("{.arg m} must be numeric, not {.obj_type_friendly {m}}.", arg = "m")
  if (anyNA(y)) abort_input("{.arg y} must not have missing values.", arg = "y")
  if (anyNA(m)) abort_input("{.arg m} must not have missing values.", arg = "m")
  sizes <- c(length(y), length(m))
  max_size <- max(sizes)
  if (max_size == 0) {
    return(numeric(0))
  }
  incompatible <- sizes != 1 & max_size %% sizes != 0
  if (any(incompatible)) {
    abort_input("{.arg y} and {.arg m} must have lengths that divide the longer of the two.", arg = "m")
  }
  y <- rep_len(y, max_size)
  m <- rep_len(m, max_size)

  if (any(m <= 0)) abort_input("{.arg m} must be positive.", arg = "m")
  if (any(y < 0)) abort_input("{.arg y} must be non-negative.", arg = "y")
  if (any(y > m)) abort_input("{.arg y} must not exceed {.arg m}.", arg = "y")

  log((y + 0.5) / (m - y + 0.5))
}


##' @title EPSG of the UTM Zone
##' @description Suggests the EPSG code for the UTM zone where the majority of the data falls.
##' @param data An object of class \code{sf} containing the coordinates.
##' @details The function determines the UTM zone and hemisphere where the majority of the data points are located and proposes the corresponding EPSG code.
##' @return An integer indicating the EPSG code of the UTM zone.
##' @examples
##' data("infect_sma")
##' propose_utm(infect_sma)
##'
##' @export
propose_utm <- function(data) {
  check_supplied(data)

  check_data(data, geometry = "any")

  # Transform to WGS84 (EPSG:4326) to ensure coordinates are in lon/lat
  data <- st_transform(data, crs = 4326)

  # Calculate UTM Zone
  utm_z <- floor((st_coordinates(data)[, 1] + 180)/6) + 1
  utm_z_u <- unique(utm_z)

  if (length(utm_z_u) > 1) {
    tab_utm <- table(utm_z)
    if (all(diff(tab_utm) == 0))
      riskmap_warn("The same number of locations falls in more than one UTM zone, so the lowest-numbered of them is used.")
    utm_z_u <- as.numeric(names(which.max(tab_utm)))
  }

  # Determine Hemisphere (fixing the latitude check)
  # latitude 0 (the Equator) is treated as northern hemisphere, per UTM convention
  ns <- ifelse(st_coordinates(data)[, 2] >= 0, 1, -1)  # Use latitude, not longitude
  ns_u <- unique(ns)

  if (length(ns_u) > 1) {
    tab_ns <- table(ns)
    if (all(diff(tab_ns) == 0))
      riskmap_warn("The same number of locations falls north and south of the Equator, so the southern zone is used.")
    ns_u <- as.numeric(names(which.max(tab_ns)))
  }

  # Construct EPSG code for UTM zone
  if (ns_u == 1) {
    out <- as.numeric(paste0(326, utm_z_u))  # Northern Hemisphere
  } else if (ns_u == -1) {
    out <- as.numeric(paste0(327, utm_z_u))  # Southern Hemisphere
  }

  return(out)
}


##' @title Matern Correlation Function
##' @description Computes the Matern correlation function.
##' @param u A vector or matrix of distances between pairs of data locations.
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @param return_sym_matrix A logical value indicating whether to return a symmetric correlation matrix. Defaults to \code{FALSE}.
##' @details The Matern correlation function is defined as
##' \deqn{\rho(u; \phi; \kappa) = (2^{\kappa-1})^{-1}(u/\phi)^\kappa K_{\kappa}(u/\phi)}
##' where \eqn{\phi} and \eqn{\kappa} are the scale and smoothness parameters, and \eqn{K_{\kappa}(\cdot)} denotes the modified Bessel function of the third kind of order \eqn{\kappa}. The parameters \eqn{\phi} and \eqn{\kappa} must be positive.
##' @return A vector or matrix of the same shape as \code{u} with the values of the Matern correlation function for the given distances, if \code{return_sym_matrix=FALSE}. If \code{return_sym_matrix=TRUE}, a symmetric correlation matrix is returned.
##' @examples
##' u <- c(10, 20, 30)
##' matern_correlation(u, phi = 1, kappa = 0.8)
##'
##' @export
matern_correlation <- function(u, phi, kappa, return_sym_matrix = FALSE) {
  check_supplied(u)
  check_supplied(phi)
  check_supplied(kappa)
  if (any(!is.finite(u) | u < 0)) {
    abort_input("{.arg u} must contain only finite, non-negative distances.", arg = "u")
  }
  check_range(phi, min = 0, allow_equal = FALSE, arg = "phi")

  input_dimensions <- if (is.matrix(u)) dim(u) else NULL
  if (is.vector(u))
    names(u) <- NULL
  if (is.matrix(u))
    dimnames(u) <- list(NULL, NULL)
  if (kappa %in% c(0.5, 1.5, 2.5)) {
    uphi <- cpp_half_integer_matern(as.numeric(u), phi, kappa)
    if (!is.null(input_dimensions)) {
      dim(uphi) <- input_dimensions
    }
  } else {
    uphi <- u / phi
    uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf,
                                                      gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)),
                   1)
    uphi[u > 600 * phi] <- 0
  }

  if(return_sym_matrix) {
    n <- (1 + sqrt(1 + 8 * length(uphi))) / 2
    varcov <- matrix(NA, n, n)
    varcov[lower.tri(varcov)] <- uphi
    varcov <- t(varcov)
    varcov[lower.tri(varcov)] <- uphi
    diag(varcov) <- 1
    out <- varcov
  } else {
    out <- uphi
  }
  return(out)
}

##' @title First Derivative with Respect to \eqn{\phi}
##' @description Computes the first derivative of the Matern correlation function with respect to \eqn{\phi}.
##' @param u A distance matrix returned from [stats::dist()]
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @return A matrix with the values of the first derivative of the Matern function with respect to \eqn{\phi} for the given distances.
##' @examples
##' x <- matrix(rnorm(20), nrow = 5)
##' u <- dist(x)
##' matern_gradient_phi(u, phi = 1, kappa = 0.8)
##'
##' @export
##'
matern_gradient_phi <- function(u, phi, kappa) {
  check_supplied(u)
  check_supplied(phi)
  check_supplied(kappa)

  if (!inherits(u, "dist")) {
    abort_input("{.arg u} must be made by {.fn stats::dist}, not {.obj_type_friendly {u}}.", arg = "u")
  }
  check_range(phi, min = 0, allow_equal = FALSE, arg = "phi")

  der.phi <- function(u, phi, kappa) {
    u <- u + 10e-16
    if(kappa == 0.5) {
      out <- (u * exp(-u / phi)) / phi^2
    } else {
      out <- ((besselK(u / phi, kappa + 1) + besselK(u / phi, kappa - 1)) *
                phi^(-kappa - 2) * u^(kappa + 1)) / (2^kappa * gamma(kappa)) -
        (kappa * 2^(1 - kappa) * besselK(u / phi, kappa) * phi^(-kappa - 1) *
           u^kappa) / gamma(kappa)
    }
    out
  }

  n <- attr(u, "Size")
  grad.phi.mat <- matrix(NA, nrow = n, ncol = n)
  ind <- lower.tri(grad.phi.mat)
  grad.phi <- der.phi(as.numeric(u), phi, kappa)
  grad.phi.mat[ind] <-  grad.phi
  grad.phi.mat <- t(grad.phi.mat)
  grad.phi.mat[ind] <-  grad.phi
  diag(grad.phi.mat) <- rep(der.phi(0, phi, kappa), n)
  grad.phi.mat
}

##' @title Second Derivative with Respect to \eqn{\phi}
##' @description Computes the second derivative of the Matern correlation
##' function with respect to \eqn{\phi}.
##' @param u A distance matrix returned from [stats::dist()]
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @return A matrix with the values of the second derivative of the
##' Matern function with respect to \eqn{\phi} for the given distances.
##' @examples
##' x <- matrix(rnorm(20), nrow = 5)
##' u <- dist(x)
##' matern_hessian_phi(u, phi = 1, kappa = 0.8)
##'
##' @export
matern_hessian_phi <- function(u, phi, kappa) {
  check_supplied(u)
  check_supplied(phi)
  check_supplied(kappa)

  if (!inherits(u, "dist")) {
    abort_input("{.arg u} must be made by {.fn stats::dist}, not {.obj_type_friendly {u}}.", arg = "u")
  }
  check_range(phi, min = 0, allow_equal = FALSE, arg = "phi")

  der2.phi <- function(u, phi, kappa) {
    u <- u + 10e-16
    if(kappa == 0.5) {
      out <- (u * (u - 2 * phi) * exp(-u / phi)) / phi^4
    } else {
      bk <- besselK(u / phi, kappa)
      bk.p1 <- besselK(u / phi, kappa + 1)
      bk.p2 <- besselK(u / phi, kappa + 2)
      bk.m1 <- besselK(u / phi, kappa - 1)
      bk.m2 <- besselK(u / phi, kappa - 2)
      out <- (2^(-kappa - 1) * phi^(-kappa - 4) * u^kappa * (bk.p2 * u^2 + 2 * bk * u^2 +
                                                               bk.m2 * u^2 - 4 * kappa * bk.p1 * phi * u - 4 *
                                                               bk.p1 * phi * u - 4 * kappa * bk.m1 * phi * u - 4 * bk.m1 * phi * u +
                                                               4 * kappa^2 * bk * phi^2 + 4 * kappa * bk * phi^2)) / (gamma(kappa))
    }
    out
  }
  n <- attr(u, "Size")
  hess.phi.mat <- matrix(NA, nrow = n, ncol = n)
  ind <- lower.tri(hess.phi.mat)
  hess.phi <- der2.phi(as.numeric(u), phi, kappa)
  hess.phi.mat[ind] <-  hess.phi
  hess.phi.mat <- t(hess.phi.mat)
  hess.phi.mat[ind] <-  hess.phi
  diag(hess.phi.mat) <- rep(der2.phi(0, phi, kappa), n)
  hess.phi.mat
}
##' @title Gaussian Process Model Specification
##' @description Specifies the terms, smoothness, and nugget effect for a Gaussian Process (GP) model.
##' @param ... Variable representing the spatial coordinates for the GP model. If left blank the
##' `geometry` column from the data is used automatically.
##' @param kappa The smoothness parameter \eqn{\kappa}. Default is `0.5`.
##' @param nugget The nugget effect, which represents the variance of the measurement error.
##' Default is `FALSE` in which case it is not estimated. If `TRUE` the value will be estimated or
##' a positive numeric value can be provided instead to fix the effect.
##' @details The function constructs a list that includes the specified terms (spatial coordinates or covariates),
##' the smoothness parameter \eqn{\kappa}, and the nugget effect. This list can be used as a specification for a Gaussian Process model.
##' A model has one spatial Gaussian process, so a formula can contain only one `gp()` term.
##' @return A list of class \code{RiskMap_gp_spec} containing the following elements:
##' \item{term}{A character vector of the specified terms.}
##' \item{kappa}{The smoothness parameter \eqn{\kappa}.}
##' \item{nugget}{The nugget effect.}
##' \item{dim}{The number of specified terms.}
##' \item{label}{A character string representing the full call for the GP model.}
##' @examples
##' # this function is not intended for use by itself
##' gp()
##' gp(nugget = TRUE)
##'
##' # use in a model formula
##' data(italy_sim)
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##' @export
gp <- function (..., kappa = 0.5, nugget = FALSE) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  term <- NULL

  check_range(kappa, min = 0, allow_equal = FALSE)

  if (!(isTRUE(nugget) || isFALSE(nugget) ||
        (is.numeric(nugget) && length(nugget) == 1L && is.finite(nugget) &&
         nugget > 0))) {
    supplied <- describe_value(nugget)
    abort_input("{.arg nugget} must be {.code TRUE}, {.code FALSE} or a positive number, not {supplied}.", arg = "nugget")
  }

  if (isFALSE(nugget)){
    nugget <- 0
  }

  if (d == 0) {
    term <- "sf"
  } else {
    if (d > 0) {
      for (i in 1:d) {
        term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      }
    }

    for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])),
                                   "term.labels")
  }
  full.call <- paste0("gp(", term[1])
  if (d > 1)
    for (i in 2:d) full.call <- paste(full.call, term[i], sep = ",")
  label <- gsub("sf", "", paste0(full.call, ")"))
  ret <- list(term = term, kappa = kappa, nugget = nugget, dim = d,
              label = label)
  class(ret) <- "RiskMap_gp_spec"
  ret
}
##' @title Random Effect Model Specification
##' @description Specifies the terms for a random effect model.
##' @param ... Variables representing the random effects in the model. Each variable adds its own independent random effect.
##' @details The function constructs a list that includes the specified terms for the random effects. This list can be used as a specification for a random effect model.
##' Several `re()` calls in a formula are combined, so `re(region) + re(province)` is the same model as `re(region, province)`.
##' @return A list of class \code{RiskMap_re_spec} containing the following elements:
##' \item{term}{A character vector of the specified terms.}
##' \item{dim}{The number of specified terms.}
##' \item{label}{A character string representing the full call for the random effect model.}
##' @note At least one variable must be provided as input.
##'
##' @examples
##' # this function is not intended for use by itself
##' re(province)
##'
##' # use in a model formula
##' data(italy_sim)
##' fit <- glgpm(
##'   formula = y ~ gp() + re(province),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' @export
re <- function (...) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  term <- NULL

  if (d == 0) {
    abort_input("{.fn re} needs at least one variable, as in {.code re(village)}.", arg = "...")
  } else {
    if (d > 0) {
      for (i in 1:d) {
        term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      }
    }
    for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])),
                                   "term.labels")
  }
  full.call <- paste0("re(", term[1])
  if (d > 1)
    for (i in 2:d) full.call <- paste0(full.call, ",", term[i])
  label <- gsub("sf", "", paste0(full.call, ")"))
  ret <- list(term = term, dim = d, label = label)
  class(ret) <- "RiskMap_re_spec"
  ret
}

##' @title Combine the re() terms of a formula
##' @description Each variable in `re()` adds its own independent random effect,
##' so several `re()` calls mean the same as one call listing all their variables.
##' @param re_terms A character vector of `re()` calls, as they appear in the formula.
##' @param formula_env The environment of the formula.
##' @return A list of class `RiskMap_re_spec`, or `NULL` if `re_terms` is empty.
##' @noRd
combine_re_terms <- function(re_terms, formula_env) {
  if (length(re_terms) == 0L) {
    return(NULL)
  }
  re_specs <- lapply(re_terms, function(re_term) {
    eval(parse(text = re_term), envir = formula_env)
  })
  if (length(re_specs) == 1L) {
    return(re_specs[[1L]])
  }
  re_variables <- unique(unlist(lapply(re_specs, function(re_spec) re_spec$term)))
  combined_call <- paste0("re(", paste(re_variables, collapse = ", "), ")")
  eval(parse(text = combined_call), envir = formula_env)
}

##' @title Read the offset column of a formula
##' @param formula_terms The output of `terms.formula()` for the formula.
##' @return The name of the offset column, or `NULL` if there is no `offset()` term.
##' @noRd
read_offset_column <- function(formula_terms) {
  offset_index <- attr(formula_terms, "offset")
  if (length(offset_index) == 0L) {
    return(NULL)
  }
  if (length(offset_index) > 1L) {
    abort_input("{.arg formula} can contain only one {.fn offset} term.", arg = "formula")
  }
  offset_call <- attr(formula_terms, "variables")[[offset_index + 1L]]
  offset_variable <- offset_call[[2L]]
  if (!is.symbol(offset_variable)) {
    offset_text <- deparse(offset_call)
    abort_input(c("{.fn offset} must name a column of {.arg data}, not {.code {offset_text}}.", "i" = "Add the values as a column and name it in {.fn offset}."), arg = "formula")
  }
  as.character(offset_variable)
}

interpret_formula <- function(formula) {
  p.env <- environment(formula)
  tf <- terms.formula(formula, specials = c("gp", "re"))
  terms <- attr(tf, "term.labels")
  nt <- length(terms)

  response <- as.character(attr(tf, "variables")[2])

  gp <- attr(tf, "specials")$gp
  re <- attr(tf, "specials")$re
  vtab <- attr(tf, "factors")

  if (length(gp) > 0) {
    for (i in 1:length(gp)) {
      ind <- (1:nt)[as.logical(vtab[gp[i], ])]
      gp[i] <- ind
    }
  }

  if (length(re) > 0) {
    for (i in 1:length(re)) {
      ind <- (1:nt)[as.logical(vtab[re[i], ])]
      re[i] <- ind
    }
  }

  len.gp <- length(gp)
  len.re <- length(re)
  if (len.gp > 1L) {
    abort_input("{.arg formula} can contain only one {.fn gp} term, because the model has one spatial Gaussian process.", arg = "formula")
  }
  gp_spec <- eval(parse(text = terms[gp]), envir = p.env)
  re_spec <- combine_re_terms(terms[re], p.env)
  offset <- read_offset_column(tf)

  if (length(terms[-c(gp, re)]) > 0) {
    pf <- paste(response, "~", paste(terms[-c(gp, re)], collapse = " + "))
  } else if (length(terms[-c(gp, re)]) == 0) {
    pf <- paste(response, "~ 1")
  }

  if (attr(tf, "intercept") == 0) {
    pf <- paste0(pf, "-1")
  }

  ret <- list(
    pf = as.formula(pf, p.env),
    gp_spec = gp_spec,
    re_spec = re_spec,
    offset = offset,
    response = response
  )
  ret
}

##' @title Extract terms from formula ignoring kappa and nugget
##' @description Recursively extract variable names from a formula/expression,
##' but for calls to gp(), only look inside unnamed (positional) arguments.
##' @param formula The formula to check
##' @return A character vector of terms
##' @noRd
get_formula_terms <- function(formula) {
  if (is.symbol(formula)) {
    return(as.character(formula))
  }

  if (is.call(formula)) {
    fn_name <- if (is.symbol(formula[[1]])) as.character(formula[[1]]) else ""

    args <- as.list(formula)[-1]
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- rep("", length(args))

    if (fn_name == "gp") {
      # only keep unnamed arguments
      args <- args[arg_names == ""]
    }

    return(unlist(lapply(args, get_formula_terms)))
  }

  NULL
}


##' @title Check that select columns of the data contain no missing values
##' @description Checks that the specified columns are complete, ignoring
##' missing values elsewhere in `data`. The error message names the argument
##' as passed by the caller (for example, `data` from `check_formula()`).
##' @param data The data to check.
##' @param columns The column names to check for missing data.
##' @param arg The argument the user gave the data as, by default the
##'   expression passed as `data`.
##' @return TRUE if there is no missing data, or raise an error if not.
##' @noRd
check_complete_data <- function(data, columns,
                                arg = deparse(substitute(data))) {
  name <- arg
  drop_coords <- st_drop_geometry(data[, columns])
  if (any(!complete.cases(drop_coords))) {
    abort_input(c("{.arg {name}} has rows with missing values in the variables of the model.", "i" = "Check or remove them."), arg = name)
  }
  invisible(TRUE)
}

##' @title Check that formula is valid and that there is no missing data
##' @description Checks that the formula object is of class formula and that all
##' the terms in the formula are present in the data
##' @param formula The formula to check
##' @param data The data to look for variables in
##' @param response_required Whether the response must be present in `data`.
##' @return TRUE if the formula is valid or raise an error if not
##' @noRd
check_formula <- function(formula, data, response_required = TRUE){

  if (!inherits(formula, "formula")) {
    abort_input("{.arg formula} must be a formula, not {.obj_type_friendly {formula}}.", arg = "formula")
  }

  if (length(formula) != 3L) {
    abort_input("{.arg formula} must have a response on the left-hand side, as in {.code y ~ gp()}.", arg = "formula")
  }

  formula_terms <- unique(get_formula_terms(formula))
  column_names <- names(data)

  contains_gp <- !is.null(attr(terms(formula, specials = "gp"), "specials")$gp)

  if (!contains_gp){
    abort_input("{.arg formula} must contain a Gaussian process term, {.fn gp}.", arg = "formula")
  }

  if (!response_required) {
    formula_terms <- setdiff(formula_terms, all.vars(formula[[2L]]))
  }
  missing_columns <- setdiff(formula_terms, column_names)
  n_missing <- length(missing_columns)
  if (n_missing > 0){

    abort_input("{.arg data} must have the column{?s} {.val {missing_columns}} named in {.arg formula}.", arg = "data")
  }

  check_complete_data(data, formula_terms)

  invisible(TRUE)
}

##' @title Extract Parameter Estimates from a "RiskMap" Model Fit
##' @description This \code{coef} method for the "RiskMap" class extracts the
##' maximum likelihood estimates from model fits obtained from the \code{\link{glgpm}} functions.
##' @param object An object of class "RiskMap" obtained as a result of a call to \code{\link{glgpm}}.
##' @param ... other parameters.
##' @return A list containing the maximum likelihood estimates. For standard models:
##' \item{beta}{A vector of coefficient estimates.}
##' \item{sigma2}{The estimate for the variance parameter \eqn{\sigma^2}.}
##' \item{phi}{The estimate for the scale parameter \eqn{\phi} of the spatial correlation function.}
##' \item{tau2}{The estimate for the nugget effect parameter \eqn{\tau^2}, if applicable.}
##' \item{sigma2_me}{The estimate for the measurement error variance \eqn{\sigma^2_{me}}, if applicable.}
##' \item{sigma2_re}{A named vector of variance estimates for the random effects, if applicable.}
##' \item{beta}{Coefficient estimates for log mean worm burden.}
##' \item{k}{Negative binomial overdispersion parameter.}
##' \item{rho}{Egg detection rate (fecundity).}
##' \item{sigma2}{Spatial process variance.}
##' \item{phi}{Spatial correlation scale.}
##' @seealso \code{\link{glgpm}}
##' @method coef RiskMap
##' @export
##'
coef.RiskMap <- function(object, ...) {

  estimate <- object$estimate

  beta_names <- colnames(as.matrix(object$D))

  res      <- list()
  res$beta <- estimate$beta
  names(res$beta) <- beta_names

  res$sigma2 <- exp(estimate$sigma2)
  res$phi    <- exp(estimate$phi)

  if (object$family == "gaussian" && !is.null(estimate$sigma2_me))
    res$sigma2_me <- exp(estimate$sigma2_me)

  if (!is.null(estimate$nu2))
    ## tau2 = nu2 * sigma2, so on the log scale their raw estimates add
    res$tau2 <- exp(estimate$nu2 + estimate$sigma2)

  if (!is.null(estimate$sigma2_re))
    res$sigma2_re <- exp(estimate$sigma2_re)

  return(res)
}

##' Extract fitted values from a RiskMap model
##'
##' @param object A fitted `RiskMap` object.
##' @param type Whether to return fitted values on the response or link scale.
##' @param ... Additional arguments, currently unused.
##' @return A numeric vector, or `NULL` when final conditional sampling was
##' skipped for a non-Gaussian model.
##' @method fitted RiskMap
##' @examples
##' data(italy_sim)
##'
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   return_samples = TRUE,
##'   messages = FALSE
##' )
##'
##' fitted(fit)
##'
##' @export
fitted.RiskMap <- function(object, type = c("response", "link"), ...) {
  type <- match_choice(type)
  if (type == "response") object$fitted_values else object$linear_predictors
}

##' @title Summarize Model Fits
##' @description Provides a \code{summary} method for the "RiskMap" class that
##'   computes standard errors and confidence intervals for likelihood-based
##'   model fits
##' @param object An object of class "RiskMap" from \code{\link{glgpm}}
##' @param ... other parameters.
##' @param level Confidence level for intervals (default 0.95).
##' @return A list of class \code{"summary.RiskMap"} with parameter estimates,
##'   standard errors, and confidence intervals. For binomial and Poisson
##'   models, its component \code{mcml} reports certified Monte Carlo maximum
##'   likelihood: \code{certificate} (an approximate bound, with probability
##'   \code{level}, on the Monte Carlo error of the estimate in statistical
##'   standard errors, \code{NA} when there is none), \code{tolerance},
##'   \code{certified}, \code{outcome}, \code{correction_refused} (whether
##'   the Newton correction of the last stage did not converge),
##'   \code{stages}, \code{fitting_draws} and \code{validation_draws} (at the
##'   last stage), \code{certification_bank} (the effective draws and Pareto
##'   shape of the certification draws) and \code{checks_to_review}, the number
##'   of checks of \code{\link{check_fit}} that need a look. \code{mcml} is
##'   \code{NULL} for Gaussian models.
##' @method summary RiskMap
##' @examples
##' data(italy_sim)
##'
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##'
##' summary(fit)
##' summary(fit, level = 0.98)
##' @export
summary.RiskMap <- function(object, ..., level = 0.95) {

  alpha  <- 1 - level
  z_crit <- qnorm(1 - alpha / 2)
  res    <- list()

  # ---------------------------------------------------------------------------
  # Helper: log-normal CI (for positive parameters)
  # ---------------------------------------------------------------------------
  lnCI <- function(est, se)
    c(Estimate      = est,
      "Lower limit" = exp(log(est) - z_crit * se / est),
      "Upper limit" = exp(log(est) + z_crit * se / est))

  # ===========================================================================
  # STANDARD RISKMAP MODELS (glgpm)
  # ===========================================================================

  link_name <- NULL
  inv_expr  <- NULL
  named_link <- if (identical(object$link_function$name, "custom")) {
    object$link_function$link
  }
  if (!is.null(named_link)) {
    link_name <- named_link
    inv_expr  <- paste("Inverse link function =", switch(
      named_link,
      probit = "pnorm(x)",
      cloglog = "1 - exp(-exp(x))",
      cauchit = "0.5 + atan(x) / pi"
    ))
  } else if (identical(object$link_function$name, "custom")) {
    link_name <- "custom"
    inv_expr  <- paste("Inverse link function =", describe_inverse_link(
      object$link_function$inv, object$call$invlink
    ))
  } else {
    if (identical(object$family, "poisson")) {
      link_name <- "canonical (log)"
      inv_expr  <- "Inverse link function = exp(x)"
    } else if (identical(object$family, "binomial")) {
      link_name <- "canonical (logit)"
      inv_expr  <- "Inverse link function = 1 / (1 + exp(-x))"
    } else if (identical(object$family, "gaussian")) {
      link_name <- "identity"
      inv_expr  <- "Inverse link function = x"
    }
  }

  n_re     <- length(object$re)
  re_names <- if (n_re > 0) names(object$re) else NULL

  beta_names <- colnames(as.matrix(object$D))

  ## `object$estimate` is a list (#92: keeps e.g. a "sigma2" covariate's name
  ## from colliding with the spatial variance parameter's own name). The
  ## delta-method covariance adjustment below needs a flat vector to line up
  ## against `covariance` (one joint matrix over every parameter); `unlist()`
  ## gives that, disambiguating the same way ("beta.sigma2" vs "sigma2").
  estimate <- unlist(object$estimate)
  nm       <- names(estimate)

  beta_flat_names <- paste0("beta.", beta_names)
  has_tau2      <- "nu2" %in% nm
  has_sigma2_me <- object$family == "gaussian" && "sigma2_me" %in% nm
  re_par_names  <- if (n_re > 0) paste0("sigma2_re.", re_names) else NULL

  ## tau2 = nu2 * sigma2, so on the log scale their raw estimates add. This is
  ## the one linear reparametrisation of the working-scale parameters that
  ## isn't just an exp(); its covariance is obtained via the delta method
  ## below (J), alongside the other, unchanged parameters.
  J <- diag(length(estimate))
  dimnames(J) <- list(nm, nm)
  if (has_tau2) {
    estimate["nu2"] <- estimate["nu2"] + estimate["sigma2"]
    J["nu2", "sigma2"] <- 1
  }

  ## The new parameters are J times the working-scale ones, so their
  ## covariance is J covariance J'. A parameter without a standard error (an
  ## NA row and column, for example a variance estimated at zero) is left
  ## out, with any new parameter that depends on it, so that the missing
  ## values do not spread to the other parameters.
  covariance <- object$covariance
  free <- !is.na(diag(covariance))
  computable <- rowSums(J[, !free, drop = FALSE] != 0) == 0
  covariance_new <- matrix(NA_real_, length(nm), length(nm),
                           dimnames = list(nm, nm))
  covariance_new[computable, computable] <-
    J[computable, free, drop = FALSE] %*%
    covariance[free, free, drop = FALSE] %*%
    t(J[computable, free, drop = FALSE])
  se_par <- sqrt(diag(covariance_new))

  non_beta <- setdiff(nm, beta_flat_names)
  estimate[non_beta] <- exp(estimate[non_beta])

  se_beta <- se_par[beta_flat_names]
  zval <- estimate[beta_flat_names] / se_beta
  res$reg_coef <- cbind(
    Estimate      = estimate[beta_flat_names],
    "Lower limit" = estimate[beta_flat_names] - se_beta * z_crit,
    "Upper limit" = estimate[beta_flat_names] + se_beta * z_crit,
    StdErr        = se_beta,
    z.value       = zval,
    p.value       = 2 * pnorm(-abs(zval))
  )
  rownames(res$reg_coef) <- beta_names

  if (object$family == "gaussian") {
    if (has_sigma2_me) {
      est_me <- estimate["sigma2_me"]
      se_me  <- se_par["sigma2_me"]
      res$me <- cbind(
        Estimate      = est_me,
        "Lower limit" = exp(log(est_me) - z_crit * se_me),
        "Upper limit" = exp(log(est_me) + z_crit * se_me)
      )
      rownames(res$me) <- unname(parameter_labels("sigma2_me",
                                                  object$distance_units))
    } else {
      res$me <- object$fix_var_me
    }
  }

  sp_names <- c("sigma2", "phi", if (has_tau2) "nu2")
  est_sp <- estimate[sp_names]
  se_sp  <- se_par[sp_names]
  res$sp <- cbind(
    Estimate      = est_sp,
    "Lower limit" = exp(log(est_sp) - z_crit * se_sp),
    "Upper limit" = exp(log(est_sp) + z_crit * se_sp)
  )
  # nu2 is reported as tau2, the variance of the nugget.
  rownames(res$sp) <- unname(parameter_labels(c("sigma2", "phi",
                                                if (has_tau2) "tau2"),
                                              object$distance_units))
  if (!is.null(object$fix_tau2)) res$tau2 <- object$fix_tau2

  if (n_re > 0) {
    est_re <- estimate[re_par_names]
    se_re  <- se_par[re_par_names]
    res$ranef <- cbind(
      Estimate      = est_re,
      "Lower limit" = exp(log(est_re) - z_crit * se_re),
      "Upper limit" = exp(log(est_re) + z_crit * se_re)
    )
    rownames(res$ranef) <- unname(parameter_labels(re_par_names,
                                                   object$distance_units))
  }

  res$level           <- level
  res$family          <- object$family
  res$kappa           <- object$kappa
  res$log_lik         <- object$log_lik
  res$cov_offset_used <- !(is.null(object$cov_offset) ||
                             all(object$cov_offset == 0))
  if (object$family == "gaussian") {
    res$aic <- 2 * length(unlist(object$estimate)) - 2 * res$log_lik
  }

  res$call               <- object$call %||% NULL
  res$link_name          <- link_name
  res$invlink_expression <- inv_expr
  res$mcml               <- summarise_mcml_fit(object)

  class(res) <- "summary.RiskMap"
  return(res)
}

##' How `summary()` writes an inverse link given in `invlink`
##'
##' The body of a function written by the user, such as `1/(1 + exp(-x))`.
##' The body of a package function such as `pnorm`, or of a primitive such
##' as `exp`, is not informative (`.Call(C_pnorm, ...)` or `NULL`), so such a
##' function is written with the name it was given in the call, as
##' `pnorm(x)`, when the call names it.
##' @param inverse The inverse link stored on the fit.
##' @param supplied The `invlink` argument of the fit's call.
##' @noRd
describe_inverse_link <- function(inverse, supplied) {
  from_package <- is.primitive(inverse) || isNamespace(environment(inverse))
  if (!from_package) {
    return(tryCatch(
      paste(deparse(body(inverse), width.cutoff = 500L), collapse = " "),
      error = function(e) "<user-supplied function>"
    ))
  }
  named <- is.name(supplied) ||
    (is.call(supplied) && identical(supplied[[1L]], as.name("::")))
  if (named) paste0(deparse(supplied), "(x)") else "<user-supplied function>"
}

##' Readable names of the parameters, as `summary()` labels them
##'
##' `names` are those of `unlist(fit$estimate)` or `unlist(coef(fit))`, with
##' `tau2` for the variance of the nugget. The prefix `beta.` of the
##' regression coefficients is dropped, and any other name is kept as it is.
##' The labels of `summary()`, `profile()` and `parametric_bootstrap()` all
##' come from here.
##' @noRd
parameter_labels <- function(names, distance_units) {
  labels <- sub("^beta\\.", "", names)
  labels[names == "sigma2"] <- "Spatial process var."
  labels[names == "phi"] <- paste0("Spatial corr. scale (",
                                   distance_units %||% "units", ")")
  labels[names == "tau2"] <- "Variance of the nugget"
  labels[names == "sigma2_me"] <- "Measurement error var."
  random <- grepl("^sigma2_re\\.", names)
  labels[random] <- paste0(sub("^sigma2_re\\.", "", names[random]),
                           " (random eff. var.)")
  stats::setNames(labels, names)
}

##' @title Print Summary of RiskMap Model
##' @description Print method for objects of class \code{"summary.RiskMap"}.
##' @param x An object of class \code{"summary.RiskMap"}.
##' @param ... other parameters.
##' @return Invisibly returns \code{x}.
##' @method print summary.RiskMap
##' @export
print.summary.RiskMap <- function(x, ...) {
  out_print(format_fit_summary(x))
  invisible(x)
}

##' The lines that print the summary of a fitted model
##' @noRd
format_fit_summary <- function(x) {
  traits <- fit_traits(x)
  title <- model_title(traits)
  coefficients <- function(table, p_values) {
    out_capture(printCoefmat(table, P.values = p_values,
                             has.Pvalue = p_values))
  }
  nugget_estimated <- "Variance of the nugget" %in% rownames(x$sp)
  monte_carlo <- if (!is.null(x$mcml)) {
    c("",
      out_heading("Monte Carlo maximum likelihood"),
      out_bullets(c(stats::setNames(describe_mcml_certificate(x$mcml),
                                    certificate_symbol(x$mcml)),
                    describe_checks_to_review(x$mcml$checks_to_review))))
  }
  c(out_rule(title),
    format_call(x$call),
    out_text(c(describe_summary_link(x$link_name, x$invlink_expression),
               paste0("The limits are those of ", x$level * 100,
                      "% confidence intervals."))),
    "",
    out_heading("Regression coefficients"),
    coefficients(x$reg_coef, TRUE),
    if (isTRUE(x$cov_offset_used)) {
      out_text("Offset included in the linear predictor")
    },
    if (isTRUE(traits$measurement_error)) {
      if (length(x$me) > 1) {
        c("", out_heading("Measurement error"), coefficients(x$me, FALSE))
      } else {
        c("", out_text(paste0("Measurement error var. fixed at ", x$me)))
      }
    },
    "",
    out_heading(paste0("Spatial Gaussian process, Matern correlation with ",
                       "kappa = ", x$kappa)),
    coefficients(x$sp, FALSE),
    if (!nugget_estimated && is.numeric(x$tau2)) {
      out_text(paste0("Variance of the nugget effect fixed at ", x$tau2))
    },
    if (!is.null(x$ranef)) {
      c("", out_heading("Unstructured random effects"),
        coefficients(x$ranef, FALSE))
    },
    "",
    out_pairs(c("Log-likelihood", if (!is.null(x$aic)) "AIC"),
              format(c(x$log_lik, x$aic))),
    monte_carlo)
}

##' One sentence on the link of a summarised fit
##'
##' The canonical links and the identity are named, and any other link is
##' given with its inverse, as `summary()` records it.
##' @noRd
describe_summary_link <- function(link_name, invlink_expression) {
  if (is.null(link_name)) return(NULL)
  canonical <- regmatches(link_name, regexec("^canonical \\((.*)\\)$",
                                             link_name))[[1L]]
  if (length(canonical)) {
    return(paste0("The canonical ", canonical[2L], " link."))
  }
  if (identical(link_name, "identity")) return("The identity link.")
  inverse <- sub("^Inverse link function = ", "", invlink_expression %||% "")
  paste0(if (identical(link_name, "custom")) "A custom" else
           paste("The", link_name),
         " link", if (nzchar(inverse)) paste0(", with inverse ", inverse), ".")
}

##' @title Format RiskMap Results as a Table
##' @description Converts a fitted "RiskMap" model, cross-validation
##' results or polygon predictions into a table that renders directly in
##' Quarto, R Markdown, HTML, LaTeX and the R console.
##' @param object An object of class "RiskMap" resulting from a call to
##' \code{\link{glgpm}}, a "summary.RiskMap" object, a
##' "summary.RiskMap_cross_validation" object, or a "RiskMap_prediction"
##' object (see \code{\link{to_table.RiskMap_prediction}}).
##' @param digits A non-negative integer giving the number of decimal places
##' used to display numeric results.
##' @param col_names Optional column names, replacing the default headers.
##' @param ... Additional arguments passed to \code{\link[knitr]{kable}}.
##' @details This function creates a presentation-ready summary table from a
##' fitted "RiskMap" model or cross-validation results for multiple models.
##' Use \code{\link{coef}} or \code{\link{summary}} when numeric results are
##' required for further analysis. Numeric values use fixed notation with
##' \code{digits} decimal places, except when scientific notation is needed to
##' represent very large or very small values clearly. Values that are
##' compared share one notation: the estimate and confidence limits of a
##' parameter, and one metric across models.
##'
##' When the input is a "RiskMap" model object, the table includes:
##' \itemize{
##'   \item Regression coefficients with their estimates, confidence intervals, and p-values.
##'   \item Parameters for the spatial process.
##'   \item Random effect variances.
##'   \item Measurement error variance, if applicable.
##' }
##'
##' When the input is a cross-validation summary object ("summary.RiskMap_cross_validation"), the table includes:
##' \itemize{
##'   \item A row for each model being compared.
##'   \item Performance metrics such as CRPS and SCRPS for each model.
##' }
##'
##' @return An object of class "knitr_kable" that can be rendered directly.
##' @importFrom knitr kable
##' @export
##' @seealso \code{\link{glgpm}}, \code{\link{summary.RiskMap_cross_validation}}
##' @examples
##' data(italy_sim)
##'
##' fit <- glgpm(
##'   formula = y ~ gp(),
##'   data = italy_sim[1:100,],
##'   family = "gaussian",
##'   messages = FALSE
##' )
##' to_table(fit)
##' to_table(fit, digits = 1)
to_table <- function(object, ...) {
  check_supplied(object)
  UseMethod("to_table")
}

##' @rdname to_table
##' @export
to_table.default <- function(object, digits = 3, col_names = NULL, ...) {
  check_positive_integer(digits, "digits", allow_zero = TRUE)

  if (inherits(object, "summary.RiskMap") ||
      inherits(object, "summary.RiskMap_cross_validation")) {
    summary_out <- object
  } else {
    summary_out <- summary(object)
  }
  if (inherits(summary_out,
               what = "summary.RiskMap", which = FALSE)) {
    tab <- rbind(summary_out$reg_coef[, 1:3], summary_out$sp, summary_out$ranef,
                 summary_out$me)
    include_row_names <- TRUE
  } else if (inherits(summary_out,
                      what = "summary.RiskMap_cross_validation", which = FALSE)) {
    n_models <- nrow(summary_out)
    n_metrics <- ncol(summary_out)
    model_names <- rownames(summary_out)
    metric_names <- toupper(colnames(summary_out))
    tab <- data.frame(Model = model_names)
    for (i in seq_len(n_metrics)) {
      tab[[paste(metric_names[i])]] <- summary_out[,i]
    }
    include_row_names <- FALSE
  } else {
    abort_input("{.arg object} must be a {.cls RiskMap} fit or a {.cls RiskMap_cross_validation} result, not {.obj_type_friendly {object}}.", arg = "object")
  }

  tab <- as.data.frame(tab, check.names = FALSE)
  numeric_columns <- which(vapply(tab, is.numeric, logical(1)))
  values <- as.matrix(tab[numeric_columns])
  needs_scientific <- is.finite(values) & values != 0 &
    (abs(values) >= 1e6 | abs(values) < 10^(-digits))
  # Values that are compared share one notation: the estimate and limits of a
  # parameter (a row of a model table), or one metric across models (a
  # column of a cross-validation table).
  scientific <- if (include_row_names) {
    matrix(rowSums(needs_scientific) > 0, nrow(values), ncol(values))
  } else {
    matrix(colSums(needs_scientific) > 0, nrow(values), ncol(values),
           byrow = TRUE)
  }
  for (j in seq_along(numeric_columns)) {
    x <- values[, j]
    out <- formatC(x, format = "f", digits = as.integer(digits))
    use_scientific <- scientific[, j] & is.finite(x)
    out[use_scientific] <- formatC(
      x[use_scientific],
      format = "e",
      digits = as.integer(digits)
    )
    tab[[numeric_columns[j]]] <- out
  }

  dots <- list(...)
  if (!is.null(col_names)) {
    if (length(col_names) != ncol(tab)) {
      abort_input("{.arg col_names} must have {ncol(tab)} elements, one for each column, not {length(col_names)}.", arg = "col_names")
    }
    dots$col.names <- col_names
  }
  if (is.null(dots$row.names))
    dots$row.names <- include_row_names
  if (is.null(dots$align))
    dots$align <- if (include_row_names) rep("r", ncol(tab)) else
      c("l", rep("r", ncol(tab) - 1L))

  do.call(kable, c(list(x = tab), dots))
}

##' @title Compute Unique Coordinate Identifiers
##'
##' @description
##' This function identifies unique coordinates from a `sf` (simple feature) object
##' and assigns an identifier to each coordinate occurrence. It returns a list
##' containing the identifiers for each row and a vector of unique identifiers.
##'
##' @param data An `sf` object containing geometrical data from which coordinates are extracted.
##'
##' @return A list with the following elements:
##' \describe{
##'   \item{ID_coords}{An integer vector where each element corresponds to a row in the input,
##'   indicating the index of the unique coordinate in the full set of unique coordinates.}
##'   \item{s_unique}{An integer vector containing the unique identifiers of all distinct coordinates.}
##' }
##'
##' @details
##' The function extracts the coordinate pairs from the `sf` object and determines the unique
##' coordinates. It then assigns each row in the input data an identifier corresponding
##' to the unique coordinate it matches.
##'
##' @examples
##' data(italy_sim)
##'
##' create_ids(italy_sim)
##'
##' @export
create_ids <- function(data) {
  check_supplied(data)
  check_data(data)
  coords_o <- st_coordinates(data)
  coords <- unique(coords_o)

  m <- nrow(coords_o)
  ID_coords <- sapply(1:m, function(i)
    which(coords_o[i,1] == coords[,1] &
            coords_o[i,2] == coords[,2]))
  out <- list()
  out$ID_coords <- ID_coords
  out$s_unique <- unique(ID_coords)
  return(out)
}


##' @title Summarize Cross-Validation Scores for Spatial RiskMap Models
##'
##' @description This function summarizes cross-validation scores for different spatial models obtained
##' from \code{\link{assess_prediction}}.
##'
##' @param object A `RiskMap_cross_validation` object containing cross-validation scores for each
##'               model, as obtained from \code{\link{assess_prediction}}.
##' @param view_all Logical. If `TRUE`, stores the average scores across test sets for each
##'                 model alongside the overall average across all models. Defaults to `TRUE`.
##' @param ... Additional arguments passed to or from other methods.
##'
##' @details
##' The function computes and returns a matrix where rows correspond to models and columns
##' correspond to performance metrics (e.g., CRPS, SCRPS). Scores are weighted by subset sizes
##' to compute averages. Attributes of the returned object include:
##' \itemize{
##'   \item `test_set_means`: A list of average scores for each test set and model.
##'   \item `overall_averages`: Overall averages for each metric across all models.
##'   \item `view_all`: Indicates whether averages across test sets are available for visualization.
##' }
##'
##' @return A matrix of summary scores with models as rows and metrics as columns, with class
##' `"summary.RiskMap_cross_validation"`.
##'
##' @seealso \code{\link{assess_prediction}}
##'
##' @export
##' @method summary RiskMap_cross_validation
summary.RiskMap_cross_validation <- function(object, view_all = TRUE, ...) {
  model_names <- names(object$model)
  n_models <- length(model_names)

  metric_names <- names(object$model[[1]]$metric)
  if (is.null(metric_names)) {
    abort_input("{.arg object} has no metrics of predictive performance.", arg = "object")
  }
  n_metrics <- length(metric_names)

  res <- matrix(NA, ncol = n_metrics, nrow = n_models)
  colnames(res) <- metric_names
  rownames(res) <- model_names

  test_set_means <- list()

  n_subs <- length(object$model[[1]]$metric[[1]])
  # Weight each test set by its number of observations, for every metric
  per_observation <- intersect(c("CRPS", "SCRPS"), metric_names)
  w <- if (length(per_observation) > 0) {
    lengths(object$model[[1]]$metric[[per_observation[1]]])
  } else {
    vapply(object$test_set, nrow, integer(1))
  }

  for (i in 1:n_models) {
    model_scores <- list()
    for (j in 1:n_metrics) {
      score_j <- rep(NA, n_subs)
      for (h in 1:n_subs) {
        score_j[h] <- mean(object$model[[i]]$metric[[j]][[h]], na.rm = TRUE)
      }
      model_scores[[j]] <- score_j
      # Test sets whose assessment failed have no scores
      ok <- is.finite(score_j)
      res[i, j] <- if (any(ok)) sum(w[ok] * score_j[ok]) / sum(w[ok]) else NA
    }
    test_set_means[[model_names[i]]] <- model_scores
  }

  overall_averages <- colMeans(res, na.rm = TRUE)

  # Attach additional attributes for printing
  attr(res, "test_set_means") <- test_set_means
  attr(res, "overall_averages") <- overall_averages
  attr(res, "view_all") <- view_all
  attr(res, "keep_par_fixed") <- object$keep_par_fixed
  attr(res, "n_failed") <- if (is.null(object$failures)) 0L else nrow(object$failures)
  class(res) <- "summary.RiskMap_cross_validation"
  return(res)
}

##' @title Print Summary of RiskMap Spatial Cross-Validation Scores
##'
##' @description This function prints the matrix of cross-validation scores produced by
##' `summary.RiskMap_cross_validation` in a readable format.
##'
##' @param x An object of class `"summary.RiskMap_cross_validation"`, typically the output of
##'          `summary.RiskMap_cross_validation`.
##' @param ... Additional arguments passed to or from other methods.
##'
##' @details
##' This method is primarily used to format and display the summary score matrix,
##' printing it to the console. It provides a clear view of the cross-validation performance
##' metrics across different spatial models.
##'
##' @return Invisibly returns \code{x}.
##' @export
##' @method print summary.RiskMap_cross_validation
print.summary.RiskMap_cross_validation <- function(x, ...) {
  out_print(format_cv_summary(x))
  invisible(x)
}

##' The lines that print the summary of a cross-validation
##' @noRd
format_cv_summary <- function(x) {
  c(out_rule("Cross-validation scores"),
    format_cv_settings(attr(x, "keep_par_fixed"), attr(x, "n_failed")),
    "",
    format_cv_averages(x),
    if (isTRUE(attr(x, "view_all"))) c("", format_cv_test_sets(x)))
}

##' How the parameters were set and how many assessments failed
##' @noRd
format_cv_settings <- function(keep_par_fixed, n_failed, before = NULL) {
  parameters <- if (!is.null(keep_par_fixed)) {
    if (keep_par_fixed) {
      "Parameters are fixed at the estimates from all the data."
    } else {
      "Parameters are re-estimated in each training set."
    }
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

##' The scores of each model averaged over the test sets
##' @noRd
format_cv_averages <- function(x) {
  scores <- data.frame(Model = rownames(x), check.names = FALSE)
  for (metric in colnames(x)) {
    scores[[metric]] <- format_cv_score(x[, metric])
  }
  c(out_heading("Averages over the test sets, weighted by their sizes"),
    out_capture(print(scores, row.names = FALSE, right = TRUE)))
}

##' Scores to four decimals, NA where a failed assessment left none
##'
##' The mean of no scores is `NaN`, which is printed as NA, as the missing
##' scores it stands for.
##' @noRd
format_cv_score <- function(values) {
  ifelse(is.finite(values), sprintf("%.4f", values), "NA")
}

##' The scores of each model in each test set
##' @noRd
format_cv_test_sets <- function(x) {
  means <- attr(x, "test_set_means")
  rows <- lapply(names(means), function(model) {
    scores <- means[[model]]
    table <- data.frame(Model = model, "Test set" = seq_along(scores[[1L]]),
                        check.names = FALSE)
    for (k in seq_along(scores)) {
      table[[colnames(x)[k]]] <- format_cv_score(scores[[k]])
    }
    table
  })
  c(out_heading("By test set"),
    out_capture(print(do.call(rbind, rows), row.names = FALSE, right = TRUE)))
}

##' @title Print the Results of a Spatial Cross-Validation
##' @description Print method for objects of class
##'   \code{"RiskMap_cross_validation"} from \code{\link{assess_prediction}}.
##'   It gives the models, the test sets, how the parameters were set and the
##'   scores of each model averaged over the test sets.
##'   \code{\link{summary.RiskMap_cross_validation}} gives the scores of each
##'   test set, and \code{\link{plot.RiskMap_cross_validation}} maps the scores
##'   at the test locations and draws the calibration curve.
##' @param x An object of class \code{"RiskMap_cross_validation"}.
##' @param ... Not used.
##' @return Invisibly returns \code{x}.
##' @method print RiskMap_cross_validation
##' @export
print.RiskMap_cross_validation <- function(x, ...) {
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
  scores <- summary(x, view_all = FALSE)
  n_models <- length(x$model)
  sizes <- vapply(x$test_set, NROW, integer(1))
  setting <- paste0(n_models, if (n_models == 1L) " model" else " models",
                    " assessed on ", describe_test_sets(sizes), ".")
  hint <- paste0(out_fun("summary", "summary.RiskMap_cross_validation"),
                 " gives the scores of each test set, and ",
                 out_fun("plot", "plot.RiskMap_cross_validation"),
                 " maps the scores at the test locations and draws the ",
                 "calibration curve.")
  c(out_rule("Cross-validation"),
    format_cv_settings(x$keep_par_fixed, NROW(x$failures), before = setting),
    "",
    format_cv_averages(scores),
    "",
    out_bullets(c(i = hint)))
}

##' Name a predictive-performance metric as it should appear on a plot
##'
##' The calibration metrics are stored as `"PIT"` and `"PIT_area"` for every
##' family, but are labelled AnPIT for discrete families.
##'
##' @param object A `RiskMap_cross_validation` object.
##' @param metric Character string naming the metric.
##' @param model Character string naming the model.
##' @return Character string.
##' @noRd
.metric_label <- function(object, metric, model) {
  pit_label <- if (is.null(object$model[[model]]$AnPIT)) "PIT" else "AnPIT"
  switch(metric,
         PIT = pit_label,
         PIT_area = paste0(pit_label, "_area"),
         metric)
}

##' Build one spatial map of a predictive-performance metric
##'
##' @noRd
.plot_metric_map <- function(object, metric, model, ...) {

  if (!model %in% names(object$model)) {
    abort_input("{.arg model} must name a model of {.arg object}, and {.val {model}} doesn't.", arg = "model")
  }

  if (!metric %in% names(object$model[[model]]$metric)) {
    abort_input("{.arg metric} must be a metric computed for model {.val {model}}, and {.val {metric}} isn't.", arg = "metric")
  }

  # Extract the test sets and number of test sets
  test_sets <- object$test_set
  n_test <- length(test_sets)

  # Combine the data and add the metric variable
  data_full <- st_as_sf(test_sets[[1]])
  data_full$value <- object$model[[model]]$metric[[metric]][[1]]

  if (n_test > 1) {
    for (i in 2:n_test) {
      test_sets[[i]]$value <- object$model[[model]]$metric[[metric]][[i]]
      data_full <- rbind(data_full, test_sets[[i]])
    }
  }

  # Check for duplicate locations and average the metric
  data_full <- data_full %>%
    mutate(geom_id = st_as_text(.data$geometry)) %>%
    group_by(.data$geom_id) %>%
    summarize(value = mean(.data$value, na.rm = TRUE),
              geometry = first(.data$geometry), .groups = "drop") %>%
    st_as_sf()

  # Create the base plot
  metric_label <- .metric_label(object, metric, model)
  out <- ggplot(data = data_full) +
    geom_sf(aes(color = .data$value), size = 2) +
    labs(title = paste0("Model ", model, ": ", metric_label),
         colour = metric_label) +
    theme_minimal()

  # Layer on any additional ggplot components passed via ...
  Reduce(`+`, list(...), out)
}

##' Build calibration-curve plots (AnPIT / PIT) for one or more models
##'
##' For Binomial or Poisson models this visualises the Aggregated
##' normalised Probability Integral Transform (AnPIT) curves stored in
##' `$AnPIT`; for Gaussian models it instead plots the empirical PIT curve
##' (ECDF of `$PIT` values) on the same grid. A 45-degree dashed red line
##' indicates perfect calibration.
##'
##' @return A named list of ggplot objects, one per `model` - or, when
##'   `pit_mode == "average"` and `combine_pit` is `TRUE`, a single-element
##'   list holding one combined plot.
##' @noRd
.plot_calibration_curve <- function(object, pit_mode, pit_test_set, model,
                                    combine_pit) {

  make_df <- function(mname) {
    m <- object$model[[mname]]
    if (!is.null(m$AnPIT)) {
      lapply(seq_along(m$AnPIT), function(j) {
        curve_vals <- m$AnPIT[[j]]
        if (length(curve_vals) == 0) return(NULL)
        data.frame(
          u_val   = seq(0, 1, length.out = length(curve_vals)),
          value   = curve_vals,
          test_set = j,
          model    = mname,
          type     = "AnPIT"
        )
      })
    } else if (!is.null(m$PIT)) {
      u_grid <- seq(0, 1, length.out = 1000)
      lapply(seq_along(m$PIT), function(j) {
        pit_vec <- m$PIT[[j]]
        if (length(pit_vec) == 0) return(NULL)
        data.frame(
          u_val   = u_grid,
          value   = ecdf(pit_vec)(u_grid),
          test_set = j,
          model    = mname,
          type     = "PIT"
        )
      })
    } else {
      NULL
    }
  }

  plot_data <- do.call(rbind, unlist(lapply(model, make_df), recursive = FALSE))

  if (is.null(plot_data) || nrow(plot_data) == 0)
    abort_input("{.arg x} has no PIT or AnPIT values to plot.", arg = "x")

  y_label <- unique(plot_data$type)
  if (length(y_label) > 1) y_label <- "Calibration curve"

  id_line <- geom_abline(intercept = 0, slope = 1,
                         linetype = "dashed", colour = "red")

  if (pit_mode == "average" && combine_pit) {
    avg <- plot_data %>%
      dplyr::group_by(.data$model, .data$u_val) %>%
      dplyr::summarize(value = mean(.data$value), .groups = "drop")

    return(list(
      ggplot(avg, aes(.data$u_val, .data$value, colour = .data$model)) +
        geom_line() + id_line +
        labs(title = "Average calibration curves",
             x = "", y = y_label) +
        theme_minimal() +
        guides(colour = guide_legend(title = "Model"))
    ))
  }

  build_plot <- function(df, title_suffix = "") {
    ggplot(df, aes(.data$u_val, .data$value,
                   colour = if (pit_mode == "all") as.factor(.data$test_set) else NULL)) +
      geom_line() + id_line +
      labs(title = title_suffix, x = "", y = unique(df$type)) +
      theme_minimal() +
      guides(colour = guide_legend(title = "Test set"))
  }

  plots <- list()
  for (mname in model) {
    df_model <- dplyr::filter(plot_data, .data$model == mname)

    title <- paste0("Model ", mname, ": ", unique(df_model$type))
    p <- switch(pit_mode,
                average = {
                  avg <- df_model %>%
                    dplyr::group_by(.data$u_val) %>%
                    dplyr::summarize(value = mean(.data$value), .groups = "drop")
                  avg$type <- unique(df_model$type)
                  build_plot(avg, paste(title, "(average)"))
                },
                single  = {
                  if (is.null(pit_test_set)) {
                    abort_input("{.arg pit_test_set} must be given when {.arg pit_mode} is {.val single}.", arg = "pit_test_set")
                  }
                  df_ts <- df_model[df_model$test_set == pit_test_set,]
                  if (nrow(df_ts) == 0) {
                    abort_input("{.arg pit_test_set} must be a test set of model {.val {mname}}, and {pit_test_set} isn't.", arg = "pit_test_set")
                  }
                  build_plot(df_ts,
                             paste0(title, " (test set ", pit_test_set, ")"))
                },
                all     = build_plot(df_model,
                                     paste(title, "(all test sets)")),
                abort_input("{.arg pit_mode} must be {.val average}, {.val single} or {.val all}.", arg = "pit_mode")
    )

    plots[[mname]] <- p
  }

  plots
}

##' Arrange one metric group's plots into a single grid, if there is more
##' than one
##' @noRd
.arrange_metric_group <- function(group_plots) {
  if (length(group_plots) == 1) {
    return(group_plots[[1]])
  }

  ncol <- 2
  do.call(
    gridExtra::grid.arrange,
    c(group_plots, list(ncol = ncol, nrow = ceiling(length(group_plots) / ncol)))
  )
}

##' @title Plot Method for RiskMap_cross_validation Objects
##'
##' @description
##' Plots whatever predictive-performance output is available in a
##' \code{RiskMap_cross_validation} object returned by
##' \code{\link{assess_prediction}}: a spatial map for each per-location
##' metric (\code{"CRPS"}, \code{"SCRPS"}, \code{"PIT_area"}) that was
##' requested via \code{metrics} in that call, and/or a calibration curve
##' (\code{"PIT"}) - the Aggregated normalised Probability Integral
##' Transform (AnPIT) curve for discrete families, or the empirical PIT curve
##' for Gaussian models - when that was requested. Plot titles and legends
##' use AnPIT or PIT according to each model's family. By default every available
##' plot is produced, one per requested model, so you don't need to know in
##' advance which metrics were computed. Each metric gets its own plot
##' (or, with more than one model, its own small grid of one panel per
##' model) rather than everything being squeezed into one combined grid.
##'
##' @param x A \code{RiskMap_cross_validation} object.
##' @param metric Character vector restricting which plot(s) to produce;
##'   one or more of \code{"CRPS"}, \code{"SCRPS"}, \code{"PIT_area"} (a
##'   spatial map of that score) and \code{"PIT"} (the calibration
##'   curve). Defaults to every metric available in \code{x}.
##' @param model Character vector of model names to include. Defaults to
##'   every model in \code{x$model}.
##' @param ... Additional \pkg{ggplot2} components (e.g.
##'   \code{scale_color_gradient()}), layered onto every spatial-map plot.
##'   Must be passed by position only after `model`, since the remaining
##'   arguments must be named.
##' @param pit_mode For the calibration curve only: one of \code{"average"}
##'   (average curve across test sets, the default), \code{"single"} (one
##'   specific test set) or \code{"all"} (every test set separately).
##' @param pit_test_set Integer; required when \code{pit_mode = "single"}.
##' @param combine_pit Logical; when \code{pit_mode = "average"}, draw the
##'   calibration curves for every model in a single panel (\code{TRUE})
##'   rather than one panel per model (\code{FALSE}, default).
##'
##' @return When exactly one metric is produced: a single \pkg{ggplot2}
##'   object (one model) or a \pkg{grid} object from \pkg{gridExtra} (more
##'   than one model). When more than one metric is produced, each metric's
##'   plot/grid is drawn in turn as a side effect (so each one renders as
##'   its own figure), and a named list of them - one per metric - is
##'   returned invisibly for later reuse.
##'
##' @seealso \code{\link{assess_prediction}}
##' @method plot RiskMap_cross_validation
##' @importFrom dplyr group_by summarize %>%
##' @export
plot.RiskMap_cross_validation <- function(x, metric = NULL, model = NULL, ...,
                                          pit_mode = "average", pit_test_set = NULL,
                                          combine_pit = FALSE) {

  all_models <- names(x$model)
  if (is.null(model)) {
    model <- all_models
  } else {
    missing_model <- setdiff(model, all_models)
    if (length(missing_model) > 0) {
      abort_input(c("{.arg model} must name models of {.arg x}.", "x" = "{.val {missing_model}} {cli::qty(length(missing_model))}{?is/are} not among them."), arg = "model")
    }
  }

  spatial_metrics <- intersect(
    c("CRPS", "SCRPS", "PIT_area"),
    names(x$model[[model[1]]]$metric)
  )
  has_calibration <- !is.null(x$model[[model[1]]]$PIT) ||
    !is.null(x$model[[model[1]]]$AnPIT)
  available_metrics <- c(if (has_calibration) "PIT", spatial_metrics)

  if (is.null(metric)) {
    metric <- available_metrics
  } else {
    unavailable <- setdiff(metric, available_metrics)
    if (length(unavailable) > 0) {
      abort_input(c("{.arg metric} must name metrics computed for {.arg x}.", "x" = "{.val {unavailable}} {cli::qty(length(unavailable))}{?was/were} not computed.", "i" = "The metrics computed are {.val {available_metrics}}."), arg = "metric")
    }
  }

  # Calibration curve first, then the spatial-map metrics in a fixed order,
  # so metrics are plotted in a stable, predictable order regardless of how
  # `metric`/the object's own field ordering happen to be arranged.
  ordered_metrics <- intersect(c("PIT", "CRPS", "SCRPS", "PIT_area"), metric)

  groups <- stats::setNames(lapply(ordered_metrics, function(m) {
    group_plots <- if (identical(m, "PIT")) {
      .plot_calibration_curve(x, pit_mode = pit_mode, pit_test_set = pit_test_set,
                              model = model, combine_pit = combine_pit)
    } else {
      stats::setNames(
        lapply(model, function(mod) .plot_metric_map(x, metric = m, model = mod, ...)),
        model
      )
    }
    .arrange_metric_group(group_plots)
  }), ordered_metrics)

  if (length(groups) == 1) {
    if (inherits(groups[[1]], "ggplot")) {
      return(groups[[1]])
    }
    return(invisible(groups[[1]]))
  }

  # draw each metric as its own figure rather than
  # squeezing every metric into a single combined grid.
  for (plot_obj in groups) {
    if (inherits(plot_obj, "ggplot")) print(plot_obj)
  }
  invisible(groups)
}

##' @title Check for valid binomial values
##'
##' @description
##' Checks that binomial data only consists of zero or positive integers and that if
##' den is provided that all values of y are less than or equal to it.
##' Some tolerance is provided for floating point errors
##'
##' @param y the data to check
##' @param den the denominator
##' @return TRUE if valid, raise an error if not
##' @noRd
check_binomial <- function(y, den){
  tolerance <- sqrt(.Machine$double.eps)
  valid <- all(y >= 0) & all(abs(y - round(y)) < tolerance)
  if (!valid) {
    abort_input("The responses of a binomial model must be non-negative whole numbers.", arg = "data")
  }

  if (!is.null(den) && !all(den >= y)) {
    abort_input("Binomial responses can't exceed their denominators.", arg = "denominator")
  }

  invisible(TRUE)
}

#' @title check_data
#' @description
#'
#' Check that the data is an sf or sfc object, with a CRS, only containing points
#' or either polygons or multipolygons.
#' If CRS == 4326 it also checks that the coordinates are possible (i.e. not
#' latitudes > 90)
#' @param data the data to check
#' @param geometry whether to check that the data contains `"point"` (default),
#' `"polygon"` (covering both polygons and multipolygons) or `"any"` to not check
#' @param type whether to check that the data is `"sf"` (default) or
#' `"sfc"` (either sf or sfc)
#' @param arg The argument the user gave the data as, by default the
#' expression passed as `data`. Internal callers name the user's argument.
#' @return TRUE if the data is valid. Raise an error if not.
#' @noRd
#'
check_data <- function(data, geometry = "point", type = "sf",
                       arg = deparse(substitute(data))) {
  if (!geometry %in% c("point", "polygon", "any")) {
    riskmap_bug("{.arg geometry} must be {.val point}, {.val polygon} or {.val any}.")
  }
  if (!type %in% c("sf", "sfc")) {
    riskmap_bug("{.arg type} must be {.val sf} or {.val sfc}.")
  }

  name <- arg

  if (type == "sf") {
    if (!inherits(data, "sf")) {
      abort_input("{.arg {name}} must be an {.cls sf} object, not {.obj_type_friendly {data}}.", arg = name)
    }
  } else if (!inherits(data, c("sf", "sfc"))) {
    abort_input("{.arg {name}} must be an {.cls sf} or {.cls sfc} object, not {.obj_type_friendly {data}}.", arg = name)
  }

  if (is.na(st_crs(data))) {
    abort_input("{.arg {name}} must have a coordinate reference system.", arg = name)
  }

  if (geometry != "any") {
    all_valid_geometry <- all(grepl(toupper(geometry), st_geometry_type(data)))
    if (!all_valid_geometry) {
      geometry_type <- switch(geometry, point = "POINT",
                              polygon = c("POLYGON", "MULTIPOLYGON"))
      abort_input("{.arg {name}} must have only {.or {.val {geometry_type}}} geometries.", arg = name)
    }
  }

  if (st_crs(data) == st_crs(4326)) {
    tryCatch(
      st_is_longlat(data$geometry),
      warning = function(w) {
        abort_input(c("{.arg {name}} has impossible longitudes or latitudes.", "i" = "Check that the coordinate columns were given in the right order when the data were converted."), arg = name)
      }
    )
  }
  invisible(TRUE)
}

#' Convert between CRS and requested distance units
#'
#' @param data An `sf` or `sfc` object with a projected CRS.
#' @param distance_units The requested coordinate units, either `"m"` or `"km"`.
#' @param arg The argument that gave the coordinate reference system, such as
#'   `model_crs`, or the data whose system it is.
#' @return The numeric factor converting one CRS unit to `distance_units`.
#' @importFrom units set_units
#' @noRd
crs_to_distance_factor <- function(data, distance_units, arg = "data") {
  crs_unit <- st_crs(data)$ud_unit

  if (is.null(crs_unit)) {
    abort_input(c("The coordinate reference system given by {.arg {arg}} has no linear units.", "i" = "Use a projected system with recognised linear units."), arg = arg)
  }

  tryCatch(
    as.numeric(
      set_units(crs_unit,
                distance_units,
                mode = "standard")
    ),
    error = function(e) {
      abort_input(c("The units of the coordinate reference system given by {.arg {arg}} can't be converted to {.val {distance_units}}.", "i" = "Use a projected system with recognised linear units."), arg = arg, parent = e)
    }
  )
}

#' Convert spatial coordinates to requested distance units
#'
#' @inheritParams crs_to_distance_factor
#' @return A numeric coordinate matrix expressed in `distance_units`.
#' @noRd
coordinates_in_units <- function(data, distance_units, arg = "data") {
  st_coordinates(data) * crs_to_distance_factor(data, distance_units, arg)
}

#' Preserve the caller's random-number state
#'
#' Capture the current random-number state and return a function that restores
#' it. If no state existed, the returned function removes any state subsequently
#' created. Callers should register the returned function with `on.exit()` before
#' calling `set.seed()`.
#'
#' @return A function that restores the captured random-number state.
#' @noRd
preserve_random_seed <- function() {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}

#' Check that an inverse link stays in the range of its family in the tails
#'
#' A draw whose mean is not a number or is outside the range has no weight
#' (the `log_likelihood()` of the family object), so a link that goes wrong only in
#' its tails would change the likelihood without any message. A mean that
#' overflows to infinity for a Poisson model is a limit and is allowed.
#' @noRd
check_invlink_tails <- function(invlink, family,
                                eta = seq(-50, 50, by = 0.5)) {
  mean <- invlink(eta)
  range <- family_traits(family)$mean_range
  outside <- is.na(mean) | mean < range[1L] | mean > range[2L]
  if (any(outside)) {
    at <- eta[which(outside)[1L]]
    abort_input("{.arg invlink} returns a value that isn't a number, or is outside the range of the family, at a linear predictor of {at}.", arg = "invlink")
  }
  invisible(TRUE)
}

#' @title check_invlink
#' @description
#'
#' Applies `invlink` to `eta` and checks that the result is a valid mean for
#' `family`: one finite numeric value per linear predictor, with
#' probabilities in `[0, 1]` for `"binomial"` and non-negative values for
#' `"poisson"`.
#' @param invlink A function mapping the linear predictor to the mean.
#' @param family One of `"gaussian"`, `"binomial"` or `"poisson"`.
#' @param eta Numeric vector of linear predictor values. Defaults to a
#'   representative spread of values, for checking `invlink` on its own
#'   before any real linear predictor is available.
#' @return The computed mean values, `invlink(eta)`. Raises an error if they
#'   are invalid.
#' @noRd
check_invlink <- function(invlink, family, eta = seq(-5, 5, length.out = 11)) {
  if (!is.function(invlink)) {
    abort_input("{.arg invlink} must be a function, not {.obj_type_friendly {invlink}}.", arg = "invlink")
  }
  mu <- invlink(as.numeric(eta))
  if (!is.numeric(mu) || length(mu) != length(eta)) {
    abort_input("{.arg invlink} must return one number for each linear predictor.", arg = "invlink")
  }
  if (any(!is.finite(mu))) {
    abort_input(c("{.arg invlink} returns means that aren't finite.", "i" = "Check the model parameters and covariates."), arg = "invlink")
  }
  traits <- family_traits(family)
  if (any(mu < traits$mean_range[1L] | mu > traits$mean_range[2L])) {
    abort_input(traits$range_message, arg = "invlink")
  }
  mu
}
