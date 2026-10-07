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
  # Check if the input is an sf object
  if (!inherits(sf_object, "sf")) {
    stop("`sf_object` must be an sf object")
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
  if (!is.numeric(y) || !is.numeric(m)) {
    stop("'y' and 'm' must be numeric")
  }
  if (anyNA(y) || anyNA(m)) {
    stop("'y' and 'm' cannot contain missing values")
  }
  sizes <- c(length(y), length(m))
  max_size <- max(sizes)
  if (max_size == 0) {
    return(numeric(0))
  }
  incompatible <- sizes != 1 & max_size %% sizes != 0
  if (any(incompatible)) {
    stop("'y' and 'm' must have compatible lengths")
  }
  y <- rep_len(y, max_size)
  m <- rep_len(m, max_size)

  if (any(m <= 0)) {
    stop("'m' must contain only positive values")
  }
  if (any(y < 0)) {
    stop("'y' must contain only non-negative values")
  }
  if (any(y > m)) {
    stop("'y' must be less than or equal to 'm'")
  }

  log((y + 0.5) / (m - y + 0.5))
}


##' @title EPSG of the UTM Zone
##' @description Suggests the EPSG code for the UTM zone where the majority of the data falls.
##' @param data An object of class \code{sf} containing the coordinates.
##' @details The function determines the UTM zone and hemisphere where the majority of the data points are located and proposes the corresponding EPSG code.
##' @return An integer indicating the EPSG code of the UTM zone.
##' @export
propose_utm <- function (data) {
  if (!inherits(data, "sf"))
    stop("'data' must be an object of class sf")
  if (is.na(st_crs(data)))
    stop("the CRS of the data is missing and must be specified; see ?st_crs")

  # Transform to WGS84 (EPSG:4326) to ensure coordinates are in lon/lat
  data <- st_transform(data, crs = 4326)

  # Calculate UTM Zone
  utm_z <- floor((st_coordinates(data)[, 1] + 180)/6) + 1
  utm_z_u <- unique(utm_z)

  if (length(utm_z_u) > 1) {
    tab_utm <- table(utm_z)
    if (all(diff(tab_utm) == 0))
      warning("An equal amount of locations falls in different UTM zones")
    utm_z_u <- as.numeric(names(which.max(tab_utm)))
  }

  # Determine Hemisphere (fixing the latitude check)
  # latitude 0 (the Equator) is treated as northern hemisphere, per UTM convention
  ns <- ifelse(st_coordinates(data)[, 2] >= 0, 1, -1)  # Use latitude, not longitude
  ns_u <- unique(ns)

  if (length(ns_u) > 1) {
    tab_ns <- table(ns_u)
    if (all(diff(tab_ns) == 0))
      warning("An equal amount of locations falls north and south of the Equator")
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
##' @param u A vector of distances between pairs of data locations.
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @param return_sym_matrix A logical value indicating whether to return a symmetric correlation matrix. Defaults to \code{FALSE}.
##' @details The Matern correlation function is defined as
##' \deqn{\rho(u; \phi; \kappa) = (2^{\kappa-1})^{-1}(u/\phi)^\kappa K_{\kappa}(u/\phi)}
##' where \eqn{\phi} and \eqn{\kappa} are the scale and smoothness parameters, and \eqn{K_{\kappa}(\cdot)} denotes the modified Bessel function of the third kind of order \eqn{\kappa}. The parameters \eqn{\phi} and \eqn{\kappa} must be positive.
##' @return A vector of the same length as \code{u} with the values of the Matern correlation function for the given distances, if \code{return_sym_matrix=FALSE}. If \code{return_sym_matrix=TRUE}, a symmetric correlation matrix is returned.
##' @export
matern_correlation <- function(u, phi, kappa, return_sym_matrix = FALSE) {
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
##' @param U A vector of distances between pairs of data locations.
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @return A matrix with the values of the first derivative of the Matern function with respect to \eqn{\phi} for the given distances.
##' @export
matern_gradient_phi <- function(U, phi, kappa) {
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

  n <- attr(U, "Size")
  grad.phi.mat <- matrix(NA, nrow = n, ncol = n)
  ind <- lower.tri(grad.phi.mat)
  grad.phi <- der.phi(as.numeric(U), phi, kappa)
  grad.phi.mat[ind] <-  grad.phi
  grad.phi.mat <- t(grad.phi.mat)
  grad.phi.mat[ind] <-  grad.phi
  diag(grad.phi.mat) <- rep(der.phi(0, phi, kappa), n)
  grad.phi.mat
}

##' @title Second Derivative with Respect to \eqn{\phi}
##' @description Computes the second derivative of the Matern correlation function with respect to \eqn{\phi}.
##' @param U A vector of distances between pairs of data locations.
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @return A matrix with the values of the second derivative of the Matern function with respect to \eqn{\phi} for the given distances.
##' @export
matern_hessian_phi <- function(U, phi, kappa) {
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
  n <- attr(U, "Size")
  hess.phi.mat <- matrix(NA, nrow = n, ncol = n)
  ind <- lower.tri(hess.phi.mat)
  hess.phi <- der2.phi(as.numeric(U), phi, kappa)
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
##' @return A list of class \code{RiskMap_gp_spec} containing the following elements:
##' \item{term}{A character vector of the specified terms.}
##' \item{kappa}{The smoothness parameter \eqn{\kappa}.}
##' \item{nugget}{The nugget effect.}
##' \item{dim}{The number of specified terms.}
##' \item{label}{A character string representing the full call for the GP model.}
##' @export
gp <- function (..., kappa = 0.5, nugget = FALSE) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  term <- NULL

  if((!is.numeric(kappa) || kappa <= 0)){
    stop("'kappa' must be positive.")
  }

  if(!(is.logical(nugget) || (is.numeric(nugget) && nugget > 0))) {
    stop("'nugget' must be either 'TRUE' or 'FALSE' or a positive real number")
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
  full.call <- paste("gp(", term[1], sep = "")
  if (d > 1)
    for (i in 2:d) full.call <- paste(full.call, ",", term[i],
                                      sep = "")
  label <- gsub("sf", "", paste(full.call, ")", sep = ""))
  ret <- list(term = term, kappa = kappa, nugget = nugget, dim = d,
              label = label)
  class(ret) <- "RiskMap_gp_spec"
  ret
}
##' @title Random Effect Model Specification
##' @description Specifies the terms for a random effect model.
##' @param ... Variables representing the random effects in the model.
##' @details The function constructs a list that includes the specified terms for the random effects. This list can be used as a specification for a random effect model.
##' @return A list of class \code{RiskMap_re_spec} containing the following elements:
##' \item{term}{A character vector of the specified terms.}
##' \item{dim}{The number of specified terms.}
##' \item{label}{A character string representing the full call for the random effect model.}
##' @note At least one variable must be provided as input.
##' @export
re <- function (...) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  term <- NULL

  if (d == 0) {
    stop("You need to provide at least one variable.")
  } else {
    if (d > 0) {
      for (i in 1:d) {
        term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      }
    }
    for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])),
                                   "term.labels")
  }
  full.call <- paste("re(", term[1], sep = "")
  if (d > 1)
    for (i in 2:d) full.call <- paste(full.call, ",", term[i],
                                      sep = "")
  label <- gsub("sf", "", paste(full.call, ")", sep = ""))
  ret <- list(term = term, dim = d, label = label)
  class(ret) <- "RiskMap_re_spec"
  ret
}

interpret_formula <- function(formula) {
  p.env <- environment(formula)
  tf <- terms.formula(formula, specials = c("gp", "re"))
  terms <- attr(tf, "term.labels")
  nt <- length(terms)

  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
  } else {
    response <- NULL
  }

  gp <- attr(tf, "specials")$gp
  re <- attr(tf, "specials")$re
  off <- attr(tf, "offset")
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
  gp_spec <- eval(parse(text = terms[gp]), envir = p.env)
  re_spec <- eval(parse(text = terms[re]), envir = p.env)

  if (length(off) > 0) {
    offset <- as.character(attr(tf, "variables")[[off[i] + 1]])[2]
  } else {
    offset <- NULL
  }

  if (length(terms[-c(gp, re)]) > 0) {
    pf <- paste(response, "~", paste(terms[-c(gp, re)], collapse = " + "))
  } else if (length(terms[-c(gp, re)]) == 0) {
    pf <- paste(response, "~ 1")
  }

  if (attr(tf, "intercept") == 0) {
    pf <- paste(pf, "-1", sep = "")
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
##' as passed by the caller (e.g. `predictors` from `setup_prediction()`,
##' `data` from `check_formula()`).
##' @param data The data to check.
##' @param columns The column names to check for missing data.
##' @return TRUE if there is no missing data, or raise an error if not.
##' @noRd
check_complete_data <- function(data, columns) {
  name <- deparse(substitute(data))
  drop_coords <- st_drop_geometry(data[, columns])
  if (any(!complete.cases(drop_coords))) {
    stop("'", name, "' contains rows with missing data - check or remove them", call. = FALSE)
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

  if(!inherits(formula, "formula")) {
    stop("'formula' must be a 'formula'
         object indicating the variables of the
         model to be fitted", call. = FALSE)
  }

  formula_terms <- unique(get_formula_terms(formula))
  column_names <- names(data)

  contains_gp <- !is.null(attr(terms(formula, specials = "gp"), "specials")$gp)

  if (!contains_gp){
    stop("The 'formula' must contain a Gaussian Process term, specified with 'gp()'", call. = FALSE)
  }

  if (!response_required) {
    formula_terms <- setdiff(formula_terms, all.vars(formula[[2L]]))
  }
  missing_columns <- setdiff(formula_terms, column_names)
  n_missing <- length(missing_columns)
  if (n_missing > 0){

    stop(paste0("The 'formula' term",
          ifelse(n_missing > 1, "s '", " '"),
          paste(missing_columns, collapse = "', '"),
          ifelse(n_missing > 1, "' are", "' is"),
         " not present in 'data'"), call. = FALSE)
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
##' \item{phi}{The estimate for the spatial range parameter \eqn{\phi}.}
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
##' @export
fitted.RiskMap <- function(object, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (type == "response") object$fitted_values else object$linear_predictors
}

##' @title Summarize Model Fits
##' @description Provides a \code{summary} method for the "RiskMap" class that
##'   computes standard errors and confidence intervals for likelihood-based
##'   model fits
##' @param object An object of class "RiskMap" from \code{\link{glgpm}}
##' @param ... other parameters.
##' @param conf_level Confidence level for intervals (default 0.95).
##' @return A list of class \code{"summary.RiskMap"} with parameter estimates,
##'   standard errors, and confidence intervals.
##' @method summary RiskMap
##' @export
summary.RiskMap <- function(object, ..., conf_level = 0.95) {

  alpha  <- 1 - conf_level
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
  if (!is.null(object$linkf) && is.list(object$linkf) &&
      is.function(object$linkf$inv)) {
    link_name <- object$linkf$name %||% "custom"
    inv_expr  <- tryCatch(
      paste0("Inverse link function = ",
             paste(deparse(body(object$linkf$inv), width.cutoff = 500L),
                   collapse = " ")),
      error = function(e) "Inverse link function = <user-supplied function>"
    )
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

  covariance <- object$covariance
  H_new          <- t(J) %*% solve(-covariance) %*% J
  covariance_new <- solve(-H_new)
  se_par         <- sqrt(diag(covariance_new))

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
      rownames(res$me) <- "Measurement error var."
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
  rownames(res$sp) <- c("Spatial process var.",
                        paste0("Spatial corr. scale (",
                               object$distance_units, ")"),
                        if (has_tau2) "Variance of the nugget")
  if (!is.null(object$fix_tau2)) res$tau2 <- object$fix_tau2

  if (n_re > 0) {
    est_re <- estimate[re_par_names]
    se_re  <- se_par[re_par_names]
    res$ranef <- cbind(
      Estimate      = est_re,
      "Lower limit" = exp(log(est_re) - z_crit * se_re),
      "Upper limit" = exp(log(est_re) + z_crit * se_re)
    )
    rownames(res$ranef) <- paste0(re_names, " (random eff. var.)")
  }

  res$conf_level      <- conf_level
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

  class(res) <- "summary.RiskMap"
  return(res)
}

##' @title Print Summary of RiskMap Model
##' @description Print method for objects of class \code{"summary.RiskMap"}.
##' @param x An object of class \code{"summary.RiskMap"}.
##' @param ... other parameters.
##' @return Invisibly returns \code{x}.
##' @method print summary.RiskMap
##' @export
print.summary.RiskMap <- function(x, ...) {

  if (!is.null(x$call)) {
    cat("Call:\n")
    cat(paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  }

  # ===========================================================================
  # STANDARD RISKMAP MODELS
  # ===========================================================================

  if (identical(x$family, "gaussian")) {
    cat("Linear geostatistical model\n")
  } else if (identical(x$family, "binomial")) {
    cat("Binomial geostatistical model\n")
  } else if (identical(x$family, "poisson")) {
    cat("Poisson geostatistical model\n")
  }

  if (!is.null(x$link_name))        cat("Link:", x$link_name, "\n")
  if (!is.null(x$invlink_expression)) cat(x$invlink_expression, "\n\n")

  cat("'Lower limit' and 'Upper limit' refer to ",
      x$conf_level * 100, "% confidence intervals\n", sep = "")

  cat("\nRegression coefficients\n")
  printCoefmat(x$reg_coef, P.values = TRUE, has.Pvalue = TRUE)
  if (isTRUE(x$cov_offset_used)) cat("Offset included in the linear predictor\n")

  if (identical(x$family, "gaussian")) {
    if (length(x$me) > 1) {
      cat("\n"); printCoefmat(x$me, P.values = FALSE, has.Pvalue = FALSE)
    } else {
      cat("\nMeasurement error var. fixed at ", x$me, "\n", sep = "")
    }
  }

  cat("\nSpatial Gaussian process\n")
  cat("Matern covariance parameters (kappa = ", x$kappa, ")\n", sep = "")

  printCoefmat(x$sp, P.values = FALSE, has.Pvalue = FALSE)
  if (!isTRUE(x$tau2))
    cat("Variance of the nugget effect fixed at ", x$tau2, "\n", sep = "")

  if (!is.null(x$ranef)) {
    cat("\nUnstructured random effects\n")
    printCoefmat(x$ranef, P.values = FALSE, has.Pvalue = FALSE)
  }

  cat("\nLog-likelihood: ", x$log_lik, "\n", sep = "")
  if (identical(x$family, "gaussian") && !is.null(x$aic))
    cat("AIC: ", x$aic, "\n", sep = "")

  return(invisible(x))
}

##' @title Format RiskMap Model and Validation Results as a Table
##' @description Converts a fitted "RiskMap" model or cross-validation
##' results into a table that renders directly in Quarto, R Markdown, HTML,
##' LaTeX and the R console.
##' @param object An object of class "RiskMap" resulting from a call to
##' \code{\link{glgpm}}, a "summary.RiskMap" object, or a
##' "summary.RiskMap_cross_validation" object.
##' @param digits A non-negative integer giving the number of decimal places
##' used to display numeric results.
##' @param ... Additional arguments passed to \code{\link[knitr]{kable}}.
##' @details This function creates a presentation-ready summary table from a
##' fitted "RiskMap" model or cross-validation results for multiple models.
##' Use \code{\link{coef}} or \code{\link{summary}} when numeric results are
##' required for further analysis. Numeric values use fixed notation with
##' \code{digits} decimal places, except when scientific notation is needed to
##' represent very large or very small values clearly.
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
##' \dontrun{
##' fit <- glgpm(y ~ x + gp(), data = example_data)
##' to_table(fit, digits = 3)
##' }
to_table <- function(object, digits = 3, ...) {
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
    stop("'object' must be a RiskMap model or RiskMap cross-validation result")
  }

  tab <- as.data.frame(tab, check.names = FALSE)
  numeric_columns <- vapply(tab, is.numeric, logical(1))
  tab[numeric_columns] <- lapply(
    tab[numeric_columns],
    function(x) {
      use_scientific <- is.finite(x) & x != 0 &
        (abs(x) >= 1e6 | abs(x) < 10^(-digits))
      out <- formatC(x, format = "f", digits = as.integer(digits))
      out[use_scientific] <- formatC(
        x[use_scientific],
        format = "e",
        digits = as.integer(digits)
      )
      out
    }
  )

  dots <- list(...)
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
##' @param data_sf An `sf` object containing geometrical data from which coordinates are extracted.
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
##' @export
##'
##'
create_ids <- function(data_sf) {
  if(!inherits(data_sf,
               what = c("sfc","sf"), which = FALSE)) {
    stop("The object passed to 'grid_pred' must be an object
         of class 'sfc'")
  }
  coords_o <- st_coordinates(data_sf)
  coords <- unique(coords_o)

  m <- nrow(coords_o)
  ID_coords <- sapply(1:m, function(i)
    which(coords_o[i,1]==coords[,1] &
            coords_o[i,2]==coords[,2]))
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
  if (is.null(metric_names)) stop("No metrics of predictive performance were computed when running 'assess_prediction'")
  n_metrics <- length(metric_names)

  res <- matrix(NA, ncol = n_metrics, nrow = n_models)
  colnames(res) <- metric_names
  rownames(res) <- model_names

  test_set_means <- list()

  n_subs <- length(object$model[[1]]$metric[[1]])
  w <- unlist(lapply(object$model[[1]]$metric[[1]], length))

  for (i in 1:n_models) {
    model_scores <- list()
    for (j in 1:n_metrics) {
      score_j <- rep(NA, n_subs)
      for (h in 1:n_subs) {
        score_j[h] <- mean(object$model[[i]]$metric[[j]][[h]])
      }
      model_scores[[j]] <- score_j
      res[i, j] <- sum(w * score_j) / sum(w)
    }
    test_set_means[[model_names[i]]] <- model_scores
  }

  overall_averages <- colMeans(res, na.rm = TRUE)

  # Attach additional attributes for printing
  attr(res, "test_set_means") <- test_set_means
  attr(res, "overall_averages") <- overall_averages
  attr(res, "view_all") <- view_all
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
##' @return This function is used for its side effect of printing to the console. It does not
##'         return a value.
##' @export
##' @method print summary.RiskMap_cross_validation
print.summary.RiskMap_cross_validation <- function(x, ...) {
  # Extract attributes
  test_set_means <- attr(x, "test_set_means")
  overall_averages <- attr(x, "overall_averages")
  view_all <- attr(x, "view_all")

  cat("Summary of Cross-Validation Scores\n")
  cat("----------------------------------\n")

  for (model_name in names(test_set_means)) {
    cat(sprintf("Model: %s\n", model_name))

    if (view_all) {
      # Print scores for each test set
      model_test_set_means <- test_set_means[[model_name]]
      n_test_sets <- length(model_test_set_means[[1]])  # Number of test sets (assumes all metrics have same length)

      for (test_set_idx in seq_len(n_test_sets)) {
        cat(sprintf("  Test Set %d:\n", test_set_idx))
        for (metric_idx in seq_along(model_test_set_means)) {
          metric_name <- colnames(x)[metric_idx]
          test_set_value <- model_test_set_means[[metric_idx]][test_set_idx]
          cat(sprintf("    %s: %.4f\n", metric_name, test_set_value))
        }
      }
    }

    # Print overall average across test sets for the model
    cat("  Overall average across test sets:\n")
    for (metric_idx in seq_along(overall_averages)) {
      metric_name <- colnames(x)[metric_idx]
      overall_avg <- x[model_name, metric_idx]
      cat(sprintf("    %s: %.4f\n", metric_name, overall_avg))
    }
    cat("\n")
  }
}

##' Build one spatial map of a predictive-performance metric
##'
##' @noRd
.plot_metric_map <- function(object, metric, model, ...) {

  if (!model %in% names(object$model)) {
    stop(paste("'model'", shQuote(model, type = "sh"), "was not found in 'object'"))
  }

  if (!metric %in% names(object$model[[model]]$metric)) {
    stop(paste("'metric'", shQuote(metric, type = "sh"), "was not computed for model", shQuote(model, type = "sh")))
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
  out <- ggplot(data = data_full) +
    geom_sf(aes(color = .data$value), size = 2) +
    ggtitle(paste("Visualizing", metric, "for model", model)) +
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
##'   `mode == "average"` and `combine_panels` is `TRUE`, a single-element
##'   list holding one combined plot.
##' @noRd
.plot_calibration_curve <- function(object, mode, test_set, model,
                                    combine_panels) {

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
    stop("No AnPIT or PIT data available for plotting.")

  y_label <- unique(plot_data$type)
  if (length(y_label) > 1) y_label <- "Calibration curve"

  id_line <- geom_abline(intercept = 0, slope = 1,
                         linetype = "dashed", colour = "red")

  if (mode == "average" && combine_panels) {
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
                   colour = if (mode == "all") as.factor(test_set) else NULL)) +
      geom_line() + id_line +
      labs(title = title_suffix, x = "", y = unique(df$type)) +
      theme_minimal() +
      guides(colour = guide_legend(title = "Test set"))
  }

  plots <- list()
  for (mname in model) {
    df_model <- dplyr::filter(plot_data, .data$model == mname)

    p <- switch(mode,
                average = {
                  avg <- df_model %>%
                    dplyr::group_by(.data$u_val) %>%
                    dplyr::summarize(value = mean(.data$value), .groups = "drop")
                  avg$type <- unique(df_model$type)
                  build_plot(avg, paste("Model", mname, ": average"))
                },
                single  = {
                  if (is.null(test_set))
                    stop("Provide `test_set` when mode = 'single'.")
                  df_ts <- dplyr::filter(df_model, test_set == test_set)
                  if (nrow(df_ts) == 0)
                    stop("No data for test_set ", test_set, " in model ", mname)
                  build_plot(df_ts,
                             paste("Model", mname, "- test set", test_set))
                },
                all     = build_plot(df_model,
                                     paste("Model", mname, "- all test sets")),
                stop("Invalid `mode`. Use 'average', 'single' or 'all'.")
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
##' metric (\code{"CRPS"}, \code{"SCRPS"}, \code{"AnPIT_area"}) that was
##' requested via \code{metrics} in that call, and/or a calibration curve
##' (\code{"AnPIT"}) - the Aggregated normalised Probability Integral
##' Transform curve for discrete families, or the empirical PIT curve for
##' Gaussian models - when that was requested. By default every available
##' plot is produced, one per requested model, so you don't need to know in
##' advance which metrics were computed. Each metric gets its own plot
##' (or, with more than one model, its own small grid of one panel per
##' model) rather than everything being squeezed into one combined grid.
##'
##' @param x A \code{RiskMap_cross_validation} object.
##' @param metric Character vector restricting which plot(s) to produce;
##'   one or more of \code{"CRPS"}, \code{"SCRPS"}, \code{"AnPIT_area"} (a
##'   spatial map of that score) and \code{"AnPIT"} (the calibration
##'   curve). Defaults to every metric available in \code{x}.
##' @param model Character vector of model names to include. Defaults to
##'   every model in \code{x$model}.
##' @param ... Additional \pkg{ggplot2} components (e.g.
##'   \code{scale_color_gradient()}), layered onto every spatial-map plot.
##'   Must be passed by position only after `model`, since the remaining
##'   arguments must be named.
##' @param mode For the calibration curve only: one of \code{"average"}
##'   (average curve across test sets, the default), \code{"single"} (one
##'   specific test set) or \code{"all"} (every test set separately).
##' @param test_set Integer; required when \code{mode = "single"}.
##' @param combine_panels Logical; when \code{mode = "average"}, draw the
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
                                          mode = "average", test_set = NULL,
                                          combine_panels = FALSE) {

  all_models <- names(x$model)
  if (is.null(model)) {
    model <- all_models
  } else {
    missing_model <- setdiff(model, all_models)
    if (length(missing_model) > 0) {
      stop(paste("'model'", shQuote(missing_model[1], type = "sh"), "was not found in 'object'"))
    }
  }

  spatial_metrics <- intersect(
    c("CRPS", "SCRPS", "AnPIT_area"),
    names(x$model[[model[1]]]$metric)
  )
  has_calibration <- !is.null(x$model[[model[1]]]$PIT) ||
    !is.null(x$model[[model[1]]]$AnPIT)
  available_metrics <- c(if (has_calibration) "AnPIT", spatial_metrics)

  if (is.null(metric)) {
    metric <- available_metrics
  } else {
    unavailable <- setdiff(metric, available_metrics)
    if (length(unavailable) > 0) {
      stop(paste("'metric'", shQuote(unavailable[1], type = "sh"),
                "was not computed for model", shQuote(model[1], type = "sh")))
    }
  }

  # Calibration curve first, then the spatial-map metrics in a fixed order,
  # so metrics are plotted in a stable, predictable order regardless of how
  # `metric`/the object's own field ordering happen to be arranged.
  ordered_metrics <- intersect(c("AnPIT", "CRPS", "SCRPS", "AnPIT_area"), metric)

  groups <- stats::setNames(lapply(ordered_metrics, function(m) {
    group_plots <- if (identical(m, "AnPIT")) {
      .plot_calibration_curve(x, mode = mode, test_set = test_set,
                              model = model, combine_panels = combine_panels)
    } else {
      stats::setNames(
        lapply(model, function(mod) .plot_metric_map(x, metric = m, model = mod, ...)),
        model
      )
    }
    .arrange_metric_group(group_plots)
  }), ordered_metrics)

  if (length(groups) == 1) {
    return(groups[[1]])
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
  stopifnot("'y' must only consist of zero or positive integers when 'family' is 'binomial'" = valid)

  if (!is.null(den)){
    valid <- all(den >= y)
    stopifnot("Values of 'den' must be greater than or equal to values of 'y'"= valid)
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
#' @param geometry whether to check that the data contains `"point"` (default) or
#' `"polygon"` (covering both polygons and multipolygons)
#' @param type whether to check that the data is `"sf"` (default) or
#' `"sfc"` (either sf or sfc)
#' @return TRUE if the data is valid. Raise an error if not.
#' @noRd
#'
check_data <- function(data, geometry = "point", type = "sf"){
  stopifnot("'geometry' must be either 'point' or 'polygon'" = geometry %in% c("point", "polygon"))
  stopifnot("'type' must be either 'sf' or 'sfc'" = type %in% c("sf", "sfc"))

  # extract name passed to function
  data_type <- paste0("'", deparse(substitute(data)), "'")

  geometry_type <- switch(geometry,
                          point = "'POINT'",
                          polygon = "'POLYGON' or 'MULTIPOLYGON'")

  if (type == "sf"){
    if (!inherits(data, "sf"))
      stop(paste(data_type, "must be of class 'sf'"), call. = FALSE)
  } else {
    if (!inherits(data, c("sf", "sfc")))
      stop(paste(data_type, "must be of class 'sf' or 'sfc'"), call. = FALSE)
  }

  if (is.na(st_crs(data)))
    stop(paste(data_type, "must contain a coordinate reference system"), call. = FALSE)

  all_valid_geometry <- all(grepl(toupper(geometry), st_geometry_type(data)))
  if (!all_valid_geometry)
    stop(paste(data_type, "can only contain", geometry_type, "geometry"), call. = FALSE)

  if (st_crs(data) == st_crs(4326)){
    tryCatch(
      st_is_longlat(data$geometry),
      warning = function(w) {
        stop(paste(data_type, "contains impossible latitude or longitude values -
             check you have specified the columns correctly when converting the data"), call. = FALSE)
      }
    )
  }
  invisible(TRUE)
}

#' Convert between CRS and requested distance units
#'
#' @param data An `sf` or `sfc` object with a projected CRS.
#' @param distance_units The requested coordinate units, either `"m"` or `"km"`.
#' @return The numeric factor converting one CRS unit to `distance_units`.
#' @importFrom units set_units
#' @noRd
crs_to_distance_factor <- function(data, distance_units) {
  crs_unit <- st_crs(data)$ud_unit

  if (is.null(crs_unit)) {
    stop(
      "The modelling CRS does not define linear coordinate units. ",
      "Use a projected CRS with recognised linear units.",
      call. = FALSE
    )
  }

  tryCatch(
    as.numeric(
      set_units(crs_unit,
                distance_units,
                mode = "standard")
    ),
    error = function(e) {
      stop(
        "The modelling CRS units cannot be converted to '",
        distance_units,
        "'. Use a projected CRS with recognised linear units.",
        call. = FALSE
      )
    }
  )
}

#' Convert spatial coordinates to requested distance units
#'
#' @inheritParams crs_to_distance_factor
#' @return A numeric coordinate matrix expressed in `distance_units`.
#' @noRd
coordinates_in_units <- function(data, distance_units) {
  st_coordinates(data) * crs_to_distance_factor(data, distance_units)
}

#' @title check_positive_integer
#' @description
#'
#' Check that a value is a single, positive integer and error if not
#' @param x the value to check
#' @param name the name of the parameter to return in error messages
#' @param allow_null whether `NULL` is permitted
#' @param allow_zero whether zero is permitted
#' @return TRUE if the data is valid. Raise an error if not.
#' @noRd
#'
check_positive_integer <- function(x, name, allow_null = FALSE,
                                   allow_zero = FALSE) {
  if (is.null(x) && allow_null) return(invisible(TRUE))
  description <- if (allow_zero) "non-negative" else "positive"
  invalid <- !is.numeric(x) || length(x) != 1L || is.na(x) ||
    !is.finite(x) || x %% 1 != 0 || x < as.integer(!allow_zero) ||
    x > .Machine$integer.max
  if (invalid) {
    stop("'", name, "' must be a single ", description, " integer",
         call. = FALSE)
  }
  invisible(TRUE)
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

#' @title check_positive_number
#' @description
#'
#' Check that a value is a single, positive number and error if not
#' @param x the value to check
#' @param type the type of value being checked. Defaults to `"starting"`.
#' @return TRUE if x is valid. Raise an error if not.
#' @noRd
#'
check_positive_number <- function(x, type = "starting ") {
  # extract name, removing any list
  name <- gsub('.*\\[\\["([^"]+)"\\]\\].*', "\\1", deparse(substitute(x)))

  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stop("The ", type, "value for '", name, "' must be a single positive number")
  }

  invisible(TRUE)
}

#' Check that a numeric value belongs to a range
#'
#' @param x The value to check.
#' @param min Lower endpoint.
#' @param max Upper endpoint.
#' @param allow_equal Whether either endpoint is allowed.
#' @param name The argument name to use in the error message. Defaults to the
#'   expression supplied as `x`.
#' @return `TRUE` invisibly when valid; otherwise raises an error.
#' @noRd
check_range <- function(x, min = -Inf, max = Inf, allow_equal = TRUE,
                        name = deparse(substitute(x))) {
  valid_number <- is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x)
  outside <- if (!valid_number) {
    TRUE
  } else if (allow_equal) {
    x < min || x > max
  } else {
    x <= min || x >= max
  }
  invalid <- !valid_number || outside
  if (invalid) {
    interval <- if (allow_equal) "between" else "strictly between"
    stop("'", name, "' must be a single finite number ", interval, " ",
         min, " and ", max,
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that a value belongs to the closed unit interval
#' @noRd
check_zero_one <- function(x, name = deparse(substitute(x))) {
  check_range(x, min = 0, max = 1, name = name)
}


#' @title check_logical
#' @description
#'
#' Check that a value is a single, non-missing logical (`TRUE` or `FALSE`)
#' and error if not
#' @param x the value to check
#' @return TRUE if x is valid. Raise an error if not.
#' @noRd
#'
check_logical <- function(x) {
  name <- deparse(substitute(x))

  if (!isTRUE(x) && !isFALSE(x)) {
    stop("'", name, "' must be either TRUE or FALSE", call. = FALSE)
  }

  invisible(TRUE)
}


#' @title check_crs
#' @description
#'
#' Check that a CRS is valid
#' @param crs the CRS to check
#' @param name the argument name to use in the error message. Defaults to the
#'   name of the variable passed as `crs`; callers wrapping this in another
#'   function should pass their own argument's name explicitly, since
#'   `substitute()` only sees the immediate call site.
#' @return TRUE if the CRS is valid. Raise an error if not.
#' @noRd
#'
check_crs <- function(crs, name = deparse(substitute(crs))){
  tryCatch(
    st_crs(crs),
    warning = function(w) {
      stop("The '", name, "' provided is not a valid CRS", call. = FALSE)
    },
    error = function(e){
      stop("The '", name, "' provided is not a valid CRS", call. = FALSE)
    }
  )
  invisible(TRUE)
}
