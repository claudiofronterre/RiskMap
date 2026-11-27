`%||%` <- function(a, b) if (!is.null(a)) a else b

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


##' @title Matern Correlation Function
##' @description Computes the Matern correlation function.
##' @param u A vector or matrix of distances between pairs of data locations.
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @param return_sym_matrix A logical value indicating whether to return a symmetric correlation matrix. Defaults to `FALSE`.
##' @details The Matern correlation function is defined as
##' \deqn{\rho(u; \phi, \kappa) = \frac{1}{2^{\kappa-1}\Gamma(\kappa)}(u/\phi)^\kappa K_{\kappa}(u/\phi)}
##' for \eqn{u > 0}, with \eqn{\rho(0; \phi, \kappa) = 1}, where \eqn{\phi} and \eqn{\kappa} are the scale and smoothness parameters, \eqn{\Gamma(\cdot)} is the gamma function, and \eqn{K_{\kappa}(\cdot)} denotes the modified Bessel function of the third kind of order \eqn{\kappa}. The parameters \eqn{\phi} and \eqn{\kappa} must be positive.
##' @return A vector or matrix of the same shape as `u` with the values of the Matern correlation function for the given distances, if `return_sym_matrix=FALSE`. If `return_sym_matrix=TRUE`, a symmetric correlation matrix is returned.
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
##' @noRd
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
##' @noRd
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
##' @param nugget Whether the model has a nugget, independent variation at
##' each location whose variance is `tau2`. Default is `FALSE`. With `TRUE` its
##' variance is estimated, unless `fixed_pars = list(tau2 = )` of [glgpm()]
##' holds it.
##' @details The function constructs a list that includes the specified terms (spatial coordinates or covariates),
##' the smoothness parameter \eqn{\kappa}, and the nugget effect. This list can be used as a specification for a Gaussian Process model.
##' A model has one spatial Gaussian process, so a formula can contain only one `gp()` term.
##' @return A list of class `RiskMap_gp_spec` containing the following elements:
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

  if (is.numeric(nugget) && length(nugget) == 1L && is.finite(nugget) &&
      nugget > 0) {
    abort_input(c("{.arg nugget} must be {.code TRUE} or {.code FALSE}, not {nugget}.", "i" = "To hold the nugget variance at {nugget}, use {.code gp(nugget = TRUE)} with {.code fixed_pars = list(tau2 = {nugget})} in {.fn glgpm}."), arg = "nugget")
  }
  if (identical(nugget, 0) || identical(nugget, 0L)) {
    abort_input(c("{.arg nugget} must be {.code TRUE} or {.code FALSE}, not 0.", "i" = "A model without a nugget has {.code nugget = FALSE}."), arg = "nugget")
  }
  if (!(isTRUE(nugget) || isFALSE(nugget))) {
    supplied <- describe_value(nugget)
    abort_input("{.arg nugget} must be {.code TRUE} or {.code FALSE}, not {supplied}.", arg = "nugget")
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
##' @return A list of class `RiskMap_re_spec` containing the following elements:
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

##' Whether the fixed effects of a model are finite numbers in each row
##'
##' A row fails when a numeric column of its model frame isn't finite, such as
##' `log(x)` at `x = 0`. When the model frame can't be built, as for
##' `ns(log(x), 2)`, the parts inside the terms that use the data are checked
##' instead, and the model frame of the rows whose parts are finite. Missing
##' values fail too.
##' @param terms The terms of the fixed effects, without a response. Those of
##'   a model keep the centre and scale of terms such as `scale(x)` in their
##'   `predvars`.
##' @param data A data frame or `sf` object with the variables of the terms.
##' @return A logical vector with an element for each row of `data`, all `TRUE`
##'   when neither can be evaluated, since other checks report that.
##' @noRd
fixed_effects_are_finite <- function(terms, data) {
  if (inherits(data, "sf")) data <- st_drop_geometry(data)
  data <- as.data.frame(data)
  every <- rep(TRUE, nrow(data))
  if (!nrow(data)) return(every)
  frame_is_finite <- function(rows) {
    frame <- tryCatch(
      suppressWarnings(model.frame(terms, data = data[rows, , drop = FALSE],
                                   na.action = na.pass)),
      error = function(error) NULL
    )
    if (is.null(frame)) return(NULL)
    finite <- rep(TRUE, sum(rows))
    for (column in frame) {
      if (is.numeric(column)) {
        finite <- finite & rowSums(!is.finite(as.matrix(column))) == 0L
      }
    }
    finite
  }
  finite <- frame_is_finite(every)
  if (!is.null(finite)) return(finite)
  parts <- formula_parts_are_finite(terms, data)
  if (!any(parts)) return(parts)
  finite <- if (!all(parts)) frame_is_finite(parts)
  if (is.null(finite)) return(every)
  parts[parts] <- finite
  parts
}

##' The terms of the fixed effects of a model formula
##'
##' @param formula A model formula, with `gp()` and perhaps `re()` terms.
##' @return The terms without a response, or `NULL` when the formula can't be
##'   read, which other checks report.
##' @noRd
fixed_effect_terms <- function(formula) {
  tryCatch(stats::terms(update(interpret_formula(formula)$pf, NULL ~ .)),
           error = function(error) NULL)
}

##' Whether the parts of the terms of a formula are finite at each row
##'
##' The parts are the calls and variables inside the variables of the terms,
##' such as `log(x)` in `ns(log(x), 2)`, that use a column of the data and give
##' a number for each row. Constants such as the breaks of `cut()` are not
##' parts.
##' @param formula A formula or terms without a response.
##' @param data A data frame.
##' @return A logical vector with an element for each row of `data`.
##' @noRd
formula_parts_are_finite <- function(formula, data) {
  terms <- stats::terms(formula)
  finite <- rep(TRUE, nrow(data))
  check_parts <- function(expression) {
    if (!is.call(expression)) return(invisible(NULL))
    parts <- as.list(expression)[-1L]
    # Without empty arguments, as in x[, 1]
    for (part in parts[!vapply(parts, identical, logical(1), quote(expr = ))]) {
      if (!any(all.vars(part) %in% names(data))) next
      values <- tryCatch(
        suppressWarnings(eval(part, data, environment(terms))),
        error = function(error) NULL
      )
      if (is.numeric(values) && length(values) == length(finite)) {
        finite <<- finite & is.finite(as.numeric(values))
      }
      check_parts(part)
    }
  }
  for (variable in as.list(attr(terms, "variables"))[-1L]) {
    check_parts(variable)
  }
  finite
}

##' Stop because rows give the model a covariate that isn't a finite number
##'
##' @param rows The numbers of the rows.
##' @param arg The argument that holds them.
##' @param use `"fit"` for data that a model is fitted to or specified with,
##'   `"simulate"` for locations of a simulation.
##' @noRd
abort_non_finite_rows <- function(rows, arg, use = c("fit", "simulate")) {
  use <- match.arg(use)
  n <- length(rows)
  # As text, whose length sets the plurals, not its value.
  rows <- as.character(rows)
  advice <- if (use == "fit") {
    "{?It is row/They are rows} {rows}. Check or remove {?it/them}."
  } else {
    "{?It is row/They are rows} {rows}. Remove {?it/them} before simulating."
  }
  abort_input(c("{n} row{?s} of {.arg {arg}} {cli::qty(n)}{?gives/give} the model a covariate that isn't a finite number, such as {.code log(0)}.", "i" = advice), arg = arg)
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
  fixed <- fixed_effect_terms(formula)
  if (!is.null(fixed)) {
    finite <- fixed_effects_are_finite(fixed, data)
    if (!all(finite)) abort_non_finite_rows(which(!finite), "data", "fit")
  }

  invisible(TRUE)
}

##' @title Extract Parameter Estimates from a "RiskMap" Model Fit
##' @description This `coef` method for the "RiskMap" class extracts the
##' maximum likelihood estimates from model fits obtained from the [glgpm()] functions.
##' @param object An object of class "RiskMap" obtained as a result of a call to [glgpm()].
##' @param ... other parameters.
##' @return A named list of every parameter of the model on its natural
##' scale, the estimates and the values held by `fixed_pars`, in the order of
##' [riskmap_parameters]:
##' \item{beta}{The regression coefficients, named by the columns of the design matrix.}
##' \item{sigma2}{The variance of the spatial process.}
##' \item{phi}{The scale of the spatial correlation, in the model's distance units.}
##' \item{tau2}{The variance of the nugget, for a model with one.}
##' \item{sigma2_me}{The measurement error variance, for a Gaussian model.}
##' \item{theta}{The size, for a negative binomial model.}
##' \item{sigma2_re}{The variances of the random effects, named by their terms, for a model with `re()` terms.}
##' The list can be given as `start_pars` of [glgpm()].
##' @seealso [glgpm()], [riskmap_parameters]
##' @method coef RiskMap
##' @export
##'
coef.RiskMap <- function(object, ...) {
  # The natural values of the parameter table, the nugget as its variance,
  # with the held parameters at their values, in the order of the table.
  parameters <- fit_parameter_table(object)
  natural_values(parameters, object$estimate, held = TRUE)[coef_order(parameters)]
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

##' @title Log-Likelihood, Covariance and Number of Observations of a Fit
##' @description The `logLik()`, `vcov()` and `nobs()` methods for fits from
##'   [glgpm()], with which [stats::AIC()] and [stats::BIC()] work too.
##' @param object A fit from [glgpm()].
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##' @details
##' `logLik()` gives the log-likelihood at the estimates. It is exact for a
##' Gaussian model. For a binomial, Poisson or negative binomial model it is a
##' Monte Carlo estimate, so that differences between such models that are
##' not well above its Monte Carlo error mean nothing, which is why
##' [summary()][summary.RiskMap()] gives no AIC for them. That of a negative
##' binomial model leaves out the terms that do not depend on the parameters,
##' as that of a Poisson model does (see [negative_binomial()]). Its degrees of
##' freedom, `df`, are the number of estimated parameters, without those that
##' `fixed_pars` holds.
##'
##' `vcov()` gives the covariance matrix of the estimates on the scale of
##' [coef()][coef.RiskMap()], by the delta method from their covariance on the
##' scale on which the model is fitted, where the variances and the scale are
##' on the log scale. Its rows and columns are named as
##' `unlist(coef(object))` names them, for the estimated parameters only. A
##' parameter without a standard error, such as a variance estimated at zero,
##' has a row and a column of `NA`, and a fit that holds every parameter has an
##' empty matrix.
##'
##' `nobs()` gives the number of observations, the rows of the data, whatever
##' their denominators.
##' @return `logLik()` gives an object of class `logLik`, `vcov()` a matrix
##'   and `nobs()` an integer.
##' @seealso [summary.RiskMap()] for standard errors and confidence intervals.
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
##' logLik(fit)
##' AIC(fit)
##' vcov(fit)
##' nobs(fit)
##' @method logLik RiskMap
##' @export
logLik.RiskMap <- function(object, ...) {
  check_dots_unused(...)
  structure(as.numeric(object$log_lik),
            df = nrow(estimated_rows(fit_parameter_table(object))),
            nobs = nobs(object), class = "logLik")
}

##' @rdname logLik.RiskMap
##' @method vcov RiskMap
##' @export
vcov.RiskMap <- function(object, ...) {
  check_dots_unused(...)
  parameters <- fit_parameter_table(object)
  rows <- estimated_rows(parameters)
  names_out <- rows$reported
  natural <- unlist(coef(object))[names_out]
  working <- unlist(object$estimate)[rows$index]
  # The derivative of each parameter on the scale of coef() in its working
  # value, and of the nugget's variance in that of the spatial variance too,
  # since its working value is the logarithm of their ratio. `depends` marks
  # the working values each parameter depends on, whatever its derivative.
  n <- nrow(rows)
  jacobian <- matrix(0, n, n, dimnames = list(names_out, names_out))
  depends <- diag(n) == 1
  for (k in seq_len(n)) {
    jacobian[k, k] <- switch(
      rows$scale[k],
      identity = 1,
      log = ,
      nugget = natural[[k]],
      family = family_parameter(rows$name[k])$natural_derivative(working[[k]])
    )
  }
  nugget <- which(rows$scale == "nugget")
  sigma2 <- which(rows$name == "sigma2")
  if (length(nugget) && length(sigma2)) {
    jacobian[nugget, sigma2] <- natural[[nugget]]
    depends[nugget, sigma2] <- TRUE
  }
  # A parameter without a standard error has an NA row and column, as does
  # any parameter that depends on it.
  covariance <- as.matrix(object$covariance)[rows$index, rows$index,
                                             drop = FALSE]
  free <- !is.na(diag(covariance))
  computable <- rowSums(depends[, !free, drop = FALSE]) == 0
  out <- matrix(NA_real_, n, n, dimnames = list(names_out, names_out))
  out[computable, computable] <-
    jacobian[computable, free, drop = FALSE] %*%
    covariance[free, free, drop = FALSE] %*%
    t(jacobian[computable, free, drop = FALSE])
  out
}

##' @rdname logLik.RiskMap
##' @method nobs RiskMap
##' @export
nobs.RiskMap <- function(object, ...) {
  check_dots_unused(...)
  length(object$y)
}

##' @title Summarize Model Fits
##' @description Provides a `summary` method for the "RiskMap" class that
##'   computes standard errors and confidence intervals for likelihood-based
##'   model fits
##' @param object An object of class "RiskMap" from [glgpm()]
##' @param level Confidence level for intervals (default 0.95).
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##' @return A list of class `"summary.RiskMap"` with parameter estimates,
##'   standard errors, and confidence intervals. For binomial, Poisson and negative binomial
##'   models, its component `mcml` reports certified Monte Carlo maximum
##'   likelihood: `certificate` (an approximate bound, with probability
##'   `level`, on the Monte Carlo error of the estimate in statistical
##'   standard errors, `NA` when there is none), `tolerance`,
##'   `certified`, `outcome`, `correction_refused` (whether
##'   the Newton correction of the last stage did not converge),
##'   `stages`, `fitting_draws` and `validation_draws` (at the
##'   last stage), `certification_bank` (the effective draws and Pareto
##'   shape of the certification draws), `checks_to_review`, the number
##'   of checks of [check_computation()] that need a look,
##'   `starting_values`, where the starting values of the parameters not
##'   given in `start_pars` came from, `"current"` (RiskMap's
##'   defaults) or `"transformed_gaussian"` (a Gaussian fit to the
##'   transformed data, which `glgpm()` uses when the Laplace
##'   approximation of the likelihood prefers it to the defaults), or
##'   `"user"` when `start_pars` gives every parameter,
##'   `starting_given`, the names of those given in
##'   `start_pars` (absent when none was), and `starting_reason`,
##'   why the defaults were kept other than by being preferred in that
##'   comparison, such as `"re_terms"` when it is not made for a model
##'   with `re()` terms, `"not_compared"` when
##'   `set_control_mcml(start_comparison = FALSE)` skips it,
##'   `"candidates_given"` when `start_pars` gives every value the
##'   comparison could change,
##'   `"transformed_failed"` when the Gaussian fit
##'   failed, or `"not_recorded"` for a fit saved by an earlier version
##'   that did not record why (`NULL` when the comparison chose the
##'   start). `mcml` is
##'   `NULL` for Gaussian models. The parameters of the family, such as
##'   the measurement error variance of a Gaussian model, are in
##'   `family_parameters`, a table of those that were estimated, and
##'   `held_family_parameters`, a list of those held at a given value.
##'   The other parameters that `fixed_pars` holds keep their rows, with
##'   no interval, and `held_parameters` says at what values they are
##'   held (absent when none is).
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
summary.RiskMap <- function(object, level = 0.95, ...) {
  check_dots_unused(..., .hints = c(
    conf_level = "The confidence level is given in {.arg level}."
  ))

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

  # NULL for a family RiskMap doesn't know.
  traits <- fit_traits(object)
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
    if (!is.null(traits)) {
      link_name <- if (traits$exact || !traits$default_is_canonical) {
        traits$default_link
      } else {
        paste0("canonical (", traits$default_link, ")")
      }
      inv_expr <- paste("Inverse link function =", traits$inverse_link_formula)
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

  parameters <- fit_parameter_table(object)
  beta_flat_names <- paste0("beta.", beta_names)
  has_tau2      <- nugget_estimated(parameters)
  family_names  <- parameters$layout$family
  re_par_names  <- if (n_re > 0) paste0("sigma2_re.", re_names) else NULL

  ## tau2 = nu2 * sigma2, so on the log scale their raw estimates add. This is
  ## the one linear reparametrisation of the working-scale parameters that
  ## isn't just an exp(); its covariance is obtained via the delta method
  ## below (J, from the parameter table), alongside the other, unchanged
  ## parameters.
  J <- parameter_jacobian(parameters, nm)
  if (has_tau2) {
    estimate["nu2"] <- estimate["nu2"] + estimate["sigma2"]
  }

  ## The new parameters are J times the working-scale ones, so their
  ## covariance is J covariance J'. A parameter without a standard error (an
  ## NA row and column, for example a variance estimated at zero) is left
  ## out, with any new parameter that depends on it, so that the missing
  ## values do not spread to the other parameters.
  ## Held parameters are known, so a new parameter that depends on one, the
  ## nugget's variance with sigma2 held, takes no variance from it.
  covariance <- object$covariance
  held <- seq_len(ncol(J)) %in% held_index(parameters)
  free <- !is.na(diag(covariance)) & !held
  computable <- rowSums(J[, !free & !held, drop = FALSE] != 0) == 0 & !held
  covariance_new <- matrix(NA_real_, length(nm), length(nm),
                           dimnames = list(nm, nm))
  covariance_new[computable, computable] <-
    J[computable, free, drop = FALSE] %*%
    covariance[free, free, drop = FALSE] %*%
    t(J[computable, free, drop = FALSE])
  se_par <- sqrt(diag(covariance_new))

  # The family's own parameters go to their natural scale through the
  # family, and the other variances and the scale through exp().
  non_beta <- setdiff(nm, c(beta_flat_names, family_names))
  estimate[non_beta] <- exp(estimate[non_beta])
  descriptors <- traits$parameters
  for (name in family_names) {
    estimate[name] <- descriptors[[name]]$natural(estimate[[name]])
  }

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

  # Each is positive and fitted on the log scale of itself or of its inverse,
  # so its interval is log-normal, as those of the variances.
  if (length(family_names)) {
    est_family <- estimate[family_names]
    se_family  <- se_par[family_names]
    res$family_parameters <- cbind(
      Estimate      = est_family,
      "Lower limit" = exp(log(est_family) - z_crit * se_family),
      "Upper limit" = exp(log(est_family) + z_crit * se_family)
    )
    rownames(res$family_parameters) <- unname(table_labels(
      parameters, family_names
    ))
  }
  res$held_family_parameters <- held_family_parameters(object)
  # The other parameters fixed_pars holds, shown as fixed.
  pinned <- describe_pinned_parameters(parameters)
  if (length(pinned)) res$held_parameters <- pinned

  sp_names <- c("sigma2", "phi", if (has_tau2) "nu2")
  est_sp <- estimate[sp_names]
  se_sp  <- se_par[sp_names]
  res$sp <- cbind(
    Estimate      = est_sp,
    "Lower limit" = exp(log(est_sp) - z_crit * se_sp),
    "Upper limit" = exp(log(est_sp) + z_crit * se_sp)
  )
  # nu2 is reported as tau2, the variance of the nugget.
  rownames(res$sp) <- unname(table_labels(parameters,
                                          c("sigma2", "phi",
                                            if (has_tau2) "tau2")))
  nugget <- fit_fix_tau2(object)
  if (!is.null(nugget)) res$tau2 <- nugget

  if (n_re > 0) {
    est_re <- estimate[re_par_names]
    se_re  <- se_par[re_par_names]
    res$ranef <- cbind(
      Estimate      = est_re,
      "Lower limit" = exp(log(est_re) - z_crit * se_re),
      "Upper limit" = exp(log(est_re) + z_crit * se_re)
    )
    rownames(res$ranef) <- unname(table_labels(parameters, re_par_names))
  }

  res$level           <- level
  res$family          <- object$family
  res$kappa           <- object$kappa
  res$log_lik         <- object$log_lik
  res$cov_offset_used <- !(is.null(object$cov_offset) ||
                             all(object$cov_offset == 0))
  # Only an exact likelihood, without Monte Carlo error, has an AIC.
  if (isTRUE(traits$exact)) {
    # The parameters fixed_pars holds are not counted.
    res$aic <- 2 * nrow(estimated_rows(parameters)) - 2 * res$log_lik
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
  for (family in riskmap_families) {
    for (name in names(family$parameters)) {
      labels[names == name] <- family$parameters[[name]]$label
    }
  }
  random <- grepl("^sigma2_re\\.", names)
  labels[random] <- paste0(sub("^sigma2_re\\.", "", names[random]),
                           " (random eff. var.)")
  stats::setNames(labels, names)
}

##' The lines of `summary()` on the family's own parameters
##'
##' The table of those the fit estimates, under the family's
##' `parameters_heading`, and a line for each that it holds at a given value.
##' Nothing for a family without parameters of its own.
##' @noRd
format_family_parameters <- function(table, held, traits) {
  held_lines <- describe_held_family_parameters(held, traits)
  c(
    if (!is.null(table)) {
      c("", out_heading(traits$parameters_heading),
        out_capture(printCoefmat(table, P.values = FALSE,
                                 has.Pvalue = FALSE)))
    },
    if (length(held_lines)) c("", out_text(held_lines))
  )
}

##' One sentence for each of the family's own parameters held at a value, as
##' in "Measurement error var. fixed at 0.1"
##' @noRd
describe_held_family_parameters <- function(held, traits) {
  vapply(names(held), function(name) {
    paste0(traits$parameters[[name]]$label, " fixed at ", held[[name]])
  }, character(1), USE.NAMES = FALSE)
}

##' One sentence for each parameter `fixed_pars` holds in the working vector,
##' as in "Spatial process var. fixed at 0.5"
##' @noRd
describe_pinned_parameters <- function(parameters) {
  pinned <- pinned_rows(parameters)
  if (!nrow(pinned)) return(character())
  paste0(pinned$label, " fixed at ", pinned$held)
}

##' @title Print Summary of RiskMap Model
##' @description Print method for objects of class `"summary.RiskMap"`.
##' @param x An object of class `"summary.RiskMap"`.
##' @param ... Not used, and an error if given, so that a misspelt argument
##'   is not ignored.
##' @return Invisibly returns `x`.
##' @method print summary.RiskMap
##' @export
print.summary.RiskMap <- function(x, ...) {
  check_dots_unused(...)
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
                    out_as_bullets(
                      describe_starting_values(x$mcml$starting_values,
                                               x$mcml$starting_reason,
                                               x$mcml$starting_given), "i"
                    ),
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
    format_family_parameters(x$family_parameters, x$held_family_parameters,
                             traits),
    "",
    out_heading(paste0("Spatial Gaussian process, Matern correlation with ",
                       "kappa = ", x$kappa)),
    coefficients(x$sp, FALSE),
    if (!nugget_estimated && is.numeric(x$tau2)) {
      out_text(paste0("Variance of the nugget effect fixed at ", x$tau2))
    },
    if (length(x$held_parameters)) c("", out_text(x$held_parameters)),
    if (!is.null(x$ranef)) {
      c("", out_heading("Unstructured random effects"),
        coefficients(x$ranef, FALSE))
    },
    "",
    out_pairs(c("Log-likelihood", if (!is.null(x$aic)) "AIC"),
              format(c(x$log_lik, x$aic))),
    monte_carlo)
}

##' Names written as code and joined, as "`phi`, `sigma2` and `tau2`"
##' @noRd
format_code_list <- function(names) {
  names <- paste0("`", names, "`")
  if (length(names) < 2L) return(names)
  paste(paste(names[-length(names)], collapse = ", "), "and",
        names[length(names)])
}

##' One sentence on where the starting values of a fit came from, and why
##'
##' `selected` and `reason` are those of `select_nongaussian_start()`, and
##' `given` the flat names of the parameters given in `start_pars`.
##' `NULL` for a fit that doesn't record them, such as one saved by an
##' earlier version, which records no reason either. A fit saved before
##' `start_pars` could be partial records `"user"` without `given`.
##' @noRd
describe_starting_values <- function(selected, reason = NULL, given = NULL) {
  if (is.null(selected)) return(NULL)
  transformed <- "a Gaussian fit to the transformed data"
  defaults <- "The starting values were RiskMap's defaults"
  laplace <- "the Laplace approximation of the likelihood"
  if (length(given) && !identical(selected, "user")) {
    others <- sub("^The starting values", "The other starting values",
                  describe_starting_values(selected, reason))
    return(paste0("The starting values of ", format_code_list(given),
                  " were given in `start_pars`. ", others))
  }
  switch(
    selected,
    user = if (length(given)) {
      "The starting values were all given in `start_pars`."
    } else {
      paste("The starting values were those given in",
            "`start_pars`, with RiskMap's defaults for any others.")
    },
    transformed_gaussian = paste0("The starting values came from ",
                                  transformed, ", which ", laplace,
                                  " preferred to RiskMap's defaults."),
    current = if (is.null(reason)) {
      paste0(defaults, ", which ", laplace, " preferred to ", transformed,
             ".")
    } else {
      paste0(defaults, switch(
        reason,
        re_terms = paste0(", which a model with `re()` terms takes without",
                          " comparing them with ", transformed, "."),
        link = paste0(", which a model with a link other than the",
                      " canonical one takes without comparing them with ",
                      transformed, "."),
        transformed_failed = paste0(", because ", transformed, " failed."),
        transformed_invalid = paste0(", because ", transformed,
                                     " gave invalid values."),
        comparison_failed = paste0(", because ", laplace,
                                   " failed at one of the two starts."),
        transformed_not_finite = paste0(", because ", laplace,
                                        " was not finite at the start from ",
                                        transformed, "."),
        neither_finite = paste0(", because ", laplace,
                                " was not finite at either start."),
        not_compared = paste0(", which `set_control_mcml(start_comparison",
                              " = FALSE)` keeps without comparing them with ",
                              transformed, "."),
        candidates_given = paste0(", since `start_pars` gives every value ",
                                  "that the comparison with ", transformed,
                                  " could change."),
        not_recorded = ".",
        "."
      ))
    },
    NULL
  )
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
##' @param object A fit from [glgpm()] or its [summary()][summary.RiskMap()],
##' an assessment of fits from [assess_fit()] or its
##' [summary()][summary.RiskMap_predictive_check()], a cross-validation from
##' [assess_prediction()] or its [summary()][summary.RiskMap_cross_validation()],
##' or a prediction over `boundaries` from [predict()][predict.RiskMap()] (see
##' [to_table.RiskMap_prediction()]).
##' @param digits A non-negative integer giving the number of decimal places
##' used to display numeric results. Small values get more (see Details).
##' @param col_names Optional column names, replacing the default headers.
##' @param reference For a cross-validation, optional name of the model that the
##' others are compared with, as in its
##' [summary()][summary.RiskMap_cross_validation()].
##' @param ... Additional arguments passed to [knitr::kable()].
##' @details This function creates a presentation-ready summary table from a
##' fitted "RiskMap" model, an assessment of fits or cross-validation results for
##' multiple models.
##' Use [coef()][coef.RiskMap()] or [summary()][summary.RiskMap()] when numeric
##' results are required for further analysis.
##'
##' Numeric values have `digits` decimal places. Values that are compared, the
##' estimate and confidence limits of a parameter or one metric across models,
##' share one notation and one number of decimal places. When one of them is
##' smaller in size than `10^(-digits)`, they get as many decimal places as
##' the smallest needs for two significant digits. They are in scientific
##' notation when that takes more than six decimal places, or when one of them
##' is a million or more in size. Zeros and missing values don't count.
##'
##' The table of a fit or its summary has
##' - the regression coefficients with their estimates and confidence
##'   intervals,
##' - the parameters of the spatial process,
##' - the variances of the random effects,
##' - the parameters of the family that were estimated, such as the
##'   measurement error variance of a Gaussian model.
##'
##' The table of a cross-validation or its summary has a row for each model,
##' with its scores, their differences from a `reference` model with their
##' standard errors in brackets, and its `PIT_area` with the limit of a correct
##' model. A model of binary outcomes has "none (binary)" as its limit, because
##' the area is not a check of calibration there (see [assess_prediction()]).
##'
##' @return An object of class "knitr_kable" that can be rendered directly.
##' @importFrom knitr kable
##' @export
##' @seealso [glgpm()], [summary.RiskMap_cross_validation()]
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
to_table.RiskMap <- function(object, digits = 3, col_names = NULL, ...) {
  to_table(summary(object), digits = digits, col_names = col_names, ...)
}

##' @rdname to_table
##' @export
to_table.summary.RiskMap <- function(object, digits = 3, col_names = NULL,
                                     ...) {
  check_positive_integer(digits, "digits", allow_zero = TRUE)
  table <- rbind(object$reg_coef[, 1:3, drop = FALSE], object$sp,
                 object$ranef, object$family_parameters)
  format_results_table(table, digits = digits, col_names = col_names, ...)
}

##' @rdname to_table
##' @export
to_table.default <- function(object, ...) {
  abort_input("{.arg object} must be a {.cls RiskMap} fit or its summary, a {.cls RiskMap_predictive_check} or {.cls RiskMap_cross_validation} result or its summary, or a {.cls RiskMap_prediction} over {.arg boundaries}, not {.obj_type_friendly {object}}.", arg = "object")
}

##' The table of a fit, formatted for `kable()`
##'
##' The values that are compared, and so formatted together by
##' `format_compared_values()`, are the estimate and limits of a parameter, a
##' row of the table. `...` goes to [knitr::kable()].
##' @noRd
format_results_table <- function(table, digits, col_names, ...) {
  tab <- as.data.frame(table, check.names = FALSE)
  values <- as.matrix(tab)
  out <- array("", dim(values))
  for (i in seq_len(nrow(values))) {
    out[i, ] <- format_compared_values(values[i, ], digits)
  }
  tab[] <- out
  render_results_table(tab, row_names = TRUE, col_names, ...)
}

##' A table of text rendered by `kable()`
##'
##' The first column, or the row names, label the rows, and the others are
##' aligned to the right. `...` goes to [knitr::kable()].
##' @noRd
render_results_table <- function(tab, row_names, col_names, ...) {
  dots <- list(...)
  if (!is.null(col_names)) {
    if (length(col_names) != ncol(tab)) {
      abort_input("{.arg col_names} must have {ncol(tab)} elements, one for each column, not {length(col_names)}.", arg = "col_names")
    }
    dots$col.names <- col_names
  }
  if (is.null(dots$row.names))
    dots$row.names <- row_names
  if (is.null(dots$align))
    dots$align <- if (row_names) rep("r", ncol(tab)) else
      c("l", rep("r", ncol(tab) - 1L))

  do.call(kable, c(list(x = tab), dots))
}

##' Values that are compared, in one notation with one number of decimals
##'
##' They have `digits` decimal places. When a value is below `10^(-digits)` in
##' size, the decimals are those that show two significant digits of the
##' smallest, and scientific notation replaces more than six. A value of a
##' million or more has scientific notation too. Zeros and values that are not
##' finite don't count.
##' @noRd
format_compared_values <- function(x, digits) {
  digits <- as.integer(digits)
  sizes <- abs(x[is.finite(x) & x != 0])
  decimals <- digits
  if (any(sizes < 10^(-digits))) {
    # The exponent of the rounded value, since 0.000996 is 0.0010
    exponent <- sub(".*e", "", formatC(min(sizes), format = "e", digits = 1L))
    decimals <- 1L - as.integer(exponent)
  }
  if (any(sizes >= 1e6) || (decimals > digits && decimals > 6L)) {
    return(formatC(x, format = "e", digits = digits))
  }
  formatC(x, format = "f", digits = decimals)
}

##' @title Check for valid binomial values
##'
##' @description
##' Checks that binomial data only consists of zero or positive integers and that if
##' den is provided that it only consists of positive integers, and that all
##' values of y are less than or equal to it. `is_whole()` gives the tolerance,
##' so both have to be rounded before they are used.
##'
##' @param y the data to check
##' @param den the denominator
##' @return TRUE if valid, raise an error if not
##' @noRd
check_binomial <- function(y, den){
  if (!all(is_whole(y)) || !all(round(y) >= 0)) {
    abort_input("The responses of a binomial model must be non-negative whole numbers.", arg = "data")
  }

  if (!is.null(den)) {
    if (!all(is_whole(den)) || !all(round(den) > 0)) {
      abort_input(c("The denominators of a binomial model must be positive whole numbers.", "i" = "Check or remove the rows where they aren't."), arg = "denominator")
    }
    if (!all(round(den) >= round(y))) {
      abort_input("Binomial responses can't exceed their denominators.", arg = "denominator")
    }
  }

  invisible(TRUE)
}

##' Which numbers are finite and whole, with a tolerance for floating point errors
##' @noRd
is_whole <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  is.finite(x) & abs(x - round(x)) < tolerance
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
#'   `distance_crs`, or the data whose system it is.
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
