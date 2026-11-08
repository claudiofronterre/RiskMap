##' Records of a fit whose nugget is integrated out
##'
##' Records at one location with the same covariates, offset and random-effect
##' levels have the same linear predictor, so they are summed into one record,
##' which leaves the likelihood unchanged up to a constant for any link of the
##' binomial and Poisson families. A location that still has several records
##' shares its nugget among them, so the likelihood does not factorise by record
##' there, and the nugget stays sampled at that location (`several`). Returns
##' the summed records (`y`, `units_m`, `D`, `cov_offset`, `ID_coords`,
##' `ID_re`), the original record of each (`record`), and `several`, one
##' logical per location.
##' @noRd
prepare_integrated_records <- function(y, units_m, D, cov_offset, ID_coords,
                                       ID_re, n_loc) {
  key <- do.call(paste, c(
    list(ID_coords), as.data.frame(D), list(cov_offset),
    if (!is.null(ID_re)) as.data.frame(ID_re), sep = "\r"
  ))
  group <- match(key, unique(key))
  first <- !duplicated(group)
  summed_y <- as.numeric(rowsum(y, group, reorder = FALSE))
  summed_units <- as.numeric(rowsum(units_m, group, reorder = FALSE))
  coords <- ID_coords[first]
  several <- tabulate(coords, nbins = n_loc) > 1L
  list(y = summed_y, units_m = summed_units, D = D[first, , drop = FALSE],
       cov_offset = cov_offset[first], ID_coords = coords,
       ID_re = if (!is.null(ID_re)) ID_re[first, , drop = FALSE],
       record = group, several = several)
}

##' The table of integrated sites for a set of records
##'
##' One compiled table object (`cpp_site_table_create()`) for the distinct
##' pairs (y, m) of the records, with `pair` giving each record's 1-based pair,
##' negative for a record whose nugget is sampled (`plain`), which takes the
##' plain kernel. Adaptive Gauss-Hermite quadrature uses 20 nodes. A custom link
##' (`linkf`) is passed to the table as the arrays of its kernel
##' (`compute_custom_kernel_arrays()`).
##' @noRd
make_site_table <- function(y, units_m, family, plain = rep(FALSE, length(y)),
                            linkf = NULL) {
  key <- paste(y, units_m)
  distinct <- !duplicated(key)
  pair <- match(key, key[distinct])
  quadrature <- build_normal_quadrature(20L)
  custom <- if (!is.null(linkf) && !identical(linkf$name, "canonical")) {
    compute_custom_kernel_arrays(linkf, family)
  }
  pointer <- cpp_site_table_create(
    as.numeric(y[distinct]), as.numeric(units_m[distinct]),
    identical(family, "poisson"), quadrature$nodes / sqrt(2),
    log(quadrature$weights) + 0.5 * log(pi), custom
  )
  list(pointer = pointer, pair = ifelse(plain, -pair, pair))
}

##' The kernel of a custom link on a lattice of linear predictors
##'
##' Each record's log-likelihood is y A(eta) + m B(eta) + c, with A and B the
##' same for every record, so one set of arrays serves all of them. For the
##' binomial family A = log p - log(1 - p) and B = log(1 - p), and for the
##' Poisson family A = log(mean) and B = -mean (with c = y log m). The
##' arrays hold A, B and their first two derivatives from the link's own
##' derivatives, on the lattice from `low` to `high` with step 0.25 / 128, so
##' the coarsest steps of the table's rule fall on lattice points. As in the
##' sampled form (`bound_log_probability()`), a probability that rounds to a
##' limit is charged at the bound of 2^-53, where that logarithm is flat. Where
##' the link is not finite, beyond the range on which it was checked, the
##' kernel is continued from the nearest point where it is.
##' @noRd
compute_custom_kernel_arrays <- function(linkf, family, low = -100,
                                         high = 80, step = 0.25 / 128) {
  eta <- low + step * (0:round((high - low) / step))
  mean <- suppressWarnings(linkf$inv(eta))
  first <- suppressWarnings(linkf$d1(eta))
  second <- suppressWarnings(linkf$d2(eta))
  invalid <- !is.finite(mean) | !is.finite(first) | !is.finite(second) |
    mean < 0 | (family == "binomial" & mean > 1)
  if (all(invalid)) stop("The custom link is not finite anywhere.", call. = FALSE)
  if (any(invalid)) {
    valid <- which(!invalid)
    nearest <- valid[pmax(1L, findInterval(seq_along(eta), valid))]
    below <- seq_along(eta) < valid[1L]
    nearest[below] <- valid[1L]
    mean[invalid] <- mean[nearest[invalid]]
    first[invalid] <- 0
    second[invalid] <- 0
  }
  if (family == "binomial") {
    bound <- 2^-53
    success <- mean > bound
    failure <- 1 - mean > bound
    log_success <- log(pmax(mean, bound))
    log_failure <- log(pmax(1 - mean, bound))
    success_first <- ifelse(success, first / mean, 0)
    success_second <- ifelse(success, second / mean - (first / mean)^2, 0)
    failure_first <- ifelse(failure, -first / (1 - mean), 0)
    failure_second <- ifelse(failure,
                             -second / (1 - mean) - (first / (1 - mean))^2, 0)
    list(low = low, step = step, A = log_success - log_failure,
         A1 = success_first - failure_first,
         A2 = success_second - failure_second, B = log_failure,
         B1 = failure_first, B2 = failure_second)
  } else {
    positive <- mean > 0
    list(low = low, step = step, A = log(pmax(mean, .Machine$double.xmin)),
         A1 = ifelse(positive, first / mean, 0),
         A2 = ifelse(positive, second / mean - (first / mean)^2, 0),
         B = -mean, B1 = -first, B2 = -second)
  }
}

##' The likelihood of records with the nugget integrated out, as a link
##'
##' A list in the place of the inverse link functions (`linkf`), with name
##' "integrated", for the nugget variance `v`. The engine's mode search, global
##' proposal and weights dispatch on that name: `terms()` gives the
##' log-likelihood, score and curvature at one linear predictor per record,
##' `nonlinear_parts()` the parts of `compute_nonlinear_parts()`, and
##' `log_likelihood()` the log-likelihood of each row of a matrix of linear
##' predictors with one column per record. With `v = 0` every record takes the
##' plain kernel.
##' @noRd
make_integrated_kernel <- function(sites, v) {
  pointer <- sites$pointer
  pair <- sites$pair
  terms <- function(eta) {
    values <- cpp_site_evaluate(pointer, v, pair, matrix(eta, ncol = 1L),
                                TRUE, 2L)
    list(log_likelihood = sum(values$h), score = as.numeric(values$h1),
         curvature = -as.numeric(values$h2))
  }
  nonlinear_parts <- function(centre, step) {
    at_centre <- cpp_site_evaluate(pointer, v, pair, matrix(centre, ncol = 1L),
                                   TRUE, 1L)
    around <- cpp_site_evaluate(pointer, v, pair, centre + step, TRUE, 0L)$h
    slope <- as.numeric(at_centre$h1)
    list(nonlinear = around - as.numeric(at_centre$h) - step * slope,
         slope = slope)
  }
  list(name = "integrated", v = v, pointer = pointer, pair = pair,
       terms = terms, nonlinear_parts = nonlinear_parts,
       log_likelihood = function(eta) {
         cpp_site_row_sums(pointer, v, pair, eta)
       })
}

##' A smooth floor at zero for the curvature of integrated sites
##'
##' epsilon log(1 + exp(C / epsilon)), which equals C to within
##' epsilon exp(-C / epsilon) for C well above epsilon and is positive
##' everywhere. A floor with continuous derivatives keeps the fixed refits a
##' smooth function of the parameters wherever an integrated site's
##' log-likelihood is not concave, which cannot be ruled out numerically.
##' @noRd
smooth_curvature_floor <- function(curvature, epsilon = 1e-3) {
  pmax(curvature, 0) + epsilon * log1p(exp(-abs(curvature) / epsilon))
}

##' Whether a link object is that of integrated sites
##' @noRd
is_integrated_link <- function(linkf) identical(linkf$name, "integrated")
