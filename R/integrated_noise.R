##' The plan of a fit whose independent effects are integrated out
##'
##' Records at one location with the same covariates, offset and random-effect
##' levels have the same linear predictor, so they are summed into one unit,
##' which leaves the likelihood unchanged up to a constant for any link of the
##' binomial and Poisson families. The nugget is treated as a term whose levels
##' are the locations, and each `re()` term as a term with its own levels.
##' `allowed` says, for the nugget and each `re()` term in that order, whether
##' it may be integrated out. Each level and location is then decided by three
##' rules.
##'
##' (a) A level that holds exactly one unit is integrated out of that unit's
##' likelihood, whose variance is the sum of the variances of all such levels.
##' At a location with a single unit this includes the nugget.
##' (b) At a location with several units, the nugget is integrated out as well
##' when the units share covariates and offset, every unit is the only unit of
##' its level in the same `re()` terms (the inner terms), and every other term
##' has one level shared by all the units. The location is then one nested
##' record, whose outer variance is the nugget's plus that of any `re()` term
##' whose level holds exactly the location's units, and whose inner variance is
##' the sum of the inner terms' variances (`site_table.h`).
##' (c) Otherwise the location's nugget stays sampled, and its levels with one
##' unit are still integrated out by (a).
##'
##' Every other level stays sampled, including the levels that span several
##' locations and those that hold several distinct units. The levels of a term
##' are independent, so a term can be integrated out at some levels and sampled
##' at others.
##'
##' Returns the engine's records, one per unit outside nested locations and one
##' per nested location (`y`, `units_m`, `D`, `cov_offset`, `ID_coords`, and
##' `ID_re` with NA for an inner term of a nested record), how each is
##' evaluated (`kind` 0 for the plain kernel, 1 for a record and 2 for a nested
##' pattern, with the classes `a` and `b` of their variances and the units'
##' `patterns`), the `classes` as a logical matrix of terms by class, the
##' latent column of each record's sampled levels (`latent_columns`, 0 where a
##' level is integrated out) with the number of sampled levels of each term
##' (`sampled_levels`), `nugget_integrated` for each location, the records in
##' which each term is an inner one (`inner`), the counts of levels integrated
##' out by term (`counts`), the terms integrated out anywhere
##' (`integrated_terms`) and whether there is any (`any`).
##' @noRd
plan_integrated_records <- function(y, units_m, D, cov_offset, ID_coords,
                                    ID_re, n_loc, allowed) {
  n_re <- if (is.null(ID_re)) 0L else ncol(ID_re)
  terms <- names(allowed)
  if (length(allowed) != n_re + 1L) {
    riskmap_bug("{.arg allowed} must have one element for the nugget and one for each random effect.")
  }
  key <- do.call(paste, c(
    list(ID_coords), as.data.frame(D), list(cov_offset),
    if (n_re) as.data.frame(ID_re), sep = "\r"
  ))
  group <- match(key, unique(key))
  first <- !duplicated(group)
  unit_y <- as.numeric(rowsum(y, group, reorder = FALSE))
  unit_m <- as.numeric(rowsum(units_m, group, reorder = FALSE))
  unit_D <- D[first, , drop = FALSE]
  unit_offset <- cov_offset[first]
  unit_location <- ID_coords[first]
  unit_re <- if (n_re) ID_re[first, , drop = FALSE]
  n_units <- length(unit_y)
  per_location <- tabulate(unit_location, nbins = n_loc)
  # For each re() term, whether a unit is the only one of its level (single),
  # and whether its level lies at one location and holds all its units
  # (whole), both only where the term may be integrated out.
  single <- whole <- matrix(FALSE, n_units, n_re)
  levels_of <- integer(n_re)
  for (j in seq_len(n_re)) {
    level <- unit_re[, j]
    levels_of[j] <- max(ID_re[, j])
    size <- tabulate(level, nbins = levels_of[j])
    places <- tabulate(level[!duplicated(cbind(level, unit_location))],
                       nbins = levels_of[j])
    single[, j] <- allowed[[j + 1L]] & size[level] == 1L
    whole[, j] <- allowed[[j + 1L]] & places[level] == 1L &
      size[level] == per_location[unit_location]
  }
  nugget_allowed <- isTRUE(allowed[[1L]])
  # Rule (b), at the locations with several units.
  predictor <- do.call(paste, c(list(unit_location), as.data.frame(unit_D),
                                list(unit_offset), sep = "\r"))
  by_location <- split(seq_len(n_units), factor(unit_location,
                                                levels = seq_len(n_loc)))
  nested <- logical(n_loc)
  for (l in which(per_location > 1L)) {
    u <- by_location[[l]]
    inner <- single[u, , drop = FALSE]
    if (!any(inner[1L, ]) || length(unique(predictor[u])) > 1L ||
        any(inner != rep(inner[1L, ], each = length(u)))) next
    shared <- vapply(which(!inner[1L, ]), function(j) {
      length(unique(unit_re[u, j])) == 1L
    }, logical(1))
    nested[l] <- all(shared) && (nugget_allowed || any(whole[u[1L], ]))
  }
  nugget_integrated <- nugget_allowed & (per_location == 1L | nested)
  # The levels integrated out: those with one unit, and those holding all the
  # units of a nested location.
  in_nested <- nested[unit_location]
  integrated_level <- lapply(seq_len(n_re), function(j) {
    out <- logical(levels_of[j])
    out[unit_re[single[, j] | (whole[, j] & in_nested), j]] <- TRUE
    out
  })
  # The engine's records, in the order of their first unit.
  emit <- !in_nested | !duplicated(unit_location)
  record_unit <- which(emit)
  nested_record <- in_nested[record_unit]
  n_records <- length(record_unit)
  records_y <- unit_y[record_unit]
  records_m <- unit_m[record_unit]
  patterns <- vector("list", n_records)
  for (r in which(nested_record)) {
    u <- by_location[[unit_location[record_unit[r]]]]
    patterns[[r]] <- list(y = unit_y[u], m = unit_m[u])
    records_y[r] <- sum(unit_y[u])
    records_m[r] <- sum(unit_m[u])
  }
  records_re <- if (n_re) unit_re[record_unit, , drop = FALSE]
  if (n_re && any(nested_record)) {
    inner_terms <- single[record_unit[nested_record], , drop = FALSE]
    records_re[nested_record, ][inner_terms] <- NA_integer_
  }
  # The classes of variances, as the terms whose variances they sum. A record
  # takes the nugget where it is its location's only unit.
  record_terms <- cbind(nugget_integrated[unit_location[record_unit]] &
                          !nested_record,
                        single[record_unit, , drop = FALSE])
  outer_terms <- cbind(nested_record & nugget_allowed,
                       whole[record_unit, , drop = FALSE] & nested_record)
  inner_terms <- cbind(FALSE, single[record_unit, , drop = FALSE] &
                         nested_record)
  record_terms[nested_record, ] <- FALSE
  candidates <- unique(rbind(record_terms[!nested_record, , drop = FALSE],
                             outer_terms[nested_record, , drop = FALSE],
                             inner_terms[nested_record, , drop = FALSE]))
  candidates <- candidates[rowSums(candidates) > 0, , drop = FALSE]
  class_key <- function(rows) apply(rows, 1L, paste, collapse = "")
  candidate_key <- class_key(candidates)
  class_of <- function(rows) match(class_key(rows), candidate_key)
  kind <- ifelse(nested_record, 2L,
                 ifelse(rowSums(record_terms) > 0, 1L, 0L))
  a <- b <- rep(1L, n_records)
  a[kind == 1L] <- class_of(record_terms[kind == 1L, , drop = FALSE])
  a[kind == 2L] <- class_of(outer_terms[kind == 2L, , drop = FALSE])
  b[kind == 2L] <- class_of(inner_terms[kind == 2L, , drop = FALSE])
  classes <- candidates
  dimnames(classes) <- list(NULL, terms)
  # Latent columns: the locations, then the sampled levels of each term.
  latent_columns <- cbind(unit_location[record_unit],
                          matrix(0L, n_records, n_re))
  sampled_levels <- integer(n_re)
  next_column <- n_loc
  for (j in seq_len(n_re)) {
    sampled <- !integrated_level[[j]]
    column <- integer(levels_of[j])
    column[sampled] <- next_column + seq_len(sum(sampled))
    sampled_levels[j] <- sum(sampled)
    next_column <- next_column + sum(sampled)
    level <- records_re[, j]
    latent_columns[!is.na(level), j + 1L] <- column[level[!is.na(level)]]
  }
  storage.mode(latent_columns) <- "integer"
  counts <- data.frame(
    term = terms,
    integrated = c(sum(nugget_integrated),
                   vapply(integrated_level, sum, integer(1))),
    levels = c(n_loc, levels_of),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(allowed[[1L]])) counts$integrated[1L] <- 0L
  list(y = records_y, units_m = records_m,
       D = unit_D[record_unit, , drop = FALSE],
       cov_offset = unit_offset[record_unit],
       ID_coords = unit_location[record_unit], ID_re = records_re,
       kind = kind, a = a, b = b, patterns = patterns, classes = classes,
       latent_columns = latent_columns, sampled_levels = sampled_levels,
       nugget_integrated = nugget_integrated,
       inner = lapply(seq_len(n_re), function(j) {
         which(nested_record & single[record_unit, j])
       }),
       counts = counts,
       integrated_terms = terms[counts$integrated > 0],
       any = any(counts$integrated > 0))
}

##' The table of integrated records for the records of an engine
##'
##' One compiled table object (`cpp_site_table_create()`) for the distinct
##' pairs (y, m) of the records, of the units of the nested patterns and of
##' those units summed, with `units` describing each record for
##' `cpp_site_evaluate()`: its `kind` (0 for the plain kernel, 1 for a record,
##' 2 for a nested pattern), its pair or pattern, and the classes `a` and `b`
##' of its variances. By default every record is one of kind 1 in the first
##' class, or of kind 0 where `plain`. `patterns` gives the units (`y` and
##' `m`) of each nested record, and identical patterns share one table.
##' Adaptive Gauss-Hermite quadrature uses 20 nodes. `family` is the family
##' object of the model's link (`make_family()`), and a custom link is passed
##' to the table as the arrays of its kernel (`kernel_arrays()`,
##' `compute_custom_kernel_arrays()`).
##' @noRd
make_site_table <- function(y, units_m, family, plain = rep(FALSE, length(y)),
                            kind = ifelse(plain, 0L, 1L),
                            a = rep(1L, length(y)), b = rep(1L, length(y)),
                            patterns = NULL) {
  nested <- which(kind == 2L)
  pattern_y <- lapply(patterns[nested], `[[`, "y")
  pattern_m <- lapply(patterns[nested], `[[`, "m")
  all_y <- c(y, unlist(pattern_y), vapply(pattern_y, sum, numeric(1)))
  all_m <- c(units_m, unlist(pattern_m), vapply(pattern_m, sum, numeric(1)))
  key <- paste(all_y, all_m)
  distinct <- !duplicated(key)
  pair_of <- function(y, m) match(paste(y, m), key[distinct])
  index <- pair_of(y, units_m)
  pattern_list <- list()
  summed <- integer()
  if (length(nested)) {
    pattern_pairs <- lapply(seq_along(nested), function(k) {
      sort(pair_of(pattern_y[[k]], pattern_m[[k]]))
    })
    pattern_key <- vapply(pattern_pairs, paste, character(1), collapse = " ")
    unique_pattern <- !duplicated(pattern_key)
    pattern_list <- pattern_pairs[unique_pattern]
    summed <- vapply(which(unique_pattern), function(k) {
      pair_of(sum(pattern_y[[k]]), sum(pattern_m[[k]]))
    }, integer(1))
    index[nested] <- match(pattern_key, pattern_key[unique_pattern])
  }
  kernel <- check_family_object(family)
  if (kernel$integrated) {
    riskmap_bug("The site table is built from the model's link, not from a kernel of integrated records.")
  }
  quadrature <- build_normal_quadrature(20L)
  custom <- kernel$kernel_arrays()
  pointer <- cpp_site_table_create(
    as.numeric(all_y[distinct]), as.numeric(all_m[distinct]),
    kernel$poisson, quadrature$nodes / sqrt(2),
    log(quadrature$weights) + 0.5 * log(pi), custom, pattern_list,
    as.integer(summed)
  )
  units <- cbind(kind = as.integer(kind), index = as.integer(index),
                 a = as.integer(a), b = as.integer(b))
  list(pointer = pointer, units = units)
}

##' The likelihood of records with effects integrated out, as a link
##'
##' A list in the place of the inverse link functions (`linkf`), with name
##' "integrated", for the variances `classes` of the classes of integrated
##' effects (`plan_integrated_records()`). The engine's mode search, global
##' proposal and weights dispatch on that name: `terms()` gives the
##' log-likelihood, score and curvature at one linear predictor per record,
##' `nonlinear_parts()` the part of each likelihood that is not linear in a
##' step from a centre, with its slope, and
##' `log_likelihood()` the log-likelihood of each row of a matrix of linear
##' predictors with one column per record. A class whose variance is zero
##' takes the plain kernel.
##' @noRd
make_integrated_kernel <- function(sites, classes) {
  pointer <- sites$pointer
  units <- sites$units
  terms <- function(eta) {
    values <- cpp_site_evaluate(pointer, classes, units,
                                matrix(eta, ncol = 1L), TRUE, 2L)
    list(log_likelihood = sum(values$h), score = as.numeric(values$h1),
         curvature = -as.numeric(values$h2))
  }
  nonlinear_parts <- function(centre, step) {
    at_centre <- cpp_site_evaluate(pointer, classes, units,
                                   matrix(centre, ncol = 1L), TRUE, 1L)
    around <- cpp_site_evaluate(pointer, classes, units, centre + step, TRUE,
                                0L)$h
    slope <- as.numeric(at_centre$h1)
    list(nonlinear = around - as.numeric(at_centre$h) - step * slope,
         slope = slope)
  }
  list(name = "integrated", classes = classes, pointer = pointer,
       units = units, terms = terms, nonlinear_parts = nonlinear_parts,
       log_likelihood = function(eta) {
         cpp_site_row_sums(pointer, classes, units, eta)
       })
}

##' A smooth floor at zero for the curvature of the global proposal's factors
##'
##' epsilon log(1 + exp(C / epsilon)), which equals C to within
##' epsilon exp(-C / epsilon) for C well above epsilon and is positive
##' everywhere. A floor with continuous derivatives keeps the fixed refits a
##' smooth function of the parameters wherever a record's log-likelihood is
##' not concave, for integrated sites, where that cannot be ruled out
##' numerically, and for links not known to be log-concave
##' (`choose_curvature_floor()`).
##' @noRd
smooth_curvature_floor <- function(curvature, epsilon = 1e-3) {
  pmax(curvature, 0) + epsilon * log1p(exp(-abs(curvature) / epsilon))
}

##' Whether a link object is that of integrated sites
##' @noRd
is_integrated_link <- function(linkf) identical(linkf$name, "integrated")
