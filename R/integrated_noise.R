##' The plan of a fit whose independent effects are integrated out
##'
##' Records at one location with the same covariates, offset and random-effect
##' levels have the same linear predictor, so they are summed into one unit,
##' which leaves the likelihood unchanged up to a constant for any link of the
##' binomial and Poisson families. The nugget is treated as a term whose levels
##' are the locations, and each `re()` term as a term with its own levels.
##' `allowed` says, for the nugget and each `re()` term in that order, whether
##' it may be integrated out. Each level and location is then decided by four
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
##' (c) At a location with several units whose linear predictors differ, the
##' units are integrated out as groups when every term but at most one has a
##' single level shared by all of them. The outer terms are the nugget and any
##' `re()` term whose level holds exactly the location's units.
##'   - If no term varies, the location is one group of all its units under
##'     its outer terms, a record of kind 3.
##'   - If one `re()` term varies and each of its levels there lies within the
##'     location, its levels are the groups, such as households whose children
##'     differ in covariates. With outer terms the location is one record of
##'     kind 4, its groups nested under the outer variance. Without them, each
##'     level with several units is a record of kind 3, and a level with one
##'     unit is integrated out by (a).
##' (d) Otherwise the location's nugget stays sampled, and its levels with one
##' unit are still integrated out by (a).
##'
##' Rule (c) applies only where the units' linear predictors differ, since
##' units that share them and differ in one term only are those of rule (b), or
##' of rule (a) without outer terms. So a fit whose units share their
##' predictors at every location has the plan it had before rule (c).
##'
##' Every other level stays sampled, including the levels that span several
##' locations and those that hold several distinct units at a location where
##' more than one term varies. The levels of a term are independent, so a term
##' can be integrated out at some levels and sampled at others.
##'
##' Returns the engine's records, one per unit outside nested locations and
##' groups, one per nested location and one per group or location of groups
##' (`y`, `units_m`, `D`, `cov_offset`, `ID_coords`, and `ID_re` with NA for an
##' inner term of a nested record or a location of groups), how each is
##' evaluated (`kind` 0 for the plain kernel, 1 for a record, 2 for a nested
##' pattern, 3 for a group and 4 for a location of groups, with the classes `a`
##' and `b` of their variances, the units' `patterns` of a nested record, and
##' the `groups` of a group or location of groups, each a list of groups whose
##' members' `y`, `m`, rows of `D` and `offset` are given relative to the
##' record's own `D` and `cov_offset`), the `classes` as a logical matrix of
##' terms by class, the latent column of each record's sampled levels
##' (`latent_columns`, 0 where a level is integrated out) with the number of
##' sampled levels of each term (`sampled_levels`), `nugget_integrated` for
##' each location, the records in which each term is an inner one (`inner`),
##' the counts of levels integrated out by term (`counts`), the terms
##' integrated out anywhere (`integrated_terms`) and whether there is any
##' (`any`).
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
  # (whole), both only where the term may be integrated out, and at how many
  # locations each level lies (places).
  single <- whole <- matrix(FALSE, n_units, n_re)
  levels_of <- integer(n_re)
  places_of <- vector("list", n_re)
  for (j in seq_len(n_re)) {
    level <- unit_re[, j]
    levels_of[j] <- max(ID_re[, j])
    size <- tabulate(level, nbins = levels_of[j])
    places <- tabulate(level[!duplicated(cbind(level, unit_location))],
                       nbins = levels_of[j])
    places_of[[j]] <- places
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
  # Rule (c), at the other locations with several units whose predictors
  # differ: 1 where the location is one group, 2 where it is a location of
  # groups, and 3 where its groups of several units are records, with the
  # term whose levels are the groups.
  grouping <- integer(n_loc)
  group_term <- integer(n_loc)
  for (l in which(per_location > 1L & !nested)) {
    u <- by_location[[l]]
    if (length(unique(predictor[u])) == 1L) next
    varying <- which(vapply(seq_len(n_re), function(j) {
      length(unique(unit_re[u, j])) > 1L
    }, logical(1)))
    has_outer <- nugget_allowed || any(whole[u[1L], ])
    if (!length(varying)) {
      if (has_outer) grouping[l] <- 1L
      next
    }
    if (length(varying) > 1L || !isTRUE(allowed[[varying + 1L]]) ||
        any(places_of[[varying]][unit_re[u, varying]] != 1L)) next
    level <- unit_re[u, varying]
    if (!has_outer && !anyDuplicated(level)) next
    grouping[l] <- if (has_outer) 2L else 3L
    group_term[l] <- varying
  }
  nugget_integrated <- nugget_allowed &
    (per_location == 1L | nested | grouping %in% c(1L, 2L))
  # Each unit's record, named by the record's first unit: its own, its
  # location's where the location is one record, or its group's where groups
  # of several units are records.
  record_of <- seq_len(n_units)
  location_record <- nested | grouping %in% c(1L, 2L)
  at_record <- location_record[unit_location]
  first_unit <- vapply(by_location, function(u) {
    if (length(u)) u[1L] else NA_integer_
  }, integer(1))
  record_of[at_record] <- first_unit[unit_location[at_record]]
  # The units whose level of a term is a group of rule (c).
  in_group <- matrix(FALSE, n_units, n_re)
  for (l in which(grouping >= 2L)) {
    u <- by_location[[l]]
    level <- unit_re[u, group_term[l]]
    lead <- u[match(level, level)]
    several <- duplicated(level) | duplicated(level, fromLast = TRUE)
    in_group[u, group_term[l]] <- grouping[l] == 2L | several
    if (grouping[l] == 3L) record_of[u[several]] <- lead[several]
  }
  # The levels integrated out: those with one unit, those holding all the
  # units of a location that is one record, and the groups of rule (c).
  in_nested <- nested[unit_location]
  in_location_record <- at_record
  integrated_level <- lapply(seq_len(n_re), function(j) {
    out <- logical(levels_of[j])
    out[unit_re[single[, j] | (whole[, j] & in_location_record) |
                  in_group[, j], j]] <- TRUE
    out
  })
  # The engine's records, in the order of their first unit.
  emit <- record_of == seq_len(n_units)
  record_unit <- which(emit)
  nested_record <- in_nested[record_unit]
  record_grouping <- integer(length(record_unit))
  record_grouping[!nested_record] <-
    grouping[unit_location[record_unit[!nested_record]]]
  # A unit of a location of rule (c) that is not in a group of several units
  # is a record of its own, under rule (a).
  record_grouping[record_grouping == 3L & !(record_unit %in% record_of[
    duplicated(record_of)])] <- 0L
  location_group_record <- record_grouping == 2L
  group_record <- record_grouping %in% c(1L, 3L)
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
  # The groups of each group and location of groups, with their members'
  # design rows and offsets relative to the record's.
  groups <- vector("list", n_records)
  members_of <- split(seq_len(n_units), record_of)
  describe <- function(u, lead) {
    list(y = unit_y[u], m = unit_m[u],
         D = sweep(unit_D[u, , drop = FALSE], 2L, unit_D[lead, ]),
         offset = unit_offset[u] - unit_offset[lead])
  }
  for (r in which(group_record | location_group_record)) {
    lead <- record_unit[r]
    u <- members_of[[as.character(lead)]]
    records_y[r] <- sum(unit_y[u])
    records_m[r] <- sum(unit_m[u])
    groups[[r]] <- if (group_record[r]) {
      list(describe(u, lead))
    } else {
      level <- unit_re[u, group_term[unit_location[lead]]]
      lapply(split(u, factor(level, levels = unique(level))), describe,
             lead = lead)
    }
  }
  records_re <- if (n_re) unit_re[record_unit, , drop = FALSE]
  if (n_re && any(nested_record)) {
    inner_terms <- single[record_unit[nested_record], , drop = FALSE]
    records_re[nested_record, ][inner_terms] <- NA_integer_
  }
  # The group term of each location of groups, as a logical matrix.
  grouped_term <- matrix(FALSE, n_records, n_re)
  if (any(location_group_record | group_record)) {
    term <- group_term[unit_location[record_unit]]
    marked <- (location_group_record | group_record) & term > 0L
    grouped_term[cbind(which(marked), term[marked])] <- TRUE
  }
  if (n_re && any(location_group_record)) {
    records_re[location_group_record, ][
      grouped_term[location_group_record, , drop = FALSE]] <- NA_integer_
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
  # A group of a location that is one group takes the outer terms, and a
  # group of a level its term. A location of groups takes the outer terms
  # outside and the groups' term inside.
  location_outer <- cbind(nugget_allowed & (record_grouping %in% c(1L, 2L)),
                          whole[record_unit, , drop = FALSE] &
                            record_grouping %in% c(1L, 2L))
  group_terms <- location_outer
  group_terms[record_grouping == 3L, ] <-
    cbind(FALSE, grouped_term)[record_grouping == 3L, , drop = FALSE]
  group_inner <- cbind(FALSE, grouped_term)
  new_kind <- group_record | location_group_record
  record_terms[new_kind, ] <- FALSE
  candidates <- unique(rbind(record_terms[!nested_record & !new_kind, ,
                                          drop = FALSE],
                             outer_terms[nested_record, , drop = FALSE],
                             inner_terms[nested_record, , drop = FALSE],
                             group_terms[group_record, , drop = FALSE],
                             location_outer[location_group_record, ,
                                            drop = FALSE],
                             group_inner[location_group_record, ,
                                         drop = FALSE]))
  candidates <- candidates[rowSums(candidates) > 0, , drop = FALSE]
  class_key <- function(rows) apply(rows, 1L, paste, collapse = "")
  candidate_key <- class_key(candidates)
  class_of <- function(rows) match(class_key(rows), candidate_key)
  kind <- ifelse(nested_record, 2L,
                 ifelse(rowSums(record_terms) > 0, 1L, 0L))
  kind[group_record] <- 3L
  kind[location_group_record] <- 4L
  a <- b <- rep(1L, n_records)
  a[kind == 1L] <- class_of(record_terms[kind == 1L, , drop = FALSE])
  a[kind == 2L] <- class_of(outer_terms[kind == 2L, , drop = FALSE])
  b[kind == 2L] <- class_of(inner_terms[kind == 2L, , drop = FALSE])
  a[kind == 3L] <- class_of(group_terms[kind == 3L, , drop = FALSE])
  a[kind == 4L] <- class_of(location_outer[kind == 4L, , drop = FALSE])
  b[kind == 4L] <- class_of(group_inner[kind == 4L, , drop = FALSE])
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
       kind = kind, a = a, b = b, patterns = patterns, groups = groups,
       classes = classes,
       latent_columns = latent_columns, sampled_levels = sampled_levels,
       nugget_integrated = nugget_integrated,
       inner = lapply(seq_len(n_re), function(j) {
         which((nested_record & single[record_unit, j]) |
                 (location_group_record & grouped_term[, j]))
       }),
       counts = counts,
       integrated_terms = terms[counts$integrated > 0],
       any = any(counts$integrated > 0))
}

##' The table of integrated records for the records of an engine
##'
##' One compiled table object (`cpp_site_table_create()`) for the distinct
##' pairs (y, m) of the records, of the units of the nested patterns and of
##' those units summed, and of the members of groups, with `units` describing
##' each record for `cpp_site_evaluate()`: its `kind` (0 for the plain kernel,
##' 1 for a record, 2 for a nested pattern, 3 for a group and 4 for a location
##' of groups), its pair, pattern, group or location, and the classes `a` and
##' `b` of its variances. By default every record is one of kind 1 in the first
##' class, or of kind 0 where `plain`. `patterns` gives the units (`y` and
##' `m`) of each nested record, and identical patterns share one table.
##' `groups` gives the groups of each group and location of groups
##' (`plan_integrated_records()`), whose members' linear predictors relative
##' to the record's are `relative$D` times the coefficients plus
##' `relative$offset`, one row per member (`relative_predictors()`), with the
##' members of each location's groups repeated at the end as the group of all
##' its members. `relative$within` and `relative$within_offset` give them
##' relative to the first member of their group, from the differences of
##' their own rows, so that a covariate shared by a group's members adds
##' exactly zero. Adaptive Gauss-Hermite quadrature uses 20 nodes. `family` is
##' the family object of the model's link (`make_family()`), and a custom link
##' is passed to the table as the arrays of its kernel (`kernel_arrays()`,
##' `compute_custom_kernel_arrays()`).
##' @noRd
make_site_table <- function(y, units_m, family, plain = rep(FALSE, length(y)),
                            kind = ifelse(plain, 0L, 1L),
                            a = rep(1L, length(y)), b = rep(1L, length(y)),
                            patterns = NULL, groups = NULL) {
  nested <- which(kind == 2L)
  pattern_y <- lapply(patterns[nested], `[[`, "y")
  pattern_m <- lapply(patterns[nested], `[[`, "m")
  # The groups of the records of kinds 3 and 4, and for each location of
  # groups its groups and the group of all its members, added after the
  # others.
  grouped <- which(kind %in% c(3L, 4L))
  group_list <- unlist(groups[grouped], recursive = FALSE)
  group_count <- lengths(groups[grouped])
  located <- kind[grouped] == 4L
  location_groups <- split(seq_along(group_list),
                           rep(seq_along(grouped), group_count))[located]
  whole_list <- lapply(location_groups, function(members) {
    parts <- group_list[members]
    list(y = unlist(lapply(parts, `[[`, "y")),
         m = unlist(lapply(parts, `[[`, "m")),
         D = do.call(rbind, lapply(parts, `[[`, "D")),
         offset = unlist(lapply(parts, `[[`, "offset")))
  })
  whole <- length(group_list) + seq_along(whole_list)
  group_list <- c(group_list, whole_list)
  member_y <- unlist(lapply(group_list, `[[`, "y"))
  member_m <- unlist(lapply(group_list, `[[`, "m"))
  all_y <- c(y, unlist(pattern_y), vapply(pattern_y, sum, numeric(1)),
             member_y)
  all_m <- c(units_m, unlist(pattern_m), vapply(pattern_m, sum, numeric(1)),
             member_m)
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
  group_pairs <- lapply(group_list, function(g) pair_of(g$y, g$m))
  if (length(grouped)) {
    index[grouped[!located]] <- cumsum(c(1L, group_count))[
      seq_along(grouped)][!located]
    index[grouped[located]] <- seq_len(sum(located))
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
    as.integer(summed), group_pairs, unname(location_groups),
    as.integer(whole)
  )
  units <- cbind(kind = as.integer(kind), index = as.integer(index),
                 a = as.integer(a), b = as.integer(b))
  relative <- if (length(group_list)) {
    within <- function(rows) sweep(rows, 2L, rows[1L, ])
    list(D = do.call(rbind, lapply(group_list, `[[`, "D")),
         offset = unlist(lapply(group_list, `[[`, "offset")),
         within = do.call(rbind, lapply(group_list, function(g) within(g$D))),
         within_offset = unlist(lapply(group_list, function(g) {
           g$offset - g$offset[1L]
         })))
  }
  list(pointer = pointer, units = units, relative = relative)
}

##' The relative linear predictors of the members of a site table's groups
##'
##' At the regression coefficients `beta`, in the order of the members
##' (`make_site_table()`), each member's relative to its record and then each
##' member's relative to the first member of its group, as
##' `cpp_site_set_relative()` takes them, or NULL for a table without groups.
##' @noRd
relative_predictors <- function(sites, beta) {
  if (is.null(sites$relative)) return(NULL)
  c(as.numeric(sites$relative$D %*% beta + sites$relative$offset),
    as.numeric(sites$relative$within %*% beta +
                 sites$relative$within_offset))
}

##' The pointer of a site table, with its groups at the coefficients `beta`
##'
##' Sets the relative linear predictors of the groups' members
##' (`cpp_site_set_relative()`), which the compiled code that takes the
##' pointer then uses, and returns the pointer. A table without groups is
##' returned unchanged.
##' @noRd
site_pointer_at <- function(sites, beta) {
  relative <- relative_predictors(sites, beta)
  if (!is.null(relative)) cpp_site_set_relative(sites$pointer, relative)
  sites$pointer
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
##' takes the plain kernel. A table with groups (`make_site_table()`) takes
##' the regression coefficients `beta`, at which each method first sets the
##' relative linear predictors of the groups' members (`site_pointer_at()`).
##' @noRd
make_integrated_kernel <- function(sites, classes, beta = NULL) {
  units <- sites$units
  if (!is.null(sites$relative) &&
      length(beta) != ncol(sites$relative$D)) {
    riskmap_bug("The kernel of a table with groups needs the regression coefficients.")
  }
  relative <- relative_predictors(sites, beta)
  pointer_at <- if (is.null(relative)) {
    function() sites$pointer
  } else {
    function() {
      cpp_site_set_relative(sites$pointer, relative)
      sites$pointer
    }
  }
  pointer <- sites$pointer
  terms <- function(eta) {
    pointer <- pointer_at()
    values <- cpp_site_evaluate(pointer, classes, units,
                                matrix(eta, ncol = 1L), TRUE, 2L)
    list(log_likelihood = sum(values$h), score = as.numeric(values$h1),
         curvature = -as.numeric(values$h2))
  }
  nonlinear_parts <- function(centre, step) {
    pointer <- pointer_at()
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
         cpp_site_row_sums(pointer_at(), classes, units, eta)
       })
}

##' A smooth floor at zero for the curvature of the global proposal's factors
##'
##' epsilon log(1 + exp(C / epsilon)), which equals C to within
##' epsilon exp(-C / epsilon) for C well above epsilon and is never negative
##' (positive in exact arithmetic, zero in double precision for C below about
##' -745 epsilon). A floor with continuous derivatives keeps the fixed refits a
##' smooth function of the parameters wherever a record's log-likelihood is
##' not concave, for integrated sites, where that cannot be ruled out
##' numerically, and, moved to the Gershgorin bound, for links not known to be
##' log-concave (`make_curvature_floor()`, `choose_curvature_floor()`).
##' @noRd
smooth_curvature_floor <- function(curvature, epsilon = 1e-3) {
  pmax(curvature, 0) + epsilon * log1p(exp(-abs(curvature) / epsilon))
}

##' Whether a link object is that of integrated sites
##' @noRd
is_integrated_link <- function(linkf) identical(linkf$name, "integrated")
