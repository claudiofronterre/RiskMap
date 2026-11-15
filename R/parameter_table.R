##' The parameters of a model, one row for each
##'
##' The one description of a model's parameters, built once from its design,
##' its nugget, its family and what it holds. Every function that names,
##' orders, converts, bounds or starts the parameters reads it, from the
##' starts of the fitting engines to `coef()`, `summary()`, `print()`, the
##' profiles and the bootstrap. The positions come from `cov_par_layout()`,
##' which the engines' internals read. The table is not stored in a fit.
##' `fit_parameter_table()` rebuilds it from a fit's fields.
##'
##' `rows` has a row for each element of the working vector, in its order,
##' and then one for each parameter held at a value, whose `index` is NA.
##'
##' - `name`, the parameter as users write it, `beta`, `sigma2`, `phi`,
##'   `tau2`, a parameter of the family such as `sigma2_me` or `theta`, or
##'   `sigma2_re`.
##' - `element`, the column of the design for a coefficient and the term for
##'   a random-effect variance, NA otherwise.
##' - `stored` and `reported`, the names of the parameter in
##'   `unlist(fit$estimate)` and in `unlist(coef(fit))`. They differ for the
##'   nugget, stored as `nu2`, the logarithm of its ratio to the spatial
##'   variance, and reported as `tau2`, its variance. NA for a held parameter.
##' - `block`, `coefficients`, `covariance`, `family` or `random_effects`.
##' - `index`, the position in the working vector.
##' - `scale`, the map between the natural and the working scale, `identity`,
##'   `log`, `nugget` (`nugget_working()` and `nugget_natural()`, which
##'   depend on the spatial variance) or `family`, the family's own map
##'   (`family_parameter()`).
##' - `lower` and `upper`, the bounds of the working value on the scale of the
##'   fitting engines, whose distances are divided by the distance scale
##'   (`scale_spatial_coordinates()`), and `shift`, whether the bounds move by
##'   the logarithm of the distance scale on the scale of the model's units
##'   (`parameter_bounds()`).
##' - `limit`, `variance` for a variance that can be held at zero
##'   (`variance_index()`), `family` for a parameter of the family that can
##'   be held at its limit (`limit_index()`), `""` otherwise.
##' - `label`, the label of `summary()` (`parameter_labels()`), and
##'   `working_label`, that of the working parameter where it differs, the
##'   ratio of the nugget to the spatial variance.
##' - `start`, the rule of the start when `start_pars` does not give one.
##'   `regression` for the coefficients of the regression without the latent
##'   effects, `family_regression` for a parameter of the family taken from
##'   the family's regression, `one`, and `distance_quantile`, the 10th
##'   percentile of the distances between locations (`default_start()`).
##' - `held`, the value of a held parameter, NA for one that is estimated. A
##'   model without a nugget holds `tau2` at zero.
##'
##' @param beta_names The names of the columns of the design, and `p` their
##'   number.
##' @param fix_tau2 `TRUE` when the nugget is estimated, its variance when it
##'   is held, 0 without a nugget (`gp()`).
##' @param estimated_family The family's parameters that are estimated, in
##'   the order of `riskmap_families`.
##' @param held_family The family's parameters held at given values, a named
##'   list (`held_family_parameters()`).
##' @param re_names The names of the random effects, and `n_re` their number.
##' @param distance_units The units of distance, for the label of the scale.
##' @return A list of class `riskmap_parameter_table` with `rows`, `layout`
##'   (from `cov_par_layout()`), `fix_tau2`, and the names `beta_names` and
##'   `re_names` as given, which name the estimate.
##' @noRd
parameter_table <- function(beta_names, fix_tau2,
                            estimated_family = character(),
                            held_family = list(), re_names = NULL,
                            p = length(beta_names), n_re = length(re_names),
                            distance_units = NULL) {
  # NULL, as estimated_family_parameters() gives for a family without
  # parameters, is no parameter, as character() is.
  estimated_family <- as.character(estimated_family)
  layout <- cov_par_layout(p, fix_tau2 = fix_tau2,
                           family_parameters = estimated_family, n_re = n_re)
  nugget <- isTRUE(fix_tau2)
  # The names that unlist() gives the estimate, built the same way.
  template <- c(
    list(beta = stats::setNames(numeric(p), beta_names), sigma2 = 0, phi = 0),
    if (nugget) list(nu2 = 0),
    stats::setNames(as.list(numeric(length(estimated_family))),
                    estimated_family),
    if (n_re > 0) list(sigma2_re = stats::setNames(numeric(n_re), re_names))
  )
  stored <- names(unlist(template))
  family_start <- function(name) {
    family_parameter(name)$start %||% NA_character_
  }
  piece <- function(name, element, block, scale, lower, upper,
                    shift = FALSE, limit = "", start = "one") {
    n <- length(element)
    data.frame(name = rep(name, n), element = element, block = rep(block, n),
               scale = rep(scale, n), lower = rep(lower, n),
               upper = rep(upper, n), shift = rep(shift, n),
               limit = rep(limit, n), start = rep(start, n),
               stringsAsFactors = FALSE)
  }
  estimated <- rbind(
    piece("beta", if (p) beta_names %||% rep(NA_character_, p) else
      character(), "coefficients", "identity", -30, 30, start = "regression"),
    piece("sigma2", NA_character_, "covariance", "log", -20, 20,
          limit = "variance"),
    piece("phi", NA_character_, "covariance", "log", -20, 10, shift = TRUE,
          start = "distance_quantile"),
    if (nugget) {
      piece("tau2", NA_character_, "covariance", "nugget", -20, 20,
            limit = "variance")
    },
    do.call(rbind, lapply(estimated_family, function(name) {
      piece(name, NA_character_, "family", "family", -20, 20,
            limit = "family", start = family_start(name))
    })),
    if (n_re > 0) {
      piece("sigma2_re", re_names %||% rep(NA_character_, n_re),
            "random_effects", "log", -20, 20, limit = "variance")
    }
  )
  estimated$stored <- stored
  estimated$reported <- ifelse(stored == "nu2" & estimated$name == "tau2",
                               "tau2", stored)
  estimated$index <- unname(unlist(layout$index[c(
    "beta", "sigma2", "phi", if (nugget) "nu2", estimated_family,
    if (n_re > 0) "sigma2_re"
  )]))
  estimated$held <- NA_real_

  held_nugget <- is.numeric(fix_tau2) && !nugget
  held <- rbind(
    if (held_nugget) {
      piece("tau2", NA_character_, "covariance", "nugget", NA_real_,
            NA_real_)
    },
    do.call(rbind, lapply(names(held_family), function(name) {
      piece(name, NA_character_, "family", "family", NA_real_, NA_real_,
            start = family_start(name))
    }))
  )
  if (!is.null(held)) {
    held$stored <- NA_character_
    held$reported <- NA_character_
    held$index <- NA_integer_
    held$held <- unname(c(if (held_nugget) as.numeric(fix_tau2),
                          vapply(held_family, as.numeric, numeric(1))))
  }
  rows <- rbind(estimated, held)
  # Labels by the names that coef() reports, and by the parameter's own name
  # for one that is held.
  rows$label <- unname(parameter_labels(
    ifelse(is.na(rows$reported), rows$name, rows$reported), distance_units
  ))
  rows$working_label <- ifelse(rows$scale == "nugget",
                               "Nugget variance ratio", rows$label)
  rows <- rows[c("name", "element", "stored", "reported", "block", "index",
                 "scale", "lower", "upper", "shift", "limit", "label",
                 "working_label", "start", "held")]
  rownames(rows) <- NULL
  structure(list(rows = rows, layout = layout, fix_tau2 = fix_tau2,
                 beta_names = beta_names, re_names = re_names),
            class = "riskmap_parameter_table")
}

##' The parameter table of a model as `glgpm()` specifies it
##'
##' The family's block is the measurement error variance of a Gaussian model
##' unless `fix_var_me` holds it (`gaussian_family_block()`), and otherwise
##' the parameters the family estimates. The held ones are those of
##' `held_family_parameters()`.
##' @noRd
model_parameter_table <- function(beta_names, fix_tau2, exact, fix_var_me,
                                  estimated_family, family_parameters,
                                  re_names, distance_units) {
  held <- family_parameters %||% list()
  if (is.numeric(fix_var_me)) held$sigma2_me <- fix_var_me
  parameter_table(
    beta_names, fix_tau2 = fix_tau2,
    estimated_family = if (exact) gaussian_family_block(fix_var_me) else
      estimated_family,
    held_family = held, re_names = re_names, distance_units = distance_units
  )
}

##' The parameter table of a fitted model
##'
##' Rebuilt from the fit's fields, the names of its design, its estimate and
##' its holds. Whether the nugget is estimated is read from the estimate, as
##' `coef()` has always read it, so that a fit returned by
##' `fit_glgpm_design()` before `glgpm()` adds its settings has its table too.
##' @noRd
fit_parameter_table <- function(fit) {
  estimate <- fit$estimate
  parameter_table(
    beta_names = colnames(as.matrix(fit$D)),
    fix_tau2 = if (!is.null(estimate$nu2)) TRUE else fit$fix_tau2,
    estimated_family = estimated_family_parameters(fit),
    held_family = held_family_parameters(fit),
    re_names = names(estimate$sigma2_re),
    n_re = length(estimate$sigma2_re),
    distance_units = fit$distance_units
  )
}

##' The rows of the parameters in the working vector
##' @noRd
estimated_rows <- function(table) {
  table$rows[!is.na(table$rows$index), , drop = FALSE]
}

##' The names of a block's estimated parameters, as the estimate names them
##'
##' One name for each parameter, `beta` and `sigma2_re` once each, in the
##' order of the working vector.
##' @noRd
estimate_names <- function(table) {
  rows <- estimated_rows(table)
  unique(ifelse(rows$scale == "nugget", "nu2", rows$name))
}

##' Whether a model estimates its nugget
##' @noRd
nugget_estimated <- function(table) {
  any(table$rows$scale == "nugget" & !is.na(table$rows$index))
}

##' The held value of a parameter, NULL when it is estimated or absent
##' @noRd
held_value <- function(table, name) {
  rows <- table$rows
  value <- rows$held[rows$name == name & is.na(rows$index)]
  if (length(value)) value
}

##' The working value of the nugget, the logarithm of its ratio to the
##' spatial variance, from both variances
##' @noRd
nugget_working <- function(tau2, sigma2) log(tau2 / sigma2)

##' The nugget variance from its working value and the working value of the
##' spatial variance, both on the log scale
##' @noRd
nugget_natural <- function(nu2, sigma2) exp(nu2 + sigma2)

##' The ratio of the nugget to the spatial variance from both variances,
##' with which the correlation matrix of the locations is augmented
##' @noRd
nugget_ratio_natural <- function(tau2, sigma2) tau2 / sigma2

##' The nugget variance of a fit, whether it is estimated or held
##'
##' From the fit's estimate through `nugget_natural()`, and NULL for a fit
##' that neither holds nor estimates one, as `fix_tau2 = FALSE` of an
##' earlier version. Reads only `fix_tau2` and the estimate, so it works on
##' a fit's essentials.
##' @noRd
fit_nugget_variance <- function(fit) {
  if (is.numeric(fit$fix_tau2)) {
    fit$fix_tau2
  } else if (!is.null(fit$estimate$nu2)) {
    nugget_natural(fit$estimate$nu2, fit$estimate$sigma2)
  }
}

##' Starting values on the working scale of the fitting engines
##'
##' `start` holds the parameters on their natural scale, named as users name
##' them (`tau2` for the nugget). Coefficients are kept, variances and the
##' scale go to the log scale, the nugget to the logarithm of its ratio to
##' the spatial variance (`nugget_working()`), and the family's parameters
##' through their own maps. The MCML engine and both Gaussian paths start
##' from it, so the start of each parameter is on the scale on which its
##' likelihood reads it.
##' @noRd
working_start <- function(table, start) {
  out <- start_value(table, start, "beta")
  for (name in setdiff(estimate_names(table), "beta")) {
    value <- start_value(table, start, name)
    out <- c(out, switch(
      name,
      sigma2 = , phi = , sigma2_re = log(value),
      nu2 = nugget_working(value, start$sigma2),
      family_parameter(name)$working(value)
    ))
  }
  unname(out)
}

##' The start of one parameter, checked to have a value for each of its rows
##'
##' `name` is as the estimate names it, so `nu2` reads the start of `tau2`. A
##' start that is missing or of the wrong length would shift every later
##' parameter in the working vector, so it stops as a bug.
##' @noRd
start_value <- function(table, start, name) {
  user <- if (name == "nu2") "tau2" else name
  value <- start[[user]]
  wanted <- sum(estimated_rows(table)$name == user)
  if (length(value) != wanted) {
    riskmap_bug("The start of {.code {user}} has {length(value)} values, not {wanted}.")
  }
  value
}

##' Split the working vector of the fitting engines into the estimate
##'
##' A named list on the working scale, `beta` named by the columns of the
##' design, then `sigma2`, `phi`, `nu2` when the nugget is estimated, the
##' family's estimated parameters and `sigma2_re` named by its terms.
##'
##' A flat named vector (as used prior to #92) can't safely be keyed by
##' parameter name: a covariate literally named e.g. `"sigma2"` collides with
##' the spatial variance parameter's own name, silently corrupting name-based
##' lookups. Structuring `estimate` as a list instead keeps those namespaces
##' separate. `summary.RiskMap()` additionally needs a *flat* vector to line
##' up against `covariance` (which is estimated as one joint matrix over
##' every parameter) - `unlist()` provides that, and disambiguates the same
##' way (e.g. `"beta.sigma2"` vs `"sigma2"`), with the names in `stored`.
##' @noRd
working_estimate <- function(table, par) {
  index <- table$layout$index
  beta <- par[index$beta]
  names(beta) <- table$beta_names
  out <- list(beta = beta,
              sigma2 = unname(par[index$sigma2]),
              phi    = unname(par[index$phi]))
  if (!is.null(index$nu2)) out$nu2 <- unname(par[index$nu2])
  for (name in table$layout$family) out[[name]] <- unname(par[index[[name]]])
  if (!is.null(index$sigma2_re)) {
    sigma2_re <- unname(par[index$sigma2_re])
    names(sigma2_re) <- table$re_names
    out$sigma2_re <- sigma2_re
  }
  out
}

##' The estimate on the natural scale of each parameter
##'
##' `estimate` is on the working scale, as a fit stores it. The coefficients
##' are named by the columns of the design, the variances and the scale are
##' `exp()` of their working values, the nugget is its variance
##' (`nugget_natural()`) and the family's parameters go through their own
##' maps. The list is in the order of the working vector, with `tau2` in the
##' place of `nu2`. With `held`, the held parameters are added at their
##' values, after the estimated ones.
##' @noRd
natural_values <- function(table, estimate, held = FALSE) {
  rows <- table$rows
  out <- list(beta = estimate$beta)
  names(out$beta) <- table$beta_names
  out$sigma2 <- exp(estimate$sigma2)
  out$phi <- exp(estimate$phi)
  for (name in setdiff(estimate_names(table), c("beta", "sigma2", "phi"))) {
    if (name == "nu2") {
      out$tau2 <- nugget_natural(estimate$nu2, estimate$sigma2)
    } else if (name == "sigma2_re") {
      out$sigma2_re <- exp(estimate$sigma2_re)
    } else {
      out[[name]] <- family_parameter(name)$natural(estimate[[name]])
    }
  }
  if (held) {
    fixed <- rows[is.na(rows$index), , drop = FALSE]
    for (k in seq_len(nrow(fixed))) out[[fixed$name[k]]] <- fixed$held[k]
  }
  out
}

##' The order of `coef()`
##'
##' The coefficients, the spatial variance and scale, the family's
##' parameters, then the nugget and the random-effect variances. It differs
##' from the working vector's, where the nugget comes before the family's
##' parameters, and changing it would change what `coef()` returns.
##' @noRd
coef_order <- function(table) {
  c("beta", "sigma2", "phi", table$layout$family,
    if (nugget_estimated(table)) "tau2",
    if (!is.null(table$layout$index$sigma2_re)) "sigma2_re")
}

##' The Jacobian of the reported parameters in the stored working ones
##'
##' On the scale on which `summary()` forms its intervals, every parameter is
##' its working value except the nugget, whose logarithm is the sum of the
##' working values of the nugget and the spatial variance. Rows and columns
##' are named by `names`, the flat names of the estimate, `stored` by
##' default, and the nugget's entry is placed by name, so that it lands on
##' the right row whatever the order of the estimate.
##' @noRd
parameter_jacobian <- function(table, names = estimated_rows(table)$stored) {
  J <- diag(length(names))
  dimnames(J) <- list(names, names)
  rows <- estimated_rows(table)
  nugget <- intersect(rows$stored[rows$scale == "nugget"], names)
  if (length(nugget)) J[nugget, "sigma2"] <- 1
  J
}

##' The bounds of the working vector
##'
##' On the scale of the fitting engines by default, and with the logarithm
##' of `distance_scale` added to those of the scale on the scale of the
##' model's distance units, as the profiles see them. `coefficients` replaces
##' the bounds of the coefficients, for the Gaussian profile that solves them
##' in closed form.
##' @return A list with `lower` and `upper`.
##' @noRd
parameter_bounds <- function(table, log_scale = NULL, coefficients = NULL) {
  rows <- estimated_rows(table)
  lower <- rows$lower
  upper <- rows$upper
  if (!is.null(coefficients)) {
    lower[rows$block == "coefficients"] <- coefficients[1L]
    upper[rows$block == "coefficients"] <- coefficients[2L]
  }
  if (!is.null(log_scale)) {
    lower[rows$shift] <- lower[rows$shift] + log_scale
    upper[rows$shift] <- upper[rows$shift] + log_scale
  }
  list(lower = lower, upper = upper)
}

##' The positions of the variances that can be held at zero
##'
##' The spatial variance, the nugget ratio and the random-effect variances.
##' @noRd
variance_index <- function(table) {
  rows <- estimated_rows(table)
  rows$index[rows$limit == "variance"]
}

##' The positions of the family's parameters that can be held at a limit
##'
##' NULL when the family estimates none.
##' @noRd
limit_index <- function(table) {
  rows <- estimated_rows(table)
  index <- rows$index[rows$limit == "family"]
  if (length(index)) index
}

##' The labels of parameters, by the names of the estimate or of `coef()`
##'
##' `names` are flat names, of `unlist(coef(fit))` or of `unlist(fit$estimate)`.
##' The nugget has the label of its variance as `tau2` and that of its ratio
##' to the spatial variance as `nu2`.
##' @noRd
table_labels <- function(table, names) {
  rows <- estimated_rows(table)
  labels <- c(stats::setNames(rows$label, rows$reported),
              stats::setNames(rows$working_label, rows$stored))
  stats::setNames(unname(labels[names]), names)
}

##' The block of each parameter, by its flat name
##'
##' For objects that keep the names of the parameters but not the fit, such
##' as a bootstrap, from the names the table gives them. `coefficients` for
##' `beta.` names, `family` for the family's parameters, `random_effects` for
##' `sigma2_re.` names and `covariance` for the others.
##' @noRd
parameter_block <- function(names) {
  block <- rep("covariance", length(names))
  block[names %in% family_parameter_names()] <- "family"
  block[grepl("^sigma2_re", names)] <- "random_effects"
  block[grepl("^beta\\.", names)] <- "coefficients"
  block
}

##' The start of a parameter by its rule, when `start_pars` does not give it
##'
##' The rules `one` and `distance_quantile` of the table. The coefficients and
##' the family's parameters start from regressions on the model's records,
##' which `glgpm()` fits.
##' @noRd
default_start <- function(table, name, coords = NULL) {
  rows <- table$rows[table$rows$name == name, , drop = FALSE]
  rule <- rows$start[1L]
  switch(
    rule,
    one = rep(1, nrow(rows)),
    distance_quantile = quantile(pairwise_distances(coords), 0.1),
    riskmap_bug("The start of {.code {name}} comes from the model's regression.")
  )
}

##' Read the starting values given in `start_pars`
##'
##' `start_pars` is NULL or a named list over any of the model's estimated
##' parameters, named as users name them (`beta`, `sigma2`, `phi`, `tau2`,
##' `sigma2_me`, `theta`, `sigma2_re`) or by the flat names of
##' `unlist(coef(fit))`, such as `"beta.x"` and `"sigma2_re.household"`, so
##' that `start_pars = coef(fit)` and a name copied from any output both work.
##' `beta` is matched by the names of the design's columns, or by position
##' when it is unnamed and has a value for every column. `sigma2_re` is
##' matched by the names of the `re()` terms, or by position when it is
##' unnamed and has a value for every term. Every value must be finite, and
##' those of the variances, the scale and the family's parameters positive.
##' A name the model lacks stops with the names it has, and a start of `phi`
##' below a hundredth of the smallest distance between the locations or above
##' ten times the largest gives a warning.
##'
##' @param table The model's parameter table.
##' @param coords The unique coordinates, in `distance_units`, for the check
##'   of `phi`.
##' @return A named list over the parameters given, `beta` and `sigma2_re` as
##'   named vectors of the elements given, the others as single values.
##' @noRd
read_start_pars <- function(table, start_pars, coords = NULL,
                            distance_units = NULL,
                            call = rlang::caller_env()) {
  if (is.null(start_pars)) return(list())
  keys <- names(start_pars)
  if (!is.list(start_pars) || (length(start_pars) &&
                               (is.null(keys) || any(!nzchar(keys))))) {
    supplied <- describe_value(start_pars)
    abort_input("{.arg start_pars} must be a named list of starting values, not {supplied}.", arg = "start_pars", call = call)
  }
  rows <- estimated_rows(table)
  held <- table$rows[is.na(table$rows$index), , drop = FALSE]
  model_names <- unique(rows$name)
  vectors <- c(beta = "coefficients", sigma2_re = "random_effects")
  given <- list()
  add <- function(name, element, value, key) {
    old <- given[[name]]
    if (name %in% names(vectors)) {
      if (any(element %in% names(old))) {
        abort_input("{.arg start_pars} gives the start of {.val {paste0(name, '.', intersect(element, names(old)))}} twice.", arg = "start_pars", call = call)
      }
      given[[name]] <<- c(old, stats::setNames(value, element))
    } else {
      if (!is.null(old)) {
        abort_input("{.arg start_pars} gives the start of {.code {name}} twice.", arg = "start_pars", call = call)
      }
      given[[name]] <<- value
    }
  }
  check_values <- function(value, key, positive) {
    bad <- !is.numeric(value) || !length(value) || any(!is.finite(value)) ||
      (positive && any(value <= 0))
    if (bad) {
      supplied <- describe_value(value)
      what <- if (positive) "positive finite numbers" else "finite numbers"
      abort_input("The starting values in {.code start_pars${key}} must be {what}, not {supplied}.", arg = "start_pars", call = call)
    }
  }
  given_keys <- keys[!vapply(start_pars, is.null, logical(1))]
  unknown <- setdiff(given_keys, c(model_names, rows$reported, held$name))
  if (length(unknown)) {
    abort_input(c("{.arg start_pars} must name parameters of the model.", "x" = "{.val {unknown}} {?is/are} not among them.", "i" = "The model's parameters are {.val {model_names}}.", "i" = "{.fn glgpm_parameters} lists them with their elements."), arg = "start_pars", call = call)
  }
  for (key in keys) {
    value <- start_pars[[key]]
    if (is.null(value)) next
    if (key %in% model_names && key %in% names(vectors)) {
      elements <- rows$element[rows$name == key]
      check_values(value, key, positive = key != "beta")
      supplied <- names(value)
      if (is.null(supplied)) {
        if (length(value) != length(elements)) {
          what <- if (key == "beta") "column of the design matrix" else
            "`re()` term"
          abort_input(c("An unnamed {.code start_pars${key}} must have a value for each {what}, {length(elements)}, not {length(value)}.", "i" = "Name its values to give some of them, among {.val {elements}}."), arg = "start_pars", call = call)
        }
        supplied <- elements
      } else if (any(!nzchar(supplied)) || any(!supplied %in% elements)) {
        unknown <- setdiff(supplied, elements)
        abort_input(c("{.code start_pars${key}} must be named by {if (key == 'beta') 'the columns of the design matrix' else 'the `re()` terms'}.", "x" = "{.val {unknown}} {?is/are} not among them.", "i" = "They are {.val {elements}}."), arg = "start_pars", call = call)
      }
      add(key, supplied, unname(value), key)
    } else if (key %in% model_names) {
      check_values(value, key, positive = TRUE)
      if (length(value) != 1L) {
        abort_input("{.code start_pars${key}} must be a single number, not {length(value)} numbers.", arg = "start_pars", call = call)
      }
      add(key, NULL, unname(value), key)
    } else if (key %in% rows$reported) {
      row <- rows[rows$reported == key, , drop = FALSE]
      check_values(value, key, positive = row$block != "coefficients")
      if (length(value) != 1L) {
        abort_input("{.code start_pars${key}} must be a single number, not {length(value)} numbers.", arg = "start_pars", call = call)
      }
      add(row$name, if (row$name %in% names(vectors)) row$element, unname(value), key)
    } else {
      value <- held$held[held$name == key]
      if (key == "tau2" && identical(value, 0)) {
        abort_input("A starting value for {.code tau2} can't be given when {.fn gp} has {.code nugget = FALSE}.", arg = "start_pars", call = call)
      }
      abort_input(c("{.code start_pars} can't give a start for {.code {key}}, which the model holds at {value}.", "i" = "Drop it from {.arg start_pars}."), arg = "start_pars", call = call)
    }
  }
  # Values of beta and sigma2_re in the order of the design and the terms.
  for (name in intersect(names(vectors), names(given))) {
    elements <- rows$element[rows$name == name]
    given[[name]] <- given[[name]][intersect(elements, names(given[[name]]))]
  }
  if (!is.null(given$phi) && !is.null(coords) && nrow(coords) > 1L) {
    distances <- as.numeric(pairwise_distances(coords))
    distances <- distances[distances > 0]
    smallest <- min(distances)
    largest <- max(distances)
    if (given$phi < smallest / 100 || given$phi > 10 * largest) {
      units <- distance_units %||% "the model's units"
      riskmap_warn(c("The start of {.code phi}, {signif(given$phi, 3)} {units}, is far from the distances between the locations.", "i" = "They range from {signif(smallest, 3)} to {signif(largest, 3)} {units}, so {.code phi} may be in other units."))
    }
  }
  given
}

##' Put the starting values given in `start_pars` in place of others
##'
##' `start` is a complete start, as `glgpm()` builds it, and `given` what
##' `read_start_pars()` read. A coefficient or random-effect variance given
##' replaces its element, found by name.
##' @noRd
apply_given_start <- function(table, start, given) {
  for (name in names(given)) {
    value <- given[[name]]
    if (name == "beta") {
      start$beta[match(names(value), table$beta_names)] <- unname(value)
    } else if (name == "sigma2_re") {
      start$sigma2_re[match(names(value), table$re_names)] <- unname(value)
    } else {
      start[[name]] <- value
    }
  }
  start
}

##' Whether `start_pars` gives every estimated parameter
##' @noRd
all_starts_given <- function(table, given) {
  rows <- estimated_rows(table)
  vectors <- rows$name %in% c("beta", "sigma2_re")
  scalars <- unique(rows$name[!vectors])
  all(scalars %in% names(given)) &&
    all(rows$element[rows$name == "beta"] %in% names(given$beta)) &&
    all(rows$element[rows$name == "sigma2_re"] %in% names(given$sigma2_re))
}

##' The flat names of the parameters given in `start_pars`, as
##' `unlist(coef(fit))` names them
##' @noRd
given_start_names <- function(given) {
  if (!length(given)) return(NULL)
  names(unlist(given))
}

##' The parameters of a geostatistical model
##'
##' @name riskmap_parameters
##' @description
##' Every parameter of a model fitted by [glgpm()] has one name, which
##' `start_pars` of [glgpm()], the result of [coef.RiskMap()], `parameters`
##' of [specify_glgpm()] and the columns of [parametric_bootstrap()] share. The names a model has follow from its
##' formula and family, so they can be listed before it is fitted, with
##' [glgpm_parameters()].
##'
##' | Name | When the model has it | Elements and scale |
##' |---|---|---|
##' | `beta` | Always | One coefficient per column of the design matrix, named as the columns, such as `"(Intercept)"`, `"log(elevation)"` or `"typeurban"` for a level of a factor |
##' | `sigma2` | Always | The variance of the spatial process |
##' | `phi` | Always | The scale of the spatial correlation, in the model's `distance_units` |
##' | `tau2` | With `gp(nugget = TRUE)` | The variance of the nugget |
##' | `sigma2_me` | Gaussian family | The variance of the measurement error |
##' | `theta` | Negative binomial family | The size, see [negative_binomial()] |
##' | `sigma2_re` | With `re()` terms | One variance per grouping variable, named by the variable as written in `re()` |
##'
##' The flat names of `unlist(coef(fit))`, such as `"beta.log(elevation)"` and
##' `"sigma2_re.household"`, name a single element and are accepted wherever
##' the names are.
##'
##' @section Starting values:
##' `start_pars` takes any of a model's estimated parameters, one, several or
##' all of them. `beta` can be a named subset of the coefficients, or unnamed
##' with one value per column. `sigma2_re` can be a named subset of the
##' variances, or unnamed with one value per term in the order of the formula.
##' Every value must be finite, and every one but the coefficients positive.
##' The others take their automatic start, given as `start` by
##' [glgpm_parameters()]. `regression` is the coefficient of the regression
##' without the latent effects, `family_regression` the size of the negative
##' binomial regression without them, `one` the value 1 and
##' `distance_quantile` the 10th percentile of the distances between the
##' locations. A name the model lacks stops with the names it has.
##' @seealso [glgpm()], [glgpm_parameters()], [coef.RiskMap()]
NULL

##' List the parameters of a model before fitting it
##'
##' The parameters [glgpm()] estimates or holds for a model, with the names
##' that `start_pars` and [coef.RiskMap()] use ([riskmap_parameters]). Only the
##' formula, the data and the family are read, so this is quick.
##'
##' @param formula,data,family,fix_var_me As in [glgpm()].
##' @param distance_units The units of distance, `"km"` (the default) or `"m"`,
##'   for the label of the scale.
##' @return A data frame of class `riskmap_parameters`, with a row per element
##'   of each parameter and the columns `name`, `element` (the column of the
##'   design matrix of a coefficient, the term of a random-effect variance),
##'   `flat` (the name in `unlist(coef(fit))`), `label` (as in
##'   [summary.RiskMap()]), `start` (the rule of its automatic start, see
##'   [riskmap_parameters]) and `held` (the value at which the model holds the
##'   parameter, `NA` when it is estimated). Its print shows the template of a
##'   `start_pars` list.
##' @seealso [riskmap_parameters], [glgpm()]
##' @examples
##' data(liberia)
##' liberia <- sf::st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
##' glgpm_parameters(npos ~ log(elevation) + gp(nugget = TRUE), data = liberia,
##'                  family = "binomial")
##' @export
glgpm_parameters <- function(formula, data, family, fix_var_me = NULL,
                             distance_units = c("km", "m")) {
  check_supplied(formula)
  check_supplied(data)
  check_supplied(family)
  check_data(data)
  check_formula(formula, data)
  family_link <- read_family(family, NULL)
  model <- family_traits(family_link$family)
  distance_units <- match_choice(distance_units)
  if (!model$exact && !is.null(fix_var_me)) {
    abort_input("{.arg fix_var_me} can be given only for a Gaussian model, not a {family_link$family} one.", arg = "fix_var_me")
  }
  inter_f <- interpret_formula(formula)
  frame <- model.frame(inter_f$pf, data = data, na.action = na.pass)
  design <- model.matrix(attr(frame, "terms"), data = data)
  parameters <- model_parameter_table(
    colnames(design), fix_tau2 = inter_f$gp_spec$nugget, exact = model$exact,
    fix_var_me = fix_var_me,
    estimated_family = family_link$estimated %||% character(),
    family_parameters = family_link$parameters,
    re_names = if (length(inter_f$re_spec)) inter_f$re_spec$term,
    distance_units = distance_units
  )
  rows <- parameters$rows
  # A model without a nugget holds it at zero in the table, which is not a
  # parameter users hold.
  rows <- rows[!(rows$name == "tau2" & rows$held %in% 0), , drop = FALSE]
  out <- data.frame(
    name = rows$name, element = rows$element,
    flat = ifelse(is.na(rows$reported), rows$name, rows$reported),
    label = rows$label, start = rows$start, held = rows$held,
    stringsAsFactors = FALSE
  )
  class(out) <- c("riskmap_parameters", "data.frame")
  out
}

##' @export
print.riskmap_parameters <- function(x, ...) {
  rows <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  estimated <- rows[is.na(rows$held), , drop = FALSE]
  template <- vapply(unique(estimated$name), function(name) {
    elements <- estimated$element[estimated$name == name]
    if (name %in% c("beta", "sigma2_re")) {
      paste0(name, " = c(", paste0("\"", elements, "\" = ", collapse = ", "),
             ")")
    } else {
      paste0(name, " = ")
    }
  }, character(1))
  shown <- rows
  shown$element[is.na(shown$element)] <- ""
  shown$held <- ifelse(is.na(shown$held), "", format(shown$held))
  out_print(c(
    out_capture(print.data.frame(shown, row.names = FALSE)),
    "",
    out_text("A template for `start_pars`:"),
    paste0("list(", paste(template, collapse = ", "), ")")
  ))
  invisible(x)
}
