# The family of a model is read through its family object (R/family.R, #253),
# so a function outside R/family.R does not compare a family string, as in
# family == "binomial" or switch(family, binomial = ...), nor look one up in a
# table keyed by family, as in c(binomial = "a", poisson = "b")[[family]]. This
# test reads every function in the installed package and fails if one that is
# not the family's own does, unless it is listed in
# family-strings-to-convert.txt, which may only shrink. It reads the code, not
# what it means, so a "gaussian" that is not a family, such as a Gaussian
# correlation, would be counted too.

family_strings <- c("gaussian", "binomial", "poisson")
comparing_calls <- c("==", "!=", "%in%", "identical", "match.arg",
                     "is.element", "startsWith", "endsWith")

# The functions of R/family.R that compare family strings, where they belong.
family_functions <- c(
  "canonical_slope_parts", "compute_custom_kernel_arrays",
  "compute_custom_link_terms", "compute_largest_record_terms",
  "compute_link_log_likelihood_rows", "family_traits", "link_nonlinear_parts",
  "make_family", "make_inverse_link_functions", "move_start_off_limits",
  "read_family", "replace_at_canonical_limits"
)

is_family_string <- function(x) {
  is.character(x) && length(x) > 0L && any(x %in% family_strings)
}

holds_family_string <- function(x) {
  if (is_family_string(x)) return(TRUE)
  if (is.call(x) && identical(x[[1L]], as.name("c"))) {
    return(any(vapply(as.list(x)[-1L], is_family_string, logical(1))))
  }
  FALSE
}

# A vector or list of several family strings, or keyed by several of them,
# the choices of a comparison or a table of the families
offers_family_strings <- function(x) {
  if (!is.call(x) || !is.name(x[[1L]]) ||
      !as.character(x[[1L]]) %in% c("c", "list")) {
    return(FALSE)
  }
  parts <- as.list(x)[-1L]
  values <- sum(vapply(parts, is_family_string, logical(1)))
  keys <- sum(names(parts) %in% family_strings)
  values > 1L || keys > 1L
}

# The name of the function a call calls, also through `::` and `:::`
called_name <- function(head) {
  if (is.name(head)) return(as.character(head))
  if (is.call(head) && length(head) == 3L &&
      as.character(head[[1L]]) %in% c("::", ":::")) {
    return(as.character(head[[3L]]))
  }
  NA_character_
}

# Whether a function, call or list of calls compares a family string, through
# a comparison or match.arg() with a family string as an argument, several
# family strings as choices or as the keys of a table, or a switch() with a
# family string as the name of a case.
compares_family_string <- function(expr) {
  if (is.function(expr)) {
    defaults <- as.list(formals(expr))
    for (i in seq_along(defaults)) {
      if (identical(defaults[[i]], quote(expr = ))) next
      if (offers_family_strings(defaults[[i]])) return(TRUE)
    }
    return(compares_family_string(defaults) ||
             compares_family_string(body(expr)))
  }
  if (is.list(expr) || is.pairlist(expr)) {
    for (i in seq_along(expr)) {
      if (identical(expr[[i]], quote(expr = ))) next
      if (compares_family_string(expr[[i]])) return(TRUE)
    }
    return(FALSE)
  }
  if (!is.call(expr)) return(FALSE)
  if (offers_family_strings(expr)) return(TRUE)
  name <- called_name(expr[[1L]])
  if (!is.na(name)) {
    arguments <- as.list(expr)[-1L]
    if (name %in% comparing_calls) {
      for (i in seq_along(arguments)) {
        if (identical(arguments[[i]], quote(expr = ))) next
        if (holds_family_string(arguments[[i]])) return(TRUE)
      }
    }
    if (name == "switch" && any(names(arguments) %in% family_strings)) {
      return(TRUE)
    }
  }
  compares_family_string(as.list(expr))
}

functions_comparing_family_strings <- function() {
  namespace <- asNamespace("RiskMap")
  names <- ls(namespace, all.names = TRUE)
  compares <- vapply(names, function(name) {
    object <- get(name, envir = namespace)
    is.function(object) && !is.primitive(object) &&
      compares_family_string(object)
  }, logical(1))
  names[compares]
}

test_that("the check finds the ways a family string is compared", {
  expect_true(compares_family_string(function(family) family == "binomial"))
  expect_true(compares_family_string(function(x) x$family != "gaussian"))
  expect_true(compares_family_string(
    function(family) identical(family, "poisson")
  ))
  expect_true(compares_family_string(
    function(family) family %in% c("binomial", "poisson")
  ))
  expect_true(compares_family_string(
    function(family) switch(family, binomial = 1, 2)
  ))
  expect_true(compares_family_string(
    function(family = c("binomial", "poisson")) match.arg(family)
  ))
  expect_true(compares_family_string(
    function(family) function(y) if (family == "poisson") y
  ))
  expect_true(compares_family_string(
    function(family) base::identical(family, "poisson")
  ))
  expect_true(compares_family_string(
    function(family) is.element(family, c("binomial", "poisson"))
  ))
  expect_true(compares_family_string(function(family) {
    engine <- c("binomial", "poisson")
    family %in% engine
  }))
  expect_true(compares_family_string(
    function(family) base::switch(family, poisson = 1, 2)
  ))
  expect_true(compares_family_string(
    function(x) c(gaussian = "linear", binomial = "binomial")[[x$family]]
  ))
  expect_false(compares_family_string(function(family) family_traits(family)))
  expect_false(compares_family_string(
    function(family = "binomial") paste("a", family)
  ))
  expect_false(compares_family_string(function(x) x[, 1L] == "a"))
})

test_that("only the family's own functions compare family strings", {
  listed <- readLines(test_path("family-strings-to-convert.txt"))
  listed <- trimws(listed[!grepl("^#", listed) & nzchar(trimws(listed))])
  comparing <- functions_comparing_family_strings()
  unlisted <- setdiff(comparing, c(family_functions, listed))
  expect(
    length(unlisted) == 0L,
    paste("These functions compare a family string. Read the family through",
          "its object (R/family.R) instead.", paste(unlisted, collapse = ", "))
  )
  converted <- setdiff(listed, comparing)
  expect(
    length(converted) == 0L,
    paste("These functions no longer compare a family string. Remove them",
          "from family-strings-to-convert.txt.",
          paste(converted, collapse = ", "))
  )
  # Where the source is kept, as with devtools::load_all(), the family's own
  # functions must be those of R/family.R.
  files <- vapply(c(family_functions, listed), function(name) {
    file <- utils::getSrcFilename(get(name, envir = asNamespace("RiskMap")))
    if (length(file)) file else NA_character_
  }, character(1))
  known <- !is.na(files)
  expect_true(all(files[intersect(names(files)[known], family_functions)] ==
                    "family.R"))
  expect_false(any(files[intersect(names(files)[known], listed)] ==
                     "family.R"))
})
