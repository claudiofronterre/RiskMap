##' Prepare Random-Effect Group Indices
##'
##' Extracts grouping variables from a data-frame-like object and converts
##' their values to the internal indices used by RiskMap model-fitting and
##' simulation functions.
##'
##' @param data A data-frame-like object containing the grouping variables.
##' @param terms A character vector naming the grouping variables.
##'
##' @return A list containing the number and names of random effects, their
##'   observation-level indices, internal level indices, and display labels.
##'
##' @noRd
prepare_random_effects <- function(data, terms) {
  if (is.null(terms)) {
    result <- list(
      n_re = 0L,
      names_re = NULL,
      ID_re = NULL,
      re_unique = NULL,
      re_unique_f = NULL
    )
    class(result) <- "RiskMap_random_effects"
    return(result)
  }

  missing_terms <- setdiff(terms, names(data))
  if (length(missing_terms) > 0L) {
    stop(
      "Random-effect variable(s) not found in `data`: ",
      paste(missing_terms, collapse = ", "),
      "."
    )
  }

  encoded_groups <- lapply(
    terms,
    function(term) encode_random_effect(data[[term]], term)
  )
  names(encoded_groups) <- terms

  ID_re <- matrix(
    unlist(
      lapply(encoded_groups, function(group) group$ids),
      use.names = FALSE
    ),
    nrow = nrow(data),
    ncol = length(terms),
    dimnames = list(NULL, terms)
  )

  result <- list(
    n_re = length(terms),
    names_re = terms,
    ID_re = ID_re,
    re_unique = lapply(encoded_groups, function(group) group$indices),
    re_unique_f = lapply(encoded_groups, function(group) group$labels)
  )
  class(result) <- "RiskMap_random_effects"
  result
}

##' Encode a Random-Effect Grouping Variable
##'
##' @param group A factor, character, integer, or numeric grouping vector.
##' @param variable_name A character string naming the grouping variable.
##'
##' @return A list containing observation-level indices, internal level
##'   indices, and original display labels.
##'
##' @noRd
encode_random_effect <- function(group, variable_name) {
  if (anyNA(group)) {
    stop(
      "Random-effect variable `", variable_name,
      "` contains missing values."
    )
  }

  if (is.factor(group)) {
    # remove empty levels
    group <- droplevels(group)
    labels <- levels(group)
    ids <- as.integer(group)
  } else if (is.character(group)) {
    labels <- unique(group)
    ids <- match(group, labels)
  } else if (is.numeric(group)) {
    if (any(!is.finite(group))) {
      stop(
        "Random-effect variable `", variable_name,
        "` contains non-finite values."
      )
    }
    labels <- unique(group)
    ids <- match(group, labels)
  } else {
    stop(
      "Random-effect variable `", variable_name,
      "` must be a factor, character, integer, or numeric vector."
    )
  }

  result <- list(
    ids = as.integer(ids),
    indices = seq_along(labels),
    labels = labels
  )
  class(result) <- "RiskMap_random_effect"
  result
}
