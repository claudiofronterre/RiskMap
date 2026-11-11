# RiskMap raises every condition through riskmap_abort(), riskmap_warn() and
# riskmap_inform(), so that each carries the RiskMap classes. Functions not yet
# converted are listed in conditions-to-convert.txt. The list may only shrink,
# because a new function that raises a condition directly fails here, and so
# does a listed function that no longer needs to be listed.

direct_condition_calls <- c("stop", "warning", "message", "stopifnot",
                            "cli_abort", "cli_warn", "cli_inform",
                            "abort", "warn", "inform")
condition_constructors <- c("riskmap_abort", "riskmap_warn", "riskmap_inform")

functions_raising_conditions_directly <- function() {
  namespace <- asNamespace("RiskMap")
  names <- ls(namespace, all.names = TRUE)
  raises <- vapply(names, function(name) {
    object <- get(name, envir = namespace)
    if (!is.function(object) || is.primitive(object)) return(FALSE)
    code <- parse(text = deparse(object), keep.source = TRUE)
    tokens <- utils::getParseData(code)
    calls <- tokens$text[tokens$token == "SYMBOL_FUNCTION_CALL"]
    any(calls %in% direct_condition_calls)
  }, logical(1))
  setdiff(names[raises], condition_constructors)
}

test_that("conditions are raised through the RiskMap constructors", {
  listed <- trimws(readLines(test_path("conditions-to-convert.txt")))
  listed <- listed[nzchar(listed) & !startsWith(listed, "#")]
  direct <- functions_raising_conditions_directly()

  unlisted <- setdiff(direct, listed)
  expect(
    length(unlisted) == 0L,
    paste("These functions raise conditions without riskmap_abort(),",
          "riskmap_warn() or riskmap_inform().",
          paste(unlisted, collapse = ", "))
  )
  converted <- setdiff(listed, direct)
  expect(
    length(converted) == 0L,
    paste("These functions no longer raise conditions directly, so remove",
          "them from conditions-to-convert.txt.",
          paste(converted, collapse = ", "))
  )
})
