# RiskMap raises every condition through riskmap_abort(), riskmap_warn() and
# riskmap_inform(), so that each carries the RiskMap classes and names the
# function the user called. This test reads every function in the installed
# package and fails if one calls stop(), warning(), message(), stopifnot(),
# signalCondition() or the cli and rlang functions that raise or print
# conditions, directly or through do.call().

direct_condition_calls <- c(
  "stop", "warning", "message", "stopifnot", "signalCondition",
  "abort", "warn", "inform", "cnd_signal",
  "cli_abort", "cli_warn", "cli_inform", "cli_text", "cli_bullets",
  "cli_alert", "cli_alert_success", "cli_alert_info", "cli_alert_warning",
  "cli_alert_danger"
)
# The constructors themselves, and the one helper that raises a caught
# condition again unchanged.
condition_constructors <- c("riskmap_abort", "riskmap_warn", "riskmap_inform",
                            "riskmap_rethrow")

# The functions called by name, and the first argument of each do.call()
called_functions <- function(object) {
  tokens <- utils::getParseData(parse(text = deparse(object),
                                      keep.source = TRUE))
  calls <- tokens$text[tokens$token == "SYMBOL_FUNCTION_CALL"]
  do_calls <- tokens$id[tokens$token == "SYMBOL_FUNCTION_CALL" &
                          tokens$text == "do.call"]
  targets <- vapply(do_calls, function(id) {
    call_id <- tokens$parent[tokens$id == tokens$parent[tokens$id == id]]
    arguments <- tokens[tokens$parent == call_id & tokens$token == "expr", ]
    arguments <- arguments[order(arguments$line1, arguments$col1), ]
    if (nrow(arguments) < 2L) return(NA_character_)
    first <- tokens$text[tokens$parent == arguments$id[2L] &
                           tokens$token %in% c("SYMBOL", "STR_CONST")]
    if (!length(first)) NA_character_ else gsub("^[\"']|[\"']$", "", first[1L])
  }, character(1))
  c(calls, targets)
}

functions_raising_conditions_directly <- function() {
  namespace <- asNamespace("RiskMap")
  names <- ls(namespace, all.names = TRUE)
  raises <- vapply(names, function(name) {
    object <- get(name, envir = namespace)
    if (!is.function(object) || is.primitive(object)) return(FALSE)
    any(called_functions(object) %in% direct_condition_calls)
  }, logical(1))
  setdiff(names[raises], condition_constructors)
}

test_that("the check finds direct calls, also through do.call()", {
  expect_true("stop" %in% called_functions(function() stop("x")))
  expect_true("stop" %in% called_functions(function() base::stop("x")))
  expect_true("cli_abort" %in% called_functions(function() cli::cli_abort("x")))
  expect_true("stop" %in% called_functions(function() do.call(stop, list("x"))))
  expect_true("warning" %in%
                called_functions(function() do.call("warning", list("x"))))
  expect_false(any(called_functions(function() do.call(paste, list("x"))) %in%
                     direct_condition_calls))
})

test_that("conditions are raised through the RiskMap constructors", {
  direct <- functions_raising_conditions_directly()
  expect(
    length(direct) == 0L,
    paste("These functions raise conditions without riskmap_abort(),",
          "riskmap_warn() or riskmap_inform().",
          paste(direct, collapse = ", "))
  )
})

test_that("leaving out a required argument gives a RiskMap input error", {
  error <- expect_input_error(check_fit(), "object",
                              "`object` is missing, with no default.")
  expect_identical(as.character(error$call[[1L]]), "check_fit")
  expect_input_error(glgpm(y ~ gp(), family = "gaussian"), "data")
  expect_input_error(variogram(gaussian_data), "variable")
  # An argument passed on without a value is missing too, and the error
  # names the function the user called.
  wrapper <- function(fit) check_fit(fit)
  error <- expect_input_error(wrapper(), "object")
  expect_identical(as.character(error$call[[1L]]), "check_fit")
})

test_that("a message without a bullet gets the information sign", {
  plain <- rlang::catch_cnd(riskmap_inform("A plain message."))
  expect_s3_class(plain, "riskmap_message")
  expect_match(rlang::cnd_message(plain),
               paste(cli::symbol$info, "A plain message."), fixed = TRUE)
  done <- rlang::catch_cnd(riskmap_inform(c(v = "Done.")))
  expect_match(rlang::cnd_message(done), paste(cli::symbol$tick, "Done."),
               fixed = TRUE)
})

test_that("a prediction warns once for the blocks that needed a jitter", {
  # Two blocks, each of two points at one place, so each block's covariance
  # is singular and needs a jitter.
  prediction <- rbind(c(0, 0), c(0, 0), c(5, 5), c(5, 5))
  conditioning <- rbind(c(50, 50), c(60, 60))
  no_weights <- function(cross) matrix(0, nrow(cross), ncol(cross))
  warnings <- testthat::capture_warnings(grouped_joint_prediction(
    prediction, conditioning, group_indices = list(1:2, 3:4),
    weight_function = no_weights, conditional_signal = c(0, 0),
    marginal_variance = 1, phi = 1, kappa = 0.5, n_samples = 3
  ))
  expect_length(warnings, 1L)
  expect_match(warnings, paste("The conditional prediction covariance of 2",
                               "of 2 blocks needed a diagonal jitter"),
               fixed = TRUE)
  # One matrix keeps the warning it always had.
  expect_warning(warn_jitter("covariance matrix", 1e-12),
                 "needed a diagonal jitter of 1e-12 to be factorised")
  expect_no_warning(warn_jitter("covariance matrix", c(0, 0)))
})
