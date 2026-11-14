# Compares two arms of the golden battery of tools/golden/run.R task by task (tools/golden/README.md). What may
# legitimately differ between runs is removed first: timings, environments and the environments of functions and
# formulas (functions are compared by their deparsed text, external pointers not at all). Printed text is compared
# line by line with times masked. Writes one line per task and, for a task that differs, the first differing paths
# with the largest difference of each, then the environments of the two arms. Exits with status 1 if a task fails.
#   Rscript --vanilla tools/golden/compare.R --a=<dir> --b=<dir> [--report=<file>] [--declared=<file or list>]
#     [--tolerance=<k> [--se-tolerance=<fraction>] [--loglik-tolerance=<units>]] [--allow-new]
# Modes.
# - Identity (the default). A task passes when its value and text are identical in the two arms.
# - Declared paths (--declared). Paths that a PR changes on purpose, such as /fit/spec, are removed before the
#   comparison and reported, so the rest must still be identical. Each entry is a regular expression matched against
#   paths as the report prints them (/fit/estimate/phi, /fit@optimizer for an attribute, /text/2 for the second
#   printed text), optionally preceded by a regular expression for the task names it applies to and a space. Entries
#   are lines of a file (# starts a comment) or a comma-separated list.
# - Tolerance (--tolerance=k), for a PR that changes the optimiser's path. A task passes when the statuses and the
#   structure of the values agree, every fit the task recorded keeps its outcome, certification and number of stages,
#   its estimates move by at most k standard errors of arm a, its standard errors by at most a fraction
#   --se-tolerance (k by default) and its log-likelihood by at most --loglik-tolerance (k by default), the printed
#   text and strings agree once numbers are masked, and its other numbers differ only if one of its fits moved.
# --allow-new lets arm b have tasks that arm a lacks, for an arm that adds tasks. They are listed, not compared. A task
# missing from arm b always fails, and so does a comparison of no task.
args <- commandArgs(TRUE)
option <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
flag <- function(name) any(args == paste0("--", name))
known <- c("a", "b", "report", "declared", "tolerance", "se-tolerance", "loglik-tolerance")
unknown_options <- args[!(sub("=.*", "", sub("^--", "", args)) %in% known & grepl("=", args)) & args != "--allow-new"]
if (length(unknown_options)) stop("unknown options: ", paste(unknown_options, collapse = " "))
suppressMessages(library(sf))
`%||%` <- function(a, b) if (is.null(a)) b else a
dir_a <- option("a"); dir_b <- option("b"); report_path <- option("report")
if (is.null(dir_a) || is.null(dir_b)) stop("give the two arms with --a=<dir> --b=<dir>")
for (dir in c(dir_a, dir_b)) if (!dir.exists(dir)) stop("no arm at ", dir)
positive <- function(name, default = NULL) {
  value <- option(name)
  if (is.null(value)) return(default)
  value <- suppressWarnings(as.numeric(value))
  if (!is.finite(value) || value <= 0) stop("--", name, " must be a positive number")
  value
}
tolerance <- positive("tolerance")
se_tolerance <- positive("se-tolerance", tolerance)
loglik_tolerance <- positive("loglik-tolerance", tolerance)
if (is.null(tolerance) && (!is.null(se_tolerance) || !is.null(loglik_tolerance))) {
  stop("--se-tolerance and --loglik-tolerance need --tolerance")
}
timing <- "seconds|timing|^started$|^elapsed|^time$|_time$|duration|^timestamp"
structural <- c("names", "dim", "dimnames", "class", "row.names", "levels", "tsp")
dropped <- new.env(); dropped$names <- character()

# Declared paths, as a table of task and path patterns.
read_declared <- function(x) {
  if (is.null(x)) return(data.frame(task = character(), path = character()))
  entries <- if (file.exists(x)) readLines(x) else strsplit(x, ",")[[1]]
  entries <- trimws(sub("#.*$", "", entries))
  entries <- entries[nzchar(entries)]
  parts <- strsplit(entries, "[[:space:]]+")
  if (any(lengths(parts) > 2)) stop("a declared entry has more than a task and a path pattern")
  table <- data.frame(task = vapply(parts, function(p) if (length(p) == 2) p[1] else ".", ""),
                      path = vapply(parts, function(p) p[length(p)], ""))
  for (pattern in c(table$task, table$path)) grepl(pattern, "", perl = TRUE)
  table
}
declared <- read_declared(option("declared"))
matches <- function(paths, patterns) {
  vapply(paths, function(p) any(vapply(patterns, grepl, logical(1), x = p, perl = TRUE)), logical(1))
}

normalise <- function(x) {
  if (is.function(x)) return(paste(deparse(x), collapse = "\n"))
  if (is.environment(x)) return("<environment>")
  if (typeof(x) == "externalptr") return("<externalptr>")
  if (inherits(x, "formula")) return(paste(deparse(x), collapse = " "))
  if (is.language(x)) return(paste(deparse(x), collapse = " "))
  at <- attributes(x)
  if (is.list(x)) {
    keep <- if (is.null(names(x))) rep(TRUE, length(x)) else !grepl(timing, names(x))
    dropped$names <- union(dropped$names, names(x)[!keep])
    x <- lapply(unclass(x)[keep], normalise)
    if (!is.null(at$names)) at$names <- at$names[keep]
  }
  if (length(at)) {
    other <- setdiff(names(at), structural)
    other_keep <- other[!grepl(timing, other)]
    dropped$names <- union(dropped$names, setdiff(other, other_keep))
    for (name in setdiff(other, other_keep)) at[[name]] <- NULL
    for (name in other_keep) at[[name]] <- normalise(at[[name]])
    attributes(x) <- at
  }
  x
}

mask_times <- function(lines) {
  lines <- gsub("[0-9]+(\\.[0-9]+)?( )?(s|sec|secs|second|seconds|minute|minutes|min|hour|hours)\\b", "<time>", lines)
  # A time whose unit print() has wrapped onto the next line, as in "The fit took 14" and "seconds.".
  gsub("took [0-9]+(\\.[0-9]+)?$", "took <time>", lines)
}
# For the tolerance mode, numbers are masked and runs of spaces collapsed, since the width of a number moves columns.
mask_numbers <- function(lines) {
  lines <- gsub("(?<![A-Za-z_])[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?", "#", lines, perl = TRUE)
  gsub("[[:space:]]+", " ", lines)
}

keys <- function(x) {
  k <- names(x) %||% rep("", length(x))
  k[!nzchar(k)] <- as.character(seq_along(x))[!nzchar(k)]
  k
}

# Removes the declared paths from a normalised value, keeping what was removed in removed$items.
prune <- function(x, path, patterns, removed) {
  at <- attributes(x)
  y <- x
  if (is.list(x)) {
    child <- if (length(x)) paste0(path, "/", keys(x)) else character()
    drop <- matches(child, patterns)
    for (i in which(drop)) removed$items[child[i]] <- list(unclass(x)[[i]])
    kept <- which(!drop)
    y <- lapply(kept, function(i) prune(unclass(x)[[i]], child[i], patterns, removed))
    if (!is.null(at$names)) at$names <- at$names[kept]
  }
  for (name in setdiff(names(at), structural)) {
    p <- paste0(path, "@", name)
    if (matches(p, patterns)) {
      removed$items[p] <- list(at[[name]])
      at[[name]] <- NULL
    } else {
      at[[name]] <- prune(at[[name]], p, patterns, removed)
    }
  }
  attributes(y) <- at
  y
}

# Paths of the leaves where two normalised objects differ, each with what differs and its kind (number, string,
# logical or structure). Lists and data frames are followed by names, attributes other than the structural ones by @.
differences <- function(a, b, path = "", found = list(), limit = 12) {
  if (length(found) >= limit || identical(a, b)) return(found)
  before <- length(found)
  add <- function(found, what, kind, where = path) {
    found[[length(found) + 1]] <- list(path = where, what = what, kind = kind)
    found
  }
  if (is.list(a) && is.list(b)) {
    names_a <- keys(a); names_b <- keys(b)
    for (name in union(names_a, names_b)) {
      ia <- match(name, names_a); ib <- match(name, names_b)
      if (is.na(ia) || is.na(ib)) {
        found <- add(found, if (is.na(ia)) "only in b" else "only in a", "structure", paste0(path, "/", name))
      } else {
        found <- differences(unclass(a)[[ia]], unclass(b)[[ib]], paste0(path, "/", name), found, limit)
      }
      if (length(found) >= limit) return(found)
    }
  } else if (is.numeric(a) && is.numeric(b) && length(a) == length(b) && !is.factor(a) && !is.factor(b)) {
    va <- as.numeric(a); vb <- as.numeric(b)
    if (!identical(is.na(va), is.na(vb))) {
      found <- add(found, sprintf("missing values differ (%d against %d)", sum(is.na(va)), sum(is.na(vb))),
                   "structure")
    } else if (!identical(va, vb)) {
      both <- !is.na(va) & is.finite(va) & is.finite(vb)
      gap <- if (any(both)) max(abs(va - vb)[both]) else NA_real_
      scale <- pmax(abs(va), abs(vb))[both]
      relative <- if (any(both)) max(ifelse(scale > 0, abs(va - vb)[both] / scale, 0)) else NA_real_
      found <- add(found, sprintf("largest difference %.3g, relative %.3g", gap, relative), "number")
      found[[length(found)]]$relative <- relative
    }
  } else if (is.character(a) && is.character(b) && length(a) == length(b)) {
    if (!identical(as.vector(a), as.vector(b))) {
      differ <- which(a != b | is.na(a) != is.na(b))
      found <- add(found, sprintf("%d of %d strings differ, first: %s | %s", length(differ), length(a),
                                  a[differ[1]], b[differ[1]]), "string")
      found[[length(found)]]$a <- as.vector(a); found[[length(found)]]$b <- as.vector(b)
    }
  } else if (is.logical(a) && is.logical(b) && length(a) == length(b)) {
    if (!identical(as.vector(a), as.vector(b))) found <- add(found, "logical values differ", "logical")
  } else {
    return(add(found, sprintf("%s of length %d against %s of length %d", class(a)[1], length(a), class(b)[1],
                              length(b)), "structure"))
  }
  # Structural attributes, without the names of a list, whose differences are reported above as paths.
  at_a <- attributes(a); at_b <- attributes(b)
  struct <- function(at) {
    keep <- sort(intersect(names(at), if (is.list(a) && is.list(b)) setdiff(structural, "names") else structural))
    at[keep]
  }
  if (!identical(struct(at_a), struct(at_b))) found <- add(found, "attributes differ", "structure", paste0(path, " (attributes)"))
  for (name in union(setdiff(names(at_a), structural), setdiff(names(at_b), structural))) {
    if (length(found) >= limit) break
    where <- paste0(path, "@", name)
    if (is.null(at_a[[name]]) || is.null(at_b[[name]])) {
      found <- add(found, if (is.null(at_a[[name]])) "only in b" else "only in a", "structure", where)
    } else {
      found <- differences(at_a[[name]], at_b[[name]], where, found, limit)
    }
  }
  if (length(found) == before) found <- add(found, sprintf("%s against %s", typeof(a), typeof(b)), "structure")
  found
}

# The fits a task recorded, record$fits for the base fit of a downstream task and the fit of a fit task in its value
# (record$fits$fit in arms of the first version of this script, which kept it twice).
fits_of <- function(record) {
  fits <- record$fits %||% list()
  if (is.null(fits$fit) && inherits(record$value$fit, "RiskMap")) fits$fit <- record$value$fit
  fits
}
summarise_fit <- function(fit) {
  o <- attr(fit, "optimizer")
  estimate <- unlist(fit$estimate)
  se <- rep(NA_real_, length(estimate)); names(se) <- names(estimate)
  if (is.matrix(fit$covariance)) {
    v <- suppressWarnings(sqrt(diag(fit$covariance)))
    names(v) <- rownames(fit$covariance) %||% colnames(fit$covariance)
    common <- intersect(names(estimate), names(v))
    se[common] <- v[common]
  }
  list(estimate = estimate, se = se, log_lik = as.numeric(fit$log_lik %||% NA)[1], outcome = o$mcml_outcome,
       certified = o$mcml_certified, stages = o$mcml_stages,
       proposals = vapply(o$mcml_history %||% list(), function(h) as.character(h$proposal %||% NA), ""))
}
# Judges the recorded fits of a task against the tolerances. A fit is pruned of the declared paths, as if it sat at
# /fit for a fit task and at /base for a downstream task, before it is called moved. A moved fit must keep its
# outcome, certification, stages and parameters, its estimates must move by at most k standard errors (of arm a, or
# of arm b where a has none), its standard errors by at most a fraction se_k, and its log-likelihood by at most
# loglik_k. Returns whether any fit moved, the reasons for a failure, notes, and the largest changes.
judge_fits <- function(fits_a, fits_b, k, se_k, loglik_k, patterns) {
  out <- list(moved = FALSE, fail = character(), note = character(), shift = 0, shift_at = "", se_change = 0,
              loglik_change = 0)
  for (name in union(names(fits_a), names(fits_b))) {
    if (is.null(fits_a[[name]]) || is.null(fits_b[[name]])) {
      out$fail <- c(out$fail, sprintf("%s only in %s", name, if (is.null(fits_a[[name]])) "b" else "a"))
      next
    }
    where <- paste0("/", name)
    pa <- normalise(fits_a[[name]]); pb <- normalise(fits_b[[name]])
    if (length(patterns)) {
      ignored <- new.env(); ignored$items <- list()
      pa <- prune(pa, where, patterns, ignored); pb <- prune(pb, where, patterns, ignored)
    }
    if (identical(pa, pb)) next
    out$moved <- TRUE
    fail <- function(...) out$fail <<- c(out$fail, paste0(name, ": ", sprintf(...)))
    sa <- summarise_fit(fits_a[[name]]); sb <- summarise_fit(fits_b[[name]])
    for (what in c("outcome", "certified", "stages")) {
      if (!identical(sa[[what]], sb[[what]])) {
        fail("%s %s against %s", what, format(sa[[what]] %||% "none"), format(sb[[what]] %||% "none"))
      }
    }
    if (!identical(sa$proposals, sb$proposals)) {
      out$note <- c(out$note, sprintf("%s: proposals by stage %s against %s", name,
                                      paste(sa$proposals, collapse = " "), paste(sb$proposals, collapse = " ")))
    }
    if (!identical(names(sa$estimate), names(sb$estimate))) {
      fail("parameters %s against %s", paste(names(sa$estimate), collapse = " "),
           paste(names(sb$estimate), collapse = " "))
      next
    }
    missing <- names(sa$estimate)[is.finite(sa$estimate) != is.finite(sb$estimate)]
    if (length(missing)) fail("%s finite in one arm only", paste(missing, collapse = ", "))
    gap <- abs(sb$estimate - sa$estimate)
    gap[!is.finite(gap)] <- 0
    se <- ifelse(is.finite(sa$se) & sa$se > 0, sa$se, sb$se)
    judged <- is.finite(se) & se > 0
    if (!any(judged) && any(gap > 0)) fail("its estimates move and no parameter has a standard error")
    shift <- gap[judged] / se[judged]
    if (length(shift) && max(shift) > out$shift) {
      out$shift <- max(shift); out$shift_at <- paste0(name, " ", names(shift)[which.max(shift)])
    }
    beyond <- names(shift)[shift > k]
    if (length(beyond)) {
      fail("%s moves by %s standard errors", paste(beyond, collapse = ", "), paste(signif(shift[beyond], 3), collapse = ", "))
    }
    loose <- names(gap)[!judged & gap > 0]
    if (length(loose)) {
      out$note <- c(out$note, sprintf("%s: %s without a standard error moves by %s", name,
                                      paste(loose, collapse = ", "), paste(signif(gap[loose], 3), collapse = ", ")))
    }
    lost <- names(sa$se)[is.finite(sa$se) != is.finite(sb$se)]
    if (length(lost)) fail("standard error of %s finite in one arm only", paste(lost, collapse = ", "))
    both <- is.finite(sa$se) & is.finite(sb$se) & sa$se > 0
    change <- abs(sb$se[both] - sa$se[both]) / sa$se[both]
    if (length(change)) {
      out$se_change <- max(out$se_change, change)
      wide <- names(change)[change > se_k]
      if (length(wide)) {
        fail("standard error of %s changes by %s", paste(wide, collapse = ", "),
             paste(sprintf("%.3g%%", 100 * change[wide]), collapse = ", "))
      }
    }
    if (is.finite(sa$log_lik) != is.finite(sb$log_lik)) {
      fail("log-likelihood finite in one arm only")
    } else if (is.finite(sa$log_lik)) {
      moved_by <- abs(sb$log_lik - sa$log_lik)
      out$loglik_change <- max(out$loglik_change, moved_by)
      if (moved_by > loglik_k) fail("log-likelihood moves by %.3g", moved_by)
    }
  }
  out
}

text_paths <- function(text) if (length(text)) paste0("/text/", seq_along(text)) else character()
compare_text <- function(text_a, text_b, mask) {
  lines <- character()
  for (i in seq_len(max(length(text_a), length(text_b)))) {
    ta <- mask(if (i <= length(text_a)) text_a[[i]] else character())
    tb <- mask(if (i <= length(text_b)) text_b[[i]] else character())
    if (!identical(ta, tb)) {
      n <- min(length(ta), length(tb))
      k <- which(ta[seq_len(n)] != tb[seq_len(n)])[1]
      if (is.na(k)) k <- n + 1
      lines <- c(lines, sprintf("    text %d, line %d: %s | %s", i, k, if (k <= length(ta)) ta[k] else "(none)",
                                if (k <= length(tb)) tb[k] else "(none)"))
    }
  }
  lines
}

tasks_a <- sub("\\.rds$", "", list.files(dir_a, "\\.rds$")); tasks_b <- sub("\\.rds$", "", list.files(dir_b, "\\.rds$"))
tasks <- sort(union(tasks_a, tasks_b))
mode <- if (is.null(tolerance)) "identity" else
  sprintf("tolerance of %g standard errors, standard errors within %g%%, log-likelihood within %g", tolerance,
          100 * se_tolerance, loglik_tolerance)
lines <- c(sprintf("Golden battery: a = %s, b = %s, %s", dir_a, dir_b, mode))
if (nrow(declared)) {
  lines <- c(lines, sprintf("Declared paths: %s", paste(ifelse(declared$task == ".", declared$path,
                                                                paste(declared$task, declared$path)), collapse = ", ")))
}
failures <- 0L; compared <- 0L; skipped <- character()
environments <- list(a = list(), b = list())
for (task in tasks) {
  path_a <- file.path(dir_a, paste0(task, ".rds")); path_b <- file.path(dir_b, paste0(task, ".rds"))
  if (!file.exists(path_a) && flag("allow-new")) { skipped <- c(skipped, task); next }
  if (!file.exists(path_a) || !file.exists(path_b)) {
    lines <- c(lines, sprintf("%-36s MISSING in %s", task, if (!file.exists(path_a)) "a" else "b"))
    failures <- failures + 1L; next
  }
  compared <- compared + 1L
  a <- readRDS(path_a); b <- readRDS(path_b)
  environments$a[[task]] <- a$environment; environments$b[[task]] <- b$environment
  if (!identical(a$status, b$status)) {
    lines <- c(lines, sprintf("%-36s STATUS %s against %s (%s | %s)", task, a$status, b$status, a$message %||% "",
                              b$message %||% ""))
    failures <- failures + 1L; next
  }
  if (identical(a$status, "failed")) {
    same <- identical(a$message, b$message) ||
      (!is.null(tolerance) && identical(mask_numbers(a$message), mask_numbers(b$message)))
    lines <- c(lines, sprintf("%-36s %s both failed: %s", task, if (same) "same" else "DIFFERENT", a$message))
    if (!same) failures <- failures + 1L
    next
  }
  value_a <- a$value; value_b <- b$value
  text_a <- value_a$text; text_b <- value_b$text
  value_a$text <- NULL; value_b$text <- NULL
  na <- normalise(value_a); nb <- normalise(value_b)
  patterns <- declared$path[vapply(declared$task, grepl, logical(1), x = task, perl = TRUE)]
  declared_lines <- character()
  if (length(patterns)) {
    removed_a <- new.env(); removed_a$items <- list(); removed_b <- new.env(); removed_b$items <- list()
    na <- prune(na, "", patterns, removed_a); nb <- prune(nb, "", patterns, removed_b)
    hidden_a <- matches(text_paths(text_a), patterns); hidden_b <- matches(text_paths(text_b), patterns)
    removed_a$items[text_paths(text_a)[hidden_a]] <- text_a[hidden_a]
    removed_b$items[text_paths(text_b)[hidden_b]] <- text_b[hidden_b]
    text_a[hidden_a] <- list(NULL); text_b[hidden_b] <- list(NULL)
    for (p in sort(union(names(removed_a$items), names(removed_b$items)))) {
      state <- if (is.null(removed_a$items[[p]])) "only in b" else if (is.null(removed_b$items[[p]])) "only in a" else
        if (identical(removed_a$items[[p]], removed_b$items[[p]])) "same" else "differs"
      declared_lines <- c(declared_lines, sprintf("    declared %s: %s", p, state))
    }
  }
  same_value <- identical(na, nb)
  # Texts are compared by position, so that a declared text, removed from one arm only, leaves no gap.
  same_text <- !length(compare_text(text_a, text_b, mask_times))
  declared_note <- if (length(declared_lines)) sprintf(", declared paths %s",
                                                       if (any(grepl("(differs|only in [ab])$", declared_lines))) "differ" else "same") else ""
  if (same_value && same_text) {
    lines <- c(lines, sprintf("%-36s identical%s (%.0f s and %.0f s)", task, declared_note, a$seconds, b$seconds),
               declared_lines[grepl("(differs|only in [ab])$", declared_lines)])
    next
  }
  if (is.null(tolerance)) {
    failures <- failures + 1L
    lines <- c(lines, sprintf("%-36s DIFFERENT%s%s", task, if (!same_value) " value" else "",
                              if (!same_text) " text" else ""))
    if (!same_value) for (d in differences(na, nb)) lines <- c(lines, sprintf("    %s  %s", d$path, d$what))
    if (!same_text) lines <- c(lines, compare_text(text_a, text_b, mask_times))
    lines <- c(lines, declared_lines[grepl("(differs|only in [ab])$", declared_lines)])
    next
  }
  # The tolerance mode.
  found <- differences(na, nb, limit = Inf)
  kinds <- vapply(found, `[[`, "", "kind")
  strings <- found[kinds == "string"]
  worded <- vapply(strings, function(d) !identical(mask_numbers(d$a), mask_numbers(d$b)), logical(1))
  judged <- judge_fits(fits_of(a), fits_of(b), tolerance, se_tolerance, loglik_tolerance, patterns)
  reasons <- c(judged$fail,
               vapply(found[kinds == "structure"], function(d) sprintf("%s  %s", d$path, d$what), ""),
               vapply(strings[worded], function(d) sprintf("%s  %s", d$path, d$what), ""))
  numbers <- found[kinds %in% c("number", "logical") | (kinds == "string")]
  if (length(numbers) && !judged$moved) {
    reasons <- c(reasons, sprintf("%d values differ with no recorded fit moved, first %s  %s", length(numbers),
                                  numbers[[1]]$path, numbers[[1]]$what))
  }
  text_lines <- compare_text(text_a, text_b, function(x) mask_numbers(mask_times(x)))
  if (length(text_lines)) reasons <- c(reasons, "printed text differs beyond its numbers", text_lines)
  relative <- vapply(found[kinds == "number"], function(d) d$relative %||% NA_real_, 0)
  worst <- if (length(relative) && any(is.finite(relative))) which.max(replace(relative, !is.finite(relative), -1))
  summary_line <- sprintf(paste0("largest shift %.3g standard errors%s, standard errors change by up to %.3g%%, the ",
                                 "log-likelihood by %.3g, %d numbers, %d logical values and %d strings differ%s"),
                          judged$shift, if (nzchar(judged$shift_at)) sprintf(" (%s)", judged$shift_at) else "",
                          100 * judged$se_change, judged$loglik_change,
                          sum(kinds == "number"), sum(kinds == "logical"), sum(kinds == "string"),
                          if (!is.null(worst)) sprintf(", largest relative difference %.3g at %s", relative[worst],
                                                       found[kinds == "number"][[worst]]$path) else "")
  if (length(reasons)) {
    failures <- failures + 1L
    lines <- c(lines, sprintf("%-36s OUTSIDE tolerance, %s", task, summary_line), paste0("    ", head(reasons, 12)))
  } else {
    lines <- c(lines, sprintf("%-36s within tolerance, %s", task, summary_line))
  }
  lines <- c(lines, if (length(judged$note)) paste0("    note: ", judged$note),
             declared_lines[grepl("(differs|only in [ab])$", declared_lines)])
}

# The environments of the two arms. Bits reproduce only on one node type with the same R, BLAS and packages.
describe <- function(records, field) {
  values <- unique(unlist(lapply(records, function(e) e[[field]])))
  if (!length(values)) "unrecorded" else paste(values, collapse = ", ")
}
lines <- c(lines, "")
for (field in c("commit", "riskmap", "r", "blas", "host")) {
  da <- describe(environments$a, field); db <- describe(environments$b, field)
  lines <- c(lines, sprintf("%-8s a: %s%s", field, da, if (identical(da, db)) ", the same in b" else sprintf(" | b: %s", db)))
}
packages_of <- function(records) {
  versions <- unlist(unname(lapply(records, function(e) e$packages)))
  if (is.null(versions)) return(NULL)
  tapply(versions, names(versions), function(v) paste(unique(v), collapse = "/"))
}
pa <- packages_of(environments$a); pb <- packages_of(environments$b)
if (!is.null(pa) && !is.null(pb)) {
  common <- intersect(names(pa), names(pb))
  changed <- common[pa[common] != pb[common]]
  lines <- c(lines, if (length(changed)) sprintf("WARNING packages differ: %s", paste(sprintf("%s %s | %s", changed,
                                                   pa[changed], pb[changed]), collapse = ", ")) else
               sprintf("packages: the %d loaded in both arms have the same versions", length(common)))
} else {
  lines <- c(lines, "packages: not recorded in both arms")
}
if (length(skipped)) lines <- c(lines, sprintf("New in b, not compared: %s", paste(skipped, collapse = ", ")))
if (compared == 0L) lines <- c(lines, "No task is in both arms, which fails the comparison.")
lines <- c(lines, sprintf("Names removed as timings: %s", paste(sort(dropped$names), collapse = ", ")),
           sprintf("Tasks: %d compared, %d %s, all %s: %s", compared, failures,
                   if (is.null(tolerance)) "differ" else "outside the tolerance",
                   if (is.null(tolerance)) "identical" else "within tolerance", failures == 0L && compared > 0L))
writeLines(lines)
if (!is.null(report_path)) writeLines(lines, report_path)
quit(save = "no", status = if (failures > 0L || compared == 0L) 1L else 0L)
