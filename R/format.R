#' Printed output of RiskMap objects
#'
#' The print methods build their lines with these helpers and write them with
#' `cat()`, through `out_print()`. They are printed output, not messages, so
#' they appear where messages are hidden, as in knitr chunks with
#' `message: false`, and in `capture.output()`. cli supplies the rules,
#' symbols, colours and inline markup, so the output follows cli's options.
#' Without colour support, in knitr and in tests, the lines are plain text,
#' and with `options(cli.unicode = FALSE)` the symbols are ASCII.
#'
#' Text is plain unless a helper is told it has cli inline markup, so that
#' braces in data, such as a formula, are printed as they are.
#' @name riskmap_output
#' @noRd
NULL

#' Write the lines of a printed object
#' @noRd
out_print <- function(lines) {
  cat(lines, sep = "\n")
  invisible(NULL)
}

#' The width of printed output
#'
#' The option `cli.width` when it is set, and otherwise R's `width`, the width
#' at which R prints the tables within the output. RStudio keeps `width` at
#' the width of its console.
#' @noRd
out_width <- function() {
  max(40L, getOption("cli.width", getOption("width", 80L)))
}

#' The first line of a printed object, a rule with its title
#' @noRd
out_rule <- function(title) {
  as.character(cli::rule(left = cli::style_bold(title), width = out_width()))
}

#' A heading within a printed object
#' @noRd
out_heading <- function(text) {
  cli::style_bold(text)
}

#' Faint text, for what matters less than the rest, such as the call
#' @noRd
out_faint <- function(lines) {
  cli::col_grey(lines)
}

#' Lines of text wrapped to the width
#'
#' Each element of `text` starts a new line. With `markup`, `text` takes cli
#' inline markup, such as `{.arg max_samples}`, interpolated in `.envir`.
#' @noRd
out_text <- function(text, markup = FALSE, indent = 0L,
                     .envir = parent.frame()) {
  if (!length(text)) return(character())
  if (markup) text <- out_markup(text, .envir)
  width <- out_width() - indent
  lines <- unlist(lapply(text, out_wrap, width = width))
  paste0(strrep(" ", indent), lines)
}

#' Wrap one paragraph to `width`, keeping code between backticks whole
#'
#' The spaces of code such as `` `gp(nugget = FALSE)` `` become non-breaking
#' while the paragraph is wrapped, as cli does with its own inline code. Code
#' longer than `width` would then stick out, so it can break as text does.
#' @noRd
out_wrap <- function(text, width) {
  code <- gregexpr("`[^`]*`", text)
  regmatches(text, code) <- lapply(regmatches(text, code), function(spans) {
    fits <- nchar(spans) + 1L <= width
    spans[fits] <- gsub(" ", "\u00a0", spans[fits], fixed = TRUE)
    spans
  })
  gsub("\u00a0", " ", cli::ansi_strwrap(text, width = width), fixed = TRUE)
}

#' Lines that start with a symbol, wrapped under their text
#'
#' The names of `bullets` choose the symbol, as in [cli::cli_bullets()]:
#' `"v"` a green tick, `"!"` a yellow exclamation mark, `"x"` a red cross,
#' `"i"` a cyan information sign, `"*"` a bullet, `">"` an arrow, and `" "`
#' no symbol, for text under the line above. `markup` is as in `out_text()`.
#' @noRd
out_bullets <- function(bullets, markup = FALSE, indent = 0L,
                        .envir = parent.frame()) {
  if (!length(bullets)) return(character())
  types <- names(bullets) %||% rep("*", length(bullets))
  if (markup) bullets <- out_markup(bullets, .envir)
  width <- out_width() - indent - 2L
  pad <- strrep(" ", indent)
  unlist(lapply(seq_along(bullets), function(k) {
    lines <- out_wrap(bullets[[k]], width = width)
    if (!length(lines)) lines <- ""
    c(paste0(pad, out_symbol(types[[k]]), " ", lines[1L]),
      if (length(lines) > 1L) paste0(pad, "  ", lines[-1L]))
  }))
}

#' Give each element of `text` the bullet `type`, for `out_bullets()`
#' @noRd
out_as_bullets <- function(text, type) {
  if (!length(text)) return(character())
  stats::setNames(text, rep(type, length(text)))
}

#' The symbol of a bullet, as cli draws it
#' @noRd
out_symbol <- function(type) {
  switch(type,
         "v" = cli::col_green(cli::symbol$tick),
         "x" = cli::col_red(cli::symbol$cross),
         "!" = cli::col_yellow("!"),
         "i" = cli::col_cyan(cli::symbol$info),
         "*" = cli::col_cyan(cli::symbol$bullet),
         ">" = cli::symbol$arrow_right,
         " ")
}

#' The symbol of the status of a check of a fit
#' @noRd
out_status_symbol <- function(status) {
  out_symbol(switch(status, OK = "v", SUSPICIOUS = "!", PROBLEM = "x", "*"))
}

#' Text with cli inline markup formatted, element by element
#' @noRd
out_markup <- function(text, .envir) {
  vapply(text, cli::format_inline, character(1), .envir = .envir,
         USE.NAMES = FALSE)
}

#' A RiskMap function named in printed output
#'
#' Formatted as cli formats a function, and linked to the help page `topic`
#' where the console supports such links, as RStudio does. Elsewhere a link
#' would be spelled out in full, so it is left out.
#' @noRd
out_fun <- function(name, topic = name) {
  links <- isTRUE(cli::ansi_has_hyperlink_support()) &&
    isTRUE(cli::ansi_hyperlink_types()$help)
  if (links) {
    cli::format_inline("{.help [{.fun {name}}](RiskMap::{topic})}")
  } else {
    cli::format_inline("{.fun {name}}")
  }
}

#' Lines of a table printed by R, such as a matrix or `printCoefmat()`
#' @noRd
out_capture <- function(expr, indent = 0L) {
  lines <- utils::capture.output(expr)
  paste0(strrep(" ", indent), lines)
}

#' Rows of labels and values, the labels padded to one width
#' @noRd
out_pairs <- function(labels, values, indent = 0L) {
  if (!length(labels)) return(character())
  width <- max(nchar(labels))
  paste0(strrep(" ", indent), formatC(labels, width = -width), "  ", values)
}
