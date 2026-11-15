# The panels that plot() of a check returns, one ggplot each.
drawn_panels <- function(drawn) {
  if (inherits(drawn, "patchwork")) c(drawn$patches$plots, list(drawn)) else
    list(drawn)
}

# The titles of the panels that plot() of a check returns, and those of the
# named panels of the check, to tell which panels it holds.
drawn_titles <- function(drawn) {
  vapply(drawn_panels(drawn), function(panel) ggplot2::get_labs(panel)$title,
         character(1))
}
panel_titles <- function(checked, names) {
  vapply(draw_check_panels(checked)[names],
         function(panel) ggplot2::get_labs(panel)$title, character(1),
         USE.NAMES = FALSE)
}
