# Printed lines joined with single spaces, so that a sentence wrapped under a
# bullet, whose lines are indented, can be matched whole.
squish <- function(lines) {
  gsub("\\s+", " ", paste(lines, collapse = " "))
}
