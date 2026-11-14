# Brand images of the website. Run from the package root after changing the
# sticker, then rebuild the site.
#
#   source("pkgdown/brand.R")
#   make_logo("sticker/base.png")  # man/figures/logo.png
#   make_card()                    # pkgdown/assets/card.png
#   make_favicons()                # pkgdown/favicon/
#
# Update the brand colours in template.bslib of _pkgdown.yml and at the top of
# pkgdown/extra.css and pkgdown/card.html at the same time.

# The logo of the site and the README, 278 pixels high. The white surround of
# the sticker is made transparent by flood-filling from the corners, so a
# sticker that already has a transparent background is left as it is.
make_logo <- function(sticker, out = "man/figures/logo.png") {
  rlang::check_installed("magick")
  img <- magick::image_read(sticker)
  info <- magick::image_info(img)
  corners <- sprintf("+%d+%d",
                     c(1, info$width - 2, 1, info$width - 2),
                     c(1, 1, info$height - 2, info$height - 2))
  for (corner in corners) {
    img <- magick::image_fill(img, "none", point = corner, fuzz = 8)
  }
  img <- magick::image_scale(magick::image_trim(img), "x278")
  magick::image_write(img, out, format = "png")
  invisible(out)
}

# The image shown when a link to the site is shared, 1200 by 630 pixels, drawn
# from pkgdown/card.html by headless Chrome.
make_card <- function(out = "pkgdown/assets/card.png") {
  chrome <- Sys.which(c("google-chrome", "chromium", "chromium-browser"))
  chrome <- chrome[nzchar(chrome)]
  if (length(chrome) == 0) {
    stop("make_card() needs Chrome or Chromium.", call. = FALSE)
  }
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  profile <- tempfile("chrome-")
  shot <- tempfile("card-", fileext = ".png")
  on.exit(unlink(c(profile, shot), recursive = TRUE), add = TRUE)
  status <- system2(chrome[[1]], shQuote(c(
    "--headless=new", "--disable-gpu", "--hide-scrollbars",
    paste0("--user-data-dir=", profile),
    "--virtual-time-budget=5000", "--window-size=1200,630",
    paste0("--screenshot=", shot),
    paste0("file://", normalizePath("pkgdown/card.html"))
  )), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(shot)) {
    stop("Chrome did not draw the card, so ", out, " is unchanged.",
         call. = FALSE)
  }
  file.copy(shot, out, overwrite = TRUE)
  invisible(out)
}

# Favicons made from the logo by pkgdown, through realfavicongenerator.net.
# The web manifest it writes points at the root of the domain and has no
# name, so the site, which lives under /RiskMap/, gets relative paths.
make_favicons <- function(manifest = "pkgdown/favicon/site.webmanifest") {
  pkgdown::build_favicons(overwrite = TRUE)
  lines <- readLines(manifest)
  lines <- sub('"src": "/', '"src": "', lines, fixed = TRUE)
  lines <- sub('"name": ""', '"name": "RiskMap"', lines, fixed = TRUE)
  lines <- sub('"short_name": ""', '"short_name": "RiskMap"', lines,
               fixed = TRUE)
  writeLines(lines, manifest)
  invisible(manifest)
}
