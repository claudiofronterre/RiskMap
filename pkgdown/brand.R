# Brand images of the website. Run from the package root after changing the
# sticker, then rebuild the site. The sticker files are those of the sticker's
# own folder: the master PNG, its night version for dark pages and the square
# favicon.
#
#   source("pkgdown/brand.R")
#   make_logo("riskmap-sticker-1200.png")                   # man/figures/logo.png
#   make_logo("riskmap-sticker-night-1200.png",
#             "pkgdown/assets/logo-night.png")              # dark pages and navbar
#   make_card("riskmap-sticker@2x.png")                     # pkgdown/assets/card.png
#   make_favicons("favicon-512.png", "favicon.svg")         # pkgdown/favicon/
#
# Update the brand colours in template.bslib of _pkgdown.yml and at the top of
# pkgdown/extra.css and pkgdown/card.html at the same time.

# The logo, 278 pixels high. A white surround is made transparent by
# flood-filling from the corners, and a transparent one is left as it is.
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
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  magick::image_write(img, out, format = "png")
  invisible(out)
}

# The image shown when a link to the site is shared, 1200 by 630 pixels, drawn
# by headless Chrome from pkgdown/card.html with the sticker beside it.
make_card <- function(sticker, out = "pkgdown/assets/card.png") {
  chrome <- Sys.which(c("google-chrome", "chromium", "chromium-browser"))
  chrome <- chrome[nzchar(chrome)]
  if (length(chrome) == 0) {
    stop("make_card() needs Chrome or Chromium.", call. = FALSE)
  }
  page <- tempfile("card-")
  dir.create(page)
  on.exit(unlink(page, recursive = TRUE), add = TRUE)
  file.copy("pkgdown/card.html", file.path(page, "card.html"))
  file.copy(sticker, file.path(page, "sticker.png"))
  shot <- file.path(page, "card.png")
  status <- system2(chrome[[1]], shQuote(c(
    "--headless=new", "--disable-gpu", "--hide-scrollbars",
    paste0("--user-data-dir=", file.path(page, "profile")),
    "--virtual-time-budget=5000", "--window-size=1200,630",
    paste0("--screenshot=", shot),
    paste0("file://", file.path(page, "card.html"))
  )), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(shot)) {
    stop("Chrome did not draw the card, so ", out, " is unchanged.",
         call. = FALSE)
  }
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  file.copy(shot, out, overwrite = TRUE)
  invisible(out)
}

# The favicons pkgdown links from every page, made from a square favicon at
# least 512 pixels wide and its SVG. The icons for phones' home screens get
# the paper background, because iOS fills transparent corners with black.
make_favicons <- function(png, svg, dir = "pkgdown/favicon",
                          paper = "#F5EFE3", ink = "#1D3557") {
  rlang::check_installed("magick")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  img <- magick::image_read(png)
  on_paper <- function(size) {
    magick::image_background(magick::image_scale(img, paste0(size, "x", size)),
                             paper, flatten = TRUE)
  }
  magick::image_write(magick::image_scale(img, "96x96"),
                      file.path(dir, "favicon-96x96.png"), format = "png")
  ico <- magick::image_join(lapply(c(16, 32, 48), function(size) {
    magick::image_scale(img, paste0(size, "x", size))
  }))
  magick::image_write(ico, file.path(dir, "favicon.ico"), format = "ico")
  magick::image_write(on_paper(180), file.path(dir, "apple-touch-icon.png"),
                      format = "png")
  for (size in c(192, 512)) {
    magick::image_write(on_paper(size),
                        file.path(dir, sprintf("web-app-manifest-%dx%d.png",
                                               size, size)),
                        format = "png")
  }
  file.copy(svg, file.path(dir, "favicon.svg"), overwrite = TRUE)
  icons <- sprintf('    {"src": "web-app-manifest-%1$dx%1$d.png", "sizes": "%1$dx%1$d", "type": "image/png"}',
                   c(192L, 512L))
  writeLines(c(
    "{",
    '  "name": "RiskMap",',
    '  "short_name": "RiskMap",',
    '  "icons": [',
    paste(icons, collapse = ",\n"),
    "  ],",
    sprintf('  "theme_color": "%s",', ink),
    sprintf('  "background_color": "%s",', paper),
    '  "display": "standalone"',
    "}"
  ), file.path(dir, "site.webmanifest"))
  invisible(dir)
}
