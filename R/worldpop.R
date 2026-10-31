##' Download gridded population counts from WorldPop
##'
##' Finds and downloads WorldPop population counts for a country and year,
##' for the whole population or for chosen age groups and sexes, and returns
##' them as a single raster of people per cell.
##'
##' @param country ISO3 code of the country, for example `"LBR"` for Liberia.
##' @param year Year of the population estimates.
##' @param ages `NULL` (the default) for the whole population, or a vector
##'   `c(from, to)` of ages in years, for example `c(20, Inf)` for adults aged
##'   20 and over or `c(0, 4)` for children under five. Both ends must match
##'   the WorldPop age groups (0, 1 to 4, 5 to 9, and so on in five-year
##'   groups, with an open-ended oldest group).
##' @param sex `"both"` (the default), `"female"` or `"male"`. Only used with
##'   `ages`.
##' @param resolution `"100m"` (the default) or `"1km"`.
##' @param constrained Logical. `FALSE` (the default) uses the unconstrained
##'   estimates, which spread people over all land; `TRUE` uses the
##'   constrained estimates, which place people only in mapped settlements.
##'   See Details.
##' @param un_adjusted Logical. If `TRUE`, national totals are adjusted to
##'   match United Nations population estimates. Available for the whole
##'   population of the unconstrained 2000 to 2020 estimates.
##' @param dest_dir Directory for the downloaded files. Files already there
##'   are not downloaded again. Defaults to a temporary directory.
##' @param messages Logical. Report the product and the downloads.
##'
##' @details
##' WorldPop distributes census-based population totals over a fine grid with
##' statistical models that use covariates such as settlement maps, roads and
##' night-time lights. The products differ in how people are placed:
##'
##' * Unconstrained estimates (`constrained = FALSE`) place people in every
##'   land cell, so sparsely populated or uninhabited areas receive small
##'   counts. They are available for every year from 2000 to 2020, at 100 m
##'   and 1 km, for the whole population and by age and sex. They suit
##'   population weights on grids of a few kilometres and analyses of past
##'   surveys.
##' * Constrained estimates (`constrained = TRUE`) place people only in cells
##'   mapped as built-up settlements. They describe where people live more
##'   sharply, particularly at fine resolution, but can miss small rural
##'   settlements. This function uses the WorldPop release covering 2015 to
##'   2030.
##'
##' WorldPop also produces bottom-up estimates from microcensus surveys for a
##' few countries, which this function does not cover.
##'
##' The function queries the WorldPop catalogue at
##' <https://hub.worldpop.org/rest/data>. The 1 km age and sex layers of the
##' 2000 to 2020 release are not listed in the catalogue, and their addresses
##' are built from WorldPop's naming pattern. When `ages` covers several age
##' groups or both sexes, their counts are added together.
##'
##' The result is in longitude and latitude. To aggregate the counts to the
##' cells of a prediction grid, project them with [terra::project()] and
##' `method = "sum"`, sum them into larger cells with [terra::aggregate()],
##' and pass the result to [create_grid()] as `raster`.
##'
##' @return A `SpatRaster` with one layer, `population`, giving the number of
##'   people per cell.
##' @references WorldPop (www.worldpop.org). Please cite the WorldPop
##'   products used, as listed on the WorldPop website.
##' @seealso [create_grid()]
##' @examples
##' \donttest{
##' # Adults aged 20 and over in Liberia in 2011, at 1 km
##' adults <- get_worldpop("LBR", year = 2011, ages = c(20, Inf),
##'                        resolution = "1km")
##' terra::global(adults, "sum", na.rm = TRUE)
##' }
##' @export
get_worldpop <- function(country, year, ages = NULL,
                         sex = c("both", "female", "male"),
                         resolution = c("100m", "1km"),
                         constrained = FALSE, un_adjusted = FALSE,
                         dest_dir = tempdir(), messages = TRUE) {
  if (!is.character(country) || length(country) != 1L ||
      !grepl("^[A-Za-z]{3}$", country)) {
    stop("'country' must be a three-letter ISO3 code, such as \"LBR\"",
         call. = FALSE)
  }
  iso3 <- toupper(country)
  if (!is.numeric(year) || length(year) != 1L || year != round(year)) {
    stop("'year' must be a single whole number", call. = FALSE)
  }
  sex <- match.arg(sex)
  resolution <- match.arg(resolution)
  check_logical(constrained)
  check_logical(un_adjusted)
  check_logical(messages)
  if (!is.null(ages) && (!is.numeric(ages) || length(ages) != 2L ||
                         anyNA(ages) || ages[1] > ages[2])) {
    stop("'ages' must be NULL or c(from, to), for example c(20, Inf)",
         call. = FALSE)
  }

  product <- worldpop_product(year, ages, resolution, constrained,
                              un_adjusted)
  urls <- worldpop_files(product, iso3, year, ages, sex)
  if (messages) {
    message("WorldPop ", product$label, ", ", iso3, " ", year,
            if (!is.null(ages)) worldpop_age_label(ages, sex), ": ",
            length(urls), " file", if (length(urls) > 1L) "s")
  }
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  paths <- worldpop_download(urls, dest_dir, messages)
  population <- terra::rast(paths)
  if (terra::nlyr(population) > 1L) population <- sum(population)
  names(population) <- "population"
  population
}

#' The WorldPop product for a request
#' @noRd
worldpop_product <- function(year, ages, resolution, constrained,
                             un_adjusted) {
  by_age <- !is.null(ages)
  if (constrained) {
    if (year < 2015 || year > 2030) {
      stop("Constrained estimates are available for 2015 to 2030; use ",
           "constrained = FALSE for earlier years", call. = FALSE)
    }
    if (un_adjusted) {
      stop("'un_adjusted' is available only for the unconstrained ",
           "2000 to 2020 estimates of the whole population", call. = FALSE)
    }
    dataset <- if (by_age) {
      paste0("age_structures/G2_CN_Age_R25A_", resolution)
    } else {
      paste0("pop/G2_CN_POP_R25A_", resolution)
    }
    return(list(dataset = dataset, release = "2015-2030",
                label = paste("constrained 2015-2030 estimates at",
                              resolution)))
  }
  if (year < 2000 || year > 2020) {
    stop("Unconstrained estimates are available for 2000 to 2020; use ",
         "constrained = TRUE for 2015 to 2030", call. = FALSE)
  }
  if (un_adjusted && by_age) {
    stop("'un_adjusted' is available only for the whole population",
         call. = FALSE)
  }
  dataset <- if (by_age) {
    if (resolution == "1km") NA_character_ else "age_structures/aswpgp"
  } else if (resolution == "1km") {
    if (un_adjusted) "pop/wpicuadj1km" else "pop/wpic1km"
  } else {
    if (un_adjusted) "pop/wpgpunadj" else "pop/wpgp"
  }
  list(dataset = dataset, release = "2000-2020",
       label = paste0("unconstrained 2000-2020 estimates",
                      if (un_adjusted) " adjusted to UN totals",
                      " at ", resolution))
}

#' Addresses of the files for a request
#' @noRd
worldpop_files <- function(product, iso3, year, ages, sex) {
  if (is.na(product$dataset)) {
    # The 1 km age and sex layers of 2000-2020 are not in the catalogue.
    groups <- c(0, 1, seq(5, 80, by = 5))
    keep <- worldpop_age_groups(groups, ages)
    sexes <- switch(sex, both = c("f", "m"), female = "f", male = "m")
    files <- sprintf("%s_%s_%d_%d_1km.tif", tolower(iso3),
                     rep(sexes, each = length(keep)), keep, year)
    return(paste0("https://data.worldpop.org/GIS/AgeSex_structures/",
                  "Global_2000_2020_1km/unconstrained/", year, "/", iso3,
                  "/", files))
  }
  catalogue <- worldpop_catalogue(product$dataset, iso3)
  entry <- catalogue$files[as.character(catalogue$popyear) ==
                             as.character(year)]
  if (!length(entry)) {
    stop("WorldPop has no ", product$label, " for ", iso3, " in ", year,
         call. = FALSE)
  }
  files <- grep("\\.tif$", entry[[1]], value = TRUE)
  if (is.null(ages)) {
    if (length(files) != 1L) {
      stop("Unexpected WorldPop catalogue entry for ", iso3, " in ", year,
           call. = FALSE)
    }
    return(files)
  }
  pattern <- "^[a-z]{3}_([a-z])_([0-9]+)_[0-9]{4}"
  names <- basename(files)
  parsed <- regmatches(names, regexec(pattern, names))
  ok <- lengths(parsed) == 3L
  files <- files[ok]
  file_sex <- vapply(parsed[ok], `[`, "", 2)
  file_age <- as.numeric(vapply(parsed[ok], `[`, "", 3))
  keep <- worldpop_age_groups(sort(unique(file_age)), ages)
  sexes <- if (sex == "both") {
    if (any(file_sex == "t")) "t" else c("f", "m")
  } else {
    substr(sex, 1, 1)
  }
  files[file_sex %in% sexes & file_age %in% keep]
}

#' Lower bounds of the age groups between two ages
#' @noRd
worldpop_age_groups <- function(groups, ages) {
  upper <- c(groups[-1] - 1, Inf)
  if (!ages[1] %in% groups || !ages[2] %in% upper) {
    stop("'ages' must start and end on WorldPop age groups: groups start ",
         "at ", paste(groups, collapse = ", "), " and the oldest is ",
         "open-ended, so use, for example, c(20, Inf) or c(15, 49)",
         call. = FALSE)
  }
  groups[groups >= ages[1] & upper <= ages[2]]
}

#' Describe the age and sex selection
#' @noRd
worldpop_age_label <- function(ages, sex) {
  paste0(", ", if (sex == "both") "" else paste0(sex, "s "), "aged ",
         ages[1], if (is.finite(ages[2])) paste0(" to ", ages[2]) else
           " and over")
}

#' Query the WorldPop catalogue
#' @noRd
worldpop_catalogue <- function(dataset, iso3) {
  url <- paste0("https://hub.worldpop.org/rest/data/", dataset, "?iso3=",
                iso3)
  result <- tryCatch(jsonlite::fromJSON(url), error = function(error) {
    stop("Could not reach the WorldPop catalogue (", conditionMessage(error),
         ")", call. = FALSE)
  })
  data <- result$data
  if (is.null(data) || !length(data)) {
    stop("WorldPop has no data for ", iso3, " in ", dataset, call. = FALSE)
  }
  data
}

#' Download files that are not already present
#' @noRd
worldpop_download <- function(urls, dest_dir, messages) {
  old_timeout <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old_timeout), add = TRUE)
  paths <- file.path(dest_dir, basename(urls))
  for (i in seq_along(urls)) {
    if (file.exists(paths[i]) && file.size(paths[i]) > 0) next
    status <- tryCatch(
      utils::download.file(urls[i], paths[i], mode = "wb", quiet = TRUE),
      error = function(error) error
    )
    if (inherits(status, "error") || !identical(status, 0L)) {
      unlink(paths[i])
      stop("Could not download ", urls[i], call. = FALSE)
    }
  }
  paths
}
