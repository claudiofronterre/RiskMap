##' @title Propose a projected CRS for the data
##'
##' @description Proposes the projected coordinate reference system (CRS) that
##' keeps distances between locations closest to their values on the Earth,
##' for computing distances and fitting geostatistical models. A WGS 84 UTM
##' zone is proposed when it is nearly as good as the best projection, and a
##' conformal projection centred on the data otherwise.
##'
##' @param data An `sf` or `sfc` object of points or polygons, in any CRS.
##' @param epsg_tolerance The largest distortion of distances, as a fraction,
##'   for which a UTM zone is proposed. The default of 0.005 accepts a UTM
##'   zone that distorts distances by at most 0.5% over the data.
##' @param messages Whether to report the CRS proposed and its distortion.
##'   Defaults to the option `RiskMap.messages`, or `TRUE` if it is unset.
##'
##' @details No planar map keeps every distance, so a projected CRS stretches
##' or shrinks them by a factor that varies over the map. The distortion of a
##' CRS over the data is the largest relative error of short distances over
##' the region the data cover, measured on the WGS 84 ellipsoid (see
##' [crs_distortion()]). It also bounds the relative error of every distance
##' between two locations of the data.
##'
##' The region is the convex hull of the data, whether points or polygons. Two
##' kinds of CRS are compared.
##' - The WGS 84 UTM zone at the centre of the region, the centre of the
##'   smallest spherical cap enclosing it, and its two neighbours. The
##'   hemisphere is that of the centre, with the Equator counted as north,
##'   and does not change distances.
##' - Four conformal projections centred on the data: oblique stereographic,
##'   transverse Mercator, oblique Mercator and Lambert conformal conic. Their
##'   parameters minimise the largest scale error over the region
##'   (Chebyshev's criterion), and their scale factor balances the smallest
##'   and largest scale about 1. Their centres are rounded to 0.5 degrees and
##'   their scale factors to four decimals, so that data with similar extents
##'   usually get the same CRS.
##'
##' Conformal projections keep angles, so a covariance that is isotropic on
##' the Earth stays isotropic locally in the plane, and only its scale changes
##' slightly over the map.
##'
##' The best UTM zone is proposed when its distortion is at most
##' `epsg_tolerance`, and the best centred projection otherwise. A region about
##' 1,000 km across or less usually gets a UTM zone. The best centred
##' projection keeps distances within 0.5% over every African country. No
##' planar CRS keeps them within 1% over a region that holds a disc about
##' 3,600 km across, such as sub-Saharan Africa, and a warning of class
##' `riskmap_warning_crs` then says so.
##'
##' Functions that compute distances, such as [glgpm()], [variogram()] and
##' [create_grid()], call `propose_crs()` for data in longitude and latitude.
##' When several analyses should share one CRS, for example analyses of
##' different surveys in one country, call `propose_crs()` once on the
##' country's boundary and pass the result as `model_crs` or `distance_crs`.
##'
##' @return The EPSG code of the UTM zone as an integer, or an `sf` `crs`
##' object for a centred projection. Its name begins with "RiskMap" and gives
##' the family, the centre (the standard parallel for the Lambert conformal
##' conic projection), the azimuth of the oblique Mercator projection and the
##' scale factor `k0`, for example "RiskMap oblique stereographic 23E 4.5S k0
##' 0.9952". Either can be passed to [sf::st_transform()] and to the
##' `model_crs` argument of [glgpm()].
##'
##' @references
##' Milnor J (1969). A problem in cartography. *American Mathematical Monthly*
##' 76, 1101-1112.
##'
##' Snyder JP (1987). *Map Projections: A Working Manual*. USGS Professional
##' Paper 1395.
##'
##' @seealso [crs_distortion()] for the distortion of any CRS over the data.
##' @examples
##' data("infect_sma")
##' propose_crs(infect_sma)
##'
##' # A region too wide for a UTM zone gets a centred projection
##' library(sf)
##' wide <- st_as_sf(data.frame(lon = c(12, 31, 12, 31), lat = c(-13, -13, 5, 5)),
##'                  coords = c("lon", "lat"), crs = 4326)
##' crs <- propose_crs(wide)
##' crs$Name
##' crs_distortion(wide, crs)
##' @export
propose_crs <- function(data, epsg_tolerance = 0.005,
                        messages = getOption("RiskMap.messages", TRUE)) {
  check_supplied(data)
  check_data(data, geometry = "any", type = "sfc")
  check_range(epsg_tolerance, min = 0, max = 1, allow_equal = FALSE)
  check_logical(messages)

  choice <- choose_crs(data, epsg_tolerance)
  if (messages) {
    pct <- format_distortion(choice$error)
    label <- choice$label
    if (choice$utm) {
      riskmap_inform("Proposed {label}, which distorts distances by at most {pct} over {.arg data}.")
    } else {
      utm_label <- choice$utm_label
      utm_note <- if (is.finite(choice$utm_error)) {
        utm_pct <- format_distortion(choice$utm_error)
        "The best UTM zone, {utm_label}, distorts them by up to {utm_pct}, more than {.arg epsg_tolerance} allows."
      } else {
        "No UTM zone can be used, because the projection fails at some locations of {.arg data}."
      }
      riskmap_inform(c(
        "Proposed a custom {choice$family} projection, {.val {label}}, which distorts distances by at most {pct} over {.arg data}.",
        "i" = utm_note
      ))
    }
  }
  if (choice$error > crs_distortion_limit) {
    # A large tolerance can accept a UTM zone when a centred projection is
    # better.
    best <- if (choice$utm) choose_crs(data, 0) else choice
    if (best$error > crs_distortion_limit) {
      warn_region_too_wide(best, "data")
    } else {
      label <- choice$label
      pct <- format_distortion(choice$error)
      best_pct <- format_distortion(best$error)
      riskmap_warn(c(
        "{label} distorts distances by up to {pct} over {.arg data}.",
        "i" = "A smaller {.arg epsg_tolerance} gives a projection centred on the data, which distorts them by at most {best_pct}."
      ), class = "riskmap_warning_crs")
    }
  }
  choice$crs
}

##' @title Distortion of distances by a projected CRS
##'
##' @description Measures how much a projected coordinate reference system
##' (CRS) stretches or shrinks distances over the region the data cover.
##'
##' @param data An `sf` or `sfc` object of points or polygons, in any CRS.
##' @param crs The projected CRS to assess, as an EPSG code, a string or an
##'   `sf` `crs` object. With `NULL`, the default, the CRS of `data`, which
##'   must then be projected.
##'
##' @details At each location of the region, the CRS maps a small circle on
##' the WGS 84 ellipsoid to an ellipse whose semi-axes, relative to the
##' circle, are the largest and smallest scale factors there (Tissot's
##' indicatrix). They are computed from the derivatives of the projection by
##' central differences over 10 m. The distortion is the largest of the
##' differences between these scale factors and 1 over the region. It also
##' bounds the relative error of the distance between any two locations of the
##' region, because the straight line between them in the plane and the
##' shortest path between them on the Earth both stay within the region, or
##' very close to it for a region thousands of kilometres across.
##'
##' The region is the convex hull of the data, whether points or polygons, so
##' that it holds the paths between them, for example across the sea between
##' islands. It is represented by up to 2,000 vertices of the data, its
##' boundary and a grid of points inside it.
##'
##' The anisotropy is the largest ratio of the two scale factors at a
##' location, minus 1. It is zero, up to rounding, for a conformal CRS, such as
##' UTM or any CRS from [propose_crs()], so that distances from a location are
##' changed by the same factor in every direction. An equal-area CRS changes
##' them more in some directions than in others.
##'
##' @return A data frame with one row and the columns
##' \describe{
##'   \item{`max_error_pct`}{the distortion of distances, in percent}
##'   \item{`max_anisotropy_pct`}{the anisotropy, in percent}
##'   \item{`min_scale`, `max_scale`}{the smallest and largest scale factors
##'   over the region}
##'   \item{`n_locations`}{the number of locations at which they were
##'   computed}
##' }
##' @seealso [propose_crs()]
##' @examples
##' data("infect_sma")
##' # The UTM zone proposed for the data
##' crs_distortion(infect_sma, crs = propose_crs(infect_sma, messages = FALSE))
##' # Pseudo-Mercator, which is not conformal
##' crs_distortion(infect_sma, crs = 3857)
##' @export
crs_distortion <- function(data, crs = NULL) {
  check_supplied(data)
  check_data(data, geometry = "any", type = "sfc")
  if (is.null(crs)) {
    if (st_is_longlat(data)) {
      abort_input("{.arg data} are in longitude and latitude, so {.arg crs} must be given.", arg = "crs")
    }
    crs <- st_crs(data)
    crs_arg <- "data"
  } else {
    check_crs(crs)
    crs <- st_crs(crs)
    if (isTRUE(crs$IsGeographic)) {
      abort_input("{.arg crs} must be a projected coordinate reference system, not longitude and latitude.", arg = "crs")
    }
    crs_arg <- "crs"
  }
  crs_to_distance_factor(crs, "m", crs_arg)

  region <- crs_region(data)
  measure <- measure_crs(region$lonlat, crs)
  if (is.null(measure)) {
    abort_input(c(
      "The distortion of the coordinate reference system given by {.arg {crs_arg}} can't be computed over {.arg data}.",
      "i" = "Some locations may lie outside the area where it is defined."
    ), arg = crs_arg)
  }
  data.frame(
    max_error_pct = 100 * measure$error,
    max_anisotropy_pct = 100 * measure$anisotropy,
    min_scale = measure$min_scale,
    max_scale = measure$max_scale,
    n_locations = measure$n_locations
  )
}

# The largest distortion of distances, as a fraction, accepted without a
# warning.
crs_distortion_limit <- 0.01

# Whether the distortion of a CRS that the user gave, or that the data carry,
# is checked. RiskMap turns it off around refits of a fitted model, whose
# CRS was checked when the model was fitted (`without_crs_check()`).
crs_check_state <- new.env(parent = emptyenv())
crs_check_state$enabled <- TRUE

#' Evaluate code without checking the distortion of the CRS of the data
#'
#' For refits inside RiskMap, such as those of `assess_prediction()`, which
#' would otherwise repeat for every refit the warning given when the model was
#' fitted.
#' @param code The code to evaluate.
#' @noRd
without_crs_check <- function(code) {
  enabled <- crs_check_state$enabled
  crs_check_state$enabled <- FALSE
  on.exit(crs_check_state$enabled <- enabled, add = TRUE)
  code
}

#' Project data for computing distances
#'
#' The convention of every function that computes distances. A CRS given by
#' the user is used, and must be projected. Data in longitude and latitude
#' without one are projected to the CRS from `propose_crs()`, with a message
#' that reports it and its distortion. Projected data are used as they are.
#' A warning of class `riskmap_warning_crs` is raised when the CRS used
#' distorts distances by more than 1% over the data.
#' @param data An `sf` object.
#' @param crs `NULL` or the CRS given by the user.
#' @param crs_arg The name of the argument that gives `crs`, or `NULL` if the
#'   function has none, as `create_grid()`.
#' @param purpose What the CRS is used for, ending the message, such as
#'   `"for fitting"`.
#' @param data_arg The name of the argument that gives `data`.
#' @param call The call named in input errors.
#' @return `data`, projected.
#' @noRd
project_for_distances <- function(data, crs, crs_arg, purpose,
                                  data_arg = "data",
                                  call = rlang::caller_env()) {
  if (!is.null(crs)) {
    check_crs(crs, arg = crs_arg, call = call)
    data <- st_transform(data, crs = crs)
    if (st_is_longlat(data)) {
      abort_input("{.arg {crs_arg}} must be a projected coordinate reference system, not longitude and latitude.", arg = crs_arg, call = call)
    }
    check_crs_distortion(data, crs_arg, data_arg)
    return(data)
  }
  if (!st_is_longlat(data)) {
    check_crs_distortion(data, NULL, data_arg)
    return(data)
  }

  choice <- choose_crs(data)
  data <- st_transform(data, crs = choice$crs)
  label <- choice$label
  pct <- format_distortion(choice$error)
  target <- if (choice$utm) {
    label
  } else {
    cli::format_inline("a custom {choice$family} projection ({.val {label}})")
  }
  reason <- if (choice$utm) {
    ""
  } else {
    tolerance <- format_distortion(choice$tolerance)
    paste0(", because no UTM zone keeps distances within ", tolerance)
  }
  hint <- if (is.null(crs_arg)) {
    "Project {.arg {data_arg}} yourself to choose the projection."
  } else {
    "Set {.arg {crs_arg}} to choose the projection."
  }
  riskmap_inform(c(
    "{.arg {data_arg}} are in longitude and latitude, so they are projected to {target} {purpose}{reason}. Distances are distorted by at most {pct} over the data.",
    "i" = hint
  ))
  if (choice$error > crs_distortion_limit) {
    warn_region_too_wide(choice, data_arg)
  }
  data
}

#' Warn when a projected CRS distorts distances by more than 1%
#'
#' Used for a CRS that the user gave or that the data carry. Does nothing when
#' the check is turned off, or when the distortion can't be computed, for
#' example for a CRS that can't be related to WGS 84.
#' @param data Projected data.
#' @param crs_arg The argument that gave the CRS, or `NULL` when it is the CRS
#'   of the data.
#' @param data_arg The argument that gave the data.
#' @noRd
check_crs_distortion <- function(data, crs_arg, data_arg) {
  if (!crs_check_state$enabled) return(invisible(NULL))
  measure <- tryCatch({
    region <- crs_region(data)
    measure_crs(region$lonlat, st_crs(data))
  }, error = function(e) NULL)
  if (is.null(measure) || measure$error <= crs_distortion_limit) {
    return(invisible(NULL))
  }

  pct <- format_distortion(measure$error)
  source <- if (is.null(crs_arg)) {
    "The coordinate reference system of {.arg {data_arg}}"
  } else {
    "{.arg {crs_arg}}"
  }
  best <- tryCatch(choose_crs(data), error = function(e) NULL)
  advice <- if (is.null(best)) {
    character()
  } else if (best$error <= crs_distortion_limit) {
    best_label <- if (best$utm) best$label else cli::format_inline("a custom {best$family} projection")
    best_pct <- format_distortion(best$error)
    c("i" = "{.fn propose_crs} gives {best_label}, which distorts them by at most {best_pct}.")
  } else {
    diameter <- format_diameter(best$diameter_km)
    best_pct <- format_distortion(best$error)
    c("i" = "The data are about {diameter} km across, too wide for any planar coordinate reference system to keep distances within 1%. The best, from {.fn propose_crs}, distorts them by up to {best_pct}.")
  }
  riskmap_warn(
    c(paste0(source, " distorts distances by up to {pct} over the data."), advice),
    class = "riskmap_warning_crs"
  )
  invisible(NULL)
}

#' Warn that no planar CRS keeps distances within 1%
#' @param choice The result of `choose_crs()`.
#' @param data_arg The argument that gave the data.
#' @noRd
warn_region_too_wide <- function(choice, data_arg) {
  diameter <- format_diameter(choice$diameter_km)
  pct <- format_distortion(choice$error)
  riskmap_warn(c(
    "{.arg {data_arg}} are about {diameter} km across, too wide for any planar coordinate reference system to keep distances within 1%.",
    "i" = "The best one distorts distances by up to {pct}, so distances, and the scale of the spatial correlation, are approximate."
  ), class = "riskmap_warning_crs")
}

#' A distortion, as a fraction, in percent to two significant digits
#' @noRd
format_distortion <- function(x) {
  paste0(trimws(formatC(signif(100 * x, 2), format = "fg", digits = 2)), "%")
}

#' A distance in kilometres, rounded to 100 km
#' @noRd
format_diameter <- function(x) {
  format(max(100, round(x, -2)), big.mark = ",", scientific = FALSE)
}

#' Choose the projected CRS for the data
#'
#' @param data An `sf` or `sfc` object in any CRS.
#' @param epsg_tolerance The largest distortion for which a UTM zone is
#'   chosen.
#' @return A list with the CRS (`crs`, an integer EPSG code or a `crs`
#'   object), whether it is a UTM zone (`utm`), its label for messages
#'   (`label`), its family (`family`), its distortion (`error`), the best UTM
#'   zone's label and distortion (`utm_label`, `utm_error`), the tolerance
#'   (`tolerance`) and the diameter of the region in kilometres
#'   (`diameter_km`).
#' @noRd
choose_crs <- function(data, epsg_tolerance = 0.005) {
  region <- crs_region(data)
  utm <- best_utm_zone(region)
  choice <- list(
    crs = utm$epsg, utm = TRUE, label = utm$label, family = "UTM",
    error = utm$error, utm_label = utm$label, utm_error = utm$error,
    tolerance = epsg_tolerance, diameter_km = region$diameter_km
  )
  if (is.finite(utm$error) && utm$error <= epsg_tolerance) return(choice)

  centred <- centred_crs_candidates(region)
  errors <- vapply(centred, function(candidate) candidate$error, numeric(1))
  if (!length(centred) || !any(is.finite(errors))) {
    if (is.finite(utm$error)) return(choice)
    abort_input(c(
      "The distortion of distances can't be computed over {.arg data} for any planar coordinate reference system.",
      "i" = "The data may cover too much of the Earth. Project them yourself to use a planar one anyway."
    ), arg = "data")
  }
  best <- centred[[which.min(errors)]]
  if (is.finite(utm$error) && best$error >= utm$error) return(choice)

  choice$crs <- best$crs
  choice$utm <- FALSE
  choice$label <- best$crs$Name
  choice$family <- best$family
  choice$error <- best$error
  choice
}

#' The WGS 84 UTM zone with least distortion over the region
#'
#' Compares the zone at the centre of the region with its two neighbours, in
#' the hemisphere of the centre.
#' @param region The result of `crs_region()`.
#' @return A list with the EPSG code (`epsg`), its label and its distortion
#'   (`error`, infinite if it can't be computed).
#' @noRd
best_utm_zone <- function(region) {
  centre <- region$centre
  zone <- floor(((centre[1] + 180) %% 360) / 6) + 1
  zones <- unique(((zone + c(0, -1, 1) - 1) %% 60) + 1)
  epsg <- as.integer((if (centre[2] >= 0) 32600L else 32700L) + zones)
  errors <- vapply(epsg, function(code) {
    measure <- measure_crs(region$lonlat, st_crs(code))
    if (is.null(measure)) Inf else measure$error
  }, numeric(1))
  best <- which.min(errors)
  code <- epsg[best]
  list(epsg = code,
       label = paste0("EPSG:", code, " (", st_crs(code)$Name, ")"),
       error = errors[best])
}

#' Conformal projections centred on the region
#'
#' Chooses the parameters of each family by Chebyshev's criterion on the
#' sphere, from closed forms of the scale factor, then measures each on the
#' ellipsoid and balances its scale factor.
#' @param region The result of `crs_region()`.
#' @return A list of candidates, each with its `crs` object, `family` and
#'   distortion (`error`). A candidate PROJ can't use is left out.
#' @noRd
centred_crs_candidates <- function(region) {
  lonlat <- region$lonlat
  X <- unit_vectors(lonlat)
  centre <- region$centre
  lon <- lonlat[, 1]
  lat <- lonlat[, 2]
  lon_rel <- ((lon - centre[1] + 180) %% 360) - 180

  candidates <- list()

  # Oblique stereographic, centred on the smallest cap enclosing the region
  candidates$stereographic <- list(
    family = "oblique stereographic",
    proj = "+proj=sterea +lat_0={lat0} +lon_0={lon0}",
    lon0 = round_half(centre[1]), lat0 = round_half(centre[2])
  )

  # Transverse Mercator, with the central meridian that minimises the largest
  # of |cos(lat) sin(lon - lon0)|
  meridian <- optimise_scalar(function(l0) {
    max(abs(cos(lat * pi / 180) * sin((lon_rel - l0) * pi / 180)))
  }, range(lon_rel))
  candidates$transverse <- list(
    family = "transverse Mercator",
    proj = "+proj=tmerc +lat_0={lat0} +lon_0={lon0}",
    lon0 = round_half(centre[1] + meridian), lat0 = round_half(centre[2])
  )

  # Oblique Mercator, along the great circle closest to the whole region
  oblique <- oblique_mercator_parameters(X, centre, meridian)
  if (!is.null(oblique)) {
    candidates$oblique <- list(
      family = "oblique Mercator",
      proj = "+proj=omerc +lat_0={lat0} +lonc={lon0} +alpha={alpha} +gamma={alpha}",
      lon0 = oblique$lon0, lat0 = oblique$lat0, alpha = oblique$alpha
    )
  }

  # Lambert conformal conic, with the standard parallel that minimises the
  # range of the log scale factor, away from the Equator
  phi <- lat * pi / 180
  parallel <- optimise_scalar(function(p1) {
    diff(range(-log(cos(phi)) - sin(p1) * log(tan(pi / 4 + phi / 2))))
  }, range(phi)) * 180 / pi
  if (abs(parallel) > 2) {
    candidates$conic <- list(
      family = "Lambert conformal conic",
      proj = "+proj=lcc +lat_1={lat0} +lat_0={lat0} +lon_0={lon0}",
      lon0 = round_half(centre[1]), lat0 = round_half(parallel)
    )
  }

  out <- lapply(candidates, build_centred_crs, lonlat = lonlat)
  Filter(Negate(is.null), out)
}

#' Build and measure one centred projection
#'
#' Measures the projection with a scale factor of 1, then sets the scale
#' factor that balances the smallest and largest scales about 1, rounded to
#' four decimals, and names the CRS.
#' @param candidate A list with the family, a PROJ string template and its
#'   parameters.
#' @param lonlat The locations of the region.
#' @return A list with the `crs` object, `family` and `error`, or `NULL`.
#' @noRd
build_centred_crs <- function(candidate, lonlat) {
  lon0 <- ((candidate$lon0 + 180) %% 360) - 180
  lat0 <- candidate$lat0
  alpha <- candidate$alpha
  base <- glue_parameters(candidate$proj, lon0, lat0, alpha)
  proj <- function(k0) {
    paste0(base, " +k_0=", sprintf("%.4f", k0),
           " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  }
  measure <- tryCatch(measure_crs(lonlat, st_crs(proj(1))),
                      error = function(e) NULL)
  if (is.null(measure)) return(NULL)
  k0 <- round(2 / (measure$min_scale + measure$max_scale), 4)
  error <- max(abs(k0 * c(measure$min_scale, measure$max_scale) - 1))

  name <- paste("RiskMap", candidate$family, format_lon(lon0),
                format_lat(lat0))
  if (!is.null(alpha)) name <- paste(name, "azimuth", format_half(alpha))
  name <- paste(name, "k0", sprintf("%.4f", k0))
  crs <- tryCatch(named_crs(proj(k0), name), error = function(e) NULL)
  if (is.null(crs)) return(NULL)
  list(crs = crs, family = candidate$family, error = error)
}

#' Fill the parameters of a PROJ string template
#' @noRd
glue_parameters <- function(template, lon0, lat0, alpha = NULL) {
  out <- sub("{lon0}", format_half(lon0), template, fixed = TRUE)
  out <- gsub("{lat0}", format_half(lat0), out, fixed = TRUE)
  if (!is.null(alpha)) {
    out <- gsub("{alpha}", format_half(alpha), out, fixed = TRUE)
  }
  out
}

#' A multiple of half a degree as text, such as 23.5 or 4, whatever the
#' option `digits`
#' @noRd
format_half <- function(x) sub("\\.0$", "", sprintf("%.1f", x))

#' A CRS from a PROJ string, with a name in its WKT
#'
#' The name gives the family, the centre (the standard parallel for the
#' Lambert conformal conic projection), the azimuth of the oblique Mercator
#' projection and the scale factor, so that it identifies the CRS.
#' @noRd
named_crs <- function(proj, name) {
  wkt <- st_crs(proj)$wkt
  wkt <- sub('^PROJCRS\\["unknown"', paste0('PROJCRS["', name, '"'), wkt)
  wkt <- sub('BASEGEOGCRS\\["unknown"', 'BASEGEOGCRS["WGS 84"', wkt)
  st_crs(wkt)
}

#' Round to the nearest half degree
#' @noRd
round_half <- function(x) round(2 * x) / 2

#' A longitude or latitude for the name of a CRS, such as 23.5E or 4.5S
#' @noRd
format_lon <- function(x) {
  x <- ((x + 180) %% 360) - 180
  paste0(format_half(abs(x)), if (x < 0) "W" else "E")
}

#' @rdname format_lon
#' @noRd
format_lat <- function(x) {
  paste0(format_half(abs(x)), if (x < 0) "S" else "N")
}

#' Minimise a function of one variable over an interval
#'
#' Returns the end of an interval of zero length.
#' @noRd
optimise_scalar <- function(f, interval) {
  if (!(diff(interval) > 1e-12)) return(interval[1])
  optimize(f, interval, tol = 1e-8)$minimum
}

#' Unit vectors of locations on the sphere
#' @param lonlat A matrix of longitudes and latitudes in degrees.
#' @noRd
unit_vectors <- function(lonlat) {
  lon <- lonlat[, 1] * pi / 180
  lat <- lonlat[, 2] * pi / 180
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

#' Longitude and latitude in degrees of a vector
#' @noRd
vector_lonlat <- function(v) {
  v <- v / sqrt(sum(v^2))
  c(atan2(v[2], v[1]), asin(max(-1, min(1, v[3])))) * 180 / pi
}

#' The unit vectors east and north at a location
#' @param lonlat Longitude and latitude in degrees.
#' @noRd
local_axes <- function(lonlat) {
  lon <- lonlat[1] * pi / 180
  lat <- lonlat[2] * pi / 180
  list(east = c(-sin(lon), cos(lon), 0),
       north = c(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat)))
}

#' Vector product
#' @noRd
cross_product <- function(u, v) {
  c(u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3],
    u[1] * v[2] - u[2] * v[1])
}

#' Centre of the smallest spherical cap enclosing the locations
#'
#' Maximises the smallest cosine of the angle to the locations, over
#' directions parametrised in the plane tangent at `start`.
#' @param X Unit vectors of the locations.
#' @param start Longitude and latitude of the starting centre.
#' @return Longitude and latitude in degrees.
#' @noRd
smallest_cap_centre <- function(X, start) {
  origin <- unit_vectors(matrix(start, 1))[1, ]
  axes <- local_axes(start)
  direction <- function(p) {
    v <- origin + p[1] * axes$east + p[2] * axes$north
    v / sqrt(sum(v^2))
  }
  objective <- function(p) -min(X %*% direction(p))
  fit <- optim(c(0, 0), objective, control = list(reltol = 1e-12))
  fit <- optim(fit$par, objective, control = list(reltol = 1e-12))
  vector_lonlat(direction(fit$par))
}

#' Parameters of the oblique Mercator projection along the region
#'
#' Finds the great circle whose largest distance to the locations is
#' smallest, from three starts: the plane of least spread of the locations,
#' the central meridian of the transverse Mercator projection and the
#' Equator. Its point closest to the mean of the locations, rounded to half a
#' degree, is the centre. The azimuth of the central line through that centre
#' is then chosen again, and rounded to half a degree, avoiding 0 and 90
#' degrees, which some versions of PROJ reject.
#' @param X Unit vectors of the locations.
#' @param centre Longitude and latitude of the centre of the region.
#' @param meridian The central meridian of the transverse Mercator
#'   projection, relative to `centre`.
#' @return A list with `lon0`, `lat0` and `alpha`, or `NULL`.
#' @noRd
oblique_mercator_parameters <- function(X, centre, meridian) {
  objective <- function(q) max(abs(X %*% unit_vectors(matrix(q, 1))[1, ]))
  least_spread <- eigen(crossprod(X), symmetric = TRUE)$vectors[, 3]
  starts <- list(vector_lonlat(least_spread),
                 c(centre[1] + meridian + 90, 0),
                 c(0, 90))
  fits <- lapply(starts, function(start) {
    optim(start, objective, control = list(reltol = 1e-12))
  })
  values <- vapply(fits, function(fit) fit$value, numeric(1))
  if (!any(is.finite(values))) return(NULL)
  pole <- unit_vectors(matrix(fits[[which.min(values)]]$par, 1))[1, ]

  mean_vector <- colMeans(X)
  on_circle <- mean_vector - sum(mean_vector * pole) * pole
  if (sqrt(sum(on_circle^2)) < 1e-8) return(NULL)
  middle <- vector_lonlat(on_circle)
  centre_lonlat <- c(round_half(middle[1]), round_half(middle[2]))
  if (abs(centre_lonlat[2]) >= 89.5) return(NULL)

  origin <- unit_vectors(matrix(centre_lonlat, 1))[1, ]
  axes <- local_axes(centre_lonlat)
  tangent <- cross_product(pole, unit_vectors(matrix(middle, 1))[1, ])
  alpha0 <- atan2(sum(tangent * axes$east), sum(tangent * axes$north))
  line_pole <- function(alpha) {
    cross_product(origin, sin(alpha) * axes$east + cos(alpha) * axes$north)
  }
  width <- function(alpha) max(abs(X %*% line_pole(alpha)))
  alpha <- optimise_scalar(width, alpha0 + c(-1, 1) * 10 * pi / 180)
  alpha <- ((alpha * 180 / pi + 90) %% 180) - 90
  alpha <- round_half(alpha)
  if (alpha == 0) alpha <- 0.5
  if (abs(alpha) == 90) alpha <- 89.5
  list(lon0 = centre_lonlat[1], lat0 = centre_lonlat[2], alpha = alpha)
}

#' Locations representing the region covered by the data
#'
#' The data are projected to an azimuthal equidistant projection centred on
#' them, where the region is the convex hull of the data, whether points or
#' polygons. It is represented by up to `max_points` of the vertices of the
#' data, the boundary of the convex hull, densified, and a grid of up to 1,600
#' points inside the hull.
#' @param data An `sf` or `sfc` object in any CRS.
#' @param max_points The largest number of vertices of the data used.
#' @return A list with the matrix of longitudes and latitudes of the
#'   locations (`lonlat`), the centre of the region (`centre`), which is the
#'   centre of the smallest spherical cap enclosing it, and its diameter in
#'   kilometres (`diameter_km`).
#' @noRd
crs_region <- function(data, max_points = 2000L) {
  axis_order <- st_axis_order(FALSE)
  on.exit(st_axis_order(axis_order), add = TRUE)
  geometry <- st_geometry(data)
  geometry <- geometry[!st_is_empty(geometry)]
  if (!length(geometry)) {
    abort_input("{.arg data} must have at least one geometry that is not empty.", arg = "data")
  }
  lonlat <- st_transform(geometry, 4326)
  vertices <- unique(st_coordinates(lonlat)[, 1:2, drop = FALSE])
  mean_vector <- colMeans(unit_vectors(vertices))
  centre <- if (sqrt(sum(mean_vector^2)) > 1e-8) {
    vector_lonlat(mean_vector)
  } else {
    vertices[1, ]
  }

  plane_crs <- st_crs(sprintf(
    "+proj=aeqd +lat_0=%.10f +lon_0=%.10f +datum=WGS84 +units=m +no_defs",
    centre[2], centre[1]
  ))
  plane <- st_transform(lonlat, plane_crs)
  xy <- unique(st_coordinates(plane)[, 1:2, drop = FALSE])
  xy <- xy[is.finite(rowSums(xy)), , drop = FALSE]
  if (nrow(xy) > max_points) {
    xy <- xy[unique(round(seq(1, nrow(xy), length.out = max_points))), ,
             drop = FALSE]
  }

  hull <- st_convex_hull(st_combine(plane))
  hull_xy <- st_coordinates(hull)[, 1:2, drop = FALSE]
  extent <- sqrt(sum(apply(hull_xy, 2, function(x) diff(range(x)))^2))
  boundary_xy <- hull_xy
  if (extent > 0) {
    boundary <- st_segmentize(hull, dfMaxLength = extent / 100)
    boundary_xy <- st_coordinates(boundary)[, 1:2, drop = FALSE]
  }

  inner_xy <- NULL
  if (any(st_dimension(hull) == 2L)) {
    inner_xy <- tryCatch({
      grid <- st_make_grid(hull, n = c(40L, 40L), what = "centers")
      inside <- lengths(st_intersects(grid, hull)) > 0
      st_coordinates(grid[inside])[, 1:2, drop = FALSE]
    }, error = function(e) NULL)
  }

  points <- rbind(xy, boundary_xy, inner_xy)
  out <- sf_project(plane_crs$wkt, "EPSG:4326", points, keep = TRUE,
                    warn = FALSE, authority_compliant = FALSE)
  out <- unique(out[is.finite(rowSums(out)), , drop = FALSE])
  out[, 1] <- ((out[, 1] + 180) %% 360) - 180
  out[, 2] <- pmin(pmax(out[, 2], -89.9), 89.9)

  centre <- smallest_cap_centre(unit_vectors(out), centre)

  corners <- sf_project(plane_crs$wkt, "EPSG:4326", unique(hull_xy),
                        keep = TRUE, warn = FALSE,
                        authority_compliant = FALSE)
  corners <- corners[is.finite(rowSums(corners)), , drop = FALSE]
  cosines <- tcrossprod(unit_vectors(corners))
  diameter_km <- 6371.0088 * acos(max(-1, min(1, min(cosines))))

  list(lonlat = out, centre = as.numeric(centre), diameter_km = diameter_km)
}

#' Scale factors of a CRS at locations
#'
#' Tissot's semi-axes from the derivatives of the projection with respect to
#' distances on the WGS 84 ellipsoid, by central differences over 10 m.
#' @param lonlat A matrix of longitudes and latitudes in degrees.
#' @param crs A projected `crs` object.
#' @return A matrix with the largest (`a`) and smallest (`b`) scale factor at
#'   each location. Longitude comes first and the projected coordinates are
#'   read as easting and northing, whatever `sf::st_axis_order()` says.
#' @noRd
tissot_scales <- function(lonlat, crs) {
  semi_major <- 6378137
  flattening <- 1 / 298.257223563
  e2 <- flattening * (2 - flattening)
  step <- 10
  lon <- lonlat[, 1]
  lat <- lonlat[, 2]
  phi <- lat * pi / 180
  w <- 1 - e2 * sin(phi)^2
  meridional <- semi_major * (1 - e2) / w^1.5
  normal <- semi_major / sqrt(w)
  dlon <- step / (normal * cos(phi)) * 180 / pi
  dlat <- step / meridional * 180 / pi
  points <- rbind(cbind(lon + dlon, lat), cbind(lon - dlon, lat),
                  cbind(lon, lat + dlat), cbind(lon, lat - dlat))
  projected <- sf_project("EPSG:4326", crs$wkt, points, keep = TRUE,
                          warn = FALSE, authority_compliant = FALSE)
  projected <- projected * crs_to_distance_factor(crs, "m", "crs")
  n <- length(lon)
  east <- (projected[seq_len(n), , drop = FALSE] -
             projected[n + seq_len(n), , drop = FALSE]) / (2 * step)
  north <- (projected[2 * n + seq_len(n), , drop = FALSE] -
              projected[3 * n + seq_len(n), , drop = FALSE]) / (2 * step)
  sum_squares <- rowSums(east^2) + rowSums(north^2)
  determinant <- abs(east[, 1] * north[, 2] - east[, 2] * north[, 1])
  plus <- sqrt(sum_squares + 2 * determinant)
  minus <- sqrt(pmax(sum_squares - 2 * determinant, 0))
  cbind(a = (plus + minus) / 2, b = (plus - minus) / 2)
}

#' Distortion of a CRS over the locations of a region
#' @param lonlat A matrix of longitudes and latitudes in degrees.
#' @param crs A projected `crs` object.
#' @return A list with the distortion (`error`), the anisotropy, the smallest
#'   and largest scale factors and the number of locations, or `NULL` if the
#'   scale factors can't be computed at every location.
#' @noRd
measure_crs <- function(lonlat, crs) {
  scales <- tissot_scales(lonlat, crs)
  if (!nrow(scales) || !all(is.finite(scales)) || any(scales[, "b"] <= 0)) {
    return(NULL)
  }
  low <- min(scales)
  high <- max(scales)
  list(error = max(abs(c(low, high) - 1)),
       anisotropy = max(scales[, "a"] / scales[, "b"] - 1),
       min_scale = low, max_scale = high,
       n_locations = nrow(scales))
}
