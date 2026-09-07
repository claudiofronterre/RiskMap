##' @title Check for valid binomial values
##'
##' @description
##' Checks that binomial data only consists of zero or positive integers and that if
##' den is provided that all values of y are less than or equal to it.
##' Some tolerance is provided for floating point errors
##'
##' @param y the data to check
##' @param den the denominator
##' @return TRUE if valid, raise an error if not
##' @noRd
check_binomial <- function(y, den){
  tolerance <- sqrt(.Machine$double.eps)
  valid <- all(y >= 0) & all(abs(y - round(y)) < tolerance)
  stopifnot("'y' must only consist of zero or positive integers when 'family' is 'binomial'" = valid)

  if (!is.null(den)){
    valid <- all(den >= y)
    stopifnot("Values of 'den' must be greater than or equal to values of 'y'"= valid)
  }

  invisible(TRUE)
}

#' @title check_data
#' @description
#'
#' Check that the data is an sf or sfc object, with a CRS, only containing points
#' or either polygons or multipolygons. If CRS == 4326 it also checks that the #
#' coordinates are possible (i.e. not latitudes > 90)
#' @param data the data to check
#' @param geometry whether to check that the data contains 'point' (default) or
#' 'polygon' (covering both polygons and multipolygons)
#' @param geometry whether to check that the data is 'sf' (default) or
#' 'sfc' (either sf or sfc)
#' @return TRUE if the data is valid. Raise an error if not.
#' @noRd
#'
check_data <- function(data, geometry = "point", type = "sf"){
  stopifnot("'geometry' must be either 'point' or 'polygon'" = geometry %in% c("point", "polygon"))
  stopifnot("'type' must be either 'sf' or 'sfc'" = type %in% c("sf", "sfc"))

  # extract name passed to function
  data_type <- paste0("'", deparse(substitute(data)), "'")

  geometry_type <- switch(geometry,
                          point = "'POINT'",
                          polygon = "'POLYGON' or 'MULTIPOLYGON'")

  if (type == "sf"){
    if (!inherits(data, "sf")){
      stop(paste(data_type, "must be of class 'sf'"))
    }
  } else {
    if (!inherits(data, c("sf", "sfc"))){
      stop(paste(data_type, "must be of class 'sf' or 'sfc'"))
    }
  }

  if (is.na(sf::st_crs(data))){
    stop(paste(data_type, "must contain a coordinate reference system"))
  }

  all_valid_geometry <- all(grepl(toupper(geometry), sf::st_geometry_type(data)))
  if (!all_valid_geometry){
    stop(paste(data_type, "can only contain", geometry_type, "geometry"))
  }

  if (sf::st_crs(data) == sf::st_crs(4326)){
    tryCatch(
      sf::st_is_longlat(data$geometry),
      warning = function(w) {
        stop(paste(data_type, "contains impossible latitude or longitude values -
             check you have specified the columns correctly when converting the data"))
      }
    )
  }
  invisible(TRUE)
}

#' @title check_positive_integer
#' @description
#'
#' Check that a value is a single, positive integer and error if not
#' @param x the value to check
#' @param name the name of the parameter to return in error messages
#' @return TRUE if the data is valid. Raise an error if not.
#' @noRd
#'
check_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    stop("'", name, "' must be a single positive integer")
  }
  if (x <= 0 || x %% 1 != 0) {
    stop("'", name, "' must be a single positive integer")
  }
  invisible(TRUE)
}

#' @title check_positive_number
#' @description
#'
#' Check that a value is a single, positive number and error if not
#' @param x the value to check
#' @param type the type of value being checked. Defaults to 'starting'
#' @return TRUE if x is valid. Raise an error if not.
#' @noRd
#'
check_positive_number <- function(x, type = "starting ") {
  # extract name, removing any list
  name <- gsub('.*\\[\\["([^"]+)"\\]\\].*', "\\1", deparse(substitute(x)))

  if (!is.numeric(x) || length(x) != 1 || x <= 0 || is.na(x)) {
    stop("The ", type, "value for '", name, "' must be a single positive number")
  }

  invisible(TRUE)
}


#' @title check_crs
#' @description
#'
#' Check that a CRS is valid
#' @param crs the CRS to check
#' @return TRUE if the CRS is valid. Raise an error if not.
#' @noRd
#'
check_crs <- function(crs){
  # extract name passed to function
  variable <- deparse(substitute(crs))
  tryCatch(
    st_crs(crs),
    warning = function(w) {
      stop("The '", variable, "' provided is not a valid CRS")
    },
    error = function(e){
      stop("The '", variable, "' provided is not a valid CRS")
    }
  )
  invisible(TRUE)
}
