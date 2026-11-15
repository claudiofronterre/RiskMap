##' The version of RiskMap that is running, for the stamp of its results
##'
##' A list with `version`, the version of the loaded namespace as a
##' `package_version`, and `sha`, the commit RiskMap was installed from when
##' it was installed from GitHub (the `RemoteSha` field that remotes and pak
##' write in its DESCRIPTION), `NA` otherwise. Every development build is
##' called 2.0.0, so the commit tells them apart. The DESCRIPTION is read from
##' the folder of the namespace, which is the installed package or, under
##' `pkgload::load_all()`, the source, and only once in a session.
##' @noRd
riskmap_version_stamp <- function() {
  if (is.null(version_cache$stamp)) {
    namespace <- asNamespace("RiskMap")
    sha <- tryCatch({
      description <- file.path(getNamespaceInfo(namespace, "path"),
                               "DESCRIPTION")
      fields <- read.dcf(description, fields = "RemoteSha")
      unname(fields[1L, "RemoteSha"])
    }, error = function(error) NA_character_)
    version_cache$stamp <- list(
      version = package_version(unname(getNamespaceVersion(namespace))),
      sha = as.character(sha)
    )
  }
  version_cache$stamp
}

version_cache <- new.env(parent = emptyenv())

##' Stamp a result with the version of RiskMap that made it
##'
##' Sets the attribute `riskmap_version` (`riskmap_version_stamp()`) on a
##' result that may be saved and read again, a fit, a prediction, a profile, a
##' bootstrap, a check, a cross-validation, a simulation assessment, a
##' simulation, a simulation model, a variogram or an extrapolation. It is not
##' printed. An object without it was made before RiskMap 2.0.0.
##' @param x The result, with its class set.
##' @return `x` with the attribute.
##' @noRd
stamp_version <- function(x) {
  attr(x, "riskmap_version") <- riskmap_version_stamp()
  x
}
