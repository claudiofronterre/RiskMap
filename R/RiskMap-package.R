#' @keywords internal
"_PACKAGE"

##' @keywords internal
##' @import stats
##' @import sf
##' @import ggplot2
##' @importFrom utils flush.console
##' @importFrom graphics points
##' @importFrom Deriv Deriv
##' @importFrom numDeriv grad hessian
##' @importFrom rlang sym
##' @importFrom gridExtra grid.arrange
##' @importFrom dplyr mutate first
##' @importFrom tibble as_tibble
NULL

utils::globalVariables(c("lo", "hi", "series", "obs", "n_obs"))
