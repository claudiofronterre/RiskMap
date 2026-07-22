##' @title Construct a unified MDA penalty specification
##' @param alpha_a description
##' @param alpha_b description
##' @param gamma_type description
##' @param gamma_mean description
##' @param gamma_sd description
##' @param rho_type description
##' @param rho_mean description
##' @param rho_sd description
##' @export
make_penalty <- function(alpha_a    = NULL,
                         alpha_b    = NULL,
                         gamma_type = "lognormal",
                         gamma_mean,
                         gamma_sd,
                         rho_type   = NULL,
                         rho_mean   = NULL,
                         rho_sd     = NULL) {

  if (!gamma_type %in% c("lognormal", "normal", "gamma"))
    stop("gamma_type must be 'lognormal', 'normal', or 'gamma'")

  if (!is.null(rho_type) && !rho_type %in% c("lognormal", "normal", "gamma"))
    stop("rho_type must be 'lognormal', 'normal', or 'gamma'")

  if (!is.null(rho_type) && (is.null(rho_mean) || is.null(rho_sd)))
    stop("rho_mean and rho_sd must be provided when rho_type is specified")

  list(
    alpha_a = alpha_a,
    alpha_b = alpha_b,
    gamma_type   = gamma_type,
    gamma_mean   = gamma_mean,
    gamma_sd     = gamma_sd,
    rho_type     = rho_type,
    rho_mean     = rho_mean,
    rho_sd       = rho_sd
  )
}

##' @title Convert unified penalty to DAST format (list of 6 functions)
##' @param p description
##' @export
penalty_to_dast <- function(p) {

  ## alpha: penalty = -log prior (Beta)
  if (!is.null(p$alpha_a) && !is.null(p$alpha_b)) {
    a <- p$alpha_a; b <- p$alpha_b

    logpi    <- function(alpha) (a-1)*log(alpha) + (b-1)*log(1-alpha)
    logpi_d1 <- function(alpha) (a-1)/alpha - (b-1)/(1-alpha)
    logpi_d2 <- function(alpha) -(a-1)/alpha^2 - (b-1)/(1-alpha)^2

    pn    <- function(alpha) -logpi(alpha)
    pn_d1 <- function(alpha) -logpi_d1(alpha)
    pn_d2 <- function(alpha) -logpi_d2(alpha)
  } else {
    pn <- function(alpha) 0
    pn_d1 <- function(alpha) 0
    pn_d2 <- function(alpha) 0
  }

  ## gamma: penalty = -log prior (lognormal)
  mu <- p$gamma_mean; sd <- p$gamma_sd

  logpi_g    <- function(g) -log(g) - (log(g) - mu)^2 / (2*sd^2)
  logpi_g_d1 <- function(g) -1/g - (log(g) - mu) / (sd^2 * g)
  logpi_g_d2 <- function(g)  1/g^2 - (1 - log(g) + mu) / (sd^2 * g^2)

  pn_g    <- function(g) -logpi_g(g)
  pn_g_d1 <- function(g) -logpi_g_d1(g)
  pn_g_d2 <- function(g) -logpi_g_d2(g)

  list(pn, pn_d1, pn_d2, pn_g, pn_g_d1, pn_g_d2)
}

