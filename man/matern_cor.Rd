% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{matern_cor}
\alias{matern_cor}
\title{Matern correlation function}
\usage{
matern_cor(u, phi, kappa, return_sym_matrix = FALSE)
}
\arguments{
\item{u}{a vector with values of the distances between pairs of data locations}

\item{phi}{value of the scale parameter \eqn{\phi}.}

\item{kappa}{value of the smoothness parameter \eqn{\kappa}.}

\item{return_sym_matrix}{a logical value which indicates whether a symmetric
correlation matrix should be return. By default \code{return_sym_matrix=FALSE}.}
}
\value{
A vector of the same length of \code{u} with the value of the Matern
correlation function for the given distances, if \code{return_sym_matrix=FALSE}.
If \code{return_sym_matrix=TRUE}, a symmetric correlation matrix is returned.
}
\description{
Computes the matern correlation function.
}
\details{
The Matern correlation function is defined as
\deqn{\rho(u; \phi; \kappa) = (2^{\kappa-1})^{-1}(u/\phi)^\kappa K_{\kappa}(u/\phi)}
where \eqn{\phi} and \eqn{\kappa} are the scale and smoothness parameters, and
 and \eqn{K_{\kappa}(\cdot)} denotes the the modified Bessel function of the third
 kind of order \eqn{\kappa}. The parameters \eqn{\phi} and \eqn{\kappa} must be
 positive.
}
\author{
Emanuele Giorgi \email{e.giorgi@lancaster.ac.uk}

Claudio Fronterre \email{c.fronterr@lancaster.ac.uk}
}
