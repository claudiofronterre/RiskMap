#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <algorithm>
#include <vector>
#include "matern_kernel.h"

#ifndef FCONE
#define FCONE
#endif

using namespace Rcpp;

// Dense linear algebra of the adaptive MCML proposal. Each function replaces a
// sequence of R calls that allocated several latent-dimension square matrices
// on the R heap at every evaluation of the Monte Carlo log-likelihood, whose
// garbage collection took about half of a fit. Working memory here is on the
// C++ heap. The arithmetic follows the R code step by step, with the same
// LAPACK and BLAS routines R calls and the same order of summation, so the
// results are identical to that code.

namespace {

// chol(): upper factor by dpotrf, lower triangle set to zero. Returns false
// if the matrix is not positive definite.
bool factor_upper(std::vector<double>& matrix, const int dimension) {
  int info = 0;
  F77_CALL(dpotrf)("U", &dimension, matrix.data(), &dimension, &info FCONE);
  if (info != 0) return false;
  for (R_xlen_t column = 0; column < dimension; ++column) {
    for (R_xlen_t row = column + 1; row < dimension; ++row) {
      matrix[row + column * dimension] = 0.0;
    }
  }
  return true;
}

// chol2inv(): the inverse from an upper factor by dpotri, with the lower
// triangle filled from the upper one.
bool invert_from_factor(std::vector<double>& matrix, const int dimension) {
  int info = 0;
  F77_CALL(dpotri)("U", &dimension, matrix.data(), &dimension, &info FCONE);
  if (info != 0) return false;
  for (R_xlen_t column = 0; column < dimension; ++column) {
    for (R_xlen_t row = column + 1; row < dimension; ++row) {
      matrix[row + column * dimension] = matrix[column + row * dimension];
    }
  }
  return true;
}

// A' diag(C) A added to the prior precision P, as latent_data_hessian() and
// `+` compute it: the curvature is accumulated pair by pair of an
// observation's latent columns, observation by observation, then added.
std::vector<double> add_data_hessian(const NumericMatrix& precision,
                                     const NumericVector& curvature,
                                     const IntegerMatrix& columns) {
  const int dimension = precision.nrow();
  const R_xlen_t size = static_cast<R_xlen_t>(dimension) * dimension;
  std::vector<double> data_hessian(size, 0.0);
  for (int first = 0; first < columns.ncol(); ++first) {
    for (int second = 0; second < columns.ncol(); ++second) {
      for (int observation = 0; observation < columns.nrow(); ++observation) {
        const R_xlen_t row = columns(observation, first) - 1;
        const R_xlen_t column = columns(observation, second) - 1;
        data_hessian[row + column * dimension] += curvature[observation];
      }
    }
  }
  for (R_xlen_t index = 0; index < size; ++index) {
    data_hessian[index] = precision[index] + data_hessian[index];
  }
  return data_hessian;
}

// latent_gradient(): A'b, one latent column of each observation at a time.
std::vector<double> compute_latent_gradient(const NumericVector& linear,
                                            const IntegerMatrix& columns,
                                            const int dimension) {
  std::vector<double> gradient(dimension, 0.0);
  std::vector<double> column_sums(dimension);
  for (int column = 0; column < columns.ncol(); ++column) {
    std::fill(column_sums.begin(), column_sums.end(), 0.0);
    for (int observation = 0; observation < columns.nrow(); ++observation) {
      column_sums[columns(observation, column) - 1] += linear[observation];
    }
    for (int latent = 0; latent < dimension; ++latent) {
      gradient[latent] += column_sums[latent];
    }
  }
  return gradient;
}

// inverse_cholesky_factor(): with J the exchange matrix, chol(J H J) = R
// gives the upper factor J R^-T J of H^-1. R^-1 is found as backsolve() does,
// by dtrsm on the identity. Writes the factor to `factor` and returns false
// if H is not positive definite.
bool compute_inverse_cholesky_factor(const double* hessian,
                                     const int dimension, double* factor) {
  const R_xlen_t size = static_cast<R_xlen_t>(dimension) * dimension;
  std::vector<double> reversed(size);
  for (R_xlen_t column = 0; column < dimension; ++column) {
    for (R_xlen_t row = 0; row < dimension; ++row) {
      reversed[row + column * dimension] =
        hessian[(dimension - 1 - row) + (dimension - 1 - column) * dimension];
    }
  }
  if (!factor_upper(reversed, dimension)) return false;
  std::vector<double> inverse_root(size, 0.0);
  for (R_xlen_t index = 0; index < dimension; ++index) {
    inverse_root[index + index * dimension] = 1.0;
  }
  const double one = 1.0;
  F77_CALL(dtrsm)("L", "U", "N", "N", &dimension, &dimension, &one,
                  reversed.data(), &dimension, inverse_root.data(), &dimension
                  FCONE FCONE FCONE FCONE);
  for (R_xlen_t column = 0; column < dimension; ++column) {
    for (R_xlen_t row = 0; row < dimension; ++row) {
      factor[row + column * dimension] =
        inverse_root[(dimension - 1 - column) + (dimension - 1 - row) * dimension];
    }
  }
  return true;
}

} // namespace

// Marginal means and variances of the observations' latent linear predictors
// under the global proposal N(Q^-1 A'b, Q^-1), Q = P + A' diag(C) A, for the
// per-observation factors b (linear) and C (curvature). `columns` holds, for
// each observation, the 1-based latent columns that A adds up. This is the
// dense part of each refit of the global iteration (chol(), chol2inv(), `%*%`
// with a vector by dgemv, latent_linear_predictor() and rowsum() of the
// inverse's entries). Returns NULL when Q is not positive definite.
// [[Rcpp::export]]
SEXP cpp_global_marginals(const NumericMatrix& precision,
                          const NumericVector& linear,
                          const NumericVector& curvature,
                          const IntegerMatrix& columns) {
  const int dimension = precision.nrow();
  const int observations = columns.nrow();
  const int width = columns.ncol();
  std::vector<double> inverse = add_data_hessian(precision, curvature,
                                                 columns);
  if (!factor_upper(inverse, dimension)) return R_NilValue;
  if (!invert_from_factor(inverse, dimension)) return R_NilValue;

  std::vector<double> gradient = compute_latent_gradient(linear, columns,
                                                         dimension);
  std::vector<double> centre(dimension);
  const double one = 1.0;
  const double zero = 0.0;
  const int increment = 1;
  F77_CALL(dgemv)("N", &dimension, &dimension, &one, inverse.data(),
                  &dimension, gradient.data(), &increment, &zero,
                  centre.data(), &increment FCONE);

  NumericVector marginal_mean(observations);
  NumericVector marginal_variance(observations);
  for (int observation = 0; observation < observations; ++observation) {
    marginal_mean[observation] = centre[columns(observation, 0) - 1];
  }
  for (int column = 1; column < width; ++column) {
    for (int observation = 0; observation < observations; ++observation) {
      marginal_mean[observation] += centre[columns(observation, column) - 1];
    }
  }
  for (int first = 0; first < width; ++first) {
    for (int second = 0; second < width; ++second) {
      for (int observation = 0; observation < observations; ++observation) {
        const R_xlen_t row = columns(observation, first) - 1;
        const R_xlen_t column = columns(observation, second) - 1;
        marginal_variance[observation] += inverse[row + column * dimension];
      }
    }
  }
  return List::create(Named("marginal_mean") = marginal_mean,
                      Named("marginal_variance") = marginal_variance);
}

// The converged global proposal: the precision Q = P + A' diag(C) A made
// exactly symmetric, the centre Q^-1 A'b by two triangular solves with the
// Cholesky factor of Q, and the upper factor of Q^-1 that maps the
// standard-normal draws (`root`, NULL if that factorisation fails). Returns
// NULL when Q is not positive definite.
// [[Rcpp::export]]
SEXP cpp_global_proposal_factors(const NumericMatrix& precision,
                                 const NumericVector& linear,
                                 const NumericVector& curvature,
                                 const IntegerMatrix& columns) {
  const int dimension = precision.nrow();
  const R_xlen_t size = static_cast<R_xlen_t>(dimension) * dimension;
  std::vector<double> sum = add_data_hessian(precision, curvature, columns);
  NumericMatrix symmetric(dimension, dimension);
  for (R_xlen_t column = 0; column < dimension; ++column) {
    for (R_xlen_t row = 0; row < dimension; ++row) {
      symmetric[row + column * dimension] = (sum[row + column * dimension] +
        sum[column + row * dimension]) / 2;
    }
  }
  std::vector<double> factor(symmetric.begin(), symmetric.end());
  if (!factor_upper(factor, dimension)) return R_NilValue;

  // backsolve(root, forwardsolve(t(root), gradient)), with the transposed
  // factor as forwardsolve() receives it.
  std::vector<double> lower(size);
  for (R_xlen_t column = 0; column < dimension; ++column) {
    for (R_xlen_t row = 0; row < dimension; ++row) {
      lower[row + column * dimension] = factor[column + row * dimension];
    }
  }
  std::vector<double> gradient = compute_latent_gradient(linear, columns,
                                                         dimension);
  NumericVector centre(gradient.begin(), gradient.end());
  const double one = 1.0;
  const int one_column = 1;
  F77_CALL(dtrsm)("L", "L", "N", "N", &dimension, &one_column, &one,
                  lower.data(), &dimension, centre.begin(), &dimension
                  FCONE FCONE FCONE FCONE);
  F77_CALL(dtrsm)("L", "U", "N", "N", &dimension, &one_column, &one,
                  factor.data(), &dimension, centre.begin(), &dimension
                  FCONE FCONE FCONE FCONE);

  NumericMatrix root(dimension, dimension);
  const bool factored = compute_inverse_cholesky_factor(
    symmetric.begin(), dimension, root.begin()
  );
  return List::create(Named("precision") = symmetric,
                      Named("centre") = centre,
                      Named("root") = factored ? SEXP(root) : R_NilValue);
}

// Upper factor of the inverse of a positive-definite matrix, as
// inverse_cholesky_factor() defines it. Returns NULL if the matrix is not
// positive definite.
// [[Rcpp::export]]
SEXP cpp_inverse_cholesky_factor(const NumericMatrix& hessian) {
  const int dimension = hessian.nrow();
  NumericMatrix factor(dimension, dimension);
  if (!compute_inverse_cholesky_factor(hessian.begin(), dimension,
                                       factor.begin())) {
    return R_NilValue;
  }
  return factor;
}

// The spatial prior at one parameter value for kappa = 0.5, 1.5 or 2.5: the
// covariance sigma2 (R(phi) + nugget I), its upper Cholesky factor (with the
// attribute jitter = 0, as factor_covariance() sets it) and the precision
// from that factor (chol2inv()), symmetrised first as factor_covariance()
// does. If the covariance is not positive definite, `root` and `precision`
// are NULL and the caller factorises the covariance with jitter.
// [[Rcpp::export]]
List cpp_spatial_prior_factors(const NumericVector& distances,
                               const int dimension, const double phi,
                               const double kappa, const double sigma2,
                               const double nugget_ratio) {
  const R_xlen_t size = static_cast<R_xlen_t>(dimension) * dimension;
  NumericMatrix covariance(dimension, dimension);
  for (R_xlen_t index = 0; index < size; ++index) {
    covariance[index] = half_integer_matern(distances[index], phi, kappa);
  }
  for (R_xlen_t index = 0; index < dimension; ++index) {
    covariance[index + index * dimension] += nugget_ratio;
  }
  for (R_xlen_t index = 0; index < size; ++index) {
    covariance[index] = sigma2 * covariance[index];
  }
  std::vector<double> factor(size);
  for (R_xlen_t column = 0; column < dimension; ++column) {
    for (R_xlen_t row = 0; row < dimension; ++row) {
      factor[row + column * dimension] = (covariance[row + column * dimension] +
        covariance[column + row * dimension]) / 2;
    }
  }
  if (!factor_upper(factor, dimension)) {
    return List::create(Named("covariance") = covariance,
                        Named("root") = R_NilValue,
                        Named("precision") = R_NilValue);
  }
  NumericMatrix root(dimension, dimension);
  std::copy(factor.begin(), factor.end(), root.begin());
  root.attr("jitter") = 0.0;
  invert_from_factor(factor, dimension);
  NumericMatrix precision(dimension, dimension);
  std::copy(factor.begin(), factor.end(), precision.begin());
  return List::create(Named("covariance") = covariance,
                      Named("root") = root,
                      Named("precision") = precision);
}
