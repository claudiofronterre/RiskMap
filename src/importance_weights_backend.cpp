#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "likelihood_kernel.h"

#ifndef FCONE
#define FCONE
#endif

using namespace Rcpp;

// Log importance weights of a bank of standard-normal draws for the adaptive
// MCML engine with a canonical link, in one pass and without temporary copies
// on the R heap. It reproduces the R code of the engine's evaluate() and
// log_mixture_proposal() step by step: the same BLAS calls R makes for `%*%`
// and forwardsolve(), and the long double accumulation R uses in rowSums(),
// colSums() and sum(), so the weights are identical to it.

namespace {

typedef long double accumulator;

struct Bank {
  const double* draws;
  int rows;
  int dimension;
};

// `%*%` as R computes it: dgemv for a single column or row, dgemm otherwise.
void multiply_like_r(const double* left, const int left_rows,
                     const int inner, const double* right,
                     const int right_columns, double* result) {
  const double one = 1.0;
  const double zero = 0.0;
  const int increment = 1;
  if (right_columns == 1) {
    F77_CALL(dgemv)("N", &left_rows, &inner, &one, left, &left_rows, right,
                    &increment, &zero, result, &increment FCONE);
  } else if (left_rows == 1) {
    F77_CALL(dgemv)("T", &inner, &right_columns, &one, right, &inner, left,
                    &increment, &zero, result, &increment FCONE);
  } else {
    F77_CALL(dgemm)("N", "N", &left_rows, &right_columns, &inner, &one, left,
                    &left_rows, right, &inner, &zero, result, &left_rows
                    FCONE FCONE);
  }
}

// forwardsolve(t(upper), right_side) for the leading `size` rows and columns
// of an upper factor, overwriting `right_side` (size x columns).
void forward_solve_transposed(const NumericMatrix& upper, const int size,
                              std::vector<double>& right_side,
                              const int columns) {
  if (columns == 0) return;
  std::vector<double> lower(static_cast<R_xlen_t>(size) * size);
  const R_xlen_t stride = upper.nrow();
  for (R_xlen_t column = 0; column < size; ++column) {
    for (R_xlen_t row = 0; row < size; ++row) {
      lower[row + column * size] = upper[column + row * stride];
    }
  }
  const double one = 1.0;
  F77_CALL(dtrsm)("L", "L", "N", "N", &size, &columns, &one, lower.data(),
                  &size, right_side.data(), &size FCONE FCONE FCONE FCONE);
}

// colSums(x^2) of a size x columns matrix.
std::vector<double> column_sums_of_squares(const std::vector<double>& matrix,
                                           const int size,
                                           const int columns) {
  std::vector<double> sums(columns);
  for (R_xlen_t column = 0; column < columns; ++column) {
    accumulator sum = 0.0;
    for (R_xlen_t row = 0; row < size; ++row) {
      const double value = matrix[row + column * size];
      sum += value * value;
    }
    sums[column] = static_cast<double>(sum);
  }
  return sums;
}

// rowSums(x[rows, first:(first + count - 1)]^2) of a column-major matrix
// with `stride` rows, accumulated column by column as rowSums() does.
std::vector<double> row_sums_of_squares(const double* matrix,
                                        const R_xlen_t stride,
                                        const std::vector<int>& rows,
                                        const int first, const int count) {
  std::vector<accumulator> sums(rows.size(), 0.0);
  for (R_xlen_t column = first; column < first + count; ++column) {
    for (std::size_t index = 0; index < rows.size(); ++index) {
      const double value = matrix[rows[index] + column * stride];
      sums[index] += value * value;
    }
  }
  return std::vector<double>(sums.begin(), sums.end());
}

// 2 * sum(log(diag(factor))) for the leading `size` diagonal entries.
double log_determinant(const NumericMatrix& factor, const int size) {
  accumulator sum = 0.0;
  for (R_xlen_t index = 0; index < size; ++index) {
    sum += std::log(factor[index + index * factor.nrow()]);
  }
  return 2 * static_cast<double>(sum);
}

// Transposed submatrix t(x[rows, first:(first + count - 1)]) of the
// column-major matrix x with `stride` rows, plus an optional shift per column
// of x (the proposal centre), as count x rows.
std::vector<double> transpose_rows(const double* matrix, const R_xlen_t stride,
                                   const std::vector<int>& rows,
                                   const int first, const int count,
                                   const double* shift) {
  std::vector<double> result(static_cast<R_xlen_t>(count) * rows.size());
  for (std::size_t index = 0; index < rows.size(); ++index) {
    for (R_xlen_t column = 0; column < count; ++column) {
      const double value = matrix[rows[index] + (first + column) * stride];
      result[column + index * count] = shift == nullptr ? value :
        value + shift[first + column];
    }
  }
  return result;
}

// Deviations from the centre: draws %*% root for every draw, then the
// defensive draws replaced by their deviations under the prior (the spatial
// block through its factor, random-effect levels scaled by their standard
// deviations), as prior_deviations() computes them.
std::vector<double> compute_deviations(const Bank& bank,
                                       const NumericMatrix& root,
                                       const std::vector<int>& defensive_rows,
                                       const NumericMatrix& spatial_root,
                                       const NumericVector& sigma2_re,
                                       const IntegerVector& random_first,
                                       const IntegerVector& random_size) {
  const R_xlen_t rows = bank.rows;
  std::vector<double> deviations(rows * bank.dimension);
  multiply_like_r(bank.draws, bank.rows, bank.dimension, root.begin(),
                  bank.dimension, deviations.data());
  const int defensive = defensive_rows.size();
  if (defensive == 0) return deviations;
  const int locations = spatial_root.nrow();
  std::vector<double> selected(static_cast<R_xlen_t>(defensive) * locations);
  for (R_xlen_t column = 0; column < locations; ++column) {
    for (R_xlen_t index = 0; index < defensive; ++index) {
      selected[index + column * defensive] =
        bank.draws[defensive_rows[index] + column * rows];
    }
  }
  std::vector<double> spatial(selected.size());
  multiply_like_r(selected.data(), defensive, locations, spatial_root.begin(),
                  locations, spatial.data());
  for (R_xlen_t column = 0; column < locations; ++column) {
    for (R_xlen_t index = 0; index < defensive; ++index) {
      deviations[defensive_rows[index] + column * rows] =
        spatial[index + column * defensive];
    }
  }
  for (R_xlen_t effect = 0; effect < random_first.size(); ++effect) {
    const double scale = std::sqrt(sigma2_re[effect]);
    for (R_xlen_t column = random_first[effect] - 1;
         column < random_first[effect] - 1 + random_size[effect]; ++column) {
      for (int index = 0; index < defensive; ++index) {
        const R_xlen_t cell = defensive_rows[index] + column * rows;
        deviations[cell] = bank.draws[cell] * scale;
      }
    }
  }
  return deviations;
}

// Log density of each draw under the target: the prior of the latent effects
// (spatial block standardised by the target factor, random effects with the
// target variances) plus the canonical log-likelihood. Returns the
// standardised spatial samples in `standardised` for reuse by the mixture.
std::vector<double> compute_log_target(
    const std::vector<double>& deviations, const NumericVector& centre,
    const int rows, const NumericMatrix& spatial_root,
    const NumericVector& sigma2_re, const IntegerVector& random_first,
    const IntegerVector& random_size, const NumericVector& linear_offset,
    const IntegerMatrix& columns, const NumericVector& y,
    const NumericVector& units, const bool poisson,
    std::vector<double>& standardised) {
  const int locations = spatial_root.nrow();
  std::vector<int> all_rows(rows);
  for (int row = 0; row < rows; ++row) all_rows[row] = row;
  standardised = transpose_rows(deviations.data(), rows, all_rows, 0,
                                locations, centre.begin());
  forward_solve_transposed(spatial_root, locations, standardised, rows);
  const std::vector<double> prior_norms =
    column_sums_of_squares(standardised, locations, rows);
  const double prior_constant = locations * std::log(2 * M_PI) +
    log_determinant(spatial_root, locations);
  std::vector<double> log_target(rows);
  for (int row = 0; row < rows; ++row) {
    log_target[row] = -0.5 * (prior_constant + prior_norms[row]);
  }
  // Samples are deviations plus the centre, recomputed where needed.
  std::vector<double> samples(deviations.size());
  for (R_xlen_t column = 0; column < centre.size(); ++column) {
    for (R_xlen_t row = 0; row < rows; ++row) {
      samples[row + column * rows] = deviations[row + column * rows] +
        centre[column];
    }
  }
  for (R_xlen_t effect = 0; effect < random_first.size(); ++effect) {
    const std::vector<double> norms = row_sums_of_squares(
      samples.data(), rows, all_rows, random_first[effect] - 1,
      random_size[effect]
    );
    const double constant = random_size[effect] *
      std::log(2 * M_PI * sigma2_re[effect]);
    for (int row = 0; row < rows; ++row) {
      log_target[row] = log_target[row] -
        0.5 * (constant + norms[row] / sigma2_re[effect]);
    }
  }
  std::vector<double> log_likelihood(rows, 0.0);
  std::vector<double> eta(rows);
  for (int observation = 0; observation < columns.nrow(); ++observation) {
    const R_xlen_t location = columns(observation, 0) - 1;
    for (int row = 0; row < rows; ++row) {
      eta[row] = samples[row + location * rows] + linear_offset[observation];
    }
    for (int effect = 1; effect < columns.ncol(); ++effect) {
      const R_xlen_t level = columns(observation, effect) - 1;
      for (int row = 0; row < rows; ++row) {
        eta[row] = eta[row] + samples[row + level * rows];
      }
    }
    for (int row = 0; row < rows; ++row) {
      log_likelihood[row] += canonical_log_likelihood(
        eta[row], y[observation], units[observation], poisson
      );
    }
  }
  for (int row = 0; row < rows; ++row) {
    if (poisson && !std::isfinite(log_likelihood[row])) {
      log_likelihood[row] = -std::numeric_limits<double>::infinity();
    }
    log_target[row] = log_target[row] + log_likelihood[row];
  }
  return log_target;
}

// Log density of each draw under the defensive mixture of the Gaussian
// proposal and the prior, as log_mixture_proposal() computes it.
std::vector<double> compute_log_proposal(
    const Bank& bank, const std::vector<double>& deviations,
    const std::vector<int>& defensive_rows,
    const std::vector<int>& proposal_rows, const double share,
    const NumericMatrix& root, const NumericVector& centre,
    const NumericMatrix& spatial_root, const NumericVector& sigma2_re,
    const IntegerVector& random_first, const IntegerVector& random_size,
    const std::vector<double>* standardised) {
  const int rows = bank.rows;
  std::vector<int> all_rows(rows);
  for (int row = 0; row < rows; ++row) all_rows[row] = row;
  const std::vector<double> draw_norms = row_sums_of_squares(
    bank.draws, rows, all_rows, 0, bank.dimension
  );
  const double log_normalising = bank.dimension * std::log(2 * M_PI);
  const double log_proposal_determinant = log_determinant(root,
                                                          bank.dimension);
  std::vector<double> result(rows);
  if (share == 0) {
    for (int row = 0; row < rows; ++row) {
      result[row] = -0.5 * (log_normalising + log_proposal_determinant +
        draw_norms[row]);
    }
    return result;
  }
  std::vector<double> proposal_norms = draw_norms;
  std::vector<double> defensive_solved = transpose_rows(
    deviations.data(), rows, defensive_rows, 0, bank.dimension, nullptr
  );
  forward_solve_transposed(root, bank.dimension, defensive_solved,
                           defensive_rows.size());
  const std::vector<double> defensive_norms = column_sums_of_squares(
    defensive_solved, bank.dimension, defensive_rows.size()
  );
  for (std::size_t index = 0; index < defensive_rows.size(); ++index) {
    proposal_norms[defensive_rows[index]] = defensive_norms[index];
  }
  const int locations = spatial_root.nrow();
  std::vector<double> spatial_norms;
  if (standardised != nullptr) {
    std::vector<double> centre_solved(centre.begin(),
                                      centre.begin() + locations);
    forward_solve_transposed(spatial_root, locations, centre_solved, 1);
    spatial_norms.resize(proposal_rows.size());
    for (std::size_t index = 0; index < proposal_rows.size(); ++index) {
      accumulator sum = 0.0;
      const R_xlen_t offset =
        static_cast<R_xlen_t>(proposal_rows[index]) * locations;
      for (R_xlen_t location = 0; location < locations; ++location) {
        const double value = (*standardised)[location + offset] -
          centre_solved[location];
        sum += value * value;
      }
      spatial_norms[index] = static_cast<double>(sum);
    }
  } else {
    std::vector<double> solved = transpose_rows(
      deviations.data(), rows, proposal_rows, 0, locations, nullptr
    );
    forward_solve_transposed(spatial_root, locations, solved,
                             proposal_rows.size());
    spatial_norms = column_sums_of_squares(solved, locations,
                                           proposal_rows.size());
  }
  accumulator random_determinant = 0.0;
  for (R_xlen_t effect = 0; effect < random_first.size(); ++effect) {
    const std::vector<double> norms = row_sums_of_squares(
      deviations.data(), rows, proposal_rows, random_first[effect] - 1,
      random_size[effect]
    );
    for (std::size_t index = 0; index < proposal_rows.size(); ++index) {
      spatial_norms[index] = spatial_norms[index] +
        norms[index] / sigma2_re[effect];
    }
    random_determinant += random_size[effect] * std::log(sigma2_re[effect]);
  }
  std::vector<double> prior_norms = draw_norms;
  for (std::size_t index = 0; index < proposal_rows.size(); ++index) {
    prior_norms[proposal_rows[index]] = spatial_norms[index];
  }
  const double log_prior_determinant = log_determinant(spatial_root,
                                                       locations) +
    (random_first.size() ? static_cast<double>(random_determinant) : 0.0);
  for (int row = 0; row < rows; ++row) {
    const double log_gaussian = std::log1p(-share) - 0.5 * (
      log_normalising + log_proposal_determinant + proposal_norms[row]);
    const double log_defensive = std::log(share) - 0.5 * (
      log_normalising + log_prior_determinant + prior_norms[row]);
    if (std::isnan(log_gaussian) || std::isnan(log_defensive)) {
      result[row] = NA_REAL;
      continue;
    }
    const double largest = std::max(log_gaussian, log_defensive);
    result[row] = largest + std::log(std::exp(log_gaussian - largest) +
      std::exp(log_defensive - largest));
  }
  return result;
}

} // namespace

// Log importance weights log p(y, u) - log q(u) of each row of `draws` for a
// canonical binomial or Poisson model. `root`, `centre`, `proposal_spatial_root`
// and `proposal_sigma2_re` define the proposal and its defensive component,
// with `defensive` marking the rows drawn from the prior and `share` their
// share of the whole bank. The target is given by `target_spatial_root`,
// `target_sigma2_re` and `linear_offset` (D beta plus the offset).
// `same_parameters` says whether the target and proposal parameters coincide,
// so that the standardised spatial samples serve both densities. `columns`
// holds each observation's latent columns (location, then random-effect
// levels) and `random_first` and `random_size` the latent columns of each
// random effect.
// [[Rcpp::export]]
NumericVector cpp_importance_log_weights(
    const NumericMatrix& draws, const NumericMatrix& root,
    const NumericVector& centre, const LogicalVector& defensive,
    const double share, const NumericMatrix& proposal_spatial_root,
    const NumericVector& proposal_sigma2_re,
    const NumericMatrix& target_spatial_root,
    const NumericVector& target_sigma2_re, const bool same_parameters,
    const NumericVector& linear_offset, const IntegerMatrix& columns,
    const IntegerVector& random_first, const IntegerVector& random_size,
    const NumericVector& y, const NumericVector& units, const bool poisson) {
  const Bank bank = {draws.begin(), draws.nrow(), draws.ncol()};
  std::vector<int> defensive_rows;
  std::vector<int> proposal_rows;
  for (int row = 0; row < bank.rows; ++row) {
    (defensive[row] ? defensive_rows : proposal_rows).push_back(row);
  }
  const std::vector<double> deviations = compute_deviations(
    bank, root, defensive_rows, proposal_spatial_root, proposal_sigma2_re,
    random_first, random_size
  );
  std::vector<double> standardised;
  const std::vector<double> log_target = compute_log_target(
    deviations, centre, bank.rows, target_spatial_root, target_sigma2_re,
    random_first, random_size, linear_offset, columns, y, units, poisson,
    standardised
  );
  const std::vector<double> log_proposal = compute_log_proposal(
    bank, deviations, defensive_rows, proposal_rows, share, root, centre,
    proposal_spatial_root, proposal_sigma2_re, random_first, random_size,
    same_parameters ? &standardised : nullptr
  );
  NumericVector log_weights(bank.rows);
  for (int row = 0; row < bank.rows; ++row) {
    log_weights[row] = log_target[row] - log_proposal[row];
  }
  return log_weights;
}
