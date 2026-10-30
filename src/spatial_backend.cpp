#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "matern_kernel.h"

using namespace Rcpp;

namespace {

double euclidean_distance(const NumericMatrix& first,
                          const int first_row,
                          const NumericMatrix& second,
                          const int second_row) {
  double squared_distance = 0.0;

  for (int column = 0; column < first.ncol(); ++column) {
    const double difference = first(first_row, column) -
      second(second_row, column);
    squared_distance += difference * difference;
  }

  return std::sqrt(squared_distance);
}

} // namespace

// [[Rcpp::export]]
NumericVector cpp_pairwise_distances(const NumericMatrix coordinates) {

  const int number_locations = coordinates.nrow();
  const R_xlen_t number_distances =
    static_cast<R_xlen_t>(number_locations) * (number_locations - 1) / 2;
  NumericVector distances(number_distances);
  R_xlen_t index = 0;

  // Match the column-wise lower-triangle ordering used by stats::dist(). The
  // two-dimensional path avoids a function call and inner loop for every
  // distance, which matters for large prediction grids.
  if (coordinates.ncol() == 2) {
    for (int column = 0; column < number_locations - 1; ++column) {
      for (int row = column + 1; row < number_locations; ++row) {
        const double x_difference = coordinates(row, 0) -
          coordinates(column, 0);
        const double y_difference = coordinates(row, 1) -
          coordinates(column, 1);
        distances[index++] = std::sqrt(x_difference * x_difference +
          y_difference * y_difference);
      }
    }
  } else {
    for (int column = 0; column < number_locations - 1; ++column) {
      for (int row = column + 1; row < number_locations; ++row) {
        distances[index++] = euclidean_distance(coordinates, row,
                                                coordinates, column);
      }
    }
  }

  return distances;
}

// [[Rcpp::export]]
NumericMatrix cpp_cross_distances(const NumericMatrix first,
                                  const NumericMatrix second) {

  NumericMatrix distances(first.nrow(), second.nrow());

  if (first.ncol() == 2) {
    for (int first_row = 0; first_row < first.nrow(); ++first_row) {
      for (int second_row = 0; second_row < second.nrow(); ++second_row) {
        const double x_difference = first(first_row, 0) -
          second(second_row, 0);
        const double y_difference = first(first_row, 1) -
          second(second_row, 1);
        distances(first_row, second_row) = std::sqrt(
          x_difference * x_difference + y_difference * y_difference
        );
      }
    }
  } else {
    for (int first_row = 0; first_row < first.nrow(); ++first_row) {
      for (int second_row = 0; second_row < second.nrow(); ++second_row) {
        distances(first_row, second_row) = euclidean_distance(
          first, first_row, second, second_row
        );
      }
    }
  }

  return distances;
}

// [[Rcpp::export]]
NumericMatrix cpp_binned_semivariances(const NumericMatrix permuted_values,
                                       const IntegerVector first_index,
                                       const IntegerVector second_index,
                                       const IntegerVector bin_index,
                                       const int number_bins) {
  // The only caller, variogram(), validates the values and builds zero-based
  // pair and bin indices of equal length, so they are not rechecked here
  const int number_values = permuted_values.nrow();
  const int number_pairs = first_index.length();

  IntegerVector bin_counts(number_bins);
  for (int pair = 0; pair < number_pairs; ++pair) {
    ++bin_counts[bin_index[pair]];
  }

  // Permutations are processed in blocks. The values of a block are copied so
  // that those of each location are contiguous and the copy fits in cache,
  // which lets each pair update the whole block at once. Pairs are still
  // summed in the same order, and halving the sums at the end is exact, so the
  // result does not depend on the block size. A fixed loop length for full
  // blocks lets the compiler vectorise the update
  const int number_curves = permuted_values.ncol();
  constexpr int block_size = 32;
  std::vector<double> block_values(static_cast<size_t>(number_values) * block_size);
  std::vector<double> block_sums(static_cast<size_t>(number_bins) * block_size);
  NumericMatrix semivariances(number_bins, number_curves);
  const double* values_pointer = permuted_values.begin();
  const int* first_pointer = first_index.begin();
  const int* second_pointer = second_index.begin();
  const int* bin_pointer = bin_index.begin();
  for (int block_start = 0; block_start < number_curves; block_start += block_size) {
    const int width = std::min(block_size, number_curves - block_start);
    for (int curve = 0; curve < width; ++curve) {
      const double* values_column = values_pointer +
        static_cast<R_xlen_t>(block_start + curve) * number_values;
      for (int location = 0; location < number_values; ++location) {
        block_values[static_cast<size_t>(location) * block_size + curve] =
          values_column[location];
      }
    }
    std::fill(block_sums.begin(), block_sums.end(), 0.0);
    for (int pair = 0; pair < number_pairs; ++pair) {
      const double* first_values = block_values.data() +
        static_cast<size_t>(first_pointer[pair]) * block_size;
      const double* second_values = block_values.data() +
        static_cast<size_t>(second_pointer[pair]) * block_size;
      double* sums = block_sums.data() +
        static_cast<size_t>(bin_pointer[pair]) * block_size;
      if (width == block_size) {
        for (int curve = 0; curve < block_size; ++curve) {
          const double difference = first_values[curve] - second_values[curve];
          sums[curve] += difference * difference;
        }
      } else {
        for (int curve = 0; curve < width; ++curve) {
          const double difference = first_values[curve] - second_values[curve];
          sums[curve] += difference * difference;
        }
      }
    }
    for (int bin = 0; bin < number_bins; ++bin) {
      for (int curve = 0; curve < width; ++curve) {
        semivariances(bin, block_start + curve) = bin_counts[bin] == 0 ? NA_REAL :
          0.5 * block_sums[static_cast<size_t>(bin) * block_size + curve] /
          bin_counts[bin];
      }
    }
  }

  return semivariances;
}

// [[Rcpp::export]]
NumericVector cpp_half_integer_matern(const NumericVector distances,
                                      const double phi,
                                      const double kappa) {
  // matern_correlation() checks that phi is positive and the distances are
  // finite and non-negative before calling this kernel
  NumericVector correlation(distances.length());

  for (R_xlen_t i = 0; i < distances.length(); ++i) {
    correlation[i] = half_integer_matern(distances[i], phi, kappa);
  }

  return correlation;
}
