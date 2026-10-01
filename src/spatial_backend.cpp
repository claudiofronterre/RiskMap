#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

namespace {

void check_coordinates(const NumericMatrix& coordinates) {
  if (coordinates.ncol() < 1) {
    stop("`coordinates` must contain at least one column.");
  }

  for (R_xlen_t i = 0; i < coordinates.length(); ++i) {
    if (!R_finite(coordinates[i])) {
      stop("`coordinates` must contain only finite values.");
    }
  }
}

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
  check_coordinates(coordinates);

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
  check_coordinates(first);
  check_coordinates(second);

  if (first.ncol() != second.ncol()) {
    stop("Coordinate matrices must have the same number of columns.");
  }

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
NumericVector cpp_half_integer_matern(const NumericVector distances,
                                      const double phi,
                                      const double kappa) {
  if (!R_finite(phi) || phi <= 0.0) {
    stop("`phi` must be a finite positive number.");
  }

  NumericVector correlation(distances.length());

  for (R_xlen_t i = 0; i < distances.length(); ++i) {
    const double distance = distances[i];

    if (!R_finite(distance) || distance < 0.0) {
      stop("Distances must be finite and non-negative.");
    }

    if (distance == 0.0) {
      correlation[i] = 1.0;
      continue;
    }

    const double scaled_distance = distance / phi;

    if (scaled_distance > 600.0) {
      correlation[i] = 0.0;
    } else if (kappa == 0.5) {
      correlation[i] = std::exp(-scaled_distance);
    } else if (kappa == 1.5) {
      correlation[i] = (1.0 + scaled_distance) *
        std::exp(-scaled_distance);
    } else if (kappa == 2.5) {
      correlation[i] = (1.0 + scaled_distance +
        scaled_distance * scaled_distance / 3.0) *
        std::exp(-scaled_distance);
    } else {
      stop("The compiled Matérn kernel only supports kappa = 0.5, 1.5 or 2.5.");
    }
  }

  return correlation;
}
