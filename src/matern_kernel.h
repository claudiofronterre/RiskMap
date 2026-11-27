#ifndef RISKMAP_MATERN_KERNEL_H
#define RISKMAP_MATERN_KERNEL_H

#include <cmath>

// Matern correlation at one distance for kappa = 0.5, 1.5 or 2.5, in closed
// form. The callers check that phi is positive, the distances are finite and
// non-negative, and kappa is one of the three values.
inline double half_integer_matern(const double distance, const double phi,
                                  const double kappa) {
  if (distance == 0.0) return 1.0;
  const double scaled_distance = distance / phi;
  if (scaled_distance > 600.0) return 0.0;
  if (kappa == 0.5) return std::exp(-scaled_distance);
  if (kappa == 1.5) {
    return (1.0 + scaled_distance) * std::exp(-scaled_distance);
  }
  return (1.0 + scaled_distance +
    scaled_distance * scaled_distance / 3.0) * std::exp(-scaled_distance);
}

#endif
