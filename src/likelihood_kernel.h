#ifndef RISKMAP_LIKELIHOOD_KERNEL_H
#define RISKMAP_LIKELIHOOD_KERNEL_H

#include <cmath>

// log(1 + exp(x)) without overflow, the binomial normalising term on the
// logit scale.
inline double softplus(const double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) :
    std::log1p(std::exp(x));
}

// One observation's canonical log-likelihood at linear predictor `eta`, up to
// the constant in y: y eta - m exp(eta) (Poisson) or y eta - m log(1 +
// exp(eta)) (binomial).
inline double canonical_log_likelihood(const double eta, const double y,
                                       const double units,
                                       const bool poisson) {
  return poisson ? y * eta - units * std::exp(eta) :
    y * eta - units * softplus(eta);
}

#endif
