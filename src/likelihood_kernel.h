#ifndef RISKMAP_LIKELIHOOD_KERNEL_H
#define RISKMAP_LIKELIHOOD_KERNEL_H

#include <cfloat>
#include <cmath>

// The distribution functions of R's mathematical library, declared here so
// that this header does not depend on the order in which Rcpp, RcppEigen and
// Rmath.h are included.
extern "C" {
double Rf_pnorm5(double, double, double, int, int);
double Rf_pcauchy(double, double, double, int, int);
}

// The links whose log-likelihood is compiled, numbered as the family object's
// `code` (make_family() in R/family.R). Every other link is computed in R.
enum LinkCode : int {
  link_logit = 0,    // canonical binomial
  link_log = 1,      // canonical Poisson
  link_probit = 2,
  link_cloglog = 3,
  link_cauchit = 4
};

inline bool is_canonical_link(const int code) {
  return code == link_logit || code == link_log;
}

// log(1 + exp(x)) without overflow, the binomial normalising term on the
// logit scale.
inline double softplus(const double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) :
    std::log1p(std::exp(x));
}

// One observation's canonical log-likelihood at linear predictor `eta`, up to
// the constant in y: y eta - m exp(eta) (Poisson) or y eta - m log(1 +
// exp(eta)) (binomial). The compiled weights, the compiled rows of the
// log-likelihood, the site tables and the Langevin sampler all take it from
// here.
inline double canonical_log_likelihood(const double eta, const double y,
                                       const double units,
                                       const bool poisson) {
  return poisson ? y * eta - units * std::exp(eta) :
    y * eta - units * softplus(eta);
}

// The inverse of a named binomial link, as named_inverse_link() writes it in
// R, with the same functions of R's library.
inline double named_inverse_link(const double eta, const int code) {
  switch (code) {
  case link_probit:
    return Rf_pnorm5(eta, 0.0, 1.0, 1, 0);
  case link_cloglog:
    return -std::expm1(-std::exp(eta));
  default:
    return Rf_pcauchy(eta, 0.0, 1.0, 1, 0);
  }
}

// One observation's log-likelihood under a named binomial link, operation for
// operation as compute_link_log_likelihood_rows() and bound_log_probability()
// compute it in R: y log max(p, DBL_MIN) + (m - y) log1p(-min(p, 1 - DBL_EPSILON
// / 2)). `outside` is set where p is not a probability, which gives the whole
// draw a log-likelihood of -Inf. The value is then not used.
inline double named_link_log_likelihood(const double eta, const double y,
                                        const double units, const int code,
                                        bool& outside) {
  const double p = named_inverse_link(eta, code);
  if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
    outside = true;
    return 0.0;
  }
  const double success = std::log(p < DBL_MIN ? DBL_MIN : p);
  const double ceiling = 1.0 - DBL_EPSILON / 2.0;
  const double failure = std::log1p(-(p > ceiling ? ceiling : p));
  return success * y + failure * (units - y);
}

#endif
