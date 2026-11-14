#ifndef RISKMAP_SITE_TABLE_H
#define RISKMAP_SITE_TABLE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <list>
#include <memory>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include "likelihood_kernel.h"

// The log-likelihood of a record with its independent effects integrated out,
//   h(x; v) = log integral of exp(k(x + e)) N(e; 0, v) de,
// for the binomial and Poisson kernels k of any link, with y events of m
// trials or exposure m. A table holds h and its first two derivatives in x
// on the grid -50, -49.75, ..., 30 for each distinct pair (y, m) and one value
// of v, and the quintic Hermite interpolant through them gives h with two
// continuous derivatives between grid points. A cell is filled the first time
// it is needed, and it is a pure function of (y, m, x, v), so the result does
// not depend on the order in which cells are filled. Beyond the grid the
// integral is evaluated directly. Up to 32 tables are kept, keyed on the exact
// value of v, and the least recently used is dropped.
//
// For the canonical links the rule is adaptive Gauss-Hermite quadrature with
// 20 nodes up to v = 0.2, and above it the trapezium rule on the lattice of
// multiples of a step of at most 0.25 (for a nested pattern with an inner
// variance vi above 16, at most sqrt(vi) / 16) and at most half the integrand's
// width at its mode (Trefethen and Weideman, 2014), walked out from the mode
// until a term is negligible, which is safe because the integrand is
// log-concave. Where the kernel is flat on one side and v is large, the rest
// of that side is taken in closed form (uniform_cell()). A custom link takes
// the trapezium rule of custom_cell(), which assumes no mode. The
// derivatives are those of Stein's identity, h' = E k' and
// h'' = E k'' + Var k' under the tilted distribution, which have no
// cancellation as v goes to zero.
//
// A nested pattern is a location whose units (records after identical ones
// are summed) share their linear predictor x and each carry an effect of
// variance vi of their own, such as households, under a shared effect of
// variance vo, such as the nugget. Its log-likelihood is
//   H(x; vo, vi) = log integral of N(e; 0, vo) prod_u g(y_u, m_u, x + e; vi) de,
// computed by the same rules applied to the outer integral, with the inner
// h_u from their record table at vi, and tabled in the same way, one table per
// exact pair (vo, vi), up to 16 of them.
//
// A group is a set of units at one location whose linear predictors differ,
// such as the children of a household who differ in sex or birth order, under
// an effect of variance v that they share. Its members are units j with pair
// (y_j, m_j) and linear predictor x + r_j, where x is the linear predictor of
// the group's record and r_j the member's relative linear predictor, which
// depends on the regression coefficients. Its log-likelihood
//   h_g(x; v) = log integral of N(u; 0, v) prod_j g(y_j, m_j, x + r_j + u) du
// is that of a record whose kernel is the sum of its members' kernels, so it
// takes the rules above, except that under a custom link known to be
// log-concave, the named probit and cloglog links, it takes the trapezium
// rule walked out from its mode (uniform_cell()) at every variance, since
// custom_cell(), which assumes no mode, is built for tables filled once per
// variance and these are filled again as the coefficients change. A location
// of groups nests them under an outer
// effect of variance vo, such as the nugget,
//   H_l(x; vo, v) = log integral of N(z; 0, vo) prod_g exp(h_g(x + z; v)) dz,
// as a nested pattern nests units that share x. Since r depends on the
// coefficients, a group's integral is that of its pattern, the variance and
// its members' pairs and linear predictors relative to a reference member,
// computed from the differences of their covariates, so that h_g(x) is the
// pattern's integral at x plus that member's r. A change of the coefficients
// of covariates that its members share, such as those of the household or the
// location, leaves a group's pattern unchanged, and every group with the same
// pattern, such as households of two children without deaths who differ in
// sex, shares its integral. Each location of
// groups keeps up to 3 integrals in the same way, keyed on its variances and
// the relative linear predictors of all its members. Each kept integral holds
// only the grid cells it has needed, in a strip that grows as the draws reach
// further. A cell is a pure function of the key, the variances and its grid
// point, as in the tables above. Within a location, a group of one member,
// such as a household of one child, takes its record's h at x + r_j from the
// record table at v, which does not depend on the coefficients and is shared
// by every record with its pair.

namespace riskmap_site {

inline double logistic(const double t) { return 1.0 / (1.0 + std::exp(-t)); }

// The quintic Hermite interpolant on [0, 1] through values f, first
// derivatives d and second derivatives s at both ends (in units of the
// interval), with its first two derivatives in t.
inline void quintic(const double t, const double f0, const double f1,
                    const double d0, const double d1, const double s0,
                    const double s1, double& value, double& first,
                    double& second) {
  const double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
  const double tail = 10 * t3 - 15 * t4 + 6 * t5;
  value = (1 - tail) * f0 + tail * f1 +
    (t - 6 * t3 + 8 * t4 - 3 * t5) * d0 + (-4 * t3 + 7 * t4 - 3 * t5) * d1 +
    0.5 * (t2 - 3 * t3 + 3 * t4 - t5) * s0 + 0.5 * (t3 - 2 * t4 + t5) * s1;
  first = (30 * t2 - 60 * t3 + 30 * t4) * (f1 - f0) +
    (1 - 18 * t2 + 32 * t3 - 15 * t4) * d0 +
    (-12 * t2 + 28 * t3 - 15 * t4) * d1 +
    0.5 * (2 * t - 9 * t2 + 12 * t3 - 5 * t4) * s0 +
    0.5 * (3 * t2 - 8 * t3 + 5 * t4) * s1;
  second = (60 * t - 180 * t2 + 120 * t3) * (f1 - f0) +
    (-36 * t + 96 * t2 - 60 * t3) * d0 + (-24 * t + 84 * t2 - 60 * t3) * d1 +
    0.5 * (2 - 18 * t + 36 * t2 - 20 * t3) * s0 +
    0.5 * (6 * t - 24 * t2 + 20 * t3) * s1;
}

// The kernel of a custom link written as k(eta) = y A(eta) + m B(eta) + c,
// with A, B and their first two derivatives on a lattice from `low` with
// `step`, computed once from the link's R functions, and c = y log(m) for the
// Poisson family. Between lattice points they are interpolated by quintic
// Hermite polynomials, and beyond the lattice they are constant.
struct CustomArrays {
  double low = 0, step = 1;
  std::vector<double> A, A1, A2, B, B1, B2;
  bool poisson = false;
};

// The log-likelihood kernel of one record in its linear predictor.
struct Kernel {
  bool poisson = false;
  const CustomArrays* custom = nullptr;

  // k, k' and k'' at eta.
  void terms(const double y, const double m, const double eta, double& k,
             double& k1, double& k2) const {
    if (custom != nullptr) {
      custom_terms(y, m, eta, k, k1, k2);
      return;
    }
    // k is the shared canonical kernel (likelihood_kernel.h). The derivatives
    // keep this table's logistic, which differs from the Langevin sampler's
    // in the last bits for negative eta.
    k = canonical_log_likelihood(eta, y, m, poisson);
    if (poisson) {
      const double mean = m * std::exp(eta);
      k1 = y - mean;
      k2 = -mean;
    } else {
      const double p = logistic(eta);
      k1 = y - m * p;
      k2 = -m * p * (1.0 - p);
    }
  }

 private:
  void custom_terms(const double y, const double m, const double eta,
                    double& k, double& k1, double& k2) const {
    const CustomArrays& c = *custom;
    const long last = static_cast<long>(c.A.size()) - 1;
    const double position = (eta - c.low) / c.step;
    double a = 0, a1 = 0, a2 = 0, b = 0, b1 = 0, b2 = 0;
    if (!(position > 0.0)) {
      a = c.A[0]; b = c.B[0];
    } else if (!(position < last)) {
      a = c.A[last]; b = c.B[last];
    } else {
      const long at = static_cast<long>(position);
      const double t = position - at;
      const double h = c.step, h2 = c.step * c.step;
      quintic(t, c.A[at], c.A[at + 1], h * c.A1[at], h * c.A1[at + 1],
              h2 * c.A2[at], h2 * c.A2[at + 1], a, a1, a2);
      quintic(t, c.B[at], c.B[at + 1], h * c.B1[at], h * c.B1[at + 1],
              h2 * c.B2[at], h2 * c.B2[at + 1], b, b1, b2);
      a1 /= h; b1 /= h; a2 /= h2; b2 /= h2;
    }
    k = y * a + m * b + (c.poisson && m > 0 ? y * std::log(m) : 0.0);
    k1 = y * a1 + m * b1;
    k2 = y * a2 + m * b2;
  }
};

// h and its first two derivatives, and for a nested pattern's slope the
// tilted mean of an extra quantity (cell_rules below).
struct Moments {
  double h, h1, h2;
  double extra = 0.0;
};

struct Rule {
  std::vector<double> z;      // nodes of Gauss-Hermite for exp(-z^2)
  std::vector<double> log_w;  // their log weights
};

// The terms of one record's kernel at eta, for the rules below.
struct RecordTerms {
  static constexpr bool has_extra = false;
  const Kernel& kernel;
  const double y, m;
  void operator()(const double eta, double& k, double& k1, double& k2) const {
    kernel.terms(y, m, eta, k, k1, k2);
  }
  double extra() const { return 0.0; }
};

// The mode of k(eta) - (eta - x)^2 / (2 v), by Newton's method from x with
// steps of at most 3, and the curvature of the integrand there. `Terms` gives
// k and its first two derivatives at eta, for a record or a nested pattern.
// Newton's method can fail to arrive: a mode beyond the 600 that 200 steps of
// 3 reach, where k is nearly linear and v large, a Poisson kernel far above
// its mean, where it moves by about 1 per step, or one that overflows, above a
// linear predictor of about 709. If a step is not finite or the last one is
// above 1e-6, the mode is found instead by bisection on the slope of the log
// integrand, which decreases since the integrand is log-concave, in a bracket
// found by doubling steps from x. Wherever Newton's method arrives, as in
// every ordinary fit, its mode is kept, since the rule of aghq_cell() places
// its nodes at the mode itself.
template <class Terms>
inline double integrand_mode(const Terms& terms, const double x,
                             const double v, double& curvature) {
  double t = x, k = 0, k1 = 0, k2 = 0, last = 0;
  bool finite = true;
  for (int iteration = 0; iteration < 200; ++iteration) {
    terms(t, k, k1, k2);
    double step = (k1 - (t - x) / v) / (-k2 + 1.0 / v);
    if (!std::isfinite(step)) {
      finite = false;
      break;
    }
    step = std::max(std::min(step, 3.0), -3.0);
    t += step;
    last = step;
    if (std::fabs(step) < 1e-11) break;
  }
  if (!finite || std::fabs(last) > 1e-6) {
    // The slope of the log integrand, minus infinity where k overflows.
    const auto slope_at = [&](const double at) {
      terms(at, k, k1, k2);
      const double slope = k1 - (at - x) / v;
      return std::isnan(slope) ? -std::numeric_limits<double>::infinity() :
        slope;
    };
    double low = x, high = x, width = 1.0;
    if (slope_at(x) > 0) {
      for (int doubling = 0; doubling < 2100; ++doubling, width *= 2.0) {
        high = x + width;
        if (!(slope_at(high) > 0)) break;
        low = high;
      }
    } else {
      for (int doubling = 0; doubling < 2100; ++doubling, width *= 2.0) {
        low = x - width;
        if (slope_at(low) > 0) break;
        high = low;
      }
    }
    for (int iteration = 0; iteration < 200 &&
           high - low > 1e-12 * std::max(1.0, std::fabs(low)); ++iteration) {
      const double middle = 0.5 * (low + high);
      if (slope_at(middle) > 0) low = middle; else high = middle;
    }
    t = 0.5 * (low + high);
  }
  terms(t, k, k1, k2);
  curvature = -k2 + 1.0 / v;
  return t;
}

// With Terms::has_extra, each rule also gives the tilted mean of
// terms.extra(), read after each call of terms() at a node.
template <class Terms>
inline Moments aghq_cell(const Terms& terms, const double x, const double v,
                         const Rule& rule) {
  double curvature = 0;
  const double t = integrand_mode(terms, x, v, curvature);
  const double a = 1.0 / std::sqrt(curvature);
  const double constant = -0.5 * std::log(2.0 * M_PI * v);
  const int nodes = rule.z.size();
  std::vector<double> log_f(nodes), score(nodes), second(nodes),
    extra(Terms::has_extra ? nodes : 0);
  double top = -std::numeric_limits<double>::infinity();
  for (int j = 0; j < nodes; ++j) {
    const double eta = t + M_SQRT2 * a * rule.z[j];
    const double e = eta - x;
    double k = 0, k1 = 0, k2 = 0;
    terms(eta, k, k1, k2);
    if constexpr (Terms::has_extra) extra[j] = terms.extra();
    log_f[j] = k - e * e / (2.0 * v) + constant + rule.log_w[j] +
      rule.z[j] * rule.z[j];
    score[j] = k1;
    second[j] = k1 * k1 + k2;
    top = std::max(top, log_f[j]);
  }
  double total = 0, first = 0, square = 0, sum_extra = 0;
  for (int j = 0; j < nodes; ++j) {
    const double w = std::exp(log_f[j] - top);
    total += w;
    first += w * score[j];
    square += w * second[j];
    if constexpr (Terms::has_extra) sum_extra += w * extra[j];
  }
  const double mean = first / total;
  return {top + std::log(M_SQRT2 * a * total), mean,
          square / total - mean * mean, sum_extra / total};
}

// log(Phi(z) exp(z^2 / 2)) for z <= 0, with Phi the standard normal
// distribution function, by the complementary error function down to -30 and
// by the asymptotic series of Mills' ratio below it, whose first omitted term
// is below 2e-12 there.
inline double log_scaled_normal_tail(const double z) {
  if (z > -30.0) return 0.5 * z * z + std::log(0.5 * std::erfc(-z / M_SQRT2));
  const double r = 1.0 / (z * z);
  return -std::log(-z) - 0.5 * std::log(2.0 * M_PI) +
    std::log1p(r * (-1.0 + r * (3.0 + r * (-15.0 + r * 105.0))));
}

// The trapezium rule walks out from the mode until a term is negligible, and
// with a large v the integrand stays near its top for about sqrt(80 v) on the
// side where the kernel flattens (no events, or as many events as trials, for
// a binomial record), so the walk would grow with the square root of v. The
// canonical kernels bend only near a linear predictor of zero (the logit) or
// above it (the log), and so do sums of their integrals. After
// `rest_after_nodes` nodes on one side, where the kernel is linear to within
// its curvature `linear_curvature`, a stretch over which it stays linear is
// taken in closed form, the integral of the Gaussian factor times the
// exponential of the kernel's tangent. Beyond zero that stretch is the rest of
// the side. Before zero it runs to the last node short of the bend, found by
// bisection on the curvature and checked against the tangent, and the walk
// resumes there. Each node where such a stretch starts or ends takes half its
// weight, with the Euler-Maclaurin term of order delta^2 for that end. A walk
// shorter than `rest_after_nodes` on both sides, which is every walk for v up
// to about 125, is unchanged.
constexpr long rest_after_nodes = 400;
constexpr double linear_curvature = 1e-12;

template <class Terms>
inline Moments uniform_cell(const Terms& terms, const double x,
                            const double v, const double step_cap = 0.25) {
  double curvature = 0;
  const double t = integrand_mode(terms, x, v, curvature);
  const double width = 1.0 / std::sqrt(curvature);
  double delta = step_cap;
  while (delta > 0.5 * width && delta > 1e-6) delta *= 0.5;
  const double constant = -0.5 * std::log(2.0 * M_PI * v);
  const double sd = std::sqrt(v);
  const double cut = 40.0;
  const long centre = std::lround(t / delta);
  std::vector<double> log_f, score, second, extra;
  double top = -std::numeric_limits<double>::infinity();
  // A stretch taken in closed form: the node where it starts and, unless it
  // is the rest of its side, the node where it ends, its log mass, the
  // kernel's slope and the extra quantity over it, and the slopes of the log
  // integrand at its ends.
  struct Stretch {
    std::size_t start = 0, end = 0;
    bool ends_at_node = false;
    double log_mass = 0, k1 = 0, extra = 0, start_slope = 0, end_slope = 0;
  };
  std::vector<Stretch> stretches;
  for (int direction = 0; direction < 2; ++direction) {
    const long sign = direction == 0 ? 1 : -1;
    long walked = 0;
    bool tried = false;
    for (long j = direction == 0 ? centre : centre - 1;; j += sign) {
      const double eta = j * delta;
      const double e = eta - x;
      double k = 0, k1 = 0, k2 = 0;
      terms(eta, k, k1, k2);
      const double value = k - e * e / (2.0 * v) + constant;
      if (!std::isfinite(value)) break;
      top = std::max(top, value);
      log_f.push_back(value);
      score.push_back(k1);
      second.push_back(k1 * k1 + k2);
      if constexpr (Terms::has_extra) extra.push_back(terms.extra());
      if (value < top - cut || log_f.size() > 4000000) break;
      // The slope of the log integrand, negative going up from the mode and
      // positive going down.
      const double slope = k1 - e / v;
      if (++walked < rest_after_nodes ||
          !(std::fabs(k2) <= linear_curvature) || !(sign * slope < 0)) {
        continue;
      }
      Stretch stretch;
      stretch.start = log_f.size() - 1;
      stretch.k1 = k1;
      if constexpr (Terms::has_extra) stretch.extra = extra.back();
      stretch.start_slope = slope;
      // With k linear beyond eta, the Gaussian factor times exp(k) beyond it
      // in the walk's direction is Phi(z) times
      // exp(value + log(2 pi v) / 2 + z^2 / 2), z = -|slope| sqrt(v).
      const double z = -std::fabs(slope) * sd;
      if (sign * eta > 0) {
        stretch.log_mass = value - constant + log_scaled_normal_tail(z);
        stretches.push_back(stretch);
        break;
      }
      if (tried) continue;
      tried = true;
      // Before zero, the last node short of the bend, by bisection between
      // this node and the one at zero, worth a stretch if 64 or more ahead.
      long low = j, high = 0;
      double end_k = 0, end_k1 = 0, end_k2 = 0;
      terms(0.0, end_k, end_k1, end_k2);
      if (std::fabs(end_k2) <= linear_curvature) {
        low = high;
      } else {
        while (sign * (high - low) > 1) {
          const long middle = low + (high - low) / 2;
          terms(middle * delta, end_k, end_k1, end_k2);
          if (std::fabs(end_k2) <= linear_curvature) {
            low = middle;
          } else {
            high = middle;
          }
        }
      }
      if (sign * (low - j) < 64) continue;
      // The tangent at this node must still hold at the end of the stretch.
      const double end_eta = low * delta;
      terms(end_eta, end_k, end_k1, end_k2);
      if (!(std::fabs(end_k - (k + k1 * (end_eta - eta))) <=
            1e-9 * std::max(1.0, std::fabs(end_k)))) {
        continue;
      }
      // From z here to z - |end_eta - eta| / sqrt(v) at the end, the mass is
      // that of the rest less that beyond the end.
      const double z_end = z - std::fabs(end_eta - eta) / sd;
      stretch.ends_at_node = true;
      stretch.end = log_f.size();
      stretch.end_slope = end_k1 - (end_eta - x) / v;
      stretch.log_mass = value - constant + log_scaled_normal_tail(z) +
        std::log1p(-std::exp(log_scaled_normal_tail(z_end) -
                             log_scaled_normal_tail(z) -
                             0.5 * (z_end * z_end - z * z)));
      stretches.push_back(stretch);
      // The loop's step lands on the end of the stretch.
      j = low - sign;
    }
  }
  double total = 0, first = 0, square = 0, sum_extra = 0;
  for (std::size_t j = 0; j < log_f.size(); ++j) {
    const double w = std::exp(log_f[j] - top);
    total += w;
    first += w * score[j];
    square += w * second[j];
    if constexpr (Terms::has_extra) sum_extra += w * extra[j];
  }
  // In units of the trapezium's weights, which are multiplied by delta below:
  // half of each end node's weight off with its end term delta^2 / 12 |f'|,
  // added where the stretch starts and taken off where it ends, and the
  // closed-form mass, over which k' is constant and k'' zero.
  const auto adjust = [&](const std::size_t node, const double weight) {
    total += weight;
    first += weight * score[node];
    square += weight * second[node];
    if constexpr (Terms::has_extra) sum_extra += weight * extra[node];
  };
  for (const Stretch& stretch : stretches) {
    const double w = std::exp(log_f[stretch.start] - top);
    adjust(stretch.start, -0.5 * w +
             delta / 12.0 * w * std::fabs(stretch.start_slope));
    if (stretch.ends_at_node && stretch.end < log_f.size()) {
      const double w_end = std::exp(log_f[stretch.end] - top);
      adjust(stretch.end, -0.5 * w_end -
               delta / 12.0 * w_end * std::fabs(stretch.end_slope));
    }
    const double mass = std::exp(stretch.log_mass - top) / delta;
    total += mass;
    first += mass * stretch.k1;
    square += mass * stretch.k1 * stretch.k1;
    if constexpr (Terms::has_extra) sum_extra += mass * stretch.extra;
  }
  const double mean = first / total;
  return {top + std::log(total * delta), mean, square / total - mean * mean,
          sum_extra / total};
}

// log Phi(z) for any z, with Phi the standard normal distribution function.
inline double log_normal_cdf(const double z) {
  if (z <= 0.0) return log_scaled_normal_tail(z) - 0.5 * z * z;
  return std::log1p(-std::exp(log_scaled_normal_tail(-z) - 0.5 * z * z));
}

// The rule of custom_cell() for a kernel that is constant below `flat_low`
// and above `flat_high`, as a custom link's is beyond its lattice, when the
// range to search, x - reach to x + reach, passes one of them. The part
// beyond each end it passes is exp(k) times the Gaussian tail there, in closed
// form. The rest, from the larger of x - reach and flat_low to the smaller of
// x + reach and flat_high, is searched by a coarse pass with steps of at most
// `step_cap`, which finds where the integrand is not negligible and the
// largest curvature of k there, and that part is summed by the trapezium rule
// on nodes that include its ends, with a step of at most `step_cap` and at
// most half the narrowest width of the integrand on it, and half weights at
// the ends. With a step cap above 0.25, for a nested pattern, the curvature is
// taken in absolute value, since its kernel bends upwards where it meets its
// constant tail. Where that part reaches flat_low or flat_high, the
// Euler-Maclaurin term of order step^2 is added there, with the kernel's slope
// from inside, since it can jump at the end of the lattice. Beyond the ends k'
// and k'' are zero, or negligible for a nested pattern.
template <class Terms>
inline Moments custom_flat_cell(const Terms& terms, const double x,
                                const double v, const double reach,
                                const double flat_low, const double flat_high,
                                const double step_cap) {
  const double cut = 40.0;
  const double sd = std::sqrt(v);
  const double constant = -0.5 * std::log(2.0 * M_PI * v);
  const double from = std::max(x - reach, flat_low);
  const double to = std::min(x + reach, flat_high);
  double k = 0, k1 = 0, k2 = 0;
  // The pieces beyond the ends: log mass, and k', k'' + k'^2 and the extra
  // quantity there.
  struct Piece {
    double log_mass = -std::numeric_limits<double>::infinity();
    double score = 0, second = 0, extra = 0;
  } pieces[2];
  if (x - reach < flat_low) {
    terms(flat_low, k, k1, k2);
    pieces[0].log_mass = k + log_normal_cdf((flat_low - x) / sd);
    pieces[0].score = k1;
    pieces[0].second = k1 * k1 + k2;
    pieces[0].extra = terms.extra();
  }
  if (x + reach > flat_high) {
    terms(flat_high, k, k1, k2);
    pieces[1].log_mass = k + log_normal_cdf((x - flat_high) / sd);
    pieces[1].score = k1;
    pieces[1].second = k1 * k1 + k2;
    pieces[1].extra = terms.extra();
  }
  double peak = std::max(pieces[0].log_mass, pieces[1].log_mass);
  std::vector<double> fine, score, second, extra;
  double step = 0, low = from, high = to;
  if (from < to) {
    const double coarse = std::min(step_cap, (to - from) / 400.0);
    const long nodes = static_cast<long>(std::ceil((to - from) / coarse));
    std::vector<double> log_f(nodes + 1), worst(nodes + 1);
    double top = -std::numeric_limits<double>::infinity();
    for (long j = 0; j <= nodes; ++j) {
      const double eta = std::min(from + j * coarse, to), e = eta - x;
      terms(eta, k, k1, k2);
      log_f[j] = k - e * e / (2.0 * v) + constant;
      worst[j] = step_cap > 0.25 ? std::fabs(k2) : -k2;
      if (std::isfinite(log_f[j])) top = std::max(top, log_f[j]);
    }
    long first = nodes + 1, last = -1;
    double largest = 0;
    for (long j = 0; j <= nodes; ++j) {
      const double e = std::fabs(std::min(from + j * coarse, to) - x);
      const double slack = (coarse * coarse + 2.0 * coarse * e) / (2.0 * v);
      if (log_f[j] > top - cut - slack) {
        first = std::min(first, j);
        last = std::max(last, j);
        largest = std::max(largest, worst[j]);
      }
    }
    if (first <= last) {
      low = std::max(from, from + (first - 1) * coarse);
      high = std::min(to, from + (last + 1) * coarse);
      const double width = 1.0 / std::sqrt(largest + 1.0 / v);
      double delta = step_cap;
      while (delta > 0.5 * width && delta > 1e-7) delta *= 0.5;
      const long count = std::max(1L, static_cast<long>(std::ceil(
        (high - low) / delta)));
      step = (high - low) / count;
      fine.resize(count + 1);
      score.resize(count + 1);
      second.resize(count + 1);
      extra.resize(Terms::has_extra ? count + 1 : 0);
      for (long j = 0; j <= count; ++j) {
        const double eta = j == count ? high : low + j * step, e = eta - x;
        terms(eta, k, k1, k2);
        fine[j] = k - e * e / (2.0 * v) + constant;
        score[j] = k1;
        second[j] = k1 * k1 + k2;
        if constexpr (Terms::has_extra) extra[j] = terms.extra();
        if (std::isfinite(fine[j])) peak = std::max(peak, fine[j]);
      }
    }
  }
  // Sums in units of exp(peak) times the step, or exp(peak) without nodes.
  const double unit = step > 0 ? step : 1.0;
  double total = 0, sum_first = 0, sum_square = 0, sum_extra = 0;
  const long last = static_cast<long>(fine.size()) - 1;
  for (long j = 0; j <= last; ++j) {
    const double w = (j == 0 || j == last ? 0.5 : 1.0) *
      std::exp(fine[j] - peak);
    total += w;
    sum_first += w * score[j];
    sum_square += w * second[j];
    if constexpr (Terms::has_extra) sum_extra += w * extra[j];
  }
  // The Euler-Maclaurin term -step^2 / 12 (f'(high) - f'(low)) at an end of
  // the flat range, with f' = f (k' - (eta - x) / v) from inside.
  const double inside = 1e-9 * std::max(1.0, std::max(std::fabs(low),
                                                      std::fabs(high)));
  if (last >= 0 && low == flat_low) {
    terms(low + inside, k, k1, k2);
    total += step / 12.0 * std::exp(fine[0] - peak) * (k1 - (low - x) / v);
  }
  if (last >= 0 && high == flat_high) {
    terms(high - inside, k, k1, k2);
    total -= step / 12.0 * std::exp(fine[last] - peak) * (k1 - (high - x) / v);
  }
  for (const Piece& piece : pieces) {
    const double mass = std::exp(piece.log_mass - peak) / unit;
    total += mass;
    sum_first += mass * piece.score;
    sum_square += mass * piece.second;
    if constexpr (Terms::has_extra) sum_extra += mass * piece.extra;
  }
  const double mean = sum_first / total;
  return {peak + std::log(total * unit), mean,
          sum_square / total - mean * mean, sum_extra / total};
}

// The rule for a custom link, which assumes no mode, since concavity cannot be
// established numerically. The integrand exp(k(eta) - (eta - x)^2 / (2 v))
// is negligible where (eta - x)^2 / (2 v) exceeds sup k - k(x) plus the cut,
// which bounds the search. A coarse pass over that range finds where the
// integrand is not negligible and the largest curvature of k there, and the
// trapezium rule then sums that range with a step of at most `step_cap` and at
// most half the narrowest width of the integrand on it. The cap is 0.25, and
// larger only for a nested pattern with an inner variance above 16
// (SiteTable::pattern_step_cap()), whose curvature is then taken in absolute
// value. A record's kernel is
// constant beyond its lattice, and a nested pattern's beyond a wider range
// (SiteTable::pattern_flat_range()), from `flat_low` to `flat_high`. Where the
// range would reach beyond `flat_reach` and past one of them, which needs a
// variance above 150 at least and above 3000 near the kernel's top, the parts
// beyond them are taken in closed form (custom_flat_cell(), with steps of at
// most `step_cap`), so the work no longer grows with v.
constexpr double flat_reach = 500.0;

template <class Terms>
inline Moments custom_cell(const Terms& terms, const double x, const double v,
                           const double kernel_sup,
                           const double flat_low =
                             -std::numeric_limits<double>::infinity(),
                           const double flat_high =
                             std::numeric_limits<double>::infinity(),
                           const double step_cap = 0.25) {
  const double cut = 40.0;
  double k = 0, k1 = 0, k2 = 0;
  terms(x, k, k1, k2);
  const double reach = std::sqrt(2.0 * v * (std::max(kernel_sup - k, 0.0) +
                                            cut + 2.0));
  if (reach > flat_reach && (x - reach < flat_low || x + reach > flat_high)) {
    return custom_flat_cell(terms, x, v, reach, flat_low, flat_high,
                            step_cap);
  }
  const double coarse = std::min(step_cap, reach / 400.0);
  const long nodes = static_cast<long>(std::ceil(reach / coarse));
  std::vector<double> log_f(2 * nodes + 1), worst(2 * nodes + 1);
  double top = -std::numeric_limits<double>::infinity();
  for (long j = -nodes; j <= nodes; ++j) {
    const double e = j * coarse;
    terms(x + e, k, k1, k2);
    log_f[j + nodes] = k - e * e / (2.0 * v);
    worst[j + nodes] = step_cap > 0.25 ? std::fabs(k2) : -k2;
    if (std::isfinite(log_f[j + nodes])) top = std::max(top, log_f[j + nodes]);
  }
  // A peak between coarse nodes can be missed by at most the change of the
  // Gaussian factor over one coarse step.
  long first = 2 * nodes + 1, last = -1;
  double largest = 0;
  for (long j = -nodes; j <= nodes; ++j) {
    const double e = std::fabs(j * coarse);
    const double slack = (coarse * coarse + 2.0 * coarse * e) / (2.0 * v);
    if (log_f[j + nodes] > top - cut - slack) {
      first = std::min(first, j);
      last = std::max(last, j);
      largest = std::max(largest, worst[j + nodes]);
    }
  }
  const double width = 1.0 / std::sqrt(largest + 1.0 / v);
  double delta = step_cap;
  while (delta > 0.5 * width && delta > 1e-7) delta *= 0.5;
  const double from = x + (first - 1) * coarse, to = x + (last + 1) * coarse;
  const long start = static_cast<long>(std::floor(from / delta));
  const long end = static_cast<long>(std::ceil(to / delta));
  const double constant = -0.5 * std::log(2.0 * M_PI * v);
  std::vector<double> fine(end - start + 1), score(end - start + 1),
    second(end - start + 1), extra(Terms::has_extra ? end - start + 1 : 0);
  double peak = -std::numeric_limits<double>::infinity();
  for (long j = start; j <= end; ++j) {
    const double eta = j * delta, e = eta - x;
    terms(eta, k, k1, k2);
    fine[j - start] = k - e * e / (2.0 * v) + constant;
    score[j - start] = k1;
    second[j - start] = k1 * k1 + k2;
    if constexpr (Terms::has_extra) extra[j - start] = terms.extra();
    if (std::isfinite(fine[j - start])) peak = std::max(peak, fine[j - start]);
  }
  double total = 0, sum_first = 0, sum_square = 0, sum_extra = 0;
  for (std::size_t j = 0; j < fine.size(); ++j) {
    const double difference = fine[j] - peak;
    if (!(difference > -cut)) continue;
    const double w = std::exp(difference);
    total += w;
    sum_first += w * score[j];
    sum_square += w * second[j];
    if constexpr (Terms::has_extra) sum_extra += w * extra[j];
  }
  const double mean = sum_first / total;
  return {peak + std::log(total * delta), mean,
          sum_square / total - mean * mean, sum_extra / total};
}

const double grid_low = -50.0, grid_step = 0.25;
const int grid_points = 321;

// The quintic interpolant of tabled values at position t of the interval from
// the grid point `left` (an index into H, H1 and H2), with its derivatives up
// to `order`.
inline Moments interpolate(const std::vector<double>& H,
                           const std::vector<double>& H1,
                           const std::vector<double>& H2,
                           const std::size_t left, const double t,
                           const int order) {
  const double h0 = H[left], h1 = H[left + 1];
  const double d0 = H1[left], d1 = H1[left + 1];
  const double s0 = H2[left], s1 = H2[left + 1];
  const double step = grid_step, step2 = grid_step * grid_step;
  const double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
  const double tail = 10 * t3 - 15 * t4 + 6 * t5;
  Moments out;
  out.h = (1 - tail) * h0 + tail * h1 +
    step * ((t - 6 * t3 + 8 * t4 - 3 * t5) * d0 +
            (-4 * t3 + 7 * t4 - 3 * t5) * d1) +
    step2 * (0.5 * (t2 - 3 * t3 + 3 * t4 - t5) * s0 +
             0.5 * (t3 - 2 * t4 + t5) * s1);
  out.h1 = 0;
  out.h2 = 0;
  if (order >= 1) {
    const double tail1 = 30 * t2 - 60 * t3 + 30 * t4;
    out.h1 = (tail1 * (h1 - h0) +
              step * ((1 - 18 * t2 + 32 * t3 - 15 * t4) * d0 +
                      (-12 * t2 + 28 * t3 - 15 * t4) * d1) +
              step2 * (0.5 * (2 * t - 9 * t2 + 12 * t3 - 5 * t4) * s0 +
                       0.5 * (3 * t2 - 8 * t3 + 5 * t4) * s1)) / step;
  }
  if (order >= 2) {
    const double tail2 = 60 * t - 180 * t2 + 120 * t3;
    out.h2 = (tail2 * (h1 - h0) +
              step * ((-36 * t + 96 * t2 - 60 * t3) * d0 +
                      (-24 * t + 84 * t2 - 60 * t3) * d1) +
              step2 * (0.5 * (2 - 18 * t + 36 * t2 - 20 * t3) * s0 +
                       0.5 * (6 * t - 24 * t2 + 20 * t3) * s1)) / step2;
  }
  return out;
}

struct Table {
  double v = 0;
  std::uint64_t used = 0;
  std::vector<double> H, H1, H2;
  std::vector<unsigned char> filled;
};

struct NestedTable {
  double vo = 0, vi = 0;
  std::uint64_t used = 0;
  std::vector<double> H, H1, H2;
  std::vector<unsigned char> filled;
};

class SiteTable;

// The sum over a pattern's units of h at vi (or of the plain kernel when
// `inner` is null) and its first two derivatives at eta, and with has_extra
// the sum of (h'' + h'^2) / 2, the derivative of the sum in vi.
template <bool Extra>
struct PatternTerms {
  static constexpr bool has_extra = Extra;
  SiteTable& sites;
  Table* inner;
  const int pattern;
  mutable double last_extra = 0.0;
  void operator()(const double eta, double& k, double& k1, double& k2) const;
  double extra() const { return last_extra; }
};

// The cells of one group or location filled so far, for the grid points from
// `first` on.
struct Strip {
  int first = 0;
  std::vector<double> H, H1, H2;
  std::vector<unsigned char> filled;
};

// One integral kept for a group or a location: the strip of its cells at the
// variance v of a group, or the outer variance v and inner variance vi of a
// location, and the key of relative linear predictors it was filled for, a
// group's relative to its first member and a location's relative to its
// record. Each group and location keeps 3, the least recently used replaced
// first.
struct Entry {
  bool kept = false;
  double v = 0, vi = 0;
  std::vector<double> key;
  std::uint64_t used = 0;
  Strip strip;
};
constexpr int entries_kept = 3;
using Entries = std::array<Entry, entries_kept>;

// The integral of a pattern of members at one variance v: the pairs of a
// group's members and their linear predictors relative to a reference member,
// in a canonical order, with the bound on their log-likelihood for a custom
// link. Every group with that pattern shares it, whatever its household or
// location, as households of two children with no deaths and different sexes
// do.
struct Pattern {
  double v = 0;
  std::vector<int> pairs;
  std::vector<double> relative;
  double sup = 0;
  std::uint64_t generation = 0;
  Strip strip;
};

struct KeyHash {
  std::size_t operator()(const std::vector<std::uint64_t>& key) const {
    std::uint64_t hash = 1469598103934665603ULL;
    for (const std::uint64_t word : key) {
      hash ^= word + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    }
    return static_cast<std::size_t>(hash);
  }
};

// A group as its location's integral reads it: its pattern at the inner
// variance, null for a group of one member, which takes the record table, or
// for an inner variance of zero, the shift of its argument, which is the
// relative linear predictor of its pattern's reference member, and its
// members' relative linear predictors.
struct GroupRef {
  int group = 0;
  Pattern* pattern = nullptr;
  double shift = 0;
  const double* relative = nullptr;
};

// The sum over a pattern's members of their kernels at eta + r_j.
struct GroupTerms {
  static constexpr bool has_extra = false;
  const SiteTable& sites;
  const Pattern& pattern;
  void operator()(const double eta, double& k, double& k1, double& k2) const;
  double extra() const { return 0.0; }
};

// The sum over a location's groups of h_g at the inner variance (or of the
// groups' kernels when `inner` is false) and its first two derivatives at
// eta, and with has_extra the sum of (h'' + h'^2) / 2, the derivative of the
// sum in the inner variance. A group of one member takes the record table
// `records` at the inner variance.
template <bool Extra>
struct LocationTerms {
  static constexpr bool has_extra = Extra;
  SiteTable& sites;
  const std::vector<GroupRef>& groups;
  const bool inner;
  Table* records;
  mutable double last_extra = 0.0;
  void operator()(const double eta, double& k, double& k1, double& k2) const;
  double extra() const { return last_extra; }
};

class SiteTable {
 public:
  // `patterns` holds the 0-based pairs of each nested pattern's units and
  // `summed` the 0-based pair of each pattern's units summed into one record.
  // `groups` holds the 0-based pairs of each group's members, whose relative
  // linear predictors are numbered in that order across the groups,
  // `locations` the 0-based groups of each location of groups, and `whole`
  // for each location the 0-based group of all its members together.
  // `concave` says that a custom link is known to be log-concave, as the
  // named probit and cloglog links are, so that the integrals of groups and
  // locations can take the trapezium rule walked out from their mode.
  SiteTable(std::vector<double> y, std::vector<double> m, const bool poisson,
            Rule rule, CustomArrays custom = CustomArrays(),
            const bool has_custom = false,
            std::vector<std::vector<int>> patterns = {},
            std::vector<int> summed = {},
            const std::vector<std::vector<int>>& groups = {},
            std::vector<std::vector<int>> locations = {},
            std::vector<int> whole = {}, const bool concave = false)
      : y_(std::move(y)), m_(std::move(m)), rule_(std::move(rule)),
        custom_(std::move(custom)), patterns_(std::move(patterns)),
        summed_(std::move(summed)), locations_(std::move(locations)),
        whole_(std::move(whole)), concave_(concave) {
    group_start_.assign(1, 0);
    for (const std::vector<int>& members : groups) {
      member_pair_.insert(member_pair_.end(), members.begin(), members.end());
      group_start_.push_back(static_cast<int>(member_pair_.size()));
    }
    kernel_.poisson = poisson;
    if (has_custom) {
      custom_.poisson = poisson;
      kernel_.custom = &custom_;
      // The largest value of each pair's kernel on the lattice, plus one for
      // a peak between lattice points.
      sup_.resize(y_.size());
      for (std::size_t pair = 0; pair < y_.size(); ++pair) {
        double largest = -std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < custom_.A.size(); ++j) {
          largest = std::max(largest, y_[pair] * custom_.A[j] +
                                        m_[pair] * custom_.B[j]);
        }
        sup_[pair] = largest + 1.0 +
          (poisson && m_[pair] > 0 ? y_[pair] * std::log(m_[pair]) : 0.0);
      }
    }
    // Each pattern's bound on its units' log-likelihood, for custom links,
    // and the constant by which the plain kernel of its summed record falls
    // short of the sum of its units' plain kernels, which is zero for the
    // canonical links and y log(m) terms for a custom Poisson link.
    pattern_sup_.assign(patterns_.size(), 0.0);
    pattern_constant_.assign(patterns_.size(), 0.0);
    for (std::size_t p = 0; p < patterns_.size(); ++p) {
      double total = 0;
      for (const int pair : patterns_[p]) {
        if (!sup_.empty()) pattern_sup_[p] += sup_[pair];
        total += plain(pair, 0.0).h;
      }
      pattern_constant_[p] = total - plain(summed_[p], 0.0).h;
    }
    // Each group's and location's bound on its members' log-likelihood, for
    // custom links.
    group_sup_.assign(group_count(), 0.0);
    if (!sup_.empty()) {
      for (int g = 0; g < group_count(); ++g) {
        for (int j = group_start_[g]; j < group_start_[g + 1]; ++j) {
          group_sup_[g] += sup_[member_pair_[j]];
        }
      }
    }
    location_sup_.assign(locations_.size(), 0.0);
    for (std::size_t l = 0; l < locations_.size(); ++l) {
      for (const int g : locations_[l]) location_sup_[l] += group_sup_[g];
    }
    location_entries_.resize(locations_.size());
  }
  SiteTable(const SiteTable&) = delete;
  SiteTable& operator=(const SiteTable&) = delete;

  int pairs() const { return y_.size(); }
  int pattern_count() const { return patterns_.size(); }
  long cells_filled() const { return cells_filled_; }
  int tables_built() const { return tables_built_; }
  long nested_cells_filled() const { return nested_cells_filled_; }
  int nested_tables_built() const { return nested_tables_built_; }
  int group_count() const { return static_cast<int>(group_start_.size()) - 1; }
  int location_count() const { return locations_.size(); }
  int member_count() const { return member_pair_.size(); }
  long group_cells_filled() const { return group_cells_filled_; }
  int group_tables_built() const { return group_tables_built_; }
  long location_cells_filled() const { return location_cells_filled_; }
  int location_tables_built() const { return location_tables_built_; }
  int whole_group(const int location) const { return whole_[location]; }

  // The linear predictors of the members, which the groups and locations
  // resolved from now on take: first each member's relative to its record,
  // then each member's relative to the first member of its group, computed
  // from the differences of their covariates so that a covariate the group's
  // members share adds exactly zero.
  void set_relative(const std::vector<double>& relative) {
    if (static_cast<int>(relative.size()) != 2 * member_count()) {
      throw std::invalid_argument(
        "Two relative linear predictors are needed for each member of a group.");
    }
    if (current_ != nullptr && *current_ == relative) return;
    current_ = std::make_shared<const std::vector<double>>(relative);
  }
  const std::shared_ptr<const std::vector<double>>& relative() const {
    return current_;
  }
  int group_start(const int group) const { return group_start_[group]; }
  int group_size(const int group) const {
    return group_start_[group + 1] - group_start_[group];
  }
  const std::vector<int>& location_groups(const int location) const {
    return locations_[location];
  }

  // Starts the resolution of an evaluator's records: patterns not used by
  // this evaluator or the two before it are dropped once more than about four
  // per group are kept, so that a pattern stays valid while an evaluator
  // that resolved it lives.
  void begin_resolution() {
    ++generation_;
    const std::size_t cap = 4 * static_cast<std::size_t>(group_count()) + 256;
    if (patterns_shared_.size() <= cap) return;
    for (auto it = patterns_shared_.begin(); it != patterns_shared_.end();) {
      if (it->second.generation + 2 < generation_) {
        it = patterns_shared_.erase(it);
      } else {
        ++it;
      }
    }
  }

  // The pattern of a group at v for the current relative linear predictors,
  // its members in the order of their linear predictors relative to the
  // first member of the group, which come from the differences of their
  // covariates, and then of their pairs, with the first in that order as the
  // reference member. In `shift` the reference member's relative linear
  // predictor, so that the group's h at x is the pattern's at x + shift. A
  // household's pattern so stays the same when only coefficients of
  // covariates shared by its members change, and households alike in their
  // members' outcomes and covariate differences share it.
  Pattern& group_pattern(const int group, const double v, double& shift) {
    const int first = group_start_[group];
    const int size = group_size(group);
    const double* within = current_->data() + member_count() + first;
    order_.resize(size);
    std::iota(order_.begin(), order_.end(), 0);
    std::sort(order_.begin(), order_.end(), [&](const int a, const int b) {
      if (within[a] != within[b]) return within[a] < within[b];
      return member_pair_[first + a] < member_pair_[first + b];
    });
    const int reference = order_[0];
    shift = (*current_)[first + reference];
    words_.assign(1, bits(v));
    for (const int j : order_) {
      words_.push_back(static_cast<std::uint64_t>(member_pair_[first + j]));
      words_.push_back(bits(within[j] - within[reference]));
    }
    auto found = patterns_shared_.find(words_);
    if (found == patterns_shared_.end()) {
      Pattern fresh;
      fresh.v = v;
      for (const int j : order_) {
        const int pair = member_pair_[first + j];
        fresh.pairs.push_back(pair);
        fresh.relative.push_back(within[j] - within[reference]);
        if (!sup_.empty()) fresh.sup += sup_[pair];
      }
      found = patterns_shared_.emplace(words_, std::move(fresh)).first;
      ++group_tables_built_;
    }
    found->second.generation = generation_;
    return found->second;
  }

  // The entry of a location of groups at (vo, vi), keyed on the relative
  // linear predictors of all its members, its groups' in order.
  Entry& location_entry(const int location, const double vo,
                        const double vi) {
    key_.clear();
    for (const int g : locations_[location]) {
      const double* relative = current_->data() + group_start_[g];
      key_.insert(key_.end(), relative, relative + group_size(g));
    }
    return resolve(location_entries_[location], vo, vi, key_,
                   location_tables_built_);
  }

  // The sum over a pattern's members of their kernels at eta + r_j.
  Moments pattern_plain(const Pattern& pattern, const double eta) const {
    Moments out{0.0, 0.0, 0.0, 0.0};
    for (std::size_t j = 0; j < pattern.pairs.size(); ++j) {
      double k = 0, k1 = 0, k2 = 0;
      const int pair = pattern.pairs[j];
      kernel_.terms(y_[pair], m_[pair], eta + pattern.relative[j], k, k1, k2);
      out.h += k;
      out.h1 += k1;
      out.h2 += k2;
    }
    return out;
  }

  // The sum over a group's members of their kernels at eta + r_j, with
  // `relative` the members' r_j in order.
  Moments group_plain(const int group, const double* relative,
                      const double eta) const {
    Moments out{0.0, 0.0, 0.0, 0.0};
    const int first = group_start_[group];
    for (int j = first; j < group_start_[group + 1]; ++j) {
      double k = 0, k1 = 0, k2 = 0;
      const int pair = member_pair_[j];
      kernel_.terms(y_[pair], m_[pair], eta + relative[j - first], k, k1, k2);
      out.h += k;
      out.h1 += k1;
      out.h2 += k2;
    }
    return out;
  }

  // The integral of a pattern at t, a group's argument shifted by its
  // shift, with derivatives up to `order`.
  Moments evaluate_group(Pattern& pattern, const double t, const int order) {
    return strip_evaluate(pattern.strip, t, order, group_cells_filled_,
                          [&](const double at) {
                            return group_cell(pattern, at);
                          });
  }

  // H_l at x from the location's entry, whose groups are `groups` at the
  // entry's vi, with the record table `records` at vi.
  Moments evaluate_location(Entry& entry, const std::vector<GroupRef>& groups,
                            Table& records, const int location,
                            const double x, const int order) {
    return strip_evaluate(entry.strip, x, order, location_cells_filled_,
                          [&](const double at) {
                            return location_cell(groups, true, &records,
                                                 location, at, entry.v,
                                                 entry.vi);
                          });
  }

  // h_g at eta as its location reads it: from the group's entry, from the
  // record table `records` for a group of one member, whose h_g is that of
  // its record at eta + r_j, or the group's kernel without `inner`.
  Moments group_moments(const GroupRef& ref, const bool inner,
                        Table* records, const double eta, const int order) {
    if (!inner) return group_plain(ref.group, ref.relative, eta);
    if (ref.pattern == nullptr) {
      return evaluate(*records, member_pair_[group_start_[ref.group]],
                      eta + ref.shift, order);
    }
    return evaluate_group(*ref.pattern, eta + ref.shift, order);
  }

  // The sum over a location's groups of h_g at the inner variance, or of the
  // groups' kernels without `inner`, with derivatives up to `order`, and for
  // order 2 in `extra` the sum of (h'' + h'^2) / 2.
  Moments location_sum(const std::vector<GroupRef>& groups, const bool inner,
                       Table* records, const double eta, const int order) {
    Moments out{0.0, 0.0, 0.0, 0.0};
    for (const GroupRef& ref : groups) {
      const Moments m = group_moments(ref, inner, records, eta, order);
      out.h += m.h;
      out.h1 += m.h1;
      out.h2 += m.h2;
      if (order >= 2) out.extra += 0.5 * (m.h2 + m.h1 * m.h1);
    }
    return out;
  }

  // The tilted mean, under the location's outer integral at vo, of the sum
  // over its groups of (h'' + h'^2) / 2 at the inner variance vi, which is
  // the derivative of H_l in the inner variance, computed directly at x.
  // Without `inner` the inner variance is zero.
  double location_inner_slope(const std::vector<GroupRef>& groups,
                              const bool inner, Table* records,
                              const int location, const double x,
                              const double vo, const double vi) {
    const LocationTerms<true> terms{*this, groups, inner, records};
    if (kernel_.custom != nullptr && !concave_) {
      double low = 0, high = 0, step_cap = 0;
      location_flat_range(groups, inner ? vi : 0.0, low, high, step_cap);
      return custom_cell(terms, x, vo, location_sup_[location], low, high,
                         step_cap).extra;
    }
    return vo <= 0.2 && kernel_.custom == nullptr ?
      aghq_cell(terms, x, vo, rule_).extra :
      uniform_cell(terms, x, vo, inner_step_cap(inner ? vi : 0.0)).extra;
  }

  // The table for v, built empty if it is not kept. Tables are kept in a
  // list, so a reference stays valid until its table is dropped, and the
  // least recently used one is dropped first.
  Table& table(const double v) {
    ++clock_;
    for (Table& kept : tables_) {
      if (kept.v == v) {
        kept.used = clock_;
        return kept;
      }
    }
    if (tables_.size() >= 32) {
      auto oldest = tables_.begin();
      for (auto it = tables_.begin(); it != tables_.end(); ++it) {
        if (it->used < oldest->used) oldest = it;
      }
      tables_.erase(oldest);
    }
    Table fresh;
    fresh.v = v;
    fresh.used = clock_;
    const std::size_t size = static_cast<std::size_t>(grid_points) * pairs();
    fresh.H.assign(size, 0.0);
    fresh.H1.assign(size, 0.0);
    fresh.H2.assign(size, 0.0);
    fresh.filled.assign(size, 0);
    tables_.push_back(std::move(fresh));
    ++tables_built_;
    return tables_.back();
  }

  // The nested table for the pair (vo, vi), built empty if it is not kept.
  NestedTable& nested(const double vo, const double vi) {
    ++clock_;
    for (NestedTable& kept : nested_) {
      if (kept.vo == vo && kept.vi == vi) {
        kept.used = clock_;
        return kept;
      }
    }
    if (nested_.size() >= 16) {
      auto oldest = nested_.begin();
      for (auto it = nested_.begin(); it != nested_.end(); ++it) {
        if (it->used < oldest->used) oldest = it;
      }
      nested_.erase(oldest);
    }
    NestedTable fresh;
    fresh.vo = vo;
    fresh.vi = vi;
    fresh.used = clock_;
    const std::size_t size = static_cast<std::size_t>(grid_points) *
      patterns_.size();
    fresh.H.assign(size, 0.0);
    fresh.H1.assign(size, 0.0);
    fresh.H2.assign(size, 0.0);
    fresh.filled.assign(size, 0);
    nested_.push_back(std::move(fresh));
    ++nested_tables_built_;
    return nested_.back();
  }

  // h and, if `order` is 1 or 2, its derivatives at x for the 0-based pair.
  Moments evaluate(Table& table, const int pair, const double x,
                   const int order) {
    const double v = table.v;
    const double position = (x - grid_low) / grid_step;
    if (!(position >= 0.0 && position <= grid_points - 1.0)) {
      if (!std::isfinite(x)) {
        return {-std::numeric_limits<double>::infinity(), 0.0, 0.0};
      }
      return cell(pair, x, v);
    }
    const int at = std::min(static_cast<int>(position), grid_points - 2);
    const double t = position - at;
    const std::size_t left = static_cast<std::size_t>(pair) * grid_points + at;
    fill(table, pair, at, left);
    fill(table, pair, at + 1, left + 1);
    return interpolate(table.H, table.H1, table.H2, left, t, order);
  }

  // H of the 0-based pattern at x from its nested table, whose inner pairs
  // take the record table `inner` at the table's vi.
  Moments evaluate_nested(NestedTable& table, Table& inner, const int pattern,
                          const double x, const int order) {
    const double position = (x - grid_low) / grid_step;
    if (!(position >= 0.0 && position <= grid_points - 1.0)) {
      if (!std::isfinite(x)) {
        return {-std::numeric_limits<double>::infinity(), 0.0, 0.0};
      }
      return nested_cell(&inner, pattern, x, table.vo);
    }
    const int at = std::min(static_cast<int>(position), grid_points - 2);
    const double t = position - at;
    const std::size_t left =
      static_cast<std::size_t>(pattern) * grid_points + at;
    fill_nested(table, inner, pattern, at, left);
    fill_nested(table, inner, pattern, at + 1, left + 1);
    return interpolate(table.H, table.H1, table.H2, left, t, order);
  }

  // The tilted mean, under the pattern's outer integral at vo, of the sum
  // over its units of (h'' + h'^2) / 2 at vi, which is the derivative of H in
  // vi, computed directly at x. `inner` is null for vi = 0.
  double nested_inner_slope(Table* inner, const int pattern, const double x,
                            const double vo) {
    const PatternTerms<true> terms{*this, inner, pattern};
    if (kernel_.custom != nullptr) {
      double low = 0, high = 0, step_cap = 0;
      pattern_flat_range(inner, pattern, low, high, step_cap);
      return custom_cell(terms, x, vo, pattern_sup_[pattern], low, high,
                         step_cap).extra;
    }
    return vo <= 0.2 ? aghq_cell(terms, x, vo, rule_).extra :
      uniform_cell(terms, x, vo, pattern_step_cap(inner)).extra;
  }

  // The sum over a pattern's units of h at the record table `inner`, or of
  // the plain kernel when it is null, with derivatives up to `order`, and for
  // order 2 in `extra` the sum of (h'' + h'^2) / 2.
  Moments pattern_sum(Table* inner, const int pattern, const double eta,
                      const int order) {
    Moments out{0.0, 0.0, 0.0, 0.0};
    for (const int pair : patterns_[pattern]) {
      const Moments m = inner == nullptr ? plain(pair, eta) :
        evaluate(*inner, pair, eta, order);
      out.h += m.h;
      out.h1 += m.h1;
      out.h2 += m.h2;
      if (order >= 2) out.extra += 0.5 * (m.h2 + m.h1 * m.h1);
    }
    return out;
  }

  // The pattern's units summed into one record with the whole variance vo
  // (the record table `outer`, or the plain kernel when it is null), which is
  // the pattern's H when vi = 0.
  Moments pattern_summed(Table* outer, const int pattern, const double eta,
                         const int order) {
    Moments out = outer == nullptr ? plain(summed_[pattern], eta) :
      evaluate(*outer, summed_[pattern], eta, order);
    out.h += pattern_constant_[pattern];
    return out;
  }

  // The plain kernel of the 0-based pair, for records whose effects are not
  // integrated and for variances held at zero.
  Moments plain(const int pair, const double eta) const {
    Moments out;
    kernel_.terms(y_[pair], m_[pair], eta, out.h, out.h1, out.h2);
    return out;
  }

 private:
  // For a custom link, the range beyond which a nested pattern's kernel, the
  // sum of its units' h at the inner variance vi (or of their kernels when
  // `inner` is null), is constant, and the largest step for its outer
  // integral. Each unit's kernel is constant beyond the lattice, and its h at
  // vi differs from that constant by at most exp(D) Phi(-c) relative, a
  // distance c sqrt(vi) past the lattice, with D the rise of its kernel above
  // its value at that end of the lattice. With D the pattern's rise, which
  // bounds every unit's, and c^2 / 2 = 40 + D, that is about exp(-40) / c.
  // The step is that of pattern_step_cap().
  void pattern_flat_range(const Table* inner, const int pattern, double& low,
                          double& high, double& step_cap) const {
    low = custom_.low;
    high = custom_.low +
      static_cast<double>(custom_.A.size() - 1) * custom_.step;
    step_cap = 0.25;
    if (inner == nullptr || !(inner->v > 0.0)) return;
    double at_low = 0, at_high = 0;
    for (const int pair : patterns_[pattern]) {
      at_low += plain(pair, low).h;
      at_high += plain(pair, high).h;
    }
    const double sd = std::sqrt(inner->v);
    low -= std::sqrt(2.0 * (40.0 + std::max(pattern_sup_[pattern] - at_low,
                                            0.0))) * sd;
    high += std::sqrt(2.0 * (40.0 + std::max(pattern_sup_[pattern] - at_high,
                                             0.0))) * sd;
    step_cap = pattern_step_cap(inner);
  }

  // The largest step of a nested pattern's outer integral. Its kernel is a sum
  // of h at the inner variance vi, and for a log-concave kernel h'' is the
  // tilted variance over vi^2 minus 1 / vi, between -1 / vi and 0, so it
  // varies on the scale of sqrt(vi) and the step may grow to sqrt(vi) / 16.
  // A custom link's kernel need not be log-concave, so custom_cell() also
  // takes its curvature in absolute value with such a step. Up to vi = 16 the
  // step stays at most 0.25.
  static double pattern_step_cap(const Table* inner) {
    return inner_step_cap(inner == nullptr ? 0.0 : inner->v);
  }
  static double inner_step_cap(const double vi) {
    if (!(vi > 0.0)) return 0.25;
    return std::max(0.25, std::sqrt(vi) / 16.0);
  }

  // For a custom link, the range beyond which a pattern's kernel, the sum of
  // its members' kernels at x + r_j, is constant, since each member's kernel
  // is constant beyond the lattice. Its members are in increasing order of
  // r_j.
  void group_flat_range(const Pattern& pattern, double& low,
                        double& high) const {
    low = custom_.low - pattern.relative.back();
    high = custom_.low +
      static_cast<double>(custom_.A.size() - 1) * custom_.step -
      pattern.relative.front();
  }

  // For a custom link, the range beyond which a location's kernel, the sum
  // over its groups of h_g at the inner variance vi (or of their kernels when
  // vi is zero), is constant, and the largest step for its outer integral, as
  // pattern_flat_range() gives them for a nested pattern. The plain sum over
  // all the location's members is constant beyond `low` and `high` below,
  // and each h_g at vi differs from its constant by at most about
  // exp(-40) / c a distance c sqrt(vi) beyond them, with c^2 / 2 = 40 + D and
  // D the location's rise, which bounds each group's.
  void location_flat_range(const std::vector<GroupRef>& groups,
                           const double vi, double& low, double& high,
                           double& step_cap) const {
    double top = -std::numeric_limits<double>::infinity();
    double bottom = std::numeric_limits<double>::infinity();
    double sup = 0;
    for (const GroupRef& ref : groups) {
      for (int j = 0; j < group_size(ref.group); ++j) {
        top = std::max(top, ref.relative[j]);
        bottom = std::min(bottom, ref.relative[j]);
      }
      sup += group_sup_[ref.group];
    }
    low = custom_.low - top;
    high = custom_.low +
      static_cast<double>(custom_.A.size() - 1) * custom_.step - bottom;
    step_cap = 0.25;
    if (!(vi > 0.0)) return;
    double at_low = 0, at_high = 0;
    for (const GroupRef& ref : groups) {
      at_low += group_plain(ref.group, ref.relative, low).h;
      at_high += group_plain(ref.group, ref.relative, high).h;
    }
    const double sd = std::sqrt(vi);
    low -= std::sqrt(2.0 * (40.0 + std::max(sup - at_low, 0.0))) * sd;
    high += std::sqrt(2.0 * (40.0 + std::max(sup - at_high, 0.0))) * sd;
    step_cap = inner_step_cap(vi);
  }

  Moments group_cell(const Pattern& pattern, const double x) const {
    const GroupTerms terms{*this, pattern};
    const double v = pattern.v;
    if (kernel_.custom != nullptr && !concave_) {
      double low = 0, high = 0;
      group_flat_range(pattern, low, high);
      return custom_cell(terms, x, v, pattern.sup, low, high);
    }
    return v <= 0.2 && kernel_.custom == nullptr ?
      aghq_cell(terms, x, v, rule_) : uniform_cell(terms, x, v);
  }

  static std::uint64_t bits(const double value) {
    std::uint64_t out = 0;
    std::memcpy(&out, &value, sizeof(out));
    return out;
  }

  Moments location_cell(const std::vector<GroupRef>& groups, const bool inner,
                        Table* records, const int location, const double x,
                        const double vo, const double vi) {
    const LocationTerms<false> terms{*this, groups, inner, records};
    if (kernel_.custom != nullptr && !concave_) {
      double low = 0, high = 0, step_cap = 0;
      location_flat_range(groups, inner ? vi : 0.0, low, high, step_cap);
      return custom_cell(terms, x, vo, location_sup_[location], low, high,
                         step_cap);
    }
    return vo <= 0.2 && kernel_.custom == nullptr ?
      aghq_cell(terms, x, vo, rule_) :
      uniform_cell(terms, x, vo, inner_step_cap(inner ? vi : 0.0));
  }

  // The entry of `kept` for (v, vi) and `key`, filled empty in the least
  // recently used slot if it is not kept.
  Entry& resolve(Entries& kept, const double v, const double vi,
                 const std::vector<double>& key, int& built) {
    ++clock_;
    for (Entry& entry : kept) {
      if (entry.kept && entry.v == v && entry.vi == vi && entry.key == key) {
        entry.used = clock_;
        return entry;
      }
    }
    Entry* slot = &kept[0];
    for (Entry& entry : kept) {
      if (!entry.kept) {
        slot = &entry;
        break;
      }
      if (entry.used < slot->used) slot = &entry;
    }
    *slot = Entry();
    slot->kept = true;
    slot->v = v;
    slot->vi = vi;
    slot->key = key;
    slot->used = clock_;
    ++built;
    return *slot;
  }

  // The strip extended, if needed, to hold the grid points from `low` to
  // `high`, with a margin of 8 points where it grows.
  static void cover(Strip& strip, const int low, const int high) {
    const int first = strip.first;
    const int end = first + static_cast<int>(strip.H.size());
    if (!strip.H.empty() && low >= first && high < end) return;
    const int from = strip.H.empty() || low < first ?
      std::max(0, low - 8) : first;
    const int to = strip.H.empty() || high >= end ?
      std::min(grid_points, high + 9) : end;
    const std::size_t size = static_cast<std::size_t>(to - from);
    std::vector<double> H(size, 0.0), H1(size, 0.0), H2(size, 0.0);
    std::vector<unsigned char> filled(size, 0);
    for (std::size_t j = 0; j < strip.H.size(); ++j) {
      const std::size_t at = j + static_cast<std::size_t>(first - from);
      H[at] = strip.H[j];
      H1[at] = strip.H1[j];
      H2[at] = strip.H2[j];
      filled[at] = strip.filled[j];
    }
    strip.first = from;
    strip.H.swap(H);
    strip.H1.swap(H1);
    strip.H2.swap(H2);
    strip.filled.swap(filled);
  }

  // The interpolant at x through the strip's cells, filled by `cell` where
  // they are not yet, or the cell itself beyond the grid. `counter` counts
  // the cells filled.
  template <class Cell>
  static Moments strip_evaluate(Strip& strip, const double x, const int order,
                                long& counter, const Cell& cell) {
    const double position = (x - grid_low) / grid_step;
    if (!(position >= 0.0 && position <= grid_points - 1.0)) {
      if (!std::isfinite(x)) {
        return {-std::numeric_limits<double>::infinity(), 0.0, 0.0};
      }
      return cell(x);
    }
    const int at = std::min(static_cast<int>(position), grid_points - 2);
    const double t = position - at;
    cover(strip, at, at + 1);
    for (int g = at; g <= at + 1; ++g) {
      const std::size_t index = static_cast<std::size_t>(g - strip.first);
      if (strip.filled[index]) continue;
      const Moments filled = cell(grid_low + g * grid_step);
      strip.H[index] = filled.h;
      strip.H1[index] = filled.h1;
      strip.H2[index] = filled.h2;
      strip.filled[index] = 1;
      ++counter;
    }
    return interpolate(strip.H, strip.H1, strip.H2,
                       static_cast<std::size_t>(at - strip.first), t, order);
  }

  Moments cell(const int pair, const double x, const double v) const {
    const RecordTerms terms{kernel_, y_[pair], m_[pair]};
    if (kernel_.custom != nullptr) {
      // Beyond the lattice of the custom kernel it is constant.
      const double high = custom_.low +
        static_cast<double>(custom_.A.size() - 1) * custom_.step;
      return custom_cell(terms, x, v, sup_[pair], custom_.low, high);
    }
    return v <= 0.2 ? aghq_cell(terms, x, v, rule_) :
      uniform_cell(terms, x, v);
  }

  Moments nested_cell(Table* inner, const int pattern, const double x,
                      const double vo) {
    const PatternTerms<false> terms{*this, inner, pattern};
    if (kernel_.custom != nullptr) {
      double low = 0, high = 0, step_cap = 0;
      pattern_flat_range(inner, pattern, low, high, step_cap);
      return custom_cell(terms, x, vo, pattern_sup_[pattern], low, high,
                         step_cap);
    }
    return vo <= 0.2 ? aghq_cell(terms, x, vo, rule_) :
      uniform_cell(terms, x, vo, pattern_step_cap(inner));
  }

  void fill(Table& table, const int pair, const int g, const std::size_t at) {
    if (table.filled[at]) return;
    const Moments filled = cell(pair, grid_low + g * grid_step, table.v);
    table.H[at] = filled.h;
    table.H1[at] = filled.h1;
    table.H2[at] = filled.h2;
    table.filled[at] = 1;
    ++cells_filled_;
  }

  void fill_nested(NestedTable& table, Table& inner, const int pattern,
                   const int g, const std::size_t at) {
    if (table.filled[at]) return;
    const Moments filled = nested_cell(&inner, pattern,
                                       grid_low + g * grid_step, table.vo);
    table.H[at] = filled.h;
    table.H1[at] = filled.h1;
    table.H2[at] = filled.h2;
    table.filled[at] = 1;
    ++nested_cells_filled_;
  }

  std::vector<double> y_, m_;
  Rule rule_;
  CustomArrays custom_;
  Kernel kernel_;
  std::vector<double> sup_;
  std::vector<std::vector<int>> patterns_;
  std::vector<int> summed_;
  std::vector<double> pattern_sup_, pattern_constant_;
  std::list<Table> tables_;
  std::list<NestedTable> nested_;
  std::uint64_t clock_ = 0;
  long cells_filled_ = 0, nested_cells_filled_ = 0;
  int tables_built_ = 0, nested_tables_built_ = 0;
  // The groups: each one's members are member_pair_[group_start_[g]] up to
  // group_start_[g + 1], numbered as their relative linear predictors are.
  std::vector<int> member_pair_, group_start_;
  std::vector<std::vector<int>> locations_;
  std::vector<int> whole_;
  std::vector<double> group_sup_, location_sup_;
  std::vector<Entries> location_entries_;
  std::vector<double> key_;
  bool concave_ = false;
  // The patterns of groups, keyed on the bits of their variance, pairs and
  // relative linear predictors, with the scratch of their keys.
  std::unordered_map<std::vector<std::uint64_t>, Pattern, KeyHash>
    patterns_shared_;
  std::vector<std::uint64_t> words_;
  std::vector<int> order_;
  std::uint64_t generation_ = 0;
  // The current relative linear predictors of the members.
  std::shared_ptr<const std::vector<double>> current_;
  long group_cells_filled_ = 0, location_cells_filled_ = 0;
  int group_tables_built_ = 0, location_tables_built_ = 0;
};

template <bool Extra>
inline void PatternTerms<Extra>::operator()(const double eta, double& k,
                                            double& k1, double& k2) const {
  const Moments sum = sites.pattern_sum(inner, pattern, eta, 2);
  k = sum.h;
  k1 = sum.h1;
  k2 = sum.h2;
  last_extra = sum.extra;
}

inline void GroupTerms::operator()(const double eta, double& k, double& k1,
                                   double& k2) const {
  const Moments sum = sites.pattern_plain(pattern, eta);
  k = sum.h;
  k1 = sum.h1;
  k2 = sum.h2;
}

template <bool Extra>
inline void LocationTerms<Extra>::operator()(const double eta, double& k,
                                             double& k1, double& k2) const {
  const Moments sum = sites.location_sum(groups, inner, records, eta, 2);
  k = sum.h;
  k1 = sum.h1;
  k2 = sum.h2;
  last_extra = sum.extra;
}

// How the engine's records are evaluated. Each record is plain (kind 0, the
// plain kernel of pair `index`), a record (kind 1, pair `index` at the
// variance of class `a`, plain if that is zero), a nested pattern (kind 2,
// pattern `index` with the outer variance of class `a` and the inner variance
// of class `b`), a group (kind 3, group `index` at the variance of class `a`,
// the plain sum of its members' kernels if that is zero) or a location of
// groups (kind 4, location `index` with the outer variance of class `a` and
// the groups' variance of class `b`). Indices are 0-based.
struct Unit {
  int kind, index, a, b;
};

// The likelihood of the engine's records at one set of class variances. The
// tables are found once, on construction, and stay valid while the evaluator
// lives, since at most a handful are found and the least recently used are
// dropped first. Groups and locations of groups take the current state of
// the relative linear predictors (SiteTable::set_relative()), which must be
// set.
class Evaluator {
 public:
  Evaluator(SiteTable& sites, const std::vector<double>& classes,
            const std::vector<Unit>& units)
      : sites_(sites), units_(units) {
    std::vector<Table*> by_class(classes.size(), nullptr);
    auto record_table = [&](const int c) -> Table* {
      if (!(classes[c] > 0)) return nullptr;
      if (by_class[c] == nullptr) by_class[c] = &sites_.table(classes[c]);
      return by_class[c];
    };
    resolved_.resize(units_.size());
    for (std::size_t r = 0; r < units_.size(); ++r) {
      const Unit& unit = units_[r];
      Resolved& out = resolved_[r];
      out.index = unit.index;
      if (unit.kind == 0) {
        out.mode = 0;
      } else if (unit.kind == 1) {
        out.table = record_table(unit.a);
        out.mode = out.table == nullptr ? 0 : 1;
      } else if (unit.kind == 2) {
        const double vo = classes[unit.a], vi = classes[unit.b];
        out.vo = vo;
        if (vo > 0 && vi > 0) {
          out.mode = 2;
          out.table = record_table(unit.b);
          out.nested = nested_table(vo, vi);
        } else if (vo > 0) {
          out.mode = 3;
          out.table = record_table(unit.a);
        } else {
          out.mode = 4;
          out.table = record_table(unit.b);
        }
      } else if (unit.kind == 3 || unit.kind == 4) {
        if (relative_ == nullptr) {
          relative_ = sites_.relative();
          if (relative_ == nullptr) {
            throw std::logic_error(
              "The relative linear predictors of the groups are not set.");
          }
          sites_.begin_resolution();
        }
        if (unit.kind == 3) {
          const double v = classes[unit.a];
          out.mode = v > 0 ? 5 : 6;
          if (v > 0) {
            out.pattern = &sites_.group_pattern(unit.index, v, out.shift);
          } else {
            out.relative = relative_->data() + sites_.group_start(unit.index);
          }
        } else {
          const double vo = classes[unit.a], vi = classes[unit.b];
          out.vo = vo;
          out.vi = vi;
          if (vo > 0 && vi == 0) {
            // With the groups' variance at zero the location is one group
            // of all its members at vo.
            out.mode = 8;
            out.pattern = &sites_.group_pattern(sites_.whole_group(unit.index),
                                                vo, out.shift);
          } else {
            out.mode = vo > 0 ? 7 : vi > 0 ? 9 : 10;
          }
          // The location's groups, with their entries at vi.
          for (const int g : sites_.location_groups(unit.index)) {
            GroupRef ref;
            ref.group = g;
            ref.relative = relative_->data() + sites_.group_start(g);
            ref.shift = ref.relative[0];
            if (vi > 0 && sites_.group_size(g) > 1) {
              ref.pattern = &sites_.group_pattern(g, vi, ref.shift);
            }
            out.groups.push_back(ref);
          }
          if (vi > 0) out.table = record_table(unit.b);
          if (out.mode == 7) {
            out.entry = &sites_.location_entry(unit.index, vo, vi);
          }
        }
      } else {
        throw std::invalid_argument("A record has an unknown kind.");
      }
    }
  }

  Moments moments(const int record, const double eta, const int order) {
    Resolved& r = resolved_[record];
    switch (r.mode) {
    case 0:
      return sites_.plain(r.index, eta);
    case 1:
      return sites_.evaluate(*r.table, r.index, eta, order);
    case 2:
      return sites_.evaluate_nested(*r.nested, *r.table, r.index, eta, order);
    case 3:
      return sites_.pattern_summed(r.table, r.index, eta, order);
    case 4:
      return sites_.pattern_sum(r.table, r.index, eta, order);
    case 5:
      return sites_.evaluate_group(*r.pattern, eta + r.shift, order);
    case 6:
      return sites_.group_plain(r.index, r.relative, eta);
    case 7:
      return sites_.evaluate_location(*r.entry, r.groups, *r.table, r.index,
                                      eta, order);
    case 8:
      return sites_.evaluate_group(*r.pattern, eta + r.shift, order);
    default:
      // Modes 9 and 10, the sum of the groups' h at their variance, or of
      // their kernels.
      return sites_.location_sum(r.groups, r.mode == 9, r.table, eta, order);
    }
  }

  double log_likelihood(const int record, const double eta) {
    return moments(record, eta, 0).h;
  }

  // The derivative of a nested pattern's H, or of a location's H_l, in its
  // inner variance at eta, and zero for the other records.
  double inner_slope(const int record, const double eta) {
    const Resolved& r = resolved_[record];
    if (r.mode == 2 || r.mode == 3) {
      // With vi = 0 the inner variance's slope is that of the plain kernels.
      return sites_.nested_inner_slope(r.mode == 2 ? r.table : nullptr,
                                       r.index, eta, r.vo);
    }
    if (r.mode == 4) return sites_.pattern_sum(r.table, r.index, eta, 2).extra;
    if (r.mode == 7 || r.mode == 8) {
      // With the groups' variance at zero the slope is that of their kernels.
      return sites_.location_inner_slope(r.groups, r.mode == 7, r.table,
                                         r.index, eta, r.vo, r.vi);
    }
    if (r.mode == 9 || r.mode == 10) {
      return sites_.location_sum(r.groups, r.mode == 9, r.table, eta,
                                 2).extra;
    }
    return 0.0;
  }

 private:
  struct Resolved {
    int mode = 0, index = 0;
    double vo = 0;
    Table* table = nullptr;
    NestedTable* nested = nullptr;
    double vi = 0, shift = 0;
    // A group's pattern, a location's entry, a plain group's relative linear
    // predictors, and a location's groups.
    Pattern* pattern = nullptr;
    Entry* entry = nullptr;
    const double* relative = nullptr;
    std::vector<GroupRef> groups;
  };

  NestedTable* nested_table(const double vo, const double vi) {
    for (std::size_t k = 0; k < nested_keys_.size(); ++k) {
      if (nested_keys_[k].first == vo && nested_keys_[k].second == vi) {
        return nested_found_[k];
      }
    }
    NestedTable* found = &sites_.nested(vo, vi);
    nested_keys_.emplace_back(vo, vi);
    nested_found_.push_back(found);
    return found;
  }

  SiteTable& sites_;
  const std::vector<Unit>& units_;
  std::vector<Resolved> resolved_;
  std::vector<std::pair<double, double>> nested_keys_;
  std::vector<NestedTable*> nested_found_;
  // The relative linear predictors of the groups' members, kept while the
  // evaluator lives.
  std::shared_ptr<const std::vector<double>> relative_;
};

}  // namespace riskmap_site

#endif
