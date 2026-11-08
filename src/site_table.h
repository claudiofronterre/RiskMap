#ifndef RISKMAP_SITE_TABLE_H
#define RISKMAP_SITE_TABLE_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include "likelihood_kernel.h"

// The log-likelihood of a record with its nugget integrated out,
//   h(x; v) = log integral of exp(k(x + e)) N(e; 0, v) de,
// for the binomial and Poisson kernels k of any link, with y events of m
// trials or exposure m. A table holds h and its first two derivatives in x
// on the grid -50, -49.75, ..., 30 for each distinct pair (y, m) and one value
// of v, and the quintic Hermite interpolant through them gives h with two
// continuous derivatives between grid points. A cell is filled the first time
// it is needed, and it is a pure function of (y, m, x, v), so the result does
// not depend on the order in which cells are filled. Beyond the grid the
// integral is evaluated directly. Up to 16 tables are kept, keyed on the exact
// value of v, and the least recently used is dropped.
//
// For the canonical links the rule is adaptive Gauss-Hermite quadrature with
// 20 nodes up to v = 0.2, and above it the trapezium rule on the lattice of
// multiples of a step of at most 0.25 and at most half the integrand's width at
// its mode (Trefethen and Weideman, 2014), walked out from the mode until a term
// is negligible, which is safe because the integrand is log-concave. A custom
// link takes the trapezium rule of custom_cell(), which assumes no mode. The derivatives are those of
// Stein's identity, h' = E k' and h'' = E k'' + Var k' under the tilted
// distribution, which have no cancellation as v goes to zero.

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
    if (poisson) {
      const double mean = m * std::exp(eta);
      k = y * eta - mean;
      k1 = y - mean;
      k2 = -mean;
    } else {
      const double p = logistic(eta);
      k = y * eta - m * softplus(eta);
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

struct Moments {
  double h, h1, h2;
};

struct Rule {
  std::vector<double> z;      // nodes of Gauss-Hermite for exp(-z^2)
  std::vector<double> log_w;  // their log weights
};

// The mode of k(eta) - (eta - x)^2 / (2 v), by Newton's method from x with
// steps of at most 3, and the curvature of the integrand there.
inline double integrand_mode(const Kernel& kernel, const double y,
                             const double m, const double x, const double v,
                             double& curvature) {
  double t = x, k = 0, k1 = 0, k2 = 0;
  for (int iteration = 0; iteration < 200; ++iteration) {
    kernel.terms(y, m, t, k, k1, k2);
    double step = (k1 - (t - x) / v) / (-k2 + 1.0 / v);
    step = std::max(std::min(step, 3.0), -3.0);
    t += step;
    if (std::fabs(step) < 1e-11) break;
  }
  kernel.terms(y, m, t, k, k1, k2);
  curvature = -k2 + 1.0 / v;
  return t;
}

inline Moments aghq_cell(const Kernel& kernel, const double y, const double m,
                         const double x, const double v, const Rule& rule) {
  double curvature = 0;
  const double t = integrand_mode(kernel, y, m, x, v, curvature);
  const double a = 1.0 / std::sqrt(curvature);
  const double constant = -0.5 * std::log(2.0 * M_PI * v);
  const int nodes = rule.z.size();
  std::vector<double> log_f(nodes), score(nodes), second(nodes);
  double top = -std::numeric_limits<double>::infinity();
  for (int j = 0; j < nodes; ++j) {
    const double eta = t + M_SQRT2 * a * rule.z[j];
    const double e = eta - x;
    double k = 0, k1 = 0, k2 = 0;
    kernel.terms(y, m, eta, k, k1, k2);
    log_f[j] = k - e * e / (2.0 * v) + constant + rule.log_w[j] +
      rule.z[j] * rule.z[j];
    score[j] = k1;
    second[j] = k1 * k1 + k2;
    top = std::max(top, log_f[j]);
  }
  double total = 0, first = 0, square = 0;
  for (int j = 0; j < nodes; ++j) {
    const double w = std::exp(log_f[j] - top);
    total += w;
    first += w * score[j];
    square += w * second[j];
  }
  const double mean = first / total;
  return {top + std::log(M_SQRT2 * a * total), mean, square / total - mean * mean};
}

inline Moments uniform_cell(const Kernel& kernel, const double y,
                            const double m, const double x, const double v) {
  double curvature = 0;
  const double t = integrand_mode(kernel, y, m, x, v, curvature);
  const double width = 1.0 / std::sqrt(curvature);
  double delta = 0.25;
  while (delta > 0.5 * width && delta > 1e-6) delta *= 0.5;
  const double constant = -0.5 * std::log(2.0 * M_PI * v);
  const double cut = 40.0;
  const long centre = std::lround(t / delta);
  std::vector<double> log_f, score, second;
  double top = -std::numeric_limits<double>::infinity();
  for (int direction = 0; direction < 2; ++direction) {
    for (long j = direction == 0 ? centre : centre - 1;; j += direction == 0 ? 1 : -1) {
      const double eta = j * delta;
      const double e = eta - x;
      double k = 0, k1 = 0, k2 = 0;
      kernel.terms(y, m, eta, k, k1, k2);
      const double value = k - e * e / (2.0 * v) + constant;
      if (!std::isfinite(value)) break;
      top = std::max(top, value);
      log_f.push_back(value);
      score.push_back(k1);
      second.push_back(k1 * k1 + k2);
      if (value < top - cut || log_f.size() > 4000000) break;
    }
  }
  double total = 0, first = 0, square = 0;
  for (std::size_t j = 0; j < log_f.size(); ++j) {
    const double w = std::exp(log_f[j] - top);
    total += w;
    first += w * score[j];
    square += w * second[j];
  }
  const double mean = first / total;
  return {top + std::log(total * delta), mean, square / total - mean * mean};
}

// The rule for a custom link, which assumes no mode, since concavity cannot be
// established numerically. The integrand exp(k(eta) - (eta - x)^2 / (2 v))
// is negligible where (eta - x)^2 / (2 v) exceeds sup k - k(x) plus the cut,
// which bounds the search. A coarse pass over that range finds where the
// integrand is not negligible and the largest curvature of k there, and the
// trapezium rule then sums that range with a step of at most 0.25 and at most
// half the narrowest width of the integrand on it.
inline Moments custom_cell(const Kernel& kernel, const double y,
                           const double m, const double x, const double v,
                           const double kernel_sup) {
  const double cut = 40.0;
  double k = 0, k1 = 0, k2 = 0;
  kernel.terms(y, m, x, k, k1, k2);
  const double reach = std::sqrt(2.0 * v * (std::max(kernel_sup - k, 0.0) +
                                            cut + 2.0));
  const double coarse = std::min(0.25, reach / 400.0);
  const long nodes = static_cast<long>(std::ceil(reach / coarse));
  std::vector<double> log_f(2 * nodes + 1), worst(2 * nodes + 1);
  double top = -std::numeric_limits<double>::infinity();
  for (long j = -nodes; j <= nodes; ++j) {
    const double e = j * coarse;
    kernel.terms(y, m, x + e, k, k1, k2);
    log_f[j + nodes] = k - e * e / (2.0 * v);
    worst[j + nodes] = -k2;
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
  double delta = 0.25;
  while (delta > 0.5 * width && delta > 1e-7) delta *= 0.5;
  const double from = x + (first - 1) * coarse, to = x + (last + 1) * coarse;
  const long start = static_cast<long>(std::floor(from / delta));
  const long end = static_cast<long>(std::ceil(to / delta));
  const double constant = -0.5 * std::log(2.0 * M_PI * v);
  std::vector<double> fine(end - start + 1), score(end - start + 1),
    second(end - start + 1);
  double peak = -std::numeric_limits<double>::infinity();
  for (long j = start; j <= end; ++j) {
    const double eta = j * delta, e = eta - x;
    kernel.terms(y, m, eta, k, k1, k2);
    fine[j - start] = k - e * e / (2.0 * v) + constant;
    score[j - start] = k1;
    second[j - start] = k1 * k1 + k2;
    if (std::isfinite(fine[j - start])) peak = std::max(peak, fine[j - start]);
  }
  double total = 0, sum_first = 0, sum_square = 0;
  for (std::size_t j = 0; j < fine.size(); ++j) {
    const double difference = fine[j] - peak;
    if (!(difference > -cut)) continue;
    const double w = std::exp(difference);
    total += w;
    sum_first += w * score[j];
    sum_square += w * second[j];
  }
  const double mean = sum_first / total;
  return {peak + std::log(total * delta), mean,
          sum_square / total - mean * mean};
}

const double grid_low = -50.0, grid_step = 0.25;
const int grid_points = 321;

struct Table {
  double v = 0;
  std::uint64_t used = 0;
  std::vector<double> H, H1, H2;
  std::vector<unsigned char> filled;
};

class SiteTable {
 public:
  SiteTable(std::vector<double> y, std::vector<double> m, const bool poisson,
            Rule rule, CustomArrays custom = CustomArrays(),
            const bool has_custom = false)
      : y_(std::move(y)), m_(std::move(m)), rule_(std::move(rule)),
        custom_(std::move(custom)) {
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
  }
  SiteTable(const SiteTable&) = delete;
  SiteTable& operator=(const SiteTable&) = delete;

  int pairs() const { return y_.size(); }
  long cells_filled() const { return cells_filled_; }
  int tables_built() const { return tables_built_; }

  // The table for v, built empty if it is not kept.
  Table& table(const double v) {
    ++clock_;
    for (Table& kept : tables_) {
      if (kept.v == v) {
        kept.used = clock_;
        return kept;
      }
    }
    if (tables_.size() >= 16) {
      std::size_t oldest = 0;
      for (std::size_t i = 1; i < tables_.size(); ++i) {
        if (tables_[i].used < tables_[oldest].used) oldest = i;
      }
      tables_.erase(tables_.begin() + oldest);
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
    const double h0 = table.H[left], h1 = table.H[left + 1];
    const double d0 = table.H1[left], d1 = table.H1[left + 1];
    const double s0 = table.H2[left], s1 = table.H2[left + 1];
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

  // The plain kernel of the 0-based pair, for records whose nugget is not
  // integrated and for a nugget held at zero.
  Moments plain(const int pair, const double eta) const {
    Moments out;
    kernel_.terms(y_[pair], m_[pair], eta, out.h, out.h1, out.h2);
    return out;
  }

 private:
  Moments cell(const int pair, const double x, const double v) const {
    if (kernel_.custom != nullptr) {
      return custom_cell(kernel_, y_[pair], m_[pair], x, v, sup_[pair]);
    }
    return v <= 0.2 ? aghq_cell(kernel_, y_[pair], m_[pair], x, v, rule_) :
      uniform_cell(kernel_, y_[pair], m_[pair], x, v);
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

  std::vector<double> y_, m_;
  Rule rule_;
  CustomArrays custom_;
  Kernel kernel_;
  std::vector<double> sup_;
  std::vector<Table> tables_;
  std::uint64_t clock_ = 0;
  long cells_filled_ = 0;
  int tables_built_ = 0;
};

// The log-likelihood of one record at eta: the table for a positive `pair`
// (1-based), the plain kernel of pair -pair for a negative one, and the plain
// kernel throughout when `table` is null (a nugget held at zero).
inline double record_log_likelihood(SiteTable& sites, Table* table,
                                    const int pair, const double eta) {
  if (table == nullptr || pair < 0) {
    return sites.plain(std::abs(pair) - 1, eta).h;
  }
  return sites.evaluate(*table, pair - 1, eta, 0).h;
}

}  // namespace riskmap_site

#endif
