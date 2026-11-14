// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppEigen.h>
#include <algorithm>
#include <limits>
#include <memory>
#include <vector>
#include "likelihood_kernel.h"

using Eigen::LLT;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::IntegerMatrix;
using Rcpp::List;
using Rcpp::LogicalVector;
using Rcpp::Named;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;


// Log-likelihood only, for every row of eta (draws x observations), under
// the link `link` (a LinkCode of likelihood_kernel.h). The link is resolved
// once and no score or curvature matrices are allocated. A canonical Poisson
// mean that overflows gives -Inf for that draw. Under a named link the
// observations of a draw are summed in long double, as rowSums() sums them in
// compute_link_log_likelihood_rows(), so that the rows are R's to the last
// bit, and a draw with a probability outside [0, 1] gives -Inf.
// [[Rcpp::export]]
NumericVector cpp_link_log_likelihood(const NumericMatrix eta_r,
                                      const NumericVector y_r,
                                      const NumericVector units_r,
                                      const int link) {
  if (!is_link_code(link)) {
    Rcpp::stop("The link code %d has no compiled kernel.", link);
  }
  const int n_samples = eta_r.nrow();
  const int n = eta_r.ncol();
  NumericVector log_likelihood(n_samples);
  const double* eta = eta_r.begin();
  double* total = log_likelihood.begin();
  if (!is_canonical_link(link)) {
    std::vector<long double> sum(n_samples, 0.0L);
    std::vector<char> outside(n_samples, 0);
    for (int observation = 0; observation < n; ++observation) {
      const double y = y_r[observation];
      const double units = units_r[observation];
      const double* column =
        eta + static_cast<R_xlen_t>(observation) * n_samples;
      for (int sample = 0; sample < n_samples; ++sample) {
        bool out = false;
        sum[sample] += named_link_log_likelihood(column[sample], y, units,
                                                 link, out);
        if (out) outside[sample] = 1;
      }
    }
    for (int sample = 0; sample < n_samples; ++sample) {
      total[sample] = outside[sample] ?
        -std::numeric_limits<double>::infinity() :
        static_cast<double>(sum[sample]);
    }
    return log_likelihood;
  }
  const bool poisson = link == link_log;
  for (int observation = 0; observation < n; ++observation) {
    const double y = y_r[observation];
    const double units = units_r[observation];
    const double* column = eta + static_cast<R_xlen_t>(observation) * n_samples;
    for (int sample = 0; sample < n_samples; ++sample) {
      total[sample] += canonical_log_likelihood(column[sample], y, units,
                                                poisson);
    }
  }
  if (poisson) {
    for (int sample = 0; sample < n_samples; ++sample) {
      if (!std::isfinite(total[sample]))
        total[sample] = -std::numeric_limits<double>::infinity();
    }
  }
  return log_likelihood;
}


namespace {

// The latent effects of the sampler in its internal order. The first `n_u`
// are the locations and the levels of every re() term but one, with a dense
// prior precision. The last `n_d` are the levels of the re() term with most
// levels, whose block of the Laplace precision is diagonal because each
// record has exactly one of its levels. `columns` gives, for each record, the
// latent effects of its linear predictor (0-based), with the column of that
// term last.
struct LatentStructure {
  int n_records;
  int width;
  int width_u;
  int n_u;
  int n_d;
  const int* columns;
  // Records grouped by their level of the diagonal term.
  std::vector<int> level_start;
  std::vector<int> level_records;

  LatentStructure(const IntegerMatrix& columns_r, const int n_u_, const int n_d_)
    : n_records(columns_r.nrow()), width(columns_r.ncol()),
      width_u(columns_r.ncol() - (n_d_ > 0 ? 1 : 0)), n_u(n_u_), n_d(n_d_),
      columns(columns_r.begin()) {
    if (n_d == 0) return;
    level_start.assign(n_d + 1, 0);
    for (int r = 0; r < n_records; ++r) ++level_start[level(r) + 1];
    for (int h = 0; h < n_d; ++h) level_start[h + 1] += level_start[h];
    level_records.resize(n_records);
    std::vector<int> next(level_start.begin(), level_start.end() - 1);
    for (int r = 0; r < n_records; ++r) level_records[next[level(r)]++] = r;
  }

  int column(const int record, const int j) const {
    return columns[static_cast<R_xlen_t>(j) * n_records + record];
  }
  int level(const int record) const {
    return column(record, width - 1) - n_u;
  }
  int dimension() const { return n_u + n_d; }
};

// The Laplace precision Q = P + A' W A of the latent effects, with prior
// precision P (dense on the first block, diagonal on the second) and record
// weights W, held as the Cholesky factor R'R of the Schur complement
// S = Q_uu - Q_ud D^-1 Q_du of its diagonal block D = Q_dd. The factor of
// Q^-1 used by the chain is
//   L = [R^-1, 0; -D^-1 Q_du R^-1, D^-1/2],   L L' = Q^-1,
// so a product with L or L' is one triangular solve with R and passes over
// the records, instead of dense algebra in all the latent effects.
class StructuredFactor {
 public:
  StructuredFactor(const LatentStructure& structure,
                   const Map<const MatrixXd>& prior_u,
                   const Map<const VectorXd>& prior_d,
                   const Map<const VectorXd>& weights)
    : s_(structure), weights_(weights) {
    MatrixXd schur = prior_u;
    for (int r = 0; r < s_.n_records; ++r) {
      for (int a = 0; a < s_.width_u; ++a) {
        for (int b = 0; b < s_.width_u; ++b) {
          schur(s_.column(r, a), s_.column(r, b)) += weights_[r];
        }
      }
    }
    if (s_.n_d > 0) {
      diagonal_ = prior_d;
      for (int r = 0; r < s_.n_records; ++r) diagonal_[s_.level(r)] += weights_[r];
      valid_ = (diagonal_.array() > 0.0).all() && diagonal_.allFinite();
      if (!valid_) return;
      // Subtract q_h q_h' / D_h, with q_h the column of Q_ud of level h,
      // accumulated sparsely over the records of that level.
      VectorXd accumulator = VectorXd::Zero(s_.n_u);
      std::vector<char> marked(s_.n_u, 0);
      std::vector<int> touched;
      for (int h = 0; h < s_.n_d; ++h) {
        touched.clear();
        for (int k = s_.level_start[h]; k < s_.level_start[h + 1]; ++k) {
          const int r = s_.level_records[k];
          for (int a = 0; a < s_.width_u; ++a) {
            const int c = s_.column(r, a);
            if (!marked[c]) {
              marked[c] = 1;
              touched.push_back(c);
            }
            accumulator[c] += weights_[r];
          }
        }
        for (const int c1 : touched) {
          for (const int c2 : touched) {
            schur(c1, c2) -= accumulator[c1] * accumulator[c2] / diagonal_[h];
          }
        }
        for (const int c : touched) {
          accumulator[c] = 0.0;
          marked[c] = 0;
        }
      }
      root_diagonal_ = diagonal_.array().sqrt();
    }
    if (!schur.allFinite()) {
      valid_ = false;
      return;
    }
    cholesky_.compute(schur);
    valid_ = cholesky_.info() == Eigen::Success &&
      cholesky_.matrixLLT().diagonal().allFinite();
  }

  bool valid() const { return valid_; }

  // out = L w.
  void apply(const VectorXd& w, VectorXd& out) const {
    out.head(s_.n_u) = cholesky_.matrixU().solve(w.head(s_.n_u));
    if (s_.n_d == 0) return;
    VectorXd coupling = VectorXd::Zero(s_.n_d);
    add_du(out.head(s_.n_u), coupling);
    out.tail(s_.n_d) = (w.tail(s_.n_d).array() / root_diagonal_.array() -
      coupling.array() / diagonal_.array()).matrix();
  }

  // out = L' v.
  void apply_transpose(const VectorXd& v, VectorXd& out) const {
    VectorXd upper = v.head(s_.n_u);
    if (s_.n_d > 0) {
      const VectorXd scaled = (v.tail(s_.n_d).array() / diagonal_.array()).matrix();
      subtract_ud(scaled, upper);
      out.tail(s_.n_d) = (v.tail(s_.n_d).array() / root_diagonal_.array()).matrix();
    }
    out.head(s_.n_u) = cholesky_.matrixL().solve(upper);
  }

  // out = L^-1 m.
  void apply_inverse(const VectorXd& m, VectorXd& out) const {
    out.head(s_.n_u) = cholesky_.matrixU() * m.head(s_.n_u);
    if (s_.n_d == 0) return;
    VectorXd coupling = VectorXd::Zero(s_.n_d);
    add_du(m.head(s_.n_u), coupling);
    out.tail(s_.n_d) = (coupling.array() / root_diagonal_.array() +
      root_diagonal_.array() * m.tail(s_.n_d).array()).matrix();
  }

 private:
  // coupling += Q_du x_u, from the records.
  void add_du(const Eigen::Ref<const VectorXd>& x_u, VectorXd& coupling) const {
    for (int r = 0; r < s_.n_records; ++r) {
      double sum = 0.0;
      for (int a = 0; a < s_.width_u; ++a) sum += x_u[s_.column(r, a)];
      coupling[s_.level(r)] += weights_[r] * sum;
    }
  }
  // x_u -= Q_ud z, from the records.
  void subtract_ud(const VectorXd& z, VectorXd& x_u) const {
    for (int r = 0; r < s_.n_records; ++r) {
      const double value = weights_[r] * z[s_.level(r)];
      for (int a = 0; a < s_.width_u; ++a) x_u[s_.column(r, a)] -= value;
    }
  }

  const LatentStructure& s_;
  const VectorXd weights_;
  VectorXd diagonal_;
  VectorXd root_diagonal_;
  LLT<MatrixXd> cholesky_;
  bool valid_ = true;
};

double logistic(const double x) {
  if (x >= 0.0) return 1.0 / (1.0 + std::exp(-x));
  const double z = std::exp(x);
  return z / (1.0 + z);
}

}  // namespace

// Q^-1 g for the Laplace precision of the latent effects at the given record
// weights, by the structured factor (L L' g). Used by the Newton steps of the
// sampler's mode search. Returns NULL when Q is not positive definite at
// these weights.
// [[Rcpp::export]]
SEXP cpp_structured_solve(const IntegerMatrix columns_r, const int n_u,
                          const int n_d, const NumericMatrix prior_u_r,
                          const NumericVector prior_d_r,
                          const NumericVector weights_r,
                          const NumericVector right_hand_side_r) {
  const LatentStructure structure(columns_r, n_u, n_d);
  Map<const MatrixXd> prior_u(prior_u_r.begin(), n_u, n_u);
  Map<const VectorXd> prior_d(prior_d_r.begin(), n_d);
  Map<const VectorXd> weights(weights_r.begin(), weights_r.size());
  const StructuredFactor factor(structure, prior_u, prior_d, weights);
  if (!factor.valid()) return R_NilValue;
  const VectorXd right_hand_side = Rcpp::as<VectorXd>(right_hand_side_r);
  VectorXd half(structure.dimension()), solution(structure.dimension());
  factor.apply_transpose(right_hand_side, half);
  factor.apply(half, solution);
  return Rcpp::wrap(solution);
}

// The structured factor L as a dense matrix, for tests on small problems.
// [[Rcpp::export]]
SEXP cpp_structured_factor_matrix(const IntegerMatrix columns_r, const int n_u,
                                  const int n_d, const NumericMatrix prior_u_r,
                                  const NumericVector prior_d_r,
                                  const NumericVector weights_r) {
  const LatentStructure structure(columns_r, n_u, n_d);
  Map<const MatrixXd> prior_u(prior_u_r.begin(), n_u, n_u);
  Map<const VectorXd> prior_d(prior_d_r.begin(), n_d);
  Map<const VectorXd> weights(weights_r.begin(), weights_r.size());
  const StructuredFactor factor(structure, prior_u, prior_d, weights);
  if (!factor.valid()) return R_NilValue;
  const int dimension = structure.dimension();
  MatrixXd out(dimension, dimension);
  VectorXd unit = VectorXd::Zero(dimension), column(dimension);
  for (int j = 0; j < dimension; ++j) {
    unit[j] = 1.0;
    factor.apply(unit, column);
    out.col(j) = column;
    unit[j] = 0.0;
  }
  return Rcpp::wrap(out);
}

// Metropolis-adjusted Langevin chain for the latent effects of a binomial or
// Poisson model, preconditioned by the Laplace approximation at the mode.
// The chain runs in whitened coordinates w, with latent effects u = m + L w
// and L the structured factor of the Laplace covariance at the record weights
// W. Since P = Q - A' W A and Q = L^-T L^-1, the prior term needs no product
// with P:
//   u' P u = |w + b|^2 - sum_r W_r (eta_r - mu_r)^2,   b = L^-1 m,
//   grad_w = L' A' (s + W (eta - mu)) - w - b,
// with s the scores of the records. Canonical links are evaluated here.
// Other links are evaluated by `link_terms`, an R function of the linear
// predictors returning NULL or a list of the log-likelihood and the scores.
// The chain starts from a draw of the Laplace approximation (or at the mode).
// The step size is tuned during burn-in towards an acceptance of 0.574 by a
// Robbins-Monro update of its logarithm, and the step kept afterwards is the
// average of that logarithm over the second half of burn-in.
// [[Rcpp::export]]
List cpp_latent_mala(const NumericVector y_r, const NumericVector units_r,
                     const NumericVector mu_r, const IntegerMatrix columns_r,
                     const int n_u, const int n_d,
                     const NumericMatrix prior_u_r,
                     const NumericVector prior_d_r,
                     const NumericVector weights_r,
                     const NumericVector mode_r, const std::string family,
                     const SEXP link_terms, const int n_sim, const int burnin,
                     const int thin, double step_size,
                     const double adaptation_scale,
                     const double adaptation_exponent,
                     const bool average_step,
                     const bool start_from_draw) {
  Rcpp::RNGScope scope;
  const LatentStructure structure(columns_r, n_u, n_d);
  const int n = structure.n_records;
  const int dimension = structure.dimension();
  const int n_samples = (n_sim - burnin) / thin;
  Map<const MatrixXd> prior_u(prior_u_r.begin(), n_u, n_u);
  Map<const VectorXd> prior_d(prior_d_r.begin(), n_d);
  Map<const VectorXd> weights(weights_r.begin(), n);
  const VectorXd mode = Rcpp::as<VectorXd>(mode_r);
  const StructuredFactor factor(structure, prior_u, prior_d, weights);
  if (!factor.valid()) {
    Rcpp::stop("The Laplace precision at the mode is not positive definite.");
  }
  const bool poisson = family == "poisson";
  const bool canonical = Rf_isNull(link_terms);
  std::unique_ptr<Rcpp::Function> link_function;
  if (!canonical) link_function.reset(new Rcpp::Function(link_terms));
  NumericVector eta_r(n);

  VectorXd offset(dimension);
  factor.apply_inverse(mode, offset);

  // Evaluates the state at whitened point w into the given vectors.
  VectorXd deviation(dimension), scatter(dimension);
  VectorXd score(n);
  auto evaluate = [&](const VectorXd& w, VectorXd& gradient,
                      VectorXd& latent_out) -> double {
    factor.apply(w, deviation);
    latent_out = mode + deviation;
    double* eta = eta_r.begin();
    for (int r = 0; r < n; ++r) {
      double value = mu_r[r];
      for (int j = 0; j < structure.width; ++j)
        value += latent_out[structure.column(r, j)];
      eta[r] = value;
    }
    double log_likelihood = 0.0;
    if (canonical) {
      for (int r = 0; r < n; ++r) {
        // The log-likelihood is the shared canonical kernel. The score keeps
        // this sampler's logistic, which differs from the site tables' in
        // the last bits for negative eta.
        log_likelihood += canonical_log_likelihood(eta[r], y_r[r],
                                                   units_r[r], poisson);
        if (poisson) {
          const double mean = std::exp(eta[r]);
          if (!std::isfinite(mean)) return R_NegInf;
          score[r] = y_r[r] - units_r[r] * mean;
        } else {
          score[r] = y_r[r] - units_r[r] * logistic(eta[r]);
        }
      }
    } else {
      const SEXP terms = (*link_function)(eta_r);
      if (Rf_isNull(terms)) return R_NegInf;
      const List terms_list(terms);
      log_likelihood = Rcpp::as<double>(terms_list[0]);
      const NumericVector scores = terms_list[1];
      for (int r = 0; r < n; ++r) score[r] = scores[r];
    }
    double data_quadratic = 0.0;
    scatter.setZero();
    for (int r = 0; r < n; ++r) {
      const double centred = eta[r] - mu_r[r];
      data_quadratic += weights[r] * centred * centred;
      const double value = score[r] + weights[r] * centred;
      for (int j = 0; j < structure.width; ++j)
        scatter[structure.column(r, j)] += value;
    }
    factor.apply_transpose(scatter, gradient);
    gradient -= w + offset;
    const double log_density = log_likelihood -
      0.5 * ((w + offset).squaredNorm() - data_quadratic);
    if (!std::isfinite(log_density) || !gradient.allFinite()) return R_NegInf;
    return log_density;
  };

  // The chain starts from a draw of the Laplace approximation, which lies in
  // its typical set. From the mode itself, in d dimensions, a proposal at a
  // step tuned for the stationary spread has log acceptance ratio about
  // -h^4 d / 8, so the chain would accept almost nothing until adaptation had
  // shrunk the step. A draw where the target is not finite is replaced, and
  // after 20 such draws the chain starts at the mode.
  VectorXd current = VectorXd::Zero(dimension);
  VectorXd current_gradient(dimension), current_latent(dimension);
  double current_log_density = R_NegInf;
  if (start_from_draw) {
    for (int attempt = 0; attempt < 20 &&
           !std::isfinite(current_log_density); ++attempt) {
      for (int j = 0; j < dimension; ++j) current[j] = R::rnorm(0.0, 1.0);
      current_log_density = evaluate(current, current_gradient,
                                     current_latent);
    }
  }
  if (!std::isfinite(current_log_density)) {
    current.setZero();
    current_log_density = evaluate(current, current_gradient, current_latent);
  }
  if (!std::isfinite(current_log_density)) {
    Rcpp::stop("The target density or its gradient is not finite at the mode.");
  }
  VectorXd current_mean = current +
    0.5 * step_size * step_size * current_gradient;
  MatrixXd samples(n_samples, dimension);
  NumericVector tuning(n_sim), cumulative_acceptance(n_sim);
  LogicalVector accepted(n_sim);
  int accept_count = 0;
  double log_step_sum = 0.0;
  int log_step_count = 0;
  VectorXd proposal(dimension), proposal_gradient(dimension),
    proposal_latent(dimension), proposal_mean(dimension);

  for (int iteration = 0; iteration < n_sim; ++iteration) {
    if (iteration % 256 == 0) Rcpp::checkUserInterrupt();
    for (int j = 0; j < dimension; ++j)
      proposal[j] = current_mean[j] + step_size * R::rnorm(0.0, 1.0);
    const double proposal_log_density = evaluate(proposal, proposal_gradient,
                                                 proposal_latent);
    bool take = false;
    if (std::isfinite(proposal_log_density)) {
      proposal_mean = proposal +
        0.5 * step_size * step_size * proposal_gradient;
      const double forward = (proposal - current_mean).squaredNorm();
      const double reverse = (current - proposal_mean).squaredNorm();
      const double log_probability = proposal_log_density -
        current_log_density + (forward - reverse) /
        (2.0 * step_size * step_size);
      take = std::isfinite(log_probability) &&
        std::log(R::runif(0.0, 1.0)) < log_probability;
    }
    if (take) {
      ++accept_count;
      accepted[iteration] = true;
      current.swap(proposal);
      current_gradient.swap(proposal_gradient);
      current_latent.swap(proposal_latent);
      current_mean.swap(proposal_mean);
      current_log_density = proposal_log_density;
    }
    if (iteration + 1 > burnin && (iteration + 1 - burnin) % thin == 0) {
      samples.row((iteration + 1 - burnin) / thin - 1) =
        current_latent.transpose();
    }
    cumulative_acceptance[iteration] =
      static_cast<double>(accept_count) / (iteration + 1);
    if (iteration + 1 <= burnin) {
      step_size = std::exp(std::log(step_size) + adaptation_scale *
        std::pow(iteration + 1.0, -adaptation_exponent) *
        (static_cast<double>(take) - 0.574));
      if (average_step && 2 * (iteration + 1) > burnin) {
        log_step_sum += std::log(step_size);
        ++log_step_count;
      }
      if (average_step && iteration + 1 == burnin && log_step_count > 0) {
        step_size = std::exp(log_step_sum / log_step_count);
      }
      current_mean = current +
        0.5 * step_size * step_size * current_gradient;
    }
    tuning[iteration] = step_size;
  }

  return List::create(
    Named("samples") = samples,
    Named("tuning") = tuning,
    Named("acceptance") = accepted,
    Named("cumulative_acceptance") = cumulative_acceptance
  );
}
