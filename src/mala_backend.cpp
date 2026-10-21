// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppEigen.h>
#include <limits>

using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Upper;
using Eigen::VectorXd;
using Rcpp::IntegerVector;
using Rcpp::IntegerMatrix;
using Rcpp::List;
using Rcpp::LogicalVector;
using Rcpp::Named;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

struct State {
  double log_density;
  VectorXd gradient;
  VectorXd latent;
  bool valid;
};

double logistic(const double x) {
  if (x >= 0.0) return 1.0 / (1.0 + std::exp(-x));
  const double z = std::exp(x);
  return z / (1.0 + z);
}

double softplus(const double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) :
    std::log1p(std::exp(x));
}

}  // namespace

// [[Rcpp::export]]
List cpp_canonical_likelihood(const NumericMatrix eta_r,
                                       const NumericVector y_r,
                                       const NumericVector units_r,
                                       const std::string family) {
  const int n_samples = eta_r.nrow();
  const int n = eta_r.ncol();
  if (y_r.size() != n || units_r.size() != n)
    Rcpp::stop("Observation vectors have inconsistent lengths.");
  NumericVector log_likelihood(n_samples);
  NumericMatrix score(n_samples, n), curvature(n_samples, n);
  for (int observation = 0; observation < n; ++observation) {
    for (int sample = 0; sample < n_samples; ++sample) {
      const double eta = eta_r(sample, observation);
      if (family == "binomial") {
        const double probability = logistic(eta);
        log_likelihood[sample] += y_r[observation] * eta -
          units_r[observation] * softplus(eta);
        score(sample, observation) = y_r[observation] -
          units_r[observation] * probability;
        curvature(sample, observation) = units_r[observation] * probability *
          (1.0 - probability);
      } else if (family == "poisson") {
        const double mean = std::exp(eta);
        log_likelihood[sample] += y_r[observation] * eta -
          units_r[observation] * mean;
        score(sample, observation) = y_r[observation] -
          units_r[observation] * mean;
        curvature(sample, observation) = units_r[observation] * mean;
      } else {
        Rcpp::stop("Unsupported family.");
      }
    }
  }
  return List::create(Named("log_likelihood") = log_likelihood,
                      Named("score") = score,
                      Named("curvature") = curvature);
}

// Log-likelihood only, for every row of eta (draws x observations). The
// family branch is resolved once, no score or curvature matrices are
// allocated, and a Poisson mean that overflows gives -Inf for that draw, as in
// the R implementation.
// [[Rcpp::export]]
NumericVector cpp_canonical_log_likelihood(const NumericMatrix eta_r,
                                           const NumericVector y_r,
                                           const NumericVector units_r,
                                           const bool poisson) {
  const int n_samples = eta_r.nrow();
  const int n = eta_r.ncol();
  if (y_r.size() != n || units_r.size() != n)
    Rcpp::stop("Observation vectors have inconsistent lengths.");
  NumericVector log_likelihood(n_samples);
  const double* eta = eta_r.begin();
  double* total = log_likelihood.begin();
  for (int observation = 0; observation < n; ++observation) {
    const double y = y_r[observation];
    const double units = units_r[observation];
    const double* column = eta + static_cast<R_xlen_t>(observation) * n_samples;
    if (poisson) {
      for (int sample = 0; sample < n_samples; ++sample) {
        total[sample] += y * column[sample] - units * std::exp(column[sample]);
      }
    } else {
      for (int sample = 0; sample < n_samples; ++sample) {
        total[sample] += y * column[sample] - units * softplus(column[sample]);
      }
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

// Canonical binomial/Poisson MALA. Random effects are represented as additional
// independent latent blocks and are mapped to observations by their level IDs.
// [[Rcpp::export]]
List cpp_canonical_mala(
    const NumericVector y_r,
    const NumericVector units_r,
    const NumericVector mu_r,
    const NumericMatrix proposal_lower_r,
    const NumericMatrix proposal_upper_r,
    const NumericMatrix prior_lower_r,
    const NumericMatrix prior_upper_r,
    const NumericVector transformed_prior_mean_r,
    const NumericVector proposal_mean_r,
    const IntegerVector location_id_r,
    const int n_locations,
    const IntegerMatrix random_effect_id_r,
    const std::string family,
    const int n_sim,
    const int burnin,
    const int thin,
    double step_size,
    const double adaptation_scale,
    const double adaptation_exponent) {
  Rcpp::RNGScope scope;
  const int n = y_r.size();
  const int dimension = proposal_mean_r.size();
  const int n_samples = (n_sim - burnin) / thin;
  if (family != "binomial" && family != "poisson")
    Rcpp::stop("Unsupported family.");
  if (units_r.size() != n || mu_r.size() != n || location_id_r.size() != n)
    Rcpp::stop("Observation vectors have inconsistent lengths.");
  if (random_effect_id_r.nrow() != n)
    Rcpp::stop("Random-effect indices must have one row per observation.");
  const int number_random_effects = random_effect_id_r.ncol();
  IntegerVector random_effect_offsets(number_random_effects);
  int expected_dimension = n_locations;
  for (int effect = 0; effect < number_random_effects; ++effect) {
    random_effect_offsets[effect] = expected_dimension;
    int maximum_level = 0;
    for (int observation = 0; observation < n; ++observation) {
      maximum_level = std::max(maximum_level,
        random_effect_id_r(observation, effect));
    }
    if (maximum_level < 1) Rcpp::stop("Invalid random-effect index.");
    expected_dimension += maximum_level;
  }
  if (expected_dimension != dimension)
    Rcpp::stop("Latent dimension does not match random-effect indices.");

  Map<const VectorXd> y(y_r.begin(), n);
  Map<const VectorXd> units(units_r.begin(), n);
  Map<const VectorXd> mu(mu_r.begin(), n);
  Map<const MatrixXd> proposal_lower(proposal_lower_r.begin(), dimension,
                                      dimension);
  Map<const MatrixXd> proposal_upper(proposal_upper_r.begin(), dimension,
                                      dimension);
  Map<const MatrixXd> prior_lower(prior_lower_r.begin(), dimension, dimension);
  Map<const MatrixXd> prior_upper(prior_upper_r.begin(), dimension, dimension);
  Map<const VectorXd> transformed_prior_mean(transformed_prior_mean_r.begin(),
                                              dimension);
  Map<const VectorXd> proposal_mean(proposal_mean_r.begin(), dimension);

  // Validate every index once, so the sampling loop needs no bounds checks.
  for (int observation = 0; observation < n; ++observation) {
    const int location = location_id_r[observation] - 1;
    if (location < 0 || location >= n_locations)
      Rcpp::stop("Location index outside latent dimension.");
    for (int effect = 0; effect < number_random_effects; ++effect) {
      if (random_effect_id_r(observation, effect) < 1)
        Rcpp::stop("Invalid random-effect index.");
    }
  }
  const bool poisson = family == "poisson";

  auto evaluate = [&](const VectorXd& w) {
    State state;
    state.latent = proposal_lower * w + proposal_mean;
    VectorXd gradient_latent = VectorXd::Zero(dimension);
    double log_likelihood = 0.0;
    for (int observation = 0; observation < n; ++observation) {
      const int location = location_id_r[observation] - 1;
      double eta = mu[observation] + state.latent[location];
      for (int effect = 0; effect < number_random_effects; ++effect) {
        eta += state.latent[random_effect_offsets[effect] +
          random_effect_id_r(observation, effect) - 1];
      }
      double residual;
      if (poisson) {
        const double mean = std::exp(eta);
        if (!std::isfinite(mean)) {
          state.valid = false;
          return state;
        }
        log_likelihood += y[observation] * eta - units[observation] * mean;
        residual = y[observation] - units[observation] * mean;
      } else {
        log_likelihood += y[observation] * eta -
          units[observation] * softplus(eta);
        residual = y[observation] - units[observation] * logistic(eta);
      }
      gradient_latent[location] += residual;
      for (int effect = 0; effect < number_random_effects; ++effect) {
        gradient_latent[random_effect_offsets[effect] +
          random_effect_id_r(observation, effect) - 1] += residual;
      }
    }
    VectorXd difference = w - transformed_prior_mean;
    VectorXd prior_solution = prior_lower.triangularView<Lower>().solve(
      difference);
    prior_solution = prior_upper.triangularView<Upper>().solve(prior_solution);
    state.log_density = -0.5 * difference.dot(prior_solution) + log_likelihood;
    state.gradient = -prior_solution + proposal_upper * gradient_latent;
    state.valid = std::isfinite(state.log_density) &&
      state.gradient.allFinite();
    return state;
  };

  VectorXd current = VectorXd::Zero(dimension);
  State current_state = evaluate(current);
  if (!current_state.valid) Rcpp::stop("Invalid initial state.");
  const double initial_log_density = current_state.log_density;
  const VectorXd initial_gradient = current_state.gradient;
  VectorXd current_mean = current +
    0.5 * step_size * step_size * current_state.gradient;
  MatrixXd samples(n_samples, dimension);
  NumericVector tuning(n_sim), cumulative_acceptance(n_sim);
  LogicalVector accepted(n_sim);
  int accept_count = 0;

  for (int iteration = 0; iteration < n_sim; ++iteration) {
    VectorXd proposal(dimension);
    for (int j = 0; j < dimension; ++j)
      proposal[j] = current_mean[j] + step_size * R::rnorm(0.0, 1.0);
    State proposal_state = evaluate(proposal);
    bool take = false;
    VectorXd proposal_mean_vector(dimension);
    if (proposal_state.valid) {
      proposal_mean_vector = proposal +
        0.5 * step_size * step_size * proposal_state.gradient;
      const double forward = (proposal - current_mean).squaredNorm();
      const double reverse = (current - proposal_mean_vector).squaredNorm();
      const double log_probability = proposal_state.log_density -
        current_state.log_density + (forward - reverse) /
        (2.0 * step_size * step_size);
      take = std::isfinite(log_probability) &&
        std::log(R::runif(0.0, 1.0)) < log_probability;
    }
    if (take) {
      ++accept_count;
      accepted[iteration] = true;
      current = proposal;
      current_state = proposal_state;
      current_mean = proposal_mean_vector;
    }
    if (iteration + 1 > burnin && (iteration + 1 - burnin) % thin == 0) {
      samples.row((iteration + 1 - burnin) / thin - 1) =
        current_state.latent.transpose();
    }
    cumulative_acceptance[iteration] =
      static_cast<double>(accept_count) / (iteration + 1);
    if (iteration + 1 <= burnin) {
      step_size = std::exp(std::log(step_size) + adaptation_scale *
        std::pow(iteration + 1.0, -adaptation_exponent) *
        (static_cast<double>(take) - 0.574));
      current_mean = current +
        0.5 * step_size * step_size * current_state.gradient;
    }
    tuning[iteration] = step_size;
  }

  return List::create(
    Named("samples") = samples,
    Named("tuning") = tuning,
    Named("acceptance") = accepted,
    Named("cumulative_acceptance") = cumulative_acceptance,
    Named("initial_log_density") = initial_log_density,
    Named("initial_gradient") = initial_gradient
  );
}
