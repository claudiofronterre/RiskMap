// =============================================================================
// dsgm_spatial.stan
//
// Stan sampler for the doubly stochastic geostatistical model for joint
// prevalence-intensity modelling of soil-transmitted helminths (STH).
// Used to draw MCMC samples of the spatial process S(x) at fixed parameter
// values theta_0, which are then passed to the TMB MCML engine (dsgm_tmb.cpp).
//
// Model:
//   W_j(x_i) | S(x_i) ~ NegBin(mu(x_i), k)
//     log mu(x_i) = eta_fixed[i] + S[ID_coords[i]]
//     mu(x_i)     = exp(eta_fixed[i] + S[ID_coords[i]]) * mda_impact[i]
//
//   Y_ij | W_j(x_i) ~ Poisson(rho * W_j(x_i))
//
// Marginal likelihoods used:
//   Zero observations:
//     P(Y = 0) = 1 - P(Y > 0)  where
//     P(Y > 0) = 1 - [k / (k + mu*(1-exp(-rho)))]^k
//
//   Positive observations (Y > 0):
//     P(Y > 0) as above, plus
//     Gamma approximation for Y | Y > 0:
//       Y - 1 ~ Gamma(kappa_C, theta_C)
//       mu_C    = rho * mu / P(Y>0)
//       sigma2_C = rho*mu*(1+rho)/P(Y>0) + rho^2*mu^2/P(Y>0) * (1/k + 1 - 1/P(Y>0))
//       kappa_C  = (mu_C - 1)^2 / sigma2_C
//       theta_C  = sigma2_C / (mu_C - 1)
//
// Spatial process:
//   S = sqrt(sigma2) * L * S_raw,   S_raw ~ N(0, I)
//   R[i,j] = exp(-D_mat[i,j] / phi),  L = cholesky(R)
//   L is computed once in transformed_data since phi is fixed.
//
// All model parameters (k, rho, sigma2, phi) are FIXED at their theta_0
// values and passed as data.  Only S_raw is sampled.
//
// NOTE: parameter names follow R convention (k, rho) not C++ convention
// (omega, alpha); both Stan and TMB use their own internal names.
// =============================================================================

data {

  int<lower=1> n;        // Total number of observations
  int<lower=1> n_loc;    // Number of unique spatial locations
  int<lower=1> n_pos;    // Number of egg-positive observations
  int<lower=1> p;        // Number of fixed-effect covariates (unused here; eta pre-computed)

  // Binary response
  array[n] int<lower=0, upper=1> y;  // 1 = egg-positive, 0 = negative

  // EPG for positive observations
  vector<lower=1>[n_pos] C_pos;      // Egg counts for positives (>= 1)

  // Indices of positive observations (1-indexed, matching y)
  array[n_pos] int<lower=1, upper=n> pos_idx;

  // Location mapping (1-indexed)
  array[n] int<lower=1, upper=n_loc> ID_coords;

  // Pairwise distance matrix between unique locations
  matrix[n_loc, n_loc] D_mat;

  // Pre-computed fixed-effects linear predictor: D * beta + offset
  vector[n] eta_fixed;

  // Pre-computed MDA decay factor per observation (1.0 = no MDA)
  vector<lower=0>[n] mda_impact;

  // Fixed parameters at theta_0
  real<lower=0> k;       // NB aggregation parameter
  real<lower=0> rho;     // Per-worm egg detection rate
  real<lower=0> sigma2;  // GP variance
  real<lower=0> phi;     // GP range

}

transformed data {

  // Correlation matrix and Cholesky factor, computed once at fixed phi.
  matrix[n_loc, n_loc] R;
  matrix[n_loc, n_loc] L;
  real c_rho;            // 1 - exp(-rho), used repeatedly

  c_rho = 1.0 - exp(-rho);

  for (i in 1:n_loc) {
    R[i, i] = 1.0;
    for (j in (i + 1):n_loc) {
      R[i, j] = exp(-D_mat[i, j] / phi);
      R[j, i] = R[i, j];
    }
  }

  L = cholesky_decompose(R);

}

parameters {

  // Standardised spatial random effects.
  // S = sqrt(sigma2) * L * S_raw ~ N(0, sigma2 * R)
  vector[n_loc] S_raw;

}

transformed parameters {

  vector[n_loc] S;
  vector[n]     mu_W;    // Expected worm burden per observation
  vector[n]     pr_pos;  // P(Y > 0) per observation

  S = sqrt(sigma2) * (L * S_raw);

  for (i in 1:n) {
    mu_W[i]   = exp(eta_fixed[i] + S[ID_coords[i]]) * mda_impact[i];
    real ratio = k / (k + mu_W[i] * c_rho);
    pr_pos[i]  = fmax(fmin(1.0 - pow(ratio, k), 1.0 - 1e-10), 1e-10);
  }

}

model {

  // Prior on standardised spatial process
  S_raw ~ std_normal();

  // ----- Prevalence likelihood -----
  for (i in 1:n) {
    if (y[i] == 1) {
      target += log(pr_pos[i]);
    } else {
      target += log1m(pr_pos[i]);
    }
  }

  // ----- Intensity likelihood (positives only, Gamma approximation) -----
  for (idx in 1:n_pos) {
    int i = pos_idx[idx];

    // Conditional moments of Y | Y > 0
    real mu_C    = (rho * mu_W[i]) / pr_pos[i];
    real sigma2_C = (rho * mu_W[i] * (1.0 + rho)) / pr_pos[i]
                  + (rho^2 * mu_W[i]^2 / pr_pos[i])
                    * (1.0/k + 1.0 - 1.0/pr_pos[i]);
    sigma2_C = fmax(sigma2_C, 1e-10);

    // Gamma parameters for Y - 1
    real mu_C1  = fmax(mu_C - 1.0, 1e-6);
    real kappa_C = fmax(mu_C1^2 / sigma2_C, 1e-10);
    real theta_C = fmax(sigma2_C / mu_C1, 1e-10);

    // log p(C_pos - 1 | Gamma(kappa_C, theta_C))
    target += gamma_lpdf(fmax(C_pos[idx] - 1.0, 0.0) | kappa_C, 1.0 / theta_C);
  }

}
