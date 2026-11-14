set.seed(1)
n <- 10L
coords <- cbind(runif(n, 0, 10000), runif(n, 0, 10000))

data <- data.frame(x = coords[,1],
                   z = coords[,2],
                   cov = rnorm(n),
                   denominator =  sample(5:20, n, replace = TRUE),
                   offset = rnorm(n),
                   i = rep(1:(n/2), each = 2),
                   j = rep(1:2, each = n/2))

sigma2 <- 1
phi <- 2
kappa <- 1
Sigma <- sigma2 * matern_correlation(dist(coords), phi = phi, kappa = kappa,
                             return_sym_matrix = TRUE)
Sigma_sroot <- t(chol(Sigma))
S <- as.numeric(Sigma_sroot %*% rnorm(n))

data$y <- 1 + 0.5 * data$cov + S + rnorm(n, sd = 0.1)
gaussian_data <- st_as_sf(data, coords = c("x", "z"), crs = 32637)
latlon_data <- st_transform(gaussian_data, crs = 4326)

# reduce iterations to speed up fits
control_mcmc <- set_control_mcmc(n_sim = 1000,
                                 burnin = 200)
fixture_seed <- 12345
fixture_control_mcml <- set_control_mcml(
  initial_samples = 500, validation_samples = 500, max_samples = 500,
  max_stages = 1
)

gaussian_model <- glgpm(y ~ cov + gp() + re(i),
                        data = gaussian_data,
                        family = "gaussian",
                        messages = FALSE)

gaussian_offset_model <- glgpm(y ~ cov + gp(nugget = TRUE) + offset(offset),
                               data = gaussian_data,
                               family = "gaussian",
                               fix_var_me = 0,
                               messages = FALSE)

gaussian_intercept_model <- glgpm(y ~ gp(),
                                  data = gaussian_data,
                                  family = "gaussian",
                                  messages = FALSE)

eta <- 0.2 + 0.3 * data$cov + S
p <- plogis(eta)
data$y <- rbinom(n, size = data$denominator, prob = p)
binomial_data <- st_as_sf(data, coords = c("x", "z"), crs = 32637)

binomial_model <- suppressWarnings(glgpm(
  y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
  denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
  control_mcml = fixture_control_mcml, messages = FALSE
))

lambda <- exp(0.1 + 0.2 * data$cov + S)
data$y <- rpois(n, lambda) + 1
poisson_data <- sf::st_as_sf(data, coords = c("x", "z"), crs = 32637)

poisson_model <- suppressWarnings(glgpm(
  y ~ cov + gp() + re(i), data = poisson_data, family = "poisson",
  denominator = denominator, control_mcmc = control_mcmc, seed = fixture_seed,
  control_mcml = fixture_control_mcml, return_samples = TRUE,
  messages = FALSE
))

hull <- create_convex_hull(gaussian_data)
grid <- create_grid(hull, 3)
squares <- st_make_grid(hull, n = c(2, 2))
areal <- st_sf(geometry = squares)

prediction_data <- grid
prediction_data$cov <- seq_len(nrow(prediction_data)) / 10
prediction_data$offset <- seq_len(nrow(prediction_data)) / 100
prediction_data$i <- rep(gaussian_model$re[[1]],
                         length.out = nrow(prediction_data))
prediction_data$population <- seq_len(nrow(prediction_data))

# model_a has a random effect, which assess_prediction() warns is not yet
# included in the predictive distributions
cross_validation <- suppressWarnings(assess_prediction(
  list(model_a = gaussian_model,
       model_b = gaussian_offset_model),
  method = "cluster",
  fold = 2,
  plot_fold = FALSE,
  messages = FALSE
))

surface <- grid
surface$cov <- rnorm(nrow(grid))
surface$i <- rep(1, nrow(grid))

simulation <- simulate_glgpm(gaussian_model,
                             what = c("data", "surface"),
                             prediction_grid = surface,
                             sample_locations = gaussian_data
                             )
