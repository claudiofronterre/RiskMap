set.seed(1)
n <- 20L
coords <- cbind(runif(n, 0, 10000), runif(n, 0, 10000))
n_groups <- n / 2L

data <- data.frame(x = coords[,1],
                   z = coords[,2],
                   cov = rnorm(n),
                   denominator =  sample(5:20, n, replace = TRUE),
                   offset = rnorm(n),
                   i = rep(1:n_groups, each = 2),
                   j = rep(1:2, each = n_groups))

sigma2 <- 1
phi <- 5500
kappa <- 1
Sigma <- sigma2 * matern_correlation(dist(coords), phi = phi, kappa = kappa,
                             return_sym_matrix = TRUE)
Sigma_sroot <- t(chol(Sigma))
S <- as.numeric(Sigma_sroot %*% rnorm(n))

sigma2_re <- 1.6
re_values <- rnorm(n_groups, sd = sqrt(sigma2_re))
tau2 <- 0.3
nugget_values <- rnorm(n, sd = sqrt(tau2))

data$y <- 1 + 0.5 * data$cov + S + re_values[data$i] + nugget_values
gaussian_data <- st_as_sf(data, coords = c("x", "z"), crs = 32637)
latlon_data <- st_transform(gaussian_data, crs = 4326)

# reduce iterations to speed up fits
control_mcmc <- set_control_mcmc(n_sim = 1000,
                                 burnin = 200,
                                 seed = 12345)
fixture_control_mcml <- set_control_mcml(
  initial_samples = 500, max_samples = 500, max_stages = 1,
  min_effective_samples = 1, max_pareto_k = 100,
  max_mcse_ratio = 0.99, score_probability = 0.999999
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

eta <- 0.2 + 0.3 * data$cov + S + re_values[data$i]
p <- plogis(eta)
data$y <- rbinom(n, size = data$denominator, prob = p)
binomial_data <- st_as_sf(data, coords = c("x", "z"), crs = 32637)

binomial_model <- glgpm(
  y ~ cov + gp() + re(i), data = binomial_data, family = "binomial",
  denominator = denominator, control_mcmc = control_mcmc,
  control_mcml = fixture_control_mcml, messages = FALSE
)

lambda <- exp(1 + 0.2 * data$cov + S + re_values[data$i])
data$y <- rpois(n, lambda) + 1
poisson_data <- sf::st_as_sf(data, coords = c("x", "z"), crs = 32637)

poisson_model <- glgpm(
  y ~ cov + gp() + re(i), data = poisson_data, family = "poisson",
  denominator = denominator, control_mcmc = control_mcmc,
  control_mcml = fixture_control_mcml, return_samples = TRUE,
  messages = FALSE
)

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

cross_validation <- assess_prediction(
  list(model_a = gaussian_model,
       model_b = gaussian_offset_model),
  method = "cluster",
  fold = 2,
  plot_fold = FALSE,
  messages = FALSE
)

surface <- grid
surface$cov <- rnorm(nrow(grid))
surface$i <- rep(1, nrow(grid))

simulation <- simulate_glgpm(gaussian_model,
                             what = c("data", "surface"),
                             prediction_grid = surface,
                             sample_locations = gaussian_data
                             )
