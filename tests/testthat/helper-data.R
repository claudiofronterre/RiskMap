set.seed(1)
n <- 10
coords <- cbind(runif(n, 0, 10000), runif(n, 0, 10000))

data <- data.frame(x = coords[,1],
                   z = coords[,2],
                   cov = rnorm(n),
                   den =  sample(5:20, n, replace = TRUE),
                   offset = rnorm(n),
                   i = rep(1:(n/2), each = 2))

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
control_mcmc <- set_control_mcmc(n_sim = 1000, burnin = 200)

gaussian_model <- glgpm(y ~ cov + gp() + re(i),
                        data = gaussian_data,
                        family = "gaussian",
                        messages = FALSE)

gaussian_offset_model <- glgpm(y ~ cov + gp(nugget = TRUE) + offset(offset),
                               data = gaussian_data,
                               family = "gaussian",
                               fix_var_me = 0,
                               messages = FALSE)

eta <- 0.2 + 0.3 * data$cov + S
p <- plogis(eta)
data$y <- rbinom(n, size = data$den, prob = p)
binomial_data <- st_as_sf(data, coords = c("x", "z"), crs = 32637)

binomial_model <- glgpm(y ~ cov + gp() + re(i),
                        data = binomial_data,
                        family = "binomial",
                        den = den,
                        control_mcmc = control_mcmc,
                        messages = FALSE)

lambda <- exp(0.1 + 0.2 * data$cov + S)
data$y <- rpois(n, lambda) + 1
poisson_data <- sf::st_as_sf(data, coords = c("x", "z"), crs = 32637)

poisson_model <- glgpm(y ~ cov + gp() + re(i),
                       data = poisson_data,
                       family = "poisson",
                       den = den,
                       control_mcmc = control_mcmc,
                       return_samples = TRUE,
                       messages = FALSE)

hull <- create_convex_hull(gaussian_data)
grid <- create_grid(hull, 3, propose_utm(hull))
squares <- st_make_grid(hull, n = c(2, 2))
areal <- st_sf(geometry = squares)

