set.seed(1)
n <- 10
coords <- cbind(runif(n, 0, 1), runif(n, 0, 1))

data <- data.frame(lon = coords[,1],
                 lat = coords[,2],
                 cov = rnorm(n),
                 den =  sample(5:20, n, replace = TRUE),
                 i = 1:n)

sigma2 <- 1
phi <- 2
kappa <- 1
Sigma <- sigma2 * matern_cor(dist(coords), phi = phi, kappa = kappa,
                             return_sym_matrix = TRUE)
Sigma_sroot <- t(chol(Sigma))
S <- as.numeric(Sigma_sroot %*% rnorm(n))

data$y <- 1 + 0.5 * data$cov + S + rnorm(n, sd = 0.1)
sf_data <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326)

gaussian_model <- glgpm(y ~ cov + gp(),
                         data = sf_data,
                         family = "gaussian",
                         scale_to_km = FALSE,
                         messages = FALSE)
