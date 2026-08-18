
n_sim <- 2
sim_beta <- 0.1
sim_sigma2 <- 1.2
sim_tau2 <- 0.01
sim_phi <- 1.1
sim_sigma2_me <- 0.3


test_that("glgpm_sim produces errors as expected", {

 expect_error(glgpm_sim(n_sim = "two",
                        formula = y ~ gp(nugget = 1),
                        data = sf_data,
                        family = "gaussian",
                        sim_pars = list(beta = sim_beta,
                                        sigma2 = sim_sigma2,
                                        tau2 = sim_tau2,
                                        phi = sim_phi,
                                        sigma2_me = sim_sigma2_me),
                        messages = FALSE),
              "'n_sim' must be a single positive integer")

  expect_error(glgpm_sim(n_sim = -1,
                         formula = y ~ gp(nugget = 1),
                         data = sf_data,
                         family = "gaussian",
                         sim_pars = list(beta = sim_beta,
                                         sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         phi = sim_phi,
                                         sigma2_me = sim_sigma2_me),
                         messages = FALSE),
               "'n_sim' must be a single positive integer")

  expect_error(glgpm_sim(n_sim = -1,
                         formula = y ~ gp(nugget = 1),
                         data = sf_data,
                         family = "gaussian",
                         sim_pars = list(beta = sim_beta,
                                         sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         phi = sim_phi,
                                         sigma2_me = sim_sigma2_me),
                         messages = FALSE),
               "'n_sim' must be a single positive integer")

})

test_that("glgpm_sim produces expected output for gaussian models", {

set.seed(1)
result <- glgpm_sim(n_sim = n_sim,
                    formula = y ~ gp(nugget = 1),
                    data = sf_data,
                    family = "gaussian",
                    sim_pars = list(beta = sim_beta,
                                    sigma2 = sim_sigma2,
                                    tau2 = sim_tau2,
                                    phi = sim_phi,
                                    sigma2_me = sim_sigma2_me),
                    messages = FALSE)

expect_setequal(names(result), c("data_sim", "S_sim", "lin_pred_sim", "beta", "sigma2", "tau2", "phi", "sigma2_me"))
expect_s3_class(result$data_sim, "sf")
expect_equal(result$data_sim$geometry, sf_data$geometry)
expect_equal(nrow(result$data_sim), nrow(sf_data))
expect_equal(ncol(result$data_sim), ncol(sf_data) + n_sim)

expect_equal(nrow(result$S_sim), n_sim)
expect_equal(ncol(result$S_sim), nrow(sf_data))

expect_equal(nrow(result$lin_pred_sim), n_sim)
expect_equal(ncol(result$lin_pred_sim), nrow(sf_data))

expect_equal(result$beta, sim_beta)
expect_equal(result$sigma2, sim_sigma2)
expect_equal(result$tau2, sim_tau2)
expect_equal(result$phi, sim_phi)
expect_equal(result$sigma2_me, sim_sigma2_me)

})

test_that("glgpm_sim produces expected output from a a gaussian model", {

  model <- glgpm(y ~ cov + gp(sf),
                 data = sf_data,
                 family = "gaussian",
                 scale_to_km = FALSE,
                 messages = FALSE)

  result <- glgpm_sim(n_sim = n_sim,
                      model_fit = model)

  expect_setequal(names(result), c("data_sim", "S_sim", "lin_pred_sim", "beta", "sigma2", "tau2", "phi", "sigma2_me"))
  expect_s3_class(result$data_sim, "sf")
  expect_equal(result$data_sim$geometry, sf_data$geometry)
  expect_equal(nrow(result$data_sim), nrow(sf_data))
  expect_equal(ncol(result$data_sim), ncol(sf_data) + n_sim)

  expect_equal(nrow(result$S_sim), n_sim)
  expect_equal(ncol(result$S_sim), nrow(sf_data))

  expect_equal(nrow(result$lin_pred_sim), n_sim)
  expect_equal(ncol(result$lin_pred_sim), nrow(sf_data))

})
