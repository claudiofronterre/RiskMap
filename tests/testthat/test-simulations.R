
n_sim <- 2
sim_beta <- 0.1
sim_sigma2 <- 1.2
sim_tau2 <- 0.01
sim_phi <- 1.1
sim_sigma2_me <- 0.3


test_that("simulate_glgpm produces errors as expected", {

 expect_error(simulate_glgpm(n_sim = "two",
                        formula = y ~ gp(nugget = 1),
                        data = gaussian_data,
                        family = "gaussian"),
              "'n_sim' must be a single positive integer")

  expect_error(simulate_glgpm(n_sim = -1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian"),
               "'n_sim' must be a single positive integer")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = "not a formula",
                         data = gaussian_data,
                         family = "gaussian"),
               "'formula' must be a 'formula'")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = "gaussian_data",
                         family = "gaussian"),
               "'data' must be of class 'sf'")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "not_gaussian"),
               "'family' must be either 'gaussian'")

  expect_error(simulate_glgpm(n_sim = 1,
                         model_fit = gaussian_model,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian"),
               "if you provide 'model_fit' you should not")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(beta = NULL,
                                         sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         phi = sim_phi,
                                         sigma2_me = sim_sigma2_me)),
               "'beta' is missing")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         phi = sim_phi,
                                         sigma2_me = sim_sigma2_me)),
               "'beta' is missing")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(beta = 1:2,
                                         sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         phi = sim_phi,
                                         sigma2_me = sim_sigma2_me)),
               "the number of values provided for 'beta'")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(beta = sim_beta,
                                         tau2 = sim_tau2,
                                         phi = sim_phi,
                                         sigma2_me = sim_sigma2_me)),
               "'sigma2' is missing")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(beta = sim_beta,
                                         sigma2 = sim_sigma2,
                                         phi = sim_phi,
                                         sigma2_me = sim_sigma2_me)),
               "'tau2' is missing")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(beta = sim_beta,
                                         sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         sigma2_me = sim_sigma2_me)),
               "'phi' is missing")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(beta = sim_beta,
                                         sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         phi = sim_phi
                                         )),
               "'sigma2_me' is missing")

  expect_warning(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(beta = sim_beta,
                                         sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         phi = sim_phi,
                                         sigma2_me = 0,
                                         sigma2_re = 1
                         )),
               "'sigma2_re' will be ignored")

  expect_error(simulate_glgpm(n_sim = 1,
                           formula = y ~ gp(nugget = 1) + re(i),
                           data = gaussian_data,
                           family = "gaussian",
                           sim_pars = list(beta = sim_beta,
                                           sigma2 = sim_sigma2,
                                           tau2 = sim_tau2,
                                           phi = sim_phi,
                                           sigma2_me = 0
                           )),
                 "'sigma2_re' is missing")

  expect_error(simulate_glgpm(n_sim = 1,
                         formula = y ~ gp(nugget = 1) + re(i),
                         data = gaussian_data,
                         family = "gaussian",
                         sim_pars = list(beta = sim_beta,
                                         sigma2 = sim_sigma2,
                                         tau2 = sim_tau2,
                                         phi = sim_phi,
                                         sigma2_me = sim_sigma2_me,
                                         sigma2_re = 1:2
                         )),
               "the values passed to 'sigma2_re' in 'sim_pars'")
})

test_that("simulate_glgpm produces expected output for gaussian models", {

  set.seed(1)
  result <- simulate_glgpm(n_sim = n_sim,
                      formula = y ~ gp(nugget = 1),
                      data = gaussian_data,
                      family = "gaussian",
                      sim_pars = list(beta = sim_beta,
                                      sigma2 = sim_sigma2,
                                      tau2 = sim_tau2,
                                      phi = sim_phi,
                                      sigma2_me = 0),
                      messages = FALSE)

  expect_setequal(names(result), c("data_sim", "S_sim", "lin_pred_sim", "beta", "sigma2", "tau2", "phi", "sigma2_me"))
  expect_s3_class(result$data_sim, "sf")
  expect_equal(result$data_sim$geometry, gaussian_data$geometry)
  expect_equal(nrow(result$data_sim), nrow(gaussian_data))
  expect_equal(ncol(result$data_sim), ncol(gaussian_data) + n_sim)

  expect_equal(nrow(result$S_sim), n_sim)
  expect_equal(ncol(result$S_sim), nrow(gaussian_data))

  expect_equal(nrow(result$lin_pred_sim), n_sim)
  expect_equal(ncol(result$lin_pred_sim), nrow(gaussian_data))

  expect_equal(result$beta, sim_beta)
  expect_equal(result$sigma2, sim_sigma2)
  expect_equal(result$tau2, sim_tau2)
  expect_equal(result$phi, sim_phi)
  expect_equal(result$sigma2_me, 0)

})

test_that("simulate_glgpm produces expected output from a a gaussian model", {

  model <- glgpm(y ~ cov + gp(),
                 data = gaussian_data,
                 family = "gaussian",
                 scale_to_km = FALSE,
                 messages = FALSE)

  result <- simulate_glgpm(n_sim = n_sim,
                      model_fit = model)

  expect_setequal(names(result), c("data_sim", "S_sim", "lin_pred_sim", "beta", "sigma2", "tau2", "phi", "sigma2_me"))
  expect_s3_class(result$data_sim, "sf")
  expect_equal(result$data_sim$geometry, gaussian_data$geometry)
  expect_equal(nrow(result$data_sim), nrow(gaussian_data))
  expect_equal(ncol(result$data_sim), ncol(gaussian_data) + n_sim)

  expect_equal(nrow(result$S_sim), n_sim)
  expect_equal(ncol(result$S_sim), nrow(gaussian_data))

  expect_equal(nrow(result$lin_pred_sim), n_sim)
  expect_equal(ncol(result$lin_pred_sim), nrow(gaussian_data))

})
