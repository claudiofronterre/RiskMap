## Temporary regression check for #92: confirms that the name-based
## coef.RiskMap()/summary.RiskMap() (see R/auxiliary.R) produce output
## numerically identical to the pre-#92 implementation, which looked up
## parameters by recomputed index arithmetic instead of by name. The
## pre-#92 versions are reproduced verbatim below (from commit 74799f0)
## purely so CI can diff their output against the current implementation
## across every fixture model. Delete this file once the team is satisfied
## the refactor is behaviour-preserving (tracked in #92).

old_coef_riskmap <- function(object) {

  n_re <- length(object$re)
  if (n_re > 0) re_names <- names(object$re)

  p        <- ncol(as.matrix(object$D))
  ind_beta <- 1:p

  if (p == 1) {
    object$D <- as.matrix(object$D)
    names(object$estimate)[ind_beta] <- "Intercept"
  } else {
    names(object$estimate)[ind_beta] <- colnames(object$D)
  }
  ind_sigma2 <- p + 1
  names(object$estimate)[ind_sigma2] <- "sigma2"
  ind_phi <- p + 2
  names(object$estimate)[ind_phi] <- "phi"

  if (isTRUE(object$fix_tau2)) {
    ind_tau2 <- p + 3
    names(object$estimate)[ind_tau2] <- "tau2"
    object$estimate[ind_tau2] <- object$estimate[ind_tau2] + object$estimate[ind_sigma2]
    if (object$family == "gaussian") {
      if (is.null(object$fix_var_me)) {
        ind_sigma2_me <- p + 4
        if (n_re > 0) ind_sigma2_re <- (p + 5):(p + 4 + n_re)
      } else {
        ind_sigma2_me <- NULL
        if (n_re > 0) ind_sigma2_re <- (p + 4):(p + 3 + n_re)
      }
    } else {
      ind_sigma2_me <- NULL
      if (n_re > 0) ind_sigma2_re <- (p + 4):(p + 3 + n_re)
    }
  } else {
    ind_tau2 <- NULL
    if (object$family == "gaussian") {
      if (is.null(object$fix_var_me)) {
        ind_sigma2_me <- p + 3
        names(object$estimate)[ind_sigma2_me] <- "sigma2_me"
        if (n_re > 0) ind_sigma2_re <- (p + 4):(p + 3 + n_re)
      } else {
        ind_sigma2_me <- NULL
        if (n_re > 0) ind_sigma2_re <- (p + 3):(p + 2 + n_re)
      }
    } else {
      if (n_re > 0) ind_sigma2_re <- (p + 3):(p + 2 + n_re)
    }
  }

  ind_sp <- c(ind_sigma2, ind_phi, ind_tau2)
  object$estimate[ind_sp] <- exp(object$estimate[ind_sp])

  if (n_re > 0) {
    for (i in seq_len(n_re))
      names(object$estimate)[ind_sigma2_re[i]] <-
        paste0(re_names[i], "_sigma2_re")
  }

  res        <- list()
  res$beta   <- object$estimate[ind_beta]
  res$sigma2 <- as.numeric(object$estimate[ind_sigma2])
  res$phi    <- as.numeric(object$estimate[ind_phi])
  if (object$family == "gaussian" && !is.null(ind_sigma2_me))
    res$sigma2_me <- exp(as.numeric(object$estimate[ind_sigma2_me]))
  if (!is.null(ind_tau2))
    res$tau2 <- object$estimate[ind_tau2]
  if (n_re > 0)
    res$sigma2_re <- exp(as.numeric(object$estimate[ind_sigma2_re]))

  return(res)
}

old_summary_riskmap <- function(object, conf_level = 0.95) {

  alpha  <- 1 - conf_level
  z_crit <- qnorm(1 - alpha / 2)
  res    <- list()

  link_name <- NULL
  inv_expr  <- NULL
  if (!is.null(object$linkf) && is.list(object$linkf) &&
      is.function(object$linkf$inv)) {
    link_name <- object$linkf$name %||% "custom"
    inv_expr  <- tryCatch(
      paste0("Inverse link function = ",
             paste(deparse(body(object$linkf$inv), width.cutoff = 500L),
                   collapse = " ")),
      error = function(e) "Inverse link function = <user-supplied function>"
    )
  } else {
    if (identical(object$family, "poisson")) {
      link_name <- "canonical (log)"
      inv_expr  <- "Inverse link function = exp(x)"
    } else if (identical(object$family, "binomial")) {
      link_name <- "canonical (logit)"
      inv_expr  <- "Inverse link function = 1 / (1 + exp(-x))"
    } else if (identical(object$family, "gaussian")) {
      link_name <- "identity"
      inv_expr  <- "Inverse link function = x"
    }
  }

  n_re <- length(object$re)
  if (n_re > 0) re_names <- names(object$re)

  p        <- ncol(object$D)
  ind_beta <- seq_len(p)

  names(object$estimate)[ind_beta] <- colnames(object$D)
  ind_sigma2 <- p + 1; names(object$estimate)[ind_sigma2] <- "Spatial process var."
  ind_phi    <- p + 2; names(object$estimate)[ind_phi]    <- "Spatial corr. scale"

  if (isTRUE(object$fix_tau2)) {
    ind_tau2 <- p + 3
    names(object$estimate)[ind_tau2] <- "Variance of the nugget"
    object$estimate[ind_tau2] <- object$estimate[ind_tau2] + object$estimate[ind_sigma2]
    if (object$family == "gaussian") {
      ind_sigma2_me <- if (is.null(object$fix_var_me)) p + 4 else NULL
      if (n_re > 0) ind_sigma2_re <- (p + 5):(p + 4 + n_re)
    } else {
      ind_sigma2_re <- (p + 4):(p + 3 + n_re)
    }
  } else {
    ind_tau2 <- NULL
    if (object$family == "gaussian") {
      if (is.null(object$fix_var_me)) {
        ind_sigma2_me <- p + 3
        names(object$estimate)[ind_sigma2_me] <- "Measurement error var."
      } else {
        ind_sigma2_me <- NULL
      }
      if (n_re > 0) ind_sigma2_re <- (p + 4):(p + 3 + n_re)
    } else {
      ind_sigma2_re <- (p + 3):(p + 2 + n_re)
    }
  }

  ind_sp <- c(ind_sigma2, ind_phi, ind_tau2)

  n_p <- length(object$estimate)
  object$estimate[-c(ind_beta)] <- exp(object$estimate[-c(ind_beta)])

  if (n_re > 0)
    for (i in seq_len(n_re))
      names(object$estimate)[ind_sigma2_re[i]] <-
    paste0(re_names[i], " (random eff. var.)")

  J <- diag(n_p)
  if (length(ind_tau2) > 0) J[ind_tau2, ind_sigma2] <- 1
  H_new          <- t(J) %*% solve(-object$covariance) %*% J
  covariance_new <- solve(-H_new)
  se_par         <- sqrt(diag(covariance_new))

  zval <- object$estimate[ind_beta] / se_par[ind_beta]
  res$reg_coef <- cbind(
    Estimate      = object$estimate[ind_beta],
    "Lower limit" = object$estimate[ind_beta] - se_par[ind_beta] * z_crit,
    "Upper limit" = object$estimate[ind_beta] + se_par[ind_beta] * z_crit,
    StdErr        = se_par[ind_beta],
    z.value       = zval,
    p.value       = 2 * pnorm(-abs(zval))
  )

  if (object$family == "gaussian") {
    if (is.null(object$fix_var_me)) {
      res$me <- cbind(
        Estimate      = object$estimate[ind_sigma2_me],
        "Lower limit" = exp(log(object$estimate[ind_sigma2_me]) -
                              z_crit * se_par[ind_sigma2_me]),
        "Upper limit" = exp(log(object$estimate[ind_sigma2_me]) +
                              z_crit * se_par[ind_sigma2_me])
      )
    } else {
      res$me <- object$fix_var_me
    }
  }

  res$sp <- cbind(
    Estimate      = object$estimate[ind_sp],
    "Lower limit" = exp(log(object$estimate[ind_sp]) - z_crit * se_par[ind_sp]),
    "Upper limit" = exp(log(object$estimate[ind_sp]) + z_crit * se_par[ind_sp])
  )
  if (!is.null(object$fix_tau2)) res$tau2 <- object$fix_tau2

  if (n_re > 0)
    res$ranef <- cbind(
      Estimate      = object$estimate[ind_sigma2_re],
      "Lower limit" = exp(log(object$estimate[ind_sigma2_re]) -
                            z_crit * se_par[ind_sigma2_re]),
      "Upper limit" = exp(log(object$estimate[ind_sigma2_re]) +
                            z_crit * se_par[ind_sigma2_re])
    )

  res$conf_level      <- conf_level
  res$family          <- object$family
  res$kappa           <- object$kappa
  res$log_lik         <- object$log_lik
  res$cov_offset_used <- !(is.null(object$cov_offset) ||
                             all(object$cov_offset == 0))
  if (object$family == "gaussian") {
    res$aic <- 2 * length(object$estimate) - 2 * res$log_lik
  }

  res$call               <- object$call %||% NULL
  res$link_name          <- link_name
  res$invlink_expression <- inv_expr

  class(res) <- "summary.RiskMap"
  return(res)
}

test_that("coef() and summary() match the pre-#92 positional implementation", {
  fits <- list(gaussian_model, gaussian_offset_model, gaussian_intercept_model,
               binomial_model, poisson_model)

  for (fit in fits) {
    expect_equal(coef(fit), old_coef_riskmap(fit), ignore_attr = TRUE)
    expect_equal(summary(fit), old_summary_riskmap(fit), ignore_attr = TRUE)
  }
})

test_that("coef() and summary() match the pre-#92 positional implementation when nugget and sigma2_me are both estimated", {
  set.seed(2)
  n_loc <- 6
  coords_u <- cbind(runif(n_loc, 0, 10000), runif(n_loc, 0, 10000))
  coords_rep <- coords_u[rep(seq_len(n_loc), each = 2), ]
  d <- data.frame(x = coords_rep[, 1], z = coords_rep[, 2], cov = rnorm(2 * n_loc))
  d$y <- 1 + 0.5 * d$cov + rnorm(2 * n_loc)
  dup_data <- sf::st_as_sf(d, coords = c("x", "z"), crs = 32637)

  fit <- glgpm(y ~ cov + gp(nugget = TRUE), data = dup_data, family = "gaussian",
              messages = FALSE)

  expect_equal(coef(fit), old_coef_riskmap(fit), ignore_attr = TRUE)
  expect_equal(summary(fit), old_summary_riskmap(fit), ignore_attr = TRUE)
})
