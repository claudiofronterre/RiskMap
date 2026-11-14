# The golden regression battery of RiskMap (tools/golden/README.md). Each task saves everything a fit or a low-level
# function produces, so that two installs of RiskMap can be compared bit for bit by tools/golden/compare.R. There is
# one output file per task in --out, written under a temporary name and renamed. Tasks whose output exists are
# skipped, so a run can be resumed, and a failed task is saved as a record with its message. Run with one BLAS and
# OpenMP thread, from an install of RiskMap in a library of its own, for example
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --no-init-file tools/golden/run.R \
#     --lib=<library> --out=<arms>/<arm> [--tasks=all | --tasks=F01_binomial,L01_low_binomial | --smoke] [--list]
# The battery moved here from the MCML paper's repository (theory/family-object-golden.R), where it checked every PR of
# RiskMap #253. Its first 57 tasks and their seeds are unchanged.
args <- commandArgs(TRUE)
option <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
flag <- function(name) any(args == paste0("--", name))
known <- c("lib", "out", "tasks", "data", "commit")
unknown_options <- args[!(sub("=.*", "", sub("^--", "", args)) %in% known & grepl("=", args)) &
                        !(args %in% c("--smoke", "--list", "--any-threads"))]
if (length(unknown_options)) stop("unknown options: ", paste(unknown_options, collapse = " "))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
here <- if (is.na(script)) file.path("tools", "golden") else dirname(normalizePath(script))

# Fits reproduce bit for bit only with one BLAS and OpenMP thread, since threaded sums change their order.
threads <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"))
if (any(threads != "1") && !flag("any-threads")) {
  stop("set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 (or pass --any-threads to run anyway)")
}
if ("RhpcBLASctl" %in% loadedNamespaces()) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}
suppressMessages({
  library(RiskMap, lib.loc = option("lib", .libPaths()[1]))
  library(sf)
})
ns <- asNamespace("RiskMap")
`%||%` <- function(a, b) if (is.null(a)) b else a
ctrl <- set_control_mcmc(n_sim = 1200, burnin = 200, thin = 2)

# Data. Simulated surveys, the survey n05-2 of the MCML paper's nugget study and the household survey r01-1 of its
# study of re() terms (saved exactly in data/surveys.rds, see the README), and the package's surveys.
sim <- function(n, target, nugget, family, seed) {
  set.seed(seed); coords <- cbind(runif(n), runif(n)) * 100; km <- as.matrix(dist(coords))
  S <- as.numeric(rnorm(n) %*% chol(2 * exp(-km / 15))) + if (nugget) rnorm(n, 0, 0.5) else 0
  x <- rnorm(n); inv <- if (family == "binomial") plogis else exp
  b <- uniroot(function(b) mean(inv(b + S + 0.3 * x)) - target, c(-20, 20))$root
  m <- if (family == "binomial") rep(20, n) else rep(1, n)
  y <- if (family == "binomial") rbinom(n, m, inv(b + S + 0.3 * x)) else rpois(n, m * inv(b + S + 0.3 * x))
  st_as_sf(data.frame(y = y, m = m, x = x, g = rep(seq_len(ceiling(n / 5)), length.out = n),
                      east = coords[, 1] * 1000, north = coords[, 2] * 1000), coords = c("east", "north"), crs = 32632)
}
saved_survey <- function(name) {
  s <- readRDS(file.path(option("data", file.path(here, "data")), "surveys.rds"))[[name]]
  st_as_sf(s$frame, coords = c("east", "north"), crs = s$epsg)
}
data_cache <- new.env()
dataset <- function(name) {
  if (!is.null(data_cache[[name]])) return(data_cache[[name]])
  value <- switch(name,
    d1 = sim(120, .3, FALSE, "binomial", 11), d2 = sim(120, .05, TRUE, "binomial", 12),
    d3 = sim(120, 1, TRUE, "poisson", 13), d4 = sim(100, .9, FALSE, "binomial", 14),
    no_positives = { d <- dataset("d1"); d$y <- 0L; d },
    separated = { d <- dataset("d1"); d$y <- ifelse(d$x > 0, d$m, 0L); d },
    sparse_small = { d <- dataset("d1")[1:40, ]; set.seed(31); d$y <- rbinom(40, 20, 0.01); d },
    poisson_sparse = { d <- dataset("d3")[1:60, ]; set.seed(32); d$y <- rpois(60, 0.02); d },
    nearly_all = { d <- dataset("d4")[1:60, ]; d$y <- d$m; d$y[1:3] <- 19L; d },
    gaussian = { d <- dataset("d1"); set.seed(33); km <- as.matrix(dist(st_coordinates(d))) / 1000
                 d$z <- 1 + 0.5 * d$x + as.numeric(rnorm(nrow(d)) %*% chol(exp(-km / 15))) + rnorm(nrow(d), 0, 0.3); d },
    liberia = { data("liberia", package = "RiskMap", envir = environment()); st_as_sf(liberia, coords = c("long", "lat"), crs = 4326) },
    loaloa = { data("loaloa", package = "RiskMap", envir = environment()); loaloa },
    tz_malaria = { data("tz_malaria", package = "RiskMap", envir = environment()); tz_malaria },
    n05_2 = saved_survey("n05_2"),
    households = saved_survey("households"),
    stop("unknown data set ", name))
  assign(name, value, envir = data_cache)
  value
}

fit_binomial <- function(formula, data, ..., seed) {
  glgpm(formula, denominator = m, data = dataset(data), family = "binomial", distance_units = "km", messages = FALSE,
        seed = seed, control_mcmc = ctrl, ...)
}
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
texts <- function(...) lapply(list(...), function(x) capture.output(print(x)))

# Every fit a task makes at its top level is kept in its record beside the value (record$fits), for the tolerance
# mode of compare.R. The value itself is what the identity check compares.
made_fits <- new.env()
keep_fit <- function(name, fit) assign(name, fit, envir = made_fits)

# The fits, each a function returning the fit. The fixed part of the call is in fit_binomial().
fits <- list(
  F01_binomial = function() fit_binomial(y ~ x + gp(), "d1", seed = 1L),
  F02_binomial_nugget = function() fit_binomial(y ~ x + gp(nugget = TRUE), "d2", seed = 2L),
  F03_binomial_high = function() fit_binomial(y ~ gp(), "d4", seed = 5L),
  F04_binomial_re = function() fit_binomial(y ~ x + gp() + re(g), "d1", seed = 3L),
  F05_poisson_nugget = function() glgpm(y ~ x + gp(nugget = TRUE), denominator = m, data = dataset("d3"),
                                        family = "poisson", distance_units = "km", messages = FALSE, seed = 4L,
                                        control_mcmc = ctrl),
  F06_liberia = function() glgpm(npos ~ log(elevation) + gp(), denominator = ntest, data = dataset("liberia"),
                                 family = "binomial", messages = FALSE, seed = 6L, control_mcmc = ctrl),
  F07_no_positives = function() fit_binomial(y ~ x + gp(), "no_positives", seed = 9L, return_samples = FALSE),
  F08_separated = function() fit_binomial(y ~ x + gp(), "separated", seed = 10L, return_samples = FALSE),
  F09_sparse_small = function() fit_binomial(y ~ gp(nugget = TRUE), "sparse_small", seed = 11L, return_samples = FALSE),
  F10_poisson_sparse = function() glgpm(y ~ gp(nugget = TRUE), denominator = m, data = dataset("poisson_sparse"),
                                        family = "poisson", distance_units = "km", messages = FALSE, seed = 12L,
                                        return_samples = FALSE),
  F11_nearly_all_positive = function() fit_binomial(y ~ gp(), "nearly_all", seed = 13L, return_samples = FALSE),
  F12_tz_malaria_nugget = function() glgpm(Pf ~ EVI + gp(nugget = TRUE), denominator = Ex, data = dataset("tz_malaria"),
                                           family = "binomial", messages = FALSE, seed = 7L, control_mcmc = ctrl),
  F13_loaloa_nugget_auto = function() glgpm(NO_INF ~ ELEVATION + gp(nugget = TRUE), denominator = NO_EXAM,
                                            data = dataset("loaloa"), family = "binomial", messages = FALSE, seed = 8L,
                                            control_mcmc = ctrl),
  F14_loaloa_nugget_sampled = function() glgpm(NO_INF ~ ELEVATION + gp(nugget = TRUE), denominator = NO_EXAM,
                                               data = dataset("loaloa"), family = "binomial", messages = FALSE,
                                               seed = 8L, control_mcmc = ctrl,
                                               control_mcml = set_control_mcml(integrate_noise = FALSE)),
  F15_loaloa_nugget_integrated = function() glgpm(NO_INF ~ ELEVATION + gp(nugget = TRUE), denominator = NO_EXAM,
                                                  data = dataset("loaloa"), family = "binomial", messages = FALSE,
                                                  seed = 8L, control_mcmc = ctrl,
                                                  control_mcml = set_control_mcml(integrate_noise = TRUE)),
  F16_link_pnorm_bare = function() fit_binomial(y ~ x + gp(), "d1", seed = 1L, invlink = pnorm),
  F17_link_pnorm_closure = function() fit_binomial(y ~ x + gp(), "d1", seed = 1L, invlink = function(e) pnorm(e)),
  F18_link_pnorm_list = function() fit_binomial(y ~ x + gp(), "d1", seed = 1L,
                                                invlink = list(inv = pnorm, d1 = dnorm, d2 = function(e) -e * dnorm(e))),
  F19_named_probit = function() glgpm(y ~ x + gp(), denominator = m, data = dataset("d1"), family = binomial("probit"),
                                      distance_units = "km", messages = FALSE, seed = 1L, control_mcmc = ctrl),
  F20_named_cloglog = function() glgpm(y ~ x + gp(), denominator = m, data = dataset("d1"), family = binomial("cloglog"),
                                       distance_units = "km", messages = FALSE, seed = 1L, control_mcmc = ctrl),
  F21_named_cauchit = function() glgpm(y ~ x + gp(), denominator = m, data = dataset("d1"), family = binomial("cauchit"),
                                       distance_units = "km", messages = FALSE, seed = 1L, control_mcmc = ctrl),
  F22_imperfect_test = function() fit_binomial(y ~ x + gp(), "d1", seed = 1L,
                                               invlink = function(e) 0.01 + 0.89 * plogis(e)),
  F23_imperfect_nugget_integrated = function() fit_binomial(y ~ x + gp(nugget = TRUE), "d2", seed = 2L,
                                                            invlink = function(e) 0.01 + 0.89 * plogis(e),
                                                            control_mcml = set_control_mcml(integrate_noise = TRUE)),
  F24_loaloa_probit = function() glgpm(NO_INF ~ ELEVATION + gp(), denominator = NO_EXAM, data = dataset("loaloa"),
                                       family = binomial("probit"), messages = FALSE, seed = 1L, control_mcmc = ctrl),
  F25_loaloa_cauchit = function() glgpm(NO_INF ~ ELEVATION + gp(), denominator = NO_EXAM, data = dataset("loaloa"),
                                        family = binomial("cauchit"), messages = FALSE, seed = 1L, control_mcmc = ctrl),
  F26_households = function() glgpm(y ~ x + gp(kappa = 0.5) + re(household), denominator = m,
                                    data = dataset("households"), family = "binomial", distance_units = "km",
                                    messages = FALSE, seed = 26L, return_samples = FALSE),
  F27_held_nugget_logit = function() fit_binomial(y ~ x + gp(kappa = 0.5, nugget = TRUE), "n05_2", seed = 27L,
                                                  return_samples = FALSE,
                                                  control_mcml = set_control_mcml(integrate_noise = FALSE)),
  F28_held_nugget_probit = function() glgpm(y ~ x + gp(kappa = 0.5, nugget = TRUE), denominator = m,
                                            data = dataset("n05_2"), family = binomial("probit"),
                                            distance_units = "km", messages = FALSE, seed = 27L,
                                            return_samples = FALSE,
                                            control_mcml = set_control_mcml(integrate_noise = FALSE)),
  F29_gaussian = function() glgpm(z ~ x + gp(), data = dataset("gaussian"), family = "gaussian",
                                  distance_units = "km", messages = FALSE)
)

# What a fit task keeps: the object itself (compared after timings and closures are removed) and its printed text.
run_fit <- function(make) {
  fit <- quiet(make())
  keep_fit("fit", fit)
  check <- quiet(check_fit(fit, messages = FALSE))
  list(fit = fit, check = check,
       text = c(texts(fit, summary(fit), check), list(capture.output(print(summary(check))))))
}

# Downstream of five fits, each refitted in its task so that tasks are independent. The Gaussian fit has a
# closed-form prediction, an exact profile and the checks that say a Gaussian fit has no Monte Carlo error.
base_fits <- list(logit = fits$F01_binomial, poisson = fits$F05_poisson_nugget, probit = fits$F19_named_probit,
                  integrated = function() fit_binomial(y ~ x + gp(nugget = TRUE), "d2", seed = 2L,
                                                       control_mcml = set_control_mcml(integrate_noise = TRUE)),
                  gaussian = fits$F29_gaussian)
base_data <- c(logit = "d1", poisson = "d3", probit = "d1", integrated = "d2", gaussian = "gaussian")
prediction_grid <- function(data) {
  points <- st_make_grid(st_as_sfc(st_bbox(data)), n = c(8, 8), what = "centers")
  st_sf(x = rep(0, length(points)), m = rep(1, length(points)), geometry = points)
}
downstream <- list(
  predict = function(fit, data) {
    grid <- prediction_grid(data)
    list(grid = predict(fit, newdata = grid, control_mcmc = ctrl, seed = 21L, messages = FALSE),
         observed = predict(fit, control_mcmc = ctrl, seed = 22L, messages = FALSE),
         uncertain = predict(fit, newdata = grid, control_mcmc = ctrl, parameter_uncertainty = TRUE,
                             parameter_draws = 10L, seed = 23L, messages = FALSE),
         simulated = simulate_glgpm(fit, nsim = 2, what = c("data", "surface"), prediction_grid = grid, seed = 24L))
  },
  bootstrap = function(fit, data) {
    boot <- parametric_bootstrap(fit, n_sim = 10, seed = 31L, messages = FALSE)
    list(boot = boot, text = texts(boot, summary(boot)))
  },
  profile = function(fit, data) {
    prof <- profile(fit, parameters = "phi", n_draws = 2000L, seed = 41L, messages = FALSE)
    list(profile = prof, text = texts(prof))
  },
  extras = function(fit, data) {
    check <- check_fit(fit, extra = c("profile", "radius"), messages = FALSE)
    list(check = check, text = c(texts(check), list(capture.output(print(summary(check))))))
  },
  assess = function(fit, data) {
    cv <- assess_prediction(list(model = fit), method = "cluster", fold = 2, iter = 1, control_mcmc = ctrl,
                            seed = 51L, plot_fold = FALSE, messages = FALSE)
    simulated <- simulate_glgpm(fit, nsim = 1, what = c("data", "surface"),
                                prediction_grid = prediction_grid(data), seed = 52L)
    assessed <- assess_simulation(simulated, models = list(candidate = fit), control_mcmc = ctrl, seed = 53L,
                                  spatial_scale = "grid", pred_objective = "mse", messages = FALSE)
    list(cv = cv, assessed = assessed)
  }
)

# Low-level results at fixed inputs, through engine internals. A call that fails is kept as its message.
attempt <- function(expr) tryCatch(expr, error = function(e) list(error = conditionMessage(e)))
links <- list(binomial = list(canonical = NULL, probit = c(ns$named_inverse_link("probit"), link = "probit"),
                              cloglog = c(ns$named_inverse_link("cloglog"), link = "cloglog"),
                              cauchit = c(ns$named_inverse_link("cauchit"), link = "cauchit"),
                              imperfect = function(e) 0.01 + 0.89 * plogis(e)),
              poisson = list(canonical = NULL, softplus = function(e) log1p(exp(e))))
# Before step 1d of #253 the engine took a family string with a link through its old entry points, and from 1d it
# takes the family object, whose methods those entry points called. Each call below goes the way of the commit
# tested, so that arms before 1d can still be run.
old_entry_points <- exists("laplace_observation_terms", envir = ns, inherits = FALSE)
low_level <- function(family) {
  eta <- c(-800, -745, -60, -37, -36, -27, -20, -5, -1, 0, 1, 5, 8.5, 20, 36, 37, 60, 300, 700)
  y <- rep(c(0, 3, 20), length.out = length(eta)); m <- rep(20, length(eta))
  step <- outer(rep(1, length(eta)), c(-3, -0.5, 0.25, 2))
  out <- list()
  for (name in names(links[[family]])) {
    linkf <- attempt(ns$make_inverse_link_functions(family, links[[family]][[name]], length(eta)))
    tag <- paste(family, name, sep = "_")
    out[[paste0("linkf_", tag)]] <- if (is.function(linkf$inv)) lapply(linkf[c("inv", "d1", "d2")], function(f) f(eta)) else linkf
    if (!is.function(linkf$inv)) next
    canonical <- identical(linkf$name, "canonical")
    out[[paste0("terms_", tag)]] <- attempt(if (old_entry_points) {
      ns$laplace_observation_terms(eta, y, m, family, linkf)
    } else {
      ns$make_family(family, linkf)$terms(eta, y, m)
    })
    out[[paste0("custom_terms_", tag)]] <- attempt(ns$compute_custom_link_terms(eta, y, m, family, linkf))
    out[[paste0("loglik_", tag)]] <- attempt(if (old_entry_points) {
      ns$compute_observation_log_likelihood(rbind(eta, eta / 2), y, m, family, linkf)
    } else {
      ns$make_family(family, linkf)$log_likelihood(rbind(eta, eta / 2), y, m)
    })
    centre <- pmax(pmin(eta, 30), -30)
    out[[paste0("nonlinear_", tag)]] <- attempt(if (old_entry_points) {
      ns$compute_nonlinear_parts(centre, step, y, m, family, linkf, canonical)
    } else {
      ns$make_family(family, linkf)$nonlinear_parts(centre, step, y, m)
    })
    out[[paste0("slope_", tag)]] <- attempt(if (old_entry_points) {
      ns$compute_variance_slope_terms(
        rbind(centre / 3, centre / 6), y, m, family, rep(1:4, length.out = length(eta)), if (!canonical) linkf)
    } else {
      ns$compute_slope_terms(rbind(centre / 3, centre / 6), y, m, ns$make_family(family, if (!canonical) linkf),
                             rep(1:4, length.out = length(eta)))
    })
    out[[paste0("arrays_", tag)]] <- if (!canonical) attempt(ns$compute_custom_kernel_arrays(linkf, family, low = -40, high = 40, step = 0.25))
  }
  # Engines at fixed parameters and a fixed bank, sampled and integrated, canonical and one other link.
  n <- 40; set.seed(9); co <- cbind(runif(n), runif(n)); dd <- as.matrix(dist(co))
  yy <- if (family == "binomial") rbinom(n, 10, .4) else rpois(n, 2)
  mm <- if (family == "binomial") rep(10, n) else rep(1, n)
  other <- if (family == "binomial") links$binomial$probit else links$poisson$softplus
  for (link_name in c("canonical", "other")) for (tau in c(FALSE, TRUE)) for (noise in c(FALSE, TRUE)) {
    if (noise && !tau) next
    tag <- paste(family, link_name, if (tau) "nugget" else "plain", if (noise) "integrated" else "sampled", sep = "_")
    out[[paste0("engine_", tag)]] <- attempt({
      en <- ns$make_adaptive_mcml_engine(D = matrix(1, n, 1), y = yy, units_m = mm, distances = dd, family = family,
                                         fix_tau2 = tau, invlink = if (link_name == "other") other,
                                         integrate_noise = noise)
      par <- c(-0.3, 0.1, -1, if (tau) -1)
      set.seed(10); bank <- matrix(rnorm(300 * en$latent_dimension), 300)
      en$set_anchor(par)
      list(global = en$proposal_at(par)[c("centre", "root", "type", "factors")],
           laplace = en$proposal_at(par, "laplace")[c("centre", "root", "type")],
           weights = en$evaluate(par, bank, log_weights_only = TRUE),
           objective = en$evaluate(par + 0.01, bank))
    })
  }
  Sig <- exp(-dd / 0.3)
  out$maxim <- attempt(maxim_integrand(y = yy, denominator = mm, mu = rep(-0.2, n), Sigma = Sig, ID_coords = 1:n,
                                       family = family, hessian = TRUE, gradient = TRUE))
  out$independent_sites <- attempt({
    likelihood <- ns$make_independent_site_likelihood(matrix(1, n, 1), yy, mm, seq_len(n), numeric(n),
                                                      if (old_entry_points) family else ns$make_family(family))
    lapply(c(-2, 0.5), function(b) likelihood$log_likelihood(b, log(0.8)))
  })
  out
}
bounded_driver <- function() {
  set.seed(21); n <- 60; co <- cbind(runif(n), runif(n)); X <- cbind(1, rnorm(n), rnorm(n)); yb <- rep(0L, n)
  en <- ns$make_adaptive_mcml_engine(D = X, y = yb, units_m = rep(20, n), distances = as.matrix(dist(co)),
                                     family = "binomial")
  out <- quiet(ns$fit_adaptive_mcml(engine = en, start = c(-3, 0, 0, log(0.5), log(0.2)), initial_budget = 500,
                                    validation_budget = 500, seed = 5, lower = c(-5, -30, -30, -20, -20),
                                    upper = c(30, 30, 30, 20, 10), variance_index = 4L, linked_index = list(`4` = 5L),
                                    coefficient_index = 1:3))
  out[c("outcome", "estimate", "certificate", "log_likelihood")]
}

tasks <- c(
  lapply(fits, function(make) function() run_fit(make)),
  unlist(lapply(names(base_fits), function(base) {
    stats::setNames(lapply(names(downstream), function(what) {
      force(what); function() {
        fit <- quiet(base_fits[[base]]())
        keep_fit("base", fit)
        quiet(downstream[[what]](fit, dataset(base_data[[base]])))
      }
    }), paste0("D_", base, "_", names(downstream)))
  }), recursive = FALSE),
  list(L01_low_binomial = function() low_level("binomial"), L02_low_poisson = function() low_level("poisson"),
       L03_bounded_driver = bounded_driver)
)

# The smoke subset, about a minute on one Wayland core and a few on a laptop: fits of each family, a nugget, an re()
# term, a named link, both saved surveys, a prediction of a sampled and of a Gaussian fit, and the low-level tasks.
smoke <- c("F01_binomial", "F02_binomial_nugget", "F04_binomial_re", "F05_poisson_nugget", "F19_named_probit",
           "F26_households", "F27_held_nugget_logit", "F29_gaussian", "D_logit_predict", "D_gaussian_predict",
           "L01_low_binomial", "L02_low_poisson", "L03_bounded_driver")

if (flag("list")) { cat(names(tasks), sep = "\n"); quit(save = "no") }
out_dir <- option("out")
if (is.null(out_dir)) stop("give the arm's output directory with --out=<dir>")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wanted <- if (flag("smoke")) "smoke" else option("tasks", "all")
wanted <- switch(wanted, all = names(tasks), smoke = smoke, strsplit(wanted, ",")[[1]])
unknown <- setdiff(wanted, names(tasks))
if (length(unknown)) stop("unknown tasks: ", paste(unknown, collapse = ", "))
environment_record <- list(riskmap = as.character(packageVersion("RiskMap")),
                           riskmap_path = find.package("RiskMap"),
                           commit = option("commit", Sys.getenv("RISKMAP_COMMIT")),
                           r = R.version.string, host = Sys.info()[["nodename"]],
                           blas = tryCatch(extSoftVersion()[["BLAS"]], error = function(e) NA_character_),
                           threads = as.list(threads))
cat(sprintf("RiskMap %s from %s, commit %s, on %s\n", environment_record$riskmap, environment_record$riskmap_path,
            environment_record$commit, environment_record$host))
counts <- c(done = 0L, failed = 0L, exists = 0L)
for (task in wanted) {
  path <- file.path(out_dir, paste0(task, ".rds"))
  if (file.exists(path)) { cat(sprintf("%-34s exists\n", task)); counts[["exists"]] <- counts[["exists"]] + 1L; next }
  rm(list = ls(made_fits), envir = made_fits)
  # Every task sets its own seeds. This one only keeps a task that draws before setting one independent of the tasks
  # run before it in the process, so that any subset of tasks, in any order, gives the same outputs.
  codes <- utf8ToInt(task)
  set.seed(sum(codes * seq_along(codes)))
  started <- proc.time()[["elapsed"]]
  gc_before <- sum(gc.time()[1:2])
  value <- tryCatch(tasks[[task]](), error = function(e) structure(list(message = conditionMessage(e)), class = "failed"))
  seconds <- proc.time()[["elapsed"]] - started
  failed <- inherits(value, "failed")
  record <- list(task = task, status = if (failed) "failed" else "done",
                 value = if (!failed) value, message = if (failed) value$message,
                 seconds = seconds, gc_seconds = sum(gc.time()[1:2]) - gc_before,
                 fits = as.list(made_fits),
                 environment = c(environment_record,
                                 list(packages = vapply(sort(loadedNamespaces()),
                                                        function(p) as.character(packageVersion(p)), ""))))
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(record, temporary)
  file.rename(temporary, path)
  counts[[record$status]] <- counts[[record$status]] + 1L
  cat(sprintf("%-34s %-6s %7.1f s %s\n", task, record$status, seconds, record$message %||% ""))
}
cat(sprintf("%d done, %d failed, %d already there, in %s\n", counts[["done"]], counts[["failed"]], counts[["exists"]],
            out_dir))
