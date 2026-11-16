# Precomputes what the articles load instead of computing while the website is
# built: downloads (boundaries and population counts) and the slower steps
# (profiles, the bootstrap and the simulation study). Each section repeats the
# code that its article shows, in the same order, so the results are those the
# reader would get. Run it from the package root, with the RiskMap that the
# articles document installed, and rerun it after any change that moves
# RiskMap's numbers:
#
#   Rscript vignettes/articles/precompute.R
#
# It needs the geobounds and terra packages and an internet connection, and
# takes about 6 minutes on 4 cores.

library(RiskMap)
library(sf)
library(geobounds)
library(terra)

out <- "vignettes/articles/precomputed"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

# Liberia: the data, projected as in every article ----------------------------

data(liberia)
liberia <- st_as_sf(liberia,
                    coords = c("long", "lat"),
                    crs = 4326)
crs <- propose_crs(liberia, messages = FALSE)
liberia <- st_transform(liberia, crs)

# Boundaries of Liberia and its counties, from geoBoundaries ------------------

border <- st_transform(gb_get_adm0("Liberia"), crs)
counties <- st_transform(gb_get_adm1("Liberia"), crs)[, "shapeName"]
saveRDS(list(border = border, counties = counties),
        file.path(out, "liberia-boundaries.rds"))

# Population counts of Liberia on a 5 km grid, from WorldPop ------------------

population <- get_worldpop("LBR", year = 2020, resolution = "1km",
                           messages = FALSE)
population <- project(population, crs(vect(border)),
                      method = "sum", res = 5000)
names(population) <- "population"
grid_population <- create_grid(border, raster = population)
saveRDS(grid_population, file.path(out, "liberia-population-grid.rds"))

# Checking a fit: profiles and the parametric bootstrap -----------------------

set.seed(2026)
fit <- glgpm(npos ~ log_elevation + gp(),
             data = liberia,
             family = "binomial",
             denominator = ntest,
             messages = FALSE)
profiles <- profile(fit,
                    parameters = c("log_elevation", "sigma2", "phi"),
                    pairs = list(c("sigma2", "phi")),
                    cores = 4,
                    messages = FALSE)
bootstrap <- parametric_bootstrap(fit,
                                  n_sim = 200,
                                  cores = 4,
                                  messages = FALSE)
saveRDS(list(profiles = profiles, bootstrap = bootstrap),
        file.path(out, "checking-a-fit.rds"))

# Assessing predictions and simulations: the simulation study -----------------

set.seed(2026)
fit_spatial <- glgpm(npos ~ gp(),
                     data = liberia,
                     family = "binomial",
                     denominator = ntest,
                     messages = FALSE)
fit_elevation <- glgpm(npos ~ log_elevation + gp(),
                       data = liberia,
                       family = "binomial",
                       denominator = ntest,
                       messages = FALSE)
truth <- specify_glgpm(npos ~ gp(),
                       data = liberia,
                       family = "binomial",
                       denominator = ntest,
                       parameters = coef(fit_spatial))
grid <- create_grid(border, spacing = 5)
simulations <- simulate_glgpm(truth,
                              nsim = 50,
                              what = c("data", "surface"),
                              prediction_grid = grid,
                              seed = 1)
fit_smooth <- glgpm(npos ~ gp(kappa = 1.5),
                    data = liberia,
                    family = "binomial",
                    denominator = ntest,
                    messages = FALSE)
assessment <- assess_simulation(simulations,
                                models = list(exponential = fit_spatial,
                                              smooth = fit_smooth),
                                spatial_scale = c("location", "area"),
                                aggregate = "mean",
                                boundaries = counties,
                                label = "shapeName",
                                objective = c("mse", "classify"),
                                breaks = c(0, 0.1, 0.2, 1),
                                seed = 1,
                                messages = FALSE)
saveRDS(list(assessment = assessment),
        file.path(out, "assessing-predictions.rds"))

# Prediction targets: boundaries of Ghana and its regions ---------------------

data(malnutrition)
ghana_crs <- propose_crs(malnutrition, messages = FALSE)
ghana <- st_transform(gb_get_adm0("Ghana"), ghana_crs)
regions <- st_transform(gb_get_adm1("Ghana"), ghana_crs)
saveRDS(list(ghana = ghana, regions = regions[, "shapeName"]),
        file.path(out, "ghana-boundaries.rds"))
