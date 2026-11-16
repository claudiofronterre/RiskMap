# Put the covariates of tz_covariates on the scales of tz_malaria.
#
# tz_malaria takes Temperature and NTL from two rasters of Tanzania for 2015,
# which tz_covariates did not use. This script replaces those two columns of
# tz_covariates with the values of the same rasters at its grid points, gives
# ITN as a percentage, as in tz_malaria, and marks missing values with NA
# instead of -9999. The other columns, EVI, Precipitation and Population, were
# already on the scales of tz_malaria and are left as they are.
#
# The rasters are in longitude and latitude (EPSG:4326).
#   - Tanzania_Annual_LST_2015.tif: MODIS daytime land surface temperature
#     (band LST_Day_1km, about 1 km), the mean over 2015, in degrees Celsius.
#   - Tanzania_Annual_NightLights_2015.tif: VIIRS Day/Night Band monthly
#     composites (band avg_rad, about 500 m), the mean over 2015, in
#     nanowatts per square centimetre per steradian.
#
# Run from the root of the package, with the folder that holds the rasters in
# the environment variable TZ_RASTERS, as in
#   TZ_RASTERS=<folder> Rscript data-raw/tz_covariates.R
# Running it again on its own output changes nothing.

library(sf)
library(terra)

raster_dir <- Sys.getenv("TZ_RASTERS")
if (!nzchar(raster_dir) || !dir.exists(raster_dir)) {
  stop("Set TZ_RASTERS to the folder that holds the rasters.")
}
rasters <- c(
  Temperature = "Tanzania_Annual_LST_2015.tif",
  NTL = "Tanzania_Annual_NightLights_2015.tif"
)

load("data/tz_covariates.rda")
points <- vect(st_transform(tz_covariates, 4326))
for (name in names(rasters)) {
  values <- extract(rast(file.path(raster_dir, rasters[[name]])), points,
                    ID = FALSE)[[1]]
  tz_covariates[[name]] <- as.numeric(values)
}

# ITN was a proportion with -9999 for missing values, and is a percentage in
# tz_malaria.
missing_itn <- which(tz_covariates$ITN == -9999)
tz_covariates$ITN[missing_itn] <- NA
if (max(tz_covariates$ITN, na.rm = TRUE) <= 1) {
  tz_covariates$ITN <- 100 * tz_covariates$ITN
}

save(tz_covariates, file = "data/tz_covariates.rda", compress = "xz")
