# =============================================================================
# PHASE 1: ENVIREM VARIABLE GENERATION
# =============================================================================
# Workflow overview:
#   1. Load twelve monthly layers for each required climate input.
#   2. Rename the layers using the naming convention expected by ENVIREM.
#   3. Derive monthly extraterrestrial solar radiation.
#   4. Generate the complete ENVIREM set and save one GeoTIFF per variable.
#
# Note: All input rasters must have matching geometry, spatial resolution, CRS, 
# and extent.

library(terra)
library(envirem)

# 1. Directory Paths
input_dir  <- "C:/VERA/Variables/Turkiye_Monthly_Bioclimatics"
output_dir <- "C:/VERA/Variables/Turkiye_Envirems"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. Reading and Loading Files
cat(">>> Loading monthly climate data...\n")
prec_files <- list.files(input_dir, pattern = "prec", full.names = TRUE)
tmin_files <- list.files(input_dir, pattern = "tmin", full.names = TRUE)
tmax_files <- list.files(input_dir, pattern = "tmax", full.names = TRUE)
tavg_files <- list.files(input_dir, pattern = "tavg", full.names = TRUE)

prec  <- rast(sort(prec_files))
tmin  <- rast(sort(tmin_files))
tmax  <- rast(sort(tmax_files))
tmean <- rast(sort(tavg_files)) 

# 3. Adapting Naming to ENVIREM Default Template
months <- sprintf("%02d", 1:12) 

names(prec)  <- paste0("precip_", months)
names(tmin)  <- paste0("tmin_", months)
names(tmax)  <- paste0("tmax_", months)
names(tmean) <- paste0("tmean_", months)

master_clim <- c(prec, tmin, tmax, tmean)

cat(">>> Verifying raster names:\n")
verifyRasterNames(master_clim)

# 4. Solar Radiation Calculation
cat(">>> Calculating monthly solar radiation (ET solrad)...\n")
# Note: Year 50 represents the year 2000 (2000 - 1950) in the palinsol package.
solar <- ETsolradRasters(tmean[[1]], 
                         year = 50, 
                         outputDir = NULL)

names(solar) <- paste0("et_solrad_", months)

# 5. Generate and Save ENVIREM Variables
cat(">>> Generating 17 ENVIREM variables (this may take a while)...\n")
allGrids <- generateEnvirem(masterstack = master_clim, solradstack = solar, var = 'all')

filenames <- paste0(output_dir, '/', names(allGrids), '.tif')
writeRaster(allGrids, filenames, overwrite = TRUE)

cat("========================================\n")
cat("PHASE 1 COMPLETED SUCCESSFULLY!\n")
cat("========================================\n")