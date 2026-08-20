# =============================================================================
# PHASE 2: DOMAIN CLEANING & PIXEL-LEVEL THINNING
# =============================================================================

library(terra)

# 1. Define File Paths
tif_path     <- "C:/VERA/Variables/Turkiye_Current_Bioclimatics_and_Envirems/Bio1.tif"
raw1_path    <- "C:/VERA/Occurrences/Raw/Sitta_krueperi_RAW_1.csv"
raw2_path    <- "C:/VERA/Occurrences/Raw/Sitta_krueperi_RAW_2.csv"
thinned_path <- "C:/VERA/Occurrences/Raw/Sitta_krueperi_THINNED.csv"

# 2. Import Data
bio1 <- rast(tif_path)
occ  <- read.csv(raw1_path)

# 3. Stage 1: Domain Cleaning (Remove points outside valid raster cells)
pts <- vect(occ, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
ext_vals <- terra::extract(bio1, pts)
occ_clean <- occ[!is.na(ext_vals[, 2]), ]
write.csv(occ_clean, raw2_path, row.names = FALSE, quote = FALSE)

# 4. Stage 2: Pixel-Level Spatial Thinning
cell_nums <- cellFromXY(bio1, as.matrix(occ_clean[, c("Longitude", "Latitude")]))
occ_thinned <- occ_clean[!duplicated(cell_nums), ]
write.csv(occ_thinned, thinned_path, row.names = FALSE, quote = FALSE)

cat("========================================\n")
cat("PHASE 2 COMPLETED\n")
cat("========================================\n")
cat("Initial (RAW_1)           :", nrow(occ), "points\n")
cat("After NA Cleaning (RAW_2):", nrow(occ_clean), "points\n")
cat("After Thinning (THINNED)  :", nrow(occ_thinned), "points\n")
cat("========================================\n")



# =============================================================================
# PHASE 3: GEOGRAPHIC OUTLIER ISOLATION (97%)
# =============================================================================

library(terra)

# 1. Define File Paths
thinned_path <- "C:/VERA/Occurrences/Raw/Sitta_krueperi_THINNED.csv"
out_dir      <- "C:/VERA/Occurrences/Edited"
final_path   <- paste0(out_dir, "/Sitta_krueperi.csv")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 2. Import Thinned Data
occ_thinned <- read.csv(thinned_path)

# 3. Stage 3: Geographic Isolation
pts <- vect(occ_thinned, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
# Project to a metric CRS for accurate distance and covariance calculations
pts_metric <- terra::project(pts, "EPSG:3857") 

coords <- crds(pts_metric)
cov_mat   <- stats::cov(coords)
center_pt <- colMeans(coords)

if (det(cov_mat) > 1e-6) {
  md <- stats::mahalanobis(coords, center = center_pt, cov = cov_mat)
  cutoff <- stats::quantile(md, probs = 0.97, na.rm = TRUE)
  valid_idx <- which(md <= cutoff)
  occ_final <- occ_thinned[valid_idx, ]
} else {
  # If spatial covariance cannot be estimated reliably
  occ_final <- occ_thinned 
}

# 4. Save Final VERA Occurrence Table
write.csv(occ_final, final_path, row.names = FALSE, quote = FALSE)

n_removed <- nrow(occ_thinned) - nrow(occ_final)
cat("========================================\n")
cat("PHASE 3 COMPLETED\n")
cat("========================================\n")
cat("Initial (THINNED)    :", nrow(occ_thinned), "points\n")
cat("Removed (3% Outlier) :", n_removed, "points\n")
cat("Final VERA Data      :", nrow(occ_final), "points\n")
cat("Saved Path           :", final_path, "\n")
cat("========================================\n")