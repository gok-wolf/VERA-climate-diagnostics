#!/usr/bin/env Rscript
# =============================================================================
# VERA -- Variable Ecological Restriction Analysis
# Sitta krueperi single-species tutorial pipeline (19-predictor version)
# -----------------------------------------------------------------------------
# Self-contained script (no external config file, no external profile loop).
# Derived from VERA_pipeline_North_America_2.R (2026-08-15 canonical revision).
# Executes the full canonical VERA workflow -- deterministic execution, explicit
# tail fallback, co-dominance diagnostics, consistent ridge geometry, integrated
# add-on diagnostics (addon01 / addon02), machine-readable manifests, code MD5,
# session info, output inventory -- for a single focal taxon under the
# classical 19-predictor (Bio1-Bio19) profile only.
#
# This 19-only variant targets the most common research setup, where the
# WorldClim/CHELSA nineteen bioclimatic variables are the sole predictors.
# A separate 36-predictor script (Bio + Envirems) is provided alongside.
#
# Focal taxon        : Sitta krueperi (Kruper's Nuthatch)
# Species code       : Skr
# Occurrences file   : C:/VERA/Occurrences/Edited/Sitta_krueperi.csv
# Predictor rasters  : C:/VERA/Variables/Turkiye_Current_Bioclimatics_and_Envirems
# Results depot      : C:/VERA/Results/19/Skr_current
#
# Tail-scale definition (asymmetric VRS denominator) is anchored at the global
# occurrence mean via a common-mu RMS distance: this matches the parameterisation
# used by the canonical multi-species pipeline.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tibble)
})

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

# =============================================================================
# SECTION 1 -- CONFIGURATION (hard-coded; no external config file)
# =============================================================================

cfg <- list(

  # Focal-taxon identity.
  species_code         = "Skr",
  species_label        = "Sitta krueperi",
  output_suffix        = "Skr_current",

  # Inputs.
  occ_csv              = "C:/VERA/Occurrences/Edited/Sitta_krueperi.csv",
  current_vars_dir     = "C:/VERA/Variables/Turkiye_Current_Bioclimatics_and_Envirems",

  # No polygon mask required: the Turkiye raster stack is already region-cropped.
  # Set to a shapefile path if you later want to apply an additional boundary.
  mask_shp             = NULL,

  # Result depot. Fixed to the 19-predictor branch.
  output_dir           = "C:/VERA/Results/19/Skr_current",
  temp_subdir          = "temp_terra",

  # Occurrence CSV column detection candidates.
  lon_candidates       = c("decimalLongitude", "Longitude", "longitude", "lon", "x"),
  lat_candidates       = c("decimalLatitude",  "Latitude",  "latitude",  "lat", "y"),

  # Fixed 19-predictor canonical set (Bio1..Bio19).
  focal_vars           = paste0("Bio", 1:19),

  min_tail_n                        = 5,

  # Tail-scale definition used by the asymmetric VRS denominator.
  #   "global_mu_rms" = sqrt(mean((x_tail - global_mu)^2)); canonical, because
  #                     numerator and denominator use the same anchor.
  #   "tail_mean_sd"  = legacy sample SD around each tail's own mean.
  tail_scale_center                 = "global_mu_rms",

  shrinkage_method                  = "n_over_n_plus_k",
  shrinkage_k                       = 10,
  z2_cap                            = 16,
  equality_tolerance                = 1e-12,
  bootstrap_asymmetry               = TRUE,
  bootstrap_iter                    = 1000,
  bootstrap_frac                    = 1.0,
  random_seed                       = 20260809L,
  asymmetry_ratio_threshold         = 1.2,
  truncation_bg_tail_prop_threshold = 0.05,

  # Canonical occupied-climate partition.
  core_occupied_rule                = "mahal_quantile",
  core_occupied_prob                = 0.9,

  # Canonical Mahalanobis geometry and empirical tier calibration.
  mahal_reference                   = "all_occurrences",
  tier_break_reference              = "core_occurrences",
  tier_method                       = "mahal_distance",
  tier_breaks_method                = "empirical_quantile",
  tier_breaks_quantiles             = c(0.80, 0.95, 0.99),
  mahal_distance_breaks             = c(1, 2, 3),
  ridge_sequence                    = c(0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1),
  write_fixed_break_comparator      = FALSE,
  min_valid_frac                    = 0.50,
  divergence_threshold              = 0.10,
  tail_balance_tolerance            = 0.10,
  tail_no_signal_epsilon            = 1e-12,
  latitude_band_width               = 2,

  multimodality_mode                = "warn",
  memfrac                           = 0.8,
  todisk                            = FALSE,
  raster_compression                = c("COMPRESS=LZW", "PREDICTOR=2", "TILED=YES",
                                        "BLOCKXSIZE=256", "BLOCKYSIZE=256",
                                        "NUM_THREADS=ALL_CPUS", "BIGTIFF=IF_SAFER"),
  cleanup_temp_at_end               = TRUE,
  gc_interval                       = 12L,
  write_output_checksums            = TRUE,
  precompute_mahal_reference        = TRUE,
  progress                          = 1,

  write_optional_diagnostic_rasters = TRUE,
  write_index_rasters               = TRUE,
  # Reporting-only attribution and cap-audit metrics.
  write_reporting_audit       = TRUE,
  reporting_contribution_coverage         = 0.80,
  allow_missing_predictors          = FALSE,
  script_file                       = env_or_default("VERA_SCRIPT_FILE", "")
)

resolve_script_file <- function(explicit_path = "") {
  candidates <- character(0)

  if (nzchar(explicit_path)) candidates <- c(candidates, explicit_path)

  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg)) {
    candidates <- c(candidates, sub("^--file=", "", script_arg[1]))
  }

  # source() records the active input file as `ofile` in one of the calling
  # frames. Reading it here lets a normal source("VERA_Sitta_krueperi_Tutorial.R")
  # call produce a reproducible script path and checksum.
  source_files <- unlist(lapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile) && length(frame$ofile)) as.character(frame$ofile[1]) else NULL
  }), use.names = FALSE)
  if (length(source_files)) candidates <- c(candidates, rev(source_files))

  candidates <- c(
    candidates,
    file.path(getwd(), "VERA_Sitta_krueperi_Tutorial.R"),
    file.path(getwd(), "VERA_pipeline.R")
  )
  candidates <- unique(candidates[nzchar(candidates)])

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

# Resolve while the source() calling frame is still available.
analysis_script_file <- resolve_script_file(cfg$script_file)

gdal_opts <- list(gdal = cfg$raster_compression)

validate_config <- function(cfg) {
  if (!is.character(cfg$species_label) || length(cfg$species_label) != 1L ||
      is.na(cfg$species_label) || !nzchar(cfg$species_label)) {
    stop("species_label must be one non-empty canonical taxon label.")
  }
  allowed_shrinkage <- c("n_over_n_plus_k", "sqrt_n_over_sqrt_n_plus_k", "capped_linear")
  if (!cfg$shrinkage_method %in% allowed_shrinkage) {
    stop("Unknown shrinkage_method: ", cfg$shrinkage_method)
  }
  allowed_tail_scale_centers <- c("global_mu_rms", "tail_mean_sd")
  if (!is.character(cfg$tail_scale_center) || length(cfg$tail_scale_center) != 1L ||
      !cfg$tail_scale_center %in% allowed_tail_scale_centers) {
    stop("tail_scale_center must be 'global_mu_rms' or 'tail_mean_sd'.")
  }
  if (!identical(cfg$core_occupied_rule, "mahal_quantile")) {
    stop("The publication pipeline requires core_occupied_rule='mahal_quantile'.")
  }
  if (!identical(cfg$mahal_reference, "all_occurrences")) {
    stop("The publication pipeline requires mahal_reference='all_occurrences'.")
  }
  if (!identical(cfg$tier_break_reference, "core_occurrences")) {
    stop("The canonical empirical-tier workflow requires tier_break_reference='core_occurrences'.")
  }
  if (!identical(cfg$tier_method, "mahal_distance")) {
    stop("The publication pipeline requires tier_method='mahal_distance'.")
  }
  if (!cfg$tier_breaks_method %in% c("empirical_quantile", "fixed", "chi_square")) {
    stop("Unknown tier_breaks_method: ", cfg$tier_breaks_method)
  }
  if (!is.finite(cfg$z2_cap) || cfg$z2_cap <= 0) stop("z2_cap must be positive.")
  if (!is.finite(cfg$equality_tolerance) || cfg$equality_tolerance < 0) stop("equality_tolerance must be non-negative.")
  if (!is.finite(cfg$random_seed)) stop("random_seed must be finite.")
  if (!is.finite(cfg$min_valid_frac) || cfg$min_valid_frac <= 0 || cfg$min_valid_frac > 1) {
    stop("min_valid_frac must be in (0, 1].")
  }
  if (length(cfg$ridge_sequence) < 1 || any(!is.finite(cfg$ridge_sequence)) ||
      any(cfg$ridge_sequence < 0) || is.unsorted(cfg$ridge_sequence)) {
    stop("ridge_sequence must be a non-decreasing vector of finite non-negative values.")
  }
  # This tutorial is fixed to the 19-predictor (Bio1..Bio19) canonical profile.
  expected_focal <- paste0("Bio", 1:19)
  if (!is.character(cfg$focal_vars) || length(cfg$focal_vars) != 19L ||
      !identical(sort(cfg$focal_vars), sort(expected_focal))) {
    stop("focal_vars must be the canonical 19-predictor set: Bio1..Bio19.")
  }
  # mask_shp is optional in this tutorial: a NULL / empty value skips masking.
  if (!is.null(cfg$mask_shp) && nzchar(cfg$mask_shp) && !file.exists(cfg$mask_shp)) {
    stop("Mask shapefile was declared but not found: ", cfg$mask_shp)
  }
  if (!file.exists(cfg$occ_csv)) {
    stop("Occurrence CSV was not found: ", cfg$occ_csv)
  }
  if (!dir.exists(cfg$current_vars_dir)) {
    stop("Predictor raster directory was not found: ", cfg$current_vars_dir)
  }
  invisible(TRUE)
}

validate_config(cfg)
set.seed(as.integer(cfg$random_seed))

# =============================================================================
# SECTION 2 -- UTILITY HELPERS
# =============================================================================

csv_out <- function(out_paths, stem) {
  file.path(out_paths$csvs, paste0(stem, ".csv"))
}

pick_col <- function(nms, candidates, label) {
  hit <- candidates[candidates %in% nms][1]
  if (is.na(hit)) stop(sprintf("Could not find %s column in occurrence CSV.", label))
  hit
}

list_tifs_named <- function(path) {
  files <- sort(list.files(path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE))
  if (!length(files)) stop(sprintf("No .tif files found in %s", path))
  tibble(variable = tools::file_path_sans_ext(basename(files)), file = files)
}

validate_raster_stack <- function(r_stack, expected_vars = NULL) {
  if (terra::nlyr(r_stack) < 1) stop("Raster stack contains no layers.")
  if (!is.null(expected_vars)) {
    missing <- setdiff(expected_vars, names(r_stack))
    if (length(missing)) stop("Missing required raster predictors: ", paste(missing, collapse = ", "))
  }
  if (!nzchar(terra::crs(r_stack))) stop("Raster stack has no coordinate reference system.")
  if (terra::nlyr(r_stack) > 1) {
    ref <- r_stack[[1]]
    for (j in 2:terra::nlyr(r_stack)) {
      ok <- terra::compareGeom(ref, r_stack[[j]], crs = TRUE, ext = TRUE,
                               rowcol = TRUE, res = TRUE,
                               stopOnError = FALSE)
      if (!isTRUE(ok)) stop("Raster geometry mismatch detected for layer: ", names(r_stack)[j])
    }
  }
  invisible(TRUE)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x)
}

safe_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(0)
  mean(((x - mean(x)) / s)^3)
}

safe_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(0)
  mean(((x - mean(x)) / s)^4) - 3
}

clamp01 <- function(x, fallback = 0) {
  if (!is.finite(x)) return(fallback)
  max(0, min(1, x))
}

write_lookup_outputs <- function(vars, out_paths) {
  lookup <- tibble(
    variable       = vars,
    stressor_index = seq_along(vars)
  )
  out_file <- file.path(out_paths$csvs,
                        paste0(cfg$species_code, "_00_predictor_lookup_dictionary.csv"))
  write.csv(lookup, out_file, row.names = FALSE)
  invisible(lookup)
}

# =============================================================================
# SECTION 3 -- PER-VARIABLE STRESS REFERENCE
# =============================================================================

calc_shrinkage_weight <- function(n_tail, k = 10, method = "n_over_n_plus_k") {
  if (!is.finite(n_tail) || n_tail <= 0) return(0)
  if (!is.finite(k) || k <= 0) return(1)
  if (identical(method, "n_over_n_plus_k"))           return(clamp01(n_tail / (n_tail + k)))
  if (identical(method, "sqrt_n_over_sqrt_n_plus_k")) return(clamp01(sqrt(n_tail) / (sqrt(n_tail) + sqrt(k))))
  if (identical(method, "capped_linear"))             return(clamp01(n_tail / k))
  stop(sprintf("Unknown shrinkage_method: %s", method))
}

calc_shrunk_tail_sd <- function(raw_sd, sigma_global, n_tail, k = 10, method = "n_over_n_plus_k") {
  if (!is.finite(sigma_global) || sigma_global <= 0) sigma_global <- 1e-6
  if (!is.finite(raw_sd) || raw_sd <= 0) {
    return(list(sd = sigma_global, weight = 0, used_global_only = TRUE))
  }
  w <- calc_shrinkage_weight(n_tail, k = k, method = method)
  shrunk_sd <- w * raw_sd + (1 - w) * sigma_global
  list(sd = max(shrunk_sd, 1e-6), weight = w, used_global_only = FALSE)
}

calc_tail_scale <- function(values, common_mu, center_method = "global_mu_rms") {
  values <- values[is.finite(values)]
  if (length(values) < 2L || !is.finite(common_mu)) return(NA_real_)

  if (identical(center_method, "tail_mean_sd")) {
    return(safe_sd(values))
  }
  if (identical(center_method, "global_mu_rms")) {
    # Population RMS distance from the same global occurrence mean used by the
    # VRS numerator. This is a tail scale, not a standard error or sample SD.
    return(sqrt(mean((values - common_mu)^2)))
  }
  stop("Unknown tail-scale centering method: ", center_method)
}

bootstrap_asymmetry_diagnostic <- function(x, iter = 100, frac = 1.0,
                                           ratio_threshold = 1.2, min_tail_n = 5,
                                           tail_scale_center = "global_mu_rms") {
  x <- x[is.finite(x)]
  n_total <- length(x)

  iter <- as.integer(iter)
  if (!is.finite(iter) || iter <= 0) iter <- 100L
  if (!is.finite(frac) || frac <= 0) frac <- 0.8
  frac <- min(frac, 1)
  if (!is.finite(ratio_threshold) || ratio_threshold <= 1) ratio_threshold <- 1.2
  min_tail_n <- as.integer(min_tail_n)
  if (!is.finite(min_tail_n) || min_tail_n < 2) min_tail_n <- 5L

  if (n_total < max(10, min_tail_n * 2)) {
    return(tibble(
      bootstrap_n = 0L, asym_ratio_raw = NA_real_, asym_ratio_boot_mean = NA_real_,
      asym_ratio_ci_low = NA_real_, asym_ratio_ci_high = NA_real_,
      asymmetry_bootstrap_stable = FALSE, asymmetry_bootstrap_note = "insufficient_total_n"
    ))
  }

  sample_n <- max(2L * min_tail_n, floor(n_total * frac))
  ratios <- rep(NA_real_, iter)

  for (i in seq_len(iter)) {
    xs <- sample(x, size = sample_n, replace = TRUE)
    mu_boot <- mean(xs)
    lower_vals <- xs[xs < mu_boot]
    upper_vals <- xs[xs > mu_boot]
    lower_sd   <- calc_tail_scale(lower_vals, mu_boot, tail_scale_center)
    upper_sd   <- calc_tail_scale(upper_vals, mu_boot, tail_scale_center)
    if (length(lower_vals) >= min_tail_n && length(upper_vals) >= min_tail_n &&
        is.finite(lower_sd) && lower_sd > 0 &&
        is.finite(upper_sd) && upper_sd > 0) {
      ratios[i] <- upper_sd / lower_sd
    }
  }

  ratios <- ratios[is.finite(ratios)]
  mu_raw <- mean(x)
  lower_sd_raw <- calc_tail_scale(x[x < mu_raw], mu_raw, tail_scale_center)
  upper_sd_raw <- calc_tail_scale(x[x > mu_raw], mu_raw, tail_scale_center)
  raw_ratio <- if (is.finite(lower_sd_raw) && lower_sd_raw > 0 &&
                   is.finite(upper_sd_raw) && upper_sd_raw > 0) {
    upper_sd_raw / lower_sd_raw
  } else NA_real_

  if (!length(ratios)) {
    return(tibble(
      bootstrap_n = 0L, asym_ratio_raw = raw_ratio, asym_ratio_boot_mean = NA_real_,
      asym_ratio_ci_low = NA_real_, asym_ratio_ci_high = NA_real_,
      asymmetry_bootstrap_stable = FALSE, asymmetry_bootstrap_note = "insufficient_tail_bootstrap"
    ))
  }

  ci <- as.numeric(stats::quantile(ratios, probs = c(0.025, 0.975),
                                   na.rm = TRUE, names = FALSE, type = 7))
  stable <- (ci[1] > ratio_threshold) || (ci[2] < 1 / ratio_threshold)

  tibble(
    bootstrap_n                = length(ratios),
    asym_ratio_raw             = raw_ratio,
    asym_ratio_boot_mean       = mean(ratios, na.rm = TRUE),
    asym_ratio_ci_low          = ci[1],
    asym_ratio_ci_high         = ci[2],
    asymmetry_bootstrap_stable = stable,
    asymmetry_bootstrap_note   = ifelse(stable, "stable_asymmetry", "unstable_or_near_symmetric")
  )
}

background_truncation_diagnostic <- function(layer, x, prop_threshold = 0.05) {
  x  <- x[is.finite(x)]

  if (!is.finite(prop_threshold) || prop_threshold < 0) prop_threshold <- 0.05
  prop_threshold <- min(prop_threshold, 1)

  if (!length(x)) {
    return(tibble(
      bg_lower_tail_prop    = NA_real_, bg_upper_tail_prop    = NA_real_,
      lower_truncation_flag = FALSE,    upper_truncation_flag = FALSE,
      truncation_flag_any   = FALSE
    ))
  }

  obs_min    <- min(x); obs_max <- max(x)
  lower_prop <- as.numeric(terra::global(layer < obs_min, "mean", na.rm = TRUE)[1, 1])
  upper_prop <- as.numeric(terra::global(layer > obs_max, "mean", na.rm = TRUE)[1, 1])
  lower_flag <- is.finite(lower_prop) && lower_prop <= prop_threshold
  upper_flag <- is.finite(upper_prop) && upper_prop <= prop_threshold

  tibble(
    bg_lower_tail_prop    = lower_prop, bg_upper_tail_prop    = upper_prop,
    lower_truncation_flag = lower_flag, upper_truncation_flag = upper_flag,
    truncation_flag_any   = lower_flag || upper_flag
  )
}

# =============================================================================
# SECTION 4 -- I/O SETUP
# =============================================================================

build_output_paths <- function(base_dir, temp_name = "temp") {
  paths <- list(
    base                 = base_dir,
    rasters              = file.path(base_dir, "rasters"),
    raster_core_occupied = file.path(base_dir, "rasters", "core_occupied"),
    raster_optional      = file.path(base_dir, "rasters", "optional_diagnostics"),
    csvs                 = file.path(base_dir, "csvs"),
    addon_rasters        = file.path(base_dir, "addons", "rasters"),
    addon_csvs           = file.path(base_dir, "addons", "csvs"),
    temp                 = file.path(base_dir, temp_name)
  )
  lapply(paths, function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE))
  paths
}

build_stats_from_occ <- function(occ, current_stack, lon_col, lat_col, min_tail_n) {
  occ[[lon_col]] <- suppressWarnings(as.numeric(occ[[lon_col]]))
  occ[[lat_col]] <- suppressWarnings(as.numeric(occ[[lat_col]]))
  bad_xy <- !is.finite(occ[[lon_col]]) | !is.finite(occ[[lat_col]]) |
    occ[[lon_col]] < -180 | occ[[lon_col]] > 180 |
    occ[[lat_col]] < -90 | occ[[lat_col]] > 90
  if (any(bad_xy)) {
    stop(sprintf("Occurrence table contains %d invalid longitude/latitude rows; clean them before VERA.",
                 sum(bad_xy)))
  }
  pts <- vect(occ, geom = c(lon_col, lat_col), crs = "EPSG:4326")
  if (!terra::same.crs(pts, current_stack)) pts <- terra::project(pts, terra::crs(current_stack))
  vals    <- terra::extract(current_stack, pts, ID = FALSE)
  occ_env <- bind_cols(occ, as.data.frame(vals))
  complete_env <- complete.cases(occ_env[, names(current_stack), drop = FALSE])
  excluded_incomplete <- sum(!complete_env)
  if (excluded_incomplete > 0) {
    message(sprintf("    Removing %d occurrence rows lacking a complete finite predictor vector.",
                    excluded_incomplete))
  }
  occ_env <- occ_env[complete_env, , drop = FALSE]
  if (!nrow(occ_env)) stop("No occurrence rows retain complete finite environmental values.")

  stats_tbl <- lapply(names(current_stack), function(v) {
    x  <- occ_env[[v]]
    x  <- x[is.finite(x)]
    mu <- mean(x)
    sigma_global <- safe_sd(x)
    lower_vals   <- x[x < mu]
    upper_vals   <- x[x > mu]
    raw_tail_scale_lower <- calc_tail_scale(lower_vals, mu, cfg$tail_scale_center)
    raw_tail_scale_upper <- calc_tail_scale(upper_vals, mu, cfg$tail_scale_center)

    fallback_lower <- length(lower_vals) < min_tail_n ||
      !is.finite(raw_tail_scale_lower) || raw_tail_scale_lower <= 0
    fallback_upper <- length(upper_vals) < min_tail_n ||
      !is.finite(raw_tail_scale_upper) || raw_tail_scale_upper <= 0
    shrink_lower <- if (fallback_lower) {
      list(sd = sigma_global, weight = 0, used_global_only = TRUE)
    } else {
      calc_shrunk_tail_sd(raw_tail_scale_lower, sigma_global, length(lower_vals),
                          k = cfg$shrinkage_k, method = cfg$shrinkage_method)
    }
    shrink_upper <- if (fallback_upper) {
      list(sd = sigma_global, weight = 0, used_global_only = TRUE)
    } else {
      calc_shrunk_tail_sd(raw_tail_scale_upper, sigma_global, length(upper_vals),
                          k = cfg$shrinkage_k, method = cfg$shrinkage_method)
    }
    sd_lower <- shrink_lower$sd
    sd_upper <- shrink_upper$sd

    boot_tbl <- if (isTRUE(cfg$bootstrap_asymmetry)) {
      bootstrap_asymmetry_diagnostic(
        x,
        iter            = cfg$bootstrap_iter,
        frac            = cfg$bootstrap_frac,
        ratio_threshold = cfg$asymmetry_ratio_threshold,
        min_tail_n      = min_tail_n,
        tail_scale_center = cfg$tail_scale_center
      )
    } else {
      tibble(
        bootstrap_n = 0L, asym_ratio_raw = NA_real_, asym_ratio_boot_mean = NA_real_,
        asym_ratio_ci_low = NA_real_, asym_ratio_ci_high = NA_real_,
        asymmetry_bootstrap_stable = FALSE, asymmetry_bootstrap_note = "bootstrap_disabled"
      )
    }

    trunc_tbl <- background_truncation_diagnostic(
      current_stack[[v]], x,
      prop_threshold = cfg$truncation_bg_tail_prop_threshold
    )

    tibble(
      variable     = v, mu = mu, sigma_global = sigma_global,
      n            = length(x), n_lower = length(lower_vals), n_upper = length(upper_vals),
      tail_scale_center = cfg$tail_scale_center,
      raw_tail_scale_lower = ifelse(is.finite(raw_tail_scale_lower), raw_tail_scale_lower, NA_real_),
      raw_tail_scale_upper = ifelse(is.finite(raw_tail_scale_upper), raw_tail_scale_upper, NA_real_),
      # Compatibility aliases retained for existing readers. Under
      # global_mu_rms these columns contain tail scales, not sample SDs.
      raw_sd_lower = ifelse(is.finite(raw_tail_scale_lower), raw_tail_scale_lower, NA_real_),
      raw_sd_upper = ifelse(is.finite(raw_tail_scale_upper), raw_tail_scale_upper, NA_real_),
      shrink_weight_lower = shrink_lower$weight,
      shrink_weight_upper = shrink_upper$weight,
      sd_lower = sd_lower, sd_upper = sd_upper,
      fallback_lower = fallback_lower, fallback_upper = fallback_upper
    ) %>% bind_cols(boot_tbl, trunc_tbl)
  }) %>% bind_rows()

  list(stats = stats_tbl, occ_env = occ_env, excluded_incomplete = excluded_incomplete)
}

# =============================================================================
# SECTION 5 -- MULTIVARIATE CLIMATIC CORE SELECTION
# =============================================================================

calc_mahal_core_occupied <- function(occ_env, vars_common, core_occupied_prob) {
  if (!length(vars_common)) stop("No common variables provided for Mahalanobis core selection.")
  missing_vars <- setdiff(vars_common, names(occ_env))
  if (length(missing_vars)) {
    stop("Missing variables in occurrence environment table: ",
         paste(missing_vars, collapse = ", "))
  }
  p     <- length(vars_common)
  valid <- complete.cases(occ_env[, vars_common, drop = FALSE])
  min_valid <- p + 1L
  if (sum(valid) < min_valid) {
    stop(sprintf("Not enough valid records. Valid: %d, required: %d, variables: %d",
                 sum(valid), min_valid, p))
  }

  env_valid <- occ_env[valid, vars_common, drop = FALSE]
  vars_sd   <- apply(env_valid, 2, stats::sd)
  bad_sd    <- !is.finite(vars_sd) | vars_sd <= 0
  if (any(bad_sd)) {
    stop("Zero-variance or invalid variables: ", paste(vars_common[bad_sd], collapse = ", "))
  }
  if (!is.finite(core_occupied_prob) || core_occupied_prob <= 0 || core_occupied_prob > 1) {
    stop(sprintf("Invalid core_occupied_prob: %s", core_occupied_prob))
  }

  env_scaled <- scale(env_valid)
  cov_mat    <- stats::cov(env_scaled)
  diag_scale <- mean(diag(cov_mat), na.rm = TRUE)
  if (!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1

  ridge_seq      <- cfg$ridge_sequence
  mahal_d2       <- NULL
  mahal_distance <- NULL
  ridge_used     <- NA_real_
  last_error     <- NULL

  for (ridge in ridge_seq) {
    cov_try <- cov_mat + diag(diag_scale * ridge, ncol(cov_mat))
    mahal_try_d2 <- tryCatch(
      stats::mahalanobis(env_scaled, center = rep(0, ncol(env_scaled)), cov = cov_try),
      error = function(e) e
    )
    if (is.numeric(mahal_try_d2) &&
        length(mahal_try_d2) == nrow(env_scaled) &&
        all(is.finite(mahal_try_d2))) {
      mahal_d2       <- mahal_try_d2
      mahal_distance <- sqrt(mahal_d2)
      ridge_used     <- ridge
      break
    } else {
      last_error <- mahal_try_d2
    }
  }

  if (is.null(mahal_d2)) {
    message("\n MAHALANOBIS DIAGNOSTIC ")
    if (inherits(last_error, "error")) message("R original error: ", last_error$message)
    stop("Mahalanobis selection failed even after covariance regularization.")
  }

  occ_with_dist <- occ_env[valid, , drop = FALSE]
  occ_with_dist$mahal_d2       <- mahal_d2
  occ_with_dist$mahal_distance <- mahal_distance
  occ_with_dist$mahal_ridge    <- ridge_used
  q_threshold_dist <- as.numeric(stats::quantile(
    occ_with_dist$mahal_distance, probs = core_occupied_prob,
    na.rm = TRUE, names = FALSE, type = 7
  ))
  if (!is.finite(q_threshold_dist)) stop("Mahalanobis quantile threshold is non-finite.")
  core_occupied_idx        <- which(occ_with_dist$mahal_distance <= q_threshold_dist)
  if (!length(core_occupied_idx)) core_occupied_idx <- which.min(occ_with_dist$mahal_distance)
  selection_value          <- core_occupied_prob
  selection_threshold_dist <- q_threshold_dist
  selection_threshold_d2   <- q_threshold_dist^2

  if (!is.finite(selection_threshold_dist) || !is.finite(selection_threshold_d2)) {
    stop("Mahalanobis core threshold is non-finite.")
  }

  core_occupied_idx <- sort(unique(core_occupied_idx))
  peripheral_idx    <- setdiff(seq_len(nrow(occ_with_dist)), core_occupied_idx)

  occ_with_dist$core_rule                <- "mahal_quantile"
  occ_with_dist$core_rule_value          <- selection_value
  occ_with_dist$core_threshold_distance  <- selection_threshold_dist
  occ_with_dist$core_threshold_d2        <- selection_threshold_d2

  list(
    core_occupied      = occ_with_dist[core_occupied_idx, , drop = FALSE],
    peripheral         = if (length(peripheral_idx) > 0) {
      occ_with_dist[peripheral_idx, , drop = FALSE]
    } else {
      occ_with_dist[integer(0), , drop = FALSE]
    },
    threshold_distance = selection_threshold_dist,
    threshold_d2       = selection_threshold_d2,
    rule               = "mahal_quantile",
    rule_value         = selection_value
  )
}

build_mahal_reference <- function(occ_set, vars, label = "mahal_reference") {
  if (!length(vars)) stop(sprintf("%s: no variables provided.", label))
  missing_vars <- setdiff(vars, names(occ_set))
  if (length(missing_vars)) {
    stop(sprintf("%s: missing variables: %s", label, paste(missing_vars, collapse = ", ")))
  }

  env   <- occ_set[, vars, drop = FALSE]
  valid <- complete.cases(env)
  env   <- env[valid, , drop = FALSE]

  min_rows <- max(3L, length(vars) + 1L)
  if (nrow(env) < min_rows) {
    stop(sprintf("Too few complete rows (%d; required %d) to build %s.",
                 nrow(env), min_rows, label))
  }

  vars_sd <- apply(env, 2, stats::sd)
  bad_sd  <- !is.finite(vars_sd) | vars_sd <= 0
  if (any(bad_sd)) {
    stop(sprintf("%s: zero-variance or invalid variables: %s",
                 label, paste(vars[bad_sd], collapse = ", ")))
  }

  center     <- colMeans(env)
  scale_vec  <- apply(env, 2, stats::sd)
  env_scaled <- sweep(sweep(as.matrix(env), 2, center, "-"), 2, scale_vec, "/")
  cov_mat    <- stats::cov(env_scaled)
  diag_scale <- mean(diag(cov_mat), na.rm = TRUE)
  if (!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1

  ridge_seq  <- cfg$ridge_sequence
  cov_used   <- NULL
  ridge_used <- NA_real_

  for (ridge in ridge_seq) {
    cov_try <- cov_mat + diag(diag_scale * ridge, ncol(cov_mat))
    d2_try  <- tryCatch(
      stats::mahalanobis(env_scaled, center = rep(0, length(vars)), cov = cov_try),
      error = function(e) NULL
    )
    if (is.numeric(d2_try) &&
        length(d2_try) == nrow(env) &&
        all(is.finite(d2_try))) {
      cov_used   <- cov_try
      ridge_used <- ridge
      break
    }
  }

  if (is.null(cov_used)) {
    stop(sprintf("%s: covariance regularization failed for all ridge values.", label))
  }

  list(
    n      = nrow(env),
    center = center,
    scale  = scale_vec,
    cov    = cov_used,
    ridge  = ridge_used
  )
}

# =============================================================================
# SECTION 6 -- STRESS SURFACES AND TIER CLASSIFICATION
# =============================================================================

calc_mahal_distance_surface <- function(r_stack, core_occupied_occ_or_ref, vars,
                                        scenario_label = "current") {
  if (!length(vars)) stop("No variables provided for Mahalanobis distance surface.")
  missing_vars <- setdiff(vars, names(r_stack))
  if (length(missing_vars)) {
    stop("Missing variables in raster stack: ", paste(missing_vars, collapse = ", "))
  }

  if (is.list(core_occupied_occ_or_ref) &&
      all(c("center", "scale", "cov", "ridge") %in% names(core_occupied_occ_or_ref))) {
    center     <- core_occupied_occ_or_ref$center
    scale_vec  <- core_occupied_occ_or_ref$scale
    cov_used   <- core_occupied_occ_or_ref$cov
    ridge_used <- core_occupied_occ_or_ref$ridge
  } else {
    ref        <- build_mahal_reference(core_occupied_occ_or_ref, vars,
                                        label = paste0("surface_", scenario_label))
    center     <- ref$center
    scale_vec  <- ref$scale
    cov_used   <- ref$cov
    ridge_used <- ref$ridge
  }

  mahal_fun <- function(vals) {
    if (is.null(dim(vals))) vals <- matrix(vals, nrow = 1)
    out <- matrix(NA_real_, nrow = nrow(vals), ncol = 2)
    ok <- complete.cases(vals)
    if (any(ok)) {
      z <- sweep(sweep(vals[ok, , drop = FALSE], 2, center, "-"), 2, scale_vec, "/")
      d2 <- stats::mahalanobis(z, center = rep(0, length(vars)), cov = cov_used)
      out[ok, 1] <- d2
      out[ok, 2] <- sqrt(d2)
    }
    out
  }
  mahal_stack <- terra::app(r_stack[[vars]], fun = mahal_fun)
  names(mahal_stack) <- c(paste0(scenario_label, "_mahal_d2"),
                          paste0(scenario_label, "_mahal_distance"))
  d2_rast <- mahal_stack[[1]]
  d_rast  <- mahal_stack[[2]]

  list(distance = d_rast, d2 = d2_rast, center = center, scale = scale_vec,
       cov = cov_used, ridge = ridge_used)
}

resolve_tier_breaks <- function(method, core_occupied_distances, quantiles, fixed_breaks, n_vars) {
  if (identical(method, "empirical_quantile")) {
    if (!length(core_occupied_distances) || !any(is.finite(core_occupied_distances))) {
      stop("Cannot resolve empirical tier breaks: no finite core_occupied distances provided.")
    }
    breaks <- as.numeric(stats::quantile(
      core_occupied_distances, probs = quantiles, na.rm = TRUE, names = FALSE, type = 7
    ))
    if (any(!is.finite(breaks)) || any(diff(breaks) <= 0)) {
      stop(sprintf(
        "Empirical tier breaks not strictly increasing (got %s); inspect core_occupied distance distribution.",
        paste(round(breaks, 4), collapse = " | ")
      ))
    }
    return(breaks)
  }

  if (identical(method, "fixed")) {
    breaks <- as.numeric(fixed_breaks)
    if (length(breaks) != 3 || any(!is.finite(breaks)) || any(diff(breaks) <= 0)) {
      stop("Fixed Mahalanobis breaks must be three finite strictly increasing values.")
    }
    return(breaks)
  }

  if (identical(method, "chi_square")) {
    if (!is.finite(n_vars) || n_vars <= 0) {
      stop("Invalid n_vars for chi_square tier breaks.")
    }
    if (any(!is.finite(quantiles)) || any(quantiles <= 0 | quantiles >= 1)) {
      stop("Chi-square tier quantiles must be finite values in (0, 1).")
    }
    return(sqrt(stats::qchisq(quantiles, df = n_vars)))
  }

  stop(sprintf("Unknown tier_breaks_method: %s", method))
}

classify_mahal_tiers <- function(mahal_distance_rast, breaks = c(1, 2, 3)) {
  rcl <- matrix(c(
    -Inf,      breaks[1], 1,
    breaks[1], breaks[2], 2,
    breaks[2], breaks[3], 3,
    breaks[3], Inf,       4
  ), ncol = 3, byrow = TRUE)
  classify(mahal_distance_rast, rcl)
}

write_factor_raster <- function(stressor_rast, lookup_tbl, out_file) {
  if (is.null(stressor_rast)) return(invisible(NULL))
  out <- as.factor(stressor_rast)
  level_df <- data.frame(value = lookup_tbl$stressor_index, label = lookup_tbl$variable)
  levels(out) <- level_df
  writeRaster(out, out_file, overwrite = TRUE, wopt = gdal_opts)
  invisible(out_file)
}

calc_asym_surfaces <- function(r_stack, stats_tbl, scenario_label) {
  cat(sprintf("      Computing Z-squared (VRS) surfaces for %s (%d vars)...\n",
              scenario_label, nrow(stats_tbl)))

  sum_r_sq       <- terra::init(r_stack[[1]], 0)
  valid_count    <- terra::init(r_stack[[1]], 0)

  primary_z_sq   <- terra::init(r_stack[[1]], NA)
  primary_idx    <- terra::init(r_stack[[1]], NA)
  primary_shared_max_count <- terra::init(r_stack[[1]], 0)
  secondary_z_sq <- terra::init(r_stack[[1]], NA)
  secondary_idx  <- terra::init(r_stack[[1]], NA)

  too_high_z_sq  <- terra::init(r_stack[[1]], NA)
  too_high_idx   <- terra::init(r_stack[[1]], NA)
  too_high_shared_max_count <- terra::init(r_stack[[1]], 0)
  too_low_z_sq   <- terra::init(r_stack[[1]], NA)
  too_low_idx    <- terra::init(r_stack[[1]], NA)
  too_low_shared_max_count <- terra::init(r_stack[[1]], 0)

  sum_z_sq_hi    <- terra::init(r_stack[[1]], 0)
  sum_z_sq_lo    <- terra::init(r_stack[[1]], 0)

  MAX_VRS_SCORE <- cfg$z2_cap
  equality_tol <- cfg$equality_tolerance
  predictor_diag_rows <- vector("list", nrow(stats_tbl))

  for (i in seq_len(nrow(stats_tbl))) {
    v     <- stats_tbl$variable[i]
    mu    <- stats_tbl$mu[i]
    sd_lo <- stats_tbl$sd_lower[i]
    sd_hi <- stats_tbl$sd_upper[i]

    if (!is.finite(sd_lo) || sd_lo <= 0 || !is.finite(sd_hi) || sd_hi <= 0) {
      warning(sprintf("Skipping variable with invalid asymmetric SD: %s", v))
      next
    }

    x    <- r_stack[[v]]
    z_sq <- ifel(x < mu, ((x - mu) / sd_lo)^2, ((x - mu) / sd_hi)^2)
    z_sq <- ifel(z_sq > MAX_VRS_SCORE, MAX_VRS_SCORE, z_sq)

    # Reporting-only metrics are streamed and immediately reduced to scalars;
    # no additional predictor raster is retained or written.
    if (isTRUE(cfg$write_reporting_audit)) {
      diagnostic_pair <- c(
        z_sq,
        z_sq >= (MAX_VRS_SCORE - equality_tol)
      )
      diagnostic_global <- suppressWarnings(
        terra::global(diagnostic_pair, "mean", na.rm = TRUE)
      )
      predictor_diag_rows[[i]] <- tibble(
        variable = v,
        mean_vrs_contribution = as.numeric(diagnostic_global[1, 1]),
        cap_saturation_fraction = as.numeric(diagnostic_global[2, 1])
      )
      rm(diagnostic_pair, diagnostic_global)
    }

    z_sq_hi_only <- ifel(x > mu, z_sq, 0)
    z_sq_lo_only <- ifel(x < mu, z_sq, 0)

    valid_z     <- !is.na(z_sq)
    sum_r_sq    <- sum_r_sq + ifel(valid_z, z_sq, 0)
    valid_count <- valid_count + ifel(valid_z, 1, 0)

    sum_z_sq_hi <- sum_z_sq_hi + ifel(valid_z, z_sq_hi_only, 0)
    sum_z_sq_lo <- sum_z_sq_lo + ifel(valid_z, z_sq_lo_only, 0)

    had_primary      <- !is.na(primary_z_sq)
    is_new_primary   <- valid_z & (!had_primary | z_sq > primary_z_sq + equality_tol)
    is_codominant_primary <- valid_z & had_primary & abs(z_sq - primary_z_sq) <= equality_tol
    is_codominant_secondary <- is_codominant_primary &
      (is.na(secondary_z_sq) | primary_z_sq > secondary_z_sq)
    is_new_secondary <- valid_z & !is_new_primary & !is_codominant_primary &
      (is.na(secondary_z_sq) | z_sq > secondary_z_sq)

    new_secondary_z_sq <- ifel(is_new_primary, primary_z_sq,
                               ifel(is_codominant_secondary, primary_z_sq,
                                    ifel(is_new_secondary, z_sq, secondary_z_sq)))
    new_secondary_idx  <- ifel(is_new_primary, primary_idx,
                               ifel(is_codominant_secondary, i,
                                    ifel(is_new_secondary, i, secondary_idx)))

    primary_shared_max_count <- ifel(
      is_new_primary, 1,
      ifel(is_codominant_primary, primary_shared_max_count + 1, primary_shared_max_count)
    )
    primary_z_sq   <- ifel(is_new_primary, z_sq, primary_z_sq)
    primary_idx    <- ifel(is_new_primary, i,    primary_idx)
    secondary_z_sq <- new_secondary_z_sq
    secondary_idx  <- new_secondary_idx

    had_high      <- !is.na(too_high_z_sq)
    is_new_high   <- valid_z & (z_sq_hi_only > 0) &
      (!had_high | z_sq_hi_only > too_high_z_sq + equality_tol)
    is_codominant_high <- valid_z & (z_sq_hi_only > 0) & had_high &
      abs(z_sq_hi_only - too_high_z_sq) <= equality_tol
    too_high_shared_max_count <- ifel(
      is_new_high, 1,
      ifel(is_codominant_high, too_high_shared_max_count + 1, too_high_shared_max_count)
    )
    too_high_z_sq <- ifel(is_new_high, z_sq_hi_only, too_high_z_sq)
    too_high_idx  <- ifel(is_new_high, i,            too_high_idx)

    had_low      <- !is.na(too_low_z_sq)
    is_new_low   <- valid_z & (z_sq_lo_only > 0) &
      (!had_low | z_sq_lo_only > too_low_z_sq + equality_tol)
    is_codominant_low <- valid_z & (z_sq_lo_only > 0) & had_low &
      abs(z_sq_lo_only - too_low_z_sq) <= equality_tol
    too_low_shared_max_count <- ifel(
      is_new_low, 1,
      ifel(is_codominant_low, too_low_shared_max_count + 1, too_low_shared_max_count)
    )
    too_low_z_sq <- ifel(is_new_low, z_sq_lo_only, too_low_z_sq)
    too_low_idx  <- ifel(is_new_low, i,            too_low_idx)

    # --- LOOP INTRA-STEP CLEANUP ---
    rm(x, z_sq, z_sq_hi_only, z_sq_lo_only, valid_z, had_primary, had_high, had_low,
       is_new_primary, is_codominant_primary, is_codominant_secondary, is_new_secondary,
       is_new_high, is_codominant_high, is_new_low, is_codominant_low)

    if (is.finite(cfg$gc_interval) && cfg$gc_interval > 0 && i %% cfg$gc_interval == 0) {
      gc(verbose = FALSE)
    }
  }

  # --- POST-LOOP CLEANUP ---
  gc(verbose = FALSE)

  min_valid <- max(1L, ceiling(nrow(stats_tbl) * cfg$min_valid_frac))
  pass_mask <- valid_count >= min_valid

  pixels_censored <- as.integer(terra::global((valid_count > 0) & (valid_count < min_valid),
                                              "sum", na.rm = TRUE)[1, 1])
  if (!is.finite(pixels_censored)) pixels_censored <- 0L

  mean_vrs_raw    <- sum_r_sq / valid_count
  mean_vrs        <- ifel(pass_mask, mean_vrs_raw, NA)
  names(mean_vrs) <- "mean_vrs"

  max_vrs        <- ifel(pass_mask, primary_z_sq, NA)
  names(max_vrs) <- "max_vrs"

  delta_vrs        <- ifel(pass_mask, primary_z_sq - secondary_z_sq, NA)
  names(delta_vrs) <- "delta_vrs"

  primary_co_dominance_count <- ifel(pass_mask, primary_shared_max_count, NA)
  names(primary_co_dominance_count) <- "primary_co_dominance_count"
  primary_unique_index <- ifel(pass_mask & primary_shared_max_count == 1, primary_idx, NA)
  names(primary_unique_index) <- "primary_unique_stressor_index"
  too_high_co_dominance_count <- ifel(pass_mask, too_high_shared_max_count, NA)
  names(too_high_co_dominance_count) <- "too_high_co_dominance_count"
  too_low_co_dominance_count <- ifel(pass_mask, too_low_shared_max_count, NA)
  names(too_low_co_dominance_count) <- "too_low_co_dominance_count"

  sum_z_hi_final <- ifel(pass_mask, sum_z_sq_hi, NA)
  sum_z_lo_final <- ifel(pass_mask, sum_z_sq_lo, NA)
  names(sum_z_hi_final) <- "too_high_zsq_sum"
  names(sum_z_lo_final) <- "too_low_zsq_sum"

  predictor_diagnostics <- bind_rows(predictor_diag_rows)
  if (nrow(predictor_diagnostics)) {
    total_contribution <- sum(predictor_diagnostics$mean_vrs_contribution,
                              na.rm = TRUE)
    predictor_diagnostics <- predictor_diagnostics %>%
      mutate(
        mean_vrs_contribution_pct = if (is.finite(total_contribution) &&
                                            total_contribution > 0) {
          100 * .data$mean_vrs_contribution / total_contribution
        } else NA_real_
      ) %>%
      arrange(desc(.data$mean_vrs_contribution_pct)) %>%
      mutate(
        cumulative_vrs_contribution_pct = cumsum(
          ifelse(is.na(.data$mean_vrs_contribution_pct),
                 0, .data$mean_vrs_contribution_pct)
        ),
        contribution_coverage_member =
          lag(.data$cumulative_vrs_contribution_pct, default = 0) <
          (100 * cfg$reporting_contribution_coverage)
      ) %>%
      arrange(match(.data$variable, stats_tbl$variable))
  }

  # --- FINAL FUNCTION CLEANUP BEFORE RETURN ---
  gc(verbose = FALSE)

  list(
    mean_vrs        = mean_vrs,
    max_vrs         = max_vrs,
    delta_vrs       = delta_vrs,
    stressor_index  = ifel(pass_mask, primary_idx, NA),
    primary_unique_index = primary_unique_index,
    primary_co_dominance_count = primary_co_dominance_count,
    secondary_index = ifel(pass_mask, secondary_idx, NA),
    too_high_index  = ifel(pass_mask, too_high_idx, NA),
    too_low_index   = ifel(pass_mask, too_low_idx, NA),
    too_high_z      = ifel(pass_mask, ifel(is.na(too_high_z_sq), 0, too_high_z_sq), NA),
    too_low_z       = ifel(pass_mask, ifel(is.na(too_low_z_sq), 0, too_low_z_sq), NA),
    too_high_co_dominance_count = too_high_co_dominance_count,
    too_low_co_dominance_count = too_low_co_dominance_count,
    predictor_diagnostics = predictor_diagnostics,
    too_high_z_sum  = sum_z_hi_final,
    too_low_z_sum   = sum_z_lo_final,
    valid_count     = valid_count,
    pixels_censored = pixels_censored
  )
}

# =============================================================================
# SECTION 7 -- OCCURRENCE CLASS FATES AND REFERENCE BOUNDS
# =============================================================================

extract_class_fate <- function(class_rast, pts, scenario_name, class_labels) {
  vals  <- terra::extract(class_rast, pts, ID = FALSE)[, 1]
  vals  <- vals[is.finite(vals)]
  total <- length(vals)
  out   <- tibble(scenario = scenario_name, total_points = total)

  for (i in seq_along(class_labels)) {
    lab     <- class_labels[i]
    n_val   <- sum(vals == i)
    pct_val <- if (total == 0) NA_real_ else round(100 * n_val / total, 2)
    out[[paste0(lab, "_n")]]   <- n_val
    out[[paste0(lab, "_pct")]] <- pct_val
  }
  out
}

build_core_occupied_centroid_reference <- function(core_occupied_occ, vars, stats_tbl) {
  max_z_limit <- sqrt(cfg$z2_cap)

  bind_rows(lapply(vars, function(v) {
    x <- core_occupied_occ[[v]]
    x <- x[is.finite(x)]
    if (!length(x)) stop(sprintf("No finite core_occupied values for variable: %s", v))

    v_stats <- stats_tbl[stats_tbl$variable == v, ]
    if (!nrow(v_stats)) stop(sprintf("Missing anchor statistics for variable: %s", v))

    centroid <- v_stats$mu[1]
    sd_low   <- v_stats$sd_lower[1]
    sd_up    <- v_stats$sd_upper[1]
    if (!is.finite(centroid) || !is.finite(sd_low) || sd_low <= 0 ||
        !is.finite(sd_up) || sd_up <= 0) {
      stop(sprintf("Invalid anchor statistics for variable: %s", v))
    }

    empirical_lower  <- min(x, na.rm = TRUE)
    empirical_upper  <- max(x, na.rm = TRUE)
    parametric_lower <- centroid - (max_z_limit * sd_low)
    parametric_upper <- centroid + (max_z_limit * sd_up)

    lower_bound <- pmin(empirical_lower, parametric_lower)
    upper_bound <- pmax(empirical_upper, parametric_upper)

    tibble(
      variable                       = v,
      core_occupied_centroid         = centroid,
      core_occupied_lower_bound      = lower_bound,
      core_occupied_upper_bound      = upper_bound,
      core_occupied_lower_radius     = centroid - lower_bound,
      core_occupied_upper_radius     = upper_bound - centroid,
      n                              = length(x)
    )
  }))
}

# =============================================================================
# SECTION 8 -- OCCURRENCE-LEVEL DIAGNOSTICS
# =============================================================================

run_multimodality_warning <- function(occ_env, vars, mode = "warn") {
  diag_tbl <- lapply(vars, function(v) {
    x <- occ_env[[v]]
    x <- x[is.finite(x)]
    skewness        <- safe_skewness(x)
    excess_kurtosis <- safe_excess_kurtosis(x)
    flag_skew <- is.finite(skewness)        && abs(skewness)        >= 1
    flag_kurt <- is.finite(excess_kurtosis) && abs(excess_kurtosis) >= 1.5
    warning_flag   <- flag_skew || flag_kurt
    warning_reason <- paste(c(
      if (flag_skew) sprintf("high_skewness=%.3f",         skewness),
      if (flag_kurt) sprintf("high_excess_kurtosis=%.3f",  excess_kurtosis)
    ), collapse = " | ")

    tibble(
      variable        = v,
      n               = length(x),
      mean            = mean(x),
      sd              = safe_sd(x),
      skewness        = skewness,
      excess_kurtosis = excess_kurtosis,
      warning_flag    = warning_flag,
      warning_reason  = ifelse(nchar(warning_reason) == 0, "", warning_reason)
    )
  }) %>% bind_rows() %>%
    arrange(desc(warning_flag), desc(abs(skewness)), desc(abs(excess_kurtosis)))

  flagged <- diag_tbl %>% filter(warning_flag)
  if (nrow(flagged) > 0) {
    cat(sprintf(
      "    Distribution-shape warning: %d variable(s) may challenge the single-centroid assumption.\n",
      nrow(flagged)
    ))
    for (i in seq_len(nrow(flagged))) {
      cat(sprintf(
        "    Warning: %s distribution may weaken the single-centroid assumption (%s).\n",
        flagged$variable[i], flagged$warning_reason[i]
      ))
    }
    if (identical(mode, "stop")) {
      stop("Potentially influential distribution shape detected. Review the occurrence-diagnostics CSV for this species.")
    }
  }
  invisible(diag_tbl)
}

safe_bbox_area <- function(df, lon_col, lat_col) {
  if (nrow(df) < 2) return(NA_real_)
  dx <- diff(range(df[[lon_col]], na.rm = TRUE))
  dy <- diff(range(df[[lat_col]], na.rm = TRUE))
  as.numeric(dx * dy)
}

mean_nn_distance <- function(df, lon_col, lat_col) {
  if (nrow(df) < 2) return(NA_real_)
  xy   <- as.matrix(df[, c(lon_col, lat_col), drop = FALSE])
  dmat <- as.matrix(stats::dist(xy))
  diag(dmat) <- Inf
  mean(apply(dmat, 1, min, na.rm = TRUE), na.rm = TRUE)
}

run_core_occupied_bias_diagnostic <- function(core_occupied_df, peripheral_df, occ_df,
                                              lon_col, lat_col, ref_rast) {
  all_bbox_area           <- safe_bbox_area(occ_df,           lon_col, lat_col)
  core_occupied_bbox_area <- safe_bbox_area(core_occupied_df, lon_col, lat_col)
  all_nn                  <- mean_nn_distance(occ_df,           lon_col, lat_col)
  core_occupied_nn        <- mean_nn_distance(core_occupied_df, lon_col, lat_col)

  all_cells <- unique(terra::cellFromXY(
    ref_rast,
    as.matrix(occ_df[, c(lon_col, lat_col), drop = FALSE])
  ))
  core_occupied_cells <- unique(terra::cellFromXY(
    ref_rast,
    as.matrix(core_occupied_df[, c(lon_col, lat_col), drop = FALSE])
  ))
  all_cells           <- all_cells[!is.na(all_cells)]
  core_occupied_cells <- core_occupied_cells[!is.na(core_occupied_cells)]

  diag_tbl <- tibble(
    metric = c(
      "all_occurrences_n",
      "core_occupied_n",
      "peripheral_n",
      "core_occupied_fraction",
      "all_bbox_area_deg2",
      "core_occupied_bbox_area_deg2",
      "core_occupied_bbox_area_ratio",
      "all_mean_nn_distance_deg",
      "core_occupied_mean_nn_distance_deg",
      "core_occupied_mean_nn_ratio",
      "all_occupied_cells",
      "core_occupied_cells",
      "core_occupied_cell_ratio"
    ),
    value = c(
      nrow(occ_df),
      nrow(core_occupied_df),
      nrow(peripheral_df),
      if (nrow(occ_df) > 0) nrow(core_occupied_df) / nrow(occ_df) else NA_real_,
      all_bbox_area,
      core_occupied_bbox_area,
      if (is.finite(all_bbox_area) && all_bbox_area > 0) core_occupied_bbox_area / all_bbox_area else NA_real_,
      all_nn,
      core_occupied_nn,
      if (is.finite(all_nn) && all_nn > 0) core_occupied_nn / all_nn else NA_real_,
      length(all_cells),
      length(core_occupied_cells),
      if (length(all_cells) > 0) length(core_occupied_cells) / length(all_cells) else NA_real_
    )
  )

  area_ratio <- diag_tbl$value[diag_tbl$metric == "core_occupied_bbox_area_ratio"]
  nn_ratio   <- diag_tbl$value[diag_tbl$metric == "core_occupied_mean_nn_ratio"]
  cell_ratio <- diag_tbl$value[diag_tbl$metric == "core_occupied_cell_ratio"]

  warn_flag <- (
    (is.finite(area_ratio) && area_ratio < 0.15) ||
      (is.finite(nn_ratio)   && nn_ratio   < 0.60) ||
      (is.finite(cell_ratio) && cell_ratio < 0.20)
  )

  cat(sprintf(
    "    Core occupied spatial diagnostic: bbox_ratio=%.3f | nn_ratio=%.3f | cell_ratio=%.3f\n",
    area_ratio, nn_ratio, cell_ratio
  ))
  if (isTRUE(warn_flag)) {
    cat("    Warning: core_occupied selection may exhibit spatial concentration / bias.\n")
  }

  invisible(diag_tbl)
}

# =============================================================================
# SECTION 9 -- VALIDATION AND CONSOLIDATED OUTPUT BUILDERS
# =============================================================================

build_validation_summary <- function(shared_tbl, lookup_tbl, current_core_occupied_fate,
                                     stats_tbl,
                                     pixels_censored = NA_real_,
                                     valid_pixels    = NA_real_) {
  current_pct_total <- sum(as.numeric(current_core_occupied_fate[1, c(
    "core_climate_pct", "moderate_departure_pct",
    "restriction_zone_pct", "high_extrapolative_stress_pct"
  )]), na.rm = TRUE)

  truncation_frac <- if (!is.null(stats_tbl) && "truncation_flag_any" %in% names(stats_tbl) && nrow(stats_tbl) > 0) {
    sum(as.logical(stats_tbl$truncation_flag_any), na.rm = TRUE) / nrow(stats_tbl)
  } else NA_real_
  truncation_pass   <- is.finite(truncation_frac)
  truncation_detail <- if (is.finite(truncation_frac)) {
    sprintf("truncated_fraction=%.3f (advisory diagnostic; no acceptance threshold)", truncation_frac)
  } else "truncation_frac=NA (stats_tbl missing or no truncation_flag_any column)"

  censored_frac <- if (is.finite(pixels_censored) && is.finite(valid_pixels) &&
                       (pixels_censored + valid_pixels) > 0) {
    pixels_censored / (pixels_censored + valid_pixels)
  } else NA_real_
  censored_pass   <- is.finite(censored_frac)
  censored_detail <- if (is.finite(censored_frac)) {
    sprintf("censored_fraction=%.4f (reported diagnostic; no acceptance threshold)", censored_frac)
  } else "censored_frac=NA"

  tibble(
    check = c(
      "shared_lookup_stats_equal",
      "current_core_occupied_pct_total_approx_100",
      "truncation_predictor_fraction_reported",
      "censored_pixel_fraction_reported"
    ),
    passed = c(
      nrow(shared_tbl) == nrow(lookup_tbl) && nrow(lookup_tbl) == nrow(stats_tbl),
      abs(current_pct_total - 100) <= 0.5,
      truncation_pass,
      censored_pass
    ),
    detail = c(
      sprintf("shared=%d | lookup=%d | stats=%d",
              nrow(shared_tbl), nrow(lookup_tbl), nrow(stats_tbl)),
      sprintf("current_total=%.2f", current_pct_total),
      truncation_detail,
      censored_detail
    )
  )
}

annotate_table <- function(df, section, description_en, columns_en) {
  df %>% mutate(
    section        = section,
    description_en = description_en,
    columns_en     = columns_en,
    .before        = 1
  )
}

to_long_metrics <- function(df, id_cols, metric_cols, value_col = "value") {
  bind_rows(lapply(metric_cols, function(metric_name) {
    out <- df[, id_cols, drop = FALSE]
    out$metric <- metric_name
    out[[value_col]] <- as.character(df[[metric_name]])
    out
  }))
}

build_run_manifest_csv <- function(run_metadata, validation_tbl) {
  meta_tbl <- run_metadata %>%
    transmute(
      record_type = "run_metadata",
      item = .data$metric,
      value = as.character(.data$value)
    ) %>% annotate_table(
    "run_metadata",
    "This row stores one run setting or summary value for the current analysis.",
    "item=run setting or summary metric; value=stored value for this run."
  )

  val_tbl <- validation_tbl %>%
    transmute(
      record_type = "validation_check",
      item        = check,
      value       = ifelse(passed, "PASS", "FAIL"),
      passed      = passed,
      detail      = detail
    ) %>% annotate_table(
      "validation_checks",
      "This row reports whether an internal consistency check passed.",
      "item=name of the validation check; value=PASS or FAIL; detail=short explanation."
    )

  dict_tbl <- tibble(
    record_type = "table_dictionary",
    item = c(
      paste0(cfg$species_code, "_01_run_manifest"),
      paste0(cfg$species_code, "_02_model_reference"),
      paste0(cfg$species_code, "_03_occurrence_diagnostics"),
      paste0(cfg$species_code, "_04_occurrence_partitions"),
      paste0(cfg$species_code, "_05_current_diagnostics"),
      paste0(cfg$species_code, "_06_proximity_summary")
    ),
    value = c(
      "Run settings and internal validation checks.",
      "Predictor reference tables, anchor statistics, centroid envelope, and thresholds.",
      "Occurrence-based distribution-shape, background-truncation, asymmetry-robustness and spatial-bias diagnostics.",
      "Occurrence records split into climatic core_occupied and peripheral partitions.",
      "Current baseline diagnostic layers, tier distribution, stressor shares (primary, secondary, too-high, too-low), and VPI summary statistics.",
      "Landscape-wide Variable Proximity Index (VPI) summary across all valid pixels."
    )
  ) %>% annotate_table(
    "table_dictionary",
    "This row explains what each exported CSV file is intended to show.",
    "item=CSV file stem; value=plain-language description of that output."
  )

  bind_rows(meta_tbl, val_tbl, dict_tbl)
}

build_model_reference_csv <- function(shared_tbl, lookup_tbl, stats_tbl,
                                      core_occupied_centroid_ref, threshold_tbl,
                                      current_core_occupied_vrs_tbl,
                                      predictor_attribution_tbl = NULL) {
  shared_long <- shared_tbl %>%
    transmute(subject = variable, metric = "retained_for_modeling", value = "TRUE") %>%
    annotate_table(
      "shared_variables",
      "This row shows that one predictor was retained for the current raster stack.",
      "subject=predictor name; metric=reported item; value=stored value."
    )

  lookup_long <- lookup_tbl %>%
    transmute(subject = variable, metric = "stressor_index", value = as.character(stressor_index)) %>%
    annotate_table(
      "primary_stressor_lookup",
      "This row maps a predictor to its primary stressor raster code.",
      "subject=predictor name; metric=reported item; value=stored value."
    )

  anchor_long <- to_long_metrics(
    stats_tbl,
    id_cols = c("variable"),
    metric_cols = c(
      "mu", "sigma_global", "n", "n_lower", "n_upper",
      "tail_scale_center", "raw_tail_scale_lower", "raw_tail_scale_upper",
      "raw_sd_lower", "raw_sd_upper", "shrink_weight_lower", "shrink_weight_upper",
      "sd_lower", "sd_upper", "fallback_lower", "fallback_upper",
      "bootstrap_n", "asym_ratio_raw", "asym_ratio_boot_mean",
      "asym_ratio_ci_low", "asym_ratio_ci_high",
      "asymmetry_bootstrap_stable", "asymmetry_bootstrap_note",
      "bg_lower_tail_prop", "bg_upper_tail_prop",
      "lower_truncation_flag", "upper_truncation_flag", "truncation_flag_any"
    )
  ) %>% rename(subject = variable) %>%
    annotate_table(
      "anchor_statistics",
      "This row stores one occurrence-based anchor statistic used to compute asymmetric VRS scores.",
      "subject=predictor name; metric=reported anchor statistic; value=stored value."
    )

  centroid_long <- to_long_metrics(
    core_occupied_centroid_ref,
    id_cols = c("variable"),
    metric_cols = c(
      "core_occupied_centroid", "core_occupied_lower_bound", "core_occupied_upper_bound",
      "core_occupied_lower_radius", "core_occupied_upper_radius", "n"
    )
  ) %>% rename(subject = variable) %>%
    annotate_table(
      "core_occupied_centroid_reference",
      "This row stores one climatic-envelope value derived from climatic-core_occupied occurrences.",
      "subject=predictor name; metric=reported core_occupied-envelope statistic; value=stored value."
    )

  threshold_long <- if (all(c("quantile", "vrs_threshold") %in% names(threshold_tbl))) {
    threshold_tbl %>%
      transmute(subject = class_boundary,
                metric  = paste0("quantile_", quantile * 100),
                value   = as.character(vrs_threshold)) %>%
      annotate_table(
        "tier_thresholds",
        "This row stores one empirical VRS threshold used to separate climatic tiers.",
        "subject=tier boundary name; metric=quantile used; value=threshold value."
      )
  } else {
    threshold_tbl %>%
      transmute(subject = class_boundary,
                metric  = as.character(metric),
                value   = as.character(value)) %>%
      annotate_table(
        "tier_thresholds",
        "This row stores one Mahalanobis-distance threshold used to separate climatic tiers.",
        "subject=tier boundary name; metric=reported threshold label; value=threshold value."
      )
  }

  current_core_occupied_metric_col <- names(current_core_occupied_vrs_tbl)[1]
  current_core_occupied_vrs_long <- current_core_occupied_vrs_tbl %>%
    mutate(subject = paste0("core_occupied_occurrence_", seq_len(nrow(current_core_occupied_vrs_tbl)))) %>%
    transmute(subject,
              metric = current_core_occupied_metric_col,
              value  = as.character(.data[[current_core_occupied_metric_col]])) %>%
    annotate_table(
      "current_core_occupied_vrs_values",
      "This row stores one tier-reference value extracted from a current core_occupied occurrence.",
      "subject=core_occupied occurrence index; metric=reported item; value=stored value."
    )

  predictor_attribution_long <- if (!is.null(predictor_attribution_tbl) &&
                                    nrow(predictor_attribution_tbl)) {
    to_long_metrics(
      predictor_attribution_tbl,
      id_cols = "variable",
      metric_cols = setdiff(names(predictor_attribution_tbl), "variable")
    ) %>%
      rename(subject = variable) %>%
      annotate_table(
        "predictor_attribution_diagnostic",
        "Reporting-only attribution and cap-audit metrics for predictor contribution, cap saturation, attribution share, truncation and asymmetry evidence; these rows do not alter VRS calculations.",
        "subject=predictor name; metric=reporting audit diagnostic; value=stored raw value."
      )
  } else NULL

  bind_rows(shared_long, lookup_long, anchor_long, centroid_long, threshold_long,
            current_core_occupied_vrs_long, predictor_attribution_long)
}

build_occurrence_diagnostics_csv <- function(multimodality_tbl, core_occupied_bias_tbl, stats_tbl) {
  shape_flagged <- multimodality_tbl %>%
    filter(.data$warning_flag %in% TRUE)

  shape_summary <- tibble(
    subject = "all_predictors",
    metric = c(
      "n_predictors_evaluated",
      "n_distribution_shape_warnings",
      "distribution_shape_warning_fraction",
      "distribution_shape_warning_predictors",
      "distribution_shape_screening_rule"
    ),
    value = as.character(c(
      nrow(multimodality_tbl),
      nrow(shape_flagged),
      if (nrow(multimodality_tbl) > 0) nrow(shape_flagged) / nrow(multimodality_tbl) else NA_real_,
      paste(shape_flagged$variable, collapse = " | "),
      "absolute_skewness>=1 OR absolute_excess_kurtosis>=1.5"
    ))
  ) %>% annotate_table(
    "distribution_shape_summary",
    "This section summarises predictor distributions that may weaken a single arithmetic-centroid representation; it is a screening diagnostic and not a formal multimodality test.",
    "subject=all predictors; metric=summary item; value=count, fraction, predictor list, or screening rule."
  )

  shape_flag_list <- shape_flagged %>%
    transmute(
      subject = .data$variable,
      metric = "distribution_shape_warning_reason",
      value = .data$warning_reason
    ) %>% annotate_table(
      "distribution_shape_flagged_predictors",
      "This section lists every predictor that triggered the distribution-shape screening diagnostic.",
      "subject=predictor name; metric=warning-reason label; value=observed skewness and/or excess kurtosis criterion."
    )

  truncation_flagged <- stats_tbl %>%
    filter(.data$truncation_flag_any %in% TRUE)

  truncation_summary <- tibble(
    subject = "all_predictors",
    metric = c(
      "n_predictors_evaluated",
      "n_background_truncation_flags",
      "background_truncation_flag_fraction",
      "background_truncation_flag_predictors"
    ),
    value = as.character(c(
      nrow(stats_tbl),
      nrow(truncation_flagged),
      if (nrow(stats_tbl) > 0) nrow(truncation_flagged) / nrow(stats_tbl) else NA_real_,
      paste(truncation_flagged$variable, collapse = " | ")
    ))
  ) %>% annotate_table(
    "background_truncation_summary",
    "This section separately summarises predictors whose occupied lower or upper tail approached the available raster-background boundary.",
    "subject=all predictors; metric=summary item; value=count, fraction, or predictor list."
  )

  multi_long <- to_long_metrics(
    multimodality_tbl,
    id_cols = c("variable"),
    metric_cols = c("n", "mean", "sd", "skewness", "excess_kurtosis", "warning_flag", "warning_reason")
  ) %>% rename(subject = variable) %>%
    annotate_table(
      "multimodality_warning",
      "This legacy-named section stores predictor-level distribution-shape diagnostics; warning flags are based on skewness and kurtosis and do not by themselves establish multimodality.",
      "subject=predictor name; metric=reported diagnostic item; value=stored value."
    )

  asym_long <- to_long_metrics(
    stats_tbl,
    id_cols = c("variable"),
    metric_cols = c(
      "n_lower", "n_upper", "tail_scale_center",
      "raw_tail_scale_lower", "raw_tail_scale_upper",
      "raw_sd_lower", "raw_sd_upper",
      "shrink_weight_lower", "shrink_weight_upper",
      "asym_ratio_raw", "asym_ratio_boot_mean",
      "asym_ratio_ci_low", "asym_ratio_ci_high",
      "asymmetry_bootstrap_stable", "asymmetry_bootstrap_note",
      "bg_lower_tail_prop", "bg_upper_tail_prop",
      "lower_truncation_flag", "upper_truncation_flag", "truncation_flag_any"
    )
  ) %>% rename(subject = variable) %>%
    annotate_table(
      "asymmetry_robustness_diagnostic",
      "This row stores one robustness diagnostic for tail asymmetry, including shrinkage weights, bootstrap stability, and conservative background-tail flags.",
      "subject=predictor name; metric=reported robustness item; value=stored value."
    )

  bias_long <- core_occupied_bias_tbl %>%
    transmute(subject = metric, metric = "diagnostic_value", value = as.character(value)) %>%
    annotate_table(
      "core_occupied_spatial_bias_diagnostic",
      "This row stores one spatial diagnostic comparing all occurrences, peripheral points, and the climatic core_occupied.",
      "subject=diagnostic name; metric=reported item; value=stored value."
    )

  bind_rows(shape_summary, shape_flag_list, truncation_summary,
            multi_long, asym_long, bias_long)
}

build_occurrence_partitions_csv <- function(occ_core_occupied, occ_peripheral) {
  core_tbl <- occ_core_occupied %>% mutate(partition = "core_occupied", .before = 1)
  periph_tbl <- if (!is.null(occ_peripheral) && nrow(occ_peripheral) > 0) {
    occ_peripheral %>% mutate(partition = "peripheral", .before = 1)
  } else NULL

  bind_rows(core_tbl, periph_tbl) %>%
    mutate(section = "occurrence_partitions", .before = 1)
}

build_current_diagnostics_csv <- function(diag_summary_tbl, tier_dist_tbl,
                                          primary_share_tbl, primary_unique_share_tbl,
                                          secondary_share_tbl,
                                          too_high_share_tbl, too_low_share_tbl,
                                          vpi_summary_tbl,
                                          reporting_audit_summary_tbl = NULL) {
  out_rows <- list()

  out_rows[[1]] <- to_long_metrics(
    diag_summary_tbl,
    id_cols = c("scenario"),
    metric_cols = c("mean_vrs_global", "max_vrs_global", "mahal_distance_global",
                    "delta_vrs_mean", "delta_vrs_median",
                    "primary_codominant_pixel_pct", "primary_unique_pixel_pct",
                    "valid_pixels", "censored_pixels")
  ) %>% annotate_table(
    "diagnostic_overview",
    "Per-pixel restriction summary statistics under current climate.",
    "scenario=label of the run; metric/value=summary value across all valid pixels."
  )

  out_rows[[2]] <- to_long_metrics(
    tier_dist_tbl,
    id_cols = c("scenario"),
    metric_cols = c("core_climate_pct", "moderate_departure_pct",
                    "restriction_zone_pct", "high_extrapolative_stress_pct",
                    "core_climate_n", "moderate_departure_n",
                    "restriction_zone_n", "high_extrapolative_stress_n")
  ) %>% annotate_table(
    "tier_distribution",
    "Pixel-area share and counts in each four-tier climatic class.",
    "Tier 1=core climate; Tier 2=moderate departure; Tier 3=restriction; Tier 4=high extrapolative stress."
  )

  out_rows[[3]] <- primary_share_tbl %>%
    mutate(scenario = "current", .before = 1) %>%
    annotate_table(
      "primary_stressor_share",
      "Fraction of valid pixels at which the named variable is the largest restriction contributor.",
      "Higher share = the variable acts as the dominant stressor across more of the landscape."
    )

  if (!is.null(primary_unique_share_tbl)) {
    out_rows[[8]] <- primary_unique_share_tbl %>%
      mutate(scenario = "current", .before = 1) %>%
      annotate_table(
        "primary_unique_stressor_share",
        "Fraction of valid pixels at which the named variable is the uniquely largest restriction signal.",
        "Pixels with two or more co-dominant Primary Stressors are excluded from the variable-specific numerator and denominator."
      )
  }

  if (!is.null(secondary_share_tbl)) {
    out_rows[[4]] <- secondary_share_tbl %>%
      mutate(scenario = "current", .before = 1) %>%
      annotate_table(
        "secondary_stressor_share",
        "Fraction of valid pixels at which the named variable is the second-largest contributor.",
        "Useful for diagnosing where pressure is shared between two variables (small Delta)."
      )
  }

  if (!is.null(too_high_share_tbl)) {
    out_rows[[5]] <- too_high_share_tbl %>%
      mutate(scenario = "current", .before = 1) %>%
      annotate_table(
        "too_high_share",
        "Fraction of valid pixels at which the named variable produces the largest upper-tail stress.",
        "TooHigh isolates restriction caused specifically by exceeding the upper-tail climatic optimum."
      )
  }

  if (!is.null(too_low_share_tbl)) {
    out_rows[[6]] <- too_low_share_tbl %>%
      mutate(scenario = "current", .before = 1) %>%
      annotate_table(
        "too_low_share",
        "Fraction of valid pixels at which the named variable produces the largest lower-tail stress.",
        "TooLow isolates restriction caused specifically by falling below the lower-tail climatic optimum."
      )
  }

  out_rows[[7]] <- to_long_metrics(
    vpi_summary_tbl,
    id_cols = c("scenario"),
    metric_cols = c("mean_vpi", "median_vpi", "sd_vpi",
                    "min_vpi", "max_vpi",
                    "p10_vpi", "p25_vpi", "p75_vpi", "p90_vpi")
  ) %>% annotate_table(
    "vpi_summary",
    "Distribution statistics of the Variable Proximity Index across all valid pixels.",
    "VPI = 1 / (1 + mean_vrs); higher values indicate lower climatic restriction and greater climatic proximity."
  )

  if (!is.null(reporting_audit_summary_tbl) && nrow(reporting_audit_summary_tbl)) {
    out_rows[[9]] <- to_long_metrics(
      reporting_audit_summary_tbl,
      id_cols = "scenario",
      metric_cols = setdiff(names(reporting_audit_summary_tbl), "scenario")
    ) %>% annotate_table(
      "reporting_audit_summary",
      "Reporting-only audit summaries for weighted tail evidence and measured cap saturation; no confidence class threshold is applied in the canonical pipeline.",
      "scenario=current; metric=raw reporting audit diagnostic; value=stored value."
    )
  }

  bind_rows(out_rows[!vapply(out_rows, is.null, logical(1))]) %>%
    mutate(section = "current_diagnostics", .before = 1)
}

build_proximity_summary_csv <- function(global_vpi_row) {
  global_vpi_row %>%
    annotate_table(
      "global_proximity",
      "Landscape-wide Variable Proximity Index summary across all valid pixels.",
      "Provides a single climatic-proximity summary number per metric for the entire raster."
    ) %>%
    mutate(section = "proximity_summary", .before = 1)
}

# =============================================================================
# SECTION 10 -- PROXIMITY DIAGNOSTICS (VPI MODULE)
# =============================================================================

summarize_vpi_globally <- function(vpi_rast, mean_vrs_rast, max_vrs_rast, mahal_d_rast,
                                   delta_rast,
                                   pixels_censored, valid_pixels) {
  vpi_vals   <- terra::values(vpi_rast,      mat = FALSE); vpi_vals   <- vpi_vals[is.finite(vpi_vals)]
  vrs_vals   <- terra::values(mean_vrs_rast, mat = FALSE); vrs_vals   <- vrs_vals[is.finite(vrs_vals)]
  max_vals   <- terra::values(max_vrs_rast,  mat = FALSE); max_vals   <- max_vals[is.finite(max_vals)]
  d_vals     <- terra::values(mahal_d_rast,  mat = FALSE); d_vals     <- d_vals[is.finite(d_vals)]
  delta_vals <- terra::values(delta_rast,    mat = FALSE); delta_vals <- delta_vals[is.finite(delta_vals)]

  mean_vpi <- if (length(vpi_vals)) mean(vpi_vals) else NA_real_
  q <- function(p) if (length(vpi_vals)) {
    as.numeric(stats::quantile(vpi_vals, p, na.rm = TRUE, names = FALSE))
  } else NA_real_

  global_diag <- tibble(
    scenario              = "current",
    mean_vrs_global       = if (length(vrs_vals))   mean(vrs_vals)   else NA_real_,
    max_vrs_global        = if (length(max_vals))   mean(max_vals)   else NA_real_,
    mahal_distance_global = if (length(d_vals))     mean(d_vals)     else NA_real_,
    delta_vrs_mean        = if (length(delta_vals)) mean(delta_vals) else NA_real_,
    delta_vrs_median      = if (length(delta_vals)) median(delta_vals) else NA_real_,
    valid_pixels          = as.integer(valid_pixels),
    censored_pixels       = as.integer(pixels_censored)
  )

  global_vpi <- tibble(
    scenario   = "current",
    mean_vpi   = mean_vpi,
    median_vpi = q(0.50),
    sd_vpi     = if (length(vpi_vals) > 1) stats::sd(vpi_vals) else NA_real_,
    min_vpi    = if (length(vpi_vals)) min(vpi_vals) else NA_real_,
    max_vpi    = if (length(vpi_vals)) max(vpi_vals) else NA_real_,
    p10_vpi    = q(0.10),
    p25_vpi    = q(0.25),
    p75_vpi    = q(0.75),
    p90_vpi    = q(0.90)
  )

  list(global_diag = global_diag, global_vpi = global_vpi)
}

summarize_tier_distribution <- function(tier_rast) {
  vals  <- terra::values(tier_rast, mat = FALSE)
  vals  <- vals[is.finite(vals)]
  total <- length(vals)
  if (total == 0) {
    return(tibble(
      scenario                      = "current",
      core_climate_n                = 0L, moderate_departure_n          = 0L,
      restriction_zone_n            = 0L, high_extrapolative_stress_n   = 0L,
      core_climate_pct              = NA_real_, moderate_departure_pct        = NA_real_,
      restriction_zone_pct          = NA_real_, high_extrapolative_stress_pct = NA_real_
    ))
  }
  n_t1 <- sum(vals == 1L); n_t2 <- sum(vals == 2L)
  n_t3 <- sum(vals == 3L); n_t4 <- sum(vals == 4L)
  tibble(
    scenario                      = "current",
    core_climate_n                = as.integer(n_t1),
    moderate_departure_n          = as.integer(n_t2),
    restriction_zone_n            = as.integer(n_t3),
    high_extrapolative_stress_n   = as.integer(n_t4),
    core_climate_pct              = 100 * n_t1 / total,
    moderate_departure_pct        = 100 * n_t2 / total,
    restriction_zone_pct          = 100 * n_t3 / total,
    high_extrapolative_stress_pct = 100 * n_t4 / total
  )
}

summarize_index_share <- function(idx_rast, lookup_tbl, label_prefix) {
  vals <- terra::values(idx_rast, mat = FALSE)
  vals <- vals[is.finite(vals)]
  total <- length(vals)
  if (total == 0) return(NULL)
  tab <- table(vals)
  tibble(
    metric = paste0(label_prefix, "_",
                    lookup_tbl$variable[match(as.integer(names(tab)), lookup_tbl$stressor_index)],
                    "_pct"),
    value  = as.character(round(100 * as.numeric(tab) / total, 4))
  )
}

join_index_share <- function(predictor_tbl, share_tbl, lookup_tbl,
                             label_prefix, output_column) {
  values <- rep(0, nrow(lookup_tbl))
  if (!is.null(share_tbl) && nrow(share_tbl)) {
    expected <- paste0(label_prefix, "_", lookup_tbl$variable, "_pct")
    hit <- match(expected, share_tbl$metric)
    values[!is.na(hit)] <- suppressWarnings(
      as.numeric(share_tbl$value[hit[!is.na(hit)]])
    )
  }
  share_lookup <- tibble(variable = lookup_tbl$variable, share_value = values)
  names(share_lookup)[2] <- output_column
  left_join(predictor_tbl, share_lookup, by = "variable")
}

build_predictor_attribution <- function(base_tbl, stats_tbl, lookup_tbl,
                                            primary_share_tbl,
                                            secondary_share_tbl) {
  if (is.null(base_tbl) || !nrow(base_tbl)) return(tibble())
  equivalence_low <- 1 / cfg$asymmetry_ratio_threshold
  out <- base_tbl %>%
    left_join(
      stats_tbl %>%
        transmute(
          variable,
          truncation_flag_any,
          asymmetry_bootstrap_stable,
          asym_ratio_ci_low,
          asym_ratio_ci_high,
          asymmetry_evidence = case_when(
            .data$asymmetry_bootstrap_stable %in% TRUE ~ "stable_asymmetric",
            is.finite(.data$asym_ratio_ci_low) &
              is.finite(.data$asym_ratio_ci_high) &
              .data$asym_ratio_ci_low >= equivalence_low &
              .data$asym_ratio_ci_high <= cfg$asymmetry_ratio_threshold ~
              "stable_near_symmetric",
            TRUE ~ "unresolved"
          )
        ),
      by = "variable"
    )
  out <- join_index_share(out, primary_share_tbl, lookup_tbl,
                          "primary", "primary_share_pct")
  out <- join_index_share(out, secondary_share_tbl, lookup_tbl,
                          "secondary", "secondary_share_pct")
  out %>%
    mutate(
      attribution_share_raw = .data$primary_share_pct +
        .data$secondary_share_pct,
      attribution_weight_pct = if (sum(.data$attribution_share_raw,
                                       na.rm = TRUE) > 0) {
        100 * .data$attribution_share_raw /
          sum(.data$attribution_share_raw, na.rm = TRUE)
      } else NA_real_
    )
}

build_reporting_audit_summary <- function(predictor_tbl, cap_summary) {
  if (is.null(predictor_tbl) || !nrow(predictor_tbl)) return(tibble())
  trunc <- ifelse(is.na(predictor_tbl$truncation_flag_any), NA,
                  predictor_tbl$truncation_flag_any %in% TRUE)
  unresolved <- ifelse(is.na(predictor_tbl$asymmetry_evidence), NA,
                       predictor_tbl$asymmetry_evidence == "unresolved")
  contribution_weight <- predictor_tbl$mean_vrs_contribution_pct / 100
  attribution_weight <- predictor_tbl$attribution_weight_pct / 100
  weighted <- function(flag, weight) {
    ok <- !is.na(flag) & is.finite(weight)
    if (!any(ok) || sum(weight[ok]) <= 0) return(NA_real_)
    sum(weight[ok] * flag[ok]) / sum(weight[ok])
  }
  tibble(
    scenario = "current",
    truncation_fraction_unweighted = mean(trunc, na.rm = TRUE),
    truncation_fraction_vrs_weighted = weighted(trunc, contribution_weight),
    truncation_fraction_attribution_weighted = weighted(trunc, attribution_weight),
    unresolved_asymmetry_fraction_unweighted = mean(unresolved, na.rm = TRUE),
    unresolved_asymmetry_fraction_vrs_weighted = weighted(unresolved, contribution_weight),
    unresolved_asymmetry_fraction_attribution_weighted = weighted(unresolved, attribution_weight),
    upper_and_lower_maxima_at_cap_pct = cap_summary$upper_and_lower_maxima_at_cap_pct,
    balanced_pixels_explained_by_dual_tail_cap_pct = cap_summary$balanced_pixels_explained_by_dual_tail_cap_pct,
    balanced_pixels_without_dual_tail_cap_pct = cap_summary$balanced_pixels_without_dual_tail_cap_pct,
    any_predictor_at_cap_pct = cap_summary$any_predictor_at_cap_pct,
    mean_predictors_at_cap_per_pixel = cap_summary$mean_predictors_at_cap_per_pixel
  )
}

# =============================================================================
# SECTION 10B -- ADD-ON DIAGNOSTICS
# =============================================================================

rank01_raster <- function(r) {
  vals <- terra::values(r, mat = FALSE)
  out <- rep(NA_real_, length(vals))
  ok <- is.finite(vals)
  n <- sum(ok)
  if (n == 1L) {
    out[ok] <- 0
  } else if (n > 1L) {
    out[ok] <- (rank(vals[ok], ties.method = "average") - 1) / (n - 1)
  }
  result <- terra::setValues(r, out)
  result
}

write_class_factor <- function(class_rast, labels, filename) {
  out <- as.factor(class_rast)
  levels(out) <- data.frame(value = seq_along(labels), label = labels)
  terra::writeRaster(out, filename, overwrite = TRUE, wopt = gdal_opts)
  invisible(out)
}

summarize_class_raster <- function(class_rast, labels, definition) {
  vals <- terra::values(class_rast, mat = FALSE)
  vals <- vals[is.finite(vals)]
  total <- length(vals)
  counts <- tabulate(as.integer(vals), nbins = length(labels))
  tibble(
    class_id = seq_along(labels),
    class_name = labels,
    pixel_count = as.integer(counts),
    definition = definition,
    percent_valid = if (total > 0) 100 * counts / total else NA_real_
  )
}

classify_tail_direction <- function(high_signal, low_signal, epsilon, tau) {
  total <- high_signal + low_signal
  imbalance <- abs(high_signal - low_signal) / total
  out <- ifel(total <= epsilon, 4,
              ifel(imbalance <= tau, 3,
                   ifel(high_signal > low_signal, 1, 2)))
  out <- ifel(is.na(high_signal) | is.na(low_signal), NA, out)
  names(out) <- "tail_direction_class"
  list(class = out, total = total, imbalance = imbalance)
}

run_addon_diagnostics <- function(mean_vrs, mahal_distance, current_res,
                                  template, out_paths, species_code) {
  cat(">>> Building integrated add-on diagnostics...\n")
  vrs_rank <- rank01_raster(mean_vrs)
  mahal_rank <- rank01_raster(mahal_distance)
  names(vrs_rank) <- "vrs_rank01"
  names(mahal_rank) <- "mahal_rank01"
  divergence <- vrs_rank - mahal_rank
  abs_disagreement <- abs(divergence)
  names(divergence) <- "mahal_vrs_divergence"
  names(abs_disagreement) <- "mahal_vrs_abs_disagreement"
  delta <- cfg$divergence_threshold
  agreement <- ifel(divergence > delta, 1,
                    ifel(divergence < -delta, 2,
                         ifel(vrs_rank >= 0.5, 3, 4)))
  agreement <- ifel(is.na(divergence), NA, agreement)
  names(agreement) <- "mahal_vrs_agreement_4class"

  prefix <- file.path(out_paths$addon_rasters, paste0(species_code, "_current"))
  terra::writeRaster(vrs_rank, paste0(prefix, "_vrs_rank01.tif"), overwrite = TRUE, wopt = gdal_opts)
  terra::writeRaster(mahal_rank, paste0(prefix, "_mahal_rank01.tif"), overwrite = TRUE, wopt = gdal_opts)
  terra::writeRaster(divergence, paste0(prefix, "_mahal_vrs_divergence.tif"), overwrite = TRUE, wopt = gdal_opts)
  terra::writeRaster(abs_disagreement, paste0(prefix, "_mahal_vrs_abs_disagreement.tif"), overwrite = TRUE, wopt = gdal_opts)
  terra::writeRaster(agreement, paste0(prefix, "_mahal_vrs_agreement_4class_idx.tif"), overwrite = TRUE, wopt = gdal_opts)

  agreement_labels <- c(
    "1 VRS-dominant_higher_relative_restriction_rank",
    "2 Mahalanobis-dominant_higher_relative_distance_rank",
    "3 Concordant_higher_rank_signals",
    "4 Concordant_lower_rank_signals"
  )
  write_class_factor(agreement, agreement_labels,
                     paste0(prefix, "_mahal_vrs_agreement_4class.tif"))

  agreement_counts <- summarize_class_raster(
    agreement, agreement_labels, "mahalanobis_vrs_rank_agreement"
  )
  write.csv(agreement_counts,
            file.path(out_paths$addon_csvs, "addon02_divergence_class_counts.csv"),
            row.names = FALSE)

  div_vals <- terra::values(divergence, mat = FALSE)
  div_vals <- div_vals[is.finite(div_vals)]
  abs_vals <- abs(div_vals)
  divergence_summary <- tibble(
    metric = c("species_code", "file_prefix", "n_valid_pixels", "mean_divergence",
               "median_divergence", "sd_divergence", "mean_abs_disagreement",
               "median_abs_disagreement", "q75_abs_disagreement",
               "percent_vrs_dominant", "percent_mahalanobis_dominant",
               "percent_within_threshold", "divergence_threshold"),
    value = as.character(c(
      species_code, species_code, length(div_vals), mean(div_vals), median(div_vals),
      stats::sd(div_vals), mean(abs_vals), median(abs_vals),
      as.numeric(stats::quantile(abs_vals, 0.75, names = FALSE, type = 7)),
      100 * mean(div_vals > delta), 100 * mean(div_vals < -delta),
      100 * mean(abs_vals <= delta), delta
    ))
  )
  write.csv(divergence_summary,
            file.path(out_paths$addon_csvs, "addon02_divergence_summary.csv"),
            row.names = FALSE)

  dominant <- classify_tail_direction(
    current_res$too_high_z, current_res$too_low_z,
    cfg$tail_no_signal_epsilon, cfg$tail_balance_tolerance
  )
  net <- classify_tail_direction(
    current_res$too_high_z_sum, current_res$too_low_z_sum,
    cfg$tail_no_signal_epsilon, cfg$tail_balance_tolerance
  )

  cap_threshold <- cfg$z2_cap - cfg$equality_tolerance
  valid_direction <- !is.na(dominant$class)
  dual_tail_cap <- valid_direction &
    current_res$too_high_z >= cap_threshold &
    current_res$too_low_z >= cap_threshold
  dominant_balanced <- valid_direction & dominant$class == 3
  balanced_due_to_cap <- dominant_balanced & dual_tail_cap
  balanced_without_dual_cap <- dominant_balanced & !dual_tail_cap
  any_predictor_at_cap <- valid_direction & current_res$max_vrs >= cap_threshold

  valid_direction_n <- suppressWarnings(as.numeric(
    terra::global(valid_direction, "sum", na.rm = TRUE)[1, 1]
  ))
  pct_valid <- function(condition) {
    n <- suppressWarnings(as.numeric(
      terra::global(condition, "sum", na.rm = TRUE)[1, 1]
    ))
    if (is.finite(n) && is.finite(valid_direction_n) && valid_direction_n > 0) {
      100 * n / valid_direction_n
    } else NA_real_
  }
  dominant_balanced_n <- suppressWarnings(as.numeric(
    terra::global(dominant_balanced, "sum", na.rm = TRUE)[1, 1]
  ))
  balanced_due_to_cap_n <- suppressWarnings(as.numeric(
    terra::global(balanced_due_to_cap, "sum", na.rm = TRUE)[1, 1]
  ))
  cap_summary <- list(
    upper_and_lower_maxima_at_cap_pct = pct_valid(dual_tail_cap),
    balanced_pixels_explained_by_dual_tail_cap_pct = pct_valid(balanced_due_to_cap),
    balanced_pixels_without_dual_tail_cap_pct = pct_valid(balanced_without_dual_cap),
    balanced_due_to_cap_share_of_balanced_pct =
      if (is.finite(dominant_balanced_n) && dominant_balanced_n > 0) {
        100 * balanced_due_to_cap_n / dominant_balanced_n
      } else NA_real_,
    any_predictor_at_cap_pct = pct_valid(any_predictor_at_cap),
    mean_predictors_at_cap_per_pixel = if (!is.null(current_res$predictor_diagnostics) &&
                                              nrow(current_res$predictor_diagnostics)) {
      sum(current_res$predictor_diagnostics$cap_saturation_fraction, na.rm = TRUE)
    } else NA_real_
  )
  tail_labels <- c(
    "1 Upper-tail_dominant_high-edge_restriction",
    "2 Lower-tail_dominant_low-edge_restriction",
    "3 Balanced_or_near-symmetric",
    "4 No_directional_signal"
  )
  terra::writeRaster(dominant$class, paste0(prefix, "_tail_direction_dominant_idx.tif"),
                     overwrite = TRUE, wopt = gdal_opts)
  terra::writeRaster(net$class, paste0(prefix, "_tail_direction_net_idx.tif"),
                     overwrite = TRUE, wopt = gdal_opts)
  write_class_factor(dominant$class, tail_labels,
                     paste0(prefix, "_tail_direction_dominant.tif"))
  write_class_factor(net$class, tail_labels,
                     paste0(prefix, "_tail_direction_net.tif"))

  dominant_counts <- summarize_class_raster(dominant$class, tail_labels,
                                             "dominant_max_tail_signal")
  net_counts <- summarize_class_raster(net$class, tail_labels,
                                       "net_summed_tail_signal")
  tail_counts <- bind_rows(dominant_counts, net_counts)
  write.csv(tail_counts,
            file.path(out_paths$addon_csvs, "addon01_tail_direction_pixel_counts.csv"),
            row.names = FALSE)

  dom_vals <- terra::values(dominant$class, mat = FALSE)
  net_vals <- terra::values(net$class, mat = FALSE)
  ok <- is.finite(dom_vals) & is.finite(net_vals)
  agreement_pct <- if (any(ok)) 100 * mean(dom_vals[ok] == net_vals[ok]) else NA_real_
  cross <- as.data.frame(table(
    dominant_class = factor(dom_vals[ok], levels = seq_along(tail_labels), labels = tail_labels),
    net_class = factor(net_vals[ok], levels = seq_along(tail_labels), labels = tail_labels)
  ))
  names(cross)[3] <- "pixel_count"
  cross$percent_valid <- if (sum(cross$pixel_count) > 0) {
    100 * cross$pixel_count / sum(cross$pixel_count)
  } else NA_real_
  cross$overall_agreement_pct <- agreement_pct
  write.csv(cross,
            file.path(out_paths$addon_csvs, "addon01_tail_direction_dominant_vs_net_agreement.csv"),
            row.names = FALSE)

  tail_summary <- tibble(
    metric = c("species_code", "n_valid_pixels", "overall_dominant_net_agreement_pct",
               "balance_tolerance", "no_signal_epsilon",
               "upper_and_lower_maxima_at_cap_pct",
               "balanced_pixels_explained_by_dual_tail_cap_pct",
               "balanced_pixels_without_dual_tail_cap_pct",
               "balanced_due_to_cap_share_of_balanced_pct",
               "any_predictor_at_cap_pct",
               "mean_predictors_at_cap_per_pixel"),
    value = as.character(c(species_code, sum(ok), agreement_pct,
                           cfg$tail_balance_tolerance, cfg$tail_no_signal_epsilon,
                           cap_summary$upper_and_lower_maxima_at_cap_pct,
                           cap_summary$balanced_pixels_explained_by_dual_tail_cap_pct,
                           cap_summary$balanced_pixels_without_dual_tail_cap_pct,
                           cap_summary$balanced_due_to_cap_share_of_balanced_pct,
                           cap_summary$any_predictor_at_cap_pct,
                           cap_summary$mean_predictors_at_cap_per_pixel))
  )
  write.csv(tail_summary,
            file.path(out_paths$addon_csvs, "addon01_tail_direction_summary.csv"),
            row.names = FALSE)

  if (terra::is.lonlat(template)) {
    net_all <- terra::values(net$class, mat = FALSE)
    cells <- which(is.finite(net_all))
    xy <- terra::xyFromCell(template, cells)
    width <- cfg$latitude_band_width
    lower <- floor(xy[, 2] / width) * width
    upper <- lower + width
    lat_profile <- tibble(
      lat_band = paste0("(", lower, ",", upper, "]"),
      tail_class = tail_labels[as.integer(net_all[cells])]
    ) %>%
      count(lat_band, tail_class, name = "n_pixels")
    write.csv(lat_profile,
              file.path(out_paths$addon_csvs, "addon01_tail_direction_latitudinal_profile.csv"),
              row.names = FALSE)
  } else {
    warning("Latitudinal tail profile was not written because the raster template is projected; provide a latitude raster or transform cell centres explicitly.")
  }

  invisible(list(
    vrs_rank = vrs_rank, mahal_rank = mahal_rank, divergence = divergence,
    abs_disagreement = abs_disagreement, agreement = agreement,
    dominant_tail = dominant$class, net_tail = net$class,
    cap_summary = cap_summary
  ))
}

# =============================================================================
# SECTION 11 -- STACK PREPARATION (single species; optional polygon mask)
# =============================================================================

detect_species_label <- function(occ_csv, fallback) {
  preview <- try(utils::read.csv(
    occ_csv,
    nrows = 25,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ), silent = TRUE)

  if (inherits(preview, "try-error") || !nrow(preview)) return(fallback)

  candidates <- c("Species", "species", "scientificName", "scientific_name")
  species_col <- candidates[candidates %in% names(preview)][1]
  if (is.na(species_col)) return(fallback)

  labels <- trimws(as.character(preview[[species_col]]))
  labels <- unique(labels[!is.na(labels) & nzchar(labels)])
  if (length(labels)) labels[1] else fallback
}

# Prepare the raster stack for one predictor profile. The Turkiye variable
# rasters are already region-cropped, so a polygon mask is optional; the mask
# step is executed only when cfg$mask_shp is provided and the file exists.
prepare_stack <- function(cfg) {
  vars_required <- cfg$focal_vars

  current_tbl <- list_tifs_named(cfg$current_vars_dir)
  missing_focal <- setdiff(vars_required, current_tbl$variable)
  if (length(missing_focal) && !isTRUE(cfg$allow_missing_predictors)) {
    stop("Required focal predictor rasters are missing: ",
         paste(missing_focal, collapse = ", "))
  }
  vars_common <- vars_required[vars_required %in% current_tbl$variable]
  if (!length(vars_common)) {
    stop("No shared variables found for the selected predictor profile.")
  }

  cat(sprintf("\n>>> Building raster stack for the %d-predictor profile...\n",
              length(cfg$focal_vars)))
  current_stack <- rast(current_tbl$file[match(vars_common, current_tbl$variable)])
  names(current_stack) <- vars_common
  validate_raster_stack(current_stack, vars_common)

  if (!is.null(cfg$mask_shp) && nzchar(cfg$mask_shp) && file.exists(cfg$mask_shp)) {
    cat(">>> Applying optional polygon mask: ", cfg$mask_shp, "\n", sep = "")
    poly_mask <- vect(cfg$mask_shp)
    if (!terra::same.crs(poly_mask, current_stack)) {
      poly_mask <- project(poly_mask, crs(current_stack))
    }
    if (nrow(poly_mask) > 1L) {
      poly_mask <- aggregate(poly_mask)
    }
    current_stack <- crop(current_stack, poly_mask, snap = "out")
    current_stack <- mask(current_stack, poly_mask)
    names(current_stack) <- vars_common
  } else {
    cat(">>> No polygon mask supplied; using region-cropped Turkiye rasters as-is.\n")
  }

  validate_raster_stack(current_stack, vars_common)
  cat(sprintf("    Predictors prepared: %d\n", terra::nlyr(current_stack)))
  current_stack
}

# =============================================================================
# SECTION 12 -- SINGLE-SPECIES EXECUTION
# =============================================================================

run_vera_current <- function(run_cfg, prepared_stack) {
  assign("cfg", run_cfg, envir = .GlobalEnv)
  assign("gdal_opts", list(gdal = run_cfg$raster_compression), envir = .GlobalEnv)
  validate_config(run_cfg)
  set.seed(as.integer(run_cfg$random_seed))

out_paths <- build_output_paths(cfg$output_dir, temp_name = cfg$temp_subdir)
terraOptions(
  tempdir  = out_paths$temp,
  memfrac  = cfg$memfrac,
  todisk   = cfg$todisk,
  progress = cfg$progress,
  verbose  = FALSE
)

# Clear incomplete terra products from an earlier interrupted execution. This
# directory contains scratch data only; persistent analysis outputs are stored
# in the sibling csvs/ and rasters/ directories.
stale_temp <- list.files(out_paths$temp, full.names = TRUE, all.files = TRUE,
                         no.. = TRUE)
if (length(stale_temp)) unlink(stale_temp, recursive = TRUE, force = TRUE)

cat("\n>>> Using the prepared raster stack...\n")
current_stack <- prepared_stack
vars_common <- names(current_stack)
if (!identical(vars_common, cfg$focal_vars)) {
  stop("Prepared raster stack does not match the selected predictor profile.")
}
lookup_tbl <- write_lookup_outputs(vars_common, out_paths)
shared_tbl <- tibble(variable = vars_common)

cat(sprintf("    Shared predictors retained: %d\n", length(vars_common)))
validate_raster_stack(current_stack, vars_common)

cat(">>> Reading occurrence data and building asymmetric anchor table...\n")
occ        <- read.csv(cfg$occ_csv)

# Standardise the taxon label written to occurrence-partition outputs. This
# corrects legacy spellings in the input table without changing coordinates or
# any environmental values used by the analysis.
species_col_candidates <- c("Species", "species", "scientificName", "scientific_name")
species_col <- species_col_candidates[species_col_candidates %in% names(occ)][1]
if (is.na(species_col)) {
  occ$Species <- cfg$species_label
} else {
  occ[[species_col]] <- cfg$species_label
}

lon_col    <- pick_col(names(occ), cfg$lon_candidates, "longitude")
lat_col    <- pick_col(names(occ), cfg$lat_candidates, "latitude")
anchor_obj <- build_stats_from_occ(occ, current_stack, lon_col, lat_col, cfg$min_tail_n)
stats_tbl  <- anchor_obj$stats
occ_env    <- anchor_obj$occ_env

multimodality_tbl <- run_multimodality_warning(occ_env, vars_common,
                                               mode = cfg$multimodality_mode)

cat(">>> Selecting Climatic Core occupied climate from current occurrences using Mahalanobis empirical quantile...\n")
mahal_split <- calc_mahal_core_occupied(
  occ_env, vars_common,
  core_occupied_prob = cfg$core_occupied_prob
)
occ_core_occupied <- mahal_split$core_occupied
occ_peripheral    <- mahal_split$peripheral

core_occupied_centroid_ref <- build_core_occupied_centroid_reference(
  occ_core_occupied, vars_common, stats_tbl
)
core_partition_ridge_values <- sort(unique(occ_core_occupied$mahal_ridge))
retained_occ_n     <- nrow(occ_env)
core_occupied_pts  <- vect(occ_core_occupied, geom = c(lon_col, lat_col), crs = "EPSG:4326")
cat(sprintf("    Core occupied climate points: %d | Peripheral points: %d\n",
            nrow(occ_core_occupied), nrow(occ_peripheral)))

core_occupied_bias_tbl <- run_core_occupied_bias_diagnostic(
  occ_core_occupied, occ_peripheral, occ_env, lon_col, lat_col,
  current_stack[[1]]
)

class_labels <- c("core_climate", "moderate_departure",
                  "restriction_zone", "high_extrapolative_stress")

mahal_occ_ref <- occ_env
if (isTRUE(cfg$precompute_mahal_reference)) {
  tier_mahal_ref <- build_mahal_reference(mahal_occ_ref, vars_common, label = "mahal_reference")
  cat(sprintf("    Mahalanobis reference cached: %s (n=%d, ridge=%g)\n",
              cfg$mahal_reference, tier_mahal_ref$n, tier_mahal_ref$ridge))
} else {
  tier_mahal_ref <- mahal_occ_ref
}

cat(">>> Building asymmetric diagnostic surfaces...\n")
current_res   <- calc_asym_surfaces(current_stack, stats_tbl, "current")
current_mahal <- calc_mahal_distance_surface(current_stack, tier_mahal_ref, vars_common, "current")

core_occupied_distances_current <- as.numeric(terra::extract(
  current_mahal$distance, core_occupied_pts, ID = FALSE
)[, 1])
core_occupied_distances_current <- core_occupied_distances_current[is.finite(core_occupied_distances_current)]

tier_breaks_resolved <- resolve_tier_breaks(
  method                  = cfg$tier_breaks_method,
  core_occupied_distances = core_occupied_distances_current,
  quantiles               = cfg$tier_breaks_quantiles,
  fixed_breaks            = cfg$mahal_distance_breaks,
  n_vars                  = length(vars_common)
)
cat("    Tier method          : Mahalanobis distance\n")
cat(sprintf("    Tier breaks method   : %s\n", cfg$tier_breaks_method))
cat(sprintf("    Tier breaks resolved : %s\n",
            paste(round(tier_breaks_resolved, 4), collapse = " | ")))

current_class <- classify_mahal_tiers(current_mahal$distance, tier_breaks_resolved)

threshold_tbl <- tibble(
  class_boundary = class_labels[-1],
  metric         = paste0("mahal_distance_break_", 1:3),
  value          = as.character(tier_breaks_resolved)
)
current_core_occupied_vrs_tbl <- tibble(
  current_core_occupied_mahal_distance = core_occupied_distances_current
)

if (isTRUE(cfg$write_fixed_break_comparator) && !identical(cfg$tier_breaks_method, "fixed")) {
  current_class_fixed <- classify_mahal_tiers(current_mahal$distance, cfg$mahal_distance_breaks)
  writeRaster(current_class_fixed,
              file.path(out_paths$raster_core_occupied,
                        paste0(cfg$species_code, "_current_four_tier_status_fixedbreaks.tif")),
              overwrite = TRUE, wopt = gdal_opts)
  rm(current_class_fixed)
}

current_class_factor <- as.factor(current_class)
levels(current_class_factor) <- data.frame(value = 1:4, label = class_labels)

writeRaster(current_class_factor,
            file.path(out_paths$raster_core_occupied,
                      paste0(cfg$species_code, "_current_four_tier_status.tif")),
            overwrite = TRUE, wopt = gdal_opts)

out_stub_core <- file.path(out_paths$raster_core_occupied, paste0(cfg$species_code, "_current"))
out_stub_diag <- file.path(out_paths$raster_optional,      paste0(cfg$species_code, "_current"))

writeRaster(current_res$mean_vrs,   paste0(out_stub_core, "_mean_vrs.tif"),       overwrite = TRUE, wopt = gdal_opts)
writeRaster(current_res$max_vrs,    paste0(out_stub_core, "_max_vrs.tif"),        overwrite = TRUE, wopt = gdal_opts)
writeRaster(current_mahal$distance, paste0(out_stub_core, "_mahal_distance.tif"), overwrite = TRUE, wopt = gdal_opts)
write_factor_raster(current_res$stressor_index, lookup_tbl,
                    paste0(out_stub_core, "_primary_stressor_factor.tif"))

if (!is.null(current_res$stressor_index) && isTRUE(cfg$write_index_rasters)) {
  writeRaster(current_res$stressor_index, paste0(out_stub_core, "_primary_stressor_idx.tif"), overwrite = TRUE, wopt = gdal_opts)
}

if (isTRUE(cfg$write_index_rasters)) {
  if (!is.null(current_res$secondary_index)) writeRaster(current_res$secondary_index, paste0(out_stub_diag, "_secondary_stressor_idx.tif"), overwrite = TRUE, wopt = gdal_opts)
  if (!is.null(current_res$too_high_index))  writeRaster(current_res$too_high_index,  paste0(out_stub_diag, "_too_high_idx.tif"),           overwrite = TRUE, wopt = gdal_opts)
  if (!is.null(current_res$too_low_index))   writeRaster(current_res$too_low_index,   paste0(out_stub_diag, "_too_low_idx.tif"),            overwrite = TRUE, wopt = gdal_opts)
}

# Publication-facing diagnostics are always written. The write_index_rasters
# switch controls only the auxiliary integer-code rasters above.
writeRaster(current_res$delta_vrs,
            paste0(out_stub_diag, "_delta_vrs.tif"),
            overwrite = TRUE, wopt = gdal_opts)
writeRaster(current_res$primary_co_dominance_count,
            paste0(out_stub_diag, "_primary_co_dominance_count.tif"),
            overwrite = TRUE, wopt = gdal_opts)
writeRaster(current_res$too_high_co_dominance_count,
            paste0(out_stub_diag, "_too_high_co_dominance_count.tif"),
            overwrite = TRUE, wopt = gdal_opts)
writeRaster(current_res$too_low_co_dominance_count,
            paste0(out_stub_diag, "_too_low_co_dominance_count.tif"),
            overwrite = TRUE, wopt = gdal_opts)

write_factor_raster(
  current_res$secondary_index, lookup_tbl,
  paste0(out_stub_diag, "_secondary_stressor_factor.tif")
)
write_factor_raster(
  current_res$primary_unique_index, lookup_tbl,
  paste0(out_stub_diag, "_primary_unique_stressor_factor.tif")
)
write_factor_raster(
  current_res$too_high_index, lookup_tbl,
  paste0(out_stub_diag, "_too_high_factor.tif")
)
write_factor_raster(
  current_res$too_low_index, lookup_tbl,
  paste0(out_stub_diag, "_too_low_factor.tif")
)

writeRaster(current_res$too_high_z,
            paste0(out_stub_diag, "_too_high_zsq_max.tif"),
            overwrite = TRUE, wopt = gdal_opts)
writeRaster(current_res$too_low_z,
            paste0(out_stub_diag, "_too_low_zsq_max.tif"),
            overwrite = TRUE, wopt = gdal_opts)
writeRaster(current_res$too_high_z_sum,
            paste0(out_stub_diag, "_too_high_zsq_sum.tif"),
            overwrite = TRUE, wopt = gdal_opts)
writeRaster(current_res$too_low_z_sum,
            paste0(out_stub_diag, "_too_low_zsq_sum.tif"),
            overwrite = TRUE, wopt = gdal_opts)

if (isTRUE(cfg$write_optional_diagnostic_rasters)) {
  if (!is.null(current_res$valid_count)) {
    writeRaster(current_res$valid_count, paste0(out_stub_diag, "_valid_var_count.tif"),
                overwrite = TRUE, wopt = gdal_opts)
  }
}

current_core_occupied_fate <- extract_class_fate(current_class, core_occupied_pts,
                                                 "current_core", class_labels)

valid_pixels_total <- as.integer(terra::global(current_class, "notNA")[1, 1])

cat(sprintf("    Pixels censored on min_valid_frac: %d\n", current_res$pixels_censored))

vpi_rast        <- 1 / (1 + current_res$mean_vrs)
names(vpi_rast) <- "vpi"

writeRaster(vpi_rast, paste0(out_stub_core, "_vpi.tif"),
            overwrite = TRUE, wopt = gdal_opts)

addon_obj <- run_addon_diagnostics(
  mean_vrs       = current_res$mean_vrs,
  mahal_distance = current_mahal$distance,
  current_res    = current_res,
  template       = current_stack[[1]],
  out_paths      = out_paths,
  species_code   = cfg$species_code
)

global_obj <- summarize_vpi_globally(
  vpi_rast        = vpi_rast,
  mean_vrs_rast   = current_res$mean_vrs,
  max_vrs_rast    = current_res$max_vrs,
  mahal_d_rast    = current_mahal$distance,
  delta_rast      = current_res$delta_vrs,
  pixels_censored = current_res$pixels_censored,
  valid_pixels    = valid_pixels_total
)

primary_co_dominance_vals <- terra::values(current_res$primary_co_dominance_count, mat = FALSE)
primary_co_dominance_vals <- primary_co_dominance_vals[is.finite(primary_co_dominance_vals)]
primary_codominant_pixel_pct <- if (length(primary_co_dominance_vals)) {
  100 * mean(primary_co_dominance_vals > 1)
} else NA_real_
primary_unique_pixel_pct <- if (length(primary_co_dominance_vals)) {
  100 * mean(primary_co_dominance_vals == 1)
} else NA_real_

cat(">>> Writing consolidated CSV outputs...\n")

tier_dist_tbl <- summarize_tier_distribution(current_class)

primary_share_tbl   <- summarize_index_share(current_res$stressor_index,  lookup_tbl, "primary")
primary_unique_share_tbl <- summarize_index_share(current_res$primary_unique_index,
                                                  lookup_tbl, "primary_unique")
secondary_share_tbl <- summarize_index_share(current_res$secondary_index, lookup_tbl, "secondary")
too_high_share_tbl  <- summarize_index_share(current_res$too_high_index,  lookup_tbl, "too_high")
too_low_share_tbl   <- summarize_index_share(current_res$too_low_index,   lookup_tbl, "too_low")

predictor_attribution_tbl <- build_predictor_attribution(
  base_tbl = current_res$predictor_diagnostics,
  stats_tbl = stats_tbl,
  lookup_tbl = lookup_tbl,
  primary_share_tbl = primary_share_tbl,
  secondary_share_tbl = secondary_share_tbl
)
reporting_audit_summary_tbl <- build_reporting_audit_summary(
  predictor_tbl = predictor_attribution_tbl,
  cap_summary = addon_obj$cap_summary
)

global_diag_tbl <- global_obj$global_diag %>%
  mutate(
    primary_codominant_pixel_pct = primary_codominant_pixel_pct,
    primary_unique_pixel_pct = primary_unique_pixel_pct
  )

current_diagnostics_csv <- build_current_diagnostics_csv(
  diag_summary_tbl    = global_diag_tbl,
  tier_dist_tbl       = tier_dist_tbl,
  primary_share_tbl   = primary_share_tbl,
  primary_unique_share_tbl = primary_unique_share_tbl,
  secondary_share_tbl = secondary_share_tbl,
  too_high_share_tbl  = too_high_share_tbl,
  too_low_share_tbl   = too_low_share_tbl,
  vpi_summary_tbl     = global_obj$global_vpi,
  reporting_audit_summary_tbl = reporting_audit_summary_tbl
)

validation_tbl <- build_validation_summary(
  shared_tbl                 = shared_tbl,
  lookup_tbl                 = lookup_tbl,
  current_core_occupied_fate = current_core_occupied_fate,
  stats_tbl                  = stats_tbl,
  pixels_censored            = current_res$pixels_censored,
  valid_pixels               = valid_pixels_total
)

run_metadata <- tibble(
  metric = c("species_code", "species_label", "output_suffix", "occ_csv", "current_vars_dir",
             "output_dir", "n_focal_vars", "n_retained_predictors",
             "n_occurrences_total", "n_occurrences_excluded_incomplete",
             "n_occurrences_core_occupied", "n_occurrences_peripheral",
             "tier_breaks_method", "tier_break_resolved_1", "tier_break_resolved_2",
             "tier_break_resolved_3", "core_occupied_rule", "core_occupied_rule_value",
             "mahal_reference", "tier_break_reference",
             "core_partition_ridge_used", "landscape_mahal_ridge_used",
             "tail_scale_center", "shrinkage_method", "shrinkage_k",
             "z2_cap", "equality_tolerance", "random_seed",
             "bootstrap_iter", "bootstrap_frac", "bootstrap_replace",
             "memfrac", "todisk", "raster_compression",
             "reporting_audit_enabled", "reporting_contribution_coverage"
             ),
  value  = c(cfg$species_code,
             cfg$species_label,
             cfg$output_suffix,
             cfg$occ_csv,
             cfg$current_vars_dir,
             cfg$output_dir,
             length(cfg$focal_vars),
             length(vars_common),
             nrow(occ_env),
             anchor_obj$excluded_incomplete,
             nrow(occ_core_occupied),
             nrow(occ_peripheral),
             cfg$tier_breaks_method,
             round(as.numeric(tier_breaks_resolved[1]), 6),
             round(as.numeric(tier_breaks_resolved[2]), 6),
             round(as.numeric(tier_breaks_resolved[3]), 6),
             mahal_split$rule,
             as.numeric(mahal_split$rule_value),
             cfg$mahal_reference,
             cfg$tier_break_reference,
             paste(core_partition_ridge_values, collapse = " | "),
             current_mahal$ridge,
             cfg$tail_scale_center,
             cfg$shrinkage_method,
             cfg$shrinkage_k,
             cfg$z2_cap,
             cfg$equality_tolerance,
             cfg$random_seed,
             cfg$bootstrap_iter,
             cfg$bootstrap_frac,
             TRUE,
             cfg$memfrac,
             isTRUE(cfg$todisk),
             paste(cfg$raster_compression, collapse = " | "),
             isTRUE(cfg$write_reporting_audit),
             cfg$reporting_contribution_coverage)

) %>%
  mutate(value      = as.character(.data$value),
         section    = "run_metadata",
         subsection = "settings",
         scenario   = NA_character_,
         .before    = 1)

run_manifest_csv <- build_run_manifest_csv(run_metadata, validation_tbl)

model_reference_csv <- build_model_reference_csv(
  shared_tbl,
  lookup_tbl,
  stats_tbl,
  core_occupied_centroid_ref,
  threshold_tbl,
  current_core_occupied_vrs_tbl,
  predictor_attribution_tbl
)
occurrence_diagnostics_csv <- build_occurrence_diagnostics_csv(
  multimodality_tbl      = multimodality_tbl,
  core_occupied_bias_tbl = core_occupied_bias_tbl,
  stats_tbl              = stats_tbl
)
occurrence_partitions_csv <- build_occurrence_partitions_csv(occ_core_occupied, occ_peripheral)

write.csv(run_manifest_csv,           csv_out(out_paths, paste0(cfg$species_code, "_01_run_manifest")),           row.names = FALSE)
write.csv(model_reference_csv,        csv_out(out_paths, paste0(cfg$species_code, "_02_model_reference")),        row.names = FALSE)
write.csv(occurrence_diagnostics_csv, csv_out(out_paths, paste0(cfg$species_code, "_03_occurrence_diagnostics")), row.names = FALSE)
write.csv(occurrence_partitions_csv,  csv_out(out_paths, paste0(cfg$species_code, "_04_occurrence_partitions")),  row.names = FALSE)
write.csv(current_diagnostics_csv,    csv_out(out_paths, paste0(cfg$species_code, "_05_current_diagnostics")),    row.names = FALSE)

proximity_summary_csv <- build_proximity_summary_csv(
  global_vpi_row = global_obj$global_vpi
)
write.csv(proximity_summary_csv,
          csv_out(out_paths, paste0(cfg$species_code, "_06_proximity_summary")),
          row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(cfg$output_dir, paste0(cfg$species_code, "_sessionInfo.txt")))
dput(cfg, file = file.path(cfg$output_dir, paste0(cfg$species_code, "_configuration.R")))

script_file <- analysis_script_file
script_md5 <- if (!is.na(script_file) && file.exists(script_file)) {
  unname(tools::md5sum(script_file))
} else {
  warning(
    "The executing script file could not be resolved; set VERA_SCRIPT_FILE to write a code checksum."
  )
  NA_character_
}
write.csv(
  tibble(script_file = script_file, script_md5 = script_md5,
         analysis_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
  file.path(cfg$output_dir, paste0(cfg$species_code, "_code_manifest.csv")),
  row.names = FALSE
)

final_pixels_censored <- current_res$pixels_censored
final_landscape_mahal_ridge <- current_mahal$ridge
if (isTRUE(cfg$cleanup_temp_at_end)) {
  rm(addon_obj, current_res, current_mahal, current_class, current_class_factor,
     current_stack, vpi_rast)
  gc(verbose = FALSE)
  try(terra::tmpFiles(remove = TRUE), silent = TRUE)
  residual_temp <- list.files(out_paths$temp, full.names = TRUE, all.files = TRUE,
                              no.. = TRUE)
  if (length(residual_temp)) unlink(residual_temp, recursive = TRUE, force = TRUE)
}

public_files <- list.files(cfg$output_dir, recursive = TRUE, full.names = TRUE)
public_files <- public_files[file.info(public_files)$isdir %in% FALSE]
root_norm <- normalizePath(cfg$output_dir, winslash = "/", mustWork = TRUE)
files_norm <- normalizePath(public_files, winslash = "/", mustWork = TRUE)
temp_norm <- normalizePath(out_paths$temp, winslash = "/", mustWork = FALSE)
keep_public <- !startsWith(files_norm, paste0(temp_norm, "/"))
public_files <- public_files[keep_public]
files_norm <- files_norm[keep_public]
inventory <- tibble(
  relative_path = substring(files_norm, nchar(root_norm) + 2),
  size_bytes = as.numeric(file.info(public_files)$size),
  md5 = if (isTRUE(cfg$write_output_checksums)) {
    unname(tools::md5sum(public_files))
  } else {
    rep(NA_character_, length(public_files))
  }
)
write.csv(inventory,
          file.path(cfg$output_dir, paste0(cfg$species_code, "_output_inventory.csv")),
          row.names = FALSE)

cat("\n>>> AUTO-DIAGNOSTIC REPORT\n")
cat(sprintf("    Retained occurrences          : %d\n", retained_occ_n))
cat(sprintf("    Climatic core_occupied points : %d\n", nrow(occ_core_occupied)))
cat(sprintf("    Shared predictors             : %d\n", nrow(shared_tbl)))
cat(sprintf("    Core occupied climate rule    : %s (%.2f)\n",
            mahal_split$rule, mahal_split$rule_value))
cat(sprintf("    Core occupied threshold (D2 / Dist) : %.4f / %.4f\n",
            mahal_split$threshold_d2, mahal_split$threshold_distance))
cat(sprintf("    Mahalanobis ridge             : %s\n",
            paste(core_partition_ridge_values, collapse = " | ")))
cat(sprintf("    Landscape Mahalanobis ridge   : %g\n", final_landscape_mahal_ridge))
cat(sprintf("    Tier method                   : %s\n", cfg$tier_method))
cat(sprintf("    Tier breaks method            : %s\n", cfg$tier_breaks_method))
cat(sprintf("    Tier breaks resolved          : %s\n",
            paste(round(tier_breaks_resolved, 4), collapse = " | ")))
cat(sprintf("    Pixels censored               : %d\n", final_pixels_censored))
cat(sprintf("    Primary co-dominant pixels    : %.2f%%\n", primary_codominant_pixel_pct))

if (!is.null(tier_dist_tbl)) {
  cat(sprintf("    Tier distribution             : core=%.2f%% / mod_dep=%.2f%% / restr=%.2f%% / high_ex=%.2f%%\n",
              tier_dist_tbl$core_climate_pct, tier_dist_tbl$moderate_departure_pct,
              tier_dist_tbl$restriction_zone_pct, tier_dist_tbl$high_extrapolative_stress_pct))
}

cat(sprintf("    VPI summary                   : mean=%.4f | median=%.4f | sd=%.4f\n",
            global_obj$global_vpi$mean_vpi,
            global_obj$global_vpi$median_vpi,
            global_obj$global_vpi$sd_vpi))

cat(sprintf("    Outputs in                    : %s\n", cfg$output_dir))
cat("========================================================================\n")

invisible(list(
  species_code = cfg$species_code,
  species_label = cfg$species_label,
  output_dir = cfg$output_dir,
  predictor_profile = 19L,
  n_predictors = length(vars_common),
  validation = validation_tbl,
  vpi_summary = global_obj$global_vpi
))
}

# =============================================================================
# SECTION 13 -- SITTA KRUEPERI SINGLE EXECUTION (19-PREDICTOR PROFILE)
# =============================================================================
# One pass. Builds the raster stack from Bio1..Bio19 and writes all outputs
# to C:/VERA/Results/19/Skr_current.
# =============================================================================

# Ensure the results-root directory exists before the run creates any files.
output_root_19 <- "C:/VERA/Results/19"
if (!dir.exists(output_root_19)) {
  dir.create(output_root_19, recursive = TRUE, showWarnings = FALSE)
}

# Refresh the species_label from the CSV in case the user has a canonical
# spelling in the first row that overrides the configured fallback.
cfg$species_label <- detect_species_label(cfg$occ_csv, fallback = cfg$species_label)

validate_config(cfg)

cat("\n========================================================================\n")
cat("VERA TUTORIAL RUN  (19-predictor profile)\n")
cat(sprintf("Focal taxon        : %s (%s)\n", cfg$species_label, cfg$species_code))
cat(sprintf("Predictor profile  : 19 variables (Bio1..Bio19)\n"))
cat(sprintf("Occurrence CSV     : %s\n", cfg$occ_csv))
cat(sprintf("Predictor rasters  : %s\n", cfg$current_vars_dir))
cat(sprintf("Output directory   : %s\n", cfg$output_dir))
cat("========================================================================\n")

prepared_stack <- prepare_stack(cfg)

single_species_result <- tryCatch(
  run_vera_current(cfg, prepared_stack),
  error = function(e) {
    message(sprintf("\n!!! VERA failed for %s (19-predictor profile): %s",
                    cfg$species_code, conditionMessage(e)))
    list(
      species_code      = cfg$species_code,
      species_label     = cfg$species_label,
      output_dir        = cfg$output_dir,
      predictor_profile = 19L,
      n_predictors      = length(cfg$focal_vars),
      error             = conditionMessage(e)
    )
  }
)

rm(prepared_stack)
gc(verbose = FALSE)
try(terra::tmpFiles(remove = TRUE), silent = TRUE)

# -----------------------------------------------------------------------------
# Single-run summary written next to the results depot.
# -----------------------------------------------------------------------------
run_summary <- tibble::tibble(
  species_code      = single_species_result$species_code,
  species_label     = single_species_result$species_label,
  predictor_profile = single_species_result$predictor_profile,
  n_predictors      = single_species_result$n_predictors,
  output_dir        = single_species_result$output_dir,
  status            = if (is.null(single_species_result$error)) "success" else "failed",
  error             = if (is.null(single_species_result$error)) NA_character_ else single_species_result$error
)

summary_file <- file.path(
  output_root_19,
  paste0(cfg$species_code, "_run_summary_19.csv")
)
write.csv(run_summary, summary_file, row.names = FALSE)
cat(sprintf("\nRun summary written to: %s\n", summary_file))

cat("\n>>> SITTA KRUEPERI RUN SUMMARY (19-PREDICTOR PROFILE)\n")
print(run_summary)

invisible(run_summary)