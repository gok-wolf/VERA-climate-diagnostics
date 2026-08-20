#!/usr/bin/env Rscript

################################################################################
## VERA CORE SUITE -- MASTER SPATIAL AND OCCURRENCE RENDERER
## Focal taxon : Sitta krueperi (Skr)
## Dual-profile: 19- and 36-predictor VERA outputs rendered in the same run.
## Sources     : C:/VERA/Results/{19,36}/Skr_current
## Outputs     : C:/VERA/Results/{19,36}/Images/core_renderer_plum
## No shapefile is used; map bounds are derived from the raster extent.
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ggplot2)
  library(tidyterra)
  library(dplyr)
  library(shadowtext)
  library(scales)
  library(grid)
})

cat(">>> [1/5] Initializing canonical VERA core renderer (Sitta krueperi)...\n")

# =============================================================================
# GLOBAL SETTINGS
# =============================================================================

profiles_to_render <- c(19L, 36L)
species_code       <- "Skr"
species_label      <- "Sitta krueperi"

results_roots <- c(
  `19` = "C:/VERA/Results/19",
  `36` = "C:/VERA/Results/36"
)
run_subdir            <- "Skr_current"
image_subfolder_name  <- "core_renderer_plum"

# =============================================================================
# PUBLICATION PALETTES (shared across profiles)
# =============================================================================

pal_plum <- c("#FFF2FA", "#F3C4EC", "#D695CC", "#CC7CBF")
pal_vrs   <- pal_plum
pal_vpi   <- rev(pal_plum)
pal_mahal <- pal_plum

pal_codominance <- c(
  "Unique maximum"      = "#FFF2FA",
  "Co-dominant maximum" = "#CC7CBF"
)

pal_4tier <- c(
  "Core Climate"              = "#FFF2FA",
  "Moderate Departure"        = "#F3C4EC",
  "Restriction Zone"          = "#D695CC",
  "High Extrapolative Stress" = "#CC7CBF"
)

pal_partition <- c(
  "Core occupied" = "#F3C4EC",
  "Peripheral"    = "#CC7CBF"
)

line_partition <- c(
  "Core occupied" = "#B76AAA",
  "Peripheral"    = "#824177"
)

# Canonical 36-predictor palette (nominal; predictor identity is multihued).
pal_36 <- c(
  "PETCoQ" = "#0077FF", "PETDrQ" = "#B2EBF2", "PETWaQ" = "#00C4A7",
  "PETWeQ" = "#FFD700", "PETs"   = "#00897B", "AP"     = "#9E9D24",
  "AIT"    = "#0040FF", "CMI"    = "#9C27B0", "CNT"    = "#2E7D32",
  "EQ"     = "#795548", "GD0"    = "#8D6E99", "GD5"    = "#FE1111",
  "MaTCM"  = "#F9A825", "MeTCM"  = "#006064", "MeTWM"  = "#1565C0",
  "MiTWM"  = "#AD1457", "TI"     = "#4E342E", "Bio1"   = "#0A237E",
  "Bio2"   = "#00C853", "Bio3"   = "#7B4FA2", "Bio4"   = "#F50059",
  "Bio5"   = "#FF4081", "Bio6"   = "#80DEEA", "Bio7"   = "#558F2F",
  "Bio8"   = "#FF9F00", "Bio9"   = "#FF6F00", "Bio10"  = "#FFAA00",
  "Bio11"  = "#FF6D00", "Bio12"  = "#607D8B", "Bio13"  = "#37474F",
  "Bio14"  = "#18FFFF", "Bio15"  = "#6200EA", "Bio16"  = "#AA00FF",
  "Bio17"  = "#E53935", "Bio18"  = "#FF00CC", "Bio19"  = "#00AFD5"
)

# Predictor lookup ordering exactly matches the VERA pipeline's
# write_lookup_outputs() sequence, so integer stressor codes map to the correct
# variable name if a factor label happens to be missing.
predictor_sets <- list(
  `19` = paste0("Bio", 1:19),
  `36` = c("AIT", "AP", paste0("Bio", 1:19), "CMI", "CNT",
           "EQ", "GD0", "GD5", "MaTCM", "MeTCM", "MeTWM",
           "MiTWM", "PETCoQ", "PETDrQ", "PETs", "PETWaQ",
           "PETWeQ", "TI")
)

# =============================================================================
# SHARED HELPERS (theme and map builders)
# =============================================================================

theme_map <- theme_void(base_family = "Arial") +
  theme(
    plot.title        = element_text(size = 20, face = "plain", hjust = 0.5,
                                     colour = "black", margin = margin(b = 10)),
    legend.title      = element_text(size = 14, face = "plain", colour = "black"),
    legend.text       = element_text(size = 12, colour = "black"),
    legend.key.size   = unit(0.9, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.94), colour = NA),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.margin       = margin(4, 4, 4, 4)
  )

# `xlim_map` and `ylim_map` are set per-profile inside the loop and shared with
# add_map_geometry() via lexical scoping through a small render environment.
make_render_env <- function(template_raster) {
  ext_r <- terra::ext(template_raster)
  env <- new.env(parent = emptyenv())
  env$xlim_map <- c(ext_r$xmin, ext_r$xmax)
  env$ylim_map <- c(ext_r$ymin, ext_r$ymax)
  env$raster_crs <- terra::crs(template_raster)
  env
}

add_map_geometry <- function(p, renv) {
  p +
    coord_sf(xlim = renv$xlim_map, ylim = renv$ylim_map,
             expand = FALSE, datum = NA) +
    theme_map
}

# =============================================================================
# CONTINUOUS, STRESSOR, TIER, CO-DOMINANCE MAP BUILDERS
# =============================================================================

make_continuous_map <- function(tif_path, palette, legend_title,
                                limits = NULL, renv) {
  r <- rast(tif_path)
  names(r) <- "value"

  add_map_geometry(
    ggplot() +
      geom_spatraster(data = r, aes(fill = value), na.rm = TRUE) +
      scale_fill_gradientn(
        colours   = palette,
        values    = scales::rescale(seq_along(palette)),
        name      = legend_title,
        limits    = limits,
        oob       = scales::squish,
        na.value  = "transparent",
        guide     = guide_colourbar(
          barheight = unit(7.5, "cm"), barwidth = unit(0.75, "cm")
        )
      ),
    renv
  )
}

make_stressor_factor <- function(tif_path, v_names_profile) {
  idx_path <- sub("_factor\\.tif$", "_idx.tif", tif_path, ignore.case = TRUE)
  use_idx <- !identical(idx_path, tif_path) && file.exists(idx_path)
  r <- rast(if (use_idx) idx_path else tif_path)
  source_levels <- if (is.factor(r)) levels(r)[[1]] else NULL

  if (use_idx) {
    vals <- sort(unique(suppressWarnings(as.integer(values(r, mat = FALSE)))))
    vals <- vals[is.finite(vals) & vals >= 1L & vals <= length(v_names_profile)]
    if (!length(vals)) stop("No valid stressor codes were found in ", idx_path)
    r_factor <- as.factor(r)
    levels(r_factor) <- data.frame(value = vals, label = v_names_profile[vals])
  } else if (!is.null(source_levels) && ncol(source_levels) >= 2L) {
    ids    <- suppressWarnings(as.integer(source_levels[[1]]))
    labels <- as.character(source_levels[[2]])
    keep   <- is.finite(ids) & labels %in% v_names_profile
    if (!any(keep)) stop("No canonical predictor categories were found in ", tif_path)
    r_factor <- as.factor(r)
    levels(r_factor) <- data.frame(value = ids[keep], label = labels[keep])
  } else {
    vals <- sort(unique(suppressWarnings(as.integer(values(r, mat = FALSE)))))
    vals <- vals[is.finite(vals) & vals >= 1L & vals <= length(v_names_profile)]
    if (!length(vals)) stop("No valid stressor codes were found in ", tif_path)
    r_factor <- as.factor(r)
    levels(r_factor) <- data.frame(value = vals, label = v_names_profile[vals])
  }

  names(r_factor) <- "stressor"
  r_factor
}

make_stressor_map <- function(tif_path, legend_title, v_names_profile, renv) {
  r <- make_stressor_factor(tif_path, v_names_profile)
  labels_present <- as.character(levels(r)[[1]][[2]])
  add_map_geometry(
    ggplot() +
      geom_spatraster(data = r, aes(fill = stressor), na.rm = TRUE) +
      scale_fill_manual(
        values = pal_36[labels_present], name = legend_title,
        drop = FALSE, na.translate = FALSE, na.value = "transparent",
        guide = guide_legend(ncol = 1, override.aes = list(colour = NA))
      ),
    renv
  )
}

make_tier_map <- function(tif_path, renv) {
  r <- as.factor(rast(tif_path))
  levels(r) <- data.frame(value = 1:4, label = names(pal_4tier))
  names(r) <- "tier"
  add_map_geometry(
    ggplot() +
      geom_spatraster(data = r, aes(fill = tier), na.rm = TRUE) +
      scale_fill_manual(
        values = pal_4tier, name = "Climatic tier", drop = FALSE,
        na.translate = FALSE, na.value = "transparent"
      ),
    renv
  )
}

make_codominance_map <- function(tif_path, renv) {
  r <- rast(tif_path)
  r <- ifel(is.na(r), NA, ifel(r <= 1, 1, 2))
  r <- as.factor(r)
  levels(r) <- data.frame(value = 1:2, label = names(pal_codominance))
  names(r) <- "co_dominance"
  add_map_geometry(
    ggplot() +
      geom_spatraster(data = r, aes(fill = co_dominance), na.rm = TRUE) +
      scale_fill_manual(
        values = pal_codominance, name = "Primary maximum", drop = FALSE,
        na.translate = FALSE, na.value = "transparent"
      ),
    renv
  )
}

# =============================================================================
# COMPACT VERTICAL PREDICTOR LEGEND (2-column boxes)
# =============================================================================

extract_stressor_frequency <- function(tif_path, v_names_profile) {
  f <- terra::freq(rast(tif_path))
  if (is.null(f) || !nrow(f)) return(NULL)

  if ("label" %in% names(f)) {
    f$Variable <- as.character(f$label)
  } else if (is.character(f$value) || is.factor(f$value)) {
    f$Variable <- as.character(f$value)
  } else {
    idx <- suppressWarnings(as.integer(f$value))
    f$Variable <- ifelse(
      is.finite(idx) & idx >= 1 & idx <= length(v_names_profile),
      v_names_profile[idx], NA_character_
    )
  }

  f <- f[f$Variable %in% v_names_profile & is.finite(f$count), , drop = FALSE]
  if (!nrow(f)) return(NULL)
  f %>%
    transmute(
      Variable   = .data$Variable,
      count      = as.numeric(.data$count),
      Percentage = 100 * .data$count / sum(.data$count)
    ) %>%
    arrange(desc(.data$Percentage))
}

draw_compact_vertical_legend <- function(frequency_df, output_file,
                                         max_variables = 10L) {
  df <- head(frequency_df, max_variables)
  if (!nrow(df)) return(invisible(NULL))

  n_rows <- ceiling(nrow(df) / 2)
  df$x <- rep(c(1, 2), length.out = nrow(df))

  half_width  <- 0.49
  half_height <- 0.28
  row_step    <- 0.60

  row_centres <- rev(seq(from = 0, by = row_step, length.out = n_rows))
  df$y <- rep(row_centres, each = 2)[seq_len(nrow(df))]

  y_min  <- min(row_centres) - half_height
  y_max  <- max(row_centres) + half_height
  y_span <- y_max - y_min
  output_height <- 2 * y_span

  p <- ggplot(df) +
    geom_rect(
      aes(xmin = x - half_width, xmax = x + half_width,
          ymin = y - half_height, ymax = y + half_height,
          fill = Variable),
      colour = "black", linewidth = 0.9
    ) +
    geom_shadowtext(
      aes(x = x, y = y + 0.09, label = Variable),
      family = "Arial", colour = "white", bg.colour = "black",
      bg.r = 0.05, size = 9.4, fontface = "plain"
    ) +
    geom_shadowtext(
      aes(x = x, y = y - 0.09, label = sprintf("%.1f%%", Percentage)),
      family = "Arial", colour = "#FFF3A6", bg.colour = "black",
      bg.r = 0.05, size = 9.4, fontface = "plain"
    ) +
    scale_fill_manual(values = pal_36) +
    coord_cartesian(xlim = c(0.5, 2.5), ylim = c(y_min, y_max), expand = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  ggsave(output_file, plot = p, width = 4, height = output_height,
         units = "in", dpi = 400, bg = "white")
}

# =============================================================================
# DUAL-PROFILE LOOP
# =============================================================================

for (profile_value in profiles_to_render) {

  profile_key       <- as.character(profile_value)
  v_names_profile   <- predictor_sets[[profile_key]]

  results_root <- unname(results_roots[[profile_key]])
  run_dir      <- file.path(results_root, run_subdir)
  path_core    <- file.path(run_dir, "rasters", "core_occupied")
  path_diag    <- file.path(run_dir, "rasters", "optional_diagnostics")
  path_csv     <- file.path(run_dir, "csvs")

  # Sub-classified output under Images/<renderer_name>/
  path_img <- file.path(results_root, "Images", image_subfolder_name)
  if (!dir.exists(path_img)) dir.create(path_img, recursive = TRUE, showWarnings = FALSE)

  cat("\n========================================================================\n")
  cat(sprintf("CORE RENDERER  ::  Sitta krueperi  ::  %s-predictor profile\n",
              profile_key))
  cat(sprintf("Source run     : %s\n", run_dir))
  cat(sprintf("Image target   : %s\n", path_img))
  cat("========================================================================\n")

  # ---- Build the render environment from a template raster (raster extent) ----
  template_tif <- file.path(path_core, paste0(species_code, "_current_mean_vrs.tif"))
  if (!file.exists(template_tif)) {
    warning("Skipping profile ", profile_key,
            ": template raster missing (", template_tif, ")")
    next
  }
  template_raster <- rast(template_tif)
  renv <- make_render_env(template_raster)

  # ---- Filename builders and save helper ------------------------------------
  save_map_pair <- function(plot_object, output_stub, title_text) {
    withlegend_file <- file.path(path_img, paste0(output_stub, "_withlegend.png"))
    nolegend_file   <- file.path(path_img, paste0(output_stub, "_nolegend.png"))
    ggsave(
      withlegend_file, plot = plot_object + labs(title = title_text),
      width = 14, height = 12, units = "in", dpi = 400, bg = "white"
    )
    ggsave(
      nolegend_file,
      plot = plot_object + labs(title = NULL) +
        theme(legend.position = "none", plot.title = element_blank()),
      width = 14, height = 12, units = "in", dpi = 400, bg = "white"
    )
    cat("    Saved:", basename(withlegend_file), "and", basename(nolegend_file), "\n")
  }

  # ---- Canonical map manifest for this profile ------------------------------
  stub <- function(name) paste0(species_code, "_current_", name, ".tif")
  map_manifest <- list(
    mean_vrs = list(
      file = file.path(path_core, stub("mean_vrs")), type = "continuous",
      palette = pal_vrs, limits = c(0, 16), legend = "Mean VRS",
      title = "Mean Variable Restriction Score"
    ),
    max_vrs = list(
      file = file.path(path_core, stub("max_vrs")), type = "continuous",
      palette = pal_vrs, limits = c(0, 16), legend = "Maximum VRS",
      title = "Maximum Predictor-Specific Restriction"
    ),
    delta_vrs = list(
      file = file.path(path_diag, stub("delta_vrs")), type = "continuous",
      palette = pal_vrs, limits = c(0, 16), legend = expression(Delta * "VRS"),
      title = "Primary-Secondary Restriction Difference"
    ),
    vpi = list(
      file = file.path(path_core, stub("vpi")), type = "continuous",
      palette = pal_vpi, limits = c(0, 1), legend = "VPI",
      title = "Variable Proximity Index"
    ),
    mahal_distance = list(
      file = file.path(path_core, stub("mahal_distance")),
      type = "continuous", palette = pal_mahal, limits = NULL,
      legend = "Mahalanobis distance", title = "Covariance-Aware Mahalanobis Distance"
    ),
    four_tier_status = list(
      file = file.path(path_core, stub("four_tier_status")),
      type = "tier", title = "Empirical Climatic Tiers"
    ),
    primary_assigned = list(
      file = file.path(path_core, stub("primary_stressor_factor")),
      type = "stressor", legend = "Assigned Primary Stressor",
      title = "Assigned Primary Stressor"
    ),
    primary_unique = list(
      file = file.path(path_diag, stub("primary_unique_stressor_factor")),
      type = "stressor", legend = "Unique Primary Stressor",
      title = "Uniquely Dominant Primary Stressor"
    ),
    primary_codominance = list(
      file = file.path(path_diag, stub("primary_co_dominance_count")),
      type = "codominance", title = "Primary-Stressor Co-Dominance"
    ),
    secondary = list(
      file = file.path(path_diag, stub("secondary_stressor_factor")),
      type = "stressor", legend = "Secondary Stressor", title = "Secondary Stressor"
    ),
    too_high = list(
      file = file.path(path_diag, stub("too_high_factor")),
      type = "stressor", legend = "TooHigh Stressor",
      title = "Upper-Tail (TooHigh) Stressor"
    ),
    too_low = list(
      file = file.path(path_diag, stub("too_low_factor")),
      type = "stressor", legend = "TooLow Stressor",
      title = "Lower-Tail (TooLow) Stressor"
    )
  )

  cat(sprintf(">>> [2/5] Rendering canonical core maps (%s vars)...\n", profile_key))

  for (output_stub in names(map_manifest)) {
    s <- map_manifest[[output_stub]]
    if (!file.exists(s$file)) {
      warning("Skipping missing raster: ", s$file)
      next
    }
    cat("    +", output_stub, "\n")
    p <- switch(
      s$type,
      continuous  = make_continuous_map(s$file, s$palette, s$legend, s$limits, renv),
      stressor    = make_stressor_map(s$file, s$legend, v_names_profile, renv),
      tier        = make_tier_map(s$file, renv),
      codominance = make_codominance_map(s$file, renv),
      stop("Unknown map type: ", s$type)
    )
    save_map_pair(p, paste0("vera_core_", output_stub), s$title)
    rm(p)
    gc(verbose = FALSE)
  }

  # ---- Pixel counts + compact vertical legends ------------------------------
  cat(sprintf(">>> [3/5] Writing stressor shares and compact vertical legends (%s)...\n",
              profile_key))

  stressor_legend_manifest <- list(
    primary_assigned = file.path(path_core, stub("primary_stressor_factor")),
    primary_unique   = file.path(path_diag, stub("primary_unique_stressor_factor")),
    secondary        = file.path(path_diag, stub("secondary_stressor_factor")),
    too_high         = file.path(path_diag, stub("too_high_factor")),
    too_low          = file.path(path_diag, stub("too_low_factor"))
  )

  for (legend_name in names(stressor_legend_manifest)) {
    tif_path <- stressor_legend_manifest[[legend_name]]
    if (!file.exists(tif_path)) next
    frequency_df <- extract_stressor_frequency(tif_path, v_names_profile)
    if (is.null(frequency_df)) next

    write.csv(
      frequency_df,
      file.path(path_img, paste0("vera_pixel_counts_", legend_name, ".csv")),
      row.names = FALSE
    )
    draw_compact_vertical_legend(
      frequency_df,
      file.path(path_img,
                paste0("vera_vertical_legend_", legend_name, "_compact.png"))
    )
  }

  # ---- Occurrence core / peripheral diagnostics -----------------------------
  cat(sprintf(">>> [4/5] Rendering occurrence partition diagnostics (%s)...\n",
              profile_key))

  occurrence_file <- file.path(path_csv,
                               paste0(species_code, "_04_occurrence_partitions.csv"))
  if (!file.exists(occurrence_file)) {
    warning("Occurrence partition CSV was not found: ", occurrence_file)
  } else {
    occurrence_df <- read.csv(occurrence_file, stringsAsFactors = FALSE) %>%
      filter(.data$partition %in% c("core_occupied", "peripheral")) %>%
      mutate(
        partition_label = factor(
          .data$partition,
          levels = c("core_occupied", "peripheral"),
          labels = c("Core occupied", "Peripheral")
        )
      )

    core_threshold <- unique(occurrence_df$core_threshold_distance)
    core_threshold <- core_threshold[is.finite(core_threshold)][1]

    p_density <- ggplot(
      occurrence_df,
      aes(x = mahal_distance, fill = partition_label, colour = partition_label)
    ) +
      geom_density(alpha = 0.70, linewidth = 0.85, adjust = 1) +
      geom_vline(xintercept = core_threshold, colour = "black",
                 linetype = "dashed", linewidth = 0.85) +
      scale_fill_manual(values = pal_partition, name = "Partition") +
      scale_colour_manual(values = line_partition, name = "Partition") +
      labs(x = "Occurrence-level Mahalanobis distance", y = "Density") +
      theme_classic(base_family = "Arial") +
      theme(
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA)
      )

    ggsave(
      file.path(path_img, "vera_occurrence_core_peripheral_density_withlegend.png"),
      p_density + labs(title = "Occurrence-Level Climatic Partition"),
      width = 8, height = 5, units = "in", dpi = 400, bg = "white"
    )
    ggsave(
      file.path(path_img, "vera_occurrence_core_peripheral_density_nolegend.png"),
      p_density + labs(title = NULL) +
        theme(legend.position = "none", plot.title = element_blank()),
      width = 8, height = 5, units = "in", dpi = 400, bg = "white"
    )

    # ---- Occurrence-point map -----------------------------------------------
    # Longitude / Latitude column names are inherited from the input CSV, so
    # detect them dynamically.
    lon_candidates <- c("decimalLongitude", "Longitude", "longitude", "lon", "x")
    lat_candidates <- c("decimalLatitude",  "Latitude",  "latitude",  "lat", "y")
    lon_col <- lon_candidates[lon_candidates %in% names(occurrence_df)][1]
    lat_col <- lat_candidates[lat_candidates %in% names(occurrence_df)][1]

    if (!is.na(lon_col) && !is.na(lat_col)) {
      points_sf <- st_as_sf(
        occurrence_df, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE
      )

      # A pale-grey rectangle at the raster extent gives the domain a soft
      # geographic footprint without an external boundary shapefile and without
      # requiring the ggnewscale package.
      p_points <- ggplot() +
        annotate(
          "rect",
          xmin = renv$xlim_map[1], xmax = renv$xlim_map[2],
          ymin = renv$ylim_map[1], ymax = renv$ylim_map[2],
          fill = "#F7F7F5", colour = "grey55", linewidth = 0.30
        ) +
        geom_sf(
          data = points_sf, aes(colour = partition_label),
          alpha = 0.75, size = 1.3, stroke = 0
        ) +
        scale_colour_manual(values = line_partition, name = "Partition") +
        coord_sf(xlim = renv$xlim_map, ylim = renv$ylim_map,
                 expand = FALSE, datum = NA) +
        theme_map + theme(legend.position = "right")

      save_map_pair(
        p_points, "vera_occurrence_partition_map",
        "Occurrence-Level Climatic Partition"
      )
    } else {
      warning("Could not detect longitude/latitude columns in ", occurrence_file)
    }
  }

  cat(sprintf(">>> [5/5] Core renderer completed for %s-predictor profile.\n",
              profile_key))
}

cat("\n========================================================================\n")
cat(">>> ALL CORE-RENDERER PROFILES COMPLETED\n")
cat(">>> Focal taxon : ", species_label, " (", species_code, ")\n", sep = "")
cat(">>> Profiles    : 19 + 36\n")
cat(">>> Image roots :\n")
for (pk in names(results_roots)) {
  cat("      ", file.path(results_roots[[pk]], "Images", image_subfolder_name), "\n", sep = "")
}
cat("========================================================================\n")