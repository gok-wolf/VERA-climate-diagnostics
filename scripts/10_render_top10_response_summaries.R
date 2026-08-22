#!/usr/bin/env Rscript

# =============================================================================
# OPTIONAL VERA TOP-10 RESPONSE SUMMARY RENDERER
# Sitta krueperi; supports the 19- and 36-predictor profiles in one script.
#
# Products per profile
#   1. Ten native-unit occurrence/background density panels in one gallery
#   2. One aligned-optimum asymmetric-standardisation ridge gallery
#
# These figures summarise completed VERA outputs. They do not modify anchor
# statistics, Primary Stressor assignments, VRS surfaces, or tier products.
# Each profile output directory also receives README_10.txt.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(patchwork)
  library(ggridges)
})

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================

cfg <- list(
  species_code = "Skr",
  species_label = "Sitta krueperi",
  profiles_to_run = c(19L, 36L),
  profile_dirs = c(
    `19` = "C:/VERA/Results/19/Skr_current",
    `36` = "C:/VERA/Results/36/Skr_current"
  ),
  raster_dirs = c(
    `19` = "C:/VERA/Variables/Turkiye_Current_Bioclimatics_and_Envirems",
    `36` = "C:/VERA/Variables/Turkiye_Current_Bioclimatics_and_Envirems"
  ),
  occurrence_file = "C:/VERA/Occurrences/Edited/Sitta_krueperi.csv",
  top_n = 10L,
  background_sample_size = 50000L,
  standardized_display_limit = 4,
  output_subdir = file.path("Images", "optional_response_summaries"),
  dpi = 400,
  overwrite = TRUE,
  random_seed = 20260809L
)

font_main <- "Arial"

palette <- list(
  bg = "#FFEBB5",
  occ = "#6AB4FF",
  lower = "#FF8587",
  upper = "#9AEBB1",
  ref = "#FF0000",
  text_dark = "#000000",
  grid = "#E0E0E0"
)

long_names <- c(
  Bio1 = "Annual Mean Temperature",
  Bio2 = "Mean Diurnal Range",
  Bio3 = "Isothermality",
  Bio4 = "Temperature Seasonality",
  Bio5 = "Max Temperature of Warmest Month",
  Bio6 = "Min Temperature of Coldest Month",
  Bio7 = "Temperature Annual Range",
  Bio8 = "Mean Temperature of Wettest Quarter",
  Bio9 = "Mean Temperature of Driest Quarter",
  Bio10 = "Mean Temperature of Warmest Quarter",
  Bio11 = "Mean Temperature of Coldest Quarter",
  Bio12 = "Annual Precipitation",
  Bio13 = "Precipitation of Wettest Month",
  Bio14 = "Precipitation of Driest Month",
  Bio15 = "Precipitation Seasonality",
  Bio16 = "Precipitation of Wettest Quarter",
  Bio17 = "Precipitation of Driest Quarter",
  Bio18 = "Precipitation of Warmest Quarter",
  Bio19 = "Precipitation of Coldest Quarter",
  AIT = "Aridity Index (Thornthwaite)",
  AP = "Annual Potential Evapotranspiration",
  CMI = "Climatic Moisture Index",
  CNT = "Continentality",
  EQ = "Emberger Pluviothermic Quotient",
  GD0 = "Growing Degree Days (base 0 C)",
  GD5 = "Growing Degree Days (base 5 C)",
  MaTCM = "Maximum Temperature of Coldest Month",
  MeTCM = "Mean Temperature of Coldest Month",
  MeTWM = "Mean Temperature of Warmest Month",
  MiTWM = "Minimum Temperature of Warmest Month",
  PETCoQ = "PET of Coldest Quarter",
  PETDrQ = "PET of Driest Quarter",
  PETs = "PET Seasonality",
  PETWaQ = "PET of Warmest Quarter",
  PETWeQ = "PET of Wettest Quarter",
  TI = "Thermicity Index"
)

dir_make <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_lines_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, useBytes = TRUE)
  invisible(path)
}

read_required_csv <- function(path) {
  if (!file.exists(path)) stop("Required file was not found: ", path)
  read_csv(path, show_col_types = FALSE)
}

resolve_primary_counts_file <- function(run_dir) {
  filename <- "vera_pixel_counts_primary_assigned.csv"
  candidates <- c(
    file.path(dirname(run_dir), "Images", "core_renderer_plum", filename),
    file.path(run_dir, "Images", "core_renderer_plum", filename)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(hit[1])
  candidates[1]
}

top_primary_variables <- function(counts_file, top_n) {
  counts <- read_required_csv(counts_file)
  if (ncol(counts) < 2L) stop("Primary-count table requires at least two columns.")
  counts <- counts %>% arrange(desc(suppressWarnings(as.numeric(.[[2]]))))
  variables <- as.character(counts[[1]][seq_len(min(top_n, nrow(counts)))])
  unique(variables[!is.na(variables) & nzchar(variables)])
}

get_anchor_value <- function(model_ref, variable, metric) {
  value <- model_ref$value[
    model_ref$section == "anchor_statistics" &
      model_ref$subject == variable &
      model_ref$metric == metric
  ]
  if (!length(value)) return(NA_real_)
  suppressWarnings(as.numeric(value[1]))
}

display_name <- function(variable, separator = "  —  ") {
  label <- unname(long_names[variable])
  if (!length(label) || is.na(label) || !nzchar(label)) return(variable)
  paste0(variable, separator, label)
}

coordinate_columns <- function(data) {
  lon <- grep(
    "(?i)^(lon|longitude|decimallongitude|x)$",
    names(data), value = TRUE, perl = TRUE
  )[1]
  lat <- grep(
    "(?i)^(lat|latitude|decimallatitude|y)$",
    names(data), value = TRUE, perl = TRUE
  )[1]
  list(lon = lon, lat = lat)
}

occurrence_values <- function(variable, occurrence_df, raster) {
  if (variable %in% names(occurrence_df)) {
    values <- suppressWarnings(as.numeric(occurrence_df[[variable]]))
    return(values[is.finite(values)])
  }

  coordinates <- coordinate_columns(occurrence_df)
  if (is.na(coordinates$lon) || is.na(coordinates$lat)) {
    stop("Occurrence file lacks both predictor values and coordinate columns.")
  }
  points <- vect(
    occurrence_df,
    geom = c(coordinates$lon, coordinates$lat),
    crs = "EPSG:4326"
  )
  if (!is.na(crs(raster)) && !same.crs(points, raster)) {
    points <- project(points, crs(raster))
  }
  values <- suppressWarnings(as.numeric(extract(raster, points)[[2]]))
  values[is.finite(values)]
}

asymmetric_z <- function(x, mu, sigma_lower, sigma_upper) {
  ifelse(
    x < mu,
    (x - mu) / sigma_lower,
    (x - mu) / sigma_upper
  )
}

prepare_profile_data <- function(profile, run_dir, raster_dir) {
  model_file <- file.path(
    run_dir, "csvs", paste0(cfg$species_code, "_02_model_reference.csv")
  )
  counts_file <- resolve_primary_counts_file(run_dir)
  required <- c(model_file, counts_file, cfg$occurrence_file)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Required files were not found:\n", paste(missing, collapse = "\n"))
  }

  model_ref <- read_required_csv(model_file)
  occurrence_df <- read_required_csv(cfg$occurrence_file)
  top_variables <- top_primary_variables(counts_file, cfg$top_n)
  background_records <- list()
  occurrence_records <- list()
  metadata_records <- list()

  for (variable in top_variables) {
    raster_file <- file.path(raster_dir, paste0(variable, ".tif"))
    if (!file.exists(raster_file)) {
      warning("Raster missing for ", variable, "; predictor skipped.",
              call. = FALSE)
      next
    }

    mu <- get_anchor_value(model_ref, variable, "mu")
    sigma_lower <- get_anchor_value(model_ref, variable, "sd_lower")
    sigma_upper <- get_anchor_value(model_ref, variable, "sd_upper")
    if (!is.finite(mu) || !is.finite(sigma_lower) || sigma_lower <= 0 ||
        !is.finite(sigma_upper) || sigma_upper <= 0) {
      warning("Invalid anchor statistics for ", variable, "; skipped.",
              call. = FALSE)
      next
    }

    cat(">>> Profile ", profile, ": preparing ", variable, "\n", sep = "")
    raster <- rast(raster_file)
    sampled <- spatSample(
      raster, size = cfg$background_sample_size, method = "regular",
      na.rm = TRUE, as.df = TRUE
    )
    background <- suppressWarnings(as.numeric(sampled[[1]]))
    background <- background[is.finite(background)]
    occupied <- occurrence_values(variable, occurrence_df, raster)
    if (length(background) < 2L || length(occupied) < 2L) {
      warning("Insufficient values for ", variable, "; skipped.",
              call. = FALSE)
      next
    }

    background_records[[variable]] <- tibble(
      variable = variable,
      value = background,
      zscore = asymmetric_z(background, mu, sigma_lower, sigma_upper)
    )
    occurrence_records[[variable]] <- tibble(
      variable = variable,
      value = occupied,
      zscore = asymmetric_z(occupied, mu, sigma_lower, sigma_upper)
    )
    metadata_records[[variable]] <- tibble(
      variable = variable,
      mu = mu,
      sigma_lower = sigma_lower,
      sigma_upper = sigma_upper
    )
  }

  background <- bind_rows(background_records)
  occupied <- bind_rows(occurrence_records)
  metadata <- bind_rows(metadata_records)
  retained <- top_variables[top_variables %in% metadata$variable]
  if (!length(retained)) stop("No predictors could be prepared for profile ", profile)

  list(
    profile = profile,
    variables = retained,
    background = background,
    occupied = occupied,
    metadata = metadata
  )
}

theme_density_panel <- function() {
  theme_minimal(base_family = font_main, base_size = 18) +
    theme(
      panel.grid.major = element_line(colour = palette$grid, linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = palette$grid, fill = NA, linewidth = 0.5
      ),
      axis.text = element_text(
        colour = palette$text_dark, size = 13,
        face = "plain", family = font_main
      ),
      axis.title = element_blank(),
      plot.title = element_text(
        face = "plain", size = 16, hjust = 0.5,
        family = font_main, colour = palette$text_dark,
        margin = margin(b = 10)
      ),
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(12, 16, 12, 16)
    )
}

render_density_gallery <- function(data, output_dir) {
  plots <- list()
  for (variable in data$variables) {
    meta <- data$metadata[data$metadata$variable == variable, , drop = FALSE]
    background <- data$background[data$background$variable == variable, , drop = FALSE]
    occupied <- data$occupied[data$occupied$variable == variable, , drop = FALSE]
    x_limits <- range(
      c(
        quantile(background$value, c(0.001, 0.999),
                 na.rm = TRUE, names = FALSE),
        occupied$value
      ),
      na.rm = TRUE
    )

    plots[[variable]] <- ggplot() +
      annotate("rect", xmin = meta$mu - meta$sigma_lower, xmax = meta$mu,
               ymin = -Inf, ymax = Inf,
               fill = palette$lower, alpha = 0.30) +
      annotate("rect", xmin = meta$mu, xmax = meta$mu + meta$sigma_upper,
               ymin = -Inf, ymax = Inf,
               fill = palette$upper, alpha = 0.30) +
      geom_density(data = background, aes(.data$value),
                   fill = palette$bg, colour = NA, alpha = 0.65) +
      geom_density(data = occupied, aes(.data$value),
                   fill = palette$occ, colour = palette$text_dark,
                   alpha = 0.72, linewidth = 0.35) +
      geom_vline(xintercept = meta$mu, colour = palette$ref,
                 linewidth = 0.5) +
      coord_cartesian(xlim = x_limits, expand = FALSE) +
      labs(title = display_name(variable)) +
      theme_density_panel()
  }

  combined <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = paste0(
        cfg$species_label, " — top ", length(plots),
        " native-unit density profiles (", data$profile, " predictors)"
      ),
      theme = theme(
        plot.title = element_text(
          family = font_main, face = "plain", size = 22,
          hjust = 0.5, colour = palette$text_dark,
          margin = margin(b = 24)
        ),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    )

  output_file <- file.path(
    output_dir,
    paste0(cfg$species_code, "_", data$profile,
           "_Top", length(plots), "_Density_Profiles_400DPI.png")
  )
  if (!cfg$overwrite && file.exists(output_file)) {
    stop("Output exists and overwrite is FALSE: ", output_file)
  }
  ggsave(output_file, combined, width = 14, height = 17,
         units = "in", dpi = cfg$dpi, bg = "white")
  output_file
}

render_aligned_ridges <- function(data, output_dir) {
  z_limit <- cfg$standardized_display_limit
  variable_order <- rev(data$variables)
  background <- data$background %>%
    filter(is.finite(.data$zscore), abs(.data$zscore) <= z_limit) %>%
    mutate(variable = factor(.data$variable, levels = variable_order))
  occupied <- data$occupied %>%
    filter(is.finite(.data$zscore), abs(.data$zscore) <= z_limit) %>%
    mutate(variable = factor(.data$variable, levels = variable_order))
  labels <- setNames(
    vapply(
      variable_order,
      function(variable) display_name(variable, separator = "\n"),
      character(1)
    ),
    variable_order
  )

  plot <- ggplot() +
    annotate("rect", xmin = -1, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = palette$lower, alpha = 0.28) +
    annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf,
             fill = palette$upper, alpha = 0.28) +
    geom_density_ridges(
      data = background,
      aes(x = .data$zscore, y = .data$variable,
          height = after_stat(ndensity)),
      stat = "density", fill = palette$bg, colour = NA,
      alpha = 0.70, scale = 1.25, rel_min_height = 0.005
    ) +
    geom_density_ridges(
      data = occupied,
      aes(x = .data$zscore, y = .data$variable,
          height = after_stat(ndensity)),
      stat = "density", fill = palette$occ, colour = palette$text_dark,
      alpha = 0.85, linewidth = 0.35,
      scale = 1.25, rel_min_height = 0.005
    ) +
    geom_vline(xintercept = 0, colour = palette$ref, linewidth = 0.5) +
    scale_x_continuous(limits = c(-z_limit, z_limit), expand = c(0, 0)) +
    scale_y_discrete(labels = labels, expand = expansion(add = c(0.4, 1.4))) +
    labs(
      x = "Asymmetrically standardized departure from the occupied reference",
      y = NULL,
      title = paste0(
        cfg$species_label, " — aligned-optimum profiles (",
        data$profile, " predictors)"
      )
    ) +
    theme_minimal(base_family = font_main, base_size = 15) +
    theme(
      plot.title = element_text(
        family = font_main, face = "plain", size = 20,
        hjust = 0.5, colour = palette$text_dark,
        margin = margin(b = 22)
      ),
      axis.title.x = element_text(
        family = font_main, face = "plain", size = 15,
        colour = palette$text_dark, margin = margin(t = 14)
      ),
      axis.text.x = element_text(
        family = font_main, face = "plain", size = 12,
        colour = palette$text_dark
      ),
      axis.text.y = element_text(
        family = font_main, face = "plain", size = 14,
        colour = palette$text_dark, lineheight = 0.95
      ),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.line.x = element_line(colour = palette$grid, linewidth = 0.4),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(24, 32, 24, 24),
      legend.position = "none"
    )

  output_file <- file.path(
    output_dir,
    paste0(cfg$species_code, "_", data$profile,
           "_Top", length(data$variables),
           "_Aligned_Optimum_Ridges_400DPI.png")
  )
  if (!cfg$overwrite && file.exists(output_file)) {
    stop("Output exists and overwrite is FALSE: ", output_file)
  }
  ggsave(output_file, plot, width = 13, height = 12,
         units = "in", dpi = cfg$dpi, bg = "white")
  output_file
}

render_profile <- function(profile) {
  key <- as.character(profile)
  run_dir <- cfg$profile_dirs[[key]]
  raster_dir <- cfg$raster_dirs[[key]]
  if (is.null(run_dir) || is.null(raster_dir)) {
    stop("Missing run or raster directory for profile ", key)
  }
  if (!dir.exists(run_dir)) {
    warning("Profile directory was not found and was skipped: ", run_dir,
            call. = FALSE)
    return(invisible(NULL))
  }

  output_dir <- dir_make(file.path(run_dir, cfg$output_subdir))
  data <- prepare_profile_data(profile, run_dir, raster_dir)
  density_file <- render_density_gallery(data, output_dir)
  ridge_file <- render_aligned_ridges(data, output_dir)

  anchor_file <- file.path(
    output_dir,
    paste0(cfg$species_code, "_", profile,
           "_Top", length(data$variables), "_Summary_Anchors.csv")
  )
  write_csv(
    data$metadata %>%
      mutate(
        species_code = cfg$species_code,
        species_label = cfg$species_label,
        profile = profile,
        display_order = match(.data$variable, data$variables),
        .before = 1
      ),
    anchor_file
  )

  readme_file <- file.path(output_dir, "README_10.txt")
  model_file <- file.path(
    run_dir, "csvs", paste0(cfg$species_code, "_02_model_reference.csv")
  )
  counts_file <- resolve_primary_counts_file(run_dir)
  readme_text <- c(
    "VERA OPTIONAL TOP RESPONSE-SUMMARY PACKAGE",
    "=============================================",
    "",
    paste0("Species: ", cfg$species_label, " (", cfg$species_code, ")"),
    paste0("Predictor profile: ", profile),
    "Renderer: 10_render_top10_response_summaries.R",
    "Status: OPTIONAL summary and publication-oriented product",
    "",
    "PURPOSE",
    "-------",
    "This directory contains compact multi-predictor summaries for the most",
    "frequently assigned Primary Stressors in the completed VERA profile.",
    "Script 10 reads existing VERA calibration outputs, occurrence data, and",
    "predictor rasters. It does not recalculate or modify canonical VERA results.",
    "",
    "INPUTS",
    "------",
    paste0("- ", basename(model_file)),
    paste0("- ", basename(counts_file)),
    paste0("- Occurrence table: ", normalizePath(
      cfg$occurrence_file, winslash = "/", mustWork = FALSE
    )),
    paste0("- Predictor GeoTIFF directory: ", normalizePath(
      raster_dir, winslash = "/", mustWork = FALSE
    )),
    "",
    "PREDICTOR SELECTION",
    "-------------------",
    paste0("Requested Top N: ", cfg$top_n),
    paste0("Retained predictors: ", paste(data$variables, collapse = ", ")),
    "Predictors are ranked using the existing Primary-assignment pixel-count",
    "table. Missing rasters or invalid anchor statistics may reduce the number",
    "retained. Assignment frequency does not establish causal importance.",
    "",
    "OUTPUT GUIDE",
    "------------",
    paste0("Native-unit density gallery: ", basename(density_file)),
    "- Compares sampled landscape-background and occurrence distributions in",
    "  each predictor's native units.",
    "- The red line is mu; shaded intervals represent the lower and upper",
    "  occurrence-derived scales used by VERA.",
    "",
    paste0("Aligned-reference ridge gallery: ", basename(ridge_file)),
    "- Places predictor distributions on their asymmetric standardized axes so",
    "  that the occurrence-derived reference centre is aligned at zero.",
    "- Negative and positive values represent lower- and upper-side departure",
    "  after scaling by sigma_L and sigma_U, respectively.",
    paste0("- The displayed standardized range is limited to +/-",
           cfg$standardized_display_limit, " for visual comparison."),
    "",
    paste0("Anchor summary CSV: ", basename(anchor_file)),
    "- Records mu, sigma_L, sigma_U, profile identity, and display order for the",
    "  predictors retained in the two galleries.",
    "",
    "READING RULES",
    "-------------",
    "1. The phrase aligned optimum is graphical shorthand for aligning the",
    "   occurrence-derived reference centres. It does not demonstrate a",
    "   physiological optimum, fitness maximum, or causal ecological optimum.",
    "2. Density overlap describes distributions in the selected data; it does",
    "   not estimate suitability, occurrence probability, or model accuracy.",
    "3. Standardized ridge shapes support visual comparison across predictors,",
    "   but do not make different variables ecologically interchangeable.",
    "4. The display limit is graphical clipping, not a VRS cap or biological",
    "   threshold.",
    "5. The 19- and 36-predictor summaries are profile-specific. Differences may",
    "   reflect predictor composition and should not be read as validation or",
    "   proof that one profile is superior.",
    "6. These products are optional and should be selected according to the",
    "   research question, Supplementary Information plan, and journal policy.",
    "",
    "REPRODUCIBILITY",
    "---------------",
    paste0("Background sample size: ", cfg$background_sample_size),
    paste0("Random seed: ", cfg$random_seed),
    paste0("Standardized display limit: +/-", cfg$standardized_display_limit),
    "",
    "REPORTING NOTE",
    "--------------",
    "Report predictor provenance, occurrence preparation, anchor statistics,",
    "sampling settings, display limits, software version, and predictor-selection",
    "rules in the manuscript or Supplementary Information.",
    ""
  )
  write_lines_utf8(readme_text, readme_file)

  cat(">>> Optional Top-10 summaries completed:\n")
  cat("    ", density_file, "\n", sep = "")
  cat("    ", ridge_file, "\n", sep = "")
  invisible(c(density_file, ridge_file))
}

set.seed(cfg$random_seed)
profiles <- unique(as.integer(cfg$profiles_to_run))
if (!length(profiles) || any(!profiles %in% c(19L, 36L))) {
  stop("profiles_to_run must contain 19L, 36L, or both.")
}
for (profile in profiles) render_profile(profile)

cat("\n========================================================================\n")
cat("OPTIONAL TOP-10 RESPONSE SUMMARY RENDERER COMPLETED\n")
cat("Species: ", cfg$species_label, "\n", sep = "")
cat("Profiles requested: ", paste(profiles, collapse = ", "), "\n", sep = "")
cat("========================================================================\n")
