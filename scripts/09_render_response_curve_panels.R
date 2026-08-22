#!/usr/bin/env Rscript

# =============================================================================
# OPTIONAL VERA RESPONSE-CURVE PANEL RENDERER
# Sitta krueperi; supports the 19- and 36-predictor profiles in one script.
#
# For each requested profile, this script identifies the most frequently
# assigned Primary Stressors, exports their anchor annotations, and renders
# native-unit occurrence/background densities and asymmetric Z2 transforms.
# These are explanatory graphics; no canonical VERA result is recalculated.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
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
  top_n = 6L,
  background_sample_size = 100000L,
  z2_cap = 16,
  styles_to_run = c("Original_Enhanced", "Earth_Tone"),
  modes_to_run = c("Labeled", "Unlabeled"),
  output_subdir = file.path("Images", "optional_response_curves"),
  dpi = 400,
  random_seed = 20260809L,
  overwrite = TRUE
)

font_main <- "Arial"

styles <- list(
  Original_Enhanced = list(
    bg = "#FFEBB5", occ = "#6AB4FF",
    lower = "#FF8587", upper = "#9AEBB1",
    ref = "#FF0000", curve = "#1C2529",
    text_dark = "#000000", grid = "#E0E0E0", cap = "#555555"
  ),
  Earth_Tone = list(
    bg = "#E8DCC8", occ = "#8CA88E",
    lower = "#C68F7B", upper = "#A4C3A2",
    ref = "#A64B39", curve = "#3E362E",
    text_dark = "#3E362E", grid = "#DCCFBF", cap = "#7D7469"
  )
)

dir_make <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_required_csv <- function(path) {
  if (!file.exists(path)) stop("Required file was not found: ", path)
  read_csv(path, show_col_types = FALSE)
}

theme_vera_response <- function(is_labeled, palette) {
  result <- theme_minimal(base_family = font_main, base_size = 13) +
    theme(
      panel.grid.major = element_line(colour = palette$grid, linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = palette$text_dark, linewidth = 0.5),
      axis.ticks = element_line(colour = palette$text_dark, linewidth = 0.5),
      axis.text = element_text(
        colour = palette$text_dark, family = font_main,
        face = "plain", size = 12
      ),
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )

  if (is_labeled) {
    result + theme(
      axis.title.x = element_text(
        family = font_main, face = "plain", size = 13,
        margin = margin(t = 10)
      ),
      axis.title.y = element_text(
        family = font_main, face = "plain", size = 13,
        margin = margin(r = 10)
      ),
      plot.title = element_text(
        family = font_main, face = "plain", size = 15,
        hjust = 0.5, margin = margin(b = 14)
      ),
      plot.subtitle = element_text(
        family = font_main, face = "plain", size = 11,
        colour = "grey35", hjust = 0.5, margin = margin(b = 8)
      )
    )
  } else {
    result + theme(axis.title = element_blank(), plot.title = element_blank(),
                   plot.subtitle = element_blank())
  }
}

get_anchor_value <- function(model_ref, variable, metric, numeric = TRUE) {
  value <- model_ref$value[
    model_ref$section == "anchor_statistics" &
      model_ref$subject == variable &
      model_ref$metric == metric
  ]
  if (length(value) != 1L) {
    stop("Expected one anchor value for ", variable, " / ", metric)
  }
  if (numeric) suppressWarnings(as.numeric(value[1])) else as.character(value[1])
}

top_primary_variables <- function(counts_file, top_n) {
  counts <- read_required_csv(counts_file)
  if (ncol(counts) < 2L) stop("Primary-count table requires at least two columns.")
  counts <- counts %>% arrange(desc(suppressWarnings(as.numeric(.[[2]]))))
  variables <- as.character(counts[[1]][seq_len(min(top_n, nrow(counts)))])
  unique(variables[!is.na(variables) & nzchar(variables)])
}

sample_background <- function(raster_file, variable, sample_size) {
  if (!file.exists(raster_file)) stop("Predictor raster was not found: ", raster_file)
  raster <- rast(raster_file)
  sampled <- spatSample(
    raster, size = sample_size, method = "regular",
    na.rm = TRUE, as.df = TRUE
  )
  values <- suppressWarnings(as.numeric(sampled[[1]]))
  values[is.finite(values)]
}

asymmetric_z2 <- function(x, mu, sigma_lower, sigma_upper, z2_cap) {
  raw <- ifelse(
    x < mu,
    ((x - mu) / sigma_lower)^2,
    ((x - mu) / sigma_upper)^2
  )
  pmin(raw, z2_cap)
}

prepare_predictor <- function(variable, occurrence_df, model_ref, raster_dir) {
  if (!variable %in% names(occurrence_df)) {
    stop("Occurrence-partition table lacks predictor column: ", variable)
  }
  occurrence_values <- suppressWarnings(as.numeric(occurrence_df[[variable]]))
  occurrence_values <- occurrence_values[is.finite(occurrence_values)]
  if (length(occurrence_values) < 2L) {
    stop("Too few finite occurrence values for ", variable)
  }

  background_values <- sample_background(
    file.path(raster_dir, paste0(variable, ".tif")),
    variable,
    cfg$background_sample_size
  )

  mu <- get_anchor_value(model_ref, variable, "mu")
  sigma_lower <- get_anchor_value(model_ref, variable, "sd_lower")
  sigma_upper <- get_anchor_value(model_ref, variable, "sd_upper")
  if (!is.finite(mu) || !is.finite(sigma_lower) || sigma_lower <= 0 ||
      !is.finite(sigma_upper) || sigma_upper <= 0) {
    stop("Invalid asymmetric anchor statistics for ", variable)
  }

  background_limits <- quantile(
    background_values, probs = c(0.001, 0.999),
    na.rm = TRUE, names = FALSE, type = 7
  )
  x_limits <- range(c(background_limits, occurrence_values), na.rm = TRUE)
  x_sequence <- seq(x_limits[1], x_limits[2], length.out = 1000L)

  list(
    variable = variable,
    occurrence = tibble(value = occurrence_values),
    background = tibble(value = background_values),
    curve = tibble(
      value = x_sequence,
      z2 = asymmetric_z2(
        x_sequence, mu, sigma_lower, sigma_upper, cfg$z2_cap
      )
    ),
    mu = mu,
    sigma_lower = sigma_lower,
    sigma_upper = sigma_upper,
    x_limits = x_limits
  )
}

make_density_panel <- function(info, profile, is_labeled, palette) {
  plot <- ggplot() +
    annotate("rect", xmin = info$mu - info$sigma_lower, xmax = info$mu,
             ymin = -Inf, ymax = Inf, fill = palette$lower, alpha = 0.30) +
    annotate("rect", xmin = info$mu, xmax = info$mu + info$sigma_upper,
             ymin = -Inf, ymax = Inf, fill = palette$upper, alpha = 0.30) +
    geom_density(data = info$background, aes(.data$value),
                 fill = palette$bg, colour = NA, alpha = 0.65,
                 linewidth = 0.3, adjust = 1) +
    geom_density(data = info$occurrence, aes(.data$value),
                 fill = palette$occ, colour = palette$text_dark,
                 alpha = 0.72, linewidth = 0.5, adjust = 1) +
    geom_vline(xintercept = info$mu, colour = palette$ref, linewidth = 0.7) +
    geom_vline(
      xintercept = c(info$mu - info$sigma_lower,
                     info$mu + info$sigma_upper),
      colour = c(palette$lower, palette$upper),
      linetype = "dashed", linewidth = 0.5
    ) +
    coord_cartesian(xlim = info$x_limits, expand = FALSE) +
    theme_vera_response(is_labeled, palette)

  if (is_labeled) {
    plot <- plot + labs(
      title = paste0(info$variable, " asymmetric density"),
      subtitle = paste0(cfg$species_label, " — ", profile, " predictors"),
      x = info$variable,
      y = "Density"
    )
  }
  plot
}

make_transform_panel <- function(info, profile, is_labeled, palette) {
  plot <- ggplot() +
    annotate("rect", xmin = info$mu - info$sigma_lower, xmax = info$mu,
             ymin = 0, ymax = cfg$z2_cap,
             fill = palette$lower, alpha = 0.20) +
    annotate("rect", xmin = info$mu, xmax = info$mu + info$sigma_upper,
             ymin = 0, ymax = cfg$z2_cap,
             fill = palette$upper, alpha = 0.20) +
    geom_line(data = info$curve, aes(.data$value, .data$z2),
              colour = palette$curve, linewidth = 1, lineend = "round") +
    geom_hline(yintercept = cfg$z2_cap, colour = palette$cap,
               linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = 4, colour = palette$cap,
               linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = info$mu, colour = palette$ref, linewidth = 0.7) +
    scale_y_continuous(
      breaks = c(0, 4, cfg$z2_cap),
      limits = c(0, cfg$z2_cap * 1.03),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = info$x_limits, clip = "on") +
    theme_vera_response(is_labeled, palette)

  if (is_labeled) {
    plot <- plot + labs(
      title = paste0(info$variable, " asymmetric Z2 transformation"),
      subtitle = paste0(cfg$species_label, " — ", profile, " predictors"),
      x = info$variable,
      y = expression(VRS == Z^2)
    )
  }
  plot
}

annotation_table <- function(variables, model_ref, profile) {
  rows <- lapply(variables, function(variable) {
    tail_scale_center <- get_anchor_value(
      model_ref, variable, "tail_scale_center", numeric = FALSE
    )
    tibble(
      species_code = cfg$species_code,
      species_label = cfg$species_label,
      profile = profile,
      variable = variable,
      mu = get_anchor_value(model_ref, variable, "mu"),
      sigma_lower = get_anchor_value(model_ref, variable, "sd_lower"),
      sigma_upper = get_anchor_value(model_ref, variable, "sd_upper"),
      z2_cap = cfg$z2_cap,
      asymmetry_ratio_raw = get_anchor_value(
        model_ref, variable, "asym_ratio_raw"
      ),
      asymmetry_ratio_ci_low = get_anchor_value(
        model_ref, variable, "asym_ratio_ci_low"
      ),
      asymmetry_ratio_ci_high = get_anchor_value(
        model_ref, variable, "asym_ratio_ci_high"
      ),
      asymmetry_bootstrap_stable = get_anchor_value(
        model_ref, variable, "asymmetry_bootstrap_stable", numeric = FALSE
      ),
      tail_scale_center = tail_scale_center
    )
  })
  bind_rows(rows)
}

render_profile <- function(profile) {
  key <- as.character(profile)
  run_dir <- cfg$profile_dirs[[key]]
  raster_dir <- cfg$raster_dirs[[key]]
  if (is.null(run_dir) || is.null(raster_dir)) {
    stop("Missing run or raster directory for profile ", key)
  }

  csv_dir <- file.path(run_dir, "csvs")
  occurrence_file <- file.path(
    csv_dir, paste0(cfg$species_code, "_04_occurrence_partitions.csv")
  )
  model_file <- file.path(
    csv_dir, paste0(cfg$species_code, "_02_model_reference.csv")
  )
  counts_file <- file.path(
    run_dir, "Images", "core_renderer_plum",
    "vera_pixel_counts_primary_assigned.csv"
  )

  required <- c(occurrence_file, model_file, counts_file)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    warning(
      "Profile ", profile, " skipped; files were not found:\n",
      paste(missing, collapse = "\n"), call. = FALSE
    )
    return(invisible(NULL))
  }

  output_dir <- dir_make(file.path(run_dir, cfg$output_subdir))
  occurrence_df <- read_required_csv(occurrence_file)
  model_ref <- read_required_csv(model_file)
  top_variables <- top_primary_variables(counts_file, cfg$top_n)
  if (!length(top_variables)) stop("No Primary Stressors found for profile ", key)

  write_csv(
    annotation_table(top_variables, model_ref, profile),
    file.path(
      output_dir,
      paste0(cfg$species_code, "_", profile,
             "_top", length(top_variables), "_response_annotations.csv")
    )
  )

  for (variable in top_variables) {
    cat(">>> Profile ", profile, ": rendering ", variable, "\n", sep = "")
    info <- prepare_predictor(variable, occurrence_df, model_ref, raster_dir)

    for (style_name in cfg$styles_to_run) {
      if (!style_name %in% names(styles)) stop("Unknown style: ", style_name)
      palette <- styles[[style_name]]
      style_dir <- dir_make(file.path(output_dir, style_name))

      for (mode_name in cfg$modes_to_run) {
        if (!mode_name %in% c("Labeled", "Unlabeled")) {
          stop("Mode must be Labeled or Unlabeled: ", mode_name)
        }
        is_labeled <- identical(mode_name, "Labeled")
        density_plot <- make_density_panel(
          info, profile, is_labeled, palette
        )
        transform_plot <- make_transform_panel(
          info, profile, is_labeled, palette
        )

        density_file <- file.path(
          style_dir,
          paste0(cfg$species_code, "_", profile, "_", variable,
                 "_Density_", style_name, "_", mode_name, "_400DPI.png")
        )
        transform_file <- file.path(
          style_dir,
          paste0(cfg$species_code, "_", profile, "_", variable,
                 "_Z2_", style_name, "_", mode_name, "_400DPI.png")
        )
        if (!cfg$overwrite &&
            (file.exists(density_file) || file.exists(transform_file))) {
          stop("Output exists and overwrite is FALSE for ", variable)
        }
        ggsave(density_file, density_plot, width = 7.2, height = 5.6,
               units = "in", dpi = cfg$dpi, bg = "white")
        ggsave(transform_file, transform_plot, width = 7.2, height = 5.6,
               units = "in", dpi = cfg$dpi, bg = "white")
      }
    }
  }

  cat(">>> Optional response curves completed: ", output_dir, "\n", sep = "")
  invisible(output_dir)
}

set.seed(cfg$random_seed)
profiles <- unique(as.integer(cfg$profiles_to_run))
if (!length(profiles) || any(!profiles %in% c(19L, 36L))) {
  stop("profiles_to_run must contain 19L, 36L, or both.")
}
for (profile in profiles) render_profile(profile)

cat("\n========================================================================\n")
cat("OPTIONAL RESPONSE-CURVE PANEL RENDERER COMPLETED\n")
cat("Species: ", cfg$species_label, "\n", sep = "")
cat("Profiles requested: ", paste(profiles, collapse = ", "), "\n", sep = "")
cat("========================================================================\n")
