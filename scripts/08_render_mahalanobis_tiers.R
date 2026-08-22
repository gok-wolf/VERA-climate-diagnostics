#!/usr/bin/env Rscript

# =============================================================================
# OPTIONAL VERA MAHALANOBIS TIER-CALIBRATION FIGURES
# Sitta krueperi; supports the 19- and 36-predictor profiles in one script.
#
# This renderer reads completed canonical VERA CSV outputs. It does not
# recalculate or modify Mahalanobis distances, occurrence partitions, tier
# thresholds, VRS surfaces, or any other canonical product.
# Each profile output directory also receives README_08.txt with a concise
# interpretation and file guide.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
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
  output_subdir = file.path("Images", "optional_mahalanobis_tiers"),
  dpi = 400,
  overwrite = TRUE
)

font_main <- "Arial"

pal_tiers <- c(
  "Core Climate" = "#FFF2FA",
  "Moderate Departure" = "#F3C4EC",
  "Restriction Zone" = "#D695CC",
  "High Extrapolative Stress" = "#CC7CBF"
)

pal_partition <- c(
  "Core occupied" = "#F3C4EC",
  "Peripheral" = "#CC7CBF"
)

line_partition <- c(
  "Core occupied" = "#B76AAA",
  "Peripheral" = "#824177"
)

theme_optional <- theme_classic(base_family = font_main) +
  theme(
    plot.title = element_text(
      family = font_main, face = "plain", size = 17,
      colour = "black", hjust = 0.5, margin = margin(b = 12)
    ),
    plot.subtitle = element_text(
      family = font_main, face = "plain", size = 12,
      colour = "grey30", hjust = 0.5, margin = margin(b = 10)
    ),
    axis.title = element_text(
      family = font_main, face = "plain", size = 15, colour = "black"
    ),
    axis.text = element_text(
      family = font_main, face = "plain", size = 13, colour = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      family = font_main, face = "plain", size = 12, colour = "black"
    ),
    legend.position = "bottom",
    axis.line = element_line(colour = "black", linewidth = 0.55),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(12, 14, 12, 14)
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

save_plot_pair <- function(plot, out_dir, stub, width = 7.4, height = 5.8) {
  with_legend <- file.path(out_dir, paste0(stub, "_with_legend_400DPI.png"))
  no_legend <- file.path(out_dir, paste0(stub, "_without_legend_400DPI.png"))

  if (!cfg$overwrite && (file.exists(with_legend) || file.exists(no_legend))) {
    stop("Output already exists and overwrite is FALSE: ", stub)
  }

  ggsave(with_legend, plot, width = width, height = height,
         units = "in", dpi = cfg$dpi, bg = "white")
  ggsave(no_legend, plot + theme(legend.position = "none"),
         width = width, height = height, units = "in",
         dpi = cfg$dpi, bg = "white")
  c(with_legend, no_legend)
}

required_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

render_profile <- function(profile) {
  key <- as.character(profile)
  run_dir <- cfg$profile_dirs[[key]]
  if (is.null(run_dir) || !nzchar(run_dir)) stop("No directory for profile ", key)

  csv_dir <- file.path(run_dir, "csvs")
  occurrence_file <- file.path(
    csv_dir, paste0(cfg$species_code, "_04_occurrence_partitions.csv")
  )
  model_file <- file.path(
    csv_dir, paste0(cfg$species_code, "_02_model_reference.csv")
  )
  required <- c(occurrence_file, model_file)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    warning(
      "Profile ", key, " skipped; required files were not found:\n",
      paste(missing, collapse = "\n"), call. = FALSE
    )
    return(invisible(NULL))
  }

  out_dir <- dir_make(file.path(run_dir, cfg$output_subdir))
  occurrence_df <- read_csv(occurrence_file, show_col_types = FALSE)
  model_ref <- read_csv(model_file, show_col_types = FALSE)

  required_columns(
    occurrence_df,
    c("partition", "mahal_distance", "core_threshold_distance"),
    basename(occurrence_file)
  )
  required_columns(
    model_ref,
    c("section", "subject", "metric", "value"),
    basename(model_file)
  )

  # ---------------------------------------------------------------------------
  # Panel A: occurrence-level climatic Mahalanobis distribution
  # ---------------------------------------------------------------------------

  panel_a_df <- occurrence_df %>%
    filter(
      .data$partition %in% c("core_occupied", "peripheral"),
      is.finite(suppressWarnings(as.numeric(.data$mahal_distance)))
    ) %>%
    mutate(
      mahal_distance = suppressWarnings(as.numeric(.data$mahal_distance)),
      partition_label = factor(
        .data$partition,
        levels = c("core_occupied", "peripheral"),
        labels = c("Core occupied", "Peripheral")
      )
    )

  if (!nrow(panel_a_df)) {
    stop("No finite occurrence Mahalanobis distances for profile ", key)
  }

  core_threshold <- suppressWarnings(
    as.numeric(unique(panel_a_df$core_threshold_distance))
  )
  core_threshold <- core_threshold[is.finite(core_threshold)][1]
  if (!is.finite(core_threshold)) {
    stop("Core/peripheral threshold unresolved for profile ", key)
  }

  panel_a_xlim <- c(
    min(panel_a_df$mahal_distance, na.rm = TRUE),
    as.numeric(quantile(panel_a_df$mahal_distance, 0.995, na.rm = TRUE,
                        names = FALSE, type = 7))
  )

  panel_a <- ggplot(
    panel_a_df,
    aes(.data$mahal_distance, fill = .data$partition_label,
        colour = .data$partition_label)
  ) +
    geom_density(alpha = 0.70, linewidth = 0.85, adjust = 1) +
    geom_vline(xintercept = core_threshold, colour = "black",
               linetype = "dashed", linewidth = 0.85) +
    scale_fill_manual(values = pal_partition) +
    scale_colour_manual(values = line_partition) +
    coord_cartesian(xlim = panel_a_xlim, expand = FALSE) +
    labs(
      title = paste0(cfg$species_label, " — occurrence partition"),
      subtitle = paste0(profile, "-predictor profile"),
      x = "Occurrence-level Mahalanobis distance",
      y = "Density"
    ) +
    theme_optional

  occurrence_plot_files <- save_plot_pair(
    panel_a, out_dir,
    paste0(cfg$species_code, "_", profile,
           "_Mahalanobis_Occurrence_Distribution")
  )

  # ---------------------------------------------------------------------------
  # Panel B: empirical tier calibration
  # ---------------------------------------------------------------------------

  tier_reference <- model_ref %>%
    filter(
      .data$section == "current_core_occupied_vrs_values",
      .data$metric == "current_core_occupied_mahal_distance"
    ) %>%
    transmute(mahal_distance = suppressWarnings(as.numeric(.data$value))) %>%
    filter(is.finite(.data$mahal_distance)) %>%
    arrange(.data$mahal_distance) %>%
    mutate(empirical_probability = row_number() / n())

  if (!nrow(tier_reference)) {
    stop("No core-occupied Mahalanobis reference values for profile ", key)
  }

  tier_thresholds <- model_ref %>%
    filter(.data$section == "tier_thresholds") %>%
    transmute(
      subject = as.character(.data$subject),
      threshold = suppressWarnings(as.numeric(.data$value))
    ) %>%
    filter(is.finite(.data$threshold))

  required_subjects <- c(
    "moderate_departure", "restriction_zone", "high_extrapolative_stress"
  )
  if (!all(required_subjects %in% tier_thresholds$subject)) {
    stop("One or more tier thresholds are missing for profile ", key)
  }

  threshold_value <- function(subject_name) {
    tier_thresholds$threshold[tier_thresholds$subject == subject_name][1]
  }
  break_1 <- threshold_value("moderate_departure")
  break_2 <- threshold_value("restriction_zone")
  break_3 <- threshold_value("high_extrapolative_stress")

  break_df <- tibble(
    boundary = factor(c("80%", "95%", "99%"),
                      levels = c("80%", "95%", "99%")),
    probability = c(0.80, 0.95, 0.99),
    threshold = c(break_1, break_2, break_3)
  )

  panel_b_max <- max(
    as.numeric(quantile(tier_reference$mahal_distance, 0.995, na.rm = TRUE,
                        names = FALSE, type = 7)),
    break_3 * 1.08
  )
  tier_rectangles <- tibble(
    tier = factor(names(pal_tiers), levels = names(pal_tiers)),
    xmin = c(0, break_1, break_2, break_3),
    xmax = c(break_1, break_2, break_3, panel_b_max)
  )

  panel_b <- ggplot() +
    geom_rect(
      data = tier_rectangles,
      aes(xmin = .data$xmin, xmax = .data$xmax,
          ymin = 0, ymax = 1, fill = .data$tier),
      alpha = 1, colour = NA
    ) +
    geom_step(
      data = tier_reference,
      aes(.data$mahal_distance, .data$empirical_probability),
      colour = "black", linewidth = 1, direction = "hv"
    ) +
    geom_segment(
      data = break_df,
      aes(x = .data$threshold, xend = .data$threshold,
          y = 0, yend = .data$probability),
      colour = "black", linetype = "dashed", linewidth = 0.75
    ) +
    geom_segment(
      data = break_df,
      aes(x = 0, xend = .data$threshold,
          y = .data$probability, yend = .data$probability),
      colour = "grey35", linetype = "dotted", linewidth = 0.65
    ) +
    geom_point(
      data = break_df,
      aes(.data$threshold, .data$probability),
      shape = 21, size = 3, stroke = 0.7,
      fill = "white", colour = "black"
    ) +
    scale_fill_manual(values = pal_tiers, drop = FALSE) +
    scale_x_continuous(limits = c(0, panel_b_max),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 0.80, 0.95, 0.99, 1),
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = paste0(cfg$species_label, " — empirical tier calibration"),
      subtitle = paste0(profile, "-predictor profile"),
      x = "Core-reference Mahalanobis distance",
      y = "Empirical cumulative probability"
    ) +
    theme_optional

  tier_plot_files <- save_plot_pair(
    panel_b, out_dir,
    paste0(cfg$species_code, "_", profile,
           "_Empirical_Tier_Calibration")
  )

  threshold_file <- file.path(
    out_dir,
    paste0(cfg$species_code, "_", profile,
           "_Resolved_Tier_Thresholds.csv")
  )
  write_csv(
    mutate(
      break_df,
      species_code = cfg$species_code,
      species_label = cfg$species_label,
      profile = profile,
      .before = 1
    ),
    threshold_file
  )

  readme_file <- file.path(out_dir, "README_08.txt")
  readme_text <- c(
    "VERA OPTIONAL MAHALANOBIS TIER-CALIBRATION FIGURES",
    "====================================================",
    "",
    paste0("Species: ", cfg$species_label, " (", cfg$species_code, ")"),
    paste0("Predictor profile: ", profile),
    "Renderer: 08_render_mahalanobis_tiers.R",
    "Status: OPTIONAL publication and Supplementary Information product",
    "",
    "PURPOSE",
    "-------",
    "This directory contains publication-oriented visualisations of existing",
    "VERA Mahalanobis occurrence partitions and empirical climatic-tier",
    "calibration. Script 08 does not recalculate or modify Mahalanobis",
    "distances, occurrence partitions, tier thresholds, VRS surfaces, or any",
    "other canonical VERA product.",
    "",
    "INPUTS",
    "------",
    paste0("- ", basename(occurrence_file)),
    paste0("- ", basename(model_file)),
    "",
    "OUTPUT GUIDE",
    "------------",
    "Occurrence Distribution figures:",
    "- Display occurrence-level climatic Mahalanobis-distance densities for",
    "  the existing Core occupied and Peripheral partitions.",
    "- The dashed vertical line is the existing core/peripheral threshold.",
    "- The panel visualises the rule-based partition; it is not independent",
    "  validation of that partition.",
    "",
    "Empirical Tier Calibration figures:",
    "- Display the empirical cumulative distribution of Mahalanobis distance",
    "  within the existing core-occupied climatic reference.",
    "- The 80%, 95%, and 99% boundaries are the existing empirical tier",
    "  thresholds exported by the completed VERA run.",
    "",
    "Resolved Tier Thresholds CSV:",
    "- Records the three plotted empirical boundaries and their resolved",
    "  Mahalanobis-distance values.",
    "",
    "FILES",
    "-----",
    paste0("- ", basename(occurrence_plot_files[1])),
    paste0("- ", basename(occurrence_plot_files[2])),
    paste0("- ", basename(tier_plot_files[1])),
    paste0("- ", basename(tier_plot_files[2])),
    paste0("- ", basename(threshold_file)),
    "",
    "READING RULES",
    "-------------",
    "1. Mahalanobis distance is a covariance-aware, directionally symmetric",
    "   measure of climatic departure from the occupied reference.",
    "2. Larger distance indicates greater multivariate climatic departure. It",
    "   does not mean poorer habitat, lower occurrence probability, or lower",
    "   physiological performance.",
    "3. The four tiers are empirical diagnostic classes, not physiological",
    "   limits, suitability classes, or causal ecological thresholds.",
    "4. The 19- and 36-predictor profiles have different dimensional geometry.",
    "   Their raw Mahalanobis-distance magnitudes should not be treated as",
    "   directly interchangeable.",
    "5. Files with and without legends contain the same plotted data. The",
    "   legend-free versions are supplied only for flexible panel assembly.",
    "6. These figures are optional. Their inclusion should follow the research",
    "   question, Supplementary Information plan, and journal requirements.",
    "",
    "REPORTING NOTE",
    "--------------",
    "Report the occurrence preparation, covariance treatment, core-partition",
    "rule, empirical tier probabilities, predictor profile, and software",
    "version in the manuscript or Supplementary Information.",
    ""
  )
  write_lines_utf8(readme_text, readme_file)

  cat(
    ">>> Optional Mahalanobis tier figures completed for profile ",
    profile, ": ", out_dir, "\n", sep = ""
  )
  invisible(out_dir)
}

profiles <- unique(as.integer(cfg$profiles_to_run))
if (!length(profiles) || any(!profiles %in% c(19L, 36L))) {
  stop("profiles_to_run must contain 19L, 36L, or both.")
}

for (profile in profiles) render_profile(profile)

cat("\n========================================================================\n")
cat("OPTIONAL MAHALANOBIS TIER-CALIBRATION RENDERER COMPLETED\n")
cat("Species: ", cfg$species_label, "\n", sep = "")
cat("Profiles requested: ", paste(profiles, collapse = ", "), "\n", sep = "")
cat("========================================================================\n")