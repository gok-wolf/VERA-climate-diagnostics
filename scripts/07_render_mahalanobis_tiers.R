#!/usr/bin/env Rscript

################################################################################
## FIGURE 4 -- MAHALANOBIS GEOMETRY AND EMPIRICAL TIER CALIBRATION
## Panel A : Core-occupied versus peripheral occurrence distribution
## Panel B : Empirical tier calibration from current core-occupied references
## Focal taxon : Sitta krueperi (Skr)
## Dual-profile: 19- and 36-predictor VERA outputs rendered in the same run.
## Sources     : C:/VERA/Results/{19,36}/Skr_current/csvs
## Outputs     : C:/VERA/Results/{19,36}/Images/figure4_tier_calibration
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

cat(">>> Figure 4 panels A and B are being prepared (Sitta krueperi)...\n")

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
run_subdir           <- "Skr_current"
image_subfolder_name <- "figure4_tier_calibration"

# =============================================================================
# PALETTES (shared across profiles)
# =============================================================================

pal_tiers <- c(
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

theme_figure4 <- theme_classic(base_family = "Arial") +
  theme(
    plot.title      = element_blank(),
    axis.title      = element_text(size = 15, face = "plain", colour = "black"),
    axis.text       = element_text(size = 13, colour = "black"),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 12, colour = "black"),
    legend.position = "bottom",
    panel.border    = element_blank(),
    axis.line       = element_line(colour = "black", linewidth = 0.55),
    panel.grid      = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(10, 12, 10, 12)
  )

# =============================================================================
# DUAL-PROFILE LOOP
# =============================================================================

for (profile_value in profiles_to_render) {

  profile_key   <- as.character(profile_value)
  results_root  <- unname(results_roots[[profile_key]])
  run_dir       <- file.path(results_root, run_subdir)
  csv_dir       <- file.path(run_dir, "csvs")

  occurrence_file <- file.path(csv_dir,
                               paste0(species_code, "_04_occurrence_partitions.csv"))
  model_file      <- file.path(csv_dir,
                               paste0(species_code, "_02_model_reference.csv"))

  figure_dir <- file.path(results_root, "Images", image_subfolder_name)
  if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n========================================================================\n")
  cat(sprintf("FIGURE 4  ::  Sitta krueperi  ::  %s-predictor profile\n", profile_key))
  cat(sprintf("Source CSVs   : %s\n", csv_dir))
  cat(sprintf("Image target  : %s\n", figure_dir))
  cat("========================================================================\n")

  # ---- Source-file check --------------------------------------------------
  missing_files <- c(occurrence_file, model_file)[
    !file.exists(c(occurrence_file, model_file))
  ]
  if (length(missing_files)) {
    warning("Skipping profile ", profile_key,
            ": required files not found:\n", paste(missing_files, collapse = "\n"))
    next
  }

  occurrence_df <- read_csv(occurrence_file, show_col_types = FALSE)
  model_ref     <- read_csv(model_file, show_col_types = FALSE)

  # ---- Save helper ---------------------------------------------------------
  save_panel_pair <- function(p, output_stub, width = 7.4, height = 5.8) {
    withlegend_file <- file.path(
      figure_dir, paste0(output_stub, "_withlegend_400DPI.png")
    )
    nolegend_file <- file.path(
      figure_dir, paste0(output_stub, "_nolegend_400DPI.png")
    )
    ggsave(withlegend_file, p, width = width, height = height,
           units = "in", dpi = 400, bg = "white")
    ggsave(nolegend_file, p + theme(legend.position = "none"),
           width = width, height = height, units = "in", dpi = 400, bg = "white")
    cat("    Saved:", basename(withlegend_file), "and", basename(nolegend_file), "\n")
  }

  # =============================================================================
  # PANEL A -- OCCURRENCE PARTITION DENSITY
  # =============================================================================

  panel_a_df <- occurrence_df %>%
    filter(
      .data$partition %in% c("core_occupied", "peripheral"),
      is.finite(.data$mahal_distance)
    ) %>%
    mutate(
      partition_label = factor(
        .data$partition,
        levels = c("core_occupied", "peripheral"),
        labels = c("Core occupied", "Peripheral")
      )
    )

  if (!nrow(panel_a_df)) {
    warning("Skipping Panel A for profile ", profile_key,
            ": no finite occurrence Mahalanobis distances were available.")
  } else {
    core_threshold <- unique(panel_a_df$core_threshold_distance)
    core_threshold <- core_threshold[is.finite(core_threshold)][1]
    if (!is.finite(core_threshold)) {
      warning("Skipping Panel A for profile ", profile_key,
              ": the core/peripheral threshold was unresolved.")
    } else {
      panel_a_xlim <- c(
        min(panel_a_df$mahal_distance, na.rm = TRUE),
        as.numeric(quantile(
          panel_a_df$mahal_distance, 0.995, na.rm = TRUE,
          names = FALSE, type = 7
        ))
      )

      pA <- ggplot(
        panel_a_df,
        aes(mahal_distance, fill = partition_label, colour = partition_label)
      ) +
        geom_density(alpha = 0.70, linewidth = 0.85, adjust = 1) +
        geom_vline(
          xintercept = core_threshold, colour = "black",
          linetype = "dashed", linewidth = 0.85
        ) +
        scale_fill_manual(values = pal_partition) +
        scale_colour_manual(values = line_partition) +
        coord_cartesian(xlim = panel_a_xlim, expand = FALSE) +
        labs(x = "Occurrence-level Mahalanobis distance", y = "Density") +
        theme_figure4

      save_panel_pair(
        pA, "Figure_4_Panel_A_Mahalanobis_Occurrence_Distribution"
      )
    }
  }

  # =============================================================================
  # PANEL B -- EMPIRICAL TIER CALIBRATION
  # =============================================================================

  tier_reference <- model_ref %>%
    filter(
      .data$section == "current_core_occupied_vrs_values",
      .data$metric  == "current_core_occupied_mahal_distance"
    ) %>%
    transmute(mahal_distance = suppressWarnings(as.numeric(.data$value))) %>%
    filter(is.finite(.data$mahal_distance)) %>%
    arrange(.data$mahal_distance) %>%
    mutate(empirical_probability = row_number() / n())

  if (!nrow(tier_reference)) {
    warning("Skipping Panel B for profile ", profile_key,
            ": no current core-occupied Mahalanobis reference values were found.")
    next
  }

  tier_thresholds <- model_ref %>%
    filter(.data$section == "tier_thresholds") %>%
    transmute(
      subject   = .data$subject,
      threshold = suppressWarnings(as.numeric(.data$value))
    ) %>%
    filter(is.finite(.data$threshold))

  required_subjects <- c("moderate_departure", "restriction_zone",
                         "high_extrapolative_stress")
  if (!all(required_subjects %in% tier_thresholds$subject)) {
    warning("Skipping Panel B for profile ", profile_key,
            ": one or more empirical Mahalanobis tier thresholds are missing.")
    next
  }

  threshold_value <- function(subject_name) {
    tier_thresholds$threshold[tier_thresholds$subject == subject_name][1]
  }

  break_1 <- threshold_value("moderate_departure")
  break_2 <- threshold_value("restriction_zone")
  break_3 <- threshold_value("high_extrapolative_stress")

  break_df <- tibble(
    boundary    = factor(c("80%", "95%", "99%"),
                         levels = c("80%", "95%", "99%")),
    probability = c(0.80, 0.95, 0.99),
    threshold   = c(break_1, break_2, break_3)
  )

  panel_b_max <- max(
    as.numeric(quantile(
      tier_reference$mahal_distance, 0.995, na.rm = TRUE,
      names = FALSE, type = 7
    )),
    break_3 * 1.08
  )

  tier_rectangles <- tibble(
    tier = factor(names(pal_tiers), levels = names(pal_tiers)),
    xmin = c(0, break_1, break_2, break_3),
    xmax = c(break_1, break_2, break_3, panel_b_max)
  )

  pB <- ggplot() +
    geom_rect(
      data = tier_rectangles,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = tier),
      # Full opacity ensures exact colour identity with spatial Panel D and with
      # the single shared tier legend used in the assembled Figure 4.
      alpha = 1, colour = NA
    ) +
    geom_step(
      data = tier_reference,
      aes(mahal_distance, empirical_probability),
      colour = "black", linewidth = 1, direction = "hv"
    ) +
    geom_segment(
      data = break_df,
      aes(x = threshold, xend = threshold, y = 0, yend = probability),
      colour = "black", linetype = "dashed", linewidth = 0.75
    ) +
    geom_segment(
      data = break_df,
      aes(x = 0, xend = threshold, y = probability, yend = probability),
      colour = "grey35", linetype = "dotted", linewidth = 0.65
    ) +
    geom_point(
      data = break_df,
      aes(threshold, probability),
      shape = 21, size = 3, stroke = 0.7,
      fill = "white", colour = "black"
    ) +
    scale_fill_manual(values = pal_tiers, name = NULL, drop = FALSE) +
    scale_x_continuous(
      limits = c(0, panel_b_max), expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 0.80, 0.95, 0.99, 1),
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Core-reference Mahalanobis distance",
      y = "Empirical cumulative probability"
    ) +
    theme_figure4

  save_panel_pair(pB, "Figure_4_Panel_B_Empirical_Tier_Calibration")

  write.csv(
    break_df,
    file.path(figure_dir, "Figure_4_Panel_B_Resolved_Tier_Thresholds.csv"),
    row.names = FALSE
  )

  cat("\n")
  cat(sprintf(">>> FIGURE 4 PANELS A AND B COMPLETED (%s-predictor profile)\n",
              profile_key))
  cat(sprintf(">>> Occurrence core threshold : %.6f\n", core_threshold))
  cat(sprintf(">>> Empirical tier thresholds : %.6f | %.6f | %.6f\n",
              break_1, break_2, break_3))
  cat(sprintf(">>> Output directory          : %s\n", figure_dir))
}

cat("\n========================================================================\n")
cat(">>> ALL FIGURE 4 PROFILES COMPLETED\n")
cat(">>> Focal taxon : ", species_label, " (", species_code, ")\n", sep = "")
cat(">>> Profiles    : 19 + 36\n")
cat(">>> Panels B and D may correctly share one tier legend within each profile.\n")
cat(">>> Image roots :\n")
for (pk in names(results_roots)) {
  cat("      ", file.path(results_roots[[pk]], "Images", image_subfolder_name), "\n",
      sep = "")
}
cat("========================================================================\n")