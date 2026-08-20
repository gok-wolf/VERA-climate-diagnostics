#!/usr/bin/env Rscript

################################################################################
## VERA ADD-ON 01-02 -- MASTER MAP RENDERER AND CSV PLOTTER
## Focal taxon : Sitta krueperi (Skr)
## Dual-profile: 19- and 36-predictor VERA outputs rendered in the same run.
## Sources     : C:/VERA/Results/{19,36}/Skr_current/addons
## Outputs     : C:/VERA/Results/{19,36}/Images/addon_renderer_plum
## No shapefile is used; map bounds are derived from the raster extent.
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(tidyterra)
  library(scales)
  library(grid)
})

cat(">>> [1/5] Initializing canonical VERA Add-On renderer (Sitta krueperi)...\n")

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
image_subfolder_name  <- "addon_renderer_plum"

# =============================================================================
# PALETTES (shared across profiles)
# =============================================================================

pal_addon         <- c("#FFF2FA", "#F3C4EC", "#D695CC", "#CC7CBF")
pal_addon_reverse <- rev(pal_addon)

# Signed divergence must use distinct hues on either side of zero.
pal_diverging <- c("#3F78A8", "#AFCBE0", "#FFF7FB", "#E5B5DC", "#A84D98")

cols_tail <- c(
  "1" = "#FF6467", "2" = "#FFDC5C",
  "3" = "#51A2FF", "4" = "#7CDE8D"
)
lbls_tail <- c(
  "1" = "Upper-tail dominant\n(high-edge restriction)",
  "2" = "Lower-tail dominant\n(low-edge restriction)",
  "3" = "Balanced / near-symmetric",
  "4" = "No directional signal"
)

# Agreement classes are nominal/relational, not an ordinal four-step scale.
cols_agreement <- c(
  "1" = "#A84D98", "2" = "#3F78A8",
  "3" = "#8073AC", "4" = "#D9D9D9"
)
lbls_agreement <- c(
  "1" = "VRS-dominant\nhigher restriction rank",
  "2" = "Mahalanobis-dominant\nhigher distance rank",
  "3" = "Concordant\nhigher-rank signals",
  "4" = "Concordant\nlower-rank signals"
)

# =============================================================================
# SHARED THEMES AND HELPERS
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

theme_csv <- theme_classic(base_family = "Arial") +
  theme(
    plot.title       = element_text(size = 16, face = "plain", colour = "black"),
    plot.subtitle    = element_text(size = 12, colour = "grey30"),
    axis.title       = element_text(size = 13, face = "plain", colour = "black"),
    axis.text        = element_text(size = 11, colour = "black"),
    legend.title     = element_text(size = 12, face = "plain"),
    legend.text      = element_text(size = 11),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

clean_class_label <- function(x) {
  trimws(gsub("_", " ", gsub("^[0-9]+[ _-]*", "", as.character(x))))
}

wrap_label <- function(x, width = 22) {
  vapply(as.character(x), function(z) paste(strwrap(z, width), collapse = "\n"),
         character(1))
}

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
# CATEGORICAL, SEQUENTIAL, AND DIVERGING MAP BUILDERS
# =============================================================================

make_categorical_map <- function(idx_file, labels, colours, legend_title, renv) {
  r <- rast(idx_file)
  vals <- sort(unique(suppressWarnings(as.integer(values(r, mat = FALSE)))))
  vals <- vals[is.finite(vals) & as.character(vals) %in% names(labels)]
  if (!length(vals)) stop("No valid categorical codes were found in ", idx_file)

  r <- as.factor(r)
  levels(r) <- data.frame(value = vals, label = unname(labels[as.character(vals)]))
  names(r) <- "class"
  scale_values <- setNames(
    unname(colours[as.character(vals)]), unname(labels[as.character(vals)])
  )

  add_map_geometry(
    ggplot() +
      geom_spatraster(data = r, aes(fill = class), na.rm = TRUE) +
      scale_fill_manual(
        values = scale_values, name = legend_title, drop = FALSE,
        na.translate = FALSE, na.value = "transparent"
      ),
    renv
  )
}

make_sequential_map <- function(tif_file, legend_title, limits = NULL,
                                palette = pal_addon, renv) {
  r <- rast(tif_file)
  names(r) <- "value"
  add_map_geometry(
    ggplot() +
      geom_spatraster(data = r, aes(fill = value), na.rm = TRUE) +
      scale_fill_gradientn(
        colours = palette, values = scales::rescale(seq_along(palette)),
        limits = limits, oob = scales::squish, name = legend_title,
        na.value = "transparent",
        guide = guide_colourbar(
          barheight = unit(7.5, "cm"), barwidth = unit(0.75, "cm")
        )
      ),
    renv
  )
}

make_divergence_map <- function(tif_file, renv) {
  r <- rast(tif_file)
  names(r) <- "value"
  max_abs <- as.numeric(global(abs(r), "max", na.rm = TRUE)[1, 1])
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1

  add_map_geometry(
    ggplot() +
      geom_spatraster(data = r, aes(fill = value), na.rm = TRUE) +
      scale_fill_gradientn(
        colours = pal_diverging,
        values  = c(0, 0.25, 0.50, 0.75, 1),
        limits  = c(-max_abs, max_abs), oob = scales::squish,
        name    = "VRS - Mahalanobis\nrank difference",
        na.value = "transparent",
        guide = guide_colourbar(
          barheight = unit(7.5, "cm"), barwidth = unit(0.75, "cm")
        )
      ),
    renv
  )
}

# =============================================================================
# DUAL-PROFILE LOOP
# =============================================================================

for (profile_value in profiles_to_render) {

  profile_key <- as.character(profile_value)

  results_root  <- unname(results_roots[[profile_key]])
  run_dir       <- file.path(results_root, run_subdir)
  path_ras      <- file.path(run_dir, "addons", "rasters")
  path_csv      <- file.path(run_dir, "addons", "csvs")
  inventory_file <- file.path(run_dir, paste0(species_code, "_output_inventory.csv"))

  path_img <- file.path(results_root, "Images", image_subfolder_name)
  if (!dir.exists(path_img)) dir.create(path_img, recursive = TRUE, showWarnings = FALSE)

  cat("\n========================================================================\n")
  cat(sprintf("ADD-ON RENDERER  ::  Sitta krueperi  ::  %s-predictor profile\n",
              profile_key))
  cat(sprintf("Source addons  : %s\n", path_ras))
  cat(sprintf("Image target   : %s\n", path_img))
  cat("========================================================================\n")

  # ---- Build the render environment from a template raster -----------------
  stub <- function(name) paste0(species_code, "_current_", name, ".tif")
  template_tif <- file.path(path_ras, stub("vrs_rank01"))
  if (!file.exists(template_tif)) {
    warning("Skipping profile ", profile_key,
            ": template add-on raster missing (", template_tif, ")")
    next
  }
  template_raster <- rast(template_tif)
  renv <- make_render_env(template_raster)

  # ---- Filename builders and save helpers -----------------------------------
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

  save_csv_plot <- function(p, filename, width, height) {
    ggsave(
      file.path(path_img, filename), plot = p, width = width, height = height,
      units = "in", dpi = 400, bg = "white"
    )
  }

  find_addon_csv <- function(filename) {
    path <- file.path(path_csv, filename)
    if (file.exists(path)) path else NA_character_
  }

  # =============================================================================
  # 1. ADD-ON RASTER MAPS
  # =============================================================================

  cat(sprintf(">>> [2/5] Rendering Add-On raster maps (%s vars)...\n", profile_key))

  map_manifest <- list(
    tail_direction_dominant = list(
      file = file.path(path_ras, stub("tail_direction_dominant_idx")),
      type = "tail", title = "Dominant Tail-Direction Diagnostic"
    ),
    tail_direction_net = list(
      file = file.path(path_ras, stub("tail_direction_net_idx")),
      type = "tail", title = "Net Tail-Direction Diagnostic"
    ),
    mahal_vrs_agreement_4class = list(
      file = file.path(path_ras, stub("mahal_vrs_agreement_4class_idx")),
      type = "agreement", title = "Mahalanobis-VRS Agreement Classes"
    ),
    vrs_rank01 = list(
      file = file.path(path_ras, stub("vrs_rank01")), type = "rank",
      legend = "Mean VRS\npercentile rank", title = "Mean VRS Percentile Rank"
    ),
    mahal_rank01 = list(
      file = file.path(path_ras, stub("mahal_rank01")), type = "rank",
      legend = "Mahalanobis\npercentile rank", title = "Mahalanobis Percentile Rank"
    ),
    mahal_vrs_divergence = list(
      file = file.path(path_ras, stub("mahal_vrs_divergence")),
      type = "divergence",
      title = "Percentile-Rank Divergence: VRS minus Mahalanobis"
    ),
    mahal_vrs_abs_disagreement = list(
      file = file.path(path_ras, stub("mahal_vrs_abs_disagreement")),
      type = "absolute", legend = "Absolute rank\ndisagreement",
      title = "Absolute Rank Disagreement"
    )
  )

  for (output_stub in names(map_manifest)) {
    s <- map_manifest[[output_stub]]
    if (!file.exists(s$file)) {
      warning("Skipping missing Add-On raster: ", s$file)
      next
    }
    cat("    +", output_stub, "\n")
    p <- switch(
      s$type,
      tail       = make_categorical_map(s$file, lbls_tail, cols_tail,
                                        "Tail direction", renv),
      agreement  = make_categorical_map(s$file, lbls_agreement, cols_agreement,
                                        "Agreement class", renv),
      rank       = make_sequential_map(s$file, s$legend, limits = c(0, 1),
                                       renv = renv),
      absolute   = make_sequential_map(s$file, s$legend, renv = renv),
      divergence = make_divergence_map(s$file, renv),
      stop("Unknown Add-On map type: ", s$type)
    )
    save_map_pair(p, paste0("vera_addon_", output_stub), s$title)
    rm(p)
    gc(verbose = FALSE)
  }

  # =============================================================================
  # 2. CSV FIGURES (Add-On 01)
  # =============================================================================

  cat(sprintf(">>> [3/5] Rendering Add-On statistical figures (%s)...\n", profile_key))

  pixel_counts_file <- find_addon_csv("addon01_tail_direction_pixel_counts.csv")
  if (!is.na(pixel_counts_file)) {
    df <- read.csv(pixel_counts_file, stringsAsFactors = FALSE) %>%
      mutate(
        class_id = as.character(.data$class_id),
        class_label = factor(
          .data$class_id, levels = names(lbls_tail),
          labels = wrap_label(lbls_tail, 21)
        ),
        definition_label = recode(
          .data$definition,
          dominant_max_tail_signal = "Dominant maximum-tail signal",
          net_summed_tail_signal   = "Net summed-tail signal"
        )
      )
    p <- ggplot(df, aes(class_label, percent_valid, fill = class_id)) +
      geom_col(colour = "black", linewidth = 0.35, width = 0.68) +
      geom_text(aes(label = sprintf("%.1f%%", percent_valid)),
                vjust = -0.55, size = 4) +
      facet_wrap(~definition_label, nrow = 1) +
      scale_fill_manual(values = cols_tail, guide = "none") +
      scale_y_continuous(
        limits = c(0, max(df$percent_valid, na.rm = TRUE) + 8),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(title = "Tail-Direction Diagnostic Classes", x = NULL,
           y = "Valid landscape area (%)") +
      theme_csv + theme(axis.text.x = element_text(lineheight = 0.9))
    save_csv_plot(p, "plot_addon01_tail_direction_classes.png", 12, 6)
  }

  agreement_file <- find_addon_csv("addon01_tail_direction_dominant_vs_net_agreement.csv")
  if (!is.na(agreement_file)) {
    df <- read.csv(agreement_file, stringsAsFactors = FALSE) %>%
      mutate(
        dominant_label = wrap_label(clean_class_label(.data$dominant_class), 18),
        net_label      = wrap_label(clean_class_label(.data$net_class), 18)
      )
    df$dominant_label <- factor(df$dominant_label, levels = unique(df$dominant_label))
    df$net_label      <- factor(df$net_label, levels = rev(unique(df$net_label)))
    agreement_pct <- unique(df$overall_agreement_pct)
    agreement_pct <- agreement_pct[is.finite(agreement_pct)][1]
    p <- ggplot(df, aes(dominant_label, net_label, fill = percent_valid)) +
      geom_tile(colour = "white", linewidth = 0.8) +
      geom_text(aes(label = sprintf("%.1f%%", percent_valid)), size = 4) +
      scale_fill_gradientn(colours = pal_addon, name = "Landscape\narea (%)") +
      labs(
        title    = "Dominant versus Net Tail-Direction Agreement",
        subtitle = sprintf("Overall agreement: %.1f%%", agreement_pct),
        x = "Dominant maximum-tail classification",
        y = "Net summed-tail classification"
      ) +
      theme_csv +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 25, hjust = 1))
    save_csv_plot(p, "plot_addon01_tail_direction_agreement_heatmap.png", 10, 8)
  }

  latitudinal_file <- find_addon_csv("addon01_tail_direction_latitudinal_profile.csv")
  if (!is.na(latitudinal_file)) {
    df <- read.csv(latitudinal_file, stringsAsFactors = FALSE) %>%
      mutate(
        latitude_midpoint = (
          suppressWarnings(as.numeric(sub("^\\(([-0-9.]+),.*$", "\\1", .data$lat_band))) +
            suppressWarnings(as.numeric(sub("^.*,([-0-9.]+)\\]$", "\\1", .data$lat_band)))
        ) / 2,
        class_id = sub("^([0-9]+).*$", "\\1", .data$tail_class),
        class_label = factor(
          .data$class_id, levels = names(lbls_tail),
          labels = unname(lbls_tail)
        )
      ) %>%
      filter(is.finite(.data$latitude_midpoint)) %>%
      group_by(.data$latitude_midpoint) %>%
      mutate(percent_band = 100 * .data$n_pixels / sum(.data$n_pixels)) %>%
      ungroup()

    p <- ggplot(df, aes(latitude_midpoint, percent_band, colour = class_label)) +
      geom_line(linewidth = 1.15) +
      geom_point(size = 2) +
      scale_colour_manual(
        values = setNames(unname(cols_tail), unname(lbls_tail)),
        name = "Tail direction", drop = FALSE
      ) +
      scale_y_continuous(limits = c(0, 100), labels = label_percent(scale = 1)) +
      labs(
        title = "Latitudinal Tail-Direction Profile",
        x = "Latitude-band midpoint (degrees north)",
        y = "Share within latitude band"
      ) +
      theme_csv + theme(legend.position = "bottom")
    save_csv_plot(p, "plot_addon01_tail_direction_latitudinal_profile.png", 11, 6.5)
  }

  # =============================================================================
  # 3. CSV FIGURES (Add-On 02)
  # =============================================================================

  cat(sprintf(">>> [4/5] Rendering Add-On 02 divergence figures (%s)...\n", profile_key))

  divergence_counts_file <- find_addon_csv("addon02_divergence_class_counts.csv")
  if (!is.na(divergence_counts_file)) {
    df <- read.csv(divergence_counts_file, stringsAsFactors = FALSE) %>%
      mutate(
        class_id = as.character(.data$class_id),
        class_label = factor(
          .data$class_id, levels = names(lbls_agreement),
          labels = wrap_label(lbls_agreement, 20)
        )
      )
    p <- ggplot(df, aes(class_label, percent_valid, fill = class_id)) +
      geom_col(colour = "black", linewidth = 0.35, width = 0.68) +
      geom_text(aes(label = sprintf("%.1f%%", percent_valid)),
                vjust = -0.65, size = 4.2) +
      scale_fill_manual(values = cols_agreement, guide = "none") +
      scale_y_continuous(
        limits = c(0, max(df$percent_valid, na.rm = TRUE) + 7),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(title = "Mahalanobis-VRS Agreement Classes", x = NULL,
           y = "Valid landscape area (%)") +
      theme_csv + theme(axis.text.x = element_text(lineheight = 0.9))
    save_csv_plot(p, "plot_addon02_mahal_vrs_agreement_classes.png", 10, 6)
  }

  divergence_summary_file <- find_addon_csv("addon02_divergence_summary.csv")
  if (!is.na(divergence_summary_file)) {
    divergence_summary <- read.csv(divergence_summary_file, stringsAsFactors = FALSE)

    summary_value <- function(metric) {
      suppressWarnings(as.numeric(
        divergence_summary$value[divergence_summary$metric == metric][1]
      ))
    }

    divergence_plot_df <- tibble(
      metric = factor(
        c("VRS-dominant", "Mahalanobis-dominant", "Within threshold"),
        levels = c("VRS-dominant", "Mahalanobis-dominant", "Within threshold")
      ),
      percent = c(
        summary_value("percent_vrs_dominant"),
        summary_value("percent_mahalanobis_dominant"),
        summary_value("percent_within_threshold")
      ),
      fill_id = factor(c("1", "2", "3"), levels = c("1", "2", "3"))
    ) %>% filter(is.finite(.data$percent))

    if (nrow(divergence_plot_df)) {
      p <- ggplot(divergence_plot_df, aes(metric, percent, fill = fill_id)) +
        geom_col(colour = "black", linewidth = 0.35, width = 0.65) +
        geom_text(aes(label = sprintf("%.1f%%", percent)),
                  vjust = -0.65, size = 4.2) +
        scale_fill_manual(
          values = c("1" = "#A84D98", "2" = "#3F78A8", "3" = "#8073AC"),
          guide = "none"
        ) +
        scale_y_continuous(
          limits = c(0, max(divergence_plot_df$percent) + 8),
          expand = expansion(mult = c(0, 0))
        ) +
        labs(title = "Mahalanobis-VRS Divergence Summary",
             x = NULL, y = "Valid landscape area (%)") +
        theme_csv
      save_csv_plot(p, "plot_addon02_mahal_vrs_divergence_summary.png", 8.5, 5.8)
    }

    disagreement_plot_df <- tibble(
      metric = factor(
        c("Mean", "Median", "75th percentile"),
        levels = c("Mean", "Median", "75th percentile")
      ),
      value = c(
        summary_value("mean_abs_disagreement"),
        summary_value("median_abs_disagreement"),
        summary_value("q75_abs_disagreement")
      ),
      fill_id = factor(c("1", "2", "3"), levels = c("1", "2", "3"))
    ) %>% filter(is.finite(.data$value))

    if (nrow(disagreement_plot_df)) {
      p <- ggplot(disagreement_plot_df, aes(metric, value, fill = fill_id)) +
        geom_col(colour = "black", linewidth = 0.35, width = 0.62) +
        geom_text(aes(label = sprintf("%.3f", value)),
                  vjust = -0.65, size = 4.2) +
        scale_fill_manual(
          values = c("1" = "#F3C4EC", "2" = "#D695CC", "3" = "#CC7CBF"),
          guide = "none"
        ) +
        scale_y_continuous(
          limits = c(0, max(disagreement_plot_df$value) + 0.04),
          expand = expansion(mult = c(0, 0))
        ) +
        labs(title = "Absolute Rank Disagreement",
             x = NULL, y = "Absolute disagreement magnitude") +
        theme_csv
      save_csv_plot(p, "plot_addon02_abs_disagreement_summary.png", 8.5, 5.8)
    }
  }

  # =============================================================================
  # 4. VERA OUTPUT INVENTORY BAR CHART
  # =============================================================================

  if (file.exists(inventory_file)) {
    df <- read.csv(inventory_file, stringsAsFactors = FALSE) %>%
      mutate(
        module = case_when(
          grepl("^addons/rasters/", .data$relative_path) ~ "Add-On rasters",
          grepl("^addons/csvs/", .data$relative_path)    ~ "Add-On CSV tables",
          grepl("^rasters/core_occupied/", .data$relative_path) ~ "Core rasters",
          grepl("^rasters/optional_diagnostics/", .data$relative_path) ~ "Optional diagnostics",
          grepl("^csvs/", .data$relative_path)           ~ "Core CSV tables",
          TRUE                                            ~ "Run metadata"
        )
      ) %>% count(.data$module, name = "n_outputs") %>% arrange(.data$n_outputs)
    df$module <- factor(df$module, levels = df$module)
    p <- ggplot(df, aes(module, n_outputs, fill = n_outputs)) +
      geom_col(colour = "black", linewidth = 0.35, width = 0.68) +
      geom_text(aes(label = n_outputs), vjust = -0.55, size = 4.1) +
      scale_fill_gradientn(colours = pal_addon, guide = "none") +
      labs(title = "VERA Output Inventory", x = NULL,
           y = "Number of inventory entries") +
      theme_csv + theme(axis.text.x = element_text(angle = 25, hjust = 1))
    save_csv_plot(p, "plot_vera_output_inventory.png", 10, 6)
  } else {
    warning("Output inventory CSV was not found: ", inventory_file)
  }

  cat(sprintf(">>> [5/5] Add-On renderer completed for %s-predictor profile.\n",
              profile_key))
}

cat("\n========================================================================\n")
cat(">>> ALL ADD-ON RENDERER PROFILES COMPLETED\n")
cat(">>> Focal taxon : ", species_label, " (", species_code, ")\n", sep = "")
cat(">>> Profiles    : 19 + 36\n")
cat(">>> Image roots :\n")
for (pk in names(results_roots)) {
  cat("      ", file.path(results_roots[[pk]], "Images", image_subfolder_name), "\n",
      sep = "")
}
cat("========================================================================\n")