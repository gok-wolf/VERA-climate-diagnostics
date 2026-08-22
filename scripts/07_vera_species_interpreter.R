#!/usr/bin/env Rscript

# =============================================================================
# VERA SPECIES INTERPRETER
# Deterministic post-analysis summaries for one taxon and paired 19/36 profiles
# =============================================================================
#
# This script reads completed Paper 1 VERA outputs. It does not recalculate or
# alter VRS, Mahalanobis distance, climatic tiers, attribution, or add-ons.
# It creates a separate interpretation layer with explicit evidence provenance.
#
# Main products, per predictor profile
#   1. Species Diagnostic Brief                    (.md)
#   2. Predictor Evidence Matrix                   (.csv; optional .xlsx)
#   3. Diagnostic Alert Register                   (.md + .csv)
#   4. Spatial Interpretation Status               (.tif + .csv)
#   5. Occurrence Review Queue                     (.csv)
#   6. Narrative Evidence Log and manifest         (.csv)
#
# Paired product
#   7. 19-versus-36 Profile Sensitivity Report     (.md + .csv + .tif)
#
# INTERPRETIVE LIMIT
# These outputs describe climatic departure, diagnostic resolution, evidence
# support, and sensitivity to predictor-profile choice. They do not estimate
# occurrence probability, habitat suitability, physiological tolerance,
# demographic performance, dispersal, causality, or conservation priority.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================

cfg <- list(
  species_code  = "Skr",
  species_label = "Sitta krueperi",

  profile_dirs = c(
    `19` = "C:/VERA/Results/19/Skr_current",
    `36` = "C:/VERA/Results/36/Skr_current"
  ),

  comparison_output_dir = "C:/VERA/Results/Skr_profile_sensitivity",

  profiles_to_run = c(19L, 36L),
  overwrite       = TRUE,
  write_xlsx      = TRUE,

  # These values mirror the canonical tutorial settings and are used only to
  # classify exported evidence. They never modify any canonical VERA output.
  asymmetry_ratio_threshold = 1.20,
  equality_tolerance        = 1e-12,

  raster_compression = c(
    "COMPRESS=LZW", "PREDICTOR=2", "TILED=YES",
    "BLOCKXSIZE=256", "BLOCKYSIZE=256", "BIGTIFF=IF_SAFER"
  ),

  # Generated evidence-bound sentences are rejected if they contain any of
  # these claims. Fixed interpretive disclaimers are written separately.
  forbidden_terms = c(
    "vulnerab", "refugi", "fatal", "impassable", "extinction",
    "adaptation", "fitness", "conservation priority", "absolute barrier",
    "absolute climatic limit", "lethal", "mortality"
  )
)

gdal_opts <- list(gdal = cfg$raster_compression)

# =============================================================================
# 2. GENERAL UTILITIES
# =============================================================================

stop_missing <- function(path, label = "Required input") {
  if (!file.exists(path)) stop(label, " was not found: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

dir_make <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_vera_csv <- function(path) {
  read.csv(stop_missing(path), stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("NA", "NaN", ""))
}

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

bool <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(z))
  out[z %in% c("TRUE", "T", "1", "YES")] <- TRUE
  out[z %in% c("FALSE", "F", "0", "NO")] <- FALSE
  out
}

fmt_num <- function(x, digits = 3) {
  if (!length(x) || !is.finite(x[1])) return("NA")
  formatC(x[1], digits = digits, format = "f")
}

fmt_pct <- function(x, digits = 2) {
  if (!length(x) || !is.finite(x[1])) return("NA")
  paste0(formatC(x[1], digits = digits, format = "f"), "%")
}

first_value <- function(x, default = NA_character_) {
  x <- x[!is.na(x)]
  if (length(x)) as.character(x[1]) else default
}

resolve_script_file <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  candidates <- character(0)
  if (length(script_arg)) candidates <- sub("^--file=", "", script_arg[1])
  source_files <- unlist(lapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile) && length(frame$ofile)) frame$ofile[1] else NULL
  }), use.names = FALSE)
  candidates <- unique(c(candidates, rev(source_files),
                         file.path(getwd(), "07_vera_species_interpreter.R")))
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[1], winslash = "/", mustWork = TRUE) else NA_character_
}

analysis_script_file <- resolve_script_file()

lint_generated_text <- function(text, context = "generated narrative") {
  if (!length(text)) return(invisible(TRUE))
  pattern <- paste(cfg$forbidden_terms, collapse = "|")
  hit <- grepl(pattern, text, ignore.case = TRUE, perl = TRUE)
  if (any(hit)) {
    stop(
      "Forbidden interpretive language detected in ", context, ":\n",
      paste(text[hit], collapse = "\n")
    )
  }
  invisible(TRUE)
}

section_wide <- function(df, section_name, id_col = "subject") {
  if (!all(c("section", id_col, "metric", "value") %in% names(df))) {
    return(tibble())
  }
  x <- df[df$section == section_name & !is.na(df[[id_col]]) &
            !is.na(df$metric), c(id_col, "metric", "value"), drop = FALSE]
  if (!nrow(x)) return(tibble())
  x[[id_col]] <- as.character(x[[id_col]])
  ids <- unique(x[[id_col]])
  out <- tibble(!!id_col := ids)
  for (metric_name in unique(x$metric)) {
    z <- x[x$metric == metric_name, c(id_col, "value"), drop = FALSE]
    z <- z[!duplicated(z[[id_col]]), , drop = FALSE]
    out[[metric_name]] <- z$value[match(ids, z[[id_col]])]
  }
  out
}

section_metric <- function(df, section_name, metric_name, default = NA_real_) {
  if (!all(c("section", "metric", "value") %in% names(df))) return(default)
  x <- df$value[df$section == section_name & df$metric == metric_name]
  # Some canonical exports place all diagnostic subsections under the common
  # section label "current_diagnostics". Exact metric names remain stable, so
  # fall back to a global metric lookup when the subsection label is absent.
  if (!length(x)) x <- df$value[df$metric == metric_name]
  if (!length(x)) default else num(x[1])
}

manifest_item <- function(run_manifest, item_name, default = NA_character_) {
  if (!all(c("section", "item", "value") %in% names(run_manifest))) return(default)
  x <- run_manifest$value[run_manifest$section == "run_metadata" &
                            run_manifest$item == item_name]
  first_value(x, default)
}

share_vector <- function(current_diag, section_name, prefix, variables) {
  out <- rep(0, length(variables))
  if (!all(c("section", "metric", "value") %in% names(current_diag))) return(out)
  expected <- paste0(prefix, variables, "_pct")
  x <- current_diag[current_diag$section == section_name,
                    c("metric", "value"), drop = FALSE]
  if (!nrow(x)) {
    x <- current_diag[current_diag$metric %in% expected,
                      c("metric", "value"), drop = FALSE]
  }
  hit <- match(expected, x$metric)
  out[!is.na(hit)] <- num(x$value[hit[!is.na(hit)]])
  out
}

write_lines_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, useBytes = TRUE)
  invisible(path)
}

rank01_values <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  if (n == 1L) out[ok] <- 0
  if (n > 1L) out[ok] <- (rank(x[ok], ties.method = "average") - 1) / (n - 1)
  out
}

rank01_raster <- function(r) {
  vals <- terra::values(r, mat = FALSE)
  terra::setValues(r, rank01_values(vals))
}

same_geometry <- function(...) {
  rasters <- list(...)
  if (length(rasters) < 2L) return(TRUE)
  all(vapply(rasters[-1], function(r) {
    isTRUE(terra::compareGeom(rasters[[1]], r, stopOnError = FALSE))
  }, logical(1)))
}

write_formatted_xlsx <- function(df, path, sheet = "Predictor evidence") {
  if (!isTRUE(cfg$write_xlsx)) return("disabled")
  if (!requireNamespace("openxlsx", quietly = TRUE)) return("openxlsx_not_installed")
  tryCatch({
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeDataTable(wb, sheet, df, tableStyle = "TableStyleMedium2")
    openxlsx::freezePane(wb, sheet, firstRow = TRUE, firstCol = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = "auto")
    header <- openxlsx::createStyle(
      fontName = "Arial",
      textDecoration = "bold",
      fgFill = "#D9EAF7",
      fontColour = "#17202A",
      halign = "center",
      border = "bottom",
      borderColour = "#7A8A99"
    )
    openxlsx::addStyle(wb, sheet, header, rows = 1,
                       cols = seq_len(ncol(df)), gridExpand = TRUE)
    body <- openxlsx::createStyle(fontName = "Arial")
    if (nrow(df) > 0L) {
      openxlsx::addStyle(
        wb, sheet, body, rows = 2:(nrow(df) + 1L),
        cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE
      )
    }
    openxlsx::saveWorkbook(wb, path, overwrite = cfg$overwrite)
    "written"
  }, error = function(e) {
    warning("XLSX export failed; CSV and all other interpreter products will continue: ",
            conditionMessage(e), call. = FALSE)
    paste0("xlsx_failed: ", conditionMessage(e))
  })
}

profile_paths <- function(profile_key, base_dir) {
  csv_dir <- file.path(base_dir, "csvs")
  core_dir <- file.path(base_dir, "rasters", "core_occupied")
  diag_dir <- file.path(base_dir, "rasters", "optional_diagnostics")
  addon_csv <- file.path(base_dir, "addons", "csvs")
  addon_ras <- file.path(base_dir, "addons", "rasters")
  prefix <- paste0(cfg$species_code, "_current")
  list(
    base = base_dir,
    lookup = file.path(csv_dir, paste0(cfg$species_code, "_00_predictor_lookup_dictionary.csv")),
    run_manifest = file.path(csv_dir, paste0(cfg$species_code, "_01_run_manifest.csv")),
    model_reference = file.path(csv_dir, paste0(cfg$species_code, "_02_model_reference.csv")),
    occurrence_diagnostics = file.path(csv_dir, paste0(cfg$species_code, "_03_occurrence_diagnostics.csv")),
    occurrence_partitions = file.path(csv_dir, paste0(cfg$species_code, "_04_occurrence_partitions.csv")),
    current_diagnostics = file.path(csv_dir, paste0(cfg$species_code, "_05_current_diagnostics.csv")),
    proximity_summary = file.path(csv_dir, paste0(cfg$species_code, "_06_proximity_summary.csv")),
    tier = file.path(core_dir, paste0(prefix, "_four_tier_status.tif")),
    mean_vrs = file.path(core_dir, paste0(prefix, "_mean_vrs.tif")),
    max_vrs = file.path(core_dir, paste0(prefix, "_max_vrs.tif")),
    mahal = file.path(core_dir, paste0(prefix, "_mahal_distance.tif")),
    primary_idx = file.path(core_dir, paste0(prefix, "_primary_stressor_idx.tif")),
    vpi = file.path(core_dir, paste0(prefix, "_vpi.tif")),
    delta = file.path(diag_dir, paste0(prefix, "_delta_vrs.tif")),
    codom = file.path(diag_dir, paste0(prefix, "_primary_co_dominance_count.tif")),
    valid_count = file.path(diag_dir, paste0(prefix, "_valid_var_count.tif")),
    divergence_summary = file.path(addon_csv, "addon02_divergence_summary.csv"),
    agreement_counts = file.path(addon_csv, "addon02_divergence_class_counts.csv"),
    tail_summary = file.path(addon_csv, "addon01_tail_direction_summary.csv"),
    divergence = file.path(addon_ras, paste0(prefix, "_mahal_vrs_divergence.tif"))
  )
}

output_paths <- function(base_dir) {
  root <- dir_make(file.path(base_dir, "interpretation"))
  list(
    root = root,
    human = dir_make(file.path(root, "human_readable")),
    tables = dir_make(file.path(root, "tables")),
    rasters = dir_make(file.path(root, "rasters")),
    machine = dir_make(file.path(root, "machine_readable"))
  )
}

# =============================================================================
# 3. LOAD ONE COMPLETED VERA PROFILE
# =============================================================================

load_profile <- function(profile_key, base_dir) {
  p <- profile_paths(profile_key, base_dir)
  required <- p[c("lookup", "run_manifest", "model_reference",
                  "occurrence_diagnostics", "occurrence_partitions",
                  "current_diagnostics", "tier", "mean_vrs", "max_vrs",
                  "mahal", "primary_idx", "vpi", "delta", "codom")]
  invisible(lapply(names(required), function(nm) stop_missing(required[[nm]], nm)))

  list(
    profile = as.character(profile_key),
    paths = p,
    out = output_paths(base_dir),
    lookup = read_vera_csv(p$lookup),
    run_manifest = read_vera_csv(p$run_manifest),
    model_reference = read_vera_csv(p$model_reference),
    occurrence_diagnostics = read_vera_csv(p$occurrence_diagnostics),
    occurrence_partitions = read_vera_csv(p$occurrence_partitions),
    current_diagnostics = read_vera_csv(p$current_diagnostics),
    proximity_summary = if (file.exists(p$proximity_summary)) read_vera_csv(p$proximity_summary) else NULL,
    divergence_summary = if (file.exists(p$divergence_summary)) read_vera_csv(p$divergence_summary) else NULL,
    agreement_counts = if (file.exists(p$agreement_counts)) read_vera_csv(p$agreement_counts) else NULL,
    tail_summary = if (file.exists(p$tail_summary)) read_vera_csv(p$tail_summary) else NULL
  )
}

# =============================================================================
# 4. PREDICTOR EVIDENCE MATRIX
# =============================================================================

build_predictor_evidence <- function(obj) {
  anchors <- section_wide(obj$model_reference, "anchor_statistics", "subject")
  if (!nrow(anchors)) stop("No anchor_statistics rows in: ", obj$paths$model_reference)
  names(anchors)[1] <- "variable"

  attribution <- section_wide(obj$model_reference,
                              "predictor_attribution_diagnostic", "subject")
  if (nrow(attribution)) names(attribution)[1] <- "variable"

  evidence <- if (nrow(attribution)) {
    full_join(anchors, attribution, by = "variable", suffix = c("", "_pilot"))
  } else anchors

  variables <- evidence$variable
  evidence$primary_share_pct <- share_vector(
    obj$current_diagnostics, "primary_stressor_share", "primary_", variables
  )
  evidence$unique_primary_share_pct <- share_vector(
    obj$current_diagnostics, "primary_unique_stressor_share",
    "primary_unique_", variables
  )
  evidence$secondary_share_pct_exported <- share_vector(
    obj$current_diagnostics, "secondary_stressor_share", "secondary_", variables
  )
  evidence$too_high_share_pct <- share_vector(
    obj$current_diagnostics, "too_high_share", "too_high_", variables
  )
  evidence$too_low_share_pct <- share_vector(
    obj$current_diagnostics, "too_low_share", "too_low_", variables
  )

  numeric_cols <- c(
    "mu", "sigma_global", "n", "n_lower", "n_upper",
    "raw_tail_scale_lower", "raw_tail_scale_upper", "raw_sd_lower", "raw_sd_upper",
    "shrink_weight_lower", "shrink_weight_upper", "sd_lower", "sd_upper",
    "bootstrap_n", "asym_ratio_raw", "asym_ratio_boot_mean",
    "asym_ratio_ci_low", "asym_ratio_ci_high", "bg_lower_tail_prop",
    "bg_upper_tail_prop", "mean_vrs_contribution", "cap_saturation_fraction",
    "mean_vrs_contribution_pct", "cumulative_vrs_contribution_pct",
    "primary_share_pct_pilot", "secondary_share_pct",
    "attribution_share_raw", "attribution_weight_pct"
  )
  logical_cols <- c(
    "fallback_lower", "fallback_upper", "asymmetry_bootstrap_stable",
    "lower_truncation_flag", "upper_truncation_flag", "truncation_flag_any",
    "contribution_coverage_member"
  )
  for (nm in intersect(numeric_cols, names(evidence))) evidence[[nm]] <- num(evidence[[nm]])
  for (nm in intersect(logical_cols, names(evidence))) evidence[[nm]] <- bool(evidence[[nm]])

  if (!"asymmetry_evidence" %in% names(evidence)) {
    low <- 1 / cfg$asymmetry_ratio_threshold
    evidence$asymmetry_evidence <- case_when(
      evidence$asymmetry_bootstrap_stable %in% TRUE ~ "stable_asymmetric",
      is.finite(evidence$asym_ratio_ci_low) & is.finite(evidence$asym_ratio_ci_high) &
        evidence$asym_ratio_ci_low >= low &
        evidence$asym_ratio_ci_high <= cfg$asymmetry_ratio_threshold ~
        "stable_near_symmetric",
      TRUE ~ "unresolved"
    )
  }

  evidence <- evidence %>%
    mutate(
      profile = obj$profile,
      directional_scaling = case_when(
        !is.finite(.data$asym_ratio_raw) ~ "not_resolved",
        .data$asym_ratio_raw > cfg$asymmetry_ratio_threshold ~
          "lower_side_scaled_more_strongly",
        .data$asym_ratio_raw < 1 / cfg$asymmetry_ratio_threshold ~
          "upper_side_scaled_more_strongly",
        TRUE ~ "near_symmetric_by_configured_ratio"
      ),
      .before = 1
    )

  order_col <- if ("unique_primary_share_pct" %in% names(evidence)) {
    evidence$unique_primary_share_pct
  } else evidence$primary_share_pct
  evidence[order(order_col, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
}

# =============================================================================
# 5. EVIDENCE-BOUND NARRATIVE AND SPECIES BRIEF
# =============================================================================

new_narrative_logger <- function(profile) {
  rows <- list()
  sentence_n <- 0L

  add <- function(text, source_file, section, subject, metric, raw_value,
                  template_id) {
    lint_generated_text(text, paste0("profile ", profile, " brief"))
    sentence_n <<- sentence_n + 1L
    sentence_id <- sprintf("%s-S%03d", profile, sentence_n)
    rows[[length(rows) + 1L]] <<- tibble(
      sentence_id = sentence_id,
      profile = profile,
      source_file = source_file,
      section = section,
      subject = subject,
      metric = metric,
      raw_value = as.character(raw_value),
      template_id = template_id,
      generated_text = text
    )
    text
  }

  get <- function() if (length(rows)) bind_rows(rows) else tibble()
  list(add = add, get = get)
}

top_row <- function(df, metric) {
  if (!metric %in% names(df)) return(NULL)
  ok <- is.finite(df[[metric]])
  if (!any(ok)) return(NULL)
  df[which.max(ifelse(ok, df[[metric]], -Inf)), , drop = FALSE]
}

build_species_brief <- function(obj, evidence) {
  log <- new_narrative_logger(obj$profile)
  add <- log$add
  cur <- obj$current_diagnostics

  n_predictors <- num(manifest_item(obj$run_manifest, "n_retained_predictors"))
  n_occ <- num(manifest_item(obj$run_manifest, "n_occurrences_total"))
  core_n <- num(manifest_item(obj$run_manifest, "n_occurrences_core_occupied"))
  periph_n <- num(manifest_item(obj$run_manifest, "n_occurrences_peripheral"))

  tier_core <- section_metric(cur, "tier_distribution", "core_climate_pct")
  tier_mod <- section_metric(cur, "tier_distribution", "moderate_departure_pct")
  tier_res <- section_metric(cur, "tier_distribution", "restriction_zone_pct")
  tier_high <- section_metric(cur, "tier_distribution", "high_extrapolative_stress_pct")
  unique_pct <- section_metric(cur, "diagnostic_overview", "primary_unique_pixel_pct")
  codom_pct <- section_metric(cur, "diagnostic_overview", "primary_codominant_pixel_pct")
  mean_vpi <- section_metric(cur, "vpi_summary", "mean_vpi")
  cap_pct <- section_metric(cur, "v10_pilot_diagnostic_summary", "any_predictor_at_cap_pct")

  top_primary <- top_row(evidence, "primary_share_pct")
  top_unique <- top_row(evidence, "unique_primary_share_pct")

  stable <- evidence %>%
    filter(.data$asymmetry_bootstrap_stable %in% TRUE,
           is.finite(.data$asym_ratio_raw), .data$asym_ratio_raw > 0) %>%
    mutate(asymmetry_strength = abs(log(.data$asym_ratio_raw))) %>%
    arrange(desc(.data$asymmetry_strength))

  trunc_vars <- evidence$variable[evidence$truncation_flag_any %in% TRUE]
  fallback_vars <- evidence$variable[
    evidence$fallback_lower %in% TRUE | evidence$fallback_upper %in% TRUE
  ]

  lines <- c(
    paste0("# VERA Species Diagnostic Brief — ", cfg$species_label),
    "",
    paste0("**Predictor profile:** ", obj$profile),
    "",
    "> This report is a deterministic interpretation aid derived from exported",
    "> VERA diagnostics. It does not estimate occurrence probability, habitat",
    "> suitability, physiological tolerance, demographic performance, dispersal,",
    "> causal range limitation, or conservation priority.",
    "",
    "## Calibration context",
    ""
  )

  if (is.finite(n_predictors) && is.finite(n_occ)) {
    lines <- c(lines, add(
      sprintf("The %s-predictor profile was calibrated from %s retained occurrence records.",
              as.integer(n_predictors), as.integer(n_occ)),
      basename(obj$paths$run_manifest), "run_metadata", cfg$species_code,
      "n_retained_predictors | n_occurrences_total",
      paste(n_predictors, n_occ, sep = " | "), "BRIEF_CALIBRATION_001"
    ))
  }
  if (is.finite(core_n) && is.finite(periph_n)) {
    lines <- c(lines, add(
      sprintf("The empirical climatic partition retained %s core occurrences and %s peripheral occurrences.",
              as.integer(core_n), as.integer(periph_n)),
      basename(obj$paths$run_manifest), "run_metadata", cfg$species_code,
      "n_occurrences_core_occupied | n_occurrences_peripheral",
      paste(core_n, periph_n, sep = " | "), "BRIEF_CALIBRATION_002"
    ))
  }

  lines <- c(lines, "", "## Landscape diagnostic summary", "")
  if (all(is.finite(c(tier_core, tier_mod, tier_res, tier_high)))) {
    lines <- c(lines, add(
      sprintf(paste0("Across valid pixels, the empirical tiers comprised %s core climate, ",
                     "%s moderate departure, %s restriction-zone, and %s high extrapolative-stress classes."),
              fmt_pct(tier_core), fmt_pct(tier_mod), fmt_pct(tier_res), fmt_pct(tier_high)),
      basename(obj$paths$current_diagnostics), "tier_distribution", "current",
      "core | moderate | restriction | high_extrapolative percentages",
      paste(tier_core, tier_mod, tier_res, tier_high, sep = " | "),
      "BRIEF_TIERS_001"
    ))
  }
  if (is.finite(mean_vpi)) {
    lines <- c(lines, add(
      sprintf("Mean Variable Proximity Index across valid pixels was %s.", fmt_num(mean_vpi, 4)),
      basename(obj$paths$current_diagnostics), "vpi_summary", "current",
      "mean_vpi", mean_vpi, "BRIEF_VPI_001"
    ))
  }
  if (is.finite(cap_pct)) {
    lines <- c(lines, add(
      sprintf("At least one predictor reached the configured VRS cap in %s of valid pixels; capped values are censored at the configured ceiling.",
              fmt_pct(cap_pct)),
      basename(obj$paths$current_diagnostics), "v10_pilot_diagnostic_summary",
      "current", "any_predictor_at_cap_pct", cap_pct, "BRIEF_CAP_001"
    ))
  }

  lines <- c(lines, "", "## Predictor identity and resolution", "")
  if (!is.null(top_primary)) {
    lines <- c(lines, add(
      sprintf("%s was the most frequent tie-broken Primary Stressor, accounting for %s of valid pixels.",
              top_primary$variable, fmt_pct(top_primary$primary_share_pct)),
      basename(obj$paths$current_diagnostics), "primary_stressor_share",
      top_primary$variable, "primary_share_pct", top_primary$primary_share_pct,
      "BRIEF_PRIMARY_001"
    ))
  }
  if (!is.null(top_unique) && is.finite(top_unique$unique_primary_share_pct)) {
    lines <- c(lines, add(
      sprintf("Within strictly unique maxima, %s had the largest mapped share (%s).",
              top_unique$variable, fmt_pct(top_unique$unique_primary_share_pct)),
      basename(obj$paths$current_diagnostics), "primary_unique_stressor_share",
      top_unique$variable, "unique_primary_share_pct",
      top_unique$unique_primary_share_pct, "BRIEF_UNIQUE_001"
    ))
  }
  if (is.finite(unique_pct) && is.finite(codom_pct)) {
    lines <- c(lines, add(
      sprintf("A strict single maximum occurred in %s of valid pixels, while %s had co-dominant maxima.",
              fmt_pct(unique_pct), fmt_pct(codom_pct)),
      basename(obj$paths$current_diagnostics), "diagnostic_overview", "current",
      "primary_unique_pixel_pct | primary_codominant_pixel_pct",
      paste(unique_pct, codom_pct, sep = " | "), "BRIEF_RESOLUTION_001"
    ))
  }

  lines <- c(lines, "", "## Directional calibration", "")
  if (nrow(stable)) {
    for (i in seq_len(min(3L, nrow(stable)))) {
      r <- stable[i, ]
      direction_text <- if (r$asym_ratio_raw > 1) {
        "lower-side departures are scaled more strongly than upper-side departures"
      } else {
        "upper-side departures are scaled more strongly than lower-side departures"
      }
      lines <- c(lines, add(
        sprintf("%s had an upper-to-lower occupied-breadth ratio of %s (bootstrap interval %s–%s); under VERA scaling, %s.",
                r$variable, fmt_num(r$asym_ratio_raw, 2),
                fmt_num(r$asym_ratio_ci_low, 2), fmt_num(r$asym_ratio_ci_high, 2),
                direction_text),
        basename(obj$paths$model_reference), "anchor_statistics", r$variable,
        "asym_ratio_raw | asym_ratio_ci_low | asym_ratio_ci_high",
        paste(r$asym_ratio_raw, r$asym_ratio_ci_low, r$asym_ratio_ci_high,
              sep = " | "), "BRIEF_ASYMMETRY_001"
      ))
    }
  } else {
    lines <- c(lines, add(
      "No predictor was exported with a bootstrap-stable asymmetric breadth flag.",
      basename(obj$paths$model_reference), "anchor_statistics", "all_predictors",
      "asymmetry_bootstrap_stable", "0 stable predictors",
      "BRIEF_ASYMMETRY_NONE_001"
    ))
  }

  lines <- c(lines, "", "## Evidence cautions", "")
  if (length(trunc_vars)) {
    lines <- c(lines, add(
      sprintf("Background-boundary contact was flagged for: %s. Tail interpretation for these predictors is bounded by the supplied raster domain.",
              paste(trunc_vars, collapse = ", ")),
      basename(obj$paths$model_reference), "anchor_statistics",
      paste(trunc_vars, collapse = " | "), "truncation_flag_any", "TRUE",
      "BRIEF_TRUNCATION_001"
    ))
  } else {
    lines <- c(lines, add(
      "No predictor had an exported background-truncation flag.",
      basename(obj$paths$model_reference), "anchor_statistics", "all_predictors",
      "truncation_flag_any", "0 flagged predictors",
      "BRIEF_TRUNCATION_NONE_001"
    ))
  }
  if (length(fallback_vars)) {
    lines <- c(lines, add(
      sprintf("Global-scale fallback was used for at least one tail of: %s. Predictor-specific tail interpretation should retain this qualification.",
              paste(fallback_vars, collapse = ", ")),
      basename(obj$paths$model_reference), "anchor_statistics",
      paste(fallback_vars, collapse = " | "), "fallback_lower | fallback_upper",
      "TRUE", "BRIEF_FALLBACK_001"
    ))
  }

  if (!is.null(obj$divergence_summary) &&
      all(c("metric", "value") %in% names(obj$divergence_summary))) {
    div <- obj$divergence_summary
    med_abs <- num(div$value[div$metric == "median_abs_disagreement"])
    within <- num(div$value[div$metric == "percent_within_threshold"])
    if (length(med_abs) && is.finite(med_abs[1])) {
      lines <- c(lines, "", "## Cross-geometry comparison", "", add(
        sprintf("Median absolute VRS–Mahalanobis percentile-rank disagreement was %s%s.",
                fmt_num(med_abs[1], 4),
                if (length(within) && is.finite(within[1])) {
                  paste0("; ", fmt_pct(within[1]), " of valid pixels were within the configured divergence threshold")
                } else ""),
        basename(obj$paths$divergence_summary), "addon02_divergence_summary",
        "current", "median_abs_disagreement | percent_within_threshold",
        paste(med_abs[1], ifelse(length(within), within[1], NA), sep = " | "),
        "BRIEF_DIVERGENCE_001"
      ))
    }
  }

  lines <- c(
    lines, "", "## Reading rule", "",
    "Interpret magnitude, identity, uniqueness, direction, covariance-aware context, and diagnostic cautions jointly. No single VERA surface is a standalone model-quality or ecological-performance score.",
    ""
  )

  list(lines = lines, evidence_log = log$get())
}

# =============================================================================
# 6. DIAGNOSTIC ALERT REGISTER
# =============================================================================

build_alerts <- function(obj, evidence) {
  alerts <- list()
  alert_n <- 0L
  add_alert <- function(severity, scope, subject, code, message,
                        source_file, section, source_subject, metric, raw_value) {
    lint_generated_text(message, paste0("profile ", obj$profile, " alert"))
    alert_n <<- alert_n + 1L
    alerts[[length(alerts) + 1L]] <<- tibble(
      alert_id = sprintf("%s-A%03d", obj$profile, alert_n),
      profile = obj$profile,
      severity = severity,
      scope = scope,
      subject = subject,
      alert_code = code,
      message = message,
      source_file = source_file,
      section = section,
      source_subject = source_subject,
      metric = metric,
      raw_value = as.character(raw_value)
    )
  }

  validations <- obj$run_manifest %>%
    filter(.data$section == "validation_checks",
           toupper(as.character(.data$value)) == "FAIL")
  if (nrow(validations)) {
    for (i in seq_len(nrow(validations))) {
      add_alert(
        "REVIEW_REQUIRED", "run", validations$item[i], "VALIDATION_FAIL",
        paste0("An exported internal validation check failed: ",
               validations$item[i], ". Review the run before interpreting outputs."),
        basename(obj$paths$run_manifest), "validation_checks",
        validations$item[i], "validation_check_detail",
        validations$detail[i] %||% validations$value[i]
      )
    }
  }

  for (i in seq_len(nrow(evidence))) {
    r <- evidence[i, ]
    if (isTRUE(r$lower_truncation_flag) || isTRUE(r$upper_truncation_flag)) {
      sides <- paste(c(if (isTRUE(r$lower_truncation_flag)) "lower" else NULL,
                       if (isTRUE(r$upper_truncation_flag)) "upper" else NULL),
                     collapse = " and ")
      add_alert(
        "CAUTION", "predictor", r$variable, "BACKGROUND_BOUNDARY_CONTACT",
        sprintf("%s-tail background-boundary contact was flagged for %s; the corresponding occupied breadth is bounded by the supplied raster domain.",
                sides, r$variable),
        basename(obj$paths$model_reference), "anchor_statistics", r$variable,
        "lower_truncation_flag | upper_truncation_flag",
        paste(r$lower_truncation_flag, r$upper_truncation_flag, sep = " | ")
      )
    }
    if (isTRUE(r$fallback_lower) || isTRUE(r$fallback_upper)) {
      sides <- paste(c(if (isTRUE(r$fallback_lower)) "lower" else NULL,
                       if (isTRUE(r$fallback_upper)) "upper" else NULL),
                     collapse = " and ")
      add_alert(
        "CAUTION", "predictor", r$variable, "TAIL_SCALE_FALLBACK",
        sprintf("The %s-tail scale for %s used the global occurrence SD fallback.",
                sides, r$variable),
        basename(obj$paths$model_reference), "anchor_statistics", r$variable,
        "fallback_lower | fallback_upper",
        paste(r$fallback_lower, r$fallback_upper, sep = " | ")
      )
    }
    extreme_ratio <- is.finite(r$asym_ratio_raw) &&
      (r$asym_ratio_raw > cfg$asymmetry_ratio_threshold ||
         r$asym_ratio_raw < 1 / cfg$asymmetry_ratio_threshold)
    if (extreme_ratio && !isTRUE(r$asymmetry_bootstrap_stable)) {
      add_alert(
        "CAUTION", "predictor", r$variable, "ASYMMETRY_UNRESOLVED",
        sprintf("%s has a directional breadth ratio outside the configured near-symmetry band, but bootstrap stability was not established.",
                r$variable),
        basename(obj$paths$model_reference), "anchor_statistics", r$variable,
        "asym_ratio_raw | asymmetry_bootstrap_stable",
        paste(r$asym_ratio_raw, r$asymmetry_bootstrap_stable, sep = " | ")
      )
    }
    if ("cap_saturation_fraction" %in% names(r) &&
        is.finite(r$cap_saturation_fraction) && r$cap_saturation_fraction > 0) {
      add_alert(
        "NOTE", "predictor", r$variable, "PREDICTOR_CAP_ENGAGEMENT",
        sprintf("%s reached the configured VRS cap in %s of valid pixels.",
                r$variable, fmt_pct(100 * r$cap_saturation_fraction)),
        basename(obj$paths$model_reference), "predictor_attribution_diagnostic",
        r$variable, "cap_saturation_fraction", r$cap_saturation_fraction
      )
    }
  }

  shape <- obj$occurrence_diagnostics %>%
    filter(.data$section == "distribution_shape_flagged_predictors")
  if (nrow(shape)) {
    for (i in seq_len(nrow(shape))) {
      add_alert(
        "NOTE", "predictor", shape$subject[i], "DISTRIBUTION_SHAPE_SCREEN",
        sprintf("%s triggered the exported skewness/kurtosis distribution-shape screen; this is not a formal multimodality test.",
                shape$subject[i]),
        basename(obj$paths$occurrence_diagnostics), shape$section[i],
        shape$subject[i], shape$metric[i], shape$value[i]
      )
    }
  }

  if (!length(alerts)) {
    add_alert(
      "NOTE", "run", cfg$species_code, "NO_EXPORTED_FLAGS",
      "No validation failure, fallback, truncation, unresolved directional ratio, cap engagement, or distribution-shape flag was detected in the exported fields inspected by this interpreter.",
      "multiple", "multiple", cfg$species_code, "multiple", "none"
    )
  }

  bind_rows(alerts) %>%
    mutate(severity = factor(.data$severity,
                             levels = c("REVIEW_REQUIRED", "CAUTION", "NOTE"))) %>%
    arrange(.data$severity, .data$scope, .data$subject) %>%
    mutate(severity = as.character(.data$severity))
}

alerts_to_markdown <- function(alerts, profile) {
  lines <- c(
    paste0("# VERA Diagnostic Alert Register — ", cfg$species_label),
    "", paste0("**Predictor profile:** ", profile), "",
    "> Alerts describe exported diagnostic evidence and data limitations. They",
    "> are not ecological risk classes or model acceptance decisions.", ""
  )
  for (severity in c("REVIEW_REQUIRED", "CAUTION", "NOTE")) {
    x <- alerts[alerts$severity == severity, , drop = FALSE]
    if (!nrow(x)) next
    lines <- c(lines, paste0("## ", severity), "")
    for (i in seq_len(nrow(x))) {
      lines <- c(lines, paste0("- **", x$subject[i], " — ", x$alert_code[i],
                               ":** ", x$message[i]))
    }
    lines <- c(lines, "")
  }
  lines
}

# =============================================================================
# 7. SPATIAL INTERPRETATION STATUS
# =============================================================================

build_spatial_status <- function(obj) {
  p <- obj$paths
  required <- c(p$tier, p$max_vrs, p$codom)
  if (!all(file.exists(required))) {
    return(list(written = FALSE, counts = tibble(), files = character(0)))
  }
  tier <- rast(p$tier)
  max_vrs <- rast(p$max_vrs)
  codom <- rast(p$codom)
  if (!same_geometry(tier, max_vrs, codom)) {
    stop("Spatial status raster geometry mismatch for profile ", obj$profile)
  }

  z2_cap <- num(manifest_item(obj$run_manifest, "z2_cap"))
  if (!is.finite(z2_cap)) z2_cap <- 16
  n_predictors <- num(manifest_item(obj$run_manifest, "n_retained_predictors"))

  status <- ifel(!is.na(tier), 1, NA)
  status <- ifel(!is.na(tier) & tier >= 3 & codom <= 1, 2, status)
  status <- ifel(!is.na(tier) & tier >= 3 & codom > 1, 3, status)
  status <- ifel(!is.na(max_vrs) & max_vrs >= z2_cap - cfg$equality_tolerance,
                 4, status)

  if (file.exists(p$valid_count) && is.finite(n_predictors)) {
    valid_count <- rast(p$valid_count)
    if (!same_geometry(tier, valid_count)) {
      stop("Valid-count raster geometry mismatch for profile ", obj$profile)
    }
    status <- ifel(!is.na(valid_count) & valid_count < n_predictors, 5, status)
  }
  names(status) <- "interpretation_status"

  labels <- c(
    "Lower empirical departure (tiers 1-2)",
    "Elevated departure with strict single maximum",
    "Elevated departure with co-dominant maximum",
    "VRS cap engaged (censored at configured ceiling)",
    "Reduced valid-predictor coverage"
  )
  factor_status <- as.factor(status)
  levels(factor_status) <- data.frame(value = seq_along(labels), label = labels)

  tif <- file.path(obj$out$rasters,
                   paste0(cfg$species_code, "_", obj$profile,
                          "_interpretation_status.tif"))
  writeRaster(factor_status, tif, overwrite = cfg$overwrite, wopt = gdal_opts)

  vals <- terra::values(status, mat = FALSE)
  vals <- vals[is.finite(vals)]
  counts <- tibble(
    profile = obj$profile,
    class_id = seq_along(labels),
    class_name = labels,
    pixel_count = as.integer(tabulate(as.integer(vals), nbins = length(labels))),
    percent_valid = if (length(vals)) {
      100 * tabulate(as.integer(vals), nbins = length(labels)) / length(vals)
    } else NA_real_,
    precedence_note = paste(
      "Reduced coverage overrides cap engagement; cap engagement overrides",
      "co-dominance; co-dominance overrides strict-maximum class."
    )
  )
  count_file <- file.path(obj$out$tables,
                          paste0(cfg$species_code, "_", obj$profile,
                                 "_interpretation_status_counts.csv"))
  write.csv(counts, count_file, row.names = FALSE)
  list(written = TRUE, counts = counts, files = c(tif, count_file))
}

# =============================================================================
# 8. OCCURRENCE REVIEW QUEUE
# =============================================================================

detect_coordinate_columns <- function(df) {
  lon_candidates <- c("decimalLongitude", "Longitude", "longitude", "lon", "x")
  lat_candidates <- c("decimalLatitude", "Latitude", "latitude", "lat", "y")
  lon <- lon_candidates[lon_candidates %in% names(df)]
  lat <- lat_candidates[lat_candidates %in% names(df)]
  if (length(lon) && length(lat)) c(lon = lon[1], lat = lat[1]) else NULL
}

build_occurrence_queue <- function(obj) {
  occ <- obj$occurrence_partitions
  if (!nrow(occ)) return(list(data = tibble(), files = character(0)))
  occ$profile <- obj$profile
  coords <- detect_coordinate_columns(occ)

  raster_candidates <- c(
    mean_vrs = obj$paths$mean_vrs,
    max_vrs = obj$paths$max_vrs,
    mahal_surface = obj$paths$mahal,
    tier = obj$paths$tier,
    primary_index = obj$paths$primary_idx,
    delta_vrs = obj$paths$delta,
    primary_co_dominance_count = obj$paths$codom
  )
  raster_candidates <- raster_candidates[file.exists(raster_candidates)]

  if (!is.null(coords) && length(raster_candidates)) {
    xy_ok <- is.finite(num(occ[[coords["lon"]]])) &
      is.finite(num(occ[[coords["lat"]]]))
    extracted <- as.data.frame(matrix(NA_real_, nrow(occ), length(raster_candidates)))
    names(extracted) <- names(raster_candidates)
    if (any(xy_ok)) {
      pts_df <- occ[xy_ok, , drop = FALSE]
      pts_df[[coords["lon"]]] <- num(pts_df[[coords["lon"]]])
      pts_df[[coords["lat"]]] <- num(pts_df[[coords["lat"]]])
      pts <- vect(pts_df, geom = unname(coords), crs = "EPSG:4326")
      rr <- rast(unname(raster_candidates))
      names(rr) <- names(raster_candidates)
      if (!same.crs(pts, rr)) pts <- project(pts, crs(rr))
      extracted[xy_ok, ] <- as.data.frame(extract(rr, pts, ID = FALSE))
    }
    occ <- bind_cols(occ, extracted)
  }

  if ("primary_index" %in% names(occ) &&
      all(c("stressor_index", "variable") %in% names(obj$lookup))) {
    occ$primary_predictor <- obj$lookup$variable[
      match(as.integer(round(num(occ$primary_index))), obj$lookup$stressor_index)
    ]
  }

  distance_col <- if ("mahal_distance" %in% names(occ)) {
    "mahal_distance"
  } else if ("mahal_surface" %in% names(occ)) {
    "mahal_surface"
  } else NULL

  if (!is.null(distance_col)) {
    d <- num(occ[[distance_col]])
    occ$mahal_distance_percentile <- rank01_values(d)
    occ$review_order <- rank(-d, ties.method = "first", na.last = "keep")
  } else {
    occ$mahal_distance_percentile <- NA_real_
    occ$review_order <- seq_len(nrow(occ))
  }

  occ$review_context <- case_when(
    occ$partition == "peripheral" ~
      "Peripheral climatic position; verify provenance before strong interpretation",
    occ$partition == "core_occupied" ~
      "Core occupied climatic position",
    TRUE ~ "Partition label unavailable"
  )
  occ <- occ[order(occ$review_order, na.last = TRUE), , drop = FALSE]

  out_file <- file.path(obj$out$tables,
                        paste0(cfg$species_code, "_", obj$profile,
                               "_occurrence_review_queue.csv"))
  write.csv(occ, out_file, row.names = FALSE)
  list(data = occ, files = out_file)
}

# =============================================================================
# 9. PROFILE-LEVEL MANIFEST
# =============================================================================

file_inventory <- function(files, role, generated_utc) {
  files <- unique(files[file.exists(files)])
  if (!length(files)) return(tibble())
  tibble(
    role = role,
    file = normalizePath(files, winslash = "/", mustWork = TRUE),
    size_bytes = as.numeric(file.info(files)$size),
    md5 = unname(tools::md5sum(files)),
    generated_utc = generated_utc
  )
}

write_interpreter_manifest <- function(obj, input_files, output_files, xlsx_status) {
  utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  script_tbl <- if (!is.na(analysis_script_file) && file.exists(analysis_script_file)) {
    file_inventory(analysis_script_file, "interpreter_script", utc)
  } else tibble(
    role = "interpreter_script", file = NA_character_, size_bytes = NA_real_,
    md5 = NA_character_, generated_utc = utc
  )
  manifest <- bind_rows(
    script_tbl,
    file_inventory(input_files, "input", utc),
    file_inventory(output_files, "output", utc)
  ) %>%
    mutate(
      species_code = cfg$species_code,
      species_label = cfg$species_label,
      profile = obj$profile,
      xlsx_status = xlsx_status,
      .before = 1
    )
  path <- file.path(obj$out$machine,
                    paste0(cfg$species_code, "_", obj$profile,
                           "_interpreter_manifest.csv"))
  write.csv(manifest, path, row.names = FALSE)
  session_path <- file.path(obj$out$machine,
                            paste0(cfg$species_code, "_", obj$profile,
                                   "_interpreter_sessionInfo.txt"))
  writeLines(capture.output(sessionInfo()), session_path)
  c(path, session_path)
}

# =============================================================================
# 10. RUN ONE PROFILE
# =============================================================================

run_profile_interpreter <- function(profile_key, base_dir) {
  cat("\n>>> VERA Species Interpreter — profile ", profile_key, "\n", sep = "")
  obj <- load_profile(profile_key, base_dir)
  evidence <- build_predictor_evidence(obj)

  evidence_csv <- file.path(obj$out$tables,
                            paste0(cfg$species_code, "_", obj$profile,
                                   "_predictor_evidence_matrix.csv"))
  write.csv(evidence, evidence_csv, row.names = FALSE)
  evidence_xlsx <- file.path(obj$out$tables,
                             paste0(cfg$species_code, "_", obj$profile,
                                    "_predictor_evidence_matrix.xlsx"))
  xlsx_status <- write_formatted_xlsx(evidence, evidence_xlsx)

  brief <- build_species_brief(obj, evidence)
  brief_file <- file.path(obj$out$human,
                          paste0(cfg$species_code, "_", obj$profile,
                                 "_species_diagnostic_brief.md"))
  write_lines_utf8(brief$lines, brief_file)

  evidence_log_file <- file.path(obj$out$machine,
                                 paste0(cfg$species_code, "_", obj$profile,
                                        "_narrative_evidence.csv"))
  write.csv(brief$evidence_log, evidence_log_file, row.names = FALSE)

  alerts <- build_alerts(obj, evidence)
  alerts_csv <- file.path(obj$out$tables,
                          paste0(cfg$species_code, "_", obj$profile,
                                 "_diagnostic_alerts.csv"))
  alerts_md <- file.path(obj$out$human,
                         paste0(cfg$species_code, "_", obj$profile,
                                "_diagnostic_alerts.md"))
  write.csv(alerts, alerts_csv, row.names = FALSE)
  write_lines_utf8(alerts_to_markdown(alerts, obj$profile), alerts_md)

  spatial <- build_spatial_status(obj)
  occurrence <- build_occurrence_queue(obj)

  output_files <- c(
    evidence_csv,
    if (identical(xlsx_status, "written")) evidence_xlsx else character(0),
    brief_file, evidence_log_file, alerts_csv, alerts_md,
    spatial$files, occurrence$files
  )
  input_files <- unlist(obj$paths[c(
    "lookup", "run_manifest", "model_reference", "occurrence_diagnostics",
    "occurrence_partitions", "current_diagnostics", "proximity_summary",
    "tier", "mean_vrs", "max_vrs", "mahal", "primary_idx", "vpi",
    "delta", "codom", "valid_count", "divergence_summary",
    "agreement_counts", "tail_summary"
  )], use.names = FALSE)
  manifest_files <- write_interpreter_manifest(
    obj, input_files, output_files, xlsx_status
  )

  cat("    Brief             : ", brief_file, "\n", sep = "")
  cat("    Evidence matrix   : ", evidence_csv, "\n", sep = "")
  cat("    XLSX status       : ", xlsx_status, "\n", sep = "")
  cat("    Alerts            : ", alerts_csv, "\n", sep = "")
  cat("    Interpretation dir: ", obj$out$root, "\n", sep = "")

  list(
    obj = obj,
    evidence = evidence,
    alerts = alerts,
    spatial = spatial,
    occurrence = occurrence$data,
    output_files = c(output_files, manifest_files),
    xlsx_status = xlsx_status
  )
}

# =============================================================================
# 11. 19-VERSUS-36 PROFILE SENSITIVITY
# =============================================================================

lookup_names <- function(values, lookup) {
  out <- rep(NA_character_, length(values))
  ok <- is.finite(values)
  out[ok] <- lookup$variable[
    match(as.integer(round(values[ok])), lookup$stressor_index)
  ]
  out
}

build_profile_comparison <- function(result19, result36) {
  out_root <- dir_make(cfg$comparison_output_dir)
  human <- dir_make(file.path(out_root, "human_readable"))
  tables <- dir_make(file.path(out_root, "tables"))
  rasters <- dir_make(file.path(out_root, "rasters"))
  machine <- dir_make(file.path(out_root, "machine_readable"))

  p19 <- result19$obj$paths
  p36 <- result36$obj$paths
  mean19 <- rast(p19$mean_vrs); mean36 <- rast(p36$mean_vrs)
  tier19 <- rast(p19$tier); tier36 <- rast(p36$tier)
  pri19 <- rast(p19$primary_idx); pri36 <- rast(p36$primary_idx)
  cod19 <- rast(p19$codom); cod36 <- rast(p36$codom)
  if (!same_geometry(mean19, mean36, tier19, tier36, pri19, pri36,
                     cod19, cod36)) {
    stop("The 19- and 36-predictor rasters do not share one geometry.")
  }

  v19 <- terra::values(mean19, mat = FALSE)
  v36 <- terra::values(mean36, mat = FALSE)
  common <- is.finite(v19) & is.finite(v36)
  r19 <- rep(NA_real_, length(v19))
  r36 <- rep(NA_real_, length(v36))
  r19[common] <- rank01_values(v19[common])
  r36[common] <- rank01_values(v36[common])
  rank_diff_vals <- r36 - r19
  abs_rank_diff <- abs(rank_diff_vals[common])

  rank_diff <- terra::setValues(mean19, rank_diff_vals)
  names(rank_diff) <- "mean_vrs_rank36_minus_rank19"
  rank_diff_file <- file.path(rasters,
                              paste0(cfg$species_code,
                                     "_mean_vrs_rank36_minus_rank19.tif"))
  writeRaster(rank_diff, rank_diff_file, overwrite = cfg$overwrite, wopt = gdal_opts)

  t19 <- terra::values(tier19, mat = FALSE)
  t36 <- terra::values(tier36, mat = FALSE)
  tier_ok <- is.finite(t19) & is.finite(t36)
  tier_diff_vals <- rep(NA_real_, length(t19))
  tier_diff_vals[tier_ok] <- t36[tier_ok] - t19[tier_ok]
  tier_diff <- terra::setValues(mean19, tier_diff_vals)
  names(tier_diff) <- "tier36_minus_tier19"
  tier_diff_file <- file.path(rasters,
                              paste0(cfg$species_code, "_tier36_minus_tier19.tif"))
  writeRaster(tier_diff, tier_diff_file, overwrite = cfg$overwrite, wopt = gdal_opts)

  tier_cross <- as.data.frame(table(
    tier_19 = factor(t19[tier_ok], levels = 1:4),
    tier_36 = factor(t36[tier_ok], levels = 1:4)
  ))
  names(tier_cross)[3] <- "pixel_count"
  tier_cross$percent_common <- if (sum(tier_cross$pixel_count) > 0) {
    100 * tier_cross$pixel_count / sum(tier_cross$pixel_count)
  } else NA_real_
  tier_cross_file <- file.path(tables,
                               paste0(cfg$species_code,
                                      "_19_vs_36_tier_cross_tabulation.csv"))
  write.csv(tier_cross, tier_cross_file, row.names = FALSE)

  pv19 <- terra::values(pri19, mat = FALSE)
  pv36 <- terra::values(pri36, mat = FALSE)
  cv19 <- terra::values(cod19, mat = FALSE)
  cv36 <- terra::values(cod36, mat = FALSE)
  primary_ok <- is.finite(pv19) & is.finite(pv36)
  name19 <- lookup_names(pv19, result19$obj$lookup)
  name36 <- lookup_names(pv36, result36$obj$lookup)
  names_ok <- primary_ok & !is.na(name19) & !is.na(name36)
  cod_ok <- is.finite(cv19) & is.finite(cv36)
  both_unique <- names_ok & cod_ok & cv19 == 1 & cv36 == 1
  unique36 <- names_ok & is.finite(cv36) & cv36 == 1
  is_bio36 <- grepl("^Bio[0-9]+$", name36)
  primary_class <- rep(NA_integer_, length(pv19))
  primary_class[names_ok & name19 == name36] <- 1L
  primary_class[names_ok & name19 != name36 & is_bio36] <- 2L
  primary_class[names_ok & !is_bio36] <- 3L
  primary_rast <- terra::setValues(mean19, primary_class)
  primary_labels <- c(
    "Same tie-broken bioclimatic Primary Stressor",
    "Different tie-broken bioclimatic Primary Stressor",
    "ENVIREM predictor is tie-broken Primary in 36-profile"
  )
  primary_factor <- as.factor(primary_rast)
  levels(primary_factor) <- data.frame(value = 1:3, label = primary_labels)
  primary_file <- file.path(rasters,
                            paste0(cfg$species_code,
                                   "_19_vs_36_primary_transition_class.tif"))
  writeRaster(primary_factor, primary_file, overwrite = cfg$overwrite,
              wopt = gdal_opts)

  transitions <- tibble(
    primary_19 = name19[names_ok],
    primary_36 = name36[names_ok]
  ) %>%
    count(.data$primary_19, .data$primary_36, name = "pixel_count") %>%
    mutate(
      transition_type = ifelse(.data$primary_19 == .data$primary_36,
                               "same", "changed"),
      percent_common = 100 * .data$pixel_count / sum(.data$pixel_count)
    ) %>%
    arrange(desc(.data$pixel_count))
  transition_file <- file.path(tables,
                               paste0(cfg$species_code,
                                      "_19_vs_36_primary_transitions.csv"))
  write.csv(transitions, transition_file, row.names = FALSE)

  strict_unique_transitions <- tibble(
    primary_19 = name19[both_unique],
    primary_36 = name36[both_unique]
  ) %>%
    count(.data$primary_19, .data$primary_36, name = "pixel_count") %>%
    mutate(
      transition_type = ifelse(.data$primary_19 == .data$primary_36,
                               "same", "changed"),
      percent_both_unique = if (sum(.data$pixel_count) > 0) {
        100 * .data$pixel_count / sum(.data$pixel_count)
      } else NA_real_
    ) %>%
    arrange(desc(.data$pixel_count))
  strict_unique_transition_file <- file.path(
    tables,
    paste0(cfg$species_code,
           "_19_vs_36_strict_unique_primary_transitions.csv")
  )
  write.csv(strict_unique_transitions, strict_unique_transition_file,
            row.names = FALSE)

  unique_agreement <- if (any(cod_ok)) {
    100 * mean((cv19[cod_ok] == 1) == (cv36[cod_ok] == 1))
  } else NA_real_

  summary <- tibble(
    species_code = cfg$species_code,
    species_label = cfg$species_label,
    metric = c(
      "n_common_valid_pixels",
      "mean_vrs_spearman_correlation",
      "mean_absolute_rank_difference",
      "p90_absolute_rank_difference",
      "tier_agreement_percent",
      "tie_broken_primary_identity_agreement_percent",
      "tie_broken_envirem_primary_share_in_36_percent",
      "tie_broken_bioclimatic_primary_share_in_36_percent",
      "unique_vs_codominant_status_agreement_percent",
      "n_both_profiles_strict_unique_pixels",
      "strict_unique_primary_identity_agreement_percent",
      "n_36_profile_strict_unique_pixels",
      "envirem_share_within_36_strict_unique_pixels_percent",
      "envirem_strict_unique_share_of_common_primary_pixels_percent"
    ),
    value = c(
      sum(common),
      if (sum(common) > 2) cor(v19[common], v36[common], method = "spearman") else NA_real_,
      if (length(abs_rank_diff)) mean(abs_rank_diff) else NA_real_,
      if (length(abs_rank_diff)) quantile(abs_rank_diff, 0.90, names = FALSE) else NA_real_,
      if (any(tier_ok)) 100 * mean(t19[tier_ok] == t36[tier_ok]) else NA_real_,
      if (any(names_ok)) 100 * mean(name19[names_ok] == name36[names_ok]) else NA_real_,
      if (any(names_ok)) 100 * mean(!is_bio36[names_ok]) else NA_real_,
      if (any(names_ok)) 100 * mean(is_bio36[names_ok]) else NA_real_,
      unique_agreement,
      sum(both_unique),
      if (any(both_unique)) {
        100 * mean(name19[both_unique] == name36[both_unique])
      } else NA_real_,
      sum(unique36),
      if (any(unique36)) 100 * mean(!is_bio36[unique36]) else NA_real_,
      if (any(names_ok)) {
        100 * sum(unique36 & !is_bio36) / sum(names_ok)
      } else NA_real_
    ),
    interpretation = c(
      "Pixels with finite mean VRS in both profiles.",
      "Rank association of raw mean-VRS values; use with rank-difference maps.",
      "Average absolute difference between landscape percentile ranks.",
      "Ninetieth percentile of absolute landscape-rank difference.",
      "Share of common pixels assigned to the same empirical climatic tier.",
      "Share of common pixels with the same deterministically tie-broken named Primary Stressor.",
      "Share of common pixels assigned an ENVIREM Primary in the 36-profile after deterministic tie breaking.",
      "Share of common pixels assigned a bioclimatic Primary in the 36-profile after deterministic tie breaking.",
      "Agreement in strict-single versus co-dominant maximum status.",
      "Common pixels with a strict single maximum in both profiles.",
      "Primary-identity agreement restricted to pixels with a strict single maximum in both profiles.",
      "Common pixels with a strict single maximum in the 36-profile.",
      "ENVIREM share among pixels with a strict single maximum in the 36-profile.",
      "Share of all common Primary pixels that have a strict-single ENVIREM Primary in the 36-profile."
    )
  )
  summary_file <- file.path(tables,
                            paste0(cfg$species_code,
                                   "_19_vs_36_profile_sensitivity_summary.csv"))
  write.csv(summary, summary_file, row.names = FALSE)

  get_metric <- function(name) summary$value[summary$metric == name][1]
  report_sentences <- c(
    sprintf("Mean-VRS ranks had a Spearman association of %s across common valid pixels.",
            fmt_num(get_metric("mean_vrs_spearman_correlation"), 3)),
    sprintf("The two profiles assigned the same empirical climatic tier to %s of common pixels.",
            fmt_pct(get_metric("tier_agreement_percent"))),
    sprintf("The deterministically tie-broken named Primary Stressor was identical in %s of common pixels.",
            fmt_pct(get_metric("tie_broken_primary_identity_agreement_percent"))),
    sprintf("After deterministic tie breaking, an ENVIREM predictor was assigned as Primary in %s of common pixels under the 36-predictor profile.",
            fmt_pct(get_metric("tie_broken_envirem_primary_share_in_36_percent"))),
    sprintf("Strict-single versus co-dominant maximum status agreed in %s of common pixels.",
            fmt_pct(get_metric("unique_vs_codominant_status_agreement_percent"))),
    sprintf("Among pixels with a strict single maximum in both profiles, the named Primary Stressor was identical in %s.",
            fmt_pct(get_metric("strict_unique_primary_identity_agreement_percent"))),
    sprintf("Among pixels with a strict single maximum in the 36-predictor profile, an ENVIREM predictor was Primary in %s.",
            fmt_pct(get_metric("envirem_share_within_36_strict_unique_pixels_percent")))
  )
  lint_generated_text(report_sentences, "19-versus-36 sensitivity report")
  report <- c(
    paste0("# VERA 19-versus-36 Profile Sensitivity — ", cfg$species_label),
    "",
    "> This report describes sensitivity to predictor-profile choice. It does",
    "> not identify one profile as ecological truth or external validation.",
    "", "## Summary", "",
    paste0("- ", report_sentences),
    "", "## Reading rule", "",
    "Use the summary table together with the rank-difference, tier-difference, and Primary-transition rasters. The full transition product records deterministic tie-broken assignments; use the strict-unique transition table when the interpretation requires an unshared maximum. Raw mean-VRS magnitudes are not treated as directly interchangeable across different predictor counts.",
    ""
  )
  report_file <- file.path(human,
                           paste0(cfg$species_code,
                                  "_19_vs_36_profile_sensitivity.md"))
  write_lines_utf8(report, report_file)

  evidence_log <- tibble(
    sentence_id = sprintf("19v36-S%03d", seq_along(report_sentences)),
    source_file = basename(summary_file),
    metric = c(
      "mean_vrs_spearman_correlation", "tier_agreement_percent",
      "tie_broken_primary_identity_agreement_percent",
      "tie_broken_envirem_primary_share_in_36_percent",
      "unique_vs_codominant_status_agreement_percent",
      "strict_unique_primary_identity_agreement_percent",
      "envirem_share_within_36_strict_unique_pixels_percent"
    ),
    raw_value = c(
      get_metric("mean_vrs_spearman_correlation"), get_metric("tier_agreement_percent"),
      get_metric("tie_broken_primary_identity_agreement_percent"),
      get_metric("tie_broken_envirem_primary_share_in_36_percent"),
      get_metric("unique_vs_codominant_status_agreement_percent"),
      get_metric("strict_unique_primary_identity_agreement_percent"),
      get_metric("envirem_share_within_36_strict_unique_pixels_percent")
    ),
    template_id = c(
      "PROFILE_CORRELATION_001", "PROFILE_TIER_001", "PROFILE_PRIMARY_001",
      "PROFILE_ENVIREM_001", "PROFILE_RESOLUTION_001",
      "PROFILE_STRICT_UNIQUE_PRIMARY_001",
      "PROFILE_STRICT_UNIQUE_ENVIREM_001"
    ),
    generated_text = report_sentences
  )
  evidence_file <- file.path(machine,
                             paste0(cfg$species_code,
                                    "_19_vs_36_narrative_evidence.csv"))
  write.csv(evidence_log, evidence_file, row.names = FALSE)

  files <- c(summary_file, tier_cross_file, transition_file,
             strict_unique_transition_file, report_file, evidence_file,
             rank_diff_file, tier_diff_file, primary_file)
  utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  profile_manifest_refs <- c(
    result19$output_files[grepl("_19_interpreter_manifest\\.csv$",
                                result19$output_files)],
    result36$output_files[grepl("_36_interpreter_manifest\\.csv$",
                                result36$output_files)]
  )
  profile_manifest_refs <- unique(profile_manifest_refs[file.exists(profile_manifest_refs)])
  manifest <- bind_rows(
    if (!is.na(analysis_script_file) && file.exists(analysis_script_file)) {
      file_inventory(analysis_script_file, "interpreter_script", utc)
    } else tibble(),
    file_inventory(c(p19$mean_vrs, p36$mean_vrs, p19$tier, p36$tier,
                     p19$primary_idx, p36$primary_idx, p19$codom, p36$codom),
                   "input", utc),
    file_inventory(profile_manifest_refs, "profile_interpreter_manifest", utc),
    file_inventory(files, "output", utc)
  )
  manifest_file <- file.path(machine,
                             paste0(cfg$species_code,
                                    "_19_vs_36_interpreter_manifest.csv"))
  write.csv(manifest, manifest_file, row.names = FALSE)

  cat("\n>>> 19-versus-36 sensitivity report written to: ", out_root, "\n", sep = "")
  invisible(list(summary = summary, files = c(files, manifest_file)))
}

# =============================================================================
# 12. EXECUTION
# =============================================================================

profiles <- unique(as.integer(cfg$profiles_to_run))
if (!length(profiles) || any(!profiles %in% c(19L, 36L))) {
  stop("profiles_to_run must contain 19L, 36L, or both.")
}

results <- list()
for (profile in profiles) {
  key <- as.character(profile)
  base_dir <- cfg$profile_dirs[[key]]
  if (is.null(base_dir) || !nzchar(base_dir)) stop("No profile directory for ", key)
  results[[key]] <- run_profile_interpreter(key, base_dir)
  gc(verbose = FALSE)
}

if (all(c("19", "36") %in% names(results))) {
  comparison <- build_profile_comparison(results[["19"]], results[["36"]])
} else {
  message("Paired sensitivity report skipped because both profiles were not run.")
}

cat("\n========================================================================\n")
cat("VERA SPECIES INTERPRETER COMPLETED\n")
cat("Species: ", cfg$species_label, " (", cfg$species_code, ")\n", sep = "")
cat("Profiles: ", paste(names(results), collapse = ", "), "\n", sep = "")
cat("========================================================================\n")

invisible(results)