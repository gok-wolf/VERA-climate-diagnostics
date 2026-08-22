# VERA Species Interpreter

`07_vera_species_interpreter.R` is a deterministic post-analysis layer for a
completed VERA run. It reads exported Paper 1 outputs and writes bounded,
traceable summaries for one focal taxon. It does not modify any canonical VRS,
Mahalanobis, tier, attribution or add-on product.

## Requirements

Required R packages:

```r
terra
dplyr
tibble
```

Optional:

```r
openxlsx
```

If `openxlsx` is installed, the Predictor Evidence Matrix is written as both
CSV and formatted XLSX. Otherwise, the canonical CSV is still produced and the
manifest records `openxlsx_not_installed`.

## Configuration

Edit only the configuration block near the top of the script:

```r
cfg$species_code
cfg$species_label
cfg$profile_dirs
cfg$comparison_output_dir
```

The default Sitta paths are:

```text
C:/VERA/Results/19/Skr_current
C:/VERA/Results/36/Skr_current
```

Run:

```r
source("scripts/07_vera_species_interpreter.R")
```

The same script supports all three execution modes:

```r
profiles_to_run = 19L
profiles_to_run = 36L
profiles_to_run = c(19L, 36L)
```

The paired sensitivity report is created only when both profiles are selected.

## Profile-level outputs

Each completed VERA profile receives an `interpretation/` directory:

```text
Skr_current/
└── interpretation/
    ├── human_readable/
    │   ├── Skr_19_species_diagnostic_brief.md
    │   └── Skr_19_diagnostic_alerts.md
    ├── tables/
    │   ├── Skr_19_predictor_evidence_matrix.csv
    │   ├── Skr_19_predictor_evidence_matrix.xlsx   [optional]
    │   ├── Skr_19_diagnostic_alerts.csv
    │   ├── Skr_19_occurrence_review_queue.csv
    │   └── Skr_19_interpretation_status_counts.csv
    ├── rasters/
    │   └── Skr_19_interpretation_status.tif
    └── machine_readable/
        ├── Skr_19_narrative_evidence.csv
        ├── Skr_19_interpreter_manifest.csv
        └── Skr_19_interpreter_sessionInfo.txt
```

The 36-predictor profile receives the same files with `_36_` in their names.

## Spatial Interpretation Status

This reporting-only raster uses existing VERA tiers and diagnostic flags. It
does not introduce a new ecological model. Its classes are:

1. Lower empirical departure (tiers 1–2)
2. Elevated departure with a strict single maximum
3. Elevated departure with a co-dominant maximum
4. VRS cap engaged; values are censored at the configured ceiling
5. Reduced valid-predictor coverage

Precedence is deterministic: reduced coverage overrides cap engagement; cap
engagement overrides co-dominance; co-dominance overrides the strict-maximum
class. The accompanying CSV records the definition and pixel share of every
class. Class 4 is assigned wherever maximum VRS reaches the configured cap,
regardless of the empirical climatic tier; it records numerical censoring, not
an ecological threshold.

## Occurrence Review Queue

The queue retains every exported occurrence and sorts records by climatic
Mahalanobis distance when that field is available. Raster values are extracted
at the occurrence coordinates. The queue supports provenance and coordinate
review; it does not recommend automatic record deletion.

## Paired profile sensitivity

When both profiles are present, the script creates:

```text
C:/VERA/Results/Skr_profile_sensitivity/
```

Products include:

- a root-level `README.txt` explaining the package structure and safe reading
  rules;
- mean-VRS percentile-rank difference raster;
- tier-difference raster and tier cross-tabulation;
- tie-broken Primary Stressor transition raster and table;
- a second Primary transition table restricted to pixels with a strict single
  maximum in both profiles;
- profile-sensitivity summary CSV;
- evidence-bound Markdown report;
- narrative evidence log and manifest.

Raw mean-VRS magnitudes are not treated as interchangeable across different
predictor counts. Spatial comparison therefore emphasises percentile ranks,
tier agreement and named attribution transitions.

## Language safety

All generated evidence sentences pass through a forbidden-claim linter. Each
sentence in the Species Diagnostic Brief and profile-sensitivity report is
linked to its source file, section, metric, raw value and template identifier.

The interpreter describes climatic departure, diagnostic resolution, data
support and predictor-profile sensitivity. It does not infer occurrence
probability, habitat suitability, physiological performance, demography,
dispersal, causal range limitation or conservation priority.