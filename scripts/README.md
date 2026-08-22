# VERA tutorial scripts

These scripts form the executable *Sitta krueperi* tutorial and post-analysis
workflow.

| Order | Script | Status | Purpose |
|---:|---|---|---|
| 1 | `01_envirem_generation.R` | Preprocessing | Generate the seventeen retained ENVIREM predictor layers. |
| 2 | `02_occurrence_preparation.R` | Preprocessing | Clean and spatially thin occurrence records. |
| 3 | `03_vera_19.R` | Canonical profile | Run the complete 19-predictor VERA profile. |
| 4 | `04_vera_19_36.R` | Primary tutorial | Run independent 19- and 36-predictor VERA profiles. |
| 5 | `05_render_core_outputs.R` | Core renderer | Render core VERA maps and plots. |
| 6 | `06_render_addons.R` | Add-on renderer | Render directional and cross-geometry add-ons. |
| 7 | `07_vera_species_interpreter.R` | Post-analysis interpreter | Generate evidence-bound species summaries, diagnostic alerts, predictor tables, interpretation rasters, occurrence-review outputs and, when both profiles are selected, the 19-versus-36 profile-sensitivity report. |
| 8 | `08_render_mahalanobis_tiers.R` | Optional renderer | Render occurrence-partition and empirical Mahalanobis tier-calibration figures for the 19- and/or 36-predictor profile. |
| 9 | `09_render_response_curve_panels.R` | Optional renderer | Export Top-6 anchor annotations and rapidly render native-unit density and asymmetric Z2 panels for the 19- and/or 36-predictor profile. |
| 10 | `10_render_top10_response_summaries.R` | Optional renderer | Render Top-10 native-unit density galleries and aligned-optimum ridge summaries for the 19- and/or 36-predictor profile. |

## Species Interpreter modes

Script 7 is one shared implementation; separate copies are not required. Set
`profiles_to_run` in its configuration block to choose the desired mode:

```r
profiles_to_run = 19L          # 19-predictor profile only
profiles_to_run = 36L          # 36-predictor profile only
profiles_to_run = c(19L, 36L)  # both profiles plus paired sensitivity
```

The interpreter reads completed VERA outputs but does not modify any canonical
VRS, Mahalanobis, tier, attribution or add-on product.

## Optional figure scripts

Scripts 8, 9 and 10 are convenience renderers. They are not required to run
VERA and do not define new diagnostic calculations.

- Script 8 preserves the earlier Mahalanobis calibration story as an optional
  two-panel figure family and writes a profile-specific `README_08.txt` beside
  the exported figures.
- Script 9 produces individual response-curve panels and a machine-readable
  annotation table for the most frequent Primary Stressors, and writes
  `README_09.txt` under `Images/optional_response_curves/`.
- Script 10 produces compact gallery figures for readers who want a rapid
  overview of native-unit and standardized predictor profiles, and writes
  `README_10.txt` under `Images/optional_response_summaries/`.

All new optional graphics use Arial throughout. Each optional renderer accepts
`19L`, `36L`, or both profiles through its `profiles_to_run` configuration.
Scripts 9 and 10 should be run after Script 5 because they use
`vera_pixel_counts_primary_assigned.csv` to identify the leading predictors.

```r
source("scripts/07_vera_species_interpreter.R")
source("scripts/08_render_mahalanobis_tiers.R")
source("scripts/09_render_response_curve_panels.R")
source("scripts/10_render_top10_response_summaries.R")
```

## Editing rule

The scripts retain explicit local path settings. Edit those paths before
running the tutorial on another computer.

Do not silently modify released scripts. Record executable changes with a new
version, checksum and changelog entry, and regenerate every affected output.
