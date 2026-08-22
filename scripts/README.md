# VERA tutorial scripts

These scripts form the executable *Sitta krueperi* tutorial and post-analysis workflow.

| Order | Script | Purpose |
|---:|---|---|
| 1 | `01_envirem_generation.R` | Generate the seventeen retained ENVIREM predictor layers. |
| 2 | `02_occurrence_preparation.R` | Clean and spatially thin occurrence records. |
| 3 | `03_vera_19.R` | Run the complete 19-predictor VERA profile. |
| 4 | `04_vera_19_36.R` | Run independent 19- and 36-predictor VERA profiles. |
| 5 | `05_render_core_outputs.R` | Render core VERA maps and plots. |
| 6 | `06_render_addons.R` | Render directional and cross-geometry add-ons. |
| 7 | `07_vera_species_interpreter.R` | Generate evidence-bound species summaries, diagnostic alerts, predictor tables, interpretation rasters, occurrence-review outputs and the 19-versus-36 profile-sensitivity report. |

Script 7 is a post-analysis interpretation layer. It reads completed VERA outputs but does not modify any canonical VRS, Mahalanobis, tier, attribution or add-on product.

The scripts retain explicit local path settings. Edit those paths before running
the tutorial on another computer.

Do not silently modify released scripts. Record executable changes with a new
version, checksum and changelog entry.
