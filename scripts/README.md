# VERA tutorial scripts

These scripts form the executable *Sitta krueperi* tutorial workflow.

| Order | Script | Purpose |
|---:|---|---|
| 1 | `01_envirem_generation.R` | Generate the seventeen ENVIREM predictor layers. |
| 2 | `02_occurrence_preparation.R` | Clean and spatially thin occurrence records. |
| 3 | `03_vera_19.R` | Run the complete 19-predictor VERA profile. |
| 4 | `04_vera_19_36.R` | Run independent 19- and 36-predictor VERA profiles. |
| 5 | `05_render_core_outputs.R` | Render core VERA maps and plots. |
| 6 | `06_render_addons.R` | Render directional and cross-geometry add-ons. |
| 7 | `07_render_mahalanobis_tiers.R` | Render Mahalanobis core and tier calibration. |

The scripts retain explicit local path settings. Edit those paths before running
the tutorial on another computer.

Do not silently modify released scripts. Record executable changes with a new
version, checksum and changelog entry.
