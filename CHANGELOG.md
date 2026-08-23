# Changelog

All notable repository changes will be documented here.

## 0.2.1-dev — 2026-08-23

- Replaced historical `v10`/`pilot` reporting labels with neutral
  reporting-audit terminology throughout the canonical 19- and 36-predictor
  scripts, regional batch pipelines, configurations and Supporting Information.
- Renamed the reporting switches to `write_reporting_audit` and
  `reporting_contribution_coverage` without changing any VRS, Mahalanobis,
  attribution, tier or add-on calculation.
- Renamed reporting-only helper functions, intermediate objects and exported
  table-section labels to match the neutral terminology.
- Updated canonical script SHA-256 checksums after the naming-only refactor.

## 0.2.0-dev — 2026-08-22

- Added the evidence-bound VERA Species Interpreter with 19-only, 36-only and
  paired-profile execution modes.
- Reclassified the Mahalanobis tier-calibration figure as an optional renderer
  and extended it to both tutorial profiles.
- Added optional response-curve, Top-10 density and aligned-optimum ridge
  renderers for rapid Sitta tutorial visualisation.
- Renumbered the Species Interpreter as Script 07 and the optional Mahalanobis
  tier-calibration renderer as Script 08.
- Standardised the new graphical outputs and optional Excel table to Arial.
- Updated script checksums, package requirements and the script catalogue.
- Deferred formal Citation File Format metadata until authorship, release and
  DOI details are finalised.
- Added an automatically generated, manifest-tracked README to the paired
  profile-sensitivity output package.
- Added a profile-specific `README_08.txt` guide to each optional Mahalanobis
  tier-calibration figure directory.
- Added profile-specific `README_09.txt` and `README_10.txt` guides to the
  optional response-panel and response-summary directories.
- Corrected Scripts 09 and 10 to resolve the Script 05 Primary-assignment count
  table from the profile-level `Images/core_renderer_plum/` directory, with a
  legacy nested-path fallback.

## 0.1.0-dev — 2026-08-20

- Created the initial repository architecture.
- Added the *Sitta krueperi* data-preparation, dual-profile VERA and renderer
  scripts without modifying their executable contents.
- Added provisional citation, data-governance and reproducibility documents.
- Added the initial Quarto website source structure.
