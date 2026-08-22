---
title: "Workflow"
---

## Overview

```text
Monthly climate data
        +--> ENVIREM generation

Occurrence records
        +--> Domain cleaning
        +--> Pixel-level thinning
        +--> Documented geographic isolation

        +--> 19-predictor VERA
        +--> 36-predictor VERA

        +--> Core renderer
        +--> Add-on renderer
        +--> Species Interpreter
```

## Environmental predictors

`01_envirem_generation.R` prepares the seventeen retained ENVIREM predictors
from monthly climate rasters. All layers must share one CRS, resolution,
extent, origin and cell alignment. The bounded discrete variable
`monthCountByTemp10` is not included in the 36-predictor tutorial profile.

## Occurrence preparation

`02_occurrence_preparation.R` removes points outside the raster domain and
retains at most one occurrence per analysis-grid cell. A documented geographic
isolation stage then prepares the final tutorial occurrence table. These steps
define the empirical occupied reference and should be reported explicitly.

## VERA analysis

`04_vera_19_36.R` runs independent 19- and 36-predictor analyses. Each profile
has its own complete-case filtering, asymmetric calibration, climatic-core
partition, covariance-aware Mahalanobis companion geometry and add-on
diagnostics. Results from the two predictor profiles are not combined into one
composite surface.

## Rendering and interpretation

`05_render_core_outputs.R` and `06_render_addons.R` translate completed VERA
outputs into the core and add-on figure families without changing the
underlying rasters or consolidated diagnostic tables.

`08_vera_species_interpreter.R` provides a separate evidence-bound
post-analysis layer. It produces species briefs, alerts, evidence tables,
interpretation-status rasters and occurrence-review products. When both
profiles are selected, it also produces the 19-versus-36 sensitivity report.

## Optional publication graphics

The following scripts are optional and do not alter canonical VERA outputs:

- `07_render_mahalanobis_tiers.R` — dedicated occurrence-partition and
  Mahalanobis tier-calibration figures;
- `09_render_response_curve_panels.R` — Top-6 native-unit density and
  asymmetric-transformation panels;
- `10_render_top10_response_summaries.R` — Top-10 density and aligned-optimum
  ridge galleries.

These graphics may be supplied according to the research question,
Supplementary Information plan and journal requirements.

