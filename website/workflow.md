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
        +--> Mahalanobis tier figure
```

## Environmental predictors

`01_envirem_generation.R` prepares seventeen ENVIREM predictors from monthly
climate rasters. All layers must share one CRS, resolution, extent, origin and
cell alignment.

## Occurrence preparation

`02_occurrence_preparation.R` removes points outside the raster domain and
retains at most one occurrence per analysis-grid cell. A documented geographic
isolation stage then prepares the final tutorial occurrence table.

## VERA analysis

`04_vera_19_36.R` runs two independent analyses. Each profile has its own
complete-case filtering, asymmetric calibration, climatic-core partition,
Mahalanobis regularisation, empirical tiers and add-on diagnostics.

## Rendering

The renderer scripts consume completed VERA outputs without changing the core
rasters or consolidated diagnostic tables.
