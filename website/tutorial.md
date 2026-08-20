---
title: "Sitta krueperi tutorial"
---

## Before running

Confirm that coordinate columns are recognised, predictor filenames match the
configuration, all rasters share one geometry, output folders are writable and
source-data licences permit the intended use.

## Step 1 — Generate ENVIREM predictors

Update the paths and run:

```r
source("01_envirem_generation.R")
```

Confirm that the naming check passes and seventeen GeoTIFFs are written.

## Step 2 — Prepare occurrences

```r
source("02_occurrence_preparation.R")
```

Record the initial, domain-cleaned, pixel-thinned and final occurrence counts.

## Step 3 — Run both profiles

```r
source("04_vera_19_36.R")
```

The 19-predictor profile finishes before the independent 36-predictor profile
begins. Each writes its own rasters, CSV tables, configuration, session
information, code manifest and output inventory.

## Step 4 — Render outputs

```r
source("05_render_core_outputs.R")
source("06_render_addons.R")
source("07_render_mahalanobis_tiers.R")
```

Inspect renderer logs for missing inputs and preserve the inventories with the
final figures.
