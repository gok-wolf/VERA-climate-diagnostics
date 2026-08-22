---
title: "Variable Ecological Restriction Analysis"
subtitle: "Magnitude, identity, direction and diagnostic resolution of climatic departure"
description: "VERA documentation, illustrated tutorial, evidence-bound interpretation tools and optional publication renderers."
image: figures/vera-brand-banner.png
toc: false
page-layout: full
title-block-banner: true
---

## What is VERA?

**Variable Ecological Restriction Analysis (VERA)** is an occurrence-calibrated environmental diagnostic framework. It maps how strongly local climate departs from a species' occupied reference, identifies the predictor carrying the leading signal, records whether that signal is uniquely resolved or shared, and distinguishes high-side from low-side departures.

The primary diagnostic layer is formed by directionally scaled **Variable Restriction Scores (VRS)**. The canonical tutorial also computes a covariance-aware **Mahalanobis companion geometry**, which provides multivariate climatic context and supports the add-on percentile-rank comparison. The two geometries remain separate and are not collapsed into one composite score.

> **Interpretive boundary:** VERA does not estimate occurrence probability, habitat suitability, physiological tolerance, demographic performance or causal range limitation.

![VERA — Variable Ecological Restriction Analysis.](figures/vera-brand-banner.png){fig-alt="VERA project banner showing the full name Variable Ecological Restriction Analysis." fig-align="center" width="70%"}

## Tutorial case study

This site demonstrates the workflow using *Sitta krueperi* (Krüper's nuthatch) under independent **19-predictor** and **36-predictor** profiles.

The worked example covers occurrence preparation, ENVIREM generation, asymmetric calibration, core VERA surfaces, tie-aware attribution, directional add-ons, VRS–Mahalanobis rank comparison, evidence-bound species interpretation and publication-oriented graphics.

## What the repository provides

### Core workflow

- predictor generation and documented occurrence preparation;
- independent 19- and 36-predictor VERA analyses;
- Mean VRS, Maximum VRS, VPI and categorical attribution surfaces;
- directional and cross-geometry add-on diagnostics;
- reproducibility manifests, checksums and output inventories.

### Post-analysis interpretation

`07_vera_species_interpreter.R` generates species diagnostic briefs, predictor evidence tables, diagnostic alerts, interpretation-status rasters, occurrence-review queues and narrative evidence logs. It can run a 19-predictor profile, a 36-predictor profile, or both profiles with a paired sensitivity report.

### Optional publication renderers

The optional scripts provide rapid graphics without changing canonical VERA outputs:

- `08_render_mahalanobis_tiers.R` — occurrence-partition and empirical Mahalanobis tier-calibration figures;
- `09_render_response_curve_panels.R` — Top-6 native-unit density, asymmetric transformation and anchor-annotation products;
- `10_render_top10_response_summaries.R` — Top-10 density galleries and aligned-optimum ridge summaries.

The dedicated Mahalanobis tier-calibration figure is optional. The Mahalanobis–VRS rank-comparison diagnostic remains available in the add-on output family.

## Start here

1. Read the [workflow](workflow.md).
2. Follow the [illustrated *Sitta* tutorial](tutorial.md).
3. Use the [minimum reporting guide](tutorial.md#minimum-reporting-set) to select essential manuscript maps.
4. Consult the [output guide](outputs.md).
5. Review the [reproducibility requirements](reproducibility.md).