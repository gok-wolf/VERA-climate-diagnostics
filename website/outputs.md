---
title: "Output guide"
---

## Read outputs jointly

No single VERA layer provides a complete ecological interpretation. Magnitude,
identity, directional attribution, climatic proximity and covariance-aware
context answer different questions and should be read together.

## Mean and Maximum VRS

Mean VRS summarises average predictor-wise climatic departure. Maximum VRS
reports the strongest single predictor-specific departure. Their difference
helps distinguish broadly distributed departure from a strong leading axis.
Neither surface establishes ecological causality.

## Delta VRS

Delta VRS records separation between the two leading predictor scores. Small
values indicate weakly resolved dominance; zero indicates a shared maximum.
Delta VRS is a diagnostic-resolution measure rather than a suitability index.

## VPI

The Variable Proximity Index is a bounded inverse transformation of Mean VRS.
Higher VPI values indicate greater climatic proximity to the occupied
reference. VPI is not habitat suitability, habitat quality or occurrence
probability.

## Primary and Unique Primary Stressors

Primary Stressor provides a deterministic leading-predictor label. Unique
Primary retains the identity only where one predictor has a strict maximum.
Unique therefore means resolved diagnostic dominance, not suitable habitat or
causal limitation.

## Secondary, TooHigh and TooLow Stressors

Secondary Stressor records the predictor associated with the second-largest
score. TooHigh and TooLow identify the leading departure signals above and
below predictor-specific occupied reference centres. Together they describe
the local hierarchy and direction of climatic departure without establishing
causation.

## Directional and cross-geometry add-ons

Add-On 1 compares **Dominant Tail Direction** with **Net Tail Direction**. The
first reports the side of the strongest individual predictor signal; the
second reports the side of the cumulative predictor-wise contribution.

Add-On 2 compares landscape percentile rankings from Mean VRS and the
covariance-aware Mahalanobis companion geometry. Its agreement classes localise
concordant rankings and areas where one geometry ranks pixels more strongly
than the other. Agreement is not independent validation because both
diagnostics use the same occurrences and predictor data.

## Interpretation and optional products

`07_vera_species_interpreter.R` converts exported metrics into evidence-bound
briefs, alerts, tables, review queues and interpretation-status rasters. Every
generated narrative sentence is linked to its source metric and template.

Dedicated Mahalanobis tier-calibration figures, response-curve panels, density
galleries and aligned-optimum ridges are optional publication products. They
can support explanation or Supplementary Information but are not required to
interpret the core VERA output family.
