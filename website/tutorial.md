# VERA — A Complete Tutorial Using Krüper's Nuthatch

**Variable Ecological Restriction Analysis (VERA)** is an asymmetric climatic-niche diagnostic framework. Unlike traditional species distribution models (SDMs), VERA does not attempt to predict *where a species is likely to occur*. It attempts to answer a different, sharper question: **for every pixel of the landscape, why is the environment restrictive — which specific variable is limiting the species, from which direction (too high or too low), and how does that restriction interact with the covariance-aware multivariate niche geometry?**

This repository contains the full public tutorial pipeline: from raw climate rasters through cleaned occurrences, through the VERA core computation, all the way to publication-grade diagnostic figures. Every script here is the exact code used to produce the figures shown below.

The focal taxon throughout is *Sitta krueperi* (Krüper's Nuthatch), an eastern-Mediterranean forest species with a distribution largely restricted to Türkiye.

---

## Table of Contents

- [Before Any Code Runs: Why Clean Occurrences Matter More in VERA Than in Any Other Framework](#before-any-code-runs)
- [Repository Layout](#repository-layout)
- [Stage 1 — Building the Predictor Stack (ENVIREM)](#stage-1)
- [Stage 2 — The Three-Step Occurrence Cascade](#stage-2)
- [Stage 3 — Running the VERA Core Pipeline](#stage-3)
- [Stage 4 — Rendering the Diagnostic Figures](#stage-4)
- [Reading the Outputs — Guided Figure Tour](#reading-the-outputs)
  - [Figure 1 — The VERA Conceptual Framework](#figure-1)
  - [Figure 2 — Per-Variable Niches in Native Units](#figure-2)
  - [Figure 3 — Aligned-Optimum Cross-Variable Comparison](#figure-3)
  - [Figure 4 — Core Spatial Outputs](#figure-4)
  - [Figure 5 — Mahalanobis Distance and Empirical Climatic Tiers](#figure-5)
  - [Figure 6 — Categorical Stressor Attribution](#figure-6)
  - [Figure 7 — Add-On Branch A: Metric-Agreement Diagnostic](#figure-7)
  - [Figure 8 — Add-On Branch B: Directional Attribution](#figure-8)
- [Which Script Produced Which Figure](#which-script-produced-which-figure)
- [Reproducibility](#reproducibility)

---

<a id="before-any-code-runs"></a>
## Before Any Code Runs: Why Clean Occurrences Matter More in VERA Than in Any Other Framework

If you take one thing away from this tutorial, take this: **VERA is unusually sensitive to occurrence quality, and understanding why is essential before running any of the code that follows.**

Most SDMs consume occurrences primarily as *presence signals* fed into a discriminative learner (MaxEnt, random forest, GLM, etc.). A handful of miscoded points, spatially-biased duplicates, or lonely geographic outliers usually add noise but rarely rewrite the model's conclusions — the learner regularises them away.

VERA works differently. Every downstream output — every diagnostic map, every tier boundary, every stressor attribution — is built on top of three per-variable **anchor statistics** derived directly from your occurrence values:

- **μ** (mu) — the occurrence-derived climatic optimum
- **σ_L** — the lower-tail scale (asymmetric distance to below-optimum values)
- **σ_U** — the upper-tail scale (asymmetric distance to above-optimum values)

There is no discriminative learner between your points and these numbers. A wildly wrong occurrence:

- shifts **μ** toward the outlier's value,
- inflates the σ on whichever tail the outlier landed,
- warps the empirical VRS surface across the entire landscape,
- distorts the covariance matrix used for Mahalanobis calibration,
- and shifts the empirical 80 / 95 / 99 % percentile thresholds that define the four climatic tiers.

There is also a subtler failure mode. If your occurrence set is spatially oversampled — many points clustered near cities, roads, or protected areas — **μ** silently shifts toward the environmental conditions of the *sampling effort* rather than the *actual centroid of use*, and the tail scales artificially tighten. The maps that come out of VERA will then look ecologically confident when they are really just tracking observer behaviour.

For these reasons, this tutorial dedicates an entire preprocessing stage — a three-step cascade — to producing a defensible occurrence table before VERA is ever invoked. The stages are not optional and they are not interchangeable; each one catches a different failure mode.

> **The rule:** VERA is only as trustworthy as the occurrence table you feed it. Preprocessing is not a cosmetic step. It is the foundation on which every downstream figure is standing.

---

<a id="repository-layout"></a>
## Repository Layout

```
scripts/
├── 01_envirem_generation.R          # Build 17 ENVIREM predictors from monthly climate
├── 02_occurrence_preparation.R      # Three-step occurrence cleaning cascade
├── 03_vera_19.R                     # Convenience: 19-predictor VERA run only
├── 04_vera_19_36.R                  # Canonical: dual 19- and 36-predictor VERA
├── 05_render_core_outputs.R         # CORE renderer — spatial and per-variable diagnostics
├── 06_render_addons.R               # ADD-ON renderer — agreement + tail-direction diagnostics
└── 07_render_mahalanobis_tiers.R    # Figure 4 renderer — Mahalanobis tier calibration
```

The scripts are numbered to reflect execution order. `04_vera_19_36.R` is the canonical implementation; `03_vera_19.R` is a convenience version for users who only need the classical 19 bioclimatic predictors. Both write outputs into `C:/VERA/Results/{19,36}/Skr_current/` using the exact same directory structure, so downstream renderers work identically on either result.

The renderer scripts (`05`, `06`, `07`) consume completed VERA outputs and do not modify any core diagnostic. They can be re-run at any time without recomputing the analysis.

---

<a id="stage-1"></a>
## Stage 1 — Building the Predictor Stack (ENVIREM)

**Script:** `01_envirem_generation.R`

VERA can run on the classical 19 WorldClim / CHELSA bioclimatic variables alone, but the canonical Sitta krueperi analysis uses a 36-predictor stack: the 19 bioclimatics plus 17 ENVIREM extended variables (aridity indices, growing degree days, seasonal PET terms, thermicity, continentality, etc.).

This first script builds the 17 ENVIREM rasters from monthly climate inputs (`prec_XX`, `tmin_XX`, `tmax_XX`, `tavg_XX` for months 01–12) using the `envirem` R package. Solar radiation is derived internally for year 2000. Every ENVIREM raster is written as a separate GeoTIFF and copied into the same folder as the 19 bioclimatics, so both predictor sets share identical geometry, resolution and CRS.

> **Requirement.** All monthly rasters must have matching geometry, resolution, CRS and extent. VERA's raster validator will reject a stack with mismatched geometries.

---

<a id="stage-2"></a>
## Stage 2 — The Three-Step Occurrence Cascade

**Script:** `02_occurrence_preparation.R`

This is where the "clean occurrences" principle from the preface is operationalised. The script implements three independent cleaning stages, in this exact order. Each stage catches a distinct failure mode. Skipping any stage leaves that failure mode intact in the downstream anchor statistics.

### Stage 2.1 — Domain cleaning

Points are converted to a `SpatVector`, projected to the predictor stack's CRS, and used to extract Bio1 values. Any point whose extraction returns `NA` (i.e., falls outside the valid raster footprint or into a masked cell) is dropped.

**What this catches.** Points in the sea, on nodata cells, or beyond the raster envelope. These would silently drop out of the anchor computation as incomplete rows anyway, but removing them here makes the point count honest from the very first step and preserves the traceability of every dropped record.

### Stage 2.2 — Pixel-level thinning

For each retained occurrence, the raster cell number is computed via `cellFromXY`. Records that share a cell are deduplicated so that at most **one occurrence per raster cell** survives.

**What this catches — and why it is critical for VERA specifically.** Pixel-level thinning removes the *sampling-effort bias* that quietly distorts anchor statistics. If a single 2.5-arc-minute cell contains twelve overlapping records — three from an eBird hotspot, four from a museum accession, five from a protected-area survey — the anchor μ will be pulled twelve times toward that cell's climate. The tail scales σ_L, σ_U will tighten around it. The species will appear artificially specialised in whatever climate that one cell happens to have. Pixel thinning enforces "one climate condition, one vote": μ represents the environmental centroid of *use*, not the centroid of *observation effort*.

Unlike distance-based spatial thinning (which forces a minimum separation between points), pixel-level thinning is exactly matched to the raster resolution being modelled — the finest granularity at which VERA can distinguish environmental conditions. It is the correct thinning granularity when the predictor pixel size is already the effective unit of analysis.

### Stage 2.3 — Geographic outlier isolation (retain 97 %)

The thinned points are projected to a metric CRS (EPSG:3857), a geographic Mahalanobis distance is computed from the point cloud's spatial centroid and covariance, and the top 3 % most-distant points are removed. The remaining 97 % becomes the final occurrence table used by VERA.

**What this catches.** Coordinate errors that survived domain cleaning (a point plotting into Türkiye's ecological envelope but genuinely misidentified), specimens with imprecise georeferences that landed far outside the species' actual range, and long-distance vagrants that do not represent stable ecological use. VERA's anchor computation is unbiased against typical range-edge variation but *is* sensitive to extreme geographic outliers because they usually come attached to extreme climatic outliers.

Why 3 %? It is a deliberately gentle threshold — aggressive enough to remove obvious geographic misfits, conservative enough to preserve legitimate peripheral populations at the range's true climatic and geographic edge. The peripheral populations you *want* to keep (climatic marginal but geographically plausible) will be caught downstream by VERA's own core / peripheral partition (see [Figure 4](#figure-4)) using the covariance-aware climatic Mahalanobis distance, which is the right tool for that job. Stage 2.3 handles only the geographic pathology.

> **The full cascade.** A raw table of ~4,000 records typically drops to ~2,800 after domain cleaning, to ~1,100 after pixel thinning, and to ~1,067 after the 3 % geographic isolation. This is not data loss — it is data hygiene. Every one of those removed records was a threat to at least one of your downstream anchor statistics.

The final table lands at `C:/VERA/Occurrences/Edited/Sitta_krueperi.csv` and is the *only* occurrence input consumed by the VERA scripts that follow.

---

<a id="stage-3"></a>
## Stage 3 — Running the VERA Core Pipeline

**Scripts:** `03_vera_19.R` (convenience, 19 predictors) and `04_vera_19_36.R` (canonical, 19 + 36 predictors in the same run).

Most users should start with `04_vera_19_36.R`. It runs the same focal taxon under both the 19-predictor and 36-predictor profiles in a single pass, writing outputs into two parallel result trees:

```
C:/VERA/Results/19/Skr_current/
C:/VERA/Results/36/Skr_current/
```

Both trees have identical internal structure, so any renderer script works on either without modification. `03_vera_19.R` exists purely for users whose research question requires only the 19 classical bioclimatics — it is a strict subset of what `04_vera_19_36.R` does.

Internally, each VERA run performs the following steps for the focal taxon:

1. Reads the cleaned occurrence table and extracts predictor values at each point.
2. For every predictor, estimates the asymmetric anchor triple **(μ, σ_L, σ_U)** with shrinkage of small-sample tail scales toward the global occurrence SD, a bootstrap-stability check on the asymmetry ratio σ_U / σ_L, and a background-truncation diagnostic.
3. Partitions occurrences into a **climatic core** and a **peripheral** set using the covariance-aware Mahalanobis distance and a 0.90 quantile rule.
4. Computes the per-pixel VRS surfaces: **mean VRS**, **max VRS**, **Δ VRS**, primary / secondary / TooHigh / TooLow stressor identity rasters, co-dominance counts.
5. Computes the covariance-aware **Mahalanobis distance** raster and derives empirical **80 / 95 / 99 % tier thresholds** from the core occurrences.
6. Runs Add-On 01 (tail direction) and Add-On 02 (VRS vs Mahalanobis rank agreement).
7. Writes a machine-readable manifest, session info, code MD5 checksum, and full output inventory.

The result is a single directory tree containing everything the renderer scripts need. No further computation happens after Stage 3.

---

<a id="stage-4"></a>
## Stage 4 — Rendering the Diagnostic Figures

**Scripts:** `05_render_core_outputs.R`, `06_render_addons.R`, `07_render_mahalanobis_tiers.R`.

The renderers consume the finished VERA outputs and produce every published figure in this tutorial. They do not perform any diagnostic calculation of their own — they only translate finished CSVs and GeoTIFFs into publication-grade PNGs.

Each renderer loops internally through both the 19- and 36-predictor result trees, so a single execution populates two parallel image directories:

```
C:/VERA/Results/19/Images/{core_renderer_plum, addon_renderer_plum, figure4_tier_calibration}/
C:/VERA/Results/36/Images/{core_renderer_plum, addon_renderer_plum, figure4_tier_calibration}/
```

Missing result trees are skipped silently — running the renderers on a 19-only setup produces only the 19-profile images.

---

<a id="reading-the-outputs"></a>
## Reading the Outputs — Guided Figure Tour

The eight figures below are arranged in the order that most efficiently teaches VERA. Figure 1 introduces the framework's asymmetric geometry conceptually; Figures 2–3 show that geometry realised across the top-10 stressors; Figures 4–5 move into space, showing where and how much the landscape restricts the species; Figures 6–8 attribute those restrictions to specific predictors, specific directions, and cross-check them against a covariance-aware baseline.

---

<a id="figure-1"></a>
### Figure 1 — The VERA Conceptual Framework

![Figure 1](Figure_1.png)

This two-panel figure is the conceptual foundation for everything else. Both panels use a single predictor — **PET of the Wettest Quarter (PETWeQ)** — to make the asymmetric restriction geometry visible in native units.

**Panel (a) — The observable niche.** The blue density is the species' realised distribution along PETWeQ (the *Occupied Niche*). The yellow density is the region-wide availability of PETWeQ values (the *Available Climate Space*). The narrow blue peak sits on the low end of the gradient; the yellow density has a broad secondary peak far to the right where PETWeQ values reach 250–350 mm/month.

That right-hand yellow peak is the *Unoccupied Potential*: climate conditions that physically exist across Türkiye but are *not* used by Sitta krueperi. Because the species had physical access to those conditions and still declined to occupy them, the restriction must be intrinsic — physiological or ecological — rather than a lack of available space.

The two shaded zones anchor the asymmetric tolerance interpretation. The **red Constrained Tolerance Zone** on the left of the optimum (spanning μ − σ_L to μ) marks the narrow buffer within which the species can tolerate downward departures. The **green Extended Tolerance Zone** on the right (μ to μ + σ_U) marks the far wider buffer within which it tolerates upward departures.

**Panel (b) — The restriction gradient.** Panel (b) transforms Panel (a) into the machinery VERA actually computes. The x-axis is the same (PETWeQ in mm/month) but the y-axis is now the per-pixel **Z²** value — the squared, asymmetric, capped restriction score contributed by this variable alone.

The curve is deliberately not symmetric. On the left side of μ, Z² rises steeply as PETWeQ decreases — a downward departure of a few tens of mm/month drives Z² toward the cap almost immediately. This is **rapid stress accumulation**, the mathematical face of low lower-tail tolerance. On the right side, Z² rises far more gently — the species tolerates an increase of hundreds of mm/month before Z² approaches the cap. This is the **high ecological buffer** on the upper tail.

The two horizontal reference lines — **Z² = 4** (moderate) and the **Z² = 16 cap** — are VERA's fixed thresholds. Any predictor whose Z² reaches 16 at a pixel is contributing the maximum possible restriction score and is called **saturated** for that pixel.

**The reported asymmetry metrics** (μ = 112.33, σ_L = 32.13, σ_U = 103.48) give an asymmetry ratio of **3.31** with a 95 % bootstrap CI of 3.08–3.55 and the *bootstrap-stable* flag = Yes. That flag matters: it certifies that the extreme asymmetry is not a small-sample artefact but a robust, reproducible feature of the species' response.

> **What Figure 1 teaches you.** Every subsequent figure in this tutorial is built out of this asymmetric geometry, one predictor at a time, then summed, projected across the landscape, and cross-attributed. Understand Panel (a) and Panel (b), and every downstream figure becomes readable.

---

<a id="figure-2"></a>
### Figure 2 — Per-Variable Niches in Native Units

![Figure 2](Figure_2.png)

Figure 2 replicates Panel (a) of Figure 1 for the **top 10 primary stressors** — the ten predictors that most frequently rank as the single most restrictive variable across Türkiye. Each panel keeps the variable's *native units* on the x-axis (mm, °C, mm/month, index units), so the ecological thresholds are directly readable.

For every panel:

- **Blue** density = species' occupied climate (occurrences),
- **Yellow** density = available climate space (background),
- **Red-shaded band** = μ − σ_L to μ (Constrained Tolerance Zone),
- **Green-shaded band** = μ to μ + σ_U (Extended Tolerance Zone),
- **Red vertical line** = climatic optimum μ.

Reading this figure means asking the same three questions of each panel. **First**, does the blue occupation match, avoid, or partially overlap the yellow availability? If the yellow extends far beyond the blue, the species is being *restricted* rather than *limited by opportunity*. **Second**, is the red band narrow and the green band wide (or vice versa)? Ratio and direction of asymmetry differ meaningfully across predictors. **Third**, does the blue density terminate sharply at one edge (a physiological wall) or taper smoothly (a gradient of decreasing suitability)?

Sitta krueperi shows several distinct patterns across the ten panels — sharp lower-edge specialisation on some predictors (Bio14, PETWeQ), broadly overlapping occupation with modest asymmetry on others (Bio3, Bio7), and a striking case where the blue occupation is pressed *against the upper edge* of the available range (AIT, Aridity Index), suggesting the species uses the aridity-tolerant end of what Türkiye's climate offers.

---

<a id="figure-3"></a>
### Figure 3 — Aligned-Optimum Cross-Variable Comparison

![Figure 3](Figure_3.png)

Figure 2 shows each predictor in its own units, which is faithful but makes cross-variable comparison hard: a 30-mm departure in Bio14 and a 100-unit departure in PETs are not directly comparable. Figure 3 fixes this by standardising every predictor against its own asymmetric anchor:

- Values below μ are rescaled by σ_L (so the lower-tolerance boundary always sits at z = −1),
- Values above μ are rescaled by σ_U (so the upper-tolerance boundary always sits at z = +1),
- The optimum itself always sits at z = 0.

This is the *same* standardisation VERA uses internally to compute Z². The x-axis is therefore not an arbitrary rescaling — it is the exact restriction-score coordinate.

With every μ aligned at zero, the **central red vertical line** is the shared optimum. The **red vertical band** [−1, 0] is the Constrained Tolerance Zone (identical across all variables in these standardised units). The **green vertical band** [0, +1] is the Extended Tolerance Zone.

Reading the ridges top-to-bottom, four ecological questions become answerable at a glance:

1. **Is the blue occupation narrower or wider than the yellow availability?** Narrow blue relative to yellow = specialisation. Wide blue similar to yellow = generalist.
2. **Where does the blue peak sit relative to the red optimum line?** Aligned with the line = optimum-centred use. Shifted left or right = tail-biased use even inside the species' tolerable range.
3. **Does the blue extend into the red band, into the green band, or beyond ±1 entirely?** The further out, the more the species is enduring measurable restriction on that variable.
4. **Does the yellow availability extend far beyond the blue occupation on one side?** That is the unoccupied potential from Figure 1, generalised across every variable.

Together Figures 2 and 3 form a matched pair: Figure 2 preserves ecological realism (native units, direct interpretability), Figure 3 enables commensurable comparison (standardised units, cross-variable pattern-reading).

---

<a id="figure-4"></a>
### Figure 4 — Core Spatial Outputs

![Figure 4](Figure_4.png)

Figure 4 is the first spatial output of the pipeline. It stacks five maps of Türkiye, each showing a different facet of the species' climatic geography.

**Panel 1 — Occurrence-Level Climatic Partition.** Every occurrence in the cleaned table is classified as **Core** (light pink) or **Peripheral** (dark purple) based on its multivariate Mahalanobis distance to the occupied climatic centroid. Points inside the 0.90 quantile are Core; the remaining 10 % are Peripheral. The critical observation is that Core and Peripheral are *not* geographically segregated. Peripheral occurrences (highlighted inside the two rectangular insets) sit intermingled with Core occurrences over short distances — sometimes within the same mountain range. This is a direct empirical demonstration that in complex topography, climatic distance and geographic distance are decoupled: microclimatic variation across a slope or an elevation gradient can create peripheral climatic conditions right next to core ones.

**Panel 2 — Mean VRS.** The pixel-wise mean of the Z² restriction contributions across all predictors. Low (light) values mean the species faces low average restriction; high (dark) values mean compounded systemic stress. Mean VRS answers *how heavy is the total restriction load* at each pixel.

**Panel 3 — Maximum VRS.** For each pixel, the largest single-predictor Z² contribution. Where Max VRS is dark, at least one predictor is severely restricting the species even if the others are tolerable. Max VRS answers *does this pixel contain a hard bottleneck?* Compared to Mean VRS, Max VRS is more sensitive to Liebig-style single-variable limitation.

**Panel 4 — VPI (Variable Proximity Index).** Defined as VPI = 1 / (1 + mean VRS). Ranges from 0 (impassable) to 1 (climatic proximity to the occupied optimum). Note that VPI's visual polarity is the *inverse* of the VRS maps: **light areas here are the best habitat, not the worst.** VPI collapses the whole VERA restriction stack into a single suitability-like number and is useful when a downstream analysis needs a scalar habitat-quality proxy.

**Panel 5 — Δ VRS (Delta VRS).** This is a *diagnostic* map, not a suitability map, and it is the one panel where the light-vs-dark rule inverts its meaning. Δ VRS is the gap between the primary and secondary stressor at each pixel. Dark pixels are where a **single variable dominates** — one predictor is much more restrictive than any other. Light pixels are where **two or more predictors are near-equally restrictive** (small Δ VRS). This distinction matters for management: a dark-Δ pixel can potentially be relieved by addressing the one dominant stressor, while a light-Δ pixel presents multiple co-limiting factors simultaneously.

> **Reading rule for Panels 2–4.** *Light = good, dark = bad.* Panel 5 (Δ VRS) is diagnostic and must be read separately — it describes the *shape* of the restriction, not its intensity.

---

<a id="figure-5"></a>
### Figure 5 — Mahalanobis Distance and Empirical Climatic Tiers

![Figure 5](Figure_5.png)

Figure 5 shows VERA's covariance-aware multivariate calibration. Where Figure 4 summarises restriction one variable at a time and then averages, Figure 5 treats all predictors *simultaneously* under their joint covariance structure.

**Panel (a) — Covariance-Aware Mahalanobis Distance.** For every pixel, the Mahalanobis distance to the occupied climatic centroid, computed with a regularised covariance matrix built from the actual occurrences. Light pixels are close to the multivariate niche core; dark pixels are far. Unlike VRS, this metric fully accounts for correlations between predictors — a pixel that is moderately warm and moderately dry is treated differently from a pixel that is moderately warm and moderately wet, even if both would have the same VRS.

**Panel (b) — Empirical Cumulative Probability.** This is the *calibration engine* that produces Panel (c). The curve is the empirical CDF of Mahalanobis distances at the core occurrences. Three horizontal reference lines mark the **80th, 95th, and 99th percentiles** of the core distribution. Their intersections with the CDF define three **empirical distance thresholds** — dashed vertical lines dropped from those intersections. These thresholds are not fixed; they are derived from the species' own occupied climate distribution, so the tier boundaries mean the same thing across any taxon: the boundary between "core climate" and "moderate departure" is always the distance beyond which only 20 % of core occurrences are found, and so on.

**Panel (c) — Empirical Climatic Tiers.** The thresholds from Panel (b) are applied to Panel (a), collapsing the continuous distance surface into four categorical zones:

- **Core climate** (lightest) — the multivariate niche core; the species' typical conditions.
- **Moderate departure** — tolerable but distinguishable from the core.
- **Restriction zone** — significant multivariate departure; local persistence is possible but marginal.
- **High extrapolative stress** (darkest) — the environment lies beyond the multivariate envelope the species is known to occupy; predictions here are extrapolative.

Because the thresholds come from the core occurrence distribution itself, the tier map is directly comparable across species that were run through the same pipeline.

**Panel (d) — Occurrence-level Mahalanobis density.** The complement to Panel 1 of Figure 4. Every cleaned occurrence is scored by its Mahalanobis distance to the core centroid, and the density is split by partition. **Core Occupied** occurrences (light pink) cluster at low distances with a clean unimodal peak. **Peripheral** occurrences (dark purple) peak *past* the dashed threshold and carry a long right tail — proof that the peripheral partition is genuinely a climatically distinct subset of the observed range, not a random-noise artefact.

> **What Figure 5 adds beyond Figure 4.** VRS is variable-additive and asymmetry-sensitive. Mahalanobis is covariance-aware and asymmetry-blind. The two are *complementary*, not redundant. Areas where the two agree (see [Figure 7](#figure-7)) are the most robust conclusions; areas where they disagree reveal where VERA's asymmetric geometry adds information over a traditional multivariate baseline.

---

<a id="figure-6"></a>
### Figure 6 — Categorical Stressor Attribution

![Figure 6](Figure_6.png)

Figures 4 and 5 tell you *how restricted* each pixel is. Figure 6 tells you *which variable is doing the restricting* at each pixel — the categorical, attribution face of VERA. Colours here do not encode stress intensity; each colour is a specific predictor. The legend at the bottom pairs each colour with a variable and its percentage share of the landscape.

Six panels cover complementary attribution questions.

**Panel (a) — Primary Stressor.** For each pixel, the single variable with the largest Z² contribution. Following Liebig's Law of the Minimum, this is the single limiting factor at that pixel. Bio4 (Temperature Seasonality) dominates the eastern half of Türkiye with 29.2 % landscape share — that region's most restrictive climatic feature for Sitta krueperi is its temperature seasonality.

**Panel (b) — Secondary Stressor.** For each pixel, the variable with the *second-largest* Z². If the primary stressor were somehow relieved, this variable would step forward as the next-binding constraint. CNT (Continentality) is the leading secondary stressor at 19.1 %, appearing in dark green wherever Bio4 is primary — signalling a coupled seasonality/continentality bottleneck in Türkiye's interior.

**Panel (c) — Alternate Primary View.** A companion primary-stressor view emphasising a different colour palette; the interpretation is identical to Panel (a).

**Panel (d) — Primary Stressor Uniqueness.** A diagnostic map with only two categories.
- **Light pixels** — a *unique* single predictor clearly dominates.
- **Dark plum pixels** — multiple predictors are *tied at the Z² = 16 cap*. In these regions the environment is simultaneously intolerable along several dimensions — a compound failure state, not a single-variable bottleneck.

This distinction matters. A dark-plum pixel cannot be "rescued" by relieving one variable, because several are already saturated at the cap.

**Panel (e) — TooHigh Stressor.** Attributes the primary restriction *only where the variable exceeds μ*. Answers: "where is the environment *too much* for the species — too warm, too wet, too seasonal?"

**Panel (f) — TooLow Stressor.** Attributes the primary restriction *only where the variable falls below μ*. Answers: "where is the environment *insufficient* — too cool, too dry, too aseasonal?"

Comparing Panels (e) and (f) directly localises the *direction* of primary restriction in space — a piece of information no traditional SDM output provides.

> **What Figure 6 does that other frameworks cannot.** A traditional SDM tells you a pixel is unsuitable. VERA tells you *which variable* made it unsuitable, whether that variable is *too high or too low*, and whether the pixel is a *single-variable bottleneck* or a *multi-variable collapse*. That is a substantive shift from prediction to diagnosis.

---

<a id="figure-7"></a>
### Figure 7 — Add-On Branch A: Metric-Agreement Diagnostic

![Figure 7](Figure_7.png)

The two add-on branches step back from the primary diagnostic and cross-check VERA's outputs against each other and against alternative interpretations. **Branch A** compares VERA's asymmetric restriction signal (VRS) against the covariance-aware multivariate baseline (Mahalanobis).

The five panels form a strict derivation pipeline — each panel is derived from the ones above it. The first four are intermediate mathematical steps and are shown without individual legends because the story they tell together is what matters; the fifth panel is the categorical synthesis.

1. **Mean VRS Percentile Rank** — VERA's asymmetric restriction expressed as a landscape-relative percentile.
2. **Mahalanobis Percentile Rank** — the covariance-aware multivariate distance expressed the same way.
3. **Percentile Rank Divergence (VRS minus Mahalanobis)** — the signed difference. Purple / magenta pixels are places where VRS ranks the pixel more restrictive than Mahalanobis does; blue pixels are the reverse. Signs matter here — this map answers *which model calls this pixel worse*.
4. **Absolute Rank Disagreement** — the same difference stripped of sign. Answers *how far apart the two models are*, regardless of which is higher.
5. **Mahalanobis–VRS Agreement Classes (bottom panel)** — the categorical synthesis, with a proper legend.
   - **VRS-dominant** (magenta) — pixels where VRS calls the pixel much more restrictive than Mahalanobis. These are the places where VERA's asymmetric geometry sees a bottleneck the traditional covariance-aware distance misses. Ecologically, these are pixels where a *single* variable is far outside the species' tolerance on the harsh tail, even though the pixel is not particularly unusual in multivariate space.
   - **Mahalanobis-dominant** (blue) — the reverse. Mahalanobis flags the pixel as far from the multivariate niche core, but no single predictor is severely asymmetric-restrictive. These are pixels with an unusual *combination* of otherwise-tolerable conditions.
   - **Concordant higher-rank** (purple) — both metrics agree the pixel is highly restrictive.
   - **Concordant lower-rank** (grey) — both metrics agree the pixel is close to the niche core.

The concordant classes together form the *robust conclusions*: the two independent geometries agree, so the assessment is unlikely to be an artefact of either. The dominant classes together form the *diagnostic surface*: they show precisely where the choice between asymmetric-univariate and symmetric-multivariate reasoning changes the answer.

---

<a id="figure-8"></a>
### Figure 8 — Add-On Branch B: Directional Attribution

![Figure 8](Figure_8.png)

**Branch B** is independent from Branch A. It asks a different question: *when the environment is restrictive, is that restriction coming from the upper tail (values above μ) or the lower tail (values below μ)?*

Two panels show two ways of aggregating the answer.

**Panel (a) — Dominant Tail Direction.** For each pixel, examines the *single largest* stressor contribution and asks which side of μ it sits on. Follows Liebig's-Law logic: the direction of the single most severe stressor wins.

**Panel (b) — Net Tail Direction.** For each pixel, sums all upper-tail Z² contributions and all lower-tail Z² contributions across every predictor, then reports which side has the larger cumulative sum. This is a macro-attribution view — the direction of overall restrictive pressure across the whole predictor stack.

Three classes appear on both maps:

- **Upper-tail dominant** (red) — restriction is driven by *exceeding* the species' optima. In this region, the environment is too warm, too wet, too seasonal, too continental — depending on which variable dominates locally.
- **Lower-tail dominant** (yellow) — restriction is driven by *falling short of* the species' optima. Too cool, too dry, too aseasonal.
- **Balanced / near-symmetric** (blue) — restriction is roughly equal from both sides. The species is squeezed on both edges: one variable is too high, another is too low, and neither wins the tail contest.

Panels (a) and (b) usually broadly agree, and where they diverge, the divergence itself is informative. Areas where the dominant map calls "upper-tail" while the net map calls "lower-tail" (or vice versa) are pixels where a single sharp upper-tail spike dominates one metric but is outweighed by many small lower-tail contributions in the other. In practice such pixels flag *dimensional imbalance* — one variable is severely restrictive from one tail, but the aggregate pressure comes from the opposite tail.

Together, Figures 6 and 8 form a matched attribution pair. Figure 6 says *which variable*. Figure 8 says *which direction*. Answering both, at every pixel, is the specific contribution VERA makes beyond suitability-only frameworks.

---

<a id="which-script-produced-which-figure"></a>
## Which Script Produced Which Figure

Every figure above is a direct output of the renderer scripts consuming the finished VERA result trees. No figure required manual assembly beyond the standard VERA outputs.

| Figure | Renderer script | Output subfolder | Notes |
|---|---|---|---|
| 1 | (bespoke educational figure) | `response_curves/` | Conceptual teaching figure — one predictor, two panels |
| 2 | `05_render_core_outputs.R` + native-unit density script | `response_curves/` | Top-10 grid in native units |
| 3 | `05_render_core_outputs.R` + aligned-optimum ridge script | `response_curves/` | Top-10 aligned-optimum ridges |
| 4 | `05_render_core_outputs.R` | `core_renderer_plum/` | Occurrence partition + Mean/Max/Δ VRS + VPI |
| 5 | `07_render_mahalanobis_tiers.R` | `figure4_tier_calibration/` | Panels a–d of the empirical tier calibration |
| 6 | `05_render_core_outputs.R` | `core_renderer_plum/` | Primary / secondary / TooHigh / TooLow attribution |
| 7 | `06_render_addons.R` | `addon_renderer_plum/` | VRS vs Mahalanobis agreement pipeline (Branch A) |
| 8 | `06_render_addons.R` | `addon_renderer_plum/` | Dominant vs Net tail direction (Branch B) |

Every subfolder above sits under `C:/VERA/Results/{19 or 36}/Images/`. Running any renderer on both result trees produces mirror-image galleries for the two predictor profiles.

---

<a id="reproducibility"></a>
## Reproducibility

Every VERA run writes four files at the root of the result directory that together guarantee reproducibility:

- `Skr_sessionInfo.txt` — R version, platform, and every attached package version at the moment of the run.
- `Skr_configuration.R` — the fully-resolved `cfg` list, sufficient to rerun the analysis without re-editing the source script.
- `Skr_code_manifest.csv` — the SHA-256 style MD5 checksum of the analysis script itself, plus a UTC timestamp.
- `Skr_output_inventory.csv` — every file the run produced, with size and MD5 checksum.

Combined with the initial script checksums in `metadata/script_manifest.csv`, this makes any published figure traceable to a specific script version, a specific configuration, and a specific set of input files. If any figure ever needs to be regenerated, the same script + same inputs will produce byte-identical outputs.

> **Editing rule for this repository.** Do not silently update a released script. Changes require a version increment, a checksum update, a changelog entry, and regeneration of every affected output. The scripts deliberately retain explicit local path blocks (`C:/VERA/...`) — this makes the required inputs visible to tutorial users but means paths must be edited before execution on another computer.

---

*Krüper's Nuthatch tutorial for the VERA framework. Full source code, cleaned occurrence table, predictor stack, and all diagnostic outputs are provided in this repository. For questions or contributions, open an issue.*
