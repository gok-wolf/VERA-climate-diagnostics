# VERA: Variable Ecological Restriction Analysis

VERA is an occurrence-calibrated framework for mapping the magnitude, identity,
direction and diagnostic resolution of climatic departure from a species'
occupied environmental reference.

This repository is under active development and currently provides a complete
single-species tutorial using *Sitta krueperi* (Krüper's nuthatch). The tutorial
contains independent 19-predictor and 36-predictor workflows, the core VERA
outputs, covariance-aware Mahalanobis diagnostics, directional add-ons and
publication-oriented renderers.

> **Interpretive boundary:** VERA is a diagnostic and hypothesis-generating
> framework. Its outputs are not occurrence probabilities, habitat-suitability
> estimates, physiological tolerance limits or demonstrations of causal range
> limitation.

## Repository status

- Development stage: `0.1.0-dev`
- Tutorial taxon: *Sitta krueperi*
- Predictor profiles: 19 bioclimatic variables and 36 bioclimatic + ENVIREM
  variables
- Manuscript relationship: methodological implementation associated with the
  VERA Paper 1 project; publication details will be added after acceptance
- Public DOI: not yet assigned

## Workflow

```text
Monthly climate rasters
        |
        v
ENVIREM generation -----> 36-predictor stack
        |
Occurrence cleaning and pixel-level thinning
        |
        +--------------------+
        |                    |
        v                    v
19-predictor VERA      36-predictor VERA
        |                    |
        +----------+---------+
                   v
      Core, add-on and tier renderers
```

## Scripts

| Order | Script | Purpose |
|---:|---|---|
| 1 | `scripts/01_envirem_generation.R` | Generates the ENVIREM predictor layers from monthly climate rasters. |
| 2 | `scripts/02_occurrence_preparation.R` | Removes records outside the raster domain, thins records by raster cell and performs the documented geographic isolation step. |
| 3 | `scripts/03_vera_19.R` | Runs the complete single-species VERA workflow with 19 bioclimatic predictors. |
| 4 | `scripts/04_vera_19_36.R` | Runs independent 19- and 36-predictor VERA workflows. This is the primary dual-profile tutorial. |
| 5 | `scripts/05_render_core_outputs.R` | Produces maps and plots for the core VERA output family. |
| 6 | `scripts/06_render_addons.R` | Produces maps and summaries for the directional and cross-geometry add-ons. |
| 7 | `scripts/07_render_mahalanobis_tiers.R` | Produces the Mahalanobis core-partition and empirical-tier calibration figure. |

The original script checksums are recorded in
`metadata/script_manifest.csv`. Do not silently edit a released script without
updating its checksum, version and changelog entry.

The expected per-profile archive is documented in `EXPECTED_OUTPUTS.md`.

## Requirements

The scripts require a recent R installation and use the following packages:

```r
terra
dplyr
tibble
envirem
ggplot2
readr
scales
sf
tidyterra
shadowtext
```

Package versions used for the first public release will be frozen in an
`renv.lock` file after the pipeline has been tested on a clean machine.

## Quick start

1. Clone or download this repository.
2. Obtain the occurrence and climatic input data described in `data/README.md`.
3. Edit the path block at the beginning of the relevant script.
4. Run the dual-profile tutorial:

```r
source("scripts/04_vera_19_36.R")
```

5. Run the renderers after both analysis profiles finish:

```r
source("scripts/05_render_core_outputs.R")
source("scripts/06_render_addons.R")
source("scripts/07_render_mahalanobis_tiers.R")
```

The current tutorial scripts contain explicit Windows paths for transparency.
A portable path interface will be added only after the canonical calculations
and output inventory pass a clean-machine reproducibility audit.

## Data policy

Full-resolution climate rasters are not stored in this Git repository. The
repository will contain only a small, openly redistributable example dataset,
download instructions, source citations, licences and checksums. Large input
and demonstration archives will be deposited separately in a DOI-bearing data
repository.

Occurrence coordinates will be published only when their source licences and
species-sensitivity requirements permit redistribution.

## Documentation

The website source is stored in `website/` and will be rendered to `docs/` with
Quarto. Until the website is published, the Markdown source remains the
authoritative documentation.

Deployment instructions are provided in `WEBSITE_DEPLOYMENT.md`. The intended
site address is `https://gok-wolf.github.io/VERA-climate-diagnostics/`.

## Citation

Citation metadata are provided in `CITATION.cff`. The author list, manuscript
reference and DOI are provisional until the first archived release.

## Licence

The final software, documentation and example-data licences have not yet been
selected. See `LICENSE_DECISION.md`. Until explicit licences are added, no
permission to redistribute or reuse repository content should be inferred.

## Contact and contributions

The repository is currently maintained by the VERA development team. Issues and
contribution guidance will be enabled after the first internal reproducibility
release.