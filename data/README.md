# Data requirements and redistribution policy

The full-resolution tutorial inputs are intentionally not committed to Git.

## Required inputs

The current workflow expects:

1. an occurrence CSV containing longitude and latitude columns;
2. a valid reference raster for domain cleaning and pixel-level thinning;
3. nineteen current bioclimatic rasters named `Bio1` through `Bio19`;
4. seventeen ENVIREM rasters for the 36-predictor profile;
5. monthly precipitation and temperature rasters when ENVIREM layers must be
   regenerated.

## Before publishing an example dataset

- verify the original data source and citation;
- verify redistribution permission;
- check whether exact coordinates require masking or generalisation;
- record the download or access date;
- record the coordinate reference system and raster resolution;
- calculate SHA-256 checksums;
- complete `metadata/data_manifest.csv`.

## Proposed archive strategy

Git should contain only a small, openly redistributable example. Full raster
inputs and reference outputs should be deposited in a DOI-bearing repository,
with the DOI and checksums recorded here.
