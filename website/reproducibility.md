---
title: "Reproducibility"
---

## Release discipline

Each archived release should connect one immutable code version to one input
manifest and one complete output inventory. Canonical changes require a new
release rather than replacement of previous files.

## Minimum archive

- scripts and checksums;
- configuration dumps for both profiles;
- session information;
- occurrence-processing counts;
- data provenance and licences;
- raster geometry metadata;
- output inventories;
- renderer versions;
- warnings, failures and protocol deviations.

## Package environment

Generate `renv.lock` only after a successful clean-machine run so it reflects
the verified release environment.

## Data separation

Large rasters and restricted coordinates should be archived outside Git when
licence, sensitivity or repository-size constraints require it. Preserve
provenance, licences, checksums and stable download identifiers in the repo.
