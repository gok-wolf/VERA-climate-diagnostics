---
title: "Reproducibility"
---

## Release discipline

Each archived release should connect one immutable code version to one input
manifest and one complete output inventory. Canonical changes require a new
version and changelog entry rather than silent replacement of previous files.

Optional renderers do not change canonical VERA calculations, but any script
used to create a published figure should still be archived with its checksum
and configuration.

## Minimum archive

A reproducible release should preserve:

- every reported script and its SHA-256 checksum;
- the input-data manifest, provenance records and licences;
- configuration dumps for every reported predictor profile;
- R session and package-version information;
- occurrence-processing counts and retained-record identifiers;
- raster CRS, resolution, extent, origin and alignment metadata;
- complete output inventories and file checksums;
- renderer and Species Interpreter versions;
- narrative evidence logs when generated text is reported;
- warnings, failures, exclusions and protocol deviations.

## Automated repository checks

`tests/smoke_test.R` confirms that Scripts 01–10 exist and can be parsed
without executing the analyses. `tests/verify_checksums.py` compares every
script with `metadata/script_manifest.csv`.

Both tests should pass in the GitHub Actions `repository-checks` workflow before
a release is tagged. Whenever a script is added, removed, renamed or edited,
update the smoke-test list, script manifest, package requirements, script
catalogue and changelog together.

## Package environment

Generate `renv.lock` only after a successful clean-machine run so that it
records the verified release environment rather than an incidental development
session. Optional packages should be labelled clearly in the package manifest.

## Data separation

Large rasters and restricted occurrence coordinates may be archived outside
Git when licence, sensitivity or repository-size constraints require it.
Preserve provenance, licences, checksums and stable download identifiers in the
repository so that external data remain traceable to the archived analysis.

