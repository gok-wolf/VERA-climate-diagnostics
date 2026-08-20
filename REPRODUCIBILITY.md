# Reproducibility policy

## Canonical-code rule

Released analysis scripts are immutable. Any executable change requires:

- a new version number;
- an updated SHA-256 checksum in `metadata/script_manifest.csv`;
- a changelog entry;
- regeneration of affected outputs;
- a documented reason for the change.

Comment-only changes should also be recorded, even when the executable-code
hash remains unchanged.

## Required run records

Each archived VERA run should preserve:

- the exact analysis script and checksum;
- the dumped configuration;
- R and package versions;
- input-data checksums and provenance;
- occurrence filtering counts;
- predictor order;
- random seeds;
- the output inventory and checksums;
- warnings and failed diagnostic checks.

## Clean-machine audit

Before the first public release:

1. install R and dependencies in a clean environment;
2. generate `renv.lock` from the verified package library;
3. run the small example dataset from start to finish;
4. compare the generated inventories and checksums with the reference archive;
5. repeat on at least one second operating system where practical.

## Numerical expectations

GeoTIFF compression metadata and graphics may differ across platforms even when
the scientific values are identical. Validation should therefore distinguish
between exact file hashes and numerical equivalence of raster values and tables.
