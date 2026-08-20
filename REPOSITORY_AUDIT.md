# Initial repository audit — 20 August 2026

## Completed

- Seven supplied R scripts copied into a numbered workflow.
- All seven copied-script SHA-256 checksums verified against the manifest.
- README and Quarto website sources drafted.
- Parse-only R smoke test and checksum-verification test added.
- GitHub Actions check workflow drafted.
- Data, citation, licensing and release checklists created.

## Intentional non-actions

- No canonical R script was edited.
- No occurrence coordinates were added.
- No full-resolution raster was added.
- No licence was selected without author approval.
- No author list, DOI or manuscript citation was invented.
- No `renv.lock` was generated from an unverified local package library.

## Remaining blockers before public release

1. Replace `OWNER` placeholders in repository and website URLs.
2. Confirm the official repository name.
3. Finalise authorship, citation and software licence.
4. Complete data provenance and redistribution review.
5. Test all scripts with R on a clean Windows machine.
6. Generate and validate `renv.lock`.
7. Render the Quarto website and inspect internal links.
8. Add verified example inputs, outputs and figures.

## Environment limitation during scaffold preparation

The preparation workspace did not provide R or Quarto executables. YAML files
were parsed successfully and script checksums were verified, but R parsing and
Quarto rendering must be completed in the clean-machine audit.
