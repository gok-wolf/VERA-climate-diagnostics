# Development roadmap

## Phase 0 — Repository foundation

- [x] Create the repository architecture.
- [x] Copy the seven tutorial scripts without executable modification.
- [x] Record script checksums.
- [x] Draft README, citation metadata and governance documents.
- [x] Draft the Quarto website source.

## Phase 1 — Internal reproducibility release

- [ ] Run the full 19- and 36-predictor workflows on a clean machine.
- [ ] Verify the complete output inventories against the expected Paper 1
      output families.
- [ ] Freeze package versions with `renv`.
- [ ] Add a small, licensed example dataset.
- [ ] Add verified tutorial figures and expected output-tree examples.
- [ ] Render and inspect the Quarto website.

## Phase 2 — Scientific and legal review

- [ ] Finalise authorship and contributor roles.
- [ ] Select code, documentation and example-data licences.
- [ ] Complete occurrence and climate-data provenance.
- [ ] Review species-coordinate sensitivity.
- [ ] Confirm alignment with the submitted Paper 1 version.

## Phase 3 — Public archival release

- [ ] Replace repository-owner URL placeholders.
- [ ] Publish GitHub release `v0.1.0`.
- [ ] Archive the release and large demonstration files in Zenodo.
- [ ] Add the DOI to README, `CITATION.cff` and the website.
- [ ] Enable GitHub Pages.

## Phase 4 — Post-publication development

- [ ] Add platform-neutral path handling in a separately versioned update.
- [ ] Evaluate conversion of reusable functions into an R package.
- [ ] Add additional taxa as independent case studies.
- [ ] Add automated numerical-regression tests for representative raster cells.
