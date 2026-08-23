# Pre-release checklist

## Scientific identity

- [ ] Confirm the official repository name.
- [ ] Confirm the complete author and contributor list.
- [ ] Confirm the Paper 1 title and citation status.
- [ ] Confirm which script is declared the canonical public implementation.
- [ ] Confirm that the 19-only script is labelled as a convenience profile.

## Code

- [x] Preserve the initial uploaded-script checksums and document subsequent
      naming-only executable changes in the changelog and script manifest.
- [x] Record initial SHA-256 checksums.
- [ ] Run all scripts on a clean Windows environment.
- [ ] Parse all R scripts in continuous integration.
- [ ] Confirm that 19- and 36-profile output inventories are complete.
- [ ] Confirm that renderer inputs match the canonical output filenames.
- [ ] Decide whether hard-coded paths remain in the tutorial release or are
      replaced in a separately versioned portability update.

## Data

- [ ] Complete source citations and access dates.
- [ ] Confirm occurrence-data redistribution permission.
- [ ] Review coordinate sensitivity for *Sitta krueperi*.
- [ ] Confirm climate-raster redistribution terms.
- [ ] Prepare a small openly redistributable example dataset.
- [ ] Deposit full demonstration inputs and outputs in a DOI-bearing archive.
- [ ] Complete all checksums in `metadata/data_manifest.csv`.

## Software environment

- [ ] Record R version.
- [ ] Freeze package versions in `renv.lock`.
- [ ] Record system dependencies required by `terra`, `sf` and renderers.
- [ ] Test on a second machine or operating system where practical.

## Documentation

- [x] Draft README.
- [x] Draft Quarto website structure.
- [ ] Add representative, verified figures.
- [ ] Add a predictor dictionary and abbreviations table.
- [ ] Add expected output-tree examples.
- [ ] Render and inspect the GitHub Pages site.
- [ ] Replace all `OWNER` URL placeholders.

## Legal and citation

- [ ] Select the software licence.
- [ ] Select the documentation licence.
- [ ] Confirm licences for example data and figures.
- [ ] Finalise `CITATION.cff`.
- [ ] Create a GitHub release.
- [ ] Connect the release to Zenodo and record the DOI.

## Public-release gate

- [ ] No unresolved validation failures.
- [ ] No restricted or unlicensed data committed.
- [ ] No personal local paths exposed unintentionally.
- [ ] Documentation matches the released code and outputs.
- [ ] Final internal approval recorded.
