#!/usr/bin/env Rscript

# Parse-only smoke test. This confirms that the repository scripts are present
# and syntactically readable without executing analyses or requiring input data.

expected_scripts <- c(
  "scripts/01_envirem_generation.R",
  "scripts/02_occurrence_preparation.R",
  "scripts/03_vera_19.R",
  "scripts/04_vera_19_36.R",
  "scripts/05_render_core_outputs.R",
  "scripts/06_render_addons.R",
  "scripts/07_render_mahalanobis_tiers.R"
)

missing_scripts <- expected_scripts[!file.exists(expected_scripts)]
if (length(missing_scripts)) {
  stop("Missing repository scripts: ", paste(missing_scripts, collapse = ", "))
}

for (script in expected_scripts) {
  parse(file = script, keep.source = TRUE)
  message("PARSE PASS: ", script)
}

required_repository_files <- c(
  "README.md",
  "CITATION.cff",
  "REPRODUCIBILITY.md",
  "metadata/script_manifest.csv",
  "metadata/data_manifest.csv",
  "website/_quarto.yml"
)

missing_repository_files <- required_repository_files[
  !file.exists(required_repository_files)
]
if (length(missing_repository_files)) {
  stop("Missing repository files: ",
       paste(missing_repository_files, collapse = ", "))
}

message("Repository smoke test completed successfully.")