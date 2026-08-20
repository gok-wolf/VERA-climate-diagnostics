# Expected VERA output structure

The dual-profile tutorial writes one complete archive under each configured
result root:

```text
Results/19/Skr_current/
Results/36/Skr_current/
```

Each archive should contain the following structure. Optional products are
identified explicitly.

```text
Skr_current/
├── Skr_code_manifest.csv
├── Skr_configuration.R
├── Skr_output_inventory.csv
├── Skr_sessionInfo.txt
├── csvs/
│   ├── Skr_00_predictor_lookup_dictionary.csv
│   ├── Skr_01_run_manifest.csv
│   ├── Skr_02_model_reference.csv
│   ├── Skr_03_occurrence_diagnostics.csv
│   ├── Skr_04_occurrence_partitions.csv
│   ├── Skr_05_current_diagnostics.csv
│   └── Skr_06_proximity_summary.csv
├── rasters/
│   ├── core_occupied/
│   │   ├── Skr_current_four_tier_status.tif
│   │   ├── Skr_current_mean_vrs.tif
│   │   ├── Skr_current_max_vrs.tif
│   │   ├── Skr_current_mahal_distance.tif
│   │   ├── Skr_current_primary_stressor_factor.tif
│   │   ├── Skr_current_primary_stressor_idx.tif
│   │   └── Skr_current_vpi.tif
│   └── optional_diagnostics/
│       ├── Skr_current_delta_vrs.tif
│       ├── Skr_current_primary_co_dominance_count.tif
│       ├── Skr_current_primary_unique_stressor_factor.tif
│       ├── Skr_current_secondary_stressor_factor.tif
│       ├── Skr_current_secondary_stressor_idx.tif
│       ├── Skr_current_too_high_factor.tif
│       ├── Skr_current_too_high_idx.tif
│       ├── Skr_current_too_high_co_dominance_count.tif
│       ├── Skr_current_too_high_zsq_max.tif
│       ├── Skr_current_too_high_zsq_sum.tif
│       ├── Skr_current_too_low_factor.tif
│       ├── Skr_current_too_low_idx.tif
│       ├── Skr_current_too_low_co_dominance_count.tif
│       ├── Skr_current_too_low_zsq_max.tif
│       ├── Skr_current_too_low_zsq_sum.tif
│       └── Skr_current_valid_var_count.tif
└── addons/
    ├── csvs/
    │   ├── addon01_tail_direction_pixel_counts.csv
    │   ├── addon01_tail_direction_dominant_vs_net_agreement.csv
    │   ├── addon01_tail_direction_summary.csv
    │   ├── addon01_tail_direction_latitudinal_profile.csv
    │   ├── addon02_divergence_class_counts.csv
    │   └── addon02_divergence_summary.csv
    └── rasters/
        ├── Skr_current_vrs_rank01.tif
        ├── Skr_current_mahal_rank01.tif
        ├── Skr_current_mahal_vrs_divergence.tif
        ├── Skr_current_mahal_vrs_abs_disagreement.tif
        ├── Skr_current_mahal_vrs_agreement_4class_idx.tif
        ├── Skr_current_mahal_vrs_agreement_4class.tif
        ├── Skr_current_tail_direction_dominant_idx.tif
        ├── Skr_current_tail_direction_dominant.tif
        ├── Skr_current_tail_direction_net_idx.tif
        └── Skr_current_tail_direction_net.tif
```

## Conditional files

- `Skr_current_four_tier_status_fixedbreaks.tif` is written only when the fixed
  break comparator is enabled.
- Integer index rasters are written only when `write_index_rasters` is enabled.
- `Skr_current_valid_var_count.tif` is written only when optional diagnostic
  rasters are enabled.
- The latitudinal add-on table requires a longitude/latitude raster template;
  projected rasters currently trigger a warning instead.

## Validation rule

The output inventory is the authoritative record of files produced by a given
run. This expected tree is a tutorial guide and must be updated whenever the
canonical pipeline deliberately changes its public output contract.
