# Directory cleanup log

Date: 2026-04-22

Goal: minimal-disruption, manuscript-oriented directory cleanup for the current SCI submission mainline.

## Path Risk Scan

Scripts contain hard-coded references to:

- `final_exports/`
- `analysis_outputs/`
- `manuscript_output/`
- `drug_repositioning/`

Decision:

- Do not rename `analysis_outputs/`.
- Do not move existing root-level files in `final_exports/`.
- Do not move existing root-level report files in `manuscript_output/`.
- Use compatibility copies and README/index files instead.

## Directories Created

- `archive/`
- `archive/legacy_geo_indexes/`
- `archive/deprecated_route_notes/`
- `archive/old_exploration_notes/`
- `final_exports/core_objects/`
- `final_exports/key_tables/`
- `final_exports/key_figures/`
- `final_exports/manuscript_ready_exports/`
- `manuscript_output/design_books/`
- `manuscript_output/manuscript_drafts/`
- `manuscript_output/figures_main/`
- `manuscript_output/figures_supplementary/`
- `manuscript_output/tables_main/`
- `manuscript_output/tables_supplementary/`
- `manuscript_output/presentation_materials/`
- `manuscript_output/presentation_materials/current_project_plan_assets/`

## Files Copied For Manuscript-Oriented Access

Original files remain in place for compatibility.

final_exports:

- `astrocyte_pilot_rna_with_pyscenic_auc.h5ad` -> `final_exports/core_objects/`
- `auc_mtx_matched_to_h5ad.csv` -> `final_exports/key_tables/`
- `regulon_auc_mean.csv` -> `final_exports/key_tables/` and `final_exports/manuscript_ready_exports/`
- `regulon_names.txt` -> `final_exports/key_tables/` and `final_exports/manuscript_ready_exports/`
- `merge_report.txt` -> `final_exports/manuscript_ready_exports/`

manuscript_output:

- `TLE_pyscenic_celloracle_round1_report.docx` -> `manuscript_output/manuscript_drafts/`
- `TLE_pyscenic_celloracle_round1_report.txt` -> `manuscript_output/manuscript_drafts/`
- Legacy round-2 report drafts with over-strong historical filenames were moved to `manuscript_output/manuscript_drafts/`; these names are not used as manuscript claims.
- `design_docs/癫痫亚型_TF调控网络_in_silicoKO_研究设计书_V4_20260422_修订版.docx` -> `manuscript_output/design_books/`
- `current_project_plan_assets/*` -> `manuscript_output/presentation_materials/current_project_plan_assets/`

## Files Moved

Root-level legacy GEO/index files moved to `archive/legacy_geo_indexes/`:

- `external_validation_round2_GSE190452_family.soft.gz`
- `external_validation_round2_GSE190452_filelist.txt`
- `external_validation_round2_GSE190452_suppl_index.html`
- `external_validation_round2_gse275302_full.txt`
- `external_validation_round2_gse275302_page.html`
- `gse190452_suppl_index_tmp.html`

## Files Not Moved Due To Compatibility Risk

- `final_exports/astrocyte_pilot_rna_with_pyscenic_auc.h5ad`
- `final_exports/merge_report.txt`
- all root-level report files under `manuscript_output/`
- the directory `analysis_outputs/`

Reason: active scripts and report builders reference these paths directly.
