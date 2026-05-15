# Project root README

Last updated: 2026-04-22

This workspace is now organized around one current manuscript-oriented mainline:

1. GSE268807 4-sample astrocyte pilot
2. pySCENIC regulon inference and AUCell integration
3. CellOracle round-1 in silico KO
4. sample-level robustness
5. functional interpretation / GO / KEGG / program convergence
6. supportive external/contextual analysis
7. exploratory drug-signature / after-BBB mechanism-direction clue module

Current scope: dry-lab only. No wet-lab validation, no SCENIC+ restart, no expanded high-memory multiome workflow, and no docking are part of the current submission plan.

## Critical Input And Core Object

Pilot input provenance is kept in:

- `input/gse268807_astrocyte_pilot/`

The core discovery object is:

- `final_exports/astrocyte_pilot_rna_with_pyscenic_auc.h5ad`

A compatibility copy is also available at:

- `final_exports/core_objects/astrocyte_pilot_rna_with_pyscenic_auc.h5ad`

The original path is intentionally preserved because multiple scripts still reference it directly.

## Key Result Directories

- `output/`: original pySCENIC output files, including GRN, regulons, and AUC matrix.
- `final_exports/`: final merged discovery object and manuscript-facing compatibility exports.
- `analysis_outputs/`: discovery audit, regulon comparison, shortlist, structural checks, and migrated GSE268807 pilot audit notes.
- `celloracle_run/`: CellOracle prepared data, network files, and round-1 KO outputs.
- `robustness_validation/`: sample-level robustness and pseudobulk robustness outputs.
- `functional_interpretation/`: TF gene sets, GO/KEGG enrichment, and program convergence outputs.
- `external_validation/`: GSE140393 single-group supportive expression analysis.
- `external_validation_round2/`: GSE190452 cross-syndrome supportive analysis.
- `drug_repositioning/`: current v2.1 exploratory drug-signature and after-BBB mechanism-direction clue results.
- `manuscript_output/`: design books, report drafts, presentation assets, and manuscript-facing material.

## Archive / Legacy

Archive and legacy content should not be used as current main conclusions unless explicitly noted:

- `archive/legacy_geo_indexes/`: root-level GEO download/index leftovers moved out of the active project root.
- `archive/deprecated_route_notes/`: reserved for deprecated route notes.
- `archive/old_exploration_notes/`: reserved for old exploration notes.
- `drug_repositioning/archive/v2_pre_v21_process/`: v1/v2/pre-BBB process audit materials. Current drug result is v2.1 after-BBB in `drug_repositioning/`.
- `scripts/drug_repositioning/archive/legacy_v1_scripts/`: old v1 drug repositioning scripts.

## Compatibility Rule

This cleanup intentionally avoids breaking existing scripts. When scripts had hard-coded paths, files were copied into clearer submission-oriented subfolders while originals were kept in place.
