# Current mainline run/order map

Last updated: 2026-04-22

This file describes the current manuscript-facing analysis order. It is not an instruction to rerun analyses.

## 1. Pilot Input Provenance

Directory:

- `input/gse268807_astrocyte_pilot/`

Role:

- Keeps the 4-sample GSE268807 astrocyte pilot provenance.
- Contains metadata, raw MTX triples, source metadata, and four QC h5ad files.

## 2. pySCENIC Raw Outputs

Directory:

- `output/`

Role:

- Stores pySCENIC outputs: `grn_adj.tsv`, `regulons.csv`, and `auc_mtx.csv`.

## 3. Final Discovery Object

Directory:

- `final_exports/`

Core object:

- `final_exports/astrocyte_pilot_rna_with_pyscenic_auc.h5ad`

Role:

- Main GSE268807 astrocyte pilot object with pySCENIC AUC merged back.
- Compatibility copy is in `final_exports/core_objects/`.

## 4. CellOracle

Directory:

- `celloracle_run/`

Role:

- Prepared CellOracle input, network files, and round-1 KO outputs for BHLHE40, NFE2L2, SOX2, and THRB.

## 5. Sample-Level Robustness

Directory:

- `robustness_validation/`

Role:

- leave-one-sample-out/sample-level robustness, pseudobulk robustness checks, shortlist sensitivity, and master summaries.

## 6. Functional Interpretation

Directory:

- `functional_interpretation/`

Role:

- TF gene set construction, GO/KEGG enrichment, program convergence, and integrated functional interpretation.

## 7. Supportive External Analysis: GSE140393

Directory:

- `external_validation/`

Role:

- Single-group lesion-oriented supportive expression analysis.
- This is supportive expression-level evidence, not formal regulon-level replication.

## 8. Supportive External Analysis: GSE190452

Directory:

- `external_validation_round2/`

Role:

- Cross-syndrome supportive analysis in TLE/non-epileptic control context.
- This should be described as supportive/cross-syndrome evidence, not formal same-disease validation.

## 9. Drug Repositioning

Directory:

- `drug_repositioning/`

Role:

- Current v2.1 after-BBB exploratory mechanism-direction module.
- Headline mechanism-direction leads are in `14_primary_leads_after_bbb_v21.csv`.
- This module is supplementary/exploratory and is not a treatment recommendation.

## 10. Manuscript Outputs

Directory:

- `manuscript_output/`

Role:

- Design books, report drafts, figure/table staging folders, and presentation assets for submission preparation.
