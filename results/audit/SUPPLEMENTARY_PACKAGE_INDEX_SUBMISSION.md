# Supplementary Package Index (Submission Refresh)

This supplementary package was generated from the current local project state without rerunning upstream analyses.

## Mainline boundaries

- Current manuscript mainline: GSE268807 astrocyte pilot + pySCENIC + CellOracle + sample-level robustness + functional interpretation/program convergence + supportive external analysis.
- GSE140393 is treated as **single-group supportive analysis**.
- GSE190452 is treated as **cross-syndrome supportive analysis**.
- Neither GSE140393 nor GSE190452 is described here as formal same-disease external validation.
- No SCENIC+ current mainline is used.
- No wet-lab validation is included in the current package.
- Drug repositioning remains **exploratory only** and is kept in supplementary material.

## Main-text tables

- `manuscript_output/tables_main/Table1_cohort_and_analytical_objects_submission.csv`: Compact Table 1 for the manuscript main text.
- `manuscript_output/tables_main/Table2_shortlist_main_submission.csv`: Main TF shortlist table supporting Fig. 2.
- `manuscript_output/tables_main/Table3_integrated_priority_main_submission.csv`: Integrated priority table supporting Fig. 3 / Fig. 4.

## Supplementary figures

- `FigS1_discovery_object_qc_and_input_overview_submission.png`: Discovery object QC, sample composition, upstream QC/subsampling trace, and current object lineage for the GSE268807 astrocyte pilot.
  - Main-text linkage: extends manuscript Fig. 1 / discovery object context.
  - Data sources: `analysis_outputs/gse268807_astrocyte_pilot/qc_metrics/`, `celloracle_run/prepared_data/`, `final_exports/merge_report.txt`, `robustness_validation/01_input_object_summary.txt`.

- `FigS2_full_differential_regulon_landscape_submission.png`: Extended differential regulon landscape with a larger heatmap and broader lesion/internal_control-associated regulon ranking view.
  - Main-text linkage: extends Fig. 2 and Table 2.
  - Data sources: `analysis_outputs/group_compare/`, `final_exports/auc_mtx_matched_to_h5ad.csv`, `celloracle_run/prepared_data/celloracle_round1_metadata.csv`.

- `FigS3_extended_celloracle_prioritization_submission.png`: Extended CellOracle prioritization with detailed metric heatmap and representative perturbation panels for all round1 TFs.
  - Main-text linkage: extends Fig. 3 and Table 3.
  - Data sources: `celloracle_run/ko_round1/`.

- `FigS4_sample_level_robustness_full_submission.png`: Expanded **sample-level robustness** view, including leave-one-sample-out expression/regulon effects, pseudobulk support, and shortlist sensitivity.
  - Main-text linkage: extends Fig. 4 robustness panels and Table 3.
  - Data sources: `robustness_validation/`.

- `FigS5_supportive_external_evidence_full_submission.png`: Full supportive external evidence package, integrating GSE140393 single-group supportive analysis and GSE190452 cross-syndrome supportive analysis.
  - Main-text linkage: extends Fig. 4 supportive external evidence panels and Table 3.
  - Data sources: `external_validation/`, `external_validation_round2/`.

- `FigS6_functional_interpretation_and_program_convergence_full_submission.png`: Real GO/KEGG enrichment dot plots and program convergence summaries generated from current functional interpretation result tables.
  - Main-text linkage: supports the mechanistic interpretation layer that feeds into the Fig. 5 conceptual synthesis.
  - Data sources: `functional_interpretation/`.

- `FigS7_drug_repositioning_v21_after_bbb_summary_submission.png`: Exploratory after-BBB drug repositioning v2.1 summary, kept only as supplementary mechanism-direction clue material.
  - Main-text linkage: corresponds to the exploratory drug clue module and should remain supplementary.
  - Data sources: `drug_repositioning/13_integrated_after_bbb_v21.csv`, `drug_repositioning/16_paper_ready_table_after_bbb_v21.csv`.

## Supplementary tables

- `manuscript_output/tables_supplementary/TableS1_full_differential_regulons_submission.csv`: Full differential regulon statistics; extends Fig. 2.
- `manuscript_output/tables_supplementary/TableS2_candidate_tf_shortlist_extended_submission.csv`: Extended candidate TF shortlist / secondary candidate metrics; extends Table 2.
- `manuscript_output/tables_supplementary/TableS3_celloracle_round1_full_metrics_submission.csv`: Full CellOracle round1 metrics; extends Fig. 3 and Table 3.
- `manuscript_output/tables_supplementary/TableS4_sample_level_robustness_full_submission.csv`: Detailed **sample-level robustness** package, including leave-one-sample-out, pseudobulk, shortlist sensitivity, and integrated summaries.
- `manuscript_output/tables_supplementary/TableS5_functional_interpretation_program_convergence_full_submission.csv`: Full GO/KEGG/program convergence package.
- `manuscript_output/tables_supplementary/TableS6_supportive_external_analysis_full_submission.csv`: Full supportive external analysis package for GSE140393 and GSE190452.
- `manuscript_output/tables_supplementary/TableS7_drug_repositioning_v21_after_bbb_full_submission.csv`: Full after-BBB drug repositioning v2.1 table from the current latest route.

## Current main-text anchors already present

- `manuscript_output/figures_main/Figure1_study_design.png`
- `manuscript_output/figures_main/Figure2_regulon_landscape.png`
- `manuscript_output/figures_main/Figure3_celloracle_perturbation.png`
- `manuscript_output/figures_main/Figure4_robustness_support.png`
- `manuscript_output/figures_main/Figure5_program_model.png`
- `manuscript_output/tables_main/Table1_cohort_and_analytical_objects_submission.csv`
- `manuscript_output/tables_main/Table2_shortlist_main_submission.csv`
- `manuscript_output/tables_main/Table3_integrated_priority_main_submission.csv`
