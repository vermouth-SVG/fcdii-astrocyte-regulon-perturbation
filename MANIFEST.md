# Repository manifest

Audited manuscript version: `Manuscript.docx`

Manuscript title: `Single-cell regulatory analysis identifies astrocyte stress response and homeostatic support in focal cortical dysplasia type II`

Audit date: 2026-09-01

## 1. Manuscript files

| File | Purpose |
|---|---|
| `Manuscript.docx` | Current BMC Genomics submission manuscript. |
| `Supplementary_Figures_S1_S7_preview_not_for_submission.docx` | Current Supplementary Figures S1-S7 review package. |
| `Supplementary_Tables_ST01_ST09_submission_ready.xlsx` | Current Supplementary Tables ST01-ST09 workbook. |
| `RELEASE_NOTES_v0.3.0-submission.md` | Submission-version release notes and interpretation boundaries. |

## 2. Scripts

| Folder | Main scripts | Supported module |
|---|---|---|
| `scripts/gse268807_pilot/` | input checks and object preparation scripts | Data preprocessing and discovery object construction. |
| `scripts/` | pySCENIC, AUCell merge, differential regulon, candidate TF, and CellOracle scripts | pySCENIC/AUCell, differential regulon analysis, candidate TF prioritization, in silico perturbation simulation. |
| `scripts/robustness/` | robustness scripts | Robustness and sensitivity analyses. |
| `scripts/functional/` | functional interpretation scripts | Functional enrichment and convergence analysis. |
| `scripts/drug_repositioning/` | drug-signature annotation scripts | Exploratory DSigDB / SwissADME annotation. |
| `scripts/manuscript_figures/` | Figure 1-5 and Fig. S1-S7 builders | Manuscript figure generation. |
| `scripts/manuscript_tables/` | Main and supplementary table builders | Manuscript table generation. |

## 3. Results / tables

| Manuscript item | Release files | Purpose |
|---|---|---|
| Table 1 | `results/tables/manuscript/main/Main_Table_1.csv` | Analytical datasets and object roles. |
| Table 2 | `results/tables/manuscript/main/Main_Table_2.csv` | Integrated prioritization of candidate TFs. |
| ST01-ST09 | `results/tables/manuscript/supplementary/Supplementary_Table_S1.csv` through `Supplementary_Table_S9.csv` | Supplementary table CSV exports matching the manuscript. |
| pySCENIC/AUCell summaries | `results/tables/pyscenic_overview/` | Regulon activity summaries and retained regulon names. |
| CellOracle perturbation summaries | `results/tables/celloracle_ko/` | Per-TF perturbation metrics and Approx. RI summaries. Legacy folder name retained. |
| CellOracle formal permutation controls | `results/tables/celloracle_ko/permutation_test/` | Fixed-seed, 10,000-permutation test results and summary. |
| Robustness outputs | `results/tables/robustness/` | LOSO, pseudobulk, and sensitivity summaries. |
| Supportive external outputs | `results/tables/supportive_*` | GSE140393 and GSE190452 support summaries. |
| Functional enrichment outputs | `results/tables/functional_interpretation/` | GO/KEGG and convergence summaries. |
| Exploratory drug-signature outputs | `results/tables/exploratory_drug_signature/` | DSigDB / SwissADME annotation summaries. |

## 4. Results / figures

| Manuscript item | Release files | Purpose |
|---|---|---|
| Figure 1 | `results/figures/manuscript/main/Figure_1_study_design.*` | Study design and analytical framework. |
| Figure 2 | `results/figures/manuscript/main/Figure_2_regulon_landscape.*` | Differential regulon landscape. |
| Figure 3 | `results/figures/manuscript/main/Figure_3_celloracle_perturbation.*` | CellOracle in silico perturbation priority. |
| Figure 4 | `results/figures/manuscript/main/Figure_4_robustness_support.*` | Robustness and supportive evidence. |
| Figure 5 | `results/figures/manuscript/main/Figure_5_program_model.*` | Functional interpretation and program-level working model. |
| Fig. S1-S7 | `results/figures/manuscript/supplementary/Supplementary_Figure_S1.*` through `Supplementary_Figure_S7.*` | Supplementary figure outputs matching the manuscript. |
| CellOracle permutation null distributions | `results/figures/celloracle_ko/permutation_test/` | Formal observed-versus-randomized and group-directionality control distributions. |

## 5. Audit / provenance

| Item | Files |
|---|---|
| pySCENIC/AUCell logs and audit | `results/audit/methods_audit_pyscenic_aucell_reproducibility.md`, `configs/pyscenic_runtime_and_resources.md` |
| CellOracle runtime notes | `configs/celloracle_runtime_and_resources.md` |
| Release audit | `results/audit/repository_release_audit.md` |
| Excluded files manifest | `excluded_files_manifest.csv` |

## 6. Configs / resources

| Resource | File |
|---|---|
| pySCENIC resources | `configs/pyscenic_runtime_and_resources.md` |
| CellOracle runtime | `configs/celloracle_runtime_and_resources.md` |
| Enrichr / DSigDB / SwissADME notes | `environment_notes.md` and `results/tables/exploratory_drug_signature/` |
| External data registry | `metadata/` |

## 7. Excluded files

The manuscript-level release intentionally excludes large raw single-cell matrices, controlled-access raw data, cisTarget ranking databases, motif databases / motif annotation files, Docker image layers, and large intermediate objects.

Legacy path note: file and folder names containing `ko` are retained only to preserve script compatibility and provenance. They do not define current manuscript terminology.
