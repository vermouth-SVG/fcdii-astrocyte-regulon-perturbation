# Repository manifest

Audited manuscript version: `Manuscript.docx` (BMC Genomics submission version, 2026-10-25)

Manuscript title: `Single-nucleus transcriptomic reanalysis prioritizes NFE2L2 and THRB regulatory programs in focal cortical dysplasia type II astrocytes`

Audit date: 2026-10-25

## 1. Manuscript and submission files

| File | Purpose |
|---|---|
| `Manuscript.docx` | Current BMC Genomics submission manuscript (2026-10-25 version). |
| `Cover_Letter.docx` | Cover letter to the editors (2026-10-25 version). |
| `Additional_file_1.zip` | Supplementary Figures S1-S6 (TIFF, 600 dpi) with `Supplementary_Figure_Legends.md`. |
| `Additional_file_2.zip` | Supplementary Tables ST01-ST09 (XLSX) with `Supplementary_Table_Legends.md`. |
| `Supplementary_Tables_ST01_ST09_submission_ready.xlsx` | Browsable copy of the ST01-ST09 workbook. |
| `RELEASE_NOTES_20261025.md` | Changes in the 2026-10-25 submission version. |
| `RELEASE_NOTES_v0.4.0-submission.md` | Archived analysis-version release notes (2026-09-10). |
| `RELEASE_NOTES_v0.3.0-submission.md` | Archived v0.3.0 release notes (superseded). |

The former `Supplementary_Figures_S1_S7_preview_not_for_submission.docx` review package was removed in the 2026-10-25 update; Supplementary Figures S1-S6 are distributed as TIFF files inside `Additional_file_1.zip` and in `results/figures/manuscript/supplementary/`.

## 2. Scripts

| Folder | Main scripts | Supported module |
|---|---|---|
| `scripts/gse268807_pilot/` | input checks and object preparation scripts | Data preprocessing and discovery object construction. |
| `scripts/` | pySCENIC, AUCell merge, differential regulon, candidate TF, and CellOracle scripts | pySCENIC/AUCell, differential regulon analysis, candidate TF prioritization, in silico perturbation simulation. |
| `scripts/robustness/` | robustness scripts | Robustness and sensitivity analyses. |
| `scripts/functional/` | functional interpretation scripts | Functional enrichment and convergence analysis. |
| `scripts/drug_repositioning/` | drug-signature annotation scripts | Exploratory DSigDB / SwissADME annotation. |
| `scripts/manuscript_figures/` | Figure 1-5 and supplementary figure builders | Manuscript figure generation. |
| `scripts/manuscript_tables/` | Main and supplementary table builders | Manuscript table generation. |

The scripts are the 2026-09-10 analysis snapshot and are retained unchanged.

## 3. Results / tables

| Manuscript item | Release files | Purpose |
|---|---|---|
| Table 1 | `results/tables/manuscript/main/Main_Table_1.csv` | Analytical datasets and object roles. |
| Table 2 | `results/tables/manuscript/main/Main_Table_2.csv` | Integrated prioritization of candidate TFs. |
| ST01-ST09 | `Additional_file_2.zip` (XLSX) and `results/tables/manuscript/supplementary/Supplementary_Table_S1.csv` through `Supplementary_Table_S9.csv` | Supplementary table exports matching the manuscript. |
| pySCENIC/AUCell summaries | `results/tables/pyscenic_overview/` | Regulon activity summaries and retained regulon names. |
| Current seven-candidate CellOracle summaries | `results/tables/celloracle_perturbation_7tf/` | Fixed-seed overall, group, sample, donor-paired, and metric-sensitivity summaries for all seven candidates. |
| Legacy CellOracle outputs | `results/tables/celloracle_ko/` | Retained for provenance; older four-candidate and recovery-index wording is superseded. |
| Robustness outputs | `results/tables/robustness/` | LOSO, pseudobulk, and sensitivity summaries. |
| Contextual external outputs | `results/tables/supportive_*` | Legacy path names retained; GSE140393 and GSE190452 are contextual, not independent validation. |
| Current submission tables | `results/tables/submission_additional_files_current/` | Submission-facing supplementary table files. |
| Revision boundary | `REVISION_NOTES_20260902.md` | Defines the current seven-candidate evidence and superseded terminology. |
| Functional enrichment outputs | `results/tables/functional_interpretation/` | GO/KEGG and convergence summaries. |
| Exploratory drug-signature outputs | `results/tables/exploratory_drug_signature/` | DSigDB / SwissADME annotation summaries. |

## 4. Results / figures

| Manuscript item | Release files | Purpose |
|---|---|---|
| Figure 1 | `results/figures/manuscript/main/Figure_1_study_design.tiff`, `.png`, `.pdf` | Study design and analytical framework. |
| Figure 2 | `results/figures/manuscript/main/Figure_2_regulon_landscape.tiff`, `.png`, `.pdf` | Differential regulon landscape. |
| Figure 3 | `results/figures/manuscript/main/Figure_3_celloracle_perturbation.tiff`, `.png`, `.pdf` | CellOracle in silico perturbation priority. |
| Figure 4 | `results/figures/manuscript/main/Figure_4_robustness_support.tiff`, `.png`, `.pdf` | Robustness and supportive evidence. |
| Figure 5 | `results/figures/manuscript/main/Figure_5_program_model.tiff`, `.png`, `.pdf` | Functional interpretation and program-level working model. |
| Fig. S1-S6 | `results/figures/manuscript/supplementary/Supplementary_Figure_S1.*` through `Supplementary_Figure_S6.*` | Supplementary figure outputs matching the manuscript, plus `Supplementary_Figure_Legends.md`. |

The main TIFF files are 170 mm wide (BMC Genomics full-page width). The PNG and PDF files are exports of the same submitted images. Vector SVG copies of the main figures and the retired Supplementary Figure S7 are not part of the 2026-10-25 release.

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
