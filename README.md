# Single-nucleus transcriptomic reanalysis prioritizes NFE2L2 and THRB regulatory programs in focal cortical dysplasia type II astrocytes

## Overview

This repository is the manuscript-level reproducibility package for the BMC Genomics submission of October 2026. It contains the submitted manuscript, the cover letter, the main and supplementary figures, the supplementary tables, the data and code availability statements, processed analysis summaries, analysis and figure-generation scripts, configuration files and provenance notes.

Current submission version: **2026-10-25**. The analysis snapshot in `scripts/` and `results/tables/` is the 2026-09-10 revision and was not changed by the October manuscript update; the manuscript, figures and supplementary files were updated to the submitted versions.

Large raw single-cell matrices, controlled-access data, cisTarget databases, motif annotation files, Docker image layers and large intermediate objects are not redistributed. Full upstream reruns require the original public data and third-party resources.

## Study design and evidence boundary

- Discovery dataset: GSE268807, preprocessed 10x Genomics Multiome astrocyte subset.
- Design: four samples from two paired donors (lesion L1, L2; patient-matched internal control C1, C2).
- Cells: 2,322 FCD II astrocytes, including 888 lesion-tissue and 1,434 internal-control cells, spanning 36,601 genes.
- Regulons: 105 motif-supported pySCENIC regulons. Seventy-six differed between lesion tissue and internal control at FDR < 0.05 (41 higher and 35 lower in lesion tissue).
- Two leading programs: NFE2L2, the leading program increased in lesion tissue, and THRB, the leading homeostatic program reduced in lesion tissue.
- Supporting candidates reported briefly: SATB2, SOX2, SMAD3, SOX10 and CUX1. BHLHE40 was excluded at the entry step.
- GSE140393 and GSE190452 are contextual datasets; they can test the direction of an expression signal but are not FCD II-specific, matched or regulon-level validation datasets.
- Functional enrichment, program convergence and DSigDB drug-signature results are contextual or hypothesis-generating annotations, not mechanism tests or treatment recommendations.
- The internal control is histologically normal, same-hemisphere, same-lobe cortex adjacent to the lesion rather than normal cortex from an unaffected brain.

## Submission files

| File | Content |
|---|---|
| `Manuscript.docx` | Submitted manuscript (BMC Genomics, 2026-10-25 version; no tracked changes, no comments). |
| `Cover_Letter.docx` | Cover letter to the editors. |
| `Additional_file_1.zip` | Supplementary Figures S1-S6 (TIFF, 600 dpi) with legends. |
| `Additional_file_2.zip` | Supplementary Tables ST01-ST09 (XLSX) with legends. |
| `Supplementary_Tables_ST01_ST09_submission_ready.xlsx` | Browsable copy of the ST01-ST09 workbook. |

## Figures

- Main figures 1-5: `results/figures/manuscript/main/` - the submitted TIFF files plus PNG and PDF exports of the same images. The TIFF files are 170 mm wide, the BMC Genomics full-page width.
- Supplementary figures S1-S6: `results/figures/manuscript/supplementary/` - the submitted TIFF files plus PNG and PDF exports, and `Supplementary_Figure_Legends.md`.
- The former Supplementary Figure S7 was retired in the 2026-10 revision and the remaining supplementary figures were renumbered, so S1-S6 is the complete current supplementary set.
- Vector SVG copies of the main figures were removed in the 2026-10 revision because the figures were rebuilt and only raster outputs exist for the submitted version.

## Reproducibility entry points

- Analysis scripts: `scripts/` (pySCENIC/AUCell, differential regulon analysis, candidate prioritisation, CellOracle in silico perturbation, robustness, functional interpretation, drug-signature annotation).
- Figure generation scripts: `scripts/manuscript_figures/`.
- Table generation scripts: `scripts/manuscript_tables/`.
- Seven-candidate CellOracle summaries: `results/tables/celloracle_perturbation_7tf/`.
- Current submission tables: `results/tables/submission_additional_files_current/`.
- Runtime documentation: `configs/` and `environment_notes.md`.
- Revision boundary and legacy notes: `REVISION_NOTES_20260902.md`.

## Original computational environment

- pySCENIC: Docker image `aertslab/pyscenic:0.12.1`.
- CellOracle: Docker image `kenjikamimoto126/celloracle_ubuntu:0.18.0`.
- Fixed CellOracle seed: `20240817`, reset for each candidate.

## Version history

- **2026-10-25** - Manuscript, cover letter, main figures (170 mm TIFF plus PNG/PDF exports) and supplementary files updated to the BMC Genomics submission version; supplementary figures renumbered to S1-S6; vector SVG figure copies removed; README, MANIFEST, CITATION and availability statements updated.
- **2026-09-10** - v0.4.0 analysis snapshot: seven-candidate CellOracle rerun, processed summaries and figure scripts. Preserved unchanged under `scripts/` and `results/`.
- **2026-09-01** - v0.3.0 submission release (archived in the git history).

## Citation

Ye X, Wang W, Li Y, Hu F, Hu L, Ruan Z, Lu Y, Li X. Single-nucleus transcriptomic reanalysis prioritizes NFE2L2 and THRB regulatory programs in focal cortical dysplasia type II astrocytes. Manuscript submitted to BMC Genomics, 2026. See `CITATION.cff` for repository citation metadata.

## License

Code is released under the MIT License. Data remain subject to the terms of their original repositories and databases.

## Contact

Correspondence: see the corresponding author details in `Manuscript.docx`.
