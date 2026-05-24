# Regulon and in silico perturbation analysis identifies asymmetric NFE2L2–THRB programs in FCD II astrocytes

Current manuscript title: `Regulon and in silico perturbation analysis identifies asymmetric NFE2L2–THRB programs in FCD II astrocytes`

Chinese title: `Regulon 活性与计算扰动分析识别 FCD II 星形胶质细胞中的非对称 NFE2L2–THRB 程序`

Running title: `NFE2L2–THRB programs in FCD II astrocytes`

## Overview

This repository is a manuscript-level reproducibility package for the associated manuscript. It contains analysis scripts, processed summary tables, figure-generation inputs/outputs, final figure/table outputs, provenance notes, and audit reports.

This repository supports manuscript-level layered reproducibility using processed summary tables, analysis scripts, figure inputs/outputs, and provenance notes. Full upstream reruns require external raw data and third-party resources.

The repository does not redistribute large raw single-cell matrices, controlled-access raw data, cisTarget ranking databases, motif annotation files, Docker image layers, or large intermediate objects.

## Study design

- Discovery object: GSE268807 astrocyte discovery object.
- Discovery design: 4 biological samples, including 2 lesion samples and 2 internal-control samples from 2 donors.
- Discovery dimensions: 2,322 astrocytes and 36,601 genes.
- pySCENIC/AUCell output: 105 final pySCENIC regulons matched to the discovery astrocyte object.
- Supportive datasets: GSE140393 and GSE190452, used only for supportive expression and gene-set analyses.
- Contextual published evidence: Fang 2025 / HRA010445, used only as published contextual evidence. Controlled-access raw matrices from that source were not incorporated into this analysis.

## Main analysis modules

The manuscript methods are organized around the following analysis modules:

1. Data preprocessing.
2. pySCENIC/AUCell regulon inference and activity scoring.
3. Differential regulon analysis.
4. Candidate transcription factor prioritization.
5. CellOracle in silico perturbation simulation.
6. Robustness and sensitivity analyses.
7. Supportive external data analyses.
8. Functional enrichment and convergence analysis.
9. Exploratory drug-signature annotation.
10. Statistical analysis and figure generation.

## Data sources

GSE268807, GSE140393, and GSE190452 were obtained from the Gene Expression Omnibus. GSE268807 provides the astrocyte discovery object. GSE140393 and GSE190452 are used for supportive expression and gene-set analyses only.

Fang 2025 / HRA010445 is cited as published contextual evidence for FCD II tissue context. Controlled-access raw matrices from that source were not included in this repository and were not incorporated into the present analysis.

Large raw data files and third-party databases are not redistributed. Public raw matrices should be obtained from GEO or the original data repositories. cisTarget ranking databases, motif annotation tables, and related pySCENIC resources should be obtained from AertsLab/cisTarget resources according to their access terms and licenses.

## How to reproduce

### Lightweight reproduction

Use the included processed summary tables and figure/table scripts to reproduce manuscript-level figures and tables. This route is intended for manuscript-level inspection and does not require redistributing large raw matrices or cisTarget databases.

Recommended starting points:

- `results/tables/manuscript/`
- `results/figures/manuscript/`
- `scripts/manuscript_figures/`
- `scripts/manuscript_tables/`
- `results/audit/repository_release_audit.md`
- `MANIFEST.md`

### Full upstream rerun

A full upstream rerun requires external downloads and local infrastructure, including public source data, cisTarget ranking databases, motif annotation files, Docker/Python environments, and sufficient disk and memory for single-cell and pySCENIC workflows.

This repository does not claim a one-command, bitwise-identical full rerun of all upstream analyses. It provides scripts, processed summaries, final outputs, and provenance notes to support layered reproducibility.

## Expected outputs

Manuscript-facing outputs include:

- Main figures: Figure 1 through Figure 5.
- Supplementary figures: Supplementary Figure S1 through S7.
- Main tables: Table 1 and Table 2.
- Supplementary Tables: ST01 through ST09.
- pySCENIC/AUCell summaries.
- CellOracle perturbation summaries.
- Robustness and sensitivity outputs.
- Supportive external analysis outputs.
- Functional enrichment and convergence outputs.
- Exploratory drug-signature annotation outputs.

Some retained folders and files use legacy names such as `celloracle_ko`, `ko_round1`, or `round1_ko_*`. These names are retained to avoid breaking script paths and do not define manuscript terminology. Current manuscript-facing terminology is `in silico perturbation`.

## Citation

Ye X. Regulon and in silico perturbation analysis identifies asymmetric NFE2L2–THRB programs in FCD II astrocytes. Manuscript in preparation/submission.

See `CITATION.cff` for repository citation metadata.

## License

Code is released under the MIT License. Data are subject to the terms of the original public repositories and databases. Large third-party databases are not redistributed in this repository.

## Contact

Correspondence information will be provided in the associated manuscript.
