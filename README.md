# FCD II astrocyte TF-regulon and in silico perturbation analysis

## Overview

This repository is a release-ready reproducibility package for the study "Transcription factor regulatory programs and in silico knockout prioritization in FCD II astrocytes using public single-cell data". The project is a computational discovery and candidate prioritization study based on public single-cell datasets. It provides analysis scripts, processed summary tables, figure-generation materials, final figure/table outputs, and provenance notes supporting the manuscript.

The repository does not contain large raw single-cell matrices, controlled-access raw data, cisTarget ranking databases, motif databases, Docker image layers, or large intermediate objects.

## Study design

Discovery object: GSE268807 astrocyte pilot.

Fixed discovery cohort summary:

- 4 samples.
- 2 lesion and 2 internal-control samples.
- 2 donors.
- 2,322 astrocytes.
- 36,601 genes.
- 105 pySCENIC regulons.

Supportive and contextual evidence:

- GSE140393.
- GSE190452.
- Fang 2025 / HRA010445 as published pathway-level contextual evidence only.

Main analysis modules:

- pySCENIC/AUCell regulon activity analysis.
- Differential regulon analysis.
- TF prioritization.
- CellOracle in silico perturbation analysis.
- Leave-one-sample-out and pseudobulk robustness analyses.
- Supportive/contextual evidence from public datasets and published pathway-level context.
- GO/KEGG functional interpretation.
- Exploratory drug-signature enrichment.

Interpretation limits:

- This repository supports computational discovery and candidate prioritization.
- The analyses are not formal validation.
- External datasets are used as supportive/contextual evidence only and are not interpreted as TF-level validation.
- The CellOracle results are not wet-lab perturbation results.
- The exploratory drug-signature analysis provides hypothesis-generating computational clues only and is not a therapeutic recommendation.

## Repository structure

- `scripts/`: analysis, post-processing, robustness, functional interpretation, exploratory drug-signature, and figure/table generation scripts.
- `configs/`: small configuration files and resource notes, including pySCENIC resource filenames and CellOracle configuration summaries.
- `metadata/`: public sample-level metadata and dataset registry files needed to interpret the processed outputs.
- `results/tables/`: processed summary tables used for manuscript figures, tables, and supplementary results.
- `results/figures/`: final and intermediate processed figure outputs used for manuscript support.
- `results/audit/`: reproducibility audits, redacted logs, merge reports, structure checks, and provenance notes.
- `docs/`: data/code availability statements, environment notes, release checklist, safety report, and other documentation.

## Data sources

This study used public data from the following sources:

- GSE268807: discovery astrocyte pilot source dataset.
- GSE140393: supportive/contextual single-cell evidence.
- GSE190452: supportive/contextual single-cell evidence.
- Fang 2025 / HRA010445: used only as published pathway-level contextual evidence; controlled-access raw matrices were not included in this study.

Large raw data files and third-party databases are not redistributed in this GitHub repository. Raw public matrices should be obtained from GEO or the original data repositories. cisTarget ranking databases, motif annotation tables, and related pySCENIC resources should be obtained from AertsLab/cisTarget resources using the filenames and notes in `environment_notes.md` and `configs/pyscenic_runtime_and_resources.md`.

## Key reproducibility notes

- pySCENIC was run in Docker using image `aertslab/pyscenic:0.12.1`.
- GRN inference used GRNBoost2 by the pySCENIC 0.12.1 CLI default; the recovered command did not explicitly pass `--method`.
- TF list: AertsLab cisTarget `allTFs_hg38.txt`, 1,892 TF entries.
- Ranking database: `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`.
- Motif annotation table: `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`.
- AUCell output was matched back to 2,322 astrocytes with no missing or extra cells.
- The final downstream comparison used 105 pySCENIC regulons.
- CellOracle outputs are in silico simulated perturbations and should be interpreted as candidate-prioritization evidence only.

## How to reproduce

### 1. Lightweight reproduction

Use the included processed summary tables and figure-generation scripts to inspect or regenerate manuscript-level tables and figures. This route does not require redistributing large raw matrices or cisTarget databases.

Recommended starting points:

- `results/tables/manuscript/` for final main and supplementary tables.
- `results/figures/manuscript/` for final main and supplementary figures.
- `scripts/manuscript_tables/` and `scripts/manuscript_figures/` for table and figure generation code.
- `results/audit/` for provenance and pySCENIC/AUCell reproducibility notes.

### 2. Full computational reproduction

A full rerun requires downloading public source data and large third-party resources, including cisTarget ranking databases and motif annotation files. It also requires a local Docker/Python/R environment and sufficient memory and disk space for single-cell and pySCENIC workflows. See `environment_notes.md` for resource requirements and known limitations.

This repository does not promise a single-command, bitwise-identical rerun of all analyses. It provides the scripts, processed summaries, and provenance needed to reproduce or audit the computational workflow in layers.

## Expected outputs

The release package contains the following manuscript-supporting outputs:

- Main figures: Figure 1 through Figure 5.
- Supplementary figures: Supplementary Figure S1 through S7.
- Main tables: Main Table 1 and Main Table 2.
- Supplementary Tables ST01-ST09.
- Key processed result tables for pySCENIC/AUCell, differential regulon analysis, TF prioritization, CellOracle KO summaries, robustness analyses, external supportive evidence, functional interpretation, and exploratory drug-signature analysis.

## Citation

If you use this repository, please cite the associated manuscript:

Ye X. Transcription factor regulatory programs and in silico knockout prioritization in FCD II astrocytes using public single-cell data. Manuscript in preparation/submission.

See `CITATION.cff` for repository citation metadata.

## License

Code is released under the MIT License. Data are subject to the terms of the original public repositories and databases. Large third-party databases are not redistributed in this repository.

## Contact

Correspondence: Correspondence information will be provided in the associated manuscript.
