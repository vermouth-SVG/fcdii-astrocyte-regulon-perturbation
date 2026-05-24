# Repository release audit

Audit date: 2026-05-24

## 1. Manuscript version audited

- File: `Manuscript.docx`
- File timestamp observed in repository: 2026-05-24
- English title: `Regulon and in silico perturbation analysis identifies asymmetric NFE2L2–THRB programs in FCD II astrocytes`
- Chinese title: `Regulon 活性与计算扰动分析识别 FCD II 星形胶质细胞中的非对称 NFE2L2–THRB 程序`
- Running title: `NFE2L2–THRB programs in FCD II astrocytes`

## 2. Repository status

- Current status: private repository for pre-submission review.
- This private status is acceptable for the current pre-submission stage.
- Intended public release: before submission or upon publication, depending on author decision.
- No Zenodo DOI is required for the current review package. Zenodo archival can remain a future submission checklist item.

## 3. Manuscript constants extracted

- Discovery object: GSE268807 astrocyte discovery object.
- Biological samples: 4 total, including 2 lesion samples and 2 internal-control samples.
- Donors: 2.
- Cells: 2,322 astrocytes.
- Genes: 36,601 genes.
- Final pySCENIC regulons: 105.
- Supportive datasets: GSE140393 and GSE190452 for supportive expression / gene-set analyses only.
- Contextual published evidence: Fang 2025 / HRA010445 only; controlled-access raw matrices from that source were not incorporated.

## 4. Abstract extraction

- Research object: FCD II astrocytes from the GSE268807 astrocyte discovery object.
- Discovery object size stated in the Abstract: 2,322 cells and 36,601 genes.
- Supportive datasets stated in the Abstract: GSE140393 and GSE190452 for supportive expression and gene-set analyses.
- Regulon-level result stated in the manuscript Results: 105 pySCENIC regulons.
- Main candidates stated in the Abstract: lesion-associated NFE2L2, internal-control-associated THRB, secondary BHLHE40, and lower-priority target-limited SOX2.
- CellOracle metrics stated in the Abstract: NFE2L2 has the highest mean shift (0.699), and THRB has the highest Approx. Recovery Index (0.689).

## 5. Methods extraction

| Module | Manuscript-matched details |
|---|---|
| Data preprocessing | GSE268807 preprocessed astrocyte object; 4 biological samples, 2 lesion and 2 internal-control samples, 2 donors, 2,322 astrocytes, 36,601 genes. |
| pySCENIC/AUCell | Docker image `aertslab/pyscenic:0.12.1`; GRN inference by GRNBoost2; 4 workers; AUCell matched to 2,322 astrocytes; 105 final pySCENIC regulons. |
| cisTarget / motif resources | Retained repository run uses `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` and `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`; see inconsistency note below because manuscript text currently states v80/v9. |
| Differential regulon analysis | Lesion versus internal-control AUCell comparison; Delta AUC = mean AUC lesion minus mean AUC internal-control; two-sided Wilcoxon rank-sum tests; Benjamini-Hochberg FDR; significance threshold FDR < 0.05. |
| Candidate TF prioritization | Integrated TF expression direction, regulon activity direction, positive-cell proportion, FDR, target count, CellOracle metrics, and robustness metrics; no additional manual weights. |
| CellOracle in silico perturbation | Docker image `kenjikamimoto126/celloracle_ubuntu:0.18.0`; CellOracle 0.18.0; Python 3.10.11 recovered from container; inputs include matched astrocyte expression object, candidate TFs, metadata, low-dimensional embedding, and inferred regulatory network. |
| CellOracle metrics | Perturbation vectors, mean shift, Approx. RI, and random perturbation controls. Approx. RI is study-defined cosine alignment to the lesion-to-internal-control reference vector, not a CellOracle official metric or clinical recovery metric. |
| Robustness / sensitivity | Leave-one-sample-out, pseudobulk sample x group summaries, and threshold sensitivity analyses. |
| Supportive external analyses | GSE140393 and GSE190452 used only for supportive expression / gene-set analyses, not formal validation or independent regulon-level replication. |
| Functional enrichment | Enrichr GO/KEGG, Benjamini-Hochberg FDR, significance threshold FDR < 0.05. |
| Exploratory drug-signature annotation | Enrichr DSigDB and SwissADME BBB annotation; hypothesis-generating only. |

## 6. Availability statements extracted

Data availability in the manuscript states that GSE268807, GSE140393, and GSE190452 are public GEO datasets; Fang 2025 / HRA010445 is used only as published contextual evidence; controlled-access raw matrices from that source are not incorporated; processed summary tables, figure-generation inputs, and reproducibility materials are provided in the GitHub repository; large raw matrices, large intermediate objects, cisTarget ranking databases, motif annotation files, and third-party resources are not redistributed.

Code availability in the manuscript states that the GitHub repository contains analysis scripts, configuration files, processed outputs, figure/table generation scripts, and reproducibility notes.

## 7. Figure and table references

Main manuscript references:

- Figure 1: Study design and analytical framework.
- Figure 2: Differential regulon landscape in the astrocyte pilot.
- Figure 3: CellOracle in silico perturbation priority.
- Figure 4: Robustness and supportive evidence.
- Figure 5: Functional interpretation and program-level working model.
- Table 1: Analytical datasets and object roles.
- Table 2: Integrated prioritization of candidate transcription factors.

Supplementary references:

- Fig. S1: Discovery object provenance and AUCell-matched input overview.
- Fig. S2: Full differential regulon landscape in the GSE268807 astrocyte pilot.
- Fig. S3: Complete CellOracle in silico perturbation metrics and randomized-control assessment.
- Fig. S4: Complete sample-level robustness and sensitivity analyses.
- Fig. S5: Supportive/contextual evidence from GSE140393 and GSE190452.
- Fig. S6: Full functional enrichment landscape and program convergence.
- Fig. S7: Exploratory drug-signature enrichment and BBB-layered annotation.
- Table ST01-ST09: dataset roles, TF shortlist, TF priority metrics, differential regulons, CellOracle perturbation metrics, robustness, supportive evidence, functional convergence, and exploratory drug-signature annotation.

## 8. Consistency checks

| Check | Status | Notes |
|---|---|---|
| Title consistency | Fixed in repository docs | `README.md`, `CITATION.cff`, and manifest files now use the current manuscript title. |
| Terminology consistency | Fixed in active docs | Active manuscript-facing docs now use `in silico perturbation`. Legacy paths containing `ko` are retained for compatibility and recorded as legacy filenames. |
| Dataset counts | Matched | Repository docs now match 4 samples, 2 lesion samples, 2 internal-control samples, 2 donors, 2,322 astrocytes, 36,601 genes, and 105 regulons. |
| Figure/table references | Matched | Manifest and README list Figure 1-5, Fig. S1-S7, Table 1-2, and ST01-ST09. |
| pySCENIC runtime | Matched to retained run | Docker image `aertslab/pyscenic:0.12.1`, GRNBoost2, 4 workers, 105 final regulons. |
| cisTarget resources | Matched | Manuscript resource sentence has been updated to the retained v10 resources. |
| CellOracle environment | Matched | Docker image `kenjikamimoto126/celloracle_ubuntu:0.18.0`, CellOracle 0.18.0, Python 3.10.11 recovered from container. |
| Approx. RI interpretation | Fixed in docs | Documented as a study-defined cosine-similarity-based alignment metric, not a CellOracle official metric or clinical recovery metric. |
| Supplementary table/figure existence | Present | Main and supplementary figure/table outputs are present in manuscript output folders and release package copies. |
| Data availability statement | Matched | Large raw matrices, controlled-access matrices, cisTarget resources, motif annotation files, and Docker layers are not redistributed. |
| Code availability statement | Matched | Repository supports layered reproducibility with scripts, processed outputs, figure/table generation materials, and audit notes. |

## 9. Known limitations

- Large raw single-cell matrices are not redistributed.
- Controlled-access data are not redistributed.
- cisTarget ranking databases are not redistributed.
- Motif annotation files are not redistributed.
- Docker image layers are not redistributed.
- Large intermediate objects are not redistributed.
- Full upstream reruns require external downloads and third-party resources.
- pySCENIC stochastic and multiprocessing steps may not be bitwise identical because no fixed seed was recovered for the retained discovery GRN or AUCell commands.
- CellOracle outputs are computational in silico perturbation results and should not be interpreted as experimental perturbation evidence.

## 10. Required manuscript changes before submission

- [x] pySCENIC Methods resource sentence updated in `Manuscript.docx` to match the retained run: cisTarget hg38 gene-based ranking database with 10 kb upstream/downstream of transcription start sites and v10 motif annotation.
- [x] No unresolved manuscript-repository inconsistency remains in the audited manuscript package.

## 11. Future submission checklist

- [ ] Decide whether to make the GitHub repository public before submission or upon publication.
- [ ] Optionally archive the release package on Zenodo after public release and add the DOI to the manuscript and `CITATION.cff`.
- [ ] Confirm final author list, affiliations, funding, and corresponding author details in the manuscript and citation metadata.
