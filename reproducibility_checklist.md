# Reproducibility checklist

- [x] Public data accessions are listed: GSE268807, GSE140393, GSE190452, and Fang 2025 / HRA010445 contextual evidence.
- [x] Active analysis scripts are included.
- [x] Processed summary tables are included.
- [x] Final manuscript figures are included.
- [x] Final main and supplementary tables are included.
- [x] Large raw files are excluded with reasons in `excluded_files_manifest.csv`.
- [x] Third-party ranking and motif database files are excluded or documented with source filenames/URLs.
- [x] pySCENIC/AUCell audit files are included under `results/audit/`.
- [x] CellOracle results are labeled as in silico simulated perturbations.
- [x] Exploratory drug-signature results are labeled as exploratory and not therapeutic recommendations.
- [x] No files larger than 50 MB are included in this release package based on the preparation audit.
- [x] No h5ad, loom, h5, raw mtx, feather, or CellOracle binary object files are included.
- [x] Public safety scan was performed and documented in `public_release_safety_report.md`.
- [x] Manuscript claims are supported by included processed tables and audit files: 4 samples, 2 lesion / 2 internal-control samples, 2 donors, 2,322 astrocytes, 36,601 genes, and 105 pySCENIC regulons.
- [ ] Manual review before public upload: confirm GitHub URL and decide whether Zenodo archival is needed.
- [ ] Manual review before public upload: confirm final author/contact information before submission.
