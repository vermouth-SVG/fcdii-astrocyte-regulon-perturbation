# Final public release review

Review date: 2026-05-15

Scope: `manuscript_output/github_release` only. No original project files were modified and no GitHub push was performed.

## Recommendation

**superseded by `final_public_release_review_v2.md`**

Repository URL and correspondence placeholders were resolved after this first review. See `final_public_release_review_v2.md` for the current public-release decision. The manuscript-facing README, data availability statement, and code availability statement use conservative language and do not claim formal validation, wet-lab perturbation, therapeutic recommendation, or TF-level validation from supportive datasets.

## High-Risk Issues

None identified in the final reviewed package.

Checks performed:

- No unredacted credential-like hits were found. The only remaining occurrences of credential-scan terminology appear in safety/review reports as descriptions of the scan, not as credentials.
- No unredacted private email addresses or submission-system terms were detected.
- No unredacted Windows user-home paths, local drive-letter project paths, Linux home paths, or WSL mount paths were detected. Redacted placeholders such as `<PROJECT_ROOT_WINDOWS>`, `<PROJECT_ROOT_WSL>`, and `<WINDOWS_USER_HOME>` remain in provenance notes where needed.
- No files larger than 50 MB were detected.
- No `.h5ad`, `.loom`, `.h5`, `.mtx`, `.feather`, `.parquet`, compressed archives, pickle/joblib/numpy binary files, database files, Docker image layers, or cache folders were detected.

## Medium-Risk Issues

- Some legacy folder names and script names retain `external_validation` because they are historical path labels. The README now states that external datasets are supportive/contextual evidence only and are not interpreted as TF-level validation.
- Some exploratory drug-signature file names and internal columns retain historical labels such as `lead_layer_v21`, and SwissADME retains the standard field `Leadlikeness`. These are internal ranking/ADME terms, not therapeutic claims. Manuscript-facing wording and paper-ready drug-signature table values were downgraded to hypothesis-generating mechanism-direction clues.
- Audit and terminology-check files intentionally contain terms such as `external validation`, `confirmed mechanism`, and `validated drug` as negated or required-absent terms. These are not positive claims.

## Wording Downgrades Applied

- `README.md`: added explicit language that external datasets are supportive/contextual evidence only and not TF-level validation; drug-signature enrichment is hypothesis-generating only.
- `docs/README_project_root_sanitized.md`: changed supportive external analysis and drug-signature module wording; replaced pseudobulk validation wording with robustness wording.
- `metadata/gse268807/dataset_registry.csv` and `README_metadata.txt`: replaced validation-role labels with discovery/supportive/contextual labels for public-facing metadata.
- `results/tables/supportive_gse190452/external_round2_regulon_summary_cn.txt`: replaced external-validation and regulon-level-validation wording with cross-syndrome supportive/contextual evidence wording.
- `results/tables/celloracle_ko/next_step_recommendations.txt`: replaced validation-first language with follow-up and consistency-check language.
- Supportive external scripts: plot titles and help text were changed from External Validation to supportive external evidence where they were visible output labels.
- Drug-signature scripts and processed exploratory tables: positive lead/candidate wording was downgraded where it represented interpretation text; results are described as exploratory mechanism-direction clues only.

## Fang 2025 / HRA010445 Review

Passed.

Relevant wording:

- `README.md`: Fang 2025 / HRA010445 is listed as published pathway-level contextual evidence only.
- `data_availability_statement.md`: Fang 2025 / HRA010445 is described only as published pathway-level contextual evidence, and controlled-access raw matrices are not incorporated.
- `metadata/source_data_accessions.md`: Fang 2025 / HRA010445 is listed only as pathway-level contextual evidence.

No TF-level validation claim for Fang 2025 / HRA010445 was found in the manuscript-facing release documents.

## Drug-Signature Review

Passed with the medium-risk legacy-label note above.

The README and figure/table-facing wording describe drug-signature enrichment as exploratory and hypothesis-generating. The package does not present compounds as validated therapy, treatment recommendation, or therapeutic lead. Remaining `lead` strings are restricted to legacy filenames/internal columns or standard SwissADME `Leadlikeness` terminology and should not be used in manuscript prose.

## BMC Data Availability Suitability

Yes.

`data_availability_statement.md` is suitable for BMC Genomics after replacing the GitHub URL placeholder. It states that:

- GSE268807, GSE140393, and GSE190452 are publicly available GEO datasets.
- Fang 2025 / HRA010445 is used only as published pathway-level contextual evidence.
- Large raw matrices, derived large intermediates, Docker image layers, cisTarget ranking databases, and motif annotation resources are not hosted in GitHub.
- Large files should be obtained from GEO, AertsLab/cisTarget, or original providers.

The statement does not say that all data are available in this repository.

## BMC Code Availability Suitability

Yes.

`code_availability_statement.md` is suitable for BMC Genomics after replacing the GitHub URL and optional Zenodo DOI placeholders. It states that:

- The repository contains scripts, configuration/environment notes, processed summaries, figure/table materials, and reproducibility audit files.
- Large raw single-cell matrices, external databases, and bulky intermediate objects are not redistributed.
- The package supports layered reproduction and does not claim a single-command full rerun.

## Large File And Raw Matrix Review

Passed.

No disallowed raw or bulky data extensions were found in `github_release`. The full pseudobulk expression matrix was excluded from the public package; only small processed summaries and the candidate-target pseudobulk matrix remain. The small TF list `configs/pyscenic/allTFs_hg38.txt` is included, but large cisTarget `.feather` ranking databases and motif database resources are not redistributed.

## Final Notes For Submission

- Replace all placeholders before making the repository public.
- Keep the current conservative language in the README and availability statements.
- Do not describe GSE140393, GSE190452, or Fang 2025 / HRA010445 as formal validation.
- Do not describe exploratory drug-signature results as therapeutic candidates, treatment leads, validated therapies, or recommendations.
