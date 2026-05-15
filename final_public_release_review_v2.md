# Final public release review v2

Review date: 2026-05-15

Scope: `manuscript_output/github_release` only. No original project files were modified, no remote repository was created, and no GitHub push was performed.

## Decision

yes after GitHub repository creation

## GitHub URL

https://github.com/vermouth-SVG/fcdii-astrocyte-regulon-perturbation

## Files modified

- `README.md`
- `data_availability_statement.md`
- `code_availability_statement.md`
- `CITATION.cff`
- `init_git_and_push_template.sh`
- `init_git_and_push_template.ps1`
- `final_public_release_review.md`
- `final_public_release_review_v2.md`
- `file_manifest.csv`
- `file_manifest.md`

Additional public-package files reviewed or previously downgraded for conservative terminology remain within `manuscript_output/github_release`.

## Placeholder check

No unresolved public-release placeholders were detected.

The requested GitHub URL, repository-name, correspondence, and Zenodo-placeholder patterns were checked recursively and no unresolved public-release placeholders remain.

## Terminology check

No positive claims of formal validation, therapeutic recommendation, validated therapy, or clinical drug candidacy were detected.

Remaining high-risk terms occur only in negated boundary statements, terminology audit files, conservative interpretation notes, or legacy/internal path and column names. The public-facing README and availability statements retain the intended boundaries:

- external datasets are supportive/contextual evidence only;
- Fang 2025 / HRA010445 is published pathway-level contextual evidence only;
- no TF-level validation is claimed;
- CellOracle results are in silico simulated perturbation results only;
- exploratory drug-signature enrichment is hypothesis-generating only;
- no therapeutic recommendation, validated therapy, or clinical candidate claim is made.

## Safety check

- files >50 MB: 0
- forbidden raw-data/database extensions: 0
- credential-like hits requiring manual review: 0
- local absolute path remnants: 0

Credential-scan terminology appears in `.gitignore`, `public_release_safety_report.md`, or review/audit text only as descriptions of the safety checks and was not counted as a real credential leak.

## Manual actions remaining

1. Create the GitHub repository under the selected account.
2. Upload or push the contents of manuscript_output/github_release.
3. Optionally archive the repository in Zenodo and add the DOI after deposition.
4. Re-check the public GitHub page after upload.
