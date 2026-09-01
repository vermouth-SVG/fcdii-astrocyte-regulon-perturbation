# Final public release review v3

Review date: 2026-09-01

Scope: `manuscript_output/github_release` after synchronization with the BMC Genomics submission package.

## Decision

Local release tree is ready for remote push as tag `v0.3.0-submission`. Remote push and visibility verification are pending because GitHub was unreachable from the current host during this revision.

## Content checks

- Current submission manuscript copied to `Manuscript.docx`.
- Main Figures 1-5 and Supplementary Figures S1-S7 synchronized.
- Main Tables 1-2 and individual Additional-file spreadsheets synchronized.
- Formal CellOracle permutation-test script and fixed-seed outputs added.
- Approx. RI definition corrected across active code, configuration, environment notes, manuscript, and result documentation.
- No file exceeds 50 MB.
- No raw `.h5ad`, `.loom`, `.h5`, `.feather`, `.mtx`, `.rds`, `.pkl`, or `.pickle` object is included.
- Text-file credential scan found no private key, GitHub token, cloud key, password, API key, or secret value.

## Intentional unresolved author inputs

- Final author list and order.
- Complete affiliations.
- Corresponding-author contact details.
- Funding statement.
- CRediT author contributions and all-author approval.
- Optional persistent Zenodo DOI.

These are explicitly marked `AUTHOR_INPUT_NEEDED` in submission materials and are not software-release integrity failures.

## Remaining external actions

1. Push commit and tag `v0.3.0-submission` to the configured GitHub remote when network access is available.
2. Confirm repository visibility is appropriate for peer review and publication.
3. Optionally archive the tagged release in Zenodo and replace the DOI placeholder.
