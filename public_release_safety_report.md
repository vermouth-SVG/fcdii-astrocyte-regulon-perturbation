# Public release safety report

## Scope

This report summarizes automated checks performed on `manuscript_output/github_release` during synchronization of submission release `v0.3.0-submission` on 2026-09-01.

## Checks performed

- Searched text-bearing files for private-key headers, GitHub tokens, cloud access keys, passwords, API keys, and secret assignments. Embedded base64 image payloads in SVG files were excluded from the text credential scan because they produce non-semantic random-string matches.
- Searched for local absolute path patterns, including Windows user-home paths, Linux home paths, WSL mount paths, and Windows drive-letter project paths.
- Checked for files larger than 50 MB.
- Checked for excluded raw-data/database extensions such as `.h5ad`, `.loom`, `.h5`, `.mtx`, `.feather`, `.rds`, `.pkl`, and `.pickle`.
- Redacted local absolute paths in copied text files where detected.

## Results

- Sensitive credential-like hits requiring review: 0.
- Local absolute project path retained in the historical upload report: 1 intentional provenance entry; no user-home or credential path was added by this release.
- Files larger than 50 MB included: 0.
- Forbidden raw-data/database extensions included: 0.

## Notes

The `.gitignore` file intentionally contains patterns such as `*.token`, `secrets*`, and `credentials*`; these are protective ignore rules rather than credentials. These policy terms were not counted as credential hits. The largest included file is below 7 MB.

## Manual review required

- Push local commit/tag `v0.3.0-submission` when GitHub connectivity is available and confirm repository visibility.
- Add a Zenodo DOI if the repository is archived.
- Confirm final author order, affiliations, corresponding-author details, funding, and CRediT contributions before submission.
