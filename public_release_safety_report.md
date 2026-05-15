# Public release safety report

## Scope

This report summarizes automated checks performed on `manuscript_output/github_release` during release package preparation.

## Checks performed

- Searched for credential-like patterns: token, password, secret, api_key, github_pat, sk-, OPENAI_API_KEY, and bearer patterns.
- Searched for local absolute path patterns, including Windows user-home paths, Linux home paths, WSL mount paths, and Windows drive-letter project paths.
- Checked for files larger than 50 MB.
- Checked for excluded raw-data/database extensions such as `.h5ad`, `.loom`, `.h5`, `.mtx`, `.feather`, `.gz`, `.tar`, `.oracle`, `.links`, and `.log`.
- Redacted local absolute paths in copied text files where detected.

## Results

- Sensitive credential-like hits requiring review: 0.
- Local absolute path files remaining after redaction: 0.
- Files larger than 50 MB included: 0.
- Forbidden raw-data/database extensions included: 0.

## Notes

The `.gitignore` file intentionally contains patterns such as `*.token`, `secrets*`, and `credentials*`; these are protective ignore rules rather than credentials. These policy terms were not counted as credential hits.

## Manual review required

- Replace GitHub URL placeholder before public release.
- Add Zenodo DOI if the repository is archived.
- Confirm final author/contact information before submission.
