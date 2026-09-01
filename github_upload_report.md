# GitHub upload report

## Current submission-release status (2026-09-01)

The current BMC Genomics submission package has been synchronized to the local
release tree and prepared for the annotated tag `v0.3.0-submission`. The remote
repository could not be reached from this host during this revision, so the new
commit/tag has not yet been pushed and repository visibility has not yet been
verified. The successful upload recorded below refers to the earlier repository
state, not to the 2026-09-01 submission release.

## Target repository

https://github.com/vermouth-SVG/fcdii-astrocyte-regulon-perturbation

## Local source directory

D:\docker_run_pyscenic\manuscript_output\github_release

## Tool used

Windows Git 2.54.0.windows.1

## Safety check before upload

- unresolved placeholders: 0
- wrong owner remnants: 0
- files >50 MB: 0
- forbidden raw-data/database files: 0
- credential-like hits requiring manual review: 0
- local absolute path remnants: 0

Credential-scan terms were detected only in `.gitignore`, `public_release_safety_report.md`, or safety-scan report text and were not counted as credential leaks.

## Earlier upload result

success

## Repository visibility

unknown

Visibility was not changed by this upload task. Please confirm in the GitHub web interface that the repository remains private.

## Commands executed

```text
Get-Location
git --version
git status
git init
git config user.name "vermouth-SVG"
git config user.email "vermouth-SVG@users.noreply.github.com"
git branch -M main
git remote remove origin
git remote add origin https://github.com/vermouth-SVG/fcdii-astrocyte-regulon-perturbation.git
git remote -v
git add .
git commit -m "Initial reproducibility release for FCD II astrocyte regulon perturbation analysis"
git push -u origin main
git remote -v
git status
```

No token or authentication secret was written or displayed.

## Next manual checks

- Open GitHub repository page.
- Confirm repository visibility is private.
- Confirm README renders correctly.
- Confirm expected files and folders are present.
- Send repository link back for final manuscript Data/Code availability check.
