# GitHub upload report

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

## Upload result

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
git commit -m "Initial reproducibility release for FCD II astrocyte TF-regulon analysis"
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
