#!/usr/bin/env bash
set -euo pipefail

git init
git add .
git commit -m "Initial reproducibility release for FCD II astrocyte regulon perturbation analysis"
git branch -M main
git remote add origin https://github.com/vermouth-SVG/fcdii-astrocyte-regulon-perturbation.git
git push -u origin main
