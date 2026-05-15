#!/usr/bin/env python3
from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path


def ensure_dependencies() -> None:
    required = {
        "anndata": "anndata",
        "numpy": "numpy",
        "pandas": "pandas",
    }
    missing: list[str] = []
    for module_name, package_name in required.items():
        try:
            importlib.import_module(module_name)
        except ModuleNotFoundError:
            missing.append(package_name)
    if missing:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet", *sorted(set(missing))])


ensure_dependencies()

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
H5AD_IN = ROOT / "external_validation_round2" / "input" / "external_validation_round2.h5ad"
AUC_CSV = ROOT / "external_validation_round2" / "pyscenic_recalc" / "output" / "auc_mtx.csv"
OUT_DIR = ROOT / "external_validation_round2" / "final"
OUT_H5AD = OUT_DIR / "external_validation_round2_with_auc.h5ad"
OUT_REPORT = OUT_DIR / "merge_report.txt"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(H5AD_IN)
    adata.obs_names = adata.obs_names.astype(str)

    auc = pd.read_csv(AUC_CSV)
    first_col = auc.columns[0]
    auc = auc.rename(columns={first_col: "Cell"}).copy()
    auc["Cell"] = auc["Cell"].astype(str)
    auc = auc.drop_duplicates(subset="Cell").set_index("Cell")

    shared = adata.obs_names.intersection(auc.index)
    missing_in_auc = [x for x in adata.obs_names if x not in auc.index]
    missing_in_h5ad = [x for x in auc.index if x not in set(adata.obs_names)]
    if len(shared) == 0:
        raise ValueError("No shared cells between round2 h5ad and AUC matrix.")

    adata2 = adata[shared].copy()
    auc2 = auc.loc[shared].copy()

    adata2.obsm["X_pyscenic_auc"] = auc2.to_numpy(dtype=np.float32)
    adata2.uns["pyscenic_regulon_names"] = auc2.columns.astype(str).tolist()
    adata2.uns["pyscenic_info"] = {
        "auc_csv": str(AUC_CSV),
        "n_cells": int(auc2.shape[0]),
        "n_regulons": int(auc2.shape[1]),
        "obsm_key": "X_pyscenic_auc",
        "source_dataset": "GSE190452",
    }

    adata2.write_h5ad(OUT_H5AD, compression="gzip")

    lines = [
        "round2 auc merge report",
        f"h5ad_in={H5AD_IN}",
        f"auc_csv={AUC_CSV}",
        f"h5ad_cells={adata.n_obs}",
        f"h5ad_genes={adata.n_vars}",
        f"auc_cells={auc.shape[0]}",
        f"auc_regulons={auc.shape[1]}",
        f"shared_cells={len(shared)}",
        f"missing_in_auc={len(missing_in_auc)}",
        f"missing_in_h5ad={len(missing_in_h5ad)}",
        f"out_h5ad={OUT_H5AD}",
    ]
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote {OUT_H5AD}")
    print(f"Wrote {OUT_REPORT}")


if __name__ == "__main__":
    main()
