#!/usr/bin/env python3
from __future__ import annotations

import importlib
import json
import math
import subprocess
import sys
from pathlib import Path


def ensure_dependencies() -> None:
    required = {
        "anndata": "anndata",
        "numpy": "numpy",
        "pandas": "pandas",
        "scipy": "scipy",
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
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
INPUT_H5AD = ROOT / "external_validation_round2" / "input" / "external_validation_round2.h5ad"
OUT_DIR = ROOT / "external_validation_round2" / "pyscenic_recalc"
TXT_PATH = OUT_DIR / "readiness_report.txt"
JSON_PATH = OUT_DIR / "readiness_report.json"


def inspect_matrix(matrix, sample_n: int = 100000) -> dict[str, object]:
    info: dict[str, object] = {
        "type": type(matrix).__name__,
        "dtype": str(getattr(matrix, "dtype", "")),
        "is_sparse": bool(sparse.issparse(matrix)),
        "shape": tuple(int(x) for x in matrix.shape),
        "nonzero": None,
        "min_nonzero": None,
        "max_nonzero": None,
        "sample_integerish": None,
    }
    if sparse.issparse(matrix):
        x = matrix.tocsr()
        data = x.data
        info["nonzero"] = int(data.size)
        if data.size:
            info["min_nonzero"] = float(data.min())
            info["max_nonzero"] = float(data.max())
            probe = data[: min(sample_n, data.size)]
            info["sample_integerish"] = bool(np.allclose(probe, np.round(probe)))
        else:
            info["min_nonzero"] = 0.0
            info["max_nonzero"] = 0.0
            info["sample_integerish"] = True
        return info

    arr = np.asarray(matrix)
    flat = arr.ravel()
    nonzero = flat[flat != 0]
    info["nonzero"] = int(nonzero.size)
    if nonzero.size:
        info["min_nonzero"] = float(nonzero.min())
        info["max_nonzero"] = float(nonzero.max())
        probe = nonzero[: min(sample_n, nonzero.size)]
        info["sample_integerish"] = bool(np.allclose(probe, np.round(probe)))
    else:
        info["min_nonzero"] = 0.0
        info["max_nonzero"] = 0.0
        info["sample_integerish"] = True
    return info


def classify_representation(info: dict[str, object]) -> str:
    integerish = bool(info.get("sample_integerish"))
    max_nonzero = float(info.get("max_nonzero") or 0.0)
    dtype = str(info.get("dtype") or "").lower()
    if integerish and ("int" in dtype or max_nonzero >= 20):
        return "raw_counts_like"
    if not integerish and max_nonzero <= 25.0:
        return "normalized_or_log_like"
    if integerish:
        return "counts_like"
    return "unknown"


def resolve_counts_source(adata: ad.AnnData) -> tuple[str | None, dict[str, object] | None]:
    if "counts" in adata.layers:
        info = inspect_matrix(adata.layers["counts"])
        return "layers['counts']", info
    if adata.raw is not None:
        info = inspect_matrix(adata.raw.X)
        if classify_representation(info).endswith("counts_like"):
            return "raw.X", info
    x_info = inspect_matrix(adata.X)
    if classify_representation(x_info).endswith("counts_like"):
        return "X", x_info
    return None, None


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(INPUT_H5AD)

    x_info = inspect_matrix(adata.X)
    x_repr = classify_representation(x_info)
    counts_source, counts_info = resolve_counts_source(adata)

    group_cols = {col: int(adata.obs[col].astype("string").nunique()) for col in ["group", "sample", "donor", "cluster", "celltype"] if col in adata.obs.columns}
    group_counts = {}
    for col in group_cols:
        group_counts[col] = adata.obs[col].astype("string").value_counts(dropna=False).to_dict()

    readiness_issues: list[str] = []
    if counts_source is None:
        readiness_issues.append("未检测到可靠的原始 counts 来源；pySCENIC 不建议直接使用当前 X。")
    if adata.n_obs < 100:
        readiness_issues.append("细胞数过少。")
    if adata.n_vars < 1000:
        readiness_issues.append("基因数过少。")
    if "group" not in adata.obs.columns:
        readiness_issues.append("缺少 group 列。")
    if "sample" not in adata.obs.columns and "donor" not in adata.obs.columns:
        readiness_issues.append("缺少 sample/donor 列。")

    suitable = counts_source is not None and adata.n_obs >= 100 and adata.n_vars >= 1000

    payload = {
        "input_h5ad": str(INPUT_H5AD),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "obsm_keys": list(adata.obsm.keys()),
        "uns_keys": sorted(list(adata.uns.keys())),
        "layers": list(adata.layers.keys()),
        "has_raw": bool(adata.raw is not None),
        "x_info": x_info,
        "x_representation": x_repr,
        "counts_source": counts_source,
        "counts_info": counts_info,
        "metadata_columns": group_cols,
        "metadata_counts": group_counts,
        "suitable_for_pyscenic": bool(suitable),
        "issues": readiness_issues,
        "recommendation": (
            "当前对象可用于 round2 支持性队列 pySCENIC 复算。建议从 layers['counts'] 导出表达矩阵，并保持 group/sample/donor/cluster 信息用于后续 projected regulon-level supportive analysis。"
            if suitable
            else "当前对象不建议直接进入 pySCENIC，请先处理 issues 字段中的问题。"
        ),
    }

    lines = [
        "GSE190452 round2 pySCENIC readiness report",
        f"input_h5ad: {INPUT_H5AD}",
        f"cells: {adata.n_obs}",
        f"genes: {adata.n_vars}",
        f"x_representation: {x_repr}",
        f"counts_source: {counts_source or 'not_found'}",
        f"suitable_for_pyscenic: {suitable}",
        "",
        "[metadata_columns]",
    ]
    lines.extend(f"{key}\t{value}" for key, value in group_cols.items())
    lines.extend(["", "[obsm_keys]"])
    lines.extend(list(adata.obsm.keys()) or ["<none>"])
    lines.extend(["", "[uns_keys]"])
    lines.extend(sorted(list(adata.uns.keys())) or ["<none>"])
    lines.extend(["", "[layers]"])
    lines.extend(list(adata.layers.keys()) or ["<none>"])
    lines.extend(["", "[x_info]"])
    lines.extend(f"{k}: {v}" for k, v in x_info.items())
    if counts_info is not None:
        lines.extend(["", "[counts_info]"])
        lines.extend(f"{k}: {v}" for k, v in counts_info.items())
    if readiness_issues:
        lines.extend(["", "[issues]"])
        lines.extend(readiness_issues)
    lines.extend(["", "[recommendation]", payload["recommendation"]])

    TXT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    JSON_PATH.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {TXT_PATH}")
    print(f"Wrote {JSON_PATH}")
    print(f"suitable_for_pyscenic={suitable}")


if __name__ == "__main__":
    main()
