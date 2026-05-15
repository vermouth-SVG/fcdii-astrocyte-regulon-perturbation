#!/usr/bin/env python3
from __future__ import annotations

import importlib
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def decode_strings(values) -> list[str]:
    output = []
    for value in values:
        if isinstance(value, bytes):
            output.append(value.decode("utf-8"))
        else:
            output.append(str(value))
    return output


def check_python_modules(module_names: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for name in module_names:
        try:
            module = importlib.import_module(name)
            version = getattr(module, "__version__", "unknown")
            rows.append({"module": name, "status": "ok", "version": str(version), "error": ""})
        except Exception as exc:
            rows.append({"module": name, "status": "missing", "version": "", "error": repr(exc)})
    return rows


def check_docker_daemon() -> dict[str, str | bool]:
    try:
        result = subprocess.run(
            ["docker", "version"],
            capture_output=True,
            text=True,
            check=False,
            timeout=20,
        )
    except Exception as exc:
        return {
            "docker_cli_found": False,
            "docker_daemon_ready": False,
            "message": repr(exc),
        }

    combined = (result.stdout or "") + "\n" + (result.stderr or "")
    return {
        "docker_cli_found": True,
        "docker_daemon_ready": bool(result.returncode == 0),
        "message": combined.strip(),
    }


def read_categorical_or_vector(group: h5py.Group | h5py.Dataset) -> list[str]:
    if isinstance(group, h5py.Group) and "categories" in group and "codes" in group:
        categories = decode_strings(group["categories"][...])
        codes = np.asarray(group["codes"][...], dtype=int)
        values = [categories[code] if code >= 0 else "NA" for code in codes]
        return values
    if isinstance(group, h5py.Dataset):
        return decode_strings(group[...])
    raise TypeError(f"Unsupported obs object type: {type(group)}")


def summarize_matrix(group: h5py.Group) -> dict[str, object]:
    shape = tuple(int(x) for x in group.attrs["shape"])
    data = np.asarray(group["data"][...])
    return {
        "shape": list(shape),
        "dtype": str(data.dtype),
        "nnz": int(data.size),
        "min_nonzero": float(data.min()) if data.size else 0.0,
        "max_nonzero": float(data.max()) if data.size else 0.0,
        "integer_like": bool(np.allclose(data, np.round(data))) if data.size else True,
    }


def build_obs_overview(obs_group: h5py.Group) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for key in obs_group.keys():
        if key == "_index":
            continue
        values = read_categorical_or_vector(obs_group[key])
        series = pd.Series(values, dtype="string")
        counts = series.value_counts(dropna=False)
        rows.append(
            {
                "column": key,
                "n_unique": int(series.nunique(dropna=False)),
                "n_missing": int((series == "NA").sum()),
                "top_values": " | ".join(f"{idx}:{int(val)}" for idx, val in counts.head(6).items()),
                "looks_like_group_or_cluster": bool(2 <= series.nunique(dropna=False) <= 30),
            }
        )
    return pd.DataFrame(rows).sort_values(by=["looks_like_group_or_cluster", "n_unique", "column"], ascending=[False, True, True]).reset_index(drop=True)


def find_embedding_candidates(obsm_keys: list[str]) -> list[str]:
    candidates = []
    for key in obsm_keys:
        low = key.lower()
        if "umap" in low or "pca" in low or "tsne" in low or "draw_graph" in low:
            candidates.append(key)
    return candidates


def main() -> None:
    root = project_root_from_file(__file__)
    h5ad_path = root / "final_exports" / "astrocyte_pilot_rna_with_pyscenic_auc.h5ad"
    output_dir = root / "celloracle_run" / "checks"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] Reading h5ad structure: {h5ad_path}")
    with h5py.File(h5ad_path, "r") as handle:
        obs_index = decode_strings(handle["obs"]["_index"][...])
        obs_overview = build_obs_overview(handle["obs"])
        x_summary = summarize_matrix(handle["X"])
        raw_summary = summarize_matrix(handle["raw"]["X"]) if "raw" in handle and "X" in handle["raw"] else None
        obsm_keys = list(handle["obsm"].keys())
        uns_keys = list(handle["uns"].keys())
        var_names = decode_strings(handle["var"]["_index"][:])

    print("[2/4] Checking environment prerequisites")
    module_status = check_python_modules(["anndata", "scanpy", "celloracle", "h5py", "pandas", "numpy", "scipy"])
    docker_status = check_docker_daemon()

    print("[3/4] Building report")
    embedding_candidates = find_embedding_candidates(obsm_keys)
    cluster_candidates = obs_overview.loc[obs_overview["looks_like_group_or_cluster"], "column"].tolist()

    missing_items: list[str] = []
    if not x_summary["integer_like"]:
        missing_items.append("X is not integer-like raw counts.")
    if raw_summary is None:
        missing_items.append("raw/X is missing.")
    elif not raw_summary["integer_like"]:
        missing_items.append("raw/X exists but is not integer-like raw counts.")
    if not embedding_candidates:
        missing_items.append("No UMAP/PCA/tSNE-like embedding found in obsm.")
    if not cluster_candidates:
        missing_items.append("No suitable group/cluster-like obs column detected.")

    for row in module_status:
        if row["module"] in {"anndata", "scanpy", "celloracle"} and row["status"] != "ok":
            missing_items.append(f"Python module missing: {row['module']}")

    if not docker_status["docker_daemon_ready"]:
        missing_items.append("Docker daemon is not ready on this machine.")

    report = {
        "h5ad_path": str(h5ad_path),
        "n_cells": len(obs_index),
        "n_genes": len(var_names),
        "counts_in_X": x_summary,
        "counts_in_raw_X": raw_summary,
        "obs_columns": obs_overview["column"].tolist(),
        "obsm_keys": obsm_keys,
        "uns_keys": uns_keys,
        "cluster_candidates": cluster_candidates,
        "preferred_cluster_column": "group" if "group" in cluster_candidates else (cluster_candidates[0] if cluster_candidates else None),
        "embedding_candidates": embedding_candidates,
        "immediately_ready_for_celloracle": bool(raw_summary is not None and raw_summary["integer_like"] and bool(cluster_candidates) and bool(embedding_candidates)),
        "ready_after_minimal_preprocessing": bool(raw_summary is not None and raw_summary["integer_like"] and bool(cluster_candidates)),
        "missing_items": missing_items,
        "python_module_status": module_status,
        "docker_status": docker_status,
    }

    obs_overview.to_csv(output_dir / "obs_overview_for_celloracle.csv", index=False)
    pd.DataFrame(module_status).to_csv(output_dir / "python_module_status.csv", index=False)
    pd.DataFrame({"obsm_key": obsm_keys}).to_csv(output_dir / "obsm_keys.csv", index=False)
    pd.DataFrame({"uns_key": uns_keys}).to_csv(output_dir / "uns_keys.csv", index=False)

    with open(output_dir / "celloracle_input_check_report.json", "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)

    text_lines = [
        f"h5ad: {h5ad_path}",
        f"cells: {report['n_cells']}",
        f"genes: {report['n_genes']}",
        f"preferred_cluster_column: {report['preferred_cluster_column']}",
        f"embedding_candidates: {', '.join(embedding_candidates) if embedding_candidates else 'NONE'}",
        f"ready_after_minimal_preprocessing: {report['ready_after_minimal_preprocessing']}",
        f"immediately_ready_for_celloracle: {report['immediately_ready_for_celloracle']}",
        "missing_items:",
    ]
    if missing_items:
        text_lines.extend(f"- {item}" for item in missing_items)
    else:
        text_lines.append("- NONE")
    with open(output_dir / "celloracle_input_check_report.txt", "w", encoding="utf-8") as handle:
        handle.write("\n".join(text_lines) + "\n")

    print("[4/4] Done")
    print(f"  preferred_cluster_column={report['preferred_cluster_column']}")
    print(f"  embedding_candidates={embedding_candidates if embedding_candidates else ['NONE']}")
    print(f"  ready_after_minimal_preprocessing={report['ready_after_minimal_preprocessing']}")
    print(f"  outputs={output_dir}")


if __name__ == "__main__":
    main()
