#!/usr/bin/env python3
from __future__ import annotations

import gzip
import importlib
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


def ensure_dependencies() -> None:
    required = {
        "anndata": "anndata",
        "scanpy": "scanpy",
        "igraph": "igraph",
        "leidenalg": "leidenalg",
        "pandas": "pandas",
        "numpy": "numpy",
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
import pandas as pd
import scanpy as sc
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_DIR = ROOT / "external_validation_round2" / "source_gse190452"
DEFAULT_OUTPUT_H5AD = ROOT / "external_validation_round2" / "input" / "external_validation_round2.h5ad"
DEFAULT_REPORT_TXT = ROOT / "external_validation_round2" / "input" / "build_report.txt"
DEFAULT_REPORT_JSON = ROOT / "external_validation_round2" / "input" / "build_report.json"


@dataclass
class SampleBundle:
    sample_id: str
    matrix_name: str
    feature_name: str
    barcode_name: str


def parse_args() -> dict[str, object]:
    return {
        "source_dir": DEFAULT_SOURCE_DIR,
        "output_h5ad": DEFAULT_OUTPUT_H5AD,
        "report_txt": DEFAULT_REPORT_TXT,
        "report_json": DEFAULT_REPORT_JSON,
        "hvg_n_top": 3000,
        "n_pcs": 30,
        "n_neighbors": 15,
        "leiden_resolution": 0.6,
        "target_sum": 1e4,
    }


def make_unique_index(values: pd.Series | list[str], fallback_values: pd.Series | list[str] | None = None) -> pd.Index:
    raw = ["" if pd.isna(v) else str(v).strip() for v in values]
    fallback = None
    if fallback_values is not None:
        fallback = ["" if pd.isna(v) else str(v).strip() for v in fallback_values]
    seen: dict[str, int] = {}
    out: list[str] = []
    for i, value in enumerate(raw):
        base = value or (fallback[i] if fallback is not None and fallback[i] else f"feature_{i+1}")
        n = seen.get(base, 0)
        out.append(base if n == 0 else f"{base}_{n}")
        seen[base] = n + 1
    return pd.Index(out)


def parse_soft_metadata(soft_path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, object] | None = None
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(finalize_soft_record(current))
                current = {
                    "sample_id": line.split("=", 1)[1].strip(),
                    "title": "",
                    "characteristics": [],
                }
                continue
            if current is None:
                continue
            if line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                current["characteristics"].append(line.split("=", 1)[1].strip())
    if current is not None:
        records.append(finalize_soft_record(current))
    return pd.DataFrame(records)


def finalize_soft_record(record: dict[str, object]) -> dict[str, str]:
    char_map: dict[str, str] = {}
    for item in record.get("characteristics", []):
        if ":" in item:
            key, value = item.split(":", 1)
            char_map[key.strip().lower()] = value.strip()

    diagnosis = char_map.get("diagnosis", "") or char_map.get("treatment", "")
    lowered = diagnosis.lower()
    group = "unknown"
    disease = diagnosis or "unknown"
    if "control" in lowered or "non-epileptic" in lowered or "untreated" in lowered or "normal" in lowered:
        group = "internal_control"
        disease = "control"
    elif "tle" in lowered or "epilep" in lowered:
        group = "lesion"
        disease = "TLE"

    title = str(record.get("title", ""))
    sample_name = title.replace("single-cell_", "").strip() or str(record.get("sample_id", ""))
    return {
        "sample_id": str(record.get("sample_id", "")),
        "title": title,
        "diagnosis_raw": diagnosis,
        "group": group,
        "disease": disease,
        "sample": sample_name,
        "donor": sample_name,
        "tissue": char_map.get("tissue", ""),
        "molecule_subtype": char_map.get("molecule subtype", ""),
    }


def discover_sample_bundles(extract_dir: Path, metadata_df: pd.DataFrame) -> list[SampleBundle]:
    names = {path.name for path in extract_dir.iterdir() if path.is_file()}
    bundles: list[SampleBundle] = []
    for sample_id in metadata_df["sample_id"].tolist():
        matrix = next((name for name in names if name.startswith(f"{sample_id}_") and name.endswith("matrix.mtx.gz")), None)
        features = next((name for name in names if name.startswith(f"{sample_id}_") and name.endswith("features.tsv.gz")), None)
        barcodes = next((name for name in names if name.startswith(f"{sample_id}_") and name.endswith("barcodes.tsv.gz")), None)
        if matrix and features and barcodes:
            bundles.append(SampleBundle(sample_id, matrix, features, barcodes))
    return bundles


def read_features(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    if df.shape[1] >= 3:
        out = df.iloc[:, :3].copy()
        out.columns = ["gene_id", "gene_symbol", "feature_type"]
    else:
        out = df.iloc[:, :2].copy()
        out.columns = ["gene_id", "gene_symbol"]
        out["feature_type"] = "Gene Expression"
    return out


def load_sample(
    extract_dir: Path,
    bundle: SampleBundle,
    sample_meta: dict[str, object],
    reference_gene_ids: list[str] | None,
) -> tuple[ad.AnnData, list[str]]:
    features = read_features(extract_dir / bundle.feature_name)
    barcodes = pd.read_csv(
        extract_dir / bundle.barcode_name,
        sep="\t",
        header=None,
        compression="gzip",
        names=["cell_barcode"],
    )
    matrix = mmread(extract_dir / bundle.matrix_name).tocsr().transpose().tocsr()
    gene_ids = features["gene_id"].astype(str).tolist()
    if reference_gene_ids is not None and gene_ids != reference_gene_ids:
        raise ValueError(f"Feature order mismatch in {bundle.sample_id}")

    obs = pd.DataFrame(index=[f"{bundle.sample_id}__{bc}" for bc in barcodes["cell_barcode"].astype(str).tolist()])
    obs["cell_barcode"] = barcodes["cell_barcode"].astype(str).tolist()
    for key, value in sample_meta.items():
        obs[key] = value
    obs["geo_series"] = "GSE190452"
    obs["source_dataset"] = "GSE190452"

    var = features.copy()
    var.index = pd.Index(var["gene_id"].astype(str), name="gene_id")

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    return adata, gene_ids


def build_counts_h5ad(source_dir: Path) -> tuple[ad.AnnData, pd.DataFrame]:
    soft_path = source_dir / "GSE190452_family.soft.gz"
    extract_dir = source_dir / "extracted"
    meta = parse_soft_metadata(soft_path)
    if meta.empty:
        raise ValueError(f"No sample metadata parsed from {soft_path}")
    meta = meta.set_index("sample_id", drop=False)
    bundles = discover_sample_bundles(extract_dir, meta.reset_index(drop=True))
    if not bundles:
        raise ValueError(f"No sample bundles found in {extract_dir}")

    adatas: list[ad.AnnData] = []
    ref_gene_ids: list[str] | None = None
    for bundle in bundles:
        adata_sample, ref_gene_ids = load_sample(
            extract_dir=extract_dir,
            bundle=bundle,
            sample_meta=meta.loc[bundle.sample_id].to_dict(),
            reference_gene_ids=ref_gene_ids,
        )
        adatas.append(adata_sample)

    adata = ad.concat(adatas, axis=0, join="outer", merge="first", label="input_batch", fill_value=0)
    adata.obs_names_make_unique()
    adata.var["gene_id"] = adata.var.index.astype(str)
    gene_symbols = adata.var["gene_symbol"] if "gene_symbol" in adata.var.columns else pd.Series(adata.var.index.astype(str), index=adata.var.index)
    adata.var_names = make_unique_index(gene_symbols, adata.var["gene_id"])
    adata.var_names.name = "gene_symbol_unique"
    return adata, meta.reset_index(drop=True)


def run_minimal_scanpy_processing(adata: ad.AnnData, hvg_n_top: int, n_pcs: int, n_neighbors: int, leiden_resolution: float, target_sum: float) -> ad.AnnData:
    sc.settings.verbosity = 0

    adata.layers["counts"] = adata.X.copy()
    sc.pp.filter_genes(adata, min_counts=1)
    adata.obs["n_counts"] = np.asarray(adata.layers["counts"].sum(axis=1)).ravel().astype(float)
    adata.obs["n_genes_by_counts"] = np.asarray((adata.layers["counts"] > 0).sum(axis=1)).ravel().astype(float)

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    cluster_adata = adata.copy()
    sc.pp.highly_variable_genes(cluster_adata, flavor="seurat", n_top_genes=hvg_n_top, subset=True)
    sc.pp.pca(cluster_adata, n_comps=n_pcs)
    sc.pp.neighbors(cluster_adata, n_neighbors=n_neighbors, n_pcs=min(n_pcs, cluster_adata.obsm["X_pca"].shape[1]))
    sc.tl.umap(cluster_adata)
    sc.tl.leiden(cluster_adata, resolution=leiden_resolution, key_added="leiden")

    adata.obs["cluster"] = "cluster_" + cluster_adata.obs["leiden"].astype(str).values
    adata.obs["celltype"] = adata.obs["cluster"].astype(str)
    adata.obs["cluster_source"] = "generated_leiden"
    adata.obsm["X_pca"] = cluster_adata.obsm["X_pca"].copy()
    adata.obsm["X_umap"] = cluster_adata.obsm["X_umap"].copy()
    adata.uns["cluster_generation"] = {
        "method": "scanpy_minimal_pipeline",
        "hvg_n_top": int(hvg_n_top),
        "n_pcs": int(n_pcs),
        "n_neighbors": int(n_neighbors),
        "leiden_resolution": float(leiden_resolution),
        "note": "No curated celltype annotation was present in source GEO files; celltype column mirrors generated Leiden clusters.",
    }
    return adata


def write_reports(adata: ad.AnnData, meta: pd.DataFrame, report_txt: Path, report_json: Path) -> None:
    group_counts = adata.obs["group"].astype("string").value_counts(dropna=False).to_dict()
    sample_counts = adata.obs["sample"].astype("string").value_counts(dropna=False).to_dict()
    cluster_counts = adata.obs["cluster"].astype("string").value_counts(dropna=False).to_dict() if "cluster" in adata.obs.columns else {}

    lines = [
        "external_validation_round2 build report",
        f"source_series: GSE190452",
        f"cells: {adata.n_obs}",
        f"genes: {adata.n_vars}",
        f"group_column: group",
        f"disease_column: disease",
        f"sample_column: sample",
        f"donor_column: donor",
        f"celltype_column: celltype",
        f"cluster_column: cluster",
        f"has_umap: {'X_umap' in adata.obsm}",
        "",
        "[group_counts]",
    ]
    lines.extend([f"{k}\t{v}" for k, v in group_counts.items()])
    lines.extend(["", "[sample_counts]"])
    lines.extend([f"{k}\t{v}" for k, v in sample_counts.items()])
    lines.extend(["", "[cluster_counts]"])
    lines.extend([f"{k}\t{v}" for k, v in cluster_counts.items()])
    report_txt.parent.mkdir(parents=True, exist_ok=True)
    report_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload = {
        "source_series": "GSE190452",
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "obs_columns": list(adata.obs.columns),
        "obsm_keys": list(adata.obsm.keys()),
        "uns_keys": sorted(list(adata.uns.keys())),
        "group_counts": group_counts,
        "sample_counts": sample_counts,
        "cluster_counts": cluster_counts,
        "sample_metadata": meta.to_dict(orient="records"),
    }
    report_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    source_dir = Path(args["source_dir"])
    output_h5ad = Path(args["output_h5ad"])
    report_txt = Path(args["report_txt"])
    report_json = Path(args["report_json"])

    adata, meta = build_counts_h5ad(source_dir)
    adata = run_minimal_scanpy_processing(
        adata=adata,
        hvg_n_top=int(args["hvg_n_top"]),
        n_pcs=int(args["n_pcs"]),
        n_neighbors=int(args["n_neighbors"]),
        leiden_resolution=float(args["leiden_resolution"]),
        target_sum=float(args["target_sum"]),
    )

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad, compression="gzip")
    write_reports(adata, meta, report_txt, report_json)

    print(f"Wrote {output_h5ad}")
    print(f"Wrote {report_txt}")
    print(f"Wrote {report_json}")
    print(f"Cells={adata.n_obs}, Genes={adata.n_vars}, Clusters={adata.obs['cluster'].astype(str).nunique()}")


if __name__ == "__main__":
    main()
