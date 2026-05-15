#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import re
import tarfile
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


DEFAULT_SERIES = "GSE140393"
DEFAULT_RAW_TAR = "external_validation/source_gse140393/GSE140393_RAW.tar"
DEFAULT_SOFT = "external_validation/source_gse140393/GSE140393_family.soft.gz"
DEFAULT_EXTRACT_DIR = "external_validation/source_gse140393/extracted"
DEFAULT_OUTPUT = "external_validation/input/external_validation.h5ad"
DEFAULT_REPORT = "external_validation/input/input_build_report.txt"


@dataclass
class SampleBundle:
    sample_id: str
    matrix_name: str
    feature_name: str
    barcode_name: str


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_file(__file__)
    parser = argparse.ArgumentParser(
        description="Build supportive external AnnData/h5ad from GEO raw matrices or an existing h5ad."
    )
    parser.add_argument(
        "--source-h5ad",
        default=None,
        help="Optional existing supportive external h5ad. If provided and found, it is standardized and copied forward.",
    )
    parser.add_argument(
        "--raw-tar",
        default=str(root / DEFAULT_RAW_TAR),
        help="GEO RAW tar archive with matrix.mtx.gz/barcodes.tsv.gz/features.tsv.gz files.",
    )
    parser.add_argument(
        "--soft-file",
        default=str(root / DEFAULT_SOFT),
        help="GEO family SOFT metadata file.",
    )
    parser.add_argument(
        "--extract-dir",
        default=str(root / DEFAULT_EXTRACT_DIR),
        help="Directory used to extract GEO supplementary files.",
    )
    parser.add_argument(
        "--output-h5ad",
        default=str(root / DEFAULT_OUTPUT),
        help="Output h5ad path.",
    )
    parser.add_argument(
        "--report-txt",
        default=str(root / DEFAULT_REPORT),
        help="Output text report path.",
    )
    parser.add_argument(
        "--series-accession",
        default=DEFAULT_SERIES,
        help="GEO series accession for report metadata.",
    )
    parser.add_argument(
        "--keep-all-cells",
        action="store_true",
        help="Disable barcode QC filtering and keep every barcode from the raw matrices.",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=500,
        help="Minimum total UMI counts per barcode for raw GEO matrices.",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        help="Minimum detected genes per barcode for raw GEO matrices.",
    )
    return parser.parse_args()


def make_unique_index(values: pd.Series | list[str], fallback_values: pd.Series | list[str] | None = None) -> pd.Index:
    raw = ["" if pd.isna(v) else str(v).strip() for v in values]
    fallback = None
    if fallback_values is not None:
        fallback = ["" if pd.isna(v) else str(v).strip() for v in fallback_values]
    seen: dict[str, int] = {}
    out: list[str] = []
    for i, value in enumerate(raw):
        base = value
        if not base:
            base = fallback[i] if fallback is not None and fallback[i] else f"feature_{i + 1}"
        count = seen.get(base, 0)
        if count == 0:
            out.append(base)
        else:
            out.append(f"{base}_{count}")
        seen[base] = count + 1
    return pd.Index(out)


def extract_sample_number(text: str | None) -> str | None:
    if not text:
        return None
    match = re.search(r"sample\s+([A-Za-z0-9_-]+)$", text.strip())
    if match:
        return match.group(1)
    numbers = re.findall(r"(\d{3,})", text)
    if numbers:
        return numbers[-1]
    return None


def normalize_condition(value: str | None) -> tuple[str, str]:
    text = (value or "").strip()
    lowered = text.lower()
    if "epilepsy" in lowered or "tle" in lowered or "lesion" in lowered:
        return "lesion", "TLE"
    if "control" in lowered or "normal" in lowered or "autopsy" in lowered:
        return "internal_control", "control"
    return "unknown", text or "unknown"


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
                    "characteristics": [],
                    "descriptions": [],
                    "supp_files": [],
                }
                continue
            if current is None:
                continue
            if line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                current["characteristics"].append(line.split("=", 1)[1].strip())
            elif line.startswith("!Sample_description = "):
                current["descriptions"].append(line.split("=", 1)[1].strip())
            elif line.startswith("!Sample_platform_id = "):
                current["platform_id"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_supplementary_file_"):
                current["supp_files"].append(line.split("=", 1)[1].strip())
    if current is not None:
        records.append(finalize_soft_record(current))
    return pd.DataFrame(records)


def finalize_soft_record(record: dict[str, object]) -> dict[str, object]:
    characteristics = record.get("characteristics", [])
    descriptions = record.get("descriptions", [])
    char_map: dict[str, str] = {}
    for item in characteristics:
        if ":" in item:
            key, value = item.split(":", 1)
            char_map[key.strip().lower()] = value.strip()
    title = str(record.get("title", ""))
    condition_raw = char_map.get("condition", "")
    group, disease = normalize_condition(condition_raw)
    sample = extract_sample_number(title) or extract_sample_number(" ".join(map(str, descriptions))) or str(record["sample_id"])
    chemistry = next((x.replace("10X Genomics, ", "").strip() for x in descriptions if x.startswith("10X Genomics")), "")
    return {
        "sample_id": record["sample_id"],
        "title": title,
        "condition": condition_raw,
        "group": group,
        "disease": disease,
        "tissue": char_map.get("tissue", ""),
        "celltype": char_map.get("cell type", ""),
        "sorting_strategy": char_map.get("sorting strategy", ""),
        "sample": sample,
        "donor": sample,
        "platform_id": record.get("platform_id", ""),
        "chemistry": chemistry,
        "supp_files": ";".join(x for x in record.get("supp_files", []) if x and x != "NONE"),
        "raw_descriptions": " | ".join(map(str, descriptions)),
    }


def extract_raw_tar(raw_tar: Path, extract_dir: Path) -> None:
    extract_dir.mkdir(parents=True, exist_ok=True)
    existing = list(extract_dir.glob("*matrix.mtx.gz"))
    if existing:
        return
    with tarfile.open(raw_tar, "r") as tar:
        tar.extractall(path=extract_dir)


def discover_sample_bundles(extract_dir: Path, metadata_df: pd.DataFrame) -> list[SampleBundle]:
    names = {path.name for path in extract_dir.iterdir() if path.is_file()}
    bundles: list[SampleBundle] = []
    for sample_id in metadata_df["sample_id"].tolist():
        matrix = next((name for name in names if name.startswith(f"{sample_id}_") and name.endswith("matrix.mtx.gz")), None)
        features = next((name for name in names if name.startswith(f"{sample_id}_") and name.endswith("features.tsv.gz")), None)
        barcodes = next((name for name in names if name.startswith(f"{sample_id}_") and name.endswith("barcodes.tsv.gz")), None)
        if matrix and features and barcodes:
            bundles.append(
                SampleBundle(
                    sample_id=sample_id,
                    matrix_name=matrix,
                    feature_name=features,
                    barcode_name=barcodes,
                )
            )
    return bundles


def read_features(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    if df.shape[1] == 2:
        df.columns = ["gene_id", "gene_symbol"]
        df["feature_type"] = "Gene Expression"
    else:
        df = df.iloc[:, :3].copy()
        df.columns = ["gene_id", "gene_symbol", "feature_type"]
    return df


def load_10x_sample(
    extract_dir: Path,
    bundle: SampleBundle,
    sample_meta: dict[str, object],
    reference_gene_ids: list[str] | None,
    keep_all_cells: bool,
    min_counts: int,
    min_genes: int,
) -> tuple[ad.AnnData, list[str]]:
    feature_path = extract_dir / bundle.feature_name
    barcode_path = extract_dir / bundle.barcode_name
    matrix_path = extract_dir / bundle.matrix_name

    features = read_features(feature_path)
    gene_ids = features["gene_id"].astype(str).tolist()
    if reference_gene_ids is not None and gene_ids != reference_gene_ids:
        raise ValueError(f"Feature order mismatch for {bundle.sample_id}.")

    barcodes = pd.read_csv(
        barcode_path,
        sep="\t",
        header=None,
        compression="gzip",
        names=["cell_barcode"],
    )
    matrix = mmread(matrix_path).tocsr().transpose().tocsr()
    if matrix.shape[0] != barcodes.shape[0]:
        raise ValueError(f"Barcode count mismatch for {bundle.sample_id}.")
    if matrix.shape[1] != features.shape[0]:
        raise ValueError(f"Feature count mismatch for {bundle.sample_id}.")

    total_counts = np.asarray(matrix.sum(axis=1)).ravel()
    detected_genes = np.asarray((matrix > 0).sum(axis=1)).ravel()
    if keep_all_cells:
        keep_mask = np.ones(matrix.shape[0], dtype=bool)
    else:
        keep_mask = (total_counts >= min_counts) & (detected_genes >= min_genes)
    matrix = matrix[keep_mask].tocsr()
    barcodes = barcodes.loc[keep_mask].reset_index(drop=True)
    total_counts = total_counts[keep_mask]
    detected_genes = detected_genes[keep_mask]

    obs = pd.DataFrame(index=[f"{bundle.sample_id}__{barcode}" for barcode in barcodes["cell_barcode"].astype(str).tolist()])
    obs["cell_barcode"] = barcodes["cell_barcode"].astype(str).tolist()
    for key, value in sample_meta.items():
        obs[key] = value
    obs["n_counts"] = total_counts
    obs["n_genes_by_counts"] = detected_genes
    obs["source_dataset"] = DEFAULT_SERIES
    obs["geo_series"] = DEFAULT_SERIES
    obs["is_external_validation"] = True

    var = features.copy()
    var.index = pd.Index(var["gene_id"].astype(str).tolist(), name="gene_id")

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    return adata, gene_ids


def standardize_existing_h5ad(input_path: Path) -> ad.AnnData:
    adata = ad.read_h5ad(input_path)
    obs = adata.obs.copy()

    def ensure_column(target: str, candidates: list[str], default: str) -> None:
        for candidate in candidates:
            if candidate in obs.columns:
                obs[target] = obs[candidate].astype("string").fillna(default)
                return
        obs[target] = default

    ensure_column("group", ["group", "disease_group", "condition", "status"], "unknown")
    ensure_column("disease", ["disease", "condition", "diagnosis"], "unknown")
    ensure_column("sample", ["sample", "sample_id", "orig_ident"], "unknown_sample")
    ensure_column("donor", ["donor", "donor_id", "patient", "patient_id", "sample"], "unknown_donor")
    ensure_column("celltype", ["celltype", "cell_type", "annotation", "cluster"], "unknown_celltype")
    if "source_dataset" not in obs.columns:
        obs["source_dataset"] = input_path.stem
    adata.obs = obs
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    return adata


def write_report(
    report_path: Path,
    adata: ad.AnnData,
    source_summary: dict[str, object],
    group_column: str = "group",
    celltype_column: str = "celltype",
) -> None:
    obs = adata.obs.copy()
    group_counts = obs[group_column].astype("string").value_counts(dropna=False) if group_column in obs.columns else pd.Series(dtype="int64")
    celltype_counts = obs[celltype_column].astype("string").value_counts(dropna=False) if celltype_column in obs.columns else pd.Series(dtype="int64")
    sample_counts = obs["sample"].astype("string").value_counts(dropna=False) if "sample" in obs.columns else pd.Series(dtype="int64")

    lines = [
        f"source_dataset: {source_summary.get('source_dataset', '')}",
        f"source_type: {source_summary.get('source_type', '')}",
        f"source_series: {source_summary.get('source_series', '')}",
        f"source_files: {source_summary.get('source_files', '')}",
        f"cells: {adata.n_obs}",
        f"genes: {adata.n_vars}",
        f"group_column: {group_column}",
        f"celltype_column: {celltype_column}",
        f"sample_column: sample",
        f"donor_column: donor",
        f"has_pyscenic_auc: {('X_pyscenic_auc' in adata.obsm and 'pyscenic_regulon_names' in adata.uns)}",
        f"barcode_qc_filter: {source_summary.get('barcode_qc_filter', '')}",
        "",
        "[group_counts]",
    ]
    lines.extend([f"{idx}\t{val}" for idx, val in group_counts.items()])
    lines.extend(["", "[celltype_counts]"])
    lines.extend([f"{idx}\t{val}" for idx, val in celltype_counts.head(20).items()])
    lines.extend(["", "[sample_counts]"])
    lines.extend([f"{idx}\t{val}" for idx, val in sample_counts.items()])
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_json_report(report_path: Path, adata: ad.AnnData, source_summary: dict[str, object]) -> None:
    obs = adata.obs
    payload = {
        "source_dataset": source_summary.get("source_dataset", ""),
        "source_type": source_summary.get("source_type", ""),
        "source_series": source_summary.get("source_series", ""),
        "source_files": source_summary.get("source_files", ""),
        "barcode_qc_filter": source_summary.get("barcode_qc_filter", ""),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "group_column": "group",
        "celltype_column": "celltype",
        "sample_column": "sample",
        "donor_column": "donor",
        "group_counts": obs["group"].astype("string").value_counts(dropna=False).to_dict() if "group" in obs.columns else {},
        "celltype_counts": obs["celltype"].astype("string").value_counts(dropna=False).head(50).to_dict() if "celltype" in obs.columns else {},
        "sample_counts": obs["sample"].astype("string").value_counts(dropna=False).to_dict() if "sample" in obs.columns else {},
        "has_pyscenic_auc": bool("X_pyscenic_auc" in adata.obsm and "pyscenic_regulon_names" in adata.uns),
        "obsm_keys": list(adata.obsm.keys()),
        "uns_keys": sorted(list(adata.uns.keys())),
    }
    report_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def build_from_gse140393(
    raw_tar: Path,
    soft_file: Path,
    extract_dir: Path,
    keep_all_cells: bool,
    min_counts: int,
    min_genes: int,
) -> tuple[ad.AnnData, dict[str, object]]:
    if not raw_tar.exists():
        raise FileNotFoundError(f"RAW tar not found: {raw_tar}")
    if not soft_file.exists():
        raise FileNotFoundError(f"SOFT file not found: {soft_file}")

    metadata = parse_soft_metadata(soft_file)
    extract_raw_tar(raw_tar, extract_dir)
    bundles = discover_sample_bundles(extract_dir, metadata)
    if not bundles:
        raise ValueError(f"No usable 10x bundles discovered in {extract_dir}.")

    metadata = metadata.set_index("sample_id", drop=False)
    adatas: list[ad.AnnData] = []
    reference_gene_ids: list[str] | None = None
    for bundle in bundles:
        sample_meta = metadata.loc[bundle.sample_id].to_dict()
        adata_sample, reference_gene_ids = load_10x_sample(
            extract_dir=extract_dir,
            bundle=bundle,
            sample_meta=sample_meta,
            reference_gene_ids=reference_gene_ids,
            keep_all_cells=keep_all_cells,
            min_counts=min_counts,
            min_genes=min_genes,
        )
        adatas.append(adata_sample)

    adata = ad.concat(adatas, axis=0, join="outer", merge="first", label="input_batch", fill_value=0)
    adata.obs_names_make_unique()
    adata.var["gene_id"] = adata.var.index.astype(str)
    gene_symbols = adata.var["gene_symbol"] if "gene_symbol" in adata.var.columns else pd.Series(adata.var.index.astype(str), index=adata.var.index)
    adata.var_names = make_unique_index(gene_symbols, adata.var["gene_id"])
    adata.var_names.name = "gene_symbol_unique"
    adata.layers["counts"] = adata.X.copy()
    adata.uns["external_validation_source"] = {
        "series": DEFAULT_SERIES,
        "raw_tar": str(raw_tar),
        "soft_file": str(soft_file),
        "samples_loaded": [bundle.sample_id for bundle in bundles],
    }
    source_summary = {
        "source_dataset": DEFAULT_SERIES,
        "source_type": "geo_raw_tar",
        "source_series": DEFAULT_SERIES,
        "source_files": f"{raw_tar.name}; {soft_file.name}",
        "barcode_qc_filter": "none" if keep_all_cells else f"n_counts>={min_counts}, n_genes_by_counts>={min_genes}",
    }
    return adata, source_summary


def main() -> None:
    args = parse_args()
    output_h5ad = Path(args.output_h5ad)
    report_txt = Path(args.report_txt)
    report_json = report_txt.with_suffix(".json")

    if args.source_h5ad and Path(args.source_h5ad).exists():
        adata = standardize_existing_h5ad(Path(args.source_h5ad))
        source_summary = {
            "source_dataset": Path(args.source_h5ad).stem,
            "source_type": "existing_h5ad",
            "source_series": args.series_accession,
            "source_files": Path(args.source_h5ad).name,
        }
    else:
        adata, source_summary = build_from_gse140393(
            raw_tar=Path(args.raw_tar),
            soft_file=Path(args.soft_file),
            extract_dir=Path(args.extract_dir),
            keep_all_cells=args.keep_all_cells,
            min_counts=args.min_counts,
            min_genes=args.min_genes,
        )

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    report_txt.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad, compression="gzip")
    write_report(report_txt, adata, source_summary)
    write_json_report(report_json, adata, source_summary)
    print(f"Wrote {output_h5ad}")
    print(f"Wrote {report_txt}")
    print(f"Cells={adata.n_obs}, Genes={adata.n_vars}")


if __name__ == "__main__":
    main()
