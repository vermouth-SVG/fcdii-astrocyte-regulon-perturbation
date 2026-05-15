#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = ROOT / "external_validation_round2" / "input" / "external_validation_round2.h5ad"
DEFAULT_GMT = ROOT / "external_validation_round2" / "pyscenic_recalc" / "output" / "discovery_regulons_4tf.gmt"
DEFAULT_OUT_CSV = ROOT / "external_validation_round2" / "pyscenic_recalc" / "input" / "expression_for_pyscenic_4tf_targets.csv"
DEFAULT_OUT_GENES = ROOT / "external_validation_round2" / "pyscenic_recalc" / "input" / "genes_4tf_targets.txt"
DEFAULT_REPORT_TXT = ROOT / "external_validation_round2" / "pyscenic_recalc" / "input" / "build_4tf_matrix_report.txt"
DEFAULT_REPORT_JSON = ROOT / "external_validation_round2" / "pyscenic_recalc" / "input" / "build_4tf_matrix_report.json"
DEFAULT_CHUNK_SIZE = 2000


def parse_gmt(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    ordered_union: list[str] = []
    seen = set()
    by_regulon: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            regulon = parts[0]
            genes = [gene.strip() for gene in parts[2:] if gene.strip()]
            by_regulon[regulon] = genes
            for gene in genes:
                if gene not in seen:
                    seen.add(gene)
                    ordered_union.append(gene)
    return ordered_union, by_regulon


def choose_matrix_source(adata: ad.AnnData) -> tuple[str, bool]:
    if "counts" in adata.layers.keys():
        return "counts", True
    return "X", False


def match_genes_to_var_names(var_names: pd.Index, target_genes: list[str]) -> tuple[list[int], list[str], list[str]]:
    exact_map = {str(name): idx for idx, name in enumerate(var_names.tolist())}
    upper_map: dict[str, list[int]] = {}
    for idx, name in enumerate(var_names.tolist()):
        upper_map.setdefault(str(name).upper(), []).append(idx)

    matched_pairs: list[tuple[int, str]] = []
    missing: list[str] = []
    seen_indices = set()
    for gene in target_genes:
        idx = exact_map.get(gene)
        if idx is None:
            candidates = upper_map.get(gene.upper(), [])
            if len(candidates) == 1:
                idx = candidates[0]
        if idx is None:
            missing.append(gene)
            continue
        if idx in seen_indices:
            continue
        seen_indices.add(idx)
        matched_pairs.append((idx, str(var_names[idx])))

    matched_pairs.sort(key=lambda x: x[0])
    matched_indices = [idx for idx, _ in matched_pairs]
    matched_names = [name for _, name in matched_pairs]
    return matched_indices, matched_names, missing


def extract_chunk(view: ad.AnnData, matrix_key: str):
    matrix = view.layers["counts"] if matrix_key == "counts" else view.X
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build minimal round2 expression matrix for 4 discovery regulons AUCell projection.")
    parser.add_argument("--h5ad", default=str(DEFAULT_H5AD))
    parser.add_argument("--gmt", default=str(DEFAULT_GMT))
    parser.add_argument("--output-csv", default=str(DEFAULT_OUT_CSV))
    parser.add_argument("--output-genes", default=str(DEFAULT_OUT_GENES))
    parser.add_argument("--report-txt", default=str(DEFAULT_REPORT_TXT))
    parser.add_argument("--report-json", default=str(DEFAULT_REPORT_JSON))
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad)
    gmt_path = Path(args.gmt)
    output_csv = Path(args.output_csv)
    output_genes = Path(args.output_genes)
    report_txt = Path(args.report_txt)
    report_json = Path(args.report_json)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_genes.parent.mkdir(parents=True, exist_ok=True)
    report_txt.parent.mkdir(parents=True, exist_ok=True)
    report_json.parent.mkdir(parents=True, exist_ok=True)

    target_genes, regulon_gene_map = parse_gmt(gmt_path)
    if not target_genes:
        raise ValueError(f"No target genes found in GMT: {gmt_path}")

    adata = ad.read_h5ad(h5ad_path, backed="r")
    matrix_key, has_counts = choose_matrix_source(adata)
    matched_indices, matched_names, missing_genes = match_genes_to_var_names(adata.var_names, target_genes)
    if not matched_indices:
        raise ValueError("No GMT target genes matched round2 h5ad var_names.")

    obs_names = adata.obs_names.astype(str).tolist()
    n_cells = adata.n_obs
    n_genes = len(matched_names)

    with output_genes.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(matched_names) + "\n")

    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["CellID", *matched_names])
        for start in range(0, n_cells, args.chunk_size):
            end = min(start + args.chunk_size, n_cells)
            chunk_view = adata[start:end, matched_indices]
            chunk_array = extract_chunk(chunk_view, matrix_key)
            if chunk_array.ndim == 1:
                chunk_array = chunk_array.reshape(-1, n_genes)
            for row_idx, cell_id in enumerate(obs_names[start:end]):
                writer.writerow([cell_id, *chunk_array[row_idx].tolist()])
            del chunk_view
            del chunk_array

    report = {
        "input_h5ad": str(h5ad_path),
        "input_gmt": str(gmt_path),
        "counts_layer_present": has_counts,
        "matrix_source_used": "layers['counts']" if has_counts else "X",
        "n_cells": int(n_cells),
        "n_genes_h5ad": int(adata.n_vars),
        "n_regulons_in_gmt": int(len(regulon_gene_map)),
        "n_target_genes_union_from_gmt": int(len(target_genes)),
        "n_target_genes_matched_in_h5ad": int(len(matched_names)),
        "n_target_genes_missing_in_h5ad": int(len(missing_genes)),
        "missing_genes_preview": missing_genes[:50],
        "output_expression_csv": str(output_csv),
        "output_expression_csv_size_bytes": int(output_csv.stat().st_size) if output_csv.exists() else 0,
        "output_genes_txt": str(output_genes),
        "chunk_size": int(args.chunk_size),
        "regulon_gene_counts": {key: len(value) for key, value in regulon_gene_map.items()},
    }
    matrix_source_text = "layers['counts']" if has_counts else "X"

    report_lines = [
        f"Input h5ad: {h5ad_path}",
        f"Input GMT: {gmt_path}",
        f"Counts layer present: {has_counts}",
        f"Matrix source used: {matrix_source_text}",
        f"Cells exported: {n_cells}",
        f"Genes in h5ad: {adata.n_vars}",
        f"Discovery regulons in GMT: {len(regulon_gene_map)}",
        f"Union target genes in GMT: {len(target_genes)}",
        f"Matched target genes in h5ad: {len(matched_names)}",
        f"Missing target genes in h5ad: {len(missing_genes)}",
        f"Output expression CSV: {output_csv}",
        f"Output CSV size bytes: {report['output_expression_csv_size_bytes']}",
        f"Output genes list: {output_genes}",
        f"Chunk size: {args.chunk_size}",
    ]
    if missing_genes:
        report_lines.append("Missing genes preview: " + ", ".join(missing_genes[:50]))

    report_txt.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"Wrote {output_csv}")
    print(f"Wrote {output_genes}")
    print(f"Wrote {report_txt}")
    print(f"Wrote {report_json}")
    print(f"Cells={n_cells}, genes={n_genes}, counts_used={has_counts}")


if __name__ == "__main__":
    main()
