#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
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
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = ROOT / "external_validation_round2" / "input" / "external_validation_round2.h5ad"
DEFAULT_OUT_DIR = ROOT / "external_validation_round2" / "pyscenic_recalc" / "input"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export round2 external h5ad into pySCENIC-ready input files.")
    parser.add_argument("--input-h5ad", default=str(DEFAULT_H5AD))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--batch-size", type=int, default=256, help="Number of cells per CSV write batch.")
    parser.add_argument(
        "--counts-source",
        choices=["auto", "layer_counts", "raw", "x"],
        default="auto",
        help="Prefer counts from layers['counts'], raw.X, or X.",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Only write barcodes/genes/report without emitting expression_for_pyscenic.csv.",
    )
    return parser.parse_args()


def inspect_integerish(matrix, sample_n: int = 100000) -> tuple[bool, float, float]:
    if sparse.issparse(matrix):
        data = matrix.data
        if data.size == 0:
            return True, 0.0, 0.0
        probe = data[: min(sample_n, data.size)]
        return bool(np.allclose(probe, np.round(probe))), float(data.min()), float(data.max())
    arr = np.asarray(matrix)
    flat = arr.ravel()
    nz = flat[flat != 0]
    if nz.size == 0:
        return True, 0.0, 0.0
    probe = nz[: min(sample_n, nz.size)]
    return bool(np.allclose(probe, np.round(probe))), float(nz.min()), float(nz.max())


def resolve_counts_matrix(adata: ad.AnnData, strategy: str):
    if strategy in {"auto", "layer_counts"} and "counts" in adata.layers:
        return adata.layers["counts"], "layers['counts']"
    if strategy in {"auto", "raw"} and adata.raw is not None:
        integerish, _, _ = inspect_integerish(adata.raw.X)
        if integerish:
            return adata.raw.X, "raw.X"
    if strategy in {"auto", "x"}:
        integerish, _, _ = inspect_integerish(adata.X)
        if integerish:
            return adata.X, "X"
    raise ValueError("Could not resolve a reliable counts source. Try --counts-source layer_counts/raw/x explicitly.")


def estimate_csv_size_gb(n_obs: int, n_vars: int, max_value: float) -> float:
    digits = max(1, int(math.log10(max(max_value, 1.0))) + 1)
    approx_chars_per_value = max(2, digits + 1)
    approx_total = n_obs * (n_vars * approx_chars_per_value + 24)
    return approx_total / (1024**3)


def main() -> None:
    args = parse_args()
    input_h5ad = Path(args.input_h5ad)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "expression_for_pyscenic.csv"
    barcodes_path = out_dir / "barcodes.tsv"
    genes_path = out_dir / "genes.tsv"
    report_path = out_dir / "export_report.txt"

    adata = ad.read_h5ad(input_h5ad)
    counts, counts_source = resolve_counts_matrix(adata, args.counts_source)
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(np.asarray(counts))

    integerish, min_nonzero, max_nonzero = inspect_integerish(counts)
    obs_names = adata.obs_names.astype(str).tolist()
    gene_names = adata.var_names.astype(str).tolist()

    barcodes_path.write_text("\n".join(obs_names) + "\n", encoding="utf-8")
    genes_path.write_text("\n".join(gene_names) + "\n", encoding="utf-8")

    wrote_csv = False
    if not args.report_only:
        if csv_path.exists():
            csv_path.unlink()
        header = pd.DataFrame({"CellID": []})
        header = pd.concat([header, pd.DataFrame(columns=gene_names)], axis=1)
        header.to_csv(csv_path, index=False)

        batch_size = max(1, int(args.batch_size))
        for start in range(0, counts.shape[0], batch_size):
            end = min(start + batch_size, counts.shape[0])
            batch = counts[start:end].toarray()
            if integerish:
                batch = np.rint(batch).astype(np.int64, copy=False)
            df = pd.DataFrame(batch, columns=gene_names)
            df.insert(0, "CellID", obs_names[start:end])
            df.to_csv(csv_path, mode="a", index=False, header=False)
        wrote_csv = True

    estimated_gb = estimate_csv_size_gb(counts.shape[0], counts.shape[1], max_nonzero)
    lines = [
        "round2 external_validation pySCENIC export report",
        f"input_h5ad: {input_h5ad}",
        f"counts_source: {counts_source}",
        f"cells: {counts.shape[0]}",
        f"genes: {counts.shape[1]}",
        f"counts_integerish: {integerish}",
        f"counts_min_nonzero: {min_nonzero}",
        f"counts_max_nonzero: {max_nonzero}",
        f"estimated_expression_csv_size_gb: {estimated_gb:.2f}",
        f"expression_csv: {csv_path}",
        f"barcodes_tsv: {barcodes_path}",
        f"genes_tsv: {genes_path}",
        f"expression_csv_written: {wrote_csv}",
    ]
    if args.report_only:
        lines.append("note: report-only mode used; expression_for_pyscenic.csv was intentionally skipped.")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote {barcodes_path}")
    print(f"Wrote {genes_path}")
    print(f"Wrote {report_path}")
    if wrote_csv:
        print(f"Wrote {csv_path}")
    else:
        print("Skipped expression_for_pyscenic.csv due to --report-only")


if __name__ == "__main__":
    main()
