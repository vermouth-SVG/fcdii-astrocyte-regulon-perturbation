# -*- coding: utf-8 -*-
"""
GSE 单文件 QC + 可选降采样。16GB 本机仍建议分 GSE、设 max_cells。

输入：h5ad（X 或 --layer 指定 counts）
输出：QC 后 h5ad；05_results/qc_metrics_<gse>.csv

不做整合、注释；不下载原始矩阵。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    try:
        import anndata as ad
    except ImportError as e:
        raise SystemExit("请先安装 anndata：pip install anndata") from e

    p = argparse.ArgumentParser(description="Scanpy 单文件 QC + 可选降采样")
    p.add_argument("--input", required=True, type=Path, help="输入 .h5ad 路径")
    p.add_argument("--output", required=True, type=Path, help="输出 .h5ad 路径")
    p.add_argument("--gse-id", default="GSE_unknown", help="用于 qc 报告文件名")
    p.add_argument("--layer", default=None, help="若 counts 不在 X 而在 layer，填 layer 名")
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--max-genes", type=int, default=8000)
    p.add_argument("--max-pct-mito", type=float, default=20.0, help="线粒体基因占比上限（%）")
    p.add_argument("--max-cells", type=int, default=None, help="随机降采样到此细胞数（可选）")
    p.add_argument("--random-seed", type=int, default=0)
    args = p.parse_args()

    root = Path(__file__).resolve().parent.parent
    args.output.parent.mkdir(parents=True, exist_ok=True)
    (root / "05_results").mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    n0 = int(adata.n_obs)

    if args.layer:
        if args.layer not in adata.layers:
            raise SystemExit(f"指定 layer 不存在: {args.layer}")
        X = adata.layers[args.layer]
    else:
        X = adata.X

    if hasattr(X, "getnnz"):
        n_genes = np.asarray(X.getnnz(axis=1)).ravel()
    else:
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel()

    names = adata.var_names.astype(str)
    mito = names.str.upper().str.startswith("MT-")
    if mito.any():
        import scipy.sparse as sp

        mito_mask = mito.to_numpy() if hasattr(mito, "to_numpy") else np.asarray(mito)
        if sp.issparse(X):
            mito_idx = np.where(mito_mask)[0]
            mt = np.asarray(X[:, mito_idx].sum(axis=1)).ravel().astype(float)
            tot = np.asarray(X.sum(axis=1)).ravel().astype(float) + 1e-9
        else:
            Xm = np.asarray(X)
            mt = Xm[:, mito_mask].sum(axis=1)
            tot = Xm.sum(axis=1) + 1e-9
        pct_mito = mt / tot * 100.0
    else:
        pct_mito = np.zeros(adata.n_obs)

    adata.obs["n_genes_by_counts"] = n_genes
    adata.obs["pct_counts_mito"] = pct_mito

    keep = (
        (adata.obs["n_genes_by_counts"] >= args.min_genes)
        & (adata.obs["n_genes_by_counts"] <= args.max_genes)
        & (adata.obs["pct_counts_mito"] < args.max_pct_mito)
    )
    adata = adata[keep].copy()
    n1 = int(adata.n_obs)

    if args.max_cells is not None and n1 > args.max_cells:
        rng = np.random.default_rng(args.random_seed)
        chosen = np.sort(rng.choice(adata.n_obs, size=args.max_cells, replace=False))
        adata = adata[chosen, :].copy()
    n2 = int(adata.n_obs)

    adata.write_h5ad(args.output)

    metrics = pd.DataFrame(
        [
            {
                "gse_id": args.gse_id,
                "n_cells_input": n0,
                "n_cells_after_qc": n1,
                "n_cells_after_subsample": n2,
                "min_genes": args.min_genes,
                "max_genes": args.max_genes,
                "max_pct_mito": args.max_pct_mito,
                "max_cells": args.max_cells,
            }
        ]
    )
    out_csv = root / "05_results" / f"qc_metrics_{args.gse_id}.csv"
    metrics.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"【完成】输出 h5ad: {args.output}")
    print(f"【完成】QC 摘要: {out_csv}")


if __name__ == "__main__":
    main()
