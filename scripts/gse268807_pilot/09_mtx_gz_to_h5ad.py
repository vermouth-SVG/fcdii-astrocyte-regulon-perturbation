# -*- coding: utf-8 -*-
"""
将 GEO 提供的单样本 Cell Ranger ARC MTX 三联（*.mtx.gz + barcodes + features，
同目录）按 feature_type 拆出 Gene Expression 或 Peaks 后写出 h5ad，用于「小样本」
本地跑通 02→03 流程，或为云端 SCENIC+ 准备 RNA/peak 输入。

目录示例（解压与否均可，视 scanpy 版本而定；若报错请先解压为 .mtx/.tsv）：
  02_raw/GSE268807/suppl/GSE268807_G150_D/

用法：
  python scripts/09_mtx_gz_to_h5ad.py --mtx-dir "02_raw/GSE268807/suppl" --prefix GSE268807_G150_D --output 03_processed/GSE268807/G150_D_raw.h5ad
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def _make_unique(values: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for val in values:
        if val not in seen:
            seen[val] = 0
            out.append(val)
            continue
        seen[val] += 1
        out.append(f"{val}-{seen[val]}")
    return out


def _infer_clinical_sample_label(prefix: str) -> str:
    label = prefix.removeprefix("GSE268807_")
    special = {
        "G120_F1_N": "G120_N.1",
        "G120_D_TL": "G120_D.2",
        "G120_D_FL": "G120_D.3",
        "G133_N_FL": "G133_N.2",
        "G133_D_FL": "G133_D.2",
    }
    return special.get(label, label)


def _read_barcodes(path: Path) -> list[str]:
    with _open_text(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def _read_features(path: Path) -> tuple[list[str], list[str], list[str], list[bool]]:
    gene_ids: list[str] = []
    gene_names: list[str] = []
    feature_types: list[str] = []
    keep_gene_expression: list[bool] = []
    with _open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                raise SystemExit(f"features 文件列数不足: {path}")
            feature_type = parts[2] if len(parts) >= 3 else "Gene Expression"
            gene_ids.append(parts[0])
            gene_names.append(parts[1])
            feature_types.append(feature_type)
            keep_gene_expression.append(feature_type == "Gene Expression")
    return gene_ids, gene_names, feature_types, keep_gene_expression


def _attach_sample_metadata(adata, metadata_csv: Path, clinical_sample_label: str) -> None:
    if not metadata_csv.is_file():
        print(f"【提示】未找到 metadata：{metadata_csv}，obs 仅写入 sample_prefix。")
        return

    try:
        import pandas as pd
    except ImportError as e:
        raise SystemExit("写入 sample metadata 需要 pandas。") from e

    sm = pd.read_csv(metadata_csv, dtype=str, keep_default_na=False)
    rows = sm.loc[sm["clinical_sample_label"].astype(str).eq(clinical_sample_label)].copy()
    if "assay" in rows.columns:
        rows = rows.loc[rows["assay"].astype(str).str.contains("GEX|snRNA", case=False, regex=True)]
        rows = rows.loc[~rows["assay"].astype(str).str.contains("ATAC", case=False, regex=False)]
    if len(rows) != 1:
        print(
            f"【提示】metadata 中 clinical_sample_label={clinical_sample_label!r} "
            f"匹配到 {len(rows)} 行，obs 未自动注入 donor/group。"
        )
        return

    row = rows.iloc[0]
    for col in rows.columns:
        adata.obs[col] = str(row[col])


def main() -> int:
    try:
        import anndata as ad
        import numpy as np
        import scipy.io
    except ImportError as e:
        raise SystemExit("需要 anndata/numpy/scipy：pip install anndata") from e

    p = argparse.ArgumentParser()
    p.add_argument("--mtx-dir", required=True, type=Path, help="含 mtx/features/barcodes 的目录")
    p.add_argument("--prefix", required=True, help="文件名前缀，如 GSE268807_G150_D")
    p.add_argument("--output", required=True, type=Path, help="输出 .h5ad")
    p.add_argument(
        "--feature-type",
        default="Gene Expression",
        help='从 features 第三列选择的类型；常用 "Gene Expression" 或 Peaks。',
    )
    p.add_argument(
        "--metadata-csv",
        type=Path,
        default=Path("01_metadata/sample_metadata.csv"),
        help="用于写入 donor_id/group/subtype 等 obs 列",
    )
    p.add_argument(
        "--clinical-sample-label",
        default=None,
        help="默认按文件前缀推断；特殊样本如 G120_D_TL -> G120_D.2",
    )
    args = p.parse_args()

    d = args.mtx_dir.resolve()
    pre = args.prefix
    mtx = d / f"{pre}_matrix.mtx.gz"
    if not mtx.is_file():
        mtx = d / f"{pre}_matrix.mtx"
    bc = d / f"{pre}_barcodes.tsv.gz"
    if not bc.is_file():
        bc = d / f"{pre}_barcodes.tsv"
    feat = d / f"{pre}_features.tsv.gz"
    if not feat.is_file():
        feat = d / f"{pre}_features.tsv"
    if not mtx.is_file() or not bc.is_file() or not feat.is_file():
        raise SystemExit(f"缺少三联文件：检查前缀 {pre} 与目录 {d}")

    barcodes = _read_barcodes(bc)
    gene_ids, gene_names, feature_types, keep_gene_expression = _read_features(feat)
    keep_selected = [ft == args.feature_type for ft in feature_types]
    keep = np.asarray(keep_selected, dtype=bool)
    if not keep.any():
        present = sorted(set(feature_types))
        raise SystemExit(f"{feat} 中未找到 feature_type == {args.feature_type!r} 的行；现有类型: {present}")

    with gzip.open(mtx, "rb") if mtx.suffix == ".gz" else mtx.open("rb") as fh:
        mat = scipy.io.mmread(fh).tocsr()
    if mat.shape[0] != len(gene_names) or mat.shape[1] != len(barcodes):
        raise SystemExit(
            f"矩阵维度与 features/barcodes 不一致: matrix={mat.shape}, "
            f"features={len(gene_names)}, barcodes={len(barcodes)}"
        )

    X = mat[keep, :].transpose().tocsr()
    kept_ids = [g for g, flag in zip(gene_ids, keep_selected) if flag]
    kept_names = [g for g, flag in zip(gene_names, keep_selected) if flag]
    kept_types = [g for g, flag in zip(feature_types, keep_selected) if flag]

    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes
    adata.var_names = _make_unique(kept_names)
    adata.var["feature_id"] = kept_ids
    adata.var["feature_name"] = kept_names
    adata.var["feature_type"] = kept_types
    if args.feature_type == "Gene Expression":
        adata.var["gene_id"] = kept_ids
        adata.var["gene_symbol"] = kept_names
    adata.obs["sample_prefix"] = pre

    label = args.clinical_sample_label or _infer_clinical_sample_label(pre)
    adata.obs["clinical_sample_label"] = label
    _attach_sample_metadata(adata, args.metadata_csv, label)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print(
        f"【完成】{args.output} n_obs={adata.n_obs} n_vars={adata.n_vars} "
        f"(feature_type={args.feature_type}；原始 feature 数={len(gene_names)})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
