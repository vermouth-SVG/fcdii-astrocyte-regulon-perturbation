# -*- coding: utf-8 -*-
"""
合并多个已 QC 的 h5ad，生成本地子集或云端输入用的 discovery h5ad。

典型用途：
  python scripts/10_concat_h5ad.py \
    --inputs 03_processed/GSE268807/G120_F1_N_qc_test500.h5ad 03_processed/GSE268807/G150_D_qc_test500.h5ad \
    --output 03_processed/GSE268807/GSE268807_smoke_control_lesion.h5ad
"""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    try:
        import anndata as ad
    except ImportError as e:
        raise SystemExit("需要 anndata：pip install anndata") from e

    p = argparse.ArgumentParser(description="合并多个 h5ad（默认 inner join genes）")
    p.add_argument("--inputs", nargs="+", required=True, type=Path, help="输入 h5ad 列表")
    p.add_argument("--output", required=True, type=Path, help="输出 h5ad")
    p.add_argument("--join", choices=["inner", "outer"], default="inner")
    p.add_argument(
        "--merge",
        choices=["same", "unique", "first", "only"],
        default="same",
        help="AnnData var/uns 元数据合并策略；默认保留各输入一致的 var 注释。",
    )
    args = p.parse_args()

    missing = [str(x) for x in args.inputs if not x.is_file()]
    if missing:
        raise SystemExit(f"输入 h5ad 不存在: {missing}")

    adatas = []
    keys = []
    for path in args.inputs:
        x = ad.read_h5ad(path)
        key = path.stem
        x.obs["source_h5ad"] = key
        adatas.append(x)
        keys.append(key)

    out = ad.concat(
        adatas,
        join=args.join,
        label="concat_batch",
        keys=keys,
        index_unique="-",
        merge=args.merge,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(args.output)
    print(f"【完成】{args.output} n_obs={out.n_obs} n_vars={out.n_vars}")


if __name__ == "__main__":
    main()
