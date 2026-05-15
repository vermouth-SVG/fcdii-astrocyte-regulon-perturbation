#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from pyscenic_analysis_common import (
    choose_group_column,
    default_h5ad_path,
    default_output_root,
    ensure_dir,
    extract_auc_dataframe,
    find_umap_key,
    read_adata,
)


def parse_args() -> argparse.Namespace:
    default_h5ad = default_h5ad_path(__file__)
    default_output_dir = default_output_root(__file__) / "structure_check"

    parser = argparse.ArgumentParser(
        description="Inspect h5ad structure and summarize metadata columns for downstream pySCENIC analysis."
    )
    parser.add_argument("--h5ad", default=str(default_h5ad), help="Input h5ad file.")
    parser.add_argument(
        "--output-dir",
        default=str(default_output_dir),
        help="Directory for CSV/JSON/TXT outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    h5ad_path = Path(args.h5ad).resolve()
    output_dir = ensure_dir(args.output_dir)

    print(f"[1/5] Reading h5ad: {h5ad_path}")
    adata = read_adata(h5ad_path)
    auc_df = extract_auc_dataframe(adata)
    selected_group, obs_summary = choose_group_column(adata.obs)
    candidate_summary = obs_summary[obs_summary["is_group_candidate"]].copy()
    umap_key = find_umap_key(adata.obsm.keys())

    print("[2/5] Writing obs/obsm/uns summaries")
    obs_summary.to_csv(output_dir / "obs_columns_overview.csv", index=False)
    candidate_summary.to_csv(output_dir / "group_candidate_columns.csv", index=False)

    pd.DataFrame({"obsm_key": list(adata.obsm.keys())}).to_csv(output_dir / "obsm_keys.csv", index=False)
    pd.DataFrame({"uns_key": list(adata.uns.keys())}).to_csv(output_dir / "uns_keys.csv", index=False)
    pd.DataFrame({"layer_key": list(adata.layers.keys())}).to_csv(output_dir / "layer_keys.csv", index=False)

    regulon_stats = pd.DataFrame(
        {
            "regulon": auc_df.columns,
            "mean_auc": auc_df.mean(axis=0).values,
            "std_auc": auc_df.std(axis=0).values,
            "nonzero_fraction": (auc_df > 0).mean(axis=0).values,
        }
    ).sort_values(by=["mean_auc", "std_auc"], ascending=[False, False])
    regulon_stats.to_csv(output_dir / "regulon_overview_stats.csv", index=False)
    regulon_stats.head(25).to_csv(output_dir / "top_regulons_overall.csv", index=False)

    group_counts = (
        adata.obs[selected_group]
        .astype("string")
        .fillna("NA")
        .value_counts(dropna=False)
        .rename_axis("group")
        .reset_index(name="n_cells")
    )
    group_counts.to_csv(output_dir / "selected_group_counts.csv", index=False)

    print("[3/5] Writing structure summary")
    summary_payload = {
        "h5ad_path": str(h5ad_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": list(adata.obs.columns),
        "obsm_keys": list(adata.obsm.keys()),
        "uns_keys": list(adata.uns.keys()),
        "layer_keys": list(adata.layers.keys()),
        "selected_group_column": selected_group,
        "selected_group_levels": group_counts.to_dict(orient="records"),
        "umap_key": umap_key,
        "n_regulons": int(auc_df.shape[1]),
    }
    with open(output_dir / "structure_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary_payload, handle, indent=2, ensure_ascii=False)

    with open(output_dir / "selected_group_column.txt", "w", encoding="utf-8") as handle:
        handle.write(selected_group + "\n")

    with open(output_dir / "umap_key.txt", "w", encoding="utf-8") as handle:
        handle.write((umap_key or "NONE") + "\n")

    print("[4/5] Console summary")
    print(f"  cells={adata.n_obs}, genes={adata.n_vars}, regulons={auc_df.shape[1]}")
    print(f"  selected_group_column={selected_group}")
    print(f"  umap_key={umap_key or 'NONE'}")
    if not candidate_summary.empty:
        print("  top_group_candidates=")
        print(candidate_summary.head(6)[["column", "n_unique_non_na", "group_candidate_score"]].to_string(index=False))

    print("[5/5] Done")
    print(f"  outputs={output_dir}")


if __name__ == "__main__":
    main()
