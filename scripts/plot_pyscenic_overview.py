#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pyscenic_analysis_common import (
    default_h5ad_path,
    default_output_root,
    ensure_dir,
    extract_auc_dataframe,
    find_umap_key,
    read_adata,
)


def parse_args() -> argparse.Namespace:
    default_h5ad = default_h5ad_path(__file__)
    default_output_dir = default_output_root(__file__) / "overview"

    parser = argparse.ArgumentParser(
        description="Create an overview of pySCENIC regulon activity from a merged h5ad."
    )
    parser.add_argument("--h5ad", default=str(default_h5ad), help="Input h5ad file.")
    parser.add_argument("--output-dir", default=str(default_output_dir), help="Output directory.")
    parser.add_argument(
        "--top-n",
        type=int,
        default=12,
        help="Number of top regulons to use for plots.",
    )
    return parser.parse_args()


def make_regulon_stats(auc_df: pd.DataFrame) -> pd.DataFrame:
    stats = pd.DataFrame(
        {
            "regulon": auc_df.columns,
            "mean_auc": auc_df.mean(axis=0).values,
            "std_auc": auc_df.std(axis=0).values,
            "variance_auc": auc_df.var(axis=0).values,
            "nonzero_fraction": (auc_df > 0).mean(axis=0).values,
        }
    )
    return stats.sort_values(by=["mean_auc", "std_auc"], ascending=[False, False]).reset_index(drop=True)


def plot_top_bar(
    values: pd.Series,
    title: str,
    ylabel: str,
    output_png: Path,
    color: str,
) -> None:
    fig, ax = plt.subplots(figsize=(max(8.0, 0.45 * len(values)), 5.5))
    ax.bar(np.arange(len(values)), values.values, color=color)
    ax.set_xticks(np.arange(len(values)))
    ax.set_xticklabels(list(values.index), rotation=45, ha="right")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_umap_panels(
    embedding: np.ndarray,
    auc_df: pd.DataFrame,
    regulons: list[str],
    output_png: Path,
) -> None:
    n_panels = len(regulons)
    ncols = min(4, n_panels)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 4.3 * nrows), squeeze=False)

    for idx, regulon in enumerate(regulons):
        ax = axes[idx // ncols][idx % ncols]
        values = auc_df[regulon].to_numpy(dtype=float)
        order = np.argsort(values)
        scatter = ax.scatter(
            embedding[order, 0],
            embedding[order, 1],
            c=values[order],
            s=6,
            cmap="viridis",
            linewidths=0,
        )
        ax.set_title(regulon)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.03)
        cbar.set_label("AUC")

    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Top regulon activity on UMAP", y=1.01)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    h5ad_path = Path(args.h5ad).resolve()
    output_dir = ensure_dir(args.output_dir)

    print(f"[1/5] Reading h5ad: {h5ad_path}")
    adata = read_adata(h5ad_path)
    auc_df = extract_auc_dataframe(adata)

    print("[2/5] Computing regulon summary statistics")
    regulon_stats = make_regulon_stats(auc_df)
    regulon_stats.to_csv(output_dir / "regulon_overview_stats.csv", index=False)
    regulon_stats.head(max(args.top_n, 25)).to_csv(output_dir / "top_regulons_overall.csv", index=False)

    top_mean = regulon_stats.head(args.top_n).set_index("regulon")["mean_auc"]
    top_std = regulon_stats.sort_values(by=["std_auc", "mean_auc"], ascending=[False, False]).head(args.top_n)
    top_std_series = top_std.set_index("regulon")["std_auc"]
    top_std.to_csv(output_dir / "top_regulons_by_variability.csv", index=False)

    print("[3/5] Writing overview plots")
    plot_top_bar(
        values=top_mean,
        title=f"Top {len(top_mean)} regulons by mean AUC",
        ylabel="Mean AUC",
        output_png=output_dir / "top_regulons_mean_auc.png",
        color="#2f6db3",
    )
    plot_top_bar(
        values=top_std_series,
        title=f"Top {len(top_std_series)} regulons by AUC variability",
        ylabel="AUC standard deviation",
        output_png=output_dir / "top_regulons_std_auc.png",
        color="#d95f02",
    )

    umap_key = find_umap_key(adata.obsm.keys())
    if umap_key is not None:
        print(f"[4/5] UMAP found: {umap_key}")
        embedding = np.asarray(adata.obsm[umap_key], dtype=float)
        umap_regulons = list(top_std["regulon"])
        pd.DataFrame({"regulon": umap_regulons}).to_csv(output_dir / "top_regulons_for_umap.csv", index=False)
        plot_umap_panels(
            embedding=embedding,
            auc_df=auc_df,
            regulons=umap_regulons,
            output_png=output_dir / "top_regulons_umap.png",
        )
    else:
        print("[4/5] No UMAP key found, skipping UMAP plot")
        with open(output_dir / "umap_not_found.txt", "w", encoding="utf-8") as handle:
            handle.write("No UMAP embedding found in adata.obsm.\n")

    print("[5/5] Done")
    print(f"  outputs={output_dir}")


if __name__ == "__main__":
    main()
