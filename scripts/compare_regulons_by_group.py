#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from pyscenic_analysis_common import (
    benjamini_hochberg,
    choose_group_column,
    default_h5ad_path,
    default_output_root,
    ensure_dir,
    extract_auc_dataframe,
    plot_heatmap,
    read_adata,
    row_zscore,
)


def parse_args() -> argparse.Namespace:
    default_h5ad = default_h5ad_path(__file__)
    default_output_dir = default_output_root(__file__) / "group_compare"

    parser = argparse.ArgumentParser(
        description="Compare pySCENIC regulon activity across automatically selected metadata groups."
    )
    parser.add_argument("--h5ad", default=str(default_h5ad), help="Input h5ad file.")
    parser.add_argument("--output-dir", default=str(default_output_dir), help="Output directory.")
    parser.add_argument(
        "--group-column",
        default=None,
        help="Optional obs column to force as grouping variable. Default: auto-select.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=12,
        help="Number of regulons to highlight in boxplots.",
    )
    parser.add_argument(
        "--min-group-size",
        type=int,
        default=20,
        help="Drop group levels smaller than this size before statistics.",
    )
    return parser.parse_args()


def sanitize_group_series(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("NA")


def two_group_stats(
    auc_df: pd.DataFrame,
    groups: pd.Series,
) -> pd.DataFrame:
    levels = list(groups.unique())
    if len(levels) != 2:
        raise ValueError("two_group_stats requires exactly two group levels.")

    g1, g2 = levels
    mask1 = groups == g1
    mask2 = groups == g2
    rows: list[dict[str, object]] = []

    for regulon in auc_df.columns:
        x = auc_df.loc[mask1, regulon].to_numpy(dtype=float)
        y = auc_df.loc[mask2, regulon].to_numpy(dtype=float)
        mean1 = float(np.mean(x))
        mean2 = float(np.mean(y))
        median1 = float(np.median(x))
        median2 = float(np.median(y))

        if np.allclose(x, x[0]) and np.allclose(y, y[0]) and np.isclose(x[0], y[0]):
            stat = 0.0
            p_value = 1.0
        else:
            stat, p_value = stats.mannwhitneyu(x, y, alternative="two-sided", method="auto")

        rows.append(
            {
                "regulon": regulon,
                "test": "mann_whitney_u",
                "group_1": g1,
                "group_2": g2,
                f"n_{g1}": int(x.shape[0]),
                f"n_{g2}": int(y.shape[0]),
                f"mean_{g1}": mean1,
                f"mean_{g2}": mean2,
                f"median_{g1}": median1,
                f"median_{g2}": median2,
                "mean_diff_group2_minus_group1": mean2 - mean1,
                "abs_mean_diff": abs(mean2 - mean1),
                "statistic": float(stat),
                "p_value": float(p_value),
            }
        )

    result = pd.DataFrame(rows)
    result["fdr_bh"] = benjamini_hochberg(result["p_value"].to_numpy(dtype=float))
    result = result.sort_values(by=["fdr_bh", "p_value", "abs_mean_diff"], ascending=[True, True, False]).reset_index(drop=True)
    return result


def multi_group_stats(
    auc_df: pd.DataFrame,
    groups: pd.Series,
) -> pd.DataFrame:
    levels = list(groups.unique())
    if len(levels) <= 2:
        raise ValueError("multi_group_stats requires more than two group levels.")

    level_masks = {level: (groups == level).to_numpy() for level in levels}
    rows: list[dict[str, object]] = []

    for regulon in auc_df.columns:
        arrays = [auc_df.loc[level_masks[level], regulon].to_numpy(dtype=float) for level in levels]
        level_means = {f"mean_{level}": float(np.mean(arr)) for level, arr in zip(levels, arrays)}

        if all(np.allclose(arr, arrays[0][0]) for arr in arrays) and all(np.isclose(np.mean(arr), np.mean(arrays[0])) for arr in arrays):
            stat = 0.0
            p_value = 1.0
        else:
            stat, p_value = stats.kruskal(*arrays, nan_policy="omit")

        mean_values = np.array(list(level_means.values()), dtype=float)
        rows.append(
            {
                "regulon": regulon,
                "test": "kruskal_wallis",
                "n_groups": len(levels),
                "statistic": float(stat),
                "p_value": float(p_value),
                "max_mean_minus_min_mean": float(mean_values.max() - mean_values.min()),
                **level_means,
            }
        )

    result = pd.DataFrame(rows)
    result["fdr_bh"] = benjamini_hochberg(result["p_value"].to_numpy(dtype=float))
    result = result.sort_values(
        by=["fdr_bh", "p_value", "max_mean_minus_min_mean"],
        ascending=[True, True, False],
    ).reset_index(drop=True)
    return result


def choose_top_regulons_for_boxplot(
    stats_df: pd.DataFrame,
    group_mean_auc: pd.DataFrame,
    top_n: int,
) -> list[str]:
    if "fdr_bh" in stats_df.columns:
        ranked = stats_df.copy()
        if "abs_mean_diff" in ranked.columns:
            ranked = ranked.sort_values(by=["fdr_bh", "abs_mean_diff"], ascending=[True, False])
        elif "max_mean_minus_min_mean" in ranked.columns:
            ranked = ranked.sort_values(by=["fdr_bh", "max_mean_minus_min_mean"], ascending=[True, False])
        return list(ranked["regulon"].head(top_n))

    variability = group_mean_auc.var(axis=1).sort_values(ascending=False)
    return list(variability.head(top_n).index)


def plot_group_boxplots(
    auc_df: pd.DataFrame,
    groups: pd.Series,
    regulons: list[str],
    output_png: Path,
) -> None:
    levels = list(groups.unique())
    n_panels = len(regulons)
    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.3 * ncols, 4.2 * nrows), squeeze=False)
    colors = plt.cm.Set2(np.linspace(0.0, 1.0, len(levels)))

    for idx, regulon in enumerate(regulons):
        ax = axes[idx // ncols][idx % ncols]
        data = [auc_df.loc[groups == level, regulon].to_numpy(dtype=float) for level in levels]
        box = ax.boxplot(data, patch_artist=True, tick_labels=levels, showfliers=False)
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.8)
        ax.set_title(regulon)
        ax.set_ylabel("AUC")
        ax.tick_params(axis="x", rotation=30)

    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Top regulons by selected grouping column", y=1.01)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    h5ad_path = Path(args.h5ad).resolve()
    output_dir = ensure_dir(args.output_dir)

    print(f"[1/6] Reading h5ad: {h5ad_path}")
    adata = read_adata(h5ad_path)
    auc_df = extract_auc_dataframe(adata)

    print("[2/6] Selecting grouping column")
    selected_group, obs_summary = choose_group_column(adata.obs, preferred_column=args.group_column)
    obs_summary.to_csv(output_dir / "obs_columns_overview.csv", index=False)
    obs_summary[obs_summary["is_group_candidate"]].to_csv(output_dir / "group_candidate_columns.csv", index=False)

    groups = sanitize_group_series(adata.obs[selected_group])
    group_counts = groups.value_counts(dropna=False)
    kept_levels = group_counts[group_counts >= args.min_group_size].index.tolist()
    if len(kept_levels) < 2:
        raise ValueError(
            f"Selected group column '{selected_group}' has fewer than 2 usable levels after min_group_size={args.min_group_size}."
        )
    groups = groups[groups.isin(kept_levels)].copy()
    auc_df = auc_df.loc[groups.index].copy()
    group_counts = groups.value_counts(dropna=False)

    group_counts.rename_axis("group").reset_index(name="n_cells").to_csv(
        output_dir / "selected_group_counts.csv", index=False
    )
    with open(output_dir / "selected_group_column.txt", "w", encoding="utf-8") as handle:
        handle.write(selected_group + "\n")

    print(f"  selected_group_column={selected_group}")
    print(f"  kept_levels={list(group_counts.index)}")

    print("[3/6] Computing group mean AUC")
    group_mean_auc = auc_df.groupby(groups, observed=False).mean().T
    group_mean_auc = group_mean_auc.loc[group_mean_auc.mean(axis=1).sort_values(ascending=False).index]
    group_mean_auc.to_csv(output_dir / "group_mean_auc.csv")

    group_mean_auc_z = row_zscore(group_mean_auc)
    group_mean_auc_z.to_csv(output_dir / "group_mean_auc_row_zscore.csv")

    plot_heatmap(
        data=group_mean_auc,
        output_png=output_dir / "group_mean_heatmap_raw.png",
        title=f"Group mean regulon AUC: {selected_group}",
        cmap="viridis",
        colorbar_label="Mean AUC",
    )
    plot_heatmap(
        data=group_mean_auc_z,
        output_png=output_dir / "group_mean_heatmap_row_zscore.png",
        title=f"Group mean regulon AUC (row z-score): {selected_group}",
        cmap="RdBu_r",
        colorbar_label="Row z-score",
    )

    print("[4/6] Running group statistics")
    if group_counts.shape[0] == 2:
        stats_df = two_group_stats(auc_df=auc_df, groups=groups)
        stats_mode = "two_group_mann_whitney_u"
    else:
        stats_df = multi_group_stats(auc_df=auc_df, groups=groups)
        stats_mode = "multi_group_kruskal_wallis"
    stats_df.to_csv(output_dir / "regulon_group_statistics.csv", index=False)

    print("[5/6] Creating boxplots")
    top_regulons = choose_top_regulons_for_boxplot(stats_df=stats_df, group_mean_auc=group_mean_auc, top_n=args.top_n)
    pd.DataFrame({"regulon": top_regulons}).to_csv(output_dir / "top_regulons_for_boxplot.csv", index=False)
    plot_group_boxplots(
        auc_df=auc_df,
        groups=groups,
        regulons=top_regulons,
        output_png=output_dir / "top_regulons_group_boxplots.png",
    )

    print("[6/6] Writing summary")
    summary = {
        "h5ad_path": str(h5ad_path),
        "selected_group_column": selected_group,
        "group_levels": [{"group": str(k), "n_cells": int(v)} for k, v in group_counts.items()],
        "stats_mode": stats_mode,
        "n_regulons": int(auc_df.shape[1]),
        "top_regulons_for_boxplot": top_regulons,
    }
    with open(output_dir / "comparison_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    print(f"  stats_mode={stats_mode}")
    print(f"  outputs={output_dir}")


if __name__ == "__main__":
    main()
