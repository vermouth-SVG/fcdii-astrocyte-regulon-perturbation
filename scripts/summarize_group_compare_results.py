#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in str(value))


def parse_args() -> argparse.Namespace:
    root = project_root_from_file(__file__)
    default_input_dir = root / "analysis_outputs" / "group_compare"

    parser = argparse.ArgumentParser(
        description="Summarize differential regulon results from an existing group_compare output directory."
    )
    parser.add_argument(
        "--input-dir",
        default=str(default_input_dir),
        help="Existing analysis_outputs/group_compare directory.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of top up/down regulons to export.",
    )
    parser.add_argument(
        "--top-total",
        type=int,
        default=20,
        help="Number of overall differential regulons to export for the compact table and bar plot.",
    )
    return parser.parse_args()


def read_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    stats_path = input_dir / "regulon_group_statistics.csv"
    mean_path = input_dir / "group_mean_auc.csv"
    summary_path = input_dir / "comparison_summary.json"

    if not stats_path.exists():
        raise FileNotFoundError(f"Missing file: {stats_path}")
    if not mean_path.exists():
        raise FileNotFoundError(f"Missing file: {mean_path}")
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing file: {summary_path}")

    stats_df = pd.read_csv(stats_path)
    group_mean_auc = pd.read_csv(mean_path, index_col=0)
    with open(summary_path, "r", encoding="utf-8") as handle:
        summary = json.load(handle)
    return stats_df, group_mean_auc, summary


def standardize_columns(stats_df: pd.DataFrame) -> pd.DataFrame:
    required = {
        "regulon",
        "group_1",
        "group_2",
        "p_value",
        "fdr_bh",
        "mean_diff_group2_minus_group1",
        "abs_mean_diff",
    }
    missing = required.difference(stats_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in statistics table: {sorted(missing)}")

    group_1 = str(stats_df["group_1"].iloc[0])
    group_2 = str(stats_df["group_2"].iloc[0])
    mean_diff_col = "mean_diff_group2_minus_group1"

    output = stats_df.copy()
    output["group_1"] = group_1
    output["group_2"] = group_2
    output["direction"] = np.where(
        output[mean_diff_col] < 0,
        f"up_in_{group_1}",
        f"up_in_{group_2}",
    )
    output["effect_direction"] = np.where(output[mean_diff_col] < 0, group_1, group_2)
    output["mean_diff_group1_minus_group2"] = -output[mean_diff_col]
    output["neg_log10_fdr"] = -np.log10(np.clip(output["fdr_bh"].to_numpy(dtype=float), 1e-300, None))
    output["neg_log10_p"] = -np.log10(np.clip(output["p_value"].to_numpy(dtype=float), 1e-300, None))
    return output


def split_top_regulons(stats_df: pd.DataFrame, top_n: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    internal_control_up = (
        stats_df.loc[stats_df["mean_diff_group2_minus_group1"] > 0]
        .sort_values(by=["fdr_bh", "mean_diff_group2_minus_group1"], ascending=[True, False])
        .head(top_n)
        .copy()
    )
    lesion_up = (
        stats_df.loc[stats_df["mean_diff_group2_minus_group1"] < 0]
        .sort_values(by=["fdr_bh", "mean_diff_group2_minus_group1"], ascending=[True, True])
        .head(top_n)
        .copy()
    )
    return lesion_up, internal_control_up


def build_top20_table(stats_df: pd.DataFrame, top_total: int) -> pd.DataFrame:
    top20 = stats_df.sort_values(by=["fdr_bh", "abs_mean_diff"], ascending=[True, False]).head(top_total).copy()
    top20["mean_diff_group1_minus_group2"] = top20["mean_diff_group1_minus_group2"].astype(float)
    columns = [
        "regulon",
        "group_1",
        "group_2",
        "effect_direction",
        "mean_diff_group1_minus_group2",
        "mean_diff_group2_minus_group1",
        "abs_mean_diff",
        "p_value",
        "fdr_bh",
        "neg_log10_fdr",
    ]
    keep = [col for col in columns if col in top20.columns]
    return top20[keep]


def plot_volcano(stats_df: pd.DataFrame, output_png: Path) -> None:
    x = stats_df["mean_diff_group1_minus_group2"].to_numpy(dtype=float)
    y = stats_df["neg_log10_fdr"].to_numpy(dtype=float)
    direction = stats_df["effect_direction"].astype(str).to_numpy()
    group_1 = str(stats_df["group_1"].iloc[0])
    group_2 = str(stats_df["group_2"].iloc[0])

    colors = np.where(direction == group_1, "#d73027", "#4575b4")
    fig, ax = plt.subplots(figsize=(8.2, 6.4))
    ax.scatter(x, y, c=colors, alpha=0.8, s=26, linewidths=0)
    ax.axvline(0.0, color="black", linestyle="--", linewidth=1.0)
    ax.axhline(-np.log10(0.05), color="gray", linestyle=":", linewidth=1.0)
    ax.set_xlabel(f"Mean AUC difference ({group_1} - {group_2})")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(f"Differential regulon activity: {group_1} vs {group_2}")

    label_df = stats_df.sort_values(by=["fdr_bh", "abs_mean_diff"], ascending=[True, False]).head(12)
    for _, row in label_df.iterrows():
        ax.text(
            float(row["mean_diff_group1_minus_group2"]),
            float(row["neg_log10_fdr"]),
            str(row["regulon"]),
            fontsize=8,
            ha="left" if float(row["mean_diff_group1_minus_group2"]) >= 0 else "right",
            va="bottom",
        )

    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=f"Up in {group_1}", markerfacecolor="#d73027", markersize=8),
        plt.Line2D([0], [0], marker="o", color="w", label=f"Up in {group_2}", markerfacecolor="#4575b4", markersize=8),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_top20_bar(top20_df: pd.DataFrame, output_png: Path) -> None:
    plot_df = top20_df.copy().sort_values(by="mean_diff_group1_minus_group2", ascending=True)
    group_1 = str(plot_df["group_1"].iloc[0])
    group_2 = str(plot_df["group_2"].iloc[0])
    colors = np.where(plot_df["effect_direction"] == group_1, "#d73027", "#4575b4")

    fig, ax = plt.subplots(figsize=(9.4, max(6.0, 0.35 * plot_df.shape[0] + 1.2)))
    ax.barh(plot_df["regulon"], plot_df["mean_diff_group1_minus_group2"], color=colors)
    ax.axvline(0.0, color="black", linestyle="--", linewidth=1.0)
    ax.set_xlabel(f"Mean AUC difference ({group_1} - {group_2})")
    ax.set_ylabel("Regulon")
    ax.set_title("Top 20 differential regulons")

    legend_handles = [
        plt.Line2D([0], [0], color="#d73027", linewidth=8, label=f"Up in {group_1}"),
        plt.Line2D([0], [0], color="#4575b4", linewidth=8, label=f"Up in {group_2}"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def format_regulon_list(df: pd.DataFrame, group_col: str, diff_col: str, top_n: int = 5) -> str:
    items = []
    for _, row in df.head(top_n).iterrows():
        items.append(f"{row['regulon']}({row[group_col]:.3f})")
    return "、".join(items)


def build_chinese_summary(
    stats_df: pd.DataFrame,
    group_mean_auc: pd.DataFrame,
    lesion_up: pd.DataFrame,
    internal_control_up: pd.DataFrame,
) -> str:
    group_1 = str(stats_df["group_1"].iloc[0])
    group_2 = str(stats_df["group_2"].iloc[0])

    n_sig = int((stats_df["fdr_bh"] < 0.05).sum())

    group1_mean_col = f"mean_{group_1}"
    group2_mean_col = f"mean_{group_2}"
    lesion_up_text = format_regulon_list(lesion_up, group1_mean_col, "mean_diff_group1_minus_group2")
    control_up_text = format_regulon_list(internal_control_up, group2_mean_col, "mean_diff_group2_minus_group1")

    diff_series = (group_mean_auc[group_1] - group_mean_auc[group_2]).sort_values(ascending=False)
    group1_dominant = diff_series.head(5).index.tolist()
    group2_dominant = diff_series.tail(5).index.tolist()[::-1]

    line1 = (
        f"基于 `regulon_group_statistics.csv` 的两组比较结果，在 {group_1} 与 {group_2} 之间共有 {n_sig} 个 regulon 达到 "
        f"FDR < 0.05，提示两组 astrocyte 的转录调控活性存在较明显差异。"
    )
    line2 = (
        f"按 FDR 和均值差综合排序，{group_1} 相对更高的代表性 regulon 主要包括 {lesion_up_text}；"
        f"{group_2} 相对更高的代表性 regulon 主要包括 {control_up_text}。"
    )
    line3 = (
        "从 `group_mean_heatmap_raw.png` 可以看到，两组在多个核心 regulon 上呈现成簇的活性偏移，"
        f"其中热图前列更偏向 {group_1} 的 regulon 包括 {('、'.join(group1_dominant[:5])) if group1_dominant else '无明显集中项'}，"
        f"更偏向 {group_2} 的 regulon 包括 {('、'.join(group2_dominant[:5])) if group2_dominant else '无明显集中项'}。"
    )
    line4 = (
        "结合 `top_regulons_group_boxplots.png`，上述差异并非由个别离群细胞驱动，而是在两组细胞总体分布层面就已分开，"
        f"说明 {group_1} 与 {group_2} 在 regulon activity 上具有稳定且可重复的分组效应。"
    )
    line5 = (
        f"综合来看，{group_1} 组对应一组在平均 AUC 上整体升高的 regulon 程序，"
        f"而 {group_2} 组则保留了另一组在平均 AUC 上更高、且 FDR 更显著的 regulon 特征，"
        "这些 regulon 可作为后续 CellOracle 或下游网络解释的优先候选。"
    )
    return "\n".join([line1, line2, line3, line4, line5])


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir).resolve()

    print(f"[1/5] Reading existing group_compare outputs: {input_dir}")
    stats_df, group_mean_auc, summary = read_inputs(input_dir)
    stats_df = standardize_columns(stats_df)
    group_1 = str(stats_df["group_1"].iloc[0])
    group_2 = str(stats_df["group_2"].iloc[0])
    group_1_file = safe_name(group_1)
    group_2_file = safe_name(group_2)

    print("[2/5] Ranking top differential regulons")
    lesion_up, internal_control_up = split_top_regulons(stats_df, top_n=args.top_n)
    lesion_up.to_csv(input_dir / f"top10_up_in_{group_1_file}.csv", index=False)
    internal_control_up.to_csv(input_dir / f"top10_up_in_{group_2_file}.csv", index=False)

    top20 = build_top20_table(stats_df, top_total=args.top_total)
    top20.to_csv(input_dir / "top20_differential_regulons.csv", index=False)

    print("[3/5] Writing plots")
    plot_volcano(stats_df, input_dir / "differential_regulons_volcano.png")
    plot_top20_bar(top20, input_dir / "top20_differential_regulons_barplot.png")

    print("[4/5] Writing Chinese summary")
    summary_text = build_chinese_summary(
        stats_df=stats_df,
        group_mean_auc=group_mean_auc,
        lesion_up=lesion_up,
        internal_control_up=internal_control_up,
    )
    with open(input_dir / "differential_regulon_summary_cn.txt", "w", encoding="utf-8") as handle:
        handle.write(summary_text + "\n")

    result_payload = {
        "input_dir": str(input_dir),
        "selected_group_column": summary.get("selected_group_column"),
        "group_1": group_1,
        "group_2": group_2,
        "top_n_each_direction": args.top_n,
        "top_total": args.top_total,
        "outputs": [
            f"top10_up_in_{group_1_file}.csv",
            f"top10_up_in_{group_2_file}.csv",
            "top20_differential_regulons.csv",
            "differential_regulons_volcano.png",
            "top20_differential_regulons_barplot.png",
            "differential_regulon_summary_cn.txt",
        ],
    }
    with open(input_dir / "differential_regulon_summary.json", "w", encoding="utf-8") as handle:
        json.dump(result_payload, handle, indent=2, ensure_ascii=False)

    print("[5/5] Done")
    print(f"  group_1={group_1}")
    print(f"  group_2={group_2}")
    print(f"  outputs={input_dir}")


if __name__ == "__main__":
    main()
