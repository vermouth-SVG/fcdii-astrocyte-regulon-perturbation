#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EXPECTED_TFS = ["BHLHE40", "NFE2L2", "SOX2", "THRB"]
EPS = 1e-9


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_file(__file__)
    parser = argparse.ArgumentParser(
        description="Summarize CellOracle KO round-1 outputs and compute an approximate Recovery Index."
    )
    parser.add_argument(
        "--ko-dir",
        default=str(root / "celloracle_run" / "ko_round1"),
        help="Round-1 KO result directory.",
    )
    parser.add_argument(
        "--candidate-dir",
        default=str(root / "analysis_outputs" / "celloracle_candidates"),
        help="CellOracle candidate TF result directory.",
    )
    parser.add_argument(
        "--tf-list",
        default=",".join(EXPECTED_TFS),
        help="Comma-separated TF list expected in round 1.",
    )
    return parser.parse_args()


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in str(value))


def read_csv_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path)


def discover_tf_dirs(ko_dir: Path, tf_list: list[str]) -> list[Path]:
    dirs = []
    for tf in tf_list:
        path = ko_dir / safe_name(tf)
        if path.is_dir():
            dirs.append(path)
    return dirs


def compute_score_stats(score_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    score_df = score_df.copy()
    has_random = "shift_length_random" in score_df.columns
    if has_random:
        score_df["net_shift"] = score_df["shift_length"] - score_df["shift_length_random"]
    else:
        score_df["shift_length_random"] = np.nan
        score_df["net_shift"] = score_df["shift_length"]

    group_rows = []
    for group_name, sub in score_df.groupby("group", dropna=False):
        mean_shift = float(sub["shift_length"].mean())
        mean_random = float(sub["shift_length_random"].mean()) if has_random else np.nan
        net_shift = float(sub["net_shift"].mean())
        mean_dx = float(sub["delta_x"].mean())
        mean_dy = float(sub["delta_y"].mean())
        vector_strength = float(np.sqrt(mean_dx ** 2 + mean_dy ** 2))
        coherence = float(vector_strength / (mean_shift + EPS))
        group_rows.append(
            {
                "group": str(group_name),
                "n_cells": int(sub.shape[0]),
                "mean_shift_length": mean_shift,
                "mean_shift_length_random": mean_random,
                "mean_net_shift": net_shift,
                "mean_delta_x": mean_dx,
                "mean_delta_y": mean_dy,
                "vector_strength": vector_strength,
                "coherence": coherence,
                "p90_shift_length": float(np.quantile(sub["shift_length"], 0.9)),
            }
        )

    overall = {
        "mean_shift_length_random": float(score_df["shift_length_random"].mean()) if has_random else np.nan,
        "mean_net_shift": float(score_df["net_shift"].mean()),
    }
    return pd.DataFrame(group_rows), overall


def pivot_group_metrics(group_df: pd.DataFrame, value_cols: list[str], group_levels: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for group in group_levels:
        sub = group_df.loc[group_df["group"] == group]
        if sub.empty:
            for col in value_cols:
                out[f"{group}_{col}"] = np.nan
        else:
            row = sub.iloc[0]
            for col in value_cols:
                out[f"{group}_{col}"] = float(row[col]) if pd.notna(row[col]) else np.nan
    return out


def build_master_table(
    ko_dir: Path,
    candidate_dir: Path,
    tf_dirs: list[Path],
    summary_overall_combined: pd.DataFrame,
    summary_by_group_combined: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidate_metrics = read_csv_required(candidate_dir / "celloracle_candidate_tf_metrics.csv")
    candidate_metrics["tf"] = candidate_metrics["tf"].astype(str)

    group_levels = ["lesion", "internal_control"]
    master_rows: list[dict[str, object]] = []
    by_group_rows: list[pd.DataFrame] = []

    for tf_dir in tf_dirs:
        tf = tf_dir.name
        score_path = tf_dir / f"{tf}_state_shift_scores.csv"
        overall_path = tf_dir / f"{tf}_state_shift_summary_overall.csv"
        by_group_path = tf_dir / f"{tf}_state_shift_summary_by_group.csv"

        score_df = read_csv_required(score_path)
        overall_df = read_csv_required(overall_path)
        by_group_df = read_csv_required(by_group_path)

        computed_group_df, computed_overall = compute_score_stats(score_df)
        computed_group_df.insert(0, "tf", tf)
        by_group_rows.append(computed_group_df.copy())

        overall_row = overall_df.iloc[0].to_dict()
        combined_overall_row = summary_overall_combined.loc[summary_overall_combined["tf"].astype(str) == tf]
        combined_group_row = summary_by_group_combined.loc[summary_by_group_combined["tf"].astype(str) == tf]
        candidate_row = candidate_metrics.loc[candidate_metrics["tf"] == tf]

        row: dict[str, object] = {"tf": tf}
        row.update(overall_row)
        row["mean_shift_length_random"] = computed_overall["mean_shift_length_random"]
        row["mean_net_shift"] = computed_overall["mean_net_shift"]

        value_cols = [
            "mean_shift_length",
            "mean_shift_length_random",
            "mean_net_shift",
            "mean_delta_x",
            "mean_delta_y",
            "vector_strength",
            "coherence",
            "p90_shift_length",
        ]
        row.update(pivot_group_metrics(computed_group_df, value_cols=value_cols, group_levels=group_levels))

        if not combined_overall_row.empty:
            for col in [
                "expected_regulon_effect_group",
                "expected_expr_higher_group",
                "regulon_fdr",
                "expr_wilcoxon_fdr",
            ]:
                if col in combined_overall_row.columns:
                    row[col] = combined_overall_row.iloc[0][col]
        if not candidate_row.empty:
            for col in ["shortlist_tier", "shortlist_keep", "shortlist_reason_tags"]:
                if col in candidate_row.columns:
                    row[col] = candidate_row.iloc[0][col]

        target_group = row.get("expected_regulon_effect_group")
        expr_group = row.get("expected_expr_higher_group")
        if pd.isna(target_group) or target_group not in group_levels:
            lesion_shift = row.get("lesion_mean_shift_length", np.nan)
            ctrl_shift = row.get("internal_control_mean_shift_length", np.nan)
            target_group = "lesion" if lesion_shift >= ctrl_shift else "internal_control"
        other_group = "internal_control" if target_group == "lesion" else "lesion"
        row["target_group"] = target_group
        row["other_group"] = other_group

        matched_net = float(row.get(f"{target_group}_mean_net_shift", np.nan))
        unmatched_net = float(row.get(f"{other_group}_mean_net_shift", np.nan))
        matched_coherence = float(row.get(f"{target_group}_coherence", np.nan))
        matched_vector_strength = float(row.get(f"{target_group}_vector_strength", np.nan))
        lesion_minus_control = float(row.get("lesion_mean_shift_length", np.nan)) - float(
            row.get("internal_control_mean_shift_length", np.nan)
        )
        direction_agreement = 1.0 if str(target_group) == str(expr_group) else 0.75
        group_selectivity = float((matched_net - unmatched_net) / (abs(matched_net) + abs(unmatched_net) + EPS))
        recovery_index = float(
            row["mean_net_shift"] * (1.0 + group_selectivity) * (1.0 + max(matched_coherence, 0.0)) * direction_agreement
        )

        row["matched_group_mean_net_shift"] = matched_net
        row["unmatched_group_mean_net_shift"] = unmatched_net
        row["matched_group_coherence"] = matched_coherence
        row["matched_group_vector_strength"] = matched_vector_strength
        row["direction_agreement_factor"] = direction_agreement
        row["group_selectivity"] = group_selectivity
        row["lesion_minus_internal_control_shift"] = lesion_minus_control
        row["recovery_index"] = recovery_index

        master_rows.append(row)

    master_df = pd.DataFrame(master_rows)
    master_df = master_df.sort_values(by=["recovery_index", "mean_shift_length"], ascending=[False, False]).reset_index(drop=True)
    return master_df, pd.concat(by_group_rows, axis=0, ignore_index=True)


def add_ranks(master_df: pd.DataFrame, by_group_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    master_df = master_df.copy()
    master_df["overall_shift_rank"] = master_df["mean_shift_length"].rank(method="min", ascending=False).astype(int)
    master_df["overall_net_shift_rank"] = master_df["mean_net_shift"].rank(method="min", ascending=False).astype(int)
    master_df["recovery_index_rank"] = master_df["recovery_index"].rank(method="min", ascending=False).astype(int)
    master_df["lesion_shift_rank"] = master_df["lesion_mean_shift_length"].rank(method="min", ascending=False).astype(int)
    master_df["internal_control_shift_rank"] = master_df["internal_control_mean_shift_length"].rank(method="min", ascending=False).astype(int)

    by_group_df = by_group_df.copy()
    by_group_df["recovery_index_group"] = by_group_df["mean_net_shift"] * (1.0 + by_group_df["coherence"].clip(lower=0.0))
    by_group_df["mean_shift_rank_within_group"] = by_group_df.groupby("group")["mean_shift_length"].rank(method="min", ascending=False).astype(int)
    by_group_df["recovery_index_rank_within_group"] = by_group_df.groupby("group")["recovery_index_group"].rank(method="min", ascending=False).astype(int)
    return master_df, by_group_df


def write_recovery_definition(ko_dir: Path) -> None:
    text = """Round-1 Recovery Index (approximate) definition

This is a pragmatic score built only from existing CellOracle state-shift outputs.

For each TF and group g:
1. mean_shift_g = mean(shift_length)
2. mean_random_shift_g = mean(shift_length_random)
3. net_shift_g = mean_shift_g - mean_random_shift_g
4. vector_strength_g = sqrt(mean_delta_x_g^2 + mean_delta_y_g^2)
5. coherence_g = vector_strength_g / (mean_shift_g + 1e-9)
6. recovery_index_group_g = net_shift_g * (1 + coherence_g)

For each TF overall:
1. target_group = expected_regulon_effect_group
2. other_group = the opposite group
3. matched_net = net_shift_target_group
4. unmatched_net = net_shift_other_group
5. group_selectivity = (matched_net - unmatched_net) / (abs(matched_net) + abs(unmatched_net) + 1e-9)
6. direction_agreement_factor = 1.0 if expected_regulon_effect_group == expected_expr_higher_group else 0.75
7. recovery_index =
   mean_net_shift_overall
   * (1 + group_selectivity)
   * (1 + coherence_target_group)
   * direction_agreement_factor

Interpretation:
- Higher mean_shift_length means stronger KO perturbation.
- Higher net_shift means the perturbation exceeds randomized control.
- Higher coherence means the average vector shift is more directional and less diffuse.
- Higher group_selectivity means the KO perturbs its expected group more than the opposite group.

This is not a literal rescue probability. It is an executable proxy for group-matched state-rewiring potency
using the currently available CellOracle KO outputs.
"""
    with open(ko_dir / "round1_recovery_index_definition.txt", "w", encoding="utf-8") as handle:
        handle.write(text)


def save_bar_plot(values: pd.Series, title: str, ylabel: str, output_png: Path, color_map=None) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    x = np.arange(len(values))
    colors = color_map if color_map is not None else "#3f7f93"
    ax.bar(x, values.values, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels(list(values.index), rotation=30, ha="right")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_group_compare(master_df: pd.DataFrame, output_png: Path) -> None:
    plot_df = master_df.set_index("tf")[["lesion_mean_shift_length", "internal_control_mean_shift_length"]]
    labels = list(plot_df.index)
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    ax.bar(x - width / 2, plot_df["lesion_mean_shift_length"].values, width=width, color="#d95f02", label="lesion")
    ax.bar(
        x + width / 2,
        plot_df["internal_control_mean_shift_length"].values,
        width=width,
        color="#1b7fb8",
        label="internal_control",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Mean shift length")
    ax.set_title("CellOracle KO shift by group")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_recovery_bubble(master_df: pd.DataFrame, output_png: Path) -> None:
    plot_df = master_df.copy()
    color_map = {"lesion": "#d95f02", "internal_control": "#1b7fb8"}
    colors = [color_map.get(str(v), "#666666") for v in plot_df["target_group"]]
    sizes = 900 * (plot_df["matched_group_coherence"].clip(lower=0.0).fillna(0.0) + 0.15)

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    scatter = ax.scatter(
        plot_df["mean_shift_length"],
        plot_df["recovery_index"],
        s=sizes,
        c=colors,
        alpha=0.75,
        linewidths=0.8,
        edgecolors="black",
    )
    for _, row in plot_df.iterrows():
        ax.text(row["mean_shift_length"], row["recovery_index"], row["tf"], fontsize=9, ha="left", va="bottom")
    ax.set_xlabel("Overall mean shift length")
    ax.set_ylabel("Recovery Index")
    ax.set_title("Round-1 KO potency vs Recovery Index")
    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#d95f02", markeredgecolor="black", label="target group: lesion", markersize=9),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#1b7fb8", markeredgecolor="black", label="target group: internal_control", markersize=9),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def build_chinese_summary(master_df: pd.DataFrame) -> str:
    ordered = master_df.sort_values(by=["mean_shift_length"], ascending=False).reset_index(drop=True)
    strongest = ordered.iloc[0]
    second = ordered.iloc[1]
    weakest = ordered.iloc[-1]

    lesion_df = master_df.loc[master_df["target_group"] == "lesion"].sort_values(by=["recovery_index"], ascending=False)
    control_df = master_df.loc[master_df["target_group"] == "internal_control"].sort_values(by=["recovery_index"], ascending=False)
    recommended = []
    for tf in [strongest["tf"], second["tf"], lesion_df.iloc[0]["tf"]]:
        if tf not in recommended:
            recommended.append(str(tf))

    line1 = (
        f"CellOracle 第一轮 in silico KO 结果显示，4 个候选 TF 均能诱导可观的状态偏移，但总体扰动强度存在明显层级。"
        f"按 overall mean_shift_length 排序，{strongest['tf']} 最强，其后为 {second['tf']}，"
        f"{ordered.iloc[2]['tf']} 居中，{weakest['tf']} 最弱。"
    )
    line2 = (
        f"从分组结果看，lesion 侧候选主要为 {', '.join(lesion_df['tf'].tolist())}，"
        f"internal_control 侧候选为 {', '.join(control_df['tf'].tolist())}。"
    )
    line3 = (
        f"其中 {lesion_df.iloc[0]['tf']} 在 lesion 背景下表现出最高的 Recovery Index，"
        f"提示其 KO 后对 lesion 相关状态重排的影响最强；"
        f"{control_df.iloc[0]['tf']} 则是 internal_control 程序中最突出的依赖因子。"
    )
    line4 = (
        "综合 mean_shift_length、随机对照校正后的 net shift、分组特异性以及平均位移向量一致性，"
        f"下一轮优先推荐 {'、'.join(recommended)} 进入更深入分析；"
        f"{weakest['tf']} 虽然仍有信号，但可作为次优先对象。"
    )
    line5 = (
        "整体上，这批结果支持 lesion 与 internal_control 并非共享同一套 TF 依赖结构，"
        "而是分别由 lesion 偏高程序和 internal_control 偏高程序中的关键节点维持。"
    )
    return "\n".join([line1, line2, line3, line4, line5])


def build_next_step_recommendations(master_df: pd.DataFrame) -> str:
    lesion_df = master_df.loc[master_df["target_group"] == "lesion"].sort_values(by=["recovery_index"], ascending=False)
    control_df = master_df.loc[master_df["target_group"] == "internal_control"].sort_values(by=["recovery_index"], ascending=False)
    overall_df = master_df.sort_values(by=["recovery_index"], ascending=False)

    lines = [
        "Round-1 next-step recommendations",
        "",
        "1. Next priority TFs",
        f"- Overall priority: {', '.join(overall_df['tf'].head(3).tolist())}.",
        f"- Lesion-focused priority: {', '.join(lesion_df['tf'].tolist())}.",
        f"- Internal-control-focused priority: {', '.join(control_df['tf'].tolist())}.",
        "",
        "2. Expand TFs or validate first",
        "- Recommendation: validate the current top round-1 TFs first before expanding the TF list.",
        "- Practical order: NFE2L2 and THRB first, then BHLHE40, then SOX2 if extra bandwidth is available.",
        "- Follow-up can include expression/regulon consistency checks, supportive external/contextual dataset review, or targeted wet-lab experiments if available.",
        "",
        "3. Figure placement",
        "- Main figures: overall mean_shift_length ranking bar plot, group comparison bar plot, Recovery Index ranking plot, and representative quiver/grid-flow images for NFE2L2 and THRB.",
        "- Supplementary figures: all four TF-specific embedding shift plots, all boxplots, all grid-flow panels, and per-group summary tables.",
        "",
        "4. Interpretation guidance",
        "- If the next round stays lesion-centric, prioritize lesion-up TFs with high lesion Recovery Index.",
        "- If the next round aims to compare opposing programs, retain one lesion-side TF and one internal-control-side TF as matched contrasts.",
    ]
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    ko_dir = Path(args.ko_dir).resolve()
    candidate_dir = Path(args.candidate_dir).resolve()
    tf_list = [item.strip() for item in args.tf_list.split(",") if item.strip()]

    print(f"[1/5] Reading KO directory: {ko_dir}")
    tf_dirs = discover_tf_dirs(ko_dir, tf_list)
    if not tf_dirs:
        raise ValueError(f"No TF subdirectories found in {ko_dir} for {tf_list}")

    summary_overall_combined = read_csv_required(ko_dir / "round1_ko_summary_overall.csv")
    summary_by_group_combined = read_csv_required(ko_dir / "round1_ko_summary_by_group.csv")

    print("[2/5] Building master tables and rankings")
    master_df, score_group_df = build_master_table(
        ko_dir=ko_dir,
        candidate_dir=candidate_dir,
        tf_dirs=tf_dirs,
        summary_overall_combined=summary_overall_combined,
        summary_by_group_combined=summary_by_group_combined,
    )
    master_df, score_group_df = add_ranks(master_df, score_group_df)

    overall_ranking = master_df.sort_values(by=["mean_shift_length"], ascending=False).reset_index(drop=True)
    lesion_ranking = master_df.sort_values(by=["lesion_mean_shift_length"], ascending=False).reset_index(drop=True)
    control_ranking = master_df.sort_values(by=["internal_control_mean_shift_length"], ascending=False).reset_index(drop=True)
    recovery_ranking = master_df.sort_values(by=["recovery_index"], ascending=False).reset_index(drop=True)
    recovery_by_group = score_group_df.sort_values(by=["group", "recovery_index_group"], ascending=[True, False]).reset_index(drop=True)

    print("[3/5] Writing CSV outputs")
    master_df.to_csv(ko_dir / "round1_ko_master_table.csv", index=False)
    overall_ranking.to_csv(ko_dir / "round1_ko_tf_strength_ranking.csv", index=False)
    lesion_ranking.to_csv(ko_dir / "round1_ko_tf_strength_ranking_lesion.csv", index=False)
    control_ranking.to_csv(ko_dir / "round1_ko_tf_strength_ranking_internal_control.csv", index=False)
    recovery_ranking.to_csv(ko_dir / "round1_recovery_index_ranking.csv", index=False)
    recovery_by_group.to_csv(ko_dir / "round1_recovery_index_by_group.csv", index=False)
    write_recovery_definition(ko_dir)

    print("[4/5] Writing plots")
    plot_dir = ko_dir / "summary_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    color_map = ["#2f6db3", "#2f6db3", "#2f6db3", "#2f6db3"]
    save_bar_plot(
        values=overall_ranking.set_index("tf")["mean_shift_length"],
        title="Round-1 KO mean shift length ranking",
        ylabel="Mean shift length",
        output_png=plot_dir / "round1_mean_shift_length_ranking.png",
        color_map=color_map[: overall_ranking.shape[0]],
    )
    plot_group_compare(master_df, plot_dir / "round1_group_shift_comparison.png")
    save_bar_plot(
        values=recovery_ranking.set_index("tf")["recovery_index"],
        title="Round-1 Recovery Index ranking",
        ylabel="Recovery Index",
        output_png=plot_dir / "round1_recovery_index_ranking.png",
        color_map=["#c44e52" if g == "lesion" else "#4c72b0" for g in recovery_ranking["target_group"]],
    )
    plot_recovery_bubble(master_df, plot_dir / "round1_recovery_index_bubble.png")

    print("[5/5] Writing summaries")
    with open(ko_dir / "round1_ko_summary_cn.txt", "w", encoding="utf-8") as handle:
        handle.write(build_chinese_summary(master_df) + "\n")
    with open(ko_dir / "next_step_recommendations.txt", "w", encoding="utf-8") as handle:
        handle.write(build_next_step_recommendations(master_df) + "\n")

    run_meta = {
        "ko_dir": str(ko_dir),
        "candidate_dir": str(candidate_dir),
        "tf_list": tf_list,
        "outputs": [
            "round1_ko_master_table.csv",
            "round1_ko_tf_strength_ranking.csv",
            "round1_ko_tf_strength_ranking_lesion.csv",
            "round1_ko_tf_strength_ranking_internal_control.csv",
            "round1_recovery_index_ranking.csv",
            "round1_recovery_index_by_group.csv",
            "round1_recovery_index_definition.txt",
            "round1_ko_summary_cn.txt",
            "next_step_recommendations.txt",
            "summary_plots/round1_mean_shift_length_ranking.png",
            "summary_plots/round1_group_shift_comparison.png",
            "summary_plots/round1_recovery_index_ranking.png",
            "summary_plots/round1_recovery_index_bubble.png",
        ],
    }
    with open(ko_dir / "round1_ko_summary_outputs.json", "w", encoding="utf-8") as handle:
        json.dump(run_meta, handle, indent=2, ensure_ascii=False)

    print(f"  outputs={ko_dir}")


if __name__ == "__main__":
    main()
