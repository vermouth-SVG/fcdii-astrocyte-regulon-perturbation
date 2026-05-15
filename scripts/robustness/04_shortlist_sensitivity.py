#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from robustness_common import CELLORACLE_CANDIDATES, CORE_TFS, ensure_out_dir, save_heatmap


REGULON_FDR_THRESHOLDS = [0.01, 0.05, 0.10]
LOG2FC_THRESHOLDS = [0.0, 0.25, 0.5]
POSITIVE_FRACTION_THRESHOLDS = [0.0, 0.05, 0.10]
RULES = [
    "regulon_significant_plus_direction",
    "regulon_significant_plus_expression_significant",
    "regulon_significant_plus_expression_significant_plus_direction",
]


def prepare_metrics() -> pd.DataFrame:
    if not CELLORACLE_CANDIDATES.exists():
        raise FileNotFoundError(f"缺少候选 TF 指标表: {CELLORACLE_CANDIDATES}")
    df = pd.read_csv(CELLORACLE_CANDIDATES)
    df["tf"] = df["tf"].astype(str).str.upper()
    for col in [
        "regulon_fdr",
        "expr_wilcoxon_fdr",
        "log2fc_lesion_vs_internal_control",
        "positive_fraction_in_higher_group",
    ]:
        if col not in df.columns:
            df[col] = np.nan
    if "direction_consistent_with_regulon" not in df.columns:
        df["direction_consistent_with_regulon"] = False
    if "expr_higher_group" not in df.columns:
        df["expr_higher_group"] = "unknown"
    if "regulon_effect_group" not in df.columns:
        df["regulon_effect_group"] = "unknown"
    df["abs_log2fc"] = df["log2fc_lesion_vs_internal_control"].abs()
    df["support_sort_score"] = (
        -np.log10(df["regulon_fdr"].astype(float).clip(lower=1e-300)).fillna(0) * 0.35
        + df["abs_log2fc"].fillna(0) * 0.35
        + df["positive_fraction_in_higher_group"].fillna(0) * 0.20
        + df["direction_consistent_with_regulon"].astype(float).fillna(0) * 0.10
    )
    return df


def retained_mask(df: pd.DataFrame, fdr_thr: float, logfc_thr: float, pos_thr: float, rule: str) -> pd.Series:
    base = (
        (df["regulon_fdr"] <= fdr_thr)
        & (df["abs_log2fc"] >= logfc_thr)
        & (df["positive_fraction_in_higher_group"] >= pos_thr)
    )
    if rule == "regulon_significant_plus_direction":
        return base & df["direction_consistent_with_regulon"].fillna(False).astype(bool)
    if rule == "regulon_significant_plus_expression_significant":
        return base & (df["expr_wilcoxon_fdr"] <= 0.05)
    if rule == "regulon_significant_plus_expression_significant_plus_direction":
        return base & (df["expr_wilcoxon_fdr"] <= 0.05) & df["direction_consistent_with_regulon"].fillna(False).astype(bool)
    raise ValueError(rule)


def main() -> None:
    out_dir = ensure_out_dir()
    print("执行候选 TF shortlist 阈值敏感性检查...")
    metrics = prepare_metrics()

    combo_rows = []
    long_rows = []
    combo_id = 0
    for fdr_thr in REGULON_FDR_THRESHOLDS:
        for logfc_thr in LOG2FC_THRESHOLDS:
            for pos_thr in POSITIVE_FRACTION_THRESHOLDS:
                for rule in RULES:
                    combo_id += 1
                    mask = retained_mask(metrics, fdr_thr, logfc_thr, pos_thr, rule)
                    retained = metrics[mask].sort_values("support_sort_score", ascending=False).copy()
                    retained["rank_in_combo"] = np.arange(1, retained.shape[0] + 1)
                    members = retained["tf"].tolist()
                    combo_rows.append(
                        {
                            "combo_id": combo_id,
                            "regulon_fdr_threshold": fdr_thr,
                            "expr_abs_log2fc_threshold": logfc_thr,
                            "positive_fraction_threshold": pos_thr,
                            "rule": rule,
                            "n_retained": len(members),
                            "shortlist_members": ";".join(members),
                            "core4_retained": ";".join([tf for tf in CORE_TFS if tf in members]),
                        }
                    )
                    rank_map = dict(zip(retained["tf"], retained["rank_in_combo"]))
                    for _, row in metrics.iterrows():
                        tf = row["tf"]
                        long_rows.append(
                            {
                                "combo_id": combo_id,
                                "tf": tf,
                                "retained": tf in members,
                                "rank_in_combo": rank_map.get(tf, np.nan),
                                "regulon_fdr_threshold": fdr_thr,
                                "expr_abs_log2fc_threshold": logfc_thr,
                                "positive_fraction_threshold": pos_thr,
                                "rule": rule,
                                "regulon_effect_group": row["regulon_effect_group"],
                                "expr_higher_group": row["expr_higher_group"],
                                "regulon_fdr": row["regulon_fdr"],
                                "expr_wilcoxon_fdr": row["expr_wilcoxon_fdr"],
                                "log2fc_lesion_vs_internal_control": row["log2fc_lesion_vs_internal_control"],
                                "positive_fraction_in_higher_group": row["positive_fraction_in_higher_group"],
                                "support_sort_score": row["support_sort_score"],
                            }
                        )

    grid = pd.DataFrame(combo_rows)
    long_df = pd.DataFrame(long_rows)
    grid.to_csv(out_dir / "04_shortlist_sensitivity_grid.csv", index=False)
    long_df.to_csv(out_dir / "04_shortlist_sensitivity_long.csv", index=False)

    freq = (
        long_df.groupby("tf", as_index=False)
        .agg(
            retained_count=("retained", "sum"),
            total_combinations=("retained", "size"),
            retention_frequency=("retained", "mean"),
            median_rank_when_retained=("rank_in_combo", "median"),
            best_rank=("rank_in_combo", "min"),
            worst_rank=("rank_in_combo", "max"),
            regulon_effect_group=("regulon_effect_group", "first"),
            expr_higher_group=("expr_higher_group", "first"),
            regulon_fdr=("regulon_fdr", "first"),
            expr_wilcoxon_fdr=("expr_wilcoxon_fdr", "first"),
            log2fc_lesion_vs_internal_control=("log2fc_lesion_vs_internal_control", "first"),
            positive_fraction_in_higher_group=("positive_fraction_in_higher_group", "first"),
            support_sort_score=("support_sort_score", "first"),
        )
        .sort_values(["retention_frequency", "support_sort_score"], ascending=[False, False])
    )
    freq.to_csv(out_dir / "04_shortlist_membership_frequency.csv", index=False)

    priority = freq[freq["tf"].isin(CORE_TFS)].copy()
    priority["rank_stability_component"] = 1 / priority["median_rank_when_retained"].fillna(99)
    priority["priority_stability_score"] = 0.75 * priority["retention_frequency"] + 0.25 * priority["rank_stability_component"]
    priority = priority.sort_values("priority_stability_score", ascending=False)
    priority.to_csv(out_dir / "04_shortlist_priority_stability.csv", index=False)

    heat = (
        long_df.assign(threshold_pair=lambda x: "FDR<=" + x["regulon_fdr_threshold"].astype(str) + "|logFC>=" + x["expr_abs_log2fc_threshold"].astype(str))
        .groupby(["tf", "threshold_pair"])["retained"]
        .mean()
        .unstack(fill_value=0)
    )
    heat = heat.reindex(freq["tf"].tolist())
    save_heatmap(heat, out_dir / "04_shortlist_sensitivity_heatmap.png", "Shortlist retention frequency across threshold pairs", cmap="YlGnBu", center_zero=False)

    def describe(tf: str) -> str:
        row = freq[freq["tf"] == tf]
        if row.empty:
            return "未进入候选指标表"
        r = row.iloc[0]
        return (
            f"保留频率={r['retention_frequency']:.2f}, "
            f"中位排名={r['median_rank_when_retained'] if pd.notna(r['median_rank_when_retained']) else 'NA'}, "
            f"regulon方向={r['regulon_effect_group']}, 表达方向={r['expr_higher_group']}"
        )

    lines = [
        "候选 TF shortlist 阈值敏感性检查总结",
        "=" * 45,
        f"共评估阈值组合: {grid.shape[0]} 套",
        "阈值: regulon FDR = 0.01/0.05/0.10；表达 |log2FC| = 0/0.25/0.5；阳性细胞比例 = 0/0.05/0.10；规则 = 3 类。",
        "",
        "重点 4 TF:",
        f"- NFE2L2: {describe('NFE2L2')}",
        f"- THRB: {describe('THRB')}",
        f"- BHLHE40: {describe('BHLHE40')}",
        f"- SOX2: {describe('SOX2')}",
        "",
        "解释原则: 保留频率越高，说明该 TF 对阈值选择越不敏感；中位排名越靠前，说明优先级更稳定。",
        "",
        "全部 TF 保留频率:",
        freq.to_string(index=False),
    ]
    (out_dir / "04_shortlist_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("完成 shortlist 阈值敏感性检查。")
    print(f"输出目录: {out_dir}")


if __name__ == "__main__":
    main()
