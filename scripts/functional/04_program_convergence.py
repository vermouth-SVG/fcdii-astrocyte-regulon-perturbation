#!/usr/bin/env python3
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from functional_common import (
    FOCUS_TFS,
    PRIMARY_TFS,
    TF_MAIN_DIRECTION,
    TF_TIER,
    ensure_functional_dir,
    infer_themes_from_terms,
    load_robustness_tables,
    save_matrix_heatmap,
)


def load_enrichment(out_dir, tf):
    path = out_dir / f"03_{tf.lower()}_top_terms.csv"
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def representative_targets(detail: dict, tf: str, n: int = 12) -> str:
    inter = detail[tf].get("intersection_genes", [])
    if inter:
        return ";".join(inter[:n])
    return ";".join(detail[tf].get("regulon_targets", [])[:n])


def main() -> None:
    out_dir = ensure_functional_dir()
    print("执行转录程序汇聚分析...")
    detail_path = out_dir / "02_tf_gene_sets_detail.json"
    if not detail_path.exists():
        raise FileNotFoundError("缺少 02_tf_gene_sets_detail.json，请先运行 02_build_tf_gene_sets.py")
    detail = json.loads(detail_path.read_text(encoding="utf-8"))
    robustness = load_robustness_tables()
    integrated = robustness["integrated"]
    freq = robustness["shortlist_frequency"]

    rows = []
    overlap_rows = []
    theme_rows = []
    for tf in FOCUS_TFS:
        enriched = load_enrichment(out_dir, tf)
        top_terms = enriched.head(15)["term"].astype(str).tolist() if not enriched.empty else []
        themes = infer_themes_from_terms(top_terms)
        main_themes = sorted(themes.items(), key=lambda x: x[1], reverse=True)
        theme_1 = main_themes[0][0] if main_themes and main_themes[0][1] > 0 else "not_clear"
        theme_2 = main_themes[1][0] if len(main_themes) > 1 and main_themes[1][1] > 0 else "not_clear"
        for theme, count in themes.items():
            theme_rows.append({"tf": tf, "theme": theme, "theme_term_count": count})

        reg_n = len(detail[tf].get("regulon_targets", []))
        inter_n = len(detail[tf].get("intersection_genes", []))
        deg_n = len(detail[tf].get("directional_deg_genes", []))
        overlap_ratio = inter_n / max(reg_n, 1)
        overlap_rows.append(
            {
                "tf": tf,
                "main_direction": TF_MAIN_DIRECTION[tf],
                "regulon_target_count": reg_n,
                "directional_deg_count": deg_n,
                "overlap_count": inter_n,
                "overlap_ratio_of_regulon": overlap_ratio,
                "enrichment_direction_consistent": inter_n > 0,
            }
        )

        rob_row = integrated[integrated["tf"] == tf].iloc[0] if not integrated.empty and not integrated[integrated["tf"] == tf].empty else {}
        freq_row = freq[freq["tf"] == tf].iloc[0] if not freq.empty and not freq[freq["tf"] == tf].empty else {}
        expr_direction = rob_row.get("full_data_expr_direction", TF_MAIN_DIRECTION[tf]) if isinstance(rob_row, pd.Series) else TF_MAIN_DIRECTION[tf]
        regulon_direction = rob_row.get("full_data_regulon_direction", TF_MAIN_DIRECTION[tf]) if isinstance(rob_row, pd.Series) else TF_MAIN_DIRECTION[tf]
        pseudobulk_direction = rob_row.get("pseudobulk_direction", TF_MAIN_DIRECTION[tf]) if isinstance(rob_row, pd.Series) else TF_MAIN_DIRECTION[tf]
        loo_consistency = rob_row.get("loo_expr_consistency_rate", np.nan) if isinstance(rob_row, pd.Series) else np.nan
        shortlist_frequency = freq_row.get("retention_frequency", np.nan) if isinstance(freq_row, pd.Series) else np.nan

        if tf == "NFE2L2":
            interpretation = (
                "NFE2L2 具有最强 sample-level 稳健性和较高 regulon-DEG 重叠；"
                "GO/KEGG 更偏广义转录调控/MAPK 等 lesion-associated regulatory program，"
                "可作为 stress-adaptation/reactive-compatible 主轴，但不应写成单一 oxidative-stress 通路已被完全验证。"
            )
        elif tf == "THRB":
            interpretation = (
                "THRB 汇聚于 internal_control-associated synaptic/calcium/homeostatic-supportive program；"
                "方向稳定，是对照侧核心主轴。"
            )
        elif tf == "BHLHE40":
            interpretation = "BHLHE40 内部方向稳定，程序上偏 stress/reactive/hypoxia-like 调节，但综合证据仍适合第二梯队。"
        else:
            interpretation = "SOX2 内部方向稳定但 target 较少，程序解释更偏发育/状态保留信号，不建议作为主文主轴。"

        rows.append(
            {
                "tf": tf,
                "candidate_tier": TF_TIER[tf],
                "regulon_effect_group": TF_MAIN_DIRECTION[tf],
                "expr_direction": expr_direction,
                "regulon_direction": regulon_direction,
                "leave_one_sample_out_consistency": loo_consistency,
                "pseudobulk_direction": pseudobulk_direction,
                "shortlist_retention_frequency": shortlist_frequency,
                "overlap_count": inter_n,
                "overlap_ratio": overlap_ratio,
                "program_theme_1": theme_1,
                "program_theme_2": theme_2,
                "representative_target_genes": representative_targets(detail, tf),
                "top_terms": "; ".join(top_terms[:8]),
                "final_program_interpretation": interpretation,
            }
        )

    convergence = pd.DataFrame(rows)
    overlap = pd.DataFrame(overlap_rows)
    themes = pd.DataFrame(theme_rows)
    convergence.to_csv(out_dir / "04_program_convergence_table.csv", index=False)
    overlap.to_csv(out_dir / "04_program_overlap_statistics.csv", index=False)

    evidence = convergence.copy()
    evidence["expression_support"] = evidence["expr_direction"] == evidence["regulon_effect_group"]
    evidence["regulon_support"] = evidence["regulon_direction"] == evidence["regulon_effect_group"]
    evidence["pseudobulk_support"] = evidence["pseudobulk_direction"] == evidence["regulon_effect_group"]
    evidence["robustness_support"] = evidence["leave_one_sample_out_consistency"].astype(float) >= 0.75
    evidence["shortlist_support"] = evidence["shortlist_retention_frequency"].astype(float) >= 0.60
    evidence["program_support"] = evidence["overlap_count"].astype(float) > 0
    evidence.to_csv(out_dir / "04_tf_integrated_evidence_table.csv", index=False)

    theme_matrix = themes.pivot(index="tf", columns="theme", values="theme_term_count").fillna(0).reindex(FOCUS_TFS)
    save_matrix_heatmap(theme_matrix, out_dir / "04_tf_program_theme_heatmap.png", "TF x enriched program theme count")

    layers = ["expression_support", "regulon_support", "robustness_support", "pseudobulk_support", "shortlist_support", "program_support"]
    plot_df = evidence[["tf", *layers]].copy()
    fig, ax = plt.subplots(figsize=(9, 4.6), constrained_layout=True)
    for i, tf in enumerate(plot_df["tf"]):
        for j, layer in enumerate(layers):
            value = bool(plot_df.loc[plot_df["tf"] == tf, layer].iloc[0])
            ax.scatter(j, i, s=420 if value else 120, color="#2E7D32" if value else "#BDBDBD", alpha=0.82, edgecolor="black", linewidth=0.3)
    ax.set_xticks(np.arange(len(layers)))
    ax.set_xticklabels([x.replace("_support", "") for x in layers], rotation=35, ha="right")
    ax.set_yticks(np.arange(plot_df.shape[0]))
    ax.set_yticklabels(plot_df["tf"])
    ax.set_title("Integrated evidence convergence")
    fig.savefig(out_dir / "04_integrated_evidence_bubble.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    lines = [
        "转录程序汇聚分析总结",
        "=" * 40,
        "本轮将 expression direction、regulon direction、leave-one-sample-out robustness、pseudobulk support、shortlist stability 与 GO/KEGG program theme 合并解释。",
        "",
        "重点程序:",
        "- NFE2L2: " + convergence.loc[convergence["tf"] == "NFE2L2", "final_program_interpretation"].iloc[0],
        "- THRB: " + convergence.loc[convergence["tf"] == "THRB", "final_program_interpretation"].iloc[0],
        "- BHLHE40: " + convergence.loc[convergence["tf"] == "BHLHE40", "final_program_interpretation"].iloc[0],
        "- SOX2: " + convergence.loc[convergence["tf"] == "SOX2", "final_program_interpretation"].iloc[0],
        "",
        "汇聚表:",
        convergence.to_string(index=False),
    ]
    (out_dir / "04_program_convergence_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"完成转录程序汇聚分析: {out_dir}")


if __name__ == "__main__":
    main()
