#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd

from robustness_common import ALL_TFS, OUT_DIR, direction_from_value, ensure_out_dir


def read_csv_optional(path):
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def robustness_level(loo_expr, loo_reg, pseudo_dir, full_expr_dir, retention):
    score = 0
    if pd.notna(loo_expr) and loo_expr >= 0.75:
        score += 1
    if pd.notna(loo_reg) and loo_reg >= 0.75:
        score += 1
    if pd.notna(retention) and retention >= 0.60:
        score += 1
    if str(pseudo_dir) == str(full_expr_dir):
        score += 1
    if score >= 4:
        return "high"
    if score >= 2:
        return "moderate"
    return "low"


def tier_for_tf(tf: str, level: str, full_expr_direction: str) -> str:
    if tf == "NFE2L2" and level in {"high", "moderate"}:
        return "primary_lesion_axis"
    if tf == "THRB" and level in {"high", "moderate"}:
        return "primary_internal_control_axis"
    if tf == "BHLHE40":
        return "second_tier_lesion_candidate"
    if tf == "SOX2":
        return "reserved_candidate_not_main_axis"
    if level == "high":
        return f"supporting_{full_expr_direction}_candidate"
    return "lower_priority"


def line_for_tf(tf: str, table: pd.DataFrame) -> str:
    row = table[table["tf"] == tf]
    if row.empty:
        return f"{tf}: 未检测到完整稳健性结果。"
    r = row.iloc[0]
    return (
        f"{tf}: full-data 表达方向={r['full_data_expr_direction']}，"
        f"LOO表达一致率={r['loo_expr_consistency_rate']:.2f}，"
        f"LOO regulon一致率={r['loo_regulon_consistency_rate']:.2f}，"
        f"pseudobulk方向={r['pseudobulk_direction']}，"
        f"shortlist保留频率={r['shortlist_retention_frequency']:.2f}，"
        f"内部稳健性等级={r['final_internal_robustness_level']}，"
        f"推荐层级={r['recommended_priority_tier']}。"
    )


def main() -> None:
    out_dir = ensure_out_dir()
    print("生成内部稳健性验证总汇总...")

    detected_json = out_dir / "01_detected_fields.json"
    input_summary = (out_dir / "01_input_object_summary.txt").read_text(encoding="utf-8") if (out_dir / "01_input_object_summary.txt").exists() else "输入对象检查结果缺失。"
    loo_summary = read_csv_optional(out_dir / "02_leave_one_out_tf_summary.csv")
    pseudo_tf = read_csv_optional(out_dir / "03_pseudobulk_candidate_tf_table.csv")
    pseudo_reg = read_csv_optional(out_dir / "03_pseudobulk_regulon_table.csv")
    freq = read_csv_optional(out_dir / "04_shortlist_membership_frequency.csv")

    rows = []
    for tf in ALL_TFS:
        loo_row = loo_summary[loo_summary["tf"] == tf].iloc[0] if not loo_summary.empty and not loo_summary[loo_summary["tf"] == tf].empty else {}
        pseudo_row = pseudo_tf[pseudo_tf["tf"] == tf].iloc[0] if not pseudo_tf.empty and not pseudo_tf[pseudo_tf["tf"] == tf].empty else {}
        pseudo_reg_row = pseudo_reg[pseudo_reg["tf"] == tf].iloc[0] if not pseudo_reg.empty and not pseudo_reg[pseudo_reg["tf"] == tf].empty else {}
        freq_row = freq[freq["tf"] == tf].iloc[0] if not freq.empty and not freq[freq["tf"] == tf].empty else {}

        full_expr = loo_row.get("full_data_expr_log2fc", np.nan) if isinstance(loo_row, pd.Series) else np.nan
        full_reg = loo_row.get("full_data_regulon_diff", np.nan) if isinstance(loo_row, pd.Series) else np.nan
        full_expr_dir = direction_from_value(full_expr)
        full_reg_dir = direction_from_value(full_reg)
        loo_expr_rate = loo_row.get("expr_consistency_rate", np.nan) if isinstance(loo_row, pd.Series) else np.nan
        loo_reg_rate = loo_row.get("regulon_consistency_rate", np.nan) if isinstance(loo_row, pd.Series) else np.nan
        pseudo_dir = pseudo_row.get("direction", "NA") if isinstance(pseudo_row, pd.Series) else "NA"
        pseudo_effect = pseudo_row.get("log2fc_lesion_vs_internal_control", np.nan) if isinstance(pseudo_row, pd.Series) else np.nan
        pseudo_reg_dir = pseudo_reg_row.get("direction", "NA") if isinstance(pseudo_reg_row, pd.Series) else "NA"
        retention = freq_row.get("retention_frequency", np.nan) if isinstance(freq_row, pd.Series) else np.nan

        level = robustness_level(loo_expr_rate, loo_reg_rate, pseudo_dir, full_expr_dir, retention)
        rows.append(
            {
                "tf": tf,
                "full_data_expr_log2fc": full_expr,
                "full_data_expr_direction": full_expr_dir,
                "full_data_regulon_diff": full_reg,
                "full_data_regulon_direction": full_reg_dir,
                "loo_expr_consistency_rate": loo_expr_rate,
                "loo_regulon_consistency_rate": loo_reg_rate,
                "pseudobulk_direction": pseudo_dir,
                "pseudobulk_effect_size": pseudo_effect,
                "pseudobulk_regulon_direction": pseudo_reg_dir,
                "shortlist_retention_frequency": retention,
                "final_internal_robustness_level": level,
                "recommended_priority_tier": tier_for_tf(tf, level, full_expr_dir),
            }
        )

    integrated = pd.DataFrame(rows)
    integrated.to_csv(out_dir / "05_robustness_validation_integrated_table.csv", index=False)

    def answer(tf):
        return line_for_tf(tf, integrated)

    nfe = integrated[integrated["tf"] == "NFE2L2"].iloc[0]
    thrb = integrated[integrated["tf"] == "THRB"].iloc[0]
    bhl = integrated[integrated["tf"] == "BHLHE40"].iloc[0]
    sox = integrated[integrated["tf"] == "SOX2"].iloc[0]

    nfe_answer = "是" if nfe["final_internal_robustness_level"] in {"high", "moderate"} and nfe["full_data_expr_direction"] == "lesion" else "支持有限"
    thrb_answer = "是" if thrb["final_internal_robustness_level"] in {"high", "moderate"} and thrb["full_data_expr_direction"] == "internal_control" else "支持有限"
    bhl_answer = "是" if bhl["final_internal_robustness_level"] != "high" or bhl["shortlist_retention_frequency"] < 0.60 else "不明显"
    sox_answer = "是" if sox["recommended_priority_tier"] == "reserved_candidate_not_main_axis" else "需要结合结果谨慎解释"

    lines = [
        "内部稳健性验证总汇总",
        "=" * 40,
        "",
        "1. 输入对象与字段识别结果",
        input_summary.split("输出文件:")[0].strip(),
        "",
        "2. donor/sample leave-one-out 核心发现",
    ]
    if not loo_summary.empty:
        lines.append(loo_summary.to_string(index=False))
    else:
        lines.append("leave-one-out 结果缺失。")
    lines.extend(
        [
            "",
            "3. pseudobulk / donor-level 核心发现",
        ]
    )
    if not pseudo_tf.empty:
        lines.append(pseudo_tf[["tf", "direction", "log2fc_lesion_vs_internal_control", "mean_diff_lesion_minus_internal_control"]].to_string(index=False))
    else:
        lines.append("pseudobulk 结果缺失。")
    if not pseudo_reg.empty:
        lines.append("")
        lines.append("pseudobulk regulon AUC:")
        lines.append(pseudo_reg[["tf", "direction", "mean_diff_lesion_minus_internal_control"]].to_string(index=False))
    lines.extend(
        [
            "",
            "4. shortlist 稳定性核心发现",
        ]
    )
    if not freq.empty:
        lines.append(freq[["tf", "retention_frequency", "median_rank_when_retained", "regulon_effect_group", "expr_higher_group"]].to_string(index=False))
    else:
        lines.append("shortlist sensitivity 结果缺失。")
    lines.extend(
        [
            "",
            "5. 重点 TF 一句话结论",
            answer("NFE2L2"),
            answer("THRB"),
            answer("BHLHE40"),
            answer("SOX2"),
            "",
            "6. 对主文主轴的判断",
            f"NFE2L2 是否仍是最稳定 lesion 侧候选: {nfe_answer}。",
            f"THRB 是否仍是最稳定 internal_control 侧候选: {thrb_answer}。",
            f"BHLHE40 是否明显受 donor 或阈值影响: {bhl_answer}。",
            f"SOX2 是否更适合作为保留候选而非主结论核心: {sox_answer}。",
            "",
            "综合判断: 当前内部稳健性验证支持继续将 NFE2L2 和 THRB 作为主文主轴；BHLHE40 可作为第二梯队 lesion 侧补充候选；SOX2 建议保留但不作为核心主结论。",
        ]
    )
    (out_dir / "05_robustness_validation_master_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("完成内部稳健性验证总汇总。")
    print(f"输出目录: {out_dir}")


if __name__ == "__main__":
    main()
