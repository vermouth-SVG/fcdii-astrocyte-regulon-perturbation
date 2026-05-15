#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd

from functional_common import FOCUS_TFS, ensure_functional_dir


def optional_text(path):
    return path.read_text(encoding="utf-8") if path.exists() else "文件缺失。"


def main() -> None:
    out_dir = ensure_functional_dir()
    print("生成功能解释主汇总...")
    resource = optional_text(out_dir / "01_resource_check_cn.txt")
    naming = optional_text(out_dir / "00_naming_correction_note_cn.txt")
    geneset = optional_text(out_dir / "02_tf_gene_sets_summary_cn.txt")
    enrichment = optional_text(out_dir / "03_functional_enrichment_summary_cn.txt")
    convergence_txt = optional_text(out_dir / "04_program_convergence_summary_cn.txt")

    convergence = pd.read_csv(out_dir / "04_program_convergence_table.csv") if (out_dir / "04_program_convergence_table.csv").exists() else pd.DataFrame()
    go = pd.read_csv(out_dir / "03_go_enrichment_all.csv") if (out_dir / "03_go_enrichment_all.csv").exists() else pd.DataFrame()
    kegg = pd.read_csv(out_dir / "03_kegg_enrichment_all.csv") if (out_dir / "03_kegg_enrichment_all.csv").exists() else pd.DataFrame()

    rows = []
    for tf in FOCUS_TFS:
        conv = convergence[convergence["tf"] == tf].iloc[0] if not convergence.empty and not convergence[convergence["tf"] == tf].empty else {}
        top_go = "; ".join(go[go["tf"] == tf].head(5)["term"].astype(str).tolist()) if not go.empty else ""
        top_kegg = "; ".join(kegg[kegg["tf"] == tf].head(5)["term"].astype(str).tolist()) if not kegg.empty else ""
        rows.append(
            {
                "tf": tf,
                "candidate_tier": conv.get("candidate_tier", ""),
                "regulon_effect_group": conv.get("regulon_effect_group", ""),
                "expr_direction": conv.get("expr_direction", ""),
                "regulon_direction": conv.get("regulon_direction", ""),
                "sample_level_robustness": conv.get("leave_one_sample_out_consistency", ""),
                "pseudobulk_support": conv.get("pseudobulk_direction", ""),
                "shortlist_stability": conv.get("shortlist_retention_frequency", ""),
                "top_go_terms": top_go,
                "top_kegg_terms": top_kegg,
                "program_theme_summary": f"{conv.get('program_theme_1', '')}; {conv.get('program_theme_2', '')}",
                "representative_targets": conv.get("representative_target_genes", ""),
                "final_functional_priority": conv.get("final_program_interpretation", ""),
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(out_dir / "05_functional_interpretation_integrated_table.csv", index=False)

    def interp(tf):
        row = table[table["tf"] == tf]
        return row.iloc[0]["final_functional_priority"] if not row.empty else "无结果。"

    lines = [
        "GO/KEGG + 转录程序汇聚分析主汇总",
        "=" * 52,
        "",
        "1. 资源检查结果",
        resource,
        "",
        "2. 术语修正说明",
        naming,
        "",
        "3. 每个重点 TF 的基因集构建情况",
        geneset,
        "",
        "4. GO/KEGG 核心结论",
        enrichment,
        "",
        "5. 转录程序汇聚核心结论",
        convergence_txt,
        "",
        "6. 四个 TF 的一句话功能解释",
        f"- NFE2L2: {interp('NFE2L2')}",
        f"- THRB: {interp('THRB')}",
        f"- BHLHE40: {interp('BHLHE40')}",
        f"- SOX2: {interp('SOX2')}",
        "",
        "7. 明确回答",
        "NFE2L2 是否可概括为 lesion-associated stress/reactive program 核心候选: 基本可以，但更稳妥的写法是 lesion-associated transcriptional/stress-adaptation program；GO/KEGG 未把结果单独收敛到纯 oxidative-stress 通路。",
        "THRB 是否可概括为 internal_control-associated homeostatic/supportive program 核心候选: 是，尤其体现为 synaptic/calcium/homeostatic-supportive program。",
        "BHLHE40 的程序层定位是否支持其作为第二梯队: 是，内部稳定且程序上有补充解释，但综合优先级低于 NFE2L2/THRB。",
        "SOX2 是否更适合作为保留候选而非主轴: 是，内部方向稳定但 target 较少，程序解释不宜作为主线核心。",
        "当前是否具备进入下一步轻量药物重定位条件: 是。已有 TF 优先级、target gene set、GO/KEGG program theme 与 sample-level robustness，可进入轻量药物重定位；本轮未实际执行药物重定位。",
    ]
    (out_dir / "05_functional_interpretation_master_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"完成主汇总: {out_dir / '05_functional_interpretation_master_summary_cn.txt'}")


if __name__ == "__main__":
    main()
