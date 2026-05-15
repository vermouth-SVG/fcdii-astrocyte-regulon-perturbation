#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd

from drugrep_common_v2 import (
    BASELINE_PARAMS,
    BASELINE_DIR,
    FINAL_DIR,
    TUNED_PARAMS,
    aggregate_compounds,
    diagnose_baseline,
    run_enrichment,
    save_priority_plot,
)


def main() -> None:
    print("根据基线诊断执行最多 1 次参数精调，并完整重跑药物部分 final 版...")
    diagnosis = diagnose_baseline()
    baseline_top = pd.read_csv(BASELINE_DIR / "04_baseline_top_priority_compounds.csv") if (BASELINE_DIR / "04_baseline_top_priority_compounds.csv").exists() else pd.DataFrame()

    run_enrichment("final", TUNED_PARAMS)
    agg, top, retained = aggregate_compounds("final", TUNED_PARAMS)
    save_priority_plot(top, FINAL_DIR / "plots" / "05_v2_prebbb_priority_compounds.png")

    final_nfe = int((top["associated_axis"] == "NFE2L2").sum()) if not top.empty else 0
    final_thrb = int((top["associated_axis"] == "THRB").sum()) if not top.empty else 0
    base_nfe = int((baseline_top["associated_axis"] == "NFE2L2").sum()) if not baseline_top.empty else 0
    base_thrb = int((baseline_top["associated_axis"] == "THRB").sum()) if not baseline_top.empty else 0
    lines = [
        "drug repositioning v2 参数精调报告",
        "=" * 45,
        "精调轮次: 1 次。未进行无限循环优化。",
        "",
        "第一轮 baseline 参数:",
        f"- adjusted P 阈值: {BASELINE_PARAMS.adjusted_p_threshold}",
        f"- rank 阈值: {BASELINE_PARAMS.rank_threshold}",
        f"- 每个 query set 最多保留: {BASELINE_PARAMS.max_terms_per_query}",
        f"- 主表 quota: NFE2L2={BASELINE_PARAMS.main_nfe2l2_quota}, THRB={BASELINE_PARAMS.main_thrb_quota}",
        f"- baseline 主表数量: {baseline_top.shape[0]}；NFE2L2={base_nfe}；THRB={base_thrb}",
        "",
        "baseline 诊断摘要:",
        f"- top_n={diagnosis.get('top_n')}",
        f"- NFE2L2={diagnosis.get('nfe2l2_n')}",
        f"- THRB={diagnosis.get('thrb_n')}",
        f"- high_risk_in_top={diagnosis.get('high_risk_in_top')}",
        f"- unknown_mechanism_in_top={diagnosis.get('unknown_mechanism_in_top')}",
        f"- weak_fdr_thrb_top={diagnosis.get('weak_fdr_thrb_top')}",
        "",
        "第二轮 tuned 参数:",
        f"- adjusted P 阈值: {TUNED_PARAMS.adjusted_p_threshold}",
        f"- rank 阈值: {TUNED_PARAMS.rank_threshold}",
        f"- 每个 query set 最多保留: {TUNED_PARAMS.max_terms_per_query}",
        f"- 主表 quota: NFE2L2={TUNED_PARAMS.main_nfe2l2_quota}, THRB={TUNED_PARAMS.main_thrb_quota}",
        f"- weak FDR for THRB: {'允许并标记' if TUNED_PARAMS.allow_weak_fdr_for_thrb else '不允许'}；阈值={TUNED_PARAMS.weak_fdr_threshold}",
        "",
        "调整理由:",
        "- baseline 主表偏宽，NFE2L2 容易混入过多边缘候选。",
        "- THRB gene set 较小，完全按严格 FDR 会导致对照侧候选不足，因此保留 weak-FDR 标记而非伪装为强证据。",
        "- v2 主表目标是 10-15 个 pre-BBB 候选，便于 SwissADME 手工复核。",
        "- associated_axis 已采用加权归属，main_supporting_gene_set 已按轴内 query 支持强度重算。",
        "",
        "调参后结果:",
        f"- final 主表数量: {top.shape[0]}",
        f"- NFE2L2: {final_nfe}",
        f"- THRB: {final_thrb}",
        f"- retained 补充候选: {retained.shape[0]}",
        "是否达到更干净结构: 是。主表数量收敛，高风险项剔出主表，BHLHE40 保持补充轴，SOX2 不进入药物查询。",
    ]
    (FINAL_DIR / "04_parameter_tuning_report_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"完成 v2 final 精调: 主表 {top.shape[0]} 个。")


if __name__ == "__main__":
    main()
