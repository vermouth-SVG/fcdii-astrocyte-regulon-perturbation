#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd

from functional_common import (
    DEG_LOG2FC_THRESHOLD,
    DEG_POSITIVE_FRACTION_THRESHOLD,
    FOCUS_TFS,
    INPUT_H5AD,
    TF_MAIN_DIRECTION,
    TF_TIER,
    build_directional_deg_sets,
    compute_light_deg,
    detect_group_column,
    discovery_targets_for_tfs,
    ensure_functional_dir,
    load_discovery_h5ad_light,
    write_json,
)


def main() -> None:
    out_dir = ensure_functional_dir()
    print("构建 NFE2L2/THRB/BHLHE40/SOX2 的 regulon target 与方向一致 DEG 基因集...")
    data = load_discovery_h5ad_light(INPUT_H5AD)
    obs = data["obs"]
    group_col, lesion_label, control_label = detect_group_column(obs)
    deg = compute_light_deg(data, group_col, lesion_label, control_label)
    deg.to_csv(out_dir / "02_discovery_light_deg_table.csv", index=False)
    directional_sets = build_directional_deg_sets(deg)
    targets = discovery_targets_for_tfs(FOCUS_TFS)

    detail = {}
    rows = []
    deg_map = deg.set_index("gene")
    for tf in FOCUS_TFS:
        direction = TF_MAIN_DIRECTION[tf]
        regulon_targets = [g for g in targets.get(tf, []) if g in set(deg["gene"])]
        direction_genes = sorted(directional_sets[direction])
        intersection = sorted(set(regulon_targets) & set(direction_genes))
        ratio = len(intersection) / max(len(regulon_targets), 1)
        if len(regulon_targets) < 15:
            size_note = "regulon target 较小，富集结果需谨慎"
        elif len(regulon_targets) > 800:
            size_note = "regulon target 较大，富集结果可能偏向广义程序"
        else:
            size_note = "基因集大小适中"
        top_intersection = []
        for gene in intersection:
            if gene in deg_map.index:
                row = deg_map.loc[gene]
                top_intersection.append(
                    {
                        "gene": gene,
                        "log2fc": float(row["log2fc_lesion_vs_internal_control"]),
                        "positive_fraction_higher_group": float(row["positive_fraction_higher_group"]),
                    }
                )
        top_intersection = sorted(top_intersection, key=lambda x: abs(x["log2fc"]), reverse=True)
        explanation = (
            f"{tf} 主方向为 {direction}；discovery regulon target={len(regulon_targets)}，"
            f"方向一致 DEG 交集={len(intersection)}，交集比例={ratio:.3f}。{size_note}。"
        )
        rows.append(
            {
                "tf": tf,
                "candidate_tier": TF_TIER[tf],
                "main_direction": direction,
                "regulon_target_count": len(regulon_targets),
                "directional_deg_count": len(direction_genes),
                "intersection_count": len(intersection),
                "intersection_ratio_of_regulon": ratio,
                "deg_log2fc_threshold": DEG_LOG2FC_THRESHOLD,
                "deg_positive_fraction_threshold": DEG_POSITIVE_FRACTION_THRESHOLD,
                "size_note": size_note,
                "summary_cn": explanation,
            }
        )
        detail[tf] = {
            "candidate_tier": TF_TIER[tf],
            "main_direction": direction,
            "regulon_targets": regulon_targets,
            "directional_deg_genes": direction_genes,
            "intersection_genes": intersection,
            "top_intersection_genes_by_abs_log2fc": top_intersection[:50],
            "summary_cn": explanation,
        }

    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "02_tf_gene_sets_summary.csv", index=False)
    write_json(out_dir / "02_tf_gene_sets_detail.json", detail)

    lines = [
        "TF gene set 构建总结",
        "=" * 36,
        f"表达对象: {INPUT_H5AD}",
        f"分组字段: {group_col}; lesion={lesion_label}; internal_control={control_label}",
        f"轻量 DEG 定义: |log2FC| >= {DEG_LOG2FC_THRESHOLD}, higher-group positive fraction >= {DEG_POSITIVE_FRACTION_THRESHOLD}",
        "",
    ]
    lines.extend(summary["summary_cn"].tolist())
    lines.append("")
    lines.append("说明: intersection gene set 表示 discovery regulon target 与主方向一致 DEG 的交集，是后续 program convergence 的核心输入。")
    (out_dir / "02_tf_gene_sets_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"完成 gene set 构建: {out_dir / '02_tf_gene_sets_summary.csv'}")


if __name__ == "__main__":
    main()
