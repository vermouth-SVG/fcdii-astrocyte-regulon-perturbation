#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from functional_common import (
    CELLORACLE_CANDIDATES,
    DISCOVERY_REGULONS,
    FUNCTIONAL_DIR,
    INPUT_H5AD,
    ROBUSTNESS_DIR,
    ROOT,
    ensure_functional_dir,
)


def status_line(label: str, path: Path, usage: str) -> str:
    exists = path.exists()
    size = path.stat().st_size if exists and path.is_file() else 0
    return (
        f"[{label}]\n"
        f"路径: {path}\n"
        f"是否存在: {'是' if exists else '否'}\n"
        f"大小(bytes): {size}\n"
        f"后续用途: {usage}\n"
    )


def main() -> None:
    out_dir = ensure_functional_dir()
    print("检查功能解释分析所需资源...")
    resources = [
        ("发现队列 h5ad", INPUT_H5AD, "读取 raw/X 表达矩阵、obs 分组字段和 pySCENIC AUC 信息。"),
        ("pySCENIC regulon 定义", DISCOVERY_REGULONS, "优先从多层表头 regulons.csv 解析 TF-target genes。"),
        (
            "内部稳健性整合表",
            ROBUSTNESS_DIR / "05_robustness_validation_integrated_table.csv",
            "提供 sample-level robustness、pseudobulk 和 shortlist 稳定性证据。",
        ),
        ("候选 TF shortlist 指标表", CELLORACLE_CANDIDATES, "提供候选 TF 方向、regulon FDR、表达方向和阳性比例。"),
        (
            "外部 round2 expression-level 支持表",
            ROOT / "external_validation_round2" / "external_round2_support_ranking.csv",
            "只作为背景参考，本轮不重新做外部验证。",
        ),
        (
            "外部 round2 regulon-level 支持表",
            ROOT / "external_validation_round2" / "round2_regulon_support_ranking.csv",
            "只作为背景参考，本轮不重新做外部验证。",
        ),
    ]
    lines = [
        "GO/KEGG + 转录程序汇聚分析资源检查",
        "=" * 48,
        f"输出目录: {FUNCTIONAL_DIR}",
        "",
    ]
    for label, path, usage in resources:
        lines.append(status_line(label, path, usage))

    if DISCOVERY_REGULONS.exists():
        lines.append("regulons.csv 解析策略: 使用 pandas 读取 3 层表头，识别 TF/MotifID/AUC/NES/Annotation/Context/TargetGenes 列，按 direct annotation + activating context + NES/AUC 选择代表性 target set。")
    else:
        lines.append("regulons.csv 不存在时的备选策略: 本轮将无法构建 discovery regulon target genes，后续脚本会报出清楚错误。")

    lines.append("")
    lines.append("富集数据库策略:")
    lines.append("本地未发现 GO/KEGG gene set 文件时，03_run_go_kegg.py 会从 Enrichr 下载 GO_Biological_Process_2023 与 KEGG_2021_Human 并缓存到 D 盘。若联网失败，将保留基因集与汇聚分析结果，并在富集报告中明确说明。")

    out = out_dir / "01_resource_check_cn.txt"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"完成资源检查: {out}")


if __name__ == "__main__":
    main()
