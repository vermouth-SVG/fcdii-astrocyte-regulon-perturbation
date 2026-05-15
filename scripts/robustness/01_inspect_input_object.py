#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd

from robustness_common import (
    ALL_TFS,
    INPUT_H5AD,
    auc_dataframe,
    detect_group_column,
    detect_sample_fields,
    ensure_out_dir,
    load_discovery_h5ad_light,
    resolve_regulon_name,
    write_json,
)


def preview_counts(series: pd.Series, n: int = 8) -> str:
    counts = series.astype("string").fillna("NA").value_counts(dropna=False).head(n)
    return " | ".join([f"{idx}:{int(val)}" for idx, val in counts.items()])


def main() -> None:
    out_dir = ensure_out_dir()
    print("读取发现队列 h5ad 对象...")
    data = load_discovery_h5ad_light(INPUT_H5AD)
    obs = data["obs"]
    var = data["var"]
    group_col, lesion_label, control_label = detect_group_column(obs)
    sample_info = detect_sample_fields(obs, group_col)

    auc_df = auc_dataframe(data["auc"], data["regulon_names"], obs.index)
    regulon_hits = []
    if auc_df is not None:
        for tf in ALL_TFS:
            regulon_hits.append({"tf": tf, "regulon_name": resolve_regulon_name(tf, data["regulon_names"])})

    obs_rows = []
    for col in obs.columns:
        values = obs[col]
        obs_rows.append(
            {
                "column": col,
                "dtype": str(values.dtype),
                "n_unique": int(values.astype("string").nunique(dropna=True)),
                "n_missing": int(values.isna().sum()),
                "top_values": preview_counts(values),
                "is_detected_group": col == group_col,
                "is_detected_donor": col == sample_info["donor_col"],
                "is_detected_sample": col == sample_info["sample_col"],
                "is_analysis_unit": col == sample_info["analysis_unit_col"],
            }
        )
    obs_columns = pd.DataFrame(obs_rows)
    obs_columns.to_csv(out_dir / "01_obs_columns.csv", index=False)

    dist = sample_info["unit_group_distribution"].copy()
    dist.index.name = sample_info["analysis_unit_col"]
    dist.reset_index().to_csv(out_dir / "01_group_sample_distribution.csv", index=False)

    missing_messages = []
    for unit, row in sample_info["unit_group_distribution"].iterrows():
        missing = [label for label in [lesion_label, control_label] if label not in row.index or int(row.get(label, 0)) == 0]
        if missing:
            missing_messages.append(f"{sample_info['analysis_unit_col']}={unit} 缺少分组: {', '.join(missing)}")

    donor_distribution = None
    if sample_info["donor_col"]:
        donor_distribution = pd.crosstab(obs[sample_info["donor_col"]].astype(str), obs[group_col].astype(str))
        donor_distribution.index.name = sample_info["donor_col"]
        donor_distribution.reset_index().to_csv(out_dir / "01_donor_group_distribution.csv", index=False)
    sample_distribution = None
    if sample_info["sample_col"]:
        sample_distribution = pd.crosstab(obs[sample_info["sample_col"]].astype(str), obs[group_col].astype(str))
        sample_distribution.index.name = sample_info["sample_col"]
        sample_distribution.reset_index().to_csv(out_dir / "01_sample_group_distribution.csv", index=False)

    detected = {
        "input_h5ad": str(INPUT_H5AD),
        "n_cells": int(obs.shape[0]),
        "n_genes": int(len(data["var_names"])),
        "matrix_source": data["matrix_source"],
        "group_col": group_col,
        "lesion_label": lesion_label,
        "control_label": control_label,
        "donor_col": sample_info["donor_col"],
        "sample_col": sample_info["sample_col"],
        "analysis_unit_col": sample_info["analysis_unit_col"],
        "donor_n": sample_info["donor_n"],
        "sample_n": sample_info["sample_n"],
        "obsm_keys": data["obsm_keys"],
        "uns_keys": data["uns_keys"],
        "layer_keys": data["layer_keys"],
        "auc_present": auc_df is not None,
        "n_regulons": 0 if auc_df is None else int(auc_df.shape[1]),
        "candidate_regulon_hits": regulon_hits,
    }
    write_json(out_dir / "01_detected_fields.json", detected)

    lines = []
    lines.append("内部稳健性验证 - 输入对象检查")
    lines.append("=" * 40)
    lines.append(f"输入对象: {INPUT_H5AD}")
    lines.append(f"维度: {obs.shape[0]} cells x {len(data['var_names'])} genes")
    lines.append(f"表达矩阵来源: {data['matrix_source']}")
    lines.append("")
    lines.append("obs 字段:")
    lines.append(", ".join(obs.columns.astype(str).tolist()))
    lines.append("")
    lines.append("var 字段:")
    lines.append(", ".join(var.columns.astype(str).tolist()))
    lines.append("var_names 前 20 个:")
    lines.append(", ".join(data["var_names"].astype(str).tolist()[:20]))
    lines.append("")
    lines.append(f"自动识别分组字段: {group_col}")
    lines.append(f"lesion 标签: {lesion_label}")
    lines.append(f"internal_control 标签: {control_label}")
    lines.append(f"分组计数: {preview_counts(obs[group_col])}")
    lines.append("")
    lines.append(f"donor 字段: {sample_info['donor_col']}，唯一值数: {sample_info['donor_n']}")
    lines.append(f"sample 字段: {sample_info['sample_col']}，唯一值数: {sample_info['sample_n']}")
    lines.append(f"本轮 leave-one-out / pseudobulk 分析单位: {sample_info['analysis_unit_col']}")
    lines.append(f"字段选择说明: {sample_info['unit_reason']}")
    lines.append("")
    lines.append("analysis unit x group 细胞数:")
    lines.append(sample_info["unit_group_distribution"].to_string())
    if donor_distribution is not None:
        lines.append("")
        lines.append("donor x group 细胞数:")
        lines.append(donor_distribution.to_string())
    if sample_distribution is not None:
        lines.append("")
        lines.append("sample x group 细胞数:")
        lines.append(sample_distribution.to_string())
    lines.append("")
    if missing_messages:
        lines.append("提示: 以下分析单位缺少某个分组，说明该字段更接近 sample-level，而不是同一 donor 内成对两组。")
        lines.extend(missing_messages)
    else:
        lines.append("未发现分析单位缺少 lesion/internal_control 分组。")
    lines.append("")
    lines.append(f"obsm keys: {', '.join(data['obsm_keys']) if data['obsm_keys'] else '无'}")
    lines.append(f"uns keys: {', '.join(data['uns_keys']) if data['uns_keys'] else '无'}")
    lines.append(f"layers keys: {', '.join(data['layer_keys']) if data['layer_keys'] else '无'}")
    lines.append(f"pySCENIC AUC 是否存在: {'是' if auc_df is not None else '否'}")
    lines.append(f"regulon 数量: {0 if auc_df is None else auc_df.shape[1]}")
    if regulon_hits:
        lines.append("候选 TF regulon 匹配:")
        for item in regulon_hits:
            lines.append(f"- {item['tf']}: {item['regulon_name'] or '未找到'}")
    lines.append("")
    lines.append("输出文件:")
    lines.append("- 01_input_object_summary.txt")
    lines.append("- 01_obs_columns.csv")
    lines.append("- 01_group_sample_distribution.csv")
    lines.append("- 01_detected_fields.json")

    (out_dir / "01_input_object_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("完成输入对象检查。")
    print(f"输出目录: {out_dir}")


if __name__ == "__main__":
    main()
