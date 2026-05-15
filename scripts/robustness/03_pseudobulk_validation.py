#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from robustness_common import (
    ALL_TFS,
    INPUT_H5AD,
    auc_dataframe,
    detect_group_column,
    detect_sample_fields,
    differential_vector,
    direction_from_value,
    discovery_targets_for_tfs,
    ensure_out_dir,
    get_gene_vector,
    load_discovery_h5ad_light,
    make_gene_index,
    resolve_regulon_name,
    save_heatmap,
    zscore_rows,
)


def pseudobulk_expression(obs, matrix, var_names, group_col, unit_col) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    meta = []
    unit_values = sorted(obs[unit_col].astype(str).unique().tolist())
    for unit in unit_values:
        mask = obs[unit_col].astype(str).to_numpy() == unit
        sub_obs = obs.loc[mask]
        group_counts = sub_obs[group_col].astype(str).value_counts()
        group_label = group_counts.index[0] if not group_counts.empty else "NA"
        summed = np.asarray(matrix[mask, :].sum(axis=0)).ravel()
        mean_expr = summed / max(int(mask.sum()), 1)
        rows.append(mean_expr)
        meta.append(
            {
                "analysis_unit": unit,
                "group": group_label,
                "n_cells": int(mask.sum()),
                "group_counts": ";".join([f"{k}:{int(v)}" for k, v in group_counts.items()]),
            }
        )
    expr = pd.DataFrame(rows, index=[m["analysis_unit"] for m in meta], columns=var_names.astype(str))
    meta_df = pd.DataFrame(meta).set_index("analysis_unit")
    return expr, meta_df


def compare_units(values_by_unit: pd.Series, meta_df: pd.DataFrame, lesion_label: str, control_label: str) -> dict:
    groups = meta_df.loc[values_by_unit.index, "group"].astype(str)
    lesion_values = values_by_unit[groups == lesion_label].astype(float).to_numpy()
    control_values = values_by_unit[groups == control_label].astype(float).to_numpy()
    if lesion_values.size and control_values.size:
        log2fc = float(np.log2((lesion_values.mean() + 1e-9) / (control_values.mean() + 1e-9)))
        diff = float(lesion_values.mean() - control_values.mean())
    else:
        log2fc = np.nan
        diff = np.nan
    if lesion_values.size >= 2 and control_values.size >= 2:
        stat, p = stats.mannwhitneyu(lesion_values, control_values, alternative="two-sided")
        stat = float(stat)
        p = float(p)
    else:
        stat, p = np.nan, np.nan
    return {
        "n_units_lesion": int(lesion_values.size),
        "n_units_internal_control": int(control_values.size),
        "mean_lesion": float(np.mean(lesion_values)) if lesion_values.size else np.nan,
        "mean_internal_control": float(np.mean(control_values)) if control_values.size else np.nan,
        "median_lesion": float(np.median(lesion_values)) if lesion_values.size else np.nan,
        "median_internal_control": float(np.median(control_values)) if control_values.size else np.nan,
        "log2fc_lesion_vs_internal_control": log2fc,
        "mean_diff_lesion_minus_internal_control": diff,
        "direction": direction_from_value(log2fc),
        "mannwhitney_u": stat,
        "p_value_descriptive_only": p,
        "lesion_unit_values": ";".join([f"{x:.6g}" for x in lesion_values]),
        "internal_control_unit_values": ";".join([f"{x:.6g}" for x in control_values]),
    }


def build_candidate_tf_table(pb_expr: pd.DataFrame, meta_df: pd.DataFrame, lesion_label: str, control_label: str) -> pd.DataFrame:
    rows = []
    for tf in ALL_TFS:
        if tf not in pb_expr.columns:
            rows.append({"tf": tf, "gene_found": False})
            continue
        row = {"tf": tf, "gene_found": True}
        row.update(compare_units(pb_expr[tf], meta_df, lesion_label, control_label))
        rows.append(row)
    return pd.DataFrame(rows)


def build_regulon_tables(obs, auc_df, regulon_names, meta_df, group_col, unit_col, lesion_label, control_label):
    if auc_df is None:
        return pd.DataFrame(), pd.DataFrame()
    unit_rows = []
    for unit in meta_df.index.astype(str):
        mask = obs[unit_col].astype(str).to_numpy() == unit
        row = {"analysis_unit": unit, "group": meta_df.loc[unit, "group"], "n_cells": int(mask.sum())}
        for tf in ALL_TFS:
            regulon = resolve_regulon_name(tf, regulon_names)
            if regulon is None:
                continue
            row[tf] = float(auc_df.loc[obs.index[mask].astype(str), regulon].mean())
        unit_rows.append(row)
    unit_df = pd.DataFrame(unit_rows).set_index("analysis_unit")
    stats_rows = []
    for tf in ALL_TFS:
        if tf not in unit_df.columns:
            stats_rows.append({"tf": tf, "regulon_found": False})
            continue
        row = {"tf": tf, "regulon": resolve_regulon_name(tf, regulon_names), "regulon_found": True}
        row.update(compare_units(unit_df[tf], meta_df, lesion_label, control_label))
        stats_rows.append(row)
    return unit_df, pd.DataFrame(stats_rows)


def plot_unit_points(table: pd.DataFrame, value_prefix: str, path, title: str, ylabel: str):
    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)
    x = np.arange(table.shape[0])
    width = 0.28
    ax.bar(x - width / 2, table["mean_internal_control"], width=width, color="#4C78A8", alpha=0.65, label="internal_control mean")
    ax.bar(x + width / 2, table["mean_lesion"], width=width, color="#E45756", alpha=0.65, label="lesion mean")
    for i, row in table.reset_index(drop=True).iterrows():
        control_vals = [float(v) for v in str(row.get("internal_control_unit_values", "")).split(";") if v]
        lesion_vals = [float(v) for v in str(row.get("lesion_unit_values", "")).split(";") if v]
        ax.scatter(np.full(len(control_vals), i - width / 2), control_vals, color="#1F4E79", s=45, zorder=3)
        ax.scatter(np.full(len(lesion_vals), i + width / 2), lesion_vals, color="#B22222", s=45, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(table["tf"].astype(str), rotation=30, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_summary(tf_table: pd.DataFrame, reg_table: pd.DataFrame, unit_col: str, out_path):
    def one(tf, table, col="direction"):
        row = table[table["tf"] == tf]
        if row.empty:
            return "未检测到"
        r = row.iloc[0]
        return f"{r.get(col, 'NA')}，log2FC={r.get('log2fc_lesion_vs_internal_control', np.nan):.3f}，diff={r.get('mean_diff_lesion_minus_internal_control', np.nan):.4g}"

    lines = [
        "pseudobulk / donor-level supporting evidence 总结",
        "=" * 50,
        f"本轮以 {unit_col} 为 donor/sample-level 分析单位，对每个单位内细胞表达取平均，作为描述性 pseudobulk 信号。",
        "由于当前 discovery pilot 约为 4 个 sample-level 单位，本部分不假装进行复杂组学统计，重点解释方向一致性、效应量和样本级支持。",
        "",
        "候选 TF 表达层面:",
        f"- NFE2L2: {one('NFE2L2', tf_table)}",
        f"- THRB: {one('THRB', tf_table)}",
        f"- BHLHE40: {one('BHLHE40', tf_table)}",
        f"- SOX2: {one('SOX2', tf_table)}",
    ]
    if not reg_table.empty:
        lines.extend(
            [
                "",
                "regulon AUC donor/sample-level 层面:",
                f"- NFE2L2: {one('NFE2L2', reg_table)}",
                f"- THRB: {one('THRB', reg_table)}",
                f"- BHLHE40: {one('BHLHE40', reg_table)}",
                f"- SOX2: {one('SOX2', reg_table)}",
            ]
        )
    lines.append("")
    lines.append("结论: 该部分应作为 donor-level supporting evidence，而不是替代 full cell-level 检验。")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    out_dir = ensure_out_dir()
    print("读取发现队列并执行 pseudobulk / donor-level 验证...")
    data = load_discovery_h5ad_light(INPUT_H5AD)
    obs = data["obs"]
    group_col, lesion_label, control_label = detect_group_column(obs)
    sample_info = detect_sample_fields(obs, group_col)
    unit_col = sample_info["analysis_unit_col"]

    pb_expr, meta_df = pseudobulk_expression(obs, data["matrix"], data["var_names"], group_col, unit_col)
    pb_out = pd.concat([meta_df, pb_expr], axis=1)
    pb_out.to_csv(out_dir / "03_pseudobulk_expression_matrix.csv")
    meta_df.to_csv(out_dir / "03_pseudobulk_sample_metadata.csv")

    gene_index = make_gene_index(data["var_names"])
    targets = discovery_targets_for_tfs(ALL_TFS)
    target_genes = []
    for tf in ALL_TFS:
        if tf in gene_index and tf not in target_genes:
            target_genes.append(tf)
        for gene in targets.get(tf, []):
            if (gene in gene_index or gene.upper() in gene_index) and gene not in target_genes:
                target_genes.append(gene)
    target_cols = [gene for gene in target_genes if gene in pb_expr.columns]
    if target_cols:
        pd.concat([meta_df, pb_expr[target_cols]], axis=1).to_csv(out_dir / "03_pseudobulk_candidate_target_expression_matrix.csv")

    tf_table = build_candidate_tf_table(pb_expr, meta_df, lesion_label, control_label)
    tf_table.to_csv(out_dir / "03_pseudobulk_candidate_tf_table.csv", index=False)

    auc_df = auc_dataframe(data["auc"], data["regulon_names"], obs.index)
    reg_unit_df, reg_table = build_regulon_tables(obs, auc_df, data["regulon_names"], meta_df, group_col, unit_col, lesion_label, control_label)
    if not reg_unit_df.empty:
        reg_unit_df.to_csv(out_dir / "03_pseudobulk_regulon_by_donor.csv")
    reg_table.to_csv(out_dir / "03_pseudobulk_regulon_table.csv", index=False)

    plot_unit_points(tf_table[tf_table["gene_found"] == True], "expr", out_dir / "03_pseudobulk_candidate_tf_plot.png", "Donor/sample-level candidate TF expression", "Mean expression per cell")
    if not reg_table.empty:
        plot_unit_points(reg_table[reg_table["regulon_found"] == True], "regulon", out_dir / "03_pseudobulk_regulon_plot.png", "Donor/sample-level regulon activity", "Mean regulon AUC")

    ordered_units = meta_df.sort_values(["group"]).index.tolist()
    heatmap_genes = [tf for tf in ALL_TFS if tf in pb_expr.columns]
    heatmap_matrix = zscore_rows(pb_expr.loc[ordered_units, heatmap_genes].T)
    save_heatmap(heatmap_matrix, out_dir / "03_pseudobulk_heatmap.png", "Pseudobulk candidate TF expression (row z-score)")

    write_summary(tf_table, reg_table, unit_col, out_dir / "03_pseudobulk_summary_cn.txt")
    print("完成 pseudobulk / donor-level 验证。")
    print(f"输出目录: {out_dir}")


if __name__ == "__main__":
    main()
