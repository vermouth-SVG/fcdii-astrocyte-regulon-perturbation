#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from robustness_common import (
    ALL_TFS,
    INPUT_H5AD,
    auc_dataframe,
    benjamini_hochberg,
    detect_group_column,
    detect_sample_fields,
    differential_vector,
    direction_from_value,
    ensure_out_dir,
    get_gene_vector,
    load_discovery_h5ad_light,
    make_gene_index,
    resolve_regulon_name,
    sign_consistent,
)


def add_fdr_by_fold(df: pd.DataFrame, p_col: str, out_col: str) -> pd.DataFrame:
    out = df.copy()
    out[out_col] = np.nan
    for dropped, idx in out.groupby("dropped_donor").groups.items():
        out.loc[idx, out_col] = benjamini_hochberg(out.loc[idx, p_col].to_numpy(dtype=float))
    return out


def make_expr_rows(obs, matrix, gene_index, group_col, lesion_label, control_label, unit_col):
    expr_vectors = {tf: get_gene_vector(matrix, gene_index, tf) for tf in ALL_TFS}
    folds = [("FULL_DATA", np.ones(obs.shape[0], dtype=bool))]
    for unit in sorted(obs[unit_col].astype(str).unique().tolist()):
        folds.append((unit, obs[unit_col].astype(str).to_numpy() != unit))

    rows = []
    for dropped, mask in folds:
        sub_obs = obs.loc[mask]
        for tf in ALL_TFS:
            vec = expr_vectors[tf]
            if vec is None:
                stats = {k: np.nan for k in [
                    "n_lesion", "n_internal_control", "mean_lesion", "mean_internal_control",
                    "median_lesion", "median_internal_control", "positive_fraction_lesion",
                    "positive_fraction_internal_control", "log2fc_lesion_vs_internal_control",
                    "mean_diff_lesion_minus_internal_control", "p_value", "statistic",
                ]}
                found = False
            else:
                stats = differential_vector(vec[mask], sub_obs[group_col], lesion_label, control_label)
                found = True
            rows.append(
                {
                    "tf": tf,
                    "dropped_donor": dropped,
                    "analysis_unit_col": unit_col,
                    "gene_found": found,
                    "n_cells_remaining": int(mask.sum()),
                    "lesion_cells_remaining": stats["n_lesion"],
                    "internal_control_cells_remaining": stats["n_internal_control"],
                    "expr_mean_lesion": stats["mean_lesion"],
                    "expr_mean_internal_control": stats["mean_internal_control"],
                    "expr_positive_fraction_lesion": stats["positive_fraction_lesion"],
                    "expr_positive_fraction_internal_control": stats["positive_fraction_internal_control"],
                    "expr_log2fc": stats["log2fc_lesion_vs_internal_control"],
                    "expr_direction": direction_from_value(stats["log2fc_lesion_vs_internal_control"]),
                    "expr_p": stats["p_value"],
                    "expr_statistic": stats["statistic"],
                }
            )
    df = add_fdr_by_fold(pd.DataFrame(rows), "expr_p", "expr_fdr")
    baseline = df[df["dropped_donor"] == "FULL_DATA"].set_index("tf")["expr_log2fc"].to_dict()
    df["full_data_expr_log2fc"] = df["tf"].map(baseline)
    df["expr_direction_consistent_with_full"] = [
        True if d == "FULL_DATA" else sign_consistent(v, baseline.get(tf, np.nan))
        for tf, d, v in zip(df["tf"], df["dropped_donor"], df["expr_log2fc"])
    ]
    df["expr_significant_fdr05"] = df["expr_fdr"] < 0.05
    return df


def make_regulon_rows(obs, auc_df, regulon_names, group_col, lesion_label, control_label, unit_col):
    folds = [("FULL_DATA", np.ones(obs.shape[0], dtype=bool))]
    for unit in sorted(obs[unit_col].astype(str).unique().tolist()):
        folds.append((unit, obs[unit_col].astype(str).to_numpy() != unit))
    rows = []
    for dropped, mask in folds:
        sub_obs = obs.loc[mask]
        for tf in ALL_TFS:
            regulon = resolve_regulon_name(tf, regulon_names)
            if auc_df is None or regulon is None:
                stats = {k: np.nan for k in [
                    "n_lesion", "n_internal_control", "mean_lesion", "mean_internal_control",
                    "median_lesion", "median_internal_control", "positive_fraction_lesion",
                    "positive_fraction_internal_control", "log2fc_lesion_vs_internal_control",
                    "mean_diff_lesion_minus_internal_control", "p_value", "statistic",
                ]}
                found = False
            else:
                values = auc_df.loc[sub_obs.index.astype(str), regulon].to_numpy(dtype=float)
                stats = differential_vector(values, sub_obs[group_col], lesion_label, control_label)
                found = True
            rows.append(
                {
                    "tf": tf,
                    "regulon": regulon or "",
                    "dropped_donor": dropped,
                    "analysis_unit_col": unit_col,
                    "regulon_found": found,
                    "n_cells_remaining": int(mask.sum()),
                    "lesion_cells_remaining": stats["n_lesion"],
                    "internal_control_cells_remaining": stats["n_internal_control"],
                    "regulon_mean_lesion": stats["mean_lesion"],
                    "regulon_mean_internal_control": stats["mean_internal_control"],
                    "regulon_diff": stats["mean_diff_lesion_minus_internal_control"],
                    "regulon_direction": direction_from_value(stats["mean_diff_lesion_minus_internal_control"]),
                    "regulon_p": stats["p_value"],
                    "regulon_statistic": stats["statistic"],
                }
            )
    df = add_fdr_by_fold(pd.DataFrame(rows), "regulon_p", "regulon_fdr")
    baseline = df[df["dropped_donor"] == "FULL_DATA"].set_index("tf")["regulon_diff"].to_dict()
    df["full_data_regulon_diff"] = df["tf"].map(baseline)
    df["regulon_direction_consistent_with_full"] = [
        True if d == "FULL_DATA" else sign_consistent(v, baseline.get(tf, np.nan))
        for tf, d, v in zip(df["tf"], df["dropped_donor"], df["regulon_diff"])
    ]
    df["regulon_significant_fdr05"] = df["regulon_fdr"] < 0.05
    return df


def make_summary(expr_df: pd.DataFrame, reg_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    master = expr_df.merge(
        reg_df[
            [
                "tf", "dropped_donor", "regulon", "regulon_found", "regulon_diff", "regulon_direction",
                "regulon_p", "regulon_fdr", "full_data_regulon_diff",
                "regulon_direction_consistent_with_full", "regulon_significant_fdr05",
            ]
        ],
        on=["tf", "dropped_donor"],
        how="left",
    )
    master["overall_consistency_flag"] = np.where(
        master["regulon_found"].fillna(False),
        master["expr_direction_consistent_with_full"] & master["regulon_direction_consistent_with_full"],
        master["expr_direction_consistent_with_full"],
    )

    rows = []
    loo = master[master["dropped_donor"] != "FULL_DATA"].copy()
    for tf, sub in loo.groupby("tf", sort=False):
        full_expr = master[(master["tf"] == tf) & (master["dropped_donor"] == "FULL_DATA")]["expr_log2fc"].iloc[0]
        full_reg = master[(master["tf"] == tf) & (master["dropped_donor"] == "FULL_DATA")]["regulon_diff"].iloc[0]
        expr_rate = float(sub["expr_direction_consistent_with_full"].mean())
        reg_rate = float(sub["regulon_direction_consistent_with_full"].mean()) if sub["regulon_found"].any() else np.nan
        overall_rate = float(sub["overall_consistency_flag"].mean())
        min_abs_expr_ratio = float((sub["expr_log2fc"].abs() / (abs(full_expr) + 1e-9)).min()) if np.isfinite(full_expr) else np.nan
        min_abs_reg_ratio = float((sub["regulon_diff"].abs() / (abs(full_reg) + 1e-9)).min()) if np.isfinite(full_reg) else np.nan
        donor_sensitive = bool(
            overall_rate < 0.75
            or expr_rate < 0.75
            or (np.isfinite(reg_rate) and reg_rate < 0.75)
            or (np.isfinite(min_abs_expr_ratio) and min_abs_expr_ratio < 0.25)
            or (np.isfinite(min_abs_reg_ratio) and min_abs_reg_ratio < 0.25)
        )
        rows.append(
            {
                "tf": tf,
                "n_leave_one_out_tests": int(sub.shape[0]),
                "full_data_expr_log2fc": full_expr,
                "full_data_expr_direction": direction_from_value(full_expr),
                "full_data_regulon_diff": full_reg,
                "full_data_regulon_direction": direction_from_value(full_reg),
                "expr_consistency_count": int(sub["expr_direction_consistent_with_full"].sum()),
                "expr_consistency_rate": expr_rate,
                "expr_significant_count_fdr05": int(sub["expr_significant_fdr05"].sum()),
                "regulon_consistency_count": int(sub["regulon_direction_consistent_with_full"].sum()),
                "regulon_consistency_rate": reg_rate,
                "regulon_significant_count_fdr05": int(sub["regulon_significant_fdr05"].sum()),
                "overall_consistency_count": int(sub["overall_consistency_flag"].sum()),
                "overall_consistency_rate": overall_rate,
                "min_abs_expr_ratio_vs_full": min_abs_expr_ratio,
                "min_abs_regulon_ratio_vs_full": min_abs_reg_ratio,
                "donor_sensitive_flag": donor_sensitive,
            }
        )
    return master, pd.DataFrame(rows)


def plot_consistency_heatmap(master: pd.DataFrame, path):
    sub = master[master["dropped_donor"] != "FULL_DATA"].copy()
    pivot = sub.pivot(index="tf", columns="dropped_donor", values="overall_consistency_flag").reindex(ALL_TFS)
    values = pivot.astype(float).to_numpy()
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    im = ax.imshow(values, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns.astype(str), rotation=45, ha="right")
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index.astype(str))
    ax.set_title("Leave-one-out overall consistency")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            label = "Y" if values[i, j] == 1 else "N"
            ax.text(j, i, label, ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_forest(df: pd.DataFrame, value_col: str, full_col: str, path, title: str, xlabel: str):
    sub = df[df["dropped_donor"] != "FULL_DATA"].copy()
    fig, ax = plt.subplots(figsize=(8.5, 6), constrained_layout=True)
    y_positions = {tf: i for i, tf in enumerate(ALL_TFS)}
    units = sorted(sub["dropped_donor"].astype(str).unique().tolist())
    cmap = plt.get_cmap("tab10")
    for k, unit in enumerate(units):
        u = sub[sub["dropped_donor"] == unit]
        ax.scatter(u[value_col], [y_positions[x] + (k - len(units) / 2) * 0.06 for x in u["tf"]], s=38, label=unit, color=cmap(k % 10), alpha=0.85)
    full = df[df["dropped_donor"] == "FULL_DATA"]
    ax.scatter(full[value_col], [y_positions[x] for x in full["tf"]], marker="*", s=150, color="black", label="FULL_DATA", zorder=5)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()))
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_summary(summary: pd.DataFrame, unit_col: str, out_path):
    def status(tf):
        row = summary[summary["tf"] == tf]
        if row.empty:
            return "未检测到"
        r = row.iloc[0]
        if r["donor_sensitive_flag"]:
            return f"存在敏感性，overall consistency={r['overall_consistency_rate']:.2f}"
        return f"稳定，overall consistency={r['overall_consistency_rate']:.2f}"

    lines = [
        "donor/sample leave-one-out 稳健性验证总结",
        "=" * 45,
        f"本轮 leave-one-out 使用的分析单位: {unit_col}",
        "说明: 当前对象有 4 个 sample-level 单位，结果用于评估候选 TF 是否被某个样本单独驱动。",
        "",
        "核心 TF 结论:",
        f"- NFE2L2: {status('NFE2L2')}",
        f"- THRB: {status('THRB')}",
        f"- BHLHE40: {status('BHLHE40')}",
        f"- SOX2: {status('SOX2')}",
        "",
        "完整 TF 稳定性表:",
        summary.to_string(index=False),
        "",
        "判定规则: 任一 leave-one-out 后表达或 regulon 方向翻转、总体一致性低于 0.75、或效应量下降到 full-data 的 25% 以下，则标记为 donor/sample 敏感。",
    ]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    out_dir = ensure_out_dir()
    print("读取发现队列对象并执行 leave-one-out 稳健性验证...")
    data = load_discovery_h5ad_light(INPUT_H5AD)
    obs = data["obs"]
    group_col, lesion_label, control_label = detect_group_column(obs)
    sample_info = detect_sample_fields(obs, group_col)
    unit_col = sample_info["analysis_unit_col"]
    gene_index = make_gene_index(data["var_names"])
    auc_df = auc_dataframe(data["auc"], data["regulon_names"], obs.index)

    expr_df = make_expr_rows(obs, data["matrix"], gene_index, group_col, lesion_label, control_label, unit_col)
    reg_df = make_regulon_rows(obs, auc_df, data["regulon_names"], group_col, lesion_label, control_label, unit_col)
    master, summary = make_summary(expr_df, reg_df)

    expr_df.to_csv(out_dir / "02_leave_one_out_expr_results.csv", index=False)
    reg_df.to_csv(out_dir / "02_leave_one_out_regulon_results.csv", index=False)
    master.to_csv(out_dir / "02_leave_one_out_master_results.csv", index=False)
    summary.to_csv(out_dir / "02_leave_one_out_tf_summary.csv", index=False)
    master[master["dropped_donor"] == "FULL_DATA"].to_csv(out_dir / "02_full_data_baseline_results.csv", index=False)

    plot_consistency_heatmap(master, out_dir / "02_leave_one_out_consistency_heatmap.png")
    plot_forest(expr_df, "expr_log2fc", "full_data_expr_log2fc", out_dir / "02_leave_one_out_expr_forest.png", "Leave-one-out TF expression log2FC", "log2FC lesion vs internal_control")
    plot_forest(reg_df, "regulon_diff", "full_data_regulon_diff", out_dir / "02_leave_one_out_regulon_forest.png", "Leave-one-out regulon AUC difference", "mean AUC difference lesion - internal_control")
    write_summary(summary, unit_col, out_dir / "02_leave_one_out_summary_cn.txt")
    print("完成 leave-one-out 稳健性验证。")
    print(f"输出目录: {out_dir}")


if __name__ == "__main__":
    main()
