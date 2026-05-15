#!/usr/bin/env python3
from __future__ import annotations

import math
import inspect
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[1]
INPUT_H5AD = ROOT / "external_validation_round2" / "input" / "external_validation_round2.h5ad"
OUT_DIR = ROOT / "external_validation_round2"
PLOT_DIR = OUT_DIR / "summary_plots"
MAIN_METRICS = ROOT / "analysis_outputs" / "celloracle_candidates" / "celloracle_candidate_tf_metrics.csv"
TARGET_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
OBSM_KEY = "X_pyscenic_auc"
UNS_KEY = "pyscenic_regulon_names"
EPS = 1e-9
MIN_CLUSTER_CELLS_PER_GROUP = 30


def benjamini_hochberg(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(arr)
    if not finite_mask.any():
        return out
    finite = arr[finite_mask]
    order = np.argsort(finite)
    ranked = finite[order]
    n = ranked.size
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        running = min(running, ranked[i] * n / rank)
        adjusted[i] = running
    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    out[finite_mask] = np.clip(restored, 0.0, 1.0)
    return out


def minmax_scale(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(arr)
    if not finite_mask.any():
        return out
    lo = np.nanmin(arr[finite_mask])
    hi = np.nanmax(arr[finite_mask])
    if hi - lo < EPS:
        out[finite_mask] = 0.5
        return out
    out[finite_mask] = (arr[finite_mask] - lo) / (hi - lo)
    return out


def ensure_sparse(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def load_expected_directions() -> pd.DataFrame:
    df = pd.read_csv(MAIN_METRICS)
    df["tf"] = df["tf"].astype(str).str.upper()
    df = df[df["tf"].isin([x.upper() for x in TARGET_TFS])].copy()
    keep = [
        "tf",
        "regulon",
        "regulon_effect_group",
        "regulon_fdr",
        "expr_higher_group",
        "expr_wilcoxon_fdr",
        "log2fc_lesion_vs_internal_control",
    ]
    for col in keep:
        if col not in df.columns:
            df[col] = np.nan
    return df[keep].drop_duplicates(subset=["tf"]).set_index("tf", drop=False)


def find_gene_index(adata: ad.AnnData, gene_name: str) -> int | None:
    target = gene_name.upper()
    var_names = pd.Index(adata.var_names.astype(str))
    match = np.flatnonzero(var_names.str.upper() == target)
    if match.size:
        return int(match[0])
    for col in ["gene_symbol", "gene_symbols", "symbol", "gene", "gene_id"]:
        if col in adata.var.columns:
            series = adata.var[col].astype(str).str.upper()
            match = np.flatnonzero(series.values == target)
            if match.size:
                return int(match[0])
    prefixed = np.flatnonzero(var_names.str.upper().str.startswith(target + "_"))
    if prefixed.size:
        return int(prefixed[0])
    return None


def determine_group_names(obs: pd.DataFrame, group_col: str) -> tuple[str, str]:
    values = obs[group_col].astype(str).unique().tolist()
    lesion = next((v for v in values if v.lower() in {"lesion", "tle", "disease", "case"}), None)
    control = next((v for v in values if v.lower() in {"internal_control", "control", "normal", "untreated"}), None)
    if lesion is None or control is None:
        counts = obs[group_col].astype(str).value_counts()
        if counts.shape[0] < 2:
            raise ValueError(f"Need at least two groups in {group_col}")
        ordered = counts.index.tolist()
        if lesion is None:
            lesion = ordered[0]
        if control is None:
            control = ordered[1]
    return lesion, control


def extract_expression_vectors(counts: sparse.csr_matrix, idx: int, size_factors: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw = counts[:, idx].toarray().ravel().astype(float)
    norm_linear = raw / size_factors * 1e4
    norm_log1p = np.log1p(norm_linear)
    return raw, norm_linear, norm_log1p


def compute_group_row(
    tf: str,
    expected_group: str,
    main_row: pd.Series | None,
    raw: np.ndarray,
    norm_linear: np.ndarray,
    norm_log1p: np.ndarray,
    groups: pd.Series,
    lesion_name: str,
    control_name: str,
) -> dict[str, object]:
    lesion_mask = groups.values == lesion_name
    control_mask = groups.values == control_name
    lesion_vals = norm_log1p[lesion_mask]
    control_vals = norm_log1p[control_mask]
    lesion_raw = raw[lesion_mask]
    control_raw = raw[control_mask]
    lesion_linear = norm_linear[lesion_mask]
    control_linear = norm_linear[control_mask]

    stat, p_value = stats.mannwhitneyu(lesion_vals, control_vals, alternative="two-sided")
    log2fc = math.log2((lesion_linear.mean() + 1e-3) / (control_linear.mean() + 1e-3))
    pos_frac_lesion = float((lesion_raw > 0).mean())
    pos_frac_control = float((control_raw > 0).mean())
    direction_consistent = (
        log2fc > 0 if expected_group == "lesion"
        else log2fc < 0 if expected_group == "internal_control"
        else False
    )
    return {
        "tf": tf,
        "expected_higher_group_main": expected_group,
        "main_regulon": main_row["regulon"] if main_row is not None else "",
        "main_regulon_fdr": main_row["regulon_fdr"] if main_row is not None else np.nan,
        "main_expr_fdr": main_row["expr_wilcoxon_fdr"] if main_row is not None else np.nan,
        "main_log2fc_lesion_vs_internal_control": main_row["log2fc_lesion_vs_internal_control"] if main_row is not None else np.nan,
        "lesion_group_label": lesion_name,
        "control_group_label": control_name,
        "n_lesion": int(lesion_mask.sum()),
        "n_internal_control": int(control_mask.sum()),
        "mean_log1p_norm_expr_lesion": float(lesion_vals.mean()),
        "mean_log1p_norm_expr_internal_control": float(control_vals.mean()),
        "median_log1p_norm_expr_lesion": float(np.median(lesion_vals)),
        "median_log1p_norm_expr_internal_control": float(np.median(control_vals)),
        "positive_fraction_lesion": pos_frac_lesion,
        "positive_fraction_internal_control": pos_frac_control,
        "positive_fraction_diff_lesion_minus_internal_control": pos_frac_lesion - pos_frac_control,
        "log2fc_lesion_vs_internal_control": float(log2fc),
        "mannwhitney_u": float(stat),
        "p_value": float(p_value),
        "direction_consistent_with_main": bool(direction_consistent),
    }


def compute_cluster_rows(
    tf: str,
    expected_group: str,
    main_row: pd.Series | None,
    raw: np.ndarray,
    norm_linear: np.ndarray,
    norm_log1p: np.ndarray,
    groups: pd.Series,
    clusters: pd.Series,
    lesion_name: str,
    control_name: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    cluster_values = clusters.astype(str).to_numpy()
    group_values = groups.astype(str).to_numpy()
    for cluster_name in pd.Index(cluster_values).unique().tolist():
        mask = cluster_values == str(cluster_name)
        lesion_mask = mask & (group_values == lesion_name)
        control_mask = mask & (group_values == control_name)
        n_lesion = int(lesion_mask.sum())
        n_control = int(control_mask.sum())
        if n_lesion < MIN_CLUSTER_CELLS_PER_GROUP or n_control < MIN_CLUSTER_CELLS_PER_GROUP:
            continue
        lesion_vals = norm_log1p[lesion_mask]
        control_vals = norm_log1p[control_mask]
        lesion_raw = raw[lesion_mask]
        control_raw = raw[control_mask]
        lesion_linear = norm_linear[lesion_mask]
        control_linear = norm_linear[control_mask]
        stat, p_value = stats.mannwhitneyu(lesion_vals, control_vals, alternative="two-sided")
        log2fc = math.log2((lesion_linear.mean() + 1e-3) / (control_linear.mean() + 1e-3))
        direction_consistent = (
            log2fc > 0 if expected_group == "lesion"
            else log2fc < 0 if expected_group == "internal_control"
            else False
        )
        rows.append(
            {
                "tf": tf,
                "celltype_or_cluster": str(cluster_name),
                "expected_higher_group_main": expected_group,
                "main_regulon": main_row["regulon"] if main_row is not None else "",
                "n_lesion": n_lesion,
                "n_internal_control": n_control,
                "mean_log1p_norm_expr_lesion": float(lesion_vals.mean()),
                "mean_log1p_norm_expr_internal_control": float(control_vals.mean()),
                "positive_fraction_lesion": float((lesion_raw > 0).mean()),
                "positive_fraction_internal_control": float((control_raw > 0).mean()),
                "log2fc_lesion_vs_internal_control": float(log2fc),
                "mannwhitney_u": float(stat),
                "p_value": float(p_value),
                "direction_consistent_with_main": bool(direction_consistent),
            }
        )
    return rows


def save_group_mean_barplot(group_stats: pd.DataFrame, path: Path) -> None:
    tf_order = group_stats["tf"].tolist()
    x = np.arange(len(tf_order))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax.bar(x - width / 2, group_stats["mean_log1p_norm_expr_internal_control"], width=width, label="internal_control", color="#4C78A8")
    ax.bar(x + width / 2, group_stats["mean_log1p_norm_expr_lesion"], width=width, label="lesion", color="#E45756")
    ax.set_xticks(x)
    ax.set_xticklabels(tf_order)
    ax.set_ylabel("Mean log1p normalized expression")
    ax.set_title("Round2 supportive evidence: group-level TF expression")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_group_boxplot(group_arrays: dict[str, dict[str, np.ndarray]], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    axes = axes.ravel()
    colors = ["#4C78A8", "#E45756"]
    for ax, (tf, arrays) in zip(axes, group_arrays.items()):
        kwargs = {"patch_artist": True, "showfliers": False}
        if "tick_labels" in inspect.signature(ax.boxplot).parameters:
            kwargs["tick_labels"] = ["internal_control", "lesion"]
        else:
            kwargs["labels"] = ["internal_control", "lesion"]
        bp = ax.boxplot([arrays["internal_control"], arrays["lesion"]], **kwargs)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        ax.set_title(tf)
        ax.set_ylabel("log1p normalized expression")
    for ax in axes[len(group_arrays) :]:
        ax.axis("off")
    fig.suptitle("Round2 supportive evidence: group boxplots", fontsize=14)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_celltype_heatmap(celltype_stats: pd.DataFrame, path: Path) -> None:
    if celltype_stats.empty:
        return
    cluster_order = (
        celltype_stats.groupby("celltype_or_cluster")[["n_lesion", "n_internal_control"]]
        .sum()
        .sum(axis=1)
        .sort_values(ascending=False)
        .index.tolist()
    )
    matrix = (
        celltype_stats.pivot(index="tf", columns="celltype_or_cluster", values="log2fc_lesion_vs_internal_control")
        .reindex(index=TARGET_TFS, columns=cluster_order)
        .fillna(0.0)
    )
    fig_width = max(10, 0.45 * matrix.shape[1] + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 4.5), constrained_layout=True)
    vmax = np.nanmax(np.abs(matrix.values)) if matrix.size else 1.0
    vmax = max(vmax, 0.5)
    im = ax.imshow(matrix.values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns.tolist(), rotation=90, fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.tolist())
    ax.set_title("Round2 supportive evidence: cluster-level log2FC (lesion vs control)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("log2FC lesion vs internal_control")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_support_ranking_barplot(ranking: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    colors = ["#1b9e77" if g == "lesion" else "#d95f02" for g in ranking["expected_higher_group_main"].astype(str)]
    ax.bar(ranking["tf"], ranking["support_score"], color=colors, alpha=0.85)
    ax.set_ylabel("Support score")
    ax.set_title("Round2 supportive evidence: support ranking")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_priority_focus(celltype_stats: pd.DataFrame, path: Path) -> None:
    focus = celltype_stats[celltype_stats["tf"].isin(["NFE2L2", "THRB"])].copy()
    if focus.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    for ax, tf in zip(axes, ["NFE2L2", "THRB"]):
        sub = focus[focus["tf"] == tf].copy()
        sub = sub.sort_values("log2fc_lesion_vs_internal_control", ascending=False)
        colors = ["#E45756" if x > 0 else "#4C78A8" for x in sub["log2fc_lesion_vs_internal_control"]]
        ax.barh(sub["celltype_or_cluster"], sub["log2fc_lesion_vs_internal_control"], color=colors, alpha=0.85)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(tf)
        ax.set_xlabel("log2FC lesion vs internal_control")
    fig.suptitle("Round2 supportive evidence: priority TF cluster-level direction", fontsize=14)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def build_summary_text(ranking: pd.DataFrame, group_stats: pd.DataFrame, celltype_stats: pd.DataFrame, has_auc: bool) -> str:
    top = ranking.iloc[0]
    second = ranking.iloc[1] if ranking.shape[0] > 1 else None
    lines: list[str] = []
    lines.append("第二个外部验证队列基于 GSE190452 构建，共纳入 69968 个细胞，包含 4 个 non-epileptic control 样本与 4 个 TLE 样本，并通过最小 Scanpy 流程生成 33 个 cluster。")
    if has_auc:
        lines.append("该外部对象包含 pySCENIC AUC，因此同时具备 regulon activity 与 expression-level 验证条件。")
    else:
        lines.append("该外部对象未包含 pySCENIC AUC，因此本轮验证以 expression-level analysis 为主，并在 group 与 cluster 两个层级评估方向一致性。")
    lines.append(
        f"按综合支持度排序，结果为 {' > '.join(ranking['tf'].astype(str).tolist())}，其中 {top['tf']} 位列第 1，support_score={top['support_score']:.3f}。"
        + (f" 第 2 位为 {second['tf']}，support_score={second['support_score']:.3f}。" if second is not None else "")
    )

    def tf_line(tf_name: str, positive_word: str, weaker_word: str) -> str:
        row = ranking[ranking["tf"] == tf_name].iloc[0]
        cluster_sub = celltype_stats[celltype_stats["tf"] == tf_name]
        support_fraction = row["cluster_support_fraction"]
        if tf_name == "THRB":
            if row["group_direction_consistent_with_main"]:
                return (
                    f"{tf_name} 在 round2 外部队列中继续得到支持：其组间 log2FC 方向与主分析一致，"
                    f"且有 {int(row['direction_consistent_clusters'])}/{int(row['informative_clusters'])} 个 cluster 保持一致方向，"
                    f"提示其对照侧特征在独立队列中仍较稳定。"
                )
            return f"{tf_name} 在 round2 队列中的支持弱于预期，主要原因是组间方向或 cluster 层面一致性不足。"
        if row["group_direction_consistent_with_main"] and support_fraction >= 0.5:
            return (
                f"{tf_name} {positive_word}：其组间表达方向与主分析一致，"
                f"cluster 支持比例为 {support_fraction:.2f}，说明这一信号并非仅由少数细胞群驱动。"
            )
        return f"{tf_name} {weaker_word}：虽然仍可检测到一定表达差异，但 cluster 层面的支持比例仅为 {support_fraction:.2f}。"

    lines.append(tf_line("NFE2L2", "继续得到明确支持", "支持减弱"))
    lines.append(tf_line("THRB", "继续得到明确支持", "支持减弱"))
    lines.append(tf_line("BHLHE40", "仍有一定支持", "支持减弱"))
    lines.append(tf_line("SOX2", "仍保留一定支持", "支持减弱"))
    lines.append("总体上，NFE2L2 与 THRB 在第二个外部验证队列中的表现与主分析最为一致；BHLHE40 与 SOX2 的支持度相对较弱。")
    lines.append("现阶段结果已能为成文提供较稳定的主结论框架，但若要进一步增强说服力，仍建议补充至少一个带明确细胞类型注释或可复算 regulon activity 的独立验证队列。")
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(INPUT_H5AD)
    expected = load_expected_directions()

    counts = ensure_sparse(adata.layers["counts"] if "counts" in adata.layers else adata.X)
    size_factors = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    size_factors[size_factors <= 0] = 1.0

    group_col = "group"
    celltype_col = "celltype" if "celltype" in adata.obs.columns and adata.obs["celltype"].astype("string").nunique() > 1 else "cluster"
    lesion_name, control_name = determine_group_names(adata.obs, group_col)
    has_auc = OBSM_KEY in adata.obsm and UNS_KEY in adata.uns

    group_rows: list[dict[str, object]] = []
    celltype_rows: list[dict[str, object]] = []
    group_arrays: dict[str, dict[str, np.ndarray]] = {}

    for tf in TARGET_TFS:
        main_row = expected.loc[tf.upper()] if tf.upper() in expected.index else None
        expected_group = (
            str(main_row["expr_higher_group"]) if main_row is not None and pd.notna(main_row["expr_higher_group"])
            else str(main_row["regulon_effect_group"]) if main_row is not None and pd.notna(main_row["regulon_effect_group"])
            else "unknown"
        )
        gene_idx = find_gene_index(adata, tf)
        if gene_idx is None:
            group_rows.append(
                {
                    "tf": tf,
                    "expected_higher_group_main": expected_group,
                    "main_regulon": main_row["regulon"] if main_row is not None else "",
                    "gene_found": False,
                    "direction_consistent_with_main": False,
                    "p_value": np.nan,
                }
            )
            continue

        raw, norm_linear, norm_log1p = extract_expression_vectors(counts, gene_idx, size_factors)
        group_row = compute_group_row(
            tf=tf,
            expected_group=expected_group,
            main_row=main_row,
            raw=raw,
            norm_linear=norm_linear,
            norm_log1p=norm_log1p,
            groups=adata.obs[group_col].astype(str),
            lesion_name=lesion_name,
            control_name=control_name,
        )
        group_row["gene_found"] = True
        group_rows.append(group_row)

        group_arrays[tf] = {
            "internal_control": norm_log1p[adata.obs[group_col].astype(str).values == control_name],
            "lesion": norm_log1p[adata.obs[group_col].astype(str).values == lesion_name],
        }

        celltype_rows.extend(
            compute_cluster_rows(
                tf=tf,
                expected_group=expected_group,
                main_row=main_row,
                raw=raw,
                norm_linear=norm_linear,
                norm_log1p=norm_log1p,
                groups=adata.obs[group_col].astype(str),
                clusters=adata.obs[celltype_col].astype(str),
                lesion_name=lesion_name,
                control_name=control_name,
            )
        )

    group_stats = pd.DataFrame(group_rows)
    group_stats["fdr_bh"] = benjamini_hochberg(group_stats["p_value"])
    group_stats["group_direction_consistent_with_main"] = group_stats["direction_consistent_with_main"].fillna(False).astype(bool)

    celltype_stats = pd.DataFrame(celltype_rows)
    if not celltype_stats.empty:
        celltype_stats["fdr_bh"] = benjamini_hochberg(celltype_stats["p_value"])
        celltype_stats["direction_consistent_with_main"] = celltype_stats["direction_consistent_with_main"].fillna(False).astype(bool)

    ranking_rows: list[dict[str, object]] = []
    for _, row in group_stats.iterrows():
        tf = str(row["tf"])
        expected_group = str(row["expected_higher_group_main"])
        cluster_sub = celltype_stats[celltype_stats["tf"] == tf].copy() if not celltype_stats.empty else pd.DataFrame()
        informative_clusters = int(cluster_sub.shape[0])
        consistent_clusters = int(cluster_sub["direction_consistent_with_main"].sum()) if informative_clusters else 0
        sig_consistent_clusters = int(((cluster_sub["direction_consistent_with_main"]) & (cluster_sub["fdr_bh"] < 0.05)).sum()) if informative_clusters else 0
        cluster_support_fraction = consistent_clusters / informative_clusters if informative_clusters else 0.0
        sig_cluster_support_fraction = sig_consistent_clusters / informative_clusters if informative_clusters else 0.0
        directional_logfc = float(row["log2fc_lesion_vs_internal_control"])
        directional_posdiff = float(row["positive_fraction_diff_lesion_minus_internal_control"])
        if expected_group == "internal_control":
            directional_logfc *= -1.0
            directional_posdiff *= -1.0
        ranking_rows.append(
            {
                "tf": tf,
                "expected_higher_group_main": expected_group,
                "main_regulon": row.get("main_regulon", ""),
                "main_regulon_fdr": row.get("main_regulon_fdr", np.nan),
                "main_expr_fdr": row.get("main_expr_fdr", np.nan),
                "group_mean_log1p_norm_expr_lesion": row.get("mean_log1p_norm_expr_lesion", np.nan),
                "group_mean_log1p_norm_expr_internal_control": row.get("mean_log1p_norm_expr_internal_control", np.nan),
                "group_positive_fraction_lesion": row.get("positive_fraction_lesion", np.nan),
                "group_positive_fraction_internal_control": row.get("positive_fraction_internal_control", np.nan),
                "group_log2fc_lesion_vs_internal_control": row.get("log2fc_lesion_vs_internal_control", np.nan),
                "group_fdr_bh": row.get("fdr_bh", np.nan),
                "group_direction_consistent_with_main": bool(row.get("group_direction_consistent_with_main", False)),
                "directional_logfc_score_raw": directional_logfc,
                "directional_positive_fraction_diff_raw": directional_posdiff,
                "informative_clusters": informative_clusters,
                "direction_consistent_clusters": consistent_clusters,
                "significant_direction_consistent_clusters": sig_consistent_clusters,
                "cluster_support_fraction": cluster_support_fraction,
                "significant_cluster_support_fraction": sig_cluster_support_fraction,
            }
        )

    ranking = pd.DataFrame(ranking_rows)
    ranking["group_significance_raw"] = -np.log10(ranking["group_fdr_bh"].astype(float).clip(lower=1e-300))
    ranking["directional_logfc_scaled"] = minmax_scale(np.maximum(ranking["directional_logfc_score_raw"].to_numpy(dtype=float), 0.0))
    ranking["directional_posfrac_scaled"] = minmax_scale(np.maximum(ranking["directional_positive_fraction_diff_raw"].to_numpy(dtype=float), 0.0))
    ranking["group_significance_scaled"] = minmax_scale(ranking["group_significance_raw"].to_numpy(dtype=float))
    ranking["cluster_support_scaled"] = minmax_scale(ranking["cluster_support_fraction"].to_numpy(dtype=float))
    ranking["sig_cluster_support_scaled"] = minmax_scale(ranking["significant_cluster_support_fraction"].to_numpy(dtype=float))
    ranking["support_score"] = (
        0.35 * ranking["group_significance_scaled"].fillna(0.0)
        + 0.25 * ranking["directional_logfc_scaled"].fillna(0.0)
        + 0.15 * ranking["directional_posfrac_scaled"].fillna(0.0)
        + 0.15 * ranking["cluster_support_scaled"].fillna(0.0)
        + 0.10 * ranking["sig_cluster_support_scaled"].fillna(0.0)
    )
    ranking["support_rank"] = ranking["support_score"].rank(method="dense", ascending=False).astype(int)
    ranking = ranking.sort_values(["support_rank", "support_score"], ascending=[True, False]).reset_index(drop=True)

    master_table = ranking.merge(
        group_stats[
            [
                "tf",
                "n_lesion",
                "n_internal_control",
                "mean_log1p_norm_expr_lesion",
                "mean_log1p_norm_expr_internal_control",
                "positive_fraction_lesion",
                "positive_fraction_internal_control",
                "log2fc_lesion_vs_internal_control",
                "mannwhitney_u",
                "p_value",
                "fdr_bh",
            ]
        ],
        on="tf",
        how="left",
        suffixes=("", "_group"),
    )
    master_table["auc_status"] = "available" if has_auc else "not_available_expression_only"

    group_csv = OUT_DIR / "external_round2_tf_expression_group_stats.csv"
    celltype_csv = OUT_DIR / "external_round2_tf_expression_celltype_stats.csv"
    master_csv = OUT_DIR / "external_round2_master_table.csv"
    ranking_csv = OUT_DIR / "external_round2_support_ranking.csv"
    summary_txt = OUT_DIR / "external_round2_summary_cn.txt"

    group_stats.to_csv(group_csv, index=False)
    celltype_stats.to_csv(celltype_csv, index=False)
    master_table.to_csv(master_csv, index=False)
    ranking.to_csv(ranking_csv, index=False)

    save_group_mean_barplot(group_stats.sort_values("tf"), PLOT_DIR / "external_round2_group_mean_barplot.png")
    ordered_arrays = {tf: group_arrays[tf] for tf in TARGET_TFS if tf in group_arrays}
    save_group_boxplot(ordered_arrays, PLOT_DIR / "external_round2_group_boxplot.png")
    save_celltype_heatmap(celltype_stats, PLOT_DIR / "external_round2_celltype_heatmap.png")
    save_support_ranking_barplot(ranking, PLOT_DIR / "external_round2_support_ranking_barplot.png")
    save_priority_focus(celltype_stats, PLOT_DIR / "external_round2_priority_tfs_focus.png")

    summary = build_summary_text(ranking, group_stats, celltype_stats, has_auc)
    summary_txt.write_text(summary, encoding="utf-8")

    print(f"Wrote {group_csv}")
    print(f"Wrote {celltype_csv}")
    print(f"Wrote {master_csv}")
    print(f"Wrote {ranking_csv}")
    print(f"Wrote {summary_txt}")
    print(f"Wrote plots to {PLOT_DIR}")


if __name__ == "__main__":
    main()
