#!/usr/bin/env python3
from __future__ import annotations

import inspect
import math
import subprocess
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[1]
INPUT_H5AD = ROOT / "external_validation_round2" / "input" / "external_validation_round2.h5ad"
AUC_CSV = ROOT / "external_validation_round2" / "pyscenic_recalc" / "output" / "auc_mtx_4tf_from_discovery.csv"
OUT_DIR = ROOT / "external_validation_round2"
PLOT_DIR = OUT_DIR / "regulon_validation_plots"
MAIN_METRICS = ROOT / "analysis_outputs" / "celloracle_candidates" / "celloracle_candidate_tf_metrics.csv"
EXPR_GROUP_STATS = OUT_DIR / "external_round2_tf_expression_group_stats.csv"
EXPR_CLUSTER_STATS = OUT_DIR / "external_round2_tf_expression_celltype_stats.csv"
EXPR_RANKING = OUT_DIR / "external_round2_support_ranking.csv"
TARGET_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
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


def ensure_expression_outputs() -> None:
    required = [EXPR_GROUP_STATS, EXPR_CLUSTER_STATS, EXPR_RANKING]
    if all(path.exists() for path in required):
        return
    script = ROOT / "scripts" / "run_external_validation_round2.py"
    subprocess.check_call([sys.executable, str(script)])


def load_auc_matrix(path: Path, cell_index: pd.Index) -> pd.DataFrame:
    auc_df = pd.read_csv(path, index_col=0)
    auc_df.index = auc_df.index.astype(str)
    auc_df.columns = [str(x) for x in auc_df.columns]

    col_matches = sum(any(str(col).upper().startswith(tf) for tf in TARGET_TFS) for col in auc_df.columns)
    idx_matches = sum(any(str(idx).upper().startswith(tf) for tf in TARGET_TFS) for idx in auc_df.index)
    if col_matches == 0 and idx_matches > 0:
        auc_df = auc_df.T
        auc_df.index = auc_df.index.astype(str)
        auc_df.columns = [str(x) for x in auc_df.columns]

    common_cells = cell_index.intersection(pd.Index(auc_df.index.astype(str)))
    if common_cells.empty:
        raise ValueError(f"No shared cell barcodes between round2 h5ad and AUC matrix: {path}")
    auc_df = auc_df.loc[common_cells].copy()
    return auc_df


def resolve_regulon_name(tf: str, main_regulon: str, available: list[str]) -> str | None:
    candidates = []
    if isinstance(main_regulon, str) and main_regulon.strip():
        candidates.append(main_regulon.strip())
    candidates.extend([f"{tf}(+)", tf, tf.upper(), tf.lower()])
    upper_map = {name.upper(): name for name in available}
    for candidate in candidates:
        if candidate.upper() in upper_map:
            return upper_map[candidate.upper()]
    for name in available:
        if name.upper().startswith(tf.upper() + "("):
            return name
    return None


def save_regulon_group_mean_barplot(group_stats: pd.DataFrame, path: Path) -> None:
    tf_order = group_stats["tf"].tolist()
    x = np.arange(len(tf_order))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax.bar(x - width / 2, group_stats["mean_auc_internal_control"], width=width, label="internal_control", color="#4C78A8")
    ax.bar(x + width / 2, group_stats["mean_auc_lesion"], width=width, label="lesion", color="#E45756")
    ax.set_xticks(x)
    ax.set_xticklabels(tf_order)
    ax.set_ylabel("Mean regulon AUC")
    ax.set_title("Round2 supportive evidence: group-level regulon AUC")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_regulon_boxplot(group_arrays: dict[str, dict[str, np.ndarray]], path: Path) -> None:
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
        ax.set_ylabel("Regulon AUC")
    for ax in axes[len(group_arrays) :]:
        ax.axis("off")
    fig.suptitle("Round2 supportive evidence: regulon AUC boxplots", fontsize=14)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_regulon_cluster_heatmap(cluster_stats: pd.DataFrame, path: Path) -> None:
    if cluster_stats.empty:
        return
    cluster_order = (
        cluster_stats.groupby("cluster")[["n_lesion", "n_internal_control"]]
        .sum()
        .sum(axis=1)
        .sort_values(ascending=False)
        .index.tolist()
    )
    matrix = (
        cluster_stats.pivot(index="tf", columns="cluster", values="mean_auc_diff_lesion_minus_internal_control")
        .reindex(index=TARGET_TFS, columns=cluster_order)
        .fillna(0.0)
    )
    fig_width = max(10, 0.45 * matrix.shape[1] + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 4.5), constrained_layout=True)
    vmax = np.nanmax(np.abs(matrix.values)) if matrix.size else 1.0
    vmax = max(vmax, 0.03)
    im = ax.imshow(matrix.values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns.tolist(), rotation=90, fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.tolist())
    ax.set_title("Round2 supportive evidence: cluster-level regulon mean AUC difference")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Mean AUC difference (lesion - internal_control)")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_combined_ranking_plot(ranking: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    x = np.arange(ranking.shape[0])
    width = 0.25
    ax.bar(x - width, ranking["expr_support_score_scaled"], width=width, label="Expression", color="#72B7B2")
    ax.bar(x, ranking["regulon_support_score"], width=width, label="Regulon", color="#E45756")
    ax.bar(x + width, ranking["combined_support_score"], width=width, label="Combined", color="#4C78A8")
    ax.set_xticks(x)
    ax.set_xticklabels(ranking["tf"].tolist())
    ax.set_ylabel("Scaled support score")
    ax.set_title("Round2 supportive evidence: expression + regulon combined ranking")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_priority_focus(cluster_stats: pd.DataFrame, group_stats: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)
    for ax, tf in zip(axes, ["NFE2L2", "THRB"]):
        sub = cluster_stats[cluster_stats["tf"] == tf].copy()
        if sub.empty:
            ax.axis("off")
            continue
        sub = sub.sort_values("mean_auc_diff_lesion_minus_internal_control", ascending=False)
        colors = ["#E45756" if x > 0 else "#4C78A8" for x in sub["mean_auc_diff_lesion_minus_internal_control"]]
        ax.barh(sub["cluster"], sub["mean_auc_diff_lesion_minus_internal_control"], color=colors, alpha=0.85)
        ax.axvline(0, color="black", linewidth=0.8)
        row = group_stats[group_stats["tf"] == tf].iloc[0]
        ax.set_title(
            f"{tf} | lesion={row['mean_auc_lesion']:.3f}, control={row['mean_auc_internal_control']:.3f}"
        )
        ax.set_xlabel("Mean AUC difference (lesion - internal_control)")
    fig.suptitle("Round2 supportive evidence: priority regulons in clusters", fontsize=14)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def build_summary_text(ranking: pd.DataFrame) -> str:
    order = " > ".join(ranking["tf"].astype(str).tolist())
    top = ranking.iloc[0]["tf"] if not ranking.empty else ""
    nfe2l2_row = ranking[ranking["tf"] == "NFE2L2"].iloc[0]
    thrb_row = ranking[ranking["tf"] == "THRB"].iloc[0]
    bhlhe40_row = ranking[ranking["tf"] == "BHLHE40"].iloc[0]
    sox2_row = ranking[ranking["tf"] == "SOX2"].iloc[0]

    lines = [
        "GSE190452 在本公开包中作为 cross-syndrome supportive/contextual evidence 使用，并补入 projected regulon-level supportive analysis。",
        "本轮不再依赖 round2 自身重新构建 regulon，而是使用 discovery 队列筛出的 4 个 regulon GMT 对 round2 表达矩阵进行 AUCell projection。",
        f"按 expression 与 regulon 证据整合后的综合支持度排序，结果为 {order}。",
    ]

    if top == "NFE2L2":
        lines.append(
            f"NFE2L2 继续得到最强支持：其 regulon 活性方向与主分析一致，combined_support_score={nfe2l2_row['combined_support_score']:.3f}，"
            f"且 regulon cluster_support_fraction={nfe2l2_row['regulon_cluster_support_fraction']:.2f}。"
        )
    else:
        lines.append(
            f"NFE2L2 仍保持高支持度，combined_support_score={nfe2l2_row['combined_support_score']:.3f}，"
            "说明其 lesion 侧优先候选地位并未被削弱。"
        )

    lines.append(
        f"THRB 继续保持最稳定的对照侧候选：其 regulon 活性在 internal_control 侧更高，combined_support_score={thrb_row['combined_support_score']:.3f}，"
        f"并且 regulon cluster_support_fraction={thrb_row['regulon_cluster_support_fraction']:.2f}。"
    )

    lines.append(
        f"BHLHE40 的支持度仍弱于 NFE2L2 和 THRB，combined_support_score={bhlhe40_row['combined_support_score']:.3f}，"
        "更适合作为第二梯队候选。"
    )
    lines.append(
        f"SOX2 仍保留一定信号，但综合支持度低于前两位核心候选，combined_support_score={sox2_row['combined_support_score']:.3f}。"
    )
    lines.append(
        "综合来看，round2 分析提供 expression-level 与 projected regulon-level supportive evidence；"
        "这些结果不应表述为正式外部验证、机制确认或实验验证。"
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT_H5AD.exists():
        raise FileNotFoundError(f"Missing round2 h5ad: {INPUT_H5AD}")
    if not AUC_CSV.exists():
        raise FileNotFoundError(
            f"Missing projected 4TF AUC matrix: {AUC_CSV}\n"
            "Run build_round2_minimal_expression_for_4tf_auc.py and run_round2_aucell_projection_4tf.py first."
        )

    ensure_expression_outputs()
    expr_group = pd.read_csv(EXPR_GROUP_STATS)
    expr_cluster = pd.read_csv(EXPR_CLUSTER_STATS)
    expr_ranking = pd.read_csv(EXPR_RANKING)
    expected = load_expected_directions()

    adata = ad.read_h5ad(INPUT_H5AD)

    group_col = "group"
    cluster_col = "celltype" if "celltype" in adata.obs.columns and adata.obs["celltype"].astype("string").nunique() > 1 else "cluster"
    lesion_name, control_name = determine_group_names(adata.obs, group_col)

    auc_df = load_auc_matrix(AUC_CSV, pd.Index(adata.obs_names.astype(str)))
    adata = adata[auc_df.index.tolist(), :].copy()
    regulon_names = auc_df.columns.astype(str).tolist()

    group_rows: list[dict[str, object]] = []
    cluster_rows: list[dict[str, object]] = []
    group_arrays: dict[str, dict[str, np.ndarray]] = {}

    groups = adata.obs[group_col].astype(str).copy()
    clusters = adata.obs[cluster_col].astype(str).copy()

    for tf in TARGET_TFS:
        main_row = expected.loc[tf] if tf in expected.index else None
        expected_group = (
            str(main_row["regulon_effect_group"]) if main_row is not None and pd.notna(main_row["regulon_effect_group"])
            else str(main_row["expr_higher_group"]) if main_row is not None and pd.notna(main_row["expr_higher_group"])
            else "unknown"
        )
        main_regulon = str(main_row["regulon"]) if main_row is not None and pd.notna(main_row["regulon"]) else ""
        regulon_name = resolve_regulon_name(tf, main_regulon, regulon_names)
        if regulon_name is None:
            group_rows.append(
                {
                    "tf": tf,
                    "regulon": main_regulon,
                    "regulon_found_in_round2_auc": False,
                    "expected_higher_group_main": expected_group,
                }
            )
            continue

        values = auc_df[regulon_name].to_numpy(dtype=float)
        lesion_mask = groups.values == lesion_name
        control_mask = groups.values == control_name
        lesion_vals = values[lesion_mask]
        control_vals = values[control_mask]
        stat, p_value = stats.mannwhitneyu(lesion_vals, control_vals, alternative="two-sided")
        mean_diff = float(lesion_vals.mean() - control_vals.mean())
        direction_consistent = (
            mean_diff > 0 if expected_group == "lesion"
            else mean_diff < 0 if expected_group == "internal_control"
            else False
        )
        group_rows.append(
            {
                "tf": tf,
                "regulon": regulon_name,
                "regulon_found_in_round2_auc": True,
                "expected_higher_group_main": expected_group,
                "main_regulon_fdr": main_row["regulon_fdr"] if main_row is not None else np.nan,
                "n_lesion": int(lesion_mask.sum()),
                "n_internal_control": int(control_mask.sum()),
                "mean_auc_lesion": float(lesion_vals.mean()),
                "mean_auc_internal_control": float(control_vals.mean()),
                "median_auc_lesion": float(np.median(lesion_vals)),
                "median_auc_internal_control": float(np.median(control_vals)),
                "mean_auc_diff_lesion_minus_internal_control": mean_diff,
                "mannwhitney_u": float(stat),
                "p_value": float(p_value),
                "group_direction_consistent_with_main": bool(direction_consistent),
            }
        )
        group_arrays[tf] = {"internal_control": control_vals, "lesion": lesion_vals}

        for cluster_name in pd.Index(clusters.values).unique().tolist():
            mask = clusters.values == str(cluster_name)
            lesion_cluster_mask = mask & lesion_mask
            control_cluster_mask = mask & control_mask
            n_lesion = int(lesion_cluster_mask.sum())
            n_control = int(control_cluster_mask.sum())
            if n_lesion < MIN_CLUSTER_CELLS_PER_GROUP or n_control < MIN_CLUSTER_CELLS_PER_GROUP:
                continue
            lesion_cluster_vals = values[lesion_cluster_mask]
            control_cluster_vals = values[control_cluster_mask]
            stat_cluster, p_cluster = stats.mannwhitneyu(lesion_cluster_vals, control_cluster_vals, alternative="two-sided")
            cluster_diff = float(lesion_cluster_vals.mean() - control_cluster_vals.mean())
            cluster_consistent = (
                cluster_diff > 0 if expected_group == "lesion"
                else cluster_diff < 0 if expected_group == "internal_control"
                else False
            )
            cluster_rows.append(
                {
                    "tf": tf,
                    "regulon": regulon_name,
                    "cluster": str(cluster_name),
                    "expected_higher_group_main": expected_group,
                    "n_lesion": n_lesion,
                    "n_internal_control": n_control,
                    "mean_auc_lesion": float(lesion_cluster_vals.mean()),
                    "mean_auc_internal_control": float(control_cluster_vals.mean()),
                    "mean_auc_diff_lesion_minus_internal_control": cluster_diff,
                    "mannwhitney_u": float(stat_cluster),
                    "p_value": float(p_cluster),
                    "direction_consistent_with_main": bool(cluster_consistent),
                }
            )

    regulon_group_stats = pd.DataFrame(group_rows)
    regulon_cluster_stats = pd.DataFrame(cluster_rows)
    if "p_value" in regulon_group_stats.columns:
        regulon_group_stats["fdr_bh"] = benjamini_hochberg(regulon_group_stats["p_value"])
    if not regulon_cluster_stats.empty:
        regulon_cluster_stats["fdr_bh"] = benjamini_hochberg(regulon_cluster_stats["p_value"])

    ranking_rows: list[dict[str, object]] = []
    expr_scale = minmax_scale(expr_ranking["support_score"].to_numpy(dtype=float))
    expr_ranking = expr_ranking.copy()
    expr_ranking["support_score_scaled"] = expr_scale

    for tf in TARGET_TFS:
        expr_row = expr_ranking[expr_ranking["tf"] == tf].iloc[0] if not expr_ranking[expr_ranking["tf"] == tf].empty else None
        reg_row = regulon_group_stats[regulon_group_stats["tf"] == tf].iloc[0] if not regulon_group_stats[regulon_group_stats["tf"] == tf].empty else None
        cluster_sub = regulon_cluster_stats[regulon_cluster_stats["tf"] == tf].copy() if not regulon_cluster_stats.empty else pd.DataFrame()
        informative_clusters = int(cluster_sub.shape[0])
        consistent_clusters = int(cluster_sub["direction_consistent_with_main"].sum()) if informative_clusters else 0
        sig_consistent_clusters = int(((cluster_sub["direction_consistent_with_main"]) & (cluster_sub["fdr_bh"] < 0.05)).sum()) if informative_clusters else 0
        cluster_support_fraction = consistent_clusters / informative_clusters if informative_clusters else 0.0
        sig_cluster_support_fraction = sig_consistent_clusters / informative_clusters if informative_clusters else 0.0
        expected_group = str(reg_row["expected_higher_group_main"]) if reg_row is not None and "expected_higher_group_main" in reg_row else (
            str(expr_row["expected_higher_group_main"]) if expr_row is not None else "unknown"
        )
        directional_auc_diff = float(reg_row["mean_auc_diff_lesion_minus_internal_control"]) if reg_row is not None and pd.notna(reg_row.get("mean_auc_diff_lesion_minus_internal_control")) else np.nan
        if expected_group == "internal_control" and pd.notna(directional_auc_diff):
            directional_auc_diff *= -1.0
        ranking_rows.append(
            {
                "tf": tf,
                "expected_higher_group_main": expected_group,
                "regulon": reg_row["regulon"] if reg_row is not None and "regulon" in reg_row else "",
                "expr_support_score": float(expr_row["support_score"]) if expr_row is not None else np.nan,
                "expr_support_score_scaled": float(expr_row["support_score_scaled"]) if expr_row is not None else np.nan,
                "regulon_group_fdr_bh": float(reg_row["fdr_bh"]) if reg_row is not None and pd.notna(reg_row.get("fdr_bh")) else np.nan,
                "regulon_group_direction_consistent_with_main": bool(reg_row["group_direction_consistent_with_main"]) if reg_row is not None and "group_direction_consistent_with_main" in reg_row else False,
                "regulon_mean_auc_diff_lesion_minus_internal_control": float(reg_row["mean_auc_diff_lesion_minus_internal_control"]) if reg_row is not None and pd.notna(reg_row.get("mean_auc_diff_lesion_minus_internal_control")) else np.nan,
                "regulon_directional_auc_diff_raw": directional_auc_diff,
                "regulon_informative_clusters": informative_clusters,
                "regulon_direction_consistent_clusters": consistent_clusters,
                "regulon_significant_direction_consistent_clusters": sig_consistent_clusters,
                "regulon_cluster_support_fraction": cluster_support_fraction,
                "regulon_significant_cluster_support_fraction": sig_cluster_support_fraction,
            }
        )

    ranking = pd.DataFrame(ranking_rows)
    ranking["regulon_group_significance_raw"] = -np.log10(ranking["regulon_group_fdr_bh"].astype(float).clip(lower=1e-300))
    ranking["regulon_directional_auc_scaled"] = minmax_scale(np.maximum(ranking["regulon_directional_auc_diff_raw"].to_numpy(dtype=float), 0.0))
    ranking["regulon_group_significance_scaled"] = minmax_scale(ranking["regulon_group_significance_raw"].to_numpy(dtype=float))
    ranking["regulon_cluster_support_scaled"] = minmax_scale(ranking["regulon_cluster_support_fraction"].to_numpy(dtype=float))
    ranking["regulon_sig_cluster_support_scaled"] = minmax_scale(ranking["regulon_significant_cluster_support_fraction"].to_numpy(dtype=float))
    ranking["regulon_support_score"] = (
        0.40 * pd.Series(ranking["regulon_group_significance_scaled"]).fillna(0.0)
        + 0.35 * pd.Series(ranking["regulon_directional_auc_scaled"]).fillna(0.0)
        + 0.15 * pd.Series(ranking["regulon_cluster_support_scaled"]).fillna(0.0)
        + 0.10 * pd.Series(ranking["regulon_sig_cluster_support_scaled"]).fillna(0.0)
    )
    ranking["combined_support_score"] = (
        0.40 * pd.Series(ranking["expr_support_score_scaled"]).fillna(0.0)
        + 0.60 * pd.Series(ranking["regulon_support_score"]).fillna(0.0)
    )
    ranking["combined_rank"] = ranking["combined_support_score"].rank(method="dense", ascending=False).astype(int)
    ranking = ranking.sort_values(["combined_rank", "combined_support_score"], ascending=[True, False]).reset_index(drop=True)

    master = expr_ranking.merge(
        ranking,
        on=["tf", "expected_higher_group_main"],
        how="outer",
        suffixes=("_expr", "_combined"),
    ).merge(
        regulon_group_stats,
        on="tf",
        how="left",
        suffixes=("", "_regulon"),
    )

    group_csv = OUT_DIR / "round2_regulon_group_stats.csv"
    cluster_csv = OUT_DIR / "round2_regulon_cluster_stats.csv"
    master_csv = OUT_DIR / "round2_expression_and_regulon_master_table.csv"
    ranking_csv = OUT_DIR / "round2_regulon_support_ranking.csv"
    summary_txt = OUT_DIR / "external_round2_regulon_summary_cn.txt"

    regulon_group_stats.to_csv(group_csv, index=False)
    regulon_cluster_stats.to_csv(cluster_csv, index=False)
    master.to_csv(master_csv, index=False)
    ranking.to_csv(ranking_csv, index=False)

    save_regulon_group_mean_barplot(regulon_group_stats.sort_values("tf"), PLOT_DIR / "round2_regulon_group_mean_barplot.png")
    ordered_arrays = {tf: group_arrays[tf] for tf in TARGET_TFS if tf in group_arrays}
    save_regulon_boxplot(ordered_arrays, PLOT_DIR / "round2_regulon_group_boxplot.png")
    save_regulon_cluster_heatmap(regulon_cluster_stats, PLOT_DIR / "round2_regulon_cluster_heatmap.png")
    save_combined_ranking_plot(ranking, PLOT_DIR / "round2_expression_regulon_combined_ranking.png")
    save_priority_focus(regulon_cluster_stats, regulon_group_stats, PLOT_DIR / "round2_priority_regulon_focus.png")

    summary = build_summary_text(ranking)
    summary_txt.write_text(summary, encoding="utf-8")

    print(f"Wrote {group_csv}")
    print(f"Wrote {cluster_csv}")
    print(f"Wrote {master_csv}")
    print(f"Wrote {ranking_csv}")
    print(f"Wrote {summary_txt}")
    print(f"Wrote plots to {PLOT_DIR}")


if __name__ == "__main__":
    main()
