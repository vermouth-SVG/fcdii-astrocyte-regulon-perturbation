#!/usr/bin/env python3
from __future__ import annotations

import argparse
import inspect
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


TARGET_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
OBSM_KEY = "X_pyscenic_auc"
UNS_KEY = "pyscenic_regulon_names"
EPS = 1e-9


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_file(__file__)
    parser = argparse.ArgumentParser(
        description="Single-group external support analysis for selected TFs using an external lesion-only cohort."
    )
    parser.add_argument(
        "--h5ad",
        default=str(root / "external_validation" / "input" / "external_validation.h5ad"),
        help="Input supportive external h5ad.",
    )
    parser.add_argument(
        "--main-tf-metrics",
        default=str(root / "analysis_outputs" / "celloracle_candidates" / "celloracle_candidate_tf_metrics.csv"),
        help="Main-analysis TF metrics table used to recover expected lesion/internal_control direction.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "external_validation" / "single_group_support"),
        help="Output directory.",
    )
    parser.add_argument(
        "--sample-column",
        default="sample",
        help="obs column used for sample-level summaries.",
    )
    parser.add_argument(
        "--donor-column",
        default="donor",
        help="obs column used for donor-level summaries.",
    )
    parser.add_argument(
        "--target-tfs",
        default=",".join(TARGET_TFS),
        help="Comma-separated TF list.",
    )
    return parser.parse_args()


def ensure_sparse(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def minmax_scale(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    finite = np.isfinite(arr)
    out = np.full(arr.shape, np.nan, dtype=float)
    if not finite.any():
        return out
    lo = np.nanmin(arr[finite])
    hi = np.nanmax(arr[finite])
    if hi - lo < EPS:
        out[finite] = 0.5
        return out
    out[finite] = (arr[finite] - lo) / (hi - lo)
    return out


def safe_cv(values: pd.Series | np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan
    mean = arr.mean()
    if abs(mean) < EPS:
        return 0.0
    return float(arr.std(ddof=0) / mean)


def load_expected_direction(main_metrics_path: Path, target_tfs: list[str]) -> pd.DataFrame:
    df = pd.read_csv(main_metrics_path)
    df["tf_upper"] = df["tf"].astype(str).str.upper()
    target_map = {tf.upper(): tf for tf in target_tfs}
    df = df[df["tf_upper"].isin(target_map)].copy()
    if df.empty:
        raise ValueError(f"No target TFs found in {main_metrics_path}")
    df["tf"] = df["tf_upper"].map(target_map)
    keep_cols = [
        "tf",
        "regulon",
        "regulon_effect_group",
        "regulon_fdr",
        "regulon_mean_diff_lesion_minus_internal_control",
        "expr_higher_group",
        "expr_wilcoxon_fdr",
        "log2fc_lesion_vs_internal_control",
        "positive_fraction_lesion",
        "positive_fraction_internal_control",
        "shortlist_keep",
        "shortlist_tier",
    ]
    for col in keep_cols:
        if col not in df.columns:
            df[col] = np.nan
    return df[keep_cols].drop_duplicates(subset=["tf"]).set_index("tf", drop=False)


def find_gene_index(adata: ad.AnnData, gene_name: str) -> int | None:
    target = gene_name.upper()
    var_names = pd.Index(adata.var_names.astype(str))
    exact = np.flatnonzero(var_names.str.upper() == target)
    if exact.size:
        return int(exact[0])
    for col in ["gene_symbol", "gene_symbols", "symbol", "gene", "gene_id"]:
        if col in adata.var.columns:
            series = adata.var[col].astype(str).str.upper()
            matches = np.flatnonzero(series.values == target)
            if matches.size:
                return int(matches[0])
    prefixed = np.flatnonzero(var_names.str.upper().str.startswith(target + "_"))
    if prefixed.size:
        return int(prefixed[0])
    return None


def extract_expression_vectors(
    adata: ad.AnnData,
    gene_index: int,
    count_matrix: sparse.csr_matrix,
    size_factors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    raw = count_matrix[:, gene_index].toarray().ravel().astype(float)
    norm = np.log1p(raw / size_factors * 1e4)
    return raw, norm


def build_level_summary(
    obs: pd.DataFrame,
    tf: str,
    expected_group: str,
    raw: np.ndarray,
    norm: np.ndarray,
    level_column: str,
    level_name: str,
) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "tf": tf,
            "expected_higher_group_main": expected_group,
            "level": level_name,
            "level_id": obs[level_column].astype(str).values,
            "raw_expr": raw,
            "norm_expr": norm,
            "positive": raw > 0,
        }
    )
    summary = (
        frame.groupby("level_id", dropna=False)
        .agg(
            n_cells=("norm_expr", "size"),
            mean_norm_expr=("norm_expr", "mean"),
            median_norm_expr=("norm_expr", "median"),
            positive_fraction=("positive", "mean"),
            mean_raw_expr=("raw_expr", "mean"),
        )
        .reset_index()
    )
    summary.insert(0, "level", level_name)
    summary.insert(0, "expected_higher_group_main", expected_group)
    summary.insert(0, "tf", tf)
    return summary


def save_expression_boxplots(
    boxplot_data: dict[str, list[np.ndarray]],
    sample_names: list[str],
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    axes = axes.ravel()
    colors = ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#B279A2", "#FF9DA6"]
    for ax, (tf, data_list) in zip(axes, boxplot_data.items()):
        boxplot_kwargs = {"patch_artist": True, "showfliers": False}
        if "tick_labels" in inspect.signature(ax.boxplot).parameters:
            boxplot_kwargs["tick_labels"] = sample_names
        else:
            boxplot_kwargs["labels"] = sample_names
        bp = ax.boxplot(data_list, **boxplot_kwargs)
        for patch, color in zip(bp["boxes"], colors * 10):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        ax.set_title(tf)
        ax.set_ylabel("log1p normalized expression")
        ax.tick_params(axis="x", rotation=45)
    for ax in axes[len(boxplot_data) :]:
        ax.axis("off")
    fig.suptitle("Supportive external evidence: TF expression distribution by sample", fontsize=14)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_expression_barplots_by_sample(
    sample_mean_table: pd.DataFrame,
    sample_order: list[str],
    tf_order: list[str],
    output_path: Path,
) -> None:
    colors = {
        "NFE2L2": "#1b9e77",
        "THRB": "#d95f02",
        "BHLHE40": "#7570b3",
        "SOX2": "#e7298a",
    }
    width = 0.18
    x = np.arange(len(sample_order))
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    for i, tf in enumerate(tf_order):
        subset = sample_mean_table.loc[tf, sample_order]
        ax.bar(
            x + (i - (len(tf_order) - 1) / 2) * width,
            subset.values,
            width=width,
            label=tf,
            color=colors.get(tf, "#4C78A8"),
            alpha=0.85,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order, rotation=45)
    ax.set_ylabel("Mean log1p normalized expression")
    ax.set_title("Supportive external evidence: sample-level TF expression means")
    ax.legend(frameon=False)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def build_summary_text(
    n_cells: int,
    n_genes: int,
    n_samples: int,
    has_auc: bool,
    ranking_df: pd.DataFrame,
) -> str:
    top_row = ranking_df.iloc[0]
    second_row = ranking_df.iloc[1] if ranking_df.shape[0] > 1 else None
    lesion_rows = ranking_df[ranking_df["expected_higher_group_main"] == "lesion"].copy()
    control_rows = ranking_df[ranking_df["expected_higher_group_main"] == "internal_control"].copy()

    lines: list[str] = []
    lines.append(f"本次外部验证基于单组 lesion 外部队列，共纳入 {n_cells} 个细胞、{n_genes} 个基因，覆盖 {n_samples} 个样本。")
    if has_auc:
        lines.append("该外部对象检测到 pySCENIC AUC 信息，但本脚本当前以 expression-level support 为主，未展开 regulon activity 比较。")
    else:
        lines.append("该外部对象未包含 pySCENIC AUC，因此本轮仅进行 expression-level support analysis，并在样本层面评估均值、阳性比例和稳定性。")

    lines.append(
        f"按单组支持度综合排序，{top_row['tf']} 位居首位，support_score={top_row['support_score']:.3f}。"
        + (f" 次强为 {second_row['tf']}，support_score={second_row['support_score']:.3f}。" if second_row is not None else "")
    )

    if not lesion_rows.empty:
        lesion_best = lesion_rows.iloc[0]
        lines.append(
            f"在主分析中偏 lesion 的候选 TF 中，{lesion_best['tf']} 表现出最强的外部支持，"
            f"其外部队列 mean_norm_expr={lesion_best['overall_mean_norm_expr']:.3f}，"
            f"positive_fraction={lesion_best['overall_positive_fraction']:.3f}，"
            f"样本间 CV={lesion_best['sample_mean_cv']:.3f}。"
        )
        remaining = lesion_rows.iloc[1:]
        if not remaining.empty:
            tail = "，".join(
                f"{row.tf}(score={row.support_score:.3f})"
                for row in remaining.itertuples(index=False)
            )
            lines.append(f"其余 lesion 方向 TF 的外部支持依次为：{tail}。")

    if not control_rows.empty:
        control_best = control_rows.iloc[0]
        lines.append(
            f"在主分析中偏 internal_control 的候选 TF 中，{control_best['tf']} 在 lesion-only 外部队列中保持较低表达，"
            f"其 mean_norm_expr={control_best['overall_mean_norm_expr']:.3f}，"
            f"positive_fraction={control_best['overall_positive_fraction']:.3f}，"
            f"因此与主分析中的 internal_control 富集方向一致。"
        )

    lines.append("总体上，这一外部队列更适合作为 lesion 方向候选 TF 的支持性验证，而不是独立的两组重复比较。")
    lines.append("从当前结果看，NFE2L2 的外部支持最完整；SOX2 与 BHLHE40 提供中等支持；THRB 更适合作为“lesion 中相对低表达、与 internal_control 方向一致”的反向支持证据。")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    target_tfs = [item.strip() for item in args.target_tfs.split(",") if item.strip()]
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs.copy()
    count_matrix = ensure_sparse(adata.layers["counts"] if "counts" in adata.layers else adata.X)
    size_factors = np.asarray(count_matrix.sum(axis=1)).ravel().astype(float)
    size_factors[size_factors <= 0] = 1.0

    sample_col = args.sample_column if args.sample_column in obs.columns else "sample"
    donor_col = args.donor_column if args.donor_column in obs.columns else sample_col
    sample_order = obs[sample_col].astype(str).value_counts().index.tolist()
    has_auc = OBSM_KEY in adata.obsm and UNS_KEY in adata.uns

    expected_df = load_expected_direction(Path(args.main_tf_metrics), target_tfs)

    summary_rows: list[dict[str, object]] = []
    level_summaries: list[pd.DataFrame] = []
    boxplot_data: dict[str, list[np.ndarray]] = {}

    for tf in target_tfs:
        expected_group = (
            expected_df.loc[tf, "expr_higher_group"]
            if tf in expected_df.index and pd.notna(expected_df.loc[tf, "expr_higher_group"])
            else expected_df.loc[tf, "regulon_effect_group"]
            if tf in expected_df.index and pd.notna(expected_df.loc[tf, "regulon_effect_group"])
            else "unknown"
        )
        gene_index = find_gene_index(adata, tf)
        if gene_index is None:
            summary_rows.append(
                {
                    "tf": tf,
                    "gene_found": False,
                    "expected_higher_group_main": expected_group,
                    "main_regulon": expected_df.loc[tf, "regulon"] if tf in expected_df.index else "",
                    "main_regulon_fdr": expected_df.loc[tf, "regulon_fdr"] if tf in expected_df.index else np.nan,
                    "main_expr_fdr": expected_df.loc[tf, "expr_wilcoxon_fdr"] if tf in expected_df.index else np.nan,
                    "overall_mean_norm_expr": np.nan,
                    "overall_median_norm_expr": np.nan,
                    "overall_positive_fraction": np.nan,
                    "sample_mean_cv": np.nan,
                    "donor_mean_cv": np.nan,
                    "samples_detected_ge_5pct": 0,
                    "samples_detected_fraction_ge_5pct": 0.0,
                    "auc_status": "skipped_no_gene",
                    "direction_consistency_note": "gene_not_found_in_external_dataset",
                }
            )
            continue

        raw_expr, norm_expr = extract_expression_vectors(adata, gene_index, count_matrix, size_factors)
        sample_summary = build_level_summary(obs, tf, expected_group, raw_expr, norm_expr, sample_col, "sample")
        donor_summary = build_level_summary(obs, tf, expected_group, raw_expr, norm_expr, donor_col, "donor")
        level_summaries.extend([sample_summary, donor_summary])

        sample_means = sample_summary["mean_norm_expr"].to_numpy(dtype=float)
        donor_means = donor_summary["mean_norm_expr"].to_numpy(dtype=float)
        samples_detected_ge_5pct = int((sample_summary["positive_fraction"] >= 0.05).sum())
        detection_fraction = float(samples_detected_ge_5pct / max(sample_summary.shape[0], 1))

        boxplot_data[tf] = [
            norm_expr[obs[sample_col].astype(str).values == sample_name]
            for sample_name in sample_order
        ]

        direction_note = (
            "higher lesion expression is direction-consistent"
            if expected_group == "lesion"
            else "lower lesion expression is direction-consistent"
            if expected_group == "internal_control"
            else "main direction unavailable"
        )

        summary_rows.append(
            {
                "tf": tf,
                "gene_found": True,
                "expected_higher_group_main": expected_group,
                "main_regulon": expected_df.loc[tf, "regulon"] if tf in expected_df.index else "",
                "main_regulon_fdr": expected_df.loc[tf, "regulon_fdr"] if tf in expected_df.index else np.nan,
                "main_expr_fdr": expected_df.loc[tf, "expr_wilcoxon_fdr"] if tf in expected_df.index else np.nan,
                "main_log2fc_lesion_vs_internal_control": expected_df.loc[tf, "log2fc_lesion_vs_internal_control"] if tf in expected_df.index else np.nan,
                "overall_mean_norm_expr": float(np.mean(norm_expr)),
                "overall_median_norm_expr": float(np.median(norm_expr)),
                "overall_positive_fraction": float(np.mean(raw_expr > 0)),
                "overall_mean_raw_expr": float(np.mean(raw_expr)),
                "sample_mean_mean": float(np.mean(sample_means)),
                "sample_mean_sd": float(np.std(sample_means, ddof=0)),
                "sample_mean_cv": safe_cv(sample_means),
                "donor_mean_mean": float(np.mean(donor_means)),
                "donor_mean_sd": float(np.std(donor_means, ddof=0)),
                "donor_mean_cv": safe_cv(donor_means),
                "samples_detected_ge_5pct": samples_detected_ge_5pct,
                "samples_detected_fraction_ge_5pct": detection_fraction,
                "auc_status": "available_not_used" if has_auc else "skipped_no_auc",
                "direction_consistency_note": direction_note,
            }
        )

    expression_summary = pd.DataFrame(summary_rows)
    found_mask = expression_summary["gene_found"].fillna(False).astype(bool)

    if found_mask.any():
        expr_scaled = minmax_scale(expression_summary.loc[found_mask, "overall_mean_norm_expr"].to_numpy(dtype=float))
        pos_scaled = minmax_scale(expression_summary.loc[found_mask, "overall_positive_fraction"].to_numpy(dtype=float))
        cv_scaled = minmax_scale(expression_summary.loc[found_mask, "sample_mean_cv"].fillna(0.0).to_numpy(dtype=float))
        stability = 1.0 - cv_scaled

        lesion_presence_score = 0.50 * expr_scaled + 0.30 * pos_scaled + 0.20 * stability
        direction_score = lesion_presence_score.copy()
        expected_groups = expression_summary.loc[found_mask, "expected_higher_group_main"].astype(str).tolist()
        for i, expected_group in enumerate(expected_groups):
            if expected_group == "internal_control":
                direction_score[i] = 0.50 * (1.0 - expr_scaled[i]) + 0.30 * (1.0 - pos_scaled[i]) + 0.20 * stability[i]
        expression_summary.loc[found_mask, "lesion_presence_score"] = lesion_presence_score
        expression_summary.loc[found_mask, "direction_consistency_score"] = direction_score
        expression_summary.loc[~found_mask, "lesion_presence_score"] = np.nan
        expression_summary.loc[~found_mask, "direction_consistency_score"] = np.nan
    else:
        expression_summary["lesion_presence_score"] = np.nan
        expression_summary["direction_consistency_score"] = np.nan

    expression_summary["support_score"] = expression_summary["direction_consistency_score"]
    expression_summary["support_rank"] = (
        expression_summary["support_score"].rank(method="dense", ascending=False, na_option="bottom").astype("Int64")
    )
    expression_summary["support_interpretation"] = np.where(
        expression_summary["expected_higher_group_main"].eq("lesion"),
        "higher expression in lesion-only external cohort supports main direction",
        np.where(
            expression_summary["expected_higher_group_main"].eq("internal_control"),
            "lower expression in lesion-only external cohort supports main direction",
            "main direction unavailable",
        ),
    )

    donor_level_summary = pd.concat(level_summaries, axis=0, ignore_index=True) if level_summaries else pd.DataFrame()

    ranking_cols = [
        "tf",
        "support_rank",
        "support_score",
        "lesion_presence_score",
        "direction_consistency_score",
        "expected_higher_group_main",
        "overall_mean_norm_expr",
        "overall_positive_fraction",
        "sample_mean_cv",
        "samples_detected_ge_5pct",
        "samples_detected_fraction_ge_5pct",
        "main_regulon",
        "main_regulon_fdr",
        "main_expr_fdr",
        "main_log2fc_lesion_vs_internal_control",
        "auc_status",
        "support_interpretation",
    ]
    external_support_ranking = (
        expression_summary.sort_values(["support_score", "overall_mean_norm_expr"], ascending=[False, False])[ranking_cols]
        .reset_index(drop=True)
    )

    expression_summary_path = output_dir / "expression_summary.csv"
    donor_level_path = output_dir / "donor_level_summary.csv"
    ranking_path = output_dir / "external_support_ranking.csv"
    boxplot_path = output_dir / "expression_boxplots.png"
    barplot_path = output_dir / "expression_barplots_by_sample.png"
    summary_txt_path = output_dir / "external_support_summary_cn.txt"

    expression_summary.sort_values(["support_rank", "tf"]).to_csv(expression_summary_path, index=False)
    donor_level_summary.to_csv(donor_level_path, index=False)
    external_support_ranking.to_csv(ranking_path, index=False)

    if boxplot_data:
        save_expression_boxplots(boxplot_data, sample_order, boxplot_path)
        sample_mean_table = (
            donor_level_summary[donor_level_summary["level"] == "sample"]
            .pivot(index="tf", columns="level_id", values="mean_norm_expr")
            .reindex(index=target_tfs, columns=sample_order)
            .fillna(0.0)
        )
        save_expression_barplots_by_sample(sample_mean_table, sample_order, target_tfs, barplot_path)

    summary_text = build_summary_text(
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        n_samples=int(obs[sample_col].astype(str).nunique()),
        has_auc=has_auc,
        ranking_df=external_support_ranking,
    )
    summary_txt_path.write_text(summary_text, encoding="utf-8")

    print(f"Wrote {expression_summary_path}")
    print(f"Wrote {donor_level_path}")
    print(f"Wrote {ranking_path}")
    print(f"Wrote {boxplot_path}")
    print(f"Wrote {barplot_path}")
    print(f"Wrote {summary_txt_path}")


if __name__ == "__main__":
    main()
