#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats


TARGET_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
PRIORITY_TFS = ["NFE2L2", "THRB"]
EPS = 1e-9


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_file(__file__)
    parser = argparse.ArgumentParser(
        description="Supportive external evidence review for priority TFs and regulon activity using existing pySCENIC + CellOracle round-1 findings."
    )
    parser.add_argument(
        "--h5ad",
        default=str(root / "external_validation" / "input" / "external_validation.h5ad"),
        help="Supportive external h5ad path.",
    )
    parser.add_argument(
        "--auc-csv",
        default=str(root / "external_validation" / "input" / "external_validation_auc.csv"),
        help="Optional external regulon AUC CSV. Used if h5ad does not contain AUC in obsm/uns.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "external_validation"),
        help="Output root directory.",
    )
    parser.add_argument(
        "--group-column",
        default=None,
        help="Optional obs column to force as grouping variable.",
    )
    parser.add_argument(
        "--obsm-key",
        default="X_pyscenic_auc",
        help="AUC matrix key in adata.obsm if present.",
    )
    parser.add_argument(
        "--regulon-name-key",
        default="pyscenic_regulon_names",
        help="Regulon names key in adata.uns if present.",
    )
    parser.add_argument(
        "--target-tfs",
        default=",".join(TARGET_TFS),
        help="Comma-separated TF list.",
    )
    parser.add_argument(
        "--priority-tfs",
        default=",".join(PRIORITY_TFS),
        help="Comma-separated priority TF list.",
    )
    parser.add_argument(
        "--expr-pseudocount",
        type=float,
        default=1e-3,
        help="Pseudocount for expression log2FC.",
    )
    return parser.parse_args()


PREFERRED_NAME_SCORES = [
    ("disease", 220),
    ("diagnosis", 220),
    ("condition", 210),
    ("status", 200),
    ("group", 200),
    ("cohort", 170),
    ("case", 165),
    ("control", 165),
    ("treat", 160),
    ("response", 150),
    ("subtype", 140),
    ("sample", 120),
    ("patient", 110),
    ("donor", 100),
    ("cell_type", 90),
    ("celltype", 90),
    ("cluster", 80),
    ("batch", 40),
]
EXCLUDED_NAME_SCORES = [
    ("barcode", 400),
    ("join_key", 400),
    ("cell_id", 250),
    ("obs_name", 250),
    ("index", 150),
]


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


def _series_value_counts(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("NA").value_counts(dropna=False)


def _name_score(column_name: str) -> tuple[int, list[str]]:
    name = column_name.lower()
    score = 0
    reasons: list[str] = []
    for pattern, value in PREFERRED_NAME_SCORES:
        if pattern in name:
            score += value
            reasons.append(f"name:+{pattern}")
    for pattern, value in EXCLUDED_NAME_SCORES:
        if pattern in name:
            score -= value
            reasons.append(f"name:-{pattern}")
    return score, reasons


def score_group_column(column_name: str, series: pd.Series) -> tuple[float, bool, list[str]]:
    total_n = int(series.shape[0])
    non_na = series.dropna()
    n_missing = int(series.isna().sum())
    n_non_na = int(non_na.shape[0])
    n_unique = int(non_na.nunique())
    counts = _series_value_counts(non_na) if n_non_na else pd.Series(dtype="int64")
    min_count = int(counts.min()) if not counts.empty else 0
    max_count = int(counts.max()) if not counts.empty else 0
    balance_ratio = (min_count / max_count) if max_count else 0.0
    score, reasons = _name_score(column_name)

    if pd.api.types.is_bool_dtype(series.dtype):
        score += 45
        reasons.append("dtype:bool")
    elif pd.api.types.is_numeric_dtype(series.dtype):
        if n_unique <= 12:
            score += 20
            reasons.append("dtype:low_card_numeric")
        else:
            score -= 60
            reasons.append("dtype:high_card_numeric")
    else:
        score += 35
        reasons.append("dtype:categorical_like")

    if n_unique < 2:
        score -= 1000
    elif n_unique == n_non_na:
        score -= 1000
    elif n_unique <= 2:
        score += 80
    elif n_unique <= 4:
        score += 65
    elif n_unique <= 8:
        score += 45
    elif n_unique <= 12:
        score += 25
    elif n_unique <= 20:
        score += 10
    else:
        score -= 200

    if min_count >= 100:
        score += 30
    elif min_count >= 30:
        score += 20
    elif min_count >= 10:
        score += 10
    elif min_count > 0:
        score -= 10

    score += 30 * balance_ratio
    if n_missing == 0:
        score += 8
    else:
        score -= min(30, n_missing / max(total_n, 1) * 100)

    eligible = (
        n_unique >= 2
        and n_unique < max(20, total_n)
        and n_unique != n_non_na
        and "name:-barcode" not in reasons
        and "name:-join_key" not in reasons
    )
    return float(score), eligible, reasons


def choose_group_column(obs: pd.DataFrame, preferred_column: str | None) -> tuple[str, pd.DataFrame]:
    rows = []
    for col in obs.columns:
        series = obs[col]
        score, eligible, reasons = score_group_column(col, series)
        rows.append(
            {
                "column": col,
                "dtype": str(series.dtype),
                "n_unique_non_na": int(series.dropna().nunique()),
                "group_candidate_score": round(score, 3),
                "is_group_candidate": bool(eligible),
                "value_preview": " | ".join(f"{idx}:{int(v)}" for idx, v in _series_value_counts(series).head(6).items()),
                "score_reasons": "; ".join(reasons),
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        by=["group_candidate_score", "n_unique_non_na", "column"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    if preferred_column:
        if preferred_column not in obs.columns:
            raise KeyError(f"Requested group column not found: {preferred_column}")
        return preferred_column, summary
    candidates = summary.loc[summary["is_group_candidate"]]
    if candidates.empty:
        return str(summary.iloc[0]["column"]), summary
    return str(candidates.iloc[0]["column"]), summary


def choose_two_group_order(groups: pd.Series) -> list[str]:
    levels = list(groups.astype("string").dropna().unique())
    level_lower = {level: str(level).lower() for level in levels}
    disease_tokens = ["lesion", "disease", "case", "treated", "fcd", "tumor"]
    control_tokens = ["internal_control", "control", "normal", "non_lesion", "healthy", "untreated"]

    disease_hits = [lvl for lvl in levels if any(token in level_lower[lvl] for token in disease_tokens)]
    control_hits = [lvl for lvl in levels if any(token in level_lower[lvl] for token in control_tokens)]
    if len(levels) == 2 and disease_hits and control_hits:
        disease = disease_hits[0]
        control = control_hits[0]
        if disease != control:
            return [disease, control]

    counts = groups.value_counts(dropna=False)
    ordered = list(counts.index.astype(str))
    if len(ordered) >= 2:
        return ordered[:2]
    return levels


def get_expected_group_map(root: Path) -> pd.DataFrame:
    ko_path = root / "celloracle_run" / "ko_round1" / "round1_ko_master_table.csv"
    cand_path = root / "analysis_outputs" / "celloracle_candidates" / "celloracle_candidate_tf_metrics.csv"
    if ko_path.exists():
        df = pd.read_csv(ko_path)
        cols = [c for c in ["tf", "target_group", "expected_regulon_effect_group", "expected_expr_higher_group", "recovery_index"] if c in df.columns]
        return df[cols].copy()
    if cand_path.exists():
        df = pd.read_csv(cand_path)
        cols = [c for c in ["tf", "regulon_effect_group", "expr_higher_group"] if c in df.columns]
        out = df[cols].copy()
        if "regulon_effect_group" in out.columns:
            out["target_group"] = out["regulon_effect_group"]
        return out
    return pd.DataFrame(columns=["tf", "target_group"])


def read_adata(h5ad_path: Path):
    try:
        import anndata as ad
    except Exception as exc:
        raise RuntimeError("This script requires anndata. Run it in the pyscenic-analysis Docker image.") from exc
    adata = ad.read_h5ad(h5ad_path)
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    return adata


def extract_expression_df(adata, genes: list[str]) -> tuple[pd.DataFrame, list[str], str]:
    available = []
    source = "X"
    if getattr(adata, "raw", None) is not None:
        raw_var_names = pd.Index(adata.raw.var_names.astype(str))
        available = [gene for gene in genes if gene in raw_var_names]
        if available:
            matrix = adata.raw[:, available].X
            source = "raw"
            if sparse.issparse(matrix):
                matrix = matrix.toarray()
            expr_df = pd.DataFrame(np.asarray(matrix), index=adata.obs_names, columns=available)
            return expr_df, available, source

    available = [gene for gene in genes if gene in pd.Index(adata.var_names.astype(str))]
    if available:
        matrix = adata[:, available].X
        if sparse.issparse(matrix):
            matrix = matrix.toarray()
        expr_df = pd.DataFrame(np.asarray(matrix), index=adata.obs_names, columns=available)
        return expr_df, available, source

    return pd.DataFrame(index=adata.obs_names), [], source


def extract_auc_df(adata, obsm_key: str, regulon_name_key: str, auc_csv: Path | None) -> tuple[pd.DataFrame, str]:
    if obsm_key in adata.obsm and regulon_name_key in adata.uns:
        matrix = adata.obsm[obsm_key]
        if hasattr(matrix, "toarray"):
            matrix = matrix.toarray()
        matrix = np.asarray(matrix)
        regulons = list(adata.uns[regulon_name_key])
        if matrix.ndim == 2 and matrix.shape[1] == len(regulons):
            return pd.DataFrame(matrix, index=adata.obs_names, columns=regulons), f"h5ad:{obsm_key}"

    if auc_csv is not None and auc_csv.exists():
        auc = pd.read_csv(auc_csv)
        auc = auc.rename(columns={auc.columns[0]: "CellID"})
        auc["CellID"] = auc["CellID"].astype(str)
        auc = auc.drop_duplicates(subset="CellID").set_index("CellID")
        shared = [cell for cell in adata.obs_names if cell in auc.index]
        if shared:
            return auc.loc[shared].copy(), f"csv:{auc_csv.name}"

    return pd.DataFrame(index=adata.obs_names), "missing"


def find_matching_regulon(tf: str, auc_df: pd.DataFrame) -> str | None:
    candidates = [f"{tf}(+)", tf]
    for name in candidates:
        if name in auc_df.columns:
            return name
    prefix = [col for col in auc_df.columns if str(col).startswith(f"{tf}(")]
    if prefix:
        return str(prefix[0])
    return None


def compare_feature(values: pd.Series, groups: pd.Series, group_order: list[str], pseudocount: float) -> dict[str, object]:
    g1, g2 = group_order
    x = values.loc[groups == g1].to_numpy(dtype=float)
    y = values.loc[groups == g2].to_numpy(dtype=float)
    mean1 = float(np.mean(x))
    mean2 = float(np.mean(y))
    median1 = float(np.median(x))
    median2 = float(np.median(y))
    pos1 = float(np.mean(x > 0))
    pos2 = float(np.mean(y > 0))
    if np.allclose(x, x[0]) and np.allclose(y, y[0]) and np.isclose(x[0], y[0]):
        stat = 0.0
        p_value = 1.0
    else:
        stat, p_value = stats.mannwhitneyu(x, y, alternative="two-sided", method="auto")

    observed_higher_group = g1 if mean1 > mean2 else g2 if mean2 > mean1 else "tie"
    return {
        "group_1": g1,
        "group_2": g2,
        "n_group_1": int(x.shape[0]),
        "n_group_2": int(y.shape[0]),
        "mean_group_1": mean1,
        "mean_group_2": mean2,
        "median_group_1": median1,
        "median_group_2": median2,
        "positive_fraction_group_1": pos1,
        "positive_fraction_group_2": pos2,
        "mean_diff_group1_minus_group2": mean1 - mean2,
        "abs_mean_diff": abs(mean1 - mean2),
        "log2fc_group1_vs_group2": float(np.log2((mean1 + pseudocount) / (mean2 + pseudocount))),
        "observed_higher_group": observed_higher_group,
        "statistic": float(stat),
        "p_value": float(p_value),
    }


def plot_group_bars(df: pd.DataFrame, value_col_1: str, value_col_2: str, title: str, ylabel: str, output_png: Path, colors: tuple[str, str]) -> None:
    labels = df["tf"].astype(str).tolist()
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    ax.bar(x - width / 2, df[value_col_1].to_numpy(dtype=float), width=width, color=colors[0], label=df["group_1"].iloc[0])
    ax.bar(x + width / 2, df[value_col_2].to_numpy(dtype=float), width=width, color=colors[1], label=df["group_2"].iloc[0])
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_boxpanels(score_frames: dict[str, pd.DataFrame], group_order: list[str], title_prefix: str, output_png: Path, value_col: str) -> None:
    n = len(score_frames)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(10.5, 4.4 * nrows), squeeze=False)
    palette = ["#d95f02", "#1b7fb8"]

    items = list(score_frames.items())
    for idx, (tf, df) in enumerate(items):
        ax = axes[idx // ncols][idx % ncols]
        data = [df.loc[df["group"] == group, value_col].to_numpy(dtype=float) for group in group_order]
        box = ax.boxplot(data, labels=group_order, patch_artist=True, showfliers=False)
        for patch, color in zip(box["boxes"], palette):
            patch.set_facecolor(color)
            patch.set_alpha(0.82)
        ax.set_title(f"{title_prefix}: {tf}")
        ax.set_ylabel(value_col)
        ax.tick_params(axis="x", rotation=20)

    for idx in range(len(items), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_support_scatter(master_df: pd.DataFrame, output_png: Path) -> None:
    colors = ["#c44e52" if tf in PRIORITY_TFS else "#4c72b0" for tf in master_df["tf"]]
    sizes = 500 * (master_df["support_score"].clip(lower=0.25).to_numpy(dtype=float) + 0.6)
    fig, ax = plt.subplots(figsize=(8.4, 6.2))
    ax.scatter(
        master_df["expr_abs_log2fc"].to_numpy(dtype=float),
        master_df["regulon_abs_mean_diff"].fillna(0.0).to_numpy(dtype=float),
        s=sizes,
        c=colors,
        alpha=0.78,
        edgecolors="black",
        linewidths=0.8,
    )
    for _, row in master_df.iterrows():
        ax.text(row["expr_abs_log2fc"], row["regulon_abs_mean_diff"] if pd.notna(row["regulon_abs_mean_diff"]) else 0.0, row["tf"], fontsize=9, ha="left", va="bottom")
    ax.set_xlabel("Expression |log2FC|")
    ax.set_ylabel("Regulon |mean diff|")
    ax.set_title("Supportive external evidence landscape")
    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#c44e52", markeredgecolor="black", label="Priority TF (NFE2L2 / THRB)", markersize=9),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#4c72b0", markeredgecolor="black", label="Other TF", markersize=9),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def build_chinese_summary(master_df: pd.DataFrame, group_order: list[str]) -> str:
    g1, g2 = group_order
    ordered = master_df.sort_values(by=["support_score", "priority_score", "expr_fdr", "regulon_fdr"], ascending=[False, False, True, True]).reset_index(drop=True)
    priority_df = ordered.loc[ordered["tf"].isin(PRIORITY_TFS)].copy()

    def support_text(row: pd.Series) -> str:
        pieces = []
        if row["expr_support_label"] != "not_supported":
            pieces.append(f"表达{row['expr_support_label']}")
        if row["regulon_support_label"] != "not_supported":
            pieces.append(f"regulon{row['regulon_support_label']}")
        if not pieces:
            return "表达和 regulon 均未形成一致支持"
        return "，".join(pieces)

    top_priority = priority_df.iloc[0] if not priority_df.empty else ordered.iloc[0]
    second_priority = priority_df.iloc[1] if priority_df.shape[0] > 1 else ordered.iloc[min(1, ordered.shape[0] - 1)]
    lesion_side = ordered.loc[ordered["expected_group"] == g1, "tf"].tolist()
    control_side = ordered.loc[ordered["expected_group"] == g2, "tf"].tolist()

    line1 = (
        f"在外部验证队列中，本轮优先检查了 {', '.join(TARGET_TFS)} 的表达与 regulon activity。"
        f"从综合支持度看，{top_priority['tf']} 的外部验证证据最强，{second_priority['tf']} 次之。"
    )
    line2 = (
        f"优先验证的 {PRIORITY_TFS[0]} 和 {PRIORITY_TFS[1]} 中，{PRIORITY_TFS[0]} 表现为 {support_text(master_df.loc[master_df['tf'] == PRIORITY_TFS[0]].iloc[0])}；"
        f"{PRIORITY_TFS[1]} 表现为 {support_text(master_df.loc[master_df['tf'] == PRIORITY_TFS[1]].iloc[0])}。"
    )
    line3 = (
        f"按当前内外部结果的方向预期，{g1} 侧候选主要包括 {', '.join(lesion_side) if lesion_side else '无'}，"
        f"{g2} 侧候选主要包括 {', '.join(control_side) if control_side else '无'}。"
    )
    line4 = (
        "如果外部队列表达与 regulon activity 同时显著且方向一致，可视为较强再现；"
        "若仅表达或仅 regulon activity 再现，则视为部分支持，后续应结合样本来源与细胞组成差异解释。"
    )
    line5 = (
        f"基于当前结果，建议优先继续围绕 {top_priority['tf']} 与 {second_priority['tf']} 做下一步外部验证和结果整合，"
        "其余 TF 可作为补充分析对象。"
    )
    return "\n".join([line1, line2, line3, line4, line5])


def build_support_labels(expr_supported: bool, reg_supported: bool) -> tuple[str, str]:
    expr_label = "supported" if expr_supported else "not_supported"
    reg_label = "supported" if reg_supported else "not_supported"
    return expr_label, reg_label


def main() -> None:
    args = parse_args()
    root = project_root_from_file(__file__)

    h5ad_path = Path(args.h5ad).resolve()
    auc_csv_path = Path(args.auc_csv).resolve() if args.auc_csv else None
    output_dir = Path(args.output_dir).resolve()
    plot_dir = output_dir / "summary_plots"
    input_dir = output_dir / "input"
    plot_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)

    target_tfs = [tf.strip() for tf in args.target_tfs.split(",") if tf.strip()]
    priority_tfs = [tf.strip() for tf in args.priority_tfs.split(",") if tf.strip()]

    if not h5ad_path.exists():
        raise FileNotFoundError(
            f"Supportive external h5ad not found: {h5ad_path}\n"
            "Place the supportive external AnnData file there or pass --h5ad explicitly."
        )

    print(f"[1/5] Reading supportive external h5ad: {h5ad_path}")
    adata = read_adata(h5ad_path)
    group_column, obs_summary = choose_group_column(adata.obs, preferred_column=args.group_column)
    groups = adata.obs[group_column].astype("string")
    group_order = choose_two_group_order(groups)
    if len(group_order) != 2:
        raise ValueError(f"Current script expects exactly 2 major groups, got: {group_order}")
    groups = groups.loc[groups.isin(group_order)].copy()
    adata = adata[groups.index].copy()
    groups = groups.loc[adata.obs_names]

    print("[2/5] Extracting expression and regulon activity")
    expr_df, available_genes, expr_source = extract_expression_df(adata, target_tfs)
    auc_df, auc_source = extract_auc_df(adata, args.obsm_key, args.regulon_name_key, auc_csv_path)
    expected_map = get_expected_group_map(root)

    obs_summary.to_csv(output_dir / "external_obs_group_candidates.csv", index=False)

    expression_rows = []
    expression_score_frames: dict[str, pd.DataFrame] = {}
    for tf in target_tfs:
        row = {"tf": tf, "feature_type": "expression", "feature_name": tf, "available": tf in expr_df.columns}
        if tf in expr_df.columns:
            stats_row = compare_feature(expr_df[tf], groups, group_order, args.expr_pseudocount)
            row.update(stats_row)
            score_df = pd.DataFrame({"value": expr_df[tf].to_numpy(dtype=float), "group": groups.to_numpy()})
            expression_score_frames[tf] = score_df
        expression_rows.append(row)
    expression_stats = pd.DataFrame(expression_rows)
    if "p_value" in expression_stats.columns:
        expression_stats["fdr_bh"] = benjamini_hochberg(expression_stats["p_value"])
    else:
        expression_stats["fdr_bh"] = np.nan

    regulon_rows = []
    regulon_score_frames: dict[str, pd.DataFrame] = {}
    for tf in target_tfs:
        regulon_name = find_matching_regulon(tf, auc_df) if not auc_df.empty else None
        row = {"tf": tf, "feature_type": "regulon", "feature_name": regulon_name, "available": regulon_name is not None}
        if regulon_name is not None:
            stats_row = compare_feature(auc_df[regulon_name], groups, group_order, args.expr_pseudocount)
            row.update(stats_row)
            score_df = pd.DataFrame({"value": auc_df[regulon_name].to_numpy(dtype=float), "group": groups.to_numpy()})
            regulon_score_frames[tf] = score_df
        regulon_rows.append(row)
    regulon_stats = pd.DataFrame(regulon_rows)
    if "p_value" in regulon_stats.columns:
        regulon_stats["fdr_bh"] = benjamini_hochberg(regulon_stats["p_value"])
    else:
        regulon_stats["fdr_bh"] = np.nan

    print("[3/5] Building master supportive evidence table")
    expr_renamed = expression_stats.rename(
        columns={
            "feature_name": "expr_feature_name",
            "available": "expr_available",
            "mean_group_1": "expr_mean_group_1",
            "mean_group_2": "expr_mean_group_2",
            "median_group_1": "expr_median_group_1",
            "median_group_2": "expr_median_group_2",
            "positive_fraction_group_1": "expr_positive_fraction_group_1",
            "positive_fraction_group_2": "expr_positive_fraction_group_2",
            "mean_diff_group1_minus_group2": "expr_mean_diff_group1_minus_group2",
            "abs_mean_diff": "expr_abs_mean_diff",
            "log2fc_group1_vs_group2": "expr_log2fc_group1_vs_group2",
            "observed_higher_group": "expr_observed_higher_group",
            "p_value": "expr_p_value",
            "fdr_bh": "expr_fdr",
        }
    )
    reg_renamed = regulon_stats.rename(
        columns={
            "feature_name": "regulon_feature_name",
            "available": "regulon_available",
            "mean_group_1": "regulon_mean_group_1",
            "mean_group_2": "regulon_mean_group_2",
            "median_group_1": "regulon_median_group_1",
            "median_group_2": "regulon_median_group_2",
            "positive_fraction_group_1": "regulon_positive_fraction_group_1",
            "positive_fraction_group_2": "regulon_positive_fraction_group_2",
            "mean_diff_group1_minus_group2": "regulon_mean_diff_group1_minus_group2",
            "abs_mean_diff": "regulon_abs_mean_diff",
            "log2fc_group1_vs_group2": "regulon_log2fc_group1_vs_group2",
            "observed_higher_group": "regulon_observed_higher_group",
            "p_value": "regulon_p_value",
            "fdr_bh": "regulon_fdr",
        }
    )
    master = expr_renamed.merge(reg_renamed.drop(columns=[c for c in ["group_1", "group_2", "n_group_1", "n_group_2", "feature_type"] if c in reg_renamed.columns]), on="tf", how="outer")

    if "group_1" not in master.columns:
        master["group_1"] = group_order[0]
        master["group_2"] = group_order[1]
    master["priority_tier"] = np.where(master["tf"].isin(priority_tfs), "priority", "secondary")
    master["priority_score"] = np.where(master["tf"].isin(priority_tfs), 2, 1)

    expected_map["tf"] = expected_map["tf"].astype(str)
    master = master.merge(expected_map, on="tf", how="left")
    if "target_group" in master.columns:
        master["expected_group"] = master["target_group"]
    elif "expected_regulon_effect_group" in master.columns:
        master["expected_group"] = master["expected_regulon_effect_group"]
    else:
        master["expected_group"] = np.nan

    master["expr_direction_matches_expected"] = master["expr_observed_higher_group"] == master["expected_group"]
    master["regulon_direction_matches_expected"] = master["regulon_observed_higher_group"] == master["expected_group"]
    master["expr_supported"] = master["expr_available"].fillna(False) & (master["expr_fdr"].fillna(1.0) < 0.05) & master["expr_direction_matches_expected"].fillna(False)
    master["regulon_supported"] = master["regulon_available"].fillna(False) & (master["regulon_fdr"].fillna(1.0) < 0.05) & master["regulon_direction_matches_expected"].fillna(False)
    master["expr_support_label"], master["regulon_support_label"] = zip(*[build_support_labels(e, r) for e, r in zip(master["expr_supported"], master["regulon_supported"])])
    master["support_score"] = (
        master["expr_supported"].astype(int) * 2
        + master["regulon_supported"].astype(int) * 2
        + master["priority_score"]
        + master["expr_direction_matches_expected"].fillna(False).astype(int) * 0.5
        + master["regulon_direction_matches_expected"].fillna(False).astype(int) * 0.5
    )
    master["expr_abs_log2fc"] = master["expr_log2fc_group1_vs_group2"].abs()
    master["regulon_abs_mean_diff"] = master["regulon_abs_mean_diff"]
    master = master.sort_values(by=["support_score", "priority_score", "expr_fdr", "regulon_fdr", "expr_abs_log2fc"], ascending=[False, False, True, True, False]).reset_index(drop=True)

    expression_stats.to_csv(output_dir / "external_tf_expression_group_stats.csv", index=False)
    regulon_stats.to_csv(output_dir / "external_regulon_activity_group_stats.csv", index=False)
    master.to_csv(output_dir / "external_validation_master_table.csv", index=False)

    support_rank = master[[
        "tf",
        "priority_tier",
        "expected_group",
        "expr_supported",
        "regulon_supported",
        "support_score",
        "expr_fdr",
        "regulon_fdr",
        "expr_observed_higher_group",
        "regulon_observed_higher_group",
    ]].copy()
    support_rank.to_csv(output_dir / "external_validation_support_ranking.csv", index=False)

    print("[4/5] Writing plots")
    expr_plot_df = expression_stats.loc[expression_stats["available"]].copy()
    reg_plot_df = regulon_stats.loc[regulon_stats["available"]].copy()
    if not expr_plot_df.empty:
        plot_group_bars(
            df=expr_plot_df,
            value_col_1="mean_group_1",
            value_col_2="mean_group_2",
            title="Supportive external evidence: TF expression group means",
            ylabel="Mean expression",
            output_png=plot_dir / "external_tf_expression_group_means.png",
            colors=("#d95f02", "#1b7fb8"),
        )
    if not reg_plot_df.empty:
        plot_group_bars(
            df=reg_plot_df,
            value_col_1="mean_group_1",
            value_col_2="mean_group_2",
            title="Supportive external evidence: regulon activity group means",
            ylabel="Mean regulon activity (AUC)",
            output_png=plot_dir / "external_regulon_activity_group_means.png",
            colors=("#d95f02", "#1b7fb8"),
        )
    if expression_score_frames:
        plot_boxpanels(
            score_frames=expression_score_frames,
            group_order=group_order,
            title_prefix="Expression",
            output_png=plot_dir / "external_tf_expression_boxplots.png",
            value_col="value",
        )
    if regulon_score_frames:
        plot_boxpanels(
            score_frames=regulon_score_frames,
            group_order=group_order,
            title_prefix="Regulon activity",
            output_png=plot_dir / "external_regulon_activity_boxplots.png",
            value_col="value",
        )
    plot_support_scatter(master, plot_dir / "external_validation_support_scatter.png")

    print("[5/5] Writing reports")
    input_report = {
        "h5ad_path": str(h5ad_path),
        "auc_source": auc_source,
        "expression_source": expr_source,
        "group_column": group_column,
        "group_order": group_order,
        "available_genes": available_genes,
        "available_regulons": regulon_stats.loc[regulon_stats["available"], "feature_name"].astype(str).tolist(),
        "obsm_keys": list(adata.obsm.keys()),
        "obs_columns": list(adata.obs.columns),
    }
    with open(output_dir / "external_validation_input_report.json", "w", encoding="utf-8") as handle:
        json.dump(input_report, handle, indent=2, ensure_ascii=False)
    with open(output_dir / "external_validation_summary_cn.txt", "w", encoding="utf-8") as handle:
        handle.write(build_chinese_summary(master, group_order) + "\n")

    outputs = {
        "main_tables": [
            "external_tf_expression_group_stats.csv",
            "external_regulon_activity_group_stats.csv",
            "external_validation_master_table.csv",
            "external_validation_support_ranking.csv",
        ],
        "plots": [
            "summary_plots/external_tf_expression_group_means.png",
            "summary_plots/external_regulon_activity_group_means.png",
            "summary_plots/external_tf_expression_boxplots.png",
            "summary_plots/external_regulon_activity_boxplots.png",
            "summary_plots/external_validation_support_scatter.png",
        ],
        "reports": [
            "external_validation_input_report.json",
            "external_validation_summary_cn.txt",
        ],
    }
    with open(output_dir / "external_validation_outputs.json", "w", encoding="utf-8") as handle:
        json.dump(outputs, handle, indent=2, ensure_ascii=False)

    print(f"  group_column={group_column}")
    print(f"  auc_source={auc_source}")
    print(f"  outputs={output_dir}")


if __name__ == "__main__":
    main()
