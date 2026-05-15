#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats


DEFAULT_TFS = [
    "SOX2",
    "CEBPD",
    "BHLHE40",
    "ATF7",
    "NFE2L2",
    "HMGA1",
    "THRB",
    "TCF4",
    "RARB",
    "NR2F2",
    "SATB2",
]


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in str(value))


def parse_args() -> argparse.Namespace:
    root = project_root_from_file(__file__)
    parser = argparse.ArgumentParser(
        description="Screen CellOracle candidate TFs using existing pySCENIC differential regulon results and RNA expression."
    )
    parser.add_argument(
        "--h5ad",
        default=str(root / "final_exports" / "astrocyte_pilot_rna_with_pyscenic_auc.h5ad"),
        help="Merged h5ad used for pySCENIC AUC integration.",
    )
    parser.add_argument(
        "--expression-csv",
        default=str(root / "input" / "expression_for_pyscenic.csv"),
        help="Cell x gene expression matrix used for pySCENIC.",
    )
    parser.add_argument(
        "--group-compare-dir",
        default=str(root / "analysis_outputs" / "group_compare"),
        help="Directory containing current pySCENIC group comparison outputs.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "analysis_outputs" / "celloracle_candidates"),
        help="Output directory for CellOracle candidate summaries.",
    )
    parser.add_argument(
        "--candidate-tfs",
        default=",".join(DEFAULT_TFS),
        help="Comma-separated TF genes to evaluate.",
    )
    parser.add_argument(
        "--positive-fraction-threshold",
        type=float,
        default=0.05,
        help="Minimum positive-cell fraction in the higher-expression group for shortlist inclusion.",
    )
    parser.add_argument(
        "--expr-fdr-threshold",
        type=float,
        default=0.05,
        help="Expression Wilcoxon FDR threshold for shortlist inclusion.",
    )
    parser.add_argument(
        "--regulon-fdr-threshold",
        type=float,
        default=0.05,
        help="Regulon FDR threshold for shortlist inclusion.",
    )
    parser.add_argument(
        "--logfc-threshold",
        type=float,
        default=0.25,
        help="Absolute log2FC threshold for high-confidence shortlist tier.",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=1e-3,
        help="Pseudocount used for log2FC calculation on mean expression.",
    )
    return parser.parse_args()


def decode_strings(values: np.ndarray) -> list[str]:
    output: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            output.append(value.decode("utf-8"))
        else:
            output.append(str(value))
    return output


def read_obs_group_from_h5ad(h5ad_path: Path, group_column: str = "group") -> pd.DataFrame:
    with h5py.File(h5ad_path, "r") as handle:
        cell_ids = decode_strings(handle["obs"]["_index"][...])
        group_obj = handle["obs"][group_column]
        if isinstance(group_obj, h5py.Group) and "codes" in group_obj and "categories" in group_obj:
            categories = decode_strings(group_obj["categories"][...])
            codes = np.asarray(group_obj["codes"][...], dtype=int)
            groups = [categories[code] if code >= 0 else "NA" for code in codes]
        else:
            groups = decode_strings(group_obj[...])
    return pd.DataFrame({"CellID": cell_ids, group_column: groups})


def benjamini_hochberg(p_values: pd.Series | np.ndarray) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(values)
    if not finite_mask.any():
        return adjusted

    finite_values = values[finite_mask]
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    n = ranked.size
    out = np.empty(n, dtype=float)
    running = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        running = min(running, ranked[idx] * n / rank)
        out[idx] = running
    restored = np.empty(n, dtype=float)
    restored[order] = out
    adjusted[finite_mask] = np.clip(restored, 0.0, 1.0)
    return adjusted


def load_regulon_stats(group_compare_dir: Path) -> pd.DataFrame:
    stats_path = group_compare_dir / "regulon_group_statistics.csv"
    if not stats_path.exists():
        raise FileNotFoundError(f"Missing file: {stats_path}")
    stats_df = pd.read_csv(stats_path)
    if "mean_diff_group1_minus_group2" not in stats_df.columns:
        stats_df["mean_diff_group1_minus_group2"] = -stats_df["mean_diff_group2_minus_group1"]
    stats_df["effect_direction"] = np.where(
        stats_df["mean_diff_group1_minus_group2"] > 0,
        stats_df["group_1"],
        stats_df["group_2"],
    )
    rank_df = stats_df.sort_values(by=["fdr_bh", "abs_mean_diff"], ascending=[True, False]).reset_index()
    rank_df["regulon_abs_rank"] = np.arange(1, rank_df.shape[0] + 1)
    stats_df = stats_df.merge(rank_df[["index", "regulon_abs_rank"]], left_index=True, right_on="index", how="left")
    stats_df = stats_df.drop(columns=["index"]).reset_index(drop=True)
    return stats_df


def load_expression_subset(expression_csv: Path, candidate_tfs: list[str]) -> pd.DataFrame:
    usecols = ["CellID"] + candidate_tfs
    return pd.read_csv(expression_csv, usecols=usecols)


def match_regulon_row(stats_df: pd.DataFrame, tf: str) -> pd.Series | None:
    exact_candidates = [f"{tf}(+)", tf]
    for name in exact_candidates:
        hit = stats_df.loc[stats_df["regulon"] == name]
        if not hit.empty:
            return hit.iloc[0]
    prefix_hit = stats_df.loc[stats_df["regulon"].astype(str).str.startswith(f"{tf}(")]
    if not prefix_hit.empty:
        return prefix_hit.iloc[0]
    return None


def pick_higher_group(row: pd.Series, group_1: str, group_2: str) -> str:
    if row["mean_expr_group1"] > row["mean_expr_group2"]:
        return group_1
    if row["mean_expr_group2"] > row["mean_expr_group1"]:
        return group_2
    return "tie"


def build_shortlist_flags(
    results_df: pd.DataFrame,
    positive_fraction_threshold: float,
    expr_fdr_threshold: float,
    regulon_fdr_threshold: float,
    logfc_threshold: float,
) -> pd.DataFrame:
    df = results_df.copy()
    df["regulon_sig_pass"] = df["regulon_fdr"].fillna(1.0) < regulon_fdr_threshold
    df["expr_sig_pass"] = df["expr_wilcoxon_fdr"].fillna(1.0) < expr_fdr_threshold
    df["positive_fraction_pass"] = df["positive_fraction_in_higher_group"].fillna(0.0) >= positive_fraction_threshold
    df["direction_pass"] = df["direction_consistent_with_regulon"].fillna(False)
    df["logfc_pass"] = df["abs_log2fc_group1_vs_group2"].fillna(0.0) >= logfc_threshold

    df["shortlist_tier"] = "discard"
    medium_mask = (
        df["regulon_sig_pass"]
        & df["expr_sig_pass"]
        & df["positive_fraction_pass"]
        & df["direction_pass"]
    )
    high_mask = medium_mask & df["logfc_pass"]
    df.loc[medium_mask, "shortlist_tier"] = "medium"
    df.loc[high_mask, "shortlist_tier"] = "high"
    df["shortlist_keep"] = df["shortlist_tier"].isin(["high", "medium"])

    reasons = []
    for _, row in df.iterrows():
        tags = []
        if row["regulon_sig_pass"]:
            tags.append("regulon_sig")
        if row["expr_sig_pass"]:
            tags.append("expr_sig")
        if row["positive_fraction_pass"]:
            tags.append("expr_detected")
        if row["direction_pass"]:
            tags.append("direction_match")
        if row["logfc_pass"]:
            tags.append("logFC_strong")
        reasons.append(";".join(tags))
    df["shortlist_reason_tags"] = reasons
    return df


def build_chinese_summary(shortlist_df: pd.DataFrame, all_df: pd.DataFrame, group_1: str, group_2: str) -> str:
    keep_df = shortlist_df.loc[shortlist_df["shortlist_keep"]].copy()
    high_df = keep_df.loc[keep_df["shortlist_tier"] == "high"].copy()

    def fmt_rows(df: pd.DataFrame, top_n: int = 6) -> str:
        if df.empty:
            return "无"
        items = []
        for _, row in df.head(top_n).iterrows():
            items.append(
                f"{row['tf']}({row['expr_higher_group']}, log2FC={row['log2fc_group1_vs_group2']:.2f}, "
                f"exprFDR={row['expr_wilcoxon_fdr']:.2e}, regulonFDR={row['regulon_fdr']:.2e})"
            )
        return "、".join(items)

    lesion_like = keep_df.loc[keep_df["expr_higher_group"] == group_1].sort_values(
        by=["shortlist_tier", "expr_wilcoxon_fdr", "regulon_fdr", "abs_log2fc_group1_vs_group2"],
        ascending=[True, True, True, False],
    )
    control_like = keep_df.loc[keep_df["expr_higher_group"] == group_2].sort_values(
        by=["shortlist_tier", "expr_wilcoxon_fdr", "regulon_fdr", "abs_log2fc_group1_vs_group2"],
        ascending=[True, True, True, False],
    )

    line1 = (
        f"本次共检查 {all_df.shape[0]} 个候选 TF，其中 {keep_df.shape[0]} 个进入最终 shortlist，"
        f"其中高置信度 {high_df.shape[0]} 个。"
    )
    line2 = (
        f"shortlist 的判定标准为：regulon FDR 显著、TF 表达 Wilcoxon FDR 显著、在高表达组中具有足够阳性细胞比例，"
        f"且 TF 表达方向与 regulon 活性方向一致。"
    )
    line3 = f"{group_1} 侧优先候选包括：{fmt_rows(lesion_like)}。"
    line4 = f"{group_2} 侧优先候选包括：{fmt_rows(control_like)}。"
    line5 = "这些 TF 可以作为后续 CellOracle GRN 方向性验证和扰动模拟的优先起点。"
    return "\n".join([line1, line2, line3, line4, line5])


def main() -> None:
    args = parse_args()

    h5ad_path = Path(args.h5ad).resolve()
    expression_csv = Path(args.expression_csv).resolve()
    group_compare_dir = Path(args.group_compare_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    candidate_tfs = [item.strip() for item in args.candidate_tfs.split(",") if item.strip()]
    print(f"[1/6] Reading metadata from h5ad: {h5ad_path}")
    meta = read_obs_group_from_h5ad(h5ad_path, group_column="group")
    groups = meta["group"].value_counts()
    if not {"lesion", "internal_control"}.issubset(set(groups.index)):
        raise ValueError(f"Expected lesion/internal_control in group labels, got: {sorted(groups.index.tolist())}")

    group_1 = "lesion"
    group_2 = "internal_control"
    meta = meta.loc[meta["group"].isin([group_1, group_2])].copy()

    print(f"[2/6] Reading regulon statistics: {group_compare_dir}")
    regulon_stats = load_regulon_stats(group_compare_dir)

    print(f"[3/6] Reading expression subset for {len(candidate_tfs)} TFs")
    expr = load_expression_subset(expression_csv, candidate_tfs)
    if not expr["CellID"].equals(meta["CellID"]):
        expr = expr.merge(meta[["CellID"]], on="CellID", how="inner")
    merged = meta.merge(expr, on="CellID", how="inner")
    if merged.shape[0] != meta.shape[0]:
        raise ValueError(f"Cell merge mismatch: metadata={meta.shape[0]}, merged={merged.shape[0]}")

    missing_tfs = [tf for tf in candidate_tfs if tf not in merged.columns]
    if missing_tfs:
        raise ValueError(f"Missing TFs in expression matrix: {missing_tfs}")

    print("[4/6] Calculating expression metrics and Wilcoxon tests")
    rows: list[dict[str, object]] = []
    mask_group1 = merged["group"] == group_1
    mask_group2 = merged["group"] == group_2

    for tf in candidate_tfs:
        x = merged.loc[mask_group1, tf].to_numpy(dtype=float)
        y = merged.loc[mask_group2, tf].to_numpy(dtype=float)

        mean_group1 = float(np.mean(x))
        mean_group2 = float(np.mean(y))
        median_group1 = float(np.median(x))
        median_group2 = float(np.median(y))
        pos_group1 = float(np.mean(x > 0))
        pos_group2 = float(np.mean(y > 0))
        log2fc = float(np.log2((mean_group1 + args.pseudocount) / (mean_group2 + args.pseudocount)))

        if np.allclose(x, x[0]) and np.allclose(y, y[0]) and np.isclose(x[0], y[0]):
            stat = 0.0
            p_value = 1.0
        else:
            stat, p_value = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")

        regulon_row = match_regulon_row(regulon_stats, tf)
        regulon_name = regulon_row["regulon"] if regulon_row is not None else None
        regulon_fdr = float(regulon_row["fdr_bh"]) if regulon_row is not None else np.nan
        regulon_p = float(regulon_row["p_value"]) if regulon_row is not None else np.nan
        regulon_effect_group = str(regulon_row["effect_direction"]) if regulon_row is not None else None
        regulon_mean_diff = (
            float(regulon_row["mean_diff_group1_minus_group2"]) if regulon_row is not None else np.nan
        )
        regulon_rank = int(regulon_row["regulon_abs_rank"]) if regulon_row is not None else np.nan
        expr_higher_group = group_1 if mean_group1 > mean_group2 else group_2 if mean_group2 > mean_group1 else "tie"
        pos_in_higher_group = pos_group1 if expr_higher_group == group_1 else pos_group2 if expr_higher_group == group_2 else max(pos_group1, pos_group2)
        direction_match = bool(regulon_effect_group == expr_higher_group) if regulon_effect_group not in [None, "tie"] else False

        rows.append(
            {
                "tf": tf,
                "regulon": regulon_name,
                "regulon_rank_by_fdr_absdiff": regulon_rank,
                "regulon_effect_group": regulon_effect_group,
                "regulon_mean_diff_lesion_minus_internal_control": regulon_mean_diff,
                "regulon_p_value": regulon_p,
                "regulon_fdr": regulon_fdr,
                "n_lesion": int(x.shape[0]),
                "n_internal_control": int(y.shape[0]),
                "mean_expr_group1": mean_group1,
                "mean_expr_group2": mean_group2,
                "median_expr_group1": median_group1,
                "median_expr_group2": median_group2,
                "mean_expr_lesion": mean_group1,
                "mean_expr_internal_control": mean_group2,
                "positive_fraction_group1": pos_group1,
                "positive_fraction_group2": pos_group2,
                "positive_fraction_lesion": pos_group1,
                "positive_fraction_internal_control": pos_group2,
                "log2fc_group1_vs_group2": log2fc,
                "log2fc_lesion_vs_internal_control": log2fc,
                "abs_log2fc_group1_vs_group2": abs(log2fc),
                "expr_wilcoxon_statistic": float(stat),
                "expr_wilcoxon_p": float(p_value),
                "expr_higher_group": expr_higher_group,
                "positive_fraction_in_higher_group": pos_in_higher_group,
                "direction_consistent_with_regulon": direction_match,
            }
        )

    results_df = pd.DataFrame(rows)
    results_df["expr_wilcoxon_fdr"] = benjamini_hochberg(results_df["expr_wilcoxon_p"].to_numpy(dtype=float))
    results_df = build_shortlist_flags(
        results_df,
        positive_fraction_threshold=args.positive_fraction_threshold,
        expr_fdr_threshold=args.expr_fdr_threshold,
        regulon_fdr_threshold=args.regulon_fdr_threshold,
        logfc_threshold=args.logfc_threshold,
    )
    tier_rank = {"high": 0, "medium": 1, "discard": 2}
    results_df["shortlist_tier_rank"] = results_df["shortlist_tier"].map(tier_rank).fillna(9).astype(int)

    results_df = results_df.sort_values(
        by=["shortlist_keep", "shortlist_tier_rank", "expr_wilcoxon_fdr", "regulon_fdr", "abs_log2fc_group1_vs_group2"],
        ascending=[False, True, True, True, False],
    ).reset_index(drop=True)
    results_df = results_df.drop(columns=["shortlist_tier_rank"])

    print("[5/6] Writing candidate tables")
    all_metrics_path = output_dir / "celloracle_candidate_tf_metrics.csv"
    shortlist_path = output_dir / "celloracle_tf_shortlist.csv"
    high_conf_path = output_dir / "celloracle_tf_shortlist_high_confidence.csv"
    summary_txt_path = output_dir / "celloracle_tf_shortlist_summary_cn.txt"
    summary_json_path = output_dir / "celloracle_tf_shortlist_summary.json"

    results_df.to_csv(all_metrics_path, index=False)
    results_df.loc[results_df["shortlist_keep"]].to_csv(shortlist_path, index=False)
    results_df.loc[results_df["shortlist_tier"] == "high"].to_csv(high_conf_path, index=False)

    summary_text = build_chinese_summary(results_df, results_df, group_1=group_1, group_2=group_2)
    with open(summary_txt_path, "w", encoding="utf-8") as handle:
        handle.write(summary_text + "\n")

    summary_payload = {
        "group_1": group_1,
        "group_2": group_2,
        "candidate_tfs": candidate_tfs,
        "positive_fraction_threshold": args.positive_fraction_threshold,
        "expr_fdr_threshold": args.expr_fdr_threshold,
        "regulon_fdr_threshold": args.regulon_fdr_threshold,
        "logfc_threshold": args.logfc_threshold,
        "n_all_candidates": int(results_df.shape[0]),
        "n_shortlist": int(results_df["shortlist_keep"].sum()),
        "n_high_confidence": int((results_df["shortlist_tier"] == "high").sum()),
        "outputs": [
            all_metrics_path.name,
            shortlist_path.name,
            high_conf_path.name,
            summary_txt_path.name,
        ],
    }
    with open(summary_json_path, "w", encoding="utf-8") as handle:
        json.dump(summary_payload, handle, indent=2, ensure_ascii=False)

    print("[6/6] Done")
    print(f"  n_candidates={results_df.shape[0]}")
    print(f"  n_shortlist={int(results_df['shortlist_keep'].sum())}")
    print(f"  n_high_confidence={int((results_df['shortlist_tier'] == 'high').sum())}")
    print(f"  outputs={output_dir}")


if __name__ == "__main__":
    main()
