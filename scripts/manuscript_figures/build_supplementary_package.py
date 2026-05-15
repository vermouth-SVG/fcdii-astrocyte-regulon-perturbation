#!/usr/bin/env python3
from __future__ import annotations

import textwrap
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "manuscript_output" / "figures_supplementary"
TABLE_DIR = ROOT / "manuscript_output" / "tables_supplementary"
MAIN_TABLE_DIR = ROOT / "manuscript_output" / "tables_main"
INDEX_PATH = ROOT / "manuscript_output" / "SUPPLEMENTARY_PACKAGE_INDEX.md"

GROUP_COLORS = {"lesion": "#C45A4D", "internal_control": "#4D78A8"}
TF_COLORS = {
    "NFE2L2": "#C45A4D",
    "THRB": "#4D78A8",
    "BHLHE40": "#C98A2E",
    "SOX2": "#8A7FA8",
    "SATB2": "#6F8C63",
    "RARB": "#5E8D87",
    "HMGA1": "#9A6F52",
}
DRUG_LAYER_COLORS = {
    "headline_cns_mechanism_direction_leads": "#B44E43",
    "supportive_after_bbb_leads": "#C79A43",
    "peripheral_program_modulating_clues": "#7E8895",
}
TF_ORDER_MAIN = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
TF_ORDER_EXT = ["NFE2L2", "THRB", "BHLHE40", "SOX2", "SATB2", "RARB", "HMGA1"]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.0,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.0,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.4,
            "legend.fontsize": 6.2,
            "figure.titlesize": 10.0,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.65,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_csv(path: str | Path, **kwargs) -> pd.DataFrame:
    full = ROOT / path if not isinstance(path, Path) else path
    return pd.read_csv(full, **kwargs)


def safe_read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8").strip() if path.exists() else ""


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")


def save_figure(fig: plt.Figure, stem: str) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png = FIG_DIR / f"{stem}.png"
    pdf = FIG_DIR / f"{stem}.pdf"
    fig.savefig(png, dpi=360, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return [png, pdf]


def style_axes(ax: plt.Axes, grid_axis: str | None = "y") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color="#E7E7E7", linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.13, 1.05, label, transform=ax.transAxes, fontsize=10.4, fontweight="bold", va="top", ha="left")


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def neg_log10(values: pd.Series | np.ndarray) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").astype(float)
    x = np.where(np.isfinite(x) & (x > 0), x, np.nan)
    finite = x[np.isfinite(x)]
    floor = max(float(np.nanmin(finite)) * 0.1, 1e-300) if finite.size else 1e-300
    x = np.where(np.isfinite(x), x, floor)
    return -np.log10(x)


def clean_first_col(df: pd.DataFrame, name: str) -> pd.DataFrame:
    first = df.columns[0]
    if first.startswith("Unnamed") or first in {"", "H1"}:
        df = df.rename(columns={first: name})
    return df


def parse_simple_kv_text(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("["):
            continue
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
        elif ":" in line:
            k, v = line.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def normalize_qc_sample_prefix(value: str) -> str:
    text = str(value)
    text = re.sub(r"_max\d+$", "", text)
    return text


def clean_mechanism_theme(value: object) -> str:
    text = str(value).strip() if pd.notna(value) else "unassigned"
    mapping = {
        "未自动识别；需人工核实": "manual_review_needed",
        "未自动识别；补充轴线索": "secondary_axis_clue",
    }
    text = mapping.get(text, text or "unassigned")
    if any(ord(ch) > 127 for ch in text):
        text = "non_english_theme"
    return text


def trim_white(img: np.ndarray, threshold: float = 0.985) -> np.ndarray:
    arr = np.asarray(img)
    if arr.ndim == 2:
        mask = arr < threshold
    else:
        mask = np.any(arr[..., :3] < threshold, axis=2)
    coords = np.argwhere(mask)
    if coords.size == 0:
        return arr
    y0, x0 = coords.min(axis=0)[:2]
    y1, x1 = coords.max(axis=0)[:2] + 1
    return arr[y0:y1, x0:x1]


def short_layer_name(layer: str) -> str:
    mapping = {
        "headline_cns_mechanism_direction_leads": "headline",
        "supportive_after_bbb_leads": "supportive",
        "peripheral_program_modulating_clues": "peripheral",
    }
    return mapping.get(str(layer), str(layer))


def short_sample_label(x: str) -> str:
    return str(x).replace("GSE268807_", "").replace("_FL", "").replace(".", "_")


def load_metadata() -> pd.DataFrame:
    meta = read_csv("celloracle_run/prepared_data/celloracle_round1_metadata.csv")
    first = meta.columns[0]
    if first.startswith("Unnamed"):
        meta = meta.rename(columns={first: "cell_id"})
    return meta


def build_table1() -> Path:
    merge = parse_simple_kv_text(ROOT / "final_exports" / "merge_report.txt")
    build140 = parse_simple_kv_text(ROOT / "external_validation" / "input" / "input_build_report.txt")
    build190 = parse_simple_kv_text(ROOT / "external_validation_round2" / "input" / "build_report.txt")
    df = pd.DataFrame(
        [
            {
                "Dataset": "GSE268807 astrocyte pilot",
                "Role": "discovery cohort",
                "Object/subset": f"4-sample astrocyte pilot with pySCENIC AUC-matched object ({merge.get('h5ad_cells', '2322')} cells; {merge.get('auc_regulons', '105')} regulons)",
                "Grouping": "lesion vs internal_control",
                "Current use": "current main discovery cohort",
                "Interpretation level / note": "Primary source for the current manuscript mainline, including differential regulon analysis, CellOracle, and sample-level robustness.",
            },
            {
                "Dataset": "GSE140393",
                "Role": "single-group supportive analysis",
                "Object/subset": f"lesion-only external supportive object ({build140.get('cells', '11399')} cells)",
                "Grouping": "lesion only",
                "Current use": "supportive external analysis",
                "Interpretation level / note": "Expression-level supportive evidence only; no matched internal_control group and not formal external validation.",
            },
            {
                "Dataset": "GSE190452",
                "Role": "cross-syndrome supportive analysis",
                "Object/subset": f"TLE vs non-epileptic control external supportive object ({build190.get('cells', '69968')} cells)",
                "Grouping": "TLE vs non-epileptic control",
                "Current use": "supportive external analysis",
                "Interpretation level / note": "Cross-syndrome supportive evidence with expression and projected regulon support; not formal same-disease external validation.",
            },
            {
                "Dataset": "GSE275302",
                "Role": "not included/private",
                "Object/subset": "planned candidate external dataset; not incorporated into current local result chain",
                "Grouping": "not available in current project package",
                "Current use": "not included",
                "Interpretation level / note": "Kept out of the current manuscript because the project documentation marks it as private/unavailable at the time of route consolidation.",
            },
            {
                "Dataset": "SCENIC+/Kaggle old exploration route",
                "Role": "archived exploration",
                "Object/subset": "legacy proxy/Kaggle/cloud exploration artifacts and historical route notes",
                "Grouping": "not applicable",
                "Current use": "archived only",
                "Interpretation level / note": "Retained for provenance and audit only; no SCENIC+ current mainline and no current manuscript conclusions depend on this route.",
            },
        ]
    )
    out = MAIN_TABLE_DIR / "Table1_cohort_and_analytical_objects.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    return out


def build_table_supplementary() -> list[Path]:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    s1 = read_csv("analysis_outputs/group_compare/regulon_group_statistics.csv").sort_values(["fdr_bh", "abs_mean_diff"], ascending=[True, False])
    p1 = TABLE_DIR / "TableS1_full_differential_regulons.csv"
    s1.to_csv(p1, index=False)
    outputs.append(p1)

    s2 = read_csv("analysis_outputs/celloracle_candidates/celloracle_candidate_tf_metrics.csv").sort_values(
        ["shortlist_keep", "regulon_fdr", "expr_wilcoxon_fdr"], ascending=[False, True, True]
    )
    p2 = TABLE_DIR / "TableS2_candidate_tf_shortlist_extended.csv"
    s2.to_csv(p2, index=False)
    outputs.append(p2)

    s3 = read_csv("celloracle_run/ko_round1/round1_ko_master_table.csv").sort_values("mean_shift_length", ascending=False)
    p3 = TABLE_DIR / "TableS3_celloracle_round1_full_metrics.csv"
    s3.to_csv(p3, index=False)
    outputs.append(p3)

    robust_parts: list[pd.DataFrame] = []
    for rel, component in [
        ("robustness_validation/02_leave_one_sample_out_master_results_revised.csv", "leave_one_sample_out_master"),
        ("robustness_validation/02_leave_one_sample_out_tf_summary_revised.csv", "leave_one_sample_out_tf_summary"),
        ("robustness_validation/03_pseudobulk_candidate_tf_table.csv", "pseudobulk_candidate_tf"),
        ("robustness_validation/03_pseudobulk_regulon_table.csv", "pseudobulk_regulon"),
        ("robustness_validation/04_shortlist_priority_stability.csv", "shortlist_priority_stability"),
        ("robustness_validation/04_shortlist_membership_frequency.csv", "shortlist_membership_frequency"),
        ("robustness_validation/05_robustness_validation_integrated_table.csv", "integrated_robustness"),
    ]:
        df = read_csv(rel)
        df.insert(0, "analysis_component", component)
        df.insert(1, "source_file", rel.replace("\\", "/"))
        robust_parts.append(df)
    s4 = pd.concat(robust_parts, ignore_index=True, sort=False)
    p4 = TABLE_DIR / "TableS4_sample_level_robustness_full.csv"
    s4.to_csv(p4, index=False)
    outputs.append(p4)

    functional_parts: list[pd.DataFrame] = []
    for rel, component in [
        ("functional_interpretation/03_go_enrichment_all.csv", "GO_enrichment"),
        ("functional_interpretation/03_kegg_enrichment_all.csv", "KEGG_enrichment"),
        ("functional_interpretation/04_program_convergence_table.csv", "program_convergence"),
        ("functional_interpretation/04_tf_integrated_evidence_table.csv", "integrated_evidence"),
        ("functional_interpretation/05_functional_interpretation_integrated_table.csv", "functional_priority_summary"),
    ]:
        df = read_csv(rel)
        df.insert(0, "analysis_component", component)
        df.insert(1, "source_file", rel.replace("\\", "/"))
        functional_parts.append(df)
    s5 = pd.concat(functional_parts, ignore_index=True, sort=False)
    p5 = TABLE_DIR / "TableS5_functional_interpretation_program_convergence_full.csv"
    s5.to_csv(p5, index=False)
    outputs.append(p5)

    external_parts: list[pd.DataFrame] = []
    for rel, dataset, component in [
        ("external_validation/single_group_support/expression_summary.csv", "GSE140393", "single_group_expression_summary"),
        ("external_validation/single_group_support/donor_level_summary.csv", "GSE140393", "single_group_sample_donor_summary"),
        ("external_validation/single_group_support/external_support_ranking.csv", "GSE140393", "single_group_support_ranking"),
        ("external_validation_round2/external_round2_tf_expression_group_stats.csv", "GSE190452", "round2_group_expression"),
        ("external_validation_round2/external_round2_tf_expression_celltype_stats.csv", "GSE190452", "round2_cluster_expression"),
        ("external_validation_round2/external_round2_support_ranking.csv", "GSE190452", "round2_expression_support_ranking"),
        ("external_validation_round2/round2_regulon_group_stats.csv", "GSE190452", "round2_regulon_group"),
        ("external_validation_round2/round2_regulon_cluster_stats.csv", "GSE190452", "round2_regulon_cluster"),
        ("external_validation_round2/round2_regulon_support_ranking.csv", "GSE190452", "round2_regulon_support_ranking"),
        ("external_validation_round2/round2_expression_and_regulon_master_table.csv", "GSE190452", "round2_expression_regulon_master"),
    ]:
        df = read_csv(rel)
        df.insert(0, "dataset", dataset)
        df.insert(1, "analysis_component", component)
        df.insert(2, "source_file", rel.replace("\\", "/"))
        external_parts.append(df)
    s6 = pd.concat(external_parts, ignore_index=True, sort=False)
    p6 = TABLE_DIR / "TableS6_supportive_external_analysis_full.csv"
    s6.to_csv(p6, index=False)
    outputs.append(p6)

    s7 = read_csv("drug_repositioning/13_integrated_after_bbb_v21.csv")
    p7 = TABLE_DIR / "TableS7_drug_repositioning_v21_after_bbb_full.csv"
    s7.to_csv(p7, index=False)
    outputs.append(p7)

    return outputs


def build_figS1() -> list[Path]:
    meta = load_metadata()
    merge = parse_simple_kv_text(ROOT / "final_exports" / "merge_report.txt")
    sample_counts = (
        meta.groupby(["clinical_sample_label", "group"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["group", "clinical_sample_label"])
    )
    sample_order = sample_counts["clinical_sample_label"].tolist()

    qc_frames = []
    for path in sorted((ROOT / "analysis_outputs" / "gse268807_astrocyte_pilot" / "qc_metrics").glob("qc_metrics_*.csv")):
        df = pd.read_csv(path)
        df["sample_prefix"] = normalize_qc_sample_prefix(path.stem.replace("qc_metrics_", ""))
        qc_frames.append(df)
    qc = pd.concat(qc_frames, ignore_index=True)
    prefix_to_label = meta[["sample_prefix", "clinical_sample_label", "group"]].drop_duplicates()
    qc = qc.merge(prefix_to_label, on="sample_prefix", how="left")
    qc["clinical_sample_label"] = qc["clinical_sample_label"].fillna(qc["sample_prefix"])
    inferred_group = pd.Series(
        np.where(qc["sample_prefix"].astype(str).str.contains("_D_"), "lesion", "internal_control"),
        index=qc.index,
    )
    qc["group"] = qc["group"].where(qc["group"].notna(), inferred_group)
    qc = qc.sort_values(["group", "clinical_sample_label"])

    heat = meta.groupby(["clinical_sample_label", "group"]).size().unstack(fill_value=0).reindex(index=sample_order, columns=["internal_control", "lesion"]).fillna(0).astype(int)

    fig = plt.figure(figsize=(12.2, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[0.9, 1.05])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    ax_a.bar(
        sample_counts["clinical_sample_label"],
        sample_counts["n_cells"],
        color=[GROUP_COLORS.get(g, "#777777") for g in sample_counts["group"]],
        width=0.62,
    )
    for i, row in sample_counts.reset_index(drop=True).iterrows():
        ax_a.text(i, row["n_cells"] + 28, str(int(row["n_cells"])), ha="center", va="bottom", fontsize=6.0)
    ax_a.set_ylabel("Cells in final analysis object")
    ax_a.set_title("Final 4-sample astrocyte pilot composition", pad=6)
    ax_a.tick_params(axis="x", rotation=26)
    style_axes(ax_a)
    ax_a.legend(
        handles=[
            Line2D([0], [0], color=GROUP_COLORS["internal_control"], linewidth=5, label="internal_control"),
            Line2D([0], [0], color=GROUP_COLORS["lesion"], linewidth=5, label="lesion"),
        ],
        frameon=False,
        loc="upper right",
    )
    panel_label(ax_a, "A")

    stages = ["n_cells_input", "n_cells_after_qc", "n_cells_after_subsample"]
    x = np.arange(len(stages))
    for _, row in qc.iterrows():
        vals = [float(row[s]) for s in stages]
        label = str(row["clinical_sample_label"])
        color = GROUP_COLORS.get(str(row["group"]), "#777777")
        ax_b.plot(x, vals, marker="o", linewidth=1.8, color=color, alpha=0.95)
        ax_b.text(x[-1] + 0.04, vals[-1], label, fontsize=6.0, va="center", color=color)
    ax_b.set_xlim(-0.1, 2.55)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(["input", "after QC", "after subsample"])
    ax_b.set_ylabel("Cells")
    ax_b.set_title("Upstream QC and subsampling trace", pad=6)
    style_axes(ax_b)
    panel_label(ax_b, "B")

    im = ax_c.imshow(heat.values.T, aspect="auto", cmap="Blues", vmin=0, vmax=max(1, int(heat.values.max())))
    ax_c.set_xticks(np.arange(heat.shape[0]))
    ax_c.set_xticklabels(heat.index, rotation=26, ha="right")
    ax_c.set_yticks(np.arange(heat.shape[1]))
    ax_c.set_yticklabels(heat.columns)
    for i in range(heat.shape[1]):
        for j in range(heat.shape[0]):
            val = int(heat.values[j, i])
            ax_c.text(j, i, str(val), ha="center", va="center", fontsize=6.1, color="white" if val > 400 else "#222222")
    for spine in ax_c.spines.values():
        spine.set_visible(False)
    ax_c.set_title("Final sample-by-group cell matrix", pad=6)
    fig.colorbar(im, ax=ax_c, fraction=0.05, pad=0.03, label="cells")
    panel_label(ax_c, "C")

    ax_d.axis("off")
    boxes = [
        (0.02, 0.62, 0.28, 0.22, "Raw pilot input", "GSE268807 4 samples\nsample-level pilot subset"),
        (0.36, 0.62, 0.28, 0.22, "QC / subsample", "max 3000 cells per source\nper-sample QC traces retained"),
        (0.70, 0.62, 0.28, 0.22, "Merged RNA object", f"{merge.get('h5ad_cells', '2322')} cells x\n{merge.get('h5ad_genes', '36601')} genes"),
        (0.36, 0.20, 0.28, 0.22, "AUC-matched export", f"{merge.get('auc_regulons', '105')} regulons\n0 missing matched cells"),
    ]
    for x0, y0, w, h, title, body in boxes:
        rect = Rectangle((x0, y0), w, h, transform=ax_d.transAxes, facecolor="white", edgecolor="#B8B8B8", linewidth=0.9)
        ax_d.add_patch(rect)
        ax_d.text(x0 + 0.02, y0 + h - 0.05, title, transform=ax_d.transAxes, fontsize=6.7, fontweight="bold", ha="left", va="top")
        ax_d.text(x0 + 0.02, y0 + 0.05, body, transform=ax_d.transAxes, fontsize=6.2, ha="left", va="bottom")
    for start, end in [((0.30, 0.73), (0.36, 0.73)), ((0.64, 0.73), (0.70, 0.73)), ((0.50, 0.62), (0.50, 0.42))]:
        arr = FancyArrowPatch(start, end, transform=ax_d.transAxes, arrowstyle="-|>", mutation_scale=10, linewidth=0.9, color="#888888")
        ax_d.add_patch(arr)
    ax_d.text(0.02, 0.06, "Input lineage and analytical object trace", transform=ax_d.transAxes, fontsize=7.4, fontweight="bold", ha="left")
    ax_d.text(0.02, 0.00, "Current mainline object: GSE268807 astrocyte pilot -> pySCENIC AUCell-matched object -> CellOracle-ready metadata", transform=ax_d.transAxes, fontsize=6.0, ha="left")
    panel_label(ax_d, "D")

    fig.suptitle("Fig. S1. Discovery object QC and input overview", y=1.01)
    return save_figure(fig, "FigS1_discovery_object_qc_and_input_overview")


def build_figS2() -> list[Path]:
    stats = read_csv("analysis_outputs/group_compare/regulon_group_statistics.csv")
    stats["lesion_minus_internal_control"] = as_num(stats["mean_lesion"]) - as_num(stats["mean_internal_control"])
    stats["neg_log10_fdr"] = neg_log10(stats["fdr_bh"])
    stats["regulon"] = stats["regulon"].astype(str)
    meta = load_metadata()
    auc = read_csv("final_exports/auc_mtx_matched_to_h5ad.csv").rename(columns=lambda c: "cell_id" if str(c).startswith("Unnamed") else c)
    merged = auc.merge(meta[["cell_id", "clinical_sample_label", "group"]], on="cell_id", how="left")

    top40 = stats.sort_values(["fdr_bh", "abs_mean_diff"], ascending=[True, False]).head(40)["regulon"].tolist()
    sample_means = merged.groupby("clinical_sample_label")[top40].mean()
    sample_order_meta = meta[["clinical_sample_label", "group"]].drop_duplicates().sort_values(["group", "clinical_sample_label"])
    sample_order = sample_order_meta["clinical_sample_label"].tolist()
    sample_groups = sample_order_meta["group"].tolist()
    heat = sample_means.loc[sample_order].T
    heat = heat.sub(heat.mean(axis=1), axis=0)
    heat = heat.div(heat.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0.0)
    heat = heat.loc[stats.set_index("regulon").loc[top40].sort_values("lesion_minus_internal_control").index.tolist()]

    top_ic = stats[stats["lesion_minus_internal_control"] < 0].nsmallest(15, "lesion_minus_internal_control").copy()
    top_ls = stats[stats["lesion_minus_internal_control"] > 0].nlargest(15, "lesion_minus_internal_control").copy().sort_values("lesion_minus_internal_control")

    fig = plt.figure(figsize=(12.6, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.08, 0.92], height_ratios=[0.9, 1.1])
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    sub = gs[1, 1].subgridspec(1, 2, wspace=0.28)
    ax_c = fig.add_subplot(sub[0, 0])
    ax_d = fig.add_subplot(sub[0, 1])

    im = ax_a.imshow(heat.values, aspect="auto", cmap="RdBu_r", vmin=-1.7, vmax=1.7)
    ax_a.set_xticks(np.arange(len(sample_order)))
    ax_a.set_xticklabels(sample_order, rotation=28, ha="right")
    ax_a.set_yticks(np.arange(len(heat.index)))
    ax_a.set_yticklabels(heat.index)
    ax_a.tick_params(length=0)
    ax_a.axvline(1.5, color="#D6D6D6", linewidth=0.9)
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    strip = ax_a.inset_axes([0.0, 1.01, 1.0, 0.06])
    strip.set_xlim(0, len(sample_order))
    strip.set_ylim(0, 1)
    for i, group in enumerate(sample_groups):
        strip.add_patch(Rectangle((i, 0), 1, 1, facecolor=GROUP_COLORS[group], edgecolor="white", linewidth=0.6))
    strip.text(1.0, 0.5, "internal_control", color="white", ha="center", va="center", fontsize=6.0, fontweight="bold")
    strip.text(3.0, 0.5, "lesion", color="white", ha="center", va="center", fontsize=6.0, fontweight="bold")
    strip.axis("off")
    ax_a.set_title("Extended differential regulon heatmap (top 40)", pad=6)
    fig.colorbar(im, ax=ax_a, fraction=0.05, pad=0.03, label="row z-score")
    panel_label(ax_a, "A")

    ax_b.scatter(
        stats["lesion_minus_internal_control"],
        stats["neg_log10_fdr"],
        s=15,
        c=np.where(stats["lesion_minus_internal_control"] >= 0, GROUP_COLORS["lesion"], GROUP_COLORS["internal_control"]),
        alpha=0.65,
        edgecolors="none",
    )
    ax_b.axvline(0, color="#9A9A9A", linewidth=0.8)
    for reg in ["NFE2L2(+)", "THRB(+)", "BHLHE40(+)", "SOX2(+)", "SATB2(+)", "RARB(+)", "HMGA1(+)"]:
        row = stats.loc[stats["regulon"].eq(reg)]
        if row.empty:
            continue
        xx = float(row["lesion_minus_internal_control"].iloc[0])
        yy = float(row["neg_log10_fdr"].iloc[0])
        tf = reg.replace("(+)", "")
        ax_b.scatter([xx], [yy], s=42, color="white", edgecolor=TF_COLORS.get(tf, "#333333"), linewidth=1.0, zorder=4)
        ax_b.annotate(reg, (xx, yy), xytext=(4 if xx >= 0 else -4, 4), textcoords="offset points", ha="left" if xx >= 0 else "right", fontsize=6.2)
    ax_b.set_xlabel("Mean AUC difference")
    ax_b.set_ylabel("-log10(FDR)")
    ax_b.set_title("All-regulon differential summary", pad=6)
    style_axes(ax_b)
    panel_label(ax_b, "B")

    ax_c.hlines(top_ic["regulon"], xmin=top_ic["lesion_minus_internal_control"], xmax=0, color="#D0D0D0", linewidth=1.0)
    ax_c.scatter(top_ic["lesion_minus_internal_control"], top_ic["regulon"], color=GROUP_COLORS["internal_control"], s=22)
    ax_c.axvline(0, color="#9A9A9A", linewidth=0.8)
    ax_c.set_xlabel("Mean AUC difference")
    ax_c.set_title("Top internal_control-associated regulons", pad=6)
    style_axes(ax_c, grid_axis="x")
    panel_label(ax_c, "C")

    ax_d.hlines(top_ls["regulon"], xmin=0, xmax=top_ls["lesion_minus_internal_control"], color="#D0D0D0", linewidth=1.0)
    ax_d.scatter(top_ls["lesion_minus_internal_control"], top_ls["regulon"], color=GROUP_COLORS["lesion"], s=22)
    ax_d.axvline(0, color="#9A9A9A", linewidth=0.8)
    ax_d.set_xlabel("Mean AUC difference")
    ax_d.set_title("Top lesion-associated regulons", pad=6)
    style_axes(ax_d, grid_axis="x")
    panel_label(ax_d, "D")

    fig.suptitle("Fig. S2. Full differential regulon landscape", y=1.01)
    return save_figure(fig, "FigS2_full_differential_regulon_landscape")


def build_figS3() -> list[Path]:
    ranking = read_csv("celloracle_run/ko_round1/round1_ko_master_table.csv").set_index("tf").loc[TF_ORDER_MAIN].copy()
    ranking["mean_shift_length"] = as_num(ranking["mean_shift_length"])
    ranking["recovery_index"] = as_num(ranking["recovery_index"])
    ranking["mean_net_shift"] = as_num(ranking["mean_net_shift"])
    ranking["lesion_mean_shift_length"] = as_num(ranking["lesion_mean_shift_length"])
    ranking["internal_control_mean_shift_length"] = as_num(ranking["internal_control_mean_shift_length"])
    ranking["direction_agreement_factor"] = as_num(ranking["direction_agreement_factor"])

    metric_cols = [
        "mean_shift_length",
        "recovery_index",
        "mean_net_shift",
        "lesion_mean_shift_length",
        "internal_control_mean_shift_length",
        "direction_agreement_factor",
    ]
    heat = ranking[metric_cols].copy()
    heat_norm = (heat - heat.min()) / (heat.max() - heat.min()).replace(0, np.nan)
    heat_norm = heat_norm.fillna(0.0)

    fig = plt.figure(figsize=(12.6, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[0.92, 1.05, 1.05], wspace=0.08)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1])

    im = ax_a.imshow(heat_norm.values, aspect="auto", cmap="OrRd", vmin=0, vmax=1)
    ax_a.set_xticks(np.arange(len(metric_cols)))
    ax_a.set_xticklabels(
        ["shift", "Recovery\nindex", "net shift", "lesion\nshift", "internal_control\nshift", "direction\nagreement"],
        rotation=25,
        ha="right",
    )
    ax_a.set_yticks(np.arange(len(TF_ORDER_MAIN)))
    ax_a.set_yticklabels(TF_ORDER_MAIN)
    for i, tf in enumerate(TF_ORDER_MAIN):
        for j, col in enumerate(metric_cols):
            ax_a.text(j, i, f"{heat.loc[tf, col]:.2f}", ha="center", va="center", fontsize=5.8, color="white" if heat_norm.loc[tf, col] > 0.55 else "#222222")
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    ax_a.set_title("Round1 CellOracle metric heatmap", pad=6)
    fig.colorbar(im, ax=ax_a, fraction=0.05, pad=0.03, label="column-scaled")
    panel_label(ax_a, "A")

    ax_b.scatter(
        ranking["mean_shift_length"],
        ranking["recovery_index"],
        s=90 + 380 * ranking["mean_net_shift"] / ranking["mean_net_shift"].max(),
        c=[TF_COLORS[tf] for tf in TF_ORDER_MAIN],
        alpha=0.85,
        edgecolor="white",
        linewidth=0.8,
    )
    for tf in TF_ORDER_MAIN:
        ax_b.annotate(tf, (float(ranking.loc[tf, "mean_shift_length"]), float(ranking.loc[tf, "recovery_index"])), xytext=(4, 4), textcoords="offset points", fontsize=6.4)
    ax_b.set_xlabel("Mean shift length")
    ax_b.set_ylabel("Recovery index")
    ax_b.set_title("Quantitative perturbation summary", pad=6)
    style_axes(ax_b)
    panel_label(ax_b, "B")

    for ax, tf, label in [(ax_c, "NFE2L2", "C"), (ax_d, "THRB", "D"), (ax_e, "BHLHE40", "E"), (ax_f, "SOX2", "F")]:
        img = trim_white(mpimg.imread(ROOT / "celloracle_run" / "ko_round1" / tf / f"{tf}_quiver.png"))
        ax.imshow(img)
        ax.set_title(f"{tf} representative perturbation", pad=6)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        panel_label(ax, label)

    fig.suptitle("Fig. S3. Extended CellOracle prioritization", y=1.01)
    return save_figure(fig, "FigS3_extended_celloracle_prioritization")


def build_figS4() -> list[Path]:
    loso = read_csv("robustness_validation/02_leave_one_sample_out_master_results_revised.csv")
    tf_summary = read_csv("robustness_validation/02_leave_one_sample_out_tf_summary_revised.csv")
    pseudo_expr = read_csv("robustness_validation/03_pseudobulk_candidate_tf_table.csv")
    pseudo_reg = read_csv("robustness_validation/03_pseudobulk_regulon_table.csv")
    sens = read_csv("robustness_validation/04_shortlist_sensitivity_grid.csv")
    meta = load_metadata()
    id_map = meta[["sample_id", "clinical_sample_label"]].drop_duplicates()

    loso = loso.merge(id_map, left_on="dropped_sample", right_on="sample_id", how="left")
    loso["drop_label"] = np.where(loso["dropped_sample"].eq("FULL_DATA"), "Full", loso["clinical_sample_label"].fillna(loso["dropped_sample"]))
    order_cols = ["Full"] + id_map.sort_values("clinical_sample_label")["clinical_sample_label"].tolist()
    tf_order = TF_ORDER_EXT

    expr_mat = loso.pivot_table(index="tf", columns="drop_label", values="expr_log2fc", aggfunc="first").reindex(index=tf_order, columns=order_cols)
    reg_mat = loso.pivot_table(index="tf", columns="drop_label", values="regulon_diff", aggfunc="first").reindex(index=tf_order, columns=order_cols)

    pseudo = pseudo_expr.merge(
        pseudo_reg[["tf", "mean_diff_lesion_minus_internal_control"]].rename(columns={"mean_diff_lesion_minus_internal_control": "regulon_mean_diff"}),
        on="tf",
        how="left",
    )
    pseudo = pseudo[pseudo["tf"].isin(tf_order)].set_index("tf").loc[tf_order].reset_index()

    sens_rule = sens[sens["rule"].eq("regulon_significant_plus_expression_significant_plus_direction")].copy()
    sens_rule["regulon_fdr_threshold"] = as_num(sens_rule["regulon_fdr_threshold"])
    sens_rule["expr_abs_log2fc_threshold"] = as_num(sens_rule["expr_abs_log2fc_threshold"])
    sens_rule["positive_fraction_threshold"] = as_num(sens_rule["positive_fraction_threshold"])

    fig = plt.figure(figsize=(13.0, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[1.0, 1.02], wspace=0.08)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    outer_d = fig.add_subplot(gs[1, 1])
    outer_d.axis("off")
    sub = gs[1, 1].subgridspec(1, 3, wspace=0.25)
    ax_d1 = fig.add_subplot(sub[0, 0])
    ax_d2 = fig.add_subplot(sub[0, 1], sharey=ax_d1)
    ax_d3 = fig.add_subplot(sub[0, 2], sharey=ax_d1)

    im1 = ax_a.imshow(expr_mat.values, aspect="auto", cmap="RdBu_r", vmin=-1.3, vmax=1.3)
    ax_a.set_xticks(np.arange(len(order_cols)))
    ax_a.set_xticklabels(order_cols, rotation=28, ha="right")
    ax_a.set_yticks(np.arange(len(tf_order)))
    ax_a.set_yticklabels(tf_order)
    ax_a.set_title("Leave-one-sample-out expression effect", pad=6)
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    fig.colorbar(im1, ax=ax_a, fraction=0.05, pad=0.03, label="log2FC")
    panel_label(ax_a, "A")

    im2 = ax_b.imshow(reg_mat.values, aspect="auto", cmap="RdBu_r", vmin=-0.28, vmax=0.28)
    ax_b.set_xticks(np.arange(len(order_cols)))
    ax_b.set_xticklabels(order_cols, rotation=28, ha="right")
    ax_b.set_yticks(np.arange(len(tf_order)))
    ax_b.set_yticklabels(tf_order)
    ax_b.set_title("Leave-one-sample-out regulon effect", pad=6)
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    fig.colorbar(im2, ax=ax_b, fraction=0.05, pad=0.03, label="AUC difference")
    panel_label(ax_b, "B")

    ax_c.axvline(0, color="#A0A0A0", linewidth=0.8)
    ax_c.axhline(0, color="#D6D6D6", linewidth=0.8)
    ax_c.scatter(
        pseudo["log2fc_lesion_vs_internal_control"],
        pseudo["regulon_mean_diff"],
        s=70,
        c=[TF_COLORS.get(tf, "#777777") for tf in pseudo["tf"]],
        edgecolor="white",
        linewidth=0.8,
    )
    for _, row in pseudo.iterrows():
        ax_c.annotate(row["tf"], (float(row["log2fc_lesion_vs_internal_control"]), float(row["regulon_mean_diff"])), xytext=(4, 4), textcoords="offset points", fontsize=6.2)
    ax_c.set_xlabel("Pseudobulk expression log2FC")
    ax_c.set_ylabel("Pseudobulk regulon diff")
    ax_c.set_title("Pseudobulk support for candidate TFs", pad=6)
    style_axes(ax_c)
    panel_label(ax_c, "C")

    for ax, thr, label in zip([ax_d1, ax_d2, ax_d3], [0.01, 0.05, 0.10], ["D", None, None]):
        mat = sens_rule[sens_rule["regulon_fdr_threshold"].eq(thr)].pivot_table(
            index="positive_fraction_threshold",
            columns="expr_abs_log2fc_threshold",
            values="n_retained",
            aggfunc="first",
        ).reindex(index=[0.10, 0.05, 0.00], columns=[0.00, 0.25, 0.50])
        im = ax.imshow(mat.values, aspect="auto", cmap="YlOrBr", vmin=float(sens_rule["n_retained"].min()), vmax=float(sens_rule["n_retained"].max()))
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_xticklabels([f"{c:.2f}" for c in mat.columns])
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_yticklabels([f"{i:.2f}" for i in mat.index] if ax is ax_d1 else [])
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                ax.text(j, i, str(int(mat.iloc[i, j])), ha="center", va="center", fontsize=6.1)
        ax.set_title(f"regulon FDR <= {thr:.2f}", pad=2)
        ax.set_xlabel("expr |log2FC|")
        if ax is ax_d1:
            ax.set_ylabel("positive fraction")
            panel_label(ax, "D")
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("Fig. S4. Sample-level robustness full view", y=1.01)
    return save_figure(fig, "FigS4_sample_level_robustness_full")


def build_figS5() -> list[Path]:
    donor = read_csv("external_validation/single_group_support/donor_level_summary.csv")
    ranking140 = read_csv("external_validation/single_group_support/external_support_ranking.csv").set_index("tf").loc[TF_ORDER_MAIN].reset_index()
    round2 = read_csv("external_validation_round2/round2_expression_and_regulon_master_table.csv").set_index("tf").loc[TF_ORDER_MAIN].reset_index()
    cluster = read_csv("external_validation_round2/round2_regulon_cluster_stats.csv")

    donor_sample = donor[(donor["level"].eq("sample")) & (donor["tf"].isin(TF_ORDER_MAIN))].copy()
    sample_order = donor_sample[["level_id", "n_cells"]].drop_duplicates().sort_values("n_cells", ascending=False)["level_id"].tolist()
    heat = donor_sample.pivot_table(index="tf", columns="level_id", values="positive_fraction", aggfunc="first").reindex(index=TF_ORDER_MAIN, columns=sample_order)

    cluster["signed_support"] = np.where(
        cluster["expected_higher_group_main"].astype(str).eq("lesion"),
        as_num(cluster["mean_auc_diff_lesion_minus_internal_control"]),
        -as_num(cluster["mean_auc_diff_lesion_minus_internal_control"]),
    )
    cluster_piv = cluster[cluster["tf"].isin(TF_ORDER_MAIN)].pivot_table(index="tf", columns="cluster", values="signed_support", aggfunc="first")
    cluster_cols = sorted(cluster_piv.columns, key=lambda x: int(str(x).split("_")[1]))
    cluster_piv = cluster_piv.reindex(index=TF_ORDER_MAIN, columns=cluster_cols)

    fig = plt.figure(figsize=(13.2, 8.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[0.95, 1.15], height_ratios=[0.92, 1.08], wspace=0.08)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[0, 1])
    ax_d = fig.add_subplot(gs[1, 1])

    im = ax_a.imshow(heat.values, aspect="auto", cmap="Reds", vmin=0, vmax=float(np.nanmax(heat.values)))
    ax_a.set_xticks(np.arange(len(sample_order)))
    ax_a.set_xticklabels(sample_order, rotation=30, ha="right")
    ax_a.set_yticks(np.arange(len(TF_ORDER_MAIN)))
    ax_a.set_yticklabels(TF_ORDER_MAIN)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            ax_a.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6.0, color="white" if val > 0.18 else "#222222")
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    ax_a.set_title("GSE140393 sample-level expression support", pad=6)
    fig.colorbar(im, ax=ax_a, fraction=0.05, pad=0.03, label="positive fraction")
    panel_label(ax_a, "A")

    y = np.arange(len(TF_ORDER_MAIN))[::-1]
    vals = ranking140.set_index("tf").loc[TF_ORDER_MAIN, "support_score"]
    ax_b.hlines(y, 0, vals.values, color="#D2D2D2", linewidth=1.0)
    ax_b.scatter(vals.values, y, s=42, c=[TF_COLORS[tf] for tf in TF_ORDER_MAIN], edgecolor="white", linewidth=0.7)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(TF_ORDER_MAIN)
    ax_b.set_xlim(0, 1.05)
    ax_b.set_xlabel("Support score")
    ax_b.set_title("GSE140393 ranked supportive summary", pad=6)
    style_axes(ax_b, grid_axis="x")
    for yi, tf in zip(y, TF_ORDER_MAIN):
        expected = ranking140.set_index("tf").loc[tf, "expected_higher_group_main"]
        ax_b.text(float(vals.loc[tf]) + 0.02, yi, str(expected), fontsize=6.0, va="center", color="#555555")
    panel_label(ax_b, "B")

    ax_c.axvline(0, color="#A0A0A0", linewidth=0.8)
    ax_c.axhline(0.4, color="#D6D6D6", linewidth=0.8, linestyle="--")
    ax_c.scatter(
        round2["expr_support_score"],
        round2["regulon_support_score"],
        s=120 + 450 * round2["combined_support_score"],
        c=[TF_COLORS[tf] for tf in round2["tf"]],
        alpha=0.84,
        edgecolor="white",
        linewidth=0.8,
    )
    for _, row in round2.iterrows():
        ax_c.annotate(row["tf"], (float(row["expr_support_score"]), float(row["regulon_support_score"])), xytext=(4, 4), textcoords="offset points", fontsize=6.2)
    ax_c.set_xlabel("Expression support score")
    ax_c.set_ylabel("Projected regulon support score")
    ax_c.set_title("GSE190452 expression-regulon support summary", pad=6)
    style_axes(ax_c)
    panel_label(ax_c, "C")

    im2 = ax_d.imshow(cluster_piv.values, aspect="auto", cmap="RdBu_r", vmin=-0.012, vmax=0.012)
    ax_d.set_xticks(np.arange(len(cluster_cols)))
    ax_d.set_xticklabels(cluster_cols, rotation=90)
    ax_d.set_yticks(np.arange(len(TF_ORDER_MAIN)))
    ax_d.set_yticklabels(TF_ORDER_MAIN)
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    ax_d.set_title("GSE190452 cluster-level projected regulon support", pad=6)
    fig.colorbar(im2, ax=ax_d, fraction=0.05, pad=0.03, label="directionally signed AUC diff")
    panel_label(ax_d, "D")

    fig.suptitle("Fig. S5. Supportive external evidence full view", y=1.01)
    return save_figure(fig, "FigS5_supportive_external_evidence_full")


def _top_enrichment_terms(df: pd.DataFrame, tf: str, top_go: int = 8, top_kegg: int = 4) -> pd.DataFrame:
    sub = df[df["tf"].astype(str).eq(tf)].copy()
    sub = sub[sub["gene_set_type"].astype(str).eq("regulon_targets")].copy()
    if sub.empty:
        sub = df[df["tf"].astype(str).eq(tf)].copy()
    sub["fdr_bh"] = as_num(sub["fdr_bh"])
    go = sub[sub["source"].astype(str).eq("go")].nsmallest(top_go, "fdr_bh")
    kegg = sub[sub["source"].astype(str).eq("kegg")].nsmallest(top_kegg, "fdr_bh")
    out = pd.concat([go, kegg], ignore_index=True)
    out["neg_log10_fdr"] = neg_log10(out["fdr_bh"])
    out["overlap_ratio"] = as_num(out["overlap_ratio"])
    out["overlap_count"] = as_num(out["overlap_count"])
    out["term_short"] = out["term"].astype(str).str.replace(r" \(GO:\d+\)$", "", regex=True)
    return out.sort_values("neg_log10_fdr", ascending=True)


def build_figS6() -> list[Path]:
    go = read_csv("functional_interpretation/03_go_enrichment_all.csv")
    kegg = read_csv("functional_interpretation/03_kegg_enrichment_all.csv")
    enrich = pd.concat([go, kegg], ignore_index=True, sort=False)
    conv = read_csv("functional_interpretation/04_program_convergence_table.csv").set_index("tf").loc[TF_ORDER_MAIN].reset_index()
    gene_summary = read_csv("functional_interpretation/02_tf_gene_sets_summary.csv").set_index("tf").loc[TF_ORDER_MAIN].reset_index()

    fig = plt.figure(figsize=(13.4, 10.2), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.72], wspace=0.12)
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[2, 0]),
        fig.add_subplot(gs[2, 1]),
    ]

    for ax, tf, label in zip(axes[:4], TF_ORDER_MAIN, ["A", "B", "C", "D"]):
        sub = _top_enrichment_terms(enrich, tf)
        scatter = ax.scatter(
            sub["neg_log10_fdr"],
            np.arange(len(sub)),
            s=18 + 11 * sub["overlap_count"],
            c=sub["overlap_ratio"],
            cmap="YlOrRd",
            edgecolor="white",
            linewidth=0.6,
        )
        ax.set_yticks(np.arange(len(sub)))
        ax.set_yticklabels(sub["term_short"])
        ax.set_xlabel("-log10(FDR)")
        ax.set_title(f"{tf} top GO/KEGG terms", pad=6)
        style_axes(ax, grid_axis="x")
        ax.text(0.98, 0.02, "GO + KEGG", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.0, color="#666666")
        panel_label(ax, label)
    cbar = fig.colorbar(scatter, ax=axes[:4], fraction=0.015, pad=0.02)
    cbar.set_label("overlap ratio")

    conv = conv.merge(gene_summary[["tf", "regulon_target_count"]], on="tf", how="left")
    axes[4].scatter(
        conv["overlap_ratio"],
        conv["overlap_count"],
        s=80 + 0.45 * conv["regulon_target_count"],
        c=[TF_COLORS[tf] for tf in conv["tf"]],
        edgecolor="white",
        linewidth=0.8,
        alpha=0.86,
    )
    for _, row in conv.iterrows():
        axes[4].annotate(
            f"{row['tf']} ({row['program_theme_1']})",
            (float(row["overlap_ratio"]), float(row["overlap_count"])),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=6.0,
        )
    axes[4].set_xlabel("Regulon-DEG overlap ratio")
    axes[4].set_ylabel("Overlap genes")
    axes[4].set_title("Program convergence summary", pad=6)
    style_axes(axes[4])
    panel_label(axes[4], "E")

    count_rows = []
    for tf in TF_ORDER_MAIN:
        sub = enrich[(enrich["tf"].astype(str).eq(tf)) & (enrich["gene_set_type"].astype(str).eq("regulon_targets"))].copy()
        if sub.empty:
            sub = enrich[enrich["tf"].astype(str).eq(tf)].copy()
        sub["fdr_bh"] = as_num(sub["fdr_bh"])
        count_rows.append(
            {
                "tf": tf,
                "GO_sig_terms": int(((sub["source"].astype(str) == "go") & (sub["fdr_bh"] <= 0.05)).sum()),
                "KEGG_sig_terms": int(((sub["source"].astype(str) == "kegg") & (sub["fdr_bh"] <= 0.05)).sum()),
            }
        )
    count_df = pd.DataFrame(count_rows).set_index("tf").loc[TF_ORDER_MAIN]
    im = axes[5].imshow(count_df.values, aspect="auto", cmap="Blues", vmin=0, vmax=float(count_df.values.max()))
    axes[5].set_xticks(np.arange(count_df.shape[1]))
    axes[5].set_xticklabels(count_df.columns, rotation=25, ha="right")
    axes[5].set_yticks(np.arange(count_df.shape[0]))
    axes[5].set_yticklabels(count_df.index)
    for i in range(count_df.shape[0]):
        for j in range(count_df.shape[1]):
            val = int(count_df.iloc[i, j])
            axes[5].text(j, i, str(val), ha="center", va="center", fontsize=6.2, color="white" if val > 8 else "#222222")
    for spine in axes[5].spines.values():
        spine.set_visible(False)
    axes[5].set_title("Significant term counts", pad=6)
    fig.colorbar(im, ax=axes[5], fraction=0.05, pad=0.03, label="count")
    panel_label(axes[5], "F")

    fig.suptitle("Fig. S6. Functional interpretation and program convergence full view", y=1.01)
    return save_figure(fig, "FigS6_functional_interpretation_and_program_convergence_full")


def build_figS7() -> list[Path]:
    paper = read_csv("drug_repositioning/16_paper_ready_table_after_bbb_v21.csv")
    full = read_csv("drug_repositioning/13_integrated_after_bbb_v21.csv")
    reviewed = full[full["after_bbb_final_layer"].notna()].copy()
    theme_col = "mechanism_theme" if "mechanism_theme" in reviewed.columns else "indicative_mechanism"

    layer_order = [
        "headline_cns_mechanism_direction_leads",
        "supportive_after_bbb_leads",
        "peripheral_program_modulating_clues",
    ]
    layer_counts = reviewed.groupby(["after_bbb_final_layer", "associated_axis"]).size().unstack(fill_value=0).reindex(layer_order).fillna(0)
    top_paper = paper.sort_values("final_score_after_bbb", ascending=True).copy()

    top_themes = (
        reviewed[theme_col]
        .fillna("unassigned")
        .astype(str)
        .map(clean_mechanism_theme)
        .value_counts()
        .head(6)
        .index.tolist()
    )
    theme_mat = (
        reviewed.assign(mechanism_theme=reviewed[theme_col].map(clean_mechanism_theme))
        .query("mechanism_theme in @top_themes")
        .groupby(["mechanism_theme", "after_bbb_final_layer"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=top_themes, columns=layer_order)
        .fillna(0)
    )

    fig = plt.figure(figsize=(12.8, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[0.95, 1.1], height_ratios=[0.88, 1.02], wspace=0.08)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[0, 1])
    ax_d = fig.add_subplot(gs[1, 1])

    left = np.zeros(len(layer_counts.index))
    for axis in ["NFE2L2", "THRB", "BHLHE40", "SOX2"]:
        vals = layer_counts[axis].values if axis in layer_counts.columns else np.zeros(len(layer_counts.index))
        ax_a.barh(
            [short_layer_name(x) for x in layer_counts.index],
            vals,
            left=left,
            color=TF_COLORS.get(axis, "#777777"),
            alpha=0.85,
            label=axis,
        )
        left += vals
    ax_a.set_xlabel("Compounds")
    ax_a.set_title("After-BBB layer counts by axis", pad=6)
    style_axes(ax_a, grid_axis="x")
    ax_a.legend(frameon=False, loc="lower right", ncol=2)
    panel_label(ax_a, "A")

    ax_b.hlines(np.arange(len(top_paper)), xmin=0, xmax=top_paper["final_score_after_bbb"], color="#D0D0D0", linewidth=1.0)
    ax_b.scatter(
        top_paper["final_score_after_bbb"],
        np.arange(len(top_paper)),
        s=32,
        c=[DRUG_LAYER_COLORS.get(layer, "#999999") for layer in top_paper["after_bbb_final_layer"]],
        edgecolor="white",
        linewidth=0.7,
    )
    ax_b.set_yticks(np.arange(len(top_paper)))
    ax_b.set_yticklabels(top_paper["compound"])
    ax_b.set_xlabel("Final score after BBB review")
    ax_b.set_title("Paper-ready after-BBB ranked compounds", pad=6)
    style_axes(ax_b, grid_axis="x")
    panel_label(ax_b, "B")

    ax_c.scatter(
        reviewed["final_score_prebbb_v21"],
        reviewed["final_score_after_bbb"],
        s=42,
        c=[DRUG_LAYER_COLORS.get(layer, "#999999") for layer in reviewed["after_bbb_final_layer"]],
        alpha=0.75,
        edgecolor="white",
        linewidth=0.6,
    )
    ax_c.plot([0, 1], [0, 1], color="#CCCCCC", linewidth=0.8, linestyle="--")
    top_labels = reviewed.sort_values("final_score_after_bbb", ascending=False).head(8)
    for _, row in top_labels.iterrows():
        ax_c.annotate(row["compound"], (float(row["final_score_prebbb_v21"]), float(row["final_score_after_bbb"])), xytext=(4, 4), textcoords="offset points", fontsize=6.0)
    ax_c.set_xlabel("Pre-BBB final score v2.1")
    ax_c.set_ylabel("After-BBB final score")
    ax_c.set_title("Pre-BBB vs after-BBB score shift", pad=6)
    style_axes(ax_c)
    panel_label(ax_c, "C")

    im = ax_d.imshow(theme_mat.values, aspect="auto", cmap="Greens", vmin=0, vmax=float(theme_mat.values.max()))
    ax_d.set_xticks(np.arange(theme_mat.shape[1]))
    ax_d.set_xticklabels([short_layer_name(c) for c in theme_mat.columns], rotation=25, ha="right")
    ax_d.set_yticks(np.arange(theme_mat.shape[0]))
    ax_d.set_yticklabels(theme_mat.index)
    for i in range(theme_mat.shape[0]):
        for j in range(theme_mat.shape[1]):
            ax_d.text(j, i, str(int(theme_mat.iloc[i, j])), ha="center", va="center", fontsize=6.0, color="white" if theme_mat.iloc[i, j] >= 2 else "#222222")
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    ax_d.set_title("Mechanism-theme distribution across after-BBB layers", pad=6)
    fig.colorbar(im, ax=ax_d, fraction=0.05, pad=0.03, label="count")
    panel_label(ax_d, "D")

    fig.suptitle("Fig. S7. Drug repositioning v2.1 after-BBB summary", y=1.01)
    return save_figure(fig, "FigS7_drug_repositioning_v21_after_bbb_summary")


def write_index(table1: Path, figures: dict[str, list[Path]], supp_tables: list[Path]) -> None:
    figure_names = [paths[0].name for _, paths in figures.items()]
    table_names = [table1.name] + [p.name for p in supp_tables]
    text = f"""
    # Supplementary Package Index

    This supplementary package was generated from the current local project state without rerunning upstream analyses.

    ## Mainline boundaries

    - Current manuscript mainline: GSE268807 astrocyte pilot + pySCENIC + CellOracle + sample-level robustness + functional interpretation/program convergence + supportive external analysis.
    - GSE140393 is treated as **single-group supportive analysis**.
    - GSE190452 is treated as **cross-syndrome supportive analysis**.
    - Neither GSE140393 nor GSE190452 is described here as formal same-disease external validation.
    - No SCENIC+ current mainline is used.
    - No wet-lab validation is included in the current package.
    - Drug repositioning remains **exploratory only** and is kept in supplementary material.

    ## Main-text table added

    - `{table1.relative_to(ROOT).as_posix()}`: Table 1 summarizing the discovery cohort, supportive external datasets, the not-included/private GSE275302 route, and the archived SCENIC+/Kaggle exploration route.

    ## Supplementary figures

    - `{figure_names[0]}`: Discovery object QC, sample composition, upstream QC/subsampling trace, and current object lineage for the GSE268807 astrocyte pilot.
      - Main-text linkage: extends manuscript Fig. 1 / discovery object context.
      - Data sources: `analysis_outputs/gse268807_astrocyte_pilot/qc_metrics/`, `celloracle_run/prepared_data/`, `final_exports/merge_report.txt`, `robustness_validation/01_input_object_summary.txt`.

    - `{figure_names[1]}`: Extended differential regulon landscape with a larger heatmap and broader lesion/internal_control-associated regulon ranking view.
      - Main-text linkage: extends Fig. 2 and Table 2.
      - Data sources: `analysis_outputs/group_compare/`, `final_exports/auc_mtx_matched_to_h5ad.csv`, `celloracle_run/prepared_data/celloracle_round1_metadata.csv`.

    - `{figure_names[2]}`: Extended CellOracle prioritization with detailed metric heatmap and representative perturbation panels for all round1 TFs.
      - Main-text linkage: extends Fig. 3 and Table 3.
      - Data sources: `celloracle_run/ko_round1/`.

    - `{figure_names[3]}`: Expanded **sample-level robustness** view, including leave-one-sample-out expression/regulon effects, pseudobulk support, and shortlist sensitivity.
      - Main-text linkage: extends Fig. 4 robustness panels and Table 3.
      - Data sources: `robustness_validation/`.

    - `{figure_names[4]}`: Full supportive external evidence package, integrating GSE140393 single-group supportive analysis and GSE190452 cross-syndrome supportive analysis.
      - Main-text linkage: extends Fig. 4 supportive external evidence panels and Table 3.
      - Data sources: `external_validation/`, `external_validation_round2/`.

    - `{figure_names[5]}`: Real GO/KEGG enrichment dot plots and program convergence summaries generated from current functional interpretation result tables.
      - Main-text linkage: supports the mechanistic interpretation layer that feeds into the Fig. 5 conceptual synthesis.
      - Data sources: `functional_interpretation/`.

    - `{figure_names[6]}`: Exploratory after-BBB drug repositioning v2.1 summary, kept only as supplementary mechanism-direction clue material.
      - Main-text linkage: corresponds to the exploratory drug clue module and should remain supplementary.
      - Data sources: `drug_repositioning/13_integrated_after_bbb_v21.csv`, `drug_repositioning/16_paper_ready_table_after_bbb_v21.csv`.

    ## Supplementary tables

    - `{supp_tables[0].relative_to(ROOT).as_posix()}`: Full differential regulon statistics; extends Fig. 2.
    - `{supp_tables[1].relative_to(ROOT).as_posix()}`: Extended candidate TF shortlist / secondary candidate metrics; extends Table 2.
    - `{supp_tables[2].relative_to(ROOT).as_posix()}`: Full CellOracle round1 metrics; extends Fig. 3 and Table 3.
    - `{supp_tables[3].relative_to(ROOT).as_posix()}`: Detailed **sample-level robustness** package, including leave-one-sample-out, pseudobulk, shortlist sensitivity, and integrated summaries.
    - `{supp_tables[4].relative_to(ROOT).as_posix()}`: Full GO/KEGG/program convergence package.
    - `{supp_tables[5].relative_to(ROOT).as_posix()}`: Full supportive external analysis package for GSE140393 and GSE190452.
    - `{supp_tables[6].relative_to(ROOT).as_posix()}`: Full after-BBB drug repositioning v2.1 table from the current latest route.

    ## Current main-text anchors already present

    - `manuscript_output/figures_main/Fig2/Fig2_main_v2.png`
    - `manuscript_output/figures_main/Fig3/Fig3_main_v2.png`
    - `manuscript_output/figures_main/Fig4/Fig4_main_v2.png`
    - `manuscript_output/tables_main/Table1_cohort_and_analytical_objects.csv`
    - `manuscript_output/tables_main/Table2_shortlist_main.csv`
    - `manuscript_output/tables_main/Table3_integrated_priority_main.csv`
    """
    write_text(INDEX_PATH, text)


def main() -> None:
    setup_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    MAIN_TABLE_DIR.mkdir(parents=True, exist_ok=True)

    table1 = build_table1()
    supp_tables = build_table_supplementary()
    figures = {
        "FigS1": build_figS1(),
        "FigS2": build_figS2(),
        "FigS3": build_figS3(),
        "FigS4": build_figS4(),
        "FigS5": build_figS5(),
        "FigS6": build_figS6(),
        "FigS7": build_figS7(),
    }
    write_index(table1, figures, supp_tables)

    print(f"Table 1 written to: {table1}")
    for p in supp_tables:
        print(f"Supplementary table written to: {p}")
    for name, paths in figures.items():
        for p in paths:
            print(f"{name} written to: {p}")
    print(f"Supplementary index written to: {INDEX_PATH}")


if __name__ == "__main__":
    main()
