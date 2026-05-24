#!/usr/bin/env python3
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[2]
FIG_ROOT = ROOT / "manuscript_output" / "figures_main"
TABLE_ROOT = ROOT / "manuscript_output" / "tables_main"

MAIN_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
SHORTLIST_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2", "SATB2", "RARB", "HMGA1"]
GROUP_ORDER = ["internal_control", "lesion"]
REGULONS_MAIN = [f"{tf}(+)" for tf in MAIN_TFS]

COLORS = {
    "lesion": "#C45A4D",
    "internal_control": "#4D78A8",
    "neutral": "#666666",
    "grid": "#E8E8E8",
    "NFE2L2": "#C45A4D",
    "THRB": "#4D78A8",
    "BHLHE40": "#C98A2E",
    "SOX2": "#8A7FA8",
    "SATB2": "#6F8C63",
    "RARB": "#5E8D87",
    "HMGA1": "#9A6F52",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.1,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.2,
            "figure.titlesize": 10.2,
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


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.png", dpi=360, bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def style_axes(ax: plt.Axes, grid_axis: str | None = "y") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=COLORS["grid"], linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.13,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=10.4,
        fontweight="bold",
        ha="left",
        va="top",
    )


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


def trim_white(img: np.ndarray, threshold: float = 0.985) -> np.ndarray:
    arr = np.asarray(img)
    if arr.ndim == 2:
        mask = arr < threshold
    else:
        rgb = arr[..., :3]
        mask = np.any(rgb < threshold, axis=2)
    coords = np.argwhere(mask)
    if coords.size == 0:
        return arr
    y0, x0 = coords.min(axis=0)[:2]
    y1, x1 = coords.max(axis=0)[:2] + 1
    return arr[y0:y1, x0:x1]


def load_auc_with_metadata() -> pd.DataFrame:
    auc = read_csv("final_exports/auc_mtx_matched_to_h5ad.csv")
    auc = auc.rename(columns={auc.columns[0]: "cell_id"})
    meta = read_csv("celloracle_run/prepared_data/celloracle_round1_metadata.csv")
    first = meta.columns[0]
    if first.startswith("Unnamed"):
        meta = meta.rename(columns={first: "cell_id"})
    if "cell_id" not in meta.columns and "join_key" in meta.columns:
        meta["cell_id"] = meta["join_key"]
    keep = [
        "cell_id",
        "group",
        "sample_id",
        "donor_id",
        "clinical_sample_label",
        "sample_prefix",
        "subtype",
    ]
    keep = [c for c in keep if c in meta.columns]
    return auc.merge(meta[keep], on="cell_id", how="left")


def unique_keep_order(items: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def build_table2_shortlist(shortlist: pd.DataFrame) -> Path:
    out = TABLE_ROOT / "Table2_shortlist_main.csv"
    tbl = shortlist.copy()
    tbl = tbl[tbl["tf"].isin(SHORTLIST_TFS)].copy()
    order_map = {tf: i for i, tf in enumerate(["NFE2L2", "THRB", "BHLHE40", "SOX2", "SATB2", "RARB", "HMGA1"])}
    tbl["dominant_axis"] = tbl["regulon_effect_group"].map(
        {"lesion": "lesion-associated", "internal_control": "internal_control-associated"}
    )
    tbl["regulon_fdr"] = as_num(tbl["regulon_fdr"])
    tbl["expr_fdr"] = as_num(tbl["expr_wilcoxon_fdr"])
    tbl["expr_log2fc_lesion_vs_internal_control"] = as_num(tbl["log2fc_lesion_vs_internal_control"])
    tbl["regulon_auc_diff_lesion_vs_internal_control"] = as_num(
        tbl["regulon_mean_diff_lesion_minus_internal_control"]
    )
    tbl["positive_fraction_in_higher_group"] = as_num(tbl["positive_fraction_in_higher_group"])
    tbl = tbl[
        [
            "tf",
            "regulon",
            "dominant_axis",
            "expr_log2fc_lesion_vs_internal_control",
            "regulon_auc_diff_lesion_vs_internal_control",
            "regulon_fdr",
            "expr_fdr",
            "positive_fraction_in_higher_group",
            "shortlist_tier",
        ]
    ].rename(
        columns={
            "tf": "TF",
            "regulon": "Regulon",
            "dominant_axis": "Dominant_axis",
            "regulon_fdr": "Regulon_FDR",
            "expr_fdr": "Expression_FDR",
            "shortlist_tier": "Shortlist_tier",
        }
    )
    tbl = tbl.sort_values("TF", key=lambda s: s.map(order_map))
    tbl.to_csv(out, index=False)
    return out


def build_table3_integrated(
    shortlist: pd.DataFrame,
    ranking: pd.DataFrame,
    robust: pd.DataFrame,
    gse140: pd.DataFrame,
    gse190: pd.DataFrame,
) -> Path:
    out = TABLE_ROOT / "Table3_integrated_priority_main.csv"
    order_map = {tf: i for i, tf in enumerate(["NFE2L2", "THRB", "BHLHE40", "SOX2"])}
    base = shortlist[shortlist["tf"].isin(MAIN_TFS)][
        [
            "tf",
            "regulon_effect_group",
            "log2fc_lesion_vs_internal_control",
            "regulon_mean_diff_lesion_minus_internal_control",
            "regulon_fdr",
        ]
    ].copy()
    base = base.merge(
        ranking[
            [
                "tf",
                "mean_shift_length",
                "overall_shift_rank",
                "recovery_index",
                "recovery_index_rank",
            ]
        ],
        on="tf",
        how="left",
    )
    base = base.merge(
        robust[
            [
                "tf",
                "loo_expr_consistency_rate",
                "loo_regulon_consistency_rate",
                "shortlist_retention_frequency",
                "final_internal_robustness_level",
                "recommended_priority_tier",
            ]
        ],
        on="tf",
        how="left",
    )
    base = base.merge(
        gse140[["tf", "support_score"]].rename(columns={"support_score": "GSE140393_support_score"}),
        on="tf",
        how="left",
    )
    base = base.merge(
        gse190[
            [
                "tf",
                "group_log2fc_lesion_vs_internal_control",
                "cluster_support_fraction",
                "support_score",
            ]
        ].rename(
            columns={
                "group_log2fc_lesion_vs_internal_control": "GSE190452_group_log2FC",
                "cluster_support_fraction": "GSE190452_cluster_support_fraction",
                "support_score": "GSE190452_support_score",
            }
        ),
        on="tf",
        how="left",
    )
    base["dominant_axis"] = base["regulon_effect_group"].map(
        {"lesion": "lesion-associated", "internal_control": "internal_control-associated"}
    )
    base = base[
        [
            "tf",
            "dominant_axis",
            "log2fc_lesion_vs_internal_control",
            "regulon_mean_diff_lesion_minus_internal_control",
            "mean_shift_length",
            "overall_shift_rank",
            "recovery_index",
            "recovery_index_rank",
            "loo_expr_consistency_rate",
            "loo_regulon_consistency_rate",
            "shortlist_retention_frequency",
            "GSE140393_support_score",
            "GSE190452_group_log2FC",
            "GSE190452_cluster_support_fraction",
            "GSE190452_support_score",
            "final_internal_robustness_level",
            "recommended_priority_tier",
        ]
    ].rename(
        columns={
            "tf": "TF",
            "log2fc_lesion_vs_internal_control": "Discovery_expression_log2FC",
            "regulon_mean_diff_lesion_minus_internal_control": "Discovery_regulon_AUC_diff",
            "mean_shift_length": "CellOracle_mean_shift_length",
            "overall_shift_rank": "KO_rank_by_shift",
            "recovery_index": "Recovery_index",
            "recovery_index_rank": "Recovery_index_rank",
            "loo_expr_consistency_rate": "LOSO_expression_consistency",
            "loo_regulon_consistency_rate": "LOSO_regulon_consistency",
            "shortlist_retention_frequency": "Shortlist_retention_frequency",
            "final_internal_robustness_level": "Internal_robustness_level",
            "recommended_priority_tier": "Final_priority_tier",
            "dominant_axis": "Dominant_axis",
        }
    )
    base = base.sort_values("TF", key=lambda s: s.map(order_map))
    base.to_csv(out, index=False)
    return out


def add_sample_group_strip(ax: plt.Axes, sample_order: list[str], sample_groups: list[str]) -> None:
    strip = ax.inset_axes([0.0, 1.01, 1.0, 0.08])
    strip.set_xlim(0, len(sample_order))
    strip.set_ylim(0, 1)
    for i, group in enumerate(sample_groups):
        strip.add_patch(
            Rectangle(
                (i, 0),
                1,
                1,
                facecolor=COLORS[group],
                alpha=0.9,
                edgecolor="white",
                linewidth=0.6,
            )
        )
    strip.text(1.0, 0.5, "internal_control", color="white", ha="center", va="center", fontsize=6.0, fontweight="bold")
    strip.text(3.0, 0.5, "lesion", color="white", ha="center", va="center", fontsize=6.0, fontweight="bold")
    strip.axis("off")


def build_fig2() -> dict[str, list[str]]:
    out_dir = FIG_ROOT / "Fig2"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = [
        "analysis_outputs/group_compare/regulon_group_statistics.csv",
        "analysis_outputs/group_compare/top20_differential_regulons.csv",
        "analysis_outputs/celloracle_candidates/celloracle_tf_shortlist.csv",
        "final_exports/auc_mtx_matched_to_h5ad.csv",
        "celloracle_run/prepared_data/celloracle_round1_metadata.csv",
    ]

    stats = read_csv(inputs[0])
    stats["lesion_minus_internal_control"] = as_num(stats["mean_lesion"]) - as_num(stats["mean_internal_control"])
    stats["neg_log10_fdr"] = neg_log10(stats["fdr_bh"])
    top20 = read_csv(inputs[1])["regulon"].astype(str).tolist()
    shortlist = read_csv(inputs[2])
    shortlist["neg_log10_regulon_fdr"] = neg_log10(shortlist["regulon_fdr"])
    shortlist["log2fc"] = as_num(shortlist["log2fc_lesion_vs_internal_control"])
    shortlist["posfrac"] = as_num(shortlist["positive_fraction_in_higher_group"])
    auc_meta = load_auc_with_metadata()

    ordered_regs = unique_keep_order(top20 + stats.sort_values(["fdr_bh", "abs_mean_diff"], ascending=[True, False])["regulon"].astype(str).tolist())
    selected_regs = ordered_regs[:14]
    for reg in REGULONS_MAIN:
        if reg not in selected_regs:
            for i in range(len(selected_regs) - 1, -1, -1):
                if selected_regs[i] not in REGULONS_MAIN:
                    selected_regs[i] = reg
                    break
    selected_regs = unique_keep_order(selected_regs)
    if len(selected_regs) < 14:
        for reg in ordered_regs:
            if reg not in selected_regs:
                selected_regs.append(reg)
            if len(selected_regs) == 14:
                break
    sel_stats = stats.set_index("regulon").loc[selected_regs].copy()
    sel_stats = sel_stats.sort_values("lesion_minus_internal_control", ascending=True)
    selected_regs = sel_stats.index.tolist()

    sample_meta = (
        auc_meta[["clinical_sample_label", "group", "donor_id"]]
        .drop_duplicates()
        .sort_values(["group", "donor_id", "clinical_sample_label"])
    )
    sample_order = sample_meta["clinical_sample_label"].tolist()
    sample_groups = sample_meta["group"].tolist()
    sample_means = auc_meta.groupby("clinical_sample_label")[selected_regs].mean().loc[sample_order]
    heat = sample_means.T
    heat = heat.sub(heat.mean(axis=1), axis=0)
    heat = heat.div(heat.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0.0)

    fig = plt.figure(figsize=(12.2, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 0.9, 1.08], height_ratios=[0.92, 1.08], wspace=0.06)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[:, 2])

    im = ax_a.imshow(heat.values, aspect="auto", cmap="RdBu_r", vmin=-1.6, vmax=1.6)
    ax_a.set_xticks(np.arange(len(sample_order)))
    ax_a.set_xticklabels(sample_order, rotation=28, ha="right")
    ax_a.set_yticks(np.arange(len(selected_regs)))
    ax_a.set_yticklabels(selected_regs)
    ax_a.set_title("Differential regulon heatmap", pad=6)
    ax_a.tick_params(length=0)
    ax_a.axvline(1.5, color="#D6D6D6", linewidth=0.9)
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    add_sample_group_strip(ax_a, sample_order, sample_groups)
    cbar = fig.colorbar(im, ax=ax_a, fraction=0.048, pad=0.03)
    cbar.set_label("row z-score")
    cbar.ax.tick_params(labelsize=6.0)
    panel_label(ax_a, "A")

    ax_b.scatter(
        stats["lesion_minus_internal_control"],
        stats["neg_log10_fdr"],
        s=15,
        c=np.where(stats["lesion_minus_internal_control"] >= 0, COLORS["lesion"], COLORS["internal_control"]),
        alpha=0.62,
        edgecolors="none",
        zorder=2,
    )
    ax_b.axvline(0, color="#9A9A9A", linewidth=0.8)
    ax_b.set_xlabel("Mean AUC difference")
    ax_b.set_ylabel("-log10(FDR)")
    ax_b.set_title("Differential regulon summary", pad=6)
    style_axes(ax_b)
    ax_b.text(0.02, 0.98, "internal_control-associated", transform=ax_b.transAxes, ha="left", va="top", fontsize=6.1, color=COLORS["internal_control"])
    ax_b.text(0.98, 0.98, "lesion-associated", transform=ax_b.transAxes, ha="right", va="top", fontsize=6.1, color=COLORS["lesion"])
    label_offsets = {
        "THRB(+)": (-4, 6, "right"),
        "SATB2(+)": (-4, 6, "right"),
        "RARB(+)": (-4, 5, "right"),
        "SOX2(+)": (4, 4, "left"),
        "NFE2L2(+)": (6, 4, "left"),
        "BHLHE40(+)": (6, 3, "left"),
        "HMGA1(+)": (6, -2, "left"),
    }
    for reg, (dx, dy, ha) in label_offsets.items():
        row = stats[stats["regulon"].eq(reg)]
        if row.empty:
            continue
        xx = float(row["lesion_minus_internal_control"].iloc[0])
        yy = float(row["neg_log10_fdr"].iloc[0])
        tf = reg.replace("(+)", "")
        ax_b.scatter([xx], [yy], s=46, color="white", edgecolor=COLORS.get(tf, "#333333"), linewidth=1.1, zorder=4)
        ax_b.annotate(
            reg,
            (xx, yy),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=6.3,
            bbox={"boxstyle": "round,pad=0.14", "fc": "white", "ec": "none", "alpha": 0.9},
            zorder=5,
        )
    panel_label(ax_b, "B")

    pos_map = {"NFE2L2(+)": 0.0, "THRB(+)": 1.15, "BHLHE40(+)": 2.75, "SOX2(+)": 3.65}
    off_map = {"internal_control": -0.14, "lesion": 0.14}
    width_map = {"NFE2L2(+)": 0.28, "THRB(+)": 0.28, "BHLHE40(+)": 0.18, "SOX2(+)": 0.18}
    for reg in REGULONS_MAIN:
        for group in GROUP_ORDER:
            vals = pd.to_numeric(auc_meta.loc[auc_meta["group"].eq(group), reg], errors="coerce").dropna().values
            vp = ax_c.violinplot(
                [vals],
                positions=[pos_map[reg] + off_map[group]],
                widths=width_map[reg],
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )
            for body in vp["bodies"]:
                body.set_facecolor(COLORS[group])
                body.set_edgecolor(COLORS[group])
                body.set_alpha(0.58 if reg in {"NFE2L2(+)", "THRB(+)"} else 0.32)
                body.set_linewidth(0.5)
            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            xpos = pos_map[reg] + off_map[group]
            ax_c.plot([xpos, xpos], [q1, q3], color="#2A2A2A", linewidth=0.8, zorder=3)
            ax_c.scatter([xpos], [med], s=13, color="white", edgecolor="#2A2A2A", linewidth=0.6, zorder=4)
    ax_c.axvline(1.95, color="#DDDDDD", linewidth=0.8)
    ax_c.set_xlim(-0.45, 4.15)
    ax_c.set_xticks([pos_map[r] for r in REGULONS_MAIN])
    ax_c.set_xticklabels(REGULONS_MAIN, rotation=28, ha="right")
    for label in ax_c.get_xticklabels():
        txt = label.get_text()
        if txt in {"NFE2L2(+)", "THRB(+)"}:
            label.set_fontweight("bold")
    ax_c.set_ylabel("Single-cell regulon AUC")
    ax_c.set_title("Key regulon activity comparison", pad=6)
    style_axes(ax_c)
    ax_c.legend(
        handles=[
            Line2D([0], [0], color=COLORS["internal_control"], linewidth=5, alpha=0.65, label="internal_control"),
            Line2D([0], [0], color=COLORS["lesion"], linewidth=5, alpha=0.65, label="lesion"),
        ],
        frameon=False,
        loc="upper right",
        handlelength=1.6,
    )
    panel_label(ax_c, "C")

    bubble = shortlist[shortlist["tf"].isin(SHORTLIST_TFS)].copy()
    bubble = bubble.set_index("tf").loc[SHORTLIST_TFS].reset_index()
    bubble_sizes = 95 + 320 * bubble["posfrac"].fillna(0)
    ax_d.scatter(
        bubble["log2fc"],
        bubble["neg_log10_regulon_fdr"],
        s=bubble_sizes,
        c=[COLORS.get(tf, COLORS["neutral"]) for tf in bubble["tf"]],
        alpha=0.86,
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
    )
    ax_d.axvline(0, color="#9A9A9A", linewidth=0.8)
    ax_d.set_xlabel("TF expression log2FC")
    ax_d.set_ylabel("-log10(regulon FDR)")
    ax_d.set_title("Candidate TF prioritization", pad=6)
    style_axes(ax_d)
    ax_d.text(0.02, 0.98, "internal_control-associated", transform=ax_d.transAxes, ha="left", va="top", fontsize=6.1, color=COLORS["internal_control"])
    ax_d.text(0.98, 0.98, "lesion-associated", transform=ax_d.transAxes, ha="right", va="top", fontsize=6.1, color=COLORS["lesion"])
    offsets = {
        "SATB2": (-12, 5, "left"),
        "THRB": (-4, 5, "right"),
        "RARB": (-12, 4, "left"),
        "BHLHE40": (5, 4, "left"),
        "NFE2L2": (6, 4, "left"),
        "SOX2": (4, 6, "left"),
        "HMGA1": (6, 4, "left"),
    }
    for _, row in bubble.iterrows():
        dx, dy, ha = offsets.get(row["tf"], (4, 3, "left"))
        ax_d.annotate(
            row["tf"],
            (float(row["log2fc"]), float(row["neg_log10_regulon_fdr"])),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=6.6,
            ha=ha,
            va="center",
            fontweight="bold" if row["tf"] in {"NFE2L2", "THRB"} else "normal",
        )
    legend_sizes = [0.2, 0.4, 0.6]
    size_handles = [
        plt.scatter([], [], s=95 + 320 * val, color="#BDBDBD", edgecolor="white", linewidth=0.6, alpha=0.75, label=f"{val:.1f}")
        for val in legend_sizes
    ]
    ax_d.legend(
        handles=size_handles,
        title="Positive fraction",
        frameon=False,
        loc="lower right",
        borderpad=0.2,
        labelspacing=0.8,
    )
    panel_label(ax_d, "D")

    fig.suptitle(
        "Fig. 2. Differential regulon landscape and candidate TF prioritization in the astrocyte pilot",
        y=1.01,
    )
    save_figure(fig, out_dir, "Fig2_main_v2")

    notes = """
    Fig2_main_v2
    ============
    Data sources:
    - analysis_outputs/group_compare/regulon_group_statistics.csv
    - analysis_outputs/group_compare/top20_differential_regulons.csv
    - analysis_outputs/celloracle_candidates/celloracle_tf_shortlist.csv
    - final_exports/auc_mtx_matched_to_h5ad.csv
    - celloracle_run/prepared_data/celloracle_round1_metadata.csv

    Script entry:
    - scripts/manuscript_figures/build_main_figures_v2.py

    Panels:
    A. A main-text-focused differential regulon heatmap using 14 representative regulons. The four key regulons NFE2L2(+), THRB(+), BHLHE40(+), and SOX2(+) are forced to remain visible. Sample-level mean AUCell values are row-z transformed across the four pilot samples to avoid the binary appearance of a 2-column group-only heatmap.
    B. Differential regulon summary plot using mean AUC difference and -log10(FDR), with only the core regulons and a few supporting regulons labeled.
    C. Compact violin-plus-median/IQR comparison of single-cell regulon AUC. NFE2L2 and THRB are visually prioritized; BHLHE40 and SOX2 are retained as smaller secondary comparisons.
    D. Candidate TF prioritization scatter using TF expression log2FC, regulon significance, and positive fraction in the higher-expression group.

    Main-text revision:
    - The original top-20 group heatmap and boxplots were tightened into a more publication-oriented layout.
    - No upstream analysis was rerun. Exported AUCell matrix plus CellOracle metadata were used directly because local anndata import remains unavailable on this machine.
    """
    write_text(out_dir / "Fig2_notes_v2.txt", notes)
    return {"Fig2_v2": inputs}


def build_fig3() -> dict[str, list[str]]:
    out_dir = FIG_ROOT / "Fig3"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = [
        "celloracle_run/ko_round1/round1_recovery_index_ranking.csv",
        "celloracle_run/ko_round1/round1_ko_summary_by_group.csv",
        "celloracle_run/ko_round1/NFE2L2/NFE2L2_quiver.png",
        "celloracle_run/ko_round1/THRB/THRB_quiver.png",
    ]
    ranking = read_csv(inputs[0]).set_index("tf").loc[MAIN_TFS].copy()
    by_group = read_csv(inputs[1])
    by_group["mean_shift_length"] = as_num(by_group["mean_shift_length"])
    pivot = by_group.pivot(index="tf", columns="group", values="mean_shift_length").reindex(MAIN_TFS)

    fig = plt.figure(figsize=(10.8, 7.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.05], wspace=0.08)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    order = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
    y = np.arange(len(order))[::-1]
    for yi, tf in zip(y, order):
        shift = float(ranking.loc[tf, "mean_shift_length"])
        recovery = float(ranking.loc[tf, "recovery_index"])
        s_rank = int(ranking.loc[tf, "overall_shift_rank"])
        r_rank = int(ranking.loc[tf, "recovery_index_rank"])
        ax_a.plot([shift, recovery], [yi, yi], color="#C9C9C9", linewidth=1.4, zorder=1)
        ax_a.scatter([shift], [yi], s=58, color=COLORS[tf], edgecolor="white", linewidth=0.8, zorder=3)
        ax_a.scatter([recovery], [yi], s=50, facecolor="white", edgecolor=COLORS[tf], marker="s", linewidth=1.1, zorder=4)
        ax_a.text(max(shift, recovery) + 0.014, yi, f"shift #{s_rank} | RI #{r_rank}", fontsize=6.0, va="center", color="#585858")
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(order)
    for tick in ax_a.get_yticklabels():
        if tick.get_text() in {"NFE2L2", "THRB"}:
            tick.set_fontweight("bold")
    ax_a.set_xlim(0.37, 0.76)
    ax_a.set_xlabel("Perturbation-derived score")
    ax_a.set_title("Integrated perturbation prioritization", pad=6)
    style_axes(ax_a)
    ax_a.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="None", markersize=5.5, markerfacecolor="#888888", markeredgecolor="white", label="Mean shift length"),
            Line2D([0], [0], marker="s", linestyle="None", markersize=5.0, markerfacecolor="white", markeredgecolor="#666666", label="Recovery index"),
        ],
        frameon=False,
        loc="lower right",
    )
    panel_label(ax_a, "A")

    x_pos = {"lesion": 0, "internal_control": 1}
    for tf in order:
        y_vals = [float(pivot.loc[tf, "lesion"]), float(pivot.loc[tf, "internal_control"])]
        ax_b.plot([0, 1], y_vals, color=COLORS[tf], linewidth=1.7 if tf in {"NFE2L2", "THRB"} else 1.2, alpha=0.95)
        ax_b.scatter([0, 1], y_vals, s=34 if tf in {"NFE2L2", "THRB"} else 28, color=COLORS[tf], edgecolor="white", linewidth=0.7, zorder=3)
        ax_b.text(1.04, y_vals[1], tf, color=COLORS[tf], fontsize=6.2, va="center", fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal")
    ax_b.set_xlim(-0.12, 1.34)
    ax_b.set_xticks([0, 1])
    ax_b.set_xticklabels(["lesion", "internal_control"])
    ax_b.get_xticklabels()[0].set_color(COLORS["lesion"])
    ax_b.get_xticklabels()[1].set_color(COLORS["internal_control"])
    ax_b.set_ylabel("Mean shift length")
    ax_b.set_title("Group-specific perturbation comparison", pad=6)
    style_axes(ax_b)
    panel_label(ax_b, "B")

    for ax, img_rel, title, label in [
        (ax_c, inputs[2], "NFE2L2 representative perturbation", "C"),
        (ax_d, inputs[3], "THRB representative perturbation", "D"),
    ]:
        img = trim_white(mpimg.imread(ROOT / img_rel))
        ax.imshow(img)
        ax.set_title(title, pad=6)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        panel_label(ax, label)

    fig.suptitle(
        "Fig. 3. CellOracle perturbation prioritization and integrated ranking of candidate TFs",
        y=1.01,
    )
    save_figure(fig, out_dir, "Fig3_main_v2")

    notes = """
    Fig3_main_v2
    ============
    Data sources:
    - celloracle_run/ko_round1/round1_recovery_index_ranking.csv
    - celloracle_run/ko_round1/round1_ko_summary_by_group.csv
    - celloracle_run/ko_round1/NFE2L2/NFE2L2_quiver.png
    - celloracle_run/ko_round1/THRB/THRB_quiver.png

    Script entry:
    - scripts/manuscript_figures/build_main_figures_v2.py

    Panels:
    A. Dual-metric perturbation prioritization panel combining mean shift length and Recovery index in one aligned view.
    B. Group-specific perturbation comparison shown as a slope graph rather than grouped bars.
    C. NFE2L2 representative CellOracle perturbation panel using the current round1 exported quiver image.
    D. THRB representative CellOracle perturbation panel using the current round1 exported quiver image.

    Main-text revision:
    - Separate ranking bar charts were merged into one compact panel.
    - The grouped bar chart was replaced by a slope-style comparison to reduce dashboard-like appearance.
    - Recovery index is still presented as an approximate pathological reversal index derived from the current state-shift workflow, not as a newly re-derived strict geometric model.
    - BHLHE40 was retained in the ranking panels but not given a dedicated vector-field panel to keep the figure focused on the dual-core NFE2L2–THRB axis.
    """
    write_text(out_dir / "Fig3_notes_v2.txt", notes)
    return {"Fig3_v2": inputs}


def build_fig4() -> dict[str, list[str]]:
    out_dir = FIG_ROOT / "Fig4"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = [
        "robustness_validation/05_robustness_validation_integrated_table.csv",
        "robustness_validation/04_shortlist_priority_stability.csv",
        "external_validation/single_group_support/external_support_ranking.csv",
        "external_validation_round2/external_round2_support_ranking.csv",
    ]
    robust = read_csv(inputs[0]).set_index("tf")
    stability = read_csv(inputs[1]).set_index("tf").loc[MAIN_TFS].copy()
    gse140 = read_csv(inputs[2]).set_index("tf").loc[MAIN_TFS].copy()
    gse190 = read_csv(inputs[3]).set_index("tf").loc[MAIN_TFS].copy()

    rob_cols = ["loo_expr_consistency_rate", "loo_regulon_consistency_rate", "shortlist_retention_frequency"]
    rob_mat = robust.loc[MAIN_TFS, rob_cols].apply(pd.to_numeric, errors="coerce")

    fig = plt.figure(figsize=(12.0, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 0.82, 1.08], wspace=0.08)
    ax_a = fig.add_subplot(gs[:, 0])
    bgs = gs[0, 1].subgridspec(1, 2, width_ratios=[1.12, 0.95], wspace=0.18)
    ax_b1 = fig.add_subplot(bgs[0, 0])
    ax_b2 = fig.add_subplot(bgs[0, 1], sharey=ax_b1)
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[:, 2])

    im = ax_a.imshow(rob_mat.values, aspect="auto", cmap="Blues", vmin=0.6, vmax=1.0)
    ax_a.set_xticks(np.arange(len(rob_cols)))
    ax_a.set_xticklabels(["LOSO expr", "LOSO regulon", "Shortlist retention"], rotation=28, ha="right")
    ax_a.set_yticks(np.arange(len(MAIN_TFS)))
    ax_a.set_yticklabels(MAIN_TFS)
    ax_a.set_title("Sample-level robustness heatmap", pad=6)
    for i in range(rob_mat.shape[0]):
        for j in range(rob_mat.shape[1]):
            val = rob_mat.iloc[i, j]
            ax_a.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6.2, color="white" if val >= 0.84 else "#1F1F1F")
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax_a, fraction=0.05, pad=0.03)
    cbar.set_label("rate / frequency")
    cbar.ax.tick_params(labelsize=6.0)
    panel_label(ax_a, "A")

    y = np.arange(len(MAIN_TFS))[::-1]
    tf_order = MAIN_TFS
    ret = as_num(stability.loc[tf_order, "retention_frequency"])
    ranks = as_num(stability.loc[tf_order, "median_rank_when_retained"])
    for yi, tf in zip(y, tf_order):
        val = float(ret.loc[tf])
        ax_b1.plot([0, val], [yi, yi], color="#D7D7D7", linewidth=1.1, zorder=1)
        ax_b1.scatter([val], [yi], s=42, color=COLORS[tf], edgecolor="white", linewidth=0.7, zorder=3)
    ax_b1.set_xlim(0, 1.05)
    ax_b1.set_yticks(y)
    ax_b1.set_yticklabels(tf_order)
    for tick in ax_b1.get_yticklabels():
        if tick.get_text() in {"NFE2L2", "THRB"}:
            tick.set_fontweight("bold")
    ax_b1.set_xlabel("Retention freq.")
    ax_b1.set_title("Shortlist robustness summary", pad=6)
    style_axes(ax_b1, grid_axis="x")
    panel_label(ax_b1, "B")

    for yi, tf in zip(y, tf_order):
        val = float(ranks.loc[tf])
        ax_b2.scatter([val], [yi], s=38, color=COLORS[tf], marker="s", edgecolor="white", linewidth=0.7, zorder=3)
    ax_b2.set_xlim(0.5, max(6.2, float(np.nanmax(ranks.values)) + 0.5))
    ax_b2.set_xlabel("Median rank")
    ax_b2.set_yticks(y)
    ax_b2.tick_params(labelleft=False)
    style_axes(ax_b2, grid_axis="x")

    order140 = gse140.sort_values("support_rank", ascending=True).index.tolist()
    y140 = np.arange(len(order140))[::-1]
    for yi, tf in zip(y140, order140):
        score = float(gse140.loc[tf, "support_score"])
        ax_c.plot([0, score], [yi, yi], color="#D7D7D7", linewidth=1.0, zorder=1)
        ax_c.scatter([score], [yi], s=46, color=COLORS[tf], edgecolor="white", linewidth=0.7, zorder=3)
        ax_c.text(score + 0.025, yi, str(gse140.loc[tf, "expected_higher_group_main"]), fontsize=6.0, va="center", color="#555555")
    ax_c.set_xlim(0, 1.05)
    ax_c.set_yticks(y140)
    ax_c.set_yticklabels(order140)
    ax_c.set_xlabel("Support score")
    ax_c.set_title("GSE140393 single-group supportive analysis", pad=6)
    style_axes(ax_c, grid_axis="x")
    panel_label(ax_c, "C")

    ax_d.axvline(0, color="#A0A0A0", linewidth=0.8)
    ax_d.axhline(0.6, color="#D5D5D5", linewidth=0.8, linestyle="--")
    for tf in ["NFE2L2", "THRB", "SOX2", "BHLHE40"]:
        row = gse190.loc[tf]
        x = float(row["group_log2fc_lesion_vs_internal_control"])
        yv = float(row["cluster_support_fraction"])
        size = 140 + 680 * float(row["support_score"])
        consistent = bool(row["group_direction_consistent_with_main"])
        marker = "o" if consistent else "X"
        ax_d.scatter([x], [yv], s=size, color=COLORS[tf], marker=marker, alpha=0.82, edgecolor="white", linewidth=0.8, zorder=3)
    label_offsets = {
        "NFE2L2": (4, 6),
        "THRB": (4, 6),
        "SOX2": (4, 6),
        "BHLHE40": (4, 4),
    }
    for tf in ["NFE2L2", "THRB", "SOX2", "BHLHE40"]:
        row = gse190.loc[tf]
        dx, dy = label_offsets[tf]
        ax_d.annotate(
            tf,
            (float(row["group_log2fc_lesion_vs_internal_control"]), float(row["cluster_support_fraction"])),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=6.6,
            fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal",
        )
    ax_d.set_xlim(-0.62, 0.45)
    ax_d.set_ylim(0.49, 0.84)
    ax_d.set_xlabel("Group log2FC (TLE - non-epileptic control)")
    ax_d.set_ylabel("Cluster support fraction")
    ax_d.set_title("GSE190452 cross-syndrome supportive analysis", pad=6)
    style_axes(ax_d)
    ax_d.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="None", markersize=5.0, markerfacecolor="#888888", markeredgecolor="white", label="Direction-consistent"),
            Line2D([0], [0], marker="X", linestyle="None", markersize=5.0, markerfacecolor="#888888", markeredgecolor="white", label="Direction-attenuated"),
        ],
        frameon=False,
        loc="lower right",
    )
    panel_label(ax_d, "D")

    fig.suptitle(
        "Fig. 4. Sample-level robustness and supportive external evidence for the NFE2L2–THRB dual-axis model",
        y=1.01,
    )
    save_figure(fig, out_dir, "Fig4_main_v2")

    notes = """
    Fig4_main_v2
    ============
    Data sources:
    - robustness_validation/05_robustness_validation_integrated_table.csv
    - robustness_validation/04_shortlist_priority_stability.csv
    - external_validation/single_group_support/external_support_ranking.csv
    - external_validation_round2/external_round2_support_ranking.csv

    Script entry:
    - scripts/manuscript_figures/build_main_figures_v2.py

    Panels:
    A. Sample-level robustness heatmap summarizing leave-one-sample-out expression consistency, leave-one-sample-out regulon consistency, and shortlist retention frequency.
    B. Compact shortlist robustness summary using retention frequency and median rank when retained in a dot-based format rather than a bar-plus-line plot.
    C. GSE140393 single-group supportive panel shown as a ranked lollipop display. This remains supportive evidence only, not formal validation.
    D. GSE190452 cross-syndrome supportive panel using group log2FC and cluster support fraction. Marker shape distinguishes direction-consistent versus direction-attenuated support; BHLHE40 therefore remains visibly weakened.

    Main-text revision:
    - The robustness panel keeps sample-level wording and does not revert to donor-level language.
    - The two supportive datasets are explicitly kept in a supportive-evidence role rather than a formal external-validation role.
    """
    write_text(out_dir / "Fig4_notes_v2.txt", notes)
    return {"Fig4_v2": inputs}


def write_revision_summary(table_paths: list[Path]) -> None:
    lines = [
        "# Figure revision summary v2",
        "",
        "Generated by `scripts/manuscript_figures/build_main_figures_v2.py`.",
        "",
        "## Fig2",
        "",
        "- Replaced the previous top-20 group-only heatmap with a tighter 14-regulon heatmap built from sample-level mean AUCell values across the four pilot samples, while retaining row-z scaling.",
        "- Reduced label density in the differential summary panel to the core 6-8 regulons.",
        "- Replaced the more basic boxplot presentation with a compact violin-plus-median/IQR panel that visually prioritizes NFE2L2 and THRB.",
        "- Kept the shortlist bubble panel as the closing panel, with cleaner legend handling and clearer lesion-associated versus internal_control-associated separation.",
        "",
        "## Fig3",
        "",
        "- Replaced two separate ranking bar charts with one dual-metric integrated prioritization panel.",
        "- Replaced the grouped perturbation bar chart with a slope-style lesion versus internal_control comparison.",
        "- Kept only the NFE2L2 and THRB representative CellOracle vector panels in the main figure; BHLHE40 remains visible in the ranking panels and can stay supplementary as a dedicated vector panel if needed.",
        "",
        "## Fig4",
        "",
        "- Kept the robustness heatmap but tightened typography and scale for double-column use.",
        "- Replaced the retention-frequency plus line overlay panel with a compact two-metric dot summary.",
        "- Replaced the GSE140393 support bar chart with a ranked lollipop panel.",
        "- Strengthened the GSE190452 supportive panel as the right-side visual anchor and preserved BHLHE40 attenuation rather than smoothing it away.",
        "",
        "## Content moved away from the main figure style",
        "",
        "- Separate ranking bars and dashboard-like grouped bars were removed from the main-figure layout.",
        "- A dedicated BHLHE40 CellOracle vector panel was left out of the main figure to avoid overloading the central dual-axis story.",
        "- Full drug-repositioning outputs, larger GO/KEGG tables, and long external-support detail tables remain better suited to supplementary material.",
        "",
        "## Mainline consistency check",
        "",
        "- The revised figures stay on the GSE268807 astrocyte pilot + pySCENIC + CellOracle + sample-level robustness + functional interpretation/program convergence + supportive external analysis mainline.",
        "- No old SCENIC+ route was used.",
        "- GSE140393 and GSE190452 are still described as supportive evidence, not formal external validation.",
        "- No upstream analyses were rerun and no result values were altered.",
        "",
        "## Manuscript-ready tables",
        "",
    ]
    for path in table_paths:
        lines.append(f"- `{path.relative_to(ROOT).as_posix()}`")
    write_text(FIG_ROOT / "FIGURE_REVISION_SUMMARY_V2.md", "\n".join(lines))


def main() -> None:
    setup_style()
    fig_inputs: dict[str, list[str]] = {}

    shortlist = read_csv("analysis_outputs/celloracle_candidates/celloracle_tf_shortlist.csv")
    ranking = read_csv("celloracle_run/ko_round1/round1_recovery_index_ranking.csv")
    robust = read_csv("robustness_validation/05_robustness_validation_integrated_table.csv")
    gse140 = read_csv("external_validation/single_group_support/external_support_ranking.csv")
    gse190 = read_csv("external_validation_round2/external_round2_support_ranking.csv")

    fig_inputs.update(build_fig2())
    fig_inputs.update(build_fig3())
    fig_inputs.update(build_fig4())

    table_paths = [
        build_table2_shortlist(shortlist),
        build_table3_integrated(shortlist, ranking, robust, gse140, gse190),
    ]
    write_revision_summary(table_paths)
    print(f"Revised figures written to: {FIG_ROOT}")
    for path in table_paths:
        print(f"Table written to: {path}")


if __name__ == "__main__":
    main()
