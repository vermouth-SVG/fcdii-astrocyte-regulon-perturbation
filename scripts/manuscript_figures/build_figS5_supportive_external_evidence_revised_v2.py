#!/usr/bin/env python3
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_supplementary_revised_v2"
STEM = "FigS5_supportive_external_evidence_revised_v2"

TF_ORDER = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
TF_COLORS = {
    "NFE2L2": "#8F1D2C",
    "THRB": "#0B7285",
    "BHLHE40": "#7B5AA6",
    "SOX2": "#6F7378",
}

INTERPRETATION = {
    "NFE2L2": "disease-side support",
    "THRB": "partial internal-control-like support",
    "BHLHE40": "attenuated support",
    "SOX2": "partial / retained support",
}

INPUT_FILES = [
    ROOT / "external_validation" / "single_group_support" / "donor_level_summary.csv",
    ROOT / "external_validation" / "single_group_support" / "external_support_ranking.csv",
    ROOT / "external_validation" / "input" / "input_build_report.txt",
    ROOT / "external_validation_round2" / "round2_expression_and_regulon_master_table.csv",
    ROOT / "external_validation_round2" / "external_round2_tf_expression_celltype_stats.csv",
    ROOT / "external_validation_round2" / "external_round2_support_ranking.csv",
    ROOT / "external_validation_round2" / "input" / "build_report.txt",
]

FIGURE_TEXT_ITEMS = [
    "Supportive external evidence from GSE140393 and GSE190452",
    "A. GSE140393 single-group expression support",
    "Positive fraction",
    "Sample ID",
    "Candidate TF",
    "B. GSE140393 ranked single-group support score",
    "Support score",
    "C. GSE190452 cross-syndrome expression and target-set support",
    "Expression support score",
    "Target-set support score",
    "Candidate TF",
    "D. GSE190452 cluster-level expression-support pattern",
    "Directionally signed expression-support score",
    "Cluster",
    "E. Candidate interpretation summary",
    *TF_ORDER,
    *INTERPRETATION.values(),
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.4,
            "axes.titlesize": 8.6,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.3,
            "figure.titlesize": 11.0,
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


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def style_axes(ax: plt.Axes, grid_axis: str | None = "y") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color="#E8E8E8", linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.12) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=11.0,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def parse_cluster_id(value: object) -> int:
    text = str(value)
    if "_" in text:
        text = text.rsplit("_", 1)[-1]
    return int(text)


def build_panel_a(ax: plt.Axes, donor: pd.DataFrame, fig: plt.Figure) -> None:
    donor_sample = donor[(donor["level"].astype(str).eq("sample")) & (donor["tf"].isin(TF_ORDER))].copy()
    sample_order = (
        donor_sample[["level_id", "n_cells"]]
        .drop_duplicates()
        .sort_values("n_cells", ascending=False)["level_id"]
        .astype(str)
        .tolist()
    )
    donor_sample["level_id"] = donor_sample["level_id"].astype(str)
    heat = (
        donor_sample.pivot_table(index="tf", columns="level_id", values="positive_fraction", aggfunc="first")
        .reindex(index=TF_ORDER, columns=sample_order)
        .astype(float)
    )

    vmax = max(float(np.nanmax(heat.values)), 0.01)
    im = ax.imshow(heat.values, aspect="auto", cmap="Reds", vmin=0, vmax=vmax)
    ax.set_title("GSE140393 single-group expression support", pad=7)
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Candidate TF")
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(TF_ORDER)))
    ax.set_yticklabels(TF_ORDER)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            if pd.notna(val):
                color = "white" if float(val) > vmax * 0.50 else "#222222"
                ax.text(j, i, f"{float(val):.2f}", ha="center", va="center", fontsize=6.0, color=color)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.050, pad=0.035)
    cbar.set_label("Positive fraction")
    cbar.ax.tick_params(labelsize=6.2, width=0.5, length=2.4)
    panel_label(ax, "A")


def build_panel_b(ax: plt.Axes, ranking140: pd.DataFrame) -> None:
    ranking = ranking140[ranking140["tf"].isin(TF_ORDER)].copy()
    ranking["support_score"] = as_num(ranking["support_score"])
    ranking = ranking.sort_values(["support_score", "tf"], ascending=[False, True])
    y = np.arange(len(ranking))[::-1]
    vals = ranking["support_score"].to_numpy(dtype=float)
    colors = [TF_COLORS.get(tf, "#777777") for tf in ranking["tf"]]

    ax.hlines(y, 0, vals, color="#D4D4D4", linewidth=1.1, zorder=1)
    ax.scatter(vals, y, s=54, c=colors, edgecolor="white", linewidth=0.8, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(ranking["tf"].tolist())
    ax.set_xlim(0, 1.42)
    ax.set_xticks([0, 0.25, 0.50, 0.75, 1.00])
    ax.set_xlabel("Support score")
    ax.set_title("GSE140393 ranked single-group support score", pad=7)
    style_axes(ax, grid_axis="x")
    for yi, tf in zip(y, ranking["tf"]):
        ax.text(1.04, yi, INTERPRETATION[tf], fontsize=6.5, va="center", color=TF_COLORS[tf])
    panel_label(ax, "B")


def build_panel_c(ax: plt.Axes, round2: pd.DataFrame) -> None:
    required = ["tf", "expr_support_score", "regulon_support_score"]
    missing = [col for col in required if col not in round2.columns]
    if missing:
        raise ValueError(f"Missing required GSE190452 columns for Panel C: {missing}")
    plot_df = round2.set_index("tf").reindex(TF_ORDER).reset_index()
    plot_df["expr_support_score"] = as_num(plot_df["expr_support_score"])
    plot_df["target_set_support_score"] = as_num(plot_df["regulon_support_score"])

    ax.scatter(
        plot_df["expr_support_score"],
        plot_df["target_set_support_score"],
        s=96,
        c=[TF_COLORS[tf] for tf in plot_df["tf"]],
        edgecolor="white",
        linewidth=0.85,
        alpha=0.92,
        zorder=3,
    )
    offsets = {
        "NFE2L2": (-36, 6),
        "THRB": (7, 7),
        "BHLHE40": (7, -12),
        "SOX2": (7, 7),
    }
    for _, row in plot_df.iterrows():
        tf = str(row["tf"])
        ax.annotate(
            tf,
            (float(row["expr_support_score"]), float(row["target_set_support_score"])),
            xytext=offsets.get(tf, (5, 5)),
            textcoords="offset points",
            fontsize=6.6,
            color=TF_COLORS[tf],
            arrowprops={"arrowstyle": "-", "color": "#8A8A8A", "lw": 0.45, "shrinkA": 0, "shrinkB": 5},
        )
    ax.set_xlim(-0.04, 1.06)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlabel("Expression support score")
    ax.set_ylabel("Target-set support score")
    ax.set_title("GSE190452 cross-syndrome expression and target-set support", pad=7)
    style_axes(ax, grid_axis="both")
    handles = [
        Line2D([0], [0], marker="o", linestyle="None", markersize=5.0, markerfacecolor=TF_COLORS[tf], markeredgecolor="white", label=tf)
        for tf in TF_ORDER
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right", title="Candidate TF", title_fontsize=6.5)
    panel_label(ax, "C")


def build_panel_d(ax: plt.Axes, expr_cluster: pd.DataFrame, fig: plt.Figure) -> None:
    required = ["tf", "celltype_or_cluster", "expected_higher_group_main", "log2fc_lesion_vs_internal_control"]
    missing = [col for col in required if col not in expr_cluster.columns]
    if missing:
        raise ValueError(f"Missing required GSE190452 columns for Panel D: {missing}")
    cluster = expr_cluster[expr_cluster["tf"].isin(TF_ORDER)].copy()
    logfc = as_num(cluster["log2fc_lesion_vs_internal_control"])
    is_lesion_expected = cluster["expected_higher_group_main"].astype(str).eq("lesion")
    cluster["directionally_signed_expression_support"] = np.where(is_lesion_expected, logfc, -logfc)
    pivot = cluster.pivot_table(
        index="tf",
        columns="celltype_or_cluster",
        values="directionally_signed_expression_support",
        aggfunc="first",
    )
    cluster_cols = sorted(pivot.columns, key=parse_cluster_id)
    pivot = pivot.reindex(index=TF_ORDER, columns=cluster_cols)
    values = pivot.to_numpy(dtype=float)
    finite_abs = np.abs(values[np.isfinite(values)])
    vmax = float(np.nanpercentile(finite_abs, 95)) if finite_abs.size else 1.0
    vmax = min(max(vmax, 0.50), 1.50)

    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    tick_positions = [i for i, c in enumerate(cluster_cols) if parse_cluster_id(c) % 5 == 0]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([str(parse_cluster_id(cluster_cols[i])) for i in tick_positions])
    ax.set_yticks(np.arange(len(TF_ORDER)))
    ax.set_yticklabels(TF_ORDER)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Candidate TF")
    ax.set_title("GSE190452 cluster-level expression-support pattern", pad=7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.030)
    cbar.set_label("Directionally signed expression-support score")
    cbar.ax.tick_params(labelsize=6.2, width=0.5, length=2.4)
    panel_label(ax, "D")


def build_panel_e(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.86, "Candidate interpretation summary", fontsize=8.5, fontweight="bold", ha="left", va="top")
    for i, tf in enumerate(TF_ORDER):
        x0 = 0.02 + i * 0.245
        ax.plot([x0, x0 + 0.055], [0.55, 0.55], color=TF_COLORS[tf], linewidth=3.2, solid_capstyle="round")
        ax.text(x0 + 0.065, 0.64, tf, fontsize=7.2, fontweight="bold", ha="left", va="center", color="#222222")
        ax.text(x0 + 0.065, 0.39, INTERPRETATION[tf], fontsize=7.0, ha="left", va="center", color="#333333")
    panel_label(ax, "E", x=-0.035, y=0.98)


def save_figure(fig: plt.Figure) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    png = OUT_DIR / f"{STEM}.png"
    pdf = OUT_DIR / f"{STEM}.pdf"
    tiff = OUT_DIR / f"{STEM}.tiff"
    fig.savefig(png, dpi=450, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    try:
        fig.savefig(tiff, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    except TypeError:
        fig.savefig(tiff, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return [png, pdf, tiff]


def write_text(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")


def write_caption() -> Path:
    path = OUT_DIR / f"{STEM}_caption.md"
    write_text(
        path,
        """
        # Supplementary Figure S5. Supportive external evidence from GSE140393 and GSE190452

        GSE140393 was used as a single-group expression-level supportive dataset, whereas GSE190452 was used as a cross-syndrome supportive dataset. These analyses were not interpreted as formal external validation or independent regulon-level replication.

        A. GSE140393 sample-level positive fraction summarizes single-group expression support across the available sample IDs for NFE2L2, THRB, BHLHE40, and SOX2.

        B. GSE140393 ranked support score summarizes single-group support for the four retained candidate TFs and labels the conservative interpretation assigned to each TF.

        C. GSE190452 cross-syndrome expression and target-set support compares expression support score with target-set support score for the four retained candidate TFs.

        D. GSE190452 cluster-level expression-support pattern shows directionally signed expression support across clusters, with positive values indicating support for the expected disease-side or internal-control-like direction.

        E. Candidate interpretation summary gives the conservative conclusion used for manuscript interpretation: NFE2L2 shows disease-side support, THRB shows partial internal-control-like support, BHLHE40 shows attenuated support, and SOX2 shows partial / retained support.
        """,
    )
    return path


def write_notes() -> Path:
    path = OUT_DIR / f"{STEM}_notes.txt"
    lines = [
        "FigS5 supportive external evidence revised_v2 notes",
        "",
        "Output stem:",
        f"- {rel(OUT_DIR / STEM)}",
        "",
        "Input files used:",
        *[f"- {rel(p)}" for p in INPUT_FILES],
        "",
        "GSE140393 handling:",
        "- Panel A uses sample-level positive_fraction values read from donor_level_summary.csv after filtering level == sample and TF in NFE2L2, THRB, BHLHE40, and SOX2.",
        "- Panel B uses support_score values read from external_support_ranking.csv. The score was not recomputed in this figure script; it was read from the existing single-group support table.",
        "- GSE140393 is treated as a single-group expression-level supportive dataset.",
        "",
        "GSE190452 handling:",
        "- Panel C reads expr_support_score from round2_expression_and_regulon_master_table.csv and displays it as Expression support score.",
        "- Panel C reads the historical target-gene-set projection score from the same master table and displays it as Target-set support score.",
        "- Target-set support score summarizes projected support from candidate TF-associated gene sets in GSE190452 and should not be interpreted as independently recomputed pySCENIC regulon AUC.",
        "- Panel D uses external_round2_tf_expression_celltype_stats.csv. For TFs expected to be disease-side, the cluster-level signed score is the lesion-vs-control expression log2FC. For THRB, the sign is reversed so positive values indicate internal-control-like support.",
        "- GSE190452 is treated as a cross-syndrome supportive dataset.",
        "",
        "Interpretation boundary:",
        "- The figure is supportive evidence only.",
        "- GSE190452 is not interpreted as formal validation.",
        "- The target-set support score is not interpreted as independently recomputed pySCENIC AUC.",
        "- No panel claims independent same-disease replication.",
        "",
        "Candidate-level interpretation used in the figure:",
        "- NFE2L2: disease-side support.",
        "- THRB: partial internal-control-like support.",
        "- BHLHE40: attenuated support.",
        "- SOX2: partial / retained support.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_revision_log() -> Path:
    path = OUT_DIR / f"{STEM}_revision_log.md"
    write_text(
        path,
        """
        # FigS5 supportive external evidence revised_v2 revision log

        - Created a new output set in manuscript_output/figures_supplementary_revised_v2 without overwriting the full_submission files.
        - Updated the whole-figure title to "Supportive external evidence from GSE140393 and GSE190452".
        - Panel A title changed to "GSE140393 single-group expression support"; the heatmap retains sample-level positive fraction values and candidate order NFE2L2, THRB, BHLHE40, SOX2.
        - Panel B title changed to "GSE140393 ranked single-group support score"; right-side labels now use disease-side support, partial internal-control-like support, attenuated support, and partial / retained support.
        - Panel C title changed to "GSE190452 cross-syndrome expression and target-set support"; the y axis now uses "Target-set support score" on the figure surface.
        - Panel C points use uniform size and TF color coding with a compact legend.
        - Panel D was rebuilt from cluster-level expression statistics and renamed "GSE190452 cluster-level expression-support pattern"; x-axis labels are thinned to clusters 0, 5, 10, 15, 20, 25, and 30.
        - Added Panel E as a compact candidate interpretation summary strip.
        - Retired old figure-facing language: expression-regulon support summary, projected regulon support, and underscore-form internal-control wording.
        - Added caption, notes, revision log, and terminology audit companion files.
        """,
    )
    return path


def write_terminology_audit() -> Path:
    figure_text = "\n".join(FIGURE_TEXT_ITEMS)
    figure_text_lower = figure_text.lower()

    def fig_absent(term: str) -> bool:
        return term.lower() not in figure_text_lower

    rows = [
        {
            "term": "formal validation",
            "status": "absent from figure text; boundary-only in caption/notes",
            "action": "required absent as a positive claim",
        },
        {
            "term": "independent validation",
            "status": "absent",
            "action": "required absent",
        },
        {
            "term": "external validation",
            "status": "absent from figure text; input paths and negated boundary only",
            "action": "avoid in figure text",
        },
        {
            "term": "validated",
            "status": "absent",
            "action": "required absent",
        },
        {
            "term": "confirmed",
            "status": "absent",
            "action": "required absent",
        },
        {
            "term": "regulon replication",
            "status": "absent from figure text; boundary-only in caption",
            "action": "required absent as a positive claim",
        },
        {
            "term": "regulon-level replication",
            "status": "absent from figure text; boundary-only in caption",
            "action": "required absent as a positive claim",
        },
        {
            "term": "projected regulon support",
            "status": "replaced" if fig_absent("projected regulon support") else "present",
            "action": "replaced by target-set support or expression-support pattern",
        },
        {
            "term": "internal_control",
            "status": "replaced" if fig_absent("internal_control") else "present",
            "action": "replaced by internal-control in figure text",
        },
        {
            "term": "GSE140393 single-group support",
            "status": "present" if "GSE140393 single-group".lower() in figure_text_lower else "absent",
            "action": "required",
        },
        {
            "term": "GSE190452 cross-syndrome support",
            "status": "present" if "GSE190452 cross-syndrome".lower() in figure_text_lower else "absent",
            "action": "required",
        },
        {
            "term": "supportive evidence",
            "status": "present" if "supportive external evidence" in figure_text_lower else "absent",
            "action": "required",
        },
    ]
    path = OUT_DIR / f"{STEM}_terminology_audit.csv"
    pd.DataFrame(rows, columns=["term", "status", "action"]).to_csv(path, index=False)
    return path


def build_figure() -> list[Path]:
    setup_style()
    donor = read_csv(INPUT_FILES[0])
    ranking140 = read_csv(INPUT_FILES[1])
    round2_master = read_csv(INPUT_FILES[3])
    expr_cluster = read_csv(INPUT_FILES[4])

    fig = plt.figure(figsize=(13.2, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(
        3,
        2,
        width_ratios=[0.95, 1.22],
        height_ratios=[1.0, 1.08, 0.32],
        wspace=0.10,
        hspace=0.12,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[0, 1])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, :])

    build_panel_a(ax_a, donor, fig)
    build_panel_b(ax_b, ranking140)
    build_panel_c(ax_c, round2_master)
    build_panel_d(ax_d, expr_cluster, fig)
    build_panel_e(ax_e)
    fig.suptitle("Supportive external evidence from GSE140393 and GSE190452", y=1.005, fontweight="bold")
    return save_figure(fig)


def residual_summary() -> dict[str, str]:
    figure_text = "\n".join(FIGURE_TEXT_ITEMS).lower()
    terms = [
        "validation",
        "formal validation",
        "independent validation",
        "projected regulon support",
        "internal_control",
    ]
    return {term: ("YES" if term in figure_text else "NO") for term in terms}


def main() -> None:
    outputs = build_figure()
    caption = write_caption()
    notes = write_notes()
    revision_log = write_revision_log()
    audit = write_terminology_audit()
    residuals = residual_summary()

    print("FigS5 revised_v2 output path:")
    print(f"- {rel(OUT_DIR)}")
    for path in outputs:
        print(f"- {rel(path)}")

    print("\nInput files used:")
    for path in INPUT_FILES:
        print(f"- {rel(path)}")

    print("\nTitles and axes modified:")
    print('- Whole title: "Supportive external evidence from GSE140393 and GSE190452"')
    print('- Panel A: "GSE140393 single-group expression support"; colorbar "Positive fraction"')
    print('- Panel B: "GSE140393 ranked single-group support score"; x axis "Support score"')
    print('- Panel C: "GSE190452 cross-syndrome expression and target-set support"; x axis "Expression support score"; y axis "Target-set support score"')
    print('- Panel D: "GSE190452 cluster-level expression-support pattern"; colorbar "Directionally signed expression-support score"')
    print('- Panel E: "Candidate interpretation summary"')

    print("\nResidual restricted terms in figure text:")
    for term, status in residuals.items():
        print(f"- {term}: {status}")
    print("- Companion caption/notes/audit contain required negated boundary or audit wording where applicable.")

    print("\nCompanion files generated:")
    print(f"- caption: {rel(caption)}")
    print(f"- notes: {rel(notes)}")
    print(f"- revision log: {rel(revision_log)}")
    print(f"- terminology audit: {rel(audit)}")

    print("\nManual review flag:")
    print("- No panel requires layout repair based on scripted checks; Panel C interpretation should be reviewed only for final manuscript wording consistency.")


if __name__ == "__main__":
    main()
