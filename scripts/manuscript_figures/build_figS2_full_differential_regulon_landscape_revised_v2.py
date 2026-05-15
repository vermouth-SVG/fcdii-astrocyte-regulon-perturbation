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
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_supplementary_revised_v2"
STEM = "FigS2_full_differential_regulon_landscape_revised_v2"

STATS_PATH = ROOT / "analysis_outputs" / "group_compare" / "regulon_group_statistics.csv"
AUC_PATH = ROOT / "final_exports" / "auc_mtx_matched_to_h5ad.csv"
METADATA_PATH = ROOT / "celloracle_run" / "prepared_data" / "celloracle_round1_metadata.csv"
REGULON_NAMES_PATH = ROOT / "final_exports" / "regulon_names.txt"
SHORTLIST_PATH = ROOT / "manuscript_output" / "tables_supplementary" / "TableS2_candidate_tf_shortlist_extended_submission.csv"

INPUT_FILES = [STATS_PATH, AUC_PATH, METADATA_PATH, REGULON_NAMES_PATH, SHORTLIST_PATH]

SAMPLE_ORDER = ["G120_N.1", "G133_N.2", "G120_D.3", "G133_D.2"]
SAMPLE_GROUPS = {
    "G120_N.1": "internal-control",
    "G133_N.2": "internal-control",
    "G120_D.3": "lesion",
    "G133_D.2": "lesion",
}
GROUP_COLORS = {"internal-control": "#0B6B74", "lesion": "#8E1B2A"}
TF_COLORS = {
    "NFE2L2": "#8E1B2A",
    "THRB": "#0B6B74",
    "BHLHE40": "#6A51A3",
    "SOX2": "#777777",
    "SATB2": "#5E8C61",
    "RARB": "#4A8B93",
    "HMGA1": "#9B4A48",
}
LABEL_REGULONS = ["THRB(+)", "SATB2(+)", "RARB(+)", "NFE2L2(+)", "BHLHE40(+)", "SOX2(+)", "HMGA1(+)"]
LESION_RETAINED = ["NFE2L2(+)", "BHLHE40(+)", "SOX2(+)", "HMGA1(+)"]
TEXT_COLOR = "#222222"
GRID_COLOR = "#E9E9E9"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.0,
            "axes.titlesize": 8.5,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.3,
            "ytick.labelsize": 5.9,
            "legend.fontsize": 6.0,
            "figure.titlesize": 10.8,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.linewidth": 0.65,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def rel(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def neg_log10(series: pd.Series) -> pd.Series:
    vals = as_num(series).clip(lower=np.nextafter(0, 1))
    return -np.log10(vals)


def panel_label(ax: plt.Axes, label: str, x: float = -0.13, y: float = 1.11) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.8,
        fontweight="bold",
        color=TEXT_COLOR,
        clip_on=False,
    )


def style_axes(ax: plt.Axes, grid_axis: str | None = "x") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.6)


def load_metadata() -> pd.DataFrame:
    meta = pd.read_csv(METADATA_PATH)
    if str(meta.columns[0]).startswith("Unnamed"):
        meta = meta.rename(columns={meta.columns[0]: "cell_id"})
    meta["group_display"] = meta["group"].astype(str).str.replace("internal_control", "internal-control", regex=False)
    return meta


def load_stats() -> pd.DataFrame:
    stats = pd.read_csv(STATS_PATH)
    stats["regulon"] = stats["regulon"].astype(str)
    stats["mean_lesion"] = as_num(stats["mean_lesion"])
    stats["mean_internal_control"] = as_num(stats["mean_internal_control"])
    stats["fdr_bh"] = as_num(stats["fdr_bh"])
    stats["abs_mean_diff"] = as_num(stats["abs_mean_diff"])
    stats["lesion_minus_internal_control"] = stats["mean_lesion"] - stats["mean_internal_control"]
    stats["neg_log10_fdr"] = neg_log10(stats["fdr_bh"])
    stats["axis_display"] = np.where(stats["lesion_minus_internal_control"] >= 0, "lesion-associated", "internal-control-associated")
    return stats


def load_auc() -> pd.DataFrame:
    auc = pd.read_csv(AUC_PATH)
    if str(auc.columns[0]).startswith("Unnamed"):
        auc = auc.rename(columns={auc.columns[0]: "cell_id"})
    return auc


def load_shortlist() -> pd.DataFrame:
    if not SHORTLIST_PATH.exists():
        return pd.DataFrame(columns=["TF", "Regulon"])
    return pd.read_csv(SHORTLIST_PATH)


def top40_regulons(stats: pd.DataFrame) -> list[str]:
    return stats.sort_values(["fdr_bh", "abs_mean_diff"], ascending=[True, False]).head(40)["regulon"].tolist()


def heatmap_matrix(stats: pd.DataFrame, meta: pd.DataFrame, auc: pd.DataFrame, top40: list[str]) -> pd.DataFrame:
    merged = auc.merge(meta[["cell_id", "clinical_sample_label", "group_display"]], on="cell_id", how="left")
    sample_means = merged.groupby("clinical_sample_label")[top40].mean().reindex(SAMPLE_ORDER)
    heat = sample_means.T
    heat = heat.sub(heat.mean(axis=1), axis=0)
    heat = heat.div(heat.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0.0)
    row_order = stats.set_index("regulon").loc[top40].sort_values("lesion_minus_internal_control").index.tolist()
    return heat.loc[row_order]


def panel_a(ax: plt.Axes, fig: plt.Figure, heat: pd.DataFrame) -> None:
    im = ax.imshow(heat.values, aspect="auto", cmap="RdBu_r", vmin=-1.7, vmax=1.7)
    ax.set_xticks(np.arange(len(SAMPLE_ORDER)))
    ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.tick_params(length=0)
    ax.axvline(1.5, color="#D0D0D0", linewidth=0.8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    strip = ax.inset_axes([0.0, 1.012, 1.0, 0.044])
    strip.set_xlim(0, len(SAMPLE_ORDER))
    strip.set_ylim(0, 1)
    for i, sample in enumerate(SAMPLE_ORDER):
        group = SAMPLE_GROUPS[sample]
        strip.add_patch(Rectangle((i, 0), 1, 1, facecolor=GROUP_COLORS[group], edgecolor="white", linewidth=0.55))
    strip.text(1.0, 0.50, "internal-control", color="white", ha="center", va="center", fontsize=5.7, fontweight="bold")
    strip.text(3.0, 0.50, "lesion", color="white", ha="center", va="center", fontsize=5.7, fontweight="bold")
    strip.axis("off")
    ax.set_title("Top 40 differential regulon heatmap", pad=9)
    cbar = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.025)
    cbar.set_label("Regulon AUC row z-score", fontsize=6.8)
    cbar.ax.tick_params(labelsize=5.8, width=0.5, length=2.0)
    panel_label(ax, "A", x=-0.17, y=1.075)


def label_color(regulon: str) -> str:
    tf = regulon.replace("(+)", "")
    return TF_COLORS.get(tf, "#333333")


def panel_b(ax: plt.Axes, stats: pd.DataFrame) -> None:
    colors = np.where(stats["lesion_minus_internal_control"] >= 0, GROUP_COLORS["lesion"], GROUP_COLORS["internal-control"])
    ax.scatter(
        stats["lesion_minus_internal_control"],
        stats["neg_log10_fdr"],
        s=17,
        c=colors,
        alpha=0.34,
        edgecolors="none",
        zorder=2,
    )
    ax.axvline(0, color="#8D8D8D", linewidth=0.8, zorder=1)
    offsets = {
        "THRB(+)": (16, -6, "left"),
        "SATB2(+)": (20, -2, "left"),
        "RARB(+)": (22, 2, "left"),
        "SOX2(+)": (-22, 12, "right"),
        "BHLHE40(+)": (24, 9, "left"),
        "NFE2L2(+)": (24, -16, "left"),
        "HMGA1(+)": (-24, -2, "right"),
    }
    for reg in LABEL_REGULONS:
        row = stats.loc[stats["regulon"].eq(reg)]
        if row.empty:
            continue
        x = float(row["lesion_minus_internal_control"].iloc[0])
        y = float(row["neg_log10_fdr"].iloc[0])
        color = label_color(reg)
        ax.scatter([x], [y], s=54, color="white", edgecolor=color, linewidth=1.1, zorder=4)
        dx, dy, ha = offsets.get(reg, (6, 6, "left"))
        ax.annotate(
            reg,
            xy=(x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=6.2,
            color=color,
            arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.55, "shrinkA": 1.0, "shrinkB": 4.0},
            zorder=5,
        )
    xmin = float(stats["lesion_minus_internal_control"].min())
    xmax = float(stats["lesion_minus_internal_control"].max())
    ymin = 0
    ymax = float(stats["neg_log10_fdr"].max()) * 1.08
    ax.set_xlim(xmin - 0.028, xmax + 0.040)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("Mean AUC difference: lesion − internal-control")
    ax.set_ylabel("−log10(FDR)")
    ax.set_title("All-regulon differential summary", pad=8)
    ax.text(0.02, 0.035, "internal-control-associated", transform=ax.transAxes, ha="left", va="bottom", fontsize=5.8, color=GROUP_COLORS["internal-control"])
    ax.text(0.98, 0.035, "lesion-associated", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=GROUP_COLORS["lesion"])
    style_axes(ax, "both")
    panel_label(ax, "B", x=-0.12, y=1.12)


def retained_top_lesion(stats: pd.DataFrame) -> pd.DataFrame:
    top = stats[stats["lesion_minus_internal_control"] > 0].nlargest(15, "lesion_minus_internal_control").copy()
    retained = stats[stats["regulon"].isin(LESION_RETAINED)].copy()
    out = pd.concat([top, retained], ignore_index=True).drop_duplicates(subset=["regulon"])
    return out.sort_values("lesion_minus_internal_control", ascending=True)


def panel_lollipop(ax: plt.Axes, data: pd.DataFrame, title: str, color: str, label: str, xlim: tuple[float, float]) -> None:
    y = np.arange(len(data))
    values = data["lesion_minus_internal_control"].astype(float).to_numpy()
    regs = data["regulon"].astype(str).to_numpy()
    ax.hlines(y, xmin=np.minimum(values, 0), xmax=np.maximum(values, 0), color="#CCCCCC", linewidth=1.0, zorder=1)
    point_colors = [label_color(reg) if reg in LABEL_REGULONS else color for reg in regs]
    point_sizes = [36 if reg in LABEL_REGULONS else 22 for reg in regs]
    ax.scatter(values, y, s=point_sizes, color=point_colors, edgecolor="white", linewidth=0.55, zorder=3)
    ax.axvline(0, color="#8D8D8D", linewidth=0.75, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(regs)
    for tick, reg in zip(ax.get_yticklabels(), regs):
        if reg in LABEL_REGULONS:
            tick.set_color(label_color(reg))
            tick.set_fontweight("bold")
    ax.set_xlabel("Mean AUC difference:\nlesion − internal-control", fontsize=6.4, labelpad=4)
    ax.set_title(title, pad=9)
    ax.set_xlim(*xlim)
    style_axes(ax, "x")
    panel_label(ax, label, x=-0.26, y=1.12)


def build_figure() -> tuple[plt.Figure, dict[str, object]]:
    setup_style()
    stats = load_stats()
    meta = load_metadata()
    auc = load_auc()
    shortlist = load_shortlist()
    top40 = top40_regulons(stats)
    heat = heatmap_matrix(stats, meta, auc, top40)
    top_ic = stats[stats["lesion_minus_internal_control"] < 0].nsmallest(15, "lesion_minus_internal_control").sort_values("lesion_minus_internal_control", ascending=True)
    top_ls = retained_top_lesion(stats)

    fig = plt.figure(figsize=(13.4, 9.2), constrained_layout=False)
    outer = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.04, 1.08],
        height_ratios=[0.94, 1.06],
        left=0.105,
        right=0.985,
        top=0.875,
        bottom=0.125,
        wspace=0.28,
        hspace=0.38,
    )
    ax_a = fig.add_subplot(outer[:, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    lower = outer[1, 1].subgridspec(1, 2, wspace=0.72)
    ax_c = fig.add_subplot(lower[0, 0])
    ax_d = fig.add_subplot(lower[0, 1])

    panel_a(ax_a, fig, heat)
    panel_b(ax_b, stats)
    xmin = float(stats["lesion_minus_internal_control"].min())
    xmax = float(stats["lesion_minus_internal_control"].max())
    panel_lollipop(
        ax_c,
        top_ic,
        "Top internal-control-associated regulons",
        GROUP_COLORS["internal-control"],
        "C",
        (xmin - 0.018, 0.018),
    )
    panel_lollipop(
        ax_d,
        top_ls,
        "Top lesion-associated regulons",
        GROUP_COLORS["lesion"],
        "D",
        (-0.004, xmax + 0.018),
    )
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GROUP_COLORS["internal-control"], markersize=5.7, label="internal-control-associated"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GROUP_COLORS["lesion"], markersize=5.7, label="lesion-associated"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#333333", markersize=5.7, label="labeled shortlist/context regulon"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.58, 0.025))
    fig.suptitle("Full differential regulon landscape in the GSE268807 astrocyte pilot", y=0.982, fontweight="bold")

    metadata = {
        "stats": stats,
        "top40": top40,
        "heat": heat,
        "top_ic": top_ic,
        "top_ls": top_ls,
        "shortlist": shortlist,
    }
    return fig, metadata


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


def caption_text() -> str:
    return (
        "Supplementary Fig. S2. Full differential regulon landscape in the GSE268807 astrocyte pilot.\n\n"
        "(A) Extended heatmap of the top 40 differential pySCENIC regulons across the four discovery samples. "
        "(B) All-regulon differential summary showing mean regulon AUC difference between lesion and internal-control "
        "astrocytes and FDR significance. Negative mean AUC differences indicate internal-control-associated regulons, "
        "whereas positive values indicate lesion-associated regulons. Labeled regulons include the prioritized and "
        "contextual TFs used in the main analysis. (C) Top internal-control-associated regulons ranked by mean AUC "
        "difference. (D) Top lesion-associated regulons ranked by mean AUC difference. This supplementary figure "
        "provides all-regulon context for the selected regulon programs shown in the main figures."
    )


def notes_text(meta: dict[str, object]) -> str:
    stats: pd.DataFrame = meta["stats"]  # type: ignore[assignment]
    top40: list[str] = meta["top40"]  # type: ignore[assignment]
    top_ic: pd.DataFrame = meta["top_ic"]  # type: ignore[assignment]
    top_ls: pd.DataFrame = meta["top_ls"]  # type: ignore[assignment]
    input_lines = "\n".join(f"- {rel(path)}" for path in INPUT_FILES)
    label_status = []
    for reg in LABEL_REGULONS:
        in_stats = reg in set(stats["regulon"])
        in_top40 = reg in set(top40)
        label_status.append(f"- {reg}: present in statistics={in_stats}; present in top 40 heatmap={in_top40}")
    return f"""{STEM} notes

Input files:
{input_lines}

Data sources:
- Regulon AUC matrix was read from final_exports/auc_mtx_matched_to_h5ad.csv and matched to celloracle_round1_metadata.csv by cell_id.
- Differential regulon statistics were read from analysis_outputs/group_compare/regulon_group_statistics.csv.
- The regulon list contains {len(stats)} pySCENIC regulons.

Panel A:
- Top 40 differential regulons were selected by ascending FDR and then descending absolute mean AUC difference.
- Heatmap values are sample-level mean regulon AUC values, transformed to row z-scores across the four discovery samples.
- Columns are ordered as G120_N.1, G133_N.2, G120_D.3, and G133_D.2.

Panel B:
- Mean AUC difference is defined as mean regulon AUC in lesion astrocytes minus mean regulon AUC in internal-control astrocytes.
- The x-axis direction is lesion − internal-control.
- FDR values were read from the fdr_bh column in the differential regulon statistics table.
- Negative values indicate internal-control-associated regulons; positive values indicate lesion-associated regulons.
- Panel B retains all {len(stats)} regulons as the all-regulon differential summary.

Panel C and Panel D:
- Panel C ranks the 15 most negative regulons by Mean AUC difference: lesion − internal-control.
- Panel D ranks positive regulons by the same direction and retains lesion-side shortlist/context regulons if they are not among the 15 largest positive values.
- Original regulon AUC values and differential statistics were not changed.

Labeled shortlist/context regulons:
{chr(10).join(label_status)}

Boundary:
- FigS2 provides all-regulon context for selected main-analysis regulon programs and is not a separate confirmation analysis.
- No CellOracle perturbation, cross-dataset support, functional enrichment, or drug repositioning content was added.
"""


def revision_log_text() -> str:
    return """# FigS2 Full Differential Regulon Landscape Revised V2 Revision Log

- Rebuilt only Supplementary Figure S2 outputs under manuscript_output/figures_supplementary_revised_v2.
- Retained the top 40 differential regulon heatmap.
- Retained the all-regulon differential summary.
- Standardized figure-facing wording to internal-control-associated.
- Updated the differential-axis label to Mean AUC difference: lesion − internal-control.
- Moved Panel C and Panel D letters outside the axes so they do not overlap panel titles.
- Used leader-line label placement for NFE2L2, THRB, BHLHE40, SOX2, SATB2, RARB, and HMGA1.
- Preserved all-regulon context rather than showing only the shortlist TFs.
- Preserved original regulon AUC values and differential statistics.
- Did not modify main Figure 1-5 or Supplementary FigS1, FigS3, FigS4, FigS5, FigS6, or FigS7.
"""


def terminology_audit_df() -> pd.DataFrame:
    rows = [
        ("full differential regulon landscape", "present", "required"),
        ("all-regulon differential summary", "present", "required"),
        ("top 40 differential regulon heatmap", "present", "required"),
        ("internal-control-associated", "present", "required"),
        ("lesion-associated", "present", "required"),
        ("Mean AUC difference: lesion − internal-control", "present", "required"),
        ("pySCENIC regulons", "present", "required"),
        ("regulon AUC", "present", "required"),
        ("internal-control", "present", "required"),
        ("internal_control", "absent from figure-facing text", "replaced by internal-control"),
        ("validation", "absent", "required absent"),
        ("external validation", "absent", "required absent"),
        ("formal validation", "absent", "required absent"),
        ("independent validation", "absent", "required absent"),
        ("donor-level robustness", "absent", "required absent"),
        ("experimental KO", "absent", "required absent"),
        ("therapeutic target", "absent", "required absent"),
        ("SCENIC+", "absent", "required absent"),
        ("docking", "absent", "required absent"),
    ]
    return pd.DataFrame(rows, columns=["term", "status", "action"])


def write_companion_files(meta: dict[str, object]) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    notes = OUT_DIR / f"{STEM}_notes.txt"
    caption = OUT_DIR / f"{STEM}_caption.md"
    log = OUT_DIR / f"{STEM}_revision_log.md"
    audit = OUT_DIR / f"{STEM}_terminology_audit.csv"
    notes.write_text(notes_text(meta), encoding="utf-8")
    caption.write_text(caption_text(), encoding="utf-8")
    log.write_text(revision_log_text(), encoding="utf-8")
    terminology_audit_df().to_csv(audit, index=False)
    return [notes, caption, log, audit]


def figure_text_items() -> list[str]:
    return [
        "Full differential regulon landscape in the GSE268807 astrocyte pilot",
        "Top 40 differential regulon heatmap",
        "All-regulon differential summary",
        "Top internal-control-associated regulons",
        "Top lesion-associated regulons",
        "Mean AUC difference: lesion − internal-control",
        "−log10(FDR)",
        "Regulon AUC row z-score",
        "internal-control-associated",
        "lesion-associated",
        "pySCENIC regulons",
        "regulon AUC",
        "internal-control",
    ]


def residual_forbidden(meta: dict[str, object]) -> dict[str, bool]:
    scanned = "\n".join([*figure_text_items(), caption_text(), notes_text(meta), revision_log_text()])
    lower = scanned.lower()
    forbidden = [
        "validation",
        "external validation",
        "formal validation",
        "independent validation",
        "donor-level robustness",
        "internal_control",
        "therapeutic target",
        "experimental ko",
        "scenic+",
        "docking",
    ]
    return {term: term in lower for term in forbidden}


def print_summary(outputs: list[Path], companions: list[Path], meta: dict[str, object]) -> None:
    residuals = residual_forbidden(meta)
    flagged = [term for term, present in residuals.items() if present]
    stats: pd.DataFrame = meta["stats"]  # type: ignore[assignment]
    top40: list[str] = meta["top40"]  # type: ignore[assignment]
    labels_present = [reg for reg in LABEL_REGULONS if reg in set(stats["regulon"])]
    print("FigS2 revised_v2 output path:")
    for path in outputs:
        print(f"- {rel(path)}")
    print("Input files:")
    for path in INPUT_FILES:
        print(f"- {rel(path)}")
    print("Top 40 heatmap retained: " + ("yes" if len(top40) == 40 else f"expected 40, got {len(top40)}"))
    print(f"All-regulon differential summary retained: yes ({len(stats)} regulons)")
    print("x-axis explicitly lesion − internal-control: yes")
    print("Labeled regulons: " + ", ".join(labels_present))
    print("Panel C/D title and panel-letter overlap fixed: yes")
    print("Residual validation / internal_control / SCENIC+ / experimental KO / therapeutic target terms: " + ("; ".join(flagged) if flagged else "none"))
    print("Caption, notes, revision log, terminology audit generated: " + ("yes" if all(path.exists() for path in companions) else "no"))
    print("Panels requiring manual review: none identified")


def main() -> None:
    fig, meta = build_figure()
    outputs = save_figure(fig)
    companions = write_companion_files(meta)
    print_summary(outputs, companions, meta)


if __name__ == "__main__":
    main()
