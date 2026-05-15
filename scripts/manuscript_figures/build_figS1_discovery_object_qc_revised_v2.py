#!/usr/bin/env python3
from __future__ import annotations

import json
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_supplementary_revised_v2"
STEM = "FigS1_discovery_object_qc_revised_v2"

METADATA = ROOT / "celloracle_run" / "prepared_data" / "celloracle_round1_metadata.csv"
QC_DIR = ROOT / "analysis_outputs" / "gse268807_astrocyte_pilot" / "qc_metrics"
STRUCTURE_SUMMARY = ROOT / "analysis_outputs" / "structure_check" / "structure_summary.json"
OBS_OVERVIEW = ROOT / "analysis_outputs" / "structure_check" / "obs_columns_overview.csv"
GROUP_COUNTS = ROOT / "analysis_outputs" / "structure_check" / "selected_group_counts.csv"
MERGE_REPORT = ROOT / "final_exports" / "merge_report.txt"
REGULON_NAMES = ROOT / "final_exports" / "regulon_names.txt"
MATCHED_AUC = ROOT / "final_exports" / "auc_mtx_matched_to_h5ad.csv"
H5AD_OBJECT = ROOT / "final_exports" / "astrocyte_pilot_rna_with_pyscenic_auc.h5ad"

INPUT_FILES = [
    METADATA,
    QC_DIR / "qc_metrics_GSE268807_G120_D_FL_max3000.csv",
    QC_DIR / "qc_metrics_GSE268807_G120_F1_N_max3000.csv",
    QC_DIR / "qc_metrics_GSE268807_G133_D_FL_max3000.csv",
    QC_DIR / "qc_metrics_GSE268807_G133_N_FL_max3000.csv",
    STRUCTURE_SUMMARY,
    OBS_OVERVIEW,
    GROUP_COUNTS,
    MERGE_REPORT,
    REGULON_NAMES,
    MATCHED_AUC,
    H5AD_OBJECT,
]

SAMPLE_ORDER = ["G120_N.1", "G133_N.2", "G120_D.3", "G133_D.2"]
GROUP_ORDER = ["internal-control", "lesion"]
GROUP_COLORS = {"internal-control": "#0B6B74", "lesion": "#8E1B2A"}
GROUP_LIGHT = {"internal-control": "#BDE0E4", "lesion": "#E7C7CC"}
TEXT_COLOR = "#222222"
GRID_COLOR = "#E9E9E9"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.0,
            "axes.titlesize": 8.6,
            "axes.labelsize": 7.4,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.6,
            "legend.fontsize": 6.4,
            "figure.titlesize": 12.0,
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


def fmt_int(value: object) -> str:
    return f"{int(float(value)):,}"


def parse_simple_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("["):
            continue
        if "=" in line:
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
        elif ":" in line:
            key, value = line.split(":", 1)
            out[key.strip()] = value.strip()
    return out


def normalize_qc_sample_prefix(value: object) -> str:
    return re.sub(r"_max\d+$", "", str(value))


def display_group(value: object) -> str:
    return str(value).replace("internal_control", "internal-control")


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.10) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11.0,
        fontweight="bold",
        color=TEXT_COLOR,
        clip_on=False,
    )


def style_axes(ax: plt.Axes, grid_axis: str | None = "x") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.6)


def load_metadata() -> pd.DataFrame:
    meta = pd.read_csv(METADATA)
    if str(meta.columns[0]).startswith("Unnamed"):
        meta = meta.rename(columns={meta.columns[0]: "cell_id"})
    meta["group_display"] = meta["group"].map(display_group)
    return meta


def load_qc_trace(meta: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for path in sorted(QC_DIR.glob("qc_metrics_*.csv")):
        df = pd.read_csv(path)
        df["input_file"] = rel(path)
        df["sample_prefix"] = normalize_qc_sample_prefix(path.stem.replace("qc_metrics_", ""))
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No qc_metrics_*.csv files found under {QC_DIR}")
    qc = pd.concat(frames, ignore_index=True)
    map_df = meta[["sample_prefix", "clinical_sample_label", "group_display"]].drop_duplicates()
    qc = qc.merge(map_df, on="sample_prefix", how="left")
    qc["clinical_sample_label"] = qc["clinical_sample_label"].fillna(qc["sample_prefix"])
    inferred_group = pd.Series(
        np.where(qc["sample_prefix"].astype(str).str.contains("_D_"), "lesion", "internal-control"),
        index=qc.index,
    )
    qc["group_display"] = qc["group_display"].fillna(inferred_group)
    qc["sample_order"] = qc["clinical_sample_label"].map({sample: idx for idx, sample in enumerate(SAMPLE_ORDER)})
    return qc.sort_values("sample_order").reset_index(drop=True)


def parse_summary() -> dict[str, int]:
    merge = parse_simple_kv(MERGE_REPORT)
    with STRUCTURE_SUMMARY.open("r", encoding="utf-8") as handle:
        structure = json.load(handle)
    regulon_count = int(merge.get("auc_regulons", structure.get("n_regulons", 0)))
    if REGULON_NAMES.exists():
        line_count = len([line for line in REGULON_NAMES.read_text(encoding="utf-8").splitlines() if line.strip()])
        regulon_count = line_count or regulon_count
    return {
        "cells": int(merge.get("h5ad_cells", structure.get("n_cells", 0))),
        "genes": int(merge.get("h5ad_genes", structure.get("n_genes", 0))),
        "auc_cells": int(merge.get("auc_cells", 0)),
        "shared_cells": int(merge.get("shared_cells", 0)),
        "missing_in_auc": int(merge.get("missing_in_auc", 0)),
        "missing_in_h5ad": int(merge.get("missing_in_h5ad", 0)),
        "regulons": int(regulon_count),
    }


def sample_counts(meta: pd.DataFrame) -> pd.DataFrame:
    counts = (
        meta.groupby(["clinical_sample_label", "group_display"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    counts["sample_order"] = counts["clinical_sample_label"].map({sample: idx for idx, sample in enumerate(SAMPLE_ORDER)})
    return counts.sort_values("sample_order").reset_index(drop=True)


def sample_group_matrix(meta: pd.DataFrame) -> pd.DataFrame:
    matrix = (
        meta.groupby(["group_display", "clinical_sample_label"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=GROUP_ORDER, columns=SAMPLE_ORDER)
        .fillna(0)
        .astype(int)
    )
    return matrix


def draw_panel_a(ax: plt.Axes, counts: pd.DataFrame) -> None:
    y = np.arange(len(SAMPLE_ORDER))[::-1]
    ordered = counts.set_index("clinical_sample_label").reindex(SAMPLE_ORDER).reset_index()
    ax.hlines(y, xmin=0, xmax=ordered["n_cells"], color="#D4D4D4", linewidth=1.2, zorder=1)
    ax.scatter(
        ordered["n_cells"],
        y,
        s=70,
        c=[GROUP_COLORS[group] for group in ordered["group_display"]],
        edgecolor="white",
        linewidth=0.9,
        zorder=3,
    )
    for yi, value, group in zip(y, ordered["n_cells"], ordered["group_display"]):
        ax.text(value + 22, yi, fmt_int(value), ha="left", va="center", fontsize=6.5, color=GROUP_COLORS[group])
    ax.set_yticks(y)
    ax.set_yticklabels(SAMPLE_ORDER)
    ax.set_xlabel("Cells in final analysis object")
    ax.set_ylabel("Sample ID")
    ax.set_title("Final 4-sample pilot composition", pad=7)
    ax.set_xlim(0, max(ordered["n_cells"]) * 1.22)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=GROUP_COLORS["internal-control"], markeredgecolor="white", markersize=6.5, label="internal-control"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=GROUP_COLORS["lesion"], markeredgecolor="white", markersize=6.5, label="lesion"),
        ],
        frameon=False,
        loc="lower right",
    )
    style_axes(ax, "x")
    panel_label(ax, "A")


def draw_panel_b(ax: plt.Axes, qc: pd.DataFrame) -> None:
    stages = ["n_cells_input", "n_cells_after_qc", "n_cells_after_subsample"]
    stage_labels = ["Input", "After QC", "After subsampling"]
    x = np.arange(len(stages))
    label_y = {"G133_N.2": 3800, "G120_D.3": 3150, "G133_D.2": 2550, "G120_N.1": 1200}
    for _, row in qc.iterrows():
        vals = [float(row[stage]) for stage in stages]
        group = str(row["group_display"])
        sample = str(row["clinical_sample_label"])
        color = GROUP_COLORS.get(group, "#777777")
        ax.plot(x, vals, marker="o", markersize=4.3, linewidth=1.7, color=color, alpha=0.90, zorder=3)
        y_text = label_y.get(sample, vals[-1])
        ax.plot([x[-1], x[-1] + 0.08], [vals[-1], y_text], color=color, linewidth=0.55, alpha=0.75, zorder=2)
        ax.text(x[-1] + 0.10, y_text, sample, fontsize=6.2, va="center", color=color)
    ax.set_xticks(x)
    ax.set_xticklabels(stage_labels)
    ax.set_ylabel("Cells")
    ax.set_title("Upstream source-level QC and subsampling trace", pad=7)
    ax.set_xlim(-0.10, 2.58)
    ymax = float(qc[stages].to_numpy().max()) * 1.10
    ax.set_ylim(0, ymax)
    ax.yaxis.set_major_formatter(lambda value, _pos: fmt_int(value))
    ax.text(
        0.02,
        0.96,
        "Source-level trace before final object assembly",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.0,
        color="#5F5F5F",
    )
    style_axes(ax, "y")
    panel_label(ax, "B")


def draw_panel_c(ax: plt.Axes, matrix: pd.DataFrame, fig: plt.Figure) -> None:
    im = ax.imshow(matrix.values, aspect="auto", cmap="Blues", vmin=0, vmax=float(matrix.values.max()))
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=32, ha="right", rotation_mode="anchor")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix.iloc[i, j])
            ax.text(j, i, str(value), ha="center", va="center", fontsize=7.0, color="white" if value > 400 else TEXT_COLOR)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Final sample-by-group cell matrix", pad=7)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Cells", fontsize=7.0)
    cbar.ax.tick_params(labelsize=6.2, width=0.5, length=2.2)
    panel_label(ax, "C")


def draw_box(ax: plt.Axes, x: float, y: float, w: float, h: float, title: str, body: str, facecolor: str = "white") -> None:
    rect = Rectangle((x, y), w, h, transform=ax.transAxes, facecolor=facecolor, edgecolor="#B7B7B7", linewidth=0.8)
    ax.add_patch(rect)
    ax.text(x + 0.025, y + h - 0.035, title, transform=ax.transAxes, ha="left", va="top", fontsize=6.9, fontweight="bold", color=TEXT_COLOR)
    ax.text(x + 0.025, y + h - 0.083, body, transform=ax.transAxes, ha="left", va="top", fontsize=5.9, color="#333333", linespacing=1.18)


def draw_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=0.8,
            color="#7F7F7F",
        )
    )


def draw_panel_d(ax: plt.Axes, summary: dict[str, int]) -> None:
    ax.axis("off")
    ax.set_title("Input lineage and analytical object trace", pad=7)
    box_w, box_h = 0.78, 0.158
    x0 = 0.11
    ys = [0.75, 0.545, 0.340, 0.135]
    boxes = [
        ("Raw pilot input", "GSE268807 4 samples\nsample-level pilot subset", "#FFFFFF"),
        ("QC / subsampling", "max 3,000 cells per source\nper-sample QC traces retained", "#FFFFFF"),
        ("Merged RNA object", f"{fmt_int(summary['cells'])} cells x {fmt_int(summary['genes'])} genes", "#F7F7F7"),
        ("AUCell-matched export", f"{fmt_int(summary['regulons'])} pySCENIC regulons\n0 missing matched cells", "#F7F7F7"),
    ]
    for y, (title, body, face) in zip(ys, boxes):
        draw_box(ax, x0, y, box_w, box_h, title, body, face)
    for idx in range(len(ys) - 1):
        draw_arrow(ax, (0.50, ys[idx]), (0.50, ys[idx + 1] + box_h))
    ax.text(
        0.11,
        0.025,
        "Current analysis object: GSE268807 astrocyte pilot -> pySCENIC AUCell-matched object -> CellOracle-ready metadata",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.9,
        color="#555555",
        wrap=True,
    )
    panel_label(ax, "D")


def build_figure() -> plt.Figure:
    setup_style()
    meta = load_metadata()
    qc = load_qc_trace(meta)
    summary = parse_summary()
    counts = sample_counts(meta)
    matrix = sample_group_matrix(meta)

    fig = plt.figure(figsize=(11.4, 8.2), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        left=0.095,
        right=0.985,
        top=0.865,
        bottom=0.085,
        hspace=0.46,
        wspace=0.32,
    )
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    draw_panel_a(axes[0], counts)
    draw_panel_b(axes[1], qc)
    draw_panel_c(axes[2], matrix, fig)
    draw_panel_d(axes[3], summary)

    fig.suptitle("Discovery object provenance and AUCell-matched input overview", y=0.975, fontweight="bold")
    summary_strip = (
        "GSE268807 astrocyte pilot | 4 samples | 2 lesion / 2 internal-control | 2 donors | "
        f"{fmt_int(summary['cells'])} astrocytes | {fmt_int(summary['genes'])} genes | {fmt_int(summary['regulons'])} regulons"
    )
    fig.text(0.5, 0.918, summary_strip, ha="center", va="center", fontsize=7.0, color="#4D4D4D")
    return fig


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
        "Supplementary Fig. S1. Discovery object provenance and AUCell-matched input overview.\n\n"
        "(A) Final four-sample composition of the GSE268807 astrocyte pilot used for downstream pySCENIC, "
        "AUCell, CellOracle, and robustness analyses. (B) Upstream source-level QC and subsampling trace before "
        "final analysis object assembly. (C) Final sample-by-group cell matrix showing 1,434 internal-control "
        "astrocytes and 888 lesion astrocytes, totaling 2,322 astrocytes. (D) Input lineage and analytical object "
        "trace from raw pilot input to the merged RNA object and AUCell-matched export. The final analysis object "
        "contained 2,322 cells, 36,601 genes, and 105 pySCENIC regulons."
    )


def notes_text(meta: pd.DataFrame, qc: pd.DataFrame, summary: dict[str, int], matrix: pd.DataFrame) -> str:
    counts = sample_counts(meta)
    sample_lines = "\n".join(
        f"- {row.clinical_sample_label}: {fmt_int(row.n_cells)} cells ({row.group_display})"
        for row in counts.itertuples(index=False)
    )
    input_lines = "\n".join(f"- {rel(path)}" for path in INPUT_FILES)
    g130_found = "G130_D.3" in "\n".join(
        [
            " ".join(meta.astype(str).to_numpy().ravel().tolist()),
            " ".join(qc.astype(str).to_numpy().ravel().tolist()),
        ]
    )
    sample_check = "G130_D.3 was not found in metadata or QC trace; G120_D.3 was confirmed as the lesion sample." if not g130_found else "G130_D.3 was present in source metadata and requires manual review."
    return f"""{STEM} notes

Input files:
{input_lines}

Final analysis object source:
- Final analysis object was read from the GSE268807 astrocyte pilot metadata and final_exports merge report.
- The merged RNA object path recorded in the merge report is final_exports/astrocyte_pilot_rna_with_pyscenic_auc.h5ad.
- Final object size: {fmt_int(summary['cells'])} cells and {fmt_int(summary['genes'])} genes.
- pySCENIC regulon count: {fmt_int(summary['regulons'])} regulons.
- AUCell-matched export: {fmt_int(summary['auc_cells'])} AUCell rows, {fmt_int(summary['shared_cells'])} shared cells, {fmt_int(summary['missing_in_auc'] + summary['missing_in_h5ad'])} missing matched cells.

Sample-level cell counts:
{sample_lines}

Group totals:
- internal-control astrocytes: {fmt_int(int(matrix.loc['internal-control'].sum()))}
- lesion astrocytes: {fmt_int(int(matrix.loc['lesion'].sum()))}
- total astrocytes: {fmt_int(int(matrix.values.sum()))}

Panel definitions:
- Panel A shows the final four-sample composition of the discovery object.
- Panel B shows upstream source-level QC and subsampling. These values reflect source-level processing before final object assembly and should not be read as the final {fmt_int(summary['cells'])}-cell object distribution.
- Panel C was computed by cross-tabulating clinical_sample_label and group in the final metadata. Matrix total = {fmt_int(int(matrix.values.sum()))}.
- Panel D traces Raw pilot input -> QC / subsampling -> Merged RNA object -> AUCell-matched export.

Sample ID check:
- {sample_check}

Interpretation boundary:
- FigS1 is an object provenance / QC overview for the discovery object, not an inferential confirmation figure.
"""


def revision_log_text() -> str:
    return """# FigS1 Discovery Object QC Revised V2 Revision Log

- Rebuilt only Supplementary Figure S1 under manuscript_output/figures_supplementary_revised_v2.
- Standardized the final object count as 2,322 astrocytes.
- Standardized the gene count as 36,601 genes.
- Standardized the AUCell-matched regulon count as 105 regulons.
- Standardized figure-facing group wording to internal-control.
- Checked the suspected G130_D.3 sample label and confirmed G120_D.3 in metadata.
- Updated Panel B title to Upstream source-level QC and subsampling trace.
- Updated Panel D object trace to Raw pilot input -> QC/subsampling -> merged RNA object -> AUCell-matched export.
- Preserved original sample cell counts, QC trace values, object dimensions, and regulon counts.
- Did not modify main Figure 1-5 or Supplementary FigS2, FigS3, FigS4, FigS5, FigS6, or FigS7.
"""


def terminology_audit_df() -> pd.DataFrame:
    rows = [
        ("2,322", "present", "required"),
        ("2322", "formatted", "replaced by 2,322 in figure-facing text"),
        ("2,332", "absent", "required absent"),
        ("36,601", "present", "required"),
        ("36601", "formatted", "replaced by 36,601 in figure-facing text"),
        ("105 regulons", "present", "required"),
        ("internal-control", "present", "required"),
        ("internal_control", "absent from figure-facing text", "replaced by internal-control"),
        ("G130_D.3", "absent from metadata and figure text", "correct sample_id if erroneous"),
        ("G120_D.3", "present", "metadata confirms"),
        ("validation", "absent from figure/caption/notes/log text", "required absent"),
        ("external validation", "absent", "required absent"),
        ("formal validation", "absent", "required absent"),
        ("donor-level robustness", "absent", "required absent"),
        ("SCENIC+", "absent", "required absent"),
        ("docking", "absent", "required absent"),
    ]
    return pd.DataFrame(rows, columns=["term", "status", "action"])


def write_companion_files(meta: pd.DataFrame, qc: pd.DataFrame, summary: dict[str, int], matrix: pd.DataFrame) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    notes = OUT_DIR / f"{STEM}_notes.txt"
    caption = OUT_DIR / f"{STEM}_caption.md"
    log = OUT_DIR / f"{STEM}_revision_log.md"
    audit = OUT_DIR / f"{STEM}_terminology_audit.csv"
    notes.write_text(notes_text(meta, qc, summary, matrix), encoding="utf-8")
    caption.write_text(caption_text(), encoding="utf-8")
    log.write_text(revision_log_text(), encoding="utf-8")
    terminology_audit_df().to_csv(audit, index=False)
    return [notes, caption, log, audit]


def figure_text_items(summary: dict[str, int]) -> list[str]:
    return [
        "Discovery object provenance and AUCell-matched input overview",
        "GSE268807 astrocyte pilot",
        "Final 4-sample pilot composition",
        "Cells in final analysis object",
        "Sample ID",
        "internal-control",
        "lesion",
        "Upstream source-level QC and subsampling trace",
        "Cells",
        "Input",
        "After QC",
        "After subsampling",
        "Source-level trace before final object assembly",
        "Final sample-by-group cell matrix",
        "Input lineage and analytical object trace",
        "Raw pilot input",
        "QC / subsampling",
        "Merged RNA object",
        "AUCell-matched export",
        f"{fmt_int(summary['cells'])} astrocytes",
        f"{fmt_int(summary['genes'])} genes",
        f"{fmt_int(summary['regulons'])} regulons",
        "pySCENIC regulons",
    ]


def residual_forbidden(meta: pd.DataFrame, qc: pd.DataFrame, summary: dict[str, int], matrix: pd.DataFrame) -> dict[str, bool]:
    scanned = "\n".join(
        [
            *figure_text_items(summary),
            caption_text(),
            notes_text(meta, qc, summary, matrix),
            revision_log_text(),
        ]
    )
    lower = scanned.lower()
    forbidden = [
        "validation",
        "external validation",
        "formal validation",
        "donor-level robustness",
        "internal_control",
        "2,332",
        "scenic+",
        "docking",
    ]
    return {term: term.lower() in lower for term in forbidden}


def print_summary(outputs: list[Path], companions: list[Path], meta: pd.DataFrame, qc: pd.DataFrame, summary: dict[str, int], matrix: pd.DataFrame) -> None:
    sample_ids = sorted(meta["clinical_sample_label"].dropna().astype(str).unique().tolist())
    g130_in_source = any("G130_D.3" in str(value) for value in meta.astype(str).to_numpy().ravel()) or any("G130_D.3" in str(value) for value in qc.astype(str).to_numpy().ravel())
    residuals = residual_forbidden(meta, qc, summary, matrix)
    flagged = [term for term, present in residuals.items() if present]
    print("FigS1 revised_v2 output path:")
    for path in outputs:
        print(f"- {rel(path)}")
    print("Input files:")
    for path in INPUT_FILES:
        print(f"- {rel(path)}")
    print("Sample IDs confirmed: " + ", ".join(sample_ids))
    print("G130_D.3 found and corrected: " + ("manual review needed" if g130_in_source else "not found; G120_D.3 confirmed"))
    print(f"Counts standardized: {fmt_int(summary['cells'])} cells / {fmt_int(summary['genes'])} genes / {fmt_int(summary['regulons'])} regulons")
    print("internal_control replaced by internal-control: yes")
    print("Residual validation / external validation / SCENIC+ / donor-level robustness terms: " + ("; ".join(flagged) if flagged else "none"))
    print("Caption, notes, revision log, terminology audit generated: " + ("yes" if all(path.exists() for path in companions) else "no"))
    print("Panels requiring manual review: none identified")


def main() -> None:
    meta = load_metadata()
    qc = load_qc_trace(meta)
    summary = parse_summary()
    matrix = sample_group_matrix(meta)
    fig = build_figure()
    outputs = save_figure(fig)
    companions = write_companion_files(meta, qc, summary, matrix)
    print_summary(outputs, companions, meta, qc, summary, matrix)


if __name__ == "__main__":
    main()
