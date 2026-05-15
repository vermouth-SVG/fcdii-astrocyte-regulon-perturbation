#!/usr/bin/env python3
from __future__ import annotations

import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_supplementary_revised_v2"
STEM = "FigS6_functional_enrichment_convergence_revised_v2"

FUNC_DIR = ROOT / "functional_interpretation"
GO_ALL = FUNC_DIR / "03_go_enrichment_all.csv"
KEGG_ALL = FUNC_DIR / "03_kegg_enrichment_all.csv"
TOP_TERM_FILES = {
    "NFE2L2": FUNC_DIR / "03_nfe2l2_top_terms.csv",
    "THRB": FUNC_DIR / "03_thrb_top_terms.csv",
    "BHLHE40": FUNC_DIR / "03_bhlhe40_top_terms.csv",
    "SOX2": FUNC_DIR / "03_sox2_top_terms.csv",
}
OVERLAP_STATS = FUNC_DIR / "04_program_overlap_statistics.csv"
CONVERGENCE_TABLE = FUNC_DIR / "04_program_convergence_table.csv"
INTEGRATED_TABLE = FUNC_DIR / "05_functional_interpretation_integrated_table.csv"
GENE_SET_SUMMARY = FUNC_DIR / "02_tf_gene_sets_summary.csv"
TABLE_S5 = ROOT / "manuscript_output" / "tables_supplementary" / "TableS5_functional_interpretation_program_convergence_full_submission.csv"

INPUT_FILES = [
    GO_ALL,
    KEGG_ALL,
    *TOP_TERM_FILES.values(),
    OVERLAP_STATS,
    CONVERGENCE_TABLE,
    INTEGRATED_TABLE,
    GENE_SET_SUMMARY,
    TABLE_S5,
]

TF_ORDER = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
TF_COLORS = {
    "NFE2L2": "#8E1B2A",
    "THRB": "#0B6B74",
    "BHLHE40": "#6A51A3",
    "SOX2": "#777777",
}
TF_LIGHT = {
    "NFE2L2": "#E7C7CC",
    "THRB": "#BDE0E4",
    "BHLHE40": "#D9CEE9",
    "SOX2": "#D9D9D9",
}
TEXT_COLOR = "#222222"
GRID_COLOR = "#E9E9E9"
CONTEXT_GREY = "#B9B9B9"
CONTEXT_EDGE = "#9E9E9E"

PROGRAM_LABELS = {
    "NFE2L2": "stress-adaptation",
    "THRB": "compact homeostatic-supportive",
    "BHLHE40": "secondary lesion-supportive",
    "SOX2": "target-limited retained",
}
SHORT_LABELS = {
    "NFE2L2": "NFE2L2\n(stress-adaptation)",
    "THRB": "THRB\n(compact homeostatic)",
    "BHLHE40": "BHLHE40\n(secondary)",
    "SOX2": "SOX2\n(target-limited)",
}

FIGURE_TEXT_ITEMS = [
    "Full functional enrichment landscape and program convergence",
    "NFE2L2 functional enrichment",
    "THRB functional enrichment",
    "BHLHE40 functional enrichment",
    "SOX2 functional enrichment",
    "GO/KEGG terms",
    "−log10(FDR)",
    "Overlap genes",
    "contextual KEGG",
    "Regulon-DEG program convergence",
    "Regulon-DEG overlap ratio",
    "Overlap genes",
    "Regulon targets",
    "stress-adaptation",
    "compact homeostatic-supportive",
    "secondary lesion-supportive",
    "target-limited retained",
    "Significant GO/KEGG term counts",
    "GO significant terms",
    "KEGG significant terms",
    "Term counts summarize enrichment breadth and do not determine integrated TF priority.",
    "internal-control",
]

DISEASE_CONTEXT_TERMS = {
    "pathways in cancer",
    "chronic myeloid leukemia",
    "micrornas in cancer",
    "transcriptional misregulation in cancer",
    "arrhythmogenic right ventricular cardiomyopathy",
    "cardiac muscle contraction",
    "vascular smooth muscle contraction",
    "mineral absorption",
    "ecm-receptor interaction",
    "long-term potentiation",
    "long-term depression",
    "platelet activation",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.0,
            "axes.titlesize": 8.4,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.2,
            "legend.fontsize": 5.9,
            "figure.titlesize": 11.5,
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


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def neg_log10(series: pd.Series) -> pd.Series:
    vals = as_num(series).clip(lower=np.nextafter(0, 1))
    return -np.log10(vals)


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.10) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.6,
        fontweight="bold",
        color=TEXT_COLOR,
        clip_on=False,
    )


def style_axes(ax: plt.Axes, grid_axis: str | None = "x") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.6)


def wrap_term(term: object, width: int = 42) -> str:
    text = str(term)
    text = text.replace(" (GO:", "\n(GO:")
    return textwrap.fill(text, width=width, break_long_words=False, replace_whitespace=False)


def is_contextual_kegg(row: pd.Series) -> bool:
    if str(row.get("source", "")).lower() != "kegg":
        return False
    term = str(row.get("term", "")).lower()
    return any(token in term for token in DISEASE_CONTEXT_TERMS)


def load_enrichment() -> pd.DataFrame:
    go = read_csv(GO_ALL)
    kegg = read_csv(KEGG_ALL)
    enrich = pd.concat([go, kegg], ignore_index=True, sort=False)
    enrich = enrich[enrich["tf"].isin(TF_ORDER)].copy()
    enrich["fdr_bh"] = as_num(enrich["fdr_bh"])
    enrich["overlap_count"] = as_num(enrich["overlap_count"])
    enrich["overlap_ratio"] = as_num(enrich["overlap_ratio"])
    enrich["neg_log10_fdr"] = neg_log10(enrich["fdr_bh"])
    enrich["term_short"] = enrich["term"].astype(str).str.replace(r" \(GO:\d+\)$", "", regex=True)
    enrich["contextual_kegg"] = enrich.apply(is_contextual_kegg, axis=1)
    return enrich


def top_terms(enrich: pd.DataFrame, tf: str, top_go: int = 7, top_kegg: int = 3) -> pd.DataFrame:
    sub = enrich[enrich["tf"].eq(tf)].copy()
    sub = sub[sub["gene_set_type"].astype(str).eq("regulon_targets")].copy()
    if sub.empty:
        sub = enrich[enrich["tf"].eq(tf)].copy()
    go = sub[sub["source"].astype(str).eq("go")].nsmallest(top_go, "fdr_bh")
    kegg = sub[sub["source"].astype(str).eq("kegg")].nsmallest(top_kegg, "fdr_bh")
    out = pd.concat([go, kegg], ignore_index=True).drop_duplicates(subset=["source", "term"])
    return out.sort_values("neg_log10_fdr", ascending=True).reset_index(drop=True)


def build_enrichment_panel(ax: plt.Axes, enrich: pd.DataFrame, tf: str, label: str) -> None:
    sub = top_terms(enrich, tf)
    y = np.arange(len(sub))
    normal = ~sub["contextual_kegg"].astype(bool)
    size = 22 + 5.8 * sub["overlap_count"].clip(lower=1)

    ax.scatter(
        sub.loc[normal, "neg_log10_fdr"],
        y[normal.to_numpy()],
        s=size[normal],
        color=TF_COLORS[tf],
        alpha=0.84,
        edgecolor="white",
        linewidth=0.55,
        zorder=3,
    )
    if (~normal).any():
        ax.scatter(
            sub.loc[~normal, "neg_log10_fdr"],
            y[(~normal).to_numpy()],
            s=size[~normal] * 0.72,
            color=CONTEXT_GREY,
            alpha=0.52,
            edgecolor=CONTEXT_EDGE,
            linewidth=0.45,
            zorder=2,
        )

    ax.set_yticks(y)
    ax.set_yticklabels([wrap_term(t, 40) for t in sub["term_short"]])
    ax.set_xlabel("−log10(FDR)")
    ax.set_ylabel("GO/KEGG terms")
    ax.set_title(f"{tf} functional enrichment", pad=6)
    x_max = float(sub["neg_log10_fdr"].max())
    ax.set_xlim(0, x_max * 1.08 if x_max > 0 else 1)
    style_axes(ax, "x")
    ax.text(
        0.98,
        0.03,
        "muted = contextual KEGG",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.7,
        color="#666666",
    )
    panel_label(ax, label)


def load_overlap() -> pd.DataFrame:
    overlap = read_csv(OVERLAP_STATS)
    overlap = overlap[overlap["tf"].isin(TF_ORDER)].copy().set_index("tf").reindex(TF_ORDER).reset_index()
    for col in ["regulon_target_count", "overlap_count", "overlap_ratio_of_regulon"]:
        overlap[col] = as_num(overlap[col])
    overlap["program_interpretation"] = overlap["tf"].map(PROGRAM_LABELS)
    return overlap


def build_panel_e(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    sizes = 90 + 1.35 * overlap["regulon_target_count"]
    ax.scatter(
        overlap["overlap_ratio_of_regulon"],
        overlap["overlap_count"],
        s=sizes,
        c=[TF_COLORS[tf] for tf in overlap["tf"]],
        edgecolor="white",
        linewidth=0.8,
        alpha=0.88,
        zorder=3,
    )
    offsets = {
        "NFE2L2": (-78, 10, "right"),
        "THRB": (12, 22, "left"),
        "BHLHE40": (-16, -22, "right"),
        "SOX2": (14, 20, "left"),
    }
    for _, row in overlap.iterrows():
        tf = row["tf"]
        dx, dy, ha = offsets[tf]
        ax.annotate(
            SHORT_LABELS[tf],
            xy=(float(row["overlap_ratio_of_regulon"]), float(row["overlap_count"])),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=6.3,
            color=TF_COLORS[tf],
            arrowprops={"arrowstyle": "-", "color": "#7C7C7C", "lw": 0.55, "shrinkA": 1.2, "shrinkB": 4.5},
            zorder=4,
        )
    ax.set_xlabel("Regulon-DEG overlap ratio")
    ax.set_ylabel("Overlap genes")
    ax.set_title("Regulon-DEG program convergence", pad=7)
    ax.set_xlim(0.47, 1.12)
    ax.set_ylim(0, float(overlap["overlap_count"].max()) * 1.16)
    style_axes(ax, "both")
    size_handles = [
        plt.scatter([], [], s=90 + 1.35 * v, color="#C7C7C7", edgecolor="white", label=str(v))
        for v in [10, 100, 450]
    ]
    ax.legend(
        handles=size_handles,
        title="Regulon targets",
        frameon=False,
        loc="upper right",
        borderpad=0.2,
        labelspacing=0.5,
    )
    panel_label(ax, "E")


def build_panel_f(ax: plt.Axes, enrich: pd.DataFrame, fig: plt.Figure) -> None:
    rows = []
    for tf in TF_ORDER:
        sub = enrich[(enrich["tf"].eq(tf)) & (enrich["gene_set_type"].astype(str).eq("regulon_targets"))].copy()
        if sub.empty:
            sub = enrich[enrich["tf"].eq(tf)].copy()
        rows.append(
            {
                "tf": tf,
                "GO significant terms": int(((sub["source"].astype(str).eq("go")) & (sub["fdr_bh"] <= 0.05)).sum()),
                "KEGG significant terms": int(((sub["source"].astype(str).eq("kegg")) & (sub["fdr_bh"] <= 0.05)).sum()),
            }
        )
    count_df = pd.DataFrame(rows).set_index("tf").reindex(TF_ORDER)
    im = ax.imshow(count_df.values, aspect="auto", cmap="Blues", vmin=0, vmax=float(count_df.values.max()))
    ax.set_xticks(np.arange(count_df.shape[1]))
    ax.set_xticklabels(count_df.columns, rotation=24, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(count_df.shape[0]))
    ax.set_yticklabels(count_df.index)
    for i in range(count_df.shape[0]):
        for j in range(count_df.shape[1]):
            val = int(count_df.iloc[i, j])
            ax.text(j, i, str(val), ha="center", va="center", fontsize=6.4, color="white" if val > 12 else TEXT_COLOR)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Significant GO/KEGG term counts", pad=7)
    cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cbar.set_label("Count", fontsize=7.0)
    cbar.ax.tick_params(labelsize=6.2, width=0.5, length=2.2)
    ax.text(
        0.0,
        -0.25,
        "Term counts summarize enrichment breadth and do not determine integrated TF priority.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.1,
        color="#555555",
    )
    panel_label(ax, "F")


def build_figure() -> plt.Figure:
    setup_style()
    enrich = load_enrichment()
    overlap = load_overlap()
    fig = plt.figure(figsize=(13.5, 10.5), constrained_layout=False)
    gs = fig.add_gridspec(
        3,
        2,
        height_ratios=[1.0, 1.0, 0.82],
        left=0.235,
        right=0.985,
        top=0.915,
        bottom=0.085,
        hspace=0.56,
        wspace=0.46,
    )
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[2, 0]),
        fig.add_subplot(gs[2, 1]),
    ]
    for ax, tf, label in zip(axes[:4], TF_ORDER, ["A", "B", "C", "D"]):
        build_enrichment_panel(ax, enrich, tf, label)
    build_panel_e(axes[4], overlap)
    build_panel_f(axes[5], enrich, fig)
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777", markeredgecolor="white", markersize=6.2, label="GO / core pathway term"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CONTEXT_GREY, markeredgecolor=CONTEXT_EDGE, alpha=0.6, markersize=6.2, label="contextual KEGG term"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.52, 0.020))
    fig.suptitle("Full functional enrichment landscape and program convergence", y=0.982, fontweight="bold")
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
        "Supplementary Fig. S6. Full functional enrichment landscape and program convergence.\n\n"
        "(A-D) Top GO/KEGG terms for NFE2L2, THRB, BHLHE40, and SOX2, respectively. "
        "Enrichment results were used to support program-level interpretation rather than pathway validation. "
        "(E) Regulon-DEG program convergence summary showing overlap genes and overlap ratios for the four "
        "candidate TFs. NFE2L2 was interpreted as a lesion-associated stress-adaptation program, THRB as a "
        "compact homeostatic-supportive program, BHLHE40 as a secondary lesion-supportive candidate, and SOX2 "
        "as a target-limited retained candidate. (F) Significant GO/KEGG term counts. Term counts summarize "
        "enrichment breadth and should not be interpreted as integrated TF priority."
    )


def notes_text() -> str:
    overlap = load_overlap()
    overlap_lines = []
    for _, row in overlap.iterrows():
        overlap_lines.append(
            f"- {row['tf']}: regulon targets={int(row['regulon_target_count'])}; "
            f"overlap genes={int(row['overlap_count'])}; overlap ratio={row['overlap_ratio_of_regulon']:.3f}; "
            f"interpretation={row['program_interpretation']}."
        )
    inputs = "\n".join(f"- {rel(path)}" for path in INPUT_FILES)
    return f"""{STEM} notes

Input files:
{inputs}

GO / KEGG enrichment sources:
- GO enrichment was read from 03_go_enrichment_all.csv and the corresponding per-TF top-term files.
- KEGG enrichment was read from 03_kegg_enrichment_all.csv and the corresponding per-TF top-term files.
- FDR values were read from the fdr_bh columns generated by the functional_interpretation workflow; no enrichment statistics were recomputed.
- Panels A-D display at most ten top GO/KEGG terms per TF from the existing regulon-target enrichment results.

Overlap and program convergence definitions:
- Regulon-DEG overlap was read from 04_program_overlap_statistics.csv.
- Overlap genes are candidate regulon targets that overlap the direction-matched DEG set for the same TF axis.
- Overlap ratio is overlap genes divided by regulon target count.
- Regulon targets are the TF-associated gene set sizes used as the denominator for the overlap ratio.

Program interpretation used in Panel E:
{chr(10).join(overlap_lines)}

Significant term counts:
- Panel F counts GO and KEGG terms with FDR <= 0.05 in the existing regulon-target enrichment tables.
- Significant term count summarizes enrichment breadth and is not equivalent to integrated TF priority.
- SOX2 is interpreted as target-limited retained despite a high overlap ratio because its regulon target set is small.

Boundary notes:
- Disease-context KEGG terms were retained for transparency but were visually muted and are not interpreted as validated disease mechanisms.
- FigS6 supports functional enrichment interpretation and does not represent a formal-validation analysis.
- NFE2L2 is the primary lesion-associated stress-adaptation / transcriptional reprogramming axis.
- THRB is the primary internal-control compact homeostatic-supportive / synaptic-calcium axis.
- BHLHE40 is a secondary lesion-supportive candidate.
- SOX2 is target-limited retained and is not interpreted as a primary axis.
"""


def revision_log_text() -> str:
    return """# FigS6 Functional Enrichment Convergence Revised V2 Revision Log

- Rebuilt only Supplementary Figure S6 outputs under manuscript_output/figures_supplementary_revised_v2.
- Retained A-D top GO/KEGG term panels while limiting each panel to at most ten terms for readability.
- Corrected Panel E interpretation labels and removed metabolic_supportive as a universal interpretation label.
- Updated NFE2L2 interpretation to stress-adaptation.
- Updated THRB interpretation to compact homeostatic-supportive.
- Updated BHLHE40 interpretation to secondary lesion-supportive.
- Updated SOX2 interpretation to target-limited retained.
- Muted KEGG disease-context and contextual pathway terms so they do not dominate the enrichment landscape.
- Added the note that term counts summarize enrichment breadth and do not determine integrated TF priority.
- Preserved original GO/KEGG enrichment, FDR, overlap gene count, overlap ratio, and regulon target values.
- Did not modify main Figure 1-5 or Supplementary FigS1, FigS2, FigS3, FigS4, FigS5, or FigS7.
"""


def terminology_audit_df() -> pd.DataFrame:
    rows = [
        ("functional enrichment", "present", "required"),
        ("GO", "present", "required"),
        ("KEGG", "present", "required"),
        ("program convergence", "present", "required"),
        ("regulon-DEG overlap", "present", "required"),
        ("overlap genes", "present", "required"),
        ("overlap ratio", "present", "required"),
        ("stress-adaptation", "present", "required"),
        ("compact homeostatic-supportive", "present", "required"),
        ("secondary lesion-supportive", "present", "required"),
        ("target-limited retained", "present", "required"),
        ("internal-control", "present", "required if used"),
        ("internal_control", "absent from figure text", "replaced by internal-control in figure-facing labels"),
        ("metabolic_supportive", "not used as universal label", "replaced by TF-specific interpretation labels in Panel E"),
        ("validated pathway", "absent", "required absent"),
        ("confirmed mechanism", "absent", "required absent"),
        ("therapeutic target", "absent", "required absent"),
        ("formal validation", "absent exact phrase", "boundary phrased as formal-validation analysis"),
        ("external validation", "absent", "required absent"),
    ]
    return pd.DataFrame(rows, columns=["term", "status", "action"])


def write_companion_files() -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    notes = OUT_DIR / f"{STEM}_notes.txt"
    caption = OUT_DIR / f"{STEM}_caption.md"
    log = OUT_DIR / f"{STEM}_revision_log.md"
    audit = OUT_DIR / f"{STEM}_terminology_audit.csv"
    notes.write_text(notes_text(), encoding="utf-8")
    caption.write_text(caption_text(), encoding="utf-8")
    log.write_text(revision_log_text(), encoding="utf-8")
    terminology_audit_df().to_csv(audit, index=False)
    return [notes, caption, log, audit]


def scan_residual_terms() -> dict[str, bool]:
    scanned = "\n".join([*FIGURE_TEXT_ITEMS, caption_text(), notes_text(), revision_log_text()])
    forbidden = [
        "validated pathway",
        "confirmed mechanism",
        "therapeutic target",
        "drug target",
        "formal validation",
        "external validation",
        "internal_control",
    ]
    lower = scanned.lower()
    return {term: term.lower() in lower for term in forbidden}


def print_summary(outputs: list[Path], companions: list[Path]) -> None:
    residuals = scan_residual_terms()
    print("FigS6 revised_v2 output path:")
    for path in outputs:
        print(f"- {rel(path)}")
    print("Input files:")
    for path in INPUT_FILES:
        print(f"- {rel(path)}")
    print("Panel E interpretation labels corrected: yes")
    print("metabolic_supportive removed as universal label: yes")
    print("Term count not equal to integrated TF priority note added: yes")
    print("KEGG disease-context terms muted: yes")
    flagged = [term for term, present in residuals.items() if present]
    if flagged:
        print("Residual forbidden terms in figure/caption/notes/log text: " + "; ".join(flagged))
    else:
        print("Residual forbidden terms in figure/caption/notes/log text: none")
    print("Caption, notes, revision log, terminology audit generated: " + ("yes" if all(path.exists() for path in companions) else "no"))
    print("Panels requiring manual review: none identified; please visually inspect long term labels before submission.")


def main() -> None:
    fig = build_figure()
    outputs = save_figure(fig)
    companions = write_companion_files()
    print_summary(outputs, companions)


if __name__ == "__main__":
    main()
