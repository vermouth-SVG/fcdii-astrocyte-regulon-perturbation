#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_supplementary_revised_v2"
STEM = "FigS7_drug_repositioning_exploratory_revised_v2"

INPUT_FILES = [
    ROOT / "drug_repositioning" / "16_paper_ready_table_after_bbb_v21.csv",
    ROOT / "drug_repositioning" / "13_integrated_after_bbb_v21.csv",
    ROOT / "drug_repositioning" / "03_integrated_rescored_v21.csv",
    ROOT / "drug_repositioning" / "11_swissadme_parsed_and_matched_v21.csv",
    ROOT / "drug_repositioning" / "12_bbb_annotation_manual_review_v21.csv",
]

LAYER_ORDER = [
    "headline_cns_mechanism_direction_leads",
    "supportive_after_bbb_leads",
    "peripheral_program_modulating_clues",
]
LAYER_LABELS = {
    "headline_cns_mechanism_direction_leads": "mechanism-direction clues",
    "supportive_after_bbb_leads": "supportive clues",
    "peripheral_program_modulating_clues": "peripheral/program-modulating clues",
}
LAYER_COLORS = {
    "headline_cns_mechanism_direction_leads": "#9B4A55",
    "supportive_after_bbb_leads": "#B98545",
    "peripheral_program_modulating_clues": "#7C828A",
    "retained_full_list_after_bbb": "#D1D5D8",
    "low_priority_not_retained_v21": "#E5E7E9",
}
AXIS_ORDER = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]

DISCLAIMER = "Exploratory drug-signature enrichment only; not therapeutic recommendation."
PANEL_E_FINAL = "No compound is interpreted as a validated therapy or treatment recommendation."

FIGURE_TEXT_ITEMS = [
    "Exploratory drug-signature enrichment and BBB-layered prioritization",
    DISCLAIMER,
    "BBB-layered exploratory clue matrix",
    "Count",
    *LAYER_LABELS.values(),
    *AXIS_ORDER,
    "Exploratory after-BBB ranked drug-signature clues",
    "After-BBB exploratory score",
    "Pre-BBB versus after-BBB exploratory score shift",
    "Pre-BBB exploratory score",
    "After-BBB exploratory score",
    "Mechanism-theme distribution across exploratory layers",
    "Conservative interpretation of exploratory drug-signature clues",
    "enrichment-derived, BBB-reviewed, hypothesis-generating only",
    "secondary/supportive signal; requires manual review",
    "program-linked but not prioritized for CNS-directed interpretation",
    PANEL_E_FINAL,
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.titlesize": 8.4,
            "axes.labelsize": 7.3,
            "xtick.labelsize": 6.2,
            "ytick.labelsize": 6.2,
            "legend.fontsize": 6.0,
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


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def panel_label(ax: plt.Axes, label: str, x: float = -0.11, y: float = 1.12) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=10.6,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def style_axes(ax: plt.Axes, grid_axis: str | None = "y") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color="#E8E8E8", linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def short_compound_name(name: object, max_len: int = 31) -> str:
    text = str(name)
    if len(text) <= max_len:
        return text
    return text[: max_len - 1] + "."


def display_theme(value: object) -> str:
    text = str(value).strip() if pd.notna(value) else "unassigned"
    mapping = {
        "未自动识别；需人工核实": "manual review needed",
        "未自动识别；补充轴线索": "secondary-axis clue",
        "environmental_toxicant": "flagged toxicant-related theme",
        "broad_cytotoxic/topoisomerase": "flagged cytotoxic/topoisomerase theme",
        "phosphodiesterase/cAMP": "phosphodiesterase / cAMP",
        "calcium/synaptic": "calcium / synaptic",
        "redox/stress-adaptation": "redox / stress adaptation",
        "epigenetic/hdac": "epigenetic / HDAC",
        "metabolic/homeostatic; epigenetic/hdac": "metabolic/homeostatic; epigenetic / HDAC",
        "metabolic/homeostatic": "metabolic / homeostatic",
        "thyroid/hormone": "thyroid / hormone",
        "mapk/kinase": "MAPK / kinase",
    }
    mapped = mapping.get(text, text or "unassigned")
    if any(ord(ch) > 127 for ch in mapped):
        return "manual review needed"
    return mapped


def is_flagged_theme(theme: str) -> bool:
    return theme.startswith("flagged ")


def selected_clue_rows(full: pd.DataFrame) -> pd.DataFrame:
    return full[full["after_bbb_final_layer"].isin(LAYER_ORDER)].copy()


def build_disclaimer(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.add_patch(plt.Rectangle((0.02, 0.22), 0.96, 0.56, facecolor="#F4EFEF", edgecolor="#C9A7AA", linewidth=0.7))
    ax.text(0.5, 0.50, DISCLAIMER, ha="center", va="center", fontsize=8.5, fontweight="bold", color="#6D3038")


def build_panel_a(ax: plt.Axes, full: pd.DataFrame, fig: plt.Figure) -> None:
    selected = selected_clue_rows(full)
    count_mat = (
        selected.groupby(["after_bbb_final_layer", "associated_axis"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=LAYER_ORDER, columns=AXIS_ORDER)
        .fillna(0)
        .astype(int)
    )
    vmax = max(1.0, float(count_mat.values.max()))
    im = ax.imshow(count_mat.values, aspect="auto", cmap="OrRd", vmin=0, vmax=vmax)
    ax.set_xticks(np.arange(count_mat.shape[1]))
    ax.set_xticklabels(count_mat.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(count_mat.shape[0]))
    ax.set_yticklabels([LAYER_LABELS[idx] for idx in count_mat.index])
    for i in range(count_mat.shape[0]):
        for j in range(count_mat.shape[1]):
            val = int(count_mat.iloc[i, j])
            ax.text(j, i, str(val), ha="center", va="center", fontsize=6.2, color="white" if val >= 2 else "#222222")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("BBB-layered exploratory clue matrix", pad=7)
    cbar = fig.colorbar(im, ax=ax, fraction=0.050, pad=0.035)
    cbar.set_label("Count")
    cbar.ax.tick_params(labelsize=6.0, width=0.5, length=2.4)
    panel_label(ax, "A")


def build_panel_b(ax: plt.Axes, paper: pd.DataFrame) -> None:
    plot_df = paper.sort_values("final_score_after_bbb", ascending=False).copy()
    y = np.arange(len(plot_df))
    scores = as_num(plot_df["final_score_after_bbb"]).to_numpy(dtype=float)
    colors = [LAYER_COLORS.get(layer, "#A8ADB2") for layer in plot_df["after_bbb_final_layer"]]
    ax.hlines(y, 0, scores, color="#D4D4D4", linewidth=1.0, zorder=1)
    ax.scatter(scores, y, s=42, c=colors, edgecolor="white", linewidth=0.7, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([short_compound_name(x) for x in plot_df["compound"]])
    ax.invert_yaxis()
    ax.set_xlim(0, max(1.02, float(np.nanmax(scores)) + 0.05))
    ax.set_xlabel("After-BBB exploratory score")
    ax.set_title("Exploratory after-BBB ranked drug-signature clues", pad=7)
    style_axes(ax, grid_axis="x")
    panel_label(ax, "B")


def build_panel_c(ax: plt.Axes, full: pd.DataFrame, paper: pd.DataFrame) -> None:
    reviewed = full[full["after_bbb_final_layer"].notna()].copy()
    reviewed["final_score_prebbb_v21"] = as_num(reviewed["final_score_prebbb_v21"])
    reviewed["final_score_after_bbb"] = as_num(reviewed["final_score_after_bbb"])
    colors = [LAYER_COLORS.get(layer, "#D7DADB") for layer in reviewed["after_bbb_final_layer"]]
    sizes = np.where(reviewed["after_bbb_final_layer"].isin(LAYER_ORDER), 43, 20)
    alphas = np.where(reviewed["after_bbb_final_layer"].isin(LAYER_ORDER), 0.82, 0.40)
    ax.scatter(
        reviewed["final_score_prebbb_v21"],
        reviewed["final_score_after_bbb"],
        s=sizes,
        c=colors,
        alpha=alphas,
        edgecolor="white",
        linewidth=0.55,
        zorder=3,
    )
    ax.plot([0, 1], [0, 1], color="#BEBEBE", linewidth=0.8, linestyle="--", zorder=2)
    paper_compounds = set(paper["compound"].astype(str))
    labels = reviewed[reviewed["compound"].astype(str).isin(paper_compounds)].sort_values("final_score_after_bbb", ascending=False).head(4)
    offsets = {
        "VALPROIC ACID": (-92, -21),
        "lobeline": (9, -18),
        "(-)-Epigallocatechin gallate": (-126, 7),
        "S-(+)-Rolipram": (13, 14),
        "papaverine": (13, -13),
        "quercetin": (-84, -22),
    }
    for _, row in labels.iterrows():
        name = str(row["compound"])
        ax.annotate(
            short_compound_name(name, 25),
            (float(row["final_score_prebbb_v21"]), float(row["final_score_after_bbb"])),
            xytext=offsets.get(name, (5, 5)),
            textcoords="offset points",
            fontsize=5.8,
            color="#444444",
            arrowprops={"arrowstyle": "-", "color": "#8A8A8A", "lw": 0.42, "shrinkA": 0, "shrinkB": 5},
        )
    ax.set_xlabel("Pre-BBB exploratory score")
    ax.set_ylabel("After-BBB exploratory score")
    ax.set_title("Pre-BBB versus after-BBB exploratory score shift", pad=7)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    style_axes(ax, grid_axis="both")
    panel_label(ax, "C")


def build_panel_d(ax: plt.Axes, full: pd.DataFrame, fig: plt.Figure) -> None:
    reviewed = full[full["after_bbb_final_layer"].notna()].copy()
    reviewed["display_theme"] = reviewed["indicative_mechanism"].map(display_theme)
    selected = reviewed[reviewed["after_bbb_final_layer"].isin(LAYER_ORDER)].copy()
    selected_counts = selected["display_theme"].value_counts()
    theme_order = selected_counts.index.tolist()
    for theme in ["flagged toxicant-related theme", "flagged cytotoxic/topoisomerase theme"]:
        if (reviewed["display_theme"] == theme).any() and theme not in theme_order:
            theme_order.append(theme)
    theme_order = theme_order[:8]
    theme_mat = (
        selected.groupby(["display_theme", "after_bbb_final_layer"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=theme_order, columns=LAYER_ORDER)
        .fillna(0)
        .astype(int)
    )
    vmax = max(1.0, float(theme_mat.values.max()))
    im = ax.imshow(theme_mat.values, aspect="auto", cmap="Greens", vmin=0, vmax=vmax)
    ax.set_xticks(np.arange(theme_mat.shape[1]))
    ax.set_xticklabels([LAYER_LABELS[c] for c in theme_mat.columns], rotation=25, ha="right")
    ax.set_yticks(np.arange(theme_mat.shape[0]))
    ax.set_yticklabels(theme_mat.index)
    for label in ax.get_yticklabels():
        if is_flagged_theme(label.get_text()):
            label.set_color("#8A8F93")
            label.set_fontstyle("italic")
    for i in range(theme_mat.shape[0]):
        flagged = is_flagged_theme(str(theme_mat.index[i]))
        for j in range(theme_mat.shape[1]):
            val = int(theme_mat.iloc[i, j])
            color = "#8A8F93" if flagged else ("white" if val >= 2 else "#222222")
            ax.text(j, i, str(val), ha="center", va="center", fontsize=6.0, color=color)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Mechanism-theme distribution across exploratory layers", pad=7)
    cbar = fig.colorbar(im, ax=ax, fraction=0.050, pad=0.035)
    cbar.set_label("Count")
    cbar.ax.tick_params(labelsize=6.0, width=0.5, length=2.4)
    panel_label(ax, "D")


def build_panel_e(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.90, "Conservative interpretation of exploratory drug-signature clues", fontsize=8.4, fontweight="bold", va="top", ha="left")
    entries = [
        ("mechanism-direction clues", "enrichment-derived, BBB-reviewed, hypothesis-generating only", LAYER_COLORS[LAYER_ORDER[0]]),
        ("supportive clues", "secondary/supportive signal; requires manual review", LAYER_COLORS[LAYER_ORDER[1]]),
        ("peripheral/program-modulating clues", "program-linked but not prioritized for CNS-directed interpretation", LAYER_COLORS[LAYER_ORDER[2]]),
    ]
    for i, (label, desc, color) in enumerate(entries):
        x0 = 0.02 + i * 0.32
        ax.plot([x0, x0 + 0.040], [0.55, 0.55], color=color, linewidth=3.0, solid_capstyle="round")
        ax.text(x0 + 0.050, 0.62, label, fontsize=6.9, fontweight="bold", ha="left", va="center", color="#222222")
        ax.text(x0 + 0.050, 0.40, desc, fontsize=6.5, ha="left", va="center", color="#333333")
    ax.text(0.02, 0.13, PANEL_E_FINAL, fontsize=6.8, ha="left", va="center", color="#5E3338", fontweight="bold")
    panel_label(ax, "E", x=-0.035, y=1.00)


def save_figure(fig: plt.Figure) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = [
        OUT_DIR / f"{STEM}.png",
        OUT_DIR / f"{STEM}.pdf",
        OUT_DIR / f"{STEM}.tiff",
    ]
    fig.savefig(outputs[0], dpi=450, bbox_inches="tight")
    fig.savefig(outputs[1], bbox_inches="tight")
    try:
        fig.savefig(outputs[2], dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    except TypeError:
        fig.savefig(outputs[2], dpi=600, bbox_inches="tight")
    plt.close(fig)
    return outputs


def build_figure() -> list[Path]:
    setup_style()
    paper = read_csv(INPUT_FILES[0])
    full = read_csv(INPUT_FILES[1])

    fig = plt.figure(figsize=(13.4, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(
        4,
        2,
        width_ratios=[0.96, 1.12],
        height_ratios=[0.15, 0.95, 1.08, 0.36],
        wspace=0.12,
        hspace=0.12,
    )
    ax_disc = fig.add_subplot(gs[0, :])
    ax_a = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_b = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])
    ax_e = fig.add_subplot(gs[3, :])

    build_disclaimer(ax_disc)
    build_panel_a(ax_a, full, fig)
    build_panel_c(ax_c, full, paper)
    build_panel_b(ax_b, paper)
    build_panel_d(ax_d, full, fig)
    build_panel_e(ax_e)

    fig.suptitle("Exploratory drug-signature enrichment and BBB-layered prioritization", y=1.006, fontweight="bold")
    return save_figure(fig)


def write_text(path: Path, lines: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).strip() + "\n", encoding="utf-8")
    return path


def write_notes() -> Path:
    lines = [
        "FigS7 drug repositioning exploratory revised_v2 notes",
        "",
        "Input files used:",
        *[f"- {rel(p)}" for p in INPUT_FILES],
        "",
        "Score sources:",
        "- Pre-BBB score is read from final_score_prebbb_v21 in drug_repositioning/13_integrated_after_bbb_v21.csv for Panel C, and from prebbb_score_v21 in drug_repositioning/16_paper_ready_table_after_bbb_v21.csv for the paper-table subset.",
        "- After-BBB score is read from final_score_after_bbb in drug_repositioning/13_integrated_after_bbb_v21.csv and drug_repositioning/16_paper_ready_table_after_bbb_v21.csv.",
        "- No scores were recomputed for this revised figure.",
        "",
        "BBB review / SwissADME integration:",
        "- SwissADME parsed and matched properties are read from drug_repositioning/11_swissadme_parsed_and_matched_v21.csv.",
        "- Manual BBB review annotations are read from drug_repositioning/12_bbb_annotation_manual_review_v21.csv.",
        "- The integrated after-BBB table combines the v2.1 pre-BBB score with the BBB_score and manual_final_bbb_layer to produce final_score_after_bbb and after_bbb_final_layer.",
        "",
        "Compound layer definitions used for figure display:",
        "- headline_cns_mechanism_direction_leads was relabeled as mechanism-direction clues.",
        "- supportive_after_bbb_leads was relabeled as supportive clues.",
        "- peripheral_program_modulating_clues was relabeled as peripheral/program-modulating clues.",
        "- These layer names are exploratory evidence categories and are not recommendation tiers.",
        "",
        "Mechanism theme definitions:",
        "- mechanism themes are read from indicative_mechanism in the integrated after-BBB table and from mechanism_theme in the paper-ready after-BBB subset.",
        "- manual_review_needed is displayed as manual review needed.",
        "- secondary_axis_clue is displayed as secondary-axis clue.",
        "- environmental_toxicant is displayed as flagged toxicant-related theme.",
        "- broad_cytotoxic/topoisomerase is displayed as flagged cytotoxic/topoisomerase theme.",
        "- phosphodiesterase/cAMP is displayed as phosphodiesterase / cAMP.",
        "- calcium/synaptic is displayed as calcium / synaptic.",
        "- flagged toxicant-related and cytotoxic/topoisomerase themes are cautionary themes and are not prioritized mechanism-direction interpretations.",
        "",
        "Interpretation boundary:",
        "- Compounds shown in this figure are enrichment-derived exploratory clues from drug-signature enrichment.",
        "- This figure does not report wet-lab validation.",
        "- This figure does not report docking.",
        "- This figure has no treatment recommendation meaning.",
        "- The module is not FCD II-specific drug discovery and is not a clinical prioritization result.",
        "- VALPROIC ACID, lobeline, EGCG, quercetin, HDAC-related compounds, and PDE-related compounds are shown as exploratory mechanism-direction or program-modulating clues only.",
    ]
    return write_text(OUT_DIR / f"{STEM}_notes.txt", lines)


def write_caption() -> Path:
    lines = [
        "# Supplementary Figure S7. Exploratory drug-signature enrichment and BBB-layered prioritization",
        "",
        "Drug-signature enrichment was analyzed as an exploratory supplementary module. The compounds shown here represent program-guided mechanism-direction clues and should not be interpreted as therapeutic recommendations, validated drugs, or clinical candidates for FCD II.",
        "",
        "A. After-BBB exploratory clue layer matrix showing counts of exploratory clue categories across TF-associated programs.",
        "",
        "B. Ranked exploratory after-BBB drug-signature clues based on the final exploratory score after BBB review.",
        "",
        "C. Pre-BBB versus after-BBB exploratory score shift, with the dashed diagonal indicating score preservation or shift reference.",
        "",
        "D. Mechanism-theme distribution across exploratory layers. Toxicant-related and broad cytotoxic/topoisomerase themes are flagged as cautionary themes rather than positive mechanism-direction interpretations.",
        "",
        "E. Conservative interpretation strip summarizing the three exploratory clue categories and the boundary that no compound is interpreted as a therapy or recommendation.",
    ]
    return write_text(OUT_DIR / f"{STEM}_caption.md", lines)


def write_revision_log() -> Path:
    lines = [
        "# FigS7 drug repositioning exploratory revised_v2 revision log",
        "",
        "- Rebuilt only Supplementary Figure S7 and wrote all outputs to manuscript_output/figures_supplementary_revised_v2.",
        "- Changed Paper-ready after-BBB ranked compounds to Exploratory after-BBB ranked drug-signature clues.",
        "- Changed headline/supportive/peripheral labels to mechanism-direction clues/supportive clues/peripheral-program-modulating clues.",
        "- Added a visible disclaimer strip: Exploratory drug-signature enrichment only; not therapeutic recommendation.",
        "- Removed or replaced strong therapeutic, treatment, validated, and clinical-candidate style wording from positive figure claims.",
        "- Marked toxicant-related and cytotoxic/topoisomerase mechanism themes as flagged/cautionary themes.",
        "- Preserved the original input data, compound list, after-BBB scores, pre-BBB scores, and ranked order by after-BBB exploratory score.",
        "- Did not modify main-text Figure 1-5 or Supplementary FigS1-FigS6.",
    ]
    return write_text(OUT_DIR / f"{STEM}_revision_log.md", lines)


def figure_text() -> str:
    text = "\n".join(FIGURE_TEXT_ITEMS)
    for theme in [
        "manual review needed",
        "secondary-axis clue",
        "flagged toxicant-related theme",
        "flagged cytotoxic/topoisomerase theme",
        "phosphodiesterase / cAMP",
        "calcium / synaptic",
        "redox / stress adaptation",
        "epigenetic / HDAC",
    ]:
        text += "\n" + theme
    return text


def audit_status() -> dict[str, str]:
    text = figure_text().lower()
    return {
        "therapeutic recommendation": "present only in required negated disclaimer; absent as positive claim",
        "treatment recommendation": "present only in required negated interpretation strip; absent as positive claim",
        "validated drug": "absent from figure text; negated boundary-only in caption",
        "clinical candidate": "absent from figure text; negated boundary-only in caption",
        "therapeutic candidate": "absent",
        "treatment lead": "absent",
        "drug validation": "absent",
        "FCD II treatment": "absent",
        "docking": "absent from figure and caption; negated boundary-only in notes",
        "target-based discovery": "absent",
        "Paper-ready": "replaced",
        "headline": "replaced",
        "supportive clues": "present" if "supportive clues" in text else "absent",
        "peripheral/program-modulating clues": "present" if "peripheral/program-modulating clues" in text else "absent",
        "exploratory": "present" if "exploratory" in text else "absent",
        "drug-signature enrichment": "present" if "drug-signature enrichment" in text else "absent",
        "mechanism-direction clues": "present" if "mechanism-direction clues" in text else "absent",
        "not therapeutic recommendation": "present" if "not therapeutic recommendation" in text else "absent",
        "not validated therapy": "present" if "validated therapy" in text else "absent",
    }


def write_terminology_audit() -> Path:
    rows = [
        ("therapeutic recommendation", audit_status()["therapeutic recommendation"], "required absent as positive claim; required as negated disclaimer"),
        ("treatment recommendation", audit_status()["treatment recommendation"], "required absent as positive claim; required as negated interpretation"),
        ("validated drug", audit_status()["validated drug"], "required absent as positive claim"),
        ("clinical candidate", audit_status()["clinical candidate"], "required absent as positive claim"),
        ("therapeutic candidate", audit_status()["therapeutic candidate"], "required absent"),
        ("treatment lead", audit_status()["treatment lead"], "required absent"),
        ("drug validation", audit_status()["drug validation"], "required absent"),
        ("FCD II treatment", audit_status()["FCD II treatment"], "required absent"),
        ("docking", audit_status()["docking"], "required absent as analysis claim"),
        ("target-based discovery", audit_status()["target-based discovery"], "required absent"),
        ("Paper-ready", audit_status()["Paper-ready"], "replaced by exploratory"),
        ("headline", audit_status()["headline"], "replaced by mechanism-direction clues"),
        ("supportive clues", audit_status()["supportive clues"], "required"),
        ("peripheral/program-modulating clues", audit_status()["peripheral/program-modulating clues"], "required"),
        ("exploratory", audit_status()["exploratory"], "required"),
        ("drug-signature enrichment", audit_status()["drug-signature enrichment"], "required"),
        ("mechanism-direction clues", audit_status()["mechanism-direction clues"], "required"),
        ("not therapeutic recommendation", audit_status()["not therapeutic recommendation"], "required"),
        ("not validated therapy", audit_status()["not validated therapy"], "required"),
    ]
    path = OUT_DIR / f"{STEM}_terminology_audit.csv"
    pd.DataFrame(rows, columns=["term", "status", "action"]).to_csv(path, index=False)
    return path


def terminal_residuals() -> dict[str, str]:
    text = figure_text().lower()
    checks = {
        "therapeutic": "boundary-only" if "therapeutic" in text else "NO",
        "treatment": "boundary-only" if "treatment" in text else "NO",
        "validated": "boundary-only" if "validated" in text else "NO",
        "clinical candidate": "YES" if "clinical candidate" in text else "NO",
        "docking": "YES" if "docking" in text else "NO",
    }
    return checks


def main() -> None:
    outputs = build_figure()
    notes = write_notes()
    caption = write_caption()
    revision_log = write_revision_log()
    audit = write_terminology_audit()

    print("FigS7 revised_v2 output path:")
    print(f"- {rel(OUT_DIR)}")
    for path in outputs:
        print(f"- {rel(path)}")

    print("\nInput files used:")
    for path in INPUT_FILES:
        print(f"- {rel(path)}")

    print("\nDisclaimer strip added:")
    print("- YES")

    print("\nPaper-ready replaced by Exploratory:")
    print("- YES")

    print("\nheadline/supportive/peripheral replaced by downgraded terms:")
    print("- YES; mechanism-direction clues / supportive clues / peripheral-program-modulating clues")

    print("\nResidual restricted wording in figure text:")
    for term, status in terminal_residuals().items():
        print(f"- {term}: {status}")
    print("- No positive therapeutic, treatment, validated, clinical-candidate, or docking claim is present.")

    print("\nCompanion files generated:")
    print(f"- caption: {rel(caption)}")
    print(f"- notes: {rel(notes)}")
    print(f"- revision log: {rel(revision_log)}")
    print(f"- terminology audit: {rel(audit)}")

    print("\nManual review flag:")
    print("- No panel requires layout repair based on scripted checks; Panel D flagged cautionary rows should be reviewed for final wording preference.")


if __name__ == "__main__":
    main()
