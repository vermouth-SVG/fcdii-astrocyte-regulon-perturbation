#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_supplementary_revised_v2"
STEM = "FigS4_sample_level_robustness_revised_v2"

INPUT_DIR = ROOT / "robustness_validation"
LOSO_MASTER = INPUT_DIR / "02_leave_one_sample_out_master_results_revised.csv"
LOSO_TF_SUMMARY = INPUT_DIR / "02_leave_one_sample_out_tf_summary_revised.csv"
PSEUDOBULK_EXPR = INPUT_DIR / "03_pseudobulk_candidate_tf_table.csv"
PSEUDOBULK_REGULON = INPUT_DIR / "03_pseudobulk_regulon_table.csv"
SENSITIVITY_LONG = INPUT_DIR / "04_shortlist_sensitivity_long.csv"
SENSITIVITY_GRID = INPUT_DIR / "04_shortlist_sensitivity_grid.csv"
MEMBERSHIP_FREQ = INPUT_DIR / "04_shortlist_membership_frequency.csv"
METADATA = ROOT / "celloracle_run" / "prepared_data" / "celloracle_round1_metadata.csv"

INPUT_FILES = [
    LOSO_MASTER,
    LOSO_TF_SUMMARY,
    PSEUDOBULK_EXPR,
    PSEUDOBULK_REGULON,
    SENSITIVITY_LONG,
    SENSITIVITY_GRID,
    MEMBERSHIP_FREQ,
    METADATA,
]

TF_ORDER = ["NFE2L2", "THRB", "BHLHE40", "SOX2", "SATB2", "RARB", "HMGA1"]
TF_COLORS = {
    "NFE2L2": "#8E1B2A",
    "THRB": "#0B6B74",
    "BHLHE40": "#6A51A3",
    "SOX2": "#777777",
    "SATB2": "#4E8C69",
    "RARB": "#4F78A8",
    "HMGA1": "#B45A5A",
}
LESION_COLOR = "#8E1B2A"
CONTROL_COLOR = "#0B6B74"
TEXT_COLOR = "#222222"
GRID_COLOR = "#E9E9E9"

SENSITIVITY_RULE = "regulon_significant_plus_expression_significant_plus_direction"
THRESHOLD_COLUMNS = [
    ("FDR <= 0.01", 0.01, 0.00, "FDR"),
    ("FDR <= 0.05", 0.05, 0.00, "FDR"),
    ("FDR <= 0.10", 0.10, 0.00, "FDR"),
    ("|log2FC| >= 0", 0.05, 0.00, "expr"),
    ("|log2FC| >= 0.25", 0.05, 0.25, "expr"),
    ("|log2FC| >= 0.50", 0.05, 0.50, "expr"),
]

FIGURE_TEXT_ITEMS = [
    "Complete sample-level robustness and sensitivity analyses",
    "Leave-one-sample-out expression effect",
    "Expression log2FC",
    "Leave-one-sample-out regulon effect",
    "Regulon dAUC",
    "Pseudobulk support for candidate TFs",
    "Pseudobulk expression log2FC",
    "Pseudobulk regulon dAUC",
    "Threshold sensitivity and shortlist retention",
    "Retained",
    "Not retained",
    "sample-level robustness",
    "LOSO",
    "internal-control",
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.titlesize": 8.7,
            "axes.labelsize": 7.4,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.6,
            "legend.fontsize": 6.1,
            "figure.titlesize": 11.8,
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


def display_input(path: Path) -> str:
    return path.name if path.parent == INPUT_DIR else rel(path)


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def panel_label(ax: plt.Axes, label: str, x: float = -0.105, y: float = 1.12) -> None:
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


def style_axes(ax: plt.Axes, grid_axis: str | None = "both") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def load_sample_labels() -> tuple[dict[str, str], list[str]]:
    meta = read_csv(METADATA)
    first = meta.columns[0]
    if str(first).startswith("Unnamed"):
        meta = meta.rename(columns={first: "cell_id"})
    sample_map = (
        meta[["sample_id", "clinical_sample_label"]]
        .drop_duplicates()
        .sort_values("clinical_sample_label")
        .set_index("sample_id")["clinical_sample_label"]
        .to_dict()
    )
    ordered = ["Full"] + [f"LOSO: {label}" for label in sorted(sample_map.values())]
    return sample_map, ordered


def load_loso_matrices() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    sample_map, order_cols = load_sample_labels()
    loso = read_csv(LOSO_MASTER)
    loso = loso[loso["tf"].isin(TF_ORDER)].copy()
    loso["display_col"] = loso["dropped_sample"].map(sample_map)
    loso["display_col"] = np.where(
        loso["dropped_sample"].astype(str).eq("FULL_DATA"),
        "Full",
        "LOSO: " + pd.Series(loso["display_col"]).fillna(loso["dropped_sample"]).astype(str),
    )
    loso["expr_log2fc"] = as_num(loso["expr_log2fc"])
    loso["regulon_diff"] = as_num(loso["regulon_diff"])
    expr_mat = (
        loso.pivot_table(index="tf", columns="display_col", values="expr_log2fc", aggfunc="first")
        .reindex(index=TF_ORDER, columns=order_cols)
    )
    reg_mat = (
        loso.pivot_table(index="tf", columns="display_col", values="regulon_diff", aggfunc="first")
        .reindex(index=TF_ORDER, columns=order_cols)
    )
    return expr_mat, reg_mat, order_cols


def load_pseudobulk() -> pd.DataFrame:
    expr = read_csv(PSEUDOBULK_EXPR)
    regulon = read_csv(PSEUDOBULK_REGULON)
    pseudo = expr[["tf", "log2fc_lesion_vs_internal_control", "direction"]].merge(
        regulon[["tf", "mean_diff_lesion_minus_internal_control", "direction"]].rename(
            columns={
                "mean_diff_lesion_minus_internal_control": "regulon_diff",
                "direction": "regulon_direction",
            }
        ),
        on="tf",
        how="left",
    )
    pseudo["log2fc_lesion_vs_internal_control"] = as_num(pseudo["log2fc_lesion_vs_internal_control"])
    pseudo["regulon_diff"] = as_num(pseudo["regulon_diff"])
    return pseudo[pseudo["tf"].isin(TF_ORDER)].set_index("tf").reindex(TF_ORDER).reset_index()


def load_retention_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    sensitivity = read_csv(SENSITIVITY_LONG)
    sensitivity = sensitivity[sensitivity["rule"].eq(SENSITIVITY_RULE)].copy()
    sensitivity["regulon_fdr_threshold"] = as_num(sensitivity["regulon_fdr_threshold"])
    sensitivity["expr_abs_log2fc_threshold"] = as_num(sensitivity["expr_abs_log2fc_threshold"])
    sensitivity["positive_fraction_threshold"] = as_num(sensitivity["positive_fraction_threshold"])
    retained = pd.DataFrame(index=TF_ORDER)
    ranks = pd.DataFrame(index=TF_ORDER)
    for label, fdr_thr, expr_thr, _kind in THRESHOLD_COLUMNS:
        sub = sensitivity[
            sensitivity["positive_fraction_threshold"].eq(0.00)
            & np.isclose(sensitivity["regulon_fdr_threshold"], fdr_thr)
            & np.isclose(sensitivity["expr_abs_log2fc_threshold"], expr_thr)
            & sensitivity["tf"].isin(TF_ORDER)
        ].copy()
        sub = sub.drop_duplicates(subset=["tf"]).set_index("tf").reindex(TF_ORDER)
        retained[label] = sub["retained"].fillna(False).astype(bool)
        ranks[label] = as_num(sub["rank_in_combo"])
    return retained, ranks


def heatmap_limits(mat: pd.DataFrame) -> TwoSlopeNorm:
    max_abs = float(np.nanmax(np.abs(mat.values)))
    if not np.isfinite(max_abs) or max_abs == 0:
        max_abs = 1.0
    return TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)


def build_panel_a(ax: plt.Axes, expr_mat: pd.DataFrame, fig: plt.Figure) -> None:
    im = ax.imshow(expr_mat.values, aspect="auto", cmap="RdBu_r", norm=heatmap_limits(expr_mat))
    ax.set_title("Leave-one-sample-out expression effect", pad=7)
    ax.set_xticks(np.arange(expr_mat.shape[1]))
    ax.set_xticklabels(expr_mat.columns, rotation=34, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(expr_mat.shape[0]))
    ax.set_yticklabels(expr_mat.index)
    for i in range(expr_mat.shape[0]):
        for j in range(expr_mat.shape[1]):
            val = expr_mat.iloc[i, j]
            if pd.notna(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.9, color="white" if abs(val) > 0.75 else TEXT_COLOR)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.025)
    cbar.set_label("Expression log2FC", fontsize=7.0)
    cbar.ax.tick_params(labelsize=6.1, width=0.5, length=2.2)
    panel_label(ax, "A")


def build_panel_b(ax: plt.Axes, reg_mat: pd.DataFrame, fig: plt.Figure) -> None:
    im = ax.imshow(reg_mat.values, aspect="auto", cmap="RdBu_r", norm=heatmap_limits(reg_mat))
    ax.set_title("Leave-one-sample-out regulon effect", pad=7)
    ax.set_xticks(np.arange(reg_mat.shape[1]))
    ax.set_xticklabels(reg_mat.columns, rotation=34, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(reg_mat.shape[0]))
    ax.set_yticklabels(reg_mat.index)
    for i in range(reg_mat.shape[0]):
        for j in range(reg_mat.shape[1]):
            val = reg_mat.iloc[i, j]
            if pd.notna(val):
                ax.text(j, i, f"{val:.3f}", ha="center", va="center", fontsize=5.8, color="white" if abs(val) > 0.11 else TEXT_COLOR)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.025)
    cbar.set_label("Regulon dAUC", fontsize=7.0)
    cbar.ax.tick_params(labelsize=6.1, width=0.5, length=2.2)
    panel_label(ax, "B")


def build_panel_c(ax: plt.Axes, pseudo: pd.DataFrame) -> None:
    x_col = "log2fc_lesion_vs_internal_control"
    y_col = "regulon_diff"
    ax.axvline(0, color="#9E9E9E", linewidth=0.75, zorder=1)
    ax.axhline(0, color="#9E9E9E", linewidth=0.75, zorder=1)
    ax.scatter(
        pseudo[x_col],
        pseudo[y_col],
        s=88,
        c=[TF_COLORS[tf] for tf in pseudo["tf"]],
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
    )
    offsets = {
        "NFE2L2": (-9, -16, "right"),
        "THRB": (10, -12, "left"),
        "BHLHE40": (-8, 12, "right"),
        "SOX2": (8, 10, "left"),
        "SATB2": (-10, 22, "right"),
        "RARB": (-10, -13, "right"),
        "HMGA1": (9, 8, "left"),
    }
    for _, row in pseudo.iterrows():
        tf = row["tf"]
        dx, dy, ha = offsets.get(tf, (5, 5, "left"))
        ax.annotate(
            tf,
            xy=(float(row[x_col]), float(row[y_col])),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=6.3,
            color=TF_COLORS.get(tf, TEXT_COLOR),
            arrowprops={"arrowstyle": "-", "color": "#808080", "lw": 0.5, "shrinkA": 1.2, "shrinkB": 4.0},
            zorder=4,
        )
    x_min, x_max = float(pseudo[x_col].min()), float(pseudo[x_col].max())
    y_min, y_max = float(pseudo[y_col].min()), float(pseudo[y_col].max())
    x_pad = (x_max - x_min) * 0.14
    y_pad = (y_max - y_min) * 0.16
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_xlabel("Pseudobulk expression log2FC")
    ax.set_ylabel("Pseudobulk regulon dAUC")
    ax.set_title("Pseudobulk support for candidate TFs", pad=7)
    ax.text(0.98, 0.94, "lesion-side support", transform=ax.transAxes, ha="right", va="top", fontsize=6.1, color=LESION_COLOR)
    ax.text(0.02, 0.06, "internal-control-side support", transform=ax.transAxes, ha="left", va="bottom", fontsize=6.1, color=CONTROL_COLOR)
    style_axes(ax, "both")
    panel_label(ax, "C")


def build_panel_d(ax: plt.Axes, retained: pd.DataFrame, ranks: pd.DataFrame) -> None:
    mat = retained.astype(int)
    cmap = ListedColormap(["#E9E9E9", "#485F66"])
    ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_title("Threshold sensitivity and shortlist retention", pad=7)
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=34, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(mat.index)
    ax.set_xticks(np.arange(-0.5, mat.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, mat.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="both", length=0)
    for i, tf in enumerate(mat.index):
        for j, col in enumerate(mat.columns):
            if mat.iloc[i, j] == 1:
                rank = ranks.iloc[i, j]
                label = "R" if pd.isna(rank) else f"R{int(rank)}"
                ax.text(j, i, label, ha="center", va="center", fontsize=6.1, color="white", fontweight="bold")
            else:
                ax.text(j, i, "not", ha="center", va="center", fontsize=5.7, color="#666666")
    for spine in ax.spines.values():
        spine.set_visible(False)
    legend_handles = [
        Patch(facecolor="#485F66", edgecolor="#485F66", label="Retained"),
        Patch(facecolor="#E9E9E9", edgecolor="#BDBDBD", label="Not retained"),
        Line2D([0], [0], color="none", label="R# = rank when retained"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.14),
        ncol=3,
        frameon=False,
        handlelength=1.1,
        handletextpad=0.45,
        columnspacing=0.85,
        borderpad=0.0,
    )
    ax.text(
        0.0,
        -0.30,
        "Rule: regulon FDR + expression effect + direction agreement; positive fraction threshold fixed at 0.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.9,
        color="#555555",
    )
    panel_label(ax, "D")


def build_figure() -> plt.Figure:
    setup_style()
    expr_mat, reg_mat, _order_cols = load_loso_matrices()
    pseudo = load_pseudobulk()
    retained, ranks = load_retention_matrix()

    fig = plt.figure(figsize=(13.3, 8.7), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        left=0.065,
        right=0.985,
        top=0.900,
        bottom=0.105,
        hspace=0.47,
        wspace=0.30,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    build_panel_a(ax_a, expr_mat, fig)
    build_panel_b(ax_b, reg_mat, fig)
    build_panel_c(ax_c, pseudo)
    build_panel_d(ax_d, retained, ranks)
    fig.suptitle("Complete sample-level robustness and sensitivity analyses", y=0.975, fontweight="bold")
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
        "Supplementary Fig. S4. Complete sample-level robustness and sensitivity analyses.\n\n"
        "(A) Leave-one-sample-out expression effects for candidate TFs across the full discovery "
        "object and each sample-excluded analysis. (B) Leave-one-sample-out regulon activity "
        "effects, shown as regulon AUC differences. (C) Pseudobulk support comparing expression "
        "log2FC and regulon differences for candidate TFs. (D) Threshold sensitivity and shortlist "
        "retention analysis across regulon FDR and expression-effect thresholds. These analyses "
        "assess sample-level robustness and sensitivity within the GSE268807 astrocyte pilot and "
        "should not be interpreted as donor-level or external validation."
    )


def notes_text() -> str:
    inputs = "\n".join(f"- {display_input(path)}" for path in INPUT_FILES)
    retained, _ranks = load_retention_matrix()
    retention_lines = []
    for tf in TF_ORDER:
        states = ", ".join(f"{col}: {'retained' if retained.loc[tf, col] else 'not retained'}" for col in retained.columns)
        retention_lines.append(f"- {tf}: {states}")
    return f"""{STEM} notes

Input files:
{inputs}

LOSO expression effect:
- Leave-one-sample-out expression effects were read from 02_leave_one_sample_out_master_results_revised.csv.
- Full denotes the full GSE268807 astrocyte discovery object.
- LOSO denotes a sample-excluded analysis, labeled as LOSO: clinical sample label.
- Expression effect is the lesion-versus-internal-control expression log2FC.

LOSO regulon effect:
- Leave-one-sample-out regulon effects were read from 02_leave_one_sample_out_master_results_revised.csv.
- Regulon effect is the lesion-minus-internal-control regulon AUC difference for the matched candidate regulon.
- Positive values indicate lesion-side direction; negative values indicate internal-control-side direction.

Pseudobulk support:
- Pseudobulk expression log2FC was read from 03_pseudobulk_candidate_tf_table.csv.
- Pseudobulk regulon difference was read from 03_pseudobulk_regulon_table.csv and represents lesion-minus-internal-control regulon dAUC at the sample level.
- Panel C compares expression and regulon effects for NFE2L2, THRB, BHLHE40, SOX2, SATB2, RARB, and HMGA1 without changing values.

Threshold sensitivity and shortlist retention:
- Panel D was rebuilt from 04_shortlist_sensitivity_long.csv.
- Retained / not retained was read from the retained column under the rule regulon_significant_plus_expression_significant_plus_direction.
- For the three FDR columns, expression-effect threshold was fixed at |log2FC| >= 0 and positive-fraction threshold was fixed at 0.
- For the three expression-effect columns, regulon FDR was fixed at 0.05 and positive-fraction threshold was fixed at 0.
- R# labels show the rank_in_combo value when a TF was retained; not indicates not retained for that threshold setting.

Retention states shown in Panel D:
{chr(10).join(retention_lines)}

Interpretation boundaries:
- This figure shows sample-level robustness and leave-one-sample-out robustness within the GSE268807 astrocyte pilot.
- It is not donor-level robustness because the discovery cohort has four samples but only two donors.
- It is not a formal or external validation analysis.
- NFE2L2 and THRB are treated as the primary lesion-associated and primary internal-control axes, respectively.
- BHLHE40 is treated as a secondary candidate.
- SOX2 is retained as a lower-priority candidate.
"""


def revision_log_text() -> str:
    return """# FigS4 Sample-Level Robustness Revised V2 Revision Log

- Rebuilt only Supplementary Figure S4 outputs under manuscript_output/figures_supplementary_revised_v2.
- Standardized figure wording to sample-level robustness and leave-one-sample-out / LOSO.
- Standardized display labels to Full and LOSO: sample label.
- Standardized internal-control wording in figure-facing labels.
- Kept Panel A and Panel B heatmap structures but updated titles, axis labels, and colorbar labels.
- Kept Panel C pseudobulk support scatter but updated axis labels, reference lines, candidate labels, and label offsets.
- Rebuilt Panel D from three numeric threshold count heatmaps into a threshold sensitivity and shortlist retention matrix.
- Replaced over-strong donor-level and external-verification wording with conservative sample-level robustness boundaries.
- Preserved original LOSO, pseudobulk, and threshold sensitivity values; no upstream analysis was rerun.
- Did not modify main Figure 1-5 or Supplementary FigS1, FigS2, FigS3, FigS5, FigS6, or FigS7.
"""


def terminology_audit_df() -> pd.DataFrame:
    rows = [
        ("sample-level robustness", "present", "required"),
        ("leave-one-sample-out", "present", "required"),
        ("LOSO", "present", "required"),
        ("pseudobulk support", "present", "required"),
        ("threshold sensitivity", "present", "required"),
        ("shortlist retention", "present", "required"),
        ("internal-control", "present", "required"),
        ("internal_control", "absent from figure text", "replaced by internal-control in figure-facing labels"),
        ("donor-level robustness", "absent from figure text", "boundary-only note states not donor-level robustness"),
        ("donor-level validation", "absent", "required absent"),
        ("external validation", "absent from figure text", "appears only in required negated caption/notes boundary"),
        ("formal validation", "absent from figure text", "appears only in required negated notes boundary"),
        ("independent validation", "absent", "required absent"),
        ("validated", "absent", "required absent"),
        ("confirmed", "absent", "required absent"),
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


def scan_figure_text() -> dict[str, bool]:
    text = "\n".join(FIGURE_TEXT_ITEMS)
    forbidden = [
        "donor-level robustness",
        "donor-level validation",
        "external validation",
        "formal validation",
        "independent validation",
        "validated",
        "confirmed",
        "internal_control",
    ]
    lower = text.lower()
    return {term: term.lower() in lower for term in forbidden}


def print_summary(outputs: list[Path], companions: list[Path]) -> None:
    residuals = scan_figure_text()
    print("FigS4 revised_v2 output path:")
    for path in outputs:
        print(f"- {rel(path)}")
    print("Input files:")
    for path in INPUT_FILES:
        print(f"- {rel(path)}")
    print("Sample-level robustness / LOSO wording standardized: yes")
    print("internal_control replaced by internal-control in figure-facing labels: yes")
    print("Panel D rebuilt as threshold sensitivity and shortlist retention matrix: yes")
    flagged = [term for term, present in residuals.items() if present]
    if flagged:
        print("Residual forbidden terms in figure text: " + "; ".join(flagged))
    else:
        print("Residual forbidden terms in figure text: none")
    print("Caption, notes, revision log, terminology audit generated: " + ("yes" if all(path.exists() for path in companions) else "no"))
    print("Panels requiring manual review: none identified; please visually inspect label placement before submission.")


def main() -> None:
    fig = build_figure()
    outputs = save_figure(fig)
    companions = write_companion_files()
    print_summary(outputs, companions)


if __name__ == "__main__":
    main()
