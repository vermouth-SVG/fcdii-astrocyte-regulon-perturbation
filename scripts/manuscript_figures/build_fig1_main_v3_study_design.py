from __future__ import annotations

import textwrap
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_main"
STEM = "Fig1_main_v3_study_design"
CAPTION_TITLE = "Fig. 1. Study design and analytical evidence framework"
H5AD_PATH = ROOT / "final_exports" / "astrocyte_pilot_rna_with_pyscenic_auc.h5ad"

LESION_COLOR = "#9A1F2D"
CONTROL_COLOR = "#0B6670"
BHLHE40_COLOR = "#6A51A3"
SOX2_COLOR = "#8C8C8C"
LESION_FILL = "#F8EDEF"
CONTROL_FILL = "#E8F3F4"
BHLHE40_FILL = "#F1EEF8"
SOX2_FILL = "#F4F4F4"
NEUTRAL_DARK = "#1F2937"
NEUTRAL_MID = "#6B7280"
NEUTRAL_LIGHT = "#E5E7EB"
BACKGROUND = "#FFFFFF"
LIGHT_PANEL = "#F8FAFC"

FIXED_DISPLAY_VALUES = {
    "dataset": "GSE268807 astrocyte pilot",
    "samples": 4,
    "lesion_samples": 2,
    "internal_control_samples": 2,
    "donors": 2,
    "astrocytes": 2322,
    "regulons": 105,
}


WORKFLOW_NODES = [
    ("cells", "Astrocyte\npilot object", ""),
    ("GRN", "pySCENIC\nregulon activity", "GRN -> motif pruning -> AUCell"),
    ("diff", "Differential\nregulon landscape", "Fig. 2"),
    ("Perturbation", "CellOracle\nin silico perturbation", "Fig. 3"),
    ("LOSO", "Sample-level\nrobustness", "LOSO; Fig. 4"),
    ("GO", "Program\ninterpretation", "GO/KEGG + convergence; Fig. 5"),
]

EVIDENCE_LAYERS = [
    ("Discovery", "regulon landscape + TF shortlist"),
    ("Perturbation", "CellOracle in silico perturbation"),
    ("Robustness/support", "sample-level LOSO + supportive datasets"),
    ("Interpretation", "program convergence + working model"),
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.titlesize": 9.0,
            "axes.labelsize": 7.4,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "figure.facecolor": BACKGROUND,
            "savefig.facecolor": BACKGROUND,
            "axes.facecolor": BACKGROUND,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def decode_value(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_h5ad_category(obs_group, field: str) -> list[str]:
    obj = obs_group[field]
    if hasattr(obj, "keys") and "categories" in obj and "codes" in obj:
        categories = [decode_value(value) for value in obj["categories"][:]]
        codes = obj["codes"][:]
        return [categories[int(code)] for code in codes if int(code) >= 0]
    return [decode_value(value) for value in obj[:]]


def read_h5ad_metadata() -> tuple[dict, bool, list[str]]:
    values = dict(FIXED_DISPLAY_VALUES)
    warnings: list[str] = []
    if not H5AD_PATH.exists():
        warnings.append(f"h5ad file not found: {rel(H5AD_PATH)}; fixed displayed values were used.")
        return values, False, warnings

    try:
        import h5py

        with h5py.File(H5AD_PATH, "r") as handle:
            obs = handle["obs"]
            cells = len(obs["_index"])
            genes = len(handle["var"]["_index"]) if "var" in handle and "_index" in handle["var"] else None
            samples = read_h5ad_category(obs, "sample_id") if "sample_id" in obs else []
            groups = read_h5ad_category(obs, "group") if "group" in obs else []
            donors = read_h5ad_category(obs, "donor_id") if "donor_id" in obs else []
            regulons = len(handle["uns"]["pyscenic_regulon_names"]) if "pyscenic_regulon_names" in handle["uns"] else None

            sample_to_groups: dict[str, Counter] = defaultdict(Counter)
            for sample, group in zip(samples, groups):
                sample_to_groups[sample][group] += 1
            sample_group = {
                sample: counts.most_common(1)[0][0]
                for sample, counts in sample_to_groups.items()
                if counts
            }
            sample_group_counts = Counter(sample_group.values())

            values.update(
                {
                    "samples": len(set(samples)) if samples else values["samples"],
                    "lesion_samples": sample_group_counts.get("lesion", values["lesion_samples"]),
                    "internal_control_samples": sample_group_counts.get("internal_control", values["internal_control_samples"]),
                    "donors": len(set(donors)) if donors else values["donors"],
                    "astrocytes": int(cells),
                    "genes": int(genes) if genes is not None else "not available",
                    "regulons": int(regulons) if regulons is not None else values["regulons"],
                    "h5ad_obs_fields": {
                        "sample": "sample_id" if "sample_id" in obs else "missing",
                        "group": "group" if "group" in obs else "missing",
                        "donor": "donor_id" if "donor_id" in obs else "missing",
                    },
                    "sample_group_counts": dict(sample_group_counts),
                    "cell_group_counts": dict(Counter(groups)),
                }
            )
        return values, True, warnings
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"h5ad metadata read failed: {exc}; fixed displayed values were used.")
        return values, False, warnings


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.02,
        1.02,
        label,
        transform=ax.transAxes,
        fontsize=13.0,
        fontweight="bold",
        ha="left",
        va="top",
        color="black",
        clip_on=False,
    )


def rounded_box(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    fc: str,
    ec: str = NEUTRAL_LIGHT,
    lw: float = 0.8,
    radius: float = 0.025,
    alpha: float = 1.0,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        alpha=alpha,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = NEUTRAL_MID, lw: float = 0.9) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=8.0,
            linewidth=lw,
            color=color,
            shrinkA=0,
            shrinkB=0,
            clip_on=False,
        )
    )


def wrap_label(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def draw_astrocyte_icon(ax: plt.Axes, cx: float, cy: float, scale: float = 1.0, color: str = NEUTRAL_DARK) -> None:
    ax.add_patch(Circle((cx, cy), 0.028 * scale, transform=ax.transAxes, facecolor="white", edgecolor=color, linewidth=0.9))
    ax.add_patch(Circle((cx + 0.006 * scale, cy - 0.001 * scale), 0.010 * scale, transform=ax.transAxes, facecolor="#CBD5E1", edgecolor="none", alpha=0.8))
    branches = [
        (-0.070, 0.025),
        (-0.055, -0.040),
        (0.060, 0.040),
        (0.070, -0.025),
        (-0.005, 0.070),
        (0.012, -0.075),
    ]
    for dx, dy in branches:
        ax.plot([cx, cx + dx * scale], [cy, cy + dy * scale], transform=ax.transAxes, color=color, linewidth=0.75, solid_capstyle="round")
        ax.plot(
            [cx + dx * 0.68 * scale, cx + dx * 0.88 * scale],
            [cy + dy * 0.68 * scale, cy + dy * 0.88 * scale + 0.018 * scale],
            transform=ax.transAxes,
            color=color,
            linewidth=0.55,
            solid_capstyle="round",
        )


def draw_panel_a(ax: plt.Axes, values: dict) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "A")
    ax.text(0.05, 0.93, "Discovery object", transform=ax.transAxes, fontsize=9.2, fontweight="bold", color=NEUTRAL_DARK)
    ax.text(0.05, 0.865, "GSE268807 astrocyte pilot", transform=ax.transAxes, fontsize=7.6, color=NEUTRAL_DARK)

    draw_astrocyte_icon(ax, 0.18, 0.72, scale=1.15, color="#475569")
    ax.text(0.32, 0.735, "astrocyte pilot object", transform=ax.transAxes, fontsize=7.3, color=NEUTRAL_DARK, va="center")
    ax.text(0.32, 0.692, f"{values['astrocytes']:,} astrocytes", transform=ax.transAxes, fontsize=6.8, color=NEUTRAL_MID, va="center")

    ax.text(0.05, 0.57, "Samples", transform=ax.transAxes, fontsize=7.0, fontweight="bold", color=NEUTRAL_DARK)
    tile_w, tile_h = 0.155, 0.095
    tiles = [
        ("L1", LESION_COLOR, LESION_FILL, 0.05, 0.435),
        ("L2", LESION_COLOR, LESION_FILL, 0.235, 0.435),
        ("C1", CONTROL_COLOR, CONTROL_FILL, 0.05, 0.300),
        ("C2", CONTROL_COLOR, CONTROL_FILL, 0.235, 0.300),
    ]
    for label, color, fill, x, y in tiles:
        rounded_box(ax, x, y, tile_w, tile_h, fc=fill, ec=color, lw=1.0, radius=0.018)
        ax.text(x + tile_w / 2, y + tile_h / 2, label, transform=ax.transAxes, ha="center", va="center", fontsize=8.0, fontweight="bold", color=color)
    ax.text(0.43, 0.482, "lesion", transform=ax.transAxes, fontsize=6.5, color=LESION_COLOR, va="center")
    ax.text(0.43, 0.347, "internal-control", transform=ax.transAxes, fontsize=6.5, color=CONTROL_COLOR, va="center")

    stat_lines = [
        f"{values['samples']} samples",
        f"{values['lesion_samples']} lesion / {values['internal_control_samples']} internal-control",
        f"{values['donors']} donors",
        f"{values['regulons']} pySCENIC regulons",
    ]
    y = 0.165
    for line in stat_lines:
        ax.add_patch(Circle((0.065, y + 0.005), 0.006, transform=ax.transAxes, facecolor=NEUTRAL_DARK, edgecolor="none"))
        ax.text(0.09, y, line, transform=ax.transAxes, fontsize=6.9, color=NEUTRAL_DARK, va="center")
        y -= 0.048


def draw_node_icon(ax: plt.Axes, x: float, y: float, text: str, edge: str = NEUTRAL_MID) -> None:
    ax.add_patch(Circle((x, y), 0.039, transform=ax.transAxes, facecolor=LIGHT_PANEL, edgecolor=edge, linewidth=0.9, zorder=3))
    ax.text(x, y, text, transform=ax.transAxes, ha="center", va="center", fontsize=5.7, color=edge, fontweight="bold", zorder=4)


def draw_panel_b(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "B")
    ax.text(0.035, 0.93, "Analytical framework", transform=ax.transAxes, fontsize=9.2, fontweight="bold", color=NEUTRAL_DARK)

    xs = [0.07, 0.235, 0.40, 0.565, 0.73, 0.895]
    y = 0.57
    for start, end in zip(xs[:-1], xs[1:]):
        arrow(ax, (start + 0.045, y), (end - 0.047, y), color="#9CA3AF", lw=0.85)

    edge_colors = [NEUTRAL_MID, NEUTRAL_MID, LESION_COLOR, CONTROL_COLOR, NEUTRAL_MID, BHLHE40_COLOR]
    for i, ((icon, label, sub), x, edge) in enumerate(zip(WORKFLOW_NODES, xs, edge_colors)):
        draw_node_icon(ax, x, y, icon, edge=edge)
        label_y = 0.445 if i % 2 == 0 else 0.725
        sub_y = label_y - 0.083 if i % 2 == 0 else label_y + 0.070
        va_label = "top" if i % 2 == 0 else "bottom"
        va_sub = "top" if i % 2 == 0 else "bottom"
        ax.text(x, label_y, label, transform=ax.transAxes, fontsize=6.3, ha="center", va=va_label, color=NEUTRAL_DARK, linespacing=1.08)
        if sub:
            ax.text(x, sub_y, wrap_label(sub, 22), transform=ax.transAxes, fontsize=5.4, ha="center", va=va_sub, color=NEUTRAL_MID, linespacing=1.05)

    ax.plot([0.035, 0.965], [0.17, 0.17], transform=ax.transAxes, color=NEUTRAL_LIGHT, linewidth=0.7)
    ax.text(
        0.035,
        0.105,
        "Dry-lab framework using existing transcriptomic, pySCENIC, CellOracle, and functional interpretation outputs.",
        transform=ax.transAxes,
        fontsize=5.8,
        color=NEUTRAL_MID,
        ha="left",
        va="center",
    )


def draw_panel_c(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "C")
    ax.text(0.05, 0.93, "Evidence hierarchy", transform=ax.transAxes, fontsize=9.2, fontweight="bold", color=NEUTRAL_DARK)

    y_positions = [0.78, 0.64, 0.50, 0.36]
    accent_colors = [LESION_COLOR, CONTROL_COLOR, NEUTRAL_MID, BHLHE40_COLOR]
    for idx, ((heading, sub), y, color) in enumerate(zip(EVIDENCE_LAYERS, y_positions, accent_colors), start=1):
        rounded_box(ax, 0.06, y, 0.88, 0.10, fc=LIGHT_PANEL, ec=NEUTRAL_LIGHT, lw=0.7, radius=0.022)
        ax.plot([0.075, 0.075], [y + 0.018, y + 0.082], transform=ax.transAxes, color=color, linewidth=2.1, solid_capstyle="round")
        ax.text(0.105, y + 0.065, f"{idx}. {heading}", transform=ax.transAxes, fontsize=6.8, fontweight="bold", color=NEUTRAL_DARK, va="center")
        ax.text(0.105, y + 0.033, sub, transform=ax.transAxes, fontsize=5.9, color=NEUTRAL_MID, va="center")
        if idx < len(EVIDENCE_LAYERS):
            arrow(ax, (0.50, y - 0.006), (0.50, y - 0.032), color="#CBD5E1", lw=0.7)

    ax.text(0.06, 0.235, "Final output", transform=ax.transAxes, fontsize=7.0, fontweight="bold", color=NEUTRAL_DARK)
    rounded_box(ax, 0.06, 0.135, 0.41, 0.078, fc=LESION_FILL, ec=LESION_COLOR, lw=0.8, radius=0.02)
    ax.text(0.09, 0.182, "NFE2L2", transform=ax.transAxes, fontsize=7.4, fontweight="bold", color=LESION_COLOR, va="center")
    ax.text(0.09, 0.153, "lesion-associated axis", transform=ax.transAxes, fontsize=5.5, color=NEUTRAL_DARK, va="center")

    rounded_box(ax, 0.53, 0.135, 0.41, 0.078, fc=CONTROL_FILL, ec=CONTROL_COLOR, lw=0.8, radius=0.02)
    ax.text(0.56, 0.182, "THRB", transform=ax.transAxes, fontsize=7.4, fontweight="bold", color=CONTROL_COLOR, va="center")
    ax.text(0.56, 0.153, "internal-control axis", transform=ax.transAxes, fontsize=5.5, color=NEUTRAL_DARK, va="center")

    ax.text(0.085, 0.088, "BHLHE40", transform=ax.transAxes, fontsize=6.3, color=BHLHE40_COLOR, va="center")
    ax.text(0.245, 0.088, "secondary", transform=ax.transAxes, fontsize=5.6, color=NEUTRAL_MID, va="center")
    ax.text(0.555, 0.088, "SOX2", transform=ax.transAxes, fontsize=6.3, color=SOX2_COLOR, va="center")
    ax.text(0.675, 0.088, "retained", transform=ax.transAxes, fontsize=5.6, color=NEUTRAL_MID, va="center")

    ax.text(0.06, 0.025, "External datasets were used as supportive evidence only.", transform=ax.transAxes, fontsize=5.4, color=NEUTRAL_MID)


def save_outputs(fig: plt.Figure) -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {
        "png": OUT_DIR / f"{STEM}.png",
        "pdf": OUT_DIR / f"{STEM}.pdf",
        "svg": OUT_DIR / f"{STEM}.svg",
        "notes": OUT_DIR / f"{STEM}_notes.txt",
    }
    fig.savefig(outputs["png"], dpi=600, bbox_inches="tight")
    fig.savefig(outputs["pdf"], bbox_inches="tight")
    fig.savefig(outputs["svg"], bbox_inches="tight")
    plt.close(fig)
    return outputs


def write_notes(outputs: dict[str, Path], values: dict, h5ad_read_success: bool, warnings: list[str]) -> None:
    panel_a_values = {
        "dataset": values["dataset"],
        "samples": values["samples"],
        "lesion_samples": values["lesion_samples"],
        "internal_control_samples": values["internal_control_samples"],
        "donors": values["donors"],
        "astrocytes": values["astrocytes"],
        "pySCENIC_regulons": values["regulons"],
    }
    notes = f"""
    {STEM}
    ==============================

    Caption title:
    {CAPTION_TITLE}

    Version note:
    v3 is a layout-level reconstruction of Figure 1 as a publication-style vector workflow.
    No old Figure1 image was modified and no generative image model was used.

    Upstream analysis:
    No pySCENIC, CellOracle, enrichment, drug repositioning, or other upstream analysis was rerun.

    h5ad metadata:
    - h5ad path: {rel(H5AD_PATH)}
    - h5ad read success: {h5ad_read_success}
    - cells from h5ad if available: {values.get("astrocytes", "not available")}
    - genes from h5ad if available: {values.get("genes", "not available")}
    - samples from h5ad if available: {values.get("samples", "not available")}
    - pySCENIC regulons from h5ad if available: {values.get("regulons", "not available")}
    - obs fields used if available: {values.get("h5ad_obs_fields", "fixed displayed values")}
    - sample group counts from h5ad if available: {values.get("sample_group_counts", "not available")}
    - cell group counts from h5ad if available: {values.get("cell_group_counts", "not available")}

    Panel A displayed cohort/object values:
    {textwrap.indent(str(panel_a_values), "    ")}
    {"Fixed displayed values from current project notes were used because h5ad metadata could not be read." if not h5ad_read_success else "Displayed values were confirmed from h5ad metadata using h5py."}

    Panel B workflow nodes:
    {textwrap.indent(chr(10).join([f"- {label.replace(chr(10), ' ')} | {sub}" for _, label, sub in WORKFLOW_NODES]), "    ")}

    Panel C evidence hierarchy labels:
    {textwrap.indent(chr(10).join([f"- {heading}: {sub}" for heading, sub in EVIDENCE_LAYERS]), "    ")}
    Final output labels:
    - NFE2L2: lesion-associated axis
    - THRB: internal-control axis
    - BHLHE40: secondary
    - SOX2: retained

    Required interpretation notes:
    Drug repositioning clues were not included in main Figure 1 and remain supplementary/exploratory.
    External datasets are presented as supportive evidence only, not formal validation.
    Robustness is presented as sample-level/LOSO robustness, not donor-level robustness.

    Warnings / missing data:
    {textwrap.indent(chr(10).join(warnings) if warnings else "None.", "    ")}

    Output files:
    - {rel(outputs["png"])}
    - {rel(outputs["pdf"])}
    - {rel(outputs["svg"])}
    - {rel(outputs["notes"])}
    """
    outputs["notes"].write_text(textwrap.dedent(notes).strip() + "\n", encoding="utf-8")


def main() -> None:
    setup_style()
    values, h5ad_read_success, warnings = read_h5ad_metadata()

    fig = plt.figure(figsize=(14.5, 7.2), constrained_layout=False)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.28, 0.44, 0.28], wspace=0.12)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    draw_panel_a(ax_a, values)
    draw_panel_b(ax_b)
    draw_panel_c(ax_c)

    outputs = save_outputs(fig)
    write_notes(outputs, values, h5ad_read_success, warnings)

    for key in ["png", "pdf", "svg", "notes"]:
        print(f"Wrote {rel(outputs[key])}")
    print(f"h5ad read success: {h5ad_read_success}")
    if not h5ad_read_success:
        print("Used fixed displayed values from current project notes.")
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"- {warning}")
    else:
        print("Warnings: none")


if __name__ == "__main__":
    main()
