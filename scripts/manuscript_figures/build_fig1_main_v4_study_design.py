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
STEM = "Fig1_main_v4_study_design"
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
NEUTRAL_GREY = "#6B7280"
NEUTRAL_LIGHT = "#E5E7EB"
PALE_BACKGROUND = "#F8FAFC"
BACKGROUND = "#FFFFFF"

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
    ("Astrocyte pilot\nobject", ""),
    ("pySCENIC regulon\nactivity", "GRN -> motif pruning -> AUCell"),
    ("Differential regulon\nlandscape / TF\nprioritization", "Fig. 2"),
    ("CellOracle\nin silico KO", "Fig. 3"),
    ("Sample-level robustness\n+ supportive evidence", "LOSO; Fig. 4"),
    ("Program interpretation\n/ working model", "GO/KEGG + convergence; Fig. 5"),
]

EVIDENCE_LAYERS = [
    ("Discovery evidence", "regulon landscape + TF shortlist"),
    ("Perturbation evidence", "CellOracle in silico KO"),
    ("Robustness/support", "sample-level LOSO + supportive datasets"),
    ("Program interpretation", "convergence + working model"),
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.1,
            "axes.titlesize": 9.0,
            "axes.labelsize": 7.2,
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
                    "internal_control_samples": sample_group_counts.get(
                        "internal_control", values["internal_control_samples"]
                    ),
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
        0.0,
        1.015,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="left",
        va="bottom",
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
    lw: float = 0.75,
    radius: float = 0.02,
    alpha: float = 1.0,
    zorder: int = 1,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.010,rounding_size={radius}",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        alpha=alpha,
        transform=ax.transAxes,
        clip_on=False,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = NEUTRAL_GREY,
    lw: float = 0.8,
    style: str = "-|>",
    scale: float = 7.5,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle=style,
            mutation_scale=scale,
            linewidth=lw,
            color=color,
            shrinkA=0,
            shrinkB=0,
            clip_on=False,
            zorder=4,
        )
    )


def draw_astrocyte_icon(ax: plt.Axes, cx: float, cy: float, scale: float = 1.0) -> None:
    color = "#475569"
    ax.add_patch(Circle((cx, cy), 0.020 * scale, transform=ax.transAxes, facecolor="white", edgecolor=color, linewidth=0.8))
    ax.add_patch(Circle((cx + 0.004 * scale, cy), 0.007 * scale, transform=ax.transAxes, facecolor="#CBD5E1", edgecolor="none"))
    for dx, dy in [(-0.050, 0.025), (-0.044, -0.034), (0.050, 0.032), (0.052, -0.025), (-0.003, 0.052), (0.010, -0.055)]:
        ax.plot([cx, cx + dx * scale], [cy, cy + dy * scale], transform=ax.transAxes, color=color, linewidth=0.58)


def draw_panel_a(ax: plt.Axes, values: dict) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "A")
    ax.text(0.05, 0.945, "Discovery object", transform=ax.transAxes, fontsize=9.2, fontweight="bold", color=NEUTRAL_DARK)

    rounded_box(ax, 0.05, 0.815, 0.90, 0.105, fc=PALE_BACKGROUND, ec=NEUTRAL_LIGHT, lw=0.7, radius=0.022)
    draw_astrocyte_icon(ax, 0.105, 0.868, scale=0.90)
    ax.text(0.165, 0.888, "GSE268807 astrocyte pilot", transform=ax.transAxes, fontsize=7.2, fontweight="bold", color=NEUTRAL_DARK, va="center")
    ax.text(0.165, 0.848, "astrocyte pilot object", transform=ax.transAxes, fontsize=6.1, color=NEUTRAL_GREY, va="center")

    ax.text(0.05, 0.740, "Samples", transform=ax.transAxes, fontsize=6.7, fontweight="bold", color=NEUTRAL_DARK)
    tile_w, tile_h = 0.145, 0.100
    tiles = [
        ("L1", LESION_COLOR, LESION_FILL, 0.08, 0.595),
        ("L2", LESION_COLOR, LESION_FILL, 0.245, 0.595),
        ("C1", CONTROL_COLOR, CONTROL_FILL, 0.08, 0.455),
        ("C2", CONTROL_COLOR, CONTROL_FILL, 0.245, 0.455),
    ]
    for label, color, fill, x, y in tiles:
        rounded_box(ax, x, y, tile_w, tile_h, fc=fill, ec=color, lw=0.95, radius=0.018)
        ax.text(x + tile_w / 2, y + tile_h / 2, label, transform=ax.transAxes, ha="center", va="center", fontsize=8.0, fontweight="bold", color=color)
    ax.plot([0.46, 0.46], [0.454, 0.695], transform=ax.transAxes, color=NEUTRAL_LIGHT, linewidth=0.8)
    ax.text(0.50, 0.648, "lesion", transform=ax.transAxes, fontsize=6.4, color=LESION_COLOR, va="center")
    ax.text(0.50, 0.508, "internal-control", transform=ax.transAxes, fontsize=6.4, color=CONTROL_COLOR, va="center")

    metrics = [
        (f"{values['samples']}", "samples"),
        (f"{values['donors']}", "donors"),
        (f"{values['lesion_samples']}", "lesion"),
        (f"{values['internal_control_samples']}", "internal-control"),
        (f"{values['astrocytes']:,}", "astrocytes"),
        (f"{values['regulons']}", "regulons"),
    ]
    x0s = [0.05, 0.355, 0.660]
    y0s = [0.260, 0.125]
    k = 0
    for y in y0s:
        for x in x0s:
            number, label = metrics[k]
            rounded_box(ax, x, y, 0.25, 0.088, fc="white", ec=NEUTRAL_LIGHT, lw=0.7, radius=0.017)
            ax.text(x + 0.018, y + 0.055, number, transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=NEUTRAL_DARK, va="center")
            ax.text(x + 0.018, y + 0.026, label, transform=ax.transAxes, fontsize=5.5, color=NEUTRAL_GREY, va="center")
            k += 1


def draw_panel_b(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "B")
    ax.text(0.035, 0.945, "Analysis-to-evidence framework", transform=ax.transAxes, fontsize=9.2, fontweight="bold", color=NEUTRAL_DARK)

    xs = [0.070, 0.238, 0.406, 0.574, 0.742, 0.910]
    y = 0.555
    node_w = 0.128
    node_h = 0.178
    for start, end in zip(xs[:-1], xs[1:]):
        arrow(ax, (start + node_w / 2 + 0.012, y), (end - node_w / 2 - 0.012, y), color="#9CA3AF", lw=0.75, scale=7.0)

    edge_colors = [NEUTRAL_GREY, NEUTRAL_GREY, LESION_COLOR, CONTROL_COLOR, NEUTRAL_GREY, BHLHE40_COLOR]
    for idx, ((label, sub), x, edge) in enumerate(zip(WORKFLOW_NODES, xs, edge_colors), start=1):
        left = x - node_w / 2
        rounded_box(ax, left, y - node_h / 2, node_w, node_h, fc=PALE_BACKGROUND, ec=NEUTRAL_LIGHT, lw=0.7, radius=0.019)
        ax.plot([left + 0.014, left + node_w - 0.014], [y + node_h / 2 - 0.014, y + node_h / 2 - 0.014], transform=ax.transAxes, color=edge, linewidth=1.7, solid_capstyle="round")
        ax.add_patch(Circle((left + 0.024, y + node_h / 2 - 0.036), 0.013, transform=ax.transAxes, facecolor="white", edgecolor=edge, linewidth=0.8))
        ax.text(left + 0.024, y + node_h / 2 - 0.036, str(idx), transform=ax.transAxes, fontsize=4.8, fontweight="bold", color=edge, ha="center", va="center")
        ax.text(x, y + 0.018, label, transform=ax.transAxes, fontsize=5.55, color=NEUTRAL_DARK, fontweight="bold", ha="center", va="center", linespacing=1.08)
        if sub:
            ax.text(x, y - node_h / 2 - 0.040, textwrap.fill(sub, width=18), transform=ax.transAxes, fontsize=5.05, color=NEUTRAL_GREY, ha="center", va="top", linespacing=1.06)

    rounded_box(ax, 0.055, 0.155, 0.89, 0.105, fc="white", ec=NEUTRAL_LIGHT, lw=0.65, radius=0.020)
    ax.text(0.083, 0.217, "Framework scope", transform=ax.transAxes, fontsize=5.9, fontweight="bold", color=NEUTRAL_DARK, va="center")
    ax.text(
        0.083,
        0.181,
        "Existing transcriptomic, pySCENIC/AUCell, CellOracle, robustness, and functional interpretation outputs only.",
        transform=ax.transAxes,
        fontsize=5.3,
        color=NEUTRAL_GREY,
        va="center",
    )


def draw_panel_c(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "C")
    ax.text(0.05, 0.945, "Integrated output", transform=ax.transAxes, fontsize=9.2, fontweight="bold", color=NEUTRAL_DARK)

    ax.text(0.06, 0.870, "Evidence layers", transform=ax.transAxes, fontsize=6.5, fontweight="bold", color=NEUTRAL_DARK)
    y_positions = [0.760, 0.650, 0.540, 0.430]
    accent_colors = [LESION_COLOR, CONTROL_COLOR, NEUTRAL_GREY, BHLHE40_COLOR]
    for idx, ((heading, sub), y, color) in enumerate(zip(EVIDENCE_LAYERS, y_positions, accent_colors), start=1):
        rounded_box(ax, 0.06, y, 0.88, 0.078, fc=PALE_BACKGROUND, ec=NEUTRAL_LIGHT, lw=0.65, radius=0.018)
        ax.add_patch(Circle((0.095, y + 0.039), 0.018, transform=ax.transAxes, facecolor="white", edgecolor=color, linewidth=0.8))
        ax.text(0.095, y + 0.039, str(idx), transform=ax.transAxes, fontsize=5.0, fontweight="bold", color=color, ha="center", va="center")
        ax.text(0.135, y + 0.051, heading, transform=ax.transAxes, fontsize=6.0, fontweight="bold", color=NEUTRAL_DARK, va="center")
        ax.text(0.135, y + 0.024, sub, transform=ax.transAxes, fontsize=5.15, color=NEUTRAL_GREY, va="center")
        if idx < len(EVIDENCE_LAYERS):
            ax.plot([0.095, 0.095], [y - 0.010, y - 0.032], transform=ax.transAxes, color="#CBD5E1", linewidth=0.8)

    arrow(ax, (0.50, 0.400), (0.50, 0.348), color="#CBD5E1", lw=0.75, scale=7.5)
    ax.text(0.06, 0.332, "Final dual-axis output", transform=ax.transAxes, fontsize=6.5, fontweight="bold", color=NEUTRAL_DARK)

    rounded_box(ax, 0.075, 0.205, 0.36, 0.105, fc=LESION_FILL, ec=LESION_COLOR, lw=0.9, radius=0.022)
    rounded_box(ax, 0.565, 0.205, 0.36, 0.105, fc=CONTROL_FILL, ec=CONTROL_COLOR, lw=0.9, radius=0.022)
    ax.text(0.105, 0.270, "NFE2L2", transform=ax.transAxes, fontsize=7.9, fontweight="bold", color=LESION_COLOR, va="center")
    ax.text(0.105, 0.235, "lesion-associated axis", transform=ax.transAxes, fontsize=5.55, color=NEUTRAL_DARK, va="center")
    ax.text(0.595, 0.270, "THRB", transform=ax.transAxes, fontsize=7.9, fontweight="bold", color=CONTROL_COLOR, va="center")
    ax.text(0.595, 0.235, "internal-control axis", transform=ax.transAxes, fontsize=5.55, color=NEUTRAL_DARK, va="center")
    arrow(ax, (0.455, 0.258), (0.545, 0.258), color=NEUTRAL_GREY, lw=0.9, style="<->", scale=8.0)

    ax.text(0.090, 0.145, "BHLHE40", transform=ax.transAxes, fontsize=5.9, fontweight="bold", color=BHLHE40_COLOR, va="center")
    ax.text(0.238, 0.145, "secondary", transform=ax.transAxes, fontsize=5.2, color=NEUTRAL_GREY, va="center")
    ax.text(0.590, 0.145, "SOX2", transform=ax.transAxes, fontsize=5.9, fontweight="bold", color=SOX2_COLOR, va="center")
    ax.text(0.698, 0.145, "retained", transform=ax.transAxes, fontsize=5.2, color=NEUTRAL_GREY, va="center")

    ax.plot([0.06, 0.94], [0.080, 0.080], transform=ax.transAxes, color=NEUTRAL_LIGHT, linewidth=0.7)
    ax.text(0.06, 0.047, "External datasets were used as supportive evidence only.", transform=ax.transAxes, fontsize=5.0, color=NEUTRAL_GREY)


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
    v4 is a layout-level reconstruction, not an upstream analysis update.
    The script was derived from {rel(ROOT / "scripts" / "manuscript_figures" / "build_fig1_main_v3_study_design.py")} and redrawn as a compact three-panel publication-style workflow.
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

    Panel A displayed object values:
    {textwrap.indent(str(panel_a_values), "    ")}
    {"Fixed displayed values from current project notes were used because h5ad metadata could not be read." if not h5ad_read_success else "Displayed values were confirmed from h5ad metadata using h5py."}

    Panel B workflow nodes:
    {textwrap.indent(chr(10).join([f"- {label.replace(chr(10), ' ')} | {sub}" for label, sub in WORKFLOW_NODES]), "    ")}

    Panel C evidence layers and final output labels:
    {textwrap.indent(chr(10).join([f"- {heading}: {sub}" for heading, sub in EVIDENCE_LAYERS]), "    ")}
    Final output labels:
    - NFE2L2: lesion-associated axis
    - THRB: internal-control axis
    - BHLHE40: secondary
    - SOX2: retained

    Required interpretation notes:
    Drug repositioning clues were not included in main Figure 1.
    External datasets were used as supportive evidence only, not formal validation.
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
    gs = fig.add_gridspec(1, 3, width_ratios=[0.27, 0.46, 0.27], wspace=0.085)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    fig.subplots_adjust(left=0.035, right=0.985, top=0.930, bottom=0.075)

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
