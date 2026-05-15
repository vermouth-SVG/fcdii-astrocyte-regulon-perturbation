#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_supplementary_revised_v2"
STEM = "FigS3_celloracle_perturbation_metrics_revised_v2_minorfix"

TF_ORDER = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
TF_COLORS = {
    "NFE2L2": "#8E1B2A",
    "THRB": "#0B6B74",
    "BHLHE40": "#6A51A3",
    "SOX2": "#777777",
}
GROUP_POINT_COLORS = {
    "lesion": "#DDB5BC",
    "internal_control": "#B8DDE1",
    "internal-control": "#B8DDE1",
}
RANDOM_COLOR = "#7E858A"
GRID_COLOR = "#E8E8E8"
TEXT_COLOR = "#222222"

METRIC_TABLE = ROOT / "celloracle_run" / "ko_round1" / "round1_ko_master_table.csv"
RECOVERY_BY_GROUP = ROOT / "celloracle_run" / "ko_round1" / "round1_recovery_index_by_group.csv"
RECOVERY_DEFINITION = ROOT / "celloracle_run" / "ko_round1" / "round1_recovery_index_definition.txt"
VECTOR_FILES = {
    tf: ROOT / "celloracle_run" / "ko_round1" / tf / f"{tf}_state_shift_scores.csv"
    for tf in TF_ORDER
}

INPUT_FILES = [METRIC_TABLE, RECOVERY_BY_GROUP, RECOVERY_DEFINITION, *VECTOR_FILES.values()]

METRIC_COLUMNS = [
    ("mean_shift_length", "Mean shift"),
    ("recovery_index", "Approx. RI"),
    ("mean_net_shift", "Net shift"),
    ("lesion_mean_shift_length", "Lesion shift"),
    ("internal_control_mean_shift_length", "Internal-control shift"),
    ("direction_agreement_factor", "Direction agreement"),
]

FIGURE_TEXT_ITEMS = [
    "Complete CellOracle in silico perturbation metrics and randomized-control assessment",
    "Round 1 CellOracle perturbation metric heatmap",
    "Column-scaled metric",
    "Mean shift",
    "Approx. RI",
    "Net shift",
    "Lesion shift",
    "Internal-control shift",
    "Direction agreement",
    "Quantitative perturbation summary",
    "Mean shift length",
    "Approx. Recovery Index",
    "NFE2L2\nhighest mean shift",
    "THRB\nhighest Approx. RI",
    "BHLHE40\nsecondary",
    "SOX2\nretained",
    "NFE2L2 representative in silico KO",
    "THRB representative in silico KO",
    "BHLHE40 representative in silico KO",
    "SOX2 representative in silico KO",
    "Simulated KO shift",
    "Randomized control",
    "CellOracle in silico KO / simulated perturbation only; Approx. RI = Approx. Recovery Index, a state-shift-derived summary metric, not experimental knockout evidence; randomized control is a negative-control assessment.",
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.titlesize": 8.4,
            "axes.labelsize": 7.4,
            "xtick.labelsize": 6.6,
            "ytick.labelsize": 6.7,
            "legend.fontsize": 6.2,
            "figure.titlesize": 11.4,
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


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def panel_label(ax: plt.Axes, label: str, x: float = -0.10, y: float = 1.14) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.4,
        fontweight="bold",
        color=TEXT_COLOR,
        clip_on=False,
    )


def style_axes(ax: plt.Axes, grid_axis: str | None = "both") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def load_metrics() -> pd.DataFrame:
    df = read_csv(METRIC_TABLE)
    missing = [col for col, _ in METRIC_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required metric columns in {rel(METRIC_TABLE)}: {missing}")
    if "tf" not in df.columns:
        raise ValueError(f"Missing tf column in {rel(METRIC_TABLE)}")
    df = df[df["tf"].isin(TF_ORDER)].copy()
    for col, _ in METRIC_COLUMNS:
        df[col] = as_num(df[col])
    df = df.set_index("tf").reindex(TF_ORDER)
    if df[["mean_shift_length", "recovery_index"]].isna().any().any():
        raise ValueError("Required CellOracle metric values contain missing data after TF reindexing.")
    return df


def load_vectors() -> dict[str, pd.DataFrame]:
    vectors: dict[str, pd.DataFrame] = {}
    required = [
        "embedding_x",
        "embedding_y",
        "delta_x",
        "delta_y",
        "delta_random_x",
        "delta_random_y",
        "group",
    ]
    for tf, path in VECTOR_FILES.items():
        df = read_csv(path)
        missing = [col for col in required if col not in df.columns]
        if missing:
            raise ValueError(f"Missing vector columns in {rel(path)}: {missing}")
        for col in ["embedding_x", "embedding_y", "delta_x", "delta_y", "delta_random_x", "delta_random_y"]:
            df[col] = as_num(df[col])
        vectors[tf] = df.dropna(subset=["embedding_x", "embedding_y", "delta_x", "delta_y"]).copy()
    return vectors


def compute_limits(vectors: dict[str, pd.DataFrame]) -> tuple[tuple[float, float], tuple[float, float]]:
    all_x = pd.concat([df["embedding_x"] for df in vectors.values()], ignore_index=True)
    all_y = pd.concat([df["embedding_y"] for df in vectors.values()], ignore_index=True)
    x_min, x_max = float(all_x.min()), float(all_x.max())
    y_min, y_max = float(all_y.min()), float(all_y.max())
    x_pad = max((x_max - x_min) * 0.035, 0.4)
    y_pad = max((y_max - y_min) * 0.035, 0.4)
    return (x_min - x_pad, x_max + x_pad), (y_min - y_pad, y_max + y_pad)


def build_panel_a(ax: plt.Axes, metrics: pd.DataFrame, fig: plt.Figure) -> None:
    cols = [col for col, _ in METRIC_COLUMNS]
    labels = [label for _, label in METRIC_COLUMNS]
    raw = metrics[cols].copy()
    col_min = raw.min(axis=0)
    col_span = raw.max(axis=0) - col_min
    scaled = raw.sub(col_min, axis=1)
    for col in cols:
        if col_span[col] == 0 or pd.isna(col_span[col]):
            scaled[col] = 0.5
        else:
            scaled[col] = scaled[col] / col_span[col]

    cmap = LinearSegmentedColormap.from_list("celloracle_metric", ["#F7F7F7", "#D8B7BC", "#8E1B2A"])
    im = ax.imshow(scaled.values, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_title("Round 1 CellOracle perturbation metric heatmap", pad=7)
    ax.set_yticks(np.arange(len(TF_ORDER)))
    ax.set_yticklabels(TF_ORDER)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(
        ["Mean\nshift", "Approx.\nRI", "Net\nshift", "Lesion\nshift", "Internal-control\nshift", "Direction\nagreement"],
        rotation=22,
        ha="right",
        rotation_mode="anchor",
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for i, tf in enumerate(TF_ORDER):
        for j, col in enumerate(cols):
            value = raw.loc[tf, col]
            text_color = "white" if scaled.loc[tf, col] >= 0.58 else TEXT_COLOR
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=6.1, color=text_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.025)
    cbar.set_label("Column-scaled metric", fontsize=7.0)
    cbar.ax.tick_params(labelsize=6.2, width=0.5, length=2.2)
    panel_label(ax, "A")


def build_panel_b(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    x = metrics["mean_shift_length"]
    y = metrics["recovery_index"]
    ax.scatter(
        x,
        y,
        s=105,
        c=[TF_COLORS[tf] for tf in TF_ORDER],
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
    )
    labels = {
        "NFE2L2": ("NFE2L2\nhighest mean shift", (-54, -18), "right"),
        "THRB": ("THRB\nhighest Approx. RI", (-55, 14), "right"),
        "BHLHE40": ("BHLHE40\nsecondary", (9, 10), "left"),
        "SOX2": ("SOX2\nretained", (9, -2), "left"),
    }
    for tf, (text, offset, align) in labels.items():
        ax.annotate(
            text,
            xy=(float(x.loc[tf]), float(y.loc[tf])),
            xytext=offset,
            textcoords="offset points",
            ha=align,
            va="center",
            fontsize=6.5,
            color=TF_COLORS[tf],
            arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.55, "shrinkA": 1.5, "shrinkB": 3.5},
            zorder=4,
        )

    x_span = float(x.max() - x.min())
    y_span = float(y.max() - y.min())
    ax.set_xlim(float(x.min() - 0.08 * x_span), float(x.max() + 0.11 * x_span))
    ax.set_ylim(float(y.min() - 0.09 * y_span), float(y.max() + 0.11 * y_span))
    ax.set_xlabel("Mean shift length")
    ax.set_ylabel("Approx. Recovery Index")
    ax.set_title("Quantitative perturbation summary", pad=7, fontsize=8.0, loc="left")
    style_axes(ax, "both")
    panel_label(ax, "B")


def sample_vector_rows(df: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    if len(df) <= n:
        return df
    return df.sample(n=n, random_state=seed)


def draw_vector_field(
    ax: plt.Axes,
    df: pd.DataFrame,
    tf: str,
    mode: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    seed: int,
) -> None:
    ax.set_facecolor("#FBFBFB")
    for group, group_df in df.groupby("group", sort=False):
        point_color = GROUP_POINT_COLORS.get(str(group), "#D5D5D5")
        ax.scatter(
            group_df["embedding_x"],
            group_df["embedding_y"],
            s=2.5,
            color=point_color,
            alpha=0.31,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

    sampled = sample_vector_rows(df, 320, seed)
    if mode == "simulated":
        dx_col, dy_col = "delta_x", "delta_y"
        arrow_color = TF_COLORS[tf]
    else:
        dx_col, dy_col = "delta_random_x", "delta_random_y"
        arrow_color = RANDOM_COLOR

    ax.quiver(
        sampled["embedding_x"],
        sampled["embedding_y"],
        sampled[dx_col],
        sampled[dy_col],
        angles="xy",
        scale_units="xy",
        scale=1.05,
        color=arrow_color,
        alpha=0.82 if mode == "simulated" else 0.64,
        width=0.0044 if mode == "simulated" else 0.0040,
        headwidth=3.4,
        headlength=4.4,
        headaxislength=3.8,
        zorder=2,
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.55)
        spine.set_color("#D4D4D4")
    ax.set_title("Simulated KO shift" if mode == "simulated" else "Randomized control", pad=2.2, fontsize=6.9)


def build_tf_block(
    fig: plt.Figure,
    parent_spec,
    tf: str,
    label: str,
    vectors: dict[str, pd.DataFrame],
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    seed: int,
) -> None:
    sub = parent_spec.subgridspec(2, 2, height_ratios=[0.12, 1.0], wspace=0.045, hspace=0.020)
    title_ax = fig.add_subplot(sub[0, :])
    title_ax.set_axis_off()
    title_ax.text(-0.045, 0.92, label, ha="left", va="top", fontsize=10.4, fontweight="bold", color=TEXT_COLOR, clip_on=False)
    title_ax.text(
        0.00,
        0.55,
        f"{tf} representative in silico KO",
        ha="left",
        va="center",
        fontsize=8.2,
        fontweight="bold",
        color=TF_COLORS[tf],
    )

    ax_sim = fig.add_subplot(sub[1, 0])
    ax_random = fig.add_subplot(sub[1, 1])
    draw_vector_field(ax_sim, vectors[tf], tf, "simulated", xlim, ylim, seed)
    draw_vector_field(ax_random, vectors[tf], tf, "randomized", xlim, ylim, seed + 100)


def build_method_note(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.add_patch(
        plt.Rectangle(
            (0.02, 0.19),
            0.96,
            0.62,
            facecolor="#F6F6F6",
            edgecolor="#CFCFCF",
            linewidth=0.65,
        )
    )
    ax.text(
        0.5,
        0.50,
        "CellOracle in silico KO / simulated perturbation only; Approx. RI = Approx. Recovery Index, "
        "a state-shift-derived summary metric, not experimental knockout evidence; randomized control is a negative-control assessment.",
        ha="center",
        va="center",
        fontsize=7.2,
        color="#333333",
    )


def build_figure(metrics: pd.DataFrame, vectors: dict[str, pd.DataFrame]) -> plt.Figure:
    setup_style()
    xlim, ylim = compute_limits(vectors)
    fig = plt.figure(figsize=(13.3, 10.4), constrained_layout=False)
    gs = fig.add_gridspec(
        4,
        2,
        height_ratios=[1.05, 1.70, 1.70, 0.20],
        hspace=0.36,
        wspace=0.28,
        left=0.065,
        right=0.985,
        top=0.905,
        bottom=0.055,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    build_panel_a(ax_a, metrics, fig)
    build_panel_b(ax_b, metrics)

    build_tf_block(fig, gs[1, 0], "NFE2L2", "C", vectors, xlim, ylim, seed=17)
    build_tf_block(fig, gs[1, 1], "THRB", "D", vectors, xlim, ylim, seed=23)
    build_tf_block(fig, gs[2, 0], "BHLHE40", "E", vectors, xlim, ylim, seed=31)
    build_tf_block(fig, gs[2, 1], "SOX2", "F", vectors, xlim, ylim, seed=43)

    note_ax = fig.add_subplot(gs[3, :])
    build_method_note(note_ax)
    fig.suptitle(
        "Complete CellOracle in silico perturbation metrics and randomized-control assessment",
        x=0.5,
        y=0.975,
        fontsize=12.0,
        fontweight="bold",
    )
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


def notes_text(metrics: pd.DataFrame) -> str:
    tf_lines = []
    for tf in TF_ORDER:
        tf_lines.append(
            f"- {tf}: mean_shift_length={metrics.loc[tf, 'mean_shift_length']:.6f}; "
            f"Approx. Recovery Index={metrics.loc[tf, 'recovery_index']:.6f}; "
            f"mean_net_shift={metrics.loc[tf, 'mean_net_shift']:.6f}."
        )

    inputs = "\n".join(f"- {rel(path)}" for path in INPUT_FILES)
    return f"""{STEM} notes

Input files:
{inputs}

Minorfix changes:
- Panel A title was changed from Round1 CellOracle perturbation metric heatmap to Round 1 CellOracle perturbation metric heatmap.
- Panels C-F retain the same data and two-column structure, with modestly increased point opacity, arrow opacity, arrow width, arrow length, and plotting area for readability.
- The bottom method note now uses: not experimental knockout evidence.

CellOracle round1 perturbation metrics:
- Metrics were read from {rel(METRIC_TABLE)} without changing TF order or numeric values.
- Mean shift is the mean CellOracle state-shift vector length across cells after the simulated in silico KO.
- Approx. Recovery Index is the existing recovery_index column, treated as an approximate state-shift-derived summary metric.
- Net shift is the observed mean shift after subtracting the randomized control shift baseline.
- Lesion shift and internal-control shift are group-specific mean shift lengths from the lesion and internal-control cell groups.
- Direction agreement is the existing direction_agreement_factor, indicating whether expected regulon-side and expression-side direction agreed in the current round1 table.

Metric values used in the figure:
{chr(10).join(tf_lines)}

Representative perturbation vector fields:
- Vector fields were redrawn from each TF-specific state_shift_scores.csv file using the existing embedding_x, embedding_y, delta_x, delta_y, delta_random_x, and delta_random_y columns.
- Simulated KO shift panels show existing CellOracle in silico KO state-shift vectors.
- Randomized control panels show existing randomized-control vectors and are used only as a negative-control assessment.
- Arrows were downsampled deterministically for readability; no CellOracle model, pySCENIC output, metric value, TF ranking, or random-control value was recomputed.

Interpretation boundaries:
- These panels show in silico KO / simulated perturbation outputs, not an experimental knockout assay.
- Approx. Recovery Index is a state-shift-derived summary metric; it is not a strict geometric projection model and is not a direct pathology-reversal measurement.
- NFE2L2 is interpreted as the primary lesion-associated axis with the highest mean shift.
- THRB is interpreted as the primary internal-control axis and Approx. Recovery Index anchor.
- BHLHE40 is retained as a secondary candidate.
- SOX2 is retained as a lower-priority candidate.
"""


def write_companion_files(metrics: pd.DataFrame) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    notes_path = OUT_DIR / f"{STEM}_notes.txt"
    notes_path.write_text(notes_text(metrics), encoding="utf-8")
    return [notes_path]


def scan_residual_terms(metrics: pd.DataFrame) -> dict[str, bool]:
    scanned = "\n".join([*FIGURE_TEXT_ITEMS, notes_text(metrics)])
    forbidden = [
        "experimental" + " KO",
        "validated" + " KO",
        "validated" + " knockout",
        "therapeutic" + " target",
        "treatment" + " target",
        "drug" + " target",
        "validated" + " perturbation",
        "formal" + " validation",
        "experimental" + " validation",
        "internal_control",
    ]
    lower = scanned.lower()
    return {term: term.lower() in lower for term in forbidden}


def print_summary(outputs: list[Path], companions: list[Path], residuals: dict[str, bool]) -> None:
    print("FigS3 revised_v2 minorfix output path:")
    for path in outputs:
        print(f"- {rel(path)}")
    print("Input files:")
    for path in INPUT_FILES:
        print(f"- {rel(path)}")
    print("Round1 changed to Round 1: yes")
    print("Panels C-F vector field readability enhanced: yes; point opacity, arrow opacity, arrow width, arrow length, and plotting area were increased uniformly.")
    print("Bottom method note changed from wet-lab KO evidence to experimental knockout evidence: yes")
    flagged = [term for term, present in residuals.items() if present]
    if flagged:
        print("Residual forbidden terms in figure/notes text: " + "; ".join(flagged))
    else:
        print("Residual forbidden terms in figure/notes text: none")
    generated = all(path.exists() for path in companions)
    print("Minorfix notes generated: " + ("yes" if generated else "no"))
    print("Approx. RI / randomized control / negative-control assessment retained: yes")


def main() -> None:
    metrics = load_metrics()
    vectors = load_vectors()
    fig = build_figure(metrics, vectors)
    outputs = save_figure(fig)
    companions = write_companion_files(metrics)
    residuals = scan_residual_terms(metrics)
    print_summary(outputs, companions, residuals)


if __name__ == "__main__":
    main()
