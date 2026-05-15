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
OUT_DIR = ROOT / "manuscript_output" / "figures_main"
STEM = "Fig3_main_v3_celloracle_perturbation"
CAPTION_TITLE = "Fig. 3. In silico perturbation prioritizes NFE2L2 and THRB as complementary regulatory anchors"

MAIN_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]

LESION_COLOR = "#9A1F2D"
CONTROL_COLOR = "#0B6670"
BHLHE40_COLOR = "#6A51A3"
SOX2_COLOR = "#8C8C8C"
LIGHT_GREY = "#D9D9D9"
MID_GREY = "#6F6F6F"
GRID_GREY = "#EEEEEE"

TF_COLORS = {
    "NFE2L2": LESION_COLOR,
    "THRB": CONTROL_COLOR,
    "BHLHE40": BHLHE40_COLOR,
    "SOX2": SOX2_COLOR,
}

ROW_BG = {
    "NFE2L2": "#F8E8EA",
    "THRB": "#E5F2F3",
    "BHLHE40": "#EFEAF7",
    "SOX2": "#F0F0F0",
}

INTERPRETATION = {
    "NFE2L2": "dominant\nperturbation",
    "THRB": "recovery\nanchor",
    "BHLHE40": "secondary",
    "SOX2": "retained",
}

SEARCH_DIRS = [
    ROOT / "celloracle_run",
    ROOT / "manuscript_output",
    ROOT / "analysis_outputs",
    ROOT / "final_exports",
    ROOT / "robustness_validation",
]

PREFERRED_METRIC_FILES = [
    ROOT / "manuscript_output" / "tables_main" / "Table3_integrated_priority_main_submission.csv",
    ROOT / "manuscript_output" / "tables_main" / "Table3_integrated_priority_main.csv",
    ROOT / "celloracle_run" / "ko_round1" / "round1_ko_master_table.csv",
    ROOT / "celloracle_run" / "ko_round1" / "round1_recovery_index_ranking.csv",
    ROOT / "manuscript_output" / "tables_supplementary" / "TableS3_celloracle_round1_full_metrics_submission.csv",
]

COL_SYNONYMS = {
    "tf": ["tf", "TF", "gene", "target_tf"],
    "overall_shift": ["mean_shift_length", "overall_mean_shift", "shift", "Mean_shift", "mean_shift"],
    "approx_ri": ["recovery_index", "Recovery Index", "Recovery_index", "approximate_recovery_index", "RI", "ri"],
    "net_shift": ["mean_net_shift", "net_shift", "Net_shift"],
    "lesion_shift": ["lesion_mean_shift_length", "lesion_shift", "lesion_mean_shift", "Lesion_shift"],
    "control_shift": [
        "internal_control_mean_shift_length",
        "internal_control_shift",
        "internal_control_mean_shift",
        "control_shift",
        "Control_shift",
    ],
}

VECTOR_COL_SYNONYMS = {
    "x": ["embedding_x", "x", "umap_x", "X_umap_1"],
    "y": ["embedding_y", "y", "umap_y", "X_umap_2"],
    "dx": ["delta_x", "shift_x", "velocity_x"],
    "dy": ["delta_y", "shift_y", "velocity_y"],
    "group": ["group", "celloracle_group", "condition"],
    "shift_length": ["shift_length", "mean_shift_length"],
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.4,
            "axes.titlesize": 8.8,
            "axes.labelsize": 7.6,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.2,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.65,
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


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t", **kwargs)
    return pd.read_csv(path, **kwargs)


def normalize_col(name: str) -> str:
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def find_col(df: pd.DataFrame, candidates: list[str], label: str, warnings: list[str], required: bool = True) -> str | None:
    by_norm = {normalize_col(col): col for col in df.columns}
    for candidate in candidates:
        key = normalize_col(candidate)
        if key in by_norm:
            return by_norm[key]
    if required:
        warnings.append(f"Missing required column for {label}; candidates={candidates}; available={list(df.columns)}")
    return None


def metric_file_score(path: Path) -> tuple[int, dict[str, str | None]]:
    try:
        df = read_table(path, nrows=5)
    except Exception:
        return (-1, {})

    cols = {
        key: find_col(df, candidates, key, [], required=False)
        for key, candidates in COL_SYNONYMS.items()
    }
    required = ["tf", "overall_shift", "approx_ri"]
    if any(cols.get(key) is None for key in required):
        return (-1, cols)

    score = 10
    score += sum(cols.get(key) is not None for key in ["net_shift", "lesion_shift", "control_shift"]) * 4
    name = path.name.lower()
    if "round1_ko_master_table" in name:
        score += 12
    if "recovery_index_ranking" in name:
        score += 8
    if "table3_integrated_priority" in name:
        score += 6
    if "tables_supplementary" in path.as_posix().lower():
        score += 2
    return (score, cols)


def find_metric_table(notes: list[str]) -> Path:
    seen: set[Path] = set()
    candidates: list[Path] = []
    for path in PREFERRED_METRIC_FILES:
        if path.exists() and path not in seen:
            candidates.append(path)
            seen.add(path)

    for root in SEARCH_DIRS:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if path.suffix.lower() not in {".csv", ".tsv"}:
                continue
            name = path.name.lower()
            if any(token in name for token in ["celloracle", "round1", "ko", "recovery", "table3", "priority", "master"]):
                if path not in seen:
                    candidates.append(path)
                    seen.add(path)

    scored: list[tuple[int, Path]] = []
    skipped_table3: list[Path] = []
    for path in candidates:
        score, cols = metric_file_score(path)
        if path.name.lower().startswith("table3_integrated_priority") and score < 0:
            skipped_table3.append(path)
        if score >= 0:
            scored.append((score, path))

    if not scored:
        raise FileNotFoundError(
            "No CellOracle metric table with TF, mean shift, and approximate Recovery Index columns was found."
        )
    scored.sort(key=lambda item: (item[0], -len(str(item[1]))), reverse=True)
    selected = scored[0][1]
    if skipped_table3:
        notes.append(
            "Checked Table3 integrated priority table(s), but they did not contain CellOracle shift/Approx. RI metric columns; "
            f"selected {rel(selected)} instead."
        )
    return selected


def load_metrics(path: Path, warnings: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    df = read_table(path)
    cols: dict[str, str] = {}
    for key, candidates in COL_SYNONYMS.items():
        col = find_col(df, candidates, key, warnings, required=key in {"tf", "overall_shift", "approx_ri"})
        if col is not None:
            cols[key] = col

    metrics = pd.DataFrame()
    metrics["tf"] = df[cols["tf"]].astype(str)
    metrics = metrics[metrics["tf"].isin(MAIN_TFS)].copy()
    missing_tfs = [tf for tf in MAIN_TFS if tf not in set(metrics["tf"])]
    if missing_tfs:
        warnings.append(f"Metric table is missing requested TF rows: {missing_tfs}")

    source = df.set_index(cols["tf"], drop=False)
    rows = []
    for tf in MAIN_TFS:
        if tf not in source.index:
            continue
        row = source.loc[tf]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        out = {"tf": tf}
        for target in ["overall_shift", "approx_ri", "net_shift", "lesion_shift", "control_shift"]:
            col = cols.get(target)
            out[target] = pd.to_numeric(row[col], errors="coerce") if col is not None else np.nan
        rows.append(out)
    metrics = pd.DataFrame(rows)

    if "net_shift" not in cols or metrics["net_shift"].isna().all():
        metrics["net_shift"] = np.nan
        warnings.append("mean_net_shift/net_shift column was not found; Panel A uses uniform point size and Panel B leaves Net shift blank.")

    if "lesion_shift" not in cols or "control_shift" not in cols:
        group_path = ROOT / "celloracle_run" / "ko_round1" / "round1_ko_summary_by_group.csv"
        if group_path.exists():
            by_group = read_table(group_path)
            tf_col = find_col(by_group, ["tf", "TF"], "group-summary TF", warnings)
            group_col = find_col(by_group, ["group", "condition"], "group-summary group", warnings)
            shift_col = find_col(
                by_group,
                ["mean_shift_length", "overall_mean_shift", "shift", "Mean_shift"],
                "group-summary mean shift",
                warnings,
            )
            pivot = by_group.pivot(index=tf_col, columns=group_col, values=shift_col)
            if "lesion_shift" not in cols and "lesion" in pivot.columns:
                metrics["lesion_shift"] = metrics["tf"].map(pivot["lesion"])
                cols["lesion_shift"] = f"{rel(group_path)}::{shift_col}[group=lesion]"
            if "control_shift" not in cols and "internal_control" in pivot.columns:
                metrics["control_shift"] = metrics["tf"].map(pivot["internal_control"])
                cols["control_shift"] = f"{rel(group_path)}::{shift_col}[group=internal_control]"
        else:
            warnings.append("Group-specific shift columns were not found and round1_ko_summary_by_group.csv is unavailable.")

    for col in ["overall_shift", "approx_ri", "net_shift", "lesion_shift", "control_shift"]:
        if col not in metrics:
            metrics[col] = np.nan
        metrics[col] = pd.to_numeric(metrics[col], errors="coerce")

    return metrics.set_index("tf").loc[[tf for tf in MAIN_TFS if tf in set(metrics["tf"])]].reset_index(), cols


def vector_path(tf: str) -> Path:
    return ROOT / "celloracle_run" / "ko_round1" / tf / f"{tf}_state_shift_scores.csv"


def load_vector_scores(tf: str, warnings: list[str]) -> tuple[pd.DataFrame | None, dict[str, str]]:
    path = vector_path(tf)
    if not path.exists():
        warnings.append(f"Vector field score file not found for {tf}: {rel(path)}")
        return None, {}
    df = read_table(path)
    cols: dict[str, str] = {}
    for key, candidates in VECTOR_COL_SYNONYMS.items():
        col = find_col(df, candidates, f"{tf} vector {key}", warnings, required=key in {"x", "y", "dx", "dy"})
        if col is not None:
            cols[key] = col
    required = {"x", "y", "dx", "dy"}
    if not required.issubset(cols):
        warnings.append(f"Vector field file for {tf} lacked required simulated shift columns; panel will show a missing-data notice.")
        return None, cols
    out = pd.DataFrame(
        {
            "x": pd.to_numeric(df[cols["x"]], errors="coerce"),
            "y": pd.to_numeric(df[cols["y"]], errors="coerce"),
            "dx": pd.to_numeric(df[cols["dx"]], errors="coerce"),
            "dy": pd.to_numeric(df[cols["dy"]], errors="coerce"),
        }
    )
    if "group" in cols:
        out["group"] = df[cols["group"]].astype(str)
    if "shift_length" in cols:
        out["shift_length"] = pd.to_numeric(df[cols["shift_length"]], errors="coerce")
    else:
        out["shift_length"] = np.sqrt(out["dx"] ** 2 + out["dy"] ** 2)
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["x", "y", "dx", "dy"])
    return out, cols


def downsample_arrows(df: pd.DataFrame, max_arrows: int = 800, seed: int = 7) -> pd.DataFrame:
    if len(df) <= max_arrows:
        return df.copy()
    rng = np.random.default_rng(seed)
    idx = rng.choice(df.index.to_numpy(), size=max_arrows, replace=False)
    return df.loc[np.sort(idx)].copy()


def style_axes(ax: plt.Axes, grid_axis: str | None = None) -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_GREY, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.115,
        1.085,
        label,
        transform=ax.transAxes,
        fontsize=12.0,
        fontweight="bold",
        ha="left",
        va="top",
        clip_on=False,
    )


def scale_point_sizes(values: pd.Series) -> np.ndarray:
    vals = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(vals).any() or np.nanmax(vals) == np.nanmin(vals):
        return np.full(len(vals), 390.0)
    low, high = np.nanmin(vals), np.nanmax(vals)
    scaled = (vals - low) / (high - low)
    return 260.0 + 430.0 * scaled


def point_size_from_range(value: float, low: float, high: float) -> float:
    if not np.isfinite(value) or not np.isfinite(low) or not np.isfinite(high) or high == low:
        return 390.0
    return float(260.0 + 430.0 * ((value - low) / (high - low)))


def plot_panel_a(ax: plt.Axes, metrics: pd.DataFrame, warnings: list[str]) -> None:
    data = metrics.set_index("tf").loc[MAIN_TFS].copy()
    sizes = scale_point_sizes(data["net_shift"])
    if data["net_shift"].isna().all():
        warnings.append("Panel A point size is uniform because mean_net_shift was unavailable.")

    ax.scatter(
        data["overall_shift"],
        data["approx_ri"],
        s=sizes,
        c=[TF_COLORS[tf] for tf in MAIN_TFS],
        edgecolor="white",
        linewidth=1.1,
        alpha=0.94,
        zorder=3,
    )

    x_ref = float(data["overall_shift"].median())
    y_ref = float(data["approx_ri"].median())
    ax.axvline(x_ref, color=GRID_GREY, linestyle=(0, (3, 3)), linewidth=0.9, zorder=1)
    ax.axhline(y_ref, color=GRID_GREY, linestyle=(0, (3, 3)), linewidth=0.9, zorder=1)

    label_offsets = {
        "NFE2L2": (8, -20),
        "THRB": (-18, 28),
        "BHLHE40": (8, 10),
        "SOX2": (8, -12),
    }
    label_ha = {
        "NFE2L2": "left",
        "THRB": "right",
        "BHLHE40": "left",
        "SOX2": "left",
    }
    label_text = {
        "NFE2L2": "NFE2L2\nhighest shift",
        "THRB": "THRB\nhighest Approx. RI",
        "BHLHE40": "BHLHE40",
        "SOX2": "SOX2",
    }
    for tf in MAIN_TFS:
        ax.annotate(
            label_text[tf],
            (float(data.loc[tf, "overall_shift"]), float(data.loc[tf, "approx_ri"])),
            xytext=label_offsets[tf],
            textcoords="offset points",
            fontsize=7.1 if tf in {"NFE2L2", "THRB"} else 6.8,
            fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal",
            color=TF_COLORS[tf],
            ha=label_ha[tf],
            va="center",
            arrowprops={
                "arrowstyle": "-",
                "color": TF_COLORS[tf],
                "linewidth": 0.7,
                "shrinkA": 2,
                "shrinkB": 3,
            },
            zorder=4,
        )

    x_pad = max((data["overall_shift"].max() - data["overall_shift"].min()) * 0.22, 0.035)
    y_pad = max((data["approx_ri"].max() - data["approx_ri"].min()) * 0.18, 0.045)
    ax.set_xlim(float(data["overall_shift"].min() - x_pad), float(data["overall_shift"].max() + x_pad))
    ax.set_ylim(float(data["approx_ri"].min() - y_pad), float(data["approx_ri"].max() + y_pad))
    ax.set_xlabel("Mean shift length")
    ax.set_ylabel("Approx. Recovery Index")
    ax.set_title("Perturbation priority map", pad=8)
    style_axes(ax, grid_axis="both")

    if data["net_shift"].notna().any():
        vals = data["net_shift"].dropna()
        legend_vals = [float(vals.min()), float(vals.max())]
        legend_low, legend_high = float(vals.min()), float(vals.max())
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                markerfacecolor="#BBBBBB",
                markeredgecolor="white",
                markersize=np.sqrt(point_size_from_range(value, legend_low, legend_high)) / 2.6,
                label=f"{value:.2f}",
            )
            for value in legend_vals
        ]
        ax.legend(handles=handles, title="Mean net shift", frameon=False, loc="lower right", borderpad=0.2)
    panel_label(ax, "A")


def column_norm(values: pd.Series) -> pd.Series:
    vals = pd.to_numeric(values, errors="coerce")
    finite = vals[np.isfinite(vals)]
    if finite.empty:
        return pd.Series(0.0, index=vals.index)
    low, high = float(finite.min()), float(finite.max())
    if high == low:
        return pd.Series(0.5, index=vals.index)
    return ((vals - low) / (high - low)).fillna(0.0)


def plot_scorecard_bar(
    ax: plt.Axes,
    x0: float,
    y: float,
    width: float,
    value: float,
    norm_value: float,
    color: str,
    alpha: float,
) -> None:
    ax.add_patch(Rectangle((x0, y - 0.075), width, 0.15, facecolor="#ECECEC", edgecolor="none", zorder=1))
    ax.add_patch(
        Rectangle(
            (x0, y - 0.075),
            width * max(0.0, min(float(norm_value), 1.0)),
            0.15,
            facecolor=color,
            edgecolor="none",
            alpha=alpha,
            zorder=2,
        )
    )
    text = "" if pd.isna(value) else f"{value:.2f}"
    ax.text(x0 + width + 0.006, y, text, fontsize=5.9, ha="left", va="center", color="#333333")


def plot_panel_b(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    data = metrics.set_index("tf").loc[MAIN_TFS].copy()
    metric_cols = [
        ("Overall\nshift", "overall_shift"),
        ("Approx.\nRI", "approx_ri"),
        ("Net\nshift", "net_shift"),
        ("Lesion\nshift", "lesion_shift"),
        ("Control\nshift", "control_shift"),
    ]
    norms = {col: column_norm(data[col]) for _, col in metric_cols}

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.65, len(MAIN_TFS) + 0.55)
    ax.axis("off")
    ax.set_title("Quantitative perturbation profile", pad=8)

    tf_x = 0.02
    col_x = [0.17, 0.31, 0.45, 0.59, 0.73]
    bar_w = 0.075
    interp_x = 0.88
    header_y = len(MAIN_TFS) + 0.15

    ax.text(tf_x, header_y, "TF", fontsize=7.0, fontweight="bold", ha="left", va="bottom")
    for (label, _), x0 in zip(metric_cols, col_x):
        ax.text(x0 + bar_w / 2, header_y, label, fontsize=6.4, fontweight="bold", ha="center", va="bottom")
    ax.text(interp_x, header_y, "Interpretation", fontsize=6.6, fontweight="bold", ha="left", va="bottom")

    for i, tf in enumerate(MAIN_TFS):
        y = len(MAIN_TFS) - 1 - i
        alpha_bg = 0.78 if tf in {"NFE2L2", "THRB"} else 0.64
        ax.add_patch(
            Rectangle(
                (0.0, y - 0.32),
                1.0,
                0.64,
                facecolor=ROW_BG[tf],
                edgecolor="white",
                linewidth=0.9,
                alpha=alpha_bg,
                zorder=0,
            )
        )
        ax.text(
            tf_x,
            y,
            tf,
            fontsize=7.2,
            fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal",
            color=TF_COLORS[tf],
            ha="left",
            va="center",
        )
        bar_alpha = 0.88 if tf in {"NFE2L2", "THRB"} else 0.66
        for (_, col), x0 in zip(metric_cols, col_x):
            plot_scorecard_bar(
                ax,
                x0,
                y,
                bar_w,
                float(data.loc[tf, col]) if pd.notna(data.loc[tf, col]) else np.nan,
                float(norms[col].loc[tf]) if tf in norms[col].index else 0.0,
                TF_COLORS[tf],
                bar_alpha,
            )
        ax.text(interp_x, y, INTERPRETATION[tf], fontsize=6.5, ha="left", va="center", color="#333333")

    panel_label(ax, "B")


def plot_vector_panel(
    ax: plt.Axes,
    tf: str,
    scores: pd.DataFrame | None,
    annotation: str,
    arrow_counts: dict[str, int],
    limits: tuple[float, float, float, float] | None,
) -> None:
    ax.set_title(f"{tf} KO simulation", pad=8, color=TF_COLORS[tf], fontweight="bold")
    if scores is None or scores.empty:
        ax.text(0.5, 0.5, "Vector field data unavailable", transform=ax.transAxes, ha="center", va="center", color=MID_GREY)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return

    arrows = downsample_arrows(scores, max_arrows=800, seed=7)
    arrow_counts[tf] = len(arrows)

    ax.scatter(scores["x"], scores["y"], s=3.2, c=LIGHT_GREY, alpha=0.42, linewidths=0, zorder=1)
    ax.quiver(
        arrows["x"],
        arrows["y"],
        arrows["dx"],
        arrows["dy"],
        angles="xy",
        scale_units="xy",
        scale=1.35,
        width=0.0022,
        headwidth=3.3,
        headlength=4.2,
        headaxislength=3.6,
        color=TF_COLORS[tf],
        alpha=0.55,
        zorder=2,
    )

    ax.text(
        0.03,
        0.96,
        annotation,
        transform=ax.transAxes,
        fontsize=7.1,
        color=TF_COLORS[tf],
        fontweight="bold",
        ha="left",
        va="top",
    )
    if limits is not None:
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
    else:
        x_pad = (scores["x"].max() - scores["x"].min()) * 0.04
        y_pad = (scores["y"].max() - scores["y"].min()) * 0.04
        ax.set_xlim(scores["x"].min() - x_pad, scores["x"].max() + x_pad)
        ax.set_ylim(scores["y"].min() - y_pad, scores["y"].max() + y_pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def vector_limits(score_tables: list[pd.DataFrame | None]) -> tuple[float, float, float, float] | None:
    valid = [df for df in score_tables if df is not None and not df.empty]
    if not valid:
        return None
    all_x = pd.concat([df["x"] for df in valid])
    all_y = pd.concat([df["y"] for df in valid])
    x_pad = (all_x.max() - all_x.min()) * 0.045
    y_pad = (all_y.max() - all_y.min()) * 0.045
    return (float(all_x.min() - x_pad), float(all_x.max() + x_pad), float(all_y.min() - y_pad), float(all_y.max() + y_pad))


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


def format_metrics_table(metrics: pd.DataFrame) -> str:
    cols = ["tf", "overall_shift", "approx_ri", "net_shift", "lesion_shift", "control_shift"]
    table = metrics[cols].copy()
    for col in cols[1:]:
        table[col] = table[col].map(lambda value: "" if pd.isna(value) else f"{float(value):.4f}")
    return table.to_csv(index=False)


def write_notes(
    path: Path,
    metric_path: Path,
    metric_cols: dict[str, str],
    metric_selection_notes: list[str],
    vector_cols: dict[str, dict[str, str]],
    arrow_counts: dict[str, int],
    metrics: pd.DataFrame,
    outputs: dict[str, Path],
    warnings: list[str],
) -> None:
    vector_paths = [vector_path(tf) for tf in ["NFE2L2", "THRB"]]
    notes = f"""
    {STEM}
    ==============================

    Caption title:
    {CAPTION_TITLE}

    Revision scope:
    - Main-text Fig. 3 v3 was rebuilt as a CellOracle / in silico KO perturbation prioritization figure.
    - Robustness, supportive external/contextual evidence, exploratory drug-signature clues, and pySCENIC regulon landscape panels were not included.
    - No complete figure suptitle is drawn inside the graphic; panel titles only.

    Input files:
    - Metric table: {rel(metric_path)}
    - NFE2L2 vector field scores: {rel(vector_paths[0])}
    - THRB vector field scores: {rel(vector_paths[1])}

    Metric columns used:
    - TF: {metric_cols.get("tf", "missing")}
    - Overall shift / mean shift length: {metric_cols.get("overall_shift", "missing")}
    - Approx. Recovery Index: {metric_cols.get("approx_ri", "missing")}
    - Mean net shift: {metric_cols.get("net_shift", "missing")}
    - Lesion shift: {metric_cols.get("lesion_shift", "missing")}
    - Control shift: {metric_cols.get("control_shift", "missing")}

    Metric table values used:
    {textwrap.indent(format_metrics_table(metrics).strip(), "    ")}

    Panel A encoding:
    - x-axis = Mean shift length.
    - y-axis = Approx. Recovery Index.
    - point size = Mean net shift when available; otherwise uniform size.
    - NFE2L2 is highlighted as the highest mean shift, and THRB as the highest Approx. RI.

    Panel B scorecard:
    - Rows: NFE2L2, THRB, BHLHE40, SOX2.
    - Metrics: Overall shift, Approx. RI, Net shift, Lesion shift, Control shift.
    - Each metric column uses min-max normalization within the four TFs for horizontal mini-bar length; numeric values are printed beside bars.
    - Interpretation labels: NFE2L2 = dominant perturbation; THRB = recovery anchor; BHLHE40 = secondary; SOX2 = retained.

    Panel C/D vector field sources:
    - Simulated KO shift columns are embedding_x/embedding_y and delta_x/delta_y from each TF state_shift_scores.csv.
    - Randomized-control columns were not plotted in the main figure.
    - Vector column mapping: {vector_cols}

    Arrow downsampling:
    - Background points include all available cells in each TF vector table.
    - Arrows are downsampled with a fixed random seed (7) to at most 800 arrows per panel.
    - Final arrow counts: {arrow_counts}

    Required interpretation notes:
    Figure 3 v3 uses existing CellOracle outputs only; no CellOracle rerun was performed.
    Approx. Recovery Index is treated as an approximate state-shift-derived metric, not a formal geometric projection model.

    Metric-selection notes:
    {textwrap.indent(chr(10).join(metric_selection_notes) if metric_selection_notes else "None.", "    ")}

    Warnings / missing data:
    {textwrap.indent(chr(10).join(warnings) if warnings else "None.", "    ")}

    Output files:
    - {rel(outputs["png"])}
    - {rel(outputs["pdf"])}
    - {rel(outputs["svg"])}
    - {rel(outputs["notes"])}
    """
    path.write_text(textwrap.dedent(notes).strip() + "\n", encoding="utf-8")


def main() -> None:
    setup_style()
    warnings: list[str] = []
    metric_selection_notes: list[str] = []
    metric_path = find_metric_table(metric_selection_notes)
    metrics, metric_cols = load_metrics(metric_path, warnings)

    nfe_scores, nfe_cols = load_vector_scores("NFE2L2", warnings)
    thrb_scores, thrb_cols = load_vector_scores("THRB", warnings)
    limits = vector_limits([nfe_scores, thrb_scores])
    vector_cols = {"NFE2L2": nfe_cols, "THRB": thrb_cols}
    arrow_counts: dict[str, int] = {}

    fig = plt.figure(figsize=(14.0, 9.1), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.03, 1.0],
        width_ratios=[1.0, 1.05],
        hspace=0.34,
        wspace=0.24,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    plot_panel_a(ax_a, metrics, warnings)
    plot_panel_b(ax_b, metrics)
    plot_vector_panel(ax_c, "NFE2L2", nfe_scores, "Highest mean shift", arrow_counts, limits)
    panel_label(ax_c, "C")
    plot_vector_panel(ax_d, "THRB", thrb_scores, "Highest Approx. RI", arrow_counts, limits)
    panel_label(ax_d, "D")

    outputs = save_outputs(fig)
    write_notes(outputs["notes"], metric_path, metric_cols, metric_selection_notes, vector_cols, arrow_counts, metrics, outputs, warnings)

    for key in ["png", "pdf", "svg", "notes"]:
        print(f"Wrote {rel(outputs[key])}")
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"- {warning}")
    else:
        print("Warnings: none")


if __name__ == "__main__":
    main()
