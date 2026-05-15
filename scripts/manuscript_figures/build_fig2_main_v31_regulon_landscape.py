from __future__ import annotations

import math
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_main"
STEM = "Fig2_main_v31_regulon_landscape"
TITLE = (
    "Fig. 2. Differential regulon landscape identifies opposing lesion- and "
    "internal-control-associated TF programs"
)

GROUP_ORDER = ["internal_control", "lesion"]
CONTROL_COLOR = "#0B6670"
LESION_COLOR = "#9A1F2D"
BHLHE40_COLOR = "#6A51A3"
SOX2_COLOR = "#8C8C8C"
LIGHT_GREY = "#D8D8D8"
MID_GREY = "#777777"
GRID_GREY = "#EFEFEF"

TF_COLORS = {
    "NFE2L2": LESION_COLOR,
    "THRB": CONTROL_COLOR,
    "BHLHE40": BHLHE40_COLOR,
    "SOX2": SOX2_COLOR,
    "SATB2": CONTROL_COLOR,
    "RARB": CONTROL_COLOR,
    "HMGA1": LESION_COLOR,
}

PANEL_A_CONTROL_PREF = [
    "THRB(+)",
    "SATB2(+)",
    "RARB(+)",
    "TCF4(+)",
    "NR2F2(+)",
    "SMAD3(+)",
    "SKI(+)",
    "HSF2(+)",
]
PANEL_A_LESION_PREF = [
    "NFE2L2(+)",
    "BHLHE40(+)",
    "SOX2(+)",
    "HMGA1(+)",
    "CEBPD(+)",
    "IRF9(+)",
    "ATF7(+)",
    "GABPB1(+)",
    "JUNB(+)",
    "DBX2(+)",
]
PANEL_B_LABELS = [
    "THRB(+)",
    "SATB2(+)",
    "RARB(+)",
    "NFE2L2(+)",
    "BHLHE40(+)",
    "SOX2(+)",
]
SCORECARD_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2", "SATB2", "RARB", "HMGA1"]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.1,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.4,
            "legend.fontsize": 6.0,
            "figure.titlesize": 10.2,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
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


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def find_input(preferred: list[str], tokens: list[str], warnings: list[str]) -> Path:
    for item in preferred:
        path = ROOT / item
        if path.exists():
            return path
    roots = [
        ROOT / "manuscript_output",
        ROOT / "analysis_outputs",
        ROOT / "final_exports",
        ROOT / "output",
    ]
    candidates: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for suffix in ("*.csv", "*.tsv"):
            candidates.extend(root.rglob(suffix))
    for path in sorted(candidates):
        name = path.name.lower()
        if all(token.lower() in name for token in tokens):
            warnings.append(f"Input fallback used for {tokens}: {rel(path)}")
            return path
    raise FileNotFoundError(f"Could not find input matching tokens: {tokens}")


def find_col(
    df: pd.DataFrame,
    aliases: list[str],
    label: str,
    warnings: list[str],
    required: bool = True,
) -> str | None:
    norm_to_col = {normalize_name(col): col for col in df.columns}
    for alias in aliases:
        key = normalize_name(alias)
        if key in norm_to_col:
            return norm_to_col[key]
    if required:
        warnings.append(f"Missing required column for {label}; aliases tried: {aliases}")
    return None


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def neg_log10(values: pd.Series | np.ndarray) -> np.ndarray:
    x = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    finite_positive = x[np.isfinite(x) & (x > 0)]
    floor = max(float(np.nanmin(finite_positive)) * 0.1, 1e-300) if finite_positive.size else 1e-300
    x = np.where(np.isfinite(x) & (x > 0), x, floor)
    return -np.log10(x)


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.sub(df.mean(axis=1), axis=0)
    out = out.div(df.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    return out.fillna(0.0)


def load_and_standardize_stats(path: Path, warnings: list[str]) -> pd.DataFrame:
    raw = read_table(path)
    reg_col = find_col(raw, ["regulon", "Regulon"], "regulon", warnings)
    fdr_col = find_col(raw, ["fdr_bh", "FDR", "fdr", "padj", "qval", "Regulon_FDR"], "regulon FDR", warnings)
    mean_lesion_col = find_col(
        raw,
        ["mean_lesion", "Mean_AUC_lesion", "lesion_mean_auc", "mean_auc_lesion"],
        "mean AUC lesion",
        warnings,
        required=False,
    )
    mean_control_col = find_col(
        raw,
        [
            "mean_internal_control",
            "Mean_AUC_control",
            "mean_control",
            "control_mean_auc",
            "mean_auc_internal_control",
        ],
        "mean AUC internal_control",
        warnings,
        required=False,
    )
    diff_col = find_col(
        raw,
        [
            "regulon_mean_diff_lesion_minus_internal_control",
            "dAUC_lesion_minus_control",
            "mean_auc_diff",
            "auc_diff",
            "mean_difference",
            "avg_diff",
            "mean_diff_group1_minus_group2",
            "Regulon_dAUC",
        ],
        "mean AUC difference",
        warnings,
        required=False,
    )

    if reg_col is None or fdr_col is None:
        raise ValueError(f"Cannot standardize regulon statistics from {path}")

    stats = pd.DataFrame()
    stats["regulon"] = raw[reg_col].astype(str)
    stats["fdr"] = numeric(raw[fdr_col])

    if mean_lesion_col and mean_control_col:
        stats["mean_lesion"] = numeric(raw[mean_lesion_col])
        stats["mean_internal_control"] = numeric(raw[mean_control_col])
        stats["dauc_lesion_minus_control"] = stats["mean_lesion"] - stats["mean_internal_control"]
    elif diff_col:
        stats["dauc_lesion_minus_control"] = numeric(raw[diff_col])
        stats["mean_lesion"] = np.nan
        stats["mean_internal_control"] = np.nan
        warnings.append(f"Mean AUC columns unavailable in {rel(path)}; using {diff_col} for dAUC.")
    else:
        raise ValueError(f"Cannot identify mean AUC difference columns in {path}")

    stats["abs_dauc"] = stats["dauc_lesion_minus_control"].abs()
    stats["neg_log10_fdr"] = neg_log10(stats["fdr"])
    stats["axis"] = np.where(stats["dauc_lesion_minus_control"] >= 0, "lesion", "internal_control")
    return stats


def load_and_standardize_shortlist(path: Path, stats: pd.DataFrame, warnings: list[str]) -> pd.DataFrame:
    raw = read_table(path)
    tf_col = find_col(raw, ["tf", "TF"], "TF", warnings)
    regulon_col = find_col(raw, ["regulon", "Regulon"], "candidate regulon", warnings, required=False)
    dauc_col = find_col(
        raw,
        [
            "regulon_mean_diff_lesion_minus_internal_control",
            "Regulon_dAUC",
            "dAUC_lesion_minus_control",
            "mean_auc_diff",
            "auc_diff",
        ],
        "candidate regulon dAUC",
        warnings,
        required=False,
    )
    fdr_col = find_col(
        raw,
        ["regulon_fdr", "Regulon_FDR", "FDR", "fdr_bh", "fdr", "padj", "qval"],
        "candidate regulon FDR",
        warnings,
        required=False,
    )
    expr_col = find_col(
        raw,
        [
            "log2fc_lesion_vs_internal_control",
            "Expr_log2FC",
            "log2FC",
            "tf_log2fc",
            "expression_log2fc",
        ],
        "TF expression log2FC",
        warnings,
        required=False,
    )
    pos_col = find_col(
        raw,
        [
            "positive_fraction_in_higher_group",
            "Pos_frac_high_group",
            "positive_fraction",
            "pct_positive",
            "pos_frac",
        ],
        "positive fraction",
        warnings,
        required=False,
    )

    if tf_col is None:
        raise ValueError(f"Cannot standardize TF shortlist from {path}")

    out = pd.DataFrame()
    out["tf"] = raw[tf_col].astype(str)
    if regulon_col:
        out["regulon"] = raw[regulon_col].astype(str)
    else:
        out["regulon"] = out["tf"] + "(+)"
        warnings.append(f"Regulon column missing in {rel(path)}; inferred regulon as TF(+).")

    stats_lookup = stats.set_index("regulon")
    if dauc_col:
        out["dauc_lesion_minus_control"] = numeric(raw[dauc_col])
    else:
        out["dauc_lesion_minus_control"] = out["regulon"].map(stats_lookup["dauc_lesion_minus_control"])
        warnings.append(f"Candidate dAUC column missing in {rel(path)}; filled from regulon statistics.")
    if fdr_col:
        out["regulon_fdr"] = numeric(raw[fdr_col])
    else:
        out["regulon_fdr"] = out["regulon"].map(stats_lookup["fdr"])
        warnings.append(f"Candidate regulon FDR column missing in {rel(path)}; filled from regulon statistics.")
    if expr_col:
        out["expr_log2fc"] = numeric(raw[expr_col])
    else:
        out["expr_log2fc"] = np.nan
        warnings.append(f"TF expression log2FC unavailable in {rel(path)}.")
    if pos_col:
        out["positive_fraction"] = numeric(raw[pos_col])
    else:
        out["positive_fraction"] = np.nan
        warnings.append(f"Positive fraction unavailable in {rel(path)}.")

    out["neg_log10_regulon_fdr"] = neg_log10(out["regulon_fdr"])
    return out


def load_auc_with_metadata(auc_path: Path, meta_path: Path, warnings: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    auc = read_table(auc_path)
    auc = auc.rename(columns={auc.columns[0]: "cell_id"})
    meta = read_table(meta_path)
    if meta.columns[0].startswith("Unnamed") or meta.columns[0] in {"", "H1"}:
        meta = meta.rename(columns={meta.columns[0]: "cell_id"})
    if "cell_id" not in meta.columns and "join_key" in meta.columns:
        meta["cell_id"] = meta["join_key"]
        warnings.append("Metadata cell_id field inferred from join_key.")

    group_col = find_col(meta, ["group", "condition", "disease_group", "celloracle_group"], "obs group", warnings)
    sample_col = find_col(
        meta,
        ["clinical_sample_label", "sample_id", "sample", "orig.ident", "sample_prefix"],
        "obs sample",
        warnings,
    )
    donor_col = find_col(meta, ["donor_id", "donor", "patient_id"], "obs donor", warnings, required=False)

    if group_col is None or sample_col is None or "cell_id" not in meta.columns:
        raise ValueError(f"Cannot load AUCell metadata from {meta_path}")

    meta_keep = meta[["cell_id", group_col, sample_col] + ([donor_col] if donor_col else [])].copy()
    meta_keep = meta_keep.rename(columns={group_col: "group", sample_col: "sample"})
    if donor_col:
        meta_keep = meta_keep.rename(columns={donor_col: "donor"})
    else:
        meta_keep["donor"] = meta_keep["sample"].astype(str).str.extract(r"(G\d+)", expand=False).fillna("unknown")
        donor_col = "inferred_from_sample"
        warnings.append("Donor field unavailable; donor inferred from sample label.")

    merged = auc.merge(meta_keep, on="cell_id", how="left")
    if merged["group"].isna().any():
        warnings.append("Some AUCell rows lacked matched group metadata after merge.")

    fields = {"group": group_col, "sample": sample_col, "donor": donor_col}
    return merged, fields


def sample_order_table(auc_meta: pd.DataFrame) -> pd.DataFrame:
    sample_meta = auc_meta[["sample", "group", "donor"]].drop_duplicates()
    group_rank = {group: idx for idx, group in enumerate(GROUP_ORDER)}
    sample_meta["group_rank"] = sample_meta["group"].map(group_rank).fillna(99)
    sample_meta = sample_meta.sort_values(["group_rank", "donor", "sample"]).drop(columns="group_rank")
    return sample_meta.reset_index(drop=True)


def unique_keep_order(items: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


def choose_panel_a_regulons(stats: pd.DataFrame, auc_cols: set[str], warnings: list[str]) -> tuple[list[str], list[str], list[str]]:
    stats_set = set(stats["regulon"])

    def present(items: list[str]) -> list[str]:
        return [item for item in items if item in stats_set and item in auc_cols]

    control = present(PANEL_A_CONTROL_PREF)
    lesion = present(PANEL_A_LESION_PREF)
    missing_pref = [item for item in PANEL_A_CONTROL_PREF + PANEL_A_LESION_PREF if item not in stats_set or item not in auc_cols]

    def supplement(axis: str, selected: list[str], target: int) -> list[str]:
        if len(selected) >= target:
            return selected
        sub = stats[(stats["axis"].eq(axis)) & (stats["regulon"].isin(auc_cols))].copy()
        sub = sub.sort_values(["fdr", "abs_dauc"], ascending=[True, False])
        for reg in sub["regulon"].astype(str):
            if reg not in selected and reg not in control and reg not in lesion:
                selected.append(reg)
            if len(selected) >= target:
                break
        return selected

    control = supplement("internal_control", control, len(PANEL_A_CONTROL_PREF))
    lesion = supplement("lesion", lesion, len(PANEL_A_LESION_PREF))
    selected = unique_keep_order(control + lesion)

    if len(selected) < 18:
        sub = stats[stats["regulon"].isin(auc_cols)].sort_values(["fdr", "abs_dauc"], ascending=[True, False])
        for reg in sub["regulon"].astype(str):
            if reg not in selected:
                selected.append(reg)
            if len(selected) >= 18:
                break

    if len(selected) > 24:
        selected = selected[:24]
        control = [reg for reg in control if reg in selected]
        lesion = [reg for reg in lesion if reg in selected]

    if missing_pref:
        warnings.append("Preferred Panel A regulons missing or unavailable in AUCell: " + ", ".join(missing_pref))

    return selected, [reg for reg in control if reg in selected], [reg for reg in lesion if reg in selected]


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=11.0,
        fontweight="bold",
        ha="left",
        va="top",
    )


def style_axis(ax: plt.Axes, grid_axis: str | None = None) -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_GREY, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def add_sample_annotation_strip(ax: plt.Axes, sample_meta: pd.DataFrame, donor_colors: dict[str, str]) -> None:
    strip = ax.inset_axes([0.0, 1.015, 1.0, 0.105])
    strip.set_xlim(-0.5, len(sample_meta) - 0.5)
    strip.set_ylim(0, 2)
    for i, row in sample_meta.iterrows():
        group_color = CONTROL_COLOR if row["group"] == "internal_control" else LESION_COLOR
        strip.add_patch(Rectangle((i - 0.5, 1), 1.0, 1.0, facecolor=group_color, edgecolor="white", linewidth=0.5))
        strip.add_patch(
            Rectangle(
                (i - 0.5, 0),
                1.0,
                1.0,
                facecolor=donor_colors.get(str(row["donor"]), "#BBBBBB"),
                edgecolor="white",
                linewidth=0.5,
            )
        )
    strip.text(-0.62, 1.5, "group", fontsize=5.6, ha="right", va="center", clip_on=False)
    strip.text(-0.62, 0.5, "donor", fontsize=5.6, ha="right", va="center", clip_on=False)
    strip.set_xticks([])
    strip.set_yticks([])
    for spine in strip.spines.values():
        spine.set_visible(False)


def plot_panel_a(
    fig: plt.Figure,
    ax: plt.Axes,
    heat: pd.DataFrame,
    sample_meta: pd.DataFrame,
    control_regs: list[str],
    lesion_regs: list[str],
    cmap: LinearSegmentedColormap,
) -> None:
    im = ax.imshow(heat.values, aspect="auto", cmap=cmap, vmin=-1.5, vmax=1.5)
    ax.set_xticks(np.arange(len(sample_meta)))
    ax.set_xticklabels(sample_meta["sample"], rotation=32, ha="right")
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    n_control = len(control_regs)
    if 0 < n_control < heat.shape[0]:
        ax.axhline(n_control - 0.5, color="#FFFFFF", linewidth=2.3)
        ax.axhline(n_control - 0.5, color="#9A9A9A", linewidth=0.45)
    ax.text(
        -0.98,
        (n_control - 1) / 2,
        "internal-\ncontrol",
        color=CONTROL_COLOR,
        fontsize=5.7,
        ha="right",
        va="center",
        clip_on=False,
    )
    ax.text(
        -0.98,
        n_control + (len(lesion_regs) - 1) / 2,
        "lesion",
        color=LESION_COLOR,
        fontsize=5.7,
        ha="right",
        va="center",
        clip_on=False,
    )

    for label in ax.get_yticklabels():
        text = label.get_text()
        if text in {"THRB(+)", "NFE2L2(+)"}:
            label.set_fontweight("bold")
            label.set_color(TF_COLORS.get(text.replace("(+)", ""), "#222222"))
        elif text in {"BHLHE40(+)", "SOX2(+)"}:
            label.set_fontweight("semibold")
            label.set_color(TF_COLORS.get(text.replace("(+)", ""), "#555555"))

    donors = sample_meta["donor"].astype(str).unique().tolist()
    donor_palette = ["#B8B8B8", "#E1E1E1", "#C9C9C9", "#F0F0F0"]
    donor_colors = {donor: donor_palette[i % len(donor_palette)] for i, donor in enumerate(donors)}
    add_sample_annotation_strip(ax, sample_meta, donor_colors)

    cbar = fig.colorbar(im, ax=ax, fraction=0.034, pad=0.024)
    cbar.set_label("row z-score of regulon AUC")
    cbar.ax.tick_params(labelsize=5.8, width=0.5, length=2.2)
    ax.legend(
        handles=[
            Line2D([0], [0], color=CONTROL_COLOR, linewidth=5, label="internal-control"),
            Line2D([0], [0], color=LESION_COLOR, linewidth=5, label="lesion"),
        ],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.38, 1.17),
        ncol=2,
        handlelength=1.4,
        columnspacing=0.9,
        borderaxespad=0.0,
    )
    ax.set_title("Sample-level regulon activity", pad=30)
    add_panel_label(ax, "A", x=-0.17, y=1.19)


def plot_panel_b(ax: plt.Axes, stats: pd.DataFrame) -> None:
    sig = stats["fdr"] <= 0.05
    point_colors = np.where(
        sig & (stats["dauc_lesion_minus_control"] < 0),
        CONTROL_COLOR,
        np.where(sig & (stats["dauc_lesion_minus_control"] > 0), LESION_COLOR, LIGHT_GREY),
    )
    point_alpha = np.where(sig, 0.74, 0.46)
    ax.scatter(
        stats["dauc_lesion_minus_control"],
        stats["neg_log10_fdr"],
        s=15,
        c=point_colors,
        alpha=point_alpha,
        edgecolors="none",
        zorder=2,
    )
    ax.axvline(0, color="#888888", linewidth=0.75)
    fdr_line = -math.log10(0.05)
    ax.axhline(fdr_line, color="#A8A8A8", linewidth=0.75, linestyle=(0, (3, 2)))
    ax.set_xlabel("Mean AUC difference: lesion - internal_control")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("All-regulon differential landscape", pad=6)
    style_axis(ax, grid_axis=None)

    ax.text(0.04, 0.87, "internal-control-associated", transform=ax.transAxes, color=CONTROL_COLOR, fontsize=6.0, ha="left")
    ax.text(0.96, 0.87, "lesion-associated", transform=ax.transAxes, color=LESION_COLOR, fontsize=6.0, ha="right")

    offsets = {
        "THRB(+)": (-8, 8, "right"),
        "SATB2(+)": (-8, -6, "right"),
        "RARB(+)": (-8, -12, "right"),
        "SOX2(+)": (8, 7, "left"),
        "NFE2L2(+)": (8, -9, "left"),
        "BHLHE40(+)": (8, 0, "left"),
    }
    stats_i = stats.set_index("regulon")
    for reg in PANEL_B_LABELS:
        if reg not in stats_i.index:
            continue
        row = stats_i.loc[reg]
        tf = reg.replace("(+)", "")
        x = float(row["dauc_lesion_minus_control"])
        y = float(row["neg_log10_fdr"])
        is_core = reg in {"THRB(+)", "NFE2L2(+)", "BHLHE40(+)", "SOX2(+)"}
        ax.scatter([x], [y], s=46 if is_core else 32, color="white", edgecolor=TF_COLORS.get(tf, "#333333"), linewidth=1.0, zorder=4)
        dx, dy, ha = offsets.get(reg, (5, 4, "left"))
        ax.annotate(
            reg,
            (x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=6.2 if is_core else 5.5,
            fontweight="bold" if reg in {"THRB(+)", "NFE2L2(+)"} else "normal",
            ha=ha,
            va="center",
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.86},
            zorder=5,
        )
    add_panel_label(ax, "B")


def plot_axis_slope(
    ax: plt.Axes,
    axis_table: pd.DataFrame,
    regulon: str,
    color: str,
    title: str,
    note: str,
    ylabel: bool,
) -> None:
    sub = axis_table[axis_table["regulon"].eq(regulon)].copy()
    donor_order = sorted(sub["donor"].astype(str).unique().tolist())
    marker_map = {"G120": "o", "G133": "s"}
    for donor in donor_order:
        dsub = sub[sub["donor"].astype(str).eq(donor)].set_index("group")
        if not all(group in dsub.index for group in GROUP_ORDER):
            continue
        xs = [0, 1]
        ys = [float(dsub.loc[group, "mean_auc"]) for group in GROUP_ORDER]
        ax.plot(xs, ys, color="#B0B0B0", linewidth=1.1, zorder=1)
        ax.scatter(
            xs,
            ys,
            s=34,
            c=[CONTROL_COLOR, LESION_COLOR],
            marker=marker_map.get(donor, "o"),
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
    ax.set_xlim(-0.18, 1.18)
    values = sub["mean_auc"].to_numpy(dtype=float)
    pad = max((np.nanmax(values) - np.nanmin(values)) * 0.22, 0.006)
    ax.set_ylim(np.nanmin(values) - pad, np.nanmax(values) + pad)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Control", "Lesion"])
    if ylabel:
        ax.set_ylabel("Sample-level mean regulon AUC")
    else:
        ax.set_yticklabels([])
    ax.set_title(title, color=color, pad=5, fontweight="bold" if regulon in {"NFE2L2(+)", "THRB(+)"} else "normal")
    ax.text(0.04, 0.93, note, transform=ax.transAxes, fontsize=6.0, color=color, ha="left", va="top")
    style_axis(ax, grid_axis="y")


def plot_panel_c(fig: plt.Figure, spec, axis_table: pd.DataFrame) -> None:
    sub = spec.subgridspec(1, 2, wspace=0.28)
    ax_c1 = fig.add_subplot(sub[0, 0])
    ax_c2 = fig.add_subplot(sub[0, 1])
    plot_axis_slope(ax_c1, axis_table, "NFE2L2(+)", LESION_COLOR, "NFE2L2(+)", "lesion-shifted", True)
    plot_axis_slope(ax_c2, axis_table, "THRB(+)", CONTROL_COLOR, "THRB(+)", "control-shifted", False)
    ax_c1.text(0.0, 1.25, "Donor-paired axis evidence", transform=ax_c1.transAxes, fontsize=8.0, ha="left", va="bottom")
    donor_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#FFFFFF", markeredgecolor=MID_GREY, markersize=4.5, label="G120"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#FFFFFF", markeredgecolor=MID_GREY, markersize=4.5, label="G133"),
    ]
    ax_c2.legend(
        handles=donor_handles,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(1.02, 1.27),
        ncol=2,
        handlelength=0.7,
        columnspacing=0.7,
        borderaxespad=0.0,
        title="donor",
        title_fontsize=5.8,
    )
    add_panel_label(ax_c1, "C", x=-0.27, y=1.32)


def scorecard_value_color(value: float, norm: TwoSlopeNorm, cmap: LinearSegmentedColormap) -> tuple:
    if not np.isfinite(value):
        return (0.92, 0.92, 0.92, 1.0)
    return cmap(norm(value))


def plot_panel_d(ax: plt.Axes, scorecard: pd.DataFrame, cmap: LinearSegmentedColormap) -> None:
    ax.set_axis_off()
    rows = scorecard.copy().set_index("tf").loc[SCORECARD_TFS].reset_index()
    n = len(rows)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.55, n + 0.98)

    x_tf = 0.06
    x_axis = 0.22
    x_reg = 0.47
    x_expr = 0.71
    x_pri = 0.91

    headers = [
        (x_tf, "TF"),
        (x_axis, "Axis"),
        (x_reg, "Regulon evidence"),
        (x_expr, "TF expression"),
        (x_pri, "Priority"),
    ]
    for x, label in headers:
        ax.text(x, n + 0.24, label, fontsize=7.0, fontweight="bold", ha="center", va="bottom")

    max_abs_dauc = max(0.01, float(np.nanmax(np.abs(rows["dauc_lesion_minus_control"]))))
    max_abs_expr = max(0.01, float(np.nanmax(np.abs(rows["expr_log2fc"]))))
    max_fdr = max(1.0, float(np.nanmax(rows["neg_log10_regulon_fdr"])))
    reg_half_width = 0.115
    expr_half_width = 0.115

    for i, row in rows.iterrows():
        y = n - 1 - i
        tf = str(row["tf"])
        alpha = 1.0
        if tf == "SOX2":
            alpha = 0.58
        if tf in {"NFE2L2", "THRB"}:
            face = "#F7ECEE" if tf == "NFE2L2" else "#EAF3F4"
            ax.add_patch(Rectangle((0.01, y - 0.38), 0.98, 0.76, facecolor=face, edgecolor="none", zorder=0))
        elif tf == "BHLHE40":
            ax.add_patch(Rectangle((0.01, y - 0.38), 0.98, 0.76, facecolor="#F2EEF7", edgecolor="none", zorder=0))
        elif tf == "SOX2":
            ax.add_patch(Rectangle((0.01, y - 0.38), 0.98, 0.76, facecolor="#F5F5F5", edgecolor="none", zorder=0))

        ax.hlines(y - 0.47, 0.01, 0.99, color="#EEEEEE", linewidth=0.55, zorder=0)
        ax.text(
            x_tf,
            y,
            tf + "(+)",
            fontsize=7.1,
            fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal",
            color=TF_COLORS.get(tf, "#333333"),
            ha="center",
            va="center",
            alpha=alpha,
        )
        axis_text = str(row["axis_display"]).replace("-\n", "-").replace("\n", " ")
        ax.text(x_axis, y, axis_text, fontsize=6.7, ha="center", va="center", color=TF_COLORS.get(tf, MID_GREY), alpha=alpha)

        dauc = float(row["dauc_lesion_minus_control"])
        reg_color = LESION_COLOR if dauc >= 0 else CONTROL_COLOR
        ax.hlines(y, x_reg - reg_half_width, x_reg + reg_half_width, color="#DCDCDC", linewidth=2.0, zorder=1)
        ax.vlines(x_reg, y - 0.18, y + 0.18, color="#B8B8B8", linewidth=0.6, zorder=2)
        reg_end = x_reg + reg_half_width * np.clip(dauc / max_abs_dauc, -1, 1)
        ax.hlines(y, min(x_reg, reg_end), max(x_reg, reg_end), color=reg_color, linewidth=4.2, alpha=0.82 * alpha, zorder=3)
        neg_fdr = float(row["neg_log10_regulon_fdr"])
        dot_size = 22 + 78 * math.sqrt(max(0.0, neg_fdr) / max_fdr)
        ax.scatter([reg_end], [y], s=dot_size, color=reg_color, edgecolor="white", linewidth=0.6, alpha=0.90 * alpha, zorder=4)

        expr = float(row["expr_log2fc"])
        expr_color = LESION_COLOR if expr >= 0 else CONTROL_COLOR
        ax.hlines(y, x_expr - expr_half_width, x_expr + expr_half_width, color="#DCDCDC", linewidth=2.0, zorder=1)
        ax.vlines(x_expr, y - 0.18, y + 0.18, color="#B8B8B8", linewidth=0.6, zorder=2)
        expr_end = x_expr + expr_half_width * np.clip(expr / max_abs_expr, -1, 1)
        ax.hlines(y, min(x_expr, expr_end), max(x_expr, expr_end), color=expr_color, linewidth=4.2, alpha=0.82 * alpha, zorder=3)

        ax.text(x_pri, y, str(row["priority_class"]), fontsize=6.7, ha="center", va="center", color="#222222", alpha=alpha)

    ax.set_title("TF prioritization scorecard", loc="left", pad=7)
    add_panel_label(ax, "D", x=-0.055, y=1.045)


def build_scorecard(shortlist: pd.DataFrame, stats: pd.DataFrame, warnings: list[str]) -> pd.DataFrame:
    stat_lookup = stats.set_index("regulon")
    rows = []
    shortlist_i = shortlist.drop_duplicates("tf").set_index("tf")
    axis_display = {
        "THRB": "internal-control",
        "SATB2": "internal-control",
        "RARB": "internal-control",
        "NFE2L2": "lesion",
        "HMGA1": "lesion",
        "BHLHE40": "lesion / secondary",
        "SOX2": "lesion / retained",
    }
    priority = {
        "NFE2L2": "dual-axis",
        "THRB": "dual-axis",
        "BHLHE40": "secondary",
        "SOX2": "retained",
        "SATB2": "supportive",
        "RARB": "supportive",
        "HMGA1": "context",
    }
    for tf in SCORECARD_TFS:
        regulon = f"{tf}(+)"
        if tf in shortlist_i.index:
            row = shortlist_i.loc[tf]
            dauc = row.get("dauc_lesion_minus_control", np.nan)
            fdr = row.get("regulon_fdr", np.nan)
            expr = row.get("expr_log2fc", np.nan)
            pos = row.get("positive_fraction", np.nan)
            source_regulon = str(row.get("regulon", regulon))
        else:
            source_regulon = regulon
            dauc = np.nan
            fdr = np.nan
            expr = np.nan
            pos = np.nan
            warnings.append(f"{tf} absent from shortlist table; using regulon statistics where possible.")

        if (not np.isfinite(pd.to_numeric(pd.Series([dauc]), errors="coerce").iloc[0])) and regulon in stat_lookup.index:
            dauc = stat_lookup.loc[regulon, "dauc_lesion_minus_control"]
        if (not np.isfinite(pd.to_numeric(pd.Series([fdr]), errors="coerce").iloc[0])) and regulon in stat_lookup.index:
            fdr = stat_lookup.loc[regulon, "fdr"]

        rows.append(
            {
                "tf": tf,
                "regulon": source_regulon,
                "axis_display": axis_display[tf],
                "dauc_lesion_minus_control": pd.to_numeric(pd.Series([dauc]), errors="coerce").iloc[0],
                "regulon_fdr": pd.to_numeric(pd.Series([fdr]), errors="coerce").iloc[0],
                "neg_log10_regulon_fdr": neg_log10(pd.Series([fdr]))[0],
                "expr_log2fc": pd.to_numeric(pd.Series([expr]), errors="coerce").iloc[0],
                "positive_fraction": pd.to_numeric(pd.Series([pos]), errors="coerce").iloc[0],
                "priority_class": priority[tf],
            }
        )
    scorecard = pd.DataFrame(rows)
    return scorecard


def format_table_for_notes(df: pd.DataFrame, max_rows: int | None = None) -> str:
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.4g}")
    if max_rows is not None:
        out = out.head(max_rows)
    return out.to_string(index=False)


def write_notes(
    notes_path: Path,
    inputs: dict[str, Path],
    obs_fields: dict[str, str],
    selected_regs: list[str],
    control_regs: list[str],
    lesion_regs: list[str],
    panel_c_table: pd.DataFrame,
    scorecard: pd.DataFrame,
    stats: pd.DataFrame,
    warnings: list[str],
) -> None:
    panel_d_notes = scorecard[
        [
            "tf",
            "regulon",
            "axis_display",
            "dauc_lesion_minus_control",
            "regulon_fdr",
            "neg_log10_regulon_fdr",
            "expr_log2fc",
            "positive_fraction",
            "priority_class",
        ]
    ].copy()
    panel_d_notes["axis_display"] = (
        panel_d_notes["axis_display"].str.replace("-\n", "-", regex=False).str.replace("\n", " ", regex=False)
    )

    notes = f"""
    {TITLE}
    ================================================================

    Script:
    - scripts/manuscript_figures/build_fig2_main_v31_regulon_landscape.py

    Revision note:
    - v3.1 is a layout/readability revision only and does not change data inputs, axis definitions, or scientific interpretation.
    - The complete in-figure suptitle was removed to avoid overlap with panel titles; the full title should be used in the caption.
    - Panel C removes donor text labels next to points and uses marker-shape legend instead.
    - Panel D removes the positive-fraction column from the main figure and uses a compact scorecard; positive fraction remains recorded below.
    - Panel C uses per-regulon y-axis ranges to preserve directionality in the compact slope plots.

    Output files:
    - {rel(OUT_DIR / (STEM + ".png"))}
    - {rel(OUT_DIR / (STEM + ".pdf"))}
    - {rel(OUT_DIR / (STEM + ".svg"))}
    - {rel(notes_path)}

    Input files:
    - regulon statistics: {rel(inputs["stats"])}
    - TF shortlist: {rel(inputs["shortlist"])}
    - AUCell matched matrix: {rel(inputs["auc"])}
    - metadata: {rel(inputs["metadata"])}
    - main Table 2 reviewed as compatible shortlist reference: manuscript_output/tables_main/Table2_shortlist_main_submission.csv
    - Supplementary Table ST01/S2 reviewed as compatible fallback references.

    Metadata / obs-equivalent fields used:
    - group: {obs_fields["group"]}
    - sample: {obs_fields["sample"]}
    - donor: {obs_fields["donor"]}
    - h5ad was not read because matched AUCell CSV plus metadata already contained the required cell/sample/group/donor mapping.

    Regulons:
    - All-regulon Panel B background: {len(stats)} regulons.
    - Panel A selected regulons: {len(selected_regs)}.
    - Panel A internal-control-associated regulons:
      {", ".join(control_regs)}
    - Panel A lesion-associated regulons:
      {", ".join(lesion_regs)}

    Panel B x-axis definition:
    - Mean AUC difference = mean AUC in lesion - mean AUC in internal_control.
    - Negative values indicate internal-control-associated regulon activity; positive values indicate lesion-associated activity.
    - Horizontal dashed line marks FDR = 0.05.

    Panel C sample-level mean AUC table:
    {textwrap.indent(format_table_for_notes(panel_c_table), "    ")}

    Panel D TF prioritization table:
    {textwrap.indent(format_table_for_notes(panel_d_notes), "    ")}

    Missing columns / replacements / warnings:
    {textwrap.indent(chr(10).join("- " + item for item in warnings) if warnings else "- None.", "    ")}

    Figure 2 v3.1 uses existing pySCENIC/AUCell outputs only; no pySCENIC or CellOracle rerun was performed.
    """
    notes_path.write_text(textwrap.dedent(notes).strip() + "\n", encoding="utf-8")


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    inputs = {
        "stats": find_input(
            ["analysis_outputs/group_compare/regulon_group_statistics.csv", "manuscript_output/tables_supplementary/TableS1_full_differential_regulons_submission.csv"],
            ["regulon", "statistics"],
            warnings,
        ),
        "shortlist": find_input(
            ["analysis_outputs/celloracle_candidates/celloracle_tf_shortlist.csv", "manuscript_output/tables_supplementary/TableS2_candidate_tf_shortlist_extended_submission.csv"],
            ["tf", "shortlist"],
            warnings,
        ),
        "auc": find_input(
            ["final_exports/auc_mtx_matched_to_h5ad.csv", "final_exports/key_tables/auc_mtx_matched_to_h5ad.csv", "output/auc_mtx.csv"],
            ["auc", "mtx"],
            warnings,
        ),
        "metadata": find_input(
            ["celloracle_run/prepared_data/celloracle_round1_metadata.csv"],
            ["metadata"],
            warnings,
        ),
    }

    stats = load_and_standardize_stats(inputs["stats"], warnings)
    shortlist = load_and_standardize_shortlist(inputs["shortlist"], stats, warnings)
    auc_meta, obs_fields = load_auc_with_metadata(inputs["auc"], inputs["metadata"], warnings)

    sample_meta = sample_order_table(auc_meta)
    selected_regs, control_regs, lesion_regs = choose_panel_a_regulons(stats, set(auc_meta.columns), warnings)
    sample_means = auc_meta.groupby("sample")[selected_regs].mean().loc[sample_meta["sample"]]
    heat = zscore_rows(sample_means.T)

    axis_regs = ["NFE2L2(+)", "THRB(+)"]
    axis_sample_means = (
        auc_meta.groupby(["sample", "donor", "group"], dropna=False)[axis_regs]
        .mean()
        .reset_index()
        .melt(id_vars=["sample", "donor", "group"], var_name="regulon", value_name="mean_auc")
    )
    axis_sample_means["group_rank"] = axis_sample_means["group"].map({"internal_control": 0, "lesion": 1}).fillna(9)
    axis_sample_means = axis_sample_means.sort_values(["regulon", "donor", "group_rank", "sample"]).drop(columns="group_rank")

    scorecard = build_scorecard(shortlist, stats, warnings)

    div_cmap = LinearSegmentedColormap.from_list("teal_white_burgundy", [CONTROL_COLOR, "#F7F7F7", LESION_COLOR])

    fig = plt.figure(figsize=(14.5, 9.5), constrained_layout=False)
    gs = fig.add_gridspec(
        3,
        2,
        left=0.055,
        right=0.985,
        top=0.955,
        bottom=0.075,
        width_ratios=[0.82, 1.0],
        height_ratios=[0.38, 0.25, 0.37],
        wspace=0.23,
        hspace=0.54,
    )
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_d = fig.add_subplot(gs[2, 1])

    plot_panel_a(fig, ax_a, heat, sample_meta, control_regs, lesion_regs, div_cmap)
    plot_panel_b(ax_b, stats)
    plot_panel_c(fig, gs[1, 1], axis_sample_means)
    plot_panel_d(ax_d, scorecard, div_cmap)

    for ext, kwargs in {
        "png": {"dpi": 600},
        "pdf": {},
        "svg": {},
    }.items():
        fig.savefig(OUT_DIR / f"{STEM}.{ext}", bbox_inches="tight", **kwargs)
    plt.close(fig)

    write_notes(
        OUT_DIR / f"{STEM}_notes.txt",
        inputs,
        obs_fields,
        selected_regs,
        control_regs,
        lesion_regs,
        axis_sample_means,
        scorecard,
        stats,
        warnings,
    )

    print(f"Wrote {rel(OUT_DIR / (STEM + '.png'))}")
    print(f"Wrote {rel(OUT_DIR / (STEM + '.pdf'))}")
    print(f"Wrote {rel(OUT_DIR / (STEM + '.svg'))}")
    print(f"Wrote {rel(OUT_DIR / (STEM + '_notes.txt'))}")
    if warnings:
        print("Warnings:")
        for item in warnings:
            print(f"- {item}")


if __name__ == "__main__":
    main()
