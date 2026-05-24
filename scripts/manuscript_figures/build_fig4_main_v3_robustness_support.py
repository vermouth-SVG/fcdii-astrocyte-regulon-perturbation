from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_main"
STEM = "Fig4_main_v3_robustness_support"
CAPTION_TITLE = "Fig. 4. Sample-level robustness and supportive external evidence for the NFE2L2–THRB dual-axis model"

TF_ORDER = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
AXIS_ORDER = ["NFE2L2", "THRB"]
GSE190_ORDER = ["NFE2L2", "THRB", "SOX2", "BHLHE40"]

LESION_COLOR = "#9A1F2D"
CONTROL_COLOR = "#0B6670"
BHLHE40_COLOR = "#6A51A3"
SOX2_COLOR = "#8C8C8C"
NEUTRAL_GREY = "#D9D9D9"
LIGHT_BG = "#F5F5F5"
GRID_GREY = "#ECECEC"
TEXT_GREY = "#333333"

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
    "SOX2": "#F1F1F1",
}

ROBUSTNESS_CLASS = {
    "NFE2L2": "lesion-axis stable",
    "THRB": "control-axis anchor",
    "BHLHE40": "secondary",
    "SOX2": "retained",
}

SUPPORT_CLASS_190 = {
    "NFE2L2": "cross-syndrome support",
    "THRB": "partial control-like support",
    "SOX2": "partial support",
    "BHLHE40": "weak/attenuated",
}

INPUTS = {
    "robust_integrated": ROOT / "robustness_validation" / "05_robustness_validation_integrated_table.csv",
    "shortlist_stability": ROOT / "robustness_validation" / "04_shortlist_priority_stability.csv",
    "loso_master": ROOT / "robustness_validation" / "02_leave_one_sample_out_master_results_revised.csv",
    "metadata": ROOT / "celloracle_run" / "prepared_data" / "celloracle_round1_metadata.csv",
    "gse140_support": ROOT / "external_validation" / "single_group_support" / "external_support_ranking.csv",
    "gse190_support": ROOT / "external_validation_round2" / "external_round2_support_ranking.csv",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.4,
            "axes.titlesize": 8.8,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.7,
            "ytick.labelsize": 6.7,
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


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def normalize_col(name: str) -> str:
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def find_col(df: pd.DataFrame, candidates: list[str], label: str, warnings: list[str], required: bool = True) -> str | None:
    cols = {normalize_col(col): col for col in df.columns}
    for candidate in candidates:
        key = normalize_col(candidate)
        if key in cols:
            return cols[key]
    if required:
        warnings.append(f"Missing column for {label}; candidates={candidates}; available={list(df.columns)}")
    return None


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.105,
        1.075,
        label,
        transform=ax.transAxes,
        fontsize=12.0,
        fontweight="bold",
        ha="left",
        va="top",
        clip_on=False,
    )


def style_axes(ax: plt.Axes, grid_axis: str | None = None) -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_GREY, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def load_metadata(warnings: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    path = INPUTS["metadata"]
    meta = read_table(path)
    first = meta.columns[0]
    if str(first).startswith("Unnamed") or str(first) == "":
        meta = meta.rename(columns={first: "cell_id"})
    sample_col = find_col(meta, ["sample_id", "sample", "orig.ident"], "metadata sample", warnings)
    label_col = find_col(meta, ["clinical_sample_label", "sample_label", "sample_id"], "metadata sample label", warnings)
    group_col = find_col(meta, ["group", "condition", "disease_group"], "metadata group", warnings, required=False)
    cols = {"sample": sample_col or "sample_id", "label": label_col or "clinical_sample_label"}
    if group_col:
        cols["group"] = group_col
    return meta, cols


def load_panel_a_data(warnings: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    robust = read_table(INPUTS["robust_integrated"])
    stability = read_table(INPUTS["shortlist_stability"])

    tf_col_r = find_col(robust, ["tf", "TF", "gene"], "robustness TF", warnings)
    expr_col = find_col(robust, ["loo_expr_consistency_rate", "LOSO_expr_rate", "expression_loso_direction"], "LOSO expression direction", warnings)
    reg_col = find_col(robust, ["loo_regulon_consistency_rate", "LOSO_regulon_rate", "regulon_loso_direction"], "LOSO regulon direction", warnings)
    ret_col = find_col(robust, ["shortlist_retention_frequency", "retention_frequency", "retention", "shortlist_retention"], "shortlist retention", warnings)

    tf_col_s = find_col(stability, ["tf", "TF", "gene"], "shortlist-stability TF", warnings)
    rank_col = find_col(
        stability,
        ["priority_stability_score", "rank_stability_component", "rank_stability", "median_rank_when_retained", "median_rank"],
        "rank stability",
        warnings,
    )
    median_rank_col = find_col(stability, ["median_rank_when_retained", "median_rank", "rank_when_retained"], "median rank", warnings, required=False)

    robust = robust.set_index(tf_col_r)
    stability = stability.set_index(tf_col_s)
    rows = []
    for tf in TF_ORDER:
        row = {
            "tf": tf,
            "expr_loso": pd.to_numeric(robust.loc[tf, expr_col], errors="coerce") if tf in robust.index else np.nan,
            "regulon_loso": pd.to_numeric(robust.loc[tf, reg_col], errors="coerce") if tf in robust.index else np.nan,
            "shortlist_retention": pd.to_numeric(robust.loc[tf, ret_col], errors="coerce") if tf in robust.index else np.nan,
            "rank_stability": pd.to_numeric(stability.loc[tf, rank_col], errors="coerce") if tf in stability.index else np.nan,
            "median_rank": pd.to_numeric(stability.loc[tf, median_rank_col], errors="coerce") if median_rank_col and tf in stability.index else np.nan,
            "robustness_class": ROBUSTNESS_CLASS[tf],
        }
        rows.append(row)
    df = pd.DataFrame(rows)
    used_cols = {
        "TF": tf_col_r,
        "Expr. LOSO direction": expr_col,
        "Regulon LOSO direction": reg_col,
        "Shortlist retention": ret_col,
        "Rank stability": f"{rel(INPUTS['shortlist_stability'])}::{rank_col}",
        "Median rank annotation": median_rank_col or "not available",
    }
    return df, used_cols


def load_panel_b_data(warnings: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict[str, str], bool]:
    loso = read_table(INPUTS["loso_master"])
    meta, meta_cols = load_metadata(warnings)
    id_map = meta[[meta_cols["sample"], meta_cols["label"]]].drop_duplicates()

    tf_col = find_col(loso, ["tf", "TF", "gene"], "LOSO TF", warnings)
    drop_col = find_col(loso, ["dropped_sample", "leave_out_sample", "sample_left_out"], "LOSO dropped sample", warnings)
    expr_col = find_col(loso, ["expr_log2fc", "expression_log2FC", "expr_log2FC", "log2FC"], "LOSO expression effect", warnings)
    reg_col = find_col(loso, ["regulon_diff", "auc_diff", "mean_auc_diff", "regulon_delta_auc"], "LOSO regulon effect", warnings)

    loso = loso.merge(id_map, left_on=drop_col, right_on=meta_cols["sample"], how="left")
    label_col = meta_cols["label"]
    loso["drop_label"] = np.where(loso[drop_col].astype(str).eq("FULL_DATA"), "Full", loso[label_col].fillna(loso[drop_col]))
    label_map = id_map.sort_values(label_col)
    order_cols = ["Full"] + [f"- {str(label).replace('.', '_')}" for label in label_map[label_col].tolist()]
    rename_map = {
        str(label): f"- {str(label).replace('.', '_')}"
        for label in label_map[label_col].tolist()
    }
    loso["plot_label"] = loso["drop_label"].map(rename_map).fillna(loso["drop_label"])

    expr = loso.pivot_table(index=tf_col, columns="plot_label", values=expr_col, aggfunc="first").reindex(index=AXIS_ORDER, columns=order_cols)
    reg = loso.pivot_table(index=tf_col, columns="plot_label", values=reg_col, aggfunc="first").reindex(index=AXIS_ORDER, columns=order_cols)

    row_index = [
        "NFE2L2 expression log2FC",
        "NFE2L2 regulon ΔAUC",
        "THRB expression log2FC",
        "THRB regulon ΔAUC",
    ]
    raw = pd.DataFrame(
        [
            expr.loc["NFE2L2"].to_numpy(dtype=float),
            reg.loc["NFE2L2"].to_numpy(dtype=float),
            expr.loc["THRB"].to_numpy(dtype=float),
            reg.loc["THRB"].to_numpy(dtype=float),
        ],
        index=row_index,
        columns=order_cols,
    )

    scaled = raw.copy()
    for idx in scaled.index:
        max_abs = np.nanmax(np.abs(scaled.loc[idx].to_numpy(dtype=float)))
        if np.isfinite(max_abs) and max_abs > 0:
            scaled.loc[idx] = scaled.loc[idx] / max_abs
        else:
            scaled.loc[idx] = 0.0
    used_cols = {
        "TF": tf_col,
        "dropped_sample": drop_col,
        "expression log2FC": expr_col,
        "regulon ΔAUC": reg_col,
        "sample label": f"{rel(INPUTS['metadata'])}::{label_col}",
        "x-axis columns": ", ".join(order_cols),
    }
    return raw, scaled, order_cols, used_cols, True


def load_panel_c_data(warnings: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    g140 = read_table(INPUTS["gse140_support"])
    tf_col = find_col(g140, ["tf", "TF", "gene"], "GSE140393 TF", warnings)
    support_col = find_col(g140, ["support_score", "external_support_score"], "GSE140393 support score", warnings)
    pos_col = find_col(g140, ["overall_positive_fraction", "positive_fraction", "pos_frac"], "GSE140393 positive fraction", warnings, required=False)
    mean_expr_col = find_col(g140, ["overall_mean_norm_expr", "mean_norm_expr", "mean_expression"], "GSE140393 mean expression", warnings, required=False)
    axis_col = find_col(g140, ["expected_higher_group_main", "Expected_axis", "expected_axis"], "GSE140393 expected axis", warnings, required=False)

    g140 = g140[g140[tf_col].isin(TF_ORDER)].copy()
    g140["tf"] = g140[tf_col].astype(str)
    g140["support_score_plot"] = as_num(g140[support_col])
    if pos_col:
        g140["point_size_metric"] = as_num(g140[pos_col])
        size_note = pos_col
    elif mean_expr_col:
        g140["point_size_metric"] = as_num(g140[mean_expr_col])
        size_note = mean_expr_col
        warnings.append("GSE140393 positive_fraction was unavailable; mean_norm_expr was used for point size.")
    else:
        g140["point_size_metric"] = 1.0
        size_note = "uniform size"
        warnings.append("GSE140393 positive_fraction/mean_norm_expr columns were unavailable; Panel C uses uniform point size.")
    if axis_col:
        g140["expected_axis"] = g140[axis_col].astype(str)
    else:
        g140["expected_axis"] = g140["tf"].map({"NFE2L2": "lesion", "THRB": "internal_control", "BHLHE40": "lesion", "SOX2": "lesion"})
    g140 = g140.sort_values("support_score_plot", ascending=False)
    used_cols = {
        "TF": tf_col,
        "Support score": support_col,
        "Point size": size_note,
        "Expected axis": axis_col or "inferred from fixed manuscript axis",
    }
    return g140, used_cols


def load_panel_d_data(warnings: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    g190 = read_table(INPUTS["gse190_support"])
    tf_col = find_col(g190, ["tf", "TF", "gene"], "GSE190452 TF", warnings)
    logfc_col = find_col(g190, ["group_log2fc_lesion_vs_internal_control", "group_log2FC", "group_log2fc"], "GSE190452 group log2FC", warnings)
    cluster_col = find_col(g190, ["cluster_support_fraction", "cluster_support"], "GSE190452 cluster support", warnings)
    direction_col = find_col(g190, ["group_direction_consistent_with_main", "direction_status", "direction_consistency"], "GSE190452 direction", warnings)
    support_col = find_col(g190, ["support_score", "external_support_score"], "GSE190452 support score", warnings, required=False)

    g190 = g190[g190[tf_col].isin(GSE190_ORDER)].copy()
    g190["tf"] = g190[tf_col].astype(str)
    g190["group_log2fc"] = as_num(g190[logfc_col])
    g190["cluster_support"] = as_num(g190[cluster_col])
    if g190[direction_col].dtype == bool:
        g190["direction_consistent"] = g190[direction_col]
    else:
        g190["direction_consistent"] = g190[direction_col].astype(str).str.lower().isin(["true", "1", "yes", "consistent", "direction-consistent"])
    g190["support_score"] = as_num(g190[support_col]) if support_col else np.nan
    g190["support_class"] = g190["tf"].map(SUPPORT_CLASS_190)
    g190 = g190.set_index("tf").loc[GSE190_ORDER].reset_index()
    used_cols = {
        "TF": tf_col,
        "Group log2FC": logfc_col,
        "Cluster support": cluster_col,
        "Direction": direction_col,
        "Support score": support_col or "not displayed",
    }
    return g190, used_cols


def bar_score(ax: plt.Axes, x: float, y: float, width: float, value: float, color: str, alpha: float = 0.85) -> None:
    ax.add_patch(Rectangle((x, y - 0.075), width, 0.15, facecolor="#E9E9E9", edgecolor="none", zorder=1))
    fill = 0.0 if pd.isna(value) else max(0.0, min(float(value), 1.0))
    ax.add_patch(Rectangle((x, y - 0.075), width * fill, 0.15, facecolor=color, edgecolor="none", alpha=alpha, zorder=2))
    text = "" if pd.isna(value) else f"{value:.2f}"
    ax.text(x + width + 0.008, y, text, fontsize=6.0, ha="left", va="center", color=TEXT_GREY)


def plot_panel_a(ax: plt.Axes, data: pd.DataFrame) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.65, len(TF_ORDER) + 0.58)
    ax.axis("off")
    ax.set_title("Internal sample-level robustness", pad=8)

    metric_cols = [
        ("Expr. LOSO\ndirection", "expr_loso"),
        ("Regulon LOSO\ndirection", "regulon_loso"),
        ("Shortlist\nretention", "shortlist_retention"),
        ("Rank\nstability", "rank_stability"),
    ]
    tf_x = 0.02
    col_x = [0.20, 0.36, 0.52, 0.68]
    bar_w = 0.08
    class_x = 0.84
    header_y = len(TF_ORDER) + 0.15

    ax.text(tf_x, header_y, "TF", fontsize=7.0, fontweight="bold", ha="left", va="bottom")
    for (label, _), x0 in zip(metric_cols, col_x):
        ax.text(x0 + bar_w / 2, header_y, label, fontsize=6.1, fontweight="bold", ha="center", va="bottom")
    ax.text(class_x, header_y, "Robustness\nclass", fontsize=6.3, fontweight="bold", ha="left", va="bottom")

    indexed = data.set_index("tf")
    for i, tf in enumerate(TF_ORDER):
        y = len(TF_ORDER) - 1 - i
        ax.add_patch(
            Rectangle(
                (0.0, y - 0.33),
                1.0,
                0.66,
                facecolor=ROW_BG[tf],
                edgecolor="white",
                linewidth=0.9,
                alpha=0.80 if tf in {"NFE2L2", "THRB"} else 0.65,
                zorder=0,
            )
        )
        ax.text(
            tf_x,
            y,
            tf,
            fontsize=7.2,
            color=TF_COLORS[tf],
            fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal",
            ha="left",
            va="center",
        )
        for (_, col), x0 in zip(metric_cols, col_x):
            bar_score(ax, x0, y, bar_w, float(indexed.loc[tf, col]), TF_COLORS[tf], alpha=0.9 if tf in {"NFE2L2", "THRB"} else 0.68)
        median_rank = indexed.loc[tf, "median_rank"]
        if pd.notna(median_rank):
            ax.text(col_x[-1] + bar_w + 0.054, y - 0.19, f"median {median_rank:.0f}", fontsize=5.4, color="#696969", ha="right", va="center")
        ax.text(class_x, y, ROBUSTNESS_CLASS[tf], fontsize=6.5, ha="left", va="center", color=TEXT_GREY)
    panel_label(ax, "A")


def signed_effect_cmap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("teal_white_burgundy", [CONTROL_COLOR, "#FAFAFA", LESION_COLOR], N=256)


def plot_panel_b(ax: plt.Axes, raw: pd.DataFrame, scaled: pd.DataFrame, order_cols: list[str], fig: plt.Figure) -> None:
    cmap = signed_effect_cmap()
    im = ax.imshow(scaled.values, aspect="auto", cmap=cmap, vmin=-1, vmax=1)
    ax.set_title("Leave-one-sample-out axis stability", pad=8)
    ax.set_xticks(np.arange(len(order_cols)))
    ax.set_xticklabels(order_cols, rotation=28, ha="right")
    ax.set_yticks(np.arange(raw.shape[0]))
    ax.set_yticklabels(raw.index)
    for tick in ax.get_yticklabels():
        text = tick.get_text()
        if text.startswith("NFE2L2"):
            tick.set_color(LESION_COLOR)
            tick.set_fontweight("bold")
        elif text.startswith("THRB"):
            tick.set_color(CONTROL_COLOR)
            tick.set_fontweight("bold")
    for i in range(raw.shape[0]):
        for j in range(raw.shape[1]):
            value = raw.iloc[i, j]
            if pd.isna(value):
                label = ""
            elif "regulon" in raw.index[i]:
                label = f"{value:.3f}"
            else:
                label = f"{value:.2f}"
            ax.text(j, i, label, ha="center", va="center", fontsize=5.8, color="white" if abs(float(scaled.iloc[i, j])) > 0.62 else "#2A2A2A")
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("row-scaled effect direction")
    cbar.ax.tick_params(labelsize=6.0)
    panel_label(ax, "B")


def scale_sizes(values: pd.Series, min_size: float = 42, max_size: float = 160) -> np.ndarray:
    vals = as_num(values).to_numpy(dtype=float)
    if not np.isfinite(vals).any() or np.nanmax(vals) == np.nanmin(vals):
        return np.full(len(vals), (min_size + max_size) / 2)
    scaled = (vals - np.nanmin(vals)) / (np.nanmax(vals) - np.nanmin(vals))
    return min_size + (max_size - min_size) * scaled


def plot_panel_c(ax: plt.Axes, g140: pd.DataFrame) -> None:
    y = np.arange(len(g140))[::-1]
    vals = g140["support_score_plot"].to_numpy(dtype=float)
    colors = [TF_COLORS[tf] for tf in g140["tf"]]
    sizes = scale_sizes(g140["point_size_metric"], 42, 150)
    ax.hlines(y, 0, vals, color=NEUTRAL_GREY, linewidth=1.15, zorder=1)
    ax.scatter(vals, y, s=sizes, c=colors, edgecolor="white", linewidth=0.8, zorder=3)
    for yi, (_, row) in zip(y, g140.iterrows()):
        tf = row["tf"]
        ax.text(
            float(row["support_score_plot"]) + 0.025,
            yi,
            tf,
            color=TF_COLORS[tf],
            fontsize=7.0,
            fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal",
            va="center",
            ha="left",
        )
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlim(0, 1.15)
    ax.set_ylim(-0.6, len(g140) - 0.35)
    ax.set_xlabel("Support score")
    ax.set_title("GSE140393 single-group support", pad=8)
    ax.text(0.02, 0.94, "lesion-only supportive dataset", transform=ax.transAxes, fontsize=6.8, color="#595959", ha="left", va="top")
    style_axes(ax, grid_axis="x")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="#9E9E9E", markeredgecolor="white", markersize=4.5, label="low positive fraction"),
            Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="#9E9E9E", markeredgecolor="white", markersize=8.0, label="high positive fraction"),
        ],
        frameon=False,
        loc="lower right",
    )
    panel_label(ax, "C")


def signed_bar(ax: plt.Axes, x0: float, y: float, width: float, value: float) -> None:
    center = x0 + width / 2
    ax.add_patch(Rectangle((x0, y - 0.075), width, 0.15, facecolor="#EEEEEE", edgecolor="none", zorder=1))
    ax.plot([center, center], [y - 0.10, y + 0.10], color="#BFBFBF", linewidth=0.7, zorder=2)
    max_abs = 0.60
    frac = 0.0 if pd.isna(value) else min(abs(float(value)) / max_abs, 1.0)
    if pd.isna(value):
        return
    if value >= 0:
        ax.add_patch(Rectangle((center, y - 0.075), (width / 2) * frac, 0.15, facecolor=LESION_COLOR, edgecolor="none", alpha=0.86, zorder=3))
    else:
        ax.add_patch(Rectangle((center - (width / 2) * frac, y - 0.075), (width / 2) * frac, 0.15, facecolor=CONTROL_COLOR, edgecolor="none", alpha=0.86, zorder=3))
    ax.text(x0 + width + 0.015, y, f"{value:.2f}", fontsize=6.0, ha="left", va="center", color=TEXT_GREY)


def cluster_bar(ax: plt.Axes, x0: float, y: float, width: float, value: float, tf: str) -> None:
    ax.add_patch(Rectangle((x0, y - 0.075), width, 0.15, facecolor="#EEEEEE", edgecolor="none", zorder=1))
    fill = 0.0 if pd.isna(value) else max(0.0, min(float(value), 1.0))
    ax.add_patch(Rectangle((x0, y - 0.075), width * fill, 0.15, facecolor=TF_COLORS[tf], edgecolor="none", alpha=0.82 if tf != "BHLHE40" else 0.55, zorder=2))
    ax.text(x0 + width + 0.015, y, f"{value:.2f}", fontsize=6.0, ha="left", va="center", color=TEXT_GREY)


def plot_panel_d(ax: plt.Axes, g190: pd.DataFrame) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.65, len(GSE190_ORDER) + 0.58)
    ax.axis("off")
    ax.set_title("GSE190452 cross-syndrome support", pad=8)

    tf_x = 0.02
    logfc_x, logfc_w = 0.19, 0.16
    cluster_x, cluster_w = 0.44, 0.14
    direction_x = 0.67
    class_x = 0.76
    header_y = len(GSE190_ORDER) + 0.15
    ax.text(tf_x, header_y, "TF", fontsize=7.0, fontweight="bold", ha="left", va="bottom")
    ax.text(logfc_x + logfc_w / 2, header_y, "Group\nlog2FC", fontsize=6.5, fontweight="bold", ha="center", va="bottom")
    ax.text(cluster_x + cluster_w / 2, header_y, "Cluster\nsupport", fontsize=6.5, fontweight="bold", ha="center", va="bottom")
    ax.text(direction_x, header_y, "Direction", fontsize=6.5, fontweight="bold", ha="center", va="bottom")
    ax.text(class_x, header_y, "Support class", fontsize=6.5, fontweight="bold", ha="left", va="bottom")

    indexed = g190.set_index("tf")
    for i, tf in enumerate(GSE190_ORDER):
        y = len(GSE190_ORDER) - 1 - i
        alpha = 0.78 if tf in {"NFE2L2", "THRB"} else 0.58
        ax.add_patch(Rectangle((0.0, y - 0.33), 1.0, 0.66, facecolor=ROW_BG[tf], edgecolor="white", linewidth=0.9, alpha=alpha, zorder=0))
        ax.text(
            tf_x,
            y,
            tf,
            color=TF_COLORS[tf],
            fontsize=7.2,
            fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal",
            ha="left",
            va="center",
        )
        signed_bar(ax, logfc_x, y, logfc_w, float(indexed.loc[tf, "group_log2fc"]))
        cluster_bar(ax, cluster_x, y, cluster_w, float(indexed.loc[tf, "cluster_support"]), tf)
        direction_symbol = "✓" if bool(indexed.loc[tf, "direction_consistent"]) else "×"
        direction_color = "#2E7D32" if bool(indexed.loc[tf, "direction_consistent"]) else "#8C8C8C"
        ax.text(direction_x, y, direction_symbol, fontsize=8.6, fontweight="bold", color=direction_color, ha="center", va="center")
        ax.text(class_x, y, SUPPORT_CLASS_190[tf], fontsize=6.3, ha="left", va="center", color=TEXT_GREY)

    ax.legend(
        handles=[
            Line2D([0], [0], marker="$✓$", linestyle="None", color="#2E7D32", markersize=7.0, label="direction-consistent"),
            Line2D([0], [0], marker="$×$", linestyle="None", color="#8C8C8C", markersize=7.0, label="direction-attenuated"),
        ],
        frameon=False,
        loc="lower right",
        bbox_to_anchor=(0.99, -0.02),
    )
    panel_label(ax, "D")


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


def compact_table(df: pd.DataFrame) -> str:
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].map(lambda value: "" if pd.isna(value) else f"{float(value):.4f}")
    return out.to_csv(index=False)


def write_notes(
    outputs: dict[str, Path],
    panel_cols: dict[str, dict[str, str]],
    panel_a: pd.DataFrame,
    panel_b_raw: pd.DataFrame,
    panel_b_cols: list[str],
    panel_c: pd.DataFrame,
    panel_d: pd.DataFrame,
    warnings: list[str],
) -> None:
    notes = f"""
    {STEM}
    ==============================

    Caption title:
    {CAPTION_TITLE}

    Input files:
    - {rel(INPUTS["robust_integrated"])}
    - {rel(INPUTS["shortlist_stability"])}
    - {rel(INPUTS["loso_master"])}
    - {rel(INPUTS["metadata"])}
    - {rel(INPUTS["gse140_support"])}
    - {rel(INPUTS["gse190_support"])}

    Panel A data and columns:
    - Table: {rel(INPUTS["robust_integrated"])} plus {rel(INPUTS["shortlist_stability"])}
    - Columns: {panel_cols["A"]}
    - Robustness metrics: Expr. LOSO direction = leave-one-sample-out expression direction consistency; Regulon LOSO direction = leave-one-sample-out regulon direction consistency; Shortlist retention = shortlist_retention_frequency; Rank stability = normalized priority/rank stability score from the shortlist stability table.
    - Raw median rank when retained was retained only as a small annotation, not as the main visual score, to avoid over-penalizing the NFE2L2 lesion-axis candidate.
    - Panel A values:
    {textwrap.indent(compact_table(panel_a[["tf", "expr_loso", "regulon_loso", "shortlist_retention", "rank_stability", "median_rank", "robustness_class"]]).strip(), "    ")}

    Panel B data and columns:
    - Table: {rel(INPUTS["loso_master"])}
    - Columns: {panel_cols["B"]}
    - LOSO effect rows: NFE2L2 expression log2FC, NFE2L2 regulon ΔAUC, THRB expression log2FC, THRB regulon ΔAUC.
    - Columns on x-axis: {", ".join(panel_b_cols)}
    - Expression log2FC and regulon ΔAUC have different numeric scales, so Panel B uses row-scaled signed effects: each row was divided by its maximum absolute value across Full and leave-one-sample-out conditions. Sign was preserved; positive values indicate lesion-associated direction and negative values indicate internal-control-associated direction.
    - Raw effect table:
    {textwrap.indent(compact_table(panel_b_raw.reset_index().rename(columns={"index": "effect"})).strip(), "    ")}

    Panel C data and columns:
    - Table: {rel(INPUTS["gse140_support"])}
    - Columns: {panel_cols["C"]}
    - GSE140393 is treated as a single-group supportive analysis, not formal validation.
    - Panel C is a ranked lollipop using Support score on the x-axis; point size uses positive fraction when available.
    - Panel C values:
    {textwrap.indent(compact_table(panel_c[["tf", "support_score_plot", "point_size_metric", "expected_axis"]]).strip(), "    ")}

    Panel D data and columns:
    - Table: {rel(INPUTS["gse190_support"])}
    - Columns: {panel_cols["D"]}
    - GSE190452 is treated as a cross-syndrome supportive analysis, not formal external validation.
    - Group log2FC uses direction color: positive = disease/TLE-associated, negative = control-like. Cluster support is shown as a horizontal bar. Direction uses ✓ for direction-consistent and × for direction-attenuated.
    - Panel D values:
    {textwrap.indent(compact_table(panel_d[["tf", "group_log2fc", "cluster_support", "direction_consistent", "support_class"]]).strip(), "    ")}

    Required interpretation notes:
    Figure 4 v3 uses existing robustness and external supportive-analysis outputs only; no upstream analysis was rerun.
    Robustness is reported as sample-level / leave-one-sample-out robustness, not donor-level robustness.
    GSE140393 is described as single-group supportive analysis.
    GSE190452 is described as cross-syndrome supportive analysis.
    No functional enrichment, drug repositioning, CellOracle vector field, pySCENIC rerun, CellOracle rerun, or newly downloaded data were used.

    Warnings / substitutions:
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
    warnings: list[str] = []

    panel_a, cols_a = load_panel_a_data(warnings)
    panel_b_raw, panel_b_scaled, panel_b_cols, cols_b, _ = load_panel_b_data(warnings)
    panel_c, cols_c = load_panel_c_data(warnings)
    panel_d, cols_d = load_panel_d_data(warnings)

    fig = plt.figure(figsize=(14.0, 9.0), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.05], hspace=0.34, wspace=0.24)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    plot_panel_a(ax_a, panel_a)
    plot_panel_b(ax_b, panel_b_raw, panel_b_scaled, panel_b_cols, fig)
    plot_panel_c(ax_c, panel_c)
    plot_panel_d(ax_d, panel_d)

    outputs = save_outputs(fig)
    write_notes(
        outputs,
        {"A": cols_a, "B": cols_b, "C": cols_c, "D": cols_d},
        panel_a,
        panel_b_raw,
        panel_b_cols,
        panel_c,
        panel_d,
        warnings,
    )

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
