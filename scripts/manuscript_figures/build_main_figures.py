#!/usr/bin/env python3
from __future__ import annotations

import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FIG_ROOT = ROOT / "manuscript_output" / "figures_main"

TFS_MAIN = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
TFS_SHORTLIST = ["NFE2L2", "THRB", "BHLHE40", "SOX2", "SATB2", "RARB", "HMGA1"]
REGULONS_MAIN = [f"{tf}(+)" for tf in TFS_MAIN]
GROUP_ORDER = ["internal_control", "lesion"]

COLORS = {
    "lesion": "#b54a43",
    "internal_control": "#3f6fa6",
    "neutral": "#666666",
    "NFE2L2": "#b54a43",
    "THRB": "#3f6fa6",
    "BHLHE40": "#d88a30",
    "SOX2": "#7b6fa8",
    "SATB2": "#789262",
    "RARB": "#5a8f8c",
    "HMGA1": "#9a6a4a",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.titlesize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_csv(path: str | Path, **kwargs) -> pd.DataFrame:
    p = ROOT / path if not isinstance(path, Path) else path
    return pd.read_csv(p, **kwargs)


def clean_first_col(df: pd.DataFrame, name: str = "regulon") -> pd.DataFrame:
    first = df.columns[0]
    if first.startswith("Unnamed") or first in {"", "H1"}:
        df = df.rename(columns={first: name})
    return df


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def neg_log10(values: pd.Series | np.ndarray) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").astype(float)
    x = np.where(np.isfinite(x) & (x > 0), x, np.nan)
    finite = x[np.isfinite(x)]
    floor = max(float(np.nanmin(finite)) * 0.1, 1e-300) if finite.size else 1e-300
    x = np.where(np.isfinite(x), x, floor)
    return -np.log10(x)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, fontweight="bold", fontsize=11, va="top", ha="left")


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.png", dpi=360, bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")


def style_axes(ax: plt.Axes, grid_axis: str | None = "y") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color="#e6e6e6", linewidth=0.6, zorder=0)
    ax.tick_params(width=0.6, length=3)


def load_auc_with_metadata() -> pd.DataFrame:
    auc = pd.read_csv(ROOT / "final_exports" / "auc_mtx_matched_to_h5ad.csv")
    auc = auc.rename(columns={auc.columns[0]: "cell_id"})
    meta = pd.read_csv(ROOT / "celloracle_run" / "prepared_data" / "celloracle_round1_metadata.csv")
    first = meta.columns[0]
    if first.startswith("Unnamed"):
        meta = meta.rename(columns={first: "cell_id"})
    if "cell_id" not in meta.columns and "join_key" in meta.columns:
        meta["cell_id"] = meta["join_key"]
    cols = ["cell_id", "group", "sample_id", "donor_id", "clinical_sample_label", "sample_prefix", "subtype"]
    cols = [c for c in cols if c in meta.columns]
    merged = auc.merge(meta[cols], on="cell_id", how="left")
    return merged


def plot_heatmap(ax: plt.Axes, data: pd.DataFrame, title: str, cmap: str = "RdBu_r", vmin: float = -1, vmax: float = 1) -> None:
    im = ax.imshow(data.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_xticklabels(data.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_yticklabels(data.index)
    ax.set_title(title, pad=6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return im


def build_fig2() -> dict[str, list[str]]:
    out_dir = FIG_ROOT / "Fig2"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = [
        "analysis_outputs/group_compare/regulon_group_statistics.csv",
        "analysis_outputs/group_compare/group_mean_auc_row_zscore.csv",
        "analysis_outputs/group_compare/top20_differential_regulons.csv",
        "analysis_outputs/celloracle_candidates/celloracle_tf_shortlist.csv",
        "final_exports/auc_mtx_matched_to_h5ad.csv",
        "celloracle_run/prepared_data/celloracle_round1_metadata.csv",
    ]

    stats = read_csv(inputs[0])
    stats["lesion_minus_internal_control"] = as_num(stats["mean_lesion"]) - as_num(stats["mean_internal_control"])
    stats["neg_log10_fdr"] = neg_log10(stats["fdr_bh"])

    heat = clean_first_col(read_csv(inputs[1]), "regulon").set_index("regulon")
    top20 = read_csv(inputs[2])
    top_regs = top20["regulon"].astype(str).tolist()
    heat = heat.loc[[r for r in top_regs if r in heat.index], [c for c in GROUP_ORDER if c in heat.columns]]

    shortlist = read_csv(inputs[3])
    shortlist["neg_log10_regulon_fdr"] = neg_log10(shortlist["regulon_fdr"])
    shortlist["log2fc"] = as_num(shortlist["log2fc_lesion_vs_internal_control"])
    shortlist["posfrac"] = as_num(shortlist["positive_fraction_in_higher_group"])

    auc_meta = load_auc_with_metadata()

    fig = plt.figure(figsize=(12.2, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[0.92, 1.08])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    im = plot_heatmap(ax_a, heat, "Differential regulon activity (top 20)")
    panel_label(ax_a, "A")
    cbar = fig.colorbar(im, ax=ax_a, fraction=0.046, pad=0.02)
    cbar.set_label("row z-score")

    x = stats["lesion_minus_internal_control"].astype(float)
    y = stats["neg_log10_fdr"]
    colors = np.where(x >= 0, COLORS["lesion"], COLORS["internal_control"])
    ax_b.scatter(x, y, c=colors, s=18, alpha=0.78, edgecolors="none")
    ax_b.axvline(0, color="#888888", lw=0.8)
    ax_b.set_xlabel("Mean AUC difference (lesion - internal_control)")
    ax_b.set_ylabel("-log10(FDR)")
    ax_b.set_title("Differential regulon summary")
    style_axes(ax_b)
    for reg in [f"{tf}(+)" for tf in TFS_SHORTLIST]:
        row = stats[stats["regulon"].astype(str).eq(reg)]
        if row.empty:
            continue
        xx = float(row["lesion_minus_internal_control"].iloc[0])
        yy = float(row["neg_log10_fdr"].iloc[0])
        ax_b.scatter([xx], [yy], s=44, facecolor="white", edgecolor=COLORS.get(reg.replace("(+)", ""), "#222222"), linewidth=1.2, zorder=5)
        ax_b.annotate(reg, (xx, yy), xytext=(4 if xx >= 0 else -4, 4), textcoords="offset points", ha="left" if xx >= 0 else "right", fontsize=7)
    panel_label(ax_b, "B")

    positions = np.arange(len(REGULONS_MAIN))
    width = 0.32
    for i, group in enumerate(GROUP_ORDER):
        vals = []
        for reg in REGULONS_MAIN:
            if reg in auc_meta.columns:
                vals.append(pd.to_numeric(auc_meta.loc[auc_meta["group"].eq(group), reg], errors="coerce").dropna().values)
            else:
                vals.append(np.array([]))
        bp = ax_c.boxplot(
            vals,
            positions=positions + (i - 0.5) * width,
            widths=0.26,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 0.8},
            boxprops={"linewidth": 0.7},
            whiskerprops={"linewidth": 0.7},
            capprops={"linewidth": 0.7},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(COLORS[group])
            patch.set_alpha(0.58)
        ax_c.scatter([], [], color=COLORS[group], label=group)
    ax_c.set_xticks(positions)
    ax_c.set_xticklabels(REGULONS_MAIN, rotation=25, ha="right")
    ax_c.set_ylabel("Regulon AUC")
    ax_c.set_title("Key regulon activity by group")
    ax_c.legend(frameon=False, loc="upper right")
    style_axes(ax_c)
    panel_label(ax_c, "C")

    effect_colors = shortlist["regulon_effect_group"].map(COLORS).fillna("#777777")
    sizes = 80 + 420 * shortlist["posfrac"].fillna(0)
    ax_d.scatter(shortlist["log2fc"], shortlist["neg_log10_regulon_fdr"], s=sizes, c=effect_colors, alpha=0.80, edgecolor="white", linewidth=0.7)
    ax_d.axvline(0, color="#888888", lw=0.8)
    ax_d.set_xlabel("TF expression log2FC (lesion - internal_control)")
    ax_d.set_ylabel("-log10(regulon FDR)")
    ax_d.set_title("Candidate TF shortlist prioritization")
    style_axes(ax_d)
    for _, r in shortlist.iterrows():
        tf = str(r["tf"])
        x0 = float(r["log2fc"])
        y0 = float(r["neg_log10_regulon_fdr"])
        ax_d.annotate(tf, (x0, y0), xytext=(4 if x0 >= 0 else -4, 3), textcoords="offset points", ha="left" if x0 >= 0 else "right", fontsize=7, fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal")
    ax_d.scatter([], [], s=80 + 420 * 0.25, color="#aaaaaa", alpha=0.7, label="size: positive fraction")
    ax_d.legend(frameon=False, loc="lower right")
    panel_label(ax_d, "D")

    fig.suptitle("Fig. 2. Differential regulon landscape and candidate TF prioritization in the astrocyte pilot", y=1.02)
    save_figure(fig, out_dir, "Fig2_main")

    notes = """
    Fig2_main
    =========
    Data sources:
    - analysis_outputs/group_compare/regulon_group_statistics.csv
    - analysis_outputs/group_compare/group_mean_auc_row_zscore.csv
    - analysis_outputs/group_compare/top20_differential_regulons.csv
    - analysis_outputs/celloracle_candidates/celloracle_tf_shortlist.csv
    - final_exports/auc_mtx_matched_to_h5ad.csv
    - celloracle_run/prepared_data/celloracle_round1_metadata.csv

    Script entry:
    - scripts/manuscript_figures/build_main_figures.py

    Panels:
    A. Row-z scored mean regulon AUC for the top 20 differential regulons, showing lesion-high and internal_control-high structure.
    B. Differential regulon summary plot using mean AUC difference and regulon FDR. NFE2L2(+), THRB(+), BHLHE40(+), SOX2(+), SATB2(+), RARB(+), and HMGA1(+) are labeled when present.
    C. Box plots of single-cell regulon AUC for NFE2L2(+), THRB(+), BHLHE40(+), and SOX2(+). Group labels come from CellOracle round1 metadata matched to exported AUCell matrix.
    D. Candidate TF shortlist scatter using TF expression log2FC, regulon FDR, and positive fraction in the higher-expression group.

    Simplification:
    - Panel C uses exported AUC CSV plus metadata instead of reading h5ad, because the local anndata import fails due to a torch DLL initialization error.
    """
    write_text(out_dir / "Fig2_notes.txt", notes)
    return {"Fig2": inputs}


def build_fig3() -> dict[str, list[str]]:
    out_dir = FIG_ROOT / "Fig3"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = [
        "celloracle_run/ko_round1/round1_recovery_index_ranking.csv",
        "celloracle_run/ko_round1/round1_ko_summary_by_group.csv",
        "celloracle_run/ko_round1/NFE2L2/NFE2L2_quiver.png",
        "celloracle_run/ko_round1/THRB/THRB_quiver.png",
    ]
    ranking = read_csv(inputs[0])
    ranking["mean_shift_length"] = as_num(ranking["mean_shift_length"])
    ranking["recovery_index"] = as_num(ranking["recovery_index"])
    by_group = read_csv(inputs[1])
    by_group["mean_shift_length"] = as_num(by_group["mean_shift_length"])
    by_group["mean_delta_x"] = as_num(by_group["mean_delta_x"])
    by_group["mean_delta_y"] = as_num(by_group["mean_delta_y"])

    fig = plt.figure(figsize=(12.4, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[0.92, 1.08], width_ratios=[1, 1, 1.1])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d1 = fig.add_subplot(gs[1, 0:2])
    ax_d2 = fig.add_subplot(gs[1, 2])

    overall = ranking.sort_values("mean_shift_length", ascending=True)
    ax_a.barh(overall["tf"], overall["mean_shift_length"], color=[COLORS.get(tf, "#777777") for tf in overall["tf"]])
    ax_a.set_xlabel("Mean shift length")
    ax_a.set_title("Overall perturbation ranking")
    style_axes(ax_a)
    panel_label(ax_a, "A")

    ri = ranking.sort_values("recovery_index", ascending=True)
    ax_b.barh(ri["tf"], ri["recovery_index"], color=[COLORS.get(tf, "#777777") for tf in ri["tf"]])
    ax_b.set_xlabel("Recovery Index")
    ax_b.set_title("Approximate Recovery Index ranking")
    style_axes(ax_b)
    panel_label(ax_b, "B")

    pivot = by_group.pivot(index="tf", columns="group", values="mean_shift_length").reindex(TFS_MAIN)
    x = np.arange(len(pivot.index))
    width = 0.34
    for i, group in enumerate(GROUP_ORDER):
        vals = pivot[group].astype(float).values if group in pivot.columns else np.zeros(len(x))
        ax_c.bar(x + (i - 0.5) * width, vals, width=width, label=group, color=COLORS[group], alpha=0.82)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(pivot.index, rotation=25, ha="right")
    ax_c.set_ylabel("Mean shift length")
    ax_c.set_title("Group-specific perturbation magnitude")
    ax_c.legend(frameon=False)
    style_axes(ax_c)
    panel_label(ax_c, "C")

    for ax, img_rel, title, label in [
        (ax_d1, inputs[2], "NFE2L2 representative state-shift vectors", "D"),
        (ax_d2, inputs[3], "THRB representative state-shift vectors", "E"),
    ]:
        img_path = ROOT / img_rel
        if img_path.exists():
            ax.imshow(mpimg.imread(img_path))
            ax.set_title(title, pad=5)
        else:
            ax.text(0.5, 0.5, f"Missing image\\n{img_rel}", ha="center", va="center")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        panel_label(ax, label)

    fig.suptitle("Fig. 3. CellOracle perturbation prioritization and integrated ranking of candidate TFs", y=1.02)
    save_figure(fig, out_dir, "Fig3_main")

    notes = """
    Fig3_main
    =========
    Data sources:
    - celloracle_run/ko_round1/round1_recovery_index_ranking.csv
    - celloracle_run/ko_round1/round1_ko_summary_by_group.csv
    - celloracle_run/ko_round1/NFE2L2/NFE2L2_quiver.png
    - celloracle_run/ko_round1/THRB/THRB_quiver.png

    Script entry:
    - scripts/manuscript_figures/build_main_figures.py

    Panels:
    A. Overall CellOracle perturbation ranking by mean_shift_length. The observed order is NFE2L2 > THRB > BHLHE40 > SOX2.
    B. Approximate Recovery Index ranking. The observed order is THRB > NFE2L2 > BHLHE40 > SOX2.
    C. Group-specific mean shift length in lesion and internal_control cells.
    D-E. Existing CellOracle quiver/state-shift outputs for NFE2L2 and THRB embedded as representative vector-field panels.

    Interpretation boundary:
    - Recovery Index is shown as the current state-shift-derived approximate pathological reversal index, not as a strict fully re-derived geometric projection model.
    """
    write_text(out_dir / "Fig3_notes.txt", notes)
    return {"Fig3": inputs}


def build_fig4() -> dict[str, list[str]]:
    out_dir = FIG_ROOT / "Fig4"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = [
        "robustness_validation/05_robustness_validation_integrated_table.csv",
        "robustness_validation/04_shortlist_priority_stability.csv",
        "external_validation/single_group_support/external_support_ranking.csv",
        "external_validation_round2/external_round2_support_ranking.csv",
    ]
    rob = read_csv(inputs[0]).set_index("tf")
    rob_cols = ["loo_expr_consistency_rate", "loo_regulon_consistency_rate", "shortlist_retention_frequency"]
    rob_mat = rob.reindex(TFS_MAIN)[rob_cols].apply(pd.to_numeric, errors="coerce")
    stability = read_csv(inputs[1]).set_index("tf").reindex(TFS_MAIN)
    gse140 = read_csv(inputs[2]).set_index("tf").reindex(TFS_MAIN)
    gse190 = read_csv(inputs[3]).set_index("tf").reindex(TFS_MAIN)

    fig = plt.figure(figsize=(12.2, 8.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    im = ax_a.imshow(rob_mat.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax_a.set_xticks(np.arange(len(rob_cols)))
    ax_a.set_xticklabels(["LOSO expr", "LOSO regulon", "shortlist retention"], rotation=25, ha="right")
    ax_a.set_yticks(np.arange(len(rob_mat.index)))
    ax_a.set_yticklabels(rob_mat.index)
    ax_a.set_title("Sample-level robustness summary")
    for i in range(rob_mat.shape[0]):
        for j in range(rob_mat.shape[1]):
            val = rob_mat.iloc[i, j]
            ax_a.text(j, i, f"{val:.2f}" if pd.notna(val) else "NA", ha="center", va="center", fontsize=7, color="white" if pd.notna(val) and val > 0.65 else "#222222")
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    fig.colorbar(im, ax=ax_a, fraction=0.046, pad=0.02, label="consistency / frequency")
    panel_label(ax_a, "A")

    retention = pd.to_numeric(stability["retention_frequency"], errors="coerce")
    median_rank = pd.to_numeric(stability["median_rank_when_retained"], errors="coerce")
    ax_b.bar(stability.index, retention, color=[COLORS.get(tf, "#777777") for tf in stability.index], alpha=0.78)
    ax_b.set_ylim(0, 1.08)
    ax_b.set_ylabel("Retention frequency")
    ax_b.set_title("Shortlist sensitivity")
    ax_b.tick_params(axis="x", rotation=25)
    ax_b2 = ax_b.twinx()
    ax_b2.plot(stability.index, median_rank, color="#333333", marker="o", lw=1.2, label="median rank")
    ax_b2.set_ylabel("Median rank when retained")
    ax_b2.invert_yaxis()
    style_axes(ax_b)
    ax_b2.grid(False)
    panel_label(ax_b, "B")

    support = pd.to_numeric(gse140["support_score"], errors="coerce")
    ax_c.barh(gse140.index[::-1], support.iloc[::-1], color=[COLORS.get(tf, "#777777") for tf in gse140.index[::-1]], alpha=0.82)
    ax_c.set_xlim(0, 1.05)
    ax_c.set_xlabel("Support score")
    ax_c.set_title("GSE140393 single-group supportive analysis")
    style_axes(ax_c)
    for tf in gse140.index:
        val = float(gse140.loc[tf, "support_score"])
        interp = str(gse140.loc[tf, "expected_higher_group_main"])
        ax_c.text(val + 0.025, list(gse140.index[::-1]).index(tf), interp, va="center", fontsize=7)
    panel_label(ax_c, "C")

    x = pd.to_numeric(gse190["group_log2fc_lesion_vs_internal_control"], errors="coerce")
    y = pd.to_numeric(gse190["cluster_support_fraction"], errors="coerce")
    s = 120 + 450 * pd.to_numeric(gse190["support_score"], errors="coerce").fillna(0)
    ax_d.scatter(x, y, s=s, c=[COLORS.get(tf, "#777777") for tf in gse190.index], alpha=0.82, edgecolor="white", linewidth=0.8)
    ax_d.axvline(0, color="#999999", lw=0.8)
    ax_d.axhline(0.5, color="#cccccc", lw=0.8, ls="--")
    ax_d.set_xlabel("Group log2FC (TLE - non-epileptic control)")
    ax_d.set_ylabel("Cluster direction support fraction")
    ax_d.set_title("GSE190452 cross-syndrome supportive analysis")
    style_axes(ax_d)
    for tf in gse190.index:
        ax_d.annotate(tf, (float(x.loc[tf]), float(y.loc[tf])), xytext=(4, 3), textcoords="offset points", fontsize=7, fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal")
    panel_label(ax_d, "D")

    fig.suptitle("Fig. 4. Sample-level robustness and supportive external evidence for the NFE2L2–THRB dual-axis model", y=1.02)
    save_figure(fig, out_dir, "Fig4_main")

    notes = """
    Fig4_main
    =========
    Data sources:
    - robustness_validation/05_robustness_validation_integrated_table.csv
    - robustness_validation/04_shortlist_priority_stability.csv
    - external_validation/single_group_support/external_support_ranking.csv
    - external_validation_round2/external_round2_support_ranking.csv

    Script entry:
    - scripts/manuscript_figures/build_main_figures.py

    Panels:
    A. Sample-level robustness matrix summarizing leave-one-sample-out expression consistency, leave-one-sample-out regulon consistency, and shortlist retention frequency.
    B. Shortlist sensitivity panel showing retention frequency and median rank when retained.
    C. GSE140393 single-group supportive expression analysis. This is lesion-direction support only, not formal validation.
    D. GSE190452 cross-syndrome supportive analysis using group direction and cluster support fraction. This is supportive evidence, not formal same-disease validation.
    """
    write_text(out_dir / "Fig4_notes.txt", notes)
    return {"Fig4": inputs}


def build_fig1_insets() -> dict[str, list[str]]:
    out_dir = FIG_ROOT / "Fig1_insets"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = ["celloracle_run/prepared_data/celloracle_round1_metadata.csv"]
    meta = pd.read_csv(ROOT / inputs[0])
    first = meta.columns[0]
    if first.startswith("Unnamed"):
        meta = meta.rename(columns={first: "cell_id"})

    sample_counts = meta.groupby(["clinical_sample_label", "group"], dropna=False).size().reset_index(name="n_cells")
    group_counts = meta.groupby("group", dropna=False).size().reindex(GROUP_ORDER).fillna(0).astype(int)

    fig = plt.figure(figsize=(8.6, 3.0), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    sample_counts = sample_counts.sort_values(["group", "clinical_sample_label"])
    ax_a.bar(sample_counts["clinical_sample_label"], sample_counts["n_cells"], color=[COLORS.get(g, "#777777") for g in sample_counts["group"]])
    ax_a.set_ylabel("Cells")
    ax_a.set_title("4-sample astrocyte pilot")
    ax_a.tick_params(axis="x", rotation=35)
    style_axes(ax_a)
    panel_label(ax_a, "A")

    ax_b.bar(group_counts.index, group_counts.values, color=[COLORS.get(g, "#777777") for g in group_counts.index])
    ax_b.set_ylabel("Cells")
    ax_b.set_title("Group counts")
    for i, v in enumerate(group_counts.values):
        ax_b.text(i, int(v) + max(group_counts.values) * 0.03, str(int(v)), ha="center", va="bottom", fontsize=8)
    style_axes(ax_b)
    panel_label(ax_b, "B")

    ax_c.axis("off")
    summary = [
        ("Cells", f"{len(meta):,}"),
        ("Samples", f"{meta['clinical_sample_label'].nunique()}"),
        ("Donors", f"{meta['donor_id'].nunique()}"),
        ("Cell type", "Astrocyte"),
    ]
    y = 0.85
    for k, v in summary:
        ax_c.text(0.05, y, k, ha="left", va="center", color="#555555")
        ax_c.text(0.95, y, v, ha="right", va="center", fontweight="bold")
        y -= 0.18
    ax_c.set_title("Cohort mini-summary")
    panel_label(ax_c, "C")

    fig.suptitle("Fig. 1 quantitative insets", y=1.06)
    save_figure(fig, out_dir, "fig1_insets")
    notes = """
    fig1_insets
    ===========
    Data source:
    - celloracle_run/prepared_data/celloracle_round1_metadata.csv

    Script entry:
    - scripts/manuscript_figures/build_main_figures.py

    Panels:
    A. Cell counts for the 4-sample GSE268807 astrocyte pilot.
    B. lesion vs internal_control cell counts.
    C. Basic cohort mini-summary for embedding into the Fig. 1 schematic.
    """
    write_text(out_dir / "Fig1_insets_notes.txt", notes)
    return {"Fig1_insets": inputs}


def build_fig5_insets() -> dict[str, list[str]]:
    out_dir = FIG_ROOT / "Fig5_insets"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = [
        "functional_interpretation/04_program_convergence_table.csv",
        "drug_repositioning/16_paper_ready_table_after_bbb_v21.csv",
    ]
    conv = read_csv(inputs[0]).set_index("tf").reindex(["NFE2L2", "THRB", "BHLHE40"])
    drug = read_csv(inputs[1])

    fig = plt.figure(figsize=(8.8, 3.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.25, 1.05])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    overlap = pd.to_numeric(conv["overlap_count"], errors="coerce")
    ax_a.bar(conv.index, overlap, color=[COLORS.get(tf, "#777777") for tf in conv.index], alpha=0.84)
    ax_a.set_ylabel("Overlap genes")
    ax_a.set_title("Regulon-DEG overlap")
    ax_a.tick_params(axis="x", rotation=25)
    style_axes(ax_a)
    panel_label(ax_a, "A")

    ratio = pd.to_numeric(conv["overlap_ratio"], errors="coerce")
    ax_b.bar(conv.index, ratio, color=[COLORS.get(tf, "#777777") for tf in conv.index], alpha=0.84)
    ax_b.set_ylim(0, 1.05)
    ax_b.set_ylabel("Overlap ratio")
    ax_b.set_title("Program convergence ratio")
    ax_b.tick_params(axis="x", rotation=25)
    style_axes(ax_b)
    for i, tf in enumerate(conv.index):
        theme = str(conv.loc[tf, "program_theme_2"] if pd.notna(conv.loc[tf, "program_theme_2"]) else conv.loc[tf, "program_theme_1"])
        ax_b.text(i, float(ratio.loc[tf]) + 0.04, theme.replace("_", "\n"), ha="center", va="bottom", fontsize=6)
    panel_label(ax_b, "B")

    counts = drug["after_bbb_final_layer"].value_counts()
    order = ["headline_cns_mechanism_direction_leads", "supportive_after_bbb_leads", "peripheral_program_modulating_clues"]
    counts = counts.reindex(order).fillna(0).astype(int)
    labels = ["headline", "supportive", "peripheral"]
    ax_c.bar(labels, counts.values, color=["#8c3f3a", "#b78b45", "#999999"], alpha=0.84)
    ax_c.set_ylabel("Candidates")
    ax_c.set_title("Exploratory drug clue tiers")
    ax_c.tick_params(axis="x", rotation=25)
    style_axes(ax_c)
    panel_label(ax_c, "C")

    fig.suptitle("Fig. 5 quantitative insets", y=1.06)
    save_figure(fig, out_dir, "fig5_insets")
    notes = """
    fig5_insets
    ===========
    Data sources:
    - functional_interpretation/04_program_convergence_table.csv
    - drug_repositioning/16_paper_ready_table_after_bbb_v21.csv

    Script entry:
    - scripts/manuscript_figures/build_main_figures.py

    Panels:
    A. Regulon-DEG overlap counts for NFE2L2, THRB, and BHLHE40.
    B. Program convergence overlap ratios with compact theme labels.
    C. Very small conservative summary of exploratory after-BBB drug clue tiers. This is not a drug discovery main figure and should be used only as a minor inset if needed.
    """
    write_text(out_dir / "Fig5_insets_notes.txt", notes)
    return {"Fig5_insets": inputs}


def write_summary(all_inputs: dict[str, list[str]]) -> None:
    summary = ["# Main figure build summary", "", "Generated by `scripts/manuscript_figures/build_main_figures.py`.", ""]
    for fig_name, inputs in all_inputs.items():
        summary.append(f"## {fig_name}")
        summary.append("")
        summary.append("Input files:")
        for item in inputs:
            summary.append(f"- `{item}`")
        summary.append("")
    summary.extend(
        [
            "## Simplifications",
            "",
            "- Existing upstream analyses were not rerun.",
            "- h5ad reading was avoided because local `anndata` import fails due to a torch DLL initialization error; exported AUC CSV and CellOracle metadata were used instead.",
            "- Existing CellOracle quiver PNGs were embedded for representative vector-field panels rather than regenerating CellOracle state-shift calculations.",
            "- External evidence is labeled as supportive analysis, not formal external validation.",
            "- Drug repositioning is kept as a small exploratory inset only; it is not used as a central main figure.",
            "",
            "## Supplementary candidates",
            "",
            "- Full drug repositioning tables should remain supplementary.",
            "- Full GO/KEGG dotplots and long enrichment tables should remain supplementary.",
            "- Full CellOracle per-TF state-shift images can be placed in supplementary material if space is limited.",
        ]
    )
    write_text(FIG_ROOT / "FIGURE_BUILD_SUMMARY_MAIN.md", "\n".join(summary))


def main() -> None:
    setup_style()
    all_inputs: dict[str, list[str]] = {}
    for builder in [build_fig2, build_fig3, build_fig4, build_fig1_insets, build_fig5_insets]:
        all_inputs.update(builder())
    write_summary(all_inputs)
    print(f"Figures written to: {FIG_ROOT}")


if __name__ == "__main__":
    main()
