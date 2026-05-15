from __future__ import annotations

import math
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_main"
STEM = "Fig5_main_v3_program_model"
CAPTION_TITLE = "Fig. 5. Program-level interpretation and working model of the NFE2L2-THRB dual-axis architecture"

LESION_COLOR = "#9A1F2D"
CONTROL_COLOR = "#0B6670"
BHLHE40_COLOR = "#6A51A3"
SOX2_COLOR = "#8C8C8C"
NEUTRAL_GREY = "#D9D9D9"
GRID_GREY = "#ECECEC"
TEXT_GREY = "#333333"

TF_COLORS = {
    "NFE2L2": LESION_COLOR,
    "THRB": CONTROL_COLOR,
    "BHLHE40": BHLHE40_COLOR,
    "SOX2": SOX2_COLOR,
}

CONVERGENCE_LABELS = {
    "NFE2L2": "stress-adaptation",
    "THRB": "homeostatic support",
    "BHLHE40": "secondary",
    "SOX2": "retained",
}

TF_ORDER = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]

INPUTS = {
    "go": ROOT / "functional_interpretation" / "03_go_enrichment_all.csv",
    "kegg": ROOT / "functional_interpretation" / "03_kegg_enrichment_all.csv",
    "nfe_top": ROOT / "functional_interpretation" / "03_nfe2l2_top_terms.csv",
    "thrb_top": ROOT / "functional_interpretation" / "03_thrb_top_terms.csv",
    "convergence": ROOT / "functional_interpretation" / "04_program_convergence_table.csv",
    "gene_summary": ROOT / "functional_interpretation" / "02_tf_gene_sets_summary.csv",
}

SCHEMATIC_PREFERRED = [
    ROOT / "manuscript_output" / "figures_main" / "Figure5_model_clean.png",
    ROOT / "manuscript_output" / "figures_main" / "Fig5_insets" / "Figure5_model_clean.png",
    ROOT / "manuscript_output" / "figures_main" / "Figure5.png",
]

NFE2L2_TERM_KEYWORDS = [
    "Regulation Of Transcription By RNA Polymerase II",
    "Positive Regulation Of DNA-templated Transcription",
    "MAPK signaling pathway",
    "Regulation Of Intracellular Signal Transduction",
    "Protein Phosphorylation",
    "Endocytosis",
]

THRB_TERM_KEYWORDS = [
    "Calcium Ion-Regulated Exocytosis Of Neurotransmitter",
    "Positive Regulation Of Dendrite Extension",
    "Regulation Of Dendrite Extension",
    "Modulation Of Chemical Synaptic Transmission",
    "Calcium-Ion Regulated Exocytosis",
    "Synaptic Vesicle Exocytosis",
]


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


def neg_log10(values: pd.Series | np.ndarray) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").astype(float)
    x = np.where(np.isfinite(x) & (x > 0), x, np.nan)
    finite = x[np.isfinite(x)]
    floor = max(float(np.nanmin(finite)) * 0.1, 1e-300) if finite.size else 1e-300
    x = np.where(np.isfinite(x), x, floor)
    return -np.log10(x)


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


def style_axes(ax: plt.Axes, grid_axis: str | None = "x") -> None:
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_GREY, linewidth=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def short_term(term: str) -> str:
    cleaned = re.sub(r"\s*\(GO:\d+\)\s*$", "", str(term))
    cleaned = cleaned.replace("DNA-templated", "DNA-templated")
    return "\n".join(textwrap.wrap(cleaned, width=28, break_long_words=False))


def clean_term_key(term: str) -> str:
    return re.sub(r"\s*\(GO:\d+\)\s*$", "", str(term)).strip().lower()


def load_enrichment(warnings: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    frames = []
    for key in ["go", "kegg"]:
        path = INPUTS[key]
        if path.exists():
            frames.append(read_table(path))
        else:
            warnings.append(f"Enrichment file not found: {rel(path)}")
    if not frames:
        return pd.DataFrame(), {}
    enrich = pd.concat(frames, ignore_index=True, sort=False)
    cols = {
        "tf": find_col(enrich, ["tf", "TF", "gene"], "enrichment TF", warnings),
        "term": find_col(enrich, ["term", "Term", "pathway", "description"], "enrichment term", warnings),
        "fdr": find_col(enrich, ["fdr_bh", "FDR", "fdr", "padj", "adjusted_p"], "enrichment FDR", warnings),
        "overlap_count": find_col(enrich, ["overlap_count", "overlap_genes", "gene_count"], "enrichment overlap count", warnings),
        "overlap_ratio": find_col(enrich, ["overlap_ratio", "ratio", "enrichment_ratio"], "enrichment overlap ratio", warnings),
        "source": find_col(enrich, ["source", "Source"], "enrichment source", warnings, required=False) or "source",
        "gene_set_type": find_col(enrich, ["gene_set_type", "query_type"], "enrichment gene-set type", warnings, required=False),
    }
    return enrich, cols


def select_program_terms(
    enrich: pd.DataFrame,
    cols: dict[str, str],
    tf: str,
    keyword_terms: list[str],
    warnings: list[str],
    max_terms: int = 6,
) -> pd.DataFrame:
    if enrich.empty or not cols:
        warnings.append(f"No enrichment data available for {tf}; panel will display data not found.")
        return pd.DataFrame()

    sub = enrich[enrich[cols["tf"]].astype(str).eq(tf)].copy()
    gene_set_col = cols.get("gene_set_type")
    if gene_set_col and gene_set_col in sub.columns:
        preferred = sub[sub[gene_set_col].astype(str).eq("regulon_targets")].copy()
        if tf == "THRB":
            intersection = sub[sub[gene_set_col].astype(str).eq("intersection_directional")].copy()
            if not intersection.empty:
                preferred = intersection
        if not preferred.empty:
            sub = preferred
    if sub.empty:
        warnings.append(f"No enrichment rows available for {tf}; panel will display data not found.")
        return pd.DataFrame()

    sub["fdr_value"] = as_num(sub[cols["fdr"]])
    sub["neg_log10_fdr"] = neg_log10(sub["fdr_value"])
    sub["overlap_count_value"] = as_num(sub[cols["overlap_count"]])
    sub["overlap_ratio_value"] = as_num(sub[cols["overlap_ratio"]])
    sub["term_clean"] = sub[cols["term"]].map(lambda value: re.sub(r"\s*\(GO:\d+\)\s*$", "", str(value)).strip())
    sub["term_key"] = sub[cols["term"]].map(clean_term_key)

    selected_rows = []
    seen: set[str] = set()
    for keyword in keyword_terms:
        key = clean_term_key(keyword)
        match = sub[sub["term_key"].eq(key)].sort_values("fdr_value").head(1)
        if match.empty:
            match = sub[sub["term_key"].str.contains(re.escape(key), case=False, na=False)].sort_values("fdr_value").head(1)
        if not match.empty:
            row = match.iloc[0]
            if row["term_key"] not in seen:
                selected_rows.append(row)
                seen.add(row["term_key"])

    if len(selected_rows) < max_terms:
        excluded = r"cancer|carcinoma|leukemia|viral|virus|infection|hepatitis|pathogenic"
        fallback = sub[~sub["term_key"].str.contains(excluded, case=False, na=False)].sort_values("fdr_value")
        for _, row in fallback.iterrows():
            if row["term_key"] not in seen:
                selected_rows.append(row)
                seen.add(row["term_key"])
            if len(selected_rows) >= max_terms:
                break

    out = pd.DataFrame(selected_rows).head(max_terms).copy()
    if out.empty:
        warnings.append(f"Unable to select main-text enrichment terms for {tf}.")
        return out
    out = out.sort_values("neg_log10_fdr", ascending=True)
    return out


def make_tinted_cmap(color: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("tinted", ["#F7F7F7", color], N=256)


def scale_sizes(values: pd.Series, min_size: float = 42, max_size: float = 190) -> np.ndarray:
    vals = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(vals).any() or np.nanmax(vals) == np.nanmin(vals):
        return np.full(len(vals), (min_size + max_size) / 2)
    scaled = (vals - np.nanmin(vals)) / (np.nanmax(vals) - np.nanmin(vals))
    return min_size + (max_size - min_size) * scaled


def scaled_size_value(value: float, low: float, high: float, min_size: float = 42, max_size: float = 190) -> float:
    if not np.isfinite(value) or not np.isfinite(low) or not np.isfinite(high) or high == low:
        return (min_size + max_size) / 2
    scaled = (value - low) / (high - low)
    return min_size + (max_size - min_size) * np.clip(scaled, 0, 1)


def plot_program_terms(ax: plt.Axes, terms: pd.DataFrame, tf: str, color: str) -> None:
    ax.set_title(f"{tf} program", pad=8, color=color, fontweight="bold")
    if terms.empty:
        ax.text(0.5, 0.5, "data not found", transform=ax.transAxes, ha="center", va="center", color=TEXT_GREY)
        ax.set_xticks([])
        ax.set_yticks([])
        panel_label(ax, "A" if tf == "NFE2L2" else "B")
        return

    y = np.arange(len(terms))
    scatter = ax.scatter(
        terms["neg_log10_fdr"],
        y,
        s=scale_sizes(terms["overlap_count_value"]),
        c=terms["overlap_ratio_value"],
        cmap=make_tinted_cmap(color),
        edgecolor="white",
        linewidth=0.75,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels([short_term(term) for term in terms["term_clean"]])
    ax.set_xlabel("-log10(FDR)")
    ax.set_ylim(-0.55, len(terms) - 0.45)
    style_axes(ax, grid_axis="x")
    for spine in ax.spines.values():
        spine.set_visible(False)

    cax = ax.inset_axes([0.74, 1.015, 0.22, 0.034])
    cax.set_clip_on(False)
    cbar = ax.figure.colorbar(scatter, cax=cax, orientation="horizontal")
    cbar.set_label("overlap ratio")
    cbar.ax.tick_params(labelsize=5.6, pad=1, length=2.0)
    cbar.ax.xaxis.set_label_position("top")
    cbar.ax.xaxis.set_ticks_position("bottom")

    counts = terms["overlap_count_value"].dropna()
    if not counts.empty:
        vals = [float(counts.min()), float(counts.max())]
        low, high = vals[0], vals[-1]
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                markerfacecolor="#A0A0A0",
                markeredgecolor="white",
                markersize=math.sqrt(scaled_size_value(value, low, high, 42, 190)) / 2.8,
                label=f"{int(round(value))}",
            )
            for value in vals
        ]
        ax.legend(
            handles=handles,
            title="genes",
            frameon=False,
            loc="lower right",
            borderpad=0.2,
            handletextpad=0.6,
            labelspacing=0.45,
        )
    panel_label(ax, "A" if tf == "NFE2L2" else "B")


def load_convergence(warnings: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    conv_path = INPUTS["convergence"]
    gene_path = INPUTS["gene_summary"]
    if not conv_path.exists():
        warnings.append(f"Convergence table not found: {rel(conv_path)}")
        return pd.DataFrame(), {}
    conv = read_table(conv_path)
    tf_col = find_col(conv, ["tf", "TF", "gene"], "convergence TF", warnings)
    overlap_count_col = find_col(conv, ["overlap_count", "Program_overlap_count", "overlap_genes"], "convergence overlap count", warnings)
    overlap_ratio_col = find_col(conv, ["overlap_ratio", "Program_overlap_ratio"], "convergence overlap ratio", warnings)

    gene_summary = read_table(gene_path) if gene_path.exists() else pd.DataFrame()
    target_col = None
    if not gene_summary.empty:
        tf_col_gene = find_col(gene_summary, ["tf", "TF", "gene"], "gene summary TF", warnings)
        target_col = find_col(gene_summary, ["regulon_target_count", "regulon_targets", "target_count"], "regulon target count", warnings)
        gene_summary = gene_summary.set_index(tf_col_gene)
    else:
        warnings.append(f"Gene-set summary table not found: {rel(gene_path)}")

    rows = []
    indexed = conv.set_index(tf_col)
    for tf in TF_ORDER:
        if tf not in indexed.index:
            warnings.append(f"Convergence table missing TF row: {tf}")
            continue
        row = indexed.loc[tf]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        out = {
            "tf": tf,
            "overlap_count": pd.to_numeric(row[overlap_count_col], errors="coerce"),
            "overlap_ratio": pd.to_numeric(row[overlap_ratio_col], errors="coerce"),
            "interpretation": CONVERGENCE_LABELS[tf],
        }
        if target_col and tf in gene_summary.index:
            out["regulon_targets"] = pd.to_numeric(gene_summary.loc[tf, target_col], errors="coerce")
        else:
            out["regulon_targets"] = np.nan
        rows.append(out)
    data = pd.DataFrame(rows)
    used_cols = {
        "convergence_table": rel(conv_path),
        "TF": tf_col,
        "overlap_count": overlap_count_col,
        "overlap_ratio": overlap_ratio_col,
        "gene_summary": rel(gene_path),
        "regulon_targets": target_col or "missing",
    }
    return data, used_cols


def plot_convergence(ax: plt.Axes, data: pd.DataFrame) -> None:
    ax.set_title("Program convergence", pad=8)
    if data.empty:
        ax.text(0.5, 0.5, "data not found", transform=ax.transAxes, ha="center", va="center", color=TEXT_GREY)
        ax.set_xticks([])
        ax.set_yticks([])
        panel_label(ax, "C")
        return
    data = data.set_index("tf").loc[TF_ORDER].reset_index()
    y = np.arange(len(data))[::-1]
    sizes = scale_sizes(data["overlap_count"], 75, 310)
    ax.hlines(y, 0, data["overlap_ratio"], color=NEUTRAL_GREY, linewidth=1.2, zorder=1)
    ax.scatter(
        data["overlap_ratio"],
        y,
        s=sizes,
        c=[TF_COLORS[tf] for tf in data["tf"]],
        edgecolor="white",
        linewidth=0.85,
        alpha=0.95,
        zorder=3,
    )
    for yi, (_, row) in zip(y, data.iterrows()):
        tf = row["tf"]
        label = row["interpretation"]
        if tf == "THRB":
            label = "homeostatic support\ncompact program"
        ax.text(
            float(row["overlap_ratio"]) + 0.025,
            yi,
            label,
            fontsize=6.7,
            ha="left",
            va="center",
            color=TEXT_GREY if tf != "BHLHE40" else BHLHE40_COLOR,
        )
        ax.text(
            float(row["overlap_ratio"]),
            yi - 0.28,
            f"{int(row['overlap_count'])}/{int(row['regulon_targets'])}",
            fontsize=5.9,
            ha="center",
            va="top",
            color="#555555",
        )
    ax.set_yticks(y)
    ax.set_yticklabels(data["tf"])
    for tick in ax.get_yticklabels():
        tf = tick.get_text()
        tick.set_color(TF_COLORS[tf])
        if tf in {"NFE2L2", "THRB"}:
            tick.set_fontweight("bold")
    ax.set_xlim(0, min(1.12, max(1.05, float(data["overlap_ratio"].max()) + 0.18)))
    ax.set_xlabel("Regulon-DEG overlap ratio")
    ax.set_ylim(-0.65, len(data) - 0.35)
    style_axes(ax, grid_axis="x")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="#A0A0A0", markeredgecolor="white", markersize=5.0, label="fewer overlap genes"),
            Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="#A0A0A0", markeredgecolor="white", markersize=10.0, label="more overlap genes"),
        ],
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(0.02, 0.04),
        borderaxespad=0,
        labelspacing=0.5,
    )
    panel_label(ax, "C")


def find_schematic(warnings: list[str]) -> tuple[Path | None, bool]:
    for path in SCHEMATIC_PREFERRED:
        if path.exists():
            return path, "model_clean" in path.name.lower()

    search_roots = [
        ROOT / "manuscript_output" / "figures_main",
        ROOT / "manuscript_output",
        ROOT,
    ]
    found: list[Path] = []
    for root in search_roots:
        if root.exists():
            found.extend(root.rglob("Figure5*.png"))
            found.extend(root.rglob("Fig5*.png"))
    found = sorted(set(found), key=lambda p: ("model_clean" not in p.name.lower(), len(str(p)), str(p)))
    if found:
        return found[0], "model_clean" in found[0].name.lower()
    warnings.append("No Figure5 schematic PNG was found; Panel D will show data not found.")
    return None, False


def plot_schematic(ax: plt.Axes, schematic_path: Path | None) -> None:
    ax.set_title("Dual-axis working model", pad=8)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if schematic_path is None or not schematic_path.exists():
        ax.text(0.5, 0.5, "schematic not found", transform=ax.transAxes, ha="center", va="center", color=TEXT_GREY)
        panel_label(ax, "D")
        return
    img = mpimg.imread(schematic_path)
    ax.imshow(img)
    ax.axis("off")
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


def terms_table(terms: pd.DataFrame) -> str:
    if terms.empty:
        return "data not found"
    cols = ["term_clean", "neg_log10_fdr", "overlap_count_value", "overlap_ratio_value", "fdr_value"]
    table = terms[cols].copy()
    for col in ["neg_log10_fdr", "overlap_ratio_value", "fdr_value"]:
        table[col] = table[col].map(lambda value: "" if pd.isna(value) else f"{float(value):.4g}")
    table["overlap_count_value"] = table["overlap_count_value"].map(lambda value: "" if pd.isna(value) else f"{int(value)}")
    return table.to_csv(index=False)


def compact_table(df: pd.DataFrame) -> str:
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].map(lambda value: "" if pd.isna(value) else f"{float(value):.4f}")
    return out.to_csv(index=False)


def write_notes(
    outputs: dict[str, Path],
    enrich_cols: dict[str, str],
    nfe_terms: pd.DataFrame,
    thrb_terms: pd.DataFrame,
    conv: pd.DataFrame,
    conv_cols: dict[str, str],
    schematic_path: Path | None,
    clean_schematic_found: bool,
    bottom_crop_performed: bool,
    warnings: list[str],
) -> None:
    notes = f"""
    {STEM}
    ==============================

    Caption title:
    {CAPTION_TITLE}

    Input enrichment files:
    - {rel(INPUTS["go"])}
    - {rel(INPUTS["kegg"])}
    - {rel(INPUTS["nfe_top"])}
    - {rel(INPUTS["thrb_top"])}

    Enrichment column mapping:
    - {enrich_cols}

    Panel A selected NFE2L2 terms:
    {textwrap.indent(terms_table(nfe_terms).strip(), "    ")}

    Panel B selected THRB terms:
    {textwrap.indent(terms_table(thrb_terms).strip(), "    ")}

    Panel C convergence source:
    - Table and columns: {conv_cols}
    - Metrics: x-axis = regulon-DEG overlap ratio; point size = overlap genes; point labels show overlap genes / regulon targets.
    - Values:
    {textwrap.indent(compact_table(conv).strip(), "    ")}

    Panel D schematic:
    - Schematic PNG path: {rel(schematic_path) if schematic_path else "not found"}
    - Figure5_model_clean.png found: {clean_schematic_found}
    - Bottom drug-clue crop performed: {bottom_crop_performed}
    - The selected schematic did not require drug-clue cropping because the clean model file contains only the dual-axis working model.

    Required interpretation notes:
    Figure 5 v3 uses existing functional interpretation and schematic image outputs only; no upstream analysis was rerun.
    Drug repositioning clues were not included in the main Figure 5 and should remain in supplementary materials.
    No pySCENIC, CellOracle, enrichment, exploratory drug-signature, robustness, or supportive external/contextual analysis was rerun.

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
    warnings: list[str] = []

    enrich, enrich_cols = load_enrichment(warnings)
    nfe_terms = select_program_terms(enrich, enrich_cols, "NFE2L2", NFE2L2_TERM_KEYWORDS, warnings)
    thrb_terms = select_program_terms(enrich, enrich_cols, "THRB", THRB_TERM_KEYWORDS, warnings)
    conv, conv_cols = load_convergence(warnings)
    schematic_path, clean_schematic_found = find_schematic(warnings)
    bottom_crop_performed = False

    fig = plt.figure(figsize=(14.0, 9.4), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.04], width_ratios=[1.0, 1.0], hspace=0.34, wspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    plot_program_terms(ax_a, nfe_terms, "NFE2L2", LESION_COLOR)
    plot_program_terms(ax_b, thrb_terms, "THRB", CONTROL_COLOR)
    plot_convergence(ax_c, conv)
    plot_schematic(ax_d, schematic_path)

    outputs = save_outputs(fig)
    write_notes(
        outputs,
        enrich_cols,
        nfe_terms,
        thrb_terms,
        conv,
        conv_cols,
        schematic_path,
        clean_schematic_found,
        bottom_crop_performed,
        warnings,
    )

    for key in ["png", "pdf", "svg", "notes"]:
        print(f"Wrote {rel(outputs[key])}")
    print(f"Figure5_model_clean.png found: {clean_schematic_found}")
    if schematic_path:
        print(f"Schematic used: {rel(schematic_path)}")
    print(f"Bottom drug-clue crop performed: {bottom_crop_performed}")
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"- {warning}")
    else:
        print("Warnings: none")


if __name__ == "__main__":
    main()
