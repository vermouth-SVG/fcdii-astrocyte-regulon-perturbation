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
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "manuscript_output" / "figures_main"
STEM = "Fig5_main_v4_program_model"
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
    "THRB": "compact homeostatic program",
    "BHLHE40": "secondary",
    "SOX2": "target-limited retained",
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

TERM_DISPLAY_LABELS = {
    term.lower(): label
    for term, label in {
        "Regulation Of Transcription By RNA Polymerase II": "RNA polymerase II transcription",
        "Positive Regulation Of DNA-templated Transcription": "DNA-templated transcription",
        "Positive Regulation Of Transcription By RNA Polymerase II": "RNA polymerase II transcription",
        "Negative Regulation Of DNA-templated Transcription": "DNA-templated transcription regulation",
        "Regulation Of DNA-templated Transcription": "DNA-templated transcription",
        "MAPK signaling pathway": "MAPK signaling",
        "Endocytosis": "Endocytosis",
        "Regulation Of Intracellular Signal Transduction": "Intracellular signaling",
        "Protein Phosphorylation": "Protein phosphorylation",
        "Positive Regulation Of Dendrite Extension": "Dendrite extension",
        "Regulation Of Dendrite Extension": "Dendrite extension",
        "Calcium Ion-Regulated Exocytosis Of Neurotransmitter": "Ca²⁺-regulated\nneurotransmitter exocytosis",
        "Calcium-Ion Regulated Exocytosis": "Ca²⁺-regulated exocytosis",
        "Modulation Of Chemical Synaptic Transmission": "Chemical synaptic transmission",
        "Synaptic Vesicle Exocytosis": "Synaptic vesicle exocytosis",
        "Regulation Of Regulated Secretory Pathway": "Regulated secretory pathway",
        "Positive Regulation Of Cell Communication": "Cell communication",
    }.items()
}

PROGRAM_TERM_ORDER = {
    "NFE2L2": [
        "RNA polymerase II transcription",
        "DNA-templated transcription",
        "MAPK signaling",
        "Intracellular signaling",
        "Protein phosphorylation",
        "Endocytosis",
    ],
    "THRB": [
        "Dendrite extension",
        "Ca²⁺-regulated\nneurotransmitter exocytosis",
        "Ca²⁺-regulated exocytosis",
        "Chemical synaptic transmission",
        "Synaptic vesicle exocytosis",
    ],
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
        ax.grid(True, axis=grid_axis, color=GRID_GREY, linewidth=0.45, alpha=0.55, zorder=0)
    ax.tick_params(width=0.6, length=2.8)


def wrap_two_lines(label: str, width: int = 27) -> str:
    if "\n" in label:
        parts = [part.strip() for part in label.splitlines() if part.strip()]
        return "\n".join(parts[:2])
    wrapped = textwrap.wrap(label, width=width, break_long_words=False)
    if len(wrapped) <= 2:
        return "\n".join(wrapped)
    return f"{wrapped[0]}\n{' '.join(wrapped[1:])}"


def compact_fallback_label(term: str) -> str:
    cleaned = re.sub(r"\s*\(GO:\d+\)\s*$", "", str(term)).strip()
    cleaned = cleaned.replace("DNA-templated", "DNA-templated")
    return cleaned[:1].upper() + cleaned[1:]


def add_display_labels(terms: pd.DataFrame, tf: str, unmatched_labels: list[str]) -> pd.DataFrame:
    out = terms.copy()
    labels = []
    for term in out.get("term_clean", pd.Series(dtype=str)).astype(str):
        label = TERM_DISPLAY_LABELS.get(clean_term_key(term))
        if label is None:
            label = compact_fallback_label(term)
            unmatched_labels.append(f"{tf}: {term} -> {label}")
        labels.append(wrap_two_lines(label))
    out["display_label"] = labels
    return out


def display_key(label: str) -> str:
    return " ".join(str(label).replace("\n", " ").split()).lower()


def curate_program_terms(terms: pd.DataFrame, tf: str, max_terms: int = 5) -> tuple[pd.DataFrame, int]:
    if terms.empty:
        return terms, 0
    out = terms.copy()
    out["display_key"] = out["display_label"].map(display_key)
    before = len(out)
    selected_rows = []
    seen: set[str] = set()

    for label in PROGRAM_TERM_ORDER.get(tf, []):
        key = display_key(label)
        match = out[out["display_key"].eq(key)].sort_values("fdr_value").head(1)
        if not match.empty and key not in seen:
            selected_rows.append(match.iloc[0])
            seen.add(key)
        if len(selected_rows) >= max_terms:
            break

    if len(selected_rows) < max_terms:
        for _, row in out.sort_values("fdr_value").iterrows():
            key = row["display_key"]
            if key not in seen:
                selected_rows.append(row)
                seen.add(key)
            if len(selected_rows) >= max_terms:
                break

    selected = pd.DataFrame(selected_rows).copy()
    if selected.empty:
        return selected, 0
    selected["display_label"] = selected["display_label"].map(lambda value: wrap_two_lines(str(value), width=26))
    removed = before - len(selected)
    return selected, max(0, removed)


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


def scale_sizes(values: pd.Series, min_size: float = 42, max_size: float = 190) -> np.ndarray:
    vals = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(vals).any() or np.nanmax(vals) == np.nanmin(vals):
        return np.full(len(vals), (min_size + max_size) / 2)
    scaled = (vals - np.nanmin(vals)) / (np.nanmax(vals) - np.nanmin(vals))
    return min_size + (max_size - min_size) * scaled


def nice_round(value: float) -> int:
    if value <= 10:
        return int(round(value))
    if value <= 100:
        return int(round(value / 10.0) * 10)
    return int(round(value / 25.0) * 25)


def size_legend_values(values: pd.Series, mode: str = "terms") -> list[int]:
    vals = pd.to_numeric(values, errors="coerce").dropna().astype(float)
    if vals.empty:
        return []
    low = int(round(vals.min()))
    high_raw = int(round(vals.max()))
    if mode == "convergence":
        candidates = [low, 50, 150, high_raw]
        return sorted({value for value in candidates if low <= value <= high_raw})
    unique = sorted({int(round(value)) for value in vals})
    if len(unique) <= 3:
        return unique
    high = max(high_raw, nice_round(vals.max()))
    mid = nice_round((low + vals.max()) / 2.0)
    return sorted({low, mid, high})


def scaled_size_value(value: float, low: float, high: float, min_size: float = 42, max_size: float = 190) -> float:
    if not np.isfinite(value) or not np.isfinite(low) or not np.isfinite(high) or high == low:
        return (min_size + max_size) / 2
    scaled = (value - low) / (high - low)
    return min_size + (max_size - min_size) * np.clip(scaled, 0, 1)


def plot_program_facet(
    ax: plt.Axes,
    terms: pd.DataFrame,
    title: str,
    color: str,
    size_low: float,
    size_high: float,
) -> None:
    ax.set_title(title, pad=4, color=color, fontweight="bold", fontsize=7.5)
    if terms.empty:
        ax.text(0.5, 0.5, "data not found", transform=ax.transAxes, ha="center", va="center", color=TEXT_GREY)
        ax.set_xticks([])
        ax.set_yticks([])
        return

    y = np.arange(len(terms))[::-1]
    sizes = np.array([scaled_size_value(float(value), size_low, size_high, 52, 210) for value in terms["overlap_count_value"]])
    ax.hlines(y, 0, terms["neg_log10_fdr"], color=NEUTRAL_GREY, linewidth=0.8, zorder=1)
    ax.scatter(
        terms["neg_log10_fdr"],
        y,
        s=sizes,
        c=color,
        edgecolor="white",
        linewidth=0.75,
        alpha=0.95,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(terms["display_label"].tolist())
    ax.set_xlabel("-log10(FDR)")
    ax.set_ylim(-0.55, len(terms) - 0.45)
    ax.set_xlim(0, max(float(terms["neg_log10_fdr"].max()) * 1.08, 1.0))
    style_axes(ax, grid_axis="x")
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_selected_programs(frame_ax: plt.Axes, ax_nfe: plt.Axes, ax_thrb: plt.Axes, nfe_terms: pd.DataFrame, thrb_terms: pd.DataFrame) -> bool:
    frame_ax.axis("off")
    frame_ax.text(0.5, 1.22, "Selected functional programs", transform=frame_ax.transAxes, ha="center", va="bottom", fontsize=8.9, fontweight="bold", clip_on=False)
    panel_label(frame_ax, "A")

    counts = pd.concat([nfe_terms.get("overlap_count_value", pd.Series(dtype=float)), thrb_terms.get("overlap_count_value", pd.Series(dtype=float))])
    size_low = float(counts.min()) if not counts.dropna().empty else 1.0
    size_high = float(counts.max()) if not counts.dropna().empty else 1.0

    plot_program_facet(ax_nfe, nfe_terms, "NFE2L2 lesion-associated program", LESION_COLOR, size_low, size_high)
    plot_program_facet(ax_thrb, thrb_terms, "THRB homeostatic-supportive program", CONTROL_COLOR, size_low, size_high)

    legend_values = [value for value in [3, 20, 50, int(round(size_high))] if size_low <= value <= size_high]
    if not legend_values:
        legend_values = size_legend_values(counts, mode="terms")
    legend_values = sorted(set(legend_values))
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor="#A0A0A0",
            markeredgecolor="white",
            markersize=math.sqrt(scaled_size_value(value, size_low, size_high, 52, 210)) / 2.8,
            label=f"{value}",
        )
        for value in legend_values
    ]
    frame_ax.legend(
        handles=handles,
        title="overlap genes",
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.52, -0.20),
        ncol=len(handles),
        handletextpad=0.6,
        columnspacing=1.0,
        borderaxespad=0,
    )
    return True


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


def draw_score_bar(ax: plt.Axes, x: float, y: float, width: float, value: float, max_value: float, color: str, alpha: float = 0.95) -> None:
    ax.add_patch(Rectangle((x, y - 0.055), width, 0.11, facecolor="#EFEFEF", edgecolor="none", zorder=1))
    frac = 0.0 if not np.isfinite(value) or max_value <= 0 else min(max(value / max_value, 0), 1)
    ax.add_patch(Rectangle((x, y - 0.055), width * frac, 0.11, facecolor=color, edgecolor="none", alpha=alpha, zorder=2))


def plot_convergence(ax: plt.Axes, data: pd.DataFrame) -> None:
    ax.set_title("Program convergence scorecard", pad=8, fontweight="bold")
    if data.empty:
        ax.text(0.5, 0.5, "data not found", transform=ax.transAxes, ha="center", va="center", color=TEXT_GREY)
        ax.set_xticks([])
        ax.set_yticks([])
        panel_label(ax, "B")
        return
    data = data.set_index("tf").loc[TF_ORDER].reset_index()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 5.0)
    ax.axis("off")

    panel_label(ax, "B")
    headers = [
        (0.02, "TF"),
        (0.23, "Regulon\ntargets"),
        (0.43, "Overlap\ngenes"),
        (0.62, "Overlap\nratio"),
        (0.81, "Interpretation"),
    ]
    for x, label in headers:
        ax.text(x, 4.62, label, fontsize=6.6, fontweight="bold", ha="left", va="top", color=TEXT_GREY)

    max_targets = float(data["regulon_targets"].max())
    max_overlap = float(data["overlap_count"].max())
    row_y = [3.82, 2.88, 1.94, 1.00]
    bg_colors = {
        "NFE2L2": "#F8EDEF",
        "THRB": "#E8F3F4",
        "BHLHE40": "#F1EEF8",
        "SOX2": "#F4F4F4",
    }

    for yi, (_, row) in zip(row_y, data.iterrows()):
        tf = row["tf"]
        alpha = 0.45 if tf == "SOX2" else 0.95
        ax.add_patch(Rectangle((0.0, yi - 0.36), 0.98, 0.70, facecolor=bg_colors[tf], edgecolor="none", alpha=0.95, zorder=0))
        ax.text(0.02, yi, tf, fontsize=7.3, fontweight="bold" if tf in {"NFE2L2", "THRB"} else "normal", color=TF_COLORS[tf], ha="left", va="center")

        draw_score_bar(ax, 0.23, yi, 0.105, float(row["regulon_targets"]), max_targets, TF_COLORS[tf], alpha=alpha)
        ax.text(0.345, yi, f"{int(row['regulon_targets'])}", fontsize=6.6, ha="left", va="center", color=TEXT_GREY if tf != "SOX2" else SOX2_COLOR)

        draw_score_bar(ax, 0.43, yi, 0.105, float(row["overlap_count"]), max_overlap, TF_COLORS[tf], alpha=alpha)
        ax.text(0.545, yi, f"{int(row['overlap_count'])}", fontsize=6.6, ha="left", va="center", color=TEXT_GREY if tf != "SOX2" else SOX2_COLOR)

        draw_score_bar(ax, 0.62, yi, 0.12, float(row["overlap_ratio"]), 1.0, TF_COLORS[tf], alpha=alpha)
        ax.text(0.75, yi, f"{float(row['overlap_ratio']):.2f}", fontsize=6.6, ha="left", va="center", color=TEXT_GREY if tf != "SOX2" else SOX2_COLOR)

        ax.text(0.81, yi, str(row["interpretation"]), fontsize=6.5, ha="left", va="center", color=TEXT_GREY if tf not in {"BHLHE40", "SOX2"} else TF_COLORS[tf])

    for yi in [3.35, 2.41, 1.47, 0.53]:
        ax.plot([0.0, 0.98], [yi, yi], color="#FFFFFF", linewidth=0.9, zorder=4)


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
    warnings.append("No Figure5 schematic PNG was found; Panel C will show data not found.")
    return None, False


def crop_white_border(img: np.ndarray, threshold: float = 0.985, padding_fraction: float = 0.04) -> tuple[np.ndarray, bool, str]:
    arr = img.astype(float)
    if arr.max() > 1.0:
        arr = arr / 255.0
    rgb = arr[..., :3]
    alpha = arr[..., 3] if arr.ndim == 3 and arr.shape[2] >= 4 else np.ones(arr.shape[:2])
    content = (alpha > 0.05) & (np.min(rgb, axis=2) < threshold)
    if not np.any(content):
        h, w = img.shape[:2]
        y_pad = int(round(h * 0.06))
        x_pad = int(round(w * 0.06))
        if h > 2 * y_pad and w > 2 * x_pad:
            return img[y_pad : h - y_pad, x_pad : w - x_pad], True, "fixed 6% crop fallback; automatic content mask was empty"
        return img, False, "automatic crop failed; image retained unchanged"

    ys, xs = np.where(content)
    y0, y1 = int(ys.min()), int(ys.max())
    x0, x1 = int(xs.min()), int(xs.max())
    h, w = img.shape[:2]
    bbox_h = max(y1 - y0 + 1, 1)
    bbox_w = max(x1 - x0 + 1, 1)
    y_pad = int(round(bbox_h * padding_fraction))
    x_pad = int(round(bbox_w * padding_fraction))
    top = max(0, y0 - y_pad)
    bottom = min(h, y1 + y_pad + 1)
    left = max(0, x0 - x_pad)
    right = min(w, x1 + x_pad + 1)
    cropped = img[top:bottom, left:right]
    performed = cropped.shape[0] < h or cropped.shape[1] < w
    if performed:
        return cropped, True, f"automatic white-border crop; bbox=({left},{top})-({right},{bottom})"
    return img, False, "automatic crop found content across the full image; image retained unchanged"


def plot_schematic(ax: plt.Axes, schematic_path: Path | None) -> tuple[bool, str]:
    ax.set_title("Dual-axis working model", pad=8, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.axis("off")
    panel_label(ax, "C")
    if schematic_path is None or not schematic_path.exists():
        ax.text(0.5, 0.5, "schematic not found", transform=ax.transAxes, ha="center", va="center", color=TEXT_GREY)
        return False, "schematic not found"
    img = mpimg.imread(schematic_path)
    img, crop_performed, crop_note = crop_white_border(img, padding_fraction=0.03)
    image_ax = ax.inset_axes([0.0, 0.17, 1.0, 0.66])
    image_ax.imshow(img)
    image_ax.set_anchor("C")
    image_ax.axis("off")
    return crop_performed, crop_note


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
    cols = ["term_clean", "display_label", "neg_log10_fdr", "overlap_count_value", "overlap_ratio_value", "fdr_value"]
    cols = [col for col in cols if col in terms.columns]
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
    white_crop_performed: bool,
    crop_note: str,
    unmatched_term_labels: list[str],
    nfe_deduped: int,
    thrb_deduped: int,
    independent_x_axis: bool,
    warnings: list[str],
) -> None:
    notes = f"""
    {STEM}
    ==============================

    Version note:
    v4 is a layout-level reconstruction, not an upstream analysis update. It keeps the same existing
    functional interpretation, program convergence, and schematic inputs.

    Caption title:
    {CAPTION_TITLE}

    Input enrichment files:
    - {rel(INPUTS["go"])}
    - {rel(INPUTS["kegg"])}
    - {rel(INPUTS["nfe_top"])}
    - {rel(INPUTS["thrb_top"])}

    Enrichment column mapping:
    - {enrich_cols}

    Panel A selected functional programs:
    - NFE2L2 terms selected for the lesion-associated program:
    {textwrap.indent(terms_table(nfe_terms).strip(), "    ")}

    - THRB terms selected for the homeostatic-supportive program:
    {textwrap.indent(terms_table(thrb_terms).strip(), "    ")}
    - NFE2L2 terms removed during selected-term curation/de-duplication: {nfe_deduped}
    - THRB terms removed during selected-term curation/de-duplication: {thrb_deduped}
    - THRB term de-duplication successful: {thrb_deduped > 0}
    - Panel A uses independent x-axis ranges for NFE2L2 and THRB facets: {independent_x_axis}

    Panel B convergence scorecard:
    - Table and columns: {conv_cols}
    - Metrics: regulon targets, overlap genes, overlap ratio, and short interpretation class.
    - Values:
    {textwrap.indent(compact_table(conv).strip(), "    ")}

    Panel C schematic:
    - Schematic PNG path: {rel(schematic_path) if schematic_path else "not found"}
    - Figure5_model_clean.png found: {clean_schematic_found}
    - Automatic white-border crop performed: {white_crop_performed}
    - Crop note: {crop_note}
    - No GPT image regeneration was performed.
    - If schematic-internal text remains too small after journal sizing, regenerate a cleaner model with larger text and fewer internal labels in a later design pass.

    Term-label mapping:
    - Unmatched selected term labels: {", ".join(unmatched_term_labels) if unmatched_term_labels else "None."}

    Required interpretation notes:
    Figure 5 v4 uses existing functional interpretation and schematic image outputs only; no upstream analysis was rerun.
    Drug repositioning clues were not included in main Figure 5 and remain supplementary.
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
    unmatched_term_labels: list[str] = []

    enrich, enrich_cols = load_enrichment(warnings)
    nfe_terms = select_program_terms(enrich, enrich_cols, "NFE2L2", NFE2L2_TERM_KEYWORDS, warnings)
    thrb_terms = select_program_terms(enrich, enrich_cols, "THRB", THRB_TERM_KEYWORDS, warnings)
    nfe_terms = add_display_labels(nfe_terms, "NFE2L2", unmatched_term_labels)
    thrb_terms = add_display_labels(thrb_terms, "THRB", unmatched_term_labels)
    nfe_terms, nfe_deduped = curate_program_terms(nfe_terms, "NFE2L2", max_terms=5)
    thrb_terms, thrb_deduped = curate_program_terms(thrb_terms, "THRB", max_terms=5)
    conv, conv_cols = load_convergence(warnings)
    schematic_path, clean_schematic_found = find_schematic(warnings)
    white_crop_performed = False
    crop_note = "not run"
    independent_x_axis = True

    fig = plt.figure(figsize=(14.5, 8.5), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.95], width_ratios=[1.08, 1.0], hspace=0.50, wspace=0.18)
    ax_a_frame = fig.add_subplot(gs[0, 0])
    ax_a_frame.patch.set_alpha(0)
    ax_a_frame.set_zorder(5)
    sub_a = gs[0, 0].subgridspec(1, 2, wspace=0.58)
    ax_a_nfe = fig.add_subplot(sub_a[0, 0])
    ax_a_thrb = fig.add_subplot(sub_a[0, 1])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[:, 1])

    independent_x_axis = plot_selected_programs(ax_a_frame, ax_a_nfe, ax_a_thrb, nfe_terms, thrb_terms)
    plot_convergence(ax_b, conv)
    white_crop_performed, crop_note = plot_schematic(ax_c, schematic_path)

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
        white_crop_performed,
        crop_note,
        unmatched_term_labels,
        nfe_deduped,
        thrb_deduped,
        independent_x_axis,
        warnings,
    )

    for key in ["png", "pdf", "svg", "notes"]:
        print(f"Wrote {rel(outputs[key])}")
    print(f"Figure5_model_clean.png found: {clean_schematic_found}")
    if schematic_path:
        print(f"Schematic used: {rel(schematic_path)}")
    print(f"White-border crop performed: {white_crop_performed}")
    print(f"Crop note: {crop_note}")
    print(f"NFE2L2 selected-term rows removed: {nfe_deduped}")
    print(f"THRB selected-term rows removed: {thrb_deduped}")
    print(f"Panel A independent x-axis ranges: {independent_x_axis}")
    if unmatched_term_labels:
        print("Unmatched term labels:")
        for item in unmatched_term_labels:
            print(f"- {item}")
    else:
        print("Unmatched term labels: none")
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"- {warning}")
    else:
        print("Warnings: none")


if __name__ == "__main__":
    main()
