#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path
from urllib.request import urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[2]
ROBUSTNESS_SCRIPTS = ROOT / "scripts" / "robustness"
if str(ROBUSTNESS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ROBUSTNESS_SCRIPTS))

from robustness_common import (  # noqa: E402
    CELLORACLE_CANDIDATES,
    DISCOVERY_REGULONS,
    INPUT_H5AD,
    OUT_DIR as ROBUSTNESS_DIR,
    benjamini_hochberg,
    detect_group_column,
    detect_sample_fields,
    discovery_targets_for_tfs,
    load_discovery_h5ad_light,
)


FUNCTIONAL_DIR = ROOT / "functional_interpretation"
RESOURCE_DIR = FUNCTIONAL_DIR / "resources"
GENESET_DIR = RESOURCE_DIR / "gene_sets"

FOCUS_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
PRIMARY_TFS = ["NFE2L2", "THRB"]
SECOND_TIER_TFS = ["BHLHE40"]
RESERVED_TFS = ["SOX2"]
TF_MAIN_DIRECTION = {
    "NFE2L2": "lesion",
    "BHLHE40": "lesion",
    "SOX2": "lesion",
    "THRB": "internal_control",
}
TF_TIER = {
    "NFE2L2": "primary_lesion_axis",
    "THRB": "primary_internal_control_axis",
    "BHLHE40": "second_tier",
    "SOX2": "reserved_lightweight",
}
ENRICHR_LIBRARIES = {
    "GO_Biological_Process_2023": "go",
    "KEGG_2021_Human": "kegg",
}
DEG_LOG2FC_THRESHOLD = 0.25
DEG_POSITIVE_FRACTION_THRESHOLD = 0.05
EPS = 1e-9


def ensure_functional_dir() -> Path:
    FUNCTIONAL_DIR.mkdir(parents=True, exist_ok=True)
    RESOURCE_DIR.mkdir(parents=True, exist_ok=True)
    GENESET_DIR.mkdir(parents=True, exist_ok=True)
    return FUNCTIONAL_DIR


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def read_csv_optional(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def load_robustness_tables() -> dict[str, pd.DataFrame]:
    return {
        "integrated": read_csv_optional(ROBUSTNESS_DIR / "05_robustness_validation_integrated_table.csv"),
        "loo_summary": read_csv_optional(ROBUSTNESS_DIR / "02_leave_one_out_tf_summary.csv"),
        "pseudobulk_tf": read_csv_optional(ROBUSTNESS_DIR / "03_pseudobulk_candidate_tf_table.csv"),
        "pseudobulk_regulon": read_csv_optional(ROBUSTNESS_DIR / "03_pseudobulk_regulon_table.csv"),
        "shortlist_frequency": read_csv_optional(ROBUSTNESS_DIR / "04_shortlist_membership_frequency.csv"),
    }


def compute_light_deg(data: dict, group_col: str, lesion_label: str, control_label: str) -> pd.DataFrame:
    obs = data["obs"]
    matrix = data["matrix"]
    var_names = pd.Index(data["var_names"].astype(str))
    group_values = obs[group_col].astype(str).to_numpy()
    lesion_mask = group_values == str(lesion_label)
    control_mask = group_values == str(control_label)
    if lesion_mask.sum() == 0 or control_mask.sum() == 0:
        raise ValueError("无法计算 DEG：lesion 或 internal_control 细胞数为 0。")

    lesion_matrix = matrix[lesion_mask, :]
    control_matrix = matrix[control_mask, :]
    mean_lesion = np.asarray(lesion_matrix.mean(axis=0)).ravel()
    mean_control = np.asarray(control_matrix.mean(axis=0)).ravel()
    if sparse.issparse(lesion_matrix):
        pos_lesion = np.asarray((lesion_matrix > 0).mean(axis=0)).ravel()
        pos_control = np.asarray((control_matrix > 0).mean(axis=0)).ravel()
    else:
        pos_lesion = np.mean(lesion_matrix > 0, axis=0)
        pos_control = np.mean(control_matrix > 0, axis=0)
    log2fc = np.log2((mean_lesion + EPS) / (mean_control + EPS))
    direction = np.where(log2fc > 0, "lesion", np.where(log2fc < 0, "internal_control", "tie"))
    return pd.DataFrame(
        {
            "gene": var_names,
            "mean_lesion": mean_lesion,
            "mean_internal_control": mean_control,
            "positive_fraction_lesion": pos_lesion,
            "positive_fraction_internal_control": pos_control,
            "log2fc_lesion_vs_internal_control": log2fc,
            "abs_log2fc": np.abs(log2fc),
            "direction": direction,
            "positive_fraction_higher_group": np.where(log2fc >= 0, pos_lesion, pos_control),
        }
    )


def build_directional_deg_sets(deg: pd.DataFrame) -> dict[str, set[str]]:
    base = deg[
        (deg["abs_log2fc"] >= DEG_LOG2FC_THRESHOLD)
        & (deg["positive_fraction_higher_group"] >= DEG_POSITIVE_FRACTION_THRESHOLD)
    ].copy()
    return {
        "lesion": set(base.loc[base["direction"] == "lesion", "gene"].astype(str)),
        "internal_control": set(base.loc[base["direction"] == "internal_control", "gene"].astype(str)),
    }


def download_enrichr_library(library_name: str, path: Path) -> tuple[bool, str]:
    if path.exists() and path.stat().st_size > 0:
        return True, "使用本地缓存"
    url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library_name}"
    try:
        with urlopen(url, timeout=60) as response:
            text = response.read().decode("utf-8", errors="replace")
        if not text.strip():
            return False, "下载结果为空"
        path.write_text(text, encoding="utf-8", newline="\n")
        return True, "已从 Enrichr 下载并缓存"
    except Exception as exc:  # noqa: BLE001
        return False, f"下载失败: {type(exc).__name__}: {exc}"


def parse_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = {g.strip().upper() for g in parts[2:] if g.strip()}
            if genes:
                gene_sets[term] = genes
    return gene_sets


def run_ora(query_genes: set[str], universe_genes: set[str], gene_sets: dict[str, set[str]], source: str, label: str) -> pd.DataFrame:
    q = {g.upper() for g in query_genes if str(g).strip()}
    universe = {g.upper() for g in universe_genes if str(g).strip()}
    q = q & universe
    rows = []
    M = len(universe)
    N = len(q)
    if M == 0 or N == 0:
        return pd.DataFrame()
    for term, genes in gene_sets.items():
        term_genes = genes & universe
        n = len(term_genes)
        overlap = sorted(q & term_genes)
        k = len(overlap)
        if k < 1:
            continue
        p_value = stats.hypergeom.sf(k - 1, M, n, N)
        fold = (k / N) / (n / M) if n and N else np.nan
        rows.append(
            {
                "query_label": label,
                "source": source,
                "term": term,
                "query_size": N,
                "term_size": n,
                "overlap_count": k,
                "overlap_ratio": k / N,
                "fold_enrichment": fold,
                "p_value": p_value,
                "overlap_genes": ";".join(overlap[:80]),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["fdr_bh"] = benjamini_hochberg(df["p_value"].to_numpy(dtype=float))
    return df.sort_values(["fdr_bh", "p_value", "overlap_count"], ascending=[True, True, False])


THEME_KEYWORDS = {
    "transcriptional_regulatory": ["transcription", "rna polymerase", "dna-templated", "gene expression", "chromatin"],
    "oxidative_stress": ["oxidative", "reactive oxygen", "glutathione", "redox", "response to stress", "detox"],
    "hypoxia_hif_mtor": ["hypoxia", "hif", "mtor", "autophagy", "pi3k", "akt"],
    "inflammatory_reactive": ["inflammatory", "cytokine", "interleukin", "immune", "nf-kappa", "tnf", "response to wounding", "gliosis"],
    "metabolic_supportive": ["metabolic", "metabolism", "mitochond", "lipid", "fatty acid", "transport", "homeostasis", "ion"],
    "thyroid_hormone_neuro": ["thyroid", "hormone", "synap", "axon", "neuron", "neuro", "glutamate", "calcium"],
    "differentiation_development": ["differentiation", "development", "morphogenesis", "stem cell", "glial"],
    "cell_cycle_growth": ["cell cycle", "proliferation", "growth", "dna replication", "mitotic"],
}


def infer_themes_from_terms(terms: list[str]) -> dict[str, int]:
    counts = {theme: 0 for theme in THEME_KEYWORDS}
    for term in terms:
        low = str(term).lower()
        for theme, keywords in THEME_KEYWORDS.items():
            if any(k in low for k in keywords):
                counts[theme] += 1
    return counts


def top_terms_string(df: pd.DataFrame, source: str, n: int = 5) -> str:
    if df.empty or "source" not in df.columns:
        return ""
    sub = df[df["source"] == source].head(n)
    return "; ".join(sub["term"].astype(str).tolist())


def save_dotplot(df: pd.DataFrame, tf: str, path: Path, top_n: int = 12) -> None:
    if df.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis("off")
        ax.text(0.02, 0.5, f"{tf}: 无可绘制富集结果", fontsize=12)
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        return
    plot_df = df.sort_values(["fdr_bh", "p_value"]).head(top_n).copy()
    plot_df["minus_log10_fdr"] = -np.log10(plot_df["fdr_bh"].astype(float).clip(lower=1e-300))
    plot_df = plot_df.iloc[::-1]
    colors = plot_df["source"].map({"go": "#4C78A8", "kegg": "#E45756"}).fillna("#777777")
    fig_h = max(4.5, 0.38 * plot_df.shape[0] + 1.5)
    fig, ax = plt.subplots(figsize=(9.5, fig_h), constrained_layout=True)
    sizes = 35 + 20 * plot_df["overlap_count"].astype(float)
    ax.scatter(plot_df["minus_log10_fdr"], np.arange(plot_df.shape[0]), s=sizes, c=colors, alpha=0.82, edgecolor="black", linewidth=0.3)
    ax.set_yticks(np.arange(plot_df.shape[0]))
    labels = [str(x)[:75] for x in plot_df["term"]]
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("-log10(FDR)")
    ax.set_title(f"{tf} GO/KEGG enrichment")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_matrix_heatmap(matrix: pd.DataFrame, path: Path, title: str, cmap: str = "YlGnBu") -> None:
    if matrix.empty:
        return
    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.45 * matrix.shape[0] + 1.5)), constrained_layout=True)
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns.astype(str), rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.astype(str))
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
