from __future__ import annotations

from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PREFERRED_NAME_SCORES = [
    ("disease", 220),
    ("diagnosis", 220),
    ("condition", 210),
    ("status", 200),
    ("group", 200),
    ("cohort", 170),
    ("treat", 160),
    ("response", 150),
    ("subtype", 140),
    ("sample", 120),
    ("patient", 110),
    ("donor", 100),
    ("cell_type", 90),
    ("celltype", 90),
    ("cluster", 80),
    ("batch", 40),
]

EXCLUDED_NAME_SCORES = [
    ("barcode", 400),
    ("join_key", 400),
    ("cell_id", 250),
    ("obs_name", 250),
    ("index", 150),
]


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def default_h5ad_path(script_file: str | Path) -> Path:
    root = project_root_from_file(script_file)
    return root / "final_exports" / "astrocyte_pilot_rna_with_pyscenic_auc.h5ad"


def default_output_root(script_file: str | Path) -> Path:
    root = project_root_from_file(script_file)
    return root / "analysis_outputs"


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_adata(h5ad_path: str | Path) -> ad.AnnData:
    adata = ad.read_h5ad(h5ad_path)
    adata.obs_names = adata.obs_names.astype(str)
    return adata


def extract_auc_dataframe(
    adata: ad.AnnData,
    obsm_key: str = "X_pyscenic_auc",
    regulon_name_key: str = "pyscenic_regulon_names",
) -> pd.DataFrame:
    if obsm_key not in adata.obsm:
        raise KeyError(f"Missing obsm key: {obsm_key}")
    if regulon_name_key not in adata.uns:
        raise KeyError(f"Missing uns key: {regulon_name_key}")

    matrix = adata.obsm[obsm_key]
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix)

    regulon_names = list(adata.uns[regulon_name_key])
    if matrix.ndim != 2:
        raise ValueError(f"AUC matrix must be 2D, got shape={matrix.shape}")
    if matrix.shape[1] != len(regulon_names):
        raise ValueError(
            "Regulon name count does not match AUC matrix width: "
            f"{matrix.shape[1]} vs {len(regulon_names)}"
        )

    return pd.DataFrame(matrix, index=adata.obs_names, columns=regulon_names)


def find_umap_key(obsm_keys: Iterable[str]) -> str | None:
    keys = list(obsm_keys)
    exact_priority = ["X_umap", "umap", "X_UMAP", "UMAP"]
    for key in exact_priority:
        if key in keys:
            return key
    for key in keys:
        if "umap" in str(key).lower():
            return key
    return None


def _series_value_counts(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("NA").value_counts(dropna=False)


def _name_score(column_name: str) -> tuple[int, list[str]]:
    name = column_name.lower()
    score = 0
    reasons: list[str] = []
    for pattern, value in PREFERRED_NAME_SCORES:
        if pattern in name:
            score += value
            reasons.append(f"name:+{pattern}")
    for pattern, value in EXCLUDED_NAME_SCORES:
        if pattern in name:
            score -= value
            reasons.append(f"name:-{pattern}")
    return score, reasons


def score_group_column(column_name: str, series: pd.Series) -> tuple[float, bool, list[str]]:
    total_n = int(series.shape[0])
    non_na = series.dropna()
    n_missing = int(series.isna().sum())
    n_non_na = int(non_na.shape[0])
    n_unique = int(non_na.nunique())
    counts = _series_value_counts(non_na) if n_non_na else pd.Series(dtype="int64")
    min_count = int(counts.min()) if not counts.empty else 0
    max_count = int(counts.max()) if not counts.empty else 0
    balance_ratio = (min_count / max_count) if max_count else 0.0

    score, reasons = _name_score(column_name)
    dtype = series.dtype

    if pd.api.types.is_bool_dtype(dtype):
        score += 45
        reasons.append("dtype:bool")
    elif pd.api.types.is_numeric_dtype(dtype):
        if n_unique <= 12:
            score += 20
            reasons.append("dtype:low_card_numeric")
        else:
            score -= 60
            reasons.append("dtype:high_card_numeric")
    else:
        score += 35
        reasons.append("dtype:categorical_like")

    if n_unique < 2:
        score -= 1000
        reasons.append("levels:single_value")
    elif n_unique == n_non_na:
        score -= 1000
        reasons.append("levels:all_unique")
    elif n_unique <= 2:
        score += 80
        reasons.append("levels:2")
    elif n_unique <= 4:
        score += 65
        reasons.append("levels:3-4")
    elif n_unique <= 8:
        score += 45
        reasons.append("levels:5-8")
    elif n_unique <= 12:
        score += 25
        reasons.append("levels:9-12")
    elif n_unique <= 20:
        score += 10
        reasons.append("levels:13-20")
    else:
        score -= 200
        reasons.append("levels:>20")

    if min_count >= 100:
        score += 30
        reasons.append("counts:min>=100")
    elif min_count >= 30:
        score += 20
        reasons.append("counts:min>=30")
    elif min_count >= 10:
        score += 10
        reasons.append("counts:min>=10")
    elif min_count > 0:
        score -= 10
        reasons.append("counts:min<10")

    score += 30 * balance_ratio
    reasons.append(f"balance:{balance_ratio:.3f}")

    if n_missing == 0:
        score += 8
        reasons.append("missing:0")
    else:
        score -= min(30, n_missing / max(total_n, 1) * 100)
        reasons.append(f"missing:{n_missing}")

    eligible = (
        n_unique >= 2
        and n_unique < max(20, total_n)
        and n_unique != n_non_na
        and "name:-barcode" not in reasons
        and "name:-join_key" not in reasons
    )
    return float(score), eligible, reasons


def preview_values(series: pd.Series, top_n: int = 6) -> str:
    counts = _series_value_counts(series)
    if counts.empty:
        return ""
    items = []
    for value, count in counts.head(top_n).items():
        items.append(f"{value}:{int(count)}")
    return " | ".join(items)


def summarize_obs_columns(obs: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for column in obs.columns:
        series = obs[column]
        score, eligible, reasons = score_group_column(column, series)
        non_na = series.dropna()
        counts = _series_value_counts(non_na) if not non_na.empty else pd.Series(dtype="int64")
        min_count = int(counts.min()) if not counts.empty else 0
        max_count = int(counts.max()) if not counts.empty else 0
        rows.append(
            {
                "column": column,
                "dtype": str(series.dtype),
                "n_missing": int(series.isna().sum()),
                "n_unique_non_na": int(non_na.nunique()),
                "n_unique_including_na": int(series.nunique(dropna=False)),
                "min_group_size": min_count,
                "max_group_size": max_count,
                "value_preview": preview_values(series),
                "group_candidate_score": round(score, 3),
                "is_group_candidate": bool(eligible),
                "score_reasons": "; ".join(reasons),
            }
        )

    summary = pd.DataFrame(rows)
    summary = summary.sort_values(
        by=["group_candidate_score", "n_unique_non_na", "column"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    return summary


def choose_group_column(
    obs: pd.DataFrame,
    preferred_column: str | None = None,
) -> tuple[str, pd.DataFrame]:
    summary = summarize_obs_columns(obs)
    if preferred_column:
        if preferred_column not in obs.columns:
            raise KeyError(f"Requested group column not found: {preferred_column}")
        return preferred_column, summary

    candidates = summary[summary["is_group_candidate"]].copy()
    if candidates.empty:
        fallback = summary.sort_values(
            by=["n_unique_non_na", "group_candidate_score"],
            ascending=[True, False],
        ).iloc[0]
        return str(fallback["column"]), summary

    chosen = candidates.iloc[0]
    return str(chosen["column"]), summary


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    output = np.full(values.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(values)
    if not finite_mask.any():
        return output

    finite_values = values[finite_mask]
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    n = ranked.size
    adjusted = np.empty(n, dtype=float)
    previous = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        current = ranked[idx] * n / rank
        previous = min(previous, current)
        adjusted[idx] = previous

    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    output[finite_mask] = np.clip(restored, 0.0, 1.0)
    return output


def row_zscore(df: pd.DataFrame) -> pd.DataFrame:
    values = df.to_numpy(dtype=float)
    row_mean = values.mean(axis=1, keepdims=True)
    row_std = values.std(axis=1, ddof=0, keepdims=True)
    row_std[row_std == 0] = 1.0
    z = (values - row_mean) / row_std
    return pd.DataFrame(z, index=df.index, columns=df.columns)


def plot_heatmap(
    data: pd.DataFrame,
    output_png: str | Path,
    title: str,
    cmap: str = "viridis",
    colorbar_label: str = "value",
    figsize: tuple[float, float] | None = None,
) -> None:
    if data.empty:
        raise ValueError("Heatmap input is empty.")

    output_png = Path(output_png)
    if figsize is None:
        width = max(6.0, 1.2 + 1.0 * data.shape[1])
        height = max(6.0, 0.22 * data.shape[0] + 1.5)
        figsize = (width, height)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(data.to_numpy(dtype=float), aspect="auto", interpolation="nearest", cmap=cmap)
    ax.set_title(title)
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_xticklabels(list(data.columns), rotation=45, ha="right")
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_yticklabels(list(data.index))
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label(colorbar_label)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
