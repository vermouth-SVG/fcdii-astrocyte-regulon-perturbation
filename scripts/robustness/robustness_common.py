#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import math
from pathlib import Path
from typing import Iterable

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[2]
INPUT_H5AD = ROOT / "final_exports" / "astrocyte_pilot_rna_with_pyscenic_auc.h5ad"
OUT_DIR = ROOT / "robustness_validation"
GROUP_COMPARE = ROOT / "analysis_outputs" / "group_compare" / "regulon_group_statistics.csv"
CELLORACLE_CANDIDATES = ROOT / "analysis_outputs" / "celloracle_candidates" / "celloracle_candidate_tf_metrics.csv"
DISCOVERY_REGULONS = ROOT / "output" / "regulons.csv"

CORE_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
OPTIONAL_TFS = ["HMGA1", "SATB2", "RARB"]
ALL_TFS = CORE_TFS + OPTIONAL_TFS
LESION_ALIASES = {"lesion", "tle", "disease", "case", "epilepsy"}
CONTROL_ALIASES = {"internal_control", "control", "normal", "non-lesion", "nonlesion", "healthy"}
EPS = 1e-9


def ensure_out_dir() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUT_DIR


def decode_value(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return value.decode("utf-8")
    return value


def decode_array(values) -> list:
    arr = np.asarray(values)
    return [decode_value(x) for x in arr.tolist()]


def _attr_to_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        return [str(decode_value(x)) for x in value.tolist()]
    return [str(decode_value(value))]


def read_h5ad_dataframe(h5: h5py.File, key: str) -> pd.DataFrame:
    group = h5[key]
    index_key = str(decode_value(group.attrs.get("_index", "_index")))
    column_order = _attr_to_list(group.attrs.get("column-order", []))
    index_values = decode_array(group[index_key][...])
    data = {}
    for col in column_order:
        if col == index_key:
            continue
        obj = group[col]
        if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
            categories = decode_array(obj["categories"][...])
            codes = np.asarray(obj["codes"][...])
            values = []
            for code in codes:
                if int(code) < 0:
                    values.append(pd.NA)
                else:
                    values.append(categories[int(code)])
            data[col] = values
        elif isinstance(obj, h5py.Dataset):
            data[col] = decode_array(obj[...])
        else:
            data[col] = [pd.NA] * len(index_values)
    return pd.DataFrame(data, index=pd.Index(index_values, name=index_key))


def read_var_names(h5: h5py.File, key: str = "var") -> pd.Index:
    group = h5[key]
    index_key = str(decode_value(group.attrs.get("_index", "_index")))
    return pd.Index(decode_array(group[index_key][...]), name=index_key)


def read_csr_matrix(h5: h5py.File, key: str) -> sparse.csr_matrix:
    group = h5[key]
    shape = tuple(int(x) for x in group.attrs["shape"])
    data = np.asarray(group["data"][...])
    indices = np.asarray(group["indices"][...])
    indptr = np.asarray(group["indptr"][...])
    return sparse.csr_matrix((data, indices, indptr), shape=shape)


def load_discovery_h5ad_light(path: Path = INPUT_H5AD) -> dict:
    with h5py.File(path, "r") as h5:
        obs = read_h5ad_dataframe(h5, "obs")
        var = read_h5ad_dataframe(h5, "var")
        var_names = read_var_names(h5, "var")
        if "raw" in h5 and "X" in h5["raw"]:
            matrix = read_csr_matrix(h5, "raw/X")
            matrix_source = "raw/X"
            if "var" in h5["raw"]:
                var_names = read_var_names(h5, "raw/var")
        else:
            matrix = read_csr_matrix(h5, "X")
            matrix_source = "X"
        auc = None
        regulon_names = []
        if "obsm" in h5 and "X_pyscenic_auc" in h5["obsm"]:
            auc = np.asarray(h5["obsm"]["X_pyscenic_auc"][...], dtype=float)
        if "uns" in h5 and "pyscenic_regulon_names" in h5["uns"]:
            regulon_names = [str(x) for x in decode_array(h5["uns"]["pyscenic_regulon_names"][...])]
        obsm_keys = list(h5["obsm"].keys()) if "obsm" in h5 else []
        uns_keys = list(h5["uns"].keys()) if "uns" in h5 else []
        layer_keys = list(h5["layers"].keys()) if "layers" in h5 else []
    return {
        "obs": obs,
        "var": var,
        "var_names": var_names,
        "matrix": matrix,
        "matrix_source": matrix_source,
        "auc": auc,
        "regulon_names": regulon_names,
        "obsm_keys": obsm_keys,
        "uns_keys": uns_keys,
        "layer_keys": layer_keys,
        "h5ad_path": path,
    }


def detect_group_column(obs: pd.DataFrame) -> tuple[str, str, str]:
    best_col = None
    best_score = -10**9
    best_labels = ("", "")
    for col in obs.columns:
        values = obs[col].astype("string").dropna().astype(str)
        unique = values.unique().tolist()
        if len(unique) < 2 or len(unique) > 20:
            continue
        lower_map = {v.lower(): v for v in unique}
        lesion = next((lower_map[x] for x in LESION_ALIASES if x in lower_map), None)
        control = next((lower_map[x] for x in CONTROL_ALIASES if x in lower_map), None)
        score = 0
        lname = col.lower()
        if "group" in lname:
            score += 100
        if "disease" in lname or "condition" in lname:
            score += 80
        if lesion and control:
            score += 500
        if len(unique) == 2:
            score += 50
        if score > best_score:
            best_score = score
            best_col = col
            best_labels = (lesion or unique[0], control or unique[1])
    if best_col is None:
        raise ValueError("无法自动识别分组字段：obs 中没有合适的二分类 group/disease/condition 字段。")
    return best_col, best_labels[0], best_labels[1]


def detect_sample_fields(obs: pd.DataFrame, group_col: str) -> dict:
    donor_candidates = ["donor_id", "donor", "patient_id", "patient", "subject_id", "subject"]
    sample_candidates = ["sample_id", "sample", "sample_prefix", "clinical_sample_label", "pilot_batch", "batch"]
    donor_col = next((c for c in donor_candidates if c in obs.columns), None)
    sample_col = next((c for c in sample_candidates if c in obs.columns), None)

    donor_n = int(obs[donor_col].astype("string").nunique()) if donor_col else 0
    sample_n = int(obs[sample_col].astype("string").nunique()) if sample_col else 0

    if sample_col and sample_n >= 4:
        unit_col = sample_col
        unit_reason = (
            f"检测到 {donor_col}={donor_n} 个 donor，同时 {sample_col}={sample_n} 个样本；"
            "本项目为 4-sample astrocyte pilot，因此 leave-one-out/pseudobulk 默认使用样本级单位。"
        )
    elif donor_col:
        unit_col = donor_col
        unit_reason = f"使用 donor 字段 {donor_col} 作为分析单位。"
    elif sample_col:
        unit_col = sample_col
        unit_reason = f"未检测到 donor 字段，使用 sample 字段 {sample_col} 作为 donor 替代。"
    else:
        raise ValueError("无法自动识别 donor/sample 字段。")

    distribution = pd.crosstab(obs[unit_col].astype(str), obs[group_col].astype(str))
    return {
        "donor_col": donor_col,
        "sample_col": sample_col,
        "analysis_unit_col": unit_col,
        "donor_n": donor_n,
        "sample_n": sample_n,
        "unit_reason": unit_reason,
        "unit_group_distribution": distribution,
    }


def make_gene_index(var_names: Iterable[str]) -> dict[str, int]:
    index = {}
    for i, name in enumerate(var_names):
        text = str(name)
        if text not in index:
            index[text] = i
        upper = text.upper()
        if upper not in index:
            index[upper] = i
    return index


def get_gene_vector(matrix: sparse.csr_matrix, gene_index: dict[str, int], gene: str) -> np.ndarray | None:
    idx = gene_index.get(gene) if gene in gene_index else gene_index.get(gene.upper())
    if idx is None:
        return None
    return np.asarray(matrix[:, idx].toarray()).ravel()


def auc_dataframe(auc: np.ndarray | None, regulon_names: list[str], obs_index: pd.Index) -> pd.DataFrame | None:
    if auc is None or not regulon_names:
        return None
    return pd.DataFrame(auc, index=obs_index.astype(str), columns=regulon_names)


def resolve_regulon_name(tf: str, regulon_names: list[str]) -> str | None:
    candidates = [f"{tf}(+)", tf, tf.upper(), tf.lower()]
    upper_map = {str(name).upper(): str(name) for name in regulon_names}
    for candidate in candidates:
        if candidate.upper() in upper_map:
            return upper_map[candidate.upper()]
    for name in regulon_names:
        if str(name).upper().startswith(tf.upper() + "("):
            return str(name)
    return None


def benjamini_hochberg(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    finite = np.isfinite(arr)
    if not finite.any():
        return out
    p = arr[finite]
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    adj = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        running = min(running, ranked[i] * n / (i + 1))
        adj[i] = running
    restored = np.empty(n, dtype=float)
    restored[order] = adj
    out[finite] = np.clip(restored, 0, 1)
    return out


def differential_vector(values: np.ndarray, groups: pd.Series, lesion_label: str, control_label: str) -> dict:
    values = np.asarray(values, dtype=float)
    group_arr = groups.astype(str).to_numpy()
    lesion_mask = group_arr == str(lesion_label)
    control_mask = group_arr == str(control_label)
    lesion_values = values[lesion_mask]
    control_values = values[control_mask]
    if lesion_values.size == 0 or control_values.size == 0:
        return {
            "n_lesion": int(lesion_values.size),
            "n_internal_control": int(control_values.size),
            "mean_lesion": np.nan,
            "mean_internal_control": np.nan,
            "median_lesion": np.nan,
            "median_internal_control": np.nan,
            "positive_fraction_lesion": np.nan,
            "positive_fraction_internal_control": np.nan,
            "log2fc_lesion_vs_internal_control": np.nan,
            "mean_diff_lesion_minus_internal_control": np.nan,
            "p_value": np.nan,
            "statistic": np.nan,
        }
    mean_lesion = float(np.mean(lesion_values))
    mean_control = float(np.mean(control_values))
    try:
        stat, p = stats.mannwhitneyu(lesion_values, control_values, alternative="two-sided")
    except ValueError:
        stat, p = np.nan, np.nan
    return {
        "n_lesion": int(lesion_values.size),
        "n_internal_control": int(control_values.size),
        "mean_lesion": mean_lesion,
        "mean_internal_control": mean_control,
        "median_lesion": float(np.median(lesion_values)),
        "median_internal_control": float(np.median(control_values)),
        "positive_fraction_lesion": float(np.mean(lesion_values > 0)),
        "positive_fraction_internal_control": float(np.mean(control_values > 0)),
        "log2fc_lesion_vs_internal_control": float(np.log2((mean_lesion + EPS) / (mean_control + EPS))),
        "mean_diff_lesion_minus_internal_control": float(mean_lesion - mean_control),
        "p_value": float(p) if np.isfinite(p) else np.nan,
        "statistic": float(stat) if np.isfinite(stat) else np.nan,
    }


def sign_consistent(value: float, baseline: float) -> bool:
    if not np.isfinite(value) or not np.isfinite(baseline):
        return False
    if abs(baseline) < EPS:
        return abs(value) < EPS
    return math.copysign(1, value) == math.copysign(1, baseline)


def direction_from_value(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value > 0:
        return "lesion"
    if value < 0:
        return "internal_control"
    return "tie"


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float)
    mean = np.nanmean(values, axis=1, keepdims=True)
    std = np.nanstd(values, axis=1, keepdims=True)
    std[std < EPS] = 1.0
    z = (values - mean) / std
    return pd.DataFrame(z, index=matrix.index, columns=matrix.columns)


def save_heatmap(matrix: pd.DataFrame, path: Path, title: str, cmap: str = "coolwarm", center_zero: bool = True) -> None:
    if matrix.empty:
        return
    fig_w = max(7, min(18, 0.38 * matrix.shape[1] + 2))
    fig_h = max(4, min(14, 0.35 * matrix.shape[0] + 2))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    values = matrix.to_numpy(dtype=float)
    if center_zero:
        vmax = np.nanmax(np.abs(values)) if values.size else 1
        vmax = max(vmax, 1e-3)
        im = ax.imshow(values, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    else:
        im = ax.imshow(values, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns.astype(str), rotation=90, fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.astype(str), fontsize=9)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _find_multi_column(columns, target: str):
    target = target.lower()
    for col in columns:
        if any(str(level).lower() == target for level in col):
            return col
    raise KeyError(target)


def parse_target_genes(raw_value: object) -> list[str]:
    if raw_value is None or (isinstance(raw_value, float) and pd.isna(raw_value)):
        return []
    text = str(raw_value).strip()
    if not text:
        return []
    parsed = ast.literal_eval(text)
    out = []
    for item in parsed if isinstance(parsed, list) else [parsed]:
        if isinstance(item, (tuple, list)) and item:
            gene = str(item[0]).strip()
        else:
            gene = str(item).strip()
        if gene and gene not in out:
            out.append(gene)
    return out


def discovery_targets_for_tfs(tfs: list[str], regulons_csv: Path = DISCOVERY_REGULONS) -> dict[str, list[str]]:
    if not regulons_csv.exists():
        return {}
    df = pd.read_csv(regulons_csv, header=[0, 1, 2])
    tf_col = _find_multi_column(df.columns, "TF")
    auc_col = _find_multi_column(df.columns, "AUC")
    nes_col = _find_multi_column(df.columns, "NES")
    ann_col = _find_multi_column(df.columns, "Annotation")
    ctx_col = _find_multi_column(df.columns, "Context")
    target_col = _find_multi_column(df.columns, "TargetGenes")
    flat = pd.DataFrame(
        {
            "TF": df[tf_col].astype(str),
            "AUC": pd.to_numeric(df[auc_col], errors="coerce"),
            "NES": pd.to_numeric(df[nes_col], errors="coerce"),
            "Annotation": df[ann_col].astype(str),
            "Context": df[ctx_col].astype(str),
            "TargetGenes": df[target_col].astype(str),
        }
    )
    result = {}
    for tf in tfs:
        sub = flat[flat["TF"].str.upper() == tf.upper()].copy()
        if sub.empty:
            continue
        sub["direct"] = sub["Annotation"].str.lower().str.contains("direct")
        sub["activating"] = sub["Context"].str.lower().str.contains("activating")
        sub = sub.sort_values(["direct", "activating", "NES", "AUC"], ascending=[False, False, False, False])
        result[tf] = parse_target_genes(sub.iloc[0]["TargetGenes"])
    return result

