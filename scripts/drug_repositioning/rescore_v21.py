#!/usr/bin/env python3
from __future__ import annotations

import itertools
import math
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from drugrep_common_v2 import DRUGREP_DIR, ROOT


SCRIPT_DIR = Path(__file__).resolve().parent
CONFIG_DIR = SCRIPT_DIR / "config"
OUT_DIR = ROOT / "drug_repositioning"
ARCHIVE_V2_DIR = OUT_DIR / "archive" / "v2_pre_v21_process"


def first_existing(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]


SOURCE_TABLE = first_existing(DRUGREP_DIR / "03_compound_aggregated_table.csv", ARCHIVE_V2_DIR / "03_compound_aggregated_table.csv")
V2_INTEGRATED_TABLE = first_existing(
    DRUGREP_DIR / "05_drug_repositioning_integrated_table.csv",
    ARCHIVE_V2_DIR / "05_drug_repositioning_integrated_table.csv",
)
V2_TOP_TABLE = first_existing(DRUGREP_DIR / "05_top_priority_compounds.csv", ARCHIVE_V2_DIR / "05_top_priority_compounds.csv")
TAG_CONFIG = CONFIG_DIR / "drugrep_v21_compound_tags.csv"
SCORING_CONFIG = CONFIG_DIR / "drugrep_v21_scoring.yaml"
PENALTY_CONFIG = CONFIG_DIR / "drugrep_v21_penalties.yaml"

PRIMARY_AXES = {"NFE2L2", "THRB"}
STRICT_FDR = 0.05
THRB_WEAK_FDR = 0.25
RETAINED_SCORE_MIN = 0.30
SUPPORTIVE_SCORE_MIN = 0.42
PRIMARY_SCORE_MIN = 0.48

TAG_COLUMNS = [
    "compound",
    "tool_compound_like",
    "broad_epigenetic_modulator_like",
    "weak_fdr_manual_review",
    "low_cns_epilepsy_relevance",
    "axis_preference",
    "notes",
]


WEIGHT_PRESETS = [
    {
        "weight_profile": "translational_high",
        "enrichment_support_score": 0.20,
        "cross_signature_score": 0.18,
        "axis_relevance_score": 0.16,
        "mechanism_interpretability_score": 0.08,
        "translational_suitability_score": 0.26,
        "axis_assignment_confidence_score": 0.12,
    },
    {
        "weight_profile": "enrichment_kept_axis_conf",
        "enrichment_support_score": 0.28,
        "cross_signature_score": 0.14,
        "axis_relevance_score": 0.16,
        "mechanism_interpretability_score": 0.06,
        "translational_suitability_score": 0.22,
        "axis_assignment_confidence_score": 0.14,
    },
    {
        "weight_profile": "cross_axis_mechanism_min",
        "enrichment_support_score": 0.24,
        "cross_signature_score": 0.20,
        "axis_relevance_score": 0.18,
        "mechanism_interpretability_score": 0.05,
        "translational_suitability_score": 0.21,
        "axis_assignment_confidence_score": 0.12,
    },
    {
        "weight_profile": "axis_conf_high",
        "enrichment_support_score": 0.22,
        "cross_signature_score": 0.17,
        "axis_relevance_score": 0.15,
        "mechanism_interpretability_score": 0.10,
        "translational_suitability_score": 0.20,
        "axis_assignment_confidence_score": 0.16,
    },
    {
        "weight_profile": "balanced_v21",
        "enrichment_support_score": 0.24,
        "cross_signature_score": 0.16,
        "axis_relevance_score": 0.14,
        "mechanism_interpretability_score": 0.07,
        "translational_suitability_score": 0.23,
        "axis_assignment_confidence_score": 0.16,
    },
]

PENALTY_PRESETS = [
    {
        "penalty_profile": "lenient_low",
        "weak_fdr_penalty": 0.05,
        "single_signature_penalty": 0.05,
        "tool_compound_penalty": 0.08,
        "low_cns_epilepsy_relevance_penalty": 0.08,
    },
    {
        "penalty_profile": "conservative_high",
        "weak_fdr_penalty": 0.12,
        "single_signature_penalty": 0.10,
        "tool_compound_penalty": 0.20,
        "low_cns_epilepsy_relevance_penalty": 0.18,
    },
]


@dataclass(frozen=True)
class V21Config:
    weights: dict[str, float]
    penalties: dict[str, float]
    primary_target_count: int
    primary_min_support: int
    thrb_supportive_min_translational: float
    weight_profile: str
    penalty_profile: str

    @property
    def config_id(self) -> str:
        return (
            f"{self.weight_profile}__{self.penalty_profile}"
            f"__primary{self.primary_target_count}"
            f"__support{self.primary_min_support}"
            f"__thrbtrans{self.thrb_supportive_min_translational:.2f}"
        )


def as_bool(value: Any) -> bool:
    text = str(value).strip().lower()
    return text in {"1", "true", "yes", "y", "是"}


def bool_text(value: bool) -> str:
    return "true" if bool(value) else "false"


def norm_name(value: Any) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def safe_num(series: pd.Series, default: float = 0.0) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(default)


def ensure_dirs() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)


def write_static_configs() -> None:
    if not SCORING_CONFIG.exists():
        lines = [
            "v21_scoring_search:",
            "  note: pre-BBB rescoring only; no Enrichr rerun, no BBB automation, no docking.",
            "  weight_ranges:",
            "    enrichment_support_score: [0.20, 0.28]",
            "    cross_signature_score: [0.14, 0.20]",
            "    axis_relevance_score: [0.14, 0.18]",
            "    mechanism_interpretability_score: [0.05, 0.10]",
            "    translational_suitability_score: [0.18, 0.26]",
            "    axis_assignment_confidence_score: [0.10, 0.16]",
            "  primary_target_count_tested: [4, 5, 6]",
            "  primary_min_support_tested: [2, 3]",
            "  thrb_supportive_min_translational_tested: [0.75, 0.90]",
        ]
        SCORING_CONFIG.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if not PENALTY_CONFIG.exists():
        lines = [
            "v21_penalty_search:",
            "  weak_fdr_penalty: [0.05, 0.12]",
            "  single_signature_penalty: [0.05, 0.10]",
            "  tool_compound_penalty: [0.08, 0.20]",
            "  low_cns_epilepsy_relevance_penalty: [0.08, 0.18]",
            "  notes:",
            "    - tool_compound_penalty is applied to tool_compound_like or broad_epigenetic_modulator_like tags.",
            "    - BBB is not used in v2.1 pre-BBB scoring.",
        ]
        PENALTY_CONFIG.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_source_table() -> pd.DataFrame:
    if not SOURCE_TABLE.exists():
        raise FileNotFoundError(f"缺少 pre-BBB 聚合表: {SOURCE_TABLE}")
    df = pd.read_csv(SOURCE_TABLE)
    required = {
        "compound",
        "associated_axis",
        "support_count",
        "best_adjusted_p",
        "enrichment_support_score",
        "cross_signature_score",
        "axis_relevance_score",
        "mechanism_interpretability_score",
        "translational_suitability_score",
        "axis_assignment_confidence_score",
        "translational_caution_category",
        "indicative_mechanism",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError("03_compound_aggregated_table.csv 缺少 v2.1 必需列: " + ", ".join(missing))
    return df.copy()


def default_tag_for_row(row: pd.Series) -> dict[str, str]:
    compound = str(row.get("compound", ""))
    lower = norm_name(compound)
    axis = str(row.get("associated_axis", ""))
    mechanism = str(row.get("indicative_mechanism", "")).lower()
    drug_class = str(row.get("indicative_drug_class", "")).lower()
    caution = str(row.get("translational_caution_category", ""))
    best_adjusted_p = float(pd.to_numeric(row.get("best_adjusted_p", 1.0), errors="coerce"))
    support_count = int(float(pd.to_numeric(row.get("support_count", 0), errors="coerce")))
    translational_score = float(pd.to_numeric(row.get("translational_suitability_score", 0.45), errors="coerce"))
    weak_reason = str(row.get("weak_fdr_retained_reason", ""))

    high_risk_categories = {
        "high_caution_toxic_or_research_chemical",
        "oncology_or_broad_cytotoxic_caution",
        "controlled_or_abuse_liability_caution",
        "identifier_or_tool_compound_manual_review",
    }
    tool_terms = {
        "mpep",
        "pma",
        "phorbol",
        "anisomycin",
        "cycloheximide",
        "staurosporine",
        "puromycin",
        "actinomycin",
    }
    broad_epigenetic_terms = {
        "trichostatin",
        "scriptaid",
        "ms 275",
        "ms275",
        "entinostat",
        "vorinostat",
        "tacedinaline",
        "panobinostat",
        "belinostat",
    }
    cns_relevant_mechanisms = [
        "calcium",
        "synaptic",
        "phosphodiesterase",
        "camp",
        "redox",
        "stress",
        "epigenetic",
    ]

    tool_like = caution in high_risk_categories or any(term in lower for term in tool_terms)
    broad_epi = any(term in lower for term in broad_epigenetic_terms)
    if lower == "valproic acid":
        broad_epi = False
        tool_like = False

    weak_fdr = bool(weak_reason and weak_reason.lower() != "nan") or (
        axis == "THRB" and best_adjusted_p > STRICT_FDR and best_adjusted_p <= THRB_WEAK_FDR
    )
    low_context = (
        translational_score <= 0.45
        or caution == "requires_manual_translational_review"
        or (drug_class == "unknown" and caution != "cns_or_neuroactive_known_drug")
        or (support_count <= 1 and not any(token in mechanism for token in cns_relevant_mechanisms))
    )

    notes: list[str] = []
    if tool_like:
        notes.append("工具化合物/高风险或需人工复核类别。")
    if broad_epi:
        notes.append("HDAC/表观遗传调节倾向，作为机制线索保守后移。")
    if weak_fdr:
        notes.append("THRB weak-FDR 或弱 FDR 手工复核候选。")
    if low_context:
        notes.append("CNS/epilepsy 转化语境不足或需人工核查。")
    if lower == "valproic acid":
        notes.append("已知 CNS/抗癫痫相关药物背景；不按 pan-HDAC 工具化合物处理。")

    return {
        "compound": compound,
        "tool_compound_like": bool_text(tool_like),
        "broad_epigenetic_modulator_like": bool_text(broad_epi),
        "weak_fdr_manual_review": bool_text(weak_fdr),
        "low_cns_epilepsy_relevance": bool_text(low_context),
        "axis_preference": axis if axis in PRIMARY_AXES else "supplementary",
        "notes": " ".join(notes) if notes else "auto_v21_default",
    }


def ensure_compound_tags(df: pd.DataFrame) -> pd.DataFrame:
    ensure_dirs()
    default_tags = pd.DataFrame([default_tag_for_row(row) for _, row in df.iterrows()])
    if TAG_CONFIG.exists():
        existing = pd.read_csv(TAG_CONFIG, dtype=str).fillna("")
        missing_cols = [c for c in TAG_COLUMNS if c not in existing.columns]
        for col in missing_cols:
            existing[col] = ""
        existing = existing[TAG_COLUMNS]
        existing_keys = set(existing["compound"].astype(str))
        to_add = default_tags[~default_tags["compound"].astype(str).isin(existing_keys)].copy()
        if not to_add.empty:
            existing = pd.concat([existing, to_add], ignore_index=True)
        tags = existing.drop_duplicates("compound", keep="first")
    else:
        tags = default_tags
    tags = tags[TAG_COLUMNS].sort_values("compound")
    tags.to_csv(TAG_CONFIG, index=False)
    return tags


def merge_tags(df: pd.DataFrame, tags: pd.DataFrame) -> pd.DataFrame:
    merged = df.merge(tags, on="compound", how="left")
    for col in TAG_COLUMNS:
        if col == "compound":
            continue
        if col not in merged.columns:
            merged[col] = ""
    for col in [
        "tool_compound_like",
        "broad_epigenetic_modulator_like",
        "weak_fdr_manual_review",
        "low_cns_epilepsy_relevance",
    ]:
        merged[col] = merged[col].map(as_bool)
    merged["axis_preference"] = merged["axis_preference"].replace("", np.nan).fillna(merged["associated_axis"])
    merged["notes"] = merged["notes"].fillna("")
    if V2_TOP_TABLE.exists():
        v2_top = pd.read_csv(V2_TOP_TABLE, usecols=["compound"])
        top_names = set(v2_top["compound"].astype(str))
    else:
        top_names = set()
    merged["was_v2_top_priority_candidate"] = merged["compound"].astype(str).isin(top_names)
    return merged


def iter_configs() -> list[V21Config]:
    configs: list[V21Config] = []
    for weights, penalties, primary_target, primary_support, thrb_trans in itertools.product(
        WEIGHT_PRESETS,
        PENALTY_PRESETS,
        [4, 5, 6],
        [2, 3],
        [0.75, 0.90],
    ):
        weight_profile = str(weights["weight_profile"])
        penalty_profile = str(penalties["penalty_profile"])
        weight_map = {k: float(v) for k, v in weights.items() if k != "weight_profile"}
        penalty_map = {k: float(v) for k, v in penalties.items() if k != "penalty_profile"}
        total = sum(weight_map.values())
        if not math.isclose(total, 1.0, abs_tol=1e-9):
            raise ValueError(f"权重总和不是 1.0: {weight_profile}={total}")
        configs.append(
            V21Config(
                weights=weight_map,
                penalties=penalty_map,
                primary_target_count=int(primary_target),
                primary_min_support=int(primary_support),
                thrb_supportive_min_translational=float(thrb_trans),
                weight_profile=weight_profile,
                penalty_profile=penalty_profile,
            )
        )
    return configs


def rescore(df: pd.DataFrame, config: V21Config) -> pd.DataFrame:
    scored = df.copy()
    for col in config.weights:
        scored[col] = safe_num(scored[col], 0.0)
    scored["support_count"] = safe_num(scored["support_count"], 0.0)
    scored["support_count_main_axes"] = safe_num(scored.get("support_count_main_axes", pd.Series(0, index=scored.index)), 0.0)
    scored["best_adjusted_p"] = safe_num(scored["best_adjusted_p"], 1.0)
    scored["translational_suitability_score"] = safe_num(scored["translational_suitability_score"], 0.45)
    scored["axis_assignment_confidence_score"] = safe_num(scored["axis_assignment_confidence_score"], 0.0)

    scored["weighted_pre_penalty_score_v21"] = 0.0
    for col, weight in config.weights.items():
        component_col = f"v21_component_{col}"
        scored[component_col] = scored[col] * weight
        scored["weighted_pre_penalty_score_v21"] += scored[component_col]

    weak_fdr_condition = scored["weak_fdr_manual_review"] | (
        (scored["associated_axis"].astype(str) == "THRB")
        & (scored["best_adjusted_p"] > STRICT_FDR)
        & (scored["best_adjusted_p"] <= THRB_WEAK_FDR)
    )
    scored["weak_fdr_penalty_value"] = np.where(weak_fdr_condition, config.penalties["weak_fdr_penalty"], 0.0)
    scored["single_signature_penalty_value"] = np.where(scored["support_count"] <= 1, config.penalties["single_signature_penalty"], 0.0)
    tool_condition = scored["tool_compound_like"] | scored["broad_epigenetic_modulator_like"]
    scored["tool_compound_penalty_value"] = np.where(tool_condition, config.penalties["tool_compound_penalty"], 0.0)
    scored["low_cns_epilepsy_relevance_penalty_value"] = np.where(
        scored["low_cns_epilepsy_relevance"],
        config.penalties["low_cns_epilepsy_relevance_penalty"],
        0.0,
    )
    scored["total_penalty_v21"] = scored[
        [
            "weak_fdr_penalty_value",
            "single_signature_penalty_value",
            "tool_compound_penalty_value",
            "low_cns_epilepsy_relevance_penalty_value",
        ]
    ].sum(axis=1)
    scored["final_score_v21"] = (scored["weighted_pre_penalty_score_v21"] - scored["total_penalty_v21"]).clip(lower=0, upper=1)
    scored["v21_penalty_reasons"] = scored.apply(penalty_reasons, axis=1)
    return scored


def penalty_reasons(row: pd.Series) -> str:
    reasons = []
    if row.get("weak_fdr_penalty_value", 0) > 0:
        reasons.append("weak_fdr")
    if row.get("single_signature_penalty_value", 0) > 0:
        reasons.append("single_signature")
    if row.get("tool_compound_penalty_value", 0) > 0:
        reasons.append("tool_or_broad_epigenetic")
    if row.get("low_cns_epilepsy_relevance_penalty_value", 0) > 0:
        reasons.append("low_cns_epilepsy_relevance")
    return ";".join(reasons) if reasons else "none"


def layer_candidates(scored: pd.DataFrame, config: V21Config) -> pd.DataFrame:
    df = scored.copy()
    df["lead_layer_v21"] = "retained_full_list"
    df["v21_evidence_tier"] = "supplementary"
    df["v21_layer_reason"] = ""

    high_risk = df["translational_caution_category"].astype(str).isin(
        {
            "high_caution_toxic_or_research_chemical",
            "oncology_or_broad_cytotoxic_caution",
            "controlled_or_abuse_liability_caution",
            "identifier_or_tool_compound_manual_review",
        }
    )
    nfe_primary_gate = (
        (df["associated_axis"].astype(str) == "NFE2L2")
        & df["was_v2_top_priority_candidate"].astype(bool)
        & (df["best_adjusted_p"] <= STRICT_FDR)
        & (df["support_count"] >= config.primary_min_support)
        & (df["axis_assignment_confidence_score"] >= 0.62)
        & (df["translational_suitability_score"] >= 0.70)
        & (~df["tool_compound_like"])
        & (~df["broad_epigenetic_modulator_like"])
        & (~df["low_cns_epilepsy_relevance"])
        & (~high_risk)
        & (df["final_score_v21"] >= PRIMARY_SCORE_MIN)
    )
    # THRB weak-FDR candidates remain supportive/manual-review by design in v2.1.
    primary_pool = df[nfe_primary_gate].sort_values(
        ["final_score_v21", "support_count", "axis_assignment_confidence_score", "translational_suitability_score"],
        ascending=[False, False, False, False],
    )
    primary_ids = set(primary_pool.head(config.primary_target_count).index)
    df.loc[list(primary_ids), "lead_layer_v21"] = "primary_mechanism_direction_leads"
    df.loc[list(primary_ids), "v21_evidence_tier"] = "primary_prebbb"
    df.loc[list(primary_ids), "v21_layer_reason"] = (
        "NFE2L2 strict-FDR、多集合/清晰轴线支持，且未被工具化合物或低转化语境标签排除。"
    )

    thrb_supportive_gate = (
        (df["associated_axis"].astype(str) == "THRB")
        & df["was_v2_top_priority_candidate"].astype(bool)
        & (df["best_adjusted_p"] > STRICT_FDR)
        & (df["best_adjusted_p"] <= THRB_WEAK_FDR)
        & (df["support_count"] >= 2)
        & (df["axis_assignment_confidence_score"] >= 0.80)
        & (df["translational_suitability_score"] >= config.thrb_supportive_min_translational)
        & (~df["tool_compound_like"])
        & (~high_risk)
        & (df["final_score_v21"] >= SUPPORTIVE_SCORE_MIN)
    )
    nfe_supportive_gate = (
        (df["associated_axis"].astype(str) == "NFE2L2")
        & df["was_v2_top_priority_candidate"].astype(bool)
        & (~df.index.isin(primary_ids))
        & (df["best_adjusted_p"] <= STRICT_FDR)
        & (df["support_count"] >= 2)
        & (df["axis_assignment_confidence_score"] >= 0.55)
        & (df["translational_suitability_score"] >= 0.70)
        & (~high_risk)
        & (df["final_score_v21"] >= SUPPORTIVE_SCORE_MIN)
    )
    supportive_gate = (thrb_supportive_gate | nfe_supportive_gate) & (~df.index.isin(primary_ids))
    df.loc[supportive_gate, "lead_layer_v21"] = "supportive_manual_review_leads"
    df.loc[supportive_gate, "v21_evidence_tier"] = np.where(
        df.loc[supportive_gate, "associated_axis"].astype(str).eq("THRB"),
        "manual_review_weak_fdr",
        "supportive_prebbb",
    )
    df.loc[supportive_gate, "v21_layer_reason"] = np.where(
        df.loc[supportive_gate, "associated_axis"].astype(str).eq("THRB"),
        "THRB weak-FDR 候选；保留为 supportive/manual-review，不与 NFE2L2 强证据候选混为同一层级。",
        "NFE2L2 候选具一定支持，但因工具/机制/语境标签或排序未进入 primary。",
    )

    retained_gate = (
        (df["final_score_v21"] >= RETAINED_SCORE_MIN)
        | (df["support_count"] >= 2)
        | df["include_in_main_axis_pool"].astype(bool)
    )
    df.loc[~retained_gate, "lead_layer_v21"] = "low_priority_not_retained_v21"
    df.loc[~retained_gate, "v21_evidence_tier"] = "not_retained"
    df.loc[df["v21_layer_reason"].eq(""), "v21_layer_reason"] = (
        "保留在 supplementary/retained full list；未达到 v2.1 primary 或 supportive 门槛。"
    )
    df = df.sort_values(
        ["lead_layer_v21", "final_score_v21", "support_count", "best_adjusted_p"],
        ascending=[True, False, False, True],
    )
    return df


def objective_for(layered: pd.DataFrame, config: V21Config) -> dict[str, Any]:
    primary = layered[layered["lead_layer_v21"] == "primary_mechanism_direction_leads"].copy()
    supportive = layered[layered["lead_layer_v21"] == "supportive_manual_review_leads"].copy()
    retained = layered[layered["lead_layer_v21"].isin(["primary_mechanism_direction_leads", "supportive_manual_review_leads", "retained_full_list"])]
    primary_n = int(primary.shape[0])

    if primary_n:
        support_mean = float(np.minimum(primary["support_count"] / 4.0, 1.0).mean())
        axis_conf_mean = float(primary["axis_assignment_confidence_score"].mean())
        trans_mean = float(primary["translational_suitability_score"].mean())
        final_score_mean = float(primary["final_score_v21"].mean())
        primary_penalty_mean = float(primary["total_penalty_v21"].mean())
        weak_primary = int(primary["weak_fdr_penalty_value"].gt(0).sum())
        tool_primary = int(primary["tool_compound_penalty_value"].gt(0).sum())
        low_primary = int(primary["low_cns_epilepsy_relevance_penalty_value"].gt(0).sum())
    else:
        support_mean = axis_conf_mean = trans_mean = final_score_mean = 0.0
        primary_penalty_mean = 1.0
        weak_primary = tool_primary = low_primary = 10

    count_score = 1.0 if 4 <= primary_n <= 6 else max(0.0, 1.0 - abs(primary_n - 5) / 5.0)
    supportive_thrb_n = int((supportive["associated_axis"].astype(str) == "THRB").sum()) if not supportive.empty else 0
    primary_nfe_n = int((primary["associated_axis"].astype(str) == "NFE2L2").sum()) if not primary.empty else 0
    penalty_primary_count = weak_primary + tool_primary + low_primary
    penalty_conservativeness = sum(config.penalties.values())
    objective = (
        3.0 * count_score
        + 1.4 * support_mean
        + 1.2 * axis_conf_mean
        + 1.4 * trans_mean
        + 1.0 * final_score_mean
        + 0.45 * min(supportive_thrb_n / 2.0, 1.0)
        + 0.20 * min(primary_nfe_n / 2.0, 1.0)
        - 1.5 * penalty_primary_count
        - 0.35 * primary_penalty_mean
    )
    return {
        "config_id": config.config_id,
        "weight_profile": config.weight_profile,
        "penalty_profile": config.penalty_profile,
        "objective_score": objective,
        "primary_count": primary_n,
        "supportive_count": int(supportive.shape[0]),
        "retained_count": int(retained.shape[0]),
        "primary_nfe2l2_count": primary_nfe_n,
        "primary_thrb_count": int((primary["associated_axis"].astype(str) == "THRB").sum()) if primary_n else 0,
        "supportive_nfe2l2_count": int((supportive["associated_axis"].astype(str) == "NFE2L2").sum()) if not supportive.empty else 0,
        "supportive_thrb_count": supportive_thrb_n,
        "primary_support_mean_capped": support_mean,
        "primary_axis_confidence_mean": axis_conf_mean,
        "primary_translational_mean": trans_mean,
        "primary_final_score_mean": final_score_mean,
        "primary_penalty_mean": primary_penalty_mean,
        "weak_fdr_primary_count": weak_primary,
        "tool_like_primary_count": tool_primary,
        "low_context_primary_count": low_primary,
        "primary_target_count": config.primary_target_count,
        "primary_min_support": config.primary_min_support,
        "thrb_supportive_min_translational": config.thrb_supportive_min_translational,
        "penalty_conservativeness_score": penalty_conservativeness,
        **config.weights,
        **config.penalties,
    }


def select_best_config(search: pd.DataFrame) -> pd.Series:
    sort_cols = [
        "objective_score",
        "primary_count",
        "supportive_count",
        "penalty_conservativeness_score",
        "primary_translational_mean",
        "primary_axis_confidence_mean",
        "tool_like_primary_count",
        "weak_fdr_primary_count",
    ]
    ascending = [False, True, True, False, False, False, True, True]
    return search.sort_values(sort_cols, ascending=ascending).iloc[0]


def config_from_search_row(row: pd.Series) -> V21Config:
    weights = {
        "enrichment_support_score": float(row["enrichment_support_score"]),
        "cross_signature_score": float(row["cross_signature_score"]),
        "axis_relevance_score": float(row["axis_relevance_score"]),
        "mechanism_interpretability_score": float(row["mechanism_interpretability_score"]),
        "translational_suitability_score": float(row["translational_suitability_score"]),
        "axis_assignment_confidence_score": float(row["axis_assignment_confidence_score"]),
    }
    penalties = {
        "weak_fdr_penalty": float(row["weak_fdr_penalty"]),
        "single_signature_penalty": float(row["single_signature_penalty"]),
        "tool_compound_penalty": float(row["tool_compound_penalty"]),
        "low_cns_epilepsy_relevance_penalty": float(row["low_cns_epilepsy_relevance_penalty"]),
    }
    return V21Config(
        weights=weights,
        penalties=penalties,
        primary_target_count=int(row["primary_target_count"]),
        primary_min_support=int(row["primary_min_support"]),
        thrb_supportive_min_translational=float(row["thrb_supportive_min_translational"]),
        weight_profile=str(row["weight_profile"]),
        penalty_profile=str(row["penalty_profile"]),
    )


def run_parameter_search(base: pd.DataFrame) -> tuple[pd.DataFrame, V21Config, pd.DataFrame]:
    rows = []
    layered_by_config: dict[str, pd.DataFrame] = {}
    for config in iter_configs():
        layered = layer_candidates(rescore(base, config), config)
        layered_by_config[config.config_id] = layered
        rows.append(objective_for(layered, config))
    search = pd.DataFrame(rows)
    best_row = select_best_config(search)
    best_config = config_from_search_row(best_row)
    best_layered = layered_by_config[str(best_row["config_id"])]
    search["selected_best"] = search["config_id"].astype(str).eq(str(best_row["config_id"]))
    search = search.sort_values(
        [
            "selected_best",
            "objective_score",
            "primary_count",
            "supportive_count",
            "penalty_conservativeness_score",
            "primary_translational_mean",
            "primary_axis_confidence_mean",
        ],
        ascending=[False, False, True, True, False, False, False],
    )
    return search, best_config, best_layered


def quote_yaml_value(value: Any) -> str:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return str(value)
    text = str(value)
    if re.fullmatch(r"[A-Za-z0-9_.\-]+", text):
        return text
    return '"' + text.replace('"', '\\"') + '"'


def write_best_config(config: V21Config, search: pd.DataFrame) -> None:
    best_row = search.iloc[0]
    lines = [
        "v21_best_config:",
        f"  config_id: {quote_yaml_value(config.config_id)}",
        f"  selected_at: {quote_yaml_value(datetime.now().isoformat(timespec='seconds'))}",
        "  source_table: " + quote_yaml_value(str(SOURCE_TABLE)),
        "  output_dir: " + quote_yaml_value(str(OUT_DIR)),
        "  objective:",
        f"    objective_score: {float(best_row['objective_score']):.6f}",
        "    priority: primary quality, conservative penalties, supportive THRB kept separate",
        "  weights:",
    ]
    for key, value in config.weights.items():
        lines.append(f"    {key}: {value}")
    lines.append("  penalties:")
    for key, value in config.penalties.items():
        lines.append(f"    {key}: {value}")
    lines.extend(
        [
            "  thresholds:",
            f"    strict_fdr: {STRICT_FDR}",
            f"    thrb_weak_fdr: {THRB_WEAK_FDR}",
            f"    primary_target_count: {config.primary_target_count}",
            f"    primary_min_support: {config.primary_min_support}",
            f"    primary_min_score: {PRIMARY_SCORE_MIN}",
            f"    supportive_min_score: {SUPPORTIVE_SCORE_MIN}",
            f"    thrb_supportive_min_translational: {config.thrb_supportive_min_translational}",
            "  constraints:",
            "    - BBB is not used in v2.1 pre-BBB scoring.",
            "    - THRB weak-FDR candidates cannot enter primary; they are supportive/manual-review only.",
            "    - BHLHE40 is supplementary support only; SOX2 is not used for drug-query mainline.",
        ]
    )
    (OUT_DIR / "02_best_config.yaml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def add_canonical_smiles(df: pd.DataFrame) -> pd.DataFrame:
    if "canonical_smiles" in df.columns:
        return df
    enriched = df.copy()
    enriched["canonical_smiles"] = ""
    if V2_INTEGRATED_TABLE.exists():
        v2 = pd.read_csv(V2_INTEGRATED_TABLE)
        if {"compound", "canonical_smiles"}.issubset(v2.columns):
            smiles = v2[["compound", "canonical_smiles"]].dropna().drop_duplicates("compound")
            enriched = enriched.drop(columns=["canonical_smiles"]).merge(smiles, on="compound", how="left")
            enriched["canonical_smiles"] = enriched["canonical_smiles"].fillna("")
    return enriched


def paper_ready(primary: pd.DataFrame, supportive: pd.DataFrame) -> pd.DataFrame:
    rows = []
    combined = pd.concat([primary, supportive], ignore_index=True) if not supportive.empty else primary.copy()
    for _, r in combined.iterrows():
        role = (
            "primary mechanism-direction lead"
            if r["lead_layer_v21"] == "primary_mechanism_direction_leads"
            else "supportive/manual-review mechanism clue"
        )
        rows.append(
            {
                "compound": r["compound"],
                "axis": r["associated_axis"],
                "lead_layer_v21": r["lead_layer_v21"],
                "evidence_tier": r["v21_evidence_tier"],
                "support_summary": (
                    f"{int(r['support_count'])} query set(s); "
                    f"best adj.P={float(r['best_adjusted_p']):.3g}; "
                    f"combined={float(r['best_combined_score']):.3g}"
                ),
                "mechanism_theme": r.get("indicative_mechanism", ""),
                "v21_score_prebbb": round(float(r["final_score_v21"]), 6),
                "penalty_reasons": r.get("v21_penalty_reasons", ""),
                "recommended_role": role,
                "conservative_interpretation": (
                    "Exploratory program-guided mechanism-direction clue only; not a clinical recommendation "
                    "and not direct target validation."
                ),
                "manual_bbb_status": "BBB not resolved in v2.1; requires manual search/review.",
            }
        )
    return pd.DataFrame(rows)


def write_manual_bbb_template(primary: pd.DataFrame, supportive: pd.DataFrame) -> None:
    candidates = pd.concat([primary, supportive], ignore_index=True) if not supportive.empty else primary.copy()
    candidates = add_canonical_smiles(candidates)
    cols = [
        "compound",
        "associated_axis",
        "lead_layer_v21",
        "v21_evidence_tier",
        "canonical_smiles",
        "manual_bbb_layer",
        "manual_bbb_evidence_source",
        "manual_bbb_note",
    ]
    template = pd.DataFrame(
        {
            "compound": candidates["compound"],
            "associated_axis": candidates["associated_axis"],
            "lead_layer_v21": candidates["lead_layer_v21"],
            "v21_evidence_tier": candidates["v21_evidence_tier"],
            "canonical_smiles": candidates.get("canonical_smiles", ""),
            "manual_bbb_layer": "NA",
            "manual_bbb_evidence_source": "NA",
            "manual_bbb_note": "NA",
        }
    )
    template[cols].to_csv(OUT_DIR / "09_manual_bbb_template_v21.tsv", sep="\t", index=False)


def write_index(best_config: V21Config, search: pd.DataFrame, primary: pd.DataFrame, supportive: pd.DataFrame) -> None:
    lines = [
        "drug repositioning v2.1 pre-BBB 轻量收敛版运行索引",
        "=" * 60,
        f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"输入 pre-BBB 聚合表: {SOURCE_TABLE}",
        f"输出目录: {OUT_DIR}",
        "",
        "本轮做了什么:",
        "- 复用现有 v2 aggregated/integrated 结果，在本地重算 pre-BBB 分数。",
        "- 新增显式 penalty: weak-FDR、single-signature、tool/broad-epigenetic、low-CNS/epilepsy-context。",
        "- 将结果收敛为 primary / supportive / retained full list 三层。",
        "- primary/supportive 层仅从原 v2 top-priority 候选中收敛产生，避免把 retained 中的新药名升入主展示层。",
        "- 运行 120 组确定性有限参数搜索，并按 objective 自动选择最佳配置。",
        "",
        "本轮没做什么:",
        "- 未重跑 pySCENIC、CellOracle、KO、外部验证、robustness、GO/KEGG。",
        "- 未重跑 Enrichr/DSigDB 原始检索，未下载新资源。",
        "- 未执行 docking，未调用任何 docking 脚本。",
        "- BBB 仍为后续手动搜索与整理，本轮只生成 pre-BBB 结果和手动 BBB 模板。",
        "",
        "最佳配置:",
        f"- config_id: {best_config.config_id}",
        f"- objective_score: {float(search.iloc[0]['objective_score']):.6f}",
        f"- primary leads: {primary.shape[0]}",
        f"- supportive leads: {supportive.shape[0]}",
        "",
        "关键输出:",
        "- 01_parameter_search_summary.csv",
        "- 02_best_config.yaml",
        "- 03_integrated_rescored_v21.csv",
        "- 04_primary_mechanism_direction_leads.csv",
        "- 05_supportive_manual_review_leads.csv",
        "- 06_retained_full_list_v21.csv",
        "- 07_paper_ready_table_v21.csv",
        "- 08_v21_summary_cn.txt",
        "- 09_manual_bbb_template_v21.tsv",
        "- 10_compound_tag_template_v21.csv",
    ]
    (OUT_DIR / "00_v21_run_index_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def summarize(best_config: V21Config, search: pd.DataFrame, layered: pd.DataFrame, primary: pd.DataFrame, supportive: pd.DataFrame, retained: pd.DataFrame) -> None:
    top_next = search[~search["selected_best"]].head(5)[
        [
            "config_id",
            "objective_score",
            "primary_count",
            "supportive_count",
            "primary_translational_mean",
            "primary_axis_confidence_mean",
            "supportive_thrb_count",
            "tool_like_primary_count",
            "weak_fdr_primary_count",
        ]
    ]
    v2_top_path = DRUGREP_DIR / "05_top_priority_compounds.csv"
    v2_top = pd.read_csv(v2_top_path) if v2_top_path.exists() else pd.DataFrame()
    v2_top_names = set(v2_top["compound"].astype(str)) if not v2_top.empty else set()
    primary_names = set(primary["compound"].astype(str))
    supportive_names = set(supportive["compound"].astype(str))
    downgraded_from_v2 = layered[layered["compound"].astype(str).isin(v2_top_names - primary_names)].copy()
    downgraded_from_v2 = downgraded_from_v2.sort_values(["lead_layer_v21", "final_score_v21"], ascending=[True, False])

    lines = [
        "drug repositioning v2.1 pre-BBB 轻量收敛版总结",
        "=" * 62,
        "科学定位:",
        "- 结果仅作为 program-guided exploratory drug repositioning / mechanism-direction clues。",
        "- 不解释为 FCD II 特异治疗药发现、直接靶点验证或临床治疗推荐。",
        "- NFE2L2 / THRB 是主轴；BHLHE40 仅作补充支持；SOX2 不进入主药物查询主线。",
        "",
        "最佳配置:",
        f"- config_id: {best_config.config_id}",
        f"- objective_score: {float(search.iloc[0]['objective_score']):.6f}",
        "- 权重:",
    ]
    for k, v in best_config.weights.items():
        lines.append(f"  - {k}: {v}")
    lines.append("- penalty:")
    for k, v in best_config.penalties.items():
        lines.append(f"  - {k}: {v}")
    lines.extend(
        [
            f"- primary_target_count: {best_config.primary_target_count}",
            f"- primary_min_support: {best_config.primary_min_support}",
            f"- THRB supportive min translational_suitability: {best_config.thrb_supportive_min_translational}",
            "",
            "为什么选择这一版:",
            "- primary 数量落在 4-6 的目标范围内。",
            "- primary 层没有 weak-FDR、tool-like 或 low-context 候选。",
            "- THRB weak-FDR 候选被保留为 supportive/manual-review，而不是和 NFE2L2 强证据候选混为同一层。",
            "- primary/supportive 均限制在原 v2 top-priority 候选内，避免 v2.1 变成扩展药名清单。",
            "- 当 objective 并列时，优先选择 conservative_high penalty profile，并保持 primary/supportive 数量更收敛。",
            "- 相比 v2 的 12 个大主表，v2.1 更强调可解释性、转化适配性、轴线清晰和显式 penalty。",
            "",
            "次优配置对比前 5:",
            top_next.to_string(index=False),
            "",
            f"primary_mechanism_direction_leads 数量: {primary.shape[0]}",
            primary[
                [
                    "compound",
                    "associated_axis",
                    "support_count",
                    "best_adjusted_p",
                    "final_score_v21",
                    "v21_penalty_reasons",
                    "v21_layer_reason",
                ]
            ].to_string(index=False)
            if not primary.empty
            else "无",
            "",
            f"supportive_manual_review_leads 数量: {supportive.shape[0]}",
            supportive[
                [
                    "compound",
                    "associated_axis",
                    "support_count",
                    "best_adjusted_p",
                    "final_score_v21",
                    "v21_evidence_tier",
                    "v21_penalty_reasons",
                ]
            ].head(40).to_string(index=False)
            if not supportive.empty
            else "无",
            "",
            "相对 v2 主表被降级/后移的候选:",
            downgraded_from_v2[
                [
                    "compound",
                    "associated_axis",
                    "lead_layer_v21",
                    "v21_evidence_tier",
                    "final_score_v21",
                    "v21_penalty_reasons",
                    "v21_layer_reason",
                ]
            ].to_string(index=False)
            if not downgraded_from_v2.empty
            else "无",
            "",
            f"retained_full_list / supplementary 总数: {retained.shape[0]}",
            "",
            "BBB 状态:",
            "- v2.1 未使用 BBB 作为评分输入。",
            "- 已生成 09_manual_bbb_template_v21.tsv，供后续人工搜索和整理。",
            "- 手动 BBB 后应另行运行 after-BBB 继续整合，不应回写本轮 pre-BBB 评分逻辑。",
        ]
    )
    (OUT_DIR / "08_v21_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_outputs(search: pd.DataFrame, best_config: V21Config, layered: pd.DataFrame) -> None:
    ensure_dirs()
    search.to_csv(OUT_DIR / "01_parameter_search_summary.csv", index=False)
    write_best_config(best_config, search)

    ordered = layered.sort_values(["lead_layer_v21", "final_score_v21", "support_count"], ascending=[True, False, False])
    ordered.to_csv(OUT_DIR / "03_integrated_rescored_v21.csv", index=False)
    primary = ordered[ordered["lead_layer_v21"] == "primary_mechanism_direction_leads"].copy()
    supportive = ordered[ordered["lead_layer_v21"] == "supportive_manual_review_leads"].copy()
    retained = ordered[ordered["lead_layer_v21"].isin(["primary_mechanism_direction_leads", "supportive_manual_review_leads", "retained_full_list"])].copy()

    primary.to_csv(OUT_DIR / "04_primary_mechanism_direction_leads.csv", index=False)
    supportive.to_csv(OUT_DIR / "05_supportive_manual_review_leads.csv", index=False)
    retained.to_csv(OUT_DIR / "06_retained_full_list_v21.csv", index=False)
    paper_ready(primary, supportive).to_csv(OUT_DIR / "07_paper_ready_table_v21.csv", index=False)
    write_manual_bbb_template(primary, supportive)

    tags = pd.read_csv(TAG_CONFIG, dtype=str).fillna("")
    tags.to_csv(OUT_DIR / "10_compound_tag_template_v21.csv", index=False)
    write_index(best_config, search, primary, supportive)
    summarize(best_config, search, ordered, primary, supportive, retained)


def run_v21() -> tuple[pd.DataFrame, V21Config, pd.DataFrame]:
    ensure_dirs()
    write_static_configs()
    source = load_source_table()
    tags = ensure_compound_tags(source)
    base = merge_tags(source, tags)
    search, best_config, best_layered = run_parameter_search(base)
    export_outputs(search, best_config, best_layered)
    return search, best_config, best_layered


def main() -> None:
    print("开始 drug repositioning v2.1 pre-BBB 轻量收敛；不重跑上游，不做 BBB 自动搜索，不做 docking。")
    search, best_config, layered = run_v21()
    primary_n = int((layered["lead_layer_v21"] == "primary_mechanism_direction_leads").sum())
    supportive_n = int((layered["lead_layer_v21"] == "supportive_manual_review_leads").sum())
    print(f"v2.1 完成: 搜索配置 {search.shape[0]} 组。")
    print(f"最佳配置: {best_config.config_id}")
    print(f"primary={primary_n}; supportive={supportive_n}; 输出目录={OUT_DIR}")


if __name__ == "__main__":
    main()
