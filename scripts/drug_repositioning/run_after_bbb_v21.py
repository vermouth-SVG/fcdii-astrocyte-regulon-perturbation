#!/usr/bin/env python3
from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from drugrep_common_v2 import ROOT


OUT_DIR = ROOT / "drug_repositioning"
SWISSADME_PATH = OUT_DIR / "swissadme.csv"
INTEGRATED_PREBBB = OUT_DIR / "03_integrated_rescored_v21.csv"
MANUAL_BBB_TEMPLATE = OUT_DIR / "09_manual_bbb_template_v21.tsv"

BBB_LAYER_SCORES = {
    "CNS-directed candidates": 1.00,
    "possible CNS-directed candidates": 0.75,
    "peripheral/program-modulating candidates": 0.40,
}

SWISSADME_KEY_COLUMNS = [
    "Molecule",
    "Canonical SMILES",
    "GI absorption",
    "BBB permeant",
    "Pgp substrate",
    "Bioavailability Score",
    "PAINS #alerts",
    "Brenk #alerts",
    "Leadlikeness #violations",
    "Synthetic Accessibility",
]


def is_blank(value: Any) -> bool:
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "na", "none", "null"}


def clean(value: Any) -> str:
    return "NA" if is_blank(value) else str(value).strip()


def norm_smiles(value: Any) -> str:
    if is_blank(value):
        return ""
    return re.sub(r"\s+", "", str(value).strip())


def norm_col(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def yes_no(value: Any) -> str:
    text = clean(value).lower()
    if text in {"yes", "y", "true", "1"}:
        return "Yes"
    if text in {"no", "n", "false", "0"}:
        return "No"
    return clean(value)


def numeric_or_none(value: Any) -> float | None:
    if is_blank(value):
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


def read_csv_flexible(path: Path) -> pd.DataFrame:
    last_exc: Exception | None = None
    for encoding in ["utf-8-sig", "utf-8", "gb18030"]:
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError as exc:
            last_exc = exc
    if last_exc is not None:
        raise last_exc
    return pd.read_csv(path)


def find_column(df: pd.DataFrame, aliases: list[str]) -> str | None:
    normalized = {norm_col(col): col for col in df.columns}
    for alias in aliases:
        found = normalized.get(norm_col(alias))
        if found is not None:
            return found
    return None


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in [SWISSADME_PATH, INTEGRATED_PREBBB, MANUAL_BBB_TEMPLATE]:
        if not path.exists():
            raise FileNotFoundError(f"缺少 v2.1 after-BBB 必需输入: {path}")
    integrated = pd.read_csv(INTEGRATED_PREBBB)
    template = pd.read_csv(MANUAL_BBB_TEMPLATE, sep="\t")
    swiss = read_csv_flexible(SWISSADME_PATH)
    if "canonical_smiles" not in template.columns:
        raise ValueError("09_manual_bbb_template_v21.tsv 缺少 canonical_smiles 列。")
    return integrated, template, swiss


def classify_bbb(record: dict[str, Any]) -> tuple[str, float, str]:
    bbb = yes_no(record.get("BBB permeant", "NA"))
    gi = clean(record.get("GI absorption", "NA"))
    pgp = yes_no(record.get("Pgp substrate", "NA"))
    brenk = numeric_or_none(record.get("Brenk #alerts", "NA"))
    pains = numeric_or_none(record.get("PAINS #alerts", "NA"))
    bioavailability = clean(record.get("Bioavailability Score", "NA"))

    reasons: list[str] = []
    if bbb == "Yes":
        if gi.lower() == "high" and pgp == "No" and (brenk is None or brenk <= 1) and (pains is None or pains == 0):
            layer = "CNS-directed candidates"
            reasons.append("BBB=Yes、GI=High、Pgp=No，且 Brenk/PAINS 未见明显异常。")
        else:
            layer = "possible CNS-directed candidates"
            reasons.append("BBB=Yes，但存在 Pgp/Brenk/PAINS 等限制，保守归为 possible CNS-directed。")
            if gi.lower() != "high":
                reasons.append(f"GI absorption={gi}。")
            if pgp == "Yes":
                reasons.append("Pgp substrate=Yes。")
            if brenk is not None and brenk > 1:
                reasons.append(f"Brenk alerts={int(brenk)}。")
            if pains is not None and pains > 0:
                reasons.append(f"PAINS alerts={int(pains)}。")
    elif bbb == "No":
        layer = "peripheral/program-modulating candidates"
        reasons.append("BBB=No；不删除候选，改按 peripheral/program-modulating clue 解释。")
    else:
        layer = "peripheral/program-modulating candidates"
        reasons.append("BBB 结果缺失或不可识别；按保守规则归入 peripheral/program-modulating clue。")

    if bioavailability != "NA":
        reasons.append(f"Bioavailability Score={bioavailability}。")
    return layer, BBB_LAYER_SCORES[layer], " ".join(reasons)


def match_swissadme(swiss: pd.DataFrame, template: pd.DataFrame) -> pd.DataFrame:
    molecule_col = find_column(swiss, ["Molecule"])
    canonical_col = find_column(swiss, ["Canonical SMILES", "CanonicalSMILES", "SMILES"])
    if molecule_col is None and canonical_col is None:
        raise ValueError("SwissADME CSV 缺少 Molecule/SMILES 可匹配列。")

    ref = template.copy().reset_index(drop=True)
    ref["_canonical_norm"] = ref["canonical_smiles"].map(norm_smiles)
    ref_by_smiles = {
        row["_canonical_norm"]: idx
        for idx, row in ref.iterrows()
        if row["_canonical_norm"] and ref["_canonical_norm"].tolist().count(row["_canonical_norm"]) == 1
    }

    used: set[int] = set()
    rows: list[dict[str, Any]] = []
    for swiss_idx, row in swiss.reset_index(drop=True).iterrows():
        smiles_candidates: list[tuple[str, str]] = []
        if molecule_col is not None:
            smiles_candidates.append(("Molecule", norm_smiles(row.get(molecule_col, ""))))
        if canonical_col is not None:
            smiles_candidates.append(("Canonical SMILES", norm_smiles(row.get(canonical_col, ""))))

        matched_idx: int | None = None
        matched_from = ""
        for source_col, smiles in smiles_candidates:
            if smiles and smiles in ref_by_smiles and ref_by_smiles[smiles] not in used:
                matched_idx = ref_by_smiles[smiles]
                matched_from = source_col
                break

        match_method = "unmatched"
        confidence = 0.0
        if matched_idx is not None:
            used.add(matched_idx)
            match_method = "smiles_exact"
            confidence = 1.0
            matched_ref = ref.loc[matched_idx]
            compound = clean(matched_ref["compound"])
            associated_axis = clean(matched_ref.get("associated_axis", "NA"))
            prebbb_layer = clean(matched_ref.get("lead_layer_v21", "NA"))
            evidence_tier = clean(matched_ref.get("v21_evidence_tier", "NA"))
            canonical_smiles = clean(matched_ref["canonical_smiles"])
        elif swiss.shape[0] == ref.shape[0] and swiss_idx < ref.shape[0] and swiss_idx not in used:
            matched_idx = int(swiss_idx)
            used.add(matched_idx)
            matched_from = "import_order"
            match_method = "import_order_fallback"
            confidence = 0.70
            matched_ref = ref.loc[matched_idx]
            compound = clean(matched_ref["compound"])
            associated_axis = clean(matched_ref.get("associated_axis", "NA"))
            prebbb_layer = clean(matched_ref.get("lead_layer_v21", "NA"))
            evidence_tier = clean(matched_ref.get("v21_evidence_tier", "NA"))
            canonical_smiles = clean(matched_ref["canonical_smiles"])
        else:
            compound = "NA"
            associated_axis = "NA"
            prebbb_layer = "NA"
            evidence_tier = "NA"
            canonical_smiles = "NA"

        record = {
            "swissadme_row_id": int(swiss_idx) + 1,
            "swissadme_molecule": clean(row.get(molecule_col, "NA")) if molecule_col else "NA",
            "swissadme_canonical_smiles": clean(row.get(canonical_col, "NA")) if canonical_col else "NA",
            "compound": compound,
            "associated_axis": associated_axis,
            "prebbb_layer_v21": prebbb_layer,
            "prebbb_evidence_tier_v21": evidence_tier,
            "canonical_smiles": canonical_smiles,
            "matching_method": match_method,
            "matched_from_column": matched_from if matched_from else "NA",
            "matching_confidence": confidence,
        }
        for col in SWISSADME_KEY_COLUMNS:
            if col in {"Molecule", "Canonical SMILES"}:
                continue
            record[col] = clean(row.get(col, "NA"))
        layer, score, note = classify_bbb(record)
        record["manual_final_bbb_layer"] = layer
        record["BBB_score"] = score
        record["bbb_review_note"] = note
        rows.append(record)
    return pd.DataFrame(rows)


def make_bbb_annotation(matched: pd.DataFrame, template: pd.DataFrame) -> pd.DataFrame:
    by_compound = matched.drop_duplicates("compound").set_index("compound")
    rows: list[dict[str, Any]] = []
    for _, t in template.iterrows():
        compound = str(t["compound"])
        if compound in by_compound.index:
            m = by_compound.loc[compound]
            rows.append(
                {
                    "compound": compound,
                    "associated_axis": clean(t.get("associated_axis", "NA")),
                    "prebbb_layer_v21": clean(t.get("lead_layer_v21", "NA")),
                    "prebbb_evidence_tier_v21": clean(t.get("v21_evidence_tier", "NA")),
                    "canonical_smiles": clean(t.get("canonical_smiles", "NA")),
                    "swissadme_bbb_prediction": clean(m.get("BBB permeant", "NA")),
                    "swissadme_gi_absorption": clean(m.get("GI absorption", "NA")),
                    "swissadme_pgp_substrate": clean(m.get("Pgp substrate", "NA")),
                    "swissadme_bioavailability_score": clean(m.get("Bioavailability Score", "NA")),
                    "swissadme_pains_alerts": clean(m.get("PAINS #alerts", "NA")),
                    "swissadme_brenk_alerts": clean(m.get("Brenk #alerts", "NA")),
                    "swissadme_leadlikeness_violations": clean(m.get("Leadlikeness #violations", "NA")),
                    "swissadme_synthetic_accessibility": clean(m.get("Synthetic Accessibility", "NA")),
                    "manual_final_bbb_layer": clean(m.get("manual_final_bbb_layer", "NA")),
                    "BBB_score": float(m.get("BBB_score", np.nan)),
                    "manual_review_note": clean(m.get("bbb_review_note", "NA")),
                    "reviewed_by": "swissadme_csv_v21_manual_input",
                    "reviewed_date": datetime.now().date().isoformat(),
                }
            )
        else:
            rows.append(
                {
                    "compound": compound,
                    "associated_axis": clean(t.get("associated_axis", "NA")),
                    "prebbb_layer_v21": clean(t.get("lead_layer_v21", "NA")),
                    "prebbb_evidence_tier_v21": clean(t.get("v21_evidence_tier", "NA")),
                    "canonical_smiles": clean(t.get("canonical_smiles", "NA")),
                    "swissadme_bbb_prediction": "NA",
                    "swissadme_gi_absorption": "NA",
                    "swissadme_pgp_substrate": "NA",
                    "swissadme_bioavailability_score": "NA",
                    "swissadme_pains_alerts": "NA",
                    "swissadme_brenk_alerts": "NA",
                    "swissadme_leadlikeness_violations": "NA",
                    "swissadme_synthetic_accessibility": "NA",
                    "manual_final_bbb_layer": "not_reviewed_bbb_pending",
                    "BBB_score": np.nan,
                    "manual_review_note": "未在 SwissADME CSV 中匹配；保留待人工核查。",
                    "reviewed_by": "swissadme_csv_v21_manual_input",
                    "reviewed_date": datetime.now().date().isoformat(),
                }
            )
    return pd.DataFrame(rows)


def integrate_after_bbb(integrated: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    df = integrated.copy()
    if "prebbb_layer_v21" not in df.columns:
        df["prebbb_layer_v21"] = df["lead_layer_v21"]
    else:
        df["prebbb_layer_v21"] = df["lead_layer_v21"]
    df["prebbb_evidence_tier_v21"] = df["v21_evidence_tier"]
    df["final_score_prebbb_v21"] = pd.to_numeric(df["final_score_v21"], errors="coerce").fillna(0)

    keep_cols = [
        "compound",
        "canonical_smiles",
        "swissadme_bbb_prediction",
        "swissadme_gi_absorption",
        "swissadme_pgp_substrate",
        "swissadme_bioavailability_score",
        "swissadme_pains_alerts",
        "swissadme_brenk_alerts",
        "swissadme_leadlikeness_violations",
        "swissadme_synthetic_accessibility",
        "manual_final_bbb_layer",
        "BBB_score",
        "manual_review_note",
        "reviewed_by",
        "reviewed_date",
    ]
    ann = annotation[keep_cols].drop_duplicates("compound")
    drop_cols = [c for c in keep_cols if c != "compound" and c in df.columns]
    df = df.drop(columns=drop_cols).merge(ann, on="compound", how="left")

    df["after_bbb_review_status"] = np.where(df["manual_final_bbb_layer"].notna(), "reviewed_in_v21_swissadme", "not_reviewed_bbb_pending")
    df["manual_final_bbb_layer"] = df["manual_final_bbb_layer"].fillna("not_reviewed_bbb_pending")
    df["BBB_score"] = pd.to_numeric(df["BBB_score"], errors="coerce")
    reviewed = df["manual_final_bbb_layer"].isin(BBB_LAYER_SCORES.keys())
    df["final_score_after_bbb"] = np.where(
        reviewed,
        0.88 * df["final_score_prebbb_v21"] + 0.12 * df["BBB_score"],
        df["final_score_prebbb_v21"],
    )

    df["after_bbb_final_layer"] = "retained_full_list_after_bbb"
    headline_gate = (
        df["prebbb_layer_v21"].eq("primary_mechanism_direction_leads")
        & df["manual_final_bbb_layer"].eq("CNS-directed candidates")
    )
    df.loc[headline_gate, "after_bbb_final_layer"] = "headline_cns_mechanism_direction_leads"

    supportive_gate = (
        df["manual_final_bbb_layer"].isin(["CNS-directed candidates", "possible CNS-directed candidates"])
        & ~headline_gate
        & df["prebbb_layer_v21"].isin(["primary_mechanism_direction_leads", "supportive_manual_review_leads"])
    )
    df.loc[supportive_gate, "after_bbb_final_layer"] = "supportive_after_bbb_leads"

    peripheral_gate = (
        df["manual_final_bbb_layer"].eq("peripheral/program-modulating candidates")
        & df["prebbb_layer_v21"].isin(["primary_mechanism_direction_leads", "supportive_manual_review_leads"])
    )
    df.loc[peripheral_gate, "after_bbb_final_layer"] = "peripheral_program_modulating_clues"

    pending_gate = df["manual_final_bbb_layer"].eq("not_reviewed_bbb_pending")
    df.loc[pending_gate & df["prebbb_layer_v21"].eq("low_priority_not_retained_v21"), "after_bbb_final_layer"] = "low_priority_not_retained_v21"
    df.loc[pending_gate & ~df["prebbb_layer_v21"].eq("low_priority_not_retained_v21"), "after_bbb_final_layer"] = "retained_full_list_after_bbb"

    df["after_bbb_interpretation_note"] = df.apply(after_bbb_note, axis=1)
    df = df.sort_values(
        ["after_bbb_final_layer", "final_score_after_bbb", "final_score_prebbb_v21", "support_count"],
        ascending=[True, False, False, False],
    )
    df["after_bbb_rank"] = np.arange(1, df.shape[0] + 1)
    return df


def after_bbb_note(row: pd.Series) -> str:
    bbb_layer = str(row.get("manual_final_bbb_layer", ""))
    pre_layer = str(row.get("prebbb_layer_v21", ""))
    if pre_layer == "primary_mechanism_direction_leads" and bbb_layer == "CNS-directed candidates":
        return "pre-BBB primary 且 SwissADME 支持 CNS-directed；适合 after-BBB headline 展示，但仍为 exploratory clue。"
    if pre_layer == "primary_mechanism_direction_leads" and bbb_layer == "peripheral/program-modulating candidates":
        return "pre-BBB 程序证据较强，但 BBB=No；after-BBB 降为 peripheral/program-modulating mechanism clue。"
    if bbb_layer == "possible CNS-directed candidates":
        return "SwissADME 支持 BBB permeant 但存在 Pgp/Brenk/PAINS 等限制；保留为 possible/supportive。"
    if bbb_layer == "peripheral/program-modulating candidates":
        return "BBB=No 或整体更适合外周/程序调节解释；不删除，只降低 CNS headline 级别。"
    if bbb_layer == "CNS-directed candidates":
        return "SwissADME 支持 CNS-directed，但根据 pre-BBB 层级保守展示。"
    return "未进行 v2.1 SwissADME/BBB 复核；保留在 retained/supplementary。"


def make_primary_after_bbb(after: pd.DataFrame) -> pd.DataFrame:
    return after[after["after_bbb_final_layer"] == "headline_cns_mechanism_direction_leads"].copy()


def make_supportive_after_bbb(after: pd.DataFrame) -> pd.DataFrame:
    return after[
        after["after_bbb_final_layer"].isin(
            [
                "supportive_after_bbb_leads",
                "peripheral_program_modulating_clues",
            ]
        )
    ].copy()


def make_paper_ready(primary: pd.DataFrame, supportive: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    combined = pd.concat([primary, supportive], ignore_index=True)
    for _, r in combined.iterrows():
        if r["after_bbb_final_layer"] == "headline_cns_mechanism_direction_leads":
            display_role = "after-BBB headline CNS-directed mechanism-direction lead"
        elif r["after_bbb_final_layer"] == "peripheral_program_modulating_clues":
            display_role = "peripheral/program-modulating mechanism clue"
        else:
            display_role = "supportive/manual-review after-BBB clue"
        rows.append(
            {
                "compound": r["compound"],
                "axis": r["associated_axis"],
                "prebbb_layer_v21": r["prebbb_layer_v21"],
                "after_bbb_final_layer": r["after_bbb_final_layer"],
                "BBB_layer": r["manual_final_bbb_layer"],
                "prebbb_score_v21": round(float(r["final_score_prebbb_v21"]), 6),
                "BBB_score": "" if pd.isna(r["BBB_score"]) else round(float(r["BBB_score"]), 3),
                "final_score_after_bbb": round(float(r["final_score_after_bbb"]), 6),
                "support_summary": (
                    f"{int(r['support_count'])} query set(s); "
                    f"best adj.P={float(r['best_adjusted_p']):.3g}; "
                    f"combined={float(r['best_combined_score']):.3g}"
                ),
                "mechanism_theme": r.get("indicative_mechanism", ""),
                "penalty_reasons": r.get("v21_penalty_reasons", ""),
                "display_role": display_role,
                "conservative_interpretation": (
                    "Exploratory program-guided mechanism-direction clue only; not a treatment recommendation "
                    "and not direct target validation."
                ),
                "bbb_note": r.get("manual_review_note", ""),
            }
        )
    return pd.DataFrame(rows)


def write_summary(matched: pd.DataFrame, annotation: pd.DataFrame, after: pd.DataFrame, primary: pd.DataFrame, supportive: pd.DataFrame) -> None:
    matched_count = int(matched["matching_method"].ne("unmatched").sum())
    total = int(matched.shape[0])
    method_counts = matched["matching_method"].value_counts(dropna=False)
    confidence_summary = matched["matching_confidence"].value_counts(dropna=False).sort_index()
    layer_counts = annotation["manual_final_bbb_layer"].value_counts(dropna=False)
    pre_primary = after[after["prebbb_layer_v21"].eq("primary_mechanism_direction_leads")].copy()
    downgraded_primary = pre_primary[pre_primary["after_bbb_final_layer"].eq("peripheral_program_modulating_clues")].copy()

    lines = [
        "drug repositioning v2.1 after-BBB 更新总结",
        "=" * 58,
        f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"SwissADME 输入: {SWISSADME_PATH}",
        f"pre-BBB 输入: {INTEGRATED_PREBBB}",
        "",
        "流程边界:",
        "- 只导入本地 SwissADME CSV，完成 BBB 分层与 after-BBB 排序。",
        "- 未重跑 pySCENIC、CellOracle、robustness、GO/KEGG、外部验证、Enrichr/DSigDB。",
        "- 未联网，未下载资源，未做 docking。",
        "- 结果仍为 exploratory / mechanism-direction module，不是治疗推荐或直接靶点验证。",
        "",
        "匹配结果:",
        f"- 匹配成功: {matched_count} / {total}",
        "- matching_method 计数:",
        method_counts.to_string(),
        "- matching_confidence 计数:",
        confidence_summary.to_string(),
    ]
    if matched_count == total and set(method_counts.index.astype(str)) == {"smiles_exact"}:
        lines.extend(
            [
                "- 结论: 12/12 均为 SMILES 精确匹配。",
                "- matching_method = smiles_exact",
                "- matching_confidence = 1.0",
            ]
        )
    lines.extend(
        [
            "",
            "BBB 分层计数:",
            layer_counts.to_string(),
            "",
            "after-BBB headline 候选:",
            primary[
                [
                    "compound",
                    "associated_axis",
                    "prebbb_layer_v21",
                    "manual_final_bbb_layer",
                    "final_score_prebbb_v21",
                    "final_score_after_bbb",
                    "after_bbb_interpretation_note",
                ]
            ].to_string(index=False)
            if not primary.empty
            else "无",
            "",
            "pre-BBB primary 中 after-BBB 降为 peripheral/program-modulating clue 的候选:",
            downgraded_primary[
                [
                    "compound",
                    "associated_axis",
                    "manual_final_bbb_layer",
                    "swissadme_bbb_prediction",
                    "swissadme_gi_absorption",
                    "after_bbb_interpretation_note",
                ]
            ].to_string(index=False)
            if not downgraded_primary.empty
            else "无",
            "",
            "supportive / peripheral after-BBB 候选:",
            supportive[
                [
                    "compound",
                    "associated_axis",
                    "prebbb_layer_v21",
                    "manual_final_bbb_layer",
                    "final_score_after_bbb",
                    "after_bbb_final_layer",
                    "v21_penalty_reasons",
                ]
            ].to_string(index=False)
            if not supportive.empty
            else "无",
            "",
            "如何进入论文/汇报:",
            "- 主文或图中只建议突出 after-BBB headline CNS-directed mechanism-direction leads。",
            "- pre-BBB primary 但 BBB=No 的候选放入 supplementary 或机制方向讨论，表述为 peripheral/program-modulating clues。",
            "- THRB weak-FDR 候选继续放在 supportive/manual-review 层，不与 NFE2L2 强证据主线混写。",
            "",
            "输出文件:",
            "- 11_swissadme_parsed_and_matched_v21.csv",
            "- 12_bbb_annotation_manual_review_v21.csv",
            "- 13_integrated_after_bbb_v21.csv",
            "- 14_primary_leads_after_bbb_v21.csv",
            "- 15_supportive_after_bbb_v21.csv",
            "- 16_paper_ready_table_after_bbb_v21.csv",
            "- 17_after_bbb_summary_cn.txt",
            "- 18_after_bbb_transition_note_cn.txt",
        ]
    )
    (OUT_DIR / "17_after_bbb_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_transition_note(after: pd.DataFrame) -> None:
    pre_primary = after[after["prebbb_layer_v21"].eq("primary_mechanism_direction_leads")].copy()
    lines = [
        "v2.1 pre-BBB primary 与 after-BBB headline 的关系说明",
        "=" * 62,
        "pre-BBB primary 表示: 程序证据、轴线清晰度、转化适配性与 penalty 后仍较强。",
        "after-BBB headline 表示: 在 pre-BBB primary 基础上，同时获得 SwissADME CNS-directed 支持。",
        "",
        "因此二者不必完全相同:",
        "- pre-BBB primary 但 BBB=No: 保留为 mechanism-direction clue，但降低为 peripheral/program-modulating 解释。",
        "- supportive/THRB weak-FDR 即使 BBB 较好，也仍保持 supportive/manual-review 层级。",
        "- BBB=No 不删除候选，只改变展示层级与 CNS 可达性解释。",
        "",
        "pre-BBB primary 的 after-BBB 迁移:",
        pre_primary[
            [
                "compound",
                "prebbb_layer_v21",
                "manual_final_bbb_layer",
                "after_bbb_final_layer",
                "final_score_prebbb_v21",
                "final_score_after_bbb",
                "after_bbb_interpretation_note",
            ]
        ].to_string(index=False)
        if not pre_primary.empty
        else "无",
    ]
    (OUT_DIR / "18_after_bbb_transition_note_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_after_bbb_v21() -> None:
    integrated, template, swiss = load_inputs()
    matched = match_swissadme(swiss, template)
    annotation = make_bbb_annotation(matched, template)
    after = integrate_after_bbb(integrated, annotation)
    primary = make_primary_after_bbb(after)
    supportive = make_supportive_after_bbb(after)
    paper = make_paper_ready(primary, supportive)

    matched.to_csv(OUT_DIR / "11_swissadme_parsed_and_matched_v21.csv", index=False)
    annotation.to_csv(OUT_DIR / "12_bbb_annotation_manual_review_v21.csv", index=False)
    after.to_csv(OUT_DIR / "13_integrated_after_bbb_v21.csv", index=False)
    primary.to_csv(OUT_DIR / "14_primary_leads_after_bbb_v21.csv", index=False)
    supportive.to_csv(OUT_DIR / "15_supportive_after_bbb_v21.csv", index=False)
    paper.to_csv(OUT_DIR / "16_paper_ready_table_after_bbb_v21.csv", index=False)
    write_summary(matched, annotation, after, primary, supportive)
    write_transition_note(after)

    matched_count = int(matched["matching_method"].ne("unmatched").sum())
    print(f"v2.1 after-BBB 完成: SwissADME 匹配 {matched_count}/{matched.shape[0]}")
    print("BBB 分层计数:")
    print(annotation["manual_final_bbb_layer"].value_counts(dropna=False).to_string())
    print(f"after-BBB headline candidates: {primary.shape[0]}")
    if not primary.empty:
        print(primary[["compound", "manual_final_bbb_layer", "final_score_after_bbb"]].to_string(index=False))


def main() -> None:
    print("开始 v2.1 after-BBB 更新；仅读取本地 SwissADME CSV，不重跑上游、不联网、不做 docking。")
    run_after_bbb_v21()


if __name__ == "__main__":
    main()
