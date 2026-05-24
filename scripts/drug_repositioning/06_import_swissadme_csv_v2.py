#!/usr/bin/env python3
from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd

from drugrep_common_v2 import DRUGREP_DIR, MANUAL_REVIEW_DIR, ROOT, ensure_dirs


MANUAL_COLUMNS = [
    "compound",
    "canonical_smiles",
    "swissadme_bbb_prediction",
    "swissadme_gi_absorption",
    "swissadme_boiled_egg_note",
    "swissadme_bioavailability_score",
    "manual_final_bbb_layer",
    "manual_review_note",
    "reviewed_by",
    "reviewed_date",
]

PARSED_MATCHED_COLUMNS = [
    "swissadme_row_id",
    "swissadme_molecule_id",
    "compound",
    "canonical_smiles",
    "match_method",
    "matched_confidence",
    "BBB permeant",
    "GI absorption",
    "Pgp substrate",
    "Bioavailability Score",
    "PAINS",
    "Brenk",
    "Leadlikeness",
    "Synthetic accessibility",
]

SMILES_ALIASES = [
    "Canonical SMILES",
    "CanonicalSMILES",
    "SMILES",
    "Isomeric SMILES",
    "IsomericSMILES",
    "canonical_smiles",
    "isomeric_smiles",
]
NAME_ALIASES = [
    "compound",
    "Compound",
    "name",
    "Name",
    "Title",
    "Molecule name",
    "Molecule Name",
    "Molecule",
]
FIELD_ALIASES = {
    "molecule": ["Molecule", "molecule", "Molecule ID", "Molecule name", "Name"],
    "bbb": ["BBB permeant", "BBB", "BBB permeation", "BBB prediction"],
    "gi": ["GI absorption", "Gastrointestinal absorption"],
    "pgp": ["Pgp substrate", "P-gp substrate", "Pgp", "P-gp"],
    "bioavailability": ["Bioavailability Score", "Bioavailability score"],
    "pains": ["PAINS #alerts", "PAINS alerts", "PAINS"],
    "brenk": ["Brenk #alerts", "Brenk alerts", "Brenk"],
    "leadlikeness": ["Leadlikeness #violations", "Leadlikeness violations", "Leadlikeness"],
    "synthetic_accessibility": [
        "Synthetic Accessibility",
        "Synthetic accessibility",
        "Synthetic Accessibility Score",
    ],
}

GENERIC_MOLECULE_RE = re.compile(r"^molecule\s*\d+$", re.IGNORECASE)


@dataclass(frozen=True)
class SwissCandidate:
    path: Path
    reason_source: str
    in_v2_manual_inputs: bool
    mtime: float


def is_blank(value: Any) -> bool:
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "na", "none", "null"}


def norm_col(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def norm_smiles(value: Any) -> str:
    if is_blank(value):
        return ""
    return re.sub(r"\s+", "", str(value).strip())


def looks_like_smiles(value: Any) -> bool:
    text = norm_smiles(value)
    if not text or GENERIC_MOLECULE_RE.match(str(value).strip()):
        return False
    if re.search(r"\s", str(value).strip()):
        return False
    has_atom = bool(re.search(r"\[[^\]]+\]|Br|Cl|[BCNOFPSIbcno]", text))
    has_structure_token = bool(re.search(r"[=#()\[\]/\\@+\-0-9]", text))
    return has_atom and has_structure_token


def norm_name(value: Any) -> str:
    if is_blank(value):
        return ""
    text = str(value).strip().lower()
    text = re.sub(r"[\u2010-\u2015]", "-", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def clean_value(value: Any) -> str:
    return "NA" if is_blank(value) else str(value).strip()


def find_column(df: pd.DataFrame, aliases: list[str]) -> str | None:
    normalized = {norm_col(col): col for col in df.columns}
    for alias in aliases:
        col = normalized.get(norm_col(alias))
        if col is not None:
            return col
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


def is_generated_swissadme_output(path: Path) -> bool:
    lower_parts = [part.lower() for part in path.parts]
    name = path.name.lower()
    if "manual_bbb_review" in lower_parts:
        return True
    return any(token in name for token in ["parsed_and_matched", "layered_candidates", "auto_import"])


def list_swissadme_candidates() -> list[SwissCandidate]:
    manual_inputs = DRUGREP_DIR / "manual_inputs"
    preferred = [
        (manual_inputs / "swissadme.csv", "指定优先路径 1"),
        (manual_inputs / "swissadme (1).csv", "指定优先路径 2"),
        (DRUGREP_DIR / "swissadme.csv", "指定优先路径 3"),
        (ROOT / "swissadme.csv", "指定优先路径 4"),
    ]
    seen: set[str] = set()
    candidates: list[SwissCandidate] = []

    def add(path: Path, reason_source: str) -> None:
        if not path.exists() or not path.is_file():
            return
        key = str(path.resolve()).lower()
        if key in seen:
            return
        seen.add(key)
        candidates.append(
            SwissCandidate(
                path=path,
                reason_source=reason_source,
                in_v2_manual_inputs=path.parent.resolve() == manual_inputs.resolve(),
                mtime=path.stat().st_mtime,
            )
        )

    for path, reason_source in preferred:
        add(path, reason_source)

    for path in ROOT.rglob("*.csv"):
        if "swissadme" not in path.name.lower():
            continue
        if is_generated_swissadme_output(path):
            continue
        add(path, "项目目录递归搜索: 文件名包含 swissadme")
    return candidates


def choose_swissadme_csv(required: bool = True) -> tuple[Path | None, str, list[SwissCandidate]]:
    candidates = list_swissadme_candidates()
    if not candidates:
        if required:
            raise FileNotFoundError("未找到 SwissADME CSV 文件。")
        return None, "未找到 SwissADME CSV，保留现有 06_bbb_annotation_manual_review.csv。", []

    manual_candidates = [c for c in candidates if c.in_v2_manual_inputs]
    if manual_candidates:
        chosen = sorted(manual_candidates, key=lambda c: c.mtime, reverse=True)[0]
        reason = (
            "选择原因: 文件位于 drug_repositioning\\manual_inputs 优先目录；"
            f"该目录候选 {len(manual_candidates)} 个，按最新修改时间选择。"
        )
        return chosen.path, reason, candidates

    chosen = sorted(candidates, key=lambda c: c.mtime, reverse=True)[0]
    reason = "选择原因: v2 manual_inputs 目录下没有候选文件；按项目内 swissadme*.csv 最新修改时间选择。"
    return chosen.path, reason, candidates


def try_rdkit_canonical(value: Any) -> str:
    if is_blank(value):
        return ""
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return ""
    mol = Chem.MolFromSmiles(str(value).strip())
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def build_unique_map(values: list[tuple[str, int]]) -> dict[str, int]:
    bucket: dict[str, list[int]] = {}
    for value, idx in values:
        if not value:
            continue
        bucket.setdefault(value, []).append(idx)
    return {value: idxs[0] for value, idxs in bucket.items() if len(set(idxs)) == 1}


def load_reference_table() -> pd.DataFrame:
    copy_path = DRUGREP_DIR / "06_swissadme_copy_paste_with_names.tsv"
    if not copy_path.exists():
        raise FileNotFoundError(f"缺少 {copy_path}")
    ref = pd.read_csv(copy_path, sep="\t")
    required = {"compound", "canonical_smiles"}
    missing = sorted(required - set(ref.columns))
    if missing:
        raise ValueError("06_swissadme_copy_paste_with_names.tsv 缺少列: " + ", ".join(missing))

    top_path = DRUGREP_DIR / "05_top_priority_compounds.csv"
    if top_path.exists():
        top = pd.read_csv(top_path)
        if {"compound", "normalized_compound_name"}.issubset(top.columns):
            ref = ref.merge(top[["compound", "normalized_compound_name"]], on="compound", how="left")
    if "normalized_compound_name" not in ref.columns:
        ref["normalized_compound_name"] = ref["compound"].map(norm_name)
    return ref.reset_index(drop=True)


def is_reliable_name(value: Any) -> bool:
    if is_blank(value):
        return False
    text = str(value).strip()
    return GENERIC_MOLECULE_RE.match(text) is None


def make_reference_maps(ref: pd.DataFrame) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    smiles_values: list[tuple[str, int]] = []
    rdkit_values: list[tuple[str, int]] = []
    name_values: list[tuple[str, int]] = []
    for idx, row in ref.iterrows():
        for col in ["canonical_smiles", "isomeric_smiles"]:
            if col in ref.columns:
                smiles = norm_smiles(row.get(col, ""))
                if smiles:
                    smiles_values.append((smiles, idx))
                    rdkit_smiles = try_rdkit_canonical(smiles)
                    if rdkit_smiles:
                        rdkit_values.append((rdkit_smiles, idx))
        for col in ["compound", "normalized_compound_name"]:
            name = norm_name(row.get(col, ""))
            if name:
                name_values.append((name, idx))
    return build_unique_map(smiles_values), build_unique_map(rdkit_values), build_unique_map(name_values)


def swiss_value(row: pd.Series, col: str | None) -> str:
    return "NA" if col is None else clean_value(row.get(col, "NA"))


def numeric_or_none(value: Any) -> float | None:
    if is_blank(value):
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


def yes_no(value: Any) -> str:
    text = clean_value(value).lower()
    if text in {"yes", "y", "true", "1"}:
        return "Yes"
    if text in {"no", "n", "false", "0"}:
        return "No"
    return clean_value(value)


def classify_bbb(row: dict[str, Any]) -> tuple[str, str]:
    bbb = yes_no(row.get("BBB permeant", "NA"))
    gi = clean_value(row.get("GI absorption", "NA"))
    pgp = yes_no(row.get("Pgp substrate", "NA"))
    brenk = numeric_or_none(row.get("Brenk", "NA"))
    pains = numeric_or_none(row.get("PAINS", "NA"))
    bioavailability = clean_value(row.get("Bioavailability Score", "NA"))

    reasons: list[str] = []
    if bbb == "Yes":
        if gi.lower() == "high" and pgp == "No" and (brenk is None or brenk <= 1) and (pains is None or pains == 0):
            layer = "CNS-directed candidates"
            reasons.append("BBB=Yes、GI=High、Pgp=No，且 Brenk/PAINS 未见明显异常。")
        else:
            layer = "possible CNS-directed candidates"
            reasons.append("BBB=Yes，但存在限制因素，保守降为 possible CNS-directed。")
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
        reasons.append("BBB=No；不删除候选，改按外围/程序调节线索解释。")
    else:
        layer = "peripheral/program-modulating candidates"
        reasons.append("BBB 结果缺失或不可识别，按保守规则归入外围/程序调节候选。")

    if is_blank(row.get("BBB permeant", "")):
        reasons.append("BBB permeant 缺失。")
    if is_blank(row.get("GI absorption", "")):
        reasons.append("GI absorption 缺失。")
    if is_blank(row.get("Pgp substrate", "")):
        reasons.append("Pgp substrate 缺失。")
    if bioavailability != "NA":
        reasons.append(f"Bioavailability Score={bioavailability}。")
    return layer, " ".join(reasons)


def match_swissadme_rows(swiss: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    molecule_col = find_column(swiss, FIELD_ALIASES["molecule"])
    smiles_col = find_column(swiss, SMILES_ALIASES)
    name_col = find_column(swiss, NAME_ALIASES)
    bbb_col = find_column(swiss, FIELD_ALIASES["bbb"])
    gi_col = find_column(swiss, FIELD_ALIASES["gi"])
    pgp_col = find_column(swiss, FIELD_ALIASES["pgp"])
    bio_col = find_column(swiss, FIELD_ALIASES["bioavailability"])
    pains_col = find_column(swiss, FIELD_ALIASES["pains"])
    brenk_col = find_column(swiss, FIELD_ALIASES["brenk"])
    lead_col = find_column(swiss, FIELD_ALIASES["leadlikeness"])
    synth_col = find_column(swiss, FIELD_ALIASES["synthetic_accessibility"])

    smiles_map, rdkit_map, name_map = make_reference_maps(ref)
    used_ref: set[int] = set()
    matched: list[dict[str, Any]] = []

    for swiss_idx, row in swiss.reset_index(drop=True).iterrows():
        matched_idx: int | None = None
        match_method = "unmatched"
        confidence = 0.0

        smiles_candidates: list[str] = []
        if smiles_col:
            smiles_candidates.append(norm_smiles(row.get(smiles_col, "")))
        if molecule_col and molecule_col != smiles_col:
            molecule_value = row.get(molecule_col, "")
            molecule_smiles = norm_smiles(molecule_value)
            if molecule_smiles and (molecule_smiles in smiles_map or looks_like_smiles(molecule_value)):
                smiles_candidates.append(molecule_smiles)
        smiles_candidates = list(dict.fromkeys([s for s in smiles_candidates if s]))

        for swiss_smiles in smiles_candidates:
            if swiss_smiles in smiles_map and smiles_map[swiss_smiles] not in used_ref:
                matched_idx = smiles_map[swiss_smiles]
                match_method = "SMILES"
                confidence = 1.0
                break

        if matched_idx is None and rdkit_map:
            for swiss_smiles in smiles_candidates:
                swiss_rdkit = try_rdkit_canonical(swiss_smiles)
                if swiss_rdkit and swiss_rdkit in rdkit_map and rdkit_map[swiss_rdkit] not in used_ref:
                    matched_idx = rdkit_map[swiss_rdkit]
                    match_method = "SMILES_rdkit"
                    confidence = 1.0
                    break

        if matched_idx is None and name_col:
            swiss_name = row.get(name_col, "")
            if is_reliable_name(swiss_name):
                key = norm_name(swiss_name)
                if key in name_map and name_map[key] not in used_ref:
                    matched_idx = name_map[key]
                    match_method = "name"
                    confidence = 0.95

        if matched_idx is None and swiss.shape[0] == ref.shape[0] and swiss_idx < ref.shape[0] and swiss_idx not in used_ref:
            matched_idx = int(swiss_idx)
            match_method = "import_order"
            confidence = 0.90

        molecule_id = swiss_value(row, molecule_col)
        if matched_idx is not None:
            used_ref.add(matched_idx)
            compound = clean_value(ref.loc[matched_idx, "compound"])
            canonical_smiles = clean_value(ref.loc[matched_idx, "canonical_smiles"])
        else:
            compound = "NA"
            canonical_smiles = "NA"

        record = {
            "swissadme_row_id": int(swiss_idx) + 1,
            "swissadme_molecule_id": molecule_id,
            "compound": compound,
            "canonical_smiles": canonical_smiles,
            "match_method": match_method,
            "matched_confidence": confidence,
            "BBB permeant": yes_no(swiss_value(row, bbb_col)),
            "GI absorption": swiss_value(row, gi_col),
            "Pgp substrate": yes_no(swiss_value(row, pgp_col)),
            "Bioavailability Score": swiss_value(row, bio_col),
            "PAINS": swiss_value(row, pains_col),
            "Brenk": swiss_value(row, brenk_col),
            "Leadlikeness": swiss_value(row, lead_col),
            "Synthetic accessibility": swiss_value(row, synth_col),
        }
        layer, note = classify_bbb(record)
        record["manual_final_bbb_layer"] = layer
        record["manual_review_note"] = note
        record["swissadme_boiled_egg_note"] = (
            f"BBB permeant={record['BBB permeant']}; GI absorption={record['GI absorption']}; "
            f"Pgp substrate={record['Pgp substrate']}; match_method={match_method}"
        )
        matched.append(record)

    return pd.DataFrame(matched)


def make_manual_review(ref: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    today = date.today().isoformat()
    by_compound = matched.drop_duplicates("compound").set_index("compound")
    rows: list[dict[str, Any]] = []
    for _, ref_row in ref.iterrows():
        compound = clean_value(ref_row["compound"])
        canonical_smiles = clean_value(ref_row["canonical_smiles"])
        if compound in by_compound.index:
            m = by_compound.loc[compound]
            rows.append(
                {
                    "compound": compound,
                    "canonical_smiles": canonical_smiles,
                    "swissadme_bbb_prediction": clean_value(m.get("BBB permeant", "NA")),
                    "swissadme_gi_absorption": clean_value(m.get("GI absorption", "NA")),
                    "swissadme_boiled_egg_note": clean_value(m.get("swissadme_boiled_egg_note", "NA")),
                    "swissadme_bioavailability_score": clean_value(m.get("Bioavailability Score", "NA")),
                    "manual_final_bbb_layer": clean_value(m.get("manual_final_bbb_layer", "NA")),
                    "manual_review_note": clean_value(m.get("manual_review_note", "NA")),
                    "reviewed_by": "swissadme_csv_auto_import",
                    "reviewed_date": today,
                }
            )
        else:
            rows.append(
                {
                    "compound": compound,
                    "canonical_smiles": canonical_smiles,
                    "swissadme_bbb_prediction": "NA",
                    "swissadme_gi_absorption": "NA",
                    "swissadme_boiled_egg_note": "未能匹配 SwissADME 结果。",
                    "swissadme_bioavailability_score": "NA",
                    "manual_final_bbb_layer": "peripheral/program-modulating candidates",
                    "manual_review_note": "SwissADME 行未匹配；按保守规则暂归入外围/程序调节候选，需人工核查。",
                    "reviewed_by": "swissadme_csv_auto_import",
                    "reviewed_date": today,
                }
            )
    return pd.DataFrame(rows, columns=MANUAL_COLUMNS)


def make_summary(
    chosen_path: Path,
    chosen_reason: str,
    candidates: list[SwissCandidate],
    matched: pd.DataFrame,
    manual: pd.DataFrame,
) -> str:
    method_counts = matched["match_method"].value_counts(dropna=False)
    layer_counts = manual["manual_final_bbb_layer"].value_counts(dropna=False)
    matched_count = int((matched["match_method"] != "unmatched").sum())
    downgraded = manual[manual["manual_final_bbb_layer"] != "CNS-directed candidates"].copy()
    candidate_lines = [
        f"- {c.path} | 来源={c.reason_source} | manual_inputs={c.in_v2_manual_inputs} | mtime={pd.Timestamp(c.mtime, unit='s')}"
        for c in sorted(candidates, key=lambda x: (not x.in_v2_manual_inputs, -x.mtime))
    ]
    downgraded_lines = [
        f"- {r['compound']}: {r['manual_final_bbb_layer']}；{r['manual_review_note']}"
        for _, r in downgraded.iterrows()
    ]
    match_lines = [
        f"- {r['compound']}: {r['match_method']}，confidence={r['matched_confidence']}"
        for _, r in matched.iterrows()
    ]
    if (matched["match_method"] == "import_order").any():
        order_note = "- 存在导入顺序匹配；该回退仅在 SwissADME 行数与 06_swissadme_copy_paste_with_names.tsv 完全一致时启用，且不会改变原始顺序。"
    else:
        order_note = "- 本次未使用导入顺序匹配。"
    return "\n".join(
        [
            "SwissADME CSV 自动导入与 BBB 分层日志",
            "=" * 52,
            f"实际采用的 SwissADME CSV: {chosen_path}",
            chosen_reason,
            "",
            "发现的候选文件:",
            "\n".join(candidate_lines) if candidate_lines else "无",
            "",
            f"SwissADME 行数: {matched.shape[0]}",
            f"v2 主表候选匹配: {matched_count} / {manual.shape[0]}",
            "",
            "匹配方法计数:",
            method_counts.to_string(),
            "",
            "逐药匹配方法:",
            "\n".join(match_lines) if match_lines else "无",
            "",
            "BBB 自动分层计数:",
            layer_counts.to_string(),
            "",
            "因 SwissADME 结果未进入 CNS-directed 的候选:",
            "\n".join(downgraded_lines) if downgraded_lines else "无",
            "",
            "匹配不确定性说明:",
            "- 已优先尝试 SMILES 文本精确匹配；当 Molecule 列包含具体 SMILES 时也纳入 SMILES 匹配。",
            "- 如本机有 RDKit，也会尝试 RDKit 规范化 SMILES 匹配；当前无 RDKit 时仅做文本级匹配。",
            "- 若 SMILES/名称无法可靠命中，才回退到导入顺序匹配。",
            order_note,
            "- 如后续替换 CSV，应重新运行本脚本并核查 parsed_and_matched 表。",
            "",
            "输出文件:",
            f"- {DRUGREP_DIR / '06_bbb_annotation_manual_review.csv'}",
            f"- {MANUAL_REVIEW_DIR / '06_swissadme_parsed_and_matched.csv'}",
            f"- {MANUAL_REVIEW_DIR / '06_bbb_layered_candidates_manual_review.csv'}",
            f"- {MANUAL_REVIEW_DIR / '06_bbb_auto_import_summary_cn.txt'}",
            "",
            "流程边界:",
            "- 本步骤只读取 SwissADME CSV 并回填 BBB 分层。",
            "- 未重跑 Enrichr、compound aggregation、pySCENIC、CellOracle、扰动、外部验证、robustness 或 GO/KEGG。",
            "- 未执行 docking，未调用任何 docking 脚本。",
        ]
    )


def auto_import_swissadme_csv(required: bool = True) -> pd.DataFrame | None:
    ensure_dirs()
    chosen_path, chosen_reason, candidates = choose_swissadme_csv(required=required)
    if chosen_path is None:
        print(chosen_reason)
        return None

    print(f"实际采用 SwissADME CSV: {chosen_path}")
    print(chosen_reason)
    ref = load_reference_table()
    swiss = read_csv_flexible(chosen_path)
    matched = match_swissadme_rows(swiss, ref)
    manual = make_manual_review(ref, matched)

    MANUAL_REVIEW_DIR.mkdir(parents=True, exist_ok=True)
    matched[PARSED_MATCHED_COLUMNS].to_csv(MANUAL_REVIEW_DIR / "06_swissadme_parsed_and_matched.csv", index=False)
    layered_cols = PARSED_MATCHED_COLUMNS + [
        "manual_final_bbb_layer",
        "manual_review_note",
        "swissadme_boiled_egg_note",
    ]
    matched[layered_cols].to_csv(MANUAL_REVIEW_DIR / "06_bbb_layered_candidates_manual_review.csv", index=False)
    manual.to_csv(DRUGREP_DIR / "06_bbb_annotation_manual_review.csv", index=False)

    summary = make_summary(chosen_path, chosen_reason, candidates, matched, manual)
    (MANUAL_REVIEW_DIR / "06_bbb_auto_import_summary_cn.txt").write_text(summary + "\n", encoding="utf-8")

    matched_count = int((matched["match_method"] != "unmatched").sum())
    print(f"SwissADME 自动导入完成: 匹配 {matched_count} / {manual.shape[0]}。")
    print("BBB 分层计数:")
    print(manual["manual_final_bbb_layer"].value_counts(dropna=False).to_string())
    return manual


def main() -> None:
    print("开始自动读取 SwissADME CSV 并回填 v2 BBB 分层；本脚本不执行 docking。")
    auto_import_swissadme_csv(required=True)
    print("SwissADME CSV 自动导入完成；未执行 docking。")


if __name__ == "__main__":
    main()
