#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
import importlib.util
from datetime import datetime
from pathlib import Path

import pandas as pd

from drugrep_common_v2 import ALLOWED_BBB_LAYERS, DRUGREP_DIR, TUNED_PARAMS


SCRIPT_DIR = Path(__file__).resolve().parent
REQUIRED_COLUMNS = {
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
}


def is_blank(value) -> bool:
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "na", "none"}


def validate_manual(manual: pd.DataFrame) -> None:
    missing = sorted(REQUIRED_COLUMNS - set(manual.columns))
    if missing:
        raise ValueError("06_bbb_annotation_manual_review.csv 缺少列: " + ", ".join(missing))
    bad = manual[~manual["manual_final_bbb_layer"].astype(str).str.strip().isin(ALLOWED_BBB_LAYERS.keys())]
    if not bad.empty:
        raise ValueError("manual_final_bbb_layer 存在非法值；允许值: " + "; ".join(ALLOWED_BBB_LAYERS.keys()))
    must_fill = [
        "swissadme_bbb_prediction",
        "swissadme_gi_absorption",
        "swissadme_boiled_egg_note",
        "swissadme_bioavailability_score",
        "manual_review_note",
        "reviewed_by",
        "reviewed_date",
    ]
    blanks = [f"{c}: {int(manual[c].map(is_blank).sum())}" for c in must_fill if int(manual[c].map(is_blank).sum()) > 0]
    if blanks:
        raise ValueError("以下手工 BBB 列仍有空值或 NA: " + "; ".join(blanks))


def run_swissadme_auto_import() -> None:
    importer_path = SCRIPT_DIR / "06_import_swissadme_csv_v2.py"
    if not importer_path.exists():
        print("未找到 SwissADME 自动导入脚本；将直接读取现有手工 BBB 回填表。")
        return
    spec = importlib.util.spec_from_file_location("import_swissadme_csv_v2", importer_path)
    if spec is None or spec.loader is None:
        print("SwissADME 自动导入脚本加载失败；将直接读取现有手工 BBB 回填表。")
        return
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.auto_import_swissadme_csv(required=False)


def apply_manual_bbb(manual: pd.DataFrame) -> None:
    integrated_path = DRUGREP_DIR / "05_drug_repositioning_integrated_table.csv"
    if not integrated_path.exists():
        raise FileNotFoundError(f"缺少 {integrated_path}")
    df = pd.read_csv(integrated_path)
    m = manual.drop_duplicates("compound")[
        [
            "compound",
            "canonical_smiles",
            "manual_final_bbb_layer",
            "manual_review_note",
            "swissadme_bbb_prediction",
            "swissadme_gi_absorption",
            "swissadme_boiled_egg_note",
            "swissadme_bioavailability_score",
            "reviewed_by",
            "reviewed_date",
        ]
    ]
    merged = df.drop(columns=[c for c in ["canonical_smiles", "manual_final_bbb_layer", "manual_review_note", "swissadme_bbb_prediction", "swissadme_gi_absorption", "swissadme_boiled_egg_note", "swissadme_bioavailability_score", "reviewed_by", "reviewed_date", "BBB_layer", "BBB_score", "BBB_note", "final_score_after_bbb"] if c in df.columns])
    merged = merged.merge(m, on="compound", how="left")
    merged["BBB_layer"] = merged["manual_final_bbb_layer"].fillna("not_reviewed")
    merged["BBB_score"] = merged["BBB_layer"].map(ALLOWED_BBB_LAYERS).fillna(0.35)
    merged["BBB_note"] = merged.apply(
        lambda r: (
            "SwissADME 手工复核: "
            f"BBB={r.get('swissadme_bbb_prediction','')}; "
            f"GI={r.get('swissadme_gi_absorption','')}; "
            f"BOILED-Egg={r.get('swissadme_boiled_egg_note','')}; "
            f"Bioavailability={r.get('swissadme_bioavailability_score','')}; "
            f"Note={r.get('manual_review_note','')}"
        )
        if str(r.get("manual_final_bbb_layer", "")).strip() in ALLOWED_BBB_LAYERS
        else "未进入 v2 SwissADME 手工复核小清单。",
        axis=1,
    )
    merged["final_score_after_bbb"] = (
        0.82 * pd.to_numeric(merged["final_score_prebbb"], errors="coerce").fillna(0)
        + 0.18 * pd.to_numeric(merged["BBB_score"], errors="coerce").fillna(0.35)
    )
    backup_path = DRUGREP_DIR / f"05_drug_repositioning_integrated_table_pre_manual_bbb_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    if integrated_path.exists():
        integrated_path.replace(backup_path)
    merged.to_csv(integrated_path, index=False)

    reviewed_order = {compound: i for i, compound in enumerate(manual["compound"].astype(str).tolist())}
    top = merged[merged["compound"].astype(str).isin(reviewed_order)].copy()
    top["_order"] = top["compound"].astype(str).map(reviewed_order)
    top = top.sort_values("_order").drop(columns=["_order"])
    top.to_csv(DRUGREP_DIR / "07_top_priority_compounds_after_manual_bbb.csv", index=False)
    top.to_csv(DRUGREP_DIR / "05_top_priority_compounds.csv", index=False)
    retained = merged[(merged["final_score_prebbb"] >= 0.30) | merged["include_in_main_axis_pool"].astype(bool)].copy()
    retained.to_csv(DRUGREP_DIR / "05_all_retained_compounds.csv", index=False)
    lines = [
        "v2 SwissADME CSV 自动导入后整合说明",
        "=" * 42,
        f"输入回填文件: {DRUGREP_DIR / '06_bbb_annotation_manual_review.csv'}",
        f"备份 pre-BBB integrated 表: {backup_path}",
        "已基于真实 SwissADME CSV 更新 integrated/top/retained 表；未执行 docking。",
        "",
        "SwissADME 自动 BBB 分层计数:",
        manual["manual_final_bbb_layer"].value_counts(dropna=False).to_string(),
    ]
    (DRUGREP_DIR / "07_after_manual_bbb_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_paper_ready() -> None:
    top = pd.read_csv(DRUGREP_DIR / "05_top_priority_compounds.csv")
    rows = []
    for _, r in top.iterrows():
        rows.append(
            {
                "compound": r["compound"],
                "axis": r["associated_axis"],
                "support_summary": f"{r['support_count']} query set(s); best adj.P={float(r['best_adjusted_p']):.3g}; combined={float(r['best_combined_score']):.3g}",
                "BBB_layer": r.get("BBB_layer", ""),
                "SwissADME_BBB": r.get("swissadme_bbb_prediction", ""),
                "SwissADME_GI_absorption": r.get("swissadme_gi_absorption", ""),
                "SwissADME_bioavailability_score": r.get("swissadme_bioavailability_score", ""),
                "mechanism_theme": r.get("indicative_mechanism", ""),
                "recommended_role": "program-guided mechanism-direction clue" if r.get("final_priority_tier", "") != "manual_review_weak_fdr" else "supportive/manual-review program clue",
                "evidence_level": r.get("final_priority_tier", ""),
                "caution_note": "Exploratory clue only; requires independent pharmacological and experimental follow-up before interpretation beyond hypothesis generation.",
            }
        )
    paper = pd.DataFrame(rows)
    paper.to_csv(DRUGREP_DIR / "07_drug_repositioning_paper_ready_table.csv", index=False)
    paper.to_csv(DRUGREP_DIR / "07_drug_repositioning_paper_ready_table_after_manual_bbb.csv", index=False)


def make_master_summary() -> None:
    top = pd.read_csv(DRUGREP_DIR / "05_top_priority_compounds.csv")
    retained = pd.read_csv(DRUGREP_DIR / "05_all_retained_compounds.csv")
    manual = pd.read_csv(DRUGREP_DIR / "06_bbb_annotation_manual_review.csv")
    parsed_path = DRUGREP_DIR / "manual_bbb_review" / "06_swissadme_parsed_and_matched.csv"
    parsed = pd.read_csv(parsed_path) if parsed_path.exists() else pd.DataFrame()
    import_summary_path = DRUGREP_DIR / "manual_bbb_review" / "06_bbb_auto_import_summary_cn.txt"
    import_summary = import_summary_path.read_text(encoding="utf-8") if import_summary_path.exists() else ""

    adopted_line = "实际采用的 SwissADME CSV: 未记录"
    for line in import_summary.splitlines():
        if line.startswith("实际采用的 SwissADME CSV:"):
            adopted_line = line
            break

    layer_counts = manual["manual_final_bbb_layer"].value_counts(dropna=False)
    method_counts = parsed["match_method"].value_counts(dropna=False) if not parsed.empty and "match_method" in parsed.columns else pd.Series(dtype=int)
    matched_count = int((parsed["match_method"] != "unmatched").sum()) if not parsed.empty and "match_method" in parsed.columns else 0
    downgraded = manual[manual["manual_final_bbb_layer"] != "CNS-directed candidates"].copy()

    top_cols = [
        "compound",
        "associated_axis",
        "final_score_prebbb",
        "final_score_after_bbb",
        "BBB_layer",
        "swissadme_bbb_prediction",
        "swissadme_gi_absorption",
        "final_priority_tier",
    ]
    available_top_cols = [c for c in top_cols if c in top.columns]
    lines = [
        "药物重定位 v2 SwissADME 自动 BBB 后主汇总",
        "=" * 56,
        "本轮仅导入 SwissADME CSV 并继续 v2 后续整合；未重跑 pySCENIC、CellOracle、扰动、外部验证、robustness、GO/KEGG、Enrichr 或 compound aggregation。",
        "本轮未执行 docking，未调用任何 docking 脚本。",
        adopted_line,
        "",
        "匹配结果:",
        f"- v2 主表候选匹配: {matched_count} / {manual.shape[0]}",
        method_counts.to_string() if not method_counts.empty else "无匹配方法统计",
        "",
        "真实 SwissADME BBB 分层计数:",
        layer_counts.to_string(),
        "",
        "更新后主表候选:",
        top[available_top_cols].to_string(index=False) if available_top_cols else top.to_string(index=False),
        "",
        "因 SwissADME 结果未进入 CNS-directed 的候选:",
        downgraded[["compound", "manual_final_bbb_layer", "swissadme_bbb_prediction", "swissadme_gi_absorption", "manual_review_note"]].to_string(index=False) if not downgraded.empty else "无",
        "",
        "输出文件:",
        f"- {DRUGREP_DIR / '06_bbb_annotation_manual_review.csv'}",
        f"- {DRUGREP_DIR / 'manual_bbb_review' / '06_swissadme_parsed_and_matched.csv'}",
        f"- {DRUGREP_DIR / 'manual_bbb_review' / '06_bbb_layered_candidates_manual_review.csv'}",
        f"- {DRUGREP_DIR / 'manual_bbb_review' / '06_bbb_auto_import_summary_cn.txt'}",
        f"- {DRUGREP_DIR / '05_top_priority_compounds.csv'}",
        f"- {DRUGREP_DIR / '05_drug_repositioning_integrated_table.csv'}",
        f"- {DRUGREP_DIR / '05_all_retained_compounds.csv'}",
        f"- {DRUGREP_DIR / '07_drug_repositioning_paper_ready_table.csv'}",
        f"- {DRUGREP_DIR / '07_drug_repositioning_master_summary_cn.txt'}",
        "",
        f"retained 补充候选总数: {retained.shape[0]}；BBB=No 的候选未被删除，仅改变解释分层。",
    ]
    (DRUGREP_DIR / "07_drug_repositioning_master_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    print("读取 SwissADME CSV 自动回填 v2 BBB 结果并继续整合。本脚本不执行 docking。")
    run_swissadme_auto_import()
    manual_path = DRUGREP_DIR / "06_bbb_annotation_manual_review.csv"
    if not manual_path.exists():
        raise FileNotFoundError(f"缺少 {manual_path}")
    manual = pd.read_csv(manual_path)
    validate_manual(manual)
    apply_manual_bbb(manual)
    make_paper_ready()
    make_master_summary()
    print("v2 SwissADME 自动 BBB 后续整合完成。未执行 docking。")


if __name__ == "__main__":
    main()
