# -*- coding: utf-8 -*-
"""
职责：校验 01_metadata/dataset_registry.csv 与 sample_metadata.csv 的列名、GSE_id 及 Phase 2 数据集字段。

输入（默认路径，相对于项目根目录）：
  - 01_metadata/dataset_registry.csv
  - 01_metadata/sample_metadata.csv

输出：中文日志打印至 stdout。

说明：不包含任何数据下载或矩阵读取逻辑。
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# --- 必需列（含 Phase 2：平台、模态、样本数、分组、对照、原始文件类型、分析定位）---
REQUIRED_DATASET_COLUMNS = [
    "GSE_id",
    "disease",
    "subtype",
    "assay",
    "tissue",
    "species",
    "platform",
    "modality_type",
    "n_samples",
    "grouping_summary",
    "has_control",
    "raw_file_types",
    "suited_for_analysis",
    "role",
    "discovery_or_validation",
    "notes",
]

REQUIRED_SAMPLE_COLUMNS = [
    "sample_id",
    "donor_id",
    "GSE_id",
    "disease",
    "subtype",
    "group",
    "tissue",
    "assay",
    "platform",
    "source_note",
]

# Phase 2：已填写的枚举取值（其余视为「待 GEO 核对」类占位，仅给提示不失败）
MODALITY_ENUM = {"scRNA", "snRNA", "multiome", "mixed", "unclear"}
HAS_CONTROL_ENUM = {"yes", "no", "partial", "unclear"}
SUITED_ENUM = {"discovery", "validation", "supplemental_validation", "dual"}

PHASE2_ENUM_COLUMNS = {
    "modality_type": MODALITY_ENUM,
    "has_control": HAS_CONTROL_ENUM,
    "suited_for_analysis": SUITED_ENUM,
}

# 非枚举、但须在 Phase 2 从 GEO/原文补齐的文本字段（占位则列入【提示】）
PHASE2_TEXT_COLUMNS = [
    "platform",
    "n_samples",
    "grouping_summary",
    "raw_file_types",
]


def _project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _missing_columns(df: pd.DataFrame, required: list[str]) -> list[str]:
    return [c for c in required if c not in df.columns]


def _blank_gse_mask(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip()
    return s.eq("") | s.str.lower().eq("nan")


def _is_placeholder(val: str) -> bool:
    """尚未从 GEO 固化的取值：空、TBD、待补充、待核对、待 GEO 等。"""
    t = val.strip()
    if t == "" or t.lower() in ("nan", "tbd", "na"):
        return True
    if t in ("待补充", "待核对", "待GEO核对", "待GEO确认"):
        return True
    if t.startswith("待补充") or t.startswith("待GEO") or t.startswith("待核对"):
        return True
    return False


def _pure_enum_for_column(col: str, val: str) -> str | None:
    """若该列取值恰好等于该列允许枚举之一，则返回该枚举；否则返回 None。"""
    allowed = PHASE2_ENUM_COLUMNS.get(col)
    if allowed is None:
        return None
    t = val.strip()
    if t in allowed:
        return t
    return None


def main() -> None:
    root = _project_root()
    path_ds = root / "01_metadata" / "dataset_registry.csv"
    path_sm = root / "01_metadata" / "sample_metadata.csv"

    print("【元数据校验】项目根目录：", root)
    ok = True

    if not path_ds.is_file():
        print("【失败】未找到文件：", path_ds)
        sys.exit(1)
    if not path_sm.is_file():
        print("【失败】未找到文件：", path_sm)
        sys.exit(1)

    df_ds = _read_csv(path_ds)
    miss_ds = _missing_columns(df_ds, REQUIRED_DATASET_COLUMNS)
    if miss_ds:
        print("【失败】dataset_registry.csv 缺少列：", miss_ds)
        ok = False
    else:
        print("【通过】dataset_registry.csv 必需列齐全（含 Phase 2 字段）。")

    df_sm = _read_csv(path_sm)
    miss_sm = _missing_columns(df_sm, REQUIRED_SAMPLE_COLUMNS)
    if miss_sm:
        print("【失败】sample_metadata.csv 缺少列：", miss_sm)
        ok = False
    else:
        print("【通过】sample_metadata.csv 必需列齐全。")

    if not ok:
        sys.exit(1)

    # dataset：GSE_id 非空
    bad_ds = df_ds.loc[_blank_gse_mask(df_ds["GSE_id"])]
    if len(bad_ds) > 0:
        print("【失败】dataset_registry.csv 存在空 GSE_id 行，行号（从 1 计数据行）：")
        print((bad_ds.index + 2).tolist())
        ok = False
    else:
        print("【通过】dataset_registry.csv 所有行 GSE_id 非空。")

    # Phase 2：枚举列 — 允许占位；自由文本且非占位则判为需修正
    phase2_pending: list[str] = []
    phase2_bad: list[str] = []
    for col in PHASE2_ENUM_COLUMNS:
        for idx, row in df_ds.iterrows():
            gse = str(row["GSE_id"]).strip()
            raw = str(row[col]).strip()
            if _pure_enum_for_column(col, raw) is not None:
                continue
            if _is_placeholder(raw):
                phase2_pending.append(f"{gse}.{col}")
                continue
            if raw.startswith("待补充") or raw.startswith("待GEO"):
                phase2_pending.append(f"{gse}.{col}")
            else:
                phase2_bad.append(
                    f"{gse}.{col}（请改为该列允许枚举值，或使用「待补充」类占位；参见 dataset_registry 表头说明）"
                )

    for col in PHASE2_TEXT_COLUMNS:
        for idx, row in df_ds.iterrows():
            gse = str(row["GSE_id"]).strip()
            raw = str(row[col]).strip()
            if _is_placeholder(raw):
                phase2_pending.append(f"{gse}.{col}")

    if phase2_bad:
        print("【失败】Phase 2 字段取值不符合约定：")
        for item in phase2_bad:
            print("   - ", item)
        ok = False
    else:
        print("【通过】Phase 2 枚举字段未发现非法纯取值。")

    if phase2_pending:
        print("【提示】以下字段仍为占位或待 GEO 核对（Phase 2 正常状态，补齐后应改为明确枚举或具体文字）：")
        by_gse: dict[str, list[str]] = defaultdict(list)
        for key in phase2_pending:
            gse, col = key.split(".", 1)
            by_gse[gse].append(col)
        for gse in sorted(by_gse.keys()):
            cols = sorted(set(by_gse[gse]))
            print(f"   - {gse}: {', '.join(cols)}")

    # suited_for_analysis 与 discovery_or_validation 一致性
    for idx, row in df_ds.iterrows():
        gse = str(row["GSE_id"]).strip()
        suited = str(row["suited_for_analysis"]).strip()
        dorv = str(row["discovery_or_validation"]).strip()
        if suited in SUITED_ENUM:
            if suited == "discovery" and dorv != "discovery":
                print(f"【失败】{gse}: suited_for_analysis=discovery 但 discovery_or_validation 应为 discovery。")
                ok = False
            if suited in ("validation", "supplemental_validation") and dorv != "validation":
                print(
                    f"【失败】{gse}: suited_for_analysis={suited} 时 discovery_or_validation 应为 validation。"
                )
                ok = False
            if suited == "dual" and dorv not in ("discovery", "validation"):
                print(f"【失败】{gse}: suited_for_analysis=dual 时 discovery_or_validation 取值异常。")
                ok = False

    # sample：若无数据行
    if len(df_sm) == 0:
        print("【警告】sample_metadata.csv 当前无任何样本行；请从 GEO 补充后再跑校验。")
    else:
        bad_sm_gse = df_sm.loc[_blank_gse_mask(df_sm["GSE_id"])]
        if len(bad_sm_gse) > 0:
            print("【失败】sample_metadata.csv 存在空 GSE_id 行，行号（从 1 计数据行）：")
            print((bad_sm_gse.index + 2).tolist())
            ok = False
        else:
            print("【通过】sample_metadata.csv 已填写行的 GSE_id 均非空。")

        bad_sm_sid = df_sm.loc[_blank_gse_mask(df_sm["sample_id"])]
        if len(bad_sm_sid) > 0:
            print("【失败】sample_metadata.csv 存在空 sample_id 行，行号（从 1 计数据行）：")
            print((bad_sm_sid.index + 2).tolist())
            ok = False

        bad_sm_donor = df_sm.loc[_blank_gse_mask(df_sm["donor_id"])]
        if len(bad_sm_donor) > 0:
            print("【失败】sample_metadata.csv 存在空 donor_id 行，行号（从 1 计数据行）：")
            print((bad_sm_donor.index + 2).tolist())
            ok = False

        sm_sid = df_sm["sample_id"].astype(str).str.strip()
        bad_gsm = ~sm_sid.str.match(r"^GSM\d+$", na=False)
        if bad_gsm.any():
            print("【失败】sample_id 应为 GSM 数字编号格式（GSMxxxxxxx），以下行不符合：")
            print((df_sm.loc[bad_gsm, "sample_id"].astype(str).tolist()))
            ok = False

        reg_ids = set(df_ds["GSE_id"].astype(str).str.strip())
        sm_ids = df_sm["GSE_id"].astype(str).str.strip()
        unknown = sorted(set(sm_ids) - reg_ids)
        if unknown:
            print("【失败】sample_metadata 中的 GSE_id 在 dataset_registry 中不存在：", unknown)
            ok = False
        else:
            print("【通过】sample_metadata 中的 GSE_id 均在 dataset_registry 中有登记。")

    if ok:
        print("【汇总】校验完成：通过。")
        sys.exit(0)
    else:
        print("【汇总】校验完成：存在失败项，请修正后重试。")
        sys.exit(1)


if __name__ == "__main__":
    main()
