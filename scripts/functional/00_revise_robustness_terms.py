#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd

from functional_common import FUNCTIONAL_DIR, ROBUSTNESS_DIR, ROOT, ensure_functional_dir


def safe_read(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def revise_text(text: str) -> str:
    replacements = {
        "donor/sample leave-one-out": "leave-one-sample-out / 逐样本剔除",
        "donor leave-one-out": "leave-one-sample-out",
        "donor-level supporting evidence": "sample-level supporting evidence",
        "donor-level": "sample-level",
        "donor/sample-level": "sample-level",
        "donor/sample": "sample-level",
        "donor/sample leave-one-out 稳健性验证": "leave-one-sample-out / 逐样本剔除稳健性分析",
        "donor 敏感": "sample 敏感",
        "donor_sensitive": "sample_sensitive",
    }
    out = text
    for old, new in replacements.items():
        out = out.replace(old, new)
    return out


def main() -> None:
    out_dir = ensure_functional_dir()
    print("修正内部稳健性验证术语：统一为 sample-level / leave-one-sample-out...")

    original_summary = ROBUSTNESS_DIR / "05_robustness_validation_master_summary_cn.txt"
    revised_summary = ROBUSTNESS_DIR / "05_robustness_validation_master_summary_cn_revised.txt"
    summary_text = safe_read(original_summary)
    preface = (
        "术语修正说明\n"
        "=" * 30
        + "\n"
        "当前 discovery 队列中 donor_id 仅 2 个 donor，而 sample_id 有 4 个 sample。\n"
        "上一轮稳健性验证实际采用 sample_id 作为剔除单位，因此应表述为 "
        "leave-one-sample-out / 逐样本剔除稳健性分析。\n"
        "后续写作与汇报中应统一使用 sample-level robustness，而不是 strict donor-level robustness。\n\n"
    )
    revised_summary.write_text(preface + revise_text(summary_text), encoding="utf-8")

    robustness_readme = ROOT / "scripts" / "robustness" / "README_robustness_validation.md"
    revised_readme = ROOT / "scripts" / "robustness" / "README_robustness_validation_revised.md"
    readme_text = safe_read(robustness_readme)
    if readme_text:
        revised_readme.write_text(preface + revise_text(readme_text), encoding="utf-8")

    note = (
        "内部稳健性验证命名修正说明\n"
        "=" * 40
        + "\n"
        "1. 当前实际完成的是 leave-one-sample-out / 逐样本剔除稳健性分析。\n"
        "2. 原因：discovery 队列 donor_id 仅 2 个 donor，但 sample_id 有 4 个 sample；"
        "为了避免严格 donor-level LOO 的样本量不足，本轮默认采用 sample_id 作为剔除单位。\n"
        "3. 后续写作中建议使用 sample-level robustness / leave-one-sample-out robustness。\n"
        "4. pseudobulk 部分也建议表述为 sample-level supporting evidence；如需 donor-level 表述，"
        "必须明确 donor_id 只有 2 个，统计解释应非常谨慎。\n"
        "5. 本脚本不覆盖原始 robustness 结果，仅生成修正版说明文件。\n"
    )
    (ROBUSTNESS_DIR / "00_naming_correction_note_cn.txt").write_text(note, encoding="utf-8")
    (out_dir / "00_naming_correction_note_cn.txt").write_text(note, encoding="utf-8")

    master = ROBUSTNESS_DIR / "02_leave_one_out_master_results.csv"
    if master.exists():
        df = pd.read_csv(master)
        rename_map = {
            "dropped_donor": "dropped_sample",
            "donor_sensitive_flag": "sample_sensitive_flag",
        }
        df = df.rename(columns=rename_map)
        df.to_csv(ROBUSTNESS_DIR / "02_leave_one_sample_out_master_results_revised.csv", index=False)
    tf_summary = ROBUSTNESS_DIR / "02_leave_one_out_tf_summary.csv"
    if tf_summary.exists():
        df = pd.read_csv(tf_summary)
        df = df.rename(columns={"donor_sensitive_flag": "sample_sensitive_flag"})
        df.to_csv(ROBUSTNESS_DIR / "02_leave_one_sample_out_tf_summary_revised.csv", index=False)

    print(f"已生成: {revised_summary}")
    print(f"已生成: {ROBUSTNESS_DIR / '00_naming_correction_note_cn.txt'}")
    print(f"已生成: {out_dir / '00_naming_correction_note_cn.txt'}")


if __name__ == "__main__":
    main()
