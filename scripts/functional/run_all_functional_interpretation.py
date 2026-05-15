#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
STEPS = [
    "00_revise_robustness_terms.py",
    "01_check_resources.py",
    "02_build_tf_gene_sets.py",
    "03_run_go_kegg.py",
    "04_program_convergence.py",
    "05_make_functional_master_summary.py",
]


def main() -> None:
    print("开始执行 GO/KEGG + 转录程序汇聚分析全流程...")
    for step in STEPS:
        print(f"\n>>> 运行 {step}")
        subprocess.check_call([sys.executable, str(SCRIPT_DIR / step)])
    print("\n功能解释全流程完成。")


if __name__ == "__main__":
    main()
