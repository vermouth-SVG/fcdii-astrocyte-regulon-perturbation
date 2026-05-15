#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
STEPS = [
    "01_inspect_input_object.py",
    "02_leave_one_donor_out.py",
    "03_pseudobulk_validation.py",
    "04_shortlist_sensitivity.py",
    "05_make_master_summary.py",
]


def main() -> None:
    print("开始执行内部稳健性验证全流程...")
    for step in STEPS:
        path = SCRIPT_DIR / step
        print(f"\n>>> 运行 {step}")
        subprocess.check_call([sys.executable, str(path)])
    print("\n内部稳健性验证全流程完成。")


if __name__ == "__main__":
    main()
