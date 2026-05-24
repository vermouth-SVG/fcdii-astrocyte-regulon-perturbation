#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
STEPS = [
    ["01_check_drugrep_inputs_v2.py"],
    ["02_run_dsigdb_enrichr_v2.py", "--run", "baseline"],
    ["03_aggregate_compounds_v2.py", "--run", "baseline"],
    ["04_refine_parameters_and_rerun_v2.py"],
    ["05_prepare_swissadme_manual_review_v2.py"],
    ["06_make_prebbb_summary_v2.py"],
]


def main() -> None:
    print("开始执行 drug repositioning v2：只重做药物部分，停在 BBB 手工复核前。")
    print("不会重跑 pySCENIC / CellOracle / perturbation / 外部验证 / 功能解释；不会执行 docking。")
    for step in STEPS:
        print("\n>>> 运行 " + " ".join(step))
        subprocess.check_call([sys.executable, str(SCRIPT_DIR / step[0]), *step[1:]])
    print("\ndrug repositioning v2 pre-BBB 流程完成。")


if __name__ == "__main__":
    main()
