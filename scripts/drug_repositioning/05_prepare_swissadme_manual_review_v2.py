#!/usr/bin/env python3
from __future__ import annotations

from drugrep_common_v2 import prepare_swissadme_inputs


def main() -> None:
    print("生成 v2 SwissADME 手工复核输入。本步骤不执行 BBB 自动分层，不执行 docking。")
    df = prepare_swissadme_inputs()
    missing = int(df["canonical_smiles"].fillna("").astype(str).str.len().eq(0).sum())
    print(f"完成 SwissADME 输入准备: 候选 {df.shape[0]} 个，SMILES 缺失 {missing} 个。")


if __name__ == "__main__":
    main()
