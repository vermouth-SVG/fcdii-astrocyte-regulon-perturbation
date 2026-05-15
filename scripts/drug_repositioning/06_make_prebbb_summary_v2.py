#!/usr/bin/env python3
from __future__ import annotations

from drugrep_common_v2 import make_prebbb_summary


def main() -> None:
    print("生成 v2 pre-BBB 主汇总...")
    make_prebbb_summary()
    print("完成 v2 pre-BBB 主汇总。")


if __name__ == "__main__":
    main()
