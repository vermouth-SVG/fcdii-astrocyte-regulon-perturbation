#!/usr/bin/env python3
from __future__ import annotations

from rescore_v21 import run_v21


def main() -> None:
    print("运行 drug repositioning v2.1 pre-BBB 收敛版。")
    print("本入口只重算现有 drug repurposing 表的 pre-BBB 分数/分层；不联网、不重跑上游、不做 docking、不做 BBB 自动搜索。")
    search, best_config, layered = run_v21()
    primary = layered[layered["lead_layer_v21"] == "primary_mechanism_direction_leads"]
    supportive = layered[layered["lead_layer_v21"] == "supportive_manual_review_leads"]
    print(f"完成 v2.1: 参数搜索 {search.shape[0]} 组。")
    print(f"最佳配置: {best_config.config_id}")
    print("primary leads:")
    print(primary[["compound", "associated_axis", "final_score_v21", "v21_penalty_reasons"]].to_string(index=False))
    print(f"supportive leads: {supportive.shape[0]} 个。")


if __name__ == "__main__":
    main()
