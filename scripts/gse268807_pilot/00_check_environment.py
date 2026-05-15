# -*- coding: utf-8 -*-
"""
Phase 1：检查本地轻量 Python 环境。

- 【核心】metadata / 表格 / 基础作图：pandas、numpy、scipy、matplotlib —— 未通过则 exit 1。
- 【扩展】统计与 h5ad I/O：seaborn、statsmodels、anndata —— 缺失仅警告，不阻止通过。
- scanpy / leidenalg 仅用于后续降维、聚类、注释探索；当前 09→02→03 小样本链路不强依赖。

不检查 SCENIC+、pySCENIC、CellOracle 等重型依赖。

说明：在无法访问 PyPI 的网络环境下，仍应能完成 Phase 1–2 的元数据与制表；扩展包请稍后
`pip install -r 07_env/requirements.txt` 或 conda 安装。
"""

from __future__ import annotations

import sys
from typing import Iterable


def _python_ok() -> tuple[bool, str]:
    v = sys.version_info
    if v.major < 3 or (v.major == 3 and v.minor < 10):
        return False, f"当前 Python {v.major}.{v.minor}.{v.micro}，建议 3.10–3.11（64 位）。"
    if v.major == 3 and v.minor > 11:
        return True, f"当前 Python {v.major}.{v.minor}.{v.micro}（高于 3.11 时请自行确认 scanpy 轮子可用）。"
    return True, f"当前 Python {v.major}.{v.minor}.{v.micro}。"


def _try_import(module: str) -> tuple[bool, str | None]:
    try:
        __import__(module)
        return True, None
    except Exception as e:  # noqa: BLE001 — 环境检查需捕获所有导入错误
        return False, f"{type(e).__name__}: {e}"


def _report_group(title: str, modules: Iterable[str], *, optional: bool) -> bool:
    print(f"\n{title}")
    all_ok = True
    for name in modules:
        ok, err = _try_import(name)
        if ok:
            print(f"  【通过】{name}")
        else:
            all_ok = False
            tag = "【警告】" if optional else "【失败】"
            print(f"  {tag}{name}：{err}")
    return all_ok


def main() -> int:
    print("【环境检查】本地轻量依赖（不含 GRN 重型包）")

    py_ok, py_msg = _python_ok()
    print("\n【Python】", py_msg)
    if not py_ok:
        print("【汇总】请先升级 Python 后再安装依赖。")
        return 1

    ok_core = _report_group(
        "【核心】metadata / 表格 / 基础作图（未通过则整体失败）",
        ("pandas", "numpy", "scipy", "matplotlib"),
        optional=False,
    )
    ok_ext_plot = _report_group(
        "【扩展】统计作图（推荐）",
        ("seaborn",),
        optional=True,
    )
    ok_stats = _report_group(
        "【扩展】统计建模（推荐，Phase 6 用）",
        ("statsmodels",),
        optional=True,
    )
    ok_h5ad = _report_group(
        "【扩展】h5ad I/O（推荐，本地 09→02→03 需要）",
        ("anndata",),
        optional=True,
    )
    ok_scanpy = _report_group(
        "【可选】Scanpy / Leiden 聚类（后续注释、降维时再装）",
        ("scanpy", "leidenalg"),
        optional=True,
    )

    print("\n【汇总】")
    if not ok_core:
        print("  核心包缺失，请安装：pip install pandas numpy scipy matplotlib")
        return 1

    print("  【通过】核心环境满足 metadata 与表格、基础作图。")
    if not ok_h5ad:
        print("  【提示】本地 h5ad 流程至少需要：pip install anndata")
    if not ok_stats:
        print("  【提示】statsmodels 缺失时暂不能跑后续 GLM/混合模型统计。")
    if not ok_ext_plot:
        print("  【提示】seaborn 缺失时部分统计作图需改用 matplotlib 或补装。")
    if not ok_scanpy:
        print("  【提示】scanpy/leidenalg 缺失时可暂不做 Scanpy 降维聚类；需要时再用 pip 或 conda-forge 安装。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
