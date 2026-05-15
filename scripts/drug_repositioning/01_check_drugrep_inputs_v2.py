#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd

from drugrep_common_v2 import (
    BASELINE_DIR,
    DRUGREP_DIR,
    EXCLUDED_AXIS,
    FUNCTIONAL_CONVERGENCE,
    FUNCTIONAL_DETAIL,
    FUNCTIONAL_INTEGRATED,
    FUNCTIONAL_SUMMARY,
    PRIMARY_AXES,
    ROBUSTNESS_INTEGRATED,
    SUPPLEMENT_AXIS,
    build_query_gene_sets,
    ensure_dirs,
    gene_set_size_note,
    write_json,
    write_query_gene_files,
)


def main() -> None:
    ensure_dirs()
    print("检查 drug repositioning v2 输入资源...")
    queries, excluded = build_query_gene_sets()
    write_query_gene_files(queries)

    resource_rows = []
    for label, path in [
        ("functional summary", FUNCTIONAL_SUMMARY),
        ("functional integrated table", FUNCTIONAL_INTEGRATED),
        ("program convergence table", FUNCTIONAL_CONVERGENCE),
        ("gene set detail json", FUNCTIONAL_DETAIL),
        ("robustness integrated table", ROBUSTNESS_INTEGRATED),
    ]:
        resource_rows.append({"resource": label, "path": str(path), "exists": path.exists(), "size_bytes": path.stat().st_size if path.exists() else 0})
    pd.DataFrame(resource_rows).to_csv(DRUGREP_DIR / "01_resource_files_status.csv", index=False)
    pd.DataFrame(resource_rows).to_csv(BASELINE_DIR / "01_resource_files_status.csv", index=False)

    query_rows = []
    for name, payload in queries.items():
        suitable, note = gene_set_size_note(len(payload["genes"]))
        query_rows.append(
            {
                "query_gene_set_name": name,
                "axis": payload["axis"],
                "query_role": payload["query_role"],
                "gene_count": len(payload["genes"]),
                "source": payload["source"],
                "include_in_main": payload["include_in_main"],
                "suitable_for_drug_enrichment": suitable,
                "size_note": note,
            }
        )
    query_df = pd.DataFrame(query_rows)
    query_df.to_csv(DRUGREP_DIR / "01_query_gene_sets_summary.csv", index=False)
    query_df.to_csv(BASELINE_DIR / "01_query_gene_sets_summary.csv", index=False)
    write_json(DRUGREP_DIR / "01_query_gene_sets_detail.json", queries)

    lines = [
        "drug repositioning v2 输入资源检查",
        "=" * 45,
        f"输出目录: {DRUGREP_DIR}",
        "",
        "本轮只重做 drug repurposing，不重跑上游分析。",
        f"主轴: {', '.join(PRIMARY_AXES)}",
        f"补充轴: {SUPPLEMENT_AXIS}",
        f"不进入主查询: {EXCLUDED_AXIS}",
        "",
        "资源文件:",
        *[f"- {r['resource']}: {'存在' if r['exists'] else '缺失'}; {r['path']}" for r in resource_rows],
        "",
        "query gene sets:",
        query_df.to_string(index=False),
        "",
        "排除说明:",
        *[f"- {k}: {v}" for k, v in excluded.items()],
    ]
    text = "\n".join(lines) + "\n"
    (DRUGREP_DIR / "01_resource_check_cn.txt").write_text(text, encoding="utf-8")
    (BASELINE_DIR / "01_resource_check_cn.txt").write_text(text, encoding="utf-8")
    print(f"完成输入检查: {DRUGREP_DIR / '01_resource_check_cn.txt'}")


if __name__ == "__main__":
    main()
