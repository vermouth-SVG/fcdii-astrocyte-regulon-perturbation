#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from functional_common import (
    ENRICHR_LIBRARIES,
    FOCUS_TFS,
    GENESET_DIR,
    INPUT_H5AD,
    TF_TIER,
    download_enrichr_library,
    ensure_functional_dir,
    infer_themes_from_terms,
    load_discovery_h5ad_light,
    parse_gmt,
    run_ora,
    save_dotplot,
)


def load_detail(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"缺少 gene set detail: {path}。请先运行 02_build_tf_gene_sets.py")
    return json.loads(path.read_text(encoding="utf-8"))


def query_plan_for_tf(tf: str, detail: dict) -> dict[str, set[str]]:
    regulon = set(detail[tf]["regulon_targets"])
    intersection = set(detail[tf]["intersection_genes"])
    if tf in {"NFE2L2", "THRB"}:
        return {f"{tf}_regulon_targets": regulon, f"{tf}_intersection_directional": intersection}
    return {f"{tf}_regulon_targets": regulon}


def main() -> None:
    out_dir = ensure_functional_dir()
    print("执行 GO BP / KEGG 富集分析...")
    detail = load_detail(out_dir / "02_tf_gene_sets_detail.json")
    data = load_discovery_h5ad_light(INPUT_H5AD)
    universe = set(data["var_names"].astype(str))

    library_status = []
    libraries = {}
    for lib_name, source in ENRICHR_LIBRARIES.items():
        path = GENESET_DIR / f"{lib_name}.gmt"
        ok, message = download_enrichr_library(lib_name, path)
        library_status.append({"library": lib_name, "source": source, "path": str(path), "available": ok, "message": message})
        if ok:
            libraries[source] = parse_gmt(path)

    all_results = []
    top_by_tf = {}
    for tf in FOCUS_TFS:
        tf_results = []
        for query_label, genes in query_plan_for_tf(tf, detail).items():
            if len(genes) == 0:
                continue
            for source, gene_sets in libraries.items():
                enriched = run_ora(genes, universe, gene_sets, source=source, label=query_label)
                if not enriched.empty:
                    enriched.insert(0, "tf", tf)
                    enriched.insert(1, "candidate_tier", TF_TIER[tf])
                    enriched.insert(2, "gene_set_type", query_label.replace(tf + "_", ""))
                    all_results.append(enriched)
                    tf_results.append(enriched)
        tf_df = pd.concat(tf_results, ignore_index=True) if tf_results else pd.DataFrame()
        if not tf_df.empty:
            tf_df = tf_df.sort_values(["fdr_bh", "p_value", "overlap_count"], ascending=[True, True, False])
        top_by_tf[tf] = tf_df
        tf_df.head(30).to_csv(out_dir / f"03_{tf.lower()}_top_terms.csv", index=False)
        save_dotplot(tf_df, tf, out_dir / f"03_{tf}_GO_KEGG_dotplot.png")

    result = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    if not result.empty:
        go = result[result["source"] == "go"].copy()
        kegg = result[result["source"] == "kegg"].copy()
    else:
        go = pd.DataFrame()
        kegg = pd.DataFrame()
    go.to_csv(out_dir / "03_go_enrichment_all.csv", index=False)
    kegg.to_csv(out_dir / "03_kegg_enrichment_all.csv", index=False)

    lines = [
        "GO/KEGG 富集分析总结",
        "=" * 35,
        "数据库状态:",
    ]
    for item in library_status:
        lines.append(f"- {item['library']}: {'可用' if item['available'] else '不可用'}；{item['message']}；{item['path']}")
    lines.append("")
    if result.empty:
        lines.append("未获得可用富集结果。可能原因包括联网失败、基因集过小或 query genes 与背景基因交集不足。")
    else:
        for tf in FOCUS_TFS:
            tf_df = top_by_tf.get(tf, pd.DataFrame())
            lines.append(f"{tf}:")
            if tf_df.empty:
                lines.append("- 未获得显著或可排序富集结果。")
                continue
            terms = tf_df.head(12)["term"].astype(str).tolist()
            themes = infer_themes_from_terms(terms)
            nonzero = [f"{k}:{v}" for k, v in themes.items() if v > 0]
            lines.append("- Top terms: " + "; ".join(terms[:8]))
            lines.append("- 程序主题计数: " + ("; ".join(nonzero) if nonzero else "未形成明确主题"))
            if tf == "NFE2L2":
                lines.append("- 解释重点: 优先观察 oxidative stress / stress response / inflammatory-reactive / metabolic adaptation 相关主题。")
            elif tf == "THRB":
                lines.append("- 解释重点: 优先观察 thyroid hormone / neuronal-supportive / metabolic-homeostatic / synaptic support 相关主题。")
            elif tf == "BHLHE40":
                lines.append("- 解释重点: 作为第二梯队，观察 hypoxia/stress/reactive 或代谢调节相关主题。")
            else:
                lines.append("- 解释重点: SOX2 target 较少，富集结果只作轻量保留解释。")
            lines.append("")

    (out_dir / "03_functional_enrichment_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    pd.DataFrame(library_status).to_csv(out_dir / "03_gene_set_library_status.csv", index=False)
    print(f"完成 GO/KEGG 富集分析: {out_dir}")


if __name__ == "__main__":
    main()
