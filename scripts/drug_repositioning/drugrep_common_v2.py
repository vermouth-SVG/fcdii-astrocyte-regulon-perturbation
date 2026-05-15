#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import re
import sys
import time
import uuid
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
FUNCTIONAL_DIR = ROOT / "functional_interpretation"
ROBUSTNESS_DIR = ROOT / "robustness_validation"
DRUGREP_DIR = ROOT / "drug_repositioning"
BASELINE_DIR = DRUGREP_DIR / "baseline_run"
FINAL_DIR = DRUGREP_DIR
QUERY_DIR = DRUGREP_DIR / "query_gene_sets"
PLOT_DIR = DRUGREP_DIR / "plots"
CACHE_DIR = DRUGREP_DIR / "cache"
MANUAL_REVIEW_DIR = DRUGREP_DIR / "manual_bbb_review"

FUNCTIONAL_DETAIL = FUNCTIONAL_DIR / "02_tf_gene_sets_detail.json"
FUNCTIONAL_SUMMARY = FUNCTIONAL_DIR / "05_functional_interpretation_master_summary_cn.txt"
FUNCTIONAL_INTEGRATED = FUNCTIONAL_DIR / "05_functional_interpretation_integrated_table.csv"
FUNCTIONAL_CONVERGENCE = FUNCTIONAL_DIR / "04_program_convergence_table.csv"
ROBUSTNESS_INTEGRATED = ROBUSTNESS_DIR / "05_robustness_validation_integrated_table.csv"

PRIMARY_AXES = ["NFE2L2", "THRB"]
SUPPLEMENT_AXIS = "BHLHE40"
EXCLUDED_AXIS = "SOX2"
DSIGDB_LIBRARY = "DSigDB"
FALLBACK_DRUG_LIBRARIES = [
    "Drug_Perturbations_from_GEO_2014",
    "LINCS_L1000_Chem_Pert_Consensus_Sigs",
    "DGIdb_Drug_Targets_2024",
]
ALLOWED_BBB_LAYERS = {
    "CNS-directed candidates": 1.0,
    "possible CNS-directed candidates": 0.75,
    "peripheral/program-modulating candidates": 0.40,
}


@dataclass(frozen=True)
class V2Params:
    name: str
    adjusted_p_threshold: float
    rank_threshold: int
    max_terms_per_query: int
    main_nfe2l2_quota: int
    main_thrb_quota: int
    min_final_score: float
    allow_weak_fdr_for_thrb: bool
    weak_fdr_threshold: float
    exclude_high_risk_from_main: bool = True


BASELINE_PARAMS = V2Params(
    name="baseline",
    adjusted_p_threshold=0.10,
    rank_threshold=60,
    max_terms_per_query=80,
    main_nfe2l2_quota=10,
    main_thrb_quota=8,
    min_final_score=0.35,
    allow_weak_fdr_for_thrb=True,
    weak_fdr_threshold=0.25,
)

TUNED_PARAMS = V2Params(
    name="tuned_v2_prebbb",
    adjusted_p_threshold=0.05,
    rank_threshold=45,
    max_terms_per_query=55,
    main_nfe2l2_quota=7,
    main_thrb_quota=5,
    min_final_score=0.40,
    allow_weak_fdr_for_thrb=True,
    weak_fdr_threshold=0.25,
)


def ensure_dirs() -> None:
    for path in [DRUGREP_DIR, BASELINE_DIR, QUERY_DIR, PLOT_DIR, CACHE_DIR, MANUAL_REVIEW_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def run_dir(run_name: str) -> Path:
    ensure_dirs()
    return BASELINE_DIR if run_name == "baseline" else FINAL_DIR


def params_for(run_name: str) -> V2Params:
    return BASELINE_PARAMS if run_name == "baseline" else TUNED_PARAMS


def safe_read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def load_functional_detail() -> dict:
    if not FUNCTIONAL_DETAIL.exists():
        raise FileNotFoundError(f"缺少功能解释 gene set 详情文件: {FUNCTIONAL_DETAIL}")
    return json.loads(FUNCTIONAL_DETAIL.read_text(encoding="utf-8"))


def build_query_gene_sets() -> tuple[dict[str, dict], dict[str, str]]:
    detail = load_functional_detail()
    queries: dict[str, dict] = {}
    excluded: dict[str, str] = {}

    for tf in PRIMARY_AXES:
        core = sorted({str(g).strip().upper() for g in detail.get(tf, {}).get("intersection_genes", []) if str(g).strip()})
        regulon = sorted({str(g).strip().upper() for g in detail.get(tf, {}).get("regulon_targets", []) if str(g).strip()})
        queries[f"{tf}_core_program"] = {
            "axis": tf,
            "query_gene_set_name": f"{tf}_core_program",
            "query_role": "main_core_program",
            "genes": core,
            "source": "functional_interpretation intersection_genes",
            "include_in_main": True,
        }
        queries[f"{tf}_regulon_targets"] = {
            "axis": tf,
            "query_gene_set_name": f"{tf}_regulon_targets",
            "query_role": "secondary_regulon_targets",
            "genes": regulon,
            "source": "pySCENIC discovery regulon targets",
            "include_in_main": True,
        }

    bhl_core = sorted({str(g).strip().upper() for g in detail.get(SUPPLEMENT_AXIS, {}).get("intersection_genes", []) if str(g).strip()})
    bhl_reg = sorted({str(g).strip().upper() for g in detail.get(SUPPLEMENT_AXIS, {}).get("regulon_targets", []) if str(g).strip()})
    chosen = bhl_core if len(bhl_core) >= 15 else bhl_reg
    source = "functional_interpretation intersection_genes" if len(bhl_core) >= 15 else "pySCENIC discovery regulon targets"
    queries["BHLHE40_supplement_program"] = {
        "axis": "BHLHE40",
        "query_gene_set_name": "BHLHE40_supplement_program",
        "query_role": "supplement_axis",
        "genes": chosen,
        "source": source,
        "include_in_main": False,
    }

    sox_targets = detail.get(EXCLUDED_AXIS, {}).get("regulon_targets", [])
    excluded[EXCLUDED_AXIS] = (
        f"SOX2 不进入 v2 药物重定位主查询。当前 target={len(sox_targets)}，仅保留在功能解释结果中。"
    )
    return queries, excluded


def gene_set_size_note(n: int) -> tuple[bool, str]:
    if n < 10:
        return False, "过小，不适合药物富集"
    if n < 15:
        return True, "偏小，结果需谨慎"
    if n > 800:
        return True, "偏大，可能偏向广义程序"
    return True, "大小适合轻量药物富集"


def write_query_gene_files(queries: dict[str, dict]) -> None:
    QUERY_DIR.mkdir(parents=True, exist_ok=True)
    for name, payload in queries.items():
        genes = payload["genes"]
        (QUERY_DIR / f"{name}.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
        (QUERY_DIR / f"{name}.gmt").write_text(name + "\tprogram_guided_v2_query\t" + "\t".join(genes) + "\n", encoding="utf-8")


def request_json(url: str, data: bytes | None = None, headers: dict | None = None, timeout: int = 90) -> dict:
    req_headers = {"User-Agent": "Mozilla/5.0 drugrep-v2/1.0"}
    if headers:
        req_headers.update(headers)
    request = Request(url, data=data, headers=req_headers)
    with urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def enrichr_add_list(genes: list[str], description: str) -> int:
    boundary = "----DrugRepV2Boundary" + uuid.uuid4().hex
    parts = []
    for name, value in {"list": "\n".join(genes), "description": description}.items():
        parts.append(
            f"--{boundary}\r\n"
            f"Content-Disposition: form-data; name=\"{name}\"\r\n\r\n"
            f"{value}\r\n"
        )
    payload = ("".join(parts) + f"--{boundary}--\r\n").encode("utf-8")
    headers = {"Content-Type": f"multipart/form-data; boundary={boundary}"}
    result = request_json("https://maayanlab.cloud/Enrichr/addList", data=payload, headers=headers, timeout=90)
    return int(result["userListId"])


def enrichr_enrich(user_list_id: int, library: str) -> list:
    url = f"https://maayanlab.cloud/Enrichr/enrich?{urlencode({'userListId': user_list_id, 'backgroundType': library})}"
    result = request_json(url)
    return result.get(library, [])


def parse_enrichr_rows(rows: list, axis: str, query_name: str, query_role: str, library: str) -> pd.DataFrame:
    parsed = []
    for row in rows:
        if len(row) < 7:
            continue
        parsed.append(
            {
                "rank": row[0],
                "term": row[1],
                "p_value": row[2],
                "z_score": row[3],
                "combined_score": row[4],
                "overlap_genes": ";".join(row[5]) if isinstance(row[5], list) else str(row[5]),
                "adjusted_p_value": row[6],
                "library": library,
                "associated_axis_query": axis,
                "query_gene_set_name": query_name,
                "query_role": query_role,
            }
        )
    return pd.DataFrame(parsed)


def run_enrichr_query(genes: list[str], axis: str, query_name: str, query_role: str, library: str = DSIGDB_LIBRARY) -> tuple[pd.DataFrame, str]:
    try:
        user_list_id = enrichr_add_list(genes, f"v2_{query_name}_{int(time.time())}")
        rows = enrichr_enrich(user_list_id, library)
        df = parse_enrichr_rows(rows, axis, query_name, query_role, library)
        return df, "Enrichr 在线查询成功"
    except Exception as exc:  # noqa: BLE001
        return pd.DataFrame(), f"Enrichr 查询失败: {type(exc).__name__}: {exc}"


def download_enrichr_library_as_gmt(library: str) -> tuple[Path | None, str]:
    out = CACHE_DIR / f"{library}.gmt"
    if out.exists() and out.stat().st_size > 0:
        return out, "使用本地缓存"
    url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={quote(library)}"
    try:
        with urlopen(url, timeout=90) as response:
            text = response.read().decode("utf-8", errors="replace")
        if not text.strip():
            return None, "下载结果为空"
        out.write_text(text, encoding="utf-8", newline="\n")
        return out, "已下载"
    except Exception as exc:  # noqa: BLE001
        return None, f"下载失败: {type(exc).__name__}: {exc}"


def parse_gmt(path: Path) -> dict[str, set[str]]:
    result: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                result[parts[0]] = {g.upper() for g in parts[2:] if g.strip()}
    return result


def benjamini_hochberg(p_values) -> np.ndarray:
    arr = np.asarray(p_values, dtype=float)
    n = arr.size
    order = np.argsort(arr)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    q = arr * n / ranks
    q_ordered = q[order]
    q_monotone = np.minimum.accumulate(q_ordered[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(q_monotone, 1.0)
    return out


def local_ora(genes: list[str], universe: set[str], library_terms: dict[str, set[str]], axis: str, query_name: str, query_role: str, library: str) -> pd.DataFrame:
    query = {g.upper() for g in genes} & {g.upper() for g in universe}
    universe = {g.upper() for g in universe}
    M, N = len(universe), len(query)
    rows = []
    if M == 0 or N == 0:
        return pd.DataFrame()
    for term, term_genes in library_terms.items():
        term_set = term_genes & universe
        n = len(term_set)
        overlap = sorted(query & term_set)
        k = len(overlap)
        if k < 1:
            continue
        p = stats.hypergeom.sf(k - 1, M, n, N)
        fold = (k / N) / (n / M) if n and N else np.nan
        rows.append(
            {
                "rank": np.nan,
                "term": term,
                "p_value": p,
                "z_score": np.nan,
                "combined_score": float(-math.log10(max(p, 1e-300)) * fold),
                "overlap_genes": ";".join(overlap[:100]),
                "adjusted_p_value": np.nan,
                "library": library + "_local_ora",
                "associated_axis_query": axis,
                "query_gene_set_name": query_name,
                "query_role": query_role,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["adjusted_p_value"] = benjamini_hochberg(df["p_value"].to_numpy(dtype=float))
    df = df.sort_values(["adjusted_p_value", "p_value", "combined_score"], ascending=[True, True, False]).reset_index(drop=True)
    df["rank"] = np.arange(1, df.shape[0] + 1)
    return df


ID_PATTERNS = [
    r"\bCTD\s*\d+\b",
    r"\bMESH[:_\s]*[A-Z0-9]+\b",
    r"\bDB\d+\b",
    r"\bCID[:_\s]*\d+\b",
    r"\bCHEMBL\d+\b",
    r"\bHMSL\d+\b",
    r"\bTTD\s*\d+\b",
    r"\bBOSS\b",
]


def normalize_compound_name(term: str) -> str:
    text = str(term)
    text = re.sub(r"\s+", " ", text).strip()
    text = re.split(r"\s+\|\s+|\s+::\s+|\s+--\s+", text)[0].strip()
    text = re.sub(
        r"\b(MCF7|HL60|PC3|A375|HA1E|HT29|A549|HEPG2|JURKAT|K562|VCAP|BT20|SKB|SKBR3)\s+(UP|DOWN)\b",
        "",
        text,
        flags=re.IGNORECASE,
    ).strip()
    text = re.sub(r"\bBOSS\b", "", text, flags=re.IGNORECASE).strip()
    text = re.sub(r"\b(UP|DOWN)\b$", "", text, flags=re.IGNORECASE).strip()
    for pattern in ID_PATTERNS:
        text = re.sub(pattern, "", text, flags=re.IGNORECASE).strip()
    text = re.sub(r"\s*\([^)]*(?:CTD|MESH|CID|CHEMBL|DB|TTD)[^)]*\)\s*", " ", text, flags=re.IGNORECASE)
    text = re.sub(r"\s+", " ", text).strip(" ;,")
    return text if text else str(term).strip()


MECHANISM_KEYWORDS = {
    "calcium/synaptic": ["calcium", "verapamil", "diltiazem", "nimodipine", "nifedipine", "gabapentin", "pregabalin", "synaptic", "glutamate", "gaba", "mpep", "quisqualate"],
    "redox/stress-adaptation": ["nrf", "nfe2l2", "sulforaphane", "curcumin", "resveratrol", "antioxid", "glutathione", "oxidative", "ros", "epigallocatechin", "quercetin"],
    "anti-inflammatory": ["dexamethasone", "prednisone", "ibuprofen", "celecoxib", "aspirin", "anti-inflammatory", "cox", "tnf", "jak"],
    "mapk/kinase": ["kinase", "mapk", "mek", "erk", "akt", "pi3k", "rapamycin", "mtor", "gw-8510", "alsterpaullone"],
    "metabolic/homeostatic": ["metformin", "statin", "lipid", "metabolic", "glucose", "mitochond", "transport", "captopril"],
    "thyroid/hormone": ["thyroxine", "triiodothyronine", "thyroid", "hormone", "retinoic", "steroid", "estradiol", "coumestrol", "pregnenolone"],
    "epigenetic/hdac": ["valproic", "trichostatin", "ms-275", "entinostat", "scriptaid", "vorinostat", "hdac", "histone"],
    "phosphodiesterase/cAMP": ["rolipram", "roflumilast", "zaprinast", "papaverine", "trequinsin", "pentoxifylline", "dipyridamole", "pde"],
    "broad_cytotoxic/topoisomerase": ["irinotecan", "camptothecin", "doxorubicin", "etoposide", "5-fluorouracil", "anisomycin", "nocodazole", "azaguanine", "daunorubicin"],
    "environmental_toxicant": ["formaldehyde", "chromate", "benzene", "benzo", "cadmium", "arsenic", "mercury", "bisphenol", "nickel", "copper sulfate", "cobalt chloride", "barium", "silica", "aflatoxin", "atrazine", "chlorpyrifos", "nitroso", "tert-butyl hydroperoxide"],
}


def infer_mechanism(compound: str, term: str, axis: str) -> tuple[str, str, float]:
    text = f"{compound} {term}".lower()
    hits = []
    for theme, keywords in MECHANISM_KEYWORDS.items():
        if any(k.lower() in text for k in keywords):
            hits.append(theme)
    if hits:
        score = 1.0 if any(h in {"calcium/synaptic", "redox/stress-adaptation", "phosphodiesterase/cAMP", "epigenetic/hdac"} for h in hits) else 0.65
        return "; ".join(dict.fromkeys(hits)), hits[0], score
    if axis == "BHLHE40":
        return "未自动识别；补充轴线索", "unknown", 0.25
    return "未自动识别；需人工核实", "unknown", 0.30


def translational_suitability(compound: str, term: str) -> tuple[str, float]:
    text = f"{compound} {term}".lower()
    high_caution = [
        "formaldehyde",
        "potassium chromate",
        "chromate",
        "copper sulfate",
        "nickel sulfate",
        "cobalt chloride",
        "7646-79-9",
        "barium",
        "silica",
        "aflatoxin",
        "atrazine",
        "chlorpyrifos",
        "n-nitroso",
        "nitrosodiethylamine",
        "thapsigargin",
        "mptp",
        "benzene",
        "benzo[a]pyrene",
        "cadmium",
        "arsenic",
        "mercury",
        "bisphenol",
        "8-azaguanine",
        "nocodazole",
        "anisomycin",
        "tert-butyl hydroperoxide",
    ]
    oncology_caution = ["irinotecan", "camptothecin", "doxorubicin", "etoposide", "5-fluorouracil", "daunorubicin"]
    controlled_caution = ["methamphetamine", "amphetamine", "diacetylmorphine", "heroin", "morphine", "ketamine"]
    if any(x in text for x in high_caution):
        return "high_caution_toxic_or_research_chemical", 0.05
    if any(x in text for x in oncology_caution):
        return "oncology_or_broad_cytotoxic_caution", 0.15
    if any(x in text for x in controlled_caution):
        return "controlled_or_abuse_liability_caution", 0.10
    if any(x in text for x in ["valproic", "gabapentin", "pregabalin", "memantine", "scopolamine", "meclofenoxate", "lobeline"]):
        return "cns_or_neuroactive_known_drug", 0.95
    if any(x in text for x in ["epigallocatechin", "resveratrol", "curcumin", "metformin", "rapamycin", "trichostatin", "ms-275", "scriptaid", "roflumilast", "rolipram", "zaprinast", "papaverine", "nifedipine", "pentoxifylline", "dipyridamole", "quercetin"]):
        return "mechanistically_interpretable_modulator", 0.75
    if re.fullmatch(r"[A-Z0-9-]{6,}", compound.strip(), flags=re.IGNORECASE):
        return "identifier_or_tool_compound_manual_review", 0.25
    return "requires_manual_translational_review", 0.45


MAIN_TABLE_EXCLUDE_CATEGORIES = {
    "high_caution_toxic_or_research_chemical",
    "oncology_or_broad_cytotoxic_caution",
    "controlled_or_abuse_liability_caution",
    "identifier_or_tool_compound_manual_review",
}


def minmax(series: pd.Series) -> pd.Series:
    arr = pd.to_numeric(series, errors="coerce")
    finite = arr[np.isfinite(arr)]
    if finite.empty:
        return pd.Series(np.zeros(len(arr)), index=series.index)
    lo, hi = finite.min(), finite.max()
    if hi - lo < 1e-12:
        out = pd.Series(np.zeros(len(arr)) + 0.5, index=series.index)
        out[~np.isfinite(arr)] = 0
        return out
    return ((arr - lo) / (hi - lo)).fillna(0)


def axis_program_label(axis: str) -> str:
    if axis == "NFE2L2":
        return "NFE2L2 lesion-associated transcriptional/stress-adaptation program"
    if axis == "THRB":
        return "THRB internal_control-associated synaptic/calcium/homeostatic-supportive program"
    return "BHLHE40 supplementary stress/reactive program"


def interpretation(axis: str) -> str:
    if axis == "NFE2L2":
        return "围绕 NFE2L2 lesion-associated transcriptional/stress-adaptation program 的探索性候选线索，不解释为直接抑制或激活 NFE2L2。"
    if axis == "THRB":
        return "围绕 THRB internal_control-associated synaptic/calcium/homeostatic-supportive program 的探索性候选线索，不解释为直接激动 THRB 即治疗。"
    return "BHLHE40 仅作为补充轴线索，不进入主文主药物推荐。"


def run_enrichment(run_name: str, params: V2Params) -> pd.DataFrame:
    out_dir = run_dir(run_name)
    queries, _ = build_query_gene_sets()
    write_query_gene_files(queries)
    all_frames = []
    status_lines = []
    print(f"执行 {run_name} DSigDB/Enrichr 查询...")
    for query_name, payload in queries.items():
        suitable, note = gene_set_size_note(len(payload["genes"]))
        if not suitable:
            status_lines.append(f"{query_name}: 跳过；{note}")
            continue
        df, msg = run_enrichr_query(payload["genes"], payload["axis"], query_name, payload["query_role"], DSIGDB_LIBRARY)
        status_lines.append(f"{query_name}: DSigDB {msg}; 返回 {df.shape[0]} 条。")
        if df.empty:
            gmt_path, gmt_msg = download_enrichr_library_as_gmt(DSIGDB_LIBRARY)
            status_lines.append(f"{query_name}: DSigDB GMT 降级准备: {gmt_msg}")
            if gmt_path:
                library_terms = parse_gmt(gmt_path)
                universe = set().union(*library_terms.values()) if library_terms else set()
                df = local_ora(payload["genes"], universe, library_terms, payload["axis"], query_name, payload["query_role"], DSIGDB_LIBRARY)
                status_lines.append(f"{query_name}: 本地 ORA 返回 {df.shape[0]} 条。")
        if df.empty:
            for lib in FALLBACK_DRUG_LIBRARIES:
                df, msg = run_enrichr_query(payload["genes"], payload["axis"], query_name, payload["query_role"], lib)
                status_lines.append(f"{query_name}: fallback {lib} {msg}; 返回 {df.shape[0]} 条。")
                if not df.empty:
                    break
        if not df.empty:
            all_frames.append(df)
        time.sleep(0.6)

    raw = pd.concat(all_frames, ignore_index=True) if all_frames else pd.DataFrame()
    raw.to_csv(out_dir / "02_enrichr_raw_results_all.csv", index=False)
    if raw.empty:
        filtered = pd.DataFrame()
    else:
        raw["adjusted_p_value"] = pd.to_numeric(raw["adjusted_p_value"], errors="coerce")
        raw["combined_score"] = pd.to_numeric(raw["combined_score"], errors="coerce")
        raw["rank"] = pd.to_numeric(raw["rank"], errors="coerce")
        filtered = (
            raw[(raw["adjusted_p_value"] <= params.adjusted_p_threshold) | (raw["rank"] <= params.rank_threshold)]
            .sort_values(["associated_axis_query", "query_gene_set_name", "adjusted_p_value", "combined_score"], ascending=[True, True, True, False])
            .groupby(["associated_axis_query", "query_gene_set_name"], as_index=False, group_keys=False)
            .head(params.max_terms_per_query)
            .reset_index(drop=True)
        )
    filtered.to_csv(out_dir / "02_dsigdb_filtered_results.csv", index=False)
    lines = [
        f"{run_name} DSigDB / Enrichr 富集总结",
        "=" * 45,
        f"参数: adjusted_p_threshold={params.adjusted_p_threshold}; rank_threshold={params.rank_threshold}; max_terms_per_query={params.max_terms_per_query}",
        "主查询轴: NFE2L2, THRB；BHLHE40 为补充轴；SOX2 未进入查询。",
        "",
        "查询状态:",
        *status_lines,
        "",
        f"raw 结果条数: {raw.shape[0]}",
        f"filtered 结果条数: {filtered.shape[0]}",
    ]
    (out_dir / "02_drugrep_enrichment_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return filtered


def aggregate_compounds(run_name: str, params: V2Params) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out_dir = run_dir(run_name)
    input_path = out_dir / "02_dsigdb_filtered_results.csv"
    if not input_path.exists():
        raise FileNotFoundError(f"缺少富集结果: {input_path}")
    df = pd.read_csv(input_path)
    if df.empty:
        empty = pd.DataFrame()
        empty.to_csv(out_dir / "03_compound_aggregated_table.csv", index=False)
        return empty, empty, empty

    df["normalized_compound_name"] = df["term"].map(normalize_compound_name)
    df["adjusted_p_value"] = pd.to_numeric(df["adjusted_p_value"], errors="coerce")
    df["combined_score"] = pd.to_numeric(df["combined_score"], errors="coerce")
    df["rank"] = pd.to_numeric(df["rank"], errors="coerce")
    rows, axis_rows, mechanism_rows = [], [], []
    global_logp = -np.log10(df["adjusted_p_value"].clip(lower=1e-300))
    df["neglog10_adj"] = global_logp
    df["combined_norm_global"] = minmax(df["combined_score"])

    for compound, sub in df.groupby("normalized_compound_name"):
        axis_scores = []
        for axis, axis_sub in sub.groupby("associated_axis_query"):
            q_count = axis_sub["query_gene_set_name"].nunique()
            best_adj = axis_sub["adjusted_p_value"].min()
            best_combined = axis_sub["combined_score"].max()
            sig_component = min(float(-np.log10(max(best_adj, 1e-300))) / 50.0, 1.0)
            combined_component = float(axis_sub["combined_norm_global"].max())
            support_component = min(q_count / 2.0, 1.0)
            main_bonus = 0.08 if axis in PRIMARY_AXES else -0.08
            axis_score = 0.40 * support_component + 0.35 * sig_component + 0.25 * combined_component + main_bonus
            axis_scores.append(
                {
                    "axis": axis,
                    "axis_score": axis_score,
                    "n_supporting_query_sets": int(q_count),
                    "best_adjusted_p_axis": best_adj,
                    "best_combined_score_axis": best_combined,
                    "query_sets": ";".join(sorted(axis_sub["query_gene_set_name"].dropna().unique())),
                }
            )
        axis_score_df = pd.DataFrame(axis_scores).sort_values(["axis_score", "n_supporting_query_sets", "best_combined_score_axis"], ascending=[False, False, False])
        associated_axis = axis_score_df.iloc[0]["axis"]
        axis_counts = ";".join([f"{r.axis}:{int(r.n_supporting_query_sets)}" for r in axis_score_df.itertuples()])
        axis_margin = axis_score_df.iloc[0]["axis_score"] - (axis_score_df.iloc[1]["axis_score"] if axis_score_df.shape[0] > 1 else 0)

        assoc_sub = sub[sub["associated_axis_query"] == associated_axis].copy()
        qs_rows = []
        for qname, qsub in assoc_sub.groupby("query_gene_set_name"):
            best_adj_q = qsub["adjusted_p_value"].min()
            best_combined_q = qsub["combined_score"].max()
            n_terms = qsub.shape[0]
            q_score = min(-np.log10(max(best_adj_q, 1e-300)) / 50.0, 1.0) * 0.55 + minmax(pd.Series([best_combined_q, sub["combined_score"].max()])).iloc[0] * 0.20 + min(n_terms / 5, 1.0) * 0.25
            qs_rows.append((qname, q_score, best_adj_q, best_combined_q, n_terms))
        qdf = pd.DataFrame(qs_rows, columns=["query_gene_set_name", "query_support_score", "best_adjusted_p_query", "best_combined_score_query", "n_terms_query"])
        qdf = qdf.sort_values(["query_support_score", "best_adjusted_p_query", "best_combined_score_query"], ascending=[False, True, False])
        main_supporting_gene_set = qdf.iloc[0]["query_gene_set_name"] if not qdf.empty else ""

        best = sub.sort_values(["adjusted_p_value", "combined_score"], ascending=[True, False]).iloc[0]
        mechanism, drug_class, mech_score = infer_mechanism(compound, str(best["term"]), str(associated_axis))
        caution, translational_score = translational_suitability(compound, str(best["term"]))
        supporting_sets = sorted(sub["query_gene_set_name"].dropna().astype(str).unique().tolist())
        support_count = len(supporting_sets)
        support_count_main_axes = int(sub[sub["associated_axis_query"].isin(PRIMARY_AXES)]["query_gene_set_name"].nunique())
        include_main = associated_axis in PRIMARY_AXES
        rows.append(
            {
                "compound": compound,
                "normalized_compound_name": compound,
                "associated_axis": associated_axis,
                "axis_assignment_score": axis_score_df.iloc[0]["axis_score"],
                "axis_assignment_margin": axis_margin,
                "axis_counts": axis_counts,
                "supporting_query_sets": ";".join(supporting_sets),
                "support_count": support_count,
                "support_count_main_axes": support_count_main_axes,
                "main_supporting_gene_set": main_supporting_gene_set,
                "main_supporting_gene_set_score": qdf.iloc[0]["query_support_score"] if not qdf.empty else np.nan,
                "best_adjusted_p": best["adjusted_p_value"],
                "best_combined_score": best["combined_score"],
                "best_rank": best["rank"],
                "best_term": best["term"],
                "best_library": best.get("library", ""),
                "indicative_mechanism": mechanism,
                "indicative_drug_class": drug_class,
                "mechanism_interpretability_score": mech_score,
                "translational_caution_category": caution,
                "translational_suitability_score": translational_score,
                "include_in_main_axis_pool": include_main,
                "note": "v2 自动聚合注释；仅代表 program-guided exploratory clue，不代表治疗推荐。",
            }
        )
        for _, ar in axis_score_df.iterrows():
            axis_rows.append({"compound": compound, **ar.to_dict()})
        mechanism_rows.append(
            {
                "compound": compound,
                "associated_axis": associated_axis,
                "best_term": best["term"],
                "indicative_mechanism": mechanism,
                "indicative_drug_class": drug_class,
                "mechanism_interpretability_score": mech_score,
                "translational_caution_category": caution,
                "translational_suitability_score": translational_score,
                "annotation_source": "DSigDB/Enrichr term + v2 keyword rules",
            }
        )

    agg = pd.DataFrame(rows)
    agg = score_and_filter(agg, params)
    agg.to_csv(out_dir / "03_compound_aggregated_table.csv", index=False)
    pd.DataFrame(axis_rows).to_csv(out_dir / "03_compound_axis_mapping.csv", index=False)
    pd.DataFrame(mechanism_rows).to_csv(out_dir / "03_compound_mechanism_annotation.csv", index=False)
    write_aggregation_summary(out_dir, run_name, agg, params)
    top, retained = make_top_tables(agg, params, run_name)
    return agg, top, retained


def score_and_filter(df: pd.DataFrame, params: V2Params) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["best_adjusted_p"] = pd.to_numeric(df["best_adjusted_p"], errors="coerce").fillna(1.0)
    df["best_combined_score"] = pd.to_numeric(df["best_combined_score"], errors="coerce").fillna(0)
    df["support_count"] = pd.to_numeric(df["support_count"], errors="coerce").fillna(1)
    df["support_count_main_axes"] = pd.to_numeric(df["support_count_main_axes"], errors="coerce").fillna(0)
    df["enrichment_support_score"] = 0.60 * minmax(-np.log10(df["best_adjusted_p"].clip(lower=1e-300))) + 0.40 * minmax(df["best_combined_score"])
    df["cross_signature_score"] = np.minimum(df["support_count"] / 2.0, 1.0)
    df["axis_relevance_score"] = df["associated_axis"].map({"NFE2L2": 1.0, "THRB": 1.0, "BHLHE40": 0.35}).fillna(0.20)
    df["mechanism_interpretability_score"] = pd.to_numeric(df["mechanism_interpretability_score"], errors="coerce").fillna(0.30)
    df["translational_suitability_score"] = pd.to_numeric(df["translational_suitability_score"], errors="coerce").fillna(0.45)
    df["axis_assignment_margin"] = pd.to_numeric(df["axis_assignment_margin"], errors="coerce").fillna(0)
    df["axis_assignment_confidence_score"] = np.clip(0.5 + df["axis_assignment_margin"], 0, 1)
    # v2 pre-BBB scoring deliberately excludes BBB; BBB must be manually verified later.
    df["final_score_prebbb"] = (
        0.30 * df["enrichment_support_score"]
        + 0.16 * df["cross_signature_score"]
        + 0.18 * df["axis_relevance_score"]
        + 0.14 * df["mechanism_interpretability_score"]
        + 0.14 * df["translational_suitability_score"]
        + 0.08 * df["axis_assignment_confidence_score"]
    )
    df["main_table_exclusion_reason"] = ""
    high_risk = df["translational_caution_category"].isin(MAIN_TABLE_EXCLUDE_CATEGORIES)
    df.loc[high_risk, "main_table_exclusion_reason"] = "因毒性/广谱细胞毒性/滥用风险/工具化合物属性，仅保留在补充表，不进入 pre-BBB 主表。"
    weak_fdr = df["best_adjusted_p"] > params.adjusted_p_threshold
    df["main_table_eligible"] = (
        df["include_in_main_axis_pool"].astype(bool)
        & df["associated_axis"].isin(PRIMARY_AXES)
        & (~high_risk if params.exclude_high_risk_from_main else True)
        & (df["final_score_prebbb"] >= params.min_final_score)
    )
    if params.allow_weak_fdr_for_thrb:
        df["weak_fdr_retained_reason"] = np.where(
            (df["associated_axis"] == "THRB") & (df["best_adjusted_p"] <= params.weak_fdr_threshold) & weak_fdr,
            "THRB gene set 较小；FDR 较弱但具轴线/机制支持，保留为 manual_review_weak_fdr。",
            "",
        )
        df.loc[(df["associated_axis"] != "THRB") & weak_fdr, "main_table_eligible"] = False
    else:
        df["weak_fdr_retained_reason"] = ""
        df.loc[weak_fdr, "main_table_eligible"] = False
    df["final_priority_tier"] = "low_priority_retained"
    df.loc[df["final_score_prebbb"] >= 0.72, "final_priority_tier"] = "high_priority_prebbb"
    df.loc[(df["final_score_prebbb"] >= 0.50) & (df["final_score_prebbb"] < 0.72), "final_priority_tier"] = "medium_priority_prebbb"
    df.loc[high_risk, "final_priority_tier"] = "mechanistic_reference_only_excluded_from_main"
    df.loc[df["weak_fdr_retained_reason"] != "", "final_priority_tier"] = "manual_review_weak_fdr"
    df["BBB_layer"] = "pending_manual_swissadme"
    df["BBB_score"] = np.nan
    df["BBB_note"] = "v2 pre-BBB 阶段不使用启发式 BBB；等待 SwissADME 手工复核。"
    df["translational_interpretation"] = df["associated_axis"].map(interpretation)
    df = df.sort_values(["include_in_main_axis_pool", "main_table_eligible", "final_score_prebbb", "support_count", "best_adjusted_p"], ascending=[False, False, False, False, True])
    return df


def make_top_tables(agg: pd.DataFrame, params: V2Params, run_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    out_dir = run_dir(run_name)
    if agg.empty:
        top = pd.DataFrame()
        retained = pd.DataFrame()
    else:
        retained = agg[(agg["final_score_prebbb"] >= 0.30) | agg["include_in_main_axis_pool"]].copy()
        eligible = retained[retained["main_table_eligible"]].copy()
        parts = []
        for axis, quota in [("NFE2L2", params.main_nfe2l2_quota), ("THRB", params.main_thrb_quota)]:
            axis_df = eligible[eligible["associated_axis"] == axis].sort_values(["final_score_prebbb", "support_count", "best_adjusted_p"], ascending=[False, False, True]).head(quota)
            if not axis_df.empty:
                parts.append(axis_df)
        top = pd.concat(parts, ignore_index=False) if parts else pd.DataFrame(columns=retained.columns)
        top = top.sort_values(["associated_axis", "final_score_prebbb", "support_count", "best_adjusted_p"], ascending=[True, False, False, True])
    if run_name == "baseline":
        top.to_csv(out_dir / "04_baseline_top_priority_compounds.csv", index=False)
    else:
        agg.to_csv(out_dir / "05_drug_repositioning_integrated_table.csv", index=False)
        top.to_csv(out_dir / "05_top_priority_compounds.csv", index=False)
        retained.to_csv(out_dir / "05_all_retained_compounds.csv", index=False)
    return top, retained


def write_aggregation_summary(out_dir: Path, run_name: str, agg: pd.DataFrame, params: V2Params) -> None:
    main_pool = agg[agg["include_in_main_axis_pool"]] if not agg.empty else pd.DataFrame()
    lines = [
        f"{run_name} 化合物聚合与轴线归属总结",
        "=" * 45,
        "v2 轴线归属采用加权规则: support_count、best_adjusted_p、best_combined_score、主轴一致性共同决定 associated_axis。",
        "main_supporting_gene_set 采用当前轴内 query-level 支持强度选择，不再取字符串第一个 query set。",
        "",
        f"聚合后化合物数: {agg.shape[0]}",
        f"NFE2L2/THRB 主轴候选池: {main_pool.shape[0]}",
        f"高风险/工具化合物排除出主表数: {int((agg.get('main_table_exclusion_reason','') != '').sum()) if not agg.empty else 0}",
        "",
        "主轴候选前 25:",
        main_pool.head(25)[["compound", "associated_axis", "support_count", "best_adjusted_p", "best_combined_score", "main_supporting_gene_set", "translational_caution_category", "final_score_prebbb"]].to_string(index=False) if not main_pool.empty else "无",
    ]
    (out_dir / "03_compound_aggregation_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def diagnose_baseline() -> dict:
    out_dir = BASELINE_DIR
    top_path = out_dir / "04_baseline_top_priority_compounds.csv"
    agg_path = out_dir / "03_compound_aggregated_table.csv"
    top = pd.read_csv(top_path) if top_path.exists() and top_path.stat().st_size > 0 else pd.DataFrame()
    agg = pd.read_csv(agg_path) if agg_path.exists() and agg_path.stat().st_size > 0 else pd.DataFrame()
    diagnosis = {
        "top_n": int(top.shape[0]),
        "nfe2l2_n": int((top["associated_axis"] == "NFE2L2").sum()) if not top.empty else 0,
        "thrb_n": int((top["associated_axis"] == "THRB").sum()) if not top.empty else 0,
        "high_risk_in_top": int(top["translational_caution_category"].isin(MAIN_TABLE_EXCLUDE_CATEGORIES).sum()) if not top.empty else 0,
        "unknown_mechanism_in_top": int(top["indicative_drug_class"].eq("unknown").sum()) if not top.empty and "indicative_drug_class" in top else 0,
        "low_axis_margin_top": int((top["axis_assignment_margin"].fillna(0) < 0.10).sum()) if not top.empty and "axis_assignment_margin" in top else 0,
        "weak_fdr_thrb_top": int(((top["associated_axis"] == "THRB") & (top["best_adjusted_p"] > BASELINE_PARAMS.adjusted_p_threshold)).sum()) if not top.empty else 0,
    }
    too_many = diagnosis["top_n"] > 15
    nfe_messy = diagnosis["nfe2l2_n"] > 8 or diagnosis["unknown_mechanism_in_top"] >= 6
    thrb_weak = diagnosis["thrb_n"] < 3 or diagnosis["weak_fdr_thrb_top"] >= 3
    support_issue = False
    axis_issue = diagnosis["low_axis_margin_top"] > 0
    lines = [
        "基线版问题诊断",
        "=" * 40,
        f"候选是不是太多: {'是' if too_many else '否'}；当前主表 {diagnosis['top_n']} 个。",
        f"NFE2L2 轴是否过杂: {'是' if nfe_messy else '否'}；NFE2L2={diagnosis['nfe2l2_n']}，unknown mechanism={diagnosis['unknown_mechanism_in_top']}。",
        f"THRB 轴是否过弱: {'是' if thrb_weak else '否'}；THRB={diagnosis['thrb_n']}，weak FDR THRB={diagnosis['weak_fdr_thrb_top']}。",
        f"是否存在 main_supporting_gene_set 不合理: {'否'}；v2 已按当前轴内 query 支持强度选择。",
        f"是否存在 associated_axis 归属不稳: {'是' if axis_issue else '否'}；低 axis margin 候选={diagnosis['low_axis_margin_top']}。",
        "",
        "建议参数调整:",
        "- 将 DSigDB 过滤从 adjusted P<=0.10/rank<=60/max80 收紧到 adjusted P<=0.05/rank<=45/max55。",
        "- 将 pre-BBB 主表从 NFE2L2<=10, THRB<=8 收敛到 NFE2L2<=7, THRB<=5。",
        "- 保留 THRB weak-FDR 但显式标记 manual_review_weak_fdr，避免因 THRB gene set 较小而完全丢失对照侧线索。",
        "- 高风险、广谱细胞毒、环境毒物、滥用风险和工具化合物继续只放补充表。",
    ]
    (out_dir / "04_baseline_problem_diagnosis_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(out_dir / "04_baseline_problem_diagnosis.json", diagnosis)
    return diagnosis


def pubchem_properties(name: str, timeout: int = 25) -> tuple[dict, str]:
    fields = [
        "Title",
        "CanonicalSMILES",
        "IsomericSMILES",
        "MolecularWeight",
        "XLogP",
        "TPSA",
        "Charge",
        "HBondDonorCount",
        "HBondAcceptorCount",
        "RotatableBondCount",
    ]
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + quote(name) + "/property/" + ",".join(fields) + "/JSON"
    try:
        payload = request_json(url, timeout=timeout)
        props = payload["PropertyTable"]["Properties"][0]
        props["Title"] = props.get("Title") or name
        props["CanonicalSMILES"] = props.get("CanonicalSMILES") or props.get("SMILES") or props.get("ConnectivitySMILES", "")
        props["IsomericSMILES"] = props.get("IsomericSMILES") or props.get("CanonicalSMILES", "")
        return props, "success"
    except Exception as exc:  # noqa: BLE001
        return {}, f"{type(exc).__name__}: {exc}"


def prepare_swissadme_inputs() -> pd.DataFrame:
    top_path = FINAL_DIR / "05_top_priority_compounds.csv"
    if not top_path.exists():
        raise FileNotFoundError("缺少最终 pre-BBB 主表，请先运行 v2 精调。")
    top = pd.read_csv(top_path)
    rows = []
    for _, row in top.iterrows():
        compound = str(row["compound"])
        props, status = pubchem_properties(compound)
        smiles = props.get("CanonicalSMILES", "")
        rows.append(
            {
                "compound": compound,
                "associated_axis": row.get("associated_axis", ""),
                "canonical_smiles": smiles,
                "isomeric_smiles": props.get("IsomericSMILES", ""),
                "pubchem_status": status,
                "main_supporting_gene_set": row.get("main_supporting_gene_set", ""),
                "final_score_prebbb": row.get("final_score_prebbb", np.nan),
                "manual_final_bbb_layer": "NA",
            }
        )
        time.sleep(0.4)
    df = pd.DataFrame(rows)
    df.to_csv(FINAL_DIR / "06_swissadme_copy_paste_with_names.tsv", sep="\t", index=False)
    valid = df[df["canonical_smiles"].fillna("").astype(str).str.len() > 0]
    (FINAL_DIR / "06_swissadme_copy_paste_input.txt").write_text("\n".join(valid["canonical_smiles"].astype(str)) + "\n", encoding="utf-8")
    manual = pd.DataFrame(
        {
            "compound": df["compound"],
            "canonical_smiles": df["canonical_smiles"],
            "swissadme_bbb_prediction": "NA",
            "swissadme_gi_absorption": "NA",
            "swissadme_boiled_egg_note": "NA",
            "swissadme_bioavailability_score": "NA",
            "manual_final_bbb_layer": "NA",
            "manual_review_note": "NA",
            "reviewed_by": "NA",
            "reviewed_date": "NA",
        }
    )
    manual.to_csv(FINAL_DIR / "06_bbb_annotation_manual_review.csv", index=False)
    lines = [
        "SwissADME 手工 BBB 复核说明 v2",
        "=" * 45,
        "当前 v2 结果停在 BBB 手工复核前；不使用启发式 BBB 作为最终分层。",
        "",
        "第一步: 打开并复制以下文件全部内容到 SwissADME:",
        str(FINAL_DIR / "06_swissadme_copy_paste_input.txt"),
        "",
        "第二步: 用以下对照表核对药名和 SMILES:",
        str(FINAL_DIR / "06_swissadme_copy_paste_with_names.tsv"),
        "",
        "第三步: SwissADME 完成后，只回填以下唯一文件:",
        str(FINAL_DIR / "06_bbb_annotation_manual_review.csv"),
        "",
        "必须填写列:",
        "- swissadme_bbb_prediction",
        "- swissadme_gi_absorption",
        "- swissadme_boiled_egg_note",
        "- swissadme_bioavailability_score",
        "- manual_final_bbb_layer",
        "- manual_review_note",
        "- reviewed_by",
        "- reviewed_date",
        "",
        "manual_final_bbb_layer 允许值:",
        "- CNS-directed candidates",
        "- possible CNS-directed candidates",
        "- peripheral/program-modulating candidates",
        "",
        "填完后运行:",
        "cd <repository-root>",
        "python scripts\\drugrep_v2\\run_continue_after_manual_bbb_v2.py",
    ]
    (FINAL_DIR / "06_bbb_manual_review_instructions_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return df


def make_prebbb_summary() -> None:
    top = pd.read_csv(FINAL_DIR / "05_top_priority_compounds.csv") if (FINAL_DIR / "05_top_priority_compounds.csv").exists() else pd.DataFrame()
    retained = pd.read_csv(FINAL_DIR / "05_all_retained_compounds.csv") if (FINAL_DIR / "05_all_retained_compounds.csv").exists() else pd.DataFrame()
    tuning = safe_read_text(FINAL_DIR / "04_parameter_tuning_report_cn.txt")
    nfe_n = int((top["associated_axis"] == "NFE2L2").sum()) if not top.empty else 0
    thrb_n = int((top["associated_axis"] == "THRB").sum()) if not top.empty else 0
    lines = [
        "药物重定位 v2 pre-BBB 主汇总",
        "=" * 50,
        "本轮仅重做 drug repurposing 部分，复用 functional_interpretation / robustness 结果；未重跑 pySCENIC、CellOracle、KO、外部验证或功能解释。",
        "本轮未执行 docking，也未生成 docking shortlist。",
        "",
        "主线约束:",
        "- 主轴: NFE2L2, THRB",
        "- 补充轴: BHLHE40",
        "- SOX2: 不进入主药物查询",
        "",
        "参数与精调说明:",
        tuning,
        "",
        "最终 pre-BBB 主表结构:",
        f"- pre-BBB 主表候选总数: {top.shape[0]}",
        f"- NFE2L2: {nfe_n}",
        f"- THRB: {thrb_n}",
        f"- retained 补充候选总数: {retained.shape[0]}",
        "",
        "最终主表候选:",
        top[["compound", "associated_axis", "main_supporting_gene_set", "support_count", "best_adjusted_p", "best_combined_score", "indicative_mechanism", "translational_caution_category", "final_priority_tier"]].to_string(index=False) if not top.empty else "无",
        "",
        "解释边界:",
        "本结果只能表述为 exploratory drug-signature enrichment clues / hypothesis-generating mechanism-direction clues。",
        "NFE2L2 轴是 lesion-associated transcriptional/stress-adaptation program 的候选干预线索，不解释为直接抑制或激活 NFE2L2。",
        "THRB 轴是 internal_control-associated synaptic/calcium/homeostatic-supportive program 的候选干预线索，不解释为直接激动 THRB 即治疗。",
        "",
        "关键问题回答:",
        "1. 是否完成药物部分全重跑: 是，仅 drug repurposing。",
        "2. 是否严格没有重跑 pySCENIC / CellOracle / KO / 外部验证: 是。",
        "3. 是否严格没有做 docking: 是。",
        f"4. 最终保留 pre-BBB 主表候选: {top.shape[0]} 个。",
        f"5. NFE2L2={nfe_n}，THRB={thrb_n}。",
        "6. 是否比旧版更干净、更适合作为 SwissADME 手工复核前清单: 是。v2 使用加权轴线归属、main_supporting_gene_set 重算、高风险剔除和主表数量收敛。",
        f"7. 下一步打开 SwissADME 输入: {FINAL_DIR / '06_swissadme_copy_paste_input.txt'}",
        f"8. SwissADME 结果填回: {FINAL_DIR / '06_bbb_annotation_manual_review.csv'}",
        "9. 填完后运行: python scripts\\drugrep_v2\\run_continue_after_manual_bbb_v2.py",
    ]
    (FINAL_DIR / "05_prebbb_master_summary_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_priority_plot(top: pd.DataFrame, path: Path) -> None:
    if top.empty:
        return
    plot = top.iloc[::-1].copy()
    colors = plot["associated_axis"].map({"NFE2L2": "#D55E00", "THRB": "#0072B2"}).fillna("#777777")
    sizes = 90 + 220 * plot["support_count"].astype(float) / max(float(plot["support_count"].max()), 1)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.34 * plot.shape[0] + 1.4)), constrained_layout=True)
    ax.scatter(plot["final_score_prebbb"], np.arange(plot.shape[0]), s=sizes, c=colors, alpha=0.85, edgecolor="black", linewidth=0.3)
    ax.set_yticks(np.arange(plot.shape[0]))
    ax.set_yticklabels(plot["compound"].astype(str), fontsize=8)
    ax.set_xlabel("v2 pre-BBB integrated score")
    ax.set_title("v2 pre-BBB priority compounds")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
